// Low-level helpers for talking to the rusqsieve wasm C ABI.
// Views are recreated on every call because wasm memory may have grown.

export async function loadModule(url) {
  if (WebAssembly.compileStreaming) {
    try {
      return await WebAssembly.compileStreaming(fetchChecked(url));
    } catch {
      /* fall through for file:// or servers without wasm MIME type */
    }
  }
  const bytes = await (await fetchChecked(url)).arrayBuffer();
  return WebAssembly.compile(bytes);
}

export async function instantiate(module) {
  if (!(module instanceof WebAssembly.Module)) {
    throw new TypeError("expected a compiled WebAssembly.Module");
  }
  const inst = await WebAssembly.instantiate(module, {});
  if (!inst?.exports?.memory) throw new Error("wasm module does not export memory");
  return inst.exports;
}

export function putBytes(ex, bytes) {
  if (!ArrayBuffer.isView(bytes)) throw new TypeError("expected an ArrayBuffer view");
  if (bytes.byteLength > 0xffff_ffff) throw new RangeError("input exceeds wasm32 limits");
  const length = bytes.byteLength;
  const ptr = ex.qs_alloc(length, 1);
  if (!ptr && length) throw new Error("wasm allocation failed");
  try {
    new Uint8Array(ex.memory.buffer, ptr, length).set(
      new Uint8Array(bytes.buffer, bytes.byteOffset, bytes.byteLength),
    );
  } catch (error) {
    if (ptr) ex.qs_dealloc(ptr, length, 1);
    throw error;
  }
  return { ptr, len: length };
}

export function putString(ex, str) {
  return putBytes(ex, new TextEncoder().encode(str));
}

// The default is the engine's SIQS range (`engine::MAX_SIQS_BITS`), not a general input limit:
// both endpoints that call this — the sieve worker and the coordinator — are handed a composite
// that the main thread has already peeled, so anything arriving here is destined for the sieve.
// Whole-input width is bounded separately, by what `Natural` can hold.
export function validateDecimalInput(text, maximumBits = 400) {
  if (typeof text !== "string" || !/^\d+$/u.test(text)) {
    throw new Error("input must be an unsigned decimal integer");
  }
  const significant = text.replace(/^0+/u, "") || "0";
  // Use the generic digit bound only as an allocation-avoidance precheck; the BigInt bit-length
  // check below is authoritative.
  const maximumDigits = Math.ceil(maximumBits * Math.LOG10E * Math.LN2);
  if (significant.length > maximumDigits) {
    throw new Error(`input exceeds the ${maximumBits}-bit limit`);
  }
  const value = BigInt(significant);
  if (value <= 0n) throw new Error("input must be positive");
  if (value.toString(2).length > maximumBits) {
    throw new Error(`input exceeds the ${maximumBits}-bit limit`);
  }
  return significant;
}

// Read a QSV1 packet buffer handle, returning a *copy* of its payload and
// freeing it. The envelope is an ABI boundary even when both endpoints ship in
// this repository, so reject truncation, stale handles, and wrong packet kinds.
export function takePacket(ex, handle, expectedKind = null) {
  if (!handle) return null;
  try {
    const ptr = ex.qs_buffer_pointer(handle);
    const len = ex.qs_buffer_length(handle);
    const memoryLength = ex.memory.buffer.byteLength;
    if (
      !Number.isInteger(ptr) ||
      !Number.isInteger(len) ||
      ptr <= 0 ||
      len < 12 ||
      ptr > memoryLength ||
      len > memoryLength - ptr
    ) {
      throw new Error("invalid or stale wasm packet handle");
    }
    const raw = new Uint8Array(ex.memory.buffer, ptr, len);
    const view = new DataView(ex.memory.buffer, ptr, len);
    if (raw[0] !== 0x51 || raw[1] !== 0x53 || raw[2] !== 0x56 || raw[3] !== 0x31) {
      throw new Error("invalid wasm packet magic");
    }
    const kind = view.getUint16(4, true);
    const version = view.getUint16(6, true);
    const payloadLen = view.getUint32(8, true);
    if (version !== 1) throw new Error(`unsupported wasm packet version ${version}`);
    if (expectedKind !== null && kind !== expectedKind) {
      throw new Error(`unexpected wasm packet kind ${kind}; expected ${expectedKind}`);
    }
    if (payloadLen !== len - 12) throw new Error("invalid wasm packet payload length");
    return raw.slice(12);
  } finally {
    ex.qs_buffer_free(handle);
  }
}

// 128-byte little-endian Natural<16> payload -> BigInt.
export function bytesToBigInt(bytes) {
  if (!(bytes instanceof Uint8Array)) throw new TypeError("expected a Uint8Array");
  let n = 0n;
  for (let i = bytes.length - 1; i >= 0; i--) n = (n << 8n) | BigInt(bytes[i]);
  return n;
}

async function fetchChecked(url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`failed to fetch ${response.url || url}: HTTP ${response.status}`);
  }
  return response;
}
