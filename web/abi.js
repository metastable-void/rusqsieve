// Low-level helpers for talking to the rusqsieve wasm C ABI.
// Views are recreated on every call because wasm memory may have grown.

export async function loadModule(url) {
  if (WebAssembly.compileStreaming) {
    try {
      return await WebAssembly.compileStreaming(fetch(url));
    } catch {
      /* fall through for file:// or servers without wasm MIME type */
    }
  }
  const bytes = await (await fetch(url)).arrayBuffer();
  return WebAssembly.compile(bytes);
}

export async function instantiate(module) {
  const inst = await WebAssembly.instantiate(module, {});
  return inst.exports;
}

export function putBytes(ex, bytes) {
  const ptr = ex.qs_alloc(bytes.length, 1);
  if (!ptr && bytes.length) throw new Error("wasm allocation failed");
  new Uint8Array(ex.memory.buffer, ptr, bytes.length).set(bytes);
  return { ptr, len: bytes.length };
}

export function putString(ex, str) {
  return putBytes(ex, new TextEncoder().encode(str));
}

// Read a QSV1 packet buffer handle, returning a *copy* of its payload and freeing it.
export function takePacket(ex, handle) {
  if (!handle) return null;
  const ptr = ex.qs_buffer_pointer(handle);
  const len = ex.qs_buffer_length(handle);
  const raw = new Uint8Array(ex.memory.buffer, ptr, len);
  const view = new DataView(ex.memory.buffer, ptr, len);
  if (len < 12) {
    ex.qs_buffer_free(handle);
    return null;
  }
  const payloadLen = view.getUint32(8, true);
  const payload = raw.slice(12, 12 + payloadLen);
  ex.qs_buffer_free(handle);
  return payload;
}

// 128-byte little-endian Natural<16> payload -> BigInt.
export function bytesToBigInt(bytes) {
  let n = 0n;
  for (let i = bytes.length - 1; i >= 0; i--) n = (n << 8n) | BigInt(bytes[i]);
  return n;
}
