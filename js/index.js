import { MESSAGE } from "./protocol.js";

/** Factor a decimal string or unsigned byte array using a coordinator Web Worker. */
export function factor(input, options = {}) {
  const decimal = typeof input === "string" ? input : bytesToDecimal(input);
  if (!/^\d+$/.test(decimal)) return Promise.reject(new TypeError("input must be an unsigned decimal string or byte array"));
  const coordinator = new Worker(new URL("./coordinator.js", import.meta.url), { type: "module" });
  return new Promise((resolve, reject) => {
    const abort = () => coordinator.postMessage({ type: MESSAGE.CANCEL });
    options.signal?.addEventListener("abort", abort, { once: true });
    coordinator.onmessage = ({ data }) => {
      if (data.type === MESSAGE.PROGRESS) options.onProgress?.(data.progress);
      if (data.type === MESSAGE.DONE) { cleanup(); resolve(data.value); }
      if (data.type === MESSAGE.ERROR) { cleanup(); reject(new Error(data.error)); }
    };
    coordinator.onerror = (event) => { cleanup(); reject(event.error || new Error(event.message)); };
    coordinator.postMessage({ type: MESSAGE.START, decimal, options: { parallelism: options.parallelism, wasmUrl: options.wasmUrl } });
    function cleanup() { options.signal?.removeEventListener("abort", abort); coordinator.terminate(); }
  });
}

function bytesToDecimal(value) {
  const bytes = value instanceof Uint8Array ? value : new Uint8Array(value);
  let n = 0n; for (const byte of bytes) n = (n << 8n) | BigInt(byte); return n.toString();
}
