export const MESSAGE = Object.freeze({ START: "start", READY: "ready", EXECUTE: "execute", RESULT: "result", PROGRESS: "progress", DONE: "done", ERROR: "error", CANCEL: "cancel" });

export function copyBuffer(instance, handle) {
  const pointer = instance.exports.qs_buffer_pointer(handle);
  const length = instance.exports.qs_buffer_length(handle);
  const bytes = new Uint8Array(instance.exports.memory.buffer, pointer, length).slice();
  instance.exports.qs_buffer_free(handle);
  if (bytes.length < 12 || new TextDecoder().decode(bytes.subarray(0, 4)) !== "QSV1") throw new Error("invalid quadratic-sieve packet");
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const payloadLength = view.getUint32(8, true);
  if (payloadLength !== bytes.length - 12) throw new Error("invalid quadratic-sieve packet length");
  return bytes.subarray(12);
}

export function putBytes(instance, bytes) {
  const pointer = instance.exports.qs_alloc(bytes.length, 1);
  if (!pointer && bytes.length) throw new Error("WebAssembly allocation failed");
  new Uint8Array(instance.exports.memory.buffer, pointer, bytes.length).set(bytes);
  return pointer;
}
