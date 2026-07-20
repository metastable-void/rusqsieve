import { MESSAGE, copyBuffer, putBytes } from "./protocol.js";
let instance;
self.onmessage = async ({ data }) => {
  if (data.type === MESSAGE.READY) { instance = await WebAssembly.instantiate(data.module, {}); self.postMessage({ type: MESSAGE.READY }); return; }
  if (data.type !== MESSAGE.EXECUTE || !instance) return;
  try { const pointer = putBytes(instance, data.packet); const handle = instance.exports.qs_worker_execute(data.context, pointer, data.packet.length); instance.exports.qs_dealloc(pointer, data.packet.length, 1); self.postMessage({ type: MESSAGE.RESULT, packet: copyBuffer(instance, handle) }); }
  catch (error) { self.postMessage({ type: MESSAGE.ERROR, error: String(error?.message || error) }); }
};
