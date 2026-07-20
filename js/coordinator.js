import { MESSAGE, copyBuffer, putBytes } from "./protocol.js";
let cancelled = false;
self.onmessage = async ({ data }) => {
  if (data.type === MESSAGE.CANCEL) { cancelled = true; return; }
  if (data.type !== MESSAGE.START) return;
  try {
    const wasmUrl = data.options.wasmUrl || new URL("../target/wasm32-unknown-unknown/release/rusqsieve.wasm", import.meta.url);
    const module = await WebAssembly.compileStreaming(fetch(wasmUrl));
    const { exports } = await WebAssembly.instantiate(module, {});
    const hint = Math.max(1, Math.min(data.options.parallelism ?? Number.MAX_SAFE_INTEGER, navigator.hardwareConcurrency || 1));
    const workers = Array.from({ length: hint }, () => { const w = new Worker(new URL("./worker.js", import.meta.url), { type: "module" }); w.postMessage({ type: MESSAGE.READY, module }); return w; });
    const bytes = new TextEncoder().encode(data.decimal); const pointer = putBytes({ exports }, bytes);
    const session = exports.qs_session_new(pointer, bytes.length, 0, 0); exports.qs_dealloc(pointer, bytes.length, 1);
    if (!session) throw new Error("invalid input or session configuration");
    while (!cancelled) {
      const status = exports.qs_session_advance_local(session, 0, 0);
      const progressHandle = exports.qs_session_progress(session);
      if (progressHandle) { const p = copyBuffer({ exports }, progressHandle); const v = new DataView(p.buffer, p.byteOffset, p.byteLength); self.postMessage({ type: MESSAGE.PROGRESS, progress: { revision: Number(v.getBigUint64(0, true)), phase: v.getUint32(8, true), completed: Number(v.getBigUint64(12, true)) } }); }
      if (status < 0) throw new Error("factorization failed"); if (status > 0) break;
      await new Promise(resolve => setTimeout(resolve, 0));
    }
    if (cancelled) throw new DOMException("Factorization aborted", "AbortError");
    const payload = copyBuffer({ exports }, exports.qs_session_take_factors(session));
    const grouped = new TextDecoder().decode(payload).trim().split("\n").filter(Boolean).map(line => { const [prime, exponent] = line.split(":"); return { prime, exponent: Number(exponent) }; });
    const factors = grouped.flatMap(({ prime, exponent }) => Array(exponent).fill(prime));
    workers.forEach(worker => worker.terminate()); self.postMessage({ type: MESSAGE.DONE, value: { factors, grouped } });
  } catch (error) { self.postMessage({ type: MESSAGE.ERROR, error: String(error?.message || error) }); }
};
