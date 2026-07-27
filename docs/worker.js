// Sieve worker: an independent wasm instance that rebuilds the deterministic sieve
// context and sieves the polynomial-family ranges the coordinator assigns to it.
import { instantiate, putString, takePacket, validateDecimalInput } from "./abi.js";

let ex = null;
let context = 0;

self.onmessage = async ({ data }) => {
  const gen = Number.isSafeInteger(data?.gen) ? data.gen : 0;
  try {
    if (!data || typeof data !== "object") throw new Error("invalid worker command");
    if (data.cmd === "init") {
      if (ex) throw new Error("worker is already initialized");
      ex = await instantiate(data.module);
      self.postMessage({ type: "ready" });
    } else if (data.cmd === "prepare") {
      if (!ex) throw new Error("worker is not initialized");
      if (!Number.isSafeInteger(data.gen) || data.gen <= 0 || typeof data.n !== "string") {
        throw new Error("invalid prepare command");
      }
      if (context) {
        ex.qs_worker_free(context);
        context = 0;
      }
      const decimal = validateDecimalInput(data.n);
      const n = putString(ex, decimal);
      try {
        context = ex.qs_worker_prepare(n.ptr, n.len);
      } finally {
        ex.qs_dealloc(n.ptr, n.len, 1);
      }
      self.postMessage({ type: "prepared", ok: context !== 0, gen: data.gen });
    } else if (data.cmd === "sieve") {
      if (!ex || !context) throw new Error("worker has no prepared sieve context");
      if (
        !Number.isSafeInteger(data.gen) ||
        data.gen <= 0 ||
        !Number.isInteger(data.family) ||
        data.family < 0 ||
        data.family > 0xffff_ffff ||
        !Number.isInteger(data.count) ||
        data.count < 1 ||
        data.count > 4096 ||
        data.family + data.count > 0x1_0000_0000
      ) {
        throw new Error("invalid sieve command");
      }
      const handle = ex.qs_worker_sieve(context, data.family, data.count);
      const payload = takePacket(ex, handle, 10); // raw [count][len,bytes]…
      if (payload) self.postMessage({ type: "relations", payload, gen: data.gen }, [payload.buffer]);
      else self.postMessage({ type: "relations", payload: null, gen: data.gen });
    } else {
      throw new Error(`unknown worker command: ${String(data.cmd)}`);
    }
  } catch (error) {
    self.postMessage({
      type: "error",
      error: String(error?.message || error),
      gen,
      phase: data?.cmd === "init" ? "boot" : "run",
    });
  }
};
