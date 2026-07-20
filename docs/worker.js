// Sieve worker: an independent wasm instance that rebuilds the deterministic sieve
// context and sieves the polynomial-family ranges the coordinator assigns to it.
import { instantiate, putString, takePacket } from "./abi.js";

let ex = null;
let context = 0;

self.onmessage = async ({ data }) => {
  try {
    if (data.cmd === "init") {
      ex = await instantiate(data.module);
      self.postMessage({ type: "ready" });
    } else if (data.cmd === "prepare") {
      if (context) {
        ex.qs_worker_free(context);
        context = 0;
      }
      const n = putString(ex, data.n);
      context = ex.qs_worker_prepare(n.ptr, n.len);
      ex.qs_dealloc(n.ptr, n.len, 1);
      self.postMessage({ type: "prepared", ok: context !== 0 });
    } else if (data.cmd === "sieve") {
      const handle = ex.qs_worker_sieve(context, data.family, data.count);
      const payload = takePacket(ex, handle); // raw [count][len,bytes]…
      if (payload) self.postMessage({ type: "relations", payload }, [payload.buffer]);
      else self.postMessage({ type: "relations", payload: null });
    }
  } catch (error) {
    self.postMessage({ type: "error", error: String(error?.message || error) });
  }
};
