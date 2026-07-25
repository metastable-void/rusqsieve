// Coordinator worker: owns relation collection and the serial GF(2) solve so
// neither Wasm linear algebra nor extraction blocks the browser main thread.
import { instantiate, putString, putBytes, takePacket } from "./abi.js";

let ex = null;
let session = 0;
let generation = 0;

self.onmessage = async ({ data }) => {
  try {
    if (data.cmd === "init") {
      ex = await instantiate(data.module);
      self.postMessage({ type: "ready", abi: ex.qs_abi_version() });
      return;
    }
    if (data.cmd === "new") {
      if (session) ex.qs_coord_free(session);
      generation = data.gen;
      const input = putString(ex, data.n);
      session = ex.qs_coord_new(input.ptr, input.len);
      ex.qs_dealloc(input.ptr, input.len, 1);
      if (!session) throw new Error("could not build a sieve for this number");
      self.postMessage({
        type: "session",
        gen: generation,
        target: ex.qs_coord_target(session),
      });
      return;
    }
    if (data.cmd === "submit" && data.gen === generation && session) {
      const bytes = putBytes(ex, data.payload);
      const relations = ex.qs_coord_submit(session, bytes.ptr, bytes.len);
      ex.qs_dealloc(bytes.ptr, bytes.len, 1);
      const target = ex.qs_coord_target(session);
      if (relations < target) {
        self.postMessage({
          type: "submitted",
          gen: generation,
          relations,
          target,
          worker: data.worker,
        });
        return;
      }
      // Notify first: the main thread can paint the phase change while this
      // worker proceeds into the serial solve.
      self.postMessage({ type: "linalg", gen: generation, relations, target });
      const handle = ex.qs_coord_extract(session);
      const factor = takePacket(ex, handle);
      ex.qs_coord_free(session);
      session = 0;
      if (!factor) throw new Error("linear algebra found no factor");
      self.postMessage({ type: "factor", gen: generation, factor }, [factor.buffer]);
    }
  } catch (error) {
    self.postMessage({
      type: "error",
      gen: data.gen,
      error: String(error?.message || error),
    });
  }
};
