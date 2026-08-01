// Coordinator worker: owns relation collection and the serial GF(2) solve so
// neither Wasm linear algebra nor extraction blocks the browser main thread.
import {
  instantiate,
  putString,
  putBytes,
  takePacket,
  validateDecimalInput,
} from "./abi.js";

let ex = null;
let session = 0;
let generation = 0;
let baseTarget = 0;
let extractionTarget = 0;
let lastRelations = 0;

self.onmessage = async ({ data }) => {
  const messageGeneration = Number.isSafeInteger(data?.gen) ? data.gen : 0;
  try {
    if (!data || typeof data !== "object") throw new Error("invalid coordinator command");
    if (data.cmd === "init") {
      if (ex) throw new Error("coordinator is already initialized");
      ex = await instantiate(data.module);
      const abi = ex.qs_abi_version();
      if (abi !== 3) throw new Error(`unsupported rusqsieve wasm ABI ${abi}`);
      // Report the engine's sieve range rather than duplicating it in the UI. The main thread
      // uses it to reject an over-wide composite with a specific message before spinning up a
      // run that `qs_coord_new` would refuse anyway.
      self.postMessage({ type: "ready", abi, maxSiqsBits: ex.qs_max_siqs_bits() });
      return;
    }
    if (data.cmd === "new") {
      if (!ex) throw new Error("coordinator is not initialized");
      if (!Number.isSafeInteger(data.gen) || data.gen <= 0 || typeof data.n !== "string") {
        throw new Error("invalid session command");
      }
      if (session) ex.qs_coord_free(session);
      generation = data.gen;
      session = 0;
      baseTarget = 0;
      extractionTarget = 0;
      lastRelations = 0;
      const decimal = validateDecimalInput(data.n);
      const input = putString(ex, decimal);
      try {
        session = ex.qs_coord_new(input.ptr, input.len);
      } finally {
        ex.qs_dealloc(input.ptr, input.len, 1);
      }
      if (!session) throw new Error("could not build a sieve for this number");
      baseTarget = ex.qs_coord_target(session);
      if (!baseTarget) throw new Error("coordinator returned an invalid relation target");
      extractionTarget = baseTarget;
      self.postMessage({
        type: "session",
        gen: generation,
        target: extractionTarget,
        // The family budget scales with input width, so the scheduler on the main thread has to
        // read it from the session it is actually driving instead of assuming a constant.
        familyBudget: ex.qs_coord_family_budget(session),
      });
      return;
    }
    if (data.cmd === "submit") {
      if (data.gen !== generation) return;
      if (!session) throw new Error("coordinator has no active session");
      if (
        !(data.payload instanceof Uint8Array) ||
        !Number.isInteger(data.worker) ||
        data.worker < 0
      ) {
        throw new Error("invalid relation submission");
      }
      validateRelationBatch(data.payload);
      const bytes = putBytes(ex, data.payload);
      let relations;
      try {
        relations = ex.qs_coord_submit(session, bytes.ptr, bytes.len);
      } finally {
        ex.qs_dealloc(bytes.ptr, bytes.len, 1);
      }
      if (relations < lastRelations) throw new Error("coordinator relation count regressed");
      lastRelations = relations;
      if (relations < extractionTarget) {
        self.postMessage({
          type: "submitted",
          gen: generation,
          relations,
          target: extractionTarget,
          worker: data.worker,
        });
        return;
      }
      // Notify first: the main thread can paint the phase change while this
      // worker proceeds into the serial solve.
      self.postMessage({
        type: "linalg",
        gen: generation,
        relations,
        target: extractionTarget,
      });
      const handle = ex.qs_coord_extract(session);
      const factor = takePacket(ex, handle, 11);
      if (factor) {
        if (factor.byteLength !== 128) throw new Error("coordinator returned a malformed factor");
        ex.qs_coord_free(session);
        session = 0;
        self.postMessage({ type: "factor", gen: generation, factor }, [factor.buffer]);
        return;
      }

      // A relation target is a heuristic. A valid matrix can still yield only
      // trivial GCDs in the bounded dependency set, so preserve the session and
      // request a measured surplus instead of destroying all collected work.
      const additional = Math.max(64, Math.ceil(baseTarget / 20));
      extractionTarget = relations + additional;
      self.postMessage({
        type: "submitted",
        gen: generation,
        relations,
        target: extractionTarget,
        worker: data.worker,
        retry: true,
      });
      return;
    }
    throw new Error(`unknown coordinator command: ${String(data.cmd)}`);
  } catch (error) {
    self.postMessage({
      type: "error",
      gen: messageGeneration,
      phase: data?.cmd === "init" ? "boot" : "run",
      error: String(error?.message || error),
    });
  }
};

function validateRelationBatch(payload) {
  if (payload.byteLength < 4) throw new Error("truncated relation batch");
  const view = new DataView(payload.buffer, payload.byteOffset, payload.byteLength);
  const count = view.getUint32(0, true);
  if (count > 4096) throw new Error("relation batch count exceeds the worker limit");
  let offset = 4;
  for (let index = 0; index < count; index++) {
    if (offset > payload.byteLength - 4) throw new Error("truncated relation length");
    const length = view.getUint32(offset, true);
    offset += 4;
    // A serialized family starts with family:u64, polynomials:u64, count:u32.
    if (length < 20 || length > payload.byteLength - offset) {
      throw new Error("invalid serialized family length");
    }
    offset += length;
  }
  if (offset !== payload.byteLength) throw new Error("trailing bytes in relation batch");
}
