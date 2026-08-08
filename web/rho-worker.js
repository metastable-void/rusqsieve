// Deep Pollard-Brent worker.
//
// The main thread peels cheap factors inline with BigInt, which is right for an opportunistic peel
// and hopeless for a real search: a composite the sieve cannot help with — one above
// `qs_max_siqs_bits`, or a cofactor an earlier split already proved unbalanced — needs tens of
// millions of iterations. Measured under Node, BigInt runs that loop at 288k iterations/s on a
// 512-bit modulus against 1.08M/s for the same algorithm in wasm over Montgomery-encoded limbs, and
// on the main thread every one of those iterations is a frame not painted. So the search runs here,
// in wasm, off the thread that owns the page.
//
// Each worker walks its own polynomial constants, so a pool runs that many independent walks over
// one modulus and the first collision wins. Independent walks collide in about
// `1.2·sqrt(p)/sqrt(T)` iterations for `T` of them — a `sqrt(T)` speedup rather than a linear one,
// which is still the difference between reaching a 2^49 factor and a 2^52 one at equal wall clock.
//
// There is no cancellation protocol because there is no need for one: `qs_rho` runs to its budget
// and returns, and the main thread cancels by terminating the worker.
import { bytesToBigInt, instantiate, putString, takePacket, validateDecimalInput } from "./abi.js";

// Matches `MAX_INPUT_BITS` in index.js: what `Natural` can hold. The default in `abi.js` is the
// sieve's 400-bit range, which is exactly the bound this worker exists to work above.
const MAX_MODULUS_BITS = 1024;

let ex = null;

self.onmessage = async ({ data }) => {
  const gen = Number.isSafeInteger(data?.gen) ? data.gen : 0;
  try {
    if (!data || typeof data !== "object") throw new Error("invalid rho worker command");
    if (data.cmd === "init") {
      if (ex) throw new Error("rho worker is already initialized");
      ex = await instantiate(data.module);
      const abi = ex.qs_abi_version();
      if (abi !== 4) throw new Error(`unsupported rusqsieve wasm ABI ${abi}`);
      if (typeof ex.qs_rho !== "function") throw new Error("wasm module exports no qs_rho");
      self.postMessage({ type: "ready", gen, abi });
      return;
    }
    if (data.cmd === "search") {
      if (!ex) throw new Error("rho worker is not initialized");
      if (!Number.isSafeInteger(data.budget) || data.budget <= 0) {
        throw new Error("rho search requires a positive iteration budget");
      }
      if (
        !Number.isSafeInteger(data.first) ||
        data.first <= 0 ||
        !Number.isSafeInteger(data.count) ||
        data.count <= 0
      ) {
        throw new Error("rho search requires a positive constant range");
      }
      const decimal = validateDecimalInput(data.n, MAX_MODULUS_BITS);
      const input = putString(ex, decimal);
      let handle;
      try {
        handle = ex.qs_rho(input.ptr, input.len, data.budget, data.first, data.count);
      } finally {
        ex.qs_dealloc(input.ptr, input.len, 1);
      }
      if (!handle) throw new Error("wasm rejected the rho modulus");
      const payload = takePacket(ex, handle, 12);
      // An empty payload is the budget being spent without a split, which is a result, not a
      // failure: the caller decides whether the composite now goes to the sieve or is refused.
      const factor = payload.length === 0 ? null : bytesToBigInt(payload);
      self.postMessage({
        type: "done",
        gen,
        factor: factor === null ? null : factor.toString(),
      });
      return;
    }
    throw new Error("invalid rho worker command");
  } catch (error) {
    self.postMessage({ type: "error", gen, message: String(error?.message ?? error) });
  }
};
