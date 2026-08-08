// Elliptic curve method worker.
//
// This is the last stage the frontend has for a composite the sieve refuses. Trial division and
// the primality test have run, Pollard-Brent has spent a deep budget without splitting it, and the
// coordinator will not accept it — so the alternative to running curves is telling the user the
// number is too large. ECM's cost depends on the size of the *factor* rather than of the input, so
// a 500-bit number with a 25-digit factor is ordinary work here and impossible for either of the
// other two stages.
//
// Each worker runs its own stretch of the σ sequence, so a pool covers that many curves at once and
// the first success wins. There is no cancellation protocol: `qs_ecm` runs to its curve count and
// returns, and the main thread cancels by terminating the worker.
import { bytesToBigInt, instantiate, putString, takePacket, validateDecimalInput } from "./abi.js";

// Matches `MAX_INPUT_BITS` in index.js: what `Natural` can hold. The default in `abi.js` is the
// sieve's 400-bit range, which is the bound this worker exists to work above.
const MAX_MODULUS_BITS = 1024;

let ex = null;

self.onmessage = async ({ data }) => {
  const gen = Number.isSafeInteger(data?.gen) ? data.gen : 0;
  try {
    if (!data || typeof data !== "object") throw new Error("invalid ecm worker command");
    if (data.cmd === "init") {
      if (ex) throw new Error("ecm worker is already initialized");
      ex = await instantiate(data.module);
      const abi = ex.qs_abi_version();
      if (abi !== 5) throw new Error(`unsupported rusqsieve wasm ABI ${abi}`);
      if (typeof ex.qs_ecm !== "function") throw new Error("wasm module exports no qs_ecm");
      self.postMessage({ type: "ready", gen, abi });
      return;
    }
    if (data.cmd === "search") {
      if (!ex) throw new Error("ecm worker is not initialized");
      if (!Number.isSafeInteger(data.curves) || data.curves <= 0) {
        throw new Error("ecm search requires a positive curve count");
      }
      const decimal = validateDecimalInput(data.n, MAX_MODULUS_BITS);
      const input = putString(ex, decimal);
      let handle;
      try {
        // Zero bounds take the module's own schedule for this width, so the glue does not
        // duplicate a table that lives in the engine.
        handle = ex.qs_ecm(input.ptr, input.len, data.b1 >>> 0, data.b2 >>> 0, data.curves, data.seed >>> 0);
      } finally {
        ex.qs_dealloc(input.ptr, input.len, 1);
      }
      if (!handle) throw new Error("wasm rejected the ecm modulus");
      const payload = takePacket(ex, handle, 13);
      // An empty payload means every curve was spent without a split, which is a result rather
      // than a failure: the caller decides what happens to the composite next.
      const factor = payload.length === 0 ? null : bytesToBigInt(payload);
      self.postMessage({
        type: "done",
        gen,
        factor: factor === null ? null : factor.toString(),
      });
      return;
    }
    throw new Error("invalid ecm worker command");
  } catch (error) {
    self.postMessage({ type: "error", gen, message: String(error?.message ?? error) });
  }
};
