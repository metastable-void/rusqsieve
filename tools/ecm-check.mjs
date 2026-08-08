// The frontend's last stage, driven through its real worker protocol.
//
// A composite the sieve refuses and rho cannot split has exactly one thing left, and it is this.
// The check is the property that matters: a wide number with a medium-size factor — the shape that
// used to be reported as "too large" — comes back factored.
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";

const WASM =
  process.env.RUSQSIEVE_WASM || "target/wasm-scalar/wasm32-unknown-unknown/release/rusqsieve.wasm";
const module = await WebAssembly.compile(await readFile(WASM));
const exports = (await WebAssembly.instantiate(module, {})).exports;
assert.equal(exports.qs_abi_version(), 5, "the ecm worker requires wasm ABI 5");
assert.equal(typeof exports.qs_ecm, "function", "wasm module exports no qs_ecm");
assert.ok(exports.qs_ecm_default_b1(512, 1) > 0, "no default B1 for a 512-bit composite");
assert.ok(exports.qs_ecm_default_curves(512, 1) > 0, "no default curve count");
// The two schedules must actually differ: the cheap one runs in front of a sieve that will finish
// anyway, and sizing it like the committed one would delay that sieve for a lottery ticket.
assert.ok(
  exports.qs_ecm_default_b1(384, 0) < exports.qs_ecm_default_b1(384, 1),
  "the uncommitted schedule is not cheaper than the committed one",
);

function search(request) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL("./node-ecm-worker-shim.mjs", import.meta.url), {
      type: "module",
    });
    const finish = (value, error) => {
      worker.terminate();
      if (error) reject(error);
      else resolve(value);
    };
    worker.on("message", (data) => {
      if (data?.type === "ready") worker.postMessage({ cmd: "search", gen: 1, ...request });
      else if (data?.type === "done") finish(data.factor === null ? null : BigInt(data.factor));
      else if (data?.type === "error") finish(null, new Error(data.message));
    });
    worker.on("error", (error) => finish(null, error));
    worker.postMessage({ cmd: "init", module });
  });
}

// 4,294,967,311 × a 128-bit prime: a ten-digit factor, which stage 1 reaches quickly at the
// default bounds for this width.
const n = 1461501642435138422017761784666902131351499047677n;
const started = performance.now();
const found = await search({ n: n.toString(), b1: 2000, b2: 200000, curves: 64, seed: 1 });
assert.equal(found, 4294967311n, `ecm returned ${found}`);
console.log(`ok: wasm ECM split a 161-bit composite -> ${found} (${(performance.now() - started).toFixed(0)}ms)`);

await assert.rejects(
  search({ n: "0", b1: 0, b2: 0, curves: 4, seed: 1 }),
  /positive|decimal|rejected/,
  "the worker accepted an unusable modulus",
);
console.log("ok: ecm worker rejects an unusable modulus and reports ABI 5");
