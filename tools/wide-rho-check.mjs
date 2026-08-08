// The frontend's half of the deep Pollard-Brent policy.
//
// A composite the sieve cannot help with — one above its ceiling, or a cofactor an earlier split
// already proved unbalanced — gets a deep rho before the frontend gives up on it or hands it over.
// That search runs in wasm across a pool of rho workers racing disjoint polynomial constants, with
// the main thread's BigInt implementation left as the fallback for a runtime without either.
//
// This drives the real worker protocol on node worker threads, exactly as browser-arch-check.mjs
// does for the sieve, and checks the properties the policy rests on: the wasm path returns the same
// factors the BigInt path does, disjoint constants each find one, the ABI guard is real, and the
// fallback still yields often enough to keep a page painting.
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";

import { pollardBrent, pollardBrentSliced, RHO_CONSTANTS } from "../web/numtheory.js";

// Same inputs as the native `wide_composites_are_split_by_rho_rather_than_refused` test: one small
// prime times one wide prime, at 512 and 1024 bits.
const CASES = [
  {
    bits: 512,
    n: 9501012405705509564680437712617447440170980081112656222237073910419870316392859702111963091481439276805995800801743430916377894473378632368751322056628119n,
    factor: 3667435003n,
    budget: 8 << 20,
    // Brent's cost varies by several-fold around its `1.2·sqrt(p)` expectation, so whether the
    // opening peel happens to reach a given 32-bit factor is luck. This one it misses.
    openingPeelReaches: false,
  },
  {
    bits: 1024,
    n: 140102229730795429799167923021188282773066057852400709530422767028102695167748110966423996618253720918241568745283524257975332015353165215374100761210034074640519308908102560817045658913006325713539052727661523689985681266686274209226917233091828080256346478868845059891172761884504662361609922427332741249613n,
    factor: 3479286313n,
    budget: 4 << 20,
    // ...and this one it happens to hit, which is the whole argument for a guaranteed tier sized
    // at hundreds of times the expected cost rather than at it.
    openingPeelReaches: true,
  },
];

const WASM =
  process.env.RUSQSIEVE_WASM || "target/wasm-scalar/wasm32-unknown-unknown/release/rusqsieve.wasm";

function spawnRhoWorker() {
  return new Worker(new URL("./node-rho-worker-shim.mjs", import.meta.url), { type: "module" });
}

// One request/response exchange against a freshly initialized rho worker.
function search(module, request) {
  return new Promise((resolve, reject) => {
    const worker = spawnRhoWorker();
    const finish = (value, error) => {
      worker.terminate();
      if (error) reject(error);
      else resolve(value);
    };
    worker.on("message", (data) => {
      if (data?.type === "ready") {
        worker.postMessage({ cmd: "search", gen: 1, ...request });
      } else if (data?.type === "done") {
        finish(data.factor === null ? null : BigInt(data.factor));
      } else if (data?.type === "error") {
        finish(null, new Error(data.message));
      }
    });
    worker.on("error", (error) => finish(null, error));
    worker.postMessage({ cmd: "init", module });
  });
}

const module = await WebAssembly.compile(await readFile(WASM));
const exports = (await WebAssembly.instantiate(module, {})).exports;
assert.equal(exports.qs_abi_version(), 4, "the rho worker requires wasm ABI 4");
assert.equal(typeof exports.qs_rho, "function", "wasm module exports no qs_rho");

// The wrapper must keep behaving exactly as it did before it was expressed in terms of the sliced
// generator: below the ceiling this is the cheap main-thread peel, deliberately small.
assert.equal(pollardBrent(10403n, 1 << 15), 101n);
assert.equal(pollardBrent(2n * 5717n, 1 << 15), 2n);

for (const { bits, n, factor, budget, openingPeelReaches } of CASES) {
  // The gap this closes: the opening peel reaches a 32-bit factor only by luck, and refusing the
  // composite on the strength of it is what made a findable factor look unfindable.
  assert.equal(
    pollardBrent(n, 1 << 15),
    openingPeelReaches ? factor : null,
    `${bits}-bit: the opening peel's reach is not what this check assumes`,
  );

  // The wasm path, one worker, through the same message protocol index.js uses.
  const started = performance.now();
  const found = await search(module, { n: n.toString(), budget, first: 1, count: 4 });
  const wasmElapsed = performance.now() - started;
  assert.equal(found, factor, `${bits}-bit: wasm rho returned ${found}`);

  // Racing: a second worker walking a disjoint constant range must find the same factor, since
  // every constant is an independent walk over the same modulus. This is what makes the pool safe.
  const raced = await search(module, { n: n.toString(), budget, first: 5, count: 4 });
  assert.equal(raced, factor, `${bits}-bit: disjoint constants returned ${raced}`);

  // The fallback path, still sliced so the main thread keeps painting.
  const fallbackStarted = performance.now();
  const run = pollardBrentSliced(n, budget);
  let slices = 0;
  let worstSlice = 0;
  let last = performance.now();
  let step = run.next();
  while (!step.done) {
    const now = performance.now();
    worstSlice = Math.max(worstSlice, now - last);
    last = now;
    slices += 1;
    step = run.next();
  }
  const fallbackElapsed = performance.now() - fallbackStarted;
  assert.equal(step.value, factor, `${bits}-bit: fallback returned ${step.value}`);
  assert.ok(slices > 0, `${bits}-bit: the fallback never yielded to the event loop`);
  // A slice that runs long enough to drop frames defeats the point of slicing. The threshold is
  // loose because this runs on unknown CI hardware; the measured worst slice is around 50 ms.
  assert.ok(worstSlice < 500, `${bits}-bit: worst slice ${worstSlice.toFixed(0)}ms`);

  console.log(
    `ok: ${bits}-bit composite -> ${factor} · wasm ${wasmElapsed.toFixed(0)}ms ` +
      `(and ${raced} from constants 5-8) · BigInt fallback ${fallbackElapsed.toFixed(0)}ms, ` +
      `${slices} slices, worst ${worstSlice.toFixed(0)}ms`,
  );
}

// Every constant the pool can hand out has to be one the worker accepts, or a pool sized to the
// host would silently lose walks.
assert.equal(RHO_CONSTANTS.length, 31);
await assert.rejects(
  search(module, { n: "0", budget: 1024, first: 1, count: 4 }),
  /positive|decimal|rejected/,
  "the worker accepted an unusable modulus",
);
console.log("ok: rho worker rejects an unusable modulus and reports ABI 4");
