// End-to-end check of the browser factorization architecture, without a browser.
//
// Drives the exact message protocol web/index.js uses — a dedicated coordinator Worker
// running web/coordinator.js plus a pool of sieve Workers running web/worker.js — over
// node worker threads, and asserts a known semiprime comes back correctly factored.
//
// This exists because the architecture is otherwise only exercised by hand in a browser:
// `docs/` once shipped without coordinator.js at all, which 404s and leaves boot()
// awaiting a "ready" that never arrives, and nothing in the test suite noticed.
//
// usage: node tools/browser-arch-check.mjs [DECIMAL] [WORKERS]
import { readFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";
import { bytesToBigInt } from "../web/abi.js";

// A balanced semiprime from the test corpus: 21293688545713669 * 31385813854515511.
// Must be a genuine sieve input — this driver hands it straight to the engine, unlike
// web/index.js, which peels small factors and easy cases off first.
const decimal = process.argv[2] || "668319744971798315493259725219859";
const nWorkers = Number(process.argv[3] || 4);
const expected = BigInt(decimal);

const wasmPath =
  process.env.RUSQSIEVE_WASM ||
  new URL("../target/wasm-scalar/wasm32-unknown-unknown/release/rusqsieve.wasm", import.meta.url);
const module = await WebAssembly.compile(await readFile(wasmPath));

const spawn = (shim) => new Worker(new URL(shim, import.meta.url), { type: "module" });
const ready = (w, cmd) =>
  new Promise((resolve, reject) => {
    w.once("message", (data) => (data.type === "error" ? reject(new Error(data.error)) : resolve(data)));
    w.postMessage(cmd);
  });

const coord = spawn("./node-coordinator-shim.mjs");
const abi = (await ready(coord, { cmd: "init", module })).abi;
const workers = [];
for (let i = 0; i < nWorkers; i++) {
  const w = spawn("./node-worker-shim.mjs");
  await ready(w, { cmd: "init", module });
  workers.push(w);
}
console.log(`ABI v${abi}, ${nWorkers} sieve workers + 1 coordinator worker`);

const BATCH = 2;
// Read from the session rather than hard-coded: the engine scales the family budget with input
// width, and a stale constant here would stop a large run before the engine would.
let familyBudget = 0;
const gen = 1;
let nextFamily = 0;
let finished = false;
let target = 0;
let sawLinalg = false;

const started = performance.now();
const factor = await new Promise((resolve, reject) => {
  const fail = (message) => {
    if (!finished) {
      finished = true;
      reject(new Error(message));
    }
  };
  const dispatch = (w) => {
    if (finished) return;
    const family = nextFamily;
    nextFamily += BATCH;
    if (familyBudget && nextFamily > familyBudget) return fail("relation budget exhausted");
    w.postMessage({ cmd: "sieve", family, count: BATCH, gen });
  };

  coord.on("message", (data) => {
    if (data.type === "error") return fail(data.error);
    if (data.gen !== gen) return;
    if (data.type === "session") {
      target = data.target;
      familyBudget = data.familyBudget;
      if (!Number.isInteger(familyBudget) || familyBudget <= 0) {
        return fail("coordinator returned an invalid family budget");
      }
      for (const w of workers) w.postMessage({ cmd: "prepare", n: decimal, gen });
    } else if (data.type === "submitted") {
      if (!finished) dispatch(workers[data.worker]);
    } else if (data.type === "linalg") {
      sawLinalg = true;
    } else if (data.type === "factor") {
      finished = true;
      resolve(bytesToBigInt(data.factor));
    }
  });

  workers.forEach((w, workerIndex) => {
    w.on("message", (data) => {
      if (data.type === "error") return fail(data.error);
      if (finished || data.gen !== gen) return;
      if (data.type === "prepared") {
        if (!data.ok) return fail("worker could not build a sieve");
        dispatch(w);
      } else if (data.type === "relations") {
        if (data.payload) {
          coord.postMessage({ cmd: "submit", payload: data.payload, worker: workerIndex, gen });
        } else {
          dispatch(w);
        }
      }
    });
  });

  coord.postMessage({ cmd: "new", n: decimal, gen });
});
const elapsed = ((performance.now() - started) / 1000).toFixed(2);

await Promise.all([coord.terminate(), ...workers.map((w) => w.terminate())]);

const cofactor = expected / factor;
if (factor <= 1n || factor >= expected || factor * cofactor !== expected) {
  console.error(`FAIL: ${factor} is not a nontrivial factor of ${expected}`);
  process.exit(1);
}
if (!sawLinalg) {
  console.error("FAIL: the coordinator never reported the linear-algebra phase");
  process.exit(1);
}
console.log(`ok: ${expected} = ${factor} * ${cofactor}`);
console.log(`    relation target ${target}, ${nextFamily} families dispatched, ${elapsed}s`);
