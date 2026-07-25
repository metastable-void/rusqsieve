// Node/V8 benchmark for the independent-instance sieve-Worker architecture.
//
// This measures sieve throughput, and runs the coordinator inline on the main thread. The browser
// demo no longer does: since the coordinator moved into its own Worker, relation collection and the
// GF(2) solve run off the main thread there. That difference does not affect what this benchmark
// measures, but it means this is not an end-to-end model of the page — see
// `tools/browser-arch-check.mjs`, which drives the real two-kinds-of-Worker protocol.
import { readFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";
import { performance } from "node:perf_hooks";
import { instantiate, putBytes, putString, takePacket, bytesToBigInt } from "../web/abi.js";

const decimal = process.argv[2];
const threadCount = Number(process.argv[3] || 8);
const batch = Number(process.argv[4] || 2);
if (
  !/^\d+$/.test(decimal || "") ||
  !Number.isInteger(threadCount) ||
  threadCount < 1 ||
  !Number.isInteger(batch) ||
  batch < 1
) {
  console.error("usage: node tools/wasm-bench.mjs DECIMAL [WORKERS] [FAMILIES_PER_JOB]");
  process.exit(2);
}

const wasmPath =
  process.env.RUSQSIEVE_WASM ||
  new URL(
    "../target/wasm-simd/wasm32-unknown-unknown/release/rusqsieve.wasm",
    import.meta.url,
  );
const module = await WebAssembly.compile(await readFile(wasmPath));
const coord = await instantiate(module);
const workers = await Promise.all(
  Array.from({ length: threadCount }, async () => {
    const worker = new Worker(new URL("./node-worker-shim.mjs", import.meta.url), {
      type: "module",
    });
    await request(worker, { cmd: "init", module }, "ready");
    return worker;
  }),
);

const input = putString(coord, decimal);
const session = coord.qs_coord_new(input.ptr, input.len);
coord.qs_dealloc(input.ptr, input.len, 1);
if (!session) throw new Error("could not prepare coordinator");
const target = coord.qs_coord_target(session);
let relations = 0;
let nextFamily = 0;
const started = performance.now();

await Promise.all(
  workers.map(async (worker) => {
    await request(worker, { cmd: "prepare", n: decimal, gen: 1 }, "prepared");
    while (relations < target) {
      const family = nextFamily;
      nextFamily += batch;
      const result = await request(
        worker,
        { cmd: "sieve", family, count: batch, gen: 1 },
        "relations",
      );
      if (!result.payload || relations >= target) continue;
      const bytes = putBytes(coord, result.payload);
      relations = coord.qs_coord_submit(session, bytes.ptr, bytes.len);
      coord.qs_dealloc(bytes.ptr, bytes.len, 1);
    }
  }),
);
const sieveDone = performance.now();

const handle = coord.qs_coord_extract(session);
const factorBytes = takePacket(coord, handle);
coord.qs_coord_free(session);
for (const worker of workers) await worker.terminate();
if (!factorBytes) throw new Error("linear algebra found no factor");
const factor = bytesToBigInt(factorBytes);
const inputBig = BigInt(decimal);
if (factor <= 1n || factor >= inputBig || inputBig % factor !== 0n) {
  throw new Error("invalid factor");
}
console.log(factor.toString());
console.log((inputBig / factor).toString());
const elapsed = (performance.now() - started) / 1000;
const sieveElapsed = (sieveDone - started) / 1000;
console.error(
  `wasm-v8 workers=${threadCount} batch=${batch} relations=${relations}/${target} ` +
    `families=${nextFamily} sieve=${sieveElapsed.toFixed(3)}s ` +
    `finish=${(elapsed - sieveElapsed).toFixed(3)}s wall=${elapsed.toFixed(3)}s`,
);

function request(worker, message, expectedType) {
  return new Promise((resolve, reject) => {
    const onMessage = (data) => {
      if (data.type === "error") {
        cleanup();
        reject(new Error(data.error));
      } else if (data.type === expectedType) {
        cleanup();
        resolve(data);
      }
    };
    const onError = (error) => {
      cleanup();
      reject(error);
    };
    const cleanup = () => {
      worker.off("message", onMessage);
      worker.off("error", onError);
    };
    worker.on("message", onMessage);
    worker.on("error", onError);
    worker.postMessage(message);
  });
}
