// Focused non-happy-path checks for the reference browser glue.
//
// This complements browser-arch-check.mjs: that script proves a factor reaches
// the coordinator through real worker modules, while this one verifies that
// malformed ABI packets and invalid worker command sequences terminate with
// explicit errors and release owned buffers.
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";

import {
  bytesToBigInt,
  putBytes,
  takePacket,
  validateDecimalInput,
} from "../web/abi.js";

checkPacketValidation();
checkAllocationRollback();
assert.equal(bytesToBigInt(Uint8Array.of(0x34, 0x12)), 0x1234n);
assert.equal(validateDecimalInput("00015"), "15");
assert.throws(() => validateDecimalInput("12x"), /unsigned decimal/);
assert.throws(() => validateDecimalInput("0"), /positive/);
assert.throws(() => validateDecimalInput((1n << 400n).toString()), /400-bit limit/);
assert.equal(validateDecimalInput((1n << 399n).toString()), (1n << 399n).toString());

const wasmPath =
  process.env.RUSQSIEVE_WASM ||
  new URL("../target/wasm-scalar/wasm32-unknown-unknown/release/rusqsieve.wasm", import.meta.url);
const module = await WebAssembly.compile(await readFile(wasmPath));

const sieve = new Worker(new URL("./node-worker-shim.mjs", import.meta.url), { type: "module" });
const coordinator = new Worker(new URL("./node-coordinator-shim.mjs", import.meta.url), {
  type: "module",
});
try {
  assert.equal((await request(sieve, { cmd: "init", module })).type, "ready");
  const unprepared = await request(sieve, {
    cmd: "sieve",
    family: 0,
    count: 1,
    gen: 51,
  });
  assert.equal(unprepared.type, "error");
  assert.equal(unprepared.gen, 51);
  assert.match(unprepared.error, /no prepared sieve context/);

  const unknownWorker = await request(sieve, { cmd: "not-a-command", gen: 52 });
  assert.equal(unknownWorker.type, "error");
  assert.equal(unknownWorker.gen, 52);
  assert.match(unknownWorker.error, /unknown worker command/);

  const oversizedWorkerInput = await request(sieve, {
    cmd: "prepare",
    n: (1n << 400n).toString(),
    gen: 53,
  });
  assert.equal(oversizedWorkerInput.type, "error");
  assert.equal(oversizedWorkerInput.gen, 53);
  assert.match(oversizedWorkerInput.error, /400-bit limit/);

  const ready = await request(coordinator, { cmd: "init", module });
  assert.equal(ready.type, "ready");
  assert.equal(ready.abi, 2);
  // The scheduler's family cap and the UI's width limit are both engine-owned now; a runtime that
  // does not report them would leave the glue silently guessing.
  assert.equal(ready.maxSiqsBits, 400);

  const decimal = "668319744971798315493259725219859";
  const session = await request(coordinator, { cmd: "new", n: decimal, gen: 61 });
  assert.equal(session.type, "session");
  assert.equal(session.gen, 61);
  assert.ok(session.target > 0);
  assert.ok(Number.isInteger(session.familyBudget) && session.familyBudget > 0);

  const malformed = await request(coordinator, {
    cmd: "submit",
    payload: Uint8Array.of(1, 2),
    worker: 0,
    gen: 61,
  });
  assert.equal(malformed.type, "error");
  assert.equal(malformed.gen, 61);
  assert.match(malformed.error, /truncated relation batch/);

  const unknownCoordinator = await request(coordinator, { cmd: "not-a-command", gen: 61 });
  assert.equal(unknownCoordinator.type, "error");
  assert.equal(unknownCoordinator.gen, 61);
  assert.match(unknownCoordinator.error, /unknown coordinator command/);

  // Last: `new` resets the session and generation before validating, so a rejected session must
  // not run ahead of the gen-61 assertions above.
  const oversizedCoordinatorInput = await request(coordinator, {
    cmd: "new",
    n: (1n << 400n).toString(),
    gen: 62,
  });
  assert.equal(oversizedCoordinatorInput.type, "error");
  assert.match(oversizedCoordinatorInput.error, /400-bit limit/);
} finally {
  await Promise.all([sieve.terminate(), coordinator.terminate()]);
}

console.log("ok: browser glue rejects malformed packets and invalid worker command sequences");

function checkPacketValidation() {
  const memory = new WebAssembly.Memory({ initial: 1 });
  const bytes = new Uint8Array(memory.buffer);
  let pointer = 32;
  let length = 0;
  let freed = 0;
  const ex = {
    memory,
    qs_buffer_pointer: () => pointer,
    qs_buffer_length: () => length,
    qs_buffer_free: () => {
      freed++;
    },
  };
  const writePacket = (kind, payload, { magic = "QSV1", version = 1, payloadLength } = {}) => {
    const encodedMagic = new TextEncoder().encode(magic);
    bytes.set(encodedMagic.subarray(0, 4), pointer);
    const view = new DataView(memory.buffer, pointer);
    view.setUint16(4, kind, true);
    view.setUint16(6, version, true);
    view.setUint32(8, payloadLength ?? payload.byteLength, true);
    bytes.set(payload, pointer + 12);
    length = 12 + payload.byteLength;
  };

  writePacket(10, Uint8Array.of(7, 8, 9));
  assert.deepEqual(takePacket(ex, 1, 10), Uint8Array.of(7, 8, 9));
  assert.equal(freed, 1);

  writePacket(10, Uint8Array.of(1), { magic: "BAD!" });
  assert.throws(() => takePacket(ex, 2, 10), /packet magic/);
  assert.equal(freed, 2);

  writePacket(11, Uint8Array.of(1));
  assert.throws(() => takePacket(ex, 3, 10), /packet kind/);
  assert.equal(freed, 3);

  writePacket(10, Uint8Array.of(1, 2), { payloadLength: 20 });
  assert.throws(() => takePacket(ex, 4, 10), /payload length/);
  assert.equal(freed, 4);

  pointer = 0;
  length = 12;
  assert.throws(() => takePacket(ex, 5, 10), /stale wasm packet handle/);
  assert.equal(freed, 5);
}

function checkAllocationRollback() {
  const memory = new WebAssembly.Memory({ initial: 1 });
  let deallocated = null;
  const ex = {
    memory,
    qs_alloc: () => memory.buffer.byteLength - 1,
    qs_dealloc: (...args) => {
      deallocated = args;
    },
  };
  assert.throws(() => putBytes(ex, Uint8Array.of(1, 2)), RangeError);
  assert.deepEqual(deallocated, [memory.buffer.byteLength - 1, 2, 1]);
}

function request(worker, message) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => finish(new Error("worker response timed out")), 10_000);
    const onMessage = (data) => finish(null, data);
    const onError = (error) => finish(error);
    const finish = (error, data) => {
      clearTimeout(timer);
      worker.off("message", onMessage);
      worker.off("error", onError);
      if (error) reject(error);
      else resolve(data);
    };
    worker.on("message", onMessage);
    worker.on("error", onError);
    worker.postMessage(message);
  });
}
