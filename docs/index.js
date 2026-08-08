// Main thread: UI + coordinator. Peels easy factors with BigInt number theory and
// hands hard composites to a pool of wasm Web Workers running the quadratic sieve,
// with the pool sized to navigator.hardwareConcurrency.
import { loadModule, bytesToBigInt } from "./abi.js";
import { trialDivide, isPrime, perfectPower, pollardBrent, pollardBrentSliced, groupFactors, rsaNumber, bitLength } from "./numtheory.js";

const SIMD_WASM_URL = new URL("./rusqsieve-simd.wasm", import.meta.url);
const SCALAR_WASM_URL = new URL("./rusqsieve.wasm", import.meta.url);
// Small jobs reduce the tail after the relation target is reached. Two
// families was consistently best in Node/V8 from 192 through 256 bits.
const BATCH = 2;
// The engine's family budget scales with input width, so it is read per session from the
// coordinator (`qs_coord_family_budget`) rather than hard-coded here. A stale constant would
// either stop a large run early or issue families the engine has already stopped accepting.
//
// Input width alone no longer limits a run: the number is peeled with BigInt trial division,
// perfect-power detection, and Pollard-Brent first, and only the hard composite that survives is
// range-limited. These two bound what `Natural` can hold (PARTS = 16, so 1024 bits); the sieve's
// own ceiling arrives from the coordinator as `maxSiqsBits`.
const MAX_INPUT_BITS = 1024;
const MAX_DECIMAL_DIGITS = 309;
// Set from the coordinator's ready message; the fallback matches `engine::MAX_SIQS_BITS` and is
// only used if a runtime somehow reports nothing.
let maxSiqsBits = 400;
const BOOT_TIMEOUT_MS = 30_000;
// A hard wall-clock limit rejects large inputs that are still making useful
// progress. Only abandon a run when the whole worker pool and coordinator have
// been silent while the page is visible. A single stuck worker may reduce
// throughput, but the remaining workers can still complete the factorization.
const STALL_TIMEOUT_MS = 10 * 60_000;
// If the watchdog callback itself was delayed this much, the browser or machine
// was suspended; that gap is not evidence that the worker runtime hung.
const WATCHDOG_LATE_GRACE_MS = 30_000;

const els = {
  input: document.getElementById("input"),
  inputMirror: document.getElementById("input-mirror"),
  inputInfo: document.getElementById("input-info"),
  go: document.getElementById("go"),
  bar: document.getElementById("bar"),
  status: document.getElementById("status"),
  result: document.getElementById("result"),
  workers: document.getElementById("workers"),
  meter: document.getElementById("meter"),
  rsaBits: document.getElementById("rsa-bits"),
  rsaBitsLabel: document.getElementById("rsa-bits-label"),
  rsaGen: document.getElementById("rsa-gen"),
};

let coord = null; // coordinator Worker (owns its own wasm instance)
let workers = []; // sieve worker pool
// The compiled module, kept so a deep rho search can spin up its own short-lived worker pool
// without recompiling. Sieve workers and the coordinator each instantiate their own copy of it.
let wasmModule = null;
let gen = 0; // generation token so stale worker messages are ignored
let wasmFlavor = "scalar";
let runtimeReady = false;
// Scaling remains positive through 32–48 workers on the 96-thread reference
// host, while 96 workers regress from startup, memory traffic, and job overshoot.
const nWorkers = Math.max(1, Math.min(48, navigator.hardwareConcurrency || 4));

async function boot() {
  runtimeReady = false;
  let module;
  try {
    module = await withTimeout(loadModule(SIMD_WASM_URL), BOOT_TIMEOUT_MS, "SIMD wasm load");
    wasmFlavor = "SIMD";
  } catch {
    // Older engines can still use the portable artifact.
    module = await withTimeout(
      loadModule(SCALAR_WASM_URL),
      BOOT_TIMEOUT_MS,
      "scalar wasm load",
    );
    wasmFlavor = "scalar";
  }
  const nextCoord = new Worker(new URL("./coordinator.js", import.meta.url), { type: "module" });
  const nextWorkers = Array.from(
    { length: nWorkers },
    () => new Worker(new URL("./worker.js", import.meta.url), { type: "module" }),
  );
  const bootAbort = new AbortController();
  try {
    const [coordinatorReady] = await Promise.all([
      waitForWorkerReady(nextCoord, module, true, bootAbort.signal),
      ...nextWorkers.map((worker) =>
        waitForWorkerReady(worker, module, false, bootAbort.signal),
      ),
    ]);
    coord = nextCoord;
    workers = nextWorkers;
    wasmModule = module;
    runtimeReady = true;
    // Take the sieve's range from the engine that will actually run, so the UI limit and the
    // coordinator's own check can never drift apart across a rebuild.
    if (Number.isInteger(coordinatorReady.maxSiqsBits) && coordinatorReady.maxSiqsBits > 0) {
      maxSiqsBits = coordinatorReady.maxSiqsBits;
    }
    els.workers.textContent =
      `${nWorkers} worker${nWorkers === 1 ? "" : "s"} · ${wasmFlavor} · ` +
      `ABI v${coordinatorReady.abi}`;
  } catch (error) {
    bootAbort.abort();
    nextCoord.terminate();
    for (const worker of nextWorkers) worker.terminate();
    throw error;
  }
  els.go.disabled = false;
  els.status.textContent = "Ready.";
}

function shutdownRuntime() {
  runtimeReady = false;
  coord?.terminate();
  for (const worker of workers) worker.terminate();
  coord = null;
  workers = [];
  wasmModule = null;
}

async function restartRuntime() {
  shutdownRuntime();
  els.go.disabled = true;
  await boot();
}

function waitForWorkerReady(worker, module, requireAbi, signal) {
  return new Promise((resolve, reject) => {
    let settled = false;
    const timer = setTimeout(
      () => finish(new Error("worker initialization timed out")),
      BOOT_TIMEOUT_MS,
    );
    const cleanup = () => {
      clearTimeout(timer);
      worker.removeEventListener("message", onMessage);
      worker.removeEventListener("error", onError);
      worker.removeEventListener("messageerror", onMessageError);
      signal?.removeEventListener("abort", onAbort);
    };
    const finish = (error, data) => {
      if (settled) return;
      settled = true;
      cleanup();
      if (error) reject(error);
      else resolve(data);
    };
    const onMessage = ({ data }) => {
      if (data?.type === "error") {
        finish(new Error(data.error || "worker initialization failed"));
      } else if (data?.type === "ready") {
        if (requireAbi && data.abi !== 5) {
          finish(new Error(`unsupported rusqsieve wasm ABI ${String(data.abi)}`));
        } else {
          finish(null, data);
        }
      }
    };
    const onError = (event) => {
      event.preventDefault?.();
      finish(new Error(event.message || "worker failed during initialization"));
    };
    const onMessageError = () => finish(new Error("worker initialization message was invalid"));
    const onAbort = () => finish(new Error("worker initialization cancelled"));
    worker.addEventListener("message", onMessage);
    worker.addEventListener("error", onError);
    worker.addEventListener("messageerror", onMessageError);
    signal?.addEventListener("abort", onAbort, { once: true });
    if (signal?.aborted) {
      onAbort();
      return;
    }
    try {
      worker.postMessage({ cmd: "init", module });
    } catch (error) {
      finish(error);
    }
  });
}

function withTimeout(promise, milliseconds, label) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => reject(new Error(`${label} timed out`)), milliseconds);
    Promise.resolve(promise).then(
      (value) => {
        clearTimeout(timer);
        resolve(value);
      },
      (error) => {
        clearTimeout(timer);
        reject(error);
      },
    );
  });
}

// Parallel quadratic sieve for one hard composite; resolves to a nontrivial factor.
function siqsParallel(decimal, bits, report) {
  return new Promise((resolve, reject) => {
    const myGen = ++gen;
    const sieveStarted = performance.now();
    let target = 0;
    let relations = 0;
    let nextFamily = 0;
    let familyBudget = 0;
    let activeJobs = 0;
    let pendingSubmissions = 0;
    let preparedWorkers = 0;
    let finished = false;
    const workerBusy = new Array(workers.length).fill(false);
    const workerPrepared = new Array(workers.length).fill(false);
    let stallTimer = null;
    let lastActivity = performance.now();

    const clearStallTimer = () => {
      if (stallTimer !== null) clearTimeout(stallTimer);
      stallTimer = null;
    };
    const armStallWatchdog = () => {
      clearStallTimer();
      if (finished || document.hidden) return;
      const remaining = Math.max(
        0,
        STALL_TIMEOUT_MS - (performance.now() - lastActivity),
      );
      const deadline = performance.now() + remaining;
      stallTimer = setTimeout(() => {
        const now = performance.now();
        if (now - deadline > WATCHDOG_LATE_GRACE_MS) {
          lastActivity = now;
          armStallWatchdog();
          return;
        }
        const idleMilliseconds = now - lastActivity;
        if (!finished && !document.hidden && idleMilliseconds >= STALL_TIMEOUT_MS) {
          fail(new Error("factorization stalled: no worker activity for 10 minutes"));
        } else {
          armStallWatchdog();
        }
      }, remaining);
    };
    const noteActivity = () => {
      lastActivity = performance.now();
      armStallWatchdog();
    };
    const onVisibilityChange = () => {
      if (document.hidden) {
        clearStallTimer();
      } else {
        // Browser suspension looks like a long wall-clock gap. Give workers a
        // fresh interval after the page becomes active before declaring a hang.
        noteActivity();
      }
    };

    const cleanup = () => {
      clearStallTimer();
      document.removeEventListener("visibilitychange", onVisibilityChange);
      coord.onmessage = null;
      coord.onerror = null;
      coord.onmessageerror = null;
      for (const worker of workers) {
        worker.onmessage = null;
        worker.onerror = null;
        worker.onmessageerror = null;
      }
    };
    const fail = (error) => {
      if (finished) return;
      finished = true;
      cleanup();
      reject(error instanceof Error ? error : new Error(String(error)));
    };
    const succeed = (factor) => {
      if (finished) return;
      finished = true;
      cleanup();
      resolve(factor);
    };
    const maybeExhausted = () => {
      if (
        !finished &&
        familyBudget > 0 &&
        nextFamily >= familyBudget &&
        activeJobs === 0 &&
        pendingSubmissions === 0 &&
        preparedWorkers === workers.length
      ) {
        fail(
          new Error(
            `relation budget exhausted after ${familyBudget} families ` +
              `(${relations}/${target} relations)`,
          ),
        );
      }
    };
    const dispatch = (worker, workerIndex) => {
      if (finished || workerBusy[workerIndex]) return false;
      if (familyBudget <= 0 || nextFamily >= familyBudget) {
        maybeExhausted();
        return false;
      }
      const family = nextFamily;
      const count = Math.min(BATCH, familyBudget - nextFamily);
      nextFamily += count;
      workerBusy[workerIndex] = true;
      activeJobs++;
      try {
        worker.postMessage({ cmd: "sieve", family, count, gen: myGen });
      } catch (error) {
        workerBusy[workerIndex] = false;
        activeJobs--;
        fail(error);
        return false;
      }
      return true;
    };
    const finishJob = (workerIndex) => {
      if (!workerBusy[workerIndex]) {
        fail(new Error(`unexpected response from idle sieve worker ${workerIndex + 1}`));
        return false;
      }
      workerBusy[workerIndex] = false;
      activeJobs--;
      return true;
    };

    coord.onmessage = ({ data }) => {
      if (data?.gen !== myGen) return;
      noteActivity();
      if (data.type === "error") {
        fail(new Error(data.error || "coordinator failed"));
      } else if (data.type === "session") {
        if (!Number.isInteger(data.target) || data.target <= 0) {
          fail(new Error("coordinator returned an invalid relation target"));
          return;
        }
        if (!Number.isInteger(data.familyBudget) || data.familyBudget <= 0) {
          fail(new Error("coordinator returned an invalid family budget"));
          return;
        }
        target = data.target;
        familyBudget = data.familyBudget;
        try {
          for (const w of workers) {
            w.postMessage({ cmd: "prepare", n: decimal, gen: myGen });
          }
        } catch (error) {
          fail(error);
        }
      } else if (data.type === "submitted") {
        if (pendingSubmissions <= 0) {
          fail(new Error("coordinator acknowledged an unknown submission"));
          return;
        }
        pendingSubmissions--;
        if (
          !Number.isInteger(data.worker) ||
          data.worker < 0 ||
          data.worker >= workers.length ||
          !Number.isInteger(data.relations) ||
          data.relations < relations ||
          !Number.isInteger(data.target) ||
          data.target <= 0
        ) {
          fail(new Error("coordinator returned invalid progress"));
          return;
        }
        relations = data.relations;
        target = data.target;
        const now = performance.now();
        const elapsedSeconds = (now - sieveStarted) / 1000;
        // The accepted relation count accelerates as the partial-relation graph
        // accumulates edges and closes more cycles. A linear rate extrapolation
        // is therefore wildly pessimistic early in 272-bit runs. The measured
        // browser curve is close to relations ∝ time^1.6; invert that curve to
        // estimate total sieve time without embedding a machine-specific rate.
        const progress = target > 0 ? Math.min(1, relations / target) : 0;
        const etaSeconds =
          progress >= 0.03 ? elapsedSeconds * (progress ** (-1 / 1.6) - 1) : null;
        report({
          phase: "sieving",
          bits,
          relations,
          target,
          elapsedSeconds,
          etaSeconds,
        });
        if (!finished) {
          dispatch(workers[data.worker], data.worker);
          maybeExhausted();
        }
      } else if (data.type === "linalg") {
        report({ phase: "linalg" });
      } else if (data.type === "factor") {
        if (!(data.factor instanceof Uint8Array)) {
          fail(new Error("coordinator returned a malformed factor"));
          return;
        }
        const factor = bytesToBigInt(data.factor);
        const composite = BigInt(decimal);
        if (factor <= 1n || factor >= composite || composite % factor !== 0n) {
          fail(new Error("coordinator returned an invalid factor"));
          return;
        }
        succeed(factor);
      } else {
        fail(new Error(`unknown coordinator response: ${String(data.type)}`));
      }
    };
    coord.onerror = (event) => {
      event.preventDefault?.();
      fail(new Error(event.message || "coordinator worker crashed"));
    };
    coord.onmessageerror = () => fail(new Error("coordinator returned an invalid message"));

    workers.forEach((w, workerIndex) => {
      w.onmessage = ({ data }) => {
        // Every run response, including errors, is generation-scoped. Old jobs
        // may finish after a successful factor was already returned.
        if (data?.gen !== myGen) return;
        noteActivity();
        if (data.type === "error") {
          fail(new Error(data.error || `sieve worker ${workerIndex + 1} failed`));
          return;
        }
        if (finished) return;
        if (data.type === "prepared") {
          if (workerPrepared[workerIndex]) {
            fail(new Error(`sieve worker ${workerIndex + 1} prepared twice`));
            return;
          }
          if (!data.ok) {
            fail(new Error(`sieve worker ${workerIndex + 1} could not build a sieve`));
            return;
          }
          workerPrepared[workerIndex] = true;
          preparedWorkers++;
          dispatch(w, workerIndex);
        } else if (data.type === "relations") {
          if (!finishJob(workerIndex)) return;
          if (data.payload) {
            if (!(data.payload instanceof Uint8Array)) {
              fail(new Error(`sieve worker ${workerIndex + 1} returned invalid relations`));
              return;
            }
            pendingSubmissions++;
            try {
              coord.postMessage(
                { cmd: "submit", payload: data.payload, worker: workerIndex, gen: myGen },
                [data.payload.buffer],
              );
            } catch (error) {
              pendingSubmissions--;
              fail(error);
            }
            return;
          }
          fail(new Error(`sieve worker ${workerIndex + 1} could not serialize relations`));
        } else {
          fail(new Error(`unknown sieve-worker response: ${String(data.type)}`));
        }
      };
      w.onerror = (event) => {
        event.preventDefault?.();
        fail(new Error(event.message || `sieve worker ${workerIndex + 1} crashed`));
      };
      w.onmessageerror = () =>
        fail(new Error(`sieve worker ${workerIndex + 1} returned an invalid message`));
    });
    document.addEventListener("visibilitychange", onVisibilityChange);
    armStallWatchdog();
    try {
      coord.postMessage({ cmd: "new", n: decimal, gen: myGen });
    } catch (error) {
      fail(error);
    }
  });
}

// How many rho workers race a deep search, and how many polynomial constants each one walks.
// Independent walks over the same modulus collide in about `1.2·sqrt(p)/sqrt(T)` iterations for `T`
// of them, so the pool buys a `sqrt(T)` speedup; eight is where that curve has flattened enough
// that further workers cost more in startup and memory than they return.
const DEEP_RHO_WORKERS = Math.max(1, Math.min(8, nWorkers));
const DEEP_RHO_CONSTANTS_PER_WORKER = 4;

// Per-worker iteration budget for a deep search. Flat, and much smaller than it was, because ECM
// runs after rho now and is cheaper than rho for every factor size past the handover: `T` racing
// walks collide in about `1.2*sqrt(p)/sqrt(T)` iterations each, so eight workers at this budget
// reach a smallest factor near 2^46 — the same place `engine::wide_rho_budget` stops — and
// everything beyond it belongs to curves.
//
// Measured under Node against the scalar wasm module on a 512-bit modulus: 1.086M iterations/s, so
// this is 3.9 s per worker, against 0.88 s for one curve at B1=50,000. The old 2^25 budget cost
// 30.9 s per worker and reached 2^52.5 — and the same wall time buys 35 curves per worker, which
// is most of a schedule that reaches 25-digit factors. Rho stays because it is the cheaper stage
// below the handover, not because a deeper rho was ever going to reach further than curves do.
const DEEP_RHO_BUDGET = 4 << 20;

// The fallback budget, for a runtime with no wasm module or no Workers. It runs on the main thread
// and so has to stay small enough to stay sliceable, and it is one walk rather than eight, so it
// keeps a larger iteration count than the pool does for a shallower reach.
const mainThreadRhoBudget = (bits) => (bits <= 512 ? 8 << 20 : 4 << 20);

// Deep Pollard-Brent in wasm across a pool of workers, racing disjoint polynomial constants.
// Resolves to a factor or null; cancellation is termination, which is why nothing here polls.
async function deepPollard(c, bits, report) {
  if (!wasmModule) return deepPollardOnMainThread(c, bits, report);
  const budget = DEEP_RHO_BUDGET;
  let pool;
  try {
    pool = Array.from(
      { length: DEEP_RHO_WORKERS },
      () => new Worker(new URL("./rho-worker.js", import.meta.url), { type: "module" }),
    );
  } catch {
    return deepPollardOnMainThread(c, bits, report);
  }
  const started = performance.now();
  const tickReport = () =>
    report({
      phase: "deepPollard",
      n: c,
      elapsedSeconds: (performance.now() - started) / 1000,
      workers: pool.length,
      budget,
    });
  tickReport();
  const ticker = setInterval(tickReport, 500);
  try {
    return await new Promise((resolve) => {
      let outstanding = pool.length;
      // A worker that fails — no module, a rejected modulus, a runtime without wasm — is not a
      // reason to fail the run: the others keep searching, and an all-failed pool resolves null,
      // which is the same answer an exhausted budget gives.
      const retire = () => {
        outstanding -= 1;
        if (outstanding <= 0) resolve(null);
      };
      pool.forEach((worker, index) => {
        worker.onmessage = ({ data }) => {
          if (data?.type === "ready") {
            worker.postMessage({
              cmd: "search",
              gen: 1,
              n: c.toString(),
              budget,
              first: 1 + index * DEEP_RHO_CONSTANTS_PER_WORKER,
              count: DEEP_RHO_CONSTANTS_PER_WORKER,
            });
            return;
          }
          if (data?.type === "done") {
            if (typeof data.factor === "string") resolve(BigInt(data.factor));
            else retire();
            return;
          }
          retire();
        };
        worker.onerror = retire;
        worker.postMessage({ cmd: "init", module: wasmModule });
      });
    });
  } finally {
    clearInterval(ticker);
    for (const worker of pool) worker.terminate();
  }
}

// Curves per worker for a composite the sieve has refused. The engine's own schedule is the source
// of the bounds — `qs_ecm` takes zero for "use the default for this width" — so only the split
// across workers is decided here.
const ECM_WORKERS = Math.max(1, Math.min(8, nWorkers));
const ECM_CURVES_PER_WORKER = 64;

// The elliptic curve method across a pool of workers, for a composite nothing else can take.
// Resolves to a factor or null; cancellation is termination.
async function ecmSearch(c, committed, report) {
  if (!wasmModule) return null;
  let pool;
  try {
    pool = Array.from(
      { length: ECM_WORKERS },
      () => new Worker(new URL("./ecm-worker.js", import.meta.url), { type: "module" }),
    );
  } catch {
    return null;
  }
  const started = performance.now();
  const tickReport = () =>
    report({
      phase: "ecm",
      n: c,
      elapsedSeconds: (performance.now() - started) / 1000,
      workers: pool.length,
      curves: pool.length * ECM_CURVES_PER_WORKER,
    });
  tickReport();
  const ticker = setInterval(tickReport, 500);
  try {
    return await new Promise((resolve) => {
      let outstanding = pool.length;
      const retire = () => {
        outstanding -= 1;
        if (outstanding <= 0) resolve(null);
      };
      pool.forEach((worker, index) => {
        worker.onmessage = ({ data }) => {
          if (data?.type === "ready") {
            worker.postMessage({
              cmd: "search",
              gen: 1,
              n: c.toString(),
              b1: 0,
              b2: 0,
              curves: ECM_CURVES_PER_WORKER,
              committed,
              // Disjoint σ stretches, so no two workers walk the same curve.
              seed: 1 + index * ECM_CURVES_PER_WORKER,
            });
            return;
          }
          if (data?.type === "done") {
            if (typeof data.factor === "string") resolve(BigInt(data.factor));
            else retire();
            return;
          }
          retire();
        };
        worker.onerror = retire;
        worker.postMessage({ cmd: "init", module: wasmModule });
      });
    });
  } finally {
    clearInterval(ticker);
    for (const worker of pool) worker.terminate();
  }
}

async function deepPollardOnMainThread(c, bits, report) {
  const budget = mainThreadRhoBudget(bits);
  const run = pollardBrentSliced(c, budget);
  let step = run.next();
  while (!step.done) {
    report({ phase: "deepPollard", n: c, steps: step.value, budget, workers: 0 });
    await tick();
    step = run.next();
  }
  return step.value;
}

// Width from which a cofactor already known to be unbalanced is worth a deep rho rather than a
// sieve, matching `engine::DEEP_RHO_MIN_BITS`. Below it the sieve returns in seconds whatever the
// input's shape, so a long rho could only lose.
const DEEP_RHO_MIN_BITS = 257;

// Stack entries carry whether the value is a cofactor of a composite rho already split. That is the
// evidence that it is not a balanced semiprime — the one shape the cheap opening peel is sized for.
async function factorize(N, report) {
  const primes = [];
  const stack = [[N, false]];
  while (stack.length) {
    let [c, afterSplit] = stack.pop();
    report({ phase: "trial", n: c });
    await tick();
    c = trialDivide(c, primes);
    if (c === 1n) continue;
    report({ phase: "primality", n: c });
    await tick();
    if (isPrime(c)) {
      primes.push(c);
      continue;
    }
    const pp = perfectPower(c);
    if (pp) {
      for (let i = 0; i < pp.k; i++) stack.push([pp.base, afterSplit]);
      continue;
    }
    // Pollard-Brent is a cheap opportunistic peel, not the primary tool at any size: it costs
    // O(sqrt p) in the smallest factor while the sieve costs by the size of `c`, so it only wins
    // where `c` is unbalanced. This used to spend a 2^21 budget below 84 bits on the theory that
    // rho owned that range. Measured here in node (BigInt, single-threaded), 2^21 against 2^15:
    // an 80-bit balanced semiprime 825 ms vs 44 ms, an 85-bit one 724 ms vs 44 ms, while the
    // unbalanced 127-bit case splits in 0.3 ms either way — and the sieve handles those sizes in
    // milliseconds. The large budget was up to 825 ms of blocked main thread for nothing.
    report({ phase: "pollard", n: c });
    await tick();
    const d = pollardBrent(c, 1 << 15);
    if (d && d > 1n && d < c) {
      stack.push([d, true], [c / d, true]);
      continue;
    }
    // Everything cheap has been tried and `c` is a hard composite. This is the only place a width
    // limit applies: trial division, the primality test, perfect powers, and Pollard-Brent above
    // all ran without one, so a wide number built from small factors never reaches here.
    const compositeBits = bitLength(c);
    // Two cases make the cheap peel above the wrong thing to have stopped at. Either the sieve will
    // refuse this composite outright, so that peel was not an opening move but the whole attempt;
    // or `c` split under rho on the way here, which proves it is not the balanced semiprime the
    // peel is sized for — a wide product of middling primes otherwise peels down to just under the
    // sieve's ceiling and then stalls there for weeks. `engine::rho_budget` makes both calls
    // natively; the numbers differ because BigInt runs rho at roughly an eighth of the native rate,
    // not because the policy does.
    const refused = compositeBits > maxSiqsBits;
    const unbalanced = refused || afterSplit;
    if (refused || (afterSplit && compositeBits >= DEEP_RHO_MIN_BITS)) {
      const deep = await deepPollard(c, compositeBits, report);
      if (deep && deep > 1n && deep < c) {
        stack.push([deep, true], [c / deep, true]);
        continue;
      }
    }
    // Curves run wherever the balanced premise has already been disproved, which is the same rule
    // `engine::factor_node` applies: above the ceiling there is no sieve to fall through to, and a
    // composite that split earlier has a small factor and may well have a medium one — which is
    // ECM's range and something the sieve would grind through by the size of `c` instead. A
    // balanced semiprime inside the sieve's range reaches neither branch, so it never pays.
    //
    // `committed` sizes the run rather than deciding it: full curve schedule where ECM is what the
    // composite is counting on, the cheap one where a sieve run is going to happen anyway and
    // finish in seconds.
    if (unbalanced) {
      const committed = refused || compositeBits >= DEEP_RHO_MIN_BITS;
      const curved = await ecmSearch(c, committed, report);
      if (curved && curved > 1n && curved < c) {
        stack.push([curved, true], [c / curved, true]);
        continue;
      }
    }
    if (refused) {
      throw new Error(
        `this number needs the quadratic sieve on a ${compositeBits}-bit composite, ` +
          `above the ${maxSiqsBits}-bit limit — numbers of any size still factor when ` +
          `their factors are small`,
      );
    }
    const factor = await siqsParallel(c.toString(), compositeBits, report);
    if (factor <= 1n || factor >= c || c % factor !== 0n) {
      throw new Error("quadratic sieve returned an invalid factor");
    }
    stack.push([factor, afterSplit], [c / factor, afterSplit]);
  }
  return groupFactors(primes);
}

const tick = () => new Promise((r) => setTimeout(r, 0));


const PHASE_TEXT = {
  trial: (s) => `Trial division on a ${digits(s.n)}-digit number…`,
  primality: (s) => `Miller–Rabin primality test (${digits(s.n)} digits)…`,
  pollard: (s) => `Pollard's rho on a ${digits(s.n)}-digit number…`,
  ecm: (s) =>
    `Elliptic curve method on a ${digits(s.n)}-digit composite — ${s.curves} curves across ` +
    `${s.workers} wasm worker${s.workers === 1 ? "" : "s"}, ${formatDuration(s.elapsedSeconds)} elapsed…`,
  deepPollard: (s) =>
    s.workers
      ? `Deep Pollard's rho on a ${digits(s.n)}-digit composite across ${s.workers} ` +
        `wasm worker${s.workers === 1 ? "" : "s"} — ${formatDuration(s.elapsedSeconds)} elapsed…`
      : `Deep Pollard's rho on a ${digits(s.n)}-digit composite: ` +
        `${Math.round((100 * s.steps) / s.budget)}% of the iteration budget…`,
  sieving: (s) => {
    const progress =
      `Quadratic sieve: ${s.relations}/${s.target} relations across ${nWorkers} workers…`;
    if (s.bits <= 256) return progress;
    const elapsed = `elapsed ${formatDuration(s.elapsedSeconds)}`;
    const eta =
      Number.isFinite(s.etaSeconds) && s.etaSeconds >= 0
        ? `ETA ≈ ${formatDuration(s.etaSeconds)}`
        : "ETA calculating…";
    return `${progress} ${elapsed} · ${eta}`;
  },
  linalg: () => `Linear algebra over GF(2) — extracting a factor…`,
};
const digits = (n) => n.toString().length;
const normalizeNumberText = (text) =>
  text
    .replace(/[０-９]/gu, (digit) => String(digit.codePointAt(0) - 0xff10))
    .replace(/[\p{White_Space}\uFEFF]/gu, "");
const formatDuration = (seconds) => {
  const rounded = Math.max(0, Math.round(seconds));
  if (rounded < 60) return `${rounded}s`;
  const minutes = Math.floor(rounded / 60);
  const remainder = rounded % 60;
  return `${minutes}m ${String(remainder).padStart(2, "0")}s`;
};

function render(grouped, original, seconds) {
  const plain = grouped
    .map(({ prime, exponent }) => (exponent === 1 ? `${prime}` : `${prime}^${exponent}`))
    .join(" * ");
  let product = 1n;
  for (const { prime, exponent } of grouped) product *= prime ** BigInt(exponent);
  const verified = product === original;
  els.result.innerHTML = "";

  // Each factor is shown with its own bit length beneath it, joined by "·".
  const big = document.createElement("div");
  big.className = "factors";
  if (!grouped.length) {
    big.textContent = "1";
  } else {
    grouped.forEach(({ prime, exponent }, i) => {
      if (i) {
        const sep = document.createElement("span");
        sep.className = "sep";
        const sepDot = document.createElement('span');
        sepDot.textContent = "·";
        const sepSpace = document.createElement('span');
        sepSpace.textContent = "\u00a0";
        sep.append(sepDot, sepSpace);
        big.append(sep);
      }
      const factor = document.createElement("span");
      factor.className = "factor";
      const value = document.createElement("span");
      value.className = "value";
      value.textContent =  `${prime}`;
      if (exponent > 1) {
        const exp = document.createElement('span');
        exp.textContent = '^';
        exp.classList.add('exp');
        const expNumber = document.createElement('sup');
        expNumber.textContent = `${exponent}`;
        value.append(exp, expNumber);
      }
      const bits = document.createElement("span");
      bits.className = "bits";
      bits.textContent = `${bitLength(prime)} bits`;
      factor.append(value, bits);
      big.append(factor);
    });
  }

  const meta = document.createElement("div");
  meta.className = "meta";
  meta.textContent =
    `${grouped.length} distinct prime${grouped.length === 1 ? "" : "s"} · ` +
    `${bitLength(original)}-bit input · ` +
    `${verified ? "✓ verified" : "✗ VERIFICATION FAILED"} · ` +
    `${seconds.toFixed(seconds < 10 ? 2 : 1)} s`;
  const copy = document.createElement("code");
  copy.className = "plain";
  copy.textContent = plain || "1";
  els.result.append(big, meta, copy);
  els.result.classList.toggle("bad", !verified);
}

// Live "N digits · M bits" readout for whatever is currently in the input box.
function updateInputInfo() {
  const text = normalizeNumberText(els.input.value);
  const significant = text.replace(/^0+/u, "") || "0";
  if (/^\d+$/.test(text) && significant.length > MAX_DECIMAL_DIGITS) {
    els.inputInfo.textContent =
      `${significant.length} significant digits · exceeds the ${MAX_INPUT_BITS}-bit limit`;
  } else if (/^\d+$/.test(text) && BigInt(text) > 0n) {
    const N = BigInt(text);
    const bits = bitLength(N);
    els.inputInfo.textContent =
      `${text.length} digit${text.length === 1 ? "" : "s"} · ${bits} bits` +
      (bits > MAX_INPUT_BITS ? ` · limit ${MAX_INPUT_BITS}` : "");
  } else {
    els.inputInfo.textContent = "";
  }
}

function resizeNumberInput() {
  // The mirror participates in layout while the textarea overlays it. Updating
  // the mirror, rather than the control's value, cannot disturb IME state.
  els.inputMirror.textContent = `${els.input.value}\u200b`;
}

function normalizeNumberInput() {
  const normalized = normalizeNumberText(els.input.value);
  if (normalized !== els.input.value) els.input.value = normalized;
  resizeNumberInput();
  updateInputInfo();
  return normalized;
}

function insertAtNumberSelection(text) {
  els.input.setRangeText(text, els.input.selectionStart, els.input.selectionEnd, "end");
  resizeNumberInput();
  updateInputInfo();
}

async function run() {
  const text = normalizeNumberInput();
  if (!/^\d+$/.test(text)) {
    els.status.textContent = "Enter a positive whole number.";
    return;
  }
  const significant = text.replace(/^0+/u, "") || "0";
  if (significant.length > MAX_DECIMAL_DIGITS) {
    els.status.textContent = `Enter a number no wider than ${MAX_INPUT_BITS} bits.`;
    return;
  }
  const N = BigInt(text);
  if (N < 1n) {
    els.status.textContent = "Enter a positive whole number.";
    return;
  }
  if (bitLength(N) > MAX_INPUT_BITS) {
    els.status.textContent = `Enter a number no wider than ${MAX_INPUT_BITS} bits.`;
    return;
  }
  els.go.disabled = true;
  els.result.innerHTML = "";
  els.result.classList.remove("bad");
  els.meter.classList.add("busy");
  setBar(0, true);
  const t0 = performance.now();
  const report = (s) => {
    els.status.textContent = (PHASE_TEXT[s.phase] || (() => s.phase))(s);
    if (s.phase === "sieving" && s.target) setBar(s.relations / s.target, false);
    else setBar(0, true);
  };
  try {
    if (N === 1n) {
      render([], 1n, 0);
      els.status.textContent = "1 has no prime factors.";
    } else {
      const grouped = await factorize(N, report);
      render(grouped, N, (performance.now() - t0) / 1000);
      els.status.textContent = "Done.";
    }
  } catch (error) {
    const message = String(error?.message || error);
    els.status.textContent = `Error: ${message} Resetting workers…`;
    try {
      await restartRuntime();
      els.status.textContent = `Error: ${message} Worker runtime was reset.`;
    } catch (restartError) {
      els.status.textContent =
        `Error: ${message} Worker reset failed: ` +
        String(restartError?.message || restartError);
    }
  } finally {
    els.meter.classList.remove("busy");
    setBar(0, false);
    els.go.disabled = !runtimeReady;
  }
}

function setBar(fraction, indeterminate) {
  els.meter.classList.toggle("indeterminate", indeterminate);
  els.bar.style.width = indeterminate ? "100%" : `${Math.min(100, Math.max(0, fraction * 100)).toFixed(1)}%`;
}

els.go.addEventListener("click", run);
els.input.addEventListener("keydown", (e) => {
  // Enter confirms many IME candidates. Never intercept it while composition
  // is active (keyCode 229 covers older engines that omit isComposing).
  if (e.isComposing || e.keyCode === 229) return;
  if (e.key === "Enter") {
    e.preventDefault();
    if (!els.go.disabled) run();
  }
});
els.input.addEventListener("beforeinput", (e) => {
  if (e.isComposing) return;
  const lineAction = e.inputType === "insertLineBreak" || e.inputType === "insertParagraph";
  const hasLine = typeof e.data === "string" && /[\n\r\u2028\u2029]/u.test(e.data);
  if (!lineAction && !hasLine) return;
  e.preventDefault();
  if (hasLine) insertAtNumberSelection(e.data.replace(/[\n\r\u2028\u2029]/gu, ""));
});
els.input.addEventListener("paste", (e) => {
  const pasted = e.clipboardData?.getData("text");
  if (pasted == null || !/[\n\r\u2028\u2029]/u.test(pasted)) return;
  e.preventDefault();
  insertAtNumberSelection(pasted.replace(/[\n\r\u2028\u2029]/gu, ""));
});
els.input.addEventListener("input", () => {
  resizeNumberInput();
  updateInputInfo();
});
els.input.addEventListener("blur", normalizeNumberInput);

// RSA-style semiprime generator (128–384 bits, in steps of 16).
els.rsaBits.addEventListener("input", () => {
  els.rsaBitsLabel.textContent = `${els.rsaBits.value} bits`;
});
els.rsaGen.addEventListener("click", () => {
  const bits = Number(els.rsaBits.value);
  els.rsaGen.disabled = true;
  els.rsaGen.textContent = "Generating…";
  // Yield one frame so the disabled/label state paints before the (synchronous,
  // but brief) prime search runs.
  requestAnimationFrame(() => {
    try {
      els.input.value = rsaNumber(bits).toString();
      resizeNumberInput();
      updateInputInfo();
      els.input.focus();
    } catch (e) {
      els.status.textContent = "Generator error: " + (e?.message || e);
    } finally {
      els.rsaGen.disabled = false;
      els.rsaGen.textContent = "Generate";
    }
  });
});

els.go.disabled = true;
els.status.textContent = "Loading WebAssembly…";
resizeNumberInput();
boot().catch((e) => {
  shutdownRuntime();
  els.status.textContent = "Failed to load: " + (e?.message || e);
});
