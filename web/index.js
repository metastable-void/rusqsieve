// Main thread: UI + coordinator. Peels easy factors with BigInt number theory and
// hands hard composites to a pool of wasm Web Workers running the quadratic sieve,
// with the pool sized to navigator.hardwareConcurrency.
import { loadModule, bytesToBigInt } from "./abi.js";
import { trialDivide, isPrime, perfectPower, pollardBrent, groupFactors, rsaNumber, bitLength } from "./numtheory.js";

const SIMD_WASM_URL = new URL("./rusqsieve-simd.wasm", import.meta.url);
const SCALAR_WASM_URL = new URL("./rusqsieve.wasm", import.meta.url);
// Small jobs reduce the tail after the relation target is reached. Two
// families was consistently best in Node/V8 from 192 through 256 bits.
const BATCH = 2;
const MAX_FAMILIES = 2_000_000;

const els = {
  input: document.getElementById("input"),
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
let gen = 0; // generation token so stale worker messages are ignored
let wasmFlavor = "scalar";
// Scaling remains positive through 32–48 workers on the 96-thread reference
// host, while 96 workers regress from startup, memory traffic, and job overshoot.
const nWorkers = Math.max(1, Math.min(48, navigator.hardwareConcurrency || 4));

async function boot() {
  let module;
  try {
    module = await loadModule(SIMD_WASM_URL);
    wasmFlavor = "SIMD";
  } catch {
    // Older engines can still use the portable artifact.
    module = await loadModule(SCALAR_WASM_URL);
  }
  coord = new Worker(new URL("./coordinator.js", import.meta.url), { type: "module" });
  const abi = await new Promise((resolve, reject) => {
    coord.onmessage = ({ data }) => {
      if (data.type === "ready") resolve(data.abi);
      else if (data.type === "error") reject(new Error(data.error));
    };
    coord.postMessage({ cmd: "init", module });
  });
  workers = await Promise.all(
    Array.from({ length: nWorkers }, () => {
      const w = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
      return new Promise((resolve) => {
        w.addEventListener("message", function ready({ data }) {
          if (data.type === "ready") {
            w.removeEventListener("message", ready);
            resolve(w);
          }
        });
        w.postMessage({ cmd: "init", module });
      });
    }),
  );
  els.workers.textContent =
    `${nWorkers} worker${nWorkers === 1 ? "" : "s"} · ${wasmFlavor} · ABI v${abi}`;
  els.go.disabled = false;
  els.status.textContent = "Ready.";
}

// Parallel quadratic sieve for one hard composite; resolves to a nontrivial factor.
function siqsParallel(decimal, bits, report) {
  return new Promise((resolve, reject) => {
    const myGen = ++gen;
    const sieveStarted = performance.now();
    let target = 0;
    let relations = 0;
    let nextFamily = 0;
    let finished = false;

    const dispatch = (w) => {
      if (nextFamily > MAX_FAMILIES) return;
      const family = nextFamily;
      nextFamily += BATCH;
      w.postMessage({ cmd: "sieve", family, count: BATCH, gen: myGen });
    };
    coord.onmessage = ({ data }) => {
      if (data.gen !== myGen && data.type !== "error") return;
      if (data.type === "error") {
        if (!finished) {
          finished = true;
          reject(new Error(data.error));
        }
      } else if (data.type === "session") {
        target = data.target;
        for (const w of workers) w.postMessage({ cmd: "prepare", n: decimal, gen: myGen });
      } else if (data.type === "submitted") {
        relations = data.relations;
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
        if (!finished) dispatch(workers[data.worker]);
      } else if (data.type === "linalg") {
        report({ phase: "linalg" });
      } else if (data.type === "factor") {
        finished = true;
        resolve(bytesToBigInt(data.factor));
      }
    };

    workers.forEach((w, workerIndex) => {
      w.onmessage = ({ data }) => {
        // Errors are always surfaced; other messages from an obsolete generation
        // (a worker's in-flight job for a previous composite) are ignored.
        if (data.type === "error") {
          if (!finished) {
            finished = true;
            reject(new Error(data.error));
          }
          return;
        }
        if (finished || data.gen !== myGen) return;
        if (data.type === "prepared") {
          if (!data.ok) {
            finished = true;
            reject(new Error("worker could not build a sieve"));
            return;
          }
          dispatch(w);
        } else if (data.type === "relations") {
          if (data.payload) {
            coord.postMessage(
              { cmd: "submit", payload: data.payload, worker: workerIndex, gen: myGen },
              [data.payload.buffer],
            );
            return;
          }
          if (nextFamily > MAX_FAMILIES) {
            finished = true;
            reject(new Error("relation budget exhausted"));
          } else dispatch(w);
        }
      };
    });
    coord.postMessage({ cmd: "new", n: decimal, gen: myGen });
  });
}

async function factorize(N, report) {
  const primes = [];
  const stack = [N];
  while (stack.length) {
    let c = stack.pop();
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
      for (let i = 0; i < pp.k; i++) stack.push(pp.base);
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
      stack.push(d, c / d);
      continue;
    }
    const factor = await siqsParallel(c.toString(), bitLength(c), report);
    stack.push(factor, c / factor);
  }
  return groupFactors(primes);
}

const tick = () => new Promise((r) => setTimeout(r, 0));
const SUP = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const sup = (n) => String(n).replace(/\d/g, (d) => SUP[d]);

const PHASE_TEXT = {
  trial: (s) => `Trial division on a ${digits(s.n)}-digit number…`,
  primality: (s) => `Miller–Rabin primality test (${digits(s.n)} digits)…`,
  pollard: (s) => `Pollard's rho on a ${digits(s.n)}-digit number…`,
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
        sep.textContent = "·";
        big.append(sep);
      }
      const factor = document.createElement("span");
      factor.className = "factor";
      const value = document.createElement("span");
      value.className = "value";
      value.textContent = exponent === 1 ? `${prime}` : `${prime}${sup(exponent)}`;
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
  const text = els.input.value.trim();
  if (/^\d+$/.test(text) && BigInt(text) > 0n) {
    const N = BigInt(text);
    els.inputInfo.textContent = `${text.length} digit${text.length === 1 ? "" : "s"} · ${bitLength(N)} bits`;
  } else {
    els.inputInfo.textContent = "";
  }
}

async function run() {
  const text = els.input.value.trim();
  if (!/^\d+$/.test(text)) {
    els.status.textContent = "Enter a positive whole number.";
    return;
  }
  const N = BigInt(text);
  if (N < 1n) {
    els.status.textContent = "Enter a positive whole number.";
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
    els.status.textContent = "Error: " + (error?.message || error);
  } finally {
    els.meter.classList.remove("busy");
    setBar(0, false);
    els.go.disabled = false;
  }
}

function setBar(fraction, indeterminate) {
  els.meter.classList.toggle("indeterminate", indeterminate);
  els.bar.style.width = indeterminate ? "100%" : `${Math.min(100, Math.max(0, fraction * 100)).toFixed(1)}%`;
}

els.go.addEventListener("click", run);
els.input.addEventListener("keydown", (e) => {
  if (e.key === "Enter" && !els.go.disabled) run();
});
els.input.addEventListener("input", updateInputInfo);

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
boot().catch((e) => {
  els.status.textContent = "Failed to load: " + (e?.message || e);
});
