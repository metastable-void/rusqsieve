// Main thread: UI + coordinator. Peels easy factors with BigInt number theory and
// hands hard composites to a pool of wasm Web Workers running the quadratic sieve,
// with the pool sized to navigator.hardwareConcurrency.
import { instantiate, loadModule, putString, putBytes, takePacket, bytesToBigInt } from "./abi.js";
import { trialDivide, isPrime, perfectPower, pollardBrent, groupFactors, rsaNumber, bitLength } from "./numtheory.js";

const WASM_URL = new URL("./rusqsieve.wasm", import.meta.url);
const BATCH = 4; // polynomial families dispatched per sieve job
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

let coord = null; // coordinator wasm instance (main thread)
let workers = []; // sieve worker pool
let gen = 0; // generation token so stale worker messages are ignored
const nWorkers = Math.max(1, navigator.hardwareConcurrency || 4);

async function boot() {
  const module = await loadModule(WASM_URL);
  coord = await instantiate(module);
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
  els.workers.textContent = `${nWorkers} worker${nWorkers === 1 ? "" : "s"} · ABI v${coord.qs_abi_version()}`;
  els.go.disabled = false;
  els.status.textContent = "Ready.";
}

// Parallel quadratic sieve for one hard composite; resolves to a nontrivial factor.
function siqsParallel(decimal, report) {
  return new Promise((resolve, reject) => {
    const myGen = ++gen;
    const s = putString(coord, decimal);
    const session = coord.qs_coord_new(s.ptr, s.len);
    coord.qs_dealloc(s.ptr, s.len, 1);
    if (!session) return reject(new Error("could not build a sieve for this number"));
    const target = coord.qs_coord_target(session);
    let relations = 0;
    let nextFamily = 0;
    let finished = false;

    const dispatch = (w) => {
      if (nextFamily > MAX_FAMILIES) return;
      const family = nextFamily;
      nextFamily += BATCH;
      w.postMessage({ cmd: "sieve", family, count: BATCH, gen: myGen });
    };
    const finish = () => {
      finished = true;
      report({ phase: "linalg" });
      const handle = coord.qs_coord_extract(session);
      const payload = takePacket(coord, handle);
      coord.qs_coord_free(session);
      if (!payload) return reject(new Error("linear algebra found no factor"));
      resolve(bytesToBigInt(payload));
    };

    for (const w of workers) {
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
            const b = putBytes(coord, data.payload);
            relations = coord.qs_coord_submit(session, b.ptr, b.len);
            coord.qs_dealloc(b.ptr, b.len, 1);
            report({ phase: "sieving", relations, target });
          }
          if (relations >= target) finish();
          else if (nextFamily > MAX_FAMILIES) {
            finished = true;
            reject(new Error("relation budget exhausted"));
          } else dispatch(w);
        }
      };
      w.postMessage({ cmd: "prepare", n: decimal, gen: myGen });
    }
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
    // Pollard's rho is only worthwhile where it can actually finish: it is the
    // primary tool below the sieve's viable range (~80 bits), and a quick cheap
    // peel above it (to catch a small factor trial division missed). It cannot
    // split a balanced large semiprime — that is the sieve's job.
    const bits = c.toString(2).length;
    if (bits <= 84) {
      report({ phase: "pollard", n: c });
      await tick();
      const d = pollardBrent(c, 1 << 21);
      if (d && d > 1n && d < c) {
        stack.push(d, c / d);
        continue;
      }
    } else {
      const d = pollardBrent(c, 1 << 15);
      if (d && d > 1n && d < c) {
        stack.push(d, c / d);
        continue;
      }
    }
    const factor = await siqsParallel(c.toString(), report);
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
  sieving: (s) => `Quadratic sieve: ${s.relations}/${s.target} relations across ${nWorkers} workers…`,
  linalg: () => `Linear algebra over GF(2) — extracting a factor…`,
};
const digits = (n) => n.toString().length;

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

// RSA-style semiprime generator (128–384 bits, in steps of 32).
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
