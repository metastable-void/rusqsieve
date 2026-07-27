# Browser SIQS readiness and performance roadmap

## Current technical position

The SIQS kernel is credible and performance-oriented. The most important
ingredients are present and implemented as real algorithms rather than labels:

- Knuth–Schroeppel multiplier selection;
- target-fitted squarefree `A` and self-initializing Gray-code polynomials;
- byte logarithmic sieving with translated/sorted roots;
- word-at-a-time survivor discovery;
- multiply-shift factor-base gating;
- single-large-prime graph combination;
- deterministic multi-worker family merging;
- structured sparse filtering with exact provenance;
- bounded scalar/M4RI residual elimination;
- original-matrix dependency verification;
- scoped Wasm SIMD for XOR-heavy work.

The audit's Node/V8 run of the documented 192-bit input, using eight workers
and the SIMD module, returned the verified factorization in:

```text
sieve=0.749 s, finish=0.189 s, wall=0.938 s
relations=4829/4822, families=216
```

This is one audit run, not a replacement benchmark series. It is consistent
with a fast implementation, though slower than the README's recorded 0.72 s.
The difference is not statistically interpretable without repeated,
interleaved runs and load controls.

## Why it is not yet world-class browser-ready

World-class readiness is broader than a fast steady-state kernel:

1. The primary raw Wasm module must build and execute; it does. The bundled
   frontend must contain all referenced glue; the omission found during this
   audit has now been remediated.
2. Every async path must terminate on success, error, timeout, cancellation, or
   exhaustion.
3. A heuristic relation target must support recovery, not destructive failure.
4. Performance must include cold load, compilation, startup, preprocessing,
   sieve, linear algebra, memory, and cleanup.
5. Results must cover the actual target devices and competitors.
6. Release artifacts, not only source-tree `docs/`, must be tested.

## Recommended sequence

Post-audit work completed Gate 0's core items: bounded job accounting,
generation-safe errors, boot/job/run timeouts, runtime reset, 512-bit input
enforcement, strict packet framing, and retryable extraction. The focused
failure test is now part of CI. The lists below remain as regression
requirements and as guidance for broader browser coverage.

### Gate 0 — Make the browser product reliable

- Retain an extracted-archive reference audit so the remediated glue omission
  cannot recur.
- Introduce a typed message protocol with generation on every per-run response.
- Add Worker error/message-error handling and timeouts.
- Add explicit cancellation and terminate or recycle failed workers.
- Centralize active-job accounting and family-budget exhaustion.
- Enforce the supported input range.
- Preserve coordinator state and resume after `NoFactor`.

Acceptance tests should inject:

- coordinator and worker initialization failures;
- stale success and stale error messages;
- out-of-order and duplicate results;
- a null/oversized/truncated relation packet;
- a family budget of two jobs;
- a first extraction that requests more relations;
- worker termination during prepare and sieve;
- two consecutive successful runs with old work still in flight.

### Gate 1 — Establish honest end-to-end baselines

Record, for every run:

- commit and exact artifact hash;
- browser/engine/OS/device and power state;
- logical workers and physical cores;
- cold module fetch, compile, and worker startup;
- preprocessing, first relation, sieve completion, filter/LA/extraction;
- relations, families, survivors if available, and peak memory;
- verified factors and full wall time;
- median, p10/p90, and all raw samples.

Use at least:

- Chromium, Firefox, and Safari/WebKit;
- x86-64 desktop, Apple silicon, midrange Android ARM, and iPhone/iPad;
- 1, 2, 4, 8, and device-appropriate maximum workers;
- cold and warm runs.

Remove the mobile claim until this matrix has a reproducible mobile row.

### Gate 2 — Tune the whole browser pipeline

High-value near-term work:

- move BigInt preprocessing off the main thread;
- use a prime table for trial division;
- check only prime exponents for perfect powers;
- include preprocessing in the sub-second tier benchmark;
- choose worker count by input tier and device evidence, not only
  `hardwareConcurrency`;
- lazily create workers or maintain a smaller warm pool for easy inputs;
- measure context-rebuild cost, packet bytes, transfer cost, overshoot, and
  per-worker memory;
- expose profiling counters in a development-only ABI rather than parsing UI
  text.

The existing two-family batch and 48-worker ceiling are defensible on the
reference host, but not universal. Mobile devices may prefer fewer workers due
to thermal limits and memory bandwidth.

### Gate 3 — Generalize performance evidence

Keep a tuning set and a separate holdout set. At every bit tier, include:

- balanced semiprimes;
- several multipliers and smoothness profiles;
- 32/48/64/80/96-bit smallest factors;
- three-prime composites;
- repeated prime powers;
- primes and near-prime composites.

Run current Alpertron, CrypTool/Msieve-Wasm, YAFU/Yaffle-derived browser tools,
and other maintained candidates in the same browser and device. Publish
failures and unsupported cases, not only wins.

### Gate 4 — Consider deeper algorithm changes only when profiles justify them

- Keep ECM out of the default balanced-RSA artifact. If general factorization
  is desired, ship it as a separate opt-in artifact and measure its download,
  compile, startup, and code-cache effects.
- Investigate bucket/cache-blocked sieving on small-L2 mobile hardware; the
  existing Xeon result cannot answer that question.
- Move to block Lanczos/Wiedemann only when residual matrices make current M4RI
  materially dominant and a verified reference implementation is available.
- Revisit double-large-prime variants only with a co-designed survivor filter;
  the recorded isolated DLP experiments correctly showed no win.

## Suggested release performance gates

The exact thresholds should be set from a clean baseline, but the gate shape
should be:

- no factorization or product-verification failures on the fixed and holdout
  corpora;
- no unresolved Promise or live Worker after test completion;
- scalar and SIMD parity on every browser that supports each artifact;
- raw scalar and SIMD Wasm modules execute, and any advertised demo boots after
  extraction;
- Wasm bytes, cold startup, 192/224/256 wall time, and peak memory do not
  regress beyond an explicitly approved tolerance;
- mobile thermal-repeat runs do not collapse after the first sample;
- competitor language is mechanically checked against the evidence document.

## Bottom line

The right next investment is reliability and measurement, not another isolated
inner-loop optimization. The kernel is already strong enough to justify serious
competitive benchmarking. The surrounding browser product and evidence chain
are what currently prevent a world-class readiness claim.
