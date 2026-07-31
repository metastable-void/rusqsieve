# rusqsieve 0.4 performance audit

This file records the current performance status, the optimizations that remain
relevant, and measured negative results worth preserving. Resolved root causes,
superseded roadmaps, stale line numbers, and pre-0.2 public-API constraints have
been removed.

The goal is to remain faster than FLINT's single-threaded QSieve on the measured
balanced-semiprime range while preserving the independently established
fastest-browser-Wasm-SIQS result for realistic 192–256-bit balanced
semiprimes.

The matching FLINT checkout is `/home/dev/flint` and is the authoritative local
implementation reference. This LXD container is effectively dedicated
bare-metal; load changes normally come from this project, Codex, VS Code, or
Claude Code. Performance work should use factor-verified saved binaries,
interleaved A/B order, and no concurrent build/test load.

## Current verified status

### Native single-thread comparison

Clean saved-binary measurements on the reference host:

| bits | rusqsieve | FLINT QSieve | result |
|-----:|-----------:|-------------:|-------:|
| 160 | 0.35 s | 0.57 s | 1.63× faster |
| 192 | 2.49–2.53 s | 2.78–2.79 s | 1.10–1.12× faster |
| 208 | 8.59 s | 9.89 s | 1.15× faster |
| 224 | 21.43 s | 26.60 s | 1.24× faster |
| 240 | 67.25 s | 77.18 s | 1.15× faster |

After the 0.2 API/cancellation work, the fixed 192-bit case measured 2.43 s
wall/user on one thread, confirming no native regression. FLINT aborts on the
reference host's fixed 256-bit case; rusqsieve remains functional there.

These figures support the README's narrow claim: rusqsieve beats FLINT's
single-threaded QSieve on every measured native tier from 160 through 240 bits.
They do not establish a general integer-factorization record.

### Browser-shaped Wasm comparison

Node 24.15/V8, one coordinator instance, independent worker instances,
serialized relation packets, no shared memory, and scoped SIMD128 matrix/root
kernels. The most recent 2026-07-31 release gate measured:

| bits | 8-worker wall | sieve | filter/LA/extract |
|-----:|--------------:|------:|------------------:|
| 192 | 0.773 s | 0.585 s | 0.188 s |
| 224 | 3.973 s | 3.533 s | 0.440 s |
| 256 | 27.293 s | 25.316 s | 1.978 s |

The fixed 256-bit case measured 22.26 s, 14.71 s, and 13.96 s with 16, 32, and
48 workers. Ninety-six workers regressed from startup, memory traffic, and
relation overshoot, so the browser pool is capped at 48.

All outputs were factor-verified. These are engineering measurements on one
host, not cross-project records. `BENCHMARKING.md` defines the reproducible
corpus and competitor protocol.

### High-digit native validation

The 0.4 release gate also exercised both new large-matrix/high-digit paths:

| input | workers | wall | peak RSS | result |
|---:|---:|---:|---:|---|
| fixed balanced 288-bit fixture | 48 | 37.81 s | 297 MiB | exact factors |
| RSA-100 (330 bits) | 96 | 320.10 s | 1.48 GiB | exact factors |

The RSA-100 run improves on the prior 339.2 s profiled result. It remains
slower than the matched portable YAFU reference at 185.94 s; profiling places
the remaining gap in relation collection rather than sparse linear algebra.

## Optimizations present in 0.4

### Arithmetic and preprocessing

- significant-limb fixed-capacity arithmetic;
- widening multiplication and wide modular reduction;
- normalized limb long division;
- binary GCD;
- cached small-prime table;
- deterministic `u64` Miller–Rabin and Pollard–Brent fast path;
- recursive probable-prime and perfect-power handling.

The former bit-by-bit division, double-and-add modular multiplication, weak rho,
and repeated prime-list construction are resolved and must not be reintroduced.

### SIQS polynomial and sieve path

- deterministic Knuth–Schroeppel multiplier;
- numerically target-fitted squarefree `A`;
- smaller signed CRT `B`;
- Gray-code self-initialization;
- translated, sorted score-position roots;
- paired root-difference stride loop;
- byte scores and word-at-a-time high-bit candidate rejection;
- tuned tiny-prime skipping and threshold margin;
- direct signed `g(x) = Q(x)/A` reconstruction;
- Lemire/Barrett-gated candidate division;
- confirmed-score early termination;
- reusable per-worker scratch;
- portable dense-prefix cache blocking at the measured high-digit gate;
- precomputed contiguous prime and score-weight streams;
- flat candidate-factor scratch and bitset report membership;
- single- and double-large-prime relation graph with allocation-free
  union-by-size cycle combination;
- optional x86-64 SSE2 and Wasm SIMD128 root advancement.

Below 289 bits, Lemire/Barrett root gating removes most big-integer divisions
without adding a second sieve pass. High-digit DLP tiers use the retained
prime-major report resieve, whose fixed-stride scratch and one-bit membership
filter avoid per-report allocation and sparse random-access maps.

### Linear algebra

The current solver performs deterministic sparse elimination through row
weight six, tracks exact original-column provenance, and then selects compact
scalar/M4RI row-echelon below 272 input bits or a true portable 64-way
Montgomery block-Lanczos recurrence from 272 bits upward. It expands at most 64
useful dependencies and verifies each against the original matrix.

At the fixed 256-bit tier, filtering reduces roughly 20,740×20,803 to
10,877×11,343. The compact residual solve reduced extraction from about 5.19 s
to 2.93 s in the recorded native eight-thread run.

The Lanczos implementation applies sparse `Bᵀ(BV)`, handles rank-deficient
blocks with the complete three-term recurrence, retries deterministic starting
blocks, and verifies terminal and lifted dependencies. On the fixed 272-bit
browser matrix it reduced LA/extraction from 5.574 s to 1.721 s. The larger
1.4M factor base plus Lanczos reduced the fixed 288-bit wall time from
420.223 s to 258.512 s during tuning.

### Wasm

- independent coordinator/worker instances;
- deterministic worker context reconstruction;
- two-family worker jobs to limit tail overshoot;
- serialized relation packets;
- SIMD128 restricted to explicit matrix-XOR and SIQS-root kernels;
- automatic scalar fallback;
- no shared memory or Rust Wasm threads.

Whole-program `+simd128`, Binaryen 120 `wasm-opt -O3`, and `-Oz` all regressed
the measured sieve and remain excluded.

## Retained negative results

These experiments were correct but slower on the reference host:

- plain modulo root-gating before precomputed Lemire constants;
- a second stride-based resieve at small/medium factor bases;
- whole-interval unconditional cache blocking on a 1 MiB-L2 Xeon;
- larger sieve intervals without a matching cache strategy;
- branchless scalar root updates;
- `step_by` and slice-iterator score loops;
- relation targets below the safe factor-base-plus-surplus target;
- whole-program Wasm SIMD;
- Binaryen post-optimization;
- `limit-to-512-bits` by itself on scalar x86;
- AVX-512 instead of AVX2 on the Cascade Lake reference host.

Cache blocking is not rejected universally: it remains plausible on small-L2
mobile hardware, where the score interval no longer fits in private cache. It
must be validated on representative hardware rather than inferred from the
Xeon result.

## Current open work

### 1. Same-browser competitor measurements

The highest-priority evidence gap is not another rusqsieve micro-optimization.
Run current Alpertron, Msieve-Wasm/CrypTool, and other maintained browser
factorizers on the same browser, hardware, fixed inputs, warm-up policy, and
worker counts. Publish medians and ranges, not isolated best runs.

Until then, use "plausibly fastest-class browser factorizer for balanced
192–256-bit semiprimes" and "candidate for fastest browser SIQS," not "fastest
general factorizer."

### 2. Multi-input parameter tuning

The factor-base/interval tables contain successful single-input and small-sweep
tuning. Build a deterministic sweeper over several balanced semiprimes per tier
and jointly tune:

- factor-base bound;
- sieve half-width;
- tiny-prime cutoff and slack;
- survivor margin;
- multiplier candidates;
- relation surplus;
- polynomial families per job.

Retain a holdout corpus to detect overfitting.

### 3. Mobile cache behavior

Measure flat versus carried-position blocked versus true bucket sieving on
representative ARM phones/tablets and smaller-cache consumer CPUs. The desktop
gate must remain the default unless an on-target win is demonstrated.

### 4. Linear algebra beyond the current range

True block Lanczos is now the measured solver from 272 bits upward. Further
work should target matrices beyond the current RSA-100 (330-bit) validation:
reduce retry/breakdown frequency, bound peak memory, and compare alternative
sparse layouts or Wiedemann only with exact dependency verification and
retained provenance. Do not move the 272-bit crossover without matched
controls below and above it.

### 5. General-factorization coverage

Unbalanced composites with medium-size factors remain the algorithmic gap.
ECM is acceptable only behind a non-default feature and in a separate
general-purpose Wasm artifact. The default balanced-RSA build must contain no
ECM code or initialization, and the 192/224/256 corpus must remain an A/B gate
for runtime, artifact size, startup, and compilation.

## Non-negotiable regression gates

- SIQS remains a true logarithmic sieve.
- Relation congruences and combined large-prime parity remain valid.
- Every matrix dependency is verified against the original matrix.
- Returned factors reconstruct the input and pass probable-prime testing.
- Portable/browser family merging is deterministic; native completion-order
  ingestion remains duplicate-safe and dependency-verified.
- Native cancellation joins workers and returns `FactorError::Cancelled`.
- The ordinary no-observer Rust path retains no progress callback/timing cost.
- Native C exports remain limited to the documented opaque eight-symbol ABI.
- The balanced-RSA Wasm artifact remains ECM-free.
- Scalar and SIMD Wasm builds both pass.
- Native 192-bit and browser 192/224/256-bit timings, Wasm size, and startup are
  compared before accepting performance-sensitive changes.
