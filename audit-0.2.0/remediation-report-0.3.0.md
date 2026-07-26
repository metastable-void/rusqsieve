# rusqsieve 0.2 audit remediation report

Status date: 2026-07-26. Working version: 0.3.0.

## Release and Phase 1–2 confirmation

- Confirmed locally: clean starting tree; `main`, `origin/main`, and signed tag
  `v0.2.1` all resolved to `b1444f204f6b453961eaa2d38ee30e820d11f503`.
- Confirmed externally: GitHub release `v0.2.1` is published (not draft or
  prerelease) with eight uploaded target tarballs; crates.io version 0.2.1 is
  published and not yanked.
- Phase 1 and Phase 2 are recorded in the 0.2.1 changelog with inputs,
  measurements, methodology, negative results, and unmet numeric acceptance
  targets. The regression, corpus, primality, cancellation, scheduler,
  arithmetic, sieve, relation, and browser tests all passed before starting
  0.3 work.

Per main-brief item:

| Item | Status | Result |
|---|---|---|
| 1.1 dead zone | Done in 0.2.1 | Size sweep, direct SIQS polynomial test, and 309-entry variable-arity corpus. |
| 1.2 coefficient famine | Done in 0.2.1/0.3.0 | Early internal diagnosis in 0.2.1; public `PolynomialSelection` error in 0.3.0. |
| 1.3 big-integer rho | Partial by evidence | Bounded cancellable Brent landed. The requested “rho handles balanced 65–100 bits” premise measured as a 19–216× regression and was rejected; balanced inputs correctly fall through to SIQS. |
| 1.4 cancellation/panics | Done in 0.2.1 | u64 cancellation and first worker panic propagation. |
| 1.5 Baillie–PSW | Done in 0.2.1 | Rust and browser paths, Lucas execution hook, extended pseudoprime tests, measured +12.6% on primality-only calls and 0.007% of a 256-bit run. |
| 2.1 resieving | Skipped after implementation/measurement | Net loss at every tested cutoff; full figures are in `CHANGELOG.md`. |
| 2.2 arithmetic wins | Done/partial | All measured material 0.2.1 wins landed; 0.3.0 additionally adds real fixed-limb Montgomery REDC for native big rho. The survivor power-list linear scan remains intentionally unchanged at roughly 15 elements. |
| 2.3 threshold retune | Partial by numeric criterion | Correct model and per-tier offsets landed; measured gains were about 4%/2%, not the requested 10%. |
| 2.4 large-prime bound | Partial by numeric criterion | Correct independent bound landed; memory fell 2.2×, not 10×. |
| 2.5 u64 rho | Done in 0.2.1 | Real Montgomery/Brent batching, 18.4–20.8× on measured cofactors. |
| 2.6 double large primes | Skipped after measurement | No wall-time win at matched threshold; genuine-double threshold exceeded 300 seconds. |
| 2.7 family setup | Done where useful | Binary xgcd landed; unhelpful remainder precomputation rejected. |
| 2.8 `nvar` | Done in 0.2.1 | Cap raised to nine with memory accounting. |
| 2.9 `choose_a` | Done in 0.2.1 | Quality window, deduplication, hoisted pool, checked conversion. |
| 2.10 prepare deduplication | Done in 0.2.1 | One preparation path plus translation regression test. |
| 2.11 linear algebra/browser | Partial | Coordinator Worker, flat sorted filtering structures, allocation reuse, and scalar parity unroll landed. Full CSR compaction, Wasm SIMD back-substitution, and real M4RI did not; implemented M4RI measured 75% slower. |
| 2.12 inner loop | Done as directed | Winning biased scan/wrapping writes retained; unreachable blocked kernel deleted; small-L2 remeasurement remains unavailable. |

## Phase 3

| Item | Status | Result |
|---|---|---|
| 3.1 environment | Done | No library `std::env` reads remain. CLI maps all tuning variables into owned config; profiling is config-driven. |
| 3.2 effective config/path | Done | `FactorConfig` contains only honored controls. Every value through 512 significant bits, including `Natural<P>` for `P > 16`, uses optimized SIQS. |
| 3.3 dead/misnamed code | Done | Removed the fake Montgomery façade and later added a distinct, real REDC context with domain conversions and differential tests. Removed fake xgcd/Lanczos/provenance/matrix config/session/work/reference-QS/metrics surfaces. Consolidated machine-word primality, powmod, mulmod, and xorshift. Kept `is_square` for BPSW. |
| 3.4 features | Done | Removed three inert features and the type-identity-changing width feature. Remaining features gate CLI, Wasm SIMD, or fuzz hooks. |
| 3.5 C/Wasm safety | Done | Unwinding release libraries, 256 explicit-thread cap, unsafe pointer exports, copied Wasm input, ABI/status helpers, enum status return, and C progress/cancellation. |
| 3.6 public hardening | Done except serde | Non-exhaustive enums, accessor errors, must-use queries, named iterators, `IntoIterator`, `TryFrom`, consistent 512-bit bound, compact worker roots. Serde deliberately skipped: the active private wire does not use it and a new public format contract was unjustified. |
| 3.7 module splits | Partial | Engine wire, extraction, SIQS initialization, and shared u64 math were split. Larger kernel/relation and `natural` pure moves are deferred to a behavior-neutral follow-up so this breaking logic change remains reviewable. |
| 3.8 seeded witnesses | Done | ChaCha8, full-width rejection sampling in `[2,n−2]`, no skipped rounds, deterministic cross-target algorithm, tests and documentation. |

## Phase 4

| Item | Status | Result |
|---|---|---|
| 4.1 CI | Done | Musl-native debug/release, feature matrix, fmt, clippy, two Wasm paths, C execution, scheduled slow corpus. |
| 4.2 build failures | Done | CLI test feature-gated; cold `cargo test --release` now links and passes. |
| 4.3 tests/fuzzing | Done | Added requested arithmetic differential families, public edge/recursive/wide tests, seeded witnesses, BPSW hook, corpus tiers, hostile/concurrent C tests, and three cargo-fuzz targets. |
| 4.4 C/Wasm build truth | Done | C smoke compiles/runs; Makefile and release builder both enable target SIMD128. |
| 4.5 documentation truth | Done | Removed reference-QS/DLP/provenance/inert-feature/width/panic/SIMD mismatches. |
| 4.6 documentation quality | Mostly done | Safety/error/panic contracts, architecture map, entry-point timing caveats, public C rustdoc, and invariant comments added. Exhaustive complexity sections for every private kernel are deferred with the remaining module-only split. |

## Addenda

- Addendum A (ECM): ECM remains deliberately absent. A later standalone
  experiment added real Montgomery REDC and retained it after measurement, but
  no ECM budget was introduced: on balanced RSA-like inputs it would remain
  unsuccessful overhead. Any future ECM stage must be separately opt-in or
  default-zero and pass the addendum's unbalanced-input acceptance tests.
  Pollard p−1 was also not inserted as unconditional overhead.
- Addendum B: all B.1/B.2 lying or stub surfaces were removed; all B.3 text
  defects and reference checks were addressed. B.4’s honestly named
  `pollard_u64` and `Forest` behavior was not changed. The web/tools naming
  sweep and its two fixes were already completed in 0.2.1.

## Semver and measurements

- 0.2.1: patch; private correctness/performance changes only.
- 0.3.0: breaking minor; public Rust/C/Wasm removals, additions, enum changes,
  feature removals, supported-range enforcement, and wire ABI 2.
- Native big-rho Montgomery arithmetic measured 3.10× faster for a 256-bit
  modular-square kernel and 2.33 s → 1.33–1.34 s over 200 fresh invocations of
  the fixed unbalanced 224-bit rho case. Interleaved balanced 224-bit runs were
  a wash (5.55–5.58 s before, 5.55–5.62 s after); no SIQS speedup is claimed.
- Real Chromium/Web Worker measurements now cover five deterministic balanced
  cases at 216, 224, 232, 240, 256, and 272 bits. Retained eight-worker SIMD
  tier retunes lower corpus means by 35.9% at 216 bits, 16.8% at 224 bits, 28.0%
  at 232 bits, 7.5% at 240 bits, and 8.0% at 256 bits. At 272 bits, sieving
  remains the frontier; interval and factor-base sweeps, a dense prime-only
  stream, and bit-sliced dependency back-substitution produced no reliable
  end-to-end win and were removed.
- All Phase 2 before/after figures, host details, A/B/A/B method, measured
  rejections, and variance qualifications remain in the 0.2.1 changelog and
  were not re-labeled as new 0.3 measurements.

## Reference drift and additional findings

- The brief’s line numbers were for v0.2.0 and no longer matched after the
  confirmed 0.2.1 work. Semantic targets were resolved by identifier; no
  ambiguous target was guessed.
- The release matrix retains both GNU and musl Linux artifacts, per the
  maintainer's explicit release policy.
- Additional defect found in this phase: Wasm SIMD compiled only after restoring
  a narrowly scoped unsafe allowance for the intrinsic kernel; the broad Wasm
  exemption was not reintroduced.

## Verification

Passing locally:

- `cargo test --locked`
- `cargo test --locked --no-default-features`
- `cargo test --locked --release`
- `cargo clippy --locked --all-targets --all-features -- -D warnings`
- musl-target clippy with the same warning policy
- `make wasm`
- `node tools/browser-arch-check.mjs` (ABI 2, real coordinator + four workers)
- `make c-api-smoke`
- documentation-reference audit
- `cargo check --manifest-path fuzz/Cargo.toml`
