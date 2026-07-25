# Changelog

## 0.2.1 — 2026-07-25

Patch release: correctness fixes and internal-only performance/resource work;
no public Rust signatures changed.

- Fixed factorization of balanced 65–100-bit semiprimes, including
  `18446744400127067027`, by making SIQS `A`-coefficient selection satisfiable
  and adding a bounded, cancellable big-integer Pollard–Brent preprocessing
  stage.
- Distinguished deterministic polynomial-coefficient famine internally,
  stopped retrying it for 100,000 families, made the small-factor path
  cancellable, recovered poisoned worker locks, and preserved the first worker
  panic message.
- Added Baillie–PSW (base-2 Miller–Rabin plus strong Selfridge Lucas) to final
  primality decisions above 64 bits in Rust and browser RSA generation. The
  proven deterministic 64-bit witness path is unchanged.
- Added the independently verified 309-entry, 7–256-bit factorization corpus
  and dead-zone size sweep.
- Replaced discarded-quotient reductions, redundant modular reductions, and
  hot division allocations; cached small primes and segmented factor-base
  generation.
- Added twice-log2 sieve weights with slack derived from the skipped primes,
  recorded sparse-tail resieve hits, a 256× factor-base single-large-prime
  limit, allocation reuse, index-only partial-relation forests, and
  Montgomery/batched-GCD cofactor splitting.
- `find_factor` now uses the shared `prepare` path. The pinned interval
  translation is `sieve_half_width % prime`; the formerly duplicated native
  expression was numerically equivalent only while the interval remained
  positive and unchanged.
- Moved browser relation coordination and GF(2) extraction off the main thread,
  flattened filtering storage, added Four-Russians residual reduction and
  unrolled dependency back-substitution, and guarded oversized dense fallback.
- Double-large-prime collection remains disabled: it did not demonstrate a net
  wall-time win with the bounded large-prime policy.

## 0.2.0 — 2026-07-25

- Replaced the broad implementation-facing Rust API with a documented,
  high-level factorization surface.
- Made SIQS factor bases and relations, sparse-matrix kernels, scheduler
  sessions, worker packets, primality policy, and limb mutation crate-private.
- Encapsulated `FactorConfig`; supported tuning now uses builder methods for
  parallelism and progress cadence.
- Replaced `Parallelism::Exact` with the nonzero
  `Parallelism::Threads`/`Parallelism::threads` interface.
- Made progress snapshots read-only and added documented accessors.
- Clarified `PrimeFactors` cardinality with `distinct_len`, `total_len`, and
  `expanded`; removed direct map extraction.
- Added cooperative cancellation to the native SIQS engine.
- Enabled `deny(missing_docs)` for the complete public Rust API.
- Kept the decimal-string C ABI and opaque result ownership model unchanged.
- Replaced the aspirational implementation specification and historical
  performance audit with documentation of the shipped 0.2 architecture,
  verified results, and remaining work.
- Added complete crates.io metadata and an explicit publish allowlist containing
  the Rust/C/Wasm sources, production browser frontend, release tooling,
  licenses, and supporting documentation.
