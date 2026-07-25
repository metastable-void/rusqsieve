# Changelog

## 0.2.0

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
