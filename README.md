# rusqsieve

`rusqsieve` is a portable Rust factorization crate built around a fixed-capacity unsigned
integer, deterministic work packets, and quadratic-sieve relation/matrix primitives.  The default
`Natural<16>` has a 1024-bit capacity; that is a storage limit, not a promise that hard 1024-bit
semiprimes are practical (use NFS for such inputs).

The crate is not constant-time and must not be used where operand-dependent timing is secret.

```rust
use rusqsieve::{Natural, factor};

let factors = factor(Natural::<16>::from_decimal("360").unwrap()).unwrap();
assert!(factors.verify_product(&Natural::from_u64(360)));
```

Custom schedulers, including Web Workers, can use `engine::prepare`,
`EngineSession::take_jobs`, `engine::execute`, and `EngineSession::submit`.
The portable job kernel never creates threads; the native blocking API schedules
the same deterministic polynomial-family work across persistent workers.

Licensed under `Apache-2.0 OR MPL-2.0`; see `LICENSE-APACHE` and `LICENSE-MPL`.
