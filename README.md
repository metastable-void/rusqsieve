# rusqsieve

`rusqsieve` is a speed-first, portable Rust/WebAssembly implementation of the
self-initializing quadratic sieve. Its primary target is balanced, RSA-style
semiprimes from **192 through 256 bits**, including browser execution across
independent Web Workers without shared memory.

> **To our knowledge, rusqsieve holds the fastest publicly documented browser-based SIQS timings for balanced 192-, 224-, and 256-bit semiprimes on a consumer mobile device.**

Version 0.3 exposes a deliberately small, safe Rust API and an opaque native C
ABI. SIQS relations, matrix kernels, worker packets, scheduler state, primality
policy, and mutable limbs remain private implementation details.

## Performance status

On the fixed, factor-verified corpus and reference host, rusqsieve currently:

- beats FLINT's single-threaded QSieve on every measured native tier from 160
  through 240 bits;
- factors the 192-, 224-, and 256-bit browser-shaped cases in 0.72 s, 5.04 s,
  and 37.86 s under Node 24.15/V8 with eight workers;
- scales the fixed 256-bit case to 13.96 s with 48 workers on the 96-thread
  reference host.

These results make rusqsieve a **plausibly fastest-class browser integer
factorizer for balanced 192–256-bit semiprimes**, and a strong candidate for the
fastest browser SIQS in that range. This is intentionally not an unqualified
world-record or "fastest general factorizer" claim: current competitor
implementations still need same-browser, same-hardware measurements.

[BENCHMARKING.md](BENCHMARKING.md) contains the inputs, factors, commands,
measurement scope, and competitor protocol. The crate is not constant-time and
must not be used where operand-dependent timing reveals a secret.

## Installation

Add the Rust library:

```sh
cargo add rusqsieve@0.3.0
```

Install the native CLI:

```sh
cargo install rusqsieve --version 0.3.0
```

Or build the optimized native library and CLI from source:

```sh
make
```

The library is emitted as an `rlib`, static library, and platform shared
library.

## Rust API

The blocking factorization functions are available on Unix and Windows:

```rust
use rusqsieve::{FactorConfig, Natural, Parallelism, factor_with};

let input = Natural::<4>::from_decimal("360").unwrap();
let config = FactorConfig::default().with_parallelism(
    Parallelism::threads(4).expect("a nonzero worker count"),
);
let factors = factor_with(input.clone(), config).unwrap();

assert_eq!(factors.distinct_len(), 3);
assert_eq!(factors.total_len(), 6);
assert_eq!(
    factors
        .expanded()
        .map(ToString::to_string)
        .collect::<Vec<_>>(),
    ["2", "2", "2", "3", "3", "5"],
);
assert!(factors.verify_product(&input));
```

Use `factor` for defaults, `factor_with` for configuration, or
`factor_with_progress` for read-only progress snapshots and cooperative
cancellation. The ordinary no-observer path is separately monomorphized so
callback and progress-clock machinery is compiled out.

The default `Natural` has a fixed 1024-bit storage capacity. The factorization
entry points consistently accept values through 512 bits; this is not a claim
that hard 512-bit semiprimes are practical, and NFS is the appropriate
algorithm at that scale. Arithmetic operators on `Natural` wrap at capacity,
while `checked_*` methods report overflow.

All supported public items are covered by rustdoc, enforced with
`deny(missing_docs)`. The complete 0.3 contract and implementation architecture
are documented in [SPEC.md](SPEC.md); breaking changes from 0.1 are summarized
in [CHANGELOG.md](CHANGELOG.md).

## Architecture

The native driver validates and converts the public const-generic `Natural`
into the engine’s fixed working width. The engine performs trial division,
primality and perfect-power checks, bounded Pollard–Brent, then schedules SIQS
polynomial families. Completed families flow through deterministic
single-large-prime relation collection, sparse matrix filtering, verified
row-echelon dependency recovery, and GCD extraction. The Wasm worker and
coordinator exports use the same family kernel and collector through
serialized, ownership-checked packets. `qs` owns factor-base construction and
tier parameters; `natural`, `smallfactor`, `primality`, and `f2` provide the
arithmetic kernels.

## Native C API

Unix and Windows builds export a minimal decimal-string interface from the
static and shared libraries:

```c
#include <stdio.h>
#include "rusqsieve.h"

int main(void) {
    rusqsieve_factors *factors = rusqsieve_factors_new();
    if (factors == NULL)
        return 1;

    enum rusqsieve_status status = rusqsieve_factor("360", 0, factors);
    if (status == RUSQSIEVE_OK) {
        for (size_t i = 0; i < rusqsieve_factors_len(factors); ++i)
            puts(rusqsieve_factors_get(factors, i));
    }

    rusqsieve_factors_free(factors);
    return status;
}
```

The result type is completely opaque and Rust-owned. Factor strings are
borrowed until the next call using that result or until
`rusqsieve_factors_free`; callers must not free individual strings.
`threads == 0` selects available parallelism capped at 48; explicit counts are
capped at 256. Callers can check `rusqsieve_abi_version()` and render statuses
with `rusqsieve_strerror()`.

Install libraries, CLI, header, and pkg-config metadata under `/usr/local`:

```sh
sudo make install
```

The prefix and staging root are overridable:

```sh
make install PREFIX=/usr DESTDIR="$pkgdir"
```

`BINDIR`, `LIBDIR`, `INCLUDEDIR`, and `PKGCONFIGDIR` may also be overridden.
See [rusqsieve.h](rusqsieve.h) for the complete ownership, status, and
thread-safety contract.

## Browser and WebAssembly

Build the self-contained browser demo:

```sh
make docs
make serve
```

The generated `docs/` directory can be deployed directly to GitHub Pages. It
contains scalar and SIMD128 Wasm modules; the frontend attempts SIMD first and
falls back to scalar on older engines.

The browser architecture uses a dedicated coordinator Worker/Wasm instance and
an independent Wasm instance in each sieve Worker. Workers rebuild the same deterministic SIQS
context and return serialized polynomial-family relations. The coordinator
merges families deterministically, filters the matrix, solves for dependencies,
and extracts a verified nontrivial factor.

The reference glue validates packet envelopes and relation-batch framing,
ignores obsolete generations, bounds initialization/jobs/runs with timeouts,
reports relation-budget exhaustion, and rebuilds the Worker runtime after a
failure. An extraction that finds only trivial dependencies retains its
relations and requests a surplus before retrying. Browser input is capped at
the same 512-bit supported limit as native entry points.

Notable performance work includes:

- target-fitted SIQS polynomials and Gray-code root updates;
- translated, sorted roots and a paired root-difference stride loop;
- biased byte logarithmic scores with word-at-a-time candidate rejection;
- multiply-shift-gated survivor division that stops on the recorded score;
- single-large-prime relation combination with a bounded 256× factor-base limit;
- deterministic low-weight sparse matrix elimination;
- compact residual row-echelon solving;
- scoped Wasm SIMD128 XOR acceleration;
- two-family jobs and a measured 48-worker cap.

The SIMD artifact enables the Wasm SIMD128 target feature consistently in both
release build paths. Binaryen `wasm-opt -O3`/`-Oz` is not used because it
regressed the measured sieve.

## Command-line interface

`qs-factor` reads one unsigned decimal integer from standard input and prints
the sorted prime factors, including repetitions, one per stdout line:

```sh
printf '%s\n' 360 | qs-factor --threads auto --progress auto
```

Progress and elapsed time are written only to stderr, keeping stdout
machine-readable. The factorization of one succeeds with no factor lines; zero
and malformed input are errors.

## RSA challenge proof-of-work

The balanced-semiprime SIQS path is the performance-critical deployment target
for sign-in proof-of-work. Challenges must use fresh, similarly sized random
primes, be bound to the intended session, expire, and be replay-protected. A
returned factor must be nontrivial, divide the challenge, and reconstruct it
with its cofactor.

Proof-of-work is resource pricing, not authentication. Retain the normal
authentication mechanism, and never use a modulus belonging to a real RSA key.

ECM is not part of the 0.3 default path. If added later, it must be opt-in
behind a non-default feature and shipped in a separate general-purpose Wasm
artifact. The balanced-RSA artifact must contain no ECM code or initialization;
the fixed 192/224/256-bit corpus remains an A/B gate for runtime, download size,
compilation, startup, and code-cache footprint.

## Scope and limitations

The current pipeline combines trial division, Pollard–Brent rho, primality and
perfect-power checks, and SIQS. It is strongest on balanced semiprimes.
Unbalanced 192–256-bit composites with medium-size factors remain the main
general-factorization gap because ECM is absent.

The current linear algebra uses structured sparse elimination followed by a
compact row-echelon solve; it is not a true block-Lanczos recurrence. That
becomes a future concern for matrices beyond the current practical range.

## Release builds

`build-release.sh` creates versioned archives for Linux GNU and musl, FreeBSD,
Windows MSVC, Apple arm64, and WebAssembly:

```sh
SDKROOT=../MacOSX15.4.sdk ./build-release.sh
```

It uses cross-rs for Linux and FreeBSD, xwin for MSVC, cargo-zigbuild for Apple,
and native Cargo for Wasm. Native archives contain the target-appropriate
libraries, CLI, C header, pkg-config metadata, licenses, changelog, and an
elevation-aware installer. The Wasm archive contains the deployable frontend
and both scalar and SIMD128 modules.

Pass target triples as arguments to build a subset. Run
`./build-release.sh --help` for the supported list and environment overrides.

## License

Licensed under `Apache-2.0 OR MPL-2.0`; see [LICENSE-APACHE](LICENSE-APACHE) and
[LICENSE-MPL](LICENSE-MPL).
