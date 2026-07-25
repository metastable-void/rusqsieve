# rusqsieve

`rusqsieve` is a speed-first, portable Rust/WebAssembly implementation of the
self-initializing quadratic sieve. Its primary target is balanced, RSA-style
semiprimes from **192 through 256 bits**, scheduled across browser Web Workers
without shared memory.

On the fixed, factor-verified corpus, rusqsieve currently:

- beats FLINT's single-threaded quadratic sieve on every measured native tier
  from 160 through 240 bits;
- factors the 192-, 224-, and 256-bit browser cases in 0.72 s, 5.04 s, and
  37.86 s under Node 24.15/V8 with eight workers;
- scales the fixed 256-bit case to 13.96 s with 48 workers on the 96-thread
  reference host.

These results make rusqsieve a **plausibly fastest-class browser integer
factorizer for balanced 192–256-bit semiprimes**, and a particularly strong
candidate for the fastest browser SIQS in that range. This is deliberately not
an unqualified world-record claim: current Alpertron, Msieve-Wasm, and other
browser implementations still need to be run on the same browser, hardware, and
fixed corpus. See [BENCHMARKING.md](BENCHMARKING.md) for inputs, factors,
commands, measurement scope, and the competitor protocol.

The default `Natural<16>` has a 1024-bit storage capacity. That is not a promise
that hard 1024-bit semiprimes are practical; NFS is the appropriate algorithm at
that scale.

The crate is not constant-time and must not be used where operand-dependent timing is secret.

## Native API

```rust
use rusqsieve::{FactorConfig, Natural, Parallelism, factor_with};

let input = Natural::<4>::from_decimal("360").unwrap();
let config = FactorConfig::default().with_parallelism(
    Parallelism::threads(4).expect("a nonzero worker count"),
);
let factors = factor_with(input.clone(), config).unwrap();

assert_eq!(factors.distinct_len(), 3);
assert_eq!(factors.total_len(), 6);
assert!(factors.verify_product(&input));
```

The 0.2 Rust API deliberately exposes only the safe blocking interface,
configuration builders, progress snapshots, fixed-capacity input type, and
owned factorization result. SIQS relations, matrix kernels, worker packets, and
scheduler state are private implementation details. This keeps invalid
relations or mismatched worker contexts from being representable through the
supported Rust API and lets those internals evolve without breaking callers.
All public items are covered by rustdoc, enforced during the build with the
`missing_docs` lint.

### C API

Native Unix and Windows builds export a minimal owning C interface from the
`cdylib` and `staticlib`. The ABI uses standard C `size_t` for Rust `usize` and
NUL-terminated decimal strings for all integers.

`make` builds the optimized native library and `qs-factor` CLI. Installation
defaults to `/usr/local`:

```sh
sudo make install
# Packaging/staging example:
make install PREFIX=/usr DESTDIR="$pkgdir"
```

`BINDIR`, `LIBDIR`, `INCLUDEDIR`, and `PKGCONFIGDIR` are independently
overridable. Installation includes the shared and static C libraries,
`rusqsieve.h`, the CLI, and `rusqsieve.pc`.

`./build-release.sh` creates versioned `.tar.gz` archives for the supported
Linux GNU/musl, FreeBSD, Windows MSVC, Apple arm64, and WebAssembly targets.
It uses cross-rs for Linux and FreeBSD, xwin for MSVC, zigbuild for Apple
(`SDKROOT` is required), and Cargo directly for Wasm. Native archives contain
the target-appropriate libraries, CLI, header, pkg-config metadata, licenses,
and an elevation-aware installer. The Wasm archive is directly deployable and
contains both scalar and SIMD128 builds. Pass target triples as arguments to
build only a subset; run `./build-release.sh --help` for the full list.

```c
#include "rusqsieve.h"

rusqsieve_factors *factors = rusqsieve_factors_new();
if (factors == NULL)
    return 1;

int status = rusqsieve_factor("360", 0, factors);
if (status == RUSQSIEVE_OK) {
    for (size_t i = 0; i < rusqsieve_factors_len(factors); ++i)
        puts(rusqsieve_factors_get(factors, i)); /* 2, 2, 2, 3, 3, 5 */
}

rusqsieve_factors_free(factors);
```

`threads == 0` selects available parallelism capped at 48. The opaque output
owns sorted decimal prime factors, repeated according to multiplicity and
borrowed through `rusqsieve_factors_get`. C callers must release the complete
result with `rusqsieve_factors_free`; individual strings must not be freed. See
[`rusqsieve.h`](rusqsieve.h) for the full lifetime and error contract.

## Why it is fast in browsers

- Numerically target-fitted SIQS polynomials and a double-large-prime relation
  graph reduce the amount of sieving needed.
- Translated, sorted roots and paired stride loops keep division out of the
  score-write hot path.
- Low-weight sparse elimination shrinks the matrix before a compact row-echelon
  dependency solve.
- A scoped Wasm `simd128` XOR kernel accelerates linear algebra without applying
  whole-program SIMD transformations that regress the sieve.
- Two-family Worker jobs limit end-of-run overshoot. The pool follows
  `navigator.hardwareConcurrency` up to a measured cap of 48 workers.
- The deployment preserves Rust/LLVM's speed-optimized Wasm. Binaryen 120
  `wasm-opt -O3` and `-Oz` were both measured and rejected because they slowed
  the 192-bit sieve by roughly 50%.

## Browser demo

`make docs` builds a self-contained WebAssembly demo into `docs/` for GitHub Pages
(enable Pages on the `docs/` folder). It factors a number you type using the same
crate compiled to `wasm32-unknown-unknown`, sieving in parallel across a pool of Web
Workers, and renders the result in power
notation. Modern browsers use the scoped `simd128` linear-algebra kernel; a portable
scalar artifact is selected automatically on older engines. `make serve` previews it
locally at <http://localhost:8000/>.

## RSA challenge proof-of-work

The balanced-semiprime SIQS path is the performance-critical deployment target
for sign-in proof-of-work. Server-generated challenges should use similarly
sized random primes and enter the SIQS coordinator directly through the
low-level engine/Worker API. The returned factor must be nontrivial, divide the
challenge, and reconstruct the original challenge with its cofactor.

ECM is not currently part of this path. If ECM is added for unbalanced,
general-purpose composites, the project policy is:

1. ECM is opt-in behind a non-default feature and a separate general-purpose
   Wasm artifact.
2. The balanced-RSA artifact contains no ECM code or initialization.
3. The fixed 192/224/256-bit corpus is an A/B regression gate for native time,
   browser time, Wasm size, and module startup.

This separation is intentional: merely skipping ECM at runtime would still let
its additional Wasm download and compilation cost slow a sign-in challenge.

Proof-of-work is resource pricing, not authentication. Bind each challenge to
the intended session, expire it, reject replay, and retain the normal
authentication mechanism. Never use a modulus belonging to a real RSA key.

## Scope

The current pipeline combines trial division, Pollard-Brent rho, primality and
perfect-power checks, and SIQS. It does not yet include ECM, so the fastest-class
claim is specifically for balanced semiprimes. Unbalanced 192–256-bit composites
with medium-size factors remain the main general-factorization gap.

Licensed under `Apache-2.0 OR MPL-2.0`; see `LICENSE-APACHE` and `LICENSE-MPL`.
