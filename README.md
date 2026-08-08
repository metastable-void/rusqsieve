# rusqsieve

`rusqsieve` is a speed-first, portable Rust/WebAssembly implementation of the
self-initializing quadratic sieve. Its primary target is balanced, RSA-style
semiprimes from **192 through 256 bits**, including browser execution across
independent Web Workers without shared memory.

> **rusqsieve is the world's fastest browser-Wasm SIQS implementation for
> realistic balanced semiprimes in its 192–256-bit target range, established
> through cross-factorizer testing on multiple mobile devices.**

Version 0.4 exposes a deliberately small, safe Rust API and an opaque native C
ABI. SIQS relations, matrix kernels, worker packets, scheduler state, primality
policy, and mutable limbs remain private implementation details.

## Performance status

On the fixed, factor-verified corpus, the 2026-07-31 release gate and subsequent
high-digit tuning measured:

- native single-thread performance that beats FLINT's QSieve on every measured
  tier from 160 through 240 bits;
- factors the 192-, 224-, and 256-bit browser-shaped cases in 0.77 s, 3.97 s,
  and 27.29 s under Node 24.15/V8 with eight workers;
- scales the fixed 256-bit case to 13.96 s with 48 workers on the 96-thread
  reference host;
- factors the checked-in 288-bit balanced fixture in 37.8 s with 48 native
  worker threads;
- factors the reference RSA-100 semiprime in 54.01 s with 192 workers and
  1.95 GiB peak resident memory, versus 50.70 s for portable YAFU;
- factors the reference 364-bit RSA-110 semiprime in 366.12 s with 192 workers
  and 2.06 GiB peak resident memory, versus YAFU's 405.06 s and 9.41 GiB on
  the same host.

Together with independent comparisons against online factorizers on multiple
mobile devices, these results establish rusqsieve as the **world's fastest
browser-Wasm SIQS for realistic balanced 192–256-bit semiprimes**. This is a
specific SIQS and balanced-semiprime claim, not a claim to be the fastest
general-purpose factorizer for every composite shape.

[BENCHMARKING.md](BENCHMARKING.md) contains the inputs, factors, commands,
measurement scope, and competitor protocol. The crate is not constant-time and
must not be used where operand-dependent timing reveals a secret.

## Installation

Add the Rust library:

```sh
cargo add rusqsieve@0.4.3
```

Install the native CLI:

```sh
cargo install rusqsieve --version 0.4.3
```

Or build the optimized native library and CLI from source:

```sh
make
```

The library is emitted as an `rlib`, static library, and platform shared
library.

`make` and `build-release.sh` add one measured tuning flag to native builds,
`-C llvm-args=-align-all-nofallthru-blocks=5`. Without it the alignment of the
hot sieve and Pollard–Brent loops is decided by luck that unrelated edits
re-roll — two inert lines in the CLI's progress closure moved a rho-dominated
input by 8% — and forcing 32-byte alignment measured 2.1–4.0% faster on balanced
192–272-bit inputs and 6.3–8.6% faster on Pollard–Brent-dominated ones, for a
binary about 15% larger. The build probes for it, so a toolchain whose LLVM does
not accept the option falls back to an unflagged build rather than failing;
`make NATIVE_TUNE_RUSTFLAGS=` or `RUSQSIEVE_TUNE_RUSTFLAGS=` opts out. A plain
`cargo build --release` is unflagged and correspondingly unpredictable at the
few-percent level. Wasm is deliberately excluded: block alignment means nothing
in a stack machine, and artifact size is a shipping gate there.

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

The default `Natural` has a fixed 1024-bit storage capacity, and the
factorization entry points accept any value that fits it. Input width is not
itself a limit: small factors are removed by trial division, perfect-power
detection, and Pollard–Brent, none of which is width-gated, so a very wide
number built from small factors factors normally.

The limit applies to the *composite that reaches the quadratic sieve*, which is
capped at 400 bits — a little over 120 decimal digits, roughly where GNFS
overtakes SIQS by margins no sieve tuning recovers. A composite above that
returns `FactorError::SiqsCompositeTooLarge`, carrying the composite's bit
length rather than the caller's. Arithmetic operators on `Natural` wrap at
capacity, while `checked_*` methods report overflow.

Because that composite has nowhere else to go, Pollard–Brent runs a much deeper
budget above the ceiling than below it — 26 to 36 s of iterations, which
reaches a smallest factor near 2^53 at 512 bits and 2^50 at 1024, and covers
factors up to 32 bits with orders of magnitude to spare. A wide input carrying a
findable factor is therefore split rather than refused. Set
`RUSQSIEVE_RHO_ITERATIONS` for the minutes-to-hours search that 56- and 64-bit
factors cost.

A cofactor that reached the recursion by splitting under rho keeps that deep
budget down to 257 bits, because splitting proves it is not the balanced
semiprime the cheap budget is sized for. This is what lets a wide product of
many middling primes finish: a 498-bit product of ten 50-bit primes peels five
factors in rho and hands a 250-bit remainder to the sieve, 65 s in total, where
it previously stalled on a 399-bit sieve for weeks. Balanced semiprimes never
reach that branch — rho does not split them.

All supported public items are covered by rustdoc, enforced with
`deny(missing_docs)`. The complete 0.4 contract and implementation architecture
are documented in [SPEC.md](SPEC.md); breaking changes from 0.1 are summarized
in [CHANGELOG.md](CHANGELOG.md).

## Architecture

The native driver validates and converts the public const-generic `Natural`
into the engine’s fixed working width. The engine performs trial division,
primality and perfect-power checks, bounded Pollard–Brent, the elliptic curve
method where its conditions are met, then schedules SIQS polynomial families. Native collection consumes completed unique families
immediately to avoid head-of-line stalls; the portable coordinator merges
serialized families deterministically for browser/distributed use. Both flow
through the same large-prime graph, sparse matrix filtering, verified
dependency recovery, and GCD extraction. `qs` owns factor-base construction
and tier parameters; `natural`, `smallfactor`, `primality`, and `f2` provide
the arithmetic kernels.

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
ignores obsolete generations, bounds initialization with a timeout, detects a
stalled factorization after ten visible minutes without worker activity,
reports relation-budget exhaustion, and rebuilds the Worker runtime after a
failure. Advancing runs have no wall-clock limit, and time spent with the page
hidden does not count as a stall. An extraction that finds only trivial
dependencies retains its relations and requests a surplus before retrying. The
browser peels small factors on the main thread and applies the same 400-bit
sieve limit as native entry points to the composite it hands the coordinator.
Above that limit — and below it for a cofactor an earlier split already proved
unbalanced — it first runs a deep Pollard–Brent in wasm (`qs_rho`) across a pool
of rho workers, each walking a disjoint range of polynomial constants so the pool
races that many independent walks. That keeps the main thread free and reaches a
smallest factor of about 2^52 at 512 bits, against 2^45 for the sliced `BigInt`
search it replaces, which remains the fallback where wasm or `Worker` is
unavailable. The scheduler reads its family budget from the session rather than
assuming a constant, so it cannot stop a large run before the engine would.

Notable performance work includes:

- target-fitted SIQS polynomials and Gray-code root updates;
- translated, sorted roots and a paired root-difference stride loop;
- biased byte logarithmic scores with word-at-a-time candidate rejection;
- multiply-shift-gated survivor division that stops on the recorded score;
- tier-bounded single- and double-large-prime graph combination;
- deterministic low-weight sparse matrix elimination;
- compact scalar/M4RI residual solving below 272 bits and verified 64-way
  Montgomery block Lanczos from 272 bits upward;
- scoped Wasm SIMD128 acceleration for matrix XOR and SIQS root advancement;
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

ECM never runs on a balanced semiprime inside the sieve's range unless the
caller asks for it, so proof-of-work pricing is unchanged by its presence: the
conditions that admit it without asking are all evidence that the composite is
some other shape. The fixed 192/224/256-bit corpus remains an A/B gate for
runtime, download size, compilation, startup, and code-cache footprint.

## Scope and limitations

The current pipeline combines trial division, Pollard–Brent rho, primality and
perfect-power checks, the elliptic curve method, and SIQS. It is strongest on
balanced semiprimes.

**That gap is now filled by ECM.** Montgomery curves with Suyama's `σ`, stage 1
along PRAC addition chains, and a standard-continuation stage 2 over a 210-wheel,
with one gcd per stage. Its cost is governed by the size of the *factor* rather
than of the input, which is what makes it the right tool where the other two are
wrong: it recovered a 20-digit factor of a 466-bit composite in 30 s, an input
that previously returned `SiqsCompositeTooLarge`.

Balanced semiprimes still never pay for it. Curves run without being asked only
where the balanced premise has already been disproved — a composite wider than
the sieve accepts, or one whose small factor trial division peeled or whose
ancestor Pollard–Brent split — and a balanced semiprime inside the sieve's range
satisfies neither. For that range it is opt-in: `FactorConfig::with_ecm(true)`,
`RUSQSIEVE_FLAG_ENABLE_ECM`, or `qs-factor --enable-ecm`. A default run of a
256-bit balanced semiprime measured 6.72 s against 6.55 s with curves enabled,
which is noise; the reason to keep it off is that no curve can succeed there, not
that it is expensive.

What remains of the original roadmap item is stage 2's asymptotics: GMP-ECM
evaluates a polynomial at many points at once where this evaluates point by
point. The gap is throughput, not reach.

Below 272 bits, linear algebra uses structured sparse elimination followed by
compact scalar/M4RI row-echelon solving. At 272 bits and above, native and Wasm
use a true 64-way Montgomery block-Lanczos recurrence. The 281–288 tier also
uses a larger factor base now that its residual matrix no longer incurs dense
elimination cost; double-large-prime graph collection begins above 288 bits.

On the fixed browser fixtures, Lanczos reduced 272-bit LA/extraction from
5.574 s to 1.721 s. On an exact 288-bit balanced fixture, the larger base plus
Lanczos reduced wall time from 420.223 s to 258.512 s. A matched 256-bit run,
which stays on M4RI, was unchanged within noise (27.146 s versus 26.986 s).
The final shipped-artifact confirmation in real Chromium completed the
256-/272-/288-bit fixtures in 23.702/69.783/226.901 s. Post-tuning endpoint
checks completed the 272-/288-bit fixtures in 67.433/222.006 s. Wasm SIMD root
advancement and the portable score-stream changes subsequently reduced the
verified endpoints to 62.237/185.830 s; the multiplier and Q/2 policy remains
gated above 288 bits.

On the 96-thread native tuning host, the original final-default fixed anchors completed in
41.8 s at 289 bits and 103.5 s at 304 bits. Full multiplier selection, Q/2
polynomials, larger self-initializing families, allocation-free report
resieving, precomputed score weights, portable dense-prefix blocking, a
102-bit/DLP collection policy, allocation-free partial-graph updates, and a
portable structure-of-arrays score stream reduced RSA-100 from 622.6 s to a
339.2 s profiled run, then to 320.1 s in the final release gate. The score
kernel reads contiguous primes, separates repeated-hit and sparse ranges, and
uses reusable padded sentinels instead of unpredictable sparse-root bounds
branches; none of those changes depend on x86-64. x86-64 builds additionally
dispatch root advancement through an SSE2 baseline or runtime-detected AVX2
kernel, while Wasm SIMD128 has an equivalent optional path and every target
retains scalar Rust fallbacks. Biased report streams use AVX2/SSE2 movemasks or
Wasm `i8x16.bitmask` to reject 16–32 positions per branch. RSA-110 uses all
4,096 nonredundant Gray variants per 13-factor A; the exact 320-bit crossover
instead uses 1,024-pattern packets and a wider, YAFU-matched sieve geometry.

On a 192-core Xeon 6975P-C, the tuned 320-bit case completed in 33.11 s versus
33.60 s for the supplied portable YAFU binary; using all 384 logical CPUs took
30.50 s. The retuned RSA-100 midpoint completed in 54.01 s versus YAFU's
50.70 s. RSA-110 completed in 366.12 s versus YAFU's 405.06 s, while peak RSS
was 2.06 GiB versus 9.41 GiB. The row-major, eight-way block-Lanczos multiply
accounts for 15.72 s of the rusqsieve RSA-110 run. These are single-input
gates, not cross-host claims; exact inputs and commands are recorded in
`BENCHMARKING.md`.

## Release builds

`build-release.sh` creates eight versioned archives: x86-64 and arm64 Linux
GNU, x86-64 and arm64 Linux musl, x86-64 FreeBSD, x86-64 Windows MSVC, Apple
arm64, and WebAssembly:

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
The 0.4 release gate built and integrity-checked all eight archives, executed
the GNU and musl x86-64 CLIs, linked the C smoke program against the packaged
static and shared libraries, and verified the shipped frontend in real
Chromium.

## License

Licensed under `Apache-2.0 OR MPL-2.0`; see [LICENSE-APACHE](LICENSE-APACHE) and
[LICENSE-MPL](LICENSE-MPL).
