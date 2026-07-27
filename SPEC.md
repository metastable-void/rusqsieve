# rusqsieve 0.3 implementation specification

This document specifies the supported behavior and current architecture of
rusqsieve 0.3. It describes the implementation that is shipped, not an
aspirational module layout or a compatibility promise for private internals.

Normative requirements use **must**. Descriptions of tuning and implementation
strategy document the 0.3 release and may change in later compatible releases
when observable behavior is preserved.

## 1. Purpose and scope

rusqsieve is a high-performance integer-factorization package written in Rust.
Its performance-critical workload is a balanced, RSA-style semiprime between
192 and 256 bits. It supports:

- a safe native Rust API;
- an opaque decimal-string C API on Unix and Windows;
- a native command-line program;
- a raw `wasm32-unknown-unknown` module;
- a browser frontend that distributes SIQS polynomial families across
  independent Web Workers without shared memory.

The factorization pipeline combines trial division, probable-prime testing,
perfect-power detection, Pollard–Brent rho, and a self-initializing quadratic
sieve (SIQS).

The following are explicitly outside the 0.3 scope:

- ECM and the General Number Field Sieve;
- constant-time or side-channel-resistant arithmetic;
- factoring hard semiprimes near the `Natural` storage limit;
- a public Rust API for relations, matrices, scheduler state, or worker
  packets;
- a stable ABI for Rust types;
- Rust threads, shared Wasm memory, or `SharedArrayBuffer` on
  `wasm32-unknown-unknown`.

The absence of ECM is intentional for the balanced-semiprime proof-of-work
artifact. Any future ECM implementation must be opt-in and must not add runtime,
download, compilation, initialization, or code-cache cost to the default
balanced-RSA path.

## 2. Release and compatibility boundaries

### 2.1 Supported public interfaces

The supported interfaces in 0.3 are:

1. the items re-exported from `src/lib.rs`;
2. the factor/result functions, ABI query, status formatter, and opaque type
   declared in `rusqsieve.h`;
3. the `qs-factor` command-line behavior documented below;
4. the Wasm exports used by `web/abi.js`, `web/index.js`, and `web/worker.js`.

Everything else in `src/` is private implementation detail even when an
internal item uses Rust's `pub` visibility inside its private module.

Minor 0.3 releases may retune parameters, change internal representations, add
non-exhaustive error/progress variants, or replace private algorithms. They
must not expose invalid relation, matrix, pointer, or scheduler states through
the safe Rust API.

### 2.2 Targets

The safe blocking Rust API and C ABI are built on Unix and Windows. The raw
browser ABI is built only for:

```text
wasm32-unknown-unknown
```

The release builder supports:

```text
x86_64-unknown-linux-gnu
x86_64-unknown-linux-musl
aarch64-unknown-linux-gnu
aarch64-unknown-linux-musl
x86_64-unknown-freebsd
x86_64-pc-windows-msvc
aarch64-apple-darwin
wasm32-unknown-unknown
```

The musl archives contain the CLI and static library. GNU and other native archives
contain the CLI, static library, shared library, header, pkg-config metadata,
and an installer. The Wasm archive contains scalar and SIMD128 modules plus the
deployable browser frontend.

## 3. Package structure

rusqsieve is one Cargo package with one library and one optional CLI:

```text
src/
├── lib.rs              public Rust surface and target selection
├── native.rs           safe blocking native driver
├── capi.rs             native C ownership/pointer boundary
├── engine.rs           optimized SIQS engine and portable jobs
├── qs/mod.rs           factor-base construction and SIQS tier parameters
├── f2/mod.rs           sparse filtering and dependency solving
├── natural/mod.rs      fixed-capacity unsigned arithmetic
├── smallfactor.rs      cached small primes and u64 Pollard–Brent
├── primality.rs        probable-prime testing
├── factor.rs           public configuration and error vocabulary
├── factors.rs          owned factor result
├── progress.rs         public progress vocabulary
├── wasm.rs             raw Wasm ABI and handle registries
└── bin/qs-factor.rs    native CLI
```

The library emits:

```toml
crate-type = ["rlib", "cdylib", "staticlib"]
```

The default Cargo feature enables the CLI. The `wasm-simd128` feature enables
the SIMD128 linear-algebra kernel. Cargo features are additive: none changes the
identity or default const-generic width of a public type.

## 4. Public Rust API

### 4.1 Exported surface

The crate root exports only:

```rust
pub use factor::{
    FactorConfig, FactorError, Parallelism, ProgressAction, ResourceLimitKind,
};
pub use factors::{ExpandedPrimeFactors, PrimeFactorIter, PrimeFactors};
pub use natural::{
    BufferTooSmall, CapacityError, InvalidDigit, Natural, ParseNaturalError,
};
pub use progress::{
    ProgressAmount, ProgressPhase, ProgressSnapshot, ProgressTotal, ProgressUnit,
};

#[cfg(any(unix, windows))]
pub use native::{factor, factor_with, factor_with_progress};
```

It also exports the `natural!` compile-time decimal literal macro.

All public items must have rustdoc. The crate enforces this with
`#![deny(missing_docs)]`.

### 4.2 Blocking factorization

On Unix and Windows:

```rust
pub fn factor<const P: usize>(
    input: Natural<P>,
) -> Result<PrimeFactors<P>, FactorError>;

pub fn factor_with<const P: usize>(
    input: Natural<P>,
    config: FactorConfig,
) -> Result<PrimeFactors<P>, FactorError>;

pub fn factor_with_progress<const P: usize, F>(
    input: Natural<P>,
    config: FactorConfig,
    observer: F,
) -> Result<PrimeFactors<P>, FactorError>
where
    F: FnMut(&ProgressSnapshot) -> ProgressAction;
```

Observable behavior:

- zero returns `FactorError::ZeroHasNoPrimeFactorization`;
- one succeeds with an empty `PrimeFactors`;
- factors are sorted in ascending order;
- multiplicities are preserved;
- every returned factor has passed the implementation's probable-prime test;
- the result can verify its product against the original input;
- callback cancellation returns `FactorError::Cancelled`.

The no-observer path must remain separately monomorphized so progress timing and
callback machinery are compiled out of ordinary `factor` and `factor_with`
calls.

### 4.3 Configuration

`FactorConfig` is an owned, encapsulated configuration. In 0.3 its supported
controls are:

- `parallelism()` / `with_parallelism(...)`;
- `progress_interval()` / `with_progress_interval(...)`;
- `with_witness_seed(...)` for a reproducible ChaCha8 Miller–Rabin witness
  stream.

SIQS parameters and relation limits remain private so the implementation can be
improved without breaking callers. A doc-hidden tuning constructor exists only
so the CLI can map benchmark environment variables into an owned configuration;
the library itself never reads ambient environment state.

`Parallelism::Auto` detects available native parallelism when factorization
begins. `Parallelism::Threads(NonZeroUsize)` requests a nonzero worker count.
`Parallelism::threads(0)` returns `None`.

### 4.4 Factor results

`PrimeFactors<P>` owns a sorted map from factor to nonzero multiplicity. Its
public operations provide:

- iteration over distinct `(factor, multiplicity)` pairs;
- multiplicity lookup;
- expanded iteration with repetitions;
- distinct and total cardinality;
- empty-result detection;
- checked product verification.

Callers cannot construct or mutate the internal map directly.

### 4.5 Progress

`ProgressSnapshot` and `ProgressAmount` are read-only values. A snapshot exposes:

- a monotonically increasing revision;
- the input bit length;
- a non-exhaustive high-level phase;
- phase-specific completed/total/unit counters.

Totals distinguish exact, estimated, and unknown values. Estimated totals may
change and fractions are informational. Progress callbacks run on the calling
thread, must not run while worker locks are held, and should return quickly.

Returning `ProgressAction::Cancel` is cooperative. During parallel sieving the
driver sets a shared atomic cancellation flag, stops queueing jobs, asks workers
to abandon queued work, joins them, and returns `FactorError::Cancelled`.

## 5. `Natural`

`Natural<const P: usize>` is a `repr(transparent)`, inline, fixed-capacity
unsigned integer containing `P` little-endian `u64` limbs.

The default capacity is fixed at `Natural<16>` / 1024 bits. Blocking
factorization consistently rejects values wider than 512 bits, independent of
the const-generic storage width selected by the caller.

This is a storage limit, not a practical factorization claim. NFS is appropriate
for hard composites substantially beyond this project's SIQS range.

Public behavior:

- decimal parsing accepts ASCII digits and leading zeroes;
- signs, whitespace, separators, and empty input are rejected;
- overflow is reported exactly;
- big- and little-endian decoding reject nonzero excess bytes;
- serialization writes the shortest unsigned encoding, with zero represented by
  zero bytes;
- arithmetic operators wrap modulo `2^(64P)`;
- `checked_add`, `checked_sub`, and `checked_mul` report overflow;
- division by zero is represented by `None` in `div_rem`;
- formatting is canonical unsigned decimal;
- `natural!` rejects malformed or overflowing literals at compile time.

The public type exposes no mutable limb slice. Internal arithmetic uses
significant-limb-aware operations, widening multiplication, normalized limb
division, binary GCD, and modular helpers. Mathematical code must use checked,
widening, or modular operations wherever wrapping would invalidate an
invariant.

Native big-integer Pollard–Brent constructs one real Montgomery context for its
odd modulus. Polynomial values and batched products remain encoded throughout
the stage; REDC uses limb multiplication and carry propagation rather than
division. The iteration budget is unchanged from the division-based version,
so this optimization reduces failed-rho overhead rather than deepening the
search on balanced semiprimes.

## 6. Native factorization pipeline

The optimized blocking engine performs:

1. reject zero;
2. divide by cached primes through 10,000;
3. use deterministic machine-word primality and Pollard–Brent for cofactors
   fitting `u64`;
4. run probable-prime testing on larger cofactors;
5. detect perfect powers and recursively factor the base;
6. run SIQS to recover a nontrivial divisor;
7. recursively factor divisor and cofactor;
8. sort the complete factor list.

For input capacities at or below the compiled engine width, the high-level API
converts into the optimized engine without changing the value. Wider
user-selected capacities retain a private reference fallback. This distinction
is not a separate public API.

The implementation is deterministic for a fixed version, input, configuration,
and relevant tuning environment. Parallel workers do not choose random
polynomials. Relation results are merged by family number rather than arrival
order.

## 7. SIQS engine

### 7.1 Parameters and multiplier

Bit-size tiers in `qs::parameters::engine_params` select the factor-base bound,
sieve half-width, and large-prime allowance. Environment overrides are
development tuning aids and are not stable public configuration.

The engine chooses a deterministic Knuth–Schroeppel multiplier `k` and sieves
against `kN`. Extraction still computes GCDs against `N`; because `kN` is zero
modulo `N`, the square congruence remains valid for factoring `N`.

### 7.2 Factor base

The factor base contains 2, primes dividing the multiplier where applicable,
and odd primes for which the sieved value is a quadratic residue. Modular
square roots are computed with Tonelli–Shanks. Each entry stores the prime,
root, and rounded logarithmic weight.

Lemire fast-mod constants and `interval mod p` are precomputed once per engine
context. They eliminate hardware division from candidate root tests and
per-polynomial root translation.

### 7.3 Polynomial families

Each deterministic family selects a numerically target-fitted squarefree `A`
from factor-base primes. `B` is represented as the smaller signed CRT
coefficient and satisfies:

```text
B² ≡ kN (mod A)
```

Related `B` values are visited in Gray-code order. Per-prime
`2 B_j A⁻¹ mod p` increments are precomputed, so roots advance with one modular
add/subtract per prime. Roots remain sorted and are represented directly as
score-array positions.

### 7.4 Logarithmic sieve

Workers reuse score, root, increment, and candidate buffers across polynomials.

The sieve uses bit-length `u8` scores and:

1. biases every score by `128 − threshold`, so reaching the threshold sets the
   byte's high bit;
2. skips selected very small primes and derives their expected threshold slack
   from that set;
3. adds `ceil(log2 p)` weights at both sorted modular roots, wrapping rather
   than saturating whenever the engine can prove from the smallest scored prime
   that no position can carry past 255;
4. uses the paired root-difference stride loop;
5. extracts candidates eight positions at a time with one masked compare;
6. trial-divides survivors.

For each survivor, `g(x) = Q(x)/A` is reconstructed directly as a signed value.
Candidate division is gated by a precomputed multiply-shift residue test and
stops as soon as the primes divided out account for the recorded score, which is
exact when the scores did not saturate.

The sieve threshold is `log2|g(x)| − log2(large-prime bound) − small-prime slack`
plus a measured per-tier offset. The large-prime term is the acceptance bound
itself, not an independent constant: a survivor whose cofactor exceeds what
`classify_cofactor` will accept costs a full trial division for no possible
relation.

### 7.5 Relations and large primes

An accepted relation represents:

```text
t² ≡ (-1)^sign × product(p_i ^ e_i) × large_parts (mod N)
```

where `t = Ax + B` reduced modulo `N`.

Full relations have no large-prime cofactor. The shipped engine accepts one
probable-prime large factor up to 256 times the factor-base bound. Double-large-
prime collection remains disabled because it did not produce a net wall-time win. The
coordinator treats partials as edges in a large-prime graph. A cycle combines
relations only when every large-prime exponent cancels to even parity; the
corresponding square-root factors are retained for extraction.

The relation collector buffers out-of-order families and merges them in
ascending family order. Native workers and Web Workers execute the same
`sieve_family` kernel through different schedulers.

## 8. Linear algebra and extraction

Matrix columns correspond to combined relations. Rows correspond to the sign
and factor-base-prime exponent parities.

The 0.3 solver is:

1. sparse structured elimination, including singleton removal and deterministic
   low-weight row elimination through weight six;
2. compact pivot records used to expand residual dependencies back into the
   original relation-column space;
3. a compact row-echelon solve on the residual matrix;
4. expansion of at most 64 useful dependencies;
5. verification of every expanded dependency against the original parity
   matrix.

The crate does **not** implement a block-Lanczos recurrence. The current
residual solver is compact dense Gaussian/row-echelon elimination.

For a verified dependency, extraction constructs:

```text
x = product(relation roots) mod N
y = square root of the combined factor-base and large-prime product mod N
```

It tries `gcd(|x-y|, N)` and `gcd(x+y, N)` and accepts only a factor strictly
between one and `N`.

## 9. Native scheduling

The native engine creates persistent worker threads for a factor attempt and
uses bounded channels of deterministic polynomial-family identifiers.

Small inputs cap the effective worker count to avoid startup overhead. The
automatic C path caps detected parallelism at 48. The browser similarly caps
its Web Worker pool at 48 based on measured scaling on the reference host.

The coordinator:

- owns relation and matrix state;
- keeps at most a bounded multiple of the worker count in flight;
- merges families deterministically;
- stops dispatching once the relation target is met;
- joins every worker before returning;
- never calls user progress code from worker threads.

## 10. C API and ABI

The native C ABI is declared by `rusqsieve.h`. It uses `size_t` and
NUL-terminated decimal strings:

```c
typedef struct rusqsieve_factors rusqsieve_factors;

uint32_t rusqsieve_abi_version(void);
const char *rusqsieve_strerror(int status);
rusqsieve_factors *rusqsieve_factors_new(void);
void rusqsieve_factors_free(rusqsieve_factors *factors);
size_t rusqsieve_factors_len(const rusqsieve_factors *factors);
const char *rusqsieve_factors_get(
    const rusqsieve_factors *factors,
    size_t index
);
enum rusqsieve_status rusqsieve_factor(
    const char *n,
    size_t threads,
    rusqsieve_factors *factors
);
enum rusqsieve_status rusqsieve_factor_with_progress(
    const char *n,
    size_t threads,
    rusqsieve_factors *factors,
    rusqsieve_progress_callback callback,
    void *context
);
```

Ownership and lifetime rules:

- only the constructor creates a result object;
- every non-null result must be freed exactly once;
- freeing null is allowed;
- returned factor strings are borrowed and must not be freed;
- strings remain valid until the next factorization into that result or its
  destruction;
- factorization clears the previous result before reporting success or failure;
- the implementation copies `n` before clearing, so a borrowed factor may be
  reused as the next input on the same result;
- operations on the same result must not overlap across threads;
- independent result objects may be used concurrently;
- Rust panics are caught inside `rusqsieve_factor` and mapped to
  `RUSQSIEVE_INTERNAL_ERROR`; shipped library profiles use unwinding so this
  contract remains active.

`threads == 0` uses available parallelism capped at 48. Positive values are
capped at 256 and remain subject to internal small-input caps.

The shared library must export only:

```text
rusqsieve_factor
rusqsieve_factor_with_progress
rusqsieve_abi_version
rusqsieve_factors_free
rusqsieve_factors_get
rusqsieve_factors_len
rusqsieve_factors_new
rusqsieve_strerror
```

## 11. WebAssembly architecture

### 11.1 Execution model

The browser frontend uses:

- one coordinator Wasm instance on the main thread;
- one independent Wasm instance per Web Worker;
- transferable serialized relation packets;
- no shared Wasm memory and no atomics.

Each worker rebuilds the immutable deterministic sieve context from the decimal
input. The coordinator assigns disjoint family ranges, buffers results by family
number, combines relations, runs linear algebra, and extracts a factor.

JavaScript performs inexpensive recursive preprocessing with `BigInt`: trial
division, Miller–Rabin, perfect powers, and bounded Pollard–Brent. Hard
composites are sent to the Wasm SIQS coordinator.

The reference frontend enforces the same 512-bit input ceiling as the native
entry points. Per-run messages, including errors, are generation-scoped.
Worker initialization, individual sieve jobs, and complete runs are bounded by
timeouts; a failed run terminates and rebuilds the Worker runtime. The browser
scheduler issues at most the engine's 100,000-family budget and reports
exhaustion once no jobs or submissions remain. If the first dependency set
yields no nontrivial factor, the coordinator retains its session, requests a
relation surplus, and retries instead of discarding collected work.

The production frontend is `web/`. The older `js/` prototype is not part of the
published crate or release archives.

### 11.2 Active raw ABI

The browser frontend relies on:

```text
qs_abi_version
qs_alloc
qs_dealloc
qs_buffer_pointer
qs_buffer_length
qs_buffer_free

qs_worker_prepare
qs_worker_sieve
qs_worker_free

qs_coord_new
qs_coord_target
qs_coord_relations
qs_coord_submit
qs_coord_extract
qs_coord_free
```

Handles contain a slot and generation, so stale handles do not alias newly
allocated objects. Incoming pointers and lengths are checked with checked
arithmetic against current Wasm memory and a 16 MiB packet limit.

Owned result packets use a `QSV1` envelope:

```text
magic[4] = "QSV1"
kind: u16 little-endian
version: u16 little-endian
payload_length: u32 little-endian
payload[payload_length]
```

JavaScript must copy a payload before freeing its buffer handle and must recreate
typed-array views after any export that can grow Wasm memory.

`qs_coord_extract` returns one engine-validated nontrivial factor as a
fixed-width little-endian `Natural` payload. The browser recursively factors the
two parts and verifies the complete product before presenting success.

### 11.3 SIMD selection

Release packaging builds:

- a scalar module with no default features;
- a module with `-C target-feature=+simd128` and `wasm-simd128`.

The frontend attempts to compile the SIMD module first and falls back to scalar.
SIMD is intentionally scoped to the XOR-heavy row-reduction kernel. Applying
whole-program Wasm SIMD or Binaryen post-optimization has measured regressions
and is not part of the release pipeline.

## 12. CLI

`qs-factor` is available with the default `cli` feature on native targets.

It:

- reads one unsigned decimal integer from standard input;
- accepts `--threads auto|N`;
- accepts `--progress auto|always|never`;
- prints prime factors in ascending order, one per stdout line, including
  repetitions;
- writes progress and elapsed time only to stderr;
- rejects zero, malformed input, unknown options, and multiple input tokens;
- prints nothing for the successful factorization of one;
- returns a nonzero exit status on failure.

Stdout must remain machine-readable.

## 13. Build, installation, and release packaging

Default native build:

```sh
make
```

Equivalent Cargo build:

```sh
cargo build --release
```

Native installation:

```sh
sudo make install
make install PREFIX=/usr DESTDIR="$staging_root"
```

`PREFIX`, `DESTDIR`, `BINDIR`, `LIBDIR`, `INCLUDEDIR`, and `PKGCONFIGDIR` are
overridable.

Browser build and preview:

```sh
make docs
make serve
```

Release archives:

```sh
SDKROOT=../MacOSX15.4.sdk ./build-release.sh
```

The release script uses cross-rs for Linux and FreeBSD, xwin for MSVC,
cargo-zigbuild for Apple arm64, and native Cargo for Wasm. `SDKROOT` is required
for Apple and may be relative or absolute. Archives are reproducible when
`SOURCE_DATE_EPOCH` is fixed.

## 14. crates.io package contents

The package manifest must explicitly include:

- Rust sources and tests;
- `README.md`, `CHANGELOG.md`, `BENCHMARKING.md`, and this specification;
- both license files;
- the C header and pkg-config template;
- the Makefile and release builder;
- release installer templates;
- the production `web/` frontend;
- the Node/V8 benchmark harness.

Generated `docs/`, release archives, local agent state, the historical audit,
and the obsolete `js/` prototype must not be published.

`cargo package --list` is the authoritative package-content check.
`cargo package` must successfully compile the packaged source before release.

## 15. Safety and invariants

Native code denies unsafe Rust except in the isolated C ABI module. Wasm unsafe
is limited to allocator calls and bounds-checked raw memory views. Every unsafe
block must have a concrete safety comment, and
`unsafe_op_in_unsafe_fn` is denied.

The following invariants are mandatory:

1. accepted relations satisfy their square congruence;
2. combined large-prime cycles have even large-prime parity;
3. matrix indices are bounds-checked during construction/deserialization;
4. every dependency is verified against the original matrix;
5. every returned nontrivial divisor divides the composite being split;
6. every final factor passes probable-prime testing;
7. multiplying factors with multiplicity reconstructs the input;
8. stale Wasm handles and obsolete worker generations are rejected or ignored;
9. malformed C/Wasm inputs produce errors rather than unwinding across an ABI.

The crate is not constant-time and must not be used where operand-dependent
timing reveals a secret.

## 16. Testing and performance gates

Required release checks:

```sh
cargo fmt --all -- --check
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps --all-features
cargo check --locked --all-targets --all-features
cargo package --list --allow-dirty
cargo package --allow-dirty
```

Tests cover:

- randomized arithmetic differential checks against dev-only `num-bigint`;
- fixed arithmetic, modular, primality, and perfect-power cases;
- relation and dependency validity;
- deterministic portable job output;
- out-of-order relation submission;
- native parallel engine factorization;
- browser Worker architecture and malformed glue-protocol handling;
- factor ordering and multiplicity;
- zero/one and invalid input behavior;
- progress completion and cooperative cancellation;
- CLI stdout/stderr behavior;
- C ownership, reuse, null, and bounds behavior.

Performance comparisons must:

- use release builds;
- use fixed, published inputs and verify returned factors;
- report worker/thread counts;
- avoid concurrent builds and benchmarks;
- compare saved binaries in interleaved order when machine load can vary;
- keep timing assertions out of correctness tests.

The fixed 192/224/256-bit balanced-semiprime corpus in `BENCHMARKING.md` is the
regression gate for native time, browser time, Wasm size, and startup. Broader
"general factorizer" claims require the unbalanced/multi-prime corpus and
same-browser competitor protocol described there.

## 17. Current limitations and future work

The 0.3 release is optimized for balanced semiprimes. Its principal known gaps
are:

- no ECM for medium factors in unbalanced composites;
- no true sparse block-Lanczos recurrence for matrices beyond the current
  practical range;
- cache-blocked/bucket tuning still needs representative small-L2 mobile
  measurements;
- parameter tables need multi-input rather than single-sample sweeps at every
  tier;
- retain the existing same-device mobile competitor comparisons as publishable
  raw results so the world-leading browser-SIQS result is independently
  reproducible.

Future optimization must preserve the safety and correctness invariants above
and must not slow the default balanced-RSA artifact.
