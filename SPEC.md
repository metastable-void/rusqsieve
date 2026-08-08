# rusqsieve 0.4 implementation specification

This document specifies the supported behavior and current architecture of
rusqsieve 0.4, release-verified on 2026-08-01. It describes the implementation
that is shipped, not an aspirational module layout or a compatibility promise
for private internals.

Normative requirements use **must**. Descriptions of tuning and implementation
strategy document the 0.4 release and may change in later compatible releases
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

The following are explicitly outside the 0.4 scope:

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

The supported interfaces in 0.4 are:

1. the items re-exported from `src/lib.rs`;
2. the factor/result functions, ABI query, status formatter, and opaque type
   declared in `rusqsieve.h`;
3. the `qs-factor` command-line behavior documented below;
4. the Wasm exports used by `web/abi.js`, `web/index.js`,
   `web/coordinator.js`, and `web/worker.js`.

Everything else in `src/` is private implementation detail even when an
internal item uses Rust's `pub` visibility inside its private module.

Patch releases in the 0.4 series may retune parameters, change internal
representations, add non-exhaustive error/progress variants, or replace private
algorithms. They must not expose invalid relation, matrix, pointer, or scheduler
states through the safe Rust API.

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
├── engine/
│   ├── siqs.rs         polynomial construction and sieve-family kernel
│   ├── extract.rs      verified dependency extraction
│   ├── wire.rs         private worker-family serialization
│   ├── root_simd.rs    x86-64 SSE2 root advancement
│   └── root_wasm.rs    Wasm SIMD128 root advancement
├── qs/mod.rs           factor-base construction and SIQS tier parameters
├── f2/mod.rs           sparse filtering and dependency solving
├── f2/block_lanczos.rs portable 64-way Montgomery block Lanczos
├── natural/mod.rs      fixed-capacity unsigned arithmetic
├── natural/montgomery.rs  Montgomery arithmetic for the rho stage
├── smallfactor.rs      cached small primes and u64 Pollard–Brent
├── u64math.rs          shared machine-word primality/factoring kernels
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
the explicit SIMD128 matrix-XOR and root-advancement kernels. Cargo features are
additive: none changes the identity or default const-generic width of a public
type.

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
- callback cancellation returns `FactorError::Cancelled`;
- input width alone never causes rejection; only a composite that reaches the
  quadratic sieve is range limited, returning
  `FactorError::SiqsCompositeTooLarge`;
- sieving that spends its family budget short of the relation target returns
  `FactorError::InsufficientRelations`, and linear algebra that finds no
  nontrivial dependency returns `FactorError::NoDependency`; neither is
  reported as `FactorError::ResourceLimit`.

The no-observer path must remain separately monomorphized so progress timing and
callback machinery are compiled out of ordinary `factor` and `factor_with`
calls.

### 4.3 Configuration

`FactorConfig` is an owned, encapsulated configuration. In 0.4 its supported
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
factorization accepts any value within the caller's const-generic storage
width; input width alone is never a rejection reason.

The range limit applies to the composite handed to the quadratic sieve, which
is capped at 400 bits. Stages that precede the sieve — trial division, the
primality test, perfect-power detection, and Pollard–Brent — are not width
limited, so a wide input whose factors are small is ordinary work. NFS is
appropriate for hard composites beyond this project's SIQS range.

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

Big-integer Pollard–Brent constructs one Montgomery context for its odd modulus
and runs the whole walk in raw limb buffers. Polynomial values and batched
products remain encoded throughout the stage; reduction uses limb multiplication
and carry propagation rather than division.

`natural/montgomery.rs` is written for that loop specifically. Multiplication is
coarsely integrated operand scanning, so no double-width product is materialized;
squaring uses the symmetric product, which removes about a quarter of the
word-multiplies; the inner loops are monomorphized over the modulus limb count
and unrolled; and every routine reads and writes only the significant limbs of
buffers the caller holds across iterations, so nothing is copied or cleared per
operation. Limb width is chosen per target: 64 bits where a widening 64×64
multiply exists, 32 bits on wasm where it does not and every 128-bit product
would otherwise be emulated. Both widths are generated from one macro and checked
against each other, so the host test suite covers the wasm arithmetic, and both
are checked against division-based modular arithmetic at every limb count.

The iteration budgets are stated in iterations rather than seconds, so this work
reduced the wall clock of a given budget rather than deepening the search: the
stage is 1.6× to 3.1× faster depending on width, and `BENCHMARKING.md` records
the comparison against GMP's mpn assembly, which is what YAFU and FLINT run on.

At and below the sieve's 400-bit ceiling the budget stays a small fraction of
the SIQS run it precedes, because on a balanced input rho contributes nothing
and its cost is pure overhead. Above the ceiling that reasoning does not apply:
the sieve refuses the composite, so rho is the entire attempt and the
alternative to spending more is `SiqsCompositeTooLarge` on an input whose
smallest factor was findable. The budget there is a wall-clock decision instead
— 26 to 36 s per attempt across the supported widths, tiered by width because
per-iteration cost grows with the square of the limb count. Since Brent finds
`p` in roughly `1.2·sqrt(p)` iterations, that reaches a smallest factor of
about 2^53 at 512 bits, 2^51.7 at 768, and 2^50.5 at 1024; factors up to 32
bits are covered by more than two orders of magnitude of margin at every width.
`RUSQSIEVE_RHO_ITERATIONS` overrides the budget for callers who want the
minutes-to-hours search that 56- and 64-bit factors cost.

The sieve-fraction sizing rests on one premise — that the node is a balanced
semiprime the sieve will finish cheaply — and a cofactor that reached the
recursion by *splitting under rho* disproves it: such a value has at least three
prime factors and at least one was small enough for rho to find. From 257 bits
upward, where a sieve run stops costing seconds, those cofactors therefore keep
the deep budget. Balanced semiprimes never reach that branch, because rho does
not split them and nothing below them inherits the mark. Without this rule a
wide product of middling primes peeled factors only while it was above the
ceiling: a 498-bit product of ten 50-bit primes crossed 400 bits after two
splits, and its 399-bit remainder went to a sieve that wanted 206,403 relations
at roughly two per second. It now peels five factors in rho and hands a 250-bit
remainder to a sieve that returns it in under three seconds.

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

The high-level API accepts any `Natural<P>` that fits the engine's fixed
optimized width and converts it without changing the value; a value beyond that
capacity returns `FactorError::InputTooLarge`. There is no private slow
fallback.

Step 6 is the only width-limited stage. A composite wider than 400 bits reaching
it returns `FactorError::SiqsCompositeTooLarge(bits)`, where `bits` is the
composite's width, not the caller's input. Steps 1–5 run regardless of input
width, so a wide input with small factors completes through them — and a
composite past the ceiling gets the deep Pollard–Brent budget described in §5
before step 6 is asked at all, so the error means the smallest factor outran a
minute of rho, not that the input was too wide to try.

Sieving that spends its whole polynomial-family budget without reaching the
relation target returns `FactorError::InsufficientRelations`. The budget scales
with input width. Linear algebra that finds no nontrivial dependency in a matrix
that did meet its target returns `FactorError::NoDependency`. Neither is
reported as a resource limit.

Polynomial-family selection is deterministic for a fixed version, input, and
configuration. Portable/browser sessions merge results by family number.
Native collection deliberately ingests completed unique families in arrival
order to avoid head-of-line stalls behind unusually expensive families.
Correctness does not depend on that order, and the public factor list is always
sorted.

## 7. SIQS engine

### 7.1 Parameters and multiplier

Bit-size tiers in `qs::parameters::engine_params` select the factor-base bound,
sieve half-width, large-prime allowance, and whether double-large-prime
collection is active. Environment overrides are development tuning aids and
are not stable public configuration.

The measured scalar/M4RI browser policy is unchanged below 272 bits. The
281–288 tier uses a 1.4M factor-base prime bound and a 262,144 half-width.
Above that boundary, 96-thread native anchors split the scale at
296/304/312/320 bits, with bounds from 1.2M through 3M. Every high-digit tier
uses a 262,144 half-width; larger per-worker arrays lost to memory pressure.
Native and Wasm both use true sparse linear algebra from 272 bits upward.

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
3. adds `ceil(log2 p)` weights through 288 bits and nearest-integer log weights
   above 288 bits, reducing high-tier survivor false positives; additions wrap
   rather than saturate whenever the engine can prove from the smallest scored
   prime that no position can carry past 255;
4. uses the paired root-difference stride loop;
5. extracts candidates eight positions at a time with one masked compare;
6. trial-divides survivors.

For each survivor, `g(x) = Q(x)/A` is reconstructed directly as a signed value.
Candidate division is gated by a precomputed multiply-shift residue test and
stops as soon as the primes divided out account for the recorded score, which is
exact when the scores did not saturate.

The sieve threshold is `log2|g(x)| − log2(cofactor bound) − small-prime slack`
plus a measured per-tier offset. The cofactor term is the single-large-prime
bound when DLP is off and the bounded DLP product when it is on; it is not an
independent constant. A survivor whose cofactor exceeds what
`classify_cofactor` will accept costs a full trial division for no possible
relation.

### 7.5 Relations and large primes

An accepted relation represents:

```text
t² ≡ (-1)^sign × product(p_i ^ e_i) × large_parts (mod N)
```

where `t = Ax + B` reduced modulo `N`.

Full relations have no large-prime cofactor. Through 280 bits the engine accepts
one prime large factor up to 256 times the factor-base bound; the 281–288
larger-base tier uses 100×. DLP remains off in those measured tiers. Above 288
bits the engine uses a 120× or 150× single-prime bound and admits a product of
two prime large factors up to 12× or 16× the factor-base-bound square. This
measured graph-density cutoff is deliberately below both the full
single-limit square and `single_limit^1.8`: relations containing two sparse
vertices near `single_limit` rarely close useful cycles.

The coordinator treats single partials as edges to a reserved unit vertex and
double partials as ordinary graph edges. A cycle combines relations only when
every large-prime exponent cancels to even parity; the corresponding
square-root factors are retained for extraction.

The portable relation collector buffers out-of-order families and merges them
in ascending family order. Native collection ingests completed unique families
immediately. Native workers and Web Workers execute the same `sieve_family`
kernel through these different schedulers.

## 8. Linear algebra and extraction

Matrix columns correspond to combined relations. Rows correspond to the sign
and factor-base-prime exponent parities.

The 0.4 solver is:

1. sparse structured elimination, including singleton removal and deterministic
   low-weight row elimination through weight six;
2. compact pivot records used to expand residual dependencies back into the
   original relation-column space;
3. a residual solve selected by tier: compact scalar/M4RI row-echelon below
   272 bits, or 64-way Montgomery block Lanczos from 272 bits upward;
4. expansion of at most 64 useful dependencies;
5. verification of every expanded dependency against the original parity
   matrix.

The block-Lanczos path represents 64 vectors in each `u64`, repeatedly applies
the sparse symmetric product `Bᵀ(BV)` without materializing it, selects and
inverts the recurrence's nonsingular 64×64 submatrices, and includes the third
recurrence term required after a rank-deficient block. Deterministic independent
starting blocks recover from the algorithm's small breakdown probability.
Terminal candidates and every lifted dependency are verified against the
original parity matrix.

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
- ingests completed unique families immediately to avoid head-of-line stalls;
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

The reference frontend applies the engine's 400-bit sieve limit to the composite
it sends the coordinator, not to the number the user typed, and reads that limit
from the runtime (`qs_max_siqs_bits`) rather than duplicating it. A composite
above that limit — or one below it that an earlier split already proved
unbalanced, from 257 bits up — is not handed over until a deep Pollard–Brent
search has been spent on it, mirroring the native policy: the opening BigInt peel
is sized for a sieve run that is about to happen, and in these cases either no
sieve run is, or the one that would run is the wrong tool.

That search runs in wasm (`qs_rho`) across a short-lived pool of rho workers, not
on the main thread. Two things follow. The page stays interactive, so the budget
is bounded by patience rather than by frame time; and each worker walks a disjoint
range of polynomial constants, so the pool runs that many independent walks and
the first collision wins — about a `sqrt(T)` speedup for `T` workers. Measured
under Node against the scalar module, wasm runs the loop at 1.08M iterations/s
on a 512-bit modulus and 315k/s at 1024, against 288k/s and 115k/s for the main
thread's `BigInt` implementation. With eight workers and a per-worker budget of
2^25 iterations up to 512 bits, the frontend reaches a smallest factor of roughly
2^52 there and 2^50 at 1024 — parity with the native CLI, against 2^29 for the
opening peel, which cannot even guarantee a 32-bit factor.

A runtime with no compiled module or no `Worker` falls back to the main thread's
sliced `BigInt` search, which yields to the event loop about every 50 ms and
carries a smaller budget for that reason. Per-run
messages, including errors, are generation-scoped. Worker initialization,
individual sieve jobs, and complete runs are bounded by timeouts; a failed run
terminates and rebuilds the Worker runtime. The browser scheduler reads the
engine's family budget from its session (`qs_coord_family_budget`), issues at
most that many families, and reports exhaustion once no jobs or submissions
remain. If the first dependency set
yields no nontrivial factor, the coordinator retains its session, requests a
relation surplus, and retries instead of discarding collected work.

The production frontend is `web/`. The older `js/` prototype is not part of the
published crate or release archives.

### 11.2 Active raw ABI

The browser frontend relies on:

```text
qs_abi_version
qs_max_siqs_bits
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
qs_coord_family_budget
qs_coord_relations
qs_coord_submit
qs_coord_extract
qs_coord_free
```

The Wasm ABI version is 4. Version 3 added `qs_max_siqs_bits` and
`qs_coord_family_budget`, and version 4 adds `qs_rho`; the reference glue
depends on all three rather than duplicating or working around them, so an older
module paired with current glue is rejected at initialization instead of
faulting on a missing export.

`qs_rho(n_pointer, n_length, budget, first_constant, constant_count)` runs a
bounded Pollard–Brent over the decimal composite `n` and returns a packet of
kind 12: the factor as `PARTS * 8` little-endian bytes, or an empty payload when
the budget was spent without a split. `first_constant` and `constant_count`
select which polynomial constants `y^2 + c` the call walks, so a pool of workers
given disjoint ranges runs that many independent walks over one modulus and the
first collision wins. There is no cancellation protocol: the call runs to its
budget and returns, and the frontend cancels by terminating the worker.

Handles contain a 16-bit slot and 16-bit generation. Generation checks reject
ordinary stale-handle reuse; a slot can alias an ancient handle after 65,535
reuse cycles, so the raw ABI does not promise unbounded temporal uniqueness.
Incoming pointers and lengths are checked with checked arithmetic against
current Wasm memory and a 16 MiB packet limit.

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

The frontend attempts to compile the SIMD module first and falls back to
scalar. Explicit SIMD kernels accelerate the XOR-heavy row-reduction path and
SIQS root advancement; all other code retains scalar Rust fallbacks. Applying
additional whole-program Wasm transforms or Binaryen post-optimization has
measured regressions and is not part of the release pipeline.

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
8. stale Wasm handles within the documented generation window and obsolete
   worker generations are rejected or ignored;
9. malformed C/Wasm inputs produce errors rather than unwinding across an ABI.

The crate is not constant-time and must not be used where operand-dependent
timing reveals a secret.

## 16. Testing and performance gates

Required release checks:

```sh
cargo fmt --all -- --check
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-features --all-targets
cargo test --locked --no-default-features
cargo test --locked --profile release-test --test factorization \
  supplied_factorization_corpus_above_128_bits -- --ignored
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps --all-features
cargo check --locked --all-targets --all-features
cargo package --list
cargo package --locked
make test
SDKROOT=/path/to/MacOSX.sdk ./build-release.sh
```

The final command builds the eight supported archives. Release verification
must integrity-check each archive, execute host-compatible CLIs, link the C
smoke program against packaged libraries, and exercise the shipped `docs/`
frontend in a real browser when the required cross toolchains and browser
harness are available.

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

The 0.4 release is optimized for balanced semiprimes. Its principal known gaps
are:

- no ECM for medium factors in unbalanced composites;
- **a world-class, completely opt-in ECM for non-RSA numbers** is the roadmap
  answer to every gap above. Pollard–Brent costs `O(sqrt p)` in the smallest
  factor, so the deep budget above the sieve's ceiling reaches roughly 2^53 at
  512 bits and 2^50 at 1024 and then stops being payable, while the sieve charges
  by the size of `N`. Stage-1/stage-2 ECM is what covers that range, and the
  target is parity with GMP-ECM and YAFU rather than a token implementation. Its
  users are composites the balanced-semiprime path is wrong for: multi-factor
  numbers, unbalanced ones, and anything past 400 bits. Opt-in is a hard
  requirement, not a preference — a non-default feature and a separate
  general-purpose Wasm artifact, with no ECM code or initialization in the
  balanced-RSA artifact, and the fixed 192/224/256-bit corpus retained as an A/B
  gate for runtime, download size, compilation, startup, and code-cache
  footprint. RSA challenge work must not pay a byte or a cycle for it;
- the big-integer rho stage is single-threaded, so the worker pool idles through
  it. Racing independent walks across `T` workers would give roughly a `sqrt(T)`
  speedup, which is the cheap interim answer for 56- and 64-bit factors that
  `RUSQSIEVE_RHO_ITERATIONS` currently buys in minutes to hours;
- a failed factorization discards the factors already peeled. A wide composite
  whose remaining cofactor is hard returns `SiqsCompositeTooLarge` with nothing
  else, even though rho may have split several primes off it first;
- high-digit SIQS tiers still need multi-input wall-time sweeps on representative
  native hosts;
- the 369..=400-bit tier is a single-input yield measurement, not a qualified
  wall-time result. It exists so the range makes steady, reportable progress up
  to the sieve's accepted limit; RSA-110 at 364 bits remains the highest tier
  with a competitive claim;
- retained partials are not memory-bounded. Each holds a fixed-width root and a
  power list, so a run that sieves long enough for the large-prime graph to
  percolate at the highest tiers can reach multiple gigabytes of resident
  partials. The scaled family budget makes such runs reachable where the former
  flat cap ended them early. There is no accounting that would let the engine
  report `FactorError::ResourceLimit` before the allocator fails;
- cache-blocked/bucket tuning still needs representative small-L2 mobile
  measurements;
- parameter tables need multi-input rather than single-sample sweeps at every
  tier;
- retain the existing same-device mobile competitor comparisons as publishable
  raw results so the world-leading browser-SIQS result is independently
  reproducible.

Future optimization must preserve the safety and correctness invariants above
and must not slow the default balanced-RSA artifact.
