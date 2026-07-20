# Quadratic Sieve Rust Crate — Complete Implementation Specification

## 1. Purpose

Implement a production-quality Rust crate named `rusqsieve` that factors positive integers using a highly optimized, portable, parallelized Self-Initializing Quadratic Sieve (SIQS) implementation.

The crate must provide:

1. A fixed-capacity custom unsigned big-integer type optimized for factorization workloads.
2. A high-level blocking factorization API on supported native targets.
3. A scheduler-independent low-level API for parallel sieving and sparse linear algebra over `F_2`.
4. A raw `wasm32-unknown-unknown` C ABI with no `wasm-bindgen` dependency.
5. Independent JavaScript glue that creates Web Workers and executes Rust worker kernels concurrently.
6. A native CLI that reads one positive decimal integer from standard input and prints prime factors, including repetitions, one per line.
7. A structured progress-reporting system covering factor-base construction, relation collection, matrix processing, and linear algebra.

The implementation must be contained in one Cargo package and one Rust library crate. Platform-specific behavior must be selected with `cfg` gates.

The result is intended to be implementation-ready, benchmarkable, testable, and suitable for later architecture-specific optimization.

---

## 2. Scope and expectations

### 2.1 Required factorization strategy

Use SIQS as the main large-composite factorization algorithm.

The high-level factorization pipeline must also include inexpensive preprocessing:

- validation of the input;
- trial division by embedded small primes;
- probable-prime testing using Miller–Rabin;
- perfect-square and perfect-power detection;
- optional bounded Pollard rho for composites with an easily discoverable factor;
- recursive factorization until every returned factor passes the configured probable-prime test.

### 2.2 Practical range

The default `Natural<16>` stores up to 1024 bits. This is a capacity limit, not a claim that arbitrary 1024-bit semiprimes are practically factorable with SIQS.

The implementation must not artificially restrict SIQS to a fixed decimal-digit range, but parameter heuristics, resource limits, and documentation must make clear that SIQS is not intended to replace the Number Field Sieve for very large hard semiprimes.

### 2.3 Non-goals for the initial implementation

The initial implementation does not need:

- General Number Field Sieve.
- ECM.
- Tokio or any asynchronous Rust runtime.
- `wasm-bindgen`.
- Rust-native multithreading inside `wasm32-unknown-unknown`.
- Shared WebAssembly memory or `SharedArrayBuffer` support.
- Arbitrary-precision signed integers as a public type.
- Stable ABI compatibility across unrelated crate major versions.
- Constant-time cryptographic behavior.

---

## 3. Package and target structure

Use one Cargo package with one library target and one CLI binary target.

Recommended layout:

```text
rusqsieve/
├── Cargo.toml
├── README.md
├── LICENSE-APACHE
├── LICENSE-MPL
├── benches/
│   ├── natural.rs
│   ├── sieve.rs
│   └── linear_algebra.rs
├── js/
│   ├── index.js
│   ├── coordinator.js
│   ├── worker.js
│   └── protocol.js
├── src/
│   ├── lib.rs
│   ├── progress.rs
│   ├── primality.rs
│   ├── factors.rs
│   ├── natural/
│   │   ├── mod.rs
│   │   ├── add_sub.rs
│   │   ├── mul.rs
│   │   ├── div.rs
│   │   ├── gcd.rs
│   │   ├── sqrt.rs
│   │   ├── modular.rs
│   │   ├── parse.rs
│   │   └── format.rs
│   ├── factor/
│   │   ├── mod.rs
│   │   ├── config.rs
│   │   ├── coordinator.rs
│   │   ├── preprocessing.rs
│   │   ├── recursive.rs
│   │   └── state.rs
│   ├── qs/
│   │   ├── mod.rs
│   │   ├── config.rs
│   │   ├── parameters.rs
│   │   ├── multiplier.rs
│   │   ├── factor_base.rs
│   │   ├── polynomial.rs
│   │   ├── sieve.rs
│   │   ├── relation.rs
│   │   ├── partials.rs
│   │   ├── matrix.rs
│   │   └── extract.rs
│   ├── f2/
│   │   ├── mod.rs
│   │   ├── dense.rs
│   │   ├── sparse.rs
│   │   ├── filter.rs
│   │   ├── provenance.rs
│   │   └── block_lanczos.rs
│   ├── work/
│   │   ├── mod.rs
│   │   ├── job.rs
│   │   ├── result.rs
│   │   ├── context.rs
│   │   └── kernel.rs
│   ├── native/
│   │   ├── mod.rs
│   │   ├── pool.rs
│   │   └── driver.rs
│   ├── wasm/
│   │   ├── mod.rs
│   │   ├── abi.rs
│   │   ├── handles.rs
│   │   └── wire.rs
│   ├── arch/
│   │   ├── mod.rs
│   │   ├── portable.rs
│   │   ├── x86_64.rs
│   │   ├── aarch64.rs
│   │   └── wasm32.rs
│   └── bin/
│       └── qs-factor.rs
└── tests/
    ├── natural_properties.rs
    ├── primality.rs
    ├── relations.rs
    ├── matrix.rs
    ├── factorization.rs
    ├── determinism.rs
    └── cli.rs
```

### 3.1 Cargo configuration

Use Rust edition 2024.

The library must emit both `rlib` and `cdylib`:

```toml
[lib]
crate-type = ["rlib", "cdylib"]
```

The CLI must be gated behind a feature:

```toml
[features]
default = ["cli"]
cli = []
relation-checks = []
arch-optimized = []
wasm-simd128 = []
reference-qs = []
```

```toml
[[bin]]
name = "qs-factor"
path = "src/bin/qs-factor.rs"
required-features = ["cli"]
```

Recommended profiles:

```toml
[profile.release]
opt-level = 3
lto = "fat"
codegen-units = 1
panic = "abort"

[profile.bench]
opt-level = 3
lto = "thin"
codegen-units = 1
```

Native build:

```sh
cargo build --release
```

WebAssembly build:

```sh
cargo build \
  --release \
  --target wasm32-unknown-unknown \
  --lib \
  --no-default-features
```

### 3.2 Conditional compilation

Use the following target groups:

```rust
#[cfg(any(unix, windows))]
mod native;

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
mod wasm;
```

Do not use only `not(target_arch = "wasm32")` to define native support. The blocking threaded API should be exposed only on explicitly supported native OS families.

The CLI source must contain:

```rust
#[cfg(target_arch = "wasm32")]
compile_error!("the qs-factor CLI is unavailable on wasm32 targets");
```

---

## 4. Crate-level design principles

1. The mathematical core must not create threads.
2. Parallel work must be represented as deterministic bounded jobs.
3. Native threads and Web Workers must execute the same Rust kernels.
4. The factorization coordinator must own global state and user-visible progress.
5. Workers must return results and metrics, not mutate coordinator state.
6. All cross-Wasm ABI data must be serialized; Rust layouts are not an ABI.
7. Wrapping big-integer operators must never silently contaminate mathematical algorithms that require exact arithmetic.
8. Every accepted QS relation must be verifiable through an explicit invariant.
9. Every matrix dependency must be verified before factor extraction.
10. Public APIs must be documented and designed for forward compatibility using `#[non_exhaustive]` where appropriate.

---

## 5. Public crate API

Recommended root exports:

```rust
pub mod f2;
pub mod low_level;
pub mod natural;
pub mod progress;
pub mod qs;

mod factor;
mod factors;
mod primality;
mod work;
mod arch;

#[cfg(any(unix, windows))]
mod native;

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
mod wasm;

pub use factor::{
    FactorConfig,
    FactorError,
    FactorLimits,
    FactorSession,
    Parallelism,
    ProgressAction,
};

pub use factors::PrimeFactors;
pub use natural::{Natural, ParseNaturalError, WideNatural};
pub use progress::*;

#[cfg(any(unix, windows))]
pub use native::{factor, factor_with, factor_with_progress};
```

The `low_level` module should re-export stable work-unit and kernel APIs intended for custom schedulers:

```rust
pub mod low_level {
    pub use crate::f2::{
        BlockLanczos,
        DependencySet,
        MatrixOperation,
        SparseBinaryMatrix,
    };

    pub use crate::qs::{
        prepare_siqs,
        sieve_job,
        FactorBase,
        RawRelation,
        SieveContext,
        SieveScratch,
    };

    pub use crate::work::{
        execute_job,
        JobHeader,
        KernelContexts,
        MatrixMultiplyJob,
        MatrixMultiplyResult,
        SieveJob,
        SieveResult,
        WorkJob,
        WorkResult,
        WorkerScratch,
    };
}
```

---

## 6. `Natural` fixed-capacity big integer

### 6.1 Representation

```rust
#[repr(transparent)]
#[derive(Clone, Eq, PartialEq, Hash)]
pub struct Natural<const PARTS_64: usize = 16> {
    parts: [u64; PARTS_64],
}
```

Representation rules:

- Limb order is little-endian.
- `parts[0]` is the least-significant limb.
- Capacity is exactly `PARTS_64 * 64` bits.
- Unused high limbs must be zero.
- The type represents unsigned integers only.
- `Natural` should not implement `Copy` by default. Explicit cloning avoids accidental copying of large values.

### 6.2 Constants and basic access

```rust
impl<const P: usize> Natural<P> {
    pub const BITS: usize = P * 64;
    pub const ZERO: Self;
    pub const ONE: Self;
    pub const MAX: Self;

    pub const fn from_u64(value: u64) -> Self;
    pub const fn as_parts(&self) -> &[u64; P];
    pub fn as_mut_parts(&mut self) -> &mut [u64; P];

    pub fn is_zero(&self) -> bool;
    pub fn is_one(&self) -> bool;
    pub fn is_even(&self) -> bool;
    pub fn is_odd(&self) -> bool;
    pub fn bit_len(&self) -> usize;
    pub fn trailing_zeros(&self) -> usize;
    pub fn bit(&self, index: usize) -> bool;
}
```

### 6.3 Parsing

Provide a const-compatible decimal parser:

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParseNaturalError {
    Empty,
    InvalidDigit { index: usize, byte: u8 },
    Overflow,
}
```

```rust
impl<const P: usize> Natural<P> {
    pub const fn from_decimal(value: &str) -> Result<Self, ParseNaturalError>;
}
```

Implement `FromStr` by delegating to `from_decimal`:

```rust
impl<const P: usize> core::str::FromStr for Natural<P> {
    type Err = ParseNaturalError;
}
```

Parsing requirements:

- Accept only ASCII decimal digits.
- Reject empty strings.
- Do not accept leading `+` or `-`.
- Leading zeroes are allowed.
- Overflow must be detected exactly.
- Parsing must not allocate.

### 6.4 Compile-time literal macro

Export:

```rust
natural!("966680312498850986629904881784491804947363701071", 16)
```

Recommended implementation:

```rust
#[macro_export]
macro_rules! natural {
    ($value:literal, $parts:literal) => {{
        const VALUE: $crate::Natural<$parts> =
            match $crate::Natural::<$parts>::from_decimal($value) {
                Ok(value) => value,
                Err(_) => panic!("invalid or overflowing Natural literal"),
            };
        VALUE
    }};
}
```

Malformed or overflowing literals must fail at compile time.

### 6.5 Formatting

Implement:

- `Display` as canonical decimal without leading zeroes.
- `Debug` as a capacity-qualified hexadecimal representation.
- `LowerHex` and `UpperHex`.

Suggested `Debug` form:

```text
Natural<16>(0x1234abcd)
```

Decimal formatting should repeatedly divide a scratch copy by `10^19`, storing `u64` chunks.

### 6.6 Byte conversion

```rust
pub fn from_be_bytes(bytes: &[u8]) -> Result<Self, CapacityError>;
pub fn from_le_bytes(bytes: &[u8]) -> Result<Self, CapacityError>;

pub fn write_be_bytes(&self, out: &mut [u8]) -> Result<usize, BufferTooSmall>;
pub fn write_le_bytes(&self, out: &mut [u8]) -> Result<usize, BufferTooSmall>;
```

Also provide fixed-width byte serialization if useful internally.

### 6.7 Operator semantics

Ordinary arithmetic operators use fixed-width wrapping semantics modulo `2^(64 * P)`:

- `Add`, `AddAssign`
- `Sub`, `SubAssign`
- `Mul`, `MulAssign`
- `BitAnd`, `BitAndAssign`
- `BitOr`, `BitOrAssign`
- `BitXor`, `BitXorAssign`
- `Not`
- `Shl`, `ShlAssign`
- `Shr`, `ShrAssign`
- `Div`, `DivAssign`
- `Rem`, `RemAssign`

Preferred operand implementations:

```rust
&Natural + &Natural -> Natural
Natural += &Natural
```

Owned forms may delegate to borrowed forms.

Division by zero must panic for operator traits, matching ordinary Rust integer operator behavior. Fallible methods must return `None` or a typed error.

### 6.8 Exact and overflow-aware arithmetic

Provide:

```rust
pub fn overflowing_add(&self, rhs: &Self) -> (Self, bool);
pub fn overflowing_sub(&self, rhs: &Self) -> (Self, bool);
pub fn overflowing_mul(&self, rhs: &Self) -> (Self, bool);

pub fn checked_add(&self, rhs: &Self) -> Option<Self>;
pub fn checked_sub(&self, rhs: &Self) -> Option<Self>;
pub fn checked_mul(&self, rhs: &Self) -> Option<Self>;

pub fn wrapping_add(&self, rhs: &Self) -> Self;
pub fn wrapping_sub(&self, rhs: &Self) -> Self;
pub fn wrapping_mul(&self, rhs: &Self) -> Self;
```

Mathematical algorithms must use checked, widening, or modular operations as appropriate. Do not rely on wrapping operator behavior internally where overflow would invalidate a proof or invariant.

### 6.9 Widening multiplication

Stable generic const expressions cannot be assumed for `[u64; 2 * P]`. Use:

```rust
#[derive(Clone, Eq, PartialEq)]
pub struct WideNatural<const P: usize> {
    low: [u64; P],
    high: [u64; P],
}
```

Provide:

```rust
pub fn widening_mul(&self, rhs: &Self) -> WideNatural<P>;
pub fn widening_square(&self) -> WideNatural<P>;

impl<const P: usize> WideNatural<P> {
    pub fn low(&self) -> Natural<P>;
    pub fn high(&self) -> Natural<P>;
    pub fn overflowing_narrow(self) -> (Natural<P>, bool);
}
```

### 6.10 Division and number-theoretic methods

Provide:

```rust
pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)>;
pub fn div_rem_u64(&self, divisor: u64) -> Option<(Self, u64)>;

pub fn gcd(&self, rhs: &Self) -> Self;
pub fn extended_gcd(&self, rhs: &Self) -> ExtendedGcdResult<P>;

pub fn sqrt_rem(&self) -> (Self, Self);
pub fn floor_sqrt(&self) -> Self;
pub fn ceil_sqrt(&self) -> Self;
pub fn is_square(&self) -> bool;

pub fn checked_pow_u32(&self, exponent: u32) -> Option<Self>;
pub fn perfect_power(&self) -> Option<(Self, u32)>;
```

A signed internal type may be used privately for extended GCD coefficients, but it need not be public.

### 6.11 Limb kernels

Portable arithmetic primitives should be slice-oriented:

```rust
fn add_n(out: &mut [u64], a: &[u64], b: &[u64]) -> u64;
fn sub_n(out: &mut [u64], a: &[u64], b: &[u64]) -> u64;
fn add_word(out: &mut [u64], word: u64) -> u64;
fn mul_word(out: &mut [u64], a: &[u64], word: u64) -> u64;
fn mul_schoolbook(out: &mut [u64], a: &[u64], b: &[u64]);
fn square_schoolbook(out: &mut [u64], a: &[u64]);
fn shl_bits(out: &mut [u64], input: &[u64], shift: u32) -> u64;
fn shr_bits(out: &mut [u64], input: &[u64], shift: u32) -> u64;
```

Use `u128` as the portable multiply-accumulate primitive.

Begin with schoolbook multiplication, specialized squaring, and normalized long division. Add Karatsuba only after benchmarks show a benefit at relevant operand sizes.

### 6.12 Modular arithmetic

Use a reusable Montgomery context for odd moduli:

```rust
pub struct Montgomery<const P: usize> {
    modulus: Natural<P>,
    n0_inverse: u64,
    r2: Natural<P>,
    one: Natural<P>,
}
```

Required operations:

```rust
impl<const P: usize> Montgomery<P> {
    pub fn new(modulus: Natural<P>) -> Result<Self, MontgomeryError>;
    pub fn encode(&self, value: &Natural<P>) -> Natural<P>;
    pub fn decode(&self, value: &Natural<P>) -> Natural<P>;
    pub fn mul(&self, lhs: &Natural<P>, rhs: &Natural<P>) -> Natural<P>;
    pub fn square(&self, value: &Natural<P>) -> Natural<P>;
    pub fn pow(&self, base: &Natural<P>, exponent: &Natural<P>) -> Natural<P>;
    pub fn inv(&self, value: &Natural<P>) -> Option<Natural<P>>;
}
```

Also provide specialized small-modulus functions for SIQS factor-base work:

```rust
pub fn jacobi_u64(a: u64, n: u64) -> i8;
pub fn legendre_u32(n_mod_p: u32, p: u32) -> i8;
pub fn tonelli_shanks_u32(n_mod_p: u32, p: u32) -> Option<u32>;
```

---

## 7. Primality testing

### 7.1 API

```rust
#[derive(Clone, Debug)]
pub struct PrimalityConfig {
    pub rounds: NonZero<u32>,
    pub witnesses: WitnessPolicy,
}

#[derive(Clone, Debug)]
pub enum WitnessPolicy {
    FirstPrimes,
    Seeded { seed: [u8; 32] },
}
```

```rust
pub fn is_probable_prime<const P: usize>(
    n: &Natural<P>,
    config: &PrimalityConfig,
) -> bool;
```

### 7.2 Requirements

- Handle all values below 4 correctly.
- Reject even composites immediately.
- Perform trial division by a small fixed prime set first.
- Use strong Miller–Rabin tests.
- Witness generation must be reproducible with a configured seed.
- Do not claim deterministic primality for arbitrary-width inputs unless a mathematically sufficient deterministic witness set is used for the exact range.
- Documentation must use the term “probable prime.”

---

## 8. Prime factor result type

```rust
pub struct PrimeFactors<const PARTS_64: usize = 16> {
    map: BTreeMap<Natural<PARTS_64>, NonZero<usize>>,
}
```

Required methods:

```rust
impl<const P: usize> PrimeFactors<P> {
    pub fn new() -> Self;

    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&Natural<P>, NonZero<usize>)>;

    pub fn get(&self, prime: &Natural<P>) -> Option<NonZero<usize>>;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn into_map(self) -> BTreeMap<Natural<P>, NonZero<usize>>;
    pub fn verify_product(&self, original: &Natural<P>) -> bool;
}
```

Mutation should remain crate-private so these invariants hold:

- every key is at least 2;
- every exponent is nonzero;
- every key passes the configured probable-prime policy;
- the represented product equals the original input.

The natural iteration order is ascending.

---

## 9. High-level factorization API

### 9.1 Native blocking API

On `unix` and `windows`:

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

The API should consume `Natural` to avoid an unnecessary initial clone.

### 9.2 Configuration

```rust
#[derive(Clone, Debug)]
pub struct FactorConfig {
    pub parallelism: Parallelism,
    pub primality: PrimalityConfig,
    pub trial_division_limit: u32,
    pub small_factor_method: SmallFactorMethod,
    pub qs: QsConfig,
    pub limits: FactorLimits,
    pub seed: [u8; 32],
    pub progress_reporting: ProgressReportingConfig,
}
```

```rust
#[derive(Clone, Copy, Debug)]
pub enum Parallelism {
    Auto,
    Exact(NonZero<usize>),
}
```

```rust
#[derive(Clone, Debug)]
pub enum SmallFactorMethod {
    None,
    PollardRho {
        attempts: u32,
        iteration_limit: u64,
    },
}
```

`Parallelism::Auto` resolves using `std::thread::available_parallelism()` on native targets.

### 9.3 Resource limits

```rust
#[derive(Clone, Debug)]
pub struct FactorLimits {
    pub max_relations: Option<usize>,
    pub max_partial_relations: Option<usize>,
    pub max_matrix_nonzeros: Option<usize>,
    pub max_memory_bytes: Option<usize>,
    pub max_polynomial_batches: Option<u64>,
    pub max_pollard_rho_iterations: Option<u64>,
}
```

Native-only drivers may additionally enforce wall-clock cancellation, but portable kernels must not depend on clocks.

### 9.4 Recursive pipeline

The coordinator must implement this conceptual algorithm:

```text
factor_node(n):
    if n == 0: error
    if n == 1: return

    divide out configured small primes

    if n == 1: return
    if is_probable_prime(n): insert n; return

    if n is a perfect power a^k:
        factor a recursively
        multiply exponents by k
        return

    try bounded Pollard rho if enabled

    if no factor found:
        factor = SIQS(n)

    recursively factor factor
    recursively factor n / factor
```

Returned factors must be sorted by `Natural` order through `BTreeMap`.

---

## 10. Factorization session and scheduler-independent coordination

### 10.1 Session type

```rust
pub struct FactorSession<const P: usize = 16> {
    // private coordinator state
}
```

Required API:

```rust
impl<const P: usize> FactorSession<P> {
    pub fn new(
        input: Natural<P>,
        config: FactorConfig,
    ) -> Result<Self, FactorError>;

    pub fn phase(&self) -> SessionPhase;
    pub fn progress(&self) -> &ProgressSnapshot;
    pub fn progress_revision(&self) -> u64;

    pub fn advance_local(
        &mut self,
        budget: LocalWorkBudget,
    ) -> Result<AdvanceOutcome, FactorError>;

    pub fn take_jobs(
        &mut self,
        maximum: usize,
    ) -> Result<Vec<WorkJob>, FactorError>;

    pub fn submit(
        &mut self,
        result: WorkResult,
    ) -> Result<SubmitOutcome, FactorError>;

    pub fn is_finished(&self) -> bool;

    pub fn take_factors(
        self,
    ) -> Result<PrimeFactors<P>, FactorError>;
}
```

### 10.2 Session phases

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum SessionPhase {
    Preprocessing,
    BuildingFactorBase,
    RelationCollection,
    CombiningRelations,
    MatrixConstruction,
    MatrixFiltering,
    MatrixSolving,
    FactorExtraction,
    PrimalityTesting,
    Complete,
    Failed,
}
```

### 10.3 Bounded local work

Coordinator-only phases must be incremental:

```rust
#[derive(Clone, Copy, Debug)]
pub struct LocalWorkBudget {
    pub factor_base_candidates: usize,
    pub relations_to_combine: usize,
    pub matrix_items: usize,
    pub elimination_steps: usize,
    pub primality_steps: usize,
}
```

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AdvanceOutcome {
    Progressed,
    NeedsWorkers,
    Complete,
}
```

The WebAssembly coordinator must be able to yield back to JavaScript frequently. No local coordinator method should perform an unbounded multi-second task without returning.

---

## 11. Progress reporting

### 11.1 Principles

- Progress state is owned by `FactorSession`.
- Worker threads and Web Workers return metrics with results.
- Only the coordinator updates `ProgressSnapshot`.
- Frontends decide when and how to render progress.
- Progress snapshots must contain counters, not preformatted strings.
- Exact, estimated, and unknown totals must be distinguished.
- Phase transitions must always increment the revision.

### 11.2 Snapshot

```rust
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct ProgressSnapshot {
    pub revision: u64,
    pub task_id: u64,
    pub input_bits: usize,
    pub phase: ProgressPhase,
    pub amount: ProgressAmount,
    pub detail: ProgressDetail,
}
```

### 11.3 Generic progress amount

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProgressAmount {
    pub completed: u64,
    pub total: ProgressTotal,
    pub unit: ProgressUnit,
}
```

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressTotal {
    Exact(u64),
    Estimated(u64),
    Unknown,
}
```

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressUnit {
    Candidates,
    Primes,
    Polynomials,
    SievePositions,
    Relations,
    MatrixRows,
    MatrixColumns,
    MatrixNonzeros,
    Iterations,
    MatrixProducts,
    Tasks,
}
```

`fraction()` must return `None` for unknown or zero totals. It must not clamp values above 1.0 when an estimate is exceeded.

### 11.4 Progress phases

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ProgressPhase {
    Preprocessing,
    BuildingFactorBase,
    Sieving,
    CombiningRelations,
    BuildingMatrix,
    FilteringMatrix,
    LinearAlgebra,
    ExtractingFactor,
    PrimalityTesting,
    Complete,
}
```

### 11.5 Factor-base progress

```rust
#[derive(Clone, Copy, Debug)]
pub struct FactorBaseProgress {
    pub bound: u32,
    pub searched_through: u32,
    pub primes_tested: u64,
    pub primes_accepted: u64,
    pub nonresidue_primes: u64,
}
```

Primary meter:

```text
completed = searched_through
total = Exact(bound)
unit = Candidates
```

### 11.6 Sieving progress

```rust
#[derive(Clone, Copy, Debug)]
pub struct SievingProgress {
    pub polynomial_families_completed: u64,
    pub polynomials_completed: u64,
    pub sieve_positions_processed: u64,
    pub candidates_tested: u64,

    pub full_relations: u64,
    pub single_large_prime_relations: u64,
    pub double_large_prime_relations: u64,
    pub usable_relations: u64,
    pub target_relations: u64,

    pub active_workers: usize,
    pub outstanding_jobs: usize,
}
```

Primary meter:

```text
completed = usable_relations
total = Estimated(target_relations)
unit = Relations
```

The target may increase after filtering or rank analysis.

### 11.7 Relation-combination progress

```rust
#[derive(Clone, Copy, Debug)]
pub struct RelationProgress {
    pub partial_relations_examined: u64,
    pub partial_relations_total: u64,
    pub graph_vertices: u64,
    pub graph_edges: u64,
    pub cycles_found: u64,
    pub combined_relations: u64,
}
```

### 11.8 Matrix progress

```rust
#[derive(Clone, Copy, Debug)]
pub struct MatrixProgress {
    pub stage: MatrixStage,

    pub original_rows: u64,
    pub original_columns: u64,
    pub original_nonzeros: u64,

    pub current_rows: u64,
    pub current_columns: u64,
    pub current_nonzeros: u64,

    pub items_processed: u64,
    pub items_total: Option<u64>,

    pub singleton_rows_removed: u64,
    pub duplicate_columns_removed: u64,
    pub structured_eliminations: u64,
}
```

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatrixStage {
    CountingParities,
    ConstructingCsr,
    ConstructingCsc,
    RemovingSingletons,
    RemovingDuplicates,
    StructuredElimination,
    Finalizing,
}
```

### 11.9 Linear-algebra progress

```rust
#[derive(Clone, Copy, Debug)]
pub struct LinearAlgebraProgress {
    pub solver: LinearAlgebraSolver,
    pub stage: LinearAlgebraStage,

    pub matrix_rows: u64,
    pub matrix_columns: u64,
    pub matrix_nonzeros: u64,

    pub iteration: u64,
    pub estimated_iterations: Option<u64>,

    pub matrix_products_completed: u64,
    pub matrix_products_estimated: Option<u64>,

    pub dependencies_found: u32,
}
```

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LinearAlgebraSolver {
    DenseGaussian,
    BlockLanczos,
}
```

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LinearAlgebraStage {
    Initializing,
    Iterating,
    RecoveringDependencies,
    VerifyingDependencies,
    Complete,
}
```

### 11.10 Progress detail enum

```rust
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum ProgressDetail {
    None,
    FactorBase(FactorBaseProgress),
    Sieving(SievingProgress),
    Relations(RelationProgress),
    Matrix(MatrixProgress),
    LinearAlgebra(LinearAlgebraProgress),
    Complete(CompleteProgress),
}
```

### 11.11 Native observer contract

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressAction {
    Continue,
    Cancel,
}
```

The observer:

- runs only on the caller/coordinator thread;
- is never invoked by worker threads;
- is never invoked while an internal mutex is held;
- is invoked immediately on phase changes;
- is throttled for repeated same-phase updates;
- receives one final `Complete` snapshot before success;
- can cancel by returning `ProgressAction::Cancel`.

```rust
#[derive(Clone, Copy, Debug)]
pub struct ProgressReportingConfig {
    pub minimum_interval: Duration,
}
```

A 100 ms default is appropriate on native targets.

Portable snapshots must not contain `Instant` or wall-clock timestamps. Frontends calculate rates and ETAs.

---

## 12. SIQS design

### 12.1 Configuration

```rust
#[derive(Clone, Debug)]
pub struct QsConfig {
    pub multiplier: MultiplierChoice,
    pub factor_base_bound: AutoOr<u32>,
    pub sieve_half_width: AutoOr<u32>,
    pub polynomial_batch_size: u32,
    pub large_primes: LargePrimeConfig,
    pub relation_surplus: RelationSurplus,
    pub sieve_score_scale: u8,
    pub candidate_slack: u8,
    pub matrix: MatrixConfig,
}
```

```rust
#[derive(Clone, Copy, Debug)]
pub enum AutoOr<T> {
    Auto,
    Value(T),
}
```

```rust
#[derive(Clone, Debug)]
pub struct LargePrimeConfig {
    pub single_limit: AutoOr<u64>,
    pub double_product_limit: AutoOr<u64>,
    pub enable_double: bool,
}
```

```rust
#[derive(Clone, Copy, Debug)]
pub struct RelationSurplus {
    pub absolute: usize,
    pub percent: u8,
}
```

Parameter heuristics must live in a versioned module such as `qs/parameters.rs`, not be scattered throughout the implementation.

### 12.2 Multiplier selection

Implement a small multiplier search that scores candidate square-free multipliers using:

- resulting residue properties modulo small primes;
- expected factor-base density;
- bit-length growth penalty;
- special handling of powers of two.

If a candidate multiplier shares a nontrivial GCD with `n`, return that GCD immediately as a factor.

### 12.3 Factor base

Do not use a fixed prime table as the complete factor base.

Requirements:

- Embed a compact small-prime table for trial division and bootstrap prime generation.
- Dynamically generate primes up to the selected factor-base bound.
- For each odd prime `p`:
  - if `n mod p == 0`, return `p` as a factor;
  - otherwise retain `p` only if the working value is a quadratic residue modulo `p`;
  - compute and store a modular square root.
- Handle `p = 2` separately.

```rust
pub struct FactorBaseEntry {
    pub prime: u32,
    pub log_prime: u8,
    pub sqrt_n: u32,
}
```

```rust
pub struct FactorBase {
    entries: Box<[FactorBaseEntry]>,
}
```

Factor-base construction must be incremental:

```rust
pub struct FactorBaseBuilder<const P: usize> {
    // private state
}
```

```rust
impl<const P: usize> FactorBaseBuilder<P> {
    pub fn new(
        n: Natural<P>,
        bound: u32,
    ) -> Result<Self, FactorBaseError>;

    pub fn step(
        &mut self,
        candidate_budget: usize,
    ) -> Result<FactorBaseBuildStatus, FactorBaseError>;

    pub fn progress(&self) -> FactorBaseProgress;
    pub fn finish(self) -> Result<FactorBase, FactorBaseError>;
}
```

### 12.4 SIQS context

```rust
pub struct SieveContext<const P: usize> {
    pub(crate) n: Natural<P>,
    pub(crate) working_n: Natural<P>,
    pub(crate) multiplier: u32,
    pub(crate) factor_base: FactorBase,
    pub(crate) polynomial_plan: PolynomialPlan,
    pub(crate) config: QsConfig,
    pub(crate) context_id: u64,
}
```

Preparation API:

```rust
pub fn prepare_siqs<const P: usize>(
    n: &Natural<P>,
    config: &QsConfig,
) -> Result<SieveContext<P>, QsError>;
```

### 12.5 Polynomial generation

Use SIQS polynomial families. Polynomial identifiers must deterministically derive all coefficients from:

- input;
- factor base;
- configured seed;
- family number;
- polynomial index.

Workers must not make independent random choices.

The polynomial representation must support efficient root updates between related polynomials.

### 12.6 Sieving

Use a byte score array:

```rust
pub struct SieveScratch {
    scores: Vec<u8>,
    candidates: Vec<u32>,
    factor_scratch: Vec<(u32, u16)>,
    bucket_storage: Vec<BucketEntry>,
}
```

For each polynomial:

1. Initialize or approximate `log |Q(x)|` scores.
2. Add scaled `log(p)` at both modular roots.
3. Include selected prime-power contributions where profitable.
4. Use direct stepping for small factor-base primes.
5. Use bucket sieving for large-stride primes.
6. Scan candidates above a configurable threshold.
7. Trial-divide candidate values.
8. Classify full, single-large-prime, and double-large-prime relations.
9. Emit metrics and relations.

Do not begin with handwritten SIMD. First implement a correct scalar path with data layouts friendly to auto-vectorization. Add target-specific optimizations only behind dispatch after benchmarks identify real hot spots.

### 12.7 Relation invariant

```rust
pub struct RawRelation<const P: usize> {
    pub square_root: Natural<P>,
    pub sign: bool,
    pub factors: Box<[PrimePower]>,
    pub large_primes: LargePrimePart,
    pub source: RelationSource,
}
```

```rust
pub struct PrimePower {
    pub factor_base_index: u32,
    pub exponent: u16,
}
```

```rust
pub enum LargePrimePart {
    None,
    One(u64),
    Two(u64, u64),
}
```

Every relation must satisfy:

```text
square_root^2 congruent to
    (-1)^sign
    * product(factor_base[index]^exponent)
    * product(large_primes)
mod n
```

Provide:

```rust
pub fn verify_relation<const P: usize>(
    context: &SieveContext<P>,
    relation: &RawRelation<P>,
) -> bool;
```

Enable verification in tests and under the `relation-checks` feature.

### 12.8 Large-prime relation combination

Use a graph:

- vertices are large primes;
- one-large-prime relations connect a distinguished root vertex to the prime;
- two-large-prime relations connect their two primes;
- cycles create combinations in which every large-prime exponent is even.

Preserve provenance rather than immediately multiplying all values.

```rust
pub struct CombinedRelation {
    pub sources: Box<[RelationId]>,
    pub parity: SparseParity,
}
```

Resource-limit graph growth and partial-relation retention.

### 12.9 Relation collection stopping rule

The collection target must be based on:

- matrix row count;
- configured absolute and percentage surplus;
- observed filtering losses;
- any rank deficiencies found later.

If filtering or linear algebra reveals insufficient independent relations, the session must return to relation collection with a higher estimated target and a new work generation.

---

## 13. Parallel work-unit API

### 13.1 Job identity

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JobHeader {
    pub job_id: u64,
    pub generation: u64,
    pub context_id: u64,
}
```

- `job_id` uniquely identifies a job within the session.
- `generation` identifies the current scheduling generation.
- `context_id` identifies the immutable worker context.
- Late results from obsolete generations must be ignored safely.

### 13.2 Jobs

```rust
#[derive(Clone, Debug)]
pub enum WorkJob {
    Sieve(SieveJob),
    MatrixMultiply(MatrixMultiplyJob),
}
```

```rust
#[derive(Clone, Debug)]
pub struct SieveJob {
    pub header: JobHeader,
    pub family: u64,
    pub first_polynomial: u32,
    pub polynomial_count: u32,
}
```

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatrixOperation {
    Matrix,
    Transpose,
}
```

```rust
#[derive(Clone, Debug)]
pub struct MatrixMultiplyJob {
    pub header: JobHeader,
    pub operation: MatrixOperation,
    pub output_start: u32,
    pub output_end: u32,
    pub input: Box<[u64]>,
}
```

### 13.3 Results and metrics

```rust
#[derive(Debug)]
pub enum WorkResult {
    Sieve(SieveResult),
    MatrixMultiply(MatrixMultiplyResult),
    Failed(WorkFailure),
}
```

```rust
#[derive(Clone, Copy, Debug, Default)]
pub struct SieveJobMetrics {
    pub polynomial_families: u64,
    pub polynomials: u64,
    pub sieve_positions: u64,
    pub candidates_tested: u64,
    pub full_relations: u64,
    pub single_large_prime_relations: u64,
    pub double_large_prime_relations: u64,
}
```

```rust
pub struct SieveResult<const P: usize = 16> {
    pub header: JobHeader,
    pub relations: Vec<RawRelation<P>>,
    pub metrics: SieveJobMetrics,
}
```

```rust
#[derive(Clone, Copy, Debug, Default)]
pub struct MatrixJobMetrics {
    pub output_items: u64,
    pub nonzeros_visited: u64,
}
```

```rust
pub struct MatrixMultiplyResult {
    pub header: JobHeader,
    pub output_start: u32,
    pub words: Box<[u64]>,
    pub metrics: MatrixJobMetrics,
}
```

### 13.4 Worker kernel

```rust
pub struct WorkerScratch {
    pub sieve: SieveScratch,
    pub matrix: MatrixScratch,
    pub arithmetic: ArithmeticScratch,
}
```

```rust
pub fn execute_job<const P: usize>(
    contexts: &KernelContexts<P>,
    job: WorkJob,
    scratch: &mut WorkerScratch,
) -> Result<WorkResult, KernelError>;
```

Native worker threads and Wasm worker exports must call this function or its direct sub-kernels.

---

## 14. Native parallel driver

### 14.1 Thread pool

Create long-lived worker threads once per factorization call or session driver.

Do not spawn a thread for each polynomial batch or matrix multiplication.

```rust
struct NativePool {
    workers: Vec<Worker>,
    command_tx: Sender<WorkerCommand>,
    result_rx: Receiver<WorkResult>,
    outstanding: usize,
}
```

Each worker owns `WorkerScratch` and reuses allocated memory.

### 14.2 Driver loop

The blocking driver must:

1. Create the session.
2. Resolve parallelism.
3. Create the worker pool.
4. Run bounded local coordinator work.
5. Drain completed results before blocking.
6. Request jobs up to available worker capacity.
7. Dispatch jobs.
8. Submit results to the session.
9. Render or publish progress on the coordinator thread.
10. Stop workers and return factors or an error.

Callbacks must never run while holding pool or session locks.

### 14.3 Cancellation

Cancellation may come from:

- `ProgressAction::Cancel`;
- resource-limit exhaustion;
- a native deadline implemented by the driver;
- explicit future session APIs.

When cancelled:

- stop creating jobs;
- signal workers to stop after the current bounded job;
- drain or discard late results;
- return `FactorError::Cancelled`.

---

## 15. Sparse linear algebra over `F_2`

### 15.1 Matrix definition

Rows:

- row 0 represents the sign `-1`;
- one row per factor-base prime.

Columns:

- one per usable full or combined relation.

A matrix bit is 1 when the corresponding exponent is odd.

### 15.2 Storage

```rust
pub struct SparseBinaryMatrix {
    rows: u32,
    columns: u32,

    csr_offsets: Box<[u32]>,
    csr_columns: Box<[u32]>,

    csc_offsets: Box<[u32]>,
    csc_rows: Box<[u32]>,

    provenance: Box<[CombinationId]>,
}
```

Maintain both CSR and CSC forms to support efficient `M*x` and `M^T*x` kernels.

Validate index sizes before converting from `usize` to `u32`.

### 15.3 Filtering

Perform, at minimum:

1. Remove zero columns.
2. Repeatedly remove singleton rows and their only columns.
3. Remove duplicate columns.
4. Optionally perform bounded structured elimination on low-weight rows.
5. Preserve provenance for every transformed column.

Use a provenance DAG:

```rust
enum CombinationNode {
    Raw(RelationId),
    Xor(CombinationId, CombinationId),
}
```

Do not store a full original-relation bit vector for every intermediate column.

### 15.4 Block vectors

Represent 64 parallel binary vectors with one `u64` per matrix position:

```rust
pub struct F2BlockVector {
    words: Box<[u64]>,
}
```

### 15.5 Sparse multiplication

```rust
impl SparseBinaryMatrix {
    pub fn mul_m_rows(
        &self,
        input_by_column: &[u64],
        row_range: Range<usize>,
        output_by_row: &mut [u64],
    );

    pub fn mul_mt_columns(
        &self,
        input_by_row: &[u64],
        column_range: Range<usize>,
        output_by_column: &mut [u64],
    );
}
```

Partition by disjoint output ranges. Workers must not need atomics for output accumulation.

### 15.6 Solvers

```rust
#[derive(Clone, Copy, Debug)]
pub enum MatrixSolver {
    Auto,
    DenseGaussian,
    BlockLanczos,
}
```

Dense Gaussian elimination is required for:

- small matrices;
- tests;
- reference behavior;
- debugging block Lanczos.

Production sparse matrices should use block Lanczos.

### 15.7 Block Lanczos state machine

```rust
pub struct BlockLanczos {
    // private state
}
```

```rust
pub enum LanczosRequest<'a> {
    MultiplyM { input: &'a [u64] },
    MultiplyMt { input: &'a [u64] },
    Complete,
}
```

```rust
impl BlockLanczos {
    pub fn begin(matrix: &SparseBinaryMatrix) -> Self;
    pub fn request(&self) -> LanczosRequest<'_>;

    pub fn submit_product(
        &mut self,
        product: &[u64],
    ) -> Result<LanczosProgress, LinearAlgebraError>;

    pub fn dependencies(&self) -> Option<&DependencySet>;
}
```

Implement a published, mathematically sound block-Lanczos recurrence. Do not invent an informal variation.

Every returned dependency must be checked by multiplying it through the original parity matrix.

### 15.8 Factor extraction

For each verified dependency:

1. Reconstruct the selected original relations through provenance.
2. Compute the congruence of squares.
3. Form `x` and `y` modulo `n`.
4. Try `gcd(|x-y|, n)` and `gcd(x+y, n)`.
5. Reject trivial factors 1 and `n`.
6. Continue with more dependencies if necessary.

Use exact or modular arithmetic that cannot silently wrap incorrectly.

---

## 16. WebAssembly architecture

### 16.1 General model

The portable baseline uses independent Wasm instances:

```text
main page
  └── coordinator worker
        ├── compute worker 0
        ├── compute worker 1
        ├── ...
        └── compute worker N
```

Each compute worker:

- instantiates the same `.wasm` module;
- has independent linear memory;
- imports an immutable phase context;
- receives serialized jobs;
- executes the common Rust kernel;
- returns serialized results.

The coordinator worker:

- owns `FactorSession<16>`;
- performs bounded local work;
- exports worker contexts;
- schedules jobs;
- merges results;
- publishes progress;
- returns factors.

### 16.2 Wasm capacity

The raw ABI uses one concrete type:

```rust
pub(crate) const WASM_PARTS_64: usize = 16;
pub(crate) type WasmNatural = Natural<WASM_PARTS_64>;
```

The Rust library API remains generic.

### 16.3 No Rust threads on Wasm

Do not call `std::thread::spawn` on `wasm32-unknown-unknown`.

JavaScript creates Web Workers and selects concurrency using:

```javascript
const parallelism = Math.max(
  1,
  Math.min(
    options.parallelism ?? Number.MAX_SAFE_INTEGER,
    navigator.hardwareConcurrency || 1,
  ),
);
```

Treat this value as a hint and allow an explicit cap.

### 16.4 Raw ABI requirements

- Use only `extern "C"` exports.
- Use Rust 2024 `#[unsafe(no_mangle)]` syntax.
- Use `u32`, `i32`, and carefully documented `u64` values.
- Do not expose Rust structs, slices, `Vec`, `BTreeMap`, or enum layouts.
- Use opaque handles and serialized packets.
- Validate every pointer, length, handle, version, and index.
- Never trust JavaScript-provided memory ranges.

### 16.5 Memory and buffers

Required exports:

```rust
#[unsafe(no_mangle)]
pub extern "C" fn qs_abi_version() -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_alloc(size: u32, align: u32) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_dealloc(
    pointer: u32,
    size: u32,
    align: u32,
);

#[unsafe(no_mangle)]
pub extern "C" fn qs_buffer_pointer(handle: u32) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_buffer_length(handle: u32) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_buffer_free(handle: u32);
```

After an export that may grow Wasm memory, JavaScript must recreate typed-array views before reading data.

### 16.6 Session ABI

```rust
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_new(
    input_pointer: u32,
    input_length: u32,
    config_pointer: u32,
    config_length: u32,
) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_free(session: u32);

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_phase(session: u32) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_advance_local(
    session: u32,
    budget_pointer: u32,
    budget_length: u32,
) -> i32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_export_context(session: u32) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_take_jobs(
    session: u32,
    maximum_jobs: u32,
) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_submit(
    session: u32,
    result_pointer: u32,
    result_length: u32,
) -> i32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_take_factors(session: u32) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_error(session: u32) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_session_progress(session: u32) -> u32;
```

A returned `u32` is an object or buffer handle unless explicitly documented otherwise.

### 16.7 Worker ABI

```rust
#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_context_import(
    pointer: u32,
    length: u32,
) -> u32;

#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_context_free(context: u32);

#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_execute(
    context: u32,
    job_pointer: u32,
    job_length: u32,
) -> u32;
```

The worker ABI must only decode, validate, dispatch, serialize, and return. It must contain no independent factorization policy.

### 16.8 Handle registry

Use generation-checked handles.

```rust
struct Slot<T> {
    generation: u16,
    value: Option<T>,
}
```

Pack generation and slot index into `u32` where practical. Keep typed registries for sessions, contexts, and buffers to prevent type confusion.

### 16.9 Wire format

Every packet begins with:

```rust
#[repr(C)]
struct PacketHeader {
    magic: [u8; 4],
    kind: u16,
    version: u16,
    payload_length: u32,
}
```

Use magic `b"QSV1"` for ABI version 1.

Rules:

- all multibyte fields are little-endian;
- no unaligned Rust pointer casting;
- decode field by field;
- validate payload length and offset arithmetic;
- reject unknown mandatory enum values;
- allow forward-compatible optional trailing data only when version rules explicitly permit it;
- define maximum packet sizes.

### 16.10 JavaScript API

`js/index.js` should expose a high-level function:

```javascript
export async function factor(input, options = {})
```

Accepted input should include decimal strings and byte arrays.

Return a structured result preserving multiplicity, for example:

```javascript
{
  factors: ["2", "2", "3", "5"],
  grouped: [
    { prime: "2", exponent: 2 },
    { prime: "3", exponent: 1 },
    { prime: "5", exponent: 1 },
  ],
}
```

Support:

- `options.parallelism`;
- `options.onProgress(snapshot)`;
- cancellation through `AbortSignal`;
- resource-limit configuration;
- deterministic seed configuration.

The JS coordinator must call progress updates:

- after bounded local work;
- after each submitted worker result;
- periodically while jobs are outstanding.

---

## 17. CLI specification

### 17.1 Behavior

Binary name:

```text
qs-factor
```

The CLI:

1. Reads all of standard input as UTF-8 text.
2. Trims surrounding whitespace.
3. Requires exactly one unsigned decimal integer.
4. Parses it as `Natural<16>`.
5. Factors it using the native blocking API.
6. Prints prime factors in ascending order, repeated according to multiplicity, one per line.
7. Writes progress only to standard error.
8. Writes errors only to standard error.
9. Returns exit status 0 on success and nonzero on failure.

Example:

Input:

```text
360
```

Output on stdout:

```text
2
2
2
3
3
5
```

### 17.2 Edge cases

- Input `0`: error and nonzero exit.
- Input `1`: success with no stdout output.
- Prime input: one line.
- Prime power: repeated identical lines.
- Leading and trailing whitespace: accepted.
- Internal whitespace or multiple numbers: rejected.
- `+123`, `-123`, and `1_000`: rejected.
- Values above 1024 bits: rejected with a clear message.

### 17.3 Progress mode

Support:

```text
--progress auto|always|never
```

Default: `auto`.

- `auto`: show progress only if stderr is a terminal.
- `always`: always show progress.
- `never`: never show progress.

Use `std::io::IsTerminal`.

Do not require a heavy CLI parsing dependency unless justified. A small manual parser is acceptable.

### 17.4 Progress rendering

Progress output must remain on stderr so stdout is machine-readable.

Suggested forms:

```text
building factor base  [================>----------]  62.4%
                      31,842 primes tested, 15,911 accepted
```

```text
sieving               [=============>-------------]  9,421/~15,250 relations
                      48,192 polynomials, 3,811 partials, 16 workers
```

```text
linear algebra        [=========>-----------------]  iteration 118/~275
                      14,208 x 15,463, 1,184,992 nonzeros
```

For unknown totals, use a spinner and counters rather than a fake percentage.

The CLI may calculate elapsed time, throughput, and ETA from snapshot deltas, but these values do not belong in the portable core snapshot.

### 17.5 Terminal hygiene

- Avoid progress control sequences when stderr is not a terminal unless `always` is selected.
- Finish or clear the current progress line before printing an error.
- Finish progress rendering before printing factors.
- Handle broken pipe on stdout gracefully where practical.

---

## 18. Errors

```rust
#[derive(Debug)]
#[non_exhaustive]
pub enum FactorError {
    ZeroHasNoPrimeFactorization,
    CapacityExceeded,
    ResourceLimit(ResourceLimitKind),
    NoNontrivialFactor,
    InsufficientRelations,
    LinearAlgebra(LinearAlgebraError),
    InvalidRelation,
    InvalidDependency,
    Cancelled,
    Stalled,
    InternalInvariant(&'static str),
}
```

All public error types should implement `Display` and `std::error::Error` when `std` is available.

Use typed errors for parsing, arithmetic capacity, SIQS setup, worker decoding, matrix construction, and ABI validation.

Do not use panics for ordinary malformed input, resource exhaustion, or unlucky factorization failure.

Panics are reserved for internal programming errors and operator semantics such as division by zero.

---

## 19. Architecture-specific optimization

### 19.1 Portable baseline

The complete crate must work correctly without architecture-specific features.

### 19.2 Dispatch structure

```rust
mod arch {
    mod portable;

    #[cfg(all(feature = "arch-optimized", target_arch = "x86_64"))]
    mod x86_64;

    #[cfg(all(feature = "arch-optimized", target_arch = "aarch64"))]
    mod aarch64;

    #[cfg(all(feature = "wasm-simd128", target_arch = "wasm32"))]
    mod wasm32;
}
```

### 19.3 Candidate optimizations

Only implement after profiling:

- x86-64 BMI2/ADX Montgomery multiplication;
- AVX2 or AVX-512 XOR and score scanning;
- AArch64 carry-chain and NEON kernels;
- Wasm `simd128` score scanning and XOR;
- cache-blocked sparse matrix traversal;
- improved bucket-sieve layouts;
- prefetching for sparse products.

Every optimized path must have randomized differential tests against the portable implementation.

Runtime feature detection is required where one native binary supports multiple CPU generations.

---

## 20. Determinism

With identical:

- input;
- configuration;
- seed;
- parallelism-independent job generation policy;

The implementation should produce deterministic mathematical results and deterministic job contents.

Relation arrival order may vary under parallel execution. Therefore:

- assign deterministic relation IDs at coordinator acceptance time or canonicalize relation ordering before matrix construction;
- sort or canonicalize partial-relation graph inputs where needed;
- ensure tests do not depend on nondeterministic thread completion order;
- ensure returned prime factors are always ordered through `BTreeMap`.

The exact progress timing need not be deterministic.

---

## 21. Memory management

### 21.1 Reuse

- Reuse sieve arrays per worker.
- Reuse candidate buffers.
- Reuse arithmetic scratch buffers.
- Reuse block-Lanczos vectors.
- Avoid per-prime heap allocation.
- Batch relation submissions to reduce synchronization overhead.

### 21.2 Limits

Before large allocations:

- compute sizes with checked arithmetic;
- enforce configured memory limits;
- reject matrix dimensions exceeding index representation;
- avoid attempting allocations that obviously exceed addressable Wasm memory.

### 21.3 Relation storage

Store parity sparsely and full exponents compactly. Avoid storing redundant large `Natural` values when a relation can be represented by polynomial identity plus a modular square-root contribution.

Preserve enough information to reconstruct the congruence of squares exactly.

---

## 22. Testing requirements

### 22.1 `Natural`

Test:

- exhaustive `Natural<1>` arithmetic against `u64` where feasible;
- randomized `Natural<2>` through `Natural<8>` against `num-bigint` in dev-dependencies;
- carry and borrow chains;
- all-`u64::MAX` limbs;
- multiplication and squaring;
- normalized division edge cases;
- decimal round trips;
- byte-order round trips;
- shifts around 0, 63, 64, capacity-1, capacity, and above capacity;
- GCD and extended GCD identities;
- square root and remainder;
- perfect-power detection;
- Montgomery multiplication and exponentiation;
- exact overflow flags.

`num-bigint` may be a dev-dependency only and must not be used by the production implementation.

### 22.2 Primality

Test:

- values 0 through a large small range against a simple sieve;
- primes and composites near limb boundaries;
- Carmichael numbers;
- strong pseudoprimes to early bases;
- deterministic seeded witness generation;
- known large primes and composites.

### 22.3 Relations

For every generated test relation, verify the relation invariant modulo `n`.

Create deterministic fixtures containing:

- `n`;
- multiplier;
- factor base;
- polynomial identifier;
- coefficients;
- modular roots;
- candidate positions;
- factorizations;
- emitted relations.

Native and Wasm worker kernels should produce byte-identical serialized results for the same fixture where architecture-independent serialization is expected.

### 22.4 Partial relation graph

Test:

- one-large-prime pairing;
- double-large-prime cycles;
- duplicate edges;
- disconnected components;
- cycle provenance;
- resource limits;
- deterministic combination results.

### 22.5 Matrix

Generate random sparse matrices and verify:

```text
M * dependency = 0
```

Compare dense Gaussian and block Lanczos on small matrices.

Include pathological matrices:

- zero rows;
- zero columns;
- duplicate columns;
- all singleton rows;
- rank zero;
- rank one;
- rectangular matrices;
- dimensions 63, 64, and 65;
- extremely sparse rows;
- dense small matrices.

### 22.6 End-to-end factorization

Test:

- prime inputs;
- products of two similarly sized primes;
- highly unbalanced semiprimes;
- repeated prime powers;
- perfect squares;
- Carmichael numbers;
- numbers with many tiny factors;
- values near capacity boundaries;
- `0`, `1`, and `2`;
- deterministic results across parallelism settings.

For every success:

1. Verify every returned factor is probable prime under the configured policy.
2. Multiply factors with exponents using checked/widening arithmetic.
3. Verify exact equality with the original input.

### 22.7 Progress tests

Test:

- revision monotonically increases;
- phase transitions occur in valid order;
- exact totals are never contradicted;
- estimated totals may be raised;
- stale job results do not affect progress;
- cancellation through observer works;
- completion snapshot is emitted exactly once;
- callbacks execute only on the coordinator thread.

### 22.8 CLI tests

Integration tests must cover:

- factorization of 360;
- prime input;
- input 1;
- input 0;
- invalid characters;
- multiple numbers;
- overflow;
- progress disabled when redirected in auto mode;
- progress never contaminates stdout;
- `--progress never` and `--progress always`.

### 22.9 Wasm tests

Provide browser or Node-compatible tests for:

- ABI version;
- allocation and buffer lifetime;
- invalid handles;
- stale handles;
- malformed packets;
- session creation;
- context import;
- worker execution;
- progress packet decoding;
- multi-worker factorization;
- cancellation;
- deterministic result grouping.

---

## 23. Benchmarks

Add benchmarks for:

### 23.1 Arithmetic

- addition and subtraction by limb count;
- schoolbook multiplication;
- squaring;
- division;
- Montgomery multiplication;
- modular exponentiation;
- GCD.

### 23.2 SIQS

- factor-base construction;
- root generation;
- direct sieve stepping;
- bucket sieve;
- candidate scanning;
- candidate trial division;
- relation serialization.

### 23.3 Linear algebra

- CSR `M*x`;
- CSC `M^T*x`;
- filtering;
- dense elimination;
- block-Lanczos iteration.

### 23.4 End-to-end

Use a curated set of composites with increasing difficulty. Record:

- total time;
- relation collection time;
- matrix time;
- memory use where measurable;
- scaling from 1 to available cores.

Benchmarks must not become correctness tests with fragile timing assertions.

---

## 24. Documentation

The crate documentation and README must include:

- what SIQS is;
- practical limitations;
- the distinction between storage capacity and feasible factorization size;
- native usage;
- progress callback usage;
- CLI usage;
- raw Wasm build instructions;
- JS worker usage;
- deterministic seed behavior;
- probable-prime semantics;
- low-level scheduler API;
- safety and resource-limit notes.

Every public item must have rustdoc.

Include at least these examples:

### 24.1 Native factorization

```rust
use rusqsieve::{factor, Natural};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: Natural<16> =
        "966680312498850986629904881784491804947363701071".parse()?;

    let factors = factor(n)?;

    for (prime, exponent) in factors.iter() {
        println!("{prime}^{exponent}");
    }

    Ok(())
}
```

### 24.2 Progress observer

```rust
use rusqsieve::{
    factor_with_progress,
    FactorConfig,
    Natural,
    ProgressAction,
};

let n: Natural<16> = "360".parse()?;

let factors = factor_with_progress(
    n,
    FactorConfig::default(),
    |snapshot| {
        eprintln!("{:?}: {:?}", snapshot.phase, snapshot.amount);
        ProgressAction::Continue
    },
)?;
```

### 24.3 Compile-time literal

```rust
use rusqsieve::{natural, Natural};

const RSA_160: Natural<16> = natural!(
    "966680312498850986629904881784491804947363701071",
    16
);
```

---

## 25. Safety requirements

- Use `#![forbid(unsafe_op_in_unsafe_fn)]`.
- Keep unsafe code isolated in small architecture or ABI modules.
- Document every unsafe block with a concrete safety argument.
- No unsafe code is needed for portable arithmetic.
- Wasm pointer validation must use checked arithmetic and memory bounds.
- Architecture-specific intrinsics must be gated by compile-time or runtime feature checks.
- Do not use `transmute` for wire decoding.
- Do not create references from untrusted Wasm pointers until bounds and alignment are validated.

---

## 26. Code-quality requirements

- Format with `rustfmt`.
- Pass `cargo clippy --all-targets --all-features` with no unjustified warnings.
- Avoid giant functions; separate policy, math kernels, scheduling, and serialization.
- Use explicit typed state machines rather than boolean flag collections.
- Avoid hidden global mutable state.
- No random behavior without an explicit deterministic seed source.
- No terminal or JS concerns inside mathematical modules.
- No threading code inside SIQS or linear-algebra kernels.
- No platform-specific code in public mathematical types.
- Keep dependencies minimal and justified.

Suggested production dependencies: ideally none initially. Development dependencies may include property-testing, benchmarking, and a reference bigint implementation.

---

## 27. Implementation sequence

Codex should implement in the following order and keep the crate compiling after each phase.

### Phase 1: crate skeleton

- Cargo configuration.
- Module layout.
- Public error and config stubs.
- Native and Wasm cfg gates.
- CLI target skeleton.

Acceptance:

- native `cargo check` passes;
- Wasm `cargo check --target wasm32-unknown-unknown --lib --no-default-features` passes.

### Phase 2: `Natural`

- representation;
- parsing;
- formatting;
- comparison;
- bit operations;
- add/subtract;
- shifts;
- widening multiply and square;
- division;
- GCD;
- square root;
- perfect powers;
- Montgomery arithmetic.

Acceptance:

- property tests against reference bigint pass;
- const literal macro works;
- no production bigint dependency.

### Phase 3: primality and result type

- Miller–Rabin;
- deterministic seeded witnesses;
- `PrimeFactors`;
- product verification.

Acceptance:

- pseudoprime and Carmichael tests pass.

### Phase 4: high-level preprocessing

- trial division;
- perfect powers;
- optional Pollard rho;
- recursive factor tree;
- session state machine skeleton;
- progress snapshots.

Acceptance:

- easy composites factor without SIQS;
- progress and cancellation tests pass.

### Phase 5: reference QS

Behind `reference-qs`:

- simple polynomial QS;
- full relations only;
- dense elimination;
- serial execution.

Acceptance:

- small semiprimes factor;
- relations verify.

### Phase 6: SIQS relation collection

- multiplier selection;
- dynamic factor base;
- SIQS polynomial families;
- logarithmic sieving;
- deterministic sieve jobs;
- single-large-prime support;
- relation metrics.

Acceptance:

- serial SIQS factors moderate fixtures;
- relation invariants pass.

### Phase 7: native parallelism

- persistent thread pool;
- job scheduling;
- reusable worker scratch;
- coordinator-only progress callback;
- cancellation.

Acceptance:

- factors are identical for parallelism 1 and N;
- no callback executes on a worker thread;
- speedup is measurable on suitable fixtures.

### Phase 8: double-large-prime graph and matrix filtering

- partial graph;
- cycle extraction;
- sparse parity matrix;
- CSR and CSC;
- singleton filtering;
- duplicate removal;
- provenance DAG.

Acceptance:

- combined relations verify;
- filtered dependency reconstruction is correct.

### Phase 9: block Lanczos

- serial state machine;
- dense oracle comparisons;
- dependency verification;
- factor extraction;
- parallel sparse products.

Acceptance:

- random sparse matrix tests pass;
- end-to-end SIQS succeeds with block Lanczos.

### Phase 10: Wasm ABI and JS glue

- raw ABI;
- generation-checked handles;
- packet codec;
- coordinator session exports;
- worker context and execution exports;
- JS worker pool;
- progress callback;
- cancellation.

Acceptance:

- browser/Node integration tests factor the same fixtures as native;
- malformed packets fail safely;
- multi-worker execution works.

### Phase 11: CLI

- stdin parser;
- factor output;
- progress modes;
- terminal rendering;
- integration tests.

Acceptance:

```sh
printf '360\n' | qs-factor --progress never
```

prints exactly:

```text
2
2
2
3
3
5
```

### Phase 12: optimization

- benchmark;
- profile;
- optimize demonstrated bottlenecks;
- add optional architecture-specific paths;
- retain portable differential tests.

---

## 28. Final acceptance criteria

The implementation is complete when all of the following hold:

1. One Cargo package builds as an `rlib` on native and a raw `cdylib` Wasm module.
2. Native builds expose blocking factorization APIs.
3. `wasm32-unknown-unknown` builds expose only the intended raw ABI and do not attempt Rust thread creation.
4. `Natural<16>` correctly stores and operates on 1024-bit unsigned integers.
5. Compile-time decimal literals reject invalid or overflowing input.
6. High-level factorization returns probable-prime factors in a `BTreeMap` wrapper.
7. Repeated factors are represented through nonzero exponents.
8. SIQS relation collection can run serially or as deterministic parallel jobs.
9. Native execution uses persistent worker threads up to configured parallelism.
10. JavaScript creates independent Web Workers and executes common Rust worker kernels.
11. Sparse matrix construction and filtering preserve relation provenance.
12. Dense elimination works for small matrices.
13. Block Lanczos works for production sparse matrices.
14. Every accepted relation and dependency can be verified.
15. Progress snapshots cover factor-base construction, sieving, relation processing, matrix construction/filtering, and linear algebra.
16. Native progress callbacks run only on the coordinator thread.
17. The CLI prints only factors on stdout and progress/errors on stderr.
18. Input `1` succeeds with empty factor output; input `0` fails.
19. Native and Wasm factorizations produce equivalent grouped prime factors.
20. Tests, clippy, formatting, and documentation checks pass.

---

## 29. Instructions to Codex

Implement the full crate described above, not a toy demonstration.

Important constraints:

- Do not replace the custom `Natural` type with `num-bigint` or another production bigint dependency.
- Do not replace SIQS with trial division or Pollard rho alone.
- Do not use Tokio.
- Do not use `wasm-bindgen`.
- Do not create Rust threads on `wasm32-unknown-unknown`.
- Do not split the implementation into multiple Cargo packages.
- Do not expose Rust object layouts through the Wasm ABI.
- Do not omit relation verification, provenance, or dependency verification.
- Do not fake progress percentages when totals are unknown.
- Do not call user progress observers from worker threads.
- Do not silently use wrapping arithmetic in exact mathematical code.
- Do not add architecture-specific unsafe code before the portable implementation and tests are complete.

When a full optimized implementation cannot be completed in one pass, preserve this architecture and implement coherent compiling phases rather than substituting simplified incompatible APIs.

At each phase:

1. Keep native and Wasm library builds compiling.
2. Add tests before moving to the next subsystem.
3. Document invariants in code.
4. Prefer correctness and verifiability over premature micro-optimization.
5. Leave clear `TODO` markers only for optimizations, not for required correctness behavior.
