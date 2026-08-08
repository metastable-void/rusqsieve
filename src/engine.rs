//! Portable SIQS engine and scheduler-facing work kernels.
mod extract;
#[cfg(target_arch = "x86_64")]
#[allow(unsafe_code)]
mod report_simd;
#[cfg(all(target_arch = "wasm32", feature = "wasm-simd128"))]
#[allow(unsafe_code)]
mod report_wasm;
#[cfg(target_arch = "x86_64")]
#[allow(unsafe_code)]
mod root_simd;
#[cfg(all(target_arch = "wasm32", feature = "wasm-simd128"))]
#[allow(unsafe_code)]
mod root_wasm;
mod siqs;
mod wire;

use crate::f2::{MatrixError, SparseBinaryMatrix};
use crate::factor::FactorTuning;
#[cfg(any(unix, windows, target_arch = "wasm32"))]
use crate::natural::MontgomeryContext;
use crate::qs::{AutoOr, FactorBaseEntry, MultiplierChoice, QsConfig, prepare_factor_base};
use crate::{Natural, PARTS, jacobi_u64};
#[cfg(any(unix, windows))]
use crate::{PrimalityConfig, WitnessPolicy, is_probable_prime};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::Arc;
#[cfg(any(unix, windows))]
use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};
#[cfg(any(unix, windows))]
use std::sync::{Mutex, mpsc};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EnginePhase {
    Preprocessing,
    BuildingFactorBase,
    Sieving,
    LinearAlgebra,
    Extracting,
}

#[derive(Clone, Copy, Debug)]
pub struct EngineProgress {
    pub phase: EnginePhase,
    pub polynomials: u64,
    pub relations: usize,
    pub target: usize,
    pub workers: usize,
}

#[derive(Debug)]
pub enum EngineError {
    Setup(String),
    /// The composite handed to SIQS is wider than [`MAX_SIQS_BITS`].
    SiqsInputTooLarge(usize),
    InsufficientRelations,
    NoFactor,
    /// The filtered matrix admitted no nontrivial dependency. Distinct from [`Self::NoFactor`],
    /// which means dependencies existed but every one produced a trivial gcd.
    NoDependency,
    Worker(String),
    PolynomialSelection(String),
    InvalidDependency,
    ResourceLimit,
    Cancelled,
}
impl fmt::Display for EngineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for EngineError {}

#[derive(Clone)]
struct Context {
    /// The number being factored. Used for root reduction and the extraction gcd.
    n: Natural,
    /// `k·n` for the Knuth-Schroeppel multiplier `k`. The factor base, polynomial roots, and
    /// `Q(x) = (a·x+b)² − k·n` are all built against this; because `k·n ≡ 0 (mod n)`, the
    /// congruence `x² ≡ y² (mod k·n)` still yields a factor of `n` via `gcd(x−y, n)`.
    sieve_n: Natural,
    /// Use `(2Ax+B)²-kN = 4A·g(x)` when `kN ≡ 1 (mod 8)`.
    /// The fixed factor four is a square and gives smaller sieve values.
    use_q2: bool,
    base: Arc<[FactorBaseEntry]>,
    /// Contiguous factor-base primes. Keeping this separate from the 12-byte
    /// setup records removes eight unused bytes from the dominant portable
    /// score stream and also permits four logical lanes per SIMD root load.
    primes: Arc<[u32]>,
    /// Rounded binary-log score per factor-base entry.
    score_weights: Arc<[u8]>,
    /// Lemire fast-mod constant `⌊2^64 / p⌋ + 1` per factor-base prime, precomputed once. Used to
    /// test `x mod p == root` (a ~3-instruction multiply-shift) in trial division without a
    /// hardware divide, so the whole factor base can be gated per survivor cheaply.
    pinv: Arc<[u64]>,
    /// Threshold give-back for the tiny primes the score pass skips, derived from that set.
    small_slack: usize,
    /// `interval mod p` per factor-base prime. Sieve roots are residues of the signed polynomial
    /// coordinate `x`, while score-array positions represent `x + interval`; precomputing this
    /// fixed translation avoids two signed divisions per prime and polynomial in the sieve pass.
    interval_mod_p: Arc<[u32]>,
    interval: i32,
    target_a: Natural,
    a_all: Arc<[usize]>,
    a_pool: Arc<[usize]>,
    a_factor_count: usize,
    /// Sieve-threshold slack for the unfactored cofactor, in bits. This is
    /// `log2(single_limit)` in the browser tiers and `log2(double_limit)` in
    /// high-digit DLP tiers, and nothing else: a survivor whose cofactor exceeds
    /// the applicable bound is discarded by `classify_cofactor`, so admitting
    /// it costs a full trial division for no possible relation.
    /// v0.2.0 used an independent per-tier `lp_allowance` here, which at the 256-bit tier admitted
    /// 34-bit cofactors against a 27-bit acceptance bound — seven bits of pure waste.
    lp_bits: usize,
    /// Measured per-tier sieve-threshold offset in bits (see [`crate::qs::parameters`]).
    thresh_adj: i32,
    /// Maximum accepted single large prime (and maximum factor of a double).
    single_limit: u64,
    /// Maximum accepted product of two large primes, or zero when the
    /// double-large-prime variant is disabled.
    double_limit: u64,
    /// Select the sparse 64-way Montgomery block-Lanczos recurrence after
    /// filtering. The measured crossover is 272 bits on browser Wasm.
    use_block_lanczos: bool,
    relation_percent: Option<usize>,
    small_skip: u32,
    threshold_margin: i32,
    profile: bool,
    /// Bounded native worker count for sparse block-Lanczos matvecs.
    la_threads: usize,
}

/// Large-prime cofactor content of a relation.
#[derive(Clone, Copy)]
enum LargePrime {
    None,
    One(u64),
    Two(u64, u64),
}

/// Fast deterministic hashing for the large-prime graph's machine-word keys.
///
/// `HashMap`'s default SipHash is appropriate for attacker-controlled strings,
/// but relation cofactors are already validated integers and millions of
/// graph insertions make its rounds coordinator-visible. SplitMix64 retains
/// good bucket dispersion without process-random state.
#[derive(Default)]
struct PrimeHasher(u64);
impl Hasher for PrimeHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        let mut value = 0xcbf2_9ce4_8422_2325u64;
        for &byte in bytes {
            value ^= u64::from(byte);
            value = value.wrapping_mul(0x1000_0000_01b3);
        }
        self.0 = value;
    }
    fn write_u64(&mut self, value: u64) {
        let mut mixed = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        self.0 = mixed ^ (mixed >> 31);
    }
}
type PrimeMap = HashMap<u64, u32, BuildHasherDefault<PrimeHasher>>;
impl LargePrime {
    #[inline]
    fn primes(self) -> ([u64; 2], usize) {
        match self {
            LargePrime::None => ([0, 0], 0),
            LargePrime::One(a) => ([a, 0], 1),
            LargePrime::Two(a, b) => ([a, b], 2),
        }
    }
}

#[derive(Clone)]
struct Relation {
    root: Natural,
    sign: bool,
    powers: Vec<(u32, u16)>,
    large: LargePrime,
}

#[derive(Clone)]
struct Column {
    root: Natural,
    sign: bool,
    powers: Vec<(u32, u32)>,
    /// Large primes that were squared out when combining partials; each
    /// contributes once to the reconstructed square root `y`.
    extra_sqrt: Vec<u64>,
}

struct FamilyResult {
    family: u64,
    polynomials: u64,
    relations: Vec<Relation>,
    /// Total sieve survivors examined (read only by the native profiling path).
    #[allow(dead_code)]
    survivors: u64,
    /// Profile-only nanoseconds: family setup, score/scan, survivor
    /// factorization, and Gray-code root advancement.
    #[allow(dead_code)]
    timing: [u64; 4],
}

/// Per-worker reusable buffers (SPEC §7.4 — reuse sieve/candidate scratch).
#[derive(Default)]
struct EngineScratch {
    scores: Vec<u8>,
    /// The two score-array-position residues per factor-base prime for the current polynomial.
    /// These include the fixed `+interval` translation from signed polynomial coordinates.
    /// `root1[i] == u32::MAX` marks a prime that is not directly sieved (2, or a
    /// prime dividing `a`, handled by the per-polynomial linear fallback).
    root1: Vec<u32>,
    root2: Vec<u32>,
    /// Carried roots for the portable cache-blocked dense score prefix.
    dense_root1: Vec<u32>,
    dense_root2: Vec<u32>,
    /// `2·Bⱼ·a⁻¹ mod p` for each varying B-value `j` and factor-base prime `p`
    /// (row-major `[j*nfb + i]`). Adding/subtracting this advances the roots to
    /// the next self-initializing polynomial in O(1) per prime (SPEC §7.3).
    bainv: Vec<u32>,
    /// Positions surviving the score threshold, reused across polynomials.
    candidates: Vec<u32>,
    /// Fixed-stride scored-factor indices for each report. Flat storage avoids
    /// one allocator object per report in the DLP path.
    candidate_factor_counts: Vec<u8>,
    candidate_factors: Vec<u32>,
    /// One-bit report-position filter used by prime-major resieving. At 1 bit
    /// per sieve byte it stays in L2; the old u32 position map was 4 MiB at
    /// RSA-100 and incurred a cache miss on almost every root visit.
    candidate_bits: Vec<u64>,
    /// Reused while trial-dividing reports. Accepted relations clone their
    /// final compact powers; rejected reports and later polynomials allocate
    /// nothing here.
    powers_scratch: Vec<(u32, u16)>,
}

/// Immutable portable SIQS worker context.
#[derive(Clone)]
pub struct EngineContext(Arc<Context>);

/// A deterministic polynomial-family work item.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct EngineJob {
    pub family: u64,
}

/// Relations and metrics produced by one portable work item.
pub struct EngineJobResult {
    inner: FamilyResult,
    pub family: u64,
    pub polynomials: u64,
    pub relations: usize,
}

#[cfg(feature = "fuzzing")]
pub(crate) fn validate_worker_packet(bytes: &[u8]) -> bool {
    wire::deserialize_family(bytes).is_some()
}

/// Pollard-Brent iteration budget for an `n`-bit input, sized as a small fraction of what SIQS on
/// the same input is expected to cost.
///
/// A fixed budget is wrong in both directions. The stage exists for inputs SIQS is bad at — an
/// unbalanced `N`, or a recursive cofactor carrying a small factor — where SIQS pays for the size of
/// `N` while rho pays for the size of its smallest factor. On a *balanced* semiprime rho contributes
/// nothing, so its cost there is pure overhead and must scale with the run it is attached to.
///
/// The brief this work follows expected rho to beat SIQS outright from 65 to 100 bits. Measured on
/// balanced corpus semiprimes (release, single-threaded, seconds — rho with a 16 M-iteration budget
/// against SIQS alone): 70-bit 0.03 / 0.03, 80-bit 0.56 / 0.03, 90-bit 2.16 / 0.01. Rho costs
/// `O(sqrt p)` in the smallest factor while SIQS at those sizes is already trivial, so the premise
/// only holds at the very bottom of that band and the unbounded stage was a 19-216× regression above
/// it. Hence a budget rather than a bit-length gate.
///
/// `factor_base_bound × sieve_half_width` from the tier table is used as the SIQS cost proxy; it
/// tracks measured sieve time within about 3× across 90-256 bits. The divisor was calibrated when
/// each iteration used division-based `Natural::mul_mod`, holding the stage near 1% of the estimated
/// sieve: 1 024 iterations (the floor) up to 128 bits, ~18 k at 192, ~131 k at 224, ~328 k at 256.
/// Real Montgomery REDC subsequently made the repeated modular arithmetic cheaper. The iteration
/// counts deliberately remain unchanged: this reduces unsuccessful-rho overhead on balanced inputs
/// without silently spending the saving on a deeper search.
///
/// Brent finds a factor `p` in roughly `1.2·sqrt(p)` iterations, so this covers factors up to about
/// 2^26 at 192 bits and 2^34 at 256 — above the 10^4 trial-division bound and below the 2^64 point
/// where the whole input would have taken the machine-word path. Raising it further buys a narrow
/// band of larger factors for a quadratically larger cost.
///
/// Above [`MAX_SIQS_BITS`] none of that reasoning holds, and [`wide_rho_budget`] takes over.
///
/// Neither does it hold for a cofactor that reached this point *by splitting under rho*, which is
/// what `after_split` marks. The whole argument for a sieve-fraction budget is that the node is a
/// balanced semiprime the sieve will finish cheaply; a node whose parent just split under rho is
/// provably not that — it has at least three prime factors and at least one of them was small
/// enough for rho to find. Such a node gets the deep budget from [`DEEP_RHO_MIN_BITS`] upward, which
/// is where the sieve stops being cheap. Balanced semiprimes never reach that branch: rho does not
/// split them, so nothing below them ever inherits the flag.
///
/// This is the difference between finishing and not. A 498-bit product of ten 50-bit primes peeled
/// two factors while it was above the ceiling and handed the 399-bit remainder to SIQS, which wanted
/// 206 403 relations at about two per second — weeks of sieving on a number whose every factor rho
/// finds in seconds. Deepening the cofactor's budget instead peels them.
///
/// `tuning.rho_iterations` overrides every arm outright.
fn rho_budget(bits: usize, tuning: &FactorTuning, after_split: bool) -> u64 {
    if let Some(iterations) = tuning.rho_iterations {
        return iterations.max(1);
    }
    if bits > MAX_SIQS_BITS {
        return wide_rho_budget(bits);
    }
    let p = crate::qs::parameters::engine_params(bits);
    let sieve_fraction =
        (p.factor_base_bound as u64 * p.sieve_half_width as u64 / 500_000).max(1_024);
    if after_split && bits >= DEEP_RHO_MIN_BITS {
        return sieve_fraction.max(wide_rho_budget(bits));
    }
    sieve_fraction
}

/// Polynomial constants the native ladder tries in turn, `y^2 + c`.
///
/// Each is an independent walk over the same modulus, which is why the browser can hand different
/// ones to different workers and race them; the native stage runs them in sequence under one shared
/// budget. Eight is where the sequence stops paying: a constant that has not collided by its share
/// of the budget is not more likely to than a fresh one, and every additional constant dilutes the
/// budget the earlier ones get.
pub(crate) const NATIVE_RHO_CONSTANTS: core::ops::RangeInclusive<u64> = 1..=8;

/// Width from which a cofactor known to be unbalanced is worth a deep rho rather than a sieve.
///
/// Below this the sieve is simply the better tool whatever the input's shape: measured native SIQS
/// wall time is 0.38 s at 192 bits, 2.9 s at 224 and 6.8 s at 256, so spending up to a minute of
/// single-threaded rho first could only lose. From here it climbs steeply — 16 s at 272, 38 s at
/// 288, minutes through the RSA-100 and RSA-110 tiers, and out of reach in the 369..=400 tier —
/// while the deep rho budget stays flat at about a minute and finds a factor of the size that was
/// just demonstrated to exist.
const DEEP_RHO_MIN_BITS: usize = 257;

/// Pollard-Brent iteration budget for a composite the sieve will refuse.
///
/// Above [`MAX_SIQS_BITS`] rho stops being a cheap peel in front of SIQS and becomes the entire
/// factoring attempt: there is no sieve run for the budget to be a fraction of, and the alternative
/// to spending more here is [`EngineError::SiqsInputTooLarge`] on an input whose smallest factor was
/// findable. So this arm is a wall-clock decision rather than a proportion.
///
/// Per-iteration cost grows with the square of the limb count, so a flat iteration count would cost
/// three times as much wall clock at the top of the supported range as at the bottom. Measured
/// single-thread rates (release, x86-64 Xeon 8259CL, `profile_wide_rho_throughput`): 4.84 M/s at 512
/// bits, 2.20 M/s at 768, 1.34 M/s at 1024. These tiers therefore spend a comparable amount of time
/// at every width instead — 26 s, 33 s and 36 s.
///
/// Brent finds `p` in about `1.2·sqrt(p)` iterations, so that buys a smallest factor of roughly 2^53
/// at 512 bits, 2^51.7 at 768 and 2^50.5 at 1024. The sieve-derived budget this replaces was 6.29 M
/// iterations at every width above the ceiling — reach 2^44.6 — which refused a 512-bit input
/// carrying a 48-bit factor after 2.7 s of work that was nearly deep enough.
///
/// The budgets are stated in iterations, not seconds, so the arithmetic work that made the stage 1.8×
/// to 5.2× faster bought back wall clock rather than reach: the same tiers cost half to two-thirds of
/// what they did when they were chosen. Reach is where the default stops, because each additional factor bit
/// doubles the cost: 2^56 is minutes and 2^64 is tens of minutes across this width range. That is a
/// real search, not a default one, and callers who want it ask with `RUSQSIEVE_RHO_ITERATIONS`.
const fn wide_rho_budget(bits: usize) -> u64 {
    if bits <= 512 {
        128_000_000
    } else if bits <= 768 {
        72_000_000
    } else {
        48_000_000
    }
}

/// Which dispatch arm produced each split, so tests can assert that a stage ran rather than merely
/// that the answer came out right. A stage that silently never executes passes every correctness
/// test, so the cheap ladder stages are counted explicitly.
/// Counters are thread-local, not global: the dispatch ladder runs entirely on the caller's thread
/// (workers are spawned only *inside* SIQS), and `cargo test` runs test functions concurrently, so
/// process-wide counters would see other tests' factorizations.
#[cfg(test)]
pub(crate) mod stage_counts {
    use std::cell::Cell;
    thread_local! {
        pub(crate) static RHO: Cell<usize> = const { Cell::new(0) };
        pub(crate) static SIQS: Cell<usize> = const { Cell::new(0) };
    }
    pub(crate) fn bump(counter: &'static std::thread::LocalKey<Cell<usize>>) {
        counter.with(|c| c.set(c.get() + 1));
    }
    pub(crate) fn reset() {
        RHO.with(|c| c.set(0));
        SIQS.with(|c| c.set(0));
    }
    pub(crate) fn rho() -> usize {
        RHO.with(Cell::get)
    }
    pub(crate) fn siqs() -> usize {
        SIQS.with(Cell::get)
    }
}

/// Widest composite this engine will attempt with SIQS.
///
/// 400 bits is a little over 120 decimal digits, which is the practical edge of the quadratic
/// sieve: past it GNFS wins by margins no amount of sieve tuning recovers. The cap is deliberately
/// applied to the composite *handed to SIQS*, not to the caller's input, so an arbitrarily wide
/// number whose factors are small still factors normally — only the hard cofactor is bounded.
pub const MAX_SIQS_BITS: usize = 400;

/// Polynomial families any one scheduler will issue before giving up.
///
/// Both the native thread scheduler and [`EngineSession`] use this single bound; v0.2.0 had two
/// uncoordinated 100 000 caps for the same job. Reaching it means the relation target was not met
/// and surfaces as [`EngineError::InsufficientRelations`] rather than a silent wrong answer.
///
/// The budget has to scale with input width. Relations arrive at a roughly constant rate per
/// family, but the *target* grows with the factor base and most of the late yield comes from
/// large-prime cycles, whose count grows superlinearly and therefore only pays off deep into a
/// run. A flat 100 000 was sized for the ≤288-bit tiers and silently truncated everything above
/// it: a 399-bit composite exhausted the whole budget at 13% of its relation target. Measured on
/// a 384-bit semiprime at the 369..=400 tier, a complete run needs on the order of 30 000
/// families; the tiers below leave roughly an order of magnitude of headroom over that.
const fn family_budget(bits: usize) -> u64 {
    match bits {
        // Performance-qualified browser and native tiers. Unchanged: every one of these reaches
        // its target in far fewer families, so raising the ceiling would only mask a regression.
        ..=288 => 100_000,
        289..=368 => 250_000,
        _ => 750_000,
    }
}

/// Prepare an immutable context without creating threads.
pub fn prepare(n: Natural, tuning: &FactorTuning) -> Result<EngineContext, EngineError> {
    prepare_with_la_threads(n, tuning, 1)
}

fn prepare_with_la_threads(
    n: Natural,
    tuning: &FactorTuning,
    la_threads: usize,
) -> Result<EngineContext, EngineError> {
    let input_bits = n.bit_len();
    // Refuse the job here rather than at the caller's input width: this is the one choke point
    // every scheduler (native, WASM coordinator, WASM worker) passes through, and it sees the
    // composite that actually reaches the sieve rather than whatever the caller started with.
    if input_bits > MAX_SIQS_BITS {
        return Err(EngineError::SiqsInputTooLarge(input_bits));
    }
    let mut p = crate::qs::parameters::engine_params(input_bits);
    if let Some(bound) = tuning.factor_base_bound {
        p.factor_base_bound = bound;
    }
    if let Some(half_width) = tuning.sieve_half_width {
        p.sieve_half_width = half_width;
    }
    if let Some(large_prime_multiplier) = tuning.large_prime_multiplier {
        p.large_prime_mult = large_prime_multiplier;
    }
    let high_digit_tier = input_bits >= 289;
    let k = if high_digit_tier {
        knuth_schroeppel(&n)
    } else {
        knuth_schroeppel_legacy(&n)
    };
    let sieve_n = n
        .checked_mul(&Natural::from_u64(k))
        .unwrap_or_else(|| n.clone());
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(p.factor_base_bound),
        multiplier: MultiplierChoice::Value(k as u32),
    };
    let prepared = prepare_factor_base(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
    let primes: Arc<[u32]> = base.iter().map(|entry| entry.prime).collect();
    let pinv: Arc<[u64]> = base.iter().map(|e| lemire_c(e.prime)).collect();
    // Expected log weight of the tiny primes that are skipped by the score pass: `Σ log(p)/(p−1)`
    // over the skipped set, which is what the threshold must give back. Derived once here rather
    // than carried as a hand-tuned constant, and cheap to keep out of the per-polynomial path.
    let small_skip = tuning
        .small_skip
        .unwrap_or(if input_bits >= 289 { 500 } else { 100 });
    let rounded_scores = input_bits >= 289;
    let small_slack = small_slack(&base, small_skip, rounded_scores);
    let score_weights: Arc<[u8]> = base
        .iter()
        .map(|entry| score_weight(entry.prime, rounded_scores))
        .collect();
    let interval_mod_p: Arc<[u32]> = base.iter().map(|e| p.sieve_half_width % e.prime).collect();
    let use_q2 = high_digit_tier && sieve_n.mod_u64(8) == 1;
    let target_source = if use_q2 {
        sieve_n
            .checked_mul(&Natural::from_u64(2))
            .unwrap_or_else(|| sieve_n.clone())
    } else {
        sieve_n.clone()
    };
    let mut target_a = target_source
        .floor_sqrt()
        .div_rem_u64(p.sieve_half_width as u64)
        .unwrap()
        .0;
    if use_q2 {
        target_a >>= 1;
    }
    let (a_all, a_pool, a_factor_count) = siqs::build_a_candidates(&base, &target_a);
    let (single_limit, mut double_limit) = large_prime_policy(
        p.factor_base_bound,
        p.large_prime_mult,
        if p.double_large_primes {
            p.double_large_prime_mult
        } else {
            0
        },
    );
    if let Some(bound) = tuning.double_large_prime_bound {
        double_limit = bound.min(single_limit.saturating_mul(single_limit));
    }
    let context = Arc::new(Context {
        n,
        sieve_n,
        use_q2,
        base,
        primes,
        score_weights,
        pinv,
        small_slack,
        interval_mod_p,
        interval: p.sieve_half_width as i32,
        target_a,
        a_all,
        a_pool,
        a_factor_count,
        // The double cutoff is deliberately below `single_limit²`: relations
        // with two primes near the full single-prime cutoff almost never close
        // useful cycles, but admitting them floods survivor trial division.
        lp_bits: (64
            - if double_limit == 0 {
                single_limit
            } else {
                double_limit
            }
            .leading_zeros()) as usize,
        thresh_adj: p.thresh_adj + tuning.threshold_adjustment.unwrap_or(0),
        single_limit,
        double_limit,
        use_block_lanczos: input_bits >= 272,
        relation_percent: tuning.relation_percent,
        small_skip,
        threshold_margin: tuning.threshold_margin.unwrap_or(0),
        profile: tuning.profile,
        // Sparse matvec bandwidth flattened at eight workers on the retained
        // 96k-column RSA-110 matrix; 16 lost to spawn and cache contention.
        la_threads: la_threads.clamp(1, 8),
    });
    // A famine here is a parameter-selection failure, not a search failure: if no `A` can be built
    // for family 0 then none can be built for any family, and every scheduler would otherwise burn
    // through its whole family budget producing nothing. Diagnose it once, at the only choke point
    // both the native and the WASM/session schedulers pass through.
    if siqs::choose_a(&context, 0).is_none() {
        let message = format!(
            "polynomial-coefficient selection has no viable A for {}-bit input \
             (factor base {}, target A {} bits, {} candidate primes)",
            context.n.bit_len(),
            context.base.len(),
            context.target_a.bit_len(),
            context.a_pool.len(),
        );
        return Err(EngineError::PolynomialSelection(message));
    }
    Ok(EngineContext(context))
}

/// Large-prime acceptance is independent from the sieve threshold slack.
///
/// Cap DLP products at a small multiple of the factor-base square. This keeps
/// both factors in the dense part of the large-prime graph; the conventional
/// `single_limit^1.8` window flooded this engine with unique vertices.
fn large_prime_policy(
    bound: u32,
    large_prime_mult: u32,
    double_large_prime_mult: u32,
) -> (u64, u64) {
    let single_limit = (bound as u64).saturating_mul(large_prime_mult as u64);
    let double_limit = (bound as u64)
        .saturating_mul(bound as u64)
        .saturating_mul(double_large_prime_mult as u64)
        .min(single_limit.saturating_mul(single_limit));
    (single_limit, double_limit)
}

/// Execute a job using only the caller's thread and owned scratch memory.
pub fn execute(context: &EngineContext, job: EngineJob) -> EngineJobResult {
    let mut scratch = EngineScratch::default();
    let inner = siqs::sieve_family(&context.0, job.family, &mut scratch);
    EngineJobResult {
        family: inner.family,
        polynomials: inner.polynomials,
        relations: inner.relations.len(),
        inner,
    }
}

/// Scheduler-independent relation collector. Jobs may finish out of order;
/// submission is merged deterministically by family number.
pub struct EngineSession {
    context: EngineContext,
    target: usize,
    budget: u64,
    next_job: u64,
    next_merge: u64,
    polynomials: u64,
    collector: RelationCollector,
    buffered: BTreeMap<u64, FamilyResult>,
    seen_a: HashSet<Natural>,
}
impl EngineSession {
    pub fn new(context: EngineContext) -> Self {
        let target = relation_target(context.0.base.len(), context.0.relation_percent);
        let budget = family_budget(context.0.n.bit_len());
        Self {
            context,
            target,
            budget,
            next_job: 0,
            next_merge: 0,
            polynomials: 0,
            collector: RelationCollector::new(),
            buffered: BTreeMap::new(),
            seen_a: HashSet::new(),
        }
    }
    /// Polynomial families this session will issue in total. Schedulers that assign family
    /// numbers themselves (the WASM coordinator) must not issue beyond this.
    pub fn family_budget(&self) -> u64 {
        self.budget
    }
    /// Whether every family has been issued without reaching the relation target. A caller that
    /// sees this should report exhaustion rather than attempt extraction.
    pub fn budget_exhausted(&self) -> bool {
        self.next_job >= self.budget && !self.is_ready()
    }
    /// Hand out up to `maximum` polynomial families to sieve.
    ///
    /// Returns fewer than `maximum` jobs — possibly none — once the family budget is spent, which a
    /// caller distinguishes from "done" with [`EngineSession::is_ready`]. Families whose `A`
    /// duplicates one already issued are dropped at ingest rather than here, so a caller that
    /// assigns family numbers itself (as the WASM coordinator does) is equally protected.
    pub fn take_jobs(&mut self, maximum: usize) -> Vec<EngineJob> {
        if self.is_ready() {
            return Vec::new();
        }
        let mut jobs = Vec::with_capacity(maximum);
        while jobs.len() < maximum && self.next_job < self.budget {
            jobs.push(EngineJob {
                family: self.next_job,
            });
            self.next_job += 1;
        }
        jobs
    }
    pub fn submit(&mut self, result: EngineJobResult) {
        self.buffered.insert(result.family, result.inner);
        self.drain_buffered();
    }
    /// Submit a worker's serialized [`EngineJobResult`] (see [`EngineJobResult::to_bytes`]).
    /// Returns whether enough relations have now been collected. Used by the WASM/Web-Worker
    /// scheduler to feed relations sieved in other threads back into the coordinator.
    pub fn submit_bytes(&mut self, bytes: &[u8]) -> bool {
        if let Some(fr) = wire::deserialize_family(bytes) {
            self.buffered.insert(fr.family, fr);
            self.drain_buffered();
        }
        self.is_ready()
    }
    fn drain_buffered(&mut self) {
        while let Some(r) = self.buffered.remove(&self.next_merge) {
            self.next_merge += 1;
            // Only the first family to produce a given `A` contributes. Two families that pick the
            // same `A` sieve identical polynomials, so their relations are identical too, and
            // ingesting both puts duplicate columns in the matrix — every dependency those form is
            // trivial (`x ≡ ±y`), so extraction reports "no factor" on an input that factors.
            //
            // This has to happen at ingest, not at dispatch: the WASM coordinator numbers families
            // itself and never calls `take_jobs`, so filtering there left the browser path
            // unprotected. A 110-bit semiprime that the native path factors in 14 ms produced 3
            // duplicate families out of 56 and failed outright in the browser.
            let unique = siqs::choose_a(&self.context.0, r.family)
                .is_some_and(|(a, _)| self.seen_a.insert(a));
            if !unique {
                continue;
            }
            self.polynomials += r.polynomials;
            let n = &self.context.0.n;
            for rel in r.relations {
                self.collector.ingest(rel, n);
            }
        }
    }
    pub fn is_ready(&self) -> bool {
        self.collector.columns.len() >= self.target
    }
    pub fn relations(&self) -> usize {
        self.collector.columns.len()
    }
    pub fn target(&self) -> usize {
        self.target
    }
    pub fn polynomials(&self) -> u64 {
        self.polynomials
    }
    pub fn extract_factor(&self) -> Result<Natural, EngineError> {
        // Running the linear algebra on a matrix with far fewer columns than factor-base rows
        // cannot produce a dependency, and the failure it does produce says nothing useful about
        // why. Report the real condition instead.
        if !self.is_ready() {
            return Err(EngineError::InsufficientRelations);
        }
        extract::extract(&self.context.0, &self.collector.columns)
    }
}

fn relation_target(base_len: usize, percent: Option<usize>) -> usize {
    if let Some(percent) = percent {
        // Fewer than one relation per factor-base row cannot reliably produce
        // dependencies after singleton filtering.
        return (base_len * percent.clamp(100, 110) / 100).max(base_len + 1);
    }
    base_len + 64
}

#[cfg(any(unix, windows))]
pub fn factor(
    mut n: Natural,
    threads: usize,
    tuning: &FactorTuning,
    witness_seed: Option<[u8; 32]>,
    mut progress: impl FnMut(EngineProgress) -> bool,
) -> Result<Vec<Natural>, EngineError> {
    if n.is_zero() {
        return Err(EngineError::Setup("zero has no prime factorization".into()));
    }
    let mut primality = PrimalityConfig::default();
    if let Some(seed) = witness_seed {
        primality.witnesses = WitnessPolicy::Seeded { seed };
    }
    let mut factors = Vec::new();
    for &p in crate::smallfactor::small_primes() {
        if p > 10_000 {
            break;
        }
        loop {
            let (q, r) = n.div_rem_u64(p as u64).unwrap();
            if r != 0 {
                break;
            }
            factors.push(Natural::from_u64(p as u64));
            n = q
        }
    }
    if n.is_one() {
        return Ok(factors);
    }
    factor_node(
        n,
        threads.max(1),
        &primality,
        tuning,
        false,
        &mut progress,
        &mut factors,
    )?;
    factors.sort();
    Ok(factors)
}

/// Factor one node of the recursion. `after_split` records that this value is a cofactor of a
/// composite rho already split — the evidence that sizes its own rho budget; see [`rho_budget`].
#[cfg(any(unix, windows))]
#[allow(clippy::too_many_arguments)]
fn factor_node(
    n: Natural,
    threads: usize,
    pc: &PrimalityConfig,
    tuning: &FactorTuning,
    after_split: bool,
    progress: &mut impl FnMut(EngineProgress) -> bool,
    out: &mut Vec<Natural>,
) -> Result<(), EngineError> {
    if !progress(EngineProgress {
        phase: EnginePhase::Preprocessing,
        polynomials: 0,
        relations: 0,
        target: 0,
        workers: threads,
    }) {
        return Err(EngineError::Cancelled);
    }
    if n.is_one() {
        return Ok(());
    }
    // Native machine-word fast path: everything up to 64 bits is factored with
    // deterministic Miller-Rabin + Pollard-Brent in `u64`/`u128`, bypassing
    // fixed-capacity big-integer arithmetic entirely.
    if let Some(v) = n.to_u64() {
        let mut small = Vec::new();
        let completed = crate::smallfactor::factor_u64_cancellable(v, &mut small, || {
            !progress(EngineProgress {
                phase: EnginePhase::Preprocessing,
                polynomials: 0,
                relations: 0,
                target: 0,
                workers: threads,
            })
        });
        if !completed {
            return Err(EngineError::Cancelled);
        }
        out.extend(small.into_iter().map(Natural::from_u64));
        return Ok(());
    }
    if is_probable_prime(&n, pc) {
        out.push(n);
        return Ok(());
    }
    if let Some((root, k)) = n.perfect_power() {
        let mut fs = Vec::new();
        factor_node(root, threads, pc, tuning, after_split, progress, &mut fs)?;
        for _ in 0..k {
            out.extend(fs.iter().cloned())
        }
        return Ok(());
    }
    let mut split_by_rho = false;
    let d = match pollard_brent_natural(
        &n,
        rho_budget(n.bit_len(), tuning, after_split),
        NATIVE_RHO_CONSTANTS,
        || {
            progress(EngineProgress {
                phase: EnginePhase::Preprocessing,
                polynomials: 0,
                relations: 0,
                target: 0,
                workers: threads,
            })
        },
    )? {
        Some(factor) => {
            #[cfg(test)]
            stage_counts::bump(&stage_counts::RHO);
            split_by_rho = true;
            if tuning.profile {
                eprintln!(
                    "PROFILE rho input_bits={} factor_bits={} siqs=false",
                    n.bit_len(),
                    factor.bit_len()
                );
            }
            factor
        }
        None => find_factor(n.clone(), threads, tuning, progress)?,
    };
    if d.is_one() || d == n {
        return Err(EngineError::NoFactor);
    }
    // A split under rho is what proves this composite unbalanced, so it is what the children
    // inherit. A split under SIQS proves nothing about factor sizes and only passes along whatever
    // this node already knew.
    let children_after_split = after_split || split_by_rho;
    let q = n.div_rem(&d).unwrap().0;
    factor_node(d, threads, pc, tuning, children_after_split, progress, out)?;
    factor_node(q, threads, pc, tuning, children_after_split, progress, out)
}

#[cfg(any(unix, windows))]
fn find_factor(
    n: Natural,
    threads: usize,
    tuning: &FactorTuning,
    progress: &mut impl FnMut(EngineProgress) -> bool,
) -> Result<Natural, EngineError> {
    // Small inputs finish faster than 96 OS threads take to spawn and join, so
    // cap worker count by problem size to avoid parallel-startup overhead.
    let threads = match n.bit_len() {
        0..=128 => threads.min(2),
        129..=160 => threads.min(16),
        161..=184 => threads.min(48),
        _ => threads,
    }
    .max(1);
    if !progress(EngineProgress {
        phase: EnginePhase::BuildingFactorBase,
        polynomials: 0,
        relations: 0,
        target: 0,
        workers: threads,
    }) {
        return Err(EngineError::Cancelled);
    }
    #[cfg(test)]
    stage_counts::bump(&stage_counts::SIQS);
    let prof = tuning.profile;
    let t_fb = std::time::Instant::now();
    let ctx = prepare_with_la_threads(n.clone(), tuning, threads)?.0;
    let target = relation_target(ctx.base.len(), ctx.relation_percent);
    let budget = family_budget(n.bit_len());
    if prof {
        let (score_cutoff, score_bias, exact_scores) = siqs::choose_a(&ctx, 0)
            .map(|(a, _)| {
                let (bias, scan, exact) = score_plan(&ctx, &a);
                (scan.saturating_sub(bias), bias, exact)
            })
            .unwrap_or((0, 0, false));
        eprintln!(
            "PROFILE fb_build={:.3}s bits={} k={} q2={} nfb={} interval={} target_a_bits={} a_factors={} variants={} score_cutoff={} bias={} exact={} target={}",
            t_fb.elapsed().as_secs_f64(),
            ctx.n.bit_len(),
            ctx.sieve_n.div_rem(&ctx.n).unwrap().0.to_u64().unwrap_or(0),
            ctx.use_q2,
            ctx.base.len(),
            ctx.interval,
            ctx.target_a.bit_len(),
            ctx.a_factor_count,
            1usize
                << match ctx.n.bit_len() {
                    ..=320 => (ctx.a_factor_count - 1).min(10),
                    321..=333 => (ctx.a_factor_count - 1).min(11),
                    _ => (ctx.a_factor_count - 1).min(12),
                },
            score_cutoff,
            score_bias,
            exact_scores,
            target,
        );
    }
    let (job_tx, job_rx) = mpsc::channel::<Option<u64>>();
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (res_tx, res_rx) = mpsc::channel();
    let mut handles = Vec::new();
    let cancellation = Arc::new(AtomicBool::new(false));
    for _ in 0..threads {
        let rx = job_rx.clone();
        let tx = res_tx.clone();
        let c = ctx.clone();
        let cancellation = cancellation.clone();
        handles.push(std::thread::spawn(move || {
            let mut scratch = EngineScratch::default();
            loop {
                if cancellation.load(AtomicOrdering::Relaxed) {
                    break;
                }
                let job = rx.lock().unwrap_or_else(|e| e.into_inner()).recv();
                match job {
                    Ok(Some(f)) => {
                        if cancellation.load(AtomicOrdering::Relaxed) {
                            break;
                        }
                        if tx
                            .send(siqs::sieve_family_cancellable(
                                &c,
                                f,
                                &mut scratch,
                                &cancellation,
                            ))
                            .is_err()
                        {
                            break;
                        }
                    }
                    _ => break,
                }
            }
        }))
    }
    drop(res_tx);
    let mut next_send = 0u64;
    let mut completed = 0u64;
    let mut outstanding = 0usize;
    for _ in 0..threads * 2 {
        job_tx
            .send(Some(next_send))
            .map_err(|_| EngineError::Worker("worker job channel disconnected".into()))?;
        next_send += 1;
        outstanding += 1
    }
    let t_sieve = std::time::Instant::now();
    let mut collector = RelationCollector::new();
    let mut polynomials = 0u64;
    let mut total_survivors = 0u64;
    let mut relation_kinds = [0u64; 3];
    let mut timing = [0u64; 4];
    let mut seen_a = HashSet::new();
    let mut cancelled = false;
    while collector.columns.len() < target && completed < budget && !cancelled {
        let r = res_rx
            .recv()
            .map_err(|_| EngineError::Worker("worker result channel disconnected".into()))?;
        outstanding -= 1;
        completed += 1;
        // The native scheduler consumes completed families immediately.  The
        // portable EngineSession keeps deterministic family-order merging for
        // browser/distributed callers, but imposing that order here creates a
        // many-second head-of-line stall behind one slow 2,048-polynomial job.
        let unique_a = siqs::choose_a(&ctx, r.family)
            .map(|(a, _)| seen_a.insert(a))
            .unwrap_or(false);
        if unique_a {
            polynomials += r.polynomials;
            total_survivors += r.survivors;
            for (total, value) in timing.iter_mut().zip(r.timing) {
                *total += value;
            }
            for rel in r.relations {
                relation_kinds[match rel.large {
                    LargePrime::None => 0,
                    LargePrime::One(_) => 1,
                    LargePrime::Two(_, _) => 2,
                }] += 1;
                collector.ingest(rel, &n);
                if collector.columns.len() >= target {
                    break;
                }
            }
        }
        if prof && completed.is_multiple_of(50) {
            eprintln!(
                "PROFILE checkpoint elapsed={:.3}s polys={} families={} survivors={} accepted={}/{}/{} partials={} cycles={} relations={} cpu={:.1}/{:.1}/{:.1}/{:.1}",
                t_sieve.elapsed().as_secs_f64(),
                polynomials,
                completed,
                total_survivors,
                relation_kinds[0],
                relation_kinds[1],
                relation_kinds[2],
                collector.forest.relations.len(),
                collector.cycles,
                collector.columns.len(),
                timing[0] as f64 * 1e-9,
                timing[1] as f64 * 1e-9,
                timing[2] as f64 * 1e-9,
                timing[3] as f64 * 1e-9,
            );
        }
        if !progress(EngineProgress {
            phase: EnginePhase::Sieving,
            polynomials,
            relations: collector.columns.len(),
            target,
            workers: threads,
        }) {
            cancelled = true;
            cancellation.store(true, AtomicOrdering::Relaxed);
        }
        if collector.columns.len() >= target {
            // Other workers may be partway through a 2,048-polynomial family.
            // Stop them between polynomials instead of paying an otherwise
            // useless full-family tail before linear algebra can start.
            cancellation.store(true, AtomicOrdering::Relaxed);
        }
        while outstanding < threads * 2 && next_send < budget && collector.columns.len() < target {
            job_tx
                .send(Some(next_send))
                .map_err(|_| EngineError::Worker("worker job channel disconnected".into()))?;
            next_send += 1;
            outstanding += 1
        }
    }
    for _ in 0..threads {
        let _ = job_tx.send(None);
    }
    drop(job_tx);
    let mut first_panic = None;
    for h in handles {
        if let Err(payload) = h.join()
            && first_panic.is_none()
        {
            first_panic = Some(
                payload
                    .downcast_ref::<&str>()
                    .map(|s| (*s).to_owned())
                    .or_else(|| payload.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "worker panicked with a non-string payload".into()),
            );
        }
    }
    if let Some(message) = first_panic {
        return Err(EngineError::Worker(message));
    }
    if cancelled {
        return Err(EngineError::Cancelled);
    }
    // The loop above exits on one of three conditions, and only one of them means success. Without
    // this check, an exhausted family budget fell through into the linear algebra with a matrix
    // far too narrow to admit a dependency, and the resulting solver failure was reported as a
    // memory limit — an error that described neither the cause nor the remedy.
    if collector.columns.len() < target {
        if prof {
            eprintln!(
                "PROFILE budget_exhausted families={completed} relations={} target={target}",
                collector.columns.len()
            );
        }
        return Err(EngineError::InsufficientRelations);
    }
    if prof {
        eprintln!(
            "PROFILE sieve+collect={:.3}s polys={} families={} survivors={} accepted={}/{}/{} partials={} cycles={} relations={}",
            t_sieve.elapsed().as_secs_f64(),
            polynomials,
            completed,
            total_survivors,
            relation_kinds[0],
            relation_kinds[1],
            relation_kinds[2],
            collector.forest.relations.len(),
            collector.cycles,
            collector.columns.len()
        );
        eprintln!(
            "PROFILE worker_cpu setup={:.3}s score={:.3}s factor={:.3}s roots={:.3}s",
            timing[0] as f64 * 1e-9,
            timing[1] as f64 * 1e-9,
            timing[2] as f64 * 1e-9,
            timing[3] as f64 * 1e-9,
        );
    }
    if !progress(EngineProgress {
        phase: EnginePhase::LinearAlgebra,
        polynomials,
        relations: collector.columns.len(),
        target,
        workers: threads,
    }) {
        return Err(EngineError::Cancelled);
    }
    let t_la = std::time::Instant::now();
    let result = extract::extract(&ctx, &collector.columns);
    if prof {
        eprintln!(
            "PROFILE extract(LA)={:.3}s columns={}",
            t_la.elapsed().as_secs_f64(),
            collector.columns.len()
        );
    }
    if !progress(EngineProgress {
        phase: EnginePhase::Extracting,
        polynomials,
        relations: collector.columns.len(),
        target,
        workers: threads,
    }) {
        return Err(EngineError::Cancelled);
    }
    result
}

/// Expected log weight lost by skipping tiny primes.
fn small_slack(base: &[FactorBaseEntry], skip: u32, rounded: bool) -> usize {
    base.iter()
        .filter(|entry| entry.prime < skip)
        .map(|entry| {
            score_weight(entry.prime, rounded) as f64
                / (entry.prime.saturating_sub(1).max(1)) as f64
        })
        .sum::<f64>()
        .round() as usize
}

/// Sieve log weight of one factor-base prime, computed once while preparing
/// the immutable context.
#[inline(always)]
fn score_weight(prime: u32, rounded: bool) -> u8 {
    let floor = 31 - prime.leading_zeros();
    if !rounded {
        return (floor + 1) as u8;
    }
    // Round log2(p) to nearest instead of always rounding upward. The latter
    // accumulates roughly half a bit of false score per factor and becomes a
    // major survivor false-positive source in 90–100 digit tiers.
    let rounds_up = (prime as u64) * (prime as u64) >= (2u64 << (2 * floor));
    (floor + u32::from(rounds_up)) as u8
}

fn score_plan(ctx: &Context, a: &Natural) -> (u8, u8, bool) {
    let g_bits = ctx
        .sieve_n
        .bit_len()
        .saturating_sub(a.bit_len())
        .saturating_sub(usize::from(ctx.use_q2) * 2);
    let threshold = (g_bits as i32 - ctx.lp_bits as i32 - ctx.small_slack as i32
        + ctx.threshold_margin
        + ctx.thresh_adj)
        .clamp(1, u8::MAX as i32) as u8;
    let bias = 128u8.saturating_sub(threshold);
    let small_end = ctx
        .base
        .partition_point(|entry| entry.prime < ctx.small_skip);
    let smallest_scored = ctx
        .base
        .get(small_end)
        .map_or(u32::MAX, |entry| entry.prime)
        .max(3);
    let scored_bits = (32 - smallest_scored.leading_zeros() - 1).max(1);
    let score_bound = g_bits as u32 + g_bits as u32 / scored_bits + 1;
    let exact = bias as u32 + score_bound <= u8::MAX as u32;
    (bias, threshold.saturating_add(bias), exact)
}
/// Add the two root strides for one factor-base prime. Interleaving the roots and unrolling two
/// hits at a time mirrors FLINT's flat sieve kernel and cuts loop-control overhead in the dominant
/// score-write pass.
///
/// `SATURATING` is chosen per polynomial by the caller, which proves the no-overflow bound from the
/// smallest scored prime instead of assuming one (see `exact_scores` in [`sieve_one_poly`]).
/// Non-saturating writes are worth the check: forcing saturating addition everywhere measured
/// +1.5 cpu-s on 12.3 at 224-bit, about 7% of total sieve time.
#[inline(always)]
fn sieve_root_pair<const SATURATING: bool>(
    scores: &mut [u8],
    root1: usize,
    root2: usize,
    step: usize,
    weight: u8,
) {
    #[inline(always)]
    fn add<const SATURATING: bool>(slot: &mut u8, weight: u8) {
        *slot = if SATURATING {
            slot.saturating_add(weight)
        } else {
            slot.wrapping_add(weight)
        };
    }

    let len = scores.len();
    debug_assert!(root1 <= root2);
    let diff = root2 - root1;
    let mut pos = root1;
    while pos + step + diff < len {
        add::<SATURATING>(&mut scores[pos], weight);
        add::<SATURATING>(&mut scores[pos + diff], weight);
        pos += step;
        add::<SATURATING>(&mut scores[pos], weight);
        add::<SATURATING>(&mut scores[pos + diff], weight);
        pos += step;
    }
    while pos + diff < len {
        add::<SATURATING>(&mut scores[pos], weight);
        add::<SATURATING>(&mut scores[pos + diff], weight);
        pos += step;
    }
    while pos < len {
        add::<SATURATING>(&mut scores[pos], weight);
        pos += step;
    }
}

#[inline(always)]
fn add_score<const SATURATING: bool>(slot: &mut u8, weight: u8) {
    *slot = if SATURATING {
        slot.saturating_add(weight)
    } else {
        slot.wrapping_add(weight)
    };
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn sieve_linear_root<const SATURATING: bool>(
    ctx: &Context,
    b: &Natural,
    bneg: bool,
    c: &Natural,
    csign: bool,
    idx: usize,
    scores: &mut [u8],
) {
    let p = ctx.primes[idx];
    let mut bp = b.mod_u64(p as u64) as u32;
    if bneg && bp != 0 {
        bp = p - bp;
    }
    let denom = if ctx.use_q2 {
        bp
    } else {
        (2 * bp as u64 % p as u64) as u32
    };
    let Some(inv) = inv_u32(denom, p) else {
        return;
    };
    let cm = c.mod_u64(p as u64) as u32;
    let signed_c = if csign && cm != 0 { p - cm } else { cm };
    let xroot = mulmod_u32(if signed_c == 0 { 0 } else { p - signed_c }, inv, p);
    let mut pos = add_mod_u32(xroot, ctx.interval_mod_p[idx], p) as usize;
    while pos < scores.len() {
        add_score::<SATURATING>(&mut scores[pos], ctx.score_weights[idx]);
        pos += p as usize;
    }
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::needless_range_loop)]
#[allow(unsafe_code)]
fn score_polynomial<const SATURATING: bool>(
    ctx: &Context,
    b: &Natural,
    bneg: bool,
    c: &Natural,
    csign: bool,
    root1: &[u32],
    root2: &[u32],
    scores: &mut [u8],
    small_skip: u32,
    score_start: usize,
) {
    let score_len = (ctx.interval as usize) * 2;
    debug_assert!(scores.len() >= score_len + SCORE_SENTINELS);
    let score_begin = ctx.primes.partition_point(|&prime| prime < small_skip);
    let repeated_end = ctx
        .primes
        .partition_point(|&prime| (prime as usize) < score_len);

    {
        let active_scores = &mut scores[..score_len];
        // Dense blocking handled the ordinary roots in this prefix. The
        // handful of primes dividing A still have a polynomial-dependent
        // linear root.
        for idx in score_begin..score_start {
            if root1[idx] == u32::MAX {
                sieve_linear_root::<SATURATING>(ctx, b, bneg, c, csign, idx, active_scores);
            }
        }

        // Score the primes that hit repeatedly without a per-prime interval
        // size branch.
        let ordinary_start = score_start.max(score_begin);
        for idx in ordinary_start..repeated_end {
            let p = ctx.primes[idx];
            let weight = ctx.score_weights[idx];
            if root1[idx] != u32::MAX {
                sieve_root_pair::<SATURATING>(
                    active_scores,
                    root1[idx] as usize,
                    root2[idx] as usize,
                    p as usize,
                    weight,
                );
            } else {
                sieve_linear_root::<SATURATING>(ctx, b, bneg, c, csign, idx, active_scores);
            }
        }
    }

    // Sparse roots hit unpredictably. Redirect misses to a striped sentinel
    // tail with arithmetic masks, avoiding two data-dependent bounds branches
    // per prime. The reusable tail is never scanned for reports.
    let sparse_start = repeated_end.max(score_start).max(score_begin);
    let pointer = scores.as_mut_ptr();
    for idx in sparse_start..ctx.primes.len() {
        if root1[idx] == u32::MAX {
            sieve_linear_root::<SATURATING>(ctx, b, bneg, c, csign, idx, &mut scores[..score_len]);
            continue;
        }
        let r1 = root1[idx] as usize;
        let r2 = root2[idx] as usize;
        let weight = ctx.score_weights[idx];
        let sentinel = score_len + (idx & (SCORE_SENTINELS - 1));
        let mask1 = 0usize.wrapping_sub(usize::from(r1 < score_len));
        let target1 = (r1 & mask1) | (sentinel & !mask1);
        let mask2 = 0usize.wrapping_sub(usize::from(r2 < score_len));
        let target2 = (r2 & mask2) | (sentinel & !mask2);
        // SAFETY: an in-range root is below score_len; every redirected root
        // is within the SCORE_SENTINELS-element tail.
        unsafe {
            add_score::<SATURATING>(&mut *pointer.add(target1), weight);
        }
        if r2 != r1 {
            unsafe {
                add_score::<SATURATING>(&mut *pointer.add(target2), weight);
            }
        }
    }
}

const DENSE_BLOCK_LEN: usize = 32 * 1024;
const DENSE_PRIME_CUTOFF: u32 = 8 * 1024;
const SCORE_SENTINELS: usize = 1024;
#[inline(always)]
#[allow(unsafe_code)]
fn sieve_dense_root<const SATURATING: bool>(
    scores: &mut [u8],
    mut position: usize,
    step: usize,
    weight: u8,
) -> u32 {
    let pointer = scores.as_mut_ptr();
    while position + step < scores.len() {
        // SAFETY: the loop condition proves both positions are in the slice.
        unsafe {
            add_score::<SATURATING>(&mut *pointer.add(position), weight);
        }
        position += step;
        unsafe {
            add_score::<SATURATING>(&mut *pointer.add(position), weight);
        }
        position += step;
    }
    while position < scores.len() {
        unsafe {
            add_score::<SATURATING>(&mut *pointer.add(position), weight);
        }
        position += step;
    }
    (position - scores.len()) as u32
}

#[allow(clippy::too_many_arguments)]
#[allow(unsafe_code)]
fn score_dense_prefix<const SATURATING: bool>(
    ctx: &Context,
    root1: &[u32],
    root2: &[u32],
    scores: &mut [u8],
    small_skip: u32,
    dense_end: usize,
    dense_root1: &mut Vec<u32>,
    dense_root2: &mut Vec<u32>,
) {
    dense_root1.clear();
    dense_root1.extend_from_slice(&root1[..dense_end]);
    dense_root2.clear();
    dense_root2.extend_from_slice(&root2[..dense_end]);
    let dense_start = ctx.base.partition_point(|entry| entry.prime < small_skip);
    for block in scores.chunks_mut(DENSE_BLOCK_LEN) {
        for idx in dense_start..dense_end {
            // SAFETY: all slices cover the shared factor-base prefix.
            let first = unsafe { *dense_root1.get_unchecked(idx) };
            if first == u32::MAX {
                continue;
            }
            let prime = unsafe { ctx.base.get_unchecked(idx).prime as usize };
            let weight = unsafe { *ctx.score_weights.get_unchecked(idx) };
            let second = unsafe { *dense_root2.get_unchecked(idx) };
            let next1 = sieve_dense_root::<SATURATING>(block, first as usize, prime, weight);
            let next2 = sieve_dense_root::<SATURATING>(block, second as usize, prime, weight);
            unsafe {
                *dense_root1.get_unchecked_mut(idx) = next1;
                *dense_root2.get_unchecked_mut(idx) = next2;
            }
        }
    }
}

/// Divide `q` by `p` as often as it goes, returning the exponent.
///
/// The remainder is taken first even though that means two passes over the limbs on a successful
/// division: `rem_u64` does not materialize a quotient, while `div_rem_u64` builds and returns a
/// whole `Natural`, and the terminating call — the one that finds a nonzero remainder — happens on
/// every invocation. Folding the two into one `div_rem_u64` per iteration measured slower.
#[inline]
fn divide_out(q: &mut Natural, p: u64) -> u16 {
    let mut count = 0;
    while q.rem_u64(p) == 0 {
        *q = q.div_rem_u64(p).unwrap().0;
        count += 1;
    }
    count
}

/// Divide by a prime already proved to hit this polynomial position.
///
/// The ordinary helper first computes a remainder because most calls miss.
/// Resieving has established the first division here, so doing that remainder
/// pass again wastes a full limb traversal for the overwhelmingly common
/// exponent-one case.
#[inline]
fn divide_out_known(q: &mut Natural, p: u64) -> u16 {
    let (quotient, remainder) = q.div_rem_u64(p).unwrap();
    debug_assert_eq!(remainder, 0);
    *q = quotient;
    1 + divide_out(q, p)
}

/// Collect high-scoring positions, ascending.
///
/// Survivors are rare — a handful of positions in a 180–640 KiB score array — so this is a pure
/// rejection loop, and the caller biases every score by `128 − threshold` precisely so that
/// "reached the threshold" becomes "high bit set". One `and` then rejects eight positions at a time.
/// The bias is applied by the caller rather than folded in here because a runtime bias would cost an
/// `add` and an `or` per word, which measured 2.7× slower than the masked test on the 224-bit case.
#[cfg(target_arch = "x86_64")]
fn collect_candidates(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    report_simd::collect_candidates(scores, threshold, candidates);
}

#[cfg(all(target_arch = "wasm32", feature = "wasm-simd128"))]
fn collect_candidates(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    report_wasm::collect_candidates(scores, threshold, candidates);
}

#[cfg(not(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", feature = "wasm-simd128")
)))]
fn collect_candidates(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    collect_candidates_scalar(scores, threshold, candidates);
}

#[cfg(not(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", feature = "wasm-simd128")
)))]
fn collect_candidates_scalar(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    const HIGH: u64 = 0x8080_8080_8080_8080;
    debug_assert!(
        threshold >= 128,
        "scores must be biased so that the threshold is the byte's high bit"
    );
    candidates.clear();
    let mut chunks = scores.chunks_exact(8);
    for (word_index, chunk) in chunks.by_ref().enumerate() {
        let word = u64::from_ne_bytes(chunk.try_into().unwrap());
        if word & HIGH == 0 {
            continue;
        }
        let start = word_index * 8;
        // `threshold == 128` is the biased case, where the high bit alone decides and the byte
        // compare folds away. A deeper threshold leaves the word test as a prefilter only.
        if threshold == 128 {
            for (offset, &score) in chunk.iter().enumerate() {
                if score & 0x80 != 0 {
                    candidates.push((start + offset) as u32);
                }
            }
        } else {
            for (offset, &score) in chunk.iter().enumerate() {
                if score >= threshold {
                    candidates.push((start + offset) as u32);
                }
            }
        }
    }
    let start = scores.len() - chunks.remainder().len();
    for (offset, &score) in chunks.remainder().iter().enumerate() {
        if score >= threshold {
            candidates.push((start + offset) as u32);
        }
    }
}

/// Recover the scored factor-base hits for a batch of sieve reports.
///
/// Candidate-major trial division is excellent when a polynomial has one or
/// two reports. Once the threshold is deep enough for DLP, however, repeating
/// the full factor-base scan for every report dominates. Below the interval
/// length we scan the smaller prefix candidate-major. Above it, every root has
/// at most one hit, so one prime-major pass recovers all report factors.
#[allow(clippy::too_many_arguments)]
fn resieve_candidate_factors(
    ctx: &Context,
    root1: &[u32],
    root2: &[u32],
    scores: &[u8],
    _threshold: u8,
    candidates: &[u32],
    small_end: usize,
    counts: &mut Vec<u8>,
    factors: &mut Vec<u32>,
    candidate_bits: &mut Vec<u64>,
) -> bool {
    const FACTOR_STRIDE: usize = 48;
    #[inline]
    fn push_factor(
        counts: &mut [u8],
        factors: &mut [u32],
        candidate_index: usize,
        factor: u32,
    ) -> bool {
        let count = counts[candidate_index] as usize;
        if count == FACTOR_STRIDE {
            return false;
        }
        factors[candidate_index * FACTOR_STRIDE + count] = factor;
        counts[candidate_index] += 1;
        true
    }
    counts.clear();
    counts.resize(candidates.len(), 0);
    factors.resize(candidates.len() * FACTOR_STRIDE, 0);
    let bit_words = scores.len().div_ceil(64);
    if candidate_bits.len() < bit_words {
        candidate_bits.resize(bit_words, 0);
    }
    for &position in candidates {
        candidate_bits[position as usize >> 6] |= 1u64 << (position & 63);
    }
    let mut overflow = false;
    // Candidate-major costs one fast remainder per report and prime.
    // Prime-major costs roughly two interval/prime root visits with a direct
    // slot lookup. Move the crossover with the actual report density.
    let crossover_prime = ((2 * scores.len()) / candidates.len().max(1))
        .max(ctx.small_skip as usize)
        .min(u32::MAX as usize) as u32;
    let sparse_start = ctx
        .base
        .partition_point(|entry| entry.prime < crossover_prime);

    for (candidate_index, &position) in candidates.iter().enumerate() {
        for idx in small_end..sparse_start {
            let first = root1[idx];
            if first == u32::MAX {
                continue;
            }
            let prime = ctx.base[idx].prime;
            let residue = fastmod(position, prime, ctx.pinv[idx]);
            if residue == first || residue == root2[idx] {
                overflow |= !push_factor(counts, factors, candidate_index, idx as u32);
            }
        }
    }

    for idx in sparse_start..ctx.base.len() {
        let first = root1[idx];
        if first == u32::MAX {
            continue;
        }
        let prime = ctx.base[idx].prime as usize;
        for (root_number, root) in [first, root2[idx]].into_iter().enumerate() {
            if root_number == 1 && root == first {
                continue;
            }
            let mut position = root as usize;
            while position < scores.len() {
                if candidate_bits[position >> 6] & (1u64 << (position & 63)) != 0 {
                    let candidate_index = candidates
                        .binary_search(&(position as u32))
                        .expect("candidate bit filter must be exact");
                    let count = counts[candidate_index] as usize;
                    let offset = candidate_index * FACTOR_STRIDE;
                    if count == 0 || factors[offset + count - 1] != idx as u32 {
                        overflow |= !push_factor(counts, factors, candidate_index, idx as u32);
                    }
                }
                position += prime;
            }
        }
    }
    for &position in candidates {
        candidate_bits[position as usize >> 6] &= !(1u64 << (position & 63));
    }
    !overflow
}

#[allow(clippy::too_many_arguments)]
fn sieve_one_poly(
    ctx: &Context,
    a: &Natural,
    b: &Natural,
    bneg: bool,
    aidx: &[u32],
    root1: &[u32],
    root2: &[u32],
    scores: &mut Vec<u8>,
    candidates: &mut Vec<u32>,
    candidate_factor_counts: &mut Vec<u8>,
    candidate_factors: &mut Vec<u32>,
    candidate_bits: &mut Vec<u64>,
    dense_root1: &mut Vec<u32>,
    dense_root2: &mut Vec<u32>,
    powers_scratch: &mut Vec<(u32, u16)>,
    out: &mut Vec<Relation>,
    timing: &mut [u64; 4],
) -> usize {
    let base = &ctx.base;
    let len = (ctx.interval as usize) * 2;
    let small_skip = ctx.small_skip;
    let bb = b.checked_mul(b).unwrap();
    let (mut c, csign) = if bb >= ctx.sieve_n {
        let (quotient, remainder) = bb.wrapping_sub(&ctx.sieve_n).div_rem(a).unwrap();
        debug_assert!(remainder.is_zero());
        (quotient, false)
    } else {
        let (quotient, remainder) = ctx.sieve_n.wrapping_sub(&bb).div_rem(a).unwrap();
        debug_assert!(remainder.is_zero());
        (quotient, true)
    };
    if ctx.use_q2 {
        debug_assert!(c.trailing_zeros() >= 2);
        c >>= 2;
    }
    let (bias, scan_threshold, exact_scores) = score_plan(ctx, a);
    // Factor-base entries with `prime < small_skip` occupy the low indices (the base is sorted
    // ascending). Those tiny primes are not sieved — gating them would waste a `fastmod` where a
    // direct divide is cheaper (they divide most survivors) — so they are divided out directly.
    let small_end = base.partition_point(|e| e.prime < small_skip);
    scores.clear();
    scores.resize(len + SCORE_SENTINELS, bias);
    let score_started = ctx.profile.then(std::time::Instant::now);
    let dense_end = if ctx.n.bit_len() >= 289 {
        base.partition_point(|entry| entry.prime < DENSE_PRIME_CUTOFF)
    } else {
        0
    };
    if exact_scores {
        if dense_end != 0 {
            score_dense_prefix::<false>(
                ctx,
                root1,
                root2,
                &mut scores[..len],
                small_skip,
                dense_end,
                dense_root1,
                dense_root2,
            );
        }
        score_polynomial::<false>(
            ctx, b, bneg, &c, csign, root1, root2, scores, small_skip, dense_end,
        );
    } else {
        if dense_end != 0 {
            score_dense_prefix::<true>(
                ctx,
                root1,
                root2,
                &mut scores[..len],
                small_skip,
                dense_end,
                dense_root1,
                dense_root2,
            );
        }
        score_polynomial::<true>(
            ctx, b, bneg, &c, csign, root1, root2, scores, small_skip, dense_end,
        );
    }
    collect_candidates(&scores[..len], scan_threshold, candidates);
    if let Some(started) = score_started {
        timing[1] += started.elapsed().as_nanos().min(u64::MAX as u128) as u64;
    }
    if candidates.is_empty() {
        return 0;
    }
    let use_resieve = ctx.n.bit_len() >= 289
        && candidates.len() >= 4
        && resieve_candidate_factors(
            ctx,
            root1,
            root2,
            scores,
            scan_threshold,
            candidates,
            small_end,
            candidate_factor_counts,
            candidate_factors,
            candidate_bits,
        );
    let survivors = candidates.len();
    let factor_started = ctx.profile.then(std::time::Instant::now);

    let two_idx = base.iter().position(|e| e.prime == 2).map(|i| i as u32);
    for (candidate_index, &posu) in candidates.iter().enumerate() {
        let pos = posu as usize;
        // The score is the sum of one log weight per sieve hit. Once the primes divided out account
        // for that whole weight, no further scored factor-base prime can divide this candidate.
        // This is FLINT's `extra_bits < sieve[i]` stopping rule; it saves scanning the tail of the
        // factor base for the ~99% of survivors that are not smooth. It is exact only when scores
        // did not saturate, so saturating polynomials disable it by demanding an unreachable score.
        let score_target = if exact_scores {
            u16::from(scores[pos].wrapping_sub(bias))
        } else {
            u16::MAX
        };
        let mut confirmed_score = 0u16;
        let mut small_score = 0usize;
        let x = pos as i64 - ctx.interval as i64;
        let xabs = x.unsigned_abs();
        let ax = a.checked_mul(&Natural::from_u64(xabs)).unwrap();
        // Ordinary SIQS uses t=Ax+B. Q2 uses t=2Ax+B and
        // t²-kN=4A·g; the fixed square factor is recorded below.
        let tax = if ctx.use_q2 {
            ax.wrapping_add(&ax)
        } else {
            ax.clone()
        };
        let (t, tneg) = siqs::signed_add(&tax, x < 0, b, bneg);
        // Compute g directly with signs, avoiding the wide t² and division.
        let ax2 = ax.checked_mul(&Natural::from_u64(xabs)).unwrap();
        let bx_coefficient = if ctx.use_q2 {
            b.clone()
        } else {
            b.wrapping_add(b)
        };
        let bx = bx_coefficient
            .checked_mul(&Natural::from_u64(xabs))
            .unwrap();
        let (gx, gxneg) = siqs::signed_add(&ax2, false, &bx, bneg ^ (x < 0));
        let (mut q, sign) = siqs::signed_add(&gx, gxneg, &c, csign);
        if q.is_zero() {
            continue;
        }
        powers_scratch.clear();
        powers_scratch.extend(aidx.iter().copied().map(|index| (index, 1)));
        let record = |i: u32, count: u16, powers: &mut Vec<(u32, u16)>| {
            if count == 0 {
                return;
            }
            if let Some(value) = powers.iter_mut().find(|value| value.0 == i) {
                value.1 += count;
            } else {
                powers.push((i, count));
            }
        };
        if ctx.use_q2
            && let Some(ti) = two_idx
        {
            record(ti, 2, powers_scratch);
        }
        // Prime 2 (not sieved): strip via trailing zeros.
        if let Some(ti) = two_idx {
            let c2 = q.trailing_zeros();
            if c2 != 0 {
                q >>= c2;
                small_score += c2;
                record(ti, c2 as u16, powers_scratch);
            }
        }
        // Small primes are not score-sieved, but still use the same cheap position-root gate as the
        // main factor base. This follows FLINT: it avoids both their disproportionately dense score
        // writes and an unconditional big-integer remainder for every survivor.
        for (i, e) in base[..small_end].iter().enumerate() {
            let p = e.prime as u64;
            if p == 2 {
                continue;
            }
            let r1 = root1[i];
            if r1 != u32::MAX {
                let posmodp = fastmod(posu, e.prime, ctx.pinv[i]);
                if posmodp != r1 && posmodp != root2[i] {
                    continue;
                }
            }
            let count = if r1 == u32::MAX {
                divide_out(&mut q, p)
            } else {
                divide_out_known(&mut q, p)
            };
            small_score += usize::from(count) * usize::from(ctx.score_weights[i]);
            record(i as u32, count, powers_scratch);
        }
        // The score threshold gives back the *average* contribution of the
        // unsieved tiny primes. Replace that estimate with the contribution
        // actually observed before paying for the large-factor divisions.
        // Extra scored bits compensate a below-average tiny part.
        let scored = usize::from(scores[pos].wrapping_sub(bias));
        let required_scored = usize::from(scan_threshold.wrapping_sub(bias));
        let surplus = scored.saturating_sub(required_scored);
        if ctx.n.bit_len() >= 289 && small_score + surplus < ctx.small_slack {
            continue;
        }
        // Primes dividing `a` (seeded at exponent 1, root1 == MAX so not gated) — divide directly.
        for &ai in aidx {
            let p = base[ai as usize].prime as u64;
            let count = divide_out(&mut q, p);
            if count != 0 && p >= small_skip as u64 {
                confirmed_score += u16::from(ctx.score_weights[ai as usize]);
            }
            record(ai, count, powers_scratch);
        }
        if use_resieve {
            const FACTOR_STRIDE: usize = 48;
            let offset = candidate_index * FACTOR_STRIDE;
            let count = candidate_factor_counts[candidate_index] as usize;
            for &index in &candidate_factors[offset..offset + count] {
                let p = base[index as usize].prime;
                let count = divide_out_known(&mut q, p as u64);
                record(index, count, powers_scratch);
            }
        } else {
            for idx in small_end..base.len() {
                if q.is_one() || confirmed_score >= score_target {
                    break;
                }
                let r1 = root1[idx];
                if r1 == u32::MAX {
                    continue;
                }
                let p = base[idx].prime;
                // For the sparse factor-base tail `p` exceeds the score-array
                // position, so the residue is the position itself. This removes
                // two wide multiplies from most survivor/prime gates at 90–100
                // decimal digits.
                let position_mod_p = if posu < p {
                    posu
                } else {
                    fastmod(posu, p, ctx.pinv[idx])
                };
                if position_mod_p != r1 && position_mod_p != root2[idx] {
                    continue;
                }
                let count = divide_out_known(&mut q, p as u64);
                if count != 0 {
                    confirmed_score += u16::from(ctx.score_weights[idx]);
                }
                record(idx as u32, count, powers_scratch);
            }
        }
        let large = if q.is_one() {
            LargePrime::None
        } else if q.bit_len() > 64 {
            continue;
        } else {
            match classify_cofactor(q.as_parts()[0], ctx.single_limit, ctx.double_limit) {
                Some(lp) => lp,
                None => continue,
            }
        };
        let mut root = t.div_rem(&ctx.n).unwrap().1;
        if tneg && !root.is_zero() {
            root = ctx.n.wrapping_sub(&root)
        }
        out.push(Relation {
            root,
            sign,
            powers: powers_scratch.clone(),
            large,
        });
    }
    if let Some(started) = factor_started {
        timing[2] += started.elapsed().as_nanos().min(u64::MAX as u128) as u64;
    }
    survivors
}

fn to_column(r: Relation) -> Column {
    Column {
        root: r.root,
        sign: r.sign,
        powers: r.powers.into_iter().map(|(i, e)| (i, e as u32)).collect(),
        extra_sqrt: Vec::new(),
    }
}

/// Combine a set of relations whose large primes all cancel (each appears an even
/// number of times) into a single full-relation column. The cancelled large primes
/// contribute (count/2) copies to the reconstructed square root.
fn combine_cycle<'a>(rels: impl IntoIterator<Item = &'a Relation>, n: &Natural) -> Column {
    let mut root = Natural::ONE;
    let mut sign = false;
    let mut powers: BTreeMap<u32, u32> = BTreeMap::new();
    let mut lp: BTreeMap<u64, u32> = BTreeMap::new();
    for r in rels {
        root = root.mul_mod(&r.root, n);
        sign ^= r.sign;
        for &(i, e) in &r.powers {
            *powers.entry(i).or_default() += e as u32;
        }
        let (ps, k) = r.large.primes();
        for &p in &ps[..k] {
            *lp.entry(p).or_default() += 1;
        }
    }
    let mut extra_sqrt = Vec::new();
    for (p, c) in lp {
        for _ in 0..c / 2 {
            extra_sqrt.push(p);
        }
    }
    Column {
        root,
        sign,
        powers: powers.into_iter().collect(),
        extra_sqrt,
    }
}

/// Classify a factored-out cofactor (>1, fits in `u64`) as a single or double
/// large prime, or reject it. Portable (no threads / native-only deps).
fn classify_cofactor(q: u64, single_limit: u64, double_limit: u64) -> Option<LargePrime> {
    // An exact residual may be used as a graph vertex without proving it
    // prime: when the same value closes a cycle, the entire residual is
    // squared out. In normal SIQS use it is prime anyway, because all possible
    // factor-base divisors have just been removed.
    if q <= single_limit {
        return Some(LargePrime::One(q));
    }
    if double_limit == 0 || q > double_limit {
        return None;
    }
    // Almost every residual in this range is prime. One base-2 Fermat test is
    // sufficient as a cheap rejection filter; base-2 pseudoprimes continue to
    // exact splitting and factor primality checks below.
    if crate::u64math::pow_mod(2, q - 1, q) == 1 {
        return None;
    }
    let d = crate::u64math::squfof(q).or_else(|| pollard_u64(q))?;
    let e = q / d;
    if d > 1
        && e > 1
        && d <= single_limit
        && e <= single_limit
        && crate::u64math::is_prime(d)
        && crate::u64math::is_prime(e)
    {
        Some(LargePrime::Two(d.min(e), d.max(e)))
    } else {
        None
    }
}

/// Bounded, cancellable Pollard-Brent over `Natural`. `iteration_limit` is the total across every
/// polynomial constant in `constants`, not per constant, so the caller's budget is the whole cost of
/// the stage; see [`rho_budget`].
///
/// `constants` selects which walks to run. The native ladder passes [`NATIVE_RHO_CONSTANTS`] and
/// runs them in sequence; the browser hands a disjoint range to each of its rho workers so that the
/// pool runs that many independent walks at once and the first collision wins.
#[cfg(any(unix, windows, target_arch = "wasm32"))]
pub(crate) fn pollard_brent_natural(
    n: &Natural,
    iteration_limit: u64,
    constants: core::ops::RangeInclusive<u64>,
    mut keep_going: impl FnMut() -> bool,
) -> Result<Option<Natural>, EngineError> {
    if n.is_even() {
        return Ok(Some(Natural::from_u64(2)));
    }
    // One conversion at each boundary amortizes across the entire rho stage.
    // All polynomial values and batched products remain Montgomery residues;
    // gcd(qR mod n, n) == gcd(q, n) because odd n makes R invertible.
    let montgomery = MontgomeryContext::new(n).expect("rho modulus is odd and engine-sized");
    let mut iterations = 0u64;
    for c_value in constants {
        if iterations >= iteration_limit {
            return Ok(None);
        }
        if !keep_going() {
            return Err(EngineError::Cancelled);
        }
        // The whole walk runs in raw limb buffers rather than `Natural`s. Every one of these starts
        // zeroed and nothing below writes above the modulus width, so the operations never touch —
        // or have to clear — the capacity the modulus does not use. At a 128-bit modulus that is 14
        // of 16 limbs.
        let mut c = [0 as crate::natural::Limb; crate::natural::LIMB_CAP];
        let mut y = [0 as crate::natural::Limb; crate::natural::LIMB_CAP];
        let mut x = [0 as crate::natural::Limb; crate::natural::LIMB_CAP];
        let mut ys = [0 as crate::natural::Limb; crate::natural::LIMB_CAP];
        let mut q = [0 as crate::natural::Limb; crate::natural::LIMB_CAP];
        let mut difference = [0 as crate::natural::Limb; crate::natural::LIMB_CAP];
        montgomery.load(&montgomery.encode(&Natural::from_u64(c_value)), &mut c);
        montgomery.load(&montgomery.encode(&Natural::from_u64(2)), &mut y);
        let mut r = 1u64;
        let mut g = Natural::ONE;
        while g.is_one() && iterations < iteration_limit {
            x.copy_from_slice(&y);
            for _ in 0..r {
                montgomery.sqr_add_assign(&mut y, &c);
                iterations += 1;
                if iterations >= iteration_limit {
                    break;
                }
            }
            let mut k = 0u64;
            while k < r && g.is_one() && iterations < iteration_limit {
                if !keep_going() {
                    return Err(EngineError::Cancelled);
                }
                ys.copy_from_slice(&y);
                montgomery.load(&montgomery.one(), &mut q);
                // Iterations per gcd. The gcd is the one part of this loop that does not get
                // cheaper with the modulus, so the right batch depends on what a gcd costs against
                // what an iteration costs: at 128 it was 16-45% of the stage before `Natural::gcd`
                // was made width-aware, and 5-9% after. Above 512 the curve is flat and the two
                // costs a larger batch carries — the iterations run past a collision before the
                // next gcd sees it, and the wider window for `q` to collect every factor of `n` and
                // force backtracking — stop being paid for. Neither risk reaches short searches at
                // all: `batch` is bounded by `r`, so it only binds once the walk is already deep.
                let batch = (r - k).min(512);
                for _ in 0..batch {
                    montgomery.sqr_add_assign(&mut y, &c);
                    montgomery.sub_raw(&x, &y, &mut difference);
                    if !montgomery.is_zero_raw(&difference) {
                        montgomery.mul_assign(&mut q, &difference);
                    }
                    iterations += 1;
                    if iterations >= iteration_limit {
                        break;
                    }
                }
                g = montgomery.store(&q).gcd(n);
                k += batch;
            }
            r = r.saturating_mul(2);
        }
        if g == *n {
            loop {
                if !keep_going() {
                    return Err(EngineError::Cancelled);
                }
                montgomery.sqr_add_assign(&mut ys, &c);
                montgomery.sub_raw(&x, &ys, &mut difference);
                g = montgomery.store(&difference).gcd(n);
                iterations += 1;
                if !g.is_one() || iterations >= iteration_limit {
                    break;
                }
            }
        }
        if !g.is_one() && g != *n {
            return Ok(Some(g));
        }
    }
    Ok(None)
}

struct Mont64 {
    modulus: u64,
    inverse: u64,
}
impl Mont64 {
    fn new(modulus: u64) -> Self {
        debug_assert!(modulus & 1 == 1);
        let mut inverse = 1u64;
        for _ in 0..6 {
            inverse = inverse.wrapping_mul(2u64.wrapping_sub(modulus.wrapping_mul(inverse)));
        }
        Self {
            modulus,
            inverse: inverse.wrapping_neg(),
        }
    }
    #[inline]
    fn reduce(&self, value: u128) -> u64 {
        let multiplier = (value as u64).wrapping_mul(self.inverse);
        let (sum, carry) = value.overflowing_add(multiplier as u128 * self.modulus as u128);
        let mut reduced = (sum >> 64) + ((carry as u128) << 64);
        if reduced >= self.modulus as u128 {
            reduced -= self.modulus as u128;
        }
        reduced as u64
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.reduce(a as u128 * b as u128)
    }
    fn encode(&self, value: u64) -> u64 {
        (((value % self.modulus) as u128) << 64)
            .checked_rem(self.modulus as u128)
            .unwrap() as u64
    }
    fn decode(&self, value: u64) -> u64 {
        self.reduce(value as u128)
    }
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        let (sum, overflow) = a.overflowing_add(b);
        if overflow || sum >= self.modulus {
            sum.wrapping_sub(self.modulus)
        } else {
            sum
        }
    }
}

/// Pollard-Brent with Montgomery multiplication and batched GCD.
fn pollard_u64(n: u64) -> Option<u64> {
    if n.is_multiple_of(2) {
        return Some(2);
    }
    let gcd = |mut a: u64, mut b: u64| {
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        a
    };
    let mont = Mont64::new(n);
    for c in 1..64 {
        let c_mont = mont.encode(c);
        let mut y = mont.encode(2);
        let mut r = 1u64;
        let mut g = 1u64;
        let mut x = 0u64;
        let mut ys = 0u64;
        while g == 1 {
            x = y;
            for _ in 0..r {
                y = mont.add(mont.mul(y, y), c_mont);
            }
            let mut k = 0;
            while k < r && g == 1 {
                ys = y;
                let batch = (r - k).min(128);
                let mut product = mont.encode(1);
                for _ in 0..batch {
                    y = mont.add(mont.mul(y, y), c_mont);
                    let difference = x.abs_diff(y);
                    if difference != 0 {
                        product = mont.mul(product, difference);
                    }
                }
                // Multiplication by the Montgomery radix is invertible modulo
                // odd `n`, so gcd(product·R mod n, n) is the desired gcd.
                g = gcd(product, n);
                k += batch;
            }
            r = r.saturating_mul(2);
        }
        if g == n {
            loop {
                ys = mont.add(mont.mul(ys, ys), c_mont);
                g = gcd(x.abs_diff(ys), n);
                if g != 1 {
                    break;
                }
            }
        }
        if g != n {
            return Some(g);
        }
    }
    None
}

/// A spanning forest over large-prime vertices. Each relation is an edge between
/// its large primes (single-large-prime relations use the reserved unit vertex
/// `1`). A relation that closes a cycle combines every relation on the cycle into
/// a full-relation column, since all large primes on a cycle cancel.
#[derive(Default)]
struct Forest {
    id_of: PrimeMap,
    parent: Vec<u32>,
    edge: Vec<Option<u32>>,
    size: Vec<u32>,
    relations: Vec<Relation>,
    reroot_vertices: Vec<u32>,
    reroot_edges: Vec<u32>,
}
impl Forest {
    fn vertex(&mut self, prime: u64) -> u32 {
        if let Some(&id) = self.id_of.get(&prime) {
            return id;
        }
        let id = self.parent.len() as u32;
        self.id_of.insert(prime, id);
        self.parent.push(id);
        self.edge.push(None);
        self.size.push(1);
        id
    }
    fn root(&self, mut v: u32) -> u32 {
        while self.parent[v as usize] != v {
            v = self.parent[v as usize];
        }
        v
    }
    fn path(&self, mut v: u32, out: &mut Vec<u32>) {
        while self.parent[v as usize] != v {
            out.push(self.edge[v as usize].unwrap());
            v = self.parent[v as usize];
        }
    }
    /// Re-root the tree containing `v` so that `v` becomes its root.
    fn reroot(&mut self, v: u32) {
        self.reroot_vertices.clear();
        self.reroot_edges.clear();
        self.reroot_vertices.push(v);
        let mut c = v;
        while self.parent[c as usize] != c {
            self.reroot_edges.push(self.edge[c as usize].unwrap());
            c = self.parent[c as usize];
            self.reroot_vertices.push(c);
        }
        self.parent[v as usize] = v;
        self.edge[v as usize] = None;
        for i in 0..self.reroot_edges.len() {
            self.parent[self.reroot_vertices[i + 1] as usize] = self.reroot_vertices[i];
            self.edge[self.reroot_vertices[i + 1] as usize] = Some(self.reroot_edges[i]);
        }
    }
    fn link(&mut self, a: u32, b: u32, rel: Relation) {
        let relation_index = self.relations.len() as u32;
        self.relations.push(rel);
        let root_a = self.root(a);
        let root_b = self.root(b);
        let size_a = self.size[root_a as usize];
        let size_b = self.size[root_b as usize];
        if size_a >= size_b {
            self.reroot(b);
            self.parent[b as usize] = a;
            self.edge[b as usize] = Some(relation_index);
            self.size[root_a as usize] = size_a.saturating_add(size_b);
            self.size[b as usize] = 0;
        } else {
            self.reroot(a);
            self.parent[a as usize] = b;
            self.edge[a as usize] = Some(relation_index);
            self.size[root_b as usize] = size_a.saturating_add(size_b);
            self.size[a as usize] = 0;
        }
    }
}

/// Deterministically accumulates relations into matrix columns, matching partial
/// relations through the large-prime graph.
struct RelationCollector {
    forest: Forest,
    columns: Vec<Column>,
    /// Partials combined into full relations by closing a large-prime cycle. Reported by the
    /// profiling path so the retained-partial count and the cycle yield can be compared directly.
    cycles: usize,
    cycle_path: Vec<u32>,
}
impl RelationCollector {
    fn new() -> Self {
        Self {
            forest: Forest::default(),
            columns: Vec::new(),
            cycles: 0,
            cycle_path: Vec::new(),
        }
    }
    fn ingest(&mut self, rel: Relation, n: &Natural) {
        match rel.large {
            LargePrime::None => self.columns.push(to_column(rel)),
            LargePrime::One(p) => self.edge(p, 1, rel, n),
            LargePrime::Two(a, b) if a == b => {
                self.cycles += 1;
                self.columns.push(combine_cycle([&rel], n))
            }
            LargePrime::Two(a, b) => self.edge(a, b, rel, n),
        }
    }
    fn edge(&mut self, pa: u64, pb: u64, rel: Relation, n: &Natural) {
        let va = self.forest.vertex(pa);
        let vb = self.forest.vertex(pb);
        if self.forest.root(va) == self.forest.root(vb) {
            self.cycle_path.clear();
            self.forest.path(va, &mut self.cycle_path);
            self.forest.path(vb, &mut self.cycle_path);
            self.cycles += 1;
            self.columns.push(combine_cycle(
                core::iter::once(&rel).chain(
                    self.cycle_path
                        .iter()
                        .map(|&index| &self.forest.relations[index as usize]),
                ),
                n,
            ));
        } else {
            self.forest.link(va, vb, rel);
        }
    }
}
fn inv_u32(a: u32, p: u32) -> Option<u32> {
    if a == 0 {
        return None;
    }
    if p == 2 {
        return Some(1);
    }
    let (mut u, mut v) = (a, p);
    let (mut x1, mut x2) = (1u64, 0u64);
    let modulus = p as u64;
    while u != 1 && v != 1 {
        while u & 1 == 0 {
            u >>= 1;
            x1 = if x1 & 1 == 0 {
                x1 >> 1
            } else {
                (x1 + modulus) >> 1
            };
        }
        while v & 1 == 0 {
            v >>= 1;
            x2 = if x2 & 1 == 0 {
                x2 >> 1
            } else {
                (x2 + modulus) >> 1
            };
        }
        if u >= v {
            u -= v;
            x1 = if x1 >= x2 { x1 - x2 } else { x1 + modulus - x2 };
        } else {
            v -= u;
            x2 = if x2 >= x1 { x2 - x1 } else { x2 + modulus - x1 };
        }
    }
    Some(if u == 1 { x1 } else { x2 } as u32)
}
fn mulmod_u32(a: u32, b: u32, p: u32) -> u32 {
    (a as u64 * b as u64 % p as u64) as u32
}
#[inline]
fn add_mod_u32(a: u32, b: u32, p: u32) -> u32 {
    debug_assert!(a < p && b < p && p <= u32::MAX / 2);
    let sum = a + b;
    if sum >= p { sum - p } else { sum }
}
#[inline]
fn sub_mod_u32(a: u32, b: u32, p: u32) -> u32 {
    if a >= b { a - b } else { a + (p - b) }
}
/// Lemire fast-mod constant for divisor `p`: `⌊2^64 / p⌋ + 1` (via `u64::MAX / p + 1`, which equals
/// it for every `p ≥ 2`). Precomputed once per factor-base prime; see [`fastmod`].
#[inline]
fn lemire_c(p: u32) -> u64 {
    (u64::MAX / p as u64) + 1
}
/// `a mod p` by Daniel Lemire's "faster remainder" (multiply-shift, no hardware divide). Exact for
/// `a, p < 2^32` with `c == lemire_c(p)`. Used to gate trial division: a factor-base prime `p`
/// divides `g(x)` at position `x` iff `x mod p` equals one of its two sieve roots.
#[inline]
fn fastmod(a: u32, p: u32, c: u64) -> u32 {
    let lowbits = c.wrapping_mul(a as u64);
    ((lowbits as u128 * p as u128) >> 64) as u32
}
/// The established browser-tier multiplier policy, retained as a performance
/// compatibility boundary through 288 bits.
fn knuth_schroeppel_legacy(n: &Natural) -> u64 {
    const MULTIPLIERS: [u64; 29] = [
        1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37,
        38, 41, 42, 43, 47,
    ];
    let nmod8 = n.mod_u64(8);
    let mut weights = [0.0f64; MULTIPLIERS.len()];
    for (weight, &k) in weights.iter_mut().zip(&MULTIPLIERS) {
        let mod8 = nmod8 * k % 8;
        let mut value = 0.346_573_59_f64;
        if mod8 == 1 {
            value *= 4.0;
        } else if mod8 == 5 {
            value *= 2.0;
        }
        *weight = value - (k as f64).ln() / 2.0;
    }
    let mut p = 3u64;
    let mut seen = 0usize;
    while seen < 500 {
        if crate::u64math::is_prime(p) {
            seen += 1;
            let nmod = n.mod_u64(p);
            if nmod != 0 {
                let contribution = (p as f64).ln() / p as f64;
                let symbol = jacobi_u64(nmod, p) as i32;
                for (weight, &k) in weights.iter_mut().zip(&MULTIPLIERS) {
                    let km = k % p;
                    if km == 0 {
                        *weight += contribution;
                    } else if symbol * jacobi_u64(km, p) as i32 == 1 {
                        *weight += 2.0 * contribution;
                    }
                }
            }
        }
        p += 2;
    }
    weights
        .iter()
        .zip(MULTIPLIERS)
        .max_by(|(left, _), (right, _)| left.total_cmp(right))
        .map_or(1, |(_, multiplier)| multiplier)
}

/// Modified Knuth-Schroeppel multiplier selection.
///
/// The complete odd-heavy table matters: RSA-100 selects 139, which the old
/// abbreviated table (ending at 47) could not even consider. Scores are the
/// expected logarithmic size left after sieving; lower is better.
fn knuth_schroeppel(n: &Natural) -> u64 {
    const MULTIPLIERS: [u64; 114] = [
        1, 2, 3, 5, 7, 9, 10, 11, 13, 14, 15, 17, 19, 21, 23, 25, 29, 31, 33, 35, 37, 39, 41, 43,
        45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 83, 85, 87, 89, 91,
        93, 95, 97, 101, 103, 105, 107, 109, 111, 113, 115, 119, 121, 123, 127, 129, 131, 133, 137,
        139, 141, 143, 145, 147, 149, 151, 155, 157, 159, 161, 163, 165, 167, 173, 177, 179, 181,
        183, 185, 187, 191, 193, 195, 197, 199, 201, 203, 205, 209, 211, 213, 215, 217, 219, 223,
        227, 229, 231, 233, 235, 237, 239, 241, 249, 251, 253, 255,
    ];
    const KS_PRIMES: usize = 300;
    let nmod8 = n.mod_u64(8);
    let mut scores = [0.0f64; MULTIPLIERS.len()];
    for (score, &k) in scores.iter_mut().zip(&MULTIPLIERS) {
        let mod8 = (nmod8 * k) % 8;
        *score = (k as f64).ln() / 2.0;
        *score -= match mod8 {
            // kN == 1 (mod 8) permits the Q/2 polynomial.
            1 => 2.625 * core::f64::consts::LN_2,
            5 => core::f64::consts::LN_2,
            3 | 7 => 0.5 * core::f64::consts::LN_2,
            _ => 0.0,
        };
    }
    // Weight small primes for which kN is a quadratic residue. A regular
    // factor-base prime has two roots; a prime dividing kN has one.
    let mut p = 3u64;
    let mut seen = 0usize;
    while seen < KS_PRIMES {
        if crate::u64math::is_prime(p) {
            seen += 1;
            let nmod = n.mod_u64(p);
            let contribution = (p as f64).ln() / (p - 1) as f64;
            for (score, &k) in scores.iter_mut().zip(&MULTIPLIERS) {
                let knmod = nmod * (k % p) % p;
                if knmod == 0 {
                    *score -= contribution;
                } else if jacobi_u64(knmod, p) == 1 {
                    *score -= 2.0 * contribution;
                }
            }
        }
        p += 2;
    }
    scores
        .iter()
        .zip(MULTIPLIERS)
        .min_by(|(left, _), (right, _)| left.total_cmp(right))
        .map_or(1, |(_, multiplier)| multiplier)
}
#[cfg(test)]
mod tests {
    use super::*;

    /// The sieve range is enforced at `prepare`, the one point every scheduler passes through, and
    /// it is enforced on the composite rather than on anything the caller started with.
    #[test]
    fn prepare_rejects_composites_above_the_siqs_ceiling() {
        // 401 bits: one over the line, and a genuine semiprime so nothing else can reject it.
        let wide = Natural::from_decimal(
            "3872498856681097288856216400786342693014438579089271057383667301012033\
             051304191351486034920193821931899864729911531368699",
        )
        .expect("401-bit test value");
        assert_eq!(wide.bit_len(), 401);
        assert!(matches!(
            prepare(wide, &FactorTuning::default()),
            Err(EngineError::SiqsInputTooLarge(401))
        ));

        // Directly below the ceiling the same call must build a real sieve, so the bound is a
        // boundary rather than a blanket refusal of large work.
        let inside = Natural::from_decimal(
            "2088215395251987988508198170866285429326513389712811066669121976156893\
             154145354354068456327350848764053227833851662160699",
        )
        .expect("400-bit test value");
        assert_eq!(inside.bit_len(), 400);
        assert!(prepare(inside, &FactorTuning::default()).is_ok());
    }

    /// Extraction below the relation target is not a solver failure, and must not be reported as
    /// one. This is the session-side half of the fix; the native scheduler stops before its own
    /// `extract` call for the same reason.
    #[test]
    fn session_extraction_below_target_reports_insufficient_relations() {
        let p = Natural::from_u64(21_293_688_545_713_669);
        let q = Natural::from_u64(31_385_813_854_515_511);
        let context = prepare(p.checked_mul(&q).unwrap(), &FactorTuning::default()).unwrap();
        let session = EngineSession::new(context);
        assert!(session.target() > 0);
        assert_eq!(session.relations(), 0);
        assert!(!session.is_ready());
        // Not yet exhausted: no family has been issued, so there is budget left to spend.
        assert!(!session.budget_exhausted());
        assert!(session.family_budget() >= 100_000);
        assert!(matches!(
            session.extract_factor(),
            Err(EngineError::InsufficientRelations)
        ));
    }

    /// Pollard-Brent is deliberately outside the sieve's range limit: it is how a wide input made
    /// of small factors gets peeled at all. Sizing its budget must therefore stay defined — and
    /// nonzero — for widths the sieve itself refuses.
    #[test]
    fn rho_budget_is_defined_above_the_siqs_ceiling() {
        let tuning = FactorTuning::default();
        for bits in [MAX_SIQS_BITS + 1, 512, 1024] {
            assert!(
                rho_budget(bits, &tuning, false) >= 1_024,
                "{bits}-bit rho budget"
            );
        }
    }

    /// Above the ceiling the budget is the whole factoring attempt rather than a fraction of a
    /// sieve run, so it has to be enough to reach the factors the stage is being asked to find.
    /// `1.2·sqrt(p)` is Brent's expected cost; the assertions below are in those terms, with margin
    /// for the spread around the expectation.
    #[test]
    fn wide_rho_budget_reaches_the_factor_sizes_it_claims() {
        let tuning = FactorTuning::default();
        let brent_cost = |factor_bits: u32| (1.2 * 2f64.powi(factor_bits as i32).sqrt()) as u64;
        for bits in [MAX_SIQS_BITS + 1, 512, 513, 768, 769, 1024] {
            let budget = rho_budget(bits, &tuning, false);
            // The contract this release commits to: any factor up to 32 bits, with the whole
            // budget being hundreds of times the expected cost of finding one.
            assert!(
                budget >= 100 * brent_cost(32),
                "{bits}-bit budget {budget} is not a comfortable 32-bit reach"
            );
            // And the reach the wall-clock target was chosen to buy.
            assert!(
                budget >= 2 * brent_cost(48),
                "{bits}-bit budget {budget} does not reach a 48-bit factor"
            );
            // Sanity in the other direction: this arm is not a licence to run for hours.
            assert!(
                budget <= brent_cost(56),
                "{bits}-bit budget {budget} is past the documented one-minute target"
            );
        }
        // Nothing at or below the ceiling changed: SIQS runs there, so a deep rho would be pure
        // overhead on every balanced input.
        assert!(
            rho_budget(MAX_SIQS_BITS, &tuning, false) < wide_rho_budget(MAX_SIQS_BITS + 1) / 10
        );
        // Wider inputs cost more per iteration, so the iteration count comes down to hold the
        // wall clock roughly flat.
        assert!(wide_rho_budget(512) > wide_rho_budget(1024));
    }

    /// Balanced semiprimes at or below the sieve ceiling are what this engine is for, and rho
    /// contributes nothing to them: it runs to its budget, finds nothing, and hands the input to
    /// SIQS, so every iteration it spends there is overhead on the main workload. Raising the budget
    /// above the ceiling must therefore leave this range byte-for-byte alone. These are the 0.4.2
    /// values, pinned so that a future change to the wide arm cannot leak downward unnoticed.
    #[test]
    fn budgets_at_and_below_the_ceiling_are_unchanged() {
        let tuning = FactorTuning::default();
        for (bits, expected) in [
            (64, 1_024),
            (128, 1_024),
            (192, 18_022),
            (224, 45_875),
            (256, 176_947),
            (288, 734_003),
            (320, 2_703_360),
            (364, 1_572_864),
            (MAX_SIQS_BITS, 6_291_456),
        ] {
            assert_eq!(
                rho_budget(bits, &tuning, false),
                expected,
                "{bits}-bit rho budget"
            );
        }
    }

    /// A cofactor that split under rho is not a balanced semiprime, which is the only shape the
    /// sieve-fraction budget is sized for. From `DEEP_RHO_MIN_BITS` up it therefore gets the deep
    /// budget — the difference between peeling a wide multi-factor composite and handing its
    /// remainder to a sieve that would take weeks. Below that line the sieve is genuinely the better
    /// tool, so nothing changes there.
    #[test]
    fn a_split_under_rho_deepens_the_budget_where_the_sieve_is_expensive() {
        let tuning = FactorTuning::default();
        for bits in [289, 320, 364, 399, MAX_SIQS_BITS] {
            let fresh = rho_budget(bits, &tuning, false);
            let after = rho_budget(bits, &tuning, true);
            assert!(
                after >= wide_rho_budget(bits) && after > fresh * 10,
                "{bits}-bit cofactor was not deepened after a split: {fresh} -> {after}"
            );
        }
        // The line itself, and the range below it where a sieve run costs seconds.
        assert!(
            rho_budget(DEEP_RHO_MIN_BITS, &tuning, true)
                > rho_budget(DEEP_RHO_MIN_BITS, &tuning, false)
        );
        for bits in [128, 192, 224, 256] {
            assert_eq!(
                rho_budget(bits, &tuning, true),
                rho_budget(bits, &tuning, false),
                "{bits}-bit cofactor deepened where the sieve is cheaper"
            );
        }
        // Above the ceiling the budget is already the deep one, split or no split.
        assert_eq!(
            rho_budget(512, &tuning, true),
            rho_budget(512, &tuning, false)
        );
    }

    /// The override exists so that a 56- or 64-bit factor above the ceiling is reachable without a
    /// rebuild. It has to win at every width, including the ones with a computed budget.
    #[test]
    fn rho_iteration_override_replaces_both_budget_arms() {
        let tuning = FactorTuning {
            rho_iterations: Some(4_000_000_000),
            ..FactorTuning::default()
        };
        for bits in [128, MAX_SIQS_BITS, MAX_SIQS_BITS + 1, 1024] {
            assert_eq!(rho_budget(bits, &tuning, false), 4_000_000_000);
        }
        // Zero would disable the stage silently; the floor keeps it a budget rather than a switch.
        let zero = FactorTuning {
            rho_iterations: Some(0),
            ..FactorTuning::default()
        };
        assert_eq!(rho_budget(512, &zero, false), 1);
    }

    /// The family budget exists to stop a run that cannot finish, not to truncate one that can.
    /// It must never shrink as inputs get harder.
    #[test]
    fn family_budget_is_monotonic_in_input_width() {
        let mut previous = 0;
        for bits in [128, 256, 288, 289, 320, 368, 369, 400] {
            let budget = family_budget(bits);
            assert!(budget >= previous, "{bits}-bit budget regressed");
            previous = budget;
        }
        assert!(family_budget(MAX_SIQS_BITS) > family_budget(288));
    }

    /// `find_factor` used to carry its own copy of the whole of `prepare`, and the two copies had
    /// already diverged on this one expression: `prepare` translated sieve roots by
    /// `sieve_half_width % prime` while the native copy used `interval as u32 % prime`. They agree
    /// only because `interval` *is* the (positive) sieve half width, which is the invariant pinned
    /// here — together with the property the sieve actually depends on, namely that adding the
    /// precomputed translation to a root mod `p` equals translating the signed coordinate first.
    #[test]
    fn interval_translation_matches_the_signed_coordinate_shift() {
        let p = Natural::from_u64(18_446_744_073_709_551_557);
        let q = Natural::from_u64(18_446_744_073_709_551_533);
        let ctx = prepare(p.checked_mul(&q).unwrap(), &FactorTuning::default())
            .unwrap()
            .0;
        let params = crate::qs::parameters::engine_params(ctx.n.bit_len());
        assert_eq!(ctx.interval, params.sieve_half_width as i32);
        assert!(ctx.interval > 0);
        for (entry, &translation) in ctx.base.iter().zip(ctx.interval_mod_p.iter()) {
            let prime = entry.prime;
            assert_eq!(translation, params.sieve_half_width % prime);
            assert_eq!(translation, ctx.interval as u32 % prime);
            for x in [
                -(ctx.interval as i64),
                -1,
                0,
                1,
                (prime as i64) - 1,
                ctx.interval as i64 - 1,
            ] {
                let position = (x + ctx.interval as i64) as u64 % prime as u64;
                let xmod = x.rem_euclid(prime as i64) as u32;
                assert_eq!(
                    position as u32,
                    add_mod_u32(xmod, translation, prime),
                    "p={prime} x={x}"
                );
            }
        }
    }

    /// A deterministic `choose_a` famine must be reported as a parameter-selection failure and must
    /// cost nothing, rather than burning the whole family budget and reporting "no factor". The
    /// 65-85-bit dead zone made this indistinguishable from a genuine search failure in v0.2.0.
    #[test]
    fn polynomial_selection_famine_is_diagnosed_without_searching() {
        let base = [FactorBaseEntry {
            prime: 3,
            log_prime: 9,
            sqrt_n: 1,
        }];
        let target = Natural::from_u64(1 << 40);
        let (all, pool, count) = siqs::build_a_candidates(&base, &target);
        assert!(
            all.len() < count || pool.is_empty(),
            "a one-prime factor base cannot supply {count} coefficient factors"
        );
    }

    #[test]
    fn precomputed_remainders_and_root_translation_are_exact() {
        for p in 2u32..=10_000 {
            let c = lemire_c(p);
            for a in [0, 1, p - 1, p, p.saturating_add(1), u32::MAX] {
                assert_eq!(fastmod(a, p, c), a % p, "p={p}, a={a}");
            }
            for a in [0, 1, p - 1] {
                for b in [0, 1, p - 1] {
                    assert_eq!(add_mod_u32(a, b, p), (a as u64 + b as u64) as u32 % p);
                }
            }
        }
    }

    #[test]
    fn portable_jobs_are_deterministic() {
        let p = Natural::from_u64(18_446_744_073_709_551_557);
        let q = Natural::from_u64(18_446_744_073_709_551_533);
        let context = prepare(p.checked_mul(&q).unwrap(), &FactorTuning::default()).unwrap();
        let a = execute(&context, EngineJob { family: 7 });
        let b = execute(&context, EngineJob { family: 7 });
        assert_eq!(a.family, b.family);
        assert_eq!(a.polynomials, b.polynomials);
        assert_eq!(a.relations, b.relations);
        assert!(a.polynomials > 0);
    }

    /// `EngineSession` is the collector the WASM coordinator drives, and that coordinator numbers
    /// families itself rather than calling `take_jobs`. Duplicate-`A` families therefore have to be
    /// dropped where relations are ingested; dropping them at dispatch left the browser path
    /// ingesting identical relations twice, which makes duplicate matrix columns whose dependencies
    /// are all trivial, and extraction then reports no factor on an input that factors natively.
    #[cfg(any(unix, windows))]
    #[test]
    fn session_drops_duplicate_a_families_however_they_are_scheduled() {
        use std::str::FromStr;
        // 110-bit semiprime that produced 3 duplicate families out of 56 on the native path.
        let n = Natural::from_str("668319744971798315493259725219859").unwrap();
        let context = prepare(n, &FactorTuning::default()).unwrap();

        // Feed every family in order, as the coordinator does, and count what is accepted.
        let mut session = EngineSession::new(context.clone());
        let mut duplicates = 0;
        for family in 0..64u64 {
            let before = session.polynomials();
            session.submit(execute(&context, EngineJob { family }));
            if session.polynomials() == before {
                duplicates += 1;
            }
        }
        assert!(
            duplicates > 0,
            "this input is supposed to generate duplicate A values; the test has gone stale"
        );

        // Every accepted family must have contributed a distinct A.
        let mut seen = HashSet::new();
        for family in 0..64u64 {
            if let Some((a, _)) = siqs::choose_a(&context.0, family) {
                seen.insert(a);
            }
        }
        assert_eq!(
            session.polynomials(),
            seen.len() as u64 * (1 << (context.0.a_factor_count - 1).min(9)),
            "accepted polynomial count does not match the number of distinct A values"
        );
    }

    #[test]
    fn collector_accepts_out_of_order_results() {
        let p = Natural::from_u64(18_446_744_073_709_551_557);
        let q = Natural::from_u64(18_446_744_073_709_551_533);
        let context = prepare(p.checked_mul(&q).unwrap(), &FactorTuning::default()).unwrap();
        let mut session = EngineSession::new(context.clone());
        let jobs = session.take_jobs(2);
        session.submit(execute(&context, jobs[1]));
        assert_eq!(session.polynomials(), 0);
        session.submit(execute(&context, jobs[0]));
        assert!(session.polynomials() > 0);
    }

    #[cfg(any(unix, windows))]
    #[test]
    fn full_parallel_engine_factors_128_bit_semiprime() {
        let p = Natural::from_u64(18_446_744_073_709_551_557);
        let q = Natural::from_u64(18_446_744_073_709_551_533);
        let n = p.checked_mul(&q).unwrap();
        let factors = factor(n.clone(), 2, &FactorTuning::default(), None, |_| true).unwrap();
        assert_eq!(factors, [q, p]);
        assert_eq!(
            factors
                .iter()
                .try_fold(Natural::ONE, |a, b| a.checked_mul(b)),
            Some(n)
        );
    }

    /// Balanced semiprimes at 65, 70, 75, 80, 85 and 90 bits, whose cofactors also cover the band.
    const DEAD_ZONE: [&str; 6] = [
        "18446744400127067027",        // 65 bits, the minimum reproducer
        "635904368119925963561",       // 70 bits
        "20988451891514649258347",     // 75 bits
        "703713894016303629914563",    // 80 bits
        "22921914745054882120472087",  // 85 bits
        "648536833001811612107041493", // 90 bits
    ];

    /// Phase 1.1: `choose_a` could not build a candidate pool anywhere in the 65-85-bit band, because
    /// it drew from primes above 1000 (10 bits and up) while accepting only primes within one bit of
    /// `ideal_bits`, which is 6 there. The pool was empty for every family, so `polys` stayed 0.
    ///
    /// This has to be asserted against the engine directly. The Phase 1.3 rho stage now splits
    /// everything in this band before SIQS is reached, so an end-to-end `factor()` on a dead-zone
    /// input would succeed even with the bug fully present — exactly the masking the brief predicted.
    #[test]
    fn siqs_builds_polynomials_across_the_dead_zone() {
        use std::str::FromStr;
        for case in DEAD_ZONE {
            let n = Natural::from_str(case).unwrap();
            let context = prepare(n.clone(), &FactorTuning::default())
                .unwrap_or_else(|e| panic!("{case} ({} bits): {e}", n.bit_len()));
            assert!(
                !context.0.a_pool.is_empty(),
                "{case} ({} bits): empty A candidate pool",
                n.bit_len()
            );
            let result = execute(&context, EngineJob { family: 0 });
            assert!(
                result.polynomials > 0,
                "{case} ({} bits): no polynomials from family 0",
                n.bit_len()
            );
        }
    }

    /// The same band, end to end through SIQS with the rho stage bypassed, so a factor really is
    /// recovered from sieved relations rather than from the preceding ladder stage.
    #[cfg(any(unix, windows))]
    #[test]
    fn siqs_alone_factors_the_dead_zone() {
        use std::str::FromStr;
        for case in DEAD_ZONE {
            let n = Natural::from_str(case).unwrap();
            let d = find_factor(n.clone(), 2, &FactorTuning::default(), &mut |_| true)
                .unwrap_or_else(|e| panic!("{case} ({} bits): {e}", n.bit_len()));
            assert!(!d.is_one() && d != n, "{case}: trivial factor {d}");
            assert!(
                n.div_rem(&d).unwrap().1.is_zero(),
                "{case}: {d} does not divide it"
            );
        }
    }

    /// Phase 1.3: the bounded Pollard-Brent stage must split an unbalanced `N` — the case SIQS is
    /// worst at, since SIQS pays for the size of `N` while rho pays for the size of its smallest
    /// factor — without entering SIQS at all. Asserted through the stage counters rather than from
    /// the answers, because SIQS returns the same factors and would hide a dead rho stage.
    ///
    /// These are the unbalanced entries of the supplied corpus whose small factor is above the 10^4
    /// trial-division bound: 16, 20 and 24-bit factors of 127, 160 and 224-bit inputs. The 224-bit
    /// one takes 0.01 s through rho against about 5 s of sieving.
    #[cfg(any(unix, windows))]
    #[test]
    fn rho_stage_splits_unbalanced_inputs_without_entering_siqs() {
        use std::str::FromStr;
        let cases = [
            "88948294177717782578521953992989251229",
            "1185123569529286501965460691005493488051524107431",
            "13695626177198106295200293487798368178679518660650179392786377544541",
        ];
        for case in cases {
            let n = Natural::from_str(case).unwrap();
            stage_counts::reset();
            let factors = factor(n.clone(), 2, &FactorTuning::default(), None, |_| true).unwrap();
            assert_eq!(
                factors
                    .iter()
                    .try_fold(Natural::ONE, |acc, f| acc.checked_mul(f)),
                Some(n.clone()),
                "{case} did not factor back"
            );
            assert!(
                stage_counts::rho() > 0,
                "{case} ({} bits) did not reach the rho stage",
                n.bit_len()
            );
            assert_eq!(
                stage_counts::siqs(),
                0,
                "{case} ({} bits) entered SIQS",
                n.bit_len()
            );
        }
    }

    /// The other half of the previous test: on a balanced semiprime the rho stage must spend its
    /// budget and hand off, so the assertion above cannot pass by the counters being wired backwards.
    /// This is also what keeps the stage honest — an unbounded rho here measured 0.56 s at 80 bits
    /// and 2.16 s at 90 bits against 0.03 s and 0.01 s for SIQS alone.
    #[cfg(any(unix, windows))]
    #[test]
    fn balanced_semiprimes_fall_through_rho_to_siqs() {
        let p = Natural::from_u64(18_446_744_073_709_551_557);
        let q = Natural::from_u64(18_446_744_073_709_551_533);
        let n = p.checked_mul(&q).unwrap();
        stage_counts::reset();
        factor(n, 2, &FactorTuning::default(), None, |_| true).unwrap();
        assert!(stage_counts::siqs() > 0, "128-bit input skipped SIQS");
    }

    /// The budget is the total across every polynomial constant tried, not per constant. It was
    /// per-constant at first, which made the stage cost 8× its nominal budget — 27 s of a 64 s
    /// 256-bit run — while reporting the same number.
    #[cfg(any(unix, windows))]
    #[test]
    fn rho_respects_its_total_iteration_budget() {
        let p = Natural::from_u64(18_446_744_073_709_551_557);
        let q = Natural::from_u64(18_446_744_073_709_551_533);
        let n = p.checked_mul(&q).unwrap();
        let mut polls = 0usize;
        let started = std::time::Instant::now();
        let poll = || {
            polls += 1;
            true
        };
        let result = pollard_brent_natural(&n, 4_096, NATIVE_RHO_CONSTANTS, poll).unwrap();
        assert!(
            result.is_none(),
            "4096 iterations should not split a 128-bit balanced semiprime"
        );
        // Eight constants each running to a 4096-iteration budget would take multiple seconds in an
        // unoptimized test build; one shared budget is milliseconds.
        assert!(
            started.elapsed() < std::time::Duration::from_secs(5),
            "budget was not shared across constants ({polls} cancellation polls)"
        );
    }

    #[test]
    fn montgomery_brent_splits_fixed_cofactor_corpus() {
        for (left, right) in [
            (1_000_003u64, 1_000_033u64),
            (15_485_863, 15_485_867),
            (4_294_967_291, 4_294_967_279),
        ] {
            let n = left * right;
            let factor = pollard_u64(n).expect("Brent failed to split cofactor");
            assert!(factor == left || factor == right, "{n} -> {factor}");
        }
    }

    #[test]
    fn double_large_prime_classification_is_bounded_and_exact() {
        let left = 1_000_003u64;
        let right = 1_000_033u64;
        let product = left * right;
        assert!(matches!(
            classify_cofactor(left, 2_000_000, product),
            Some(LargePrime::One(value)) if value == left
        ));
        assert!(classify_cofactor(product, 2_000_000, 0).is_none());
        assert!(classify_cofactor(product, 2_000_000, product - 1).is_none());
        assert!(matches!(
            classify_cofactor(product, 2_000_000, product),
            Some(LargePrime::Two(a, b)) if a == left && b == right
        ));
        assert!(
            classify_cofactor(product, left, product).is_none(),
            "a factor above the single-prime bound was accepted"
        );
    }

    #[test]
    fn dlp_policy_is_integer_bounded() {
        let (single, double) = large_prime_policy(3_000_000, 150, 16);
        assert_eq!(single, 450_000_000);
        assert_eq!(double, 144_000_000_000_000);
        assert!(double < single.saturating_mul(single));
        assert_eq!(large_prime_policy(3_000_000, 150, 0), (single, 0));
    }

    #[test]
    fn double_large_prime_triangle_closes_and_retains_square_roots() {
        let n = Natural::from_u64(1_000_000_007);
        let relation = |a, b| Relation {
            root: Natural::ONE,
            sign: false,
            powers: Vec::new(),
            large: LargePrime::Two(a, b),
        };
        let mut collector = RelationCollector::new();
        collector.ingest(relation(101, 103), &n);
        collector.ingest(relation(103, 107), &n);
        assert!(collector.columns.is_empty());
        collector.ingest(relation(101, 107), &n);
        assert_eq!(collector.columns.len(), 1);
        assert_eq!(collector.cycles, 1);
        let mut square_roots = collector.columns[0].extra_sqrt.clone();
        square_roots.sort_unstable();
        assert_eq!(square_roots, [101, 103, 107]);
    }

    /// Iteration throughput of the big-integer rho stage at and above the sieve ceiling. This is
    /// what [`rho_budget`]'s above-ceiling arm is sized from: the budget there is a wall-clock
    /// decision, so it has to be stated in seconds, and seconds come from this measurement.
    #[cfg(any(unix, windows))]
    #[test]
    #[ignore = "manual wide-rho throughput measurement"]
    fn profile_wide_rho_throughput() {
        for bits in [128usize, 192, 256, 320, 400, 512, 768, 1024] {
            // A prime of the right width: every iteration costs what it costs on a composite of
            // the same size, and no split can cut the run short, so the whole budget is spent and
            // the timing is the loop's rather than luck's.
            let mut bytes = vec![0xa5u8; bits.div_ceil(8)];
            bytes[0] |= 1;
            let top = bytes.len() - 1;
            bytes[top] |= 0x80;
            let mut n = Natural::from_le_bytes(&bytes).expect("width fits the engine capacity");
            let two = Natural::from_u64(2);
            while !is_probable_prime(&n, &PrimalityConfig::default()) {
                n = n.checked_add(&two).expect("prime search stays in width");
            }
            assert_eq!(n.bit_len(), bits);
            const ITERATIONS: u64 = 200_000;
            let started = std::time::Instant::now();
            let split =
                pollard_brent_natural(&n, ITERATIONS, NATIVE_RHO_CONSTANTS, || true).unwrap();
            let elapsed = started.elapsed().as_secs_f64();
            eprintln!(
                "BENCH wide_rho bits={bits} iterations={ITERATIONS} elapsed={elapsed:.3}s \
                 rate={:.0}/s split={}",
                ITERATIONS as f64 / elapsed,
                split.is_some()
            );
        }
    }

    #[test]
    #[ignore = "manual cofactor-split performance measurement"]
    fn profile_pollard_u64() {
        let n = 134_217_689u64 * 134_217_757u64;
        let started = std::time::Instant::now();
        for _ in 0..1_000 {
            assert!(
                pollard_u64(std::hint::black_box(n))
                    .map(std::hint::black_box)
                    .is_some()
            );
        }
        eprintln!(
            "BENCH pollard_u64_1000={:.6}s",
            started.elapsed().as_secs_f64()
        );
    }
}
