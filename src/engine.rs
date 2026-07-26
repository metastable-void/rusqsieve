//! Portable SIQS engine and scheduler-facing work kernels.
mod extract;
mod siqs;
mod wire;

use crate::f2::SparseBinaryMatrix;
use crate::factor::FactorTuning;
#[cfg(any(unix, windows))]
use crate::natural::MontgomeryContext;
use crate::qs::{AutoOr, FactorBaseEntry, MultiplierChoice, QsConfig, prepare_factor_base};
use crate::{Natural, PARTS, jacobi_u64};
#[cfg(any(unix, windows))]
use crate::{PrimalityConfig, WitnessPolicy, is_probable_prime};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
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
    InsufficientRelations,
    NoFactor,
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
    base: Arc<[FactorBaseEntry]>,
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
    /// Sieve-threshold slack for the unfactored cofactor, in bits. This is `log2(single_limit)` and
    /// nothing else: a survivor whose cofactor exceeds the large-prime bound is discarded by
    /// `classify_cofactor`, so admitting it costs a full trial division for no possible relation.
    /// v0.2.0 used an independent per-tier `lp_allowance` here, which at the 256-bit tier admitted
    /// 34-bit cofactors against a 27-bit acceptance bound — seven bits of pure waste.
    lp_bits: usize,
    /// Measured per-tier sieve-threshold offset in bits (see [`crate::qs::parameters`]).
    thresh_adj: i32,
    /// Maximum accepted single large prime (and maximum factor of a double).
    single_limit: u64,
    /// Whether double-large-prime cofactors are captured and combined.
    double_enabled: bool,
    relation_percent: Option<usize>,
    small_skip: u32,
    threshold_margin: i32,
    profile: bool,
}

/// Large-prime cofactor content of a relation.
#[derive(Clone, Copy)]
enum LargePrime {
    None,
    One(u64),
    Two(u64, u64),
}
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
    /// `2·Bⱼ·a⁻¹ mod p` for each varying B-value `j` and factor-base prime `p`
    /// (row-major `[j*nfb + i]`). Adding/subtracting this advances the roots to
    /// the next self-initializing polynomial in O(1) per prime (SPEC §7.3).
    bainv: Vec<u32>,
    /// Positions surviving the score threshold, reused across polynomials.
    candidates: Vec<u32>,
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
fn rho_budget(bits: usize) -> u64 {
    let p = crate::qs::parameters::engine_params(bits);
    (p.factor_base_bound as u64 * p.sieve_half_width as u64 / 500_000).max(1_024)
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

/// Polynomial families any one scheduler will issue before giving up.
///
/// Both the native thread scheduler and [`EngineSession`] use this single bound; v0.2.0 had two
/// uncoordinated 100 000 caps for the same job. Reaching it means the relation target was not met
/// and surfaces as [`EngineError::InsufficientRelations`] rather than a silent wrong answer.
const MAX_FAMILIES: u64 = 100_000;

/// Prepare an immutable context without creating threads.
pub fn prepare(n: Natural, tuning: &FactorTuning) -> Result<EngineContext, EngineError> {
    let mut p = crate::qs::parameters::engine_params(n.bit_len());
    if let Some(bound) = tuning.factor_base_bound {
        p.factor_base_bound = bound;
    }
    if let Some(half_width) = tuning.sieve_half_width {
        p.sieve_half_width = half_width;
    }
    let k = knuth_schroeppel(&n);
    let sieve_n = n
        .checked_mul(&Natural::from_u64(k))
        .unwrap_or_else(|| n.clone());
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(p.factor_base_bound),
        multiplier: MultiplierChoice::Value(k as u32),
    };
    let prepared = prepare_factor_base(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
    let pinv: Arc<[u64]> = base.iter().map(|e| lemire_c(e.prime)).collect();
    // Expected log weight of the tiny primes that are skipped by the score pass: `Σ log(p)/(p−1)`
    // over the skipped set, which is what the threshold must give back. Derived once here rather
    // than carried as a hand-tuned constant, and cheap to keep out of the per-polynomial path.
    let small_skip = tuning.small_skip.unwrap_or(100);
    let small_slack = small_slack(&base, small_skip);
    let interval_mod_p: Arc<[u32]> = base.iter().map(|e| p.sieve_half_width % e.prime).collect();
    let target_a = sieve_n
        .floor_sqrt()
        .div_rem_u64(p.sieve_half_width as u64)
        .unwrap()
        .0;
    let (a_all, a_pool, a_factor_count) = siqs::build_a_candidates(&base, &target_a);
    let (single_limit, double_enabled) =
        large_prime_policy(p.factor_base_bound, p.large_prime_mult);
    let context = Arc::new(Context {
        n,
        sieve_n,
        base,
        pinv,
        small_slack,
        interval_mod_p,
        interval: p.sieve_half_width as i32,
        target_a,
        a_all,
        a_pool,
        a_factor_count,
        // A double needs room for two acceptable large primes, so the threshold has to admit twice
        // the cofactor width or no double can ever survive it — the second, independent blocker that
        // made `LargePrime::Two` unreachable in v0.2.0 even with its gate forced open.
        lp_bits: (64 - single_limit.leading_zeros()) as usize * if double_enabled { 2 } else { 1 },
        thresh_adj: p.thresh_adj + tuning.threshold_adjustment.unwrap_or(0),
        single_limit,
        double_enabled,
        relation_percent: tuning.relation_percent,
        small_skip,
        threshold_margin: tuning.threshold_margin.unwrap_or(0),
        profile: tuning.profile,
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

/// Whether to capture and combine double-large-prime cofactors.
///
/// Off, on measurement. `LargePrime::Two`, `classify_cofactor` and cycle combination are all
/// implemented and reachable — flipping this constant is the whole switch — but at 192-bit, 4
/// threads, sieve+collect seconds:
///
/// | configuration                                    | polys | survivors | cycles | time  |
/// |--------------------------------------------------|------:|----------:|-------:|------:|
/// | off (shipped)                                    | 5 632 |    45 194 |  1 945 | 0.605 |
/// | on, threshold matched to the single-prime bound   | 5 408 |    51 513 |  2 040 | 0.607 |
/// | on, threshold 4 bits deeper                       | 4 672 |    85 829 |  2 391 | 0.680 |
/// | on, threshold widened to admit genuine doubles    |     — |         — |      — | >300  |
///
/// At a matched threshold it buys 4% fewer polynomials and the extra cofactor splits eat exactly
/// that, so it is a wash. Widening the threshold to `2 · log2(large-prime bound)` — which a genuine
/// double *needs*, and which was the second, independent reason `LargePrime::Two` was unreachable in
/// v0.2.0 — floods the survivor path and did not finish in 300 s against 0.605 s.
const DOUBLE_LARGE_PRIMES: bool = false;

/// Large-prime acceptance is independent from the sieve threshold slack.
fn large_prime_policy(bound: u32, large_prime_mult: u32) -> (u64, bool) {
    (
        (bound as u64).saturating_mul(large_prime_mult as u64),
        DOUBLE_LARGE_PRIMES,
    )
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
        Self {
            context,
            target,
            next_job: 0,
            next_merge: 0,
            polynomials: 0,
            collector: RelationCollector::new(),
            buffered: BTreeMap::new(),
            seen_a: HashSet::new(),
        }
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
        while jobs.len() < maximum && self.next_job < MAX_FAMILIES {
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
        extract::extract(&self.context.0, &self.collector.columns)
    }
}

fn relation_target(base_len: usize, percent: Option<usize>) -> usize {
    if let Some(percent) = percent {
        // Fewer than one relation per factor-base row cannot produce the
        // dependencies SIQS needs, so reject misleading under-target tuning.
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
        &mut progress,
        &mut factors,
    )?;
    factors.sort();
    Ok(factors)
}

#[cfg(any(unix, windows))]
fn factor_node(
    n: Natural,
    threads: usize,
    pc: &PrimalityConfig,
    tuning: &FactorTuning,
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
        factor_node(root, threads, pc, tuning, progress, &mut fs)?;
        for _ in 0..k {
            out.extend(fs.iter().cloned())
        }
        return Ok(());
    }
    let d = match pollard_brent_natural(&n, rho_budget(n.bit_len()), || {
        progress(EngineProgress {
            phase: EnginePhase::Preprocessing,
            polynomials: 0,
            relations: 0,
            target: 0,
            workers: threads,
        })
    })? {
        Some(factor) => {
            #[cfg(test)]
            stage_counts::bump(&stage_counts::RHO);
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
    let q = n.div_rem(&d).unwrap().0;
    factor_node(d, threads, pc, tuning, progress, out)?;
    factor_node(q, threads, pc, tuning, progress, out)
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
    let ctx = prepare(n.clone(), tuning)?.0;
    let target = relation_target(ctx.base.len(), ctx.relation_percent);
    if prof {
        eprintln!(
            "PROFILE fb_build={:.3}s nfb={} interval={} target={}",
            t_fb.elapsed().as_secs_f64(),
            ctx.base.len(),
            ctx.interval,
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
                        if tx.send(siqs::sieve_family(&c, f, &mut scratch)).is_err() {
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
    let mut next_merge = 0u64;
    let mut outstanding = 0usize;
    for _ in 0..threads * 2 {
        job_tx
            .send(Some(next_send))
            .map_err(|_| EngineError::Worker("worker job channel disconnected".into()))?;
        next_send += 1;
        outstanding += 1
    }
    let t_sieve = std::time::Instant::now();
    let mut buffered = BTreeMap::new();
    let mut collector = RelationCollector::new();
    let mut polynomials = 0u64;
    let mut total_survivors = 0u64;
    let mut seen_a = HashSet::new();
    let mut cancelled = false;
    while collector.columns.len() < target && next_merge < MAX_FAMILIES && !cancelled {
        let result = res_rx
            .recv()
            .map_err(|_| EngineError::Worker("worker result channel disconnected".into()))?;
        outstanding -= 1;
        buffered.insert(result.family, result);
        while let Some(r) = buffered.remove(&next_merge) {
            next_merge += 1;
            let unique_a = siqs::choose_a(&ctx, r.family)
                .map(|(a, _)| seen_a.insert(a))
                .unwrap_or(false);
            if !unique_a {
                continue;
            }
            polynomials += r.polynomials;
            total_survivors += r.survivors;
            for rel in r.relations {
                collector.ingest(rel, &n);
                if collector.columns.len() >= target {
                    break;
                }
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
                break;
            }
        }
        while outstanding < threads * 2
            && next_send < MAX_FAMILIES
            && collector.columns.len() < target
        {
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
    if prof {
        eprintln!(
            "PROFILE sieve+collect={:.3}s polys={} families={} survivors={} partials={} cycles={} relations={}",
            t_sieve.elapsed().as_secs_f64(),
            polynomials,
            next_merge,
            total_survivors,
            collector.forest.relations.len(),
            collector.cycles,
            collector.columns.len()
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
fn small_slack(base: &[FactorBaseEntry], skip: u32) -> usize {
    base.iter()
        .filter(|entry| entry.prime < skip)
        .map(|entry| {
            score_weight(entry.prime) as f64 / (entry.prime.saturating_sub(1).max(1)) as f64
        })
        .sum::<f64>()
        .round() as usize
}

/// Sieve log weight of one factor-base prime: its bit length, i.e. `ceil(log2 p)`. Computed from the
/// prime rather than read from a parallel array, so the score pass streams one fewer array over the
/// whole factor base for every polynomial — which matters on the small-L2 targets this crate aims at.
#[inline(always)]
fn score_weight(prime: u32) -> u8 {
    (32 - prime.leading_zeros()) as u8
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

#[allow(clippy::too_many_arguments)]
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
) {
    for (idx, e) in ctx.base.iter().enumerate() {
        let p = e.prime;
        if p == 2 || p < small_skip {
            continue;
        }
        let pu = p as usize;
        let weight = score_weight(p);
        if root1[idx] != u32::MAX {
            sieve_root_pair::<SATURATING>(
                scores,
                root1[idx] as usize,
                root2[idx] as usize,
                pu,
                weight,
            );
        } else {
            // p | a: the polynomial is linear (2bx + c) mod p — one root, per poly.
            let mut bp = b.mod_u64(p as u64) as u32;
            if bneg && bp != 0 {
                bp = p - bp;
            }
            let denom = (2 * bp as u64 % p as u64) as u32;
            let Some(inv) = inv_u32(denom, p) else {
                continue;
            };
            let cm = c.mod_u64(p as u64) as u32;
            let signed_c = if csign && cm != 0 { p - cm } else { cm };
            let xroot = mulmod_u32(if signed_c == 0 { 0 } else { p - signed_c }, inv, p);
            let mut pos = add_mod_u32(xroot, ctx.interval_mod_p[idx], p) as usize;
            while pos < scores.len() {
                scores[pos] = if SATURATING {
                    scores[pos].saturating_add(weight)
                } else {
                    scores[pos].wrapping_add(weight)
                };
                pos += pu;
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

/// Collect high-scoring positions, ascending.
///
/// Survivors are rare — a handful of positions in a 180–640 KiB score array — so this is a pure
/// rejection loop, and the caller biases every score by `128 − threshold` precisely so that
/// "reached the threshold" becomes "high bit set". One `and` then rejects eight positions at a time.
/// The bias is applied by the caller rather than folded in here because a runtime bias would cost an
/// `add` and an `or` per word, which measured 2.7× slower than the masked test on the 224-bit case.
fn collect_candidates(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
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
    out: &mut Vec<Relation>,
) -> usize {
    let base = &ctx.base;
    let len = (ctx.interval as usize) * 2;
    let small_skip = ctx.small_skip;
    let bb = b.checked_mul(b).unwrap();
    let (c, csign) = if bb >= ctx.sieve_n {
        (bb.wrapping_sub(&ctx.sieve_n).div_rem(a).unwrap().0, false)
    } else {
        (ctx.sieve_n.wrapping_sub(&bb).div_rem(a).unwrap().0, true)
    };
    let g_bits = ctx.sieve_n.bit_len().saturating_sub(a.bit_len());
    // Score threshold: a survivor's scored-prime log weight must come within `lp_bits` of g(x), the
    // point past which the unfactored cofactor could no longer be an acceptable large prime.
    // `small_slack` gives back the tiny primes that are not scored, and `thresh_adj` is the measured
    // per-tier offset trading survivors against polynomials (see `qs::parameters`).
    let threshold = (g_bits as i32 - ctx.lp_bits as i32 - ctx.small_slack as i32
        + ctx.threshold_margin
        + ctx.thresh_adj)
        .clamp(1, u8::MAX as i32) as u8;
    // Factor-base entries with `prime < small_skip` occupy the low indices (the base is sorted
    // ascending). Those tiny primes are not sieved — gating them would waste a `fastmod` where a
    // direct divide is cheaper (they divide most survivors) — so they are divided out directly.
    let small_end = base.partition_point(|e| e.prime < small_skip);
    // Bias every score so that reaching `threshold` sets the byte's high bit; the candidate scan is
    // then one masked compare per eight positions (see `collect_candidates`).
    let bias = 128u8.saturating_sub(threshold);
    let scan_threshold = threshold.saturating_add(bias);
    // Wrapping score writes are one instruction cheaper than saturating ones, but they are only
    // sound when no position can carry past 255. Every scored prime is at least
    // `base[small_end].prime`, so at most `g_bits / floor(log2 p_min)` distinct scored primes can
    // divide one `g(x)`, and each adds at most `log2(p) + 1` — hence the bound below. When it does
    // not fit, fall back to saturating writes, which lose the exact-score trial-division shortcut
    // but cannot produce a false negative. This is the invariant that `RUSQSIEVE_SMALL_SKIP` could
    // previously break silently into wrapping overflow and false negatives.
    let smallest_scored = base.get(small_end).map_or(u32::MAX, |e| e.prime).max(3);
    let scored_bits = (32 - smallest_scored.leading_zeros() - 1).max(1);
    let score_bound = g_bits as u32 + g_bits as u32 / scored_bits + 1;
    let exact_scores = bias as u32 + score_bound <= u8::MAX as u32;
    scores.clear();
    scores.resize(len, bias);
    // Flat logarithmic sieve. The formerly gated cache-blocked kernel was
    // unreachable for every shipped interval and was slower when forced on.
    if exact_scores {
        score_polynomial::<false>(ctx, b, bneg, &c, csign, root1, root2, scores, small_skip);
    } else {
        score_polynomial::<true>(ctx, b, bneg, &c, csign, root1, root2, scores, small_skip);
    }
    collect_candidates(scores, scan_threshold, candidates);
    if candidates.is_empty() {
        return 0;
    }
    let survivors = candidates.len();

    let two_idx = base.iter().position(|e| e.prime == 2).map(|i| i as u32);
    let mut powers_scratch = Vec::new();

    for &posu in candidates.iter() {
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
        let x = pos as i64 - ctx.interval as i64;
        let xabs = x.unsigned_abs();
        let ax = a.checked_mul(&Natural::from_u64(xabs)).unwrap();
        // t = a·x + b, needed for the relation's square root.
        let (t, tneg) = siqs::signed_add(&ax, x < 0, b, bneg);
        // Value to factor: g(x) = Q(x)/a = a·x² + 2b·x + c, computed directly with
        // signs (c_math = ∓c per csign). This avoids the wide t² squaring and the
        // division by a — a is guaranteed to divide Q since b² ≡ n (mod a).
        let ax2 = ax.checked_mul(&Natural::from_u64(xabs)).unwrap();
        let two_bx = b
            .wrapping_add(b)
            .checked_mul(&Natural::from_u64(xabs))
            .unwrap();
        let (gx, gxneg) = siqs::signed_add(&ax2, false, &two_bx, bneg ^ (x < 0));
        let (mut q, sign) = siqs::signed_add(&gx, gxneg, &c, csign);
        if q.is_zero() {
            continue;
        }
        powers_scratch.clear();
        powers_scratch.extend(aidx.iter().copied().map(|i| (i, 1)));
        // Merge a divided-out exponent for factor-base index `i` into `powers`.
        let record = |i: u32, count: u16, powers: &mut Vec<(u32, u16)>| {
            if count == 0 {
                return;
            }
            if let Some(v) = powers.iter_mut().find(|v| v.0 == i) {
                v.1 += count;
            } else {
                powers.push((i, count));
            }
        };
        // Prime 2 (not sieved): strip via trailing zeros.
        if let Some(ti) = two_idx {
            let c2 = q.trailing_zeros();
            if c2 != 0 {
                q >>= c2;
                record(ti, c2 as u16, &mut powers_scratch);
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
            let count = divide_out(&mut q, p);
            record(i as u32, count, &mut powers_scratch);
        }
        // Primes dividing `a` (seeded at exponent 1, root1 == MAX so not gated) — divide directly.
        for &ai in aidx {
            let p = base[ai as usize].prime as u64;
            let count = divide_out(&mut q, p);
            if count != 0 && p >= small_skip as u64 {
                confirmed_score += u16::from(score_weight(p as u32));
            }
            record(ai, count, &mut powers_scratch);
        }
        // Scored factor base: gate each prime with the precomputed multiply-shift residue test and
        // divide only on a hit (FLINT `qsieve_evaluate_candidate`), stopping as soon as the recorded
        // score is accounted for or nothing is left to divide.
        //
        // msieve-style resieving was implemented here — replaying the root progressions of the
        // sparse tail over a candidate bitmap so those primes cost O(number of factors) instead of a
        // gated scan — and measured a net LOSS at every size and every cutoff tried (224-bit: gate
        // scan 2.40 -> 1.49 cpu-s but the resieve pass cost 1.59; 256-bit: 26.5 -> 12.0 against a
        // 17.1 cpu-s pass). The reason is structural: the gate below is only ~9 cycles per prime, so
        // resieving competes against an already-cheap test, and its cost falls per *polynomial*
        // while the saving accrues per *survivor* — and this engine sees only ~5 survivors per
        // polynomial. See CHANGELOG 0.2.1.
        for idx in small_end..base.len() {
            if q.is_one() || confirmed_score >= score_target {
                break;
            }
            let r1 = root1[idx];
            if r1 == u32::MAX {
                continue; // prime divides `a`; handled above
            }
            let p = base[idx].prime;
            let position_mod_p = fastmod(posu, p, ctx.pinv[idx]);
            if position_mod_p != r1 && position_mod_p != root2[idx] {
                continue;
            }
            let count = divide_out(&mut q, p as u64);
            if count != 0 {
                confirmed_score += u16::from(score_weight(p));
            }
            record(idx as u32, count, &mut powers_scratch);
        }
        let large = if q.is_one() {
            LargePrime::None
        } else if q.bit_len() > 64 {
            continue;
        } else {
            match classify_cofactor(q.as_parts()[0], ctx.single_limit, ctx.double_enabled) {
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
            powers: core::mem::take(&mut powers_scratch),
            large,
        });
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
fn classify_cofactor(q: u64, single_limit: u64, double_enabled: bool) -> Option<LargePrime> {
    if crate::u64math::is_prime(q) {
        return (q <= single_limit).then_some(LargePrime::One(q));
    }
    if !double_enabled {
        return None;
    }
    let d = pollard_u64(q)?;
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
/// polynomial constant tried, not per constant, so the caller's budget is the whole cost of the
/// stage; see [`rho_budget`].
#[cfg(any(unix, windows))]
fn pollard_brent_natural(
    n: &Natural,
    iteration_limit: u64,
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
    for c_value in 1..=8u64 {
        if iterations >= iteration_limit {
            return Ok(None);
        }
        if !keep_going() {
            return Err(EngineError::Cancelled);
        }
        let c = montgomery.encode(&Natural::from_u64(c_value));
        let mut y = montgomery.encode(&Natural::from_u64(2));
        let mut r = 1u64;
        let mut g = Natural::ONE;
        let mut x = Natural::ZERO;
        let mut ys = Natural::ZERO;
        while g.is_one() && iterations < iteration_limit {
            x = y.clone();
            for _ in 0..r {
                y = montgomery.add(&montgomery.square(&y), &c);
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
                ys = y.clone();
                let mut q = montgomery.one();
                let batch = (r - k).min(128);
                for _ in 0..batch {
                    y = montgomery.add(&montgomery.square(&y), &c);
                    let difference = if x >= y {
                        x.wrapping_sub(&y)
                    } else {
                        y.wrapping_sub(&x)
                    };
                    if !difference.is_zero() {
                        q = montgomery.multiply(&q, &difference);
                    }
                    iterations += 1;
                    if iterations >= iteration_limit {
                        break;
                    }
                }
                g = q.gcd(n);
                k += batch;
            }
            r = r.saturating_mul(2);
        }
        if g == *n {
            loop {
                if !keep_going() {
                    return Err(EngineError::Cancelled);
                }
                ys = montgomery.add(&montgomery.square(&ys), &c);
                let difference = if x >= ys {
                    x.wrapping_sub(&ys)
                } else {
                    ys.wrapping_sub(&x)
                };
                g = difference.gcd(n);
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
    id_of: HashMap<u64, u32>,
    parent: Vec<u32>,
    edge: Vec<Option<u32>>,
    relations: Vec<Relation>,
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
        let mut chain = vec![v];
        let mut edges: Vec<u32> = Vec::new();
        let mut c = v;
        while self.parent[c as usize] != c {
            edges.push(self.edge[c as usize].unwrap());
            c = self.parent[c as usize];
            chain.push(c);
        }
        self.parent[v as usize] = v;
        self.edge[v as usize] = None;
        for (i, e) in edges.into_iter().enumerate() {
            self.parent[chain[i + 1] as usize] = chain[i];
            self.edge[chain[i + 1] as usize] = Some(e);
        }
    }
    fn link(&mut self, a: u32, b: u32, rel: Relation) {
        self.reroot(b);
        let relation_index = self.relations.len() as u32;
        self.relations.push(rel);
        self.parent[b as usize] = a;
        self.edge[b as usize] = Some(relation_index);
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
}
impl RelationCollector {
    fn new() -> Self {
        Self {
            forest: Forest::default(),
            columns: Vec::new(),
            cycles: 0,
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
            let mut path = Vec::new();
            self.forest.path(va, &mut path);
            self.forest.path(vb, &mut path);
            self.cycles += 1;
            self.columns.push(combine_cycle(
                core::iter::once(&rel).chain(
                    path.iter()
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
    let sum = a as u64 + b as u64;
    if sum >= p as u64 {
        (sum - p as u64) as u32
    } else {
        sum as u32
    }
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
/// Knuth-Schroeppel multiplier selection. Chooses a small `k` such that `k·n` is a quadratic
/// residue modulo many small primes, raising the density of smooth `Q(x)` values (a standard
/// 2–3× QS speed-up). Ported from FLINT's `qsieve_knuth_schroeppel`. Returns `k` (>= 1).
fn knuth_schroeppel(n: &Natural) -> u64 {
    const MULTIPLIERS: [u64; 29] = [
        1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37,
        38, 41, 42, 43, 47,
    ];
    const KS_PRIMES: usize = 500;
    let nmod8 = n.mod_u64(8);
    let mut weights = [0.0f64; MULTIPLIERS.len()];
    for (w, &k) in weights.iter_mut().zip(&MULTIPLIERS) {
        let mod8 = (nmod8 * k) % 8;
        let mut v = 0.346_573_59_f64; // ln2 / 2
        if mod8 == 1 {
            v *= 4.0;
        } else if mod8 == 5 {
            v *= 2.0;
        }
        *w = v - (k as f64).ln() / 2.0;
    }
    // Weight each multiplier by the small primes for which `k·n` is a quadratic residue.
    let mut p = 3u64;
    let mut seen = 0usize;
    while seen < KS_PRIMES {
        if crate::u64math::is_prime(p) {
            seen += 1;
            let nmod = n.mod_u64(p);
            if nmod != 0 {
                let logpdivp = (p as f64).ln() / p as f64;
                let kron = jacobi_u64(nmod, p) as i32; // (n / p), handles even nmod
                for (w, &k) in weights.iter_mut().zip(&MULTIPLIERS) {
                    let km = k % p;
                    if km == 0 {
                        *w += logpdivp; // p | k → k·n ≡ 0 (mod p)
                    } else if kron * jacobi_u64(km, p) as i32 == 1 {
                        *w += 2.0 * logpdivp; // k·n is a QR mod p
                    }
                }
            }
        }
        p += 2;
    }
    let mut best = f64::NEG_INFINITY;
    let mut k = 1u64;
    for (&w, &m) in weights.iter().zip(&MULTIPLIERS) {
        if w > best {
            best = w;
            k = m;
        }
    }
    k
}
#[cfg(test)]
mod tests {
    use super::*;

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
        let result = pollard_brent_natural(&n, 4_096, || {
            polls += 1;
            true
        })
        .unwrap();
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
