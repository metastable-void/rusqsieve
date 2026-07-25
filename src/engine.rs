//! Portable SIQS engine and scheduler-facing work kernels.
use crate::f2::SparseBinaryMatrix;
use crate::qs::{AutoOr, FactorBaseEntry, MultiplierChoice, QsConfig, prepare_siqs};
use crate::{Natural, PARTS, jacobi_u64};
#[cfg(any(unix, windows))]
use crate::{PrimalityConfig, is_probable_prime};
use std::collections::{BTreeMap, HashMap};
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
    Worker,
    InvalidDependency,
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
    /// `interval mod p` per factor-base prime. Sieve roots are residues of the signed polynomial
    /// coordinate `x`, while score-array positions represent `x + interval`; precomputing this
    /// fixed translation avoids two signed divisions per prime and polynomial in the sieve pass.
    interval_mod_p: Arc<[u32]>,
    interval: i32,
    target_a: Natural,
    lp_allowance: usize,
    /// Maximum accepted single large prime (and maximum factor of a double).
    single_limit: u64,
    /// Whether double-large-prime cofactors are captured and combined.
    double_enabled: bool,
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

/// Per-worker reusable buffers (SPEC §21.1 — reuse sieve/candidate scratch).
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
    /// the next self-initializing polynomial in O(1) per prime (SPEC §12.5).
    bainv: Vec<u32>,
    /// Positions surviving the score threshold, reused across polynomials.
    candidates: Vec<u32>,
    /// Absolute carried hit positions for the cache-blocked sieve. Allocated only when the score
    /// interval is at least two blocks and reused across polynomials.
    block_pos1: Vec<usize>,
    block_pos2: Vec<usize>,
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

impl EngineJobResult {
    /// Serialize this family's relations for transport to a coordinator (e.g. from a
    /// Web Worker back to the main thread). Format is little-endian:
    /// `family:u64, polynomials:u64, count:u32`, then per relation
    /// `root:PARTS×u64, sign:u8, large:{tag:u8, 0/1/2 × u64}, powers_len:u32, [index:u32, exp:u16]…`.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut v = Vec::new();
        v.extend_from_slice(&self.inner.family.to_le_bytes());
        v.extend_from_slice(&self.inner.polynomials.to_le_bytes());
        v.extend_from_slice(&(self.inner.relations.len() as u32).to_le_bytes());
        for r in &self.inner.relations {
            for limb in r.root.as_parts() {
                v.extend_from_slice(&limb.to_le_bytes());
            }
            v.push(r.sign as u8);
            match r.large {
                LargePrime::None => v.push(0),
                LargePrime::One(a) => {
                    v.push(1);
                    v.extend_from_slice(&a.to_le_bytes());
                }
                LargePrime::Two(a, b) => {
                    v.push(2);
                    v.extend_from_slice(&a.to_le_bytes());
                    v.extend_from_slice(&b.to_le_bytes());
                }
            }
            v.extend_from_slice(&(r.powers.len() as u32).to_le_bytes());
            for &(i, e) in &r.powers {
                v.extend_from_slice(&i.to_le_bytes());
                v.extend_from_slice(&e.to_le_bytes());
            }
        }
        v
    }
}

/// Inverse of [`EngineJobResult::to_bytes`].
fn deserialize_family(b: &[u8]) -> Option<FamilyResult> {
    struct Cur<'a> {
        b: &'a [u8],
        o: usize,
    }
    impl Cur<'_> {
        fn take(&mut self, n: usize) -> Option<&[u8]> {
            let s = self.b.get(self.o..self.o + n)?;
            self.o += n;
            Some(s)
        }
        fn u8(&mut self) -> Option<u8> {
            Some(self.take(1)?[0])
        }
        fn u16(&mut self) -> Option<u16> {
            Some(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
        }
        fn u32(&mut self) -> Option<u32> {
            Some(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
        }
        fn u64(&mut self) -> Option<u64> {
            Some(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
        }
    }
    let mut c = Cur { b, o: 0 };
    let family = c.u64()?;
    let polynomials = c.u64()?;
    let count = c.u32()? as usize;
    let mut relations = Vec::with_capacity(count.min(1 << 20));
    for _ in 0..count {
        let root = Natural::from_le_bytes(c.take(PARTS * 8)?).ok()?;
        let sign = c.u8()? != 0;
        let large = match c.u8()? {
            0 => LargePrime::None,
            1 => LargePrime::One(c.u64()?),
            2 => LargePrime::Two(c.u64()?, c.u64()?),
            _ => return None,
        };
        let plen = c.u32()? as usize;
        let mut powers = Vec::with_capacity(plen.min(1 << 16));
        for _ in 0..plen {
            let i = c.u32()?;
            let e = c.u16()?;
            powers.push((i, e));
        }
        relations.push(Relation {
            root,
            sign,
            powers,
            large,
        });
    }
    Some(FamilyResult {
        family,
        polynomials,
        relations,
        survivors: 0,
    })
}

/// Prepare an immutable context without creating threads.
pub fn prepare(n: Natural) -> Result<EngineContext, EngineError> {
    let p = crate::qs::parameters::engine_params(n.bit_len());
    let k = knuth_schroeppel(&n);
    let sieve_n = n
        .checked_mul(&Natural::from_u64(k))
        .unwrap_or_else(|| n.clone());
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(p.factor_base_bound),
        multiplier: MultiplierChoice::Value(k as u32),
        ..QsConfig::default()
    };
    let prepared = prepare_siqs(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
    let pinv: Arc<[u64]> = base.iter().map(|e| lemire_c(e.prime)).collect();
    let interval_mod_p: Arc<[u32]> = base.iter().map(|e| p.sieve_half_width % e.prime).collect();
    let target_a = sieve_n
        .floor_sqrt()
        .div_rem_u64(p.sieve_half_width as u64)
        .unwrap()
        .0;
    let (single_limit, double_enabled) = large_prime_policy(p.factor_base_bound, p.lp_allowance);
    Ok(EngineContext(Arc::new(Context {
        n,
        sieve_n,
        base,
        pinv,
        interval_mod_p,
        interval: p.sieve_half_width as i32,
        target_a,
        lp_allowance: p.lp_allowance,
        single_limit,
        double_enabled,
    })))
}

/// Large-prime acceptance policy derived from the cofactor budget `lp_allowance`
/// (bits) and the factor-base bound. Doubles are only enabled when the budget can
/// hold two primes each above the factor base.
fn large_prime_policy(bound: u32, lp_allowance: usize) -> (u64, bool) {
    let single_limit = 1u64 << lp_allowance.min(62);
    let bound_bits = 64 - (bound as u64).max(1).leading_zeros();
    let double_enabled = lp_allowance as u32 >= 2 * bound_bits + 2;
    (single_limit, double_enabled)
}

/// Execute a job using only the caller's thread and owned scratch memory.
pub fn execute(context: &EngineContext, job: EngineJob) -> EngineJobResult {
    let mut scratch = EngineScratch::default();
    let inner = sieve_family(&context.0, job.family, &mut scratch);
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
}
impl EngineSession {
    pub fn new(context: EngineContext) -> Self {
        let target = relation_target(context.0.base.len());
        Self {
            context,
            target,
            next_job: 0,
            next_merge: 0,
            polynomials: 0,
            collector: RelationCollector::new(),
            buffered: BTreeMap::new(),
        }
    }
    pub fn take_jobs(&mut self, maximum: usize) -> Vec<EngineJob> {
        if self.is_ready() {
            return Vec::new();
        }
        (0..maximum)
            .map(|_| {
                let j = EngineJob {
                    family: self.next_job,
                };
                self.next_job += 1;
                j
            })
            .collect()
    }
    pub fn submit(&mut self, result: EngineJobResult) {
        self.buffered.insert(result.family, result.inner);
        self.drain_buffered();
    }
    /// Submit a worker's serialized [`EngineJobResult`] (see [`EngineJobResult::to_bytes`]).
    /// Returns whether enough relations have now been collected. Used by the WASM/Web-Worker
    /// scheduler to feed relations sieved in other threads back into the coordinator.
    pub fn submit_bytes(&mut self, bytes: &[u8]) -> bool {
        if let Some(fr) = deserialize_family(bytes) {
            self.buffered.insert(fr.family, fr);
            self.drain_buffered();
        }
        self.is_ready()
    }
    fn drain_buffered(&mut self) {
        while let Some(r) = self.buffered.remove(&self.next_merge) {
            self.next_merge += 1;
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
        extract(&self.context.0, &self.collector.columns)
    }
}

fn extract(ctx: &Context, columns: &[Column]) -> Result<Natural, EngineError> {
    let matrix_cols: Vec<Vec<u32>> = columns
        .iter()
        .map(|c| {
            let mut v = Vec::new();
            if c.sign {
                v.push(0)
            }
            for &(i, e) in &c.powers {
                if e & 1 != 0 {
                    v.push(i + 1)
                }
            }
            v
        })
        .collect();
    let matrix = SparseBinaryMatrix::from_columns(ctx.base.len() + 1, &matrix_cols)
        .map_err(|_| EngineError::InvalidDependency)?;
    for dep in matrix.filtered_dependencies().iter() {
        if !matrix.verify_dependency(dep) {
            return Err(EngineError::InvalidDependency);
        }
        let mut x = Natural::ONE;
        let mut y = Natural::ONE;
        let mut sums = vec![0u32; ctx.base.len()];
        for (j, c) in columns.iter().enumerate() {
            if (dep[j / 64] >> (j % 64)) & 1 == 0 {
                continue;
            }
            x = x.mul_mod(&c.root, &ctx.n);
            for &lp in &c.extra_sqrt {
                y = y.mul_mod(&Natural::from_u64(lp), &ctx.n);
            }
            for &(i, e) in &c.powers {
                sums[i as usize] += e
            }
        }
        for (e, &s) in ctx.base.iter().zip(&sums) {
            for _ in 0..s / 2 {
                y = y.mul_mod(&Natural::from_u64(e.prime as u64), &ctx.n)
            }
        }
        let d = if x >= y {
            x.wrapping_sub(&y)
        } else {
            y.wrapping_sub(&x)
        };
        let g = d.gcd(&ctx.n);
        if !g.is_one() && g != ctx.n {
            return Ok(g);
        }
        let g = x.add_mod(&y, &ctx.n).gcd(&ctx.n);
        if !g.is_one() && g != ctx.n {
            return Ok(g);
        }
    }
    Err(EngineError::NoFactor)
}

fn relation_target(base_len: usize) -> usize {
    #[cfg(any(unix, windows))]
    if let Some(percent) = std::env::var("RUSQSIEVE_REL_PERCENT")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
    {
        return (base_len * percent.clamp(50, 110) / 100).max(64);
    }
    base_len + 64
}

#[cfg(any(unix, windows))]
pub fn factor(
    mut n: Natural,
    threads: usize,
    mut progress: impl FnMut(EngineProgress) -> bool,
) -> Result<Vec<Natural>, EngineError> {
    if n.is_zero() {
        return Err(EngineError::Setup("zero has no prime factorization".into()));
    }
    let primality = PrimalityConfig::default();
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
    factor_node(n, threads.max(1), &primality, &mut progress, &mut factors)?;
    factors.sort();
    Ok(factors)
}

#[cfg(any(unix, windows))]
fn factor_node(
    n: Natural,
    threads: usize,
    pc: &PrimalityConfig,
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
        crate::smallfactor::factor_u64(v, &mut small);
        out.extend(small.into_iter().map(Natural::from_u64));
        return Ok(());
    }
    if is_probable_prime(&n, pc) {
        out.push(n);
        return Ok(());
    }
    if let Some((root, k)) = n.perfect_power() {
        let mut fs = Vec::new();
        factor_node(root, threads, pc, progress, &mut fs)?;
        for _ in 0..k {
            out.extend(fs.iter().cloned())
        }
        return Ok(());
    }
    let d = find_factor(n.clone(), threads, progress)?;
    if d.is_one() || d == n {
        return Err(EngineError::NoFactor);
    }
    let q = n.div_rem(&d).unwrap().0;
    factor_node(d, threads, pc, progress, out)?;
    factor_node(q, threads, pc, progress, out)
}

#[cfg(any(unix, windows))]
fn find_factor(
    n: Natural,
    threads: usize,
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
    let p = crate::qs::parameters::engine_params(n.bit_len());
    let bound = p.factor_base_bound;
    let interval = p.sieve_half_width as i32;
    let prof = std::env::var_os("RUSQSIEVE_PROFILE").is_some();
    let t_fb = std::time::Instant::now();
    let k = knuth_schroeppel(&n);
    let sieve_n = n
        .checked_mul(&Natural::from_u64(k))
        .unwrap_or_else(|| n.clone());
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(bound),
        multiplier: MultiplierChoice::Value(k as u32),
        ..QsConfig::default()
    };
    let prepared = prepare_siqs(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
    let pinv: Arc<[u64]> = base.iter().map(|e| lemire_c(e.prime)).collect();
    let interval_mod_p: Arc<[u32]> = base.iter().map(|e| interval as u32 % e.prime).collect();
    let target = relation_target(base.len());
    if prof {
        eprintln!(
            "PROFILE fb_build={:.3}s nfb={} interval={} target={} k={}",
            t_fb.elapsed().as_secs_f64(),
            base.len(),
            interval,
            target,
            k
        );
    }
    let target_a = sieve_n.floor_sqrt().div_rem_u64(interval as u64).unwrap().0;
    let (single_limit, double_enabled) = large_prime_policy(bound, p.lp_allowance);
    let ctx = Arc::new(Context {
        n: n.clone(),
        sieve_n,
        base: base.clone(),
        pinv,
        interval_mod_p,
        interval,
        target_a,
        lp_allowance: p.lp_allowance,
        single_limit,
        double_enabled,
    });
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
                let job = rx.lock().unwrap().recv();
                match job {
                    Ok(Some(f)) => {
                        if cancellation.load(AtomicOrdering::Relaxed) {
                            break;
                        }
                        if tx.send(sieve_family(&c, f, &mut scratch)).is_err() {
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
            .map_err(|_| EngineError::Worker)?;
        next_send += 1;
        outstanding += 1
    }
    let t_sieve = std::time::Instant::now();
    let mut buffered = BTreeMap::new();
    let mut collector = RelationCollector::new();
    let mut polynomials = 0u64;
    let mut total_survivors = 0u64;
    let mut cancelled = false;
    let max_families = 100_000u64;
    while collector.columns.len() < target && next_merge < max_families && !cancelled {
        let result = res_rx.recv().map_err(|_| EngineError::Worker)?;
        outstanding -= 1;
        buffered.insert(result.family, result);
        while let Some(r) = buffered.remove(&next_merge) {
            next_merge += 1;
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
            && next_send < max_families
            && collector.columns.len() < target
        {
            job_tx
                .send(Some(next_send))
                .map_err(|_| EngineError::Worker)?;
            next_send += 1;
            outstanding += 1
        }
    }
    for _ in 0..threads {
        let _ = job_tx.send(None);
    }
    drop(job_tx);
    for h in handles {
        let _ = h.join();
    }
    if cancelled {
        return Err(EngineError::Cancelled);
    }
    if prof {
        eprintln!(
            "PROFILE sieve+collect={:.3}s polys={} families={} survivors={} relations={}",
            t_sieve.elapsed().as_secs_f64(),
            polynomials,
            next_merge,
            total_survivors,
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
    let result = extract(&ctx, &collector.columns);
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

fn sieve_family(ctx: &Context, family: u64, scratch: &mut EngineScratch) -> FamilyResult {
    let empty = |family| FamilyResult {
        family,
        polynomials: 0,
        relations: Vec::new(),
        survivors: 0,
    };
    let Some((a, aidx)) = choose_a(ctx, family) else {
        return empty(family);
    };
    let base = &ctx.base;
    let nfb = base.len();
    let s = aidx.len();
    let nvar = (s - 1).min(6); // number of sign bits varied per family
    let variants = 1u64 << nvar;

    // SIQS B-values: b = Σ ±Bⱼ, with Bⱼ ≡ ±sqrt(n) (mod qⱼ), 0 (mod other q).
    // Keep the true signed B instead of reducing it modulo A. This is the standard self-init
    // representation used by FLINT and makes every Gray-code root update one conditional
    // add/subtract, without a per-prime correction for modular-A wraps.
    let mut bvals: Vec<Natural> = Vec::with_capacity(s);
    for &i in &aidx {
        let q = base[i as usize].prime;
        let Some((ap, _)) = a.div_rem_u64(q as u64) else {
            return empty(family);
        };
        let Some(apinv) = inv_u32(ap.mod_u64(q as u64) as u32, q) else {
            return empty(family);
        };
        let mut coeff = (base[i as usize].sqrt_n as u64 * apinv as u64) % q as u64;
        coeff = coeff.min(q as u64 - coeff);
        bvals.push(ap.checked_mul(&Natural::from_u64(coeff)).unwrap());
    }
    let mut b = Natural::ZERO;
    for bj in &bvals {
        b = b.checked_add(bj).unwrap();
    }
    let mut bneg = false;
    let two_full: Vec<Natural> = bvals[..nvar].iter().map(|bj| bj.wrapping_add(bj)).collect();

    // Per-prime precompute for the initial polynomial: both roots and, for each
    // varying B-value, the O(1) root advance `2·Bⱼ·a⁻¹ mod p`.
    scratch.root1.clear();
    scratch.root1.resize(nfb, u32::MAX);
    scratch.root2.clear();
    scratch.root2.resize(nfb, 0);
    scratch.bainv.clear();
    scratch.bainv.resize(nvar * nfb, 0);
    for (idx, e) in base.iter().enumerate() {
        let p = e.prime;
        if p == 2 {
            continue;
        }
        let ap = a.mod_u64(p as u64) as u32;
        if ap == 0 {
            continue; // p | a: linear fallback per polynomial (root1 stays MAX)
        }
        let Some(ainvp) = inv_u32(ap, p) else {
            continue;
        };
        let mut bp = b.mod_u64(p as u64) as u32;
        if bneg && bp != 0 {
            bp = p - bp;
        }
        let xroot1 = mulmod_u32((e.sqrt_n + p - bp) % p, ainvp, p);
        let xroot2 = mulmod_u32(((p - e.sqrt_n) % p + p - bp) % p, ainvp, p);
        let r1 = add_mod_u32(xroot1, ctx.interval_mod_p[idx], p);
        let r2 = add_mod_u32(xroot2, ctx.interval_mod_p[idx], p);
        scratch.root1[idx] = r1.min(r2);
        scratch.root2[idx] = r1.max(r2);
        for (j, bj) in bvals.iter().take(nvar).enumerate() {
            let bjp = bj.mod_u64(p as u64) as u32;
            let two_bjp = (2 * bjp as u64 % p as u64) as u32;
            scratch.bainv[j * nfb + idx] = mulmod_u32(two_bjp, ainvp, p);
        }
    }

    // Sieve every polynomial in Gray-code order, advancing the roots in O(1) per
    // prime between consecutive polynomials instead of recomputing them.
    let mut relations = Vec::new();
    let mut survivors = 0u64;
    for v in 0..variants {
        survivors += sieve_one_poly(
            ctx,
            &a,
            &b,
            bneg,
            &aidx,
            &scratch.root1,
            &scratch.root2,
            &mut scratch.scores,
            &mut scratch.candidates,
            &mut scratch.block_pos1,
            &mut scratch.block_pos2,
            &mut relations,
        ) as u64;
        if v + 1 >= variants {
            break;
        }
        let j = (v + 1).trailing_zeros() as usize;
        let gray = v ^ (v >> 1);
        let flip_to_one = (gray >> j) & 1 == 0;
        let add_bainv = if flip_to_one {
            (b, bneg) = signed_add(&b, bneg, &two_full[j], true);
            true
        } else {
            (b, bneg) = signed_add(&b, bneg, &two_full[j], false);
            false
        };
        let off = j * nfb;
        if add_bainv {
            for idx in 0..nfb {
                if scratch.root1[idx] == u32::MAX {
                    continue;
                }
                let p = base[idx].prime;
                let d = scratch.bainv[off + idx];
                let r1 = add_mod_u32(scratch.root1[idx], d, p);
                let r2 = add_mod_u32(scratch.root2[idx], d, p);
                scratch.root1[idx] = r1.min(r2);
                scratch.root2[idx] = r1.max(r2);
            }
        } else {
            for idx in 0..nfb {
                if scratch.root1[idx] == u32::MAX {
                    continue;
                }
                let p = base[idx].prime;
                let d = scratch.bainv[off + idx];
                let r1 = sub_mod_u32(scratch.root1[idx], d, p);
                let r2 = sub_mod_u32(scratch.root2[idx], d, p);
                scratch.root1[idx] = r1.min(r2);
                scratch.root2[idx] = r1.max(r2);
            }
        }
    }
    FamilyResult {
        family,
        polynomials: variants,
        relations,
        survivors,
    }
}

fn signed_add(a: &Natural, aneg: bool, b: &Natural, bneg: bool) -> (Natural, bool) {
    if aneg == bneg {
        let sum = a.checked_add(b).expect("signed SIQS coefficient overflow");
        let neg = aneg && !sum.is_zero();
        (sum, neg)
    } else if a >= b {
        let diff = a.wrapping_sub(b);
        let neg = aneg && !diff.is_zero();
        (diff, neg)
    } else {
        let diff = b.wrapping_sub(a);
        let neg = bneg && !diff.is_zero();
        (diff, neg)
    }
}

fn choose_a(ctx: &Context, family: u64) -> Option<(Natural, Vec<u32>)> {
    let all: Vec<usize> = ctx
        .base
        .iter()
        .enumerate()
        .filter(|(_, e)| e.prime > 1000)
        .map(|(i, _)| i)
        .collect();
    if all.len() < 8 {
        return None;
    }
    let target_bits = ctx.target_a.bit_len();
    let factor_count = target_bits.div_ceil(14).clamp(3, 10);
    let ideal_bits = target_bits.div_ceil(factor_count);
    let pool: Vec<usize> = all
        .iter()
        .copied()
        .filter(|&i| {
            let bits = (32 - ctx.base[i].prime.leading_zeros()) as usize;
            bits.abs_diff(ideal_bits) <= 1
        })
        .collect();
    if pool.len() < factor_count * 2 {
        return None;
    }
    let mut state = family ^ 0x9e3779b97f4a7c15;
    let mut a = Natural::ONE;
    let mut idx = Vec::with_capacity(factor_count);
    while idx.len() + 1 < factor_count {
        state = xorshift(state);
        let i = pool[state as usize % pool.len()];
        if idx.contains(&(i as u32)) {
            continue;
        }
        let next = a.checked_mul(&Natural::from_u64(ctx.base[i].prime as u64))?;
        a = next;
        idx.push(i as u32)
    }
    let desired = ctx.target_a.div_rem(&a)?.0;
    let desired_u64 = desired.as_parts()[0];
    let last = all
        .iter()
        .copied()
        .filter(|&i| !idx.contains(&(i as u32)))
        .min_by_key(|&i| (ctx.base[i].prime as u64).abs_diff(desired_u64))?;
    a = a.checked_mul(&Natural::from_u64(ctx.base[last].prime as u64))?;
    idx.push(last as u32);
    Some((a, idx))
}

/// Tiny-prime skipping (audit frontier #3): primes below `small_skip()` are not added to the byte
/// scores. They account for a large share of the score-write traffic (∑ 2·len/p) but contribute
/// little log weight, and they are still divided out during factoring, so skipping them only removes
/// sieve work. The score threshold is lowered by `small_slack()` to make up for their absent
/// contribution to a smooth `g(x)`. Both are read once per polynomial and cached locally.
/// `RUSQSIEVE_SMALL_SKIP` / `RUSQSIEVE_SMALL_SLACK` override them for tuning.
fn small_skip() -> u32 {
    static V: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_default("RUSQSIEVE_SMALL_SKIP", 100) as u32)
}
fn small_slack() -> usize {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_default("RUSQSIEVE_SMALL_SLACK", 8))
}
/// Extra score bits required above the smooth threshold. Raising the bar a few bits sharply cuts
/// false-positive survivors (≈99% of survivors are non-smooth) at the cost of a few more
/// polynomials. `RUSQSIEVE_THRESH_MARGIN` overrides it for tuning.
fn thresh_margin() -> i32 {
    static V: std::sync::OnceLock<i32> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_default("RUSQSIEVE_THRESH_MARGIN", 4) as i32)
}
/// Read an unsigned tuning override, defaulting when unset or non-Unix. Callers cache the result
/// in a per-knob `OnceLock` so the hot path never touches the environment or a lock.
fn env_default(name: &str, default: usize) -> usize {
    #[cfg(any(unix, windows))]
    {
        std::env::var(name)
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(default)
    }
    #[cfg(not(any(unix, windows)))]
    {
        let _ = name;
        default
    }
}

/// Add the two root strides for one factor-base prime. Interleaving the roots and unrolling two
/// hits at a time mirrors FLINT's flat sieve kernel and cuts loop-control overhead in the dominant
/// score-write pass. For the practical range where `g_bits <= 192`, scores cannot overflow: every
/// scored prime is at least 23, so the sum of rounded log weights is below
/// `g_bits * (1 + 1/log2(23))`, with ample room in a byte. Wider inputs retain saturating addition.
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
        let weight = (32 - p.leading_zeros()) as u8;
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

/// FLINT-style carried-position cache blocking. The factor base is replayed for each 256 KiB
/// block, but every prime resumes at its first unprocessed absolute hit, so no modular arithmetic
/// or re-striding setup is repeated. Score writes for large intervals remain local to one L2-sized
/// block.
#[allow(clippy::too_many_arguments)]
fn score_polynomial_blocked<const SATURATING: bool>(
    ctx: &Context,
    b: &Natural,
    bneg: bool,
    c: &Natural,
    csign: bool,
    root1: &[u32],
    root2: &[u32],
    scores: &mut [u8],
    small_skip: u32,
    pos1: &mut Vec<usize>,
    pos2: &mut Vec<usize>,
) {
    const BLOCK: usize = 4 * 65_536;
    let nfb = ctx.base.len();
    pos1.clear();
    pos1.resize(nfb, usize::MAX);
    pos2.clear();
    pos2.resize(nfb, usize::MAX);

    for (idx, e) in ctx.base.iter().enumerate() {
        let p = e.prime;
        if p == 2 || p < small_skip {
            continue;
        }
        if root1[idx] != u32::MAX {
            pos1[idx] = root1[idx] as usize;
            pos2[idx] = root2[idx] as usize;
        } else {
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
            pos1[idx] = add_mod_u32(xroot, ctx.interval_mod_p[idx], p) as usize;
        }
    }

    let mut block_end = BLOCK.min(scores.len());
    while block_end != 0 {
        for (idx, e) in ctx.base.iter().enumerate() {
            let p = e.prime;
            if p == 2 || p < small_skip {
                continue;
            }
            let step = p as usize;
            let weight = (32 - p.leading_zeros()) as u8;
            let mut a = pos1[idx];
            while a < block_end {
                scores[a] = if SATURATING {
                    scores[a].saturating_add(weight)
                } else {
                    scores[a].wrapping_add(weight)
                };
                a += step;
            }
            pos1[idx] = a;
            let mut b = pos2[idx];
            while b < block_end {
                scores[b] = if SATURATING {
                    scores[b].saturating_add(weight)
                } else {
                    scores[b].wrapping_add(weight)
                };
                b += step;
            }
            pos2[idx] = b;
        }
        if block_end == scores.len() {
            break;
        }
        block_end = (block_end + BLOCK).min(scores.len());
    }
}

/// Collect high-scoring positions. When the score array was initialized with `128 - threshold`,
/// the high bit is the threshold comparison; test eight bytes at once before examining a word's
/// individual positions. This is the portable form of FLINT's word-at-a-time candidate scan.
fn collect_candidates(
    scores: &[u8],
    threshold: u8,
    high_bit_biased: bool,
    candidates: &mut Vec<u32>,
) {
    candidates.clear();
    if !high_bit_biased {
        for (pos, &score) in scores.iter().enumerate() {
            if score >= threshold {
                candidates.push(pos as u32);
            }
        }
        return;
    }

    const HIGH_BITS: u64 = 0x8080_8080_8080_8080;
    let mut chunks = scores.chunks_exact(8);
    for (word_index, chunk) in chunks.by_ref().enumerate() {
        let word = u64::from_ne_bytes(chunk.try_into().unwrap());
        if word & HIGH_BITS == 0 {
            continue;
        }
        let base = word_index * 8;
        for (offset, &score) in chunk.iter().enumerate() {
            if score & 0x80 != 0 {
                candidates.push((base + offset) as u32);
            }
        }
    }
    let base = scores.len() - chunks.remainder().len();
    for (offset, &score) in chunks.remainder().iter().enumerate() {
        if score & 0x80 != 0 {
            candidates.push((base + offset) as u32);
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
    block_pos1: &mut Vec<usize>,
    block_pos2: &mut Vec<usize>,
    out: &mut Vec<Relation>,
) -> usize {
    let base = &ctx.base;
    let len = (ctx.interval as usize) * 2;
    // Tuning knobs read once per polynomial (cached), never in the per-prime hot loops.
    let small_skip = small_skip();
    let bb = b.checked_mul(b).unwrap();
    let (c, csign) = if bb >= ctx.sieve_n {
        (bb.wrapping_sub(&ctx.sieve_n).div_rem(a).unwrap().0, false)
    } else {
        (ctx.sieve_n.wrapping_sub(&bb).div_rem(a).unwrap().0, true)
    };
    let g_bits = ctx.sieve_n.bit_len().saturating_sub(a.bit_len());
    // Score threshold: a survivor's sieved-prime log-weight must come within `lp_allowance` bits
    // of g(x). SMALL_SLACK compensates for the tiny primes we no longer score; THRESH_MARGIN
    // raises the bar to suppress false-positive survivors. Bias practical-range score bytes so
    // the candidate comparison becomes a high-bit test, enabling a word-at-a-time scan.
    static THRESH_ADJ: std::sync::OnceLock<i32> = std::sync::OnceLock::new();
    let adj = *THRESH_ADJ.get_or_init(|| {
        std::env::var("RUSQSIEVE_THRESH_ADJ")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0)
    });
    let threshold =
        (g_bits as i32 - ctx.lp_allowance as i32 - small_slack() as i32 + thresh_margin() + adj)
            .clamp(1, 255) as u8;
    let high_bit_biased = g_bits <= 192 && threshold <= 128;
    let initial_score = if high_bit_biased { 128 - threshold } else { 0 };
    scores.clear();
    scores.resize(len, initial_score);
    // Logarithmic sieve using the self-initialized roots. Byte scores keep the
    // whole array resident in cache (SPEC §12.6). (Frontier #2 blocked/bucket sieving was
    // implemented and measured here: a naive per-block re-stride regressed 13–41% on this
    // large-L2 host because it fragments each prime's tight strided loop; see CLAUDE-AUDIT.md.)
    // The reference host has 1 MiB private L2. Keep arrays through 768 KiB on the faster flat
    // kernel; switch only once the score array itself reaches L2 capacity. Cache-constrained
    // targets can lower this in a future target-specific tuning table.
    const BLOCK_GATE: usize = 16 * 65_536;
    if scores.len() >= BLOCK_GATE {
        if g_bits <= 192 {
            score_polynomial_blocked::<false>(
                ctx, b, bneg, &c, csign, root1, root2, scores, small_skip, block_pos1, block_pos2,
            );
        } else {
            score_polynomial_blocked::<true>(
                ctx, b, bneg, &c, csign, root1, root2, scores, small_skip, block_pos1, block_pos2,
            );
        }
    } else if g_bits <= 192 {
        score_polynomial::<false>(ctx, b, bneg, &c, csign, root1, root2, scores, small_skip);
    } else {
        score_polynomial::<true>(ctx, b, bneg, &c, csign, root1, root2, scores, small_skip);
    }
    collect_candidates(scores, threshold, high_bit_biased, candidates);
    if candidates.is_empty() {
        return 0;
    }
    let survivors = candidates.len();

    let two_idx = base.iter().position(|e| e.prime == 2).map(|i| i as u32);
    // Factor-base entries with `prime < small_skip` occupy the low indices (the base is sorted
    // ascending). Those tiny primes are not sieved — gating them would waste a `fastmod` where a
    // direct divide is cheaper (they divide most survivors) — so they are divided out directly.
    let small_end = base.partition_point(|e| e.prime < small_skip);

    for &posu in candidates.iter() {
        let pos = posu as usize;
        // The score is the sum of one rounded log weight per sieve hit. Once confirmed factors
        // account for that weight, no later normal factor-base prime can divide this candidate.
        // This is FLINT's `extra_bits < sieve[i]` stopping rule and avoids scanning the tail of the
        // factor base for partial relations. It is exact on the non-saturating practical-range
        // kernel; wider saturated scores conservatively disable the shortcut.
        let score_target = if g_bits <= 192 {
            scores[pos].wrapping_sub(initial_score) as u16
        } else {
            u16::MAX
        };
        let mut confirmed_score = 0u16;
        let x = pos as i64 - ctx.interval as i64;
        let xabs = x.unsigned_abs();
        let ax = a.checked_mul(&Natural::from_u64(xabs)).unwrap();
        // t = a·x + b, needed for the relation's square root.
        let (t, tneg) = signed_add(&ax, x < 0, b, bneg);
        // Value to factor: g(x) = Q(x)/a = a·x² + 2b·x + c, computed directly with
        // signs (c_math = ∓c per csign). This avoids the wide t² squaring and the
        // division by a — a is guaranteed to divide Q since b² ≡ n (mod a).
        let ax2 = ax.checked_mul(&Natural::from_u64(xabs)).unwrap();
        let two_bx = b
            .wrapping_add(b)
            .checked_mul(&Natural::from_u64(xabs))
            .unwrap();
        let (gx, gxneg) = signed_add(&ax2, false, &two_bx, bneg ^ (x < 0));
        let (mut q, sign) = signed_add(&gx, gxneg, &c, csign);
        if q.is_zero() {
            continue;
        }
        let mut powers: Vec<(u32, u16)> = aidx.iter().copied().map(|i| (i, 1)).collect();
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
                record(ti, c2 as u16, &mut powers);
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
            let mut count = 0;
            while q.rem_u64(p) == 0 {
                q = q.div_rem_u64(p).unwrap().0;
                count += 1;
            }
            record(i as u32, count, &mut powers);
        }
        // Primes dividing `a` (seeded at exponent 1, root1 == MAX so not gated) — divide directly.
        for &ai in aidx {
            let p = base[ai as usize].prime as u64;
            let mut count = 0;
            while q.rem_u64(p) == 0 {
                q = q.div_rem_u64(p).unwrap().0;
                count += 1;
            }
            if count != 0 && p >= small_skip as u64 {
                confirmed_score += (64 - p.leading_zeros()) as u16;
            }
            record(ai, count, &mut powers);
        }
        // Barrett-gated trial division (FLINT `qsieve_evaluate_candidate` style): a normal prime `p`
        // divides `g(x)` iff the candidate's score-array position matches one of the translated
        // roots modulo `p`. Compute that residue with a precomputed multiply-shift (`fastmod`) — no
        // hardware divide — and bignum-divide only on a hit. This makes trial division O(nfb) cheap
        // tests + O(#factors) divides, with no second sieve pass, so it stays cheap at every
        // factor-base size (unlike the old O(nfb) bignum-remainder loop).
        for idx in small_end..base.len() {
            if q.is_one() || confirmed_score >= score_target {
                break;
            }
            let r1 = root1[idx];
            if r1 == u32::MAX {
                continue; // prime divides `a`; handled above
            }
            let p = base[idx].prime;
            let posmodp = fastmod(posu, p, ctx.pinv[idx]);
            debug_assert_eq!(posmodp, posu % p, "fastmod mismatch p={p}");
            if posmodp != r1 && posmodp != root2[idx] {
                continue;
            }
            let pu = p as u64;
            let mut count = 0;
            while q.rem_u64(pu) == 0 {
                q = q.div_rem_u64(pu).unwrap().0;
                count += 1;
            }
            if count != 0 {
                confirmed_score += (32 - p.leading_zeros()) as u16;
            }
            record(idx as u32, count, &mut powers);
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
            powers,
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
fn combine_cycle(rels: &[Relation], n: &Natural) -> Column {
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
    if is_prime64(q) {
        return (q <= single_limit).then_some(LargePrime::One(q));
    }
    if !double_enabled {
        return None;
    }
    let d = pollard_u64(q)?;
    let e = q / d;
    if d > 1 && e > 1 && d <= single_limit && e <= single_limit && is_prime64(d) && is_prime64(e) {
        Some(LargePrime::Two(d.min(e), d.max(e)))
    } else {
        None
    }
}

/// Pollard's rho (Floyd) for a small composite `u64`; returns a nontrivial factor.
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
    let mut c = 1u64;
    while c < 64 {
        let f = |v: u64| ((v as u128 * v as u128 + c as u128) % n as u128) as u64;
        let (mut x, mut y, mut d) = (2u64, 2u64, 1u64);
        while d == 1 {
            x = f(x);
            y = f(f(y));
            d = gcd(x.abs_diff(y), n);
        }
        if d != n {
            return Some(d);
        }
        c += 1;
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
    edge: Vec<Option<Relation>>,
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
    fn path(&self, mut v: u32, out: &mut Vec<Relation>) {
        while self.parent[v as usize] != v {
            out.push(self.edge[v as usize].clone().unwrap());
            v = self.parent[v as usize];
        }
    }
    /// Re-root the tree containing `v` so that `v` becomes its root.
    fn reroot(&mut self, v: u32) {
        let mut chain = vec![v];
        let mut edges: Vec<Relation> = Vec::new();
        let mut c = v;
        while self.parent[c as usize] != c {
            edges.push(self.edge[c as usize].clone().unwrap());
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
        self.parent[b as usize] = a;
        self.edge[b as usize] = Some(rel);
    }
}

/// Deterministically accumulates relations into matrix columns, matching partial
/// relations through the large-prime graph.
struct RelationCollector {
    forest: Forest,
    columns: Vec<Column>,
}
impl RelationCollector {
    fn new() -> Self {
        Self {
            forest: Forest::default(),
            columns: Vec::new(),
        }
    }
    fn ingest(&mut self, rel: Relation, n: &Natural) {
        match rel.large {
            LargePrime::None => self.columns.push(to_column(rel)),
            LargePrime::One(p) => self.edge(p, 1, rel, n),
            LargePrime::Two(a, b) if a == b => self.columns.push(combine_cycle(&[rel], n)),
            LargePrime::Two(a, b) => self.edge(a, b, rel, n),
        }
    }
    fn edge(&mut self, pa: u64, pb: u64, rel: Relation, n: &Natural) {
        let va = self.forest.vertex(pa);
        let vb = self.forest.vertex(pb);
        if self.forest.root(va) == self.forest.root(vb) {
            let mut cyc = vec![rel];
            self.forest.path(va, &mut cyc);
            self.forest.path(vb, &mut cyc);
            self.columns.push(combine_cycle(&cyc, n));
        } else {
            self.forest.link(va, vb, rel);
        }
    }
}
fn inv_u32(a: u32, p: u32) -> Option<u32> {
    if a == 0 {
        return None;
    }
    let (mut t, mut nt) = (0i64, 1i64);
    let (mut r, mut nr) = (p as i64, a as i64);
    while nr != 0 {
        let q = r / nr;
        (t, nt) = (nt, t - q * nt);
        (r, nr) = (nr, r - q * nr)
    }
    (r == 1).then_some(t.rem_euclid(p as i64) as u32)
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
fn xorshift(mut x: u64) -> u64 {
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    x
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
        if is_prime64(p) {
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
fn is_prime64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1
    }
    for a in [2u64, 325, 9375, 28178, 450775, 9780504, 1795265022] {
        if a % n == 0 {
            continue;
        }
        let mut x = powmod64(a % n, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut ok = false;
        for _ in 1..s {
            x = (x as u128 * x as u128 % n as u128) as u64;
            if x == n - 1 {
                ok = true;
                break;
            }
        }
        if !ok {
            return false;
        }
    }
    true
}
fn powmod64(mut a: u64, mut e: u64, n: u64) -> u64 {
    let mut r = 1;
    while e != 0 {
        if e & 1 != 0 {
            r = (r as u128 * a as u128 % n as u128) as u64
        }
        a = (a as u128 * a as u128 % n as u128) as u64;
        e >>= 1
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let context = prepare(p.checked_mul(&q).unwrap()).unwrap();
        let a = execute(&context, EngineJob { family: 7 });
        let b = execute(&context, EngineJob { family: 7 });
        assert_eq!(a.family, b.family);
        assert_eq!(a.polynomials, b.polynomials);
        assert_eq!(a.relations, b.relations);
        assert!(a.polynomials > 0);
    }

    #[test]
    fn collector_accepts_out_of_order_results() {
        let p = Natural::from_u64(18_446_744_073_709_551_557);
        let q = Natural::from_u64(18_446_744_073_709_551_533);
        let context = prepare(p.checked_mul(&q).unwrap()).unwrap();
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
        let factors = factor(n.clone(), 2, |_| true).unwrap();
        assert_eq!(factors, [q, p]);
        assert_eq!(
            factors
                .iter()
                .try_fold(Natural::ONE, |a, b| a.checked_mul(b)),
            Some(n)
        );
    }
}
