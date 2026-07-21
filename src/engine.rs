//! Portable SIQS engine and scheduler-facing work kernels.
use crate::{Natural, PARTS, jacobi_u64};
use crate::f2::SparseBinaryMatrix;
use crate::qs::{AutoOr, FactorBaseEntry, MultiplierChoice, QsConfig, prepare_siqs};
#[cfg(any(unix, windows))]
use crate::{PrimalityConfig, is_probable_prime};
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::sync::Arc;
#[cfg(any(unix, windows))]
use std::sync::{Mutex, mpsc};

#[derive(Clone, Copy, Debug)]
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
    interval: i32,
    target_a_bits: usize,
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
    /// The two sieve roots per factor-base prime for the current polynomial.
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
    /// Resieve (audit frontier #1): `cand_at[pos]` is the index into `candidates`
    /// of the survivor at `pos`, or `u32::MAX`. Kept clean between polynomials by
    /// clearing only survivor slots, so it never needs a full memset.
    cand_at: Vec<u32>,
    /// Resieve: per-survivor list of factor-base indices whose sieve stride lands
    /// on it. Only these primes are trial-divided out of the survivor's `g(x)`.
    resieve_fac: Vec<Vec<u32>>,
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
    let sieve_n = n.checked_mul(&Natural::from_u64(k)).unwrap_or_else(|| n.clone());
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(p.factor_base_bound),
        multiplier: MultiplierChoice::Value(k as u32),
        ..QsConfig::default()
    };
    let prepared = prepare_siqs(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
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
        interval: p.sieve_half_width as i32,
        target_a_bits: target_a.bit_len(),
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
        let target = context.0.base.len() + 64;
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
    if columns.len() <= ctx.base.len() {
        return Err(EngineError::InsufficientRelations);
    }
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

#[cfg(any(unix, windows))]
pub fn factor(
    mut n: Natural,
    threads: usize,
    mut progress: impl FnMut(EngineProgress),
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
    progress: &mut impl FnMut(EngineProgress),
    out: &mut Vec<Natural>,
) -> Result<(), EngineError> {
    progress(EngineProgress {
        phase: EnginePhase::Preprocessing,
        polynomials: 0,
        relations: 0,
        target: 0,
        workers: threads,
    });
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
    progress: &mut impl FnMut(EngineProgress),
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
    progress(EngineProgress {
        phase: EnginePhase::BuildingFactorBase,
        polynomials: 0,
        relations: 0,
        target: 0,
        workers: threads,
    });
    let p = crate::qs::parameters::engine_params(n.bit_len());
    let bound = p.factor_base_bound;
    let interval = p.sieve_half_width as i32;
    let prof = std::env::var_os("RUSQSIEVE_PROFILE").is_some();
    let t_fb = std::time::Instant::now();
    let k = knuth_schroeppel(&n);
    let sieve_n = n.checked_mul(&Natural::from_u64(k)).unwrap_or_else(|| n.clone());
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(bound),
        multiplier: MultiplierChoice::Value(k as u32),
        ..QsConfig::default()
    };
    let prepared = prepare_siqs(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
    let target = base.len() + 64;
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
        interval,
        target_a_bits: target_a.bit_len(),
        lp_allowance: p.lp_allowance,
        single_limit,
        double_enabled,
    });
    let (job_tx, job_rx) = mpsc::channel::<Option<u64>>();
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (res_tx, res_rx) = mpsc::channel();
    let mut handles = Vec::new();
    for _ in 0..threads {
        let rx = job_rx.clone();
        let tx = res_tx.clone();
        let c = ctx.clone();
        handles.push(std::thread::spawn(move || {
            let mut scratch = EngineScratch::default();
            loop {
                let job = rx.lock().unwrap().recv();
                match job {
                    Ok(Some(f)) => {
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
    let max_families = 100_000u64;
    while collector.columns.len() < target && next_merge < max_families {
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
            progress(EngineProgress {
                phase: EnginePhase::Sieving,
                polynomials,
                relations: collector.columns.len(),
                target,
                workers: threads,
            });
        }
        while outstanding < threads * 2 && next_send < max_families && collector.columns.len() < target
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
    if collector.columns.len() <= base.len() {
        return Err(EngineError::InsufficientRelations);
    }
    progress(EngineProgress {
        phase: EnginePhase::LinearAlgebra,
        polynomials,
        relations: collector.columns.len(),
        target,
        workers: threads,
    });
    let t_la = std::time::Instant::now();
    let result = extract(&ctx, &collector.columns);
    if prof {
        eprintln!(
            "PROFILE extract(LA)={:.3}s columns={}",
            t_la.elapsed().as_secs_f64(),
            collector.columns.len()
        );
    }
    progress(EngineProgress {
        phase: EnginePhase::Extracting,
        polynomials,
        relations: collector.columns.len(),
        target,
        workers: threads,
    });
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

    // SIQS B-values: b = Σ ±Bⱼ (mod a), with Bⱼ ≡ sqrt(n) (mod qⱼ), 0 (mod other q).
    let mut bvals: Vec<Natural> = Vec::with_capacity(s);
    for &i in &aidx {
        let q = base[i as usize].prime;
        let Some((ap, _)) = a.div_rem_u64(q as u64) else {
            return empty(family);
        };
        let Some(apinv) = inv_u32(ap.mod_u64(q as u64) as u32, q) else {
            return empty(family);
        };
        let coeff = (base[i as usize].sqrt_n as u64 * apinv as u64) % q as u64;
        bvals.push(ap.mul_mod(&Natural::from_u64(coeff), &a));
    }
    let mut b = Natural::ZERO;
    for bj in &bvals {
        b = b.add_mod(bj, &a);
    }
    // True (unreduced) 2·Bⱼ, each < 2a. Kept unreduced so the O(1) root advance can
    // account for the mod-a wrap uniformly.
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
        let bp = b.mod_u64(p as u64) as u32;
        scratch.root1[idx] = mulmod_u32((e.sqrt_n + p - bp) % p, ainvp, p);
        scratch.root2[idx] = mulmod_u32(((p - e.sqrt_n) % p + p - bp) % p, ainvp, p);
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
            &aidx,
            &scratch.root1,
            &scratch.root2,
            &mut scratch.scores,
            &mut scratch.candidates,
            &mut scratch.cand_at,
            &mut scratch.resieve_fac,
            &mut relations,
        ) as u64;
        if v + 1 >= variants {
            break;
        }
        let j = (v + 1).trailing_zeros() as usize;
        let gray = v ^ (v >> 1);
        let flip_to_one = (gray >> j) & 1 == 0;
        // Advance b to the next polynomial (kept reduced in [0, a)) and record the
        // number of a-wraps: because a·a⁻¹ ≡ 1 (mod p), each wrap shifts every
        // prime's root by the same amount, so `shift` is applied uniformly below.
        let (add_bainv, shift): (bool, i64) = if flip_to_one {
            // b_new = (b - 2Bⱼ) mod a; raw = b + 2a - 2Bⱼ ∈ (0, 3a).
            let mut raw = b.wrapping_add(&a).wrapping_add(&a).wrapping_sub(&two_full[j]);
            let mut kp = 0i64;
            while raw >= a {
                raw = raw.wrapping_sub(&a);
                kp += 1;
            }
            b = raw;
            (true, -(2 - kp))
        } else {
            let mut raw = b.wrapping_add(&two_full[j]);
            let mut k = 0i64;
            while raw >= a {
                raw = raw.wrapping_sub(&a);
                k += 1;
            }
            b = raw;
            (false, k)
        };
        let off = j * nfb;
        for idx in 0..nfb {
            if scratch.root1[idx] == u32::MAX {
                continue;
            }
            let p = base[idx].prime as i64;
            let d = scratch.bainv[off + idx] as i64;
            let delta = (if add_bainv { d } else { -d } + shift).rem_euclid(p);
            scratch.root1[idx] = ((scratch.root1[idx] as i64 + delta) % p) as u32;
            scratch.root2[idx] = ((scratch.root2[idx] as i64 + delta) % p) as u32;
        }
    }
    FamilyResult {
        family,
        polynomials: variants,
        relations,
        survivors,
    }
}

fn choose_a(ctx: &Context, family: u64) -> Option<(Natural, Vec<u32>)> {
    let pool: Vec<usize> = ctx
        .base
        .iter()
        .enumerate()
        .filter(|(_, e)| e.prime > 1000)
        .map(|(i, _)| i)
        .collect();
    if pool.len() < 8 {
        return None;
    }
    let mut state = family ^ 0x9e3779b97f4a7c15;
    let mut a = Natural::ONE;
    let mut idx = Vec::new();
    while a.bit_len() < ctx.target_a_bits && idx.len() < 12 {
        state = xorshift(state);
        let i = pool[state as usize % pool.len()];
        if idx.contains(&(i as u32)) {
            continue;
        }
        let next = a.checked_mul(&Natural::from_u64(ctx.base[i].prime as u64))?;
        a = next;
        idx.push(i as u32)
    }
    (idx.len() >= 3).then_some((a, idx))
}

/// Minimum factor-base size at which resieving (audit frontier #1) beats full trial division.
/// Below this the extra strided resieve pass costs more than the divisions it removes; above it
/// the trial-division cost dominates and resieving wins. Measured crossover on the reference host
/// (192/224/256-bit → nfb 1550/3027/12904): loss at ~3k, win at ~13k.
const RESIEVE_MIN_FB: usize = 7000;

/// Tiny-prime skipping (audit frontier #3): primes below this are not added to the byte scores.
/// They account for a large share of the score-write traffic (∑ 2·len/p) but contribute little
/// log weight, and they are still divided out during factoring, so skipping them only removes
/// sieve work. The score threshold is lowered by `SMALL_SLACK` to make up for their absent
/// contribution to a smooth `g(x)`.
const SMALL_SKIP: u32 = 20;
const SMALL_SLACK: usize = 3;
/// Extra score bits required above the smooth threshold. Tuning found that raising the bar a few
/// bits sharply cuts false-positive survivors (99% of survivors were non-smooth) at the cost of a
/// few more polynomials — a net win across 192/224/256-bit.
const THRESH_MARGIN: i32 = 4;

#[allow(clippy::too_many_arguments)]
fn sieve_one_poly(
    ctx: &Context,
    a: &Natural,
    b: &Natural,
    aidx: &[u32],
    root1: &[u32],
    root2: &[u32],
    scores: &mut Vec<u8>,
    candidates: &mut Vec<u32>,
    cand_at: &mut Vec<u32>,
    resieve_fac: &mut Vec<Vec<u32>>,
    out: &mut Vec<Relation>,
) -> usize {
    let base = &ctx.base;
    let len = (ctx.interval as usize) * 2;
    scores.clear();
    scores.resize(len, 0);
    let bb = b.checked_mul(b).unwrap();
    let (c, csign) = if bb >= ctx.sieve_n {
        (bb.wrapping_sub(&ctx.sieve_n).div_rem(a).unwrap().0, false)
    } else {
        (ctx.sieve_n.wrapping_sub(&bb).div_rem(a).unwrap().0, true)
    };
    let neg_interval = -(ctx.interval as i64);
    // Logarithmic sieve using the self-initialized roots. Byte scores keep the
    // whole array resident in cache (SPEC §12.6). (Frontier #2 blocked/bucket sieving was
    // implemented and measured here: a naive per-block re-stride regressed 13–41% on this
    // large-L2 host because it fragments each prime's tight strided loop; see CLAUDE-AUDIT.md.)
    for (idx, e) in base.iter().enumerate() {
        let p = e.prime;
        // Skip 2 (special) and tiny primes (frontier #3): both are recovered during factoring.
        if p == 2 || p < SMALL_SKIP {
            continue;
        }
        let pu = p as usize;
        let weight = (32 - p.leading_zeros()) as u8;
        if root1[idx] != u32::MAX {
            for &root in &[root1[idx], root2[idx]] {
                let start = (root as i64 - neg_interval.rem_euclid(p as i64)).rem_euclid(p as i64)
                    as usize;
                let mut pos = start;
                while pos < len {
                    scores[pos] = scores[pos].saturating_add(weight);
                    pos += pu;
                }
            }
        } else {
            // p | a: the polynomial is linear (2bx + c) mod p — one root, per poly.
            let bp = b.mod_u64(p as u64) as u32;
            let denom = (2 * bp as u64 % p as u64) as u32;
            let Some(inv) = inv_u32(denom, p) else {
                continue;
            };
            let cm = c.mod_u64(p as u64) as u32;
            let signed_c = if csign && cm != 0 { p - cm } else { cm };
            let root = mulmod_u32(if signed_c == 0 { 0 } else { p - signed_c }, inv, p);
            let start =
                (root as i64 - neg_interval.rem_euclid(p as i64)).rem_euclid(p as i64) as usize;
            let mut pos = start;
            while pos < len {
                scores[pos] = scores[pos].saturating_add(weight);
                pos += pu;
            }
        }
    }
    let g_bits = ctx.sieve_n.bit_len().saturating_sub(a.bit_len());
    // Score threshold: a survivor's sieved-prime log-weight must come within `lp_allowance` bits
    // of g(x). SMALL_SLACK compensates for the tiny primes we no longer score; THRESH_MARGIN
    // raises the bar to suppress false-positive survivors (measured to cut wasted trial division
    // for a small polynomial-count increase). RUSQSIEVE_THRESH_ADJ (read once) tunes it further.
    static THRESH_ADJ: std::sync::OnceLock<i32> = std::sync::OnceLock::new();
    let adj = *THRESH_ADJ.get_or_init(|| {
        std::env::var("RUSQSIEVE_THRESH_ADJ")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0)
    });
    let threshold = (g_bits as i32 - ctx.lp_allowance as i32 - SMALL_SLACK as i32 + THRESH_MARGIN
        + adj)
        .clamp(1, 255) as u8;
    candidates.clear();
    for (pos, &score) in scores.iter().enumerate() {
        if score >= threshold {
            candidates.push(pos as u32);
        }
    }
    if candidates.is_empty() {
        return 0;
    }
    let survivors = candidates.len();

    // Resieve (audit frontier #1) pays off only once the factor base is large enough that
    // trial-dividing every survivor by all of it outweighs one extra strided pass. Below the
    // gate we keep the original full trial division (see `RESIEVE_MIN_FB`).
    let use_resieve = base.len() >= RESIEVE_MIN_FB;
    let two_idx = base.iter().position(|e| e.prime == 2).map(|i| i as u32);
    if use_resieve {
        // Map position → survivor index (u32::MAX = none). Only survivor slots are ever written,
        // so `cand_at` stays clean between polynomials without a full memset.
        if cand_at.len() != len {
            cand_at.clear();
            cand_at.resize(len, u32::MAX);
        }
        for (si, &posu) in candidates.iter().enumerate() {
            cand_at[posu as usize] = si as u32;
        }
        if resieve_fac.len() < candidates.len() {
            resieve_fac.resize_with(candidates.len(), Vec::new);
        }
        for f in resieve_fac[..candidates.len()].iter_mut() {
            f.clear();
        }
        // Record, per survivor, which normal primes' strides land on it. Because
        // `p | g(x) ⟺ x ≡ root (mod p)`, this misses no normal prime. Prime 2 and the primes
        // dividing `a` (root1 == MAX) are divided out directly per survivor below — both are
        // cheap and sidestep the linear-root edge cases.
        for (idx, e) in base.iter().enumerate() {
            let p = e.prime;
            if p == 2 || root1[idx] == u32::MAX {
                continue;
            }
            let pu = p as usize;
            for &root in &[root1[idx], root2[idx]] {
                let start = (root as i64 - neg_interval.rem_euclid(p as i64)).rem_euclid(p as i64)
                    as usize;
                let mut pos = start;
                while pos < len {
                    let si = cand_at[pos];
                    if si != u32::MAX {
                        resieve_fac[si as usize].push(idx as u32);
                    }
                    pos += pu;
                }
            }
        }
    }

    for (si, &posu) in candidates.iter().enumerate() {
        let pos = posu as usize;
        if use_resieve {
            cand_at[pos] = u32::MAX; // restore the clean state for the next polynomial
        }
        let x = pos as i64 - ctx.interval as i64;
        let xabs = x.unsigned_abs();
        let ax = a.checked_mul(&Natural::from_u64(xabs)).unwrap();
        // t = a·x + b, needed for the relation's square root.
        let (t, tneg) = if x >= 0 {
            (ax.checked_add(b).unwrap(), false)
        } else if ax >= *b {
            (ax.wrapping_sub(b), true)
        } else {
            (b.wrapping_sub(&ax), false)
        };
        // Value to factor: g(x) = Q(x)/a = a·x² + 2b·x + c, computed directly with
        // signs (c_math = ∓c per csign). This avoids the wide t² squaring and the
        // division by a — a is guaranteed to divide Q since b² ≡ n (mod a).
        let ax2 = ax.checked_mul(&Natural::from_u64(xabs)).unwrap();
        let two_bx = b.wrapping_add(b).checked_mul(&Natural::from_u64(xabs)).unwrap();
        let mut pos_sum = ax2;
        let mut neg_sum = Natural::ZERO;
        if x >= 0 {
            pos_sum = pos_sum.checked_add(&two_bx).unwrap();
        } else {
            neg_sum = two_bx;
        }
        if csign {
            neg_sum = neg_sum.checked_add(&c).unwrap();
        } else {
            pos_sum = pos_sum.checked_add(&c).unwrap();
        }
        let (mut q, sign) = if pos_sum >= neg_sum {
            (pos_sum.wrapping_sub(&neg_sum), false)
        } else {
            (neg_sum.wrapping_sub(&pos_sum), true)
        };
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
        if use_resieve {
            // Prime 2 (not sieved): strip via trailing zeros.
            if let Some(ti) = two_idx {
                let c2 = q.trailing_zeros();
                if c2 != 0 {
                    q >>= c2;
                    record(ti, c2 as u16, &mut powers);
                }
            }
            // Primes dividing `a` (seeded at exponent 1) — a handful, divide directly.
            for &ai in aidx {
                let p = base[ai as usize].prime as u64;
                let mut count = 0;
                while q.rem_u64(p) == 0 {
                    q = q.div_rem_u64(p).unwrap().0;
                    count += 1;
                }
                record(ai, count, &mut powers);
            }
            // Exactly the normal primes the resieve pass landed on this survivor.
            for &idx in resieve_fac[si].iter() {
                let p = base[idx as usize].prime as u64;
                let mut count = 0;
                while q.rem_u64(p) == 0 {
                    q = q.div_rem_u64(p).unwrap().0;
                    count += 1;
                }
                record(idx, count, &mut powers);
            }
        } else {
            // Small factor base: full trial division. Root-gating (a cheap `pos ≡ root (mod p)`
            // pre-test) was measured ~2–4% slower here — the early `q.is_one()` break plus `q`
            // shrinking to a single limb already make these divisions cheaper than the gate's
            // per-prime modulo — so we keep the straightforward loop.
            for (i, e) in base.iter().enumerate() {
                if q.is_one() {
                    break;
                }
                let p = e.prime as u64;
                let mut count = 0;
                while q.rem_u64(p) == 0 {
                    q = q.div_rem_u64(p).unwrap().0;
                    count += 1;
                }
                record(i as u32, count, &mut powers);
            }
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
    if d > 1
        && e > 1
        && d <= single_limit
        && e <= single_limit
        && is_prime64(d)
        && is_prime64(e)
    {
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
        let factors = factor(n.clone(), 2, |_| {}).unwrap();
        assert_eq!(factors, [q, p]);
        assert_eq!(
            factors
                .iter()
                .try_fold(Natural::ONE, |a, b| a.checked_mul(b)),
            Some(n)
        );
    }
}
