//! Portable SIQS engine and scheduler-facing work kernels.
use crate::Natural;
use crate::f2::SparseBinaryMatrix;
use crate::qs::{AutoOr, FactorBaseEntry, QsConfig, prepare_siqs};
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
    n: Natural<16>,
    base: Arc<[FactorBaseEntry]>,
    interval: i32,
    target_a_bits: usize,
    large_limit: u64,
}

#[derive(Clone)]
struct Relation {
    root: Natural<16>,
    sign: bool,
    powers: Vec<(u32, u16)>,
    large: Option<u64>,
}

#[derive(Clone)]
struct Column {
    root: Natural<16>,
    sign: bool,
    powers: Vec<(u32, u32)>,
    extra_sqrt: u64,
}

struct FamilyResult {
    family: u64,
    polynomials: u64,
    relations: Vec<Relation>,
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

/// Prepare an immutable context without creating threads.
pub fn prepare(n: Natural<16>) -> Result<EngineContext, EngineError> {
    let bits = n.bit_len();
    let bound = match bits {
        0..=128 => 4_000,
        129..=180 => 20_000,
        181..=220 => 60_000,
        221..=270 => 140_000,
        271..=330 => 300_000,
        _ => 600_000,
    };
    let interval = match bits {
        0..=180 => 32_768,
        181..=230 => 65_536,
        _ => 131_072,
    };
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(bound),
        ..QsConfig::default()
    };
    let prepared = prepare_siqs(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
    let target_a = n.floor_sqrt().div_rem_u64(interval as u64).unwrap().0;
    Ok(EngineContext(Arc::new(Context {
        n,
        base,
        interval,
        target_a_bits: target_a.bit_len(),
        large_limit: (bound as u64).saturating_mul(bound as u64),
    })))
}

/// Execute a job using only the caller's thread and owned scratch memory.
pub fn execute(context: &EngineContext, job: EngineJob) -> EngineJobResult {
    let inner = sieve_family(&context.0, job.family);
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
    pending: HashMap<u64, Relation>,
    columns: Vec<Column>,
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
            pending: HashMap::new(),
            columns: Vec::new(),
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
        while let Some(r) = self.buffered.remove(&self.next_merge) {
            self.next_merge += 1;
            self.polynomials += r.polynomials;
            for rel in r.relations {
                if let Some(lp) = rel.large {
                    if let Some(first) = self.pending.remove(&lp) {
                        self.columns
                            .push(combine(first, rel, lp, &self.context.0.n))
                    } else {
                        self.pending.insert(lp, rel);
                    }
                } else {
                    self.columns.push(to_column(rel))
                }
            }
        }
    }
    pub fn is_ready(&self) -> bool {
        self.columns.len() >= self.target
    }
    pub fn relations(&self) -> usize {
        self.columns.len()
    }
    pub fn target(&self) -> usize {
        self.target
    }
    pub fn polynomials(&self) -> u64 {
        self.polynomials
    }
    pub fn extract_factor(&self) -> Result<Natural<16>, EngineError> {
        extract(&self.context.0, &self.columns)
    }
}

fn extract(ctx: &Context, columns: &[Column]) -> Result<Natural<16>, EngineError> {
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
    for dep in matrix.dense_dependencies().iter() {
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
            y = y.mul_mod(&Natural::from_u64(c.extra_sqrt), &ctx.n);
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
    mut n: Natural<16>,
    threads: usize,
    mut progress: impl FnMut(EngineProgress),
) -> Result<Vec<Natural<16>>, EngineError> {
    if n.is_zero() {
        return Err(EngineError::Setup("zero has no prime factorization".into()));
    }
    let primality = PrimalityConfig::default();
    let mut factors = Vec::new();
    for p in small_primes(10_000) {
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
    n: Natural<16>,
    threads: usize,
    pc: &PrimalityConfig,
    progress: &mut impl FnMut(EngineProgress),
    out: &mut Vec<Natural<16>>,
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
    if is_probable_prime(&n, pc) {
        out.push(n);
        return Ok(());
    }
    if n.bit_len() < 120 {
        let factors = crate::factor::factor_complete(n, &crate::FactorConfig::default())
            .map_err(|error| EngineError::Setup(error.to_string()))?;
        for (prime, exponent) in factors.iter() {
            for _ in 0..exponent.get() {
                out.push(prime.clone());
            }
        }
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
    n: Natural<16>,
    threads: usize,
    progress: &mut impl FnMut(EngineProgress),
) -> Result<Natural<16>, EngineError> {
    progress(EngineProgress {
        phase: EnginePhase::BuildingFactorBase,
        polynomials: 0,
        relations: 0,
        target: 0,
        workers: threads,
    });
    let bits = n.bit_len();
    let bound = match bits {
        0..=128 => 4_000,
        129..=180 => 20_000,
        181..=220 => 60_000,
        221..=270 => 140_000,
        271..=330 => 300_000,
        _ => 600_000,
    };
    let interval = match bits {
        0..=180 => 32_768,
        181..=230 => 65_536,
        _ => 131_072,
    };
    let qcfg = QsConfig {
        factor_base_bound: AutoOr::Value(bound),
        ..QsConfig::default()
    };
    let prepared = prepare_siqs(&n, &qcfg).map_err(|e| EngineError::Setup(e.to_string()))?;
    let base: Arc<[FactorBaseEntry]> = prepared.factor_base().entries().to_vec().into();
    let target = base.len() + 64;
    let target_a = n.floor_sqrt().div_rem_u64(interval as u64).unwrap().0;
    let ctx = Arc::new(Context {
        n: n.clone(),
        base: base.clone(),
        interval,
        target_a_bits: target_a.bit_len(),
        large_limit: (bound as u64).saturating_mul(bound as u64),
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
            loop {
                let job = rx.lock().unwrap().recv();
                match job {
                    Ok(Some(f)) => {
                        if tx.send(sieve_family(&c, f)).is_err() {
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
    let mut buffered = BTreeMap::new();
    let mut pending = HashMap::<u64, Relation>::new();
    let mut columns = Vec::new();
    let mut polynomials = 0u64;
    let max_families = 100_000u64;
    while columns.len() < target && next_merge < max_families {
        let result = res_rx.recv().map_err(|_| EngineError::Worker)?;
        outstanding -= 1;
        buffered.insert(result.family, result);
        while let Some(r) = buffered.remove(&next_merge) {
            next_merge += 1;
            polynomials += r.polynomials;
            for rel in r.relations {
                if let Some(lp) = rel.large {
                    if let Some(first) = pending.remove(&lp) {
                        columns.push(combine(first, rel, lp, &n))
                    } else {
                        pending.insert(lp, rel);
                    }
                } else {
                    columns.push(to_column(rel))
                }
                if columns.len() >= target {
                    break;
                }
            }
            progress(EngineProgress {
                phase: EnginePhase::Sieving,
                polynomials,
                relations: columns.len(),
                target,
                workers: threads,
            });
        }
        while outstanding < threads * 2 && next_send < max_families && columns.len() < target {
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
    if columns.len() <= base.len() {
        return Err(EngineError::InsufficientRelations);
    }
    progress(EngineProgress {
        phase: EnginePhase::LinearAlgebra,
        polynomials,
        relations: columns.len(),
        target,
        workers: threads,
    });
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
    let matrix = SparseBinaryMatrix::from_columns(base.len() + 1, &matrix_cols)
        .map_err(|_| EngineError::InvalidDependency)?;
    let deps = matrix.dense_dependencies();
    progress(EngineProgress {
        phase: EnginePhase::Extracting,
        polynomials,
        relations: columns.len(),
        target,
        workers: threads,
    });
    for dep in deps.iter() {
        if !matrix.verify_dependency(dep) {
            return Err(EngineError::InvalidDependency);
        }
        let mut x = Natural::ONE;
        let mut sums = vec![0u32; base.len()];
        let mut y = Natural::ONE;
        for (j, c) in columns.iter().enumerate() {
            if (dep[j / 64] >> (j % 64)) & 1 == 0 {
                continue;
            }
            x = x.mul_mod(&c.root, &n);
            y = y.mul_mod(&Natural::from_u64(c.extra_sqrt), &n);
            for &(i, e) in &c.powers {
                sums[i as usize] += e
            }
        }
        for (e, &s) in base.iter().zip(&sums) {
            for _ in 0..s / 2 {
                y = y.mul_mod(&Natural::from_u64(e.prime as u64), &n)
            }
        }
        let delta = if x >= y {
            x.wrapping_sub(&y)
        } else {
            y.wrapping_sub(&x)
        };
        let g = delta.gcd(&n);
        if !g.is_one() && g != n {
            return Ok(g);
        }
        let g = x.add_mod(&y, &n).gcd(&n);
        if !g.is_one() && g != n {
            return Ok(g);
        }
    }
    Err(EngineError::NoFactor)
}

fn sieve_family(ctx: &Context, family: u64) -> FamilyResult {
    let Some((a, aidx)) = choose_a(ctx, family) else {
        return FamilyResult {
            family,
            polynomials: 0,
            relations: Vec::new(),
        };
    };
    let variants = 1u64 << aidx.len().saturating_sub(1).min(6);
    let mut relations = Vec::new();
    for variant in 0..variants {
        let Some(b) = crt_root(ctx, &a, &aidx, variant) else {
            continue;
        };
        sieve_polynomial(ctx, &a, &b, &aidx, family, variant, &mut relations)
    }
    FamilyResult {
        family,
        polynomials: variants,
        relations,
    }
}

fn choose_a(ctx: &Context, family: u64) -> Option<(Natural<16>, Vec<u32>)> {
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

fn crt_root(ctx: &Context, a: &Natural<16>, idx: &[u32], variant: u64) -> Option<Natural<16>> {
    let mut b = Natural::ZERO;
    for (j, &i) in idx.iter().enumerate() {
        let e = ctx.base[i as usize];
        let ap = a.div_rem_u64(e.prime as u64)?.0;
        let inv = inv_u32(ap.mod_u64(e.prime as u64) as u32, e.prime)?;
        let root = if j + 1 == idx.len() || variant >> (j.min(63)) & 1 == 0 {
            e.sqrt_n
        } else {
            if e.sqrt_n == 0 { 0 } else { e.prime - e.sqrt_n }
        };
        let coeff = (root as u64 * inv as u64) % e.prime as u64;
        let term = ap.mul_mod(&Natural::from_u64(coeff), a);
        b = b.add_mod(&term, a)
    }
    Some(b)
}

fn sieve_polynomial(
    ctx: &Context,
    a: &Natural<16>,
    b: &Natural<16>,
    aidx: &[u32],
    _family: u64,
    _variant: u64,
    out: &mut Vec<Relation>,
) {
    let len = (ctx.interval as usize) * 2;
    let mut scores = vec![0u16; len];
    let bb = b.checked_mul(b).unwrap();
    let (c, csign) = if bb >= ctx.n {
        (bb.wrapping_sub(&ctx.n).div_rem(a).unwrap().0, false)
    } else {
        (ctx.n.wrapping_sub(&bb).div_rem(a).unwrap().0, true)
    };
    for (index, e) in ctx.base.iter().enumerate() {
        let p = e.prime;
        if p == 2 {
            continue;
        }
        let ap = a.mod_u64(p as u64) as u32;
        let bp = b.mod_u64(p as u64) as u32;
        let mut roots = [0u32; 2];
        let count = if ap != 0 {
            let Some(inv) = inv_u32(ap, p) else { continue };
            roots[0] = mulmod_u32((e.sqrt_n + p - bp) % p, inv, p);
            roots[1] = mulmod_u32(((p - e.sqrt_n) % p + p - bp) % p, inv, p);
            if roots[0] == roots[1] { 1 } else { 2 }
        } else {
            let denom = (2 * bp as u64 % p as u64) as u32;
            let Some(inv) = inv_u32(denom, p) else {
                continue;
            };
            let cm = c.mod_u64(p as u64) as u32;
            let signed_c = if csign && cm != 0 { p - cm } else { cm };
            roots[0] = mulmod_u32(if signed_c == 0 { 0 } else { p - signed_c }, inv, p);
            1
        };
        let weight = ((32 - p.leading_zeros()) * 2) as u16;
        for &root in &roots[..count] {
            let start = (root as i64 - (-(ctx.interval as i64)).rem_euclid(p as i64))
                .rem_euclid(p as i64) as usize;
            for pos in (start..len).step_by(p as usize) {
                scores[pos] = scores[pos].saturating_add(weight)
            }
        }
        let _ = index;
    }
    let threshold = ((ctx.n.bit_len().saturating_sub(a.bit_len() + 55)) * 2).clamp(80, 230) as u16;
    for (pos, &score) in scores.iter().enumerate() {
        if score < threshold {
            continue;
        }
        let x = pos as i64 - ctx.interval as i64;
        let ax = a.checked_mul(&Natural::from_u64(x.unsigned_abs())).unwrap();
        let (t, tneg) = if x >= 0 {
            (ax.checked_add(b).unwrap(), false)
        } else if ax >= *b {
            (ax.wrapping_sub(b), true)
        } else {
            (b.wrapping_sub(&ax), false)
        };
        let tt = t.checked_mul(&t).unwrap();
        let (mut q, sign) = if tt >= ctx.n {
            (tt.wrapping_sub(&ctx.n), false)
        } else {
            (ctx.n.wrapping_sub(&tt), true)
        };
        let (qdiv, rem) = q.div_rem(a).unwrap();
        if !rem.is_zero() {
            continue;
        }
        q = qdiv;
        if q.is_zero() {
            continue;
        }
        let mut powers: Vec<(u32, u16)> = aidx.iter().copied().map(|i| (i, 1)).collect();
        for (i, e) in ctx.base.iter().enumerate() {
            let mut count = 0;
            loop {
                let (d, r) = q.div_rem_u64(e.prime as u64).unwrap();
                if r != 0 {
                    break;
                }
                q = d;
                count += 1
            }
            if count != 0 {
                if let Some(v) = powers.iter_mut().find(|v| v.0 == i as u32) {
                    v.1 += count
                } else {
                    powers.push((i as u32, count))
                }
            }
        }
        let large = if q.is_one() {
            None
        } else if q.bit_len() <= 64
            && q.as_parts()[0] <= ctx.large_limit
            && is_prime64(q.as_parts()[0])
        {
            Some(q.as_parts()[0])
        } else {
            continue;
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
}

fn to_column(r: Relation) -> Column {
    Column {
        root: r.root,
        sign: r.sign,
        powers: r.powers.into_iter().map(|(i, e)| (i, e as u32)).collect(),
        extra_sqrt: 1,
    }
}
fn combine(a: Relation, b: Relation, lp: u64, n: &Natural<16>) -> Column {
    let mut map = BTreeMap::<u32, u32>::new();
    for (i, e) in a.powers.into_iter().chain(b.powers) {
        *map.entry(i).or_default() += e as u32
    }
    Column {
        root: a.root.mul_mod(&b.root, n),
        sign: a.sign ^ b.sign,
        powers: map.into_iter().collect(),
        extra_sqrt: lp,
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
#[cfg(any(unix, windows))]
fn small_primes(limit: u32) -> Vec<u32> {
    let mut ps = Vec::new();
    for n in 2..=limit {
        if n == 2 || n % 2 == 1 && ps.iter().take_while(|&&p| p <= n / p).all(|&p| n % p != 0) {
            ps.push(n)
        }
    }
    ps
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
        let p = Natural::<16>::from_u64(18_446_744_073_709_551_557);
        let q = Natural::<16>::from_u64(18_446_744_073_709_551_533);
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
        let p = Natural::<16>::from_u64(18_446_744_073_709_551_557);
        let q = Natural::<16>::from_u64(18_446_744_073_709_551_533);
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
        let p = Natural::<16>::from_u64(18_446_744_073_709_551_557);
        let q = Natural::<16>::from_u64(18_446_744_073_709_551_533);
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
