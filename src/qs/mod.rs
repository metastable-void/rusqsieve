//! Portable quadratic-sieve setup, relation generation, and verification.
use crate::f2::MatrixConfig;
use crate::progress::FactorBaseProgress;
use crate::{Natural, legendre_u32, tonelli_shanks_u32};
use core::fmt;

#[derive(Clone, Copy, Debug)]
pub enum AutoOr<T> {
    Auto,
    Value(T),
}
#[derive(Clone, Copy, Debug)]
pub enum MultiplierChoice {
    Auto,
    Value(u32),
}
#[derive(Clone, Debug)]
pub struct LargePrimeConfig {
    pub single_limit: AutoOr<u64>,
    pub double_product_limit: AutoOr<u64>,
    pub enable_double: bool,
}
#[derive(Clone, Copy, Debug)]
pub struct RelationSurplus {
    pub absolute: usize,
    pub percent: u8,
}
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
impl Default for QsConfig {
    fn default() -> Self {
        Self {
            multiplier: MultiplierChoice::Auto,
            factor_base_bound: AutoOr::Auto,
            sieve_half_width: AutoOr::Auto,
            polynomial_batch_size: 64,
            large_primes: LargePrimeConfig {
                single_limit: AutoOr::Auto,
                double_product_limit: AutoOr::Auto,
                enable_double: true,
            },
            relation_surplus: RelationSurplus {
                absolute: 32,
                percent: 10,
            },
            sieve_score_scale: 8,
            candidate_slack: 16,
            matrix: MatrixConfig::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FactorBaseEntry {
    pub prime: u32,
    pub log_prime: u8,
    pub sqrt_n: u32,
}
#[derive(Clone, Debug)]
pub struct FactorBase {
    entries: Box<[FactorBaseEntry]>,
}
impl FactorBase {
    pub fn entries(&self) -> &[FactorBaseEntry] {
        &self.entries
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
#[derive(Clone, Debug)]
pub enum FactorBaseError {
    InvalidBound,
    FoundFactor(u32),
    NotFinished,
}
impl fmt::Display for FactorBaseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "factor-base error: {self:?}")
    }
}
impl std::error::Error for FactorBaseError {}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FactorBaseBuildStatus {
    InProgress,
    Complete,
}
pub struct FactorBaseBuilder<const P: usize> {
    n: Natural<P>,
    bound: u32,
    candidate: u32,
    entries: Vec<FactorBaseEntry>,
    tested: u64,
    nonresidue: u64,
    finished: bool,
    /// Knuth-Schroeppel multiplier `k` (1 if none). Primes dividing `k` divide the working
    /// modulus `k·n` but are not factors of `n`; they are added as ramified factor-base entries
    /// (`sqrt_n = 0`) rather than reported as `FoundFactor`.
    multiplier: u64,
}
impl<const P: usize> FactorBaseBuilder<P> {
    pub fn new(n: Natural<P>, bound: u32) -> Result<Self, FactorBaseError> {
        if bound < 2 {
            return Err(FactorBaseError::InvalidBound);
        }
        Ok(Self {
            n,
            bound,
            candidate: 2,
            multiplier: 1,
            entries: Vec::new(),
            tested: 0,
            nonresidue: 0,
            finished: false,
        })
    }
    pub fn step(&mut self, budget: usize) -> Result<FactorBaseBuildStatus, FactorBaseError> {
        for _ in 0..budget {
            if self.candidate > self.bound {
                self.finished = true;
                return Ok(FactorBaseBuildStatus::Complete);
            }
            let p = self.candidate;
            self.candidate = if p == 2 { 3 } else { p.saturating_add(2) };
            if !prime_u32(p) {
                continue;
            }
            self.tested += 1;
            let r = self.n.mod_u64(p as u64) as u32;
            if r == 0 && self.n != Natural::from_u64(p as u64) {
                // `p | working`. If `p | k` it only divides the multiplier, not `n` — fall through
                // and add it as a ramified prime (`r == 0` ⇒ `sqrt_n = 0`). Otherwise it is a real
                // factor of `n`.
                if !self.multiplier.is_multiple_of(p as u64) {
                    return Err(FactorBaseError::FoundFactor(p));
                }
            }
            if p == 2 || legendre_u32(r, p) >= 0 {
                self.entries.push(FactorBaseEntry {
                    prime: p,
                    log_prime: ((p as f64).ln() * 8.0).round().min(255.0) as u8,
                    sqrt_n: tonelli_shanks_u32(r, p).unwrap_or(r & 1),
                })
            } else {
                self.nonresidue += 1
            }
        }
        Ok(FactorBaseBuildStatus::InProgress)
    }
    pub fn progress(&self) -> FactorBaseProgress {
        FactorBaseProgress {
            bound: self.bound,
            searched_through: self.candidate.min(self.bound),
            primes_tested: self.tested,
            primes_accepted: self.entries.len() as u64,
            nonresidue_primes: self.nonresidue,
        }
    }
    pub fn finish(self) -> Result<FactorBase, FactorBaseError> {
        if !self.finished {
            return Err(FactorBaseError::NotFinished);
        }
        Ok(FactorBase {
            entries: self.entries.into_boxed_slice(),
        })
    }
}
fn prime_u32(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    if n.is_multiple_of(2) {
        return n == 2;
    }
    let mut d = 3;
    while d <= n / d {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 2
    }
    true
}

#[derive(Clone, Debug)]
pub struct PolynomialPlan {
    pub first_x_offset: u64,
}
#[derive(Clone, Debug)]
pub struct SieveContext<const P: usize> {
    pub(crate) n: Natural<P>,
    pub(crate) working_n: Natural<P>,
    pub(crate) multiplier: u32,
    pub(crate) factor_base: FactorBase,
    pub(crate) polynomial_plan: PolynomialPlan,
    pub(crate) config: QsConfig,
    pub(crate) context_id: u64,
}
impl<const P: usize> SieveContext<P> {
    pub fn n(&self) -> &Natural<P> {
        &self.n
    }
    pub fn factor_base(&self) -> &FactorBase {
        &self.factor_base
    }
    pub fn context_id(&self) -> u64 {
        self.context_id
    }
    pub fn multiplier(&self) -> u32 {
        self.multiplier
    }
    pub fn polynomial_plan(&self) -> &PolynomialPlan {
        &self.polynomial_plan
    }
    pub fn config(&self) -> &QsConfig {
        &self.config
    }
}
#[derive(Clone, Debug)]
pub enum QsError {
    InputTooSmall,
    Capacity,
    FactorFound(u32),
    FactorBase(FactorBaseError),
}
impl fmt::Display for QsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "quadratic-sieve setup error: {self:?}")
    }
}
impl std::error::Error for QsError {}
pub fn prepare_siqs<const P: usize>(
    n: &Natural<P>,
    config: &QsConfig,
) -> Result<SieveContext<P>, QsError> {
    if *n < Natural::from_u64(2) {
        return Err(QsError::InputTooSmall);
    }
    let multiplier = match config.multiplier {
        MultiplierChoice::Auto => 1,
        MultiplierChoice::Value(k) => k.max(1),
    };
    let working = n
        .checked_mul(&Natural::from_u64(multiplier as u64))
        .ok_or(QsError::Capacity)?;
    let bound = match config.factor_base_bound {
        AutoOr::Value(v) => v,
        AutoOr::Auto => parameters::factor_base_bound(n.bit_len()),
    };
    let mut b = FactorBaseBuilder::new(working.clone(), bound).map_err(QsError::FactorBase)?;
    b.multiplier = multiplier as u64;
    loop {
        match b.step(4096) {
            Ok(FactorBaseBuildStatus::Complete) => break,
            Ok(FactorBaseBuildStatus::InProgress) => {}
            Err(FactorBaseError::FoundFactor(p)) => return Err(QsError::FactorFound(p)),
            Err(e) => return Err(QsError::FactorBase(e)),
        }
    }
    let base = b.finish().map_err(QsError::FactorBase)?;
    let mut id = 0xcbf29ce484222325;
    for &x in n.as_parts() {
        id ^= x;
        id = id.wrapping_mul(0x100000001b3)
    }
    Ok(SieveContext {
        n: n.clone(),
        working_n: working,
        multiplier,
        factor_base: base,
        polynomial_plan: PolynomialPlan { first_x_offset: 0 },
        config: config.clone(),
        context_id: id,
    })
}
pub mod parameters {
    pub fn factor_base_bound(bits: usize) -> u32 {
        match bits {
            0..=40 => 200,
            41..=64 => 2_000,
            65..=96 => 10_000,
            97..=128 => 30_000,
            _ => 100_000,
        }
    }
    pub fn sieve_half_width(bits: usize) -> u32 {
        match bits {
            0..=64 => 4096,
            65..=128 => 32768,
            _ => 131072,
        }
    }

    /// Tuned SIQS engine parameters per input bit-length. These were selected by
    /// benchmarking balanced semiprimes against `flintqs`. `lp_allowance` is the
    /// large-prime cofactor budget (bits) used to derive the sieve threshold
    /// `2 * (log2|g(x)| - lp_allowance)`.
    #[derive(Clone, Copy, Debug)]
    pub struct EngineParams {
        pub factor_base_bound: u32,
        pub sieve_half_width: u32,
        pub lp_allowance: usize,
    }
    pub fn engine_params(bits: usize) -> EngineParams {
        #[allow(unused_mut)]
        let (mut factor_base_bound, mut sieve_half_width, lp_allowance) = match bits {
            0..=100 => (3_000, 32_768, 16),
            101..=128 => (6_000, 32_768, 18),
            129..=160 => (40_000, 65_536, 22),
            161..=192 => (60_000, 65_536, 22),
            193..=224 => (60_000, 131_072, 26),
            225..=248 => (150_000, 131_072, 30),
            _ => (300_000, 131_072, 34),
        };
        // Tuning overrides (experimentation only; unset in production builds).
        #[cfg(any(unix, windows))]
        {
            if let Some(v) = std::env::var("RUSQSIEVE_FB_BOUND").ok().and_then(|s| s.parse().ok()) {
                factor_base_bound = v;
            }
            if let Some(v) = std::env::var("RUSQSIEVE_HALFW").ok().and_then(|s| s.parse().ok()) {
                sieve_half_width = v;
            }
        }
        EngineParams {
            factor_base_bound,
            sieve_half_width,
            lp_allowance,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PrimePower {
    pub factor_base_index: u32,
    pub exponent: u16,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LargePrimePart {
    None,
    One(u64),
    Two(u64, u64),
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RelationSource {
    pub family: u64,
    pub polynomial: u32,
    pub position: i32,
}
#[derive(Clone, Debug)]
pub struct RawRelation<const P: usize> {
    pub square_root: Natural<P>,
    pub sign: bool,
    pub factors: Box<[PrimePower]>,
    pub large_primes: LargePrimePart,
    pub source: RelationSource,
}
pub type RelationId = u64;
pub type SparseParity = Box<[u32]>;
#[derive(Clone, Debug)]
pub struct CombinedRelation {
    pub sources: Box<[RelationId]>,
    pub parity: SparseParity,
}
pub fn verify_relation<const P: usize>(ctx: &SieveContext<P>, r: &RawRelation<P>) -> bool {
    let n = &ctx.n;
    if n.is_zero() {
        return false;
    }
    let mut rhs = Natural::ONE.div_rem(n).unwrap().1;
    if r.sign && !rhs.is_zero() {
        rhs = n.wrapping_sub(&rhs)
    }
    for pp in &r.factors {
        let Some(e) = ctx.factor_base.entries.get(pp.factor_base_index as usize) else {
            return false;
        };
        let base = Natural::from_u64(e.prime as u64);
        for _ in 0..pp.exponent {
            rhs = rhs.mul_mod(&base, n)
        }
    }
    match r.large_primes {
        LargePrimePart::None => {}
        LargePrimePart::One(p) => rhs = rhs.mul_mod(&Natural::from_u64(p), n),
        LargePrimePart::Two(p, q) => {
            rhs = rhs.mul_mod(&Natural::from_u64(p), n);
            rhs = rhs.mul_mod(&Natural::from_u64(q), n)
        }
    }
    r.square_root.mul_mod(&r.square_root, n) == rhs
}

#[derive(Default)]
pub struct SieveScratch {
    pub(crate) scores: Vec<u8>,
    pub(crate) candidates: Vec<u32>,
    pub(crate) factor_scratch: Vec<(u32, u16)>,
    /// `ceil_sqrt(working_n) mod p` per factor-base prime, precomputed per job.
    pub(crate) roots_mod_p: Vec<u32>,
}

/// Portable deterministic Miller–Rabin for `u64` (exact for all `n < 2^64`).
fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &p in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let mulmod = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let powmod = |mut a: u64, mut e: u64| {
        let mut r = 1u64;
        while e != 0 {
            if e & 1 == 1 {
                r = mulmod(r, a);
            }
            a = mulmod(a, a);
            e >>= 1;
        }
        r
    };
    let mut d = n - 1;
    let mut s = 0u32;
    while d & 1 == 0 {
        d >>= 1;
        s += 1;
    }
    'witness: for &a in &[2u64, 325, 9375, 28178, 450775, 9780504, 1795265022] {
        let a = a % n;
        if a == 0 {
            continue;
        }
        let mut x = powmod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// Portable logarithmic quadratic sieve (SPEC §12.6).
///
/// Each job sieves `polynomial_count` contiguous segments of the `Q(x) = x² - N`
/// line (`x = ceil_sqrt(N) + offset`, `N = working_n`). For every factor-base
/// prime `p`, `log(p)` is added to the byte score array at the two roots
/// `x ≡ ±sqrt_n (mod p)`; only positions whose accumulated score clears a
/// per-segment threshold (`≈ log|Q| − large-prime allowance − slack`) are
/// trial-divided and classified. This replaces the previous approach of
/// trial-dividing every candidate.
pub fn sieve_job<const P: usize>(
    ctx: &SieveContext<P>,
    job: &crate::work::SieveJob,
    scratch: &mut SieveScratch,
) -> crate::work::SieveResult<P> {
    let SieveScratch {
        scores,
        candidates,
        factor_scratch,
        roots_mod_p,
    } = scratch;
    factor_scratch.clear();
    let fb = &ctx.factor_base.entries;
    let root = ctx.working_n.ceil_sqrt();
    let mut relations = Vec::new();
    let mut metrics = crate::work::SieveJobMetrics {
        polynomial_families: 1,
        polynomials: job.polynomial_count as u64,
        ..Default::default()
    };

    let seg = match ctx.config.sieve_half_width {
        AutoOr::Value(v) => v.max(64) as usize,
        AutoOr::Auto => parameters::sieve_half_width(ctx.n.bit_len()) as usize,
    };
    // Log scale must match `FactorBaseEntry::log_prime`, which is `round(ln(p) * 8)`.
    const LOG_SCALE: f64 = 8.0;
    let last_prime = fb.last().map(|e| e.prime as u64).unwrap_or(2);
    let lp_limit = match ctx.config.large_primes.single_limit {
        AutoOr::Value(v) => v,
        AutoOr::Auto => last_prime.saturating_mul(last_prime),
    };
    let lp_bits = 64 - lp_limit.max(1).leading_zeros();
    let lp_allow = (lp_bits as f64 * core::f64::consts::LN_2 * LOG_SCALE) as i32;
    let slack = ctx.config.candidate_slack as i32;

    // ceil_sqrt(working_n) mod p is invariant across all segments of this job.
    roots_mod_p.clear();
    roots_mod_p.reserve(fb.len());
    for e in fb.iter() {
        roots_mod_p.push(root.mod_u64(e.prime as u64) as u32);
    }

    for k in 0..job.polynomial_count {
        let seg_index = job.first_polynomial as u64
            + k as u64
            + job.family.saturating_mul(job.polynomial_count as u64);
        let seg_start = seg_index.saturating_mul(seg as u64);
        let Some(x0) = root.checked_add(&Natural::from_u64(seg_start)) else {
            continue;
        };
        let Some(x0sq) = x0.checked_mul(&x0) else {
            continue;
        };
        let Some(q0) = x0sq.checked_sub(&ctx.working_n) else {
            continue;
        };
        metrics.sieve_positions += seg as u64;
        // Threshold from the smallest |Q| in the segment (its start), so smooth
        // low-offset positions are never scored out.
        let target = (q0.bit_len() as f64 * core::f64::consts::LN_2 * LOG_SCALE) as i32;
        let threshold = (target - lp_allow - slack).clamp(1, 255) as u8;

        scores.clear();
        scores.resize(seg, 0);
        for (i, e) in fb.iter().enumerate() {
            let p = e.prime;
            let lg = e.log_prime;
            let rp = roots_mod_p[i];
            let ss = (seg_start % p as u64) as u32;
            // Solution offsets: offset ≡ ±sqrt_n − root (mod p).
            let off1 = (e.sqrt_n + p - rp) % p;
            let neg = if e.sqrt_n == 0 { 0 } else { p - e.sqrt_n };
            let off2 = (neg + p - rp) % p;
            let step = p as usize;
            let mut j = ((off1 + p - ss) % p) as usize;
            while j < seg {
                scores[j] = scores[j].saturating_add(lg);
                j += step;
            }
            if off2 != off1 {
                let mut j = ((off2 + p - ss) % p) as usize;
                while j < seg {
                    scores[j] = scores[j].saturating_add(lg);
                    j += step;
                }
            }
        }

        // Scan survivors above threshold, then trial-divide only those.
        candidates.clear();
        for (j, &score) in scores.iter().enumerate() {
            if score >= threshold {
                candidates.push(j as u32);
            }
        }
        for &j in candidates.iter() {
            let offset = seg_start + j as u64;
            let Some(x) = root.checked_add(&Natural::from_u64(offset)) else {
                continue;
            };
            let Some(xx) = x.checked_mul(&x) else {
                continue;
            };
            let Some(mut q) = xx.checked_sub(&ctx.working_n) else {
                continue;
            };
            if q.is_zero() {
                continue;
            }
            metrics.candidates_tested += 1;
            let mut powers = Vec::new();
            for (index, e) in fb.iter().enumerate() {
                if q.is_one() {
                    break;
                }
                let mut count = 0u16;
                loop {
                    let (d, r) = q.div_rem_u64(e.prime as u64).unwrap();
                    if r != 0 {
                        break;
                    }
                    q = d;
                    count += 1
                }
                if count != 0 {
                    powers.push(PrimePower {
                        factor_base_index: index as u32,
                        exponent: count,
                    })
                }
            }
            let lp = if q.is_one() {
                metrics.full_relations += 1;
                LargePrimePart::None
            } else if q.bit_len() <= 64 {
                let v = q.as_parts()[0];
                if v <= lp_limit && is_prime_u64(v) {
                    metrics.single_large_prime_relations += 1;
                    LargePrimePart::One(v)
                } else {
                    continue;
                }
            } else {
                continue;
            };
            let rel = RawRelation {
                square_root: x.div_rem(&ctx.n).unwrap().1,
                sign: false,
                factors: powers.into_boxed_slice(),
                large_primes: lp,
                source: RelationSource {
                    family: job.family,
                    polynomial: seg_index as u32,
                    position: offset as i32,
                },
            };
            if verify_relation(ctx, &rel) {
                relations.push(rel)
            }
        }
    }
    crate::work::SieveResult {
        header: job.header,
        relations,
        metrics,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sieve_job_logarithmic_sieve() {
        // A real log-sieve must (a) emit only verified relations and (b) trial-divide
        // far fewer positions than it sieves — not every candidate like the old kernel.
        let n = Natural::<2>::from_u64(1_000_003 * 1_000_033);
        let cfg = QsConfig {
            factor_base_bound: AutoOr::Value(500),
            ..QsConfig::default()
        };
        let ctx = prepare_siqs(&n, &cfg).unwrap();
        let mut scratch = SieveScratch::default();
        let job = crate::work::SieveJob {
            header: crate::work::JobHeader {
                job_id: 0,
                generation: 0,
                context_id: ctx.context_id(),
            },
            family: 0,
            first_polynomial: 0,
            polynomial_count: 16,
        };
        let result = sieve_job(&ctx, &job, &mut scratch);
        assert!(result.metrics.sieve_positions >= 16 * 64);
        // The sieve filtered: only a small fraction of positions were trial-divided.
        assert!(
            result.metrics.candidates_tested < result.metrics.sieve_positions,
            "candidates {} !< positions {}",
            result.metrics.candidates_tested,
            result.metrics.sieve_positions
        );
        // Some relations were found, and every emitted relation satisfies the invariant.
        assert!(!result.relations.is_empty(), "no relations found");
        for r in &result.relations {
            assert!(verify_relation(&ctx, r), "relation failed verification");
        }
        // Determinism: the same job yields the same relation count.
        let again = sieve_job(&ctx, &job, &mut scratch);
        assert_eq!(again.relations.len(), result.relations.len());
    }
}
