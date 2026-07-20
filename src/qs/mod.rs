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
                return Err(FactorBaseError::FoundFactor(p));
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
        let (factor_base_bound, sieve_half_width, lp_allowance) = match bits {
            0..=100 => (3_000, 32_768, 16),
            101..=128 => (6_000, 32_768, 18),
            129..=160 => (20_000, 65_536, 22),
            161..=192 => (28_000, 65_536, 22),
            193..=224 => (60_000, 131_072, 26),
            225..=248 => (150_000, 131_072, 30),
            _ => (300_000, 131_072, 34),
        };
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
}
pub fn sieve_job<const P: usize>(
    ctx: &SieveContext<P>,
    job: &crate::work::SieveJob,
    scratch: &mut SieveScratch,
) -> crate::work::SieveResult<P> {
    scratch.scores.clear();
    scratch.candidates.clear();
    scratch.factor_scratch.clear();
    let mut relations = Vec::new();
    let root = ctx.working_n.ceil_sqrt();
    let mut metrics = crate::work::SieveJobMetrics {
        polynomial_families: 1,
        polynomials: job.polynomial_count as u64,
        ..Default::default()
    };
    for k in 0..job.polynomial_count {
        let offset = job.first_polynomial as u64
            + k as u64
            + job.family.saturating_mul(job.polynomial_count as u64);
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
        let mut powers = Vec::new();
        for (index, e) in ctx.factor_base.entries.iter().enumerate() {
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
                powers.push(PrimePower {
                    factor_base_index: index as u32,
                    exponent: count,
                })
            }
        }
        metrics.candidates_tested += 1;
        let lp = if q.is_one() {
            metrics.full_relations += 1;
            LargePrimePart::None
        } else if q.bit_len() <= 64 {
            let v = q.as_parts()[0];
            metrics.single_large_prime_relations += 1;
            LargePrimePart::One(v)
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
                polynomial: job.first_polynomial + k,
                position: offset as i32,
            },
        };
        if verify_relation(ctx, &rel) {
            relations.push(rel)
        }
    }
    metrics.sieve_positions = job.polynomial_count as u64;
    crate::work::SieveResult {
        header: job.header,
        relations,
        metrics,
    }
}
