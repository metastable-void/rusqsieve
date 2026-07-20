//! Recursive factorization policy and portable session state.
use crate::f2::LinearAlgebraError;
use crate::progress::*;
use crate::qs::{QsConfig, prepare_siqs};
use crate::work::{WorkJob, WorkResult};
use crate::{Natural, PrimalityConfig, PrimeFactors, is_probable_prime};
use core::fmt;
use core::num::NonZero;

#[derive(Clone, Copy, Debug)]
pub enum Parallelism {
    Auto,
    Exact(NonZero<usize>),
}
#[derive(Clone, Debug)]
pub enum SmallFactorMethod {
    None,
    PollardRho { attempts: u32, iteration_limit: u64 },
}
#[derive(Clone, Debug, Default)]
pub struct FactorLimits {
    pub max_relations: Option<usize>,
    pub max_partial_relations: Option<usize>,
    pub max_matrix_nonzeros: Option<usize>,
    pub max_memory_bytes: Option<usize>,
    pub max_polynomial_batches: Option<u64>,
    pub max_pollard_rho_iterations: Option<u64>,
}
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
impl Default for FactorConfig {
    fn default() -> Self {
        Self {
            parallelism: Parallelism::Auto,
            primality: PrimalityConfig::default(),
            trial_division_limit: 10_000,
            small_factor_method: SmallFactorMethod::PollardRho {
                attempts: 24,
                iteration_limit: 200_000,
            },
            qs: QsConfig::default(),
            limits: FactorLimits::default(),
            seed: [0x51; 32],
            progress_reporting: ProgressReportingConfig::default(),
        }
    }
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResourceLimitKind {
    Relations,
    PartialRelations,
    MatrixNonzeros,
    Memory,
    PolynomialBatches,
    PollardRhoIterations,
}
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
impl fmt::Display for FactorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroHasNoPrimeFactorization => f.write_str("zero has no prime factorization"),
            Self::CapacityExceeded => f.write_str("integer capacity exceeded"),
            Self::ResourceLimit(k) => write!(f, "resource limit exceeded: {k:?}"),
            Self::NoNontrivialFactor => f.write_str("no nontrivial factor found"),
            Self::InsufficientRelations => f.write_str("insufficient quadratic-sieve relations"),
            Self::LinearAlgebra(e) => write!(f, "linear algebra failed: {e}"),
            Self::InvalidRelation => f.write_str("invalid quadratic-sieve relation"),
            Self::InvalidDependency => f.write_str("invalid matrix dependency"),
            Self::Cancelled => f.write_str("factorization cancelled"),
            Self::Stalled => f.write_str("factorization stalled"),
            Self::InternalInvariant(s) => write!(f, "internal invariant failed: {s}"),
        }
    }
}
impl std::error::Error for FactorError {}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressAction {
    Continue,
    Cancel,
}
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
#[derive(Clone, Copy, Debug)]
pub struct LocalWorkBudget {
    pub factor_base_candidates: usize,
    pub relations_to_combine: usize,
    pub matrix_items: usize,
    pub elimination_steps: usize,
    pub primality_steps: usize,
}
impl Default for LocalWorkBudget {
    fn default() -> Self {
        Self {
            factor_base_candidates: 4096,
            relations_to_combine: 4096,
            matrix_items: 4096,
            elimination_steps: 4096,
            primality_steps: 64,
        }
    }
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AdvanceOutcome {
    Progressed,
    NeedsWorkers,
    Complete,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SubmitOutcome {
    Accepted,
    IgnoredObsolete,
}

pub struct FactorSession<const P: usize = 16> {
    input: Natural<P>,
    config: FactorConfig,
    phase: SessionPhase,
    progress: ProgressSnapshot,
    factors: Option<PrimeFactors<P>>,
    error: Option<FactorError>,
    generation: u64,
}
impl<const P: usize> FactorSession<P> {
    pub fn new(input: Natural<P>, config: FactorConfig) -> Result<Self, FactorError> {
        if input.is_zero() {
            return Err(FactorError::ZeroHasNoPrimeFactorization);
        }
        let bits = input.bit_len();
        Ok(Self {
            input,
            config,
            phase: SessionPhase::Preprocessing,
            progress: ProgressSnapshot::initial(bits),
            factors: None,
            error: None,
            generation: 0,
        })
    }
    pub fn phase(&self) -> SessionPhase {
        self.phase
    }
    pub fn progress(&self) -> &ProgressSnapshot {
        &self.progress
    }
    pub fn progress_revision(&self) -> u64 {
        self.progress.revision
    }
    pub fn advance_local(&mut self, _: LocalWorkBudget) -> Result<AdvanceOutcome, FactorError> {
        if self.phase == SessionPhase::Complete {
            return Ok(AdvanceOutcome::Complete);
        }
        if self.phase == SessionPhase::Failed {
            return Err(self.error.take().unwrap_or(FactorError::Stalled));
        }
        match factor_complete(self.input.clone(), &self.config) {
            Ok(f) => {
                self.factors = Some(f);
                self.phase = SessionPhase::Complete;
                self.progress.revision += 1;
                self.progress.phase = ProgressPhase::Complete;
                let factors = self.factors.as_ref().unwrap();
                self.progress.amount = ProgressAmount {
                    completed: 1,
                    total: ProgressTotal::Exact(1),
                    unit: ProgressUnit::Tasks,
                };
                self.progress.detail = ProgressDetail::Complete(CompleteProgress {
                    prime_factors: factors.len() as u64,
                    factors_with_multiplicity: factors.iter().map(|(_, e)| e.get() as u64).sum(),
                });
                Ok(AdvanceOutcome::Complete)
            }
            Err(e) => {
                self.phase = SessionPhase::Failed;
                Err(e)
            }
        }
    }
    pub fn take_jobs(&mut self, _: usize) -> Result<Vec<WorkJob>, FactorError> {
        Ok(Vec::new())
    }
    pub fn submit(&mut self, result: WorkResult<P>) -> Result<SubmitOutcome, FactorError> {
        let generation = match result {
            WorkResult::Sieve(r) => r.header.generation,
            WorkResult::MatrixMultiply(r) => r.header.generation,
            WorkResult::Failed(r) => r.header.generation,
        };
        Ok(if generation == self.generation {
            SubmitOutcome::Accepted
        } else {
            SubmitOutcome::IgnoredObsolete
        })
    }
    pub fn is_finished(&self) -> bool {
        matches!(self.phase, SessionPhase::Complete | SessionPhase::Failed)
    }
    pub fn take_factors(self) -> Result<PrimeFactors<P>, FactorError> {
        self.factors
            .ok_or(self.error.unwrap_or(FactorError::Stalled))
    }
}

pub(crate) fn factor_complete<const P: usize>(
    input: Natural<P>,
    config: &FactorConfig,
) -> Result<PrimeFactors<P>, FactorError> {
    if input.is_zero() {
        return Err(FactorError::ZeroHasNoPrimeFactorization);
    }
    let mut out = PrimeFactors::new();
    factor_node(input, 1, config, &mut out, 0)?;
    Ok(out)
}
fn factor_node<const P: usize>(
    mut n: Natural<P>,
    multiplicity: usize,
    cfg: &FactorConfig,
    out: &mut PrimeFactors<P>,
    depth: usize,
) -> Result<(), FactorError> {
    if depth > Natural::<P>::BITS {
        return Err(FactorError::Stalled);
    }
    if n.is_one() {
        return Ok(());
    }
    for p in primes_to(cfg.trial_division_limit) {
        let mut count = 0;
        loop {
            let (q, r) = n.div_rem_u64(p as u64).unwrap();
            if r != 0 {
                break;
            }
            n = q;
            count += 1
        }
        if count != 0 {
            out.insert_count(Natural::from_u64(p as u64), count * multiplicity)
        }
        if n.is_one() {
            return Ok(());
        }
    }
    if is_probable_prime(&n, &cfg.primality) {
        out.insert_count(n, multiplicity);
        return Ok(());
    }
    if let Some((root, power)) = n.perfect_power() {
        return factor_node(
            root,
            multiplicity
                .checked_mul(power as usize)
                .ok_or(FactorError::CapacityExceeded)?,
            cfg,
            out,
            depth + 1,
        );
    }
    let factor = match &cfg.small_factor_method {
        SmallFactorMethod::None => None,
        SmallFactorMethod::PollardRho {
            attempts,
            iteration_limit,
        } => {
            let cap = cfg
                .limits
                .max_pollard_rho_iterations
                .map_or(*iteration_limit, |v| v.min(*iteration_limit));
            pollard_rho(&n, *attempts, cap, &cfg.seed)
        }
    }
    .or_else(|| reference_qs_factor(&n, cfg).ok())
    .ok_or(FactorError::NoNontrivialFactor)?;
    if factor.is_one() || factor == n {
        return Err(FactorError::NoNontrivialFactor);
    }
    let other = n
        .div_rem(&factor)
        .ok_or(FactorError::InternalInvariant("factor division"))?
        .0;
    factor_node(factor, multiplicity, cfg, out, depth + 1)?;
    factor_node(other, multiplicity, cfg, out, depth + 1)
}
fn primes_to(limit: u32) -> Vec<u32> {
    let mut ps = Vec::new();
    for n in 2..=limit {
        if n == 2 || n % 2 != 0 && ps.iter().take_while(|&&p| p <= n / p).all(|&p| n % p != 0) {
            ps.push(n)
        }
    }
    ps
}
fn pollard_rho<const P: usize>(
    n: &Natural<P>,
    attempts: u32,
    limit: u64,
    seed: &[u8; 32],
) -> Option<Natural<P>> {
    if n.is_even() {
        return Some(Natural::from_u64(2));
    }
    let mut state = u64::from_le_bytes(seed[..8].try_into().unwrap())
        ^ n.as_parts().first().copied().unwrap_or(0);
    for attempt in 0..attempts {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let c = Natural::from_u64(1 + (state | 1).wrapping_add(attempt as u64))
            .div_rem(n)?
            .1;
        let mut x = Natural::from_u64(2 + state).div_rem(n)?.1;
        let mut y = x.clone();
        for _ in 0..limit {
            x = x.mul_mod(&x, n).add_mod(&c, n);
            y = y.mul_mod(&y, n).add_mod(&c, n);
            y = y.mul_mod(&y, n).add_mod(&c, n);
            let d = if x >= y {
                x.wrapping_sub(&y)
            } else {
                y.wrapping_sub(&x)
            }
            .gcd(n);
            if d.is_one() {
                continue;
            }
            if d != *n {
                return Some(d);
            }
            break;
        }
    }
    None
}

fn reference_qs_factor<const P: usize>(
    n: &Natural<P>,
    cfg: &FactorConfig,
) -> Result<Natural<P>, FactorError> {
    let ctx = match prepare_siqs(n, &cfg.qs) {
        Ok(c) => c,
        Err(crate::qs::QsError::FactorFound(p)) => return Ok(Natural::from_u64(p as u64)),
        Err(_) => return Err(FactorError::NoNontrivialFactor),
    };
    let needed = ctx.factor_base().len() + cfg.qs.relation_surplus.absolute + 1;
    let max = cfg.limits.max_relations.unwrap_or(needed.saturating_mul(8));
    if max < needed {
        return Err(FactorError::ResourceLimit(ResourceLimitKind::Relations));
    }
    let max_batches = cfg.limits.max_polynomial_batches.unwrap_or(100_000);
    let mut relations = Vec::new();
    let mut scratch = crate::qs::SieveScratch::default();
    let mut first = 0u32;
    let mut batches = 0u64;
    while relations.len() < needed && batches < max_batches {
        let job = crate::work::SieveJob {
            header: crate::work::JobHeader {
                job_id: first as u64,
                generation: 0,
                context_id: ctx.context_id(),
            },
            family: 0,
            first_polynomial: first,
            polynomial_count: cfg.qs.polynomial_batch_size.max(1),
        };
        batches += 1;
        first = first.saturating_add(job.polynomial_count);
        for r in crate::qs::sieve_job(&ctx, &job, &mut scratch).relations {
            if matches!(r.large_primes, crate::qs::LargePrimePart::None) {
                relations.push(r);
                if relations.len() == max {
                    break;
                }
            }
        }
        if first == u32::MAX || relations.len() == max {
            break;
        }
    }
    if relations.len() <= ctx.factor_base().len() {
        if batches == max_batches && cfg.limits.max_polynomial_batches.is_some() {
            return Err(FactorError::ResourceLimit(
                ResourceLimitKind::PolynomialBatches,
            ));
        }
        return Err(FactorError::InsufficientRelations);
    }
    let columns: Vec<Vec<u32>> = relations
        .iter()
        .map(|r| {
            let mut rows = Vec::new();
            if r.sign {
                rows.push(0)
            }
            for p in &r.factors {
                if p.exponent & 1 != 0 {
                    rows.push(p.factor_base_index + 1)
                }
            }
            rows
        })
        .collect();
    let matrix = crate::f2::SparseBinaryMatrix::from_columns(ctx.factor_base().len() + 1, &columns)
        .map_err(|_| FactorError::InvalidDependency)?;
    let deps = matrix.dense_dependencies();
    for dep in deps.iter() {
        if !matrix.verify_dependency(dep) {
            return Err(FactorError::InvalidDependency);
        }
        let mut x = Natural::ONE;
        let mut sums = vec![0u32; ctx.factor_base().len()];
        for (col, r) in relations.iter().enumerate() {
            if (dep[col / 64] >> (col % 64)) & 1 == 0 {
                continue;
            }
            x = x.mul_mod(&r.square_root, n);
            for pp in &r.factors {
                sums[pp.factor_base_index as usize] += pp.exponent as u32
            }
        }
        let mut y = Natural::ONE;
        for (e, &sum) in ctx.factor_base().entries().iter().zip(&sums) {
            for _ in 0..sum / 2 {
                y = y.mul_mod(&Natural::from_u64(e.prime as u64), n)
            }
        }
        let delta = if x >= y {
            x.wrapping_sub(&y)
        } else {
            y.wrapping_sub(&x)
        };
        let g = delta.gcd(n);
        if !g.is_one() && g != *n {
            return Ok(g);
        }
        let sum = x.add_mod(&y, n);
        let g = sum.gcd(n);
        if !g.is_one() && g != *n {
            return Ok(g);
        }
    }
    Err(FactorError::NoNontrivialFactor)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn factors_edge_cases_and_semiprime() {
        assert!(
            factor_complete(Natural::<2>::ONE, &FactorConfig::default())
                .unwrap()
                .is_empty()
        );
        let n = Natural::<2>::from_u64(1_000_003 * 1_000_033);
        let factors = factor_complete(n.clone(), &FactorConfig::default()).unwrap();
        assert!(factors.verify_product(&n));
        assert_eq!(factors.len(), 2);
    }
    #[test]
    fn prime_power_multiplicity() {
        let n = Natural::<2>::from_u64(3u64.pow(20));
        let factors = factor_complete(n.clone(), &FactorConfig::default()).unwrap();
        assert_eq!(factors.get(&Natural::from_u64(3)).unwrap().get(), 20);
        assert!(factors.verify_product(&n));
    }
    #[test]
    fn quadratic_sieve_fallback() {
        let n = Natural::<2>::from_u64(1009 * 1013);
        let qs = QsConfig {
            factor_base_bound: crate::qs::AutoOr::Value(200),
            polynomial_batch_size: 128,
            ..QsConfig::default()
        };
        let limits = FactorLimits {
            max_relations: Some(400),
            ..FactorLimits::default()
        };
        let config = FactorConfig {
            trial_division_limit: 97,
            small_factor_method: SmallFactorMethod::None,
            qs,
            limits,
            ..FactorConfig::default()
        };
        let factors = factor_complete(n.clone(), &config).unwrap();
        assert!(factors.verify_product(&n));
        assert_eq!(factors.len(), 2);
    }
}
