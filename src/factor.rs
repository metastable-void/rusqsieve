//! Configuration and errors for blocking SIQS factorization.

use core::fmt;
use core::num::NonZeroUsize;
use core::time::Duration;

use crate::progress::ProgressReportingConfig;

/// Selects the worker count used by native blocking factorization.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum Parallelism {
    /// Detect available parallelism when factorization starts.
    Auto,
    /// Use exactly this many worker threads.
    Threads(NonZeroUsize),
}

impl Parallelism {
    /// Constructs an explicit worker count, returning `None` for zero.
    #[must_use]
    pub const fn threads(count: usize) -> Option<Self> {
        match NonZeroUsize::new(count) {
            Some(count) => Some(Self::Threads(count)),
            None => None,
        }
    }
}

/// Internal SIQS tuning overrides used by the CLI's reproducible benchmark
/// environment. Library callers should normally leave these at their defaults.
#[derive(Clone, Debug, Default)]
pub(crate) struct FactorTuning {
    pub(crate) relation_percent: Option<usize>,
    pub(crate) small_skip: Option<u32>,
    pub(crate) threshold_margin: Option<i32>,
    pub(crate) threshold_adjustment: Option<i32>,
    pub(crate) factor_base_bound: Option<u32>,
    pub(crate) sieve_half_width: Option<u32>,
    pub(crate) large_prime_multiplier: Option<u32>,
    pub(crate) double_large_prime_bound: Option<u64>,
    pub(crate) profile: bool,
}

/// Configuration for native blocking factorization.
///
/// The default uses [`Parallelism::Auto`], deterministic primality witnesses,
/// and reports progress at most every 100 milliseconds within a phase.
#[derive(Clone, Debug)]
pub struct FactorConfig {
    pub(crate) parallelism: Parallelism,
    pub(crate) witness_seed: Option<[u8; 32]>,
    pub(crate) tuning: FactorTuning,
    pub(crate) progress_reporting: ProgressReportingConfig,
}

impl FactorConfig {
    /// Returns the configured parallelism policy.
    #[must_use]
    pub const fn parallelism(&self) -> Parallelism {
        self.parallelism
    }

    /// Sets the parallelism policy.
    #[must_use]
    pub fn with_parallelism(mut self, parallelism: Parallelism) -> Self {
        self.parallelism = parallelism;
        self
    }

    /// Returns the minimum interval between progress callbacks.
    #[must_use]
    pub const fn progress_interval(&self) -> Duration {
        self.progress_reporting.minimum_interval
    }

    /// Sets the minimum interval between progress callbacks.
    ///
    /// A zero duration requests every available update and can materially slow
    /// short factorizations.
    #[must_use]
    pub fn with_progress_interval(mut self, interval: Duration) -> Self {
        self.progress_reporting.minimum_interval = interval;
        self
    }

    /// Selects deterministic seeded Miller–Rabin witnesses.
    ///
    /// The same seed and input produce the same witnesses on native and Wasm
    /// targets. This is reproducibility and adversarial diversification, not
    /// entropy: callers that need an unpredictable seed must obtain it
    /// themselves. Baillie–PSW remains the primary primality safeguard.
    #[must_use]
    pub fn with_witness_seed(mut self, seed: [u8; 32]) -> Self {
        self.witness_seed = Some(seed);
        self
    }

    /// Applies CLI-only tuning overrides used to reproduce benchmark sweeps.
    ///
    /// This interface is hidden because the values are implementation details
    /// and may change between releases.
    #[doc(hidden)]
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn with_tuning_overrides(
        mut self,
        relation_percent: Option<usize>,
        small_skip: Option<u32>,
        threshold_margin: Option<i32>,
        threshold_adjustment: Option<i32>,
        factor_base_bound: Option<u32>,
        sieve_half_width: Option<u32>,
        large_prime_multiplier: Option<u32>,
        double_large_prime_bound: Option<u64>,
        profile: bool,
    ) -> Self {
        self.tuning = FactorTuning {
            relation_percent,
            small_skip,
            threshold_margin,
            threshold_adjustment,
            factor_base_bound,
            sieve_half_width,
            large_prime_multiplier,
            double_large_prime_bound,
            profile,
        };
        self
    }
}

impl Default for FactorConfig {
    fn default() -> Self {
        Self {
            parallelism: Parallelism::Auto,
            witness_seed: None,
            tuning: FactorTuning::default(),
            progress_reporting: ProgressReportingConfig::default(),
        }
    }
}

/// The kind of internal resource limit that was exceeded.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ResourceLimitKind {
    /// Maximum safe memory or worker count.
    Memory,
}

/// An error returned while factoring an integer.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum FactorError {
    /// Zero has no finite prime factorization.
    ZeroHasNoPrimeFactorization,
    /// A value exceeded the selected [`Natural`](crate::Natural) capacity.
    CapacityExceeded,
    /// Inputs above the supported 512-bit SIQS range are rejected.
    InputTooLarge,
    /// A configured internal resource limit was exceeded.
    ResourceLimit(ResourceLimitKind),
    /// No nontrivial divisor was recovered.
    NoNontrivialFactor,
    /// Sieving ended before enough usable relations were available.
    InsufficientRelations,
    /// SIQS polynomial coefficients could not be selected.
    PolynomialSelection(String),
    /// A matrix dependency failed an internal consistency check.
    InvalidDependency,
    /// A worker thread failed.
    WorkerFailure(String),
    /// The progress callback requested cancellation.
    Cancelled,
    /// An internal consistency check failed.
    InternalFailure,
}

impl fmt::Display for FactorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroHasNoPrimeFactorization => f.write_str("zero has no prime factorization"),
            Self::CapacityExceeded => f.write_str("integer capacity exceeded"),
            Self::InputTooLarge => f.write_str("SIQS supports inputs of at most 512 bits"),
            Self::ResourceLimit(kind) => write!(f, "resource limit exceeded: {kind:?}"),
            Self::NoNontrivialFactor => f.write_str("no nontrivial factor found"),
            Self::InsufficientRelations => f.write_str("insufficient quadratic-sieve relations"),
            Self::PolynomialSelection(message) => f.write_str(message),
            Self::InvalidDependency => f.write_str("invalid matrix dependency"),
            Self::WorkerFailure(message) => write!(f, "worker failed: {message}"),
            Self::Cancelled => f.write_str("factorization cancelled"),
            Self::InternalFailure => f.write_str("internal factorization invariant failed"),
        }
    }
}

impl std::error::Error for FactorError {}

/// Controls factorization from a progress callback.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressAction {
    /// Continue factorization.
    Continue,
    /// Stop at the next cancellation point.
    Cancel,
}
