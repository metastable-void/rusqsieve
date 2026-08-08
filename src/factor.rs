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
    /// Overrides the computed Pollard-Brent iteration budget at every width. Its reason to exist is
    /// the range above the sieve ceiling, where rho is the whole attempt and the default stops at
    /// about a minute: a caller who wants a 56- or 64-bit factor out of a 500-bit input spends the
    /// minutes or hours here deliberately.
    pub(crate) rho_iterations: Option<u64>,
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
    pub(crate) ecm: bool,
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

    /// Whether the elliptic curve method may run on composites the sieve can handle.
    #[must_use]
    pub const fn ecm(&self) -> bool {
        self.ecm
    }

    /// Enables the elliptic curve method for composites within the sieve's range.
    ///
    /// ECM finds a medium-size factor — 20 to 30 digits — in time governed by the size of that
    /// factor rather than of the input, which is the one shape neither Pollard–Brent nor the
    /// quadratic sieve handles well. It is off by default because this crate's workload is balanced
    /// semiprimes, where a curve cannot succeed and is pure overhead ahead of the sieve.
    ///
    /// This switch governs a narrower range than it looks like. Curves already run without it in
    /// two cases, and neither can arise for a balanced semiprime:
    ///
    /// - a composite wider than the sieve accepts, where there is no fallback and the alternative
    ///   is [`FactorError::SiqsCompositeTooLarge`] on a number whose factor ECM would have found;
    /// - a composite already known to be unbalanced, because trial division peeled a factor or
    ///   Pollard–Brent split an ancestor. Such a number has a small factor, so it may well have a
    ///   medium one, and the sieve would charge for the size of the input instead.
    ///
    /// What is left — inside the sieve's range with no evidence either way — is where a balanced
    /// semiprime lives, and that is what this switch turns on.
    #[must_use]
    pub fn with_ecm(mut self, enabled: bool) -> Self {
        self.ecm = enabled;
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
        rho_iterations: Option<u64>,
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
            rho_iterations,
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
            ecm: false,
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
    /// Inputs wider than the selected [`Natural`](crate::Natural) capacity are rejected.
    InputTooLarge,
    /// A composite requiring the quadratic sieve exceeded the supported 400-bit range.
    ///
    /// This bounds the hard cofactor, not the caller's input: a number of any width whose
    /// factors are small enough to peel off by trial division, perfect-power detection, or
    /// Pollard–Brent never reaches the sieve and never produces this error. The payload is the
    /// bit length of the composite that did.
    ///
    /// Because such a composite has nowhere else to go, two searches run before this is
    /// returned: a deeper Pollard–Brent budget than the sieve's range gets, reaching a smallest
    /// factor near 2^46 at every width, and then a committed run of the elliptic curve method,
    /// which reaches 25-digit factors and costs by the size of the factor rather than of the
    /// input. The error therefore means the smallest factor outran both, not that the input was
    /// too wide to attempt. `RUSQSIEVE_RHO_ITERATIONS` overrides the rho budget.
    SiqsCompositeTooLarge(usize),
    /// A configured internal resource limit was exceeded.
    ResourceLimit(ResourceLimitKind),
    /// No nontrivial divisor was recovered.
    NoNontrivialFactor,
    /// The linear algebra found no nontrivial dependency in the collected relations.
    NoDependency,
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
            Self::InputTooLarge => f.write_str("input exceeds the selected integer capacity"),
            Self::SiqsCompositeTooLarge(bits) => write!(
                f,
                "a {bits}-bit composite requires the quadratic sieve, which supports at most \
                 {} bits; inputs of any width factor normally when their factors are small",
                crate::engine::MAX_SIQS_BITS
            ),
            Self::ResourceLimit(kind) => write!(f, "resource limit exceeded: {kind:?}"),
            Self::NoNontrivialFactor => f.write_str("no nontrivial factor found"),
            Self::NoDependency => {
                f.write_str("linear algebra found no nontrivial dependency; collect more relations")
            }
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
