//! Progress reporting for blocking factorization.

use core::time::Duration;

/// A point-in-time view of a factorization.
///
/// Snapshots are passed to the callback supplied to
/// [`factor_with_progress`](crate::factor_with_progress). Values are
/// monotonically revised, but the amount and its unit can change when the
/// factorization enters a new phase.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct ProgressSnapshot {
    pub(crate) revision: u64,
    pub(crate) input_bits: usize,
    pub(crate) phase: ProgressPhase,
    pub(crate) amount: ProgressAmount,
}

impl ProgressSnapshot {
    /// Returns the monotonically increasing revision number.
    pub const fn revision(&self) -> u64 {
        self.revision
    }

    /// Returns the significant width of the original input.
    pub const fn input_bits(&self) -> usize {
        self.input_bits
    }

    /// Returns the current factorization phase.
    pub const fn phase(&self) -> ProgressPhase {
        self.phase
    }

    /// Returns phase-specific completed and total work.
    pub const fn amount(&self) -> ProgressAmount {
        self.amount
    }

    pub(crate) fn initial(bits: usize) -> Self {
        Self {
            revision: 0,
            input_bits: bits,
            phase: ProgressPhase::Preprocessing,
            amount: ProgressAmount {
                completed: 0,
                total: ProgressTotal::Unknown,
                unit: ProgressUnit::Tasks,
            },
        }
    }
}

/// Phase-specific progress counters.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProgressAmount {
    pub(crate) completed: u64,
    pub(crate) total: ProgressTotal,
    pub(crate) unit: ProgressUnit,
}

impl ProgressAmount {
    /// Returns completed work in [`Self::unit`].
    pub const fn completed(self) -> u64 {
        self.completed
    }

    /// Returns the known or estimated total.
    pub const fn total(self) -> ProgressTotal {
        self.total
    }

    /// Returns the unit used by this amount.
    pub const fn unit(self) -> ProgressUnit {
        self.unit
    }

    /// Returns a completion fraction when a nonzero total is available.
    ///
    /// Estimated totals can change, so this value is informational and is not
    /// guaranteed to increase monotonically.
    pub fn fraction(self) -> Option<f64> {
        match self.total {
            ProgressTotal::Exact(0) | ProgressTotal::Estimated(0) | ProgressTotal::Unknown => None,
            ProgressTotal::Exact(n) | ProgressTotal::Estimated(n) => {
                Some(self.completed as f64 / n as f64)
            }
        }
    }
}

/// Availability and confidence of a progress total.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressTotal {
    /// An exact total.
    Exact(u64),
    /// A total that may be revised as work proceeds.
    Estimated(u64),
    /// No useful total is currently known.
    Unknown,
}

/// Unit used by a [`ProgressAmount`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ProgressUnit {
    /// Candidate values examined.
    Candidates,
    /// Prime numbers examined or accepted.
    Primes,
    /// SIQS polynomials processed.
    Polynomials,
    /// Sieve positions processed.
    SievePositions,
    /// Relations collected.
    Relations,
    /// Matrix rows processed.
    MatrixRows,
    /// Matrix columns processed.
    MatrixColumns,
    /// Matrix nonzero entries processed.
    MatrixNonzeros,
    /// Algorithm iterations completed.
    Iterations,
    /// Matrix-vector products completed.
    MatrixProducts,
    /// Coarse-grained tasks completed.
    Tasks,
}

/// Current high-level stage of factorization.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ProgressPhase {
    /// Validating and reducing the input.
    Preprocessing,
    /// Constructing the SIQS factor base.
    BuildingFactorBase,
    /// Collecting smooth and partial relations.
    Sieving,
    /// Combining partial relations.
    CombiningRelations,
    /// Constructing the sparse parity matrix.
    BuildingMatrix,
    /// Reducing the sparse matrix before solving.
    FilteringMatrix,
    /// Solving for linear dependencies.
    LinearAlgebra,
    /// Recovering a nontrivial divisor from dependencies.
    ExtractingFactor,
    /// Testing candidate factors for probable primality.
    PrimalityTesting,
    /// Factorization completed successfully.
    Complete,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct FactorBaseProgress {
    pub(crate) bound: u32,
    pub(crate) searched_through: u32,
    pub(crate) primes_tested: u64,
    pub(crate) primes_accepted: u64,
    pub(crate) nonresidue_primes: u64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ProgressReportingConfig {
    pub(crate) minimum_interval: Duration,
}

impl Default for ProgressReportingConfig {
    fn default() -> Self {
        Self {
            minimum_interval: Duration::from_millis(100),
        }
    }
}
