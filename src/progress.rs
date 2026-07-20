//! Scheduler-neutral progress snapshots.

use core::time::Duration;

#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct ProgressSnapshot {
    pub revision: u64,
    pub task_id: u64,
    pub input_bits: usize,
    pub phase: ProgressPhase,
    pub amount: ProgressAmount,
    pub detail: ProgressDetail,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProgressAmount {
    pub completed: u64,
    pub total: ProgressTotal,
    pub unit: ProgressUnit,
}
impl ProgressAmount {
    pub fn fraction(self) -> Option<f64> {
        match self.total {
            ProgressTotal::Exact(0) | ProgressTotal::Estimated(0) | ProgressTotal::Unknown => None,
            ProgressTotal::Exact(n) | ProgressTotal::Estimated(n) => {
                Some(self.completed as f64 / n as f64)
            }
        }
    }
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressTotal {
    Exact(u64),
    Estimated(u64),
    Unknown,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressUnit {
    Candidates,
    Primes,
    Polynomials,
    SievePositions,
    Relations,
    MatrixRows,
    MatrixColumns,
    MatrixNonzeros,
    Iterations,
    MatrixProducts,
    Tasks,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ProgressPhase {
    Preprocessing,
    BuildingFactorBase,
    Sieving,
    CombiningRelations,
    BuildingMatrix,
    FilteringMatrix,
    LinearAlgebra,
    ExtractingFactor,
    PrimalityTesting,
    Complete,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct FactorBaseProgress {
    pub bound: u32,
    pub searched_through: u32,
    pub primes_tested: u64,
    pub primes_accepted: u64,
    pub nonresidue_primes: u64,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct SievingProgress {
    pub polynomial_families_completed: u64,
    pub polynomials_completed: u64,
    pub sieve_positions_processed: u64,
    pub candidates_tested: u64,
    pub full_relations: u64,
    pub single_large_prime_relations: u64,
    pub double_large_prime_relations: u64,
    pub usable_relations: u64,
    pub target_relations: u64,
    pub active_workers: usize,
    pub outstanding_jobs: usize,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct RelationProgress {
    pub partial_relations_examined: u64,
    pub partial_relations_total: u64,
    pub graph_vertices: u64,
    pub graph_edges: u64,
    pub cycles_found: u64,
    pub combined_relations: u64,
}
#[derive(Clone, Copy, Debug)]
pub struct MatrixProgress {
    pub stage: MatrixStage,
    pub original_rows: u64,
    pub original_columns: u64,
    pub original_nonzeros: u64,
    pub current_rows: u64,
    pub current_columns: u64,
    pub current_nonzeros: u64,
    pub items_processed: u64,
    pub items_total: Option<u64>,
    pub singleton_rows_removed: u64,
    pub duplicate_columns_removed: u64,
    pub structured_eliminations: u64,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatrixStage {
    CountingParities,
    ConstructingCsr,
    ConstructingCsc,
    RemovingSingletons,
    RemovingDuplicates,
    StructuredElimination,
    Finalizing,
}
#[derive(Clone, Copy, Debug)]
pub struct LinearAlgebraProgress {
    pub solver: LinearAlgebraSolver,
    pub stage: LinearAlgebraStage,
    pub matrix_rows: u64,
    pub matrix_columns: u64,
    pub matrix_nonzeros: u64,
    pub iteration: u64,
    pub estimated_iterations: Option<u64>,
    pub matrix_products_completed: u64,
    pub matrix_products_estimated: Option<u64>,
    pub dependencies_found: u32,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LinearAlgebraSolver {
    DenseGaussian,
    BlockLanczos,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LinearAlgebraStage {
    Initializing,
    Iterating,
    RecoveringDependencies,
    VerifyingDependencies,
    Complete,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct CompleteProgress {
    pub prime_factors: u64,
    pub factors_with_multiplicity: u64,
}
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum ProgressDetail {
    None,
    FactorBase(FactorBaseProgress),
    Sieving(SievingProgress),
    Relations(RelationProgress),
    Matrix(MatrixProgress),
    LinearAlgebra(LinearAlgebraProgress),
    Complete(CompleteProgress),
}
#[derive(Clone, Copy, Debug)]
pub struct ProgressReportingConfig {
    pub minimum_interval: Duration,
}
impl Default for ProgressReportingConfig {
    fn default() -> Self {
        Self {
            minimum_interval: Duration::from_millis(100),
        }
    }
}
impl ProgressSnapshot {
    pub(crate) fn initial(bits: usize) -> Self {
        Self {
            revision: 0,
            task_id: 0,
            input_bits: bits,
            phase: ProgressPhase::Preprocessing,
            amount: ProgressAmount {
                completed: 0,
                total: ProgressTotal::Unknown,
                unit: ProgressUnit::Tasks,
            },
            detail: ProgressDetail::None,
        }
    }
}
