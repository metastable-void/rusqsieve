//! Deterministic bounded work packets shared by all schedulers.
use crate::f2::{MatrixOperation, SparseBinaryMatrix};
use crate::qs::{RawRelation, SieveContext, SieveScratch};
use core::fmt;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JobHeader {
    pub job_id: u64,
    pub generation: u64,
    pub context_id: u64,
}
#[derive(Clone, Debug)]
pub enum WorkJob {
    Sieve(SieveJob),
    MatrixMultiply(MatrixMultiplyJob),
}
#[derive(Clone, Debug)]
pub struct SieveJob {
    pub header: JobHeader,
    pub family: u64,
    pub first_polynomial: u32,
    pub polynomial_count: u32,
}
#[derive(Clone, Debug)]
pub struct MatrixMultiplyJob {
    pub header: JobHeader,
    pub operation: MatrixOperation,
    pub output_start: u32,
    pub output_end: u32,
    pub input: Box<[u64]>,
}
#[derive(Debug)]
pub enum WorkResult<const P: usize = 16> {
    Sieve(SieveResult<P>),
    MatrixMultiply(MatrixMultiplyResult),
    Failed(WorkFailure),
}
#[derive(Clone, Debug)]
pub struct WorkFailure {
    pub header: JobHeader,
    pub message: String,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct SieveJobMetrics {
    pub polynomial_families: u64,
    pub polynomials: u64,
    pub sieve_positions: u64,
    pub candidates_tested: u64,
    pub full_relations: u64,
    pub single_large_prime_relations: u64,
    pub double_large_prime_relations: u64,
}
#[derive(Debug)]
pub struct SieveResult<const P: usize = 16> {
    pub header: JobHeader,
    pub relations: Vec<RawRelation<P>>,
    pub metrics: SieveJobMetrics,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct MatrixJobMetrics {
    pub output_items: u64,
    pub nonzeros_visited: u64,
}
#[derive(Debug)]
pub struct MatrixMultiplyResult {
    pub header: JobHeader,
    pub output_start: u32,
    pub words: Box<[u64]>,
    pub metrics: MatrixJobMetrics,
}
#[derive(Default)]
pub struct MatrixScratch {
    pub words: Vec<u64>,
}
#[derive(Default)]
pub struct ArithmeticScratch {
    pub words: Vec<u64>,
}
#[derive(Default)]
pub struct WorkerScratch {
    pub sieve: SieveScratch,
    pub matrix: MatrixScratch,
    pub arithmetic: ArithmeticScratch,
}
pub struct KernelContexts<const P: usize> {
    pub sieve: Option<SieveContext<P>>,
    pub matrix: Option<SparseBinaryMatrix>,
}
impl<const P: usize> Default for KernelContexts<P> {
    fn default() -> Self {
        Self {
            sieve: None,
            matrix: None,
        }
    }
}
#[derive(Clone, Debug)]
pub enum KernelError {
    MissingContext,
    ContextMismatch,
    Range,
    Worker(String),
}
impl fmt::Display for KernelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "worker kernel error: {self:?}")
    }
}
impl std::error::Error for KernelError {}
pub fn execute_job<const P: usize>(
    contexts: &KernelContexts<P>,
    job: WorkJob,
    scratch: &mut WorkerScratch,
) -> Result<WorkResult<P>, KernelError> {
    match job {
        WorkJob::Sieve(job) => {
            let ctx = contexts.sieve.as_ref().ok_or(KernelError::MissingContext)?;
            if job.header.context_id != ctx.context_id() {
                return Err(KernelError::ContextMismatch);
            }
            Ok(WorkResult::Sieve(crate::qs::sieve_job(
                ctx,
                &job,
                &mut scratch.sieve,
            )))
        }
        WorkJob::MatrixMultiply(job) => {
            let matrix = contexts
                .matrix
                .as_ref()
                .ok_or(KernelError::MissingContext)?;
            let start = job.output_start as usize;
            let end = job.output_end as usize;
            if start > end {
                return Err(KernelError::Range);
            }
            let mut words = vec![0; end - start];
            match job.operation {
                MatrixOperation::Matrix => matrix.mul_m_rows(&job.input, start..end, &mut words),
                MatrixOperation::Transpose => {
                    matrix.mul_mt_columns(&job.input, start..end, &mut words)
                }
            }
            Ok(WorkResult::MatrixMultiply(MatrixMultiplyResult {
                header: job.header,
                output_start: job.output_start,
                words: words.into_boxed_slice(),
                metrics: MatrixJobMetrics {
                    output_items: (end - start) as u64,
                    nonzeros_visited: 0,
                },
            }))
        }
    }
}
