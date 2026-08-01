use crate::{FactorConfig, FactorError, Natural, PARTS, PrimeFactors, ProgressAction};
use std::time::Instant;

/// Factors a positive integer with the default configuration.
///
/// The returned factors are sorted and grouped by multiplicity. The
/// factorization of one is empty; zero returns
/// [`FactorError::ZeroHasNoPrimeFactorization`].
/// Inputs wider than 512 bits return [`FactorError::InputTooLarge`].
///
/// This function is not constant-time and must not process secret values.
///
/// # Errors
///
/// Returns [`FactorError`] for zero, oversized input, cancellation, resource
/// exhaustion, or an internal factorization failure.
///
/// # Examples
///
/// ```
/// use rusqsieve::{Natural, factor};
///
/// let factors = factor(Natural::<4>::from_u64(360))?;
/// assert_eq!(factors.total_len(), 6);
/// # Ok::<(), rusqsieve::FactorError>(())
/// ```
pub fn factor<const P: usize>(input: Natural<P>) -> Result<PrimeFactors<P>, FactorError> {
    factor_with(input, FactorConfig::default())
}

/// Factors a positive integer with an explicit configuration.
///
/// This function is not constant-time and must not process secret values.
///
/// # Errors
///
/// Returns [`FactorError`] for zero, oversized input, cancellation, resource
/// exhaustion, or an internal factorization failure.
pub fn factor_with<const P: usize>(
    input: Natural<P>,
    config: FactorConfig,
) -> Result<PrimeFactors<P>, FactorError> {
    factor_impl::<P, false, _>(input, config, |_| ProgressAction::Continue)
}

/// Factors a positive integer while reporting cooperative progress.
///
/// Returning [`ProgressAction::Cancel`] from `observer` stops the operation at
/// the next cancellation point. Callbacks execute on the calling thread and
/// should return quickly.
///
/// This function is not constant-time and must not process secret values.
///
/// # Errors
///
/// Returns [`FactorError`] for zero, oversized input, callback cancellation,
/// resource exhaustion, or an internal factorization failure.
pub fn factor_with_progress<const P: usize, F>(
    input: Natural<P>,
    config: FactorConfig,
    observer: F,
) -> Result<PrimeFactors<P>, FactorError>
where
    F: FnMut(&crate::ProgressSnapshot) -> ProgressAction,
{
    factor_impl::<P, true, _>(input, config, observer)
}

fn factor_impl<const P: usize, const REPORT_PROGRESS: bool, F>(
    input: Natural<P>,
    config: FactorConfig,
    mut observer: F,
) -> Result<PrimeFactors<P>, FactorError>
where
    F: FnMut(&crate::ProgressSnapshot) -> ProgressAction,
{
    if input.is_zero() {
        return Err(FactorError::ZeroHasNoPrimeFactorization);
    }
    // No blanket width cap. The sieve's range limit is enforced against the composite that
    // actually reaches it (see `engine::MAX_SIQS_BITS`), so a wide input built from small factors
    // still factors; only a wide *hard* cofactor is refused. Capacity is the sole limit here.
    let mut bytes = vec![0u8; P * 8];
    let written = input
        .write_le_bytes(&mut bytes)
        .map_err(|_| FactorError::CapacityExceeded)?;
    let fast_input =
        Natural::from_le_bytes(&bytes[..written]).map_err(|_| FactorError::CapacityExceeded)?;
    const MAX_WORKERS: usize = 256;
    let workers = match config.parallelism() {
        crate::Parallelism::Auto => std::thread::available_parallelism().map_or(1, usize::from),
        crate::Parallelism::Threads(n) => n.get(),
    }
    .min(MAX_WORKERS);
    let input_bits = input.bit_len();
    let mut revision = 1u64;
    if REPORT_PROGRESS {
        let initial = crate::ProgressSnapshot {
            revision,
            input_bits,
            phase: crate::ProgressPhase::Preprocessing,
            amount: crate::ProgressAmount {
                completed: 0,
                total: crate::ProgressTotal::Unknown,
                unit: crate::ProgressUnit::Tasks,
            },
        };
        if observer(&initial) == ProgressAction::Cancel {
            return Err(FactorError::Cancelled);
        }
    }
    let interval = config.progress_interval();
    let mut last_report = REPORT_PROGRESS.then(Instant::now);
    let mut last_phase = Some(crate::ProgressPhase::Preprocessing);
    let fast = crate::engine::factor(
        fast_input,
        workers,
        &config.tuning,
        config.witness_seed,
        |state| {
            if !REPORT_PROGRESS {
                return true;
            }
            let phase = match state.phase {
                crate::engine::EnginePhase::Preprocessing => crate::ProgressPhase::Preprocessing,
                crate::engine::EnginePhase::BuildingFactorBase => {
                    crate::ProgressPhase::BuildingFactorBase
                }
                crate::engine::EnginePhase::Sieving => crate::ProgressPhase::Sieving,
                crate::engine::EnginePhase::LinearAlgebra => crate::ProgressPhase::LinearAlgebra,
                crate::engine::EnginePhase::Extracting => crate::ProgressPhase::ExtractingFactor,
            };
            if last_phase == Some(phase)
                && last_report
                    .as_ref()
                    .is_some_and(|last| last.elapsed() < interval)
            {
                return true;
            }
            revision += 1;
            let snapshot = crate::ProgressSnapshot {
                revision,
                input_bits,
                phase,
                amount: crate::ProgressAmount {
                    completed: state.relations as u64,
                    total: if state.target == 0 {
                        crate::ProgressTotal::Unknown
                    } else {
                        crate::ProgressTotal::Estimated(state.target as u64)
                    },
                    unit: crate::ProgressUnit::Relations,
                },
            };
            let keep_going = observer(&snapshot) == ProgressAction::Continue;
            last_phase = Some(phase);
            last_report = Some(Instant::now());
            keep_going
        },
    )
    .map_err(|error| match error {
        crate::engine::EngineError::Cancelled => FactorError::Cancelled,
        crate::engine::EngineError::ResourceLimit => {
            FactorError::ResourceLimit(crate::ResourceLimitKind::Memory)
        }
        crate::engine::EngineError::SiqsInputTooLarge(bits) => {
            FactorError::SiqsCompositeTooLarge(bits)
        }
        crate::engine::EngineError::InsufficientRelations => FactorError::InsufficientRelations,
        crate::engine::EngineError::NoDependency => FactorError::NoDependency,
        crate::engine::EngineError::PolynomialSelection(message) => {
            FactorError::PolynomialSelection(message)
        }
        crate::engine::EngineError::InvalidDependency => FactorError::InvalidDependency,
        crate::engine::EngineError::Worker(message) => FactorError::WorkerFailure(message),
        crate::engine::EngineError::NoFactor => FactorError::NoNontrivialFactor,
        crate::engine::EngineError::Setup(_) => FactorError::InternalFailure,
    })?;
    let mut output = PrimeFactors::new();
    for value in fast {
        let mut raw = [0u8; PARTS * 8];
        let n = value
            .write_le_bytes(&mut raw)
            .map_err(|_| FactorError::CapacityExceeded)?;
        let converted =
            Natural::<P>::from_le_bytes(&raw[..n]).map_err(|_| FactorError::CapacityExceeded)?;
        output.insert_count(converted, 1)
    }
    if REPORT_PROGRESS {
        revision += 1;
        let complete = crate::ProgressSnapshot {
            revision,
            input_bits,
            phase: crate::ProgressPhase::Complete,
            amount: crate::ProgressAmount {
                completed: 1,
                total: crate::ProgressTotal::Exact(1),
                unit: crate::ProgressUnit::Tasks,
            },
        };
        let _ = observer(&complete);
    }
    Ok(output)
}
