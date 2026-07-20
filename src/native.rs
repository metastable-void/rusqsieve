use crate::{
    FactorConfig, FactorError, FactorSession, LocalWorkBudget, Natural, PrimeFactors,
    ProgressAction,
};
use std::time::Instant;

pub fn factor<const P: usize>(input: Natural<P>) -> Result<PrimeFactors<P>, FactorError> {
    factor_with(input, FactorConfig::default())
}
pub fn factor_with<const P: usize>(
    input: Natural<P>,
    config: FactorConfig,
) -> Result<PrimeFactors<P>, FactorError> {
    factor_with_progress(input, config, |_| ProgressAction::Continue)
}
pub fn factor_with_progress<const P: usize, F>(
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
    if P <= 16 {
        let mut bytes = vec![0u8; P * 8];
        let written = input
            .write_le_bytes(&mut bytes)
            .map_err(|_| FactorError::CapacityExceeded)?;
        let fast_input = Natural::<16>::from_le_bytes(&bytes[..written])
            .map_err(|_| FactorError::CapacityExceeded)?;
        let workers = match config.parallelism {
            crate::Parallelism::Auto => std::thread::available_parallelism().map_or(1, usize::from),
            crate::Parallelism::Exact(n) => n.get(),
        };
        let mut revision = 0u64;
        let fast = crate::engine::factor(fast_input, workers, |state| {
            revision += 1;
            let phase = match state.phase {
                crate::engine::EnginePhase::Preprocessing => crate::ProgressPhase::Preprocessing,
                crate::engine::EnginePhase::BuildingFactorBase => {
                    crate::ProgressPhase::BuildingFactorBase
                }
                crate::engine::EnginePhase::Sieving => crate::ProgressPhase::Sieving,
                crate::engine::EnginePhase::LinearAlgebra => crate::ProgressPhase::LinearAlgebra,
                crate::engine::EnginePhase::Extracting => crate::ProgressPhase::ExtractingFactor,
            };
            let snapshot = crate::ProgressSnapshot {
                revision,
                task_id: 0,
                input_bits: input.bit_len(),
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
                detail: crate::ProgressDetail::None,
            };
            let _ = observer(&snapshot);
        })
        .map_err(|_| FactorError::NoNontrivialFactor)?;
        let mut output = PrimeFactors::new();
        for value in fast {
            let mut raw = [0u8; 128];
            let n = value
                .write_le_bytes(&mut raw)
                .map_err(|_| FactorError::CapacityExceeded)?;
            let converted = Natural::<P>::from_le_bytes(&raw[..n])
                .map_err(|_| FactorError::CapacityExceeded)?;
            output.insert_count(converted, 1)
        }
        return Ok(output);
    }
    let interval = config.progress_reporting.minimum_interval;
    let mut session = FactorSession::new(input, config)?;
    if observer(session.progress()) == ProgressAction::Cancel {
        return Err(FactorError::Cancelled);
    }
    let mut last = Instant::now();
    loop {
        let outcome = session.advance_local(LocalWorkBudget::default())?;
        if last.elapsed() >= interval || matches!(outcome, crate::AdvanceOutcome::Complete) {
            if observer(session.progress()) == ProgressAction::Cancel {
                return Err(FactorError::Cancelled);
            }
            last = Instant::now()
        }
        if matches!(outcome, crate::AdvanceOutcome::Complete) {
            return session.take_factors();
        }
    }
}
