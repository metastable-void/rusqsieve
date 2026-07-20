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
