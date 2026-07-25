#![cfg(any(unix, windows))]
use rusqsieve::{
    FactorConfig, FactorError, Natural, Parallelism, ProgressAction, ProgressPhase, factor,
    factor_with_progress,
};

#[test]
fn factors_sorted_with_repetitions() {
    let n: Natural = Natural::from_u64(360);
    let factors = factor(n.clone()).unwrap();
    assert_eq!(
        factors
            .iter()
            .map(|(p, e)| (p.to_string(), e.get()))
            .collect::<Vec<_>>(),
        [("2".into(), 3), ("3".into(), 2), ("5".into(), 1)]
    );
    assert!(factors.verify_product(&n));
    assert_eq!(factors.distinct_len(), 3);
    assert_eq!(factors.total_len(), 6);
    assert_eq!(
        factors
            .expanded()
            .map(ToString::to_string)
            .collect::<Vec<_>>(),
        ["2", "2", "2", "3", "3", "5"]
    );
}

#[test]
fn zero_and_one() {
    let zero: Natural = Natural::ZERO;
    let one: Natural = Natural::ONE;
    assert!(matches!(
        factor(zero),
        Err(FactorError::ZeroHasNoPrimeFactorization)
    ));
    assert!(factor(one).unwrap().is_empty());
}

#[test]
fn configuration_is_validated_and_encapsulated() {
    assert_eq!(Parallelism::threads(0), None);
    let threads = Parallelism::threads(2).unwrap();
    let config = FactorConfig::default().with_parallelism(threads);
    assert_eq!(config.parallelism(), threads);
}

#[test]
fn progress_can_cancel_before_work_starts() {
    let input: Natural = Natural::from_u64(360);
    let result = factor_with_progress(input, FactorConfig::default(), |snapshot| {
        assert_eq!(snapshot.phase(), ProgressPhase::Preprocessing);
        ProgressAction::Cancel
    });
    assert!(matches!(result, Err(FactorError::Cancelled)));
}

#[test]
fn progress_finishes_with_a_complete_snapshot() {
    let mut phases = Vec::new();
    let input: Natural = Natural::from_u64(360);
    let factors = factor_with_progress(input, FactorConfig::default(), |snapshot| {
        phases.push(snapshot.phase());
        ProgressAction::Continue
    })
    .unwrap();
    assert_eq!(factors.total_len(), 6);
    assert_eq!(phases.first(), Some(&ProgressPhase::Preprocessing));
    assert_eq!(phases.last(), Some(&ProgressPhase::Complete));
}

#[test]
fn progress_cancellation_stops_parallel_sieving() {
    let p: Natural = Natural::from_u64(18_446_744_073_709_551_557);
    let q: Natural = Natural::from_u64(18_446_744_073_709_551_533);
    let input = p.checked_mul(&q).unwrap();
    let mut reached_sieving = false;
    let result = factor_with_progress(input, FactorConfig::default(), |snapshot| {
        if snapshot.phase() == ProgressPhase::Sieving {
            reached_sieving = true;
            ProgressAction::Cancel
        } else {
            ProgressAction::Continue
        }
    });
    assert!(reached_sieving);
    assert!(matches!(result, Err(FactorError::Cancelled)));
}
