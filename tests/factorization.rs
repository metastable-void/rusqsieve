#![cfg(any(unix, windows))]
use rusqsieve::{FactorError, Natural, factor};

#[test]
fn factors_sorted_with_repetitions() {
    let n = Natural::<16>::from_u64(360);
    let factors = factor(n.clone()).unwrap();
    assert_eq!(
        factors
            .iter()
            .map(|(p, e)| (p.to_string(), e.get()))
            .collect::<Vec<_>>(),
        [("2".into(), 3), ("3".into(), 2), ("5".into(), 1)]
    );
    assert!(factors.verify_product(&n));
}

#[test]
fn zero_and_one() {
    assert!(matches!(
        factor(Natural::<16>::ZERO),
        Err(FactorError::ZeroHasNoPrimeFactorization)
    ));
    assert!(factor(Natural::<16>::ONE).unwrap().is_empty());
}
