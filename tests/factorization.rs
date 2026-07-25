#![cfg(any(unix, windows))]
use rusqsieve::{
    FactorConfig, FactorError, Natural, Parallelism, ProgressAction, ProgressPhase, factor,
    factor_with_progress,
};
use std::str::FromStr;

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

fn assert_factorization(input: &str, expected: &[&str]) {
    let n = Natural::<16>::from_str(input).unwrap();
    let factors = factor(n.clone()).unwrap_or_else(|error| panic!("failed to factor {n}: {error}"));
    assert!(
        factors.verify_product(&n),
        "factor product mismatch for {n}"
    );
    let actual = factors
        .expanded()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    assert_eq!(actual, expected, "prime factorization mismatch for {n}");

    // Each listed factor must itself terminate as one prime factor. This exercises
    // the library's final primality decision rather than merely trusting the product.
    for expected_prime in expected {
        let prime = Natural::<16>::from_str(expected_prime).unwrap();
        let prime_factors = factor(prime.clone()).unwrap();
        assert_eq!(
            prime_factors.total_len(),
            1,
            "{prime} was not accepted as prime"
        );
        assert!(prime_factors.verify_product(&prime));
    }
}

#[test]
fn balanced_semiprimes_across_the_choose_a_dead_zone() {
    let cases: &[(&str, &[&str])] = &[
        // Required minimum reproducer (65 bits).
        ("18446744400127067027", &["4294967311", "4294967357"]),
        ("27072011721716628587", &["4273765633", "6334463339"]),
        ("635904368119925963561", &["19860882047", "32017931863"]),
        ("20988451891514649258347", &["91807517561", "228613652227"]),
        (
            "703713894016303629914563",
            &["815250411389", "863187413567"],
        ),
        (
            "22921914745054882120472087",
            &["2921865453193", "7844959020959"],
        ),
        (
            "648536833001811612107041493",
            &["24295067057461", "26694177524513"],
        ),
    ];
    for &(n, expected) in cases {
        assert_factorization(n, expected);
    }
}

#[test]
fn supplied_factorization_corpus() {
    let corpus = include_str!("data/rusqsieve-factorization-corpus.txt");
    let mut tested = 0usize;
    for (line_index, line) in corpus.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split_ascii_whitespace().collect::<Vec<_>>();
        assert!(
            fields.len() >= 2,
            "corpus line {} has no factorization",
            line_index + 1
        );
        assert_factorization(fields[0], &fields[1..]);
        tested += 1;
    }
    assert_eq!(tested, 309, "unexpected corpus entry count");
}

#[test]
#[ignore = "manual interleaved performance measurement"]
fn profile_single_input() {
    let input = std::env::var("RUSQSIEVE_BENCH_INPUT").unwrap();
    let n = Natural::<16>::from_str(&input).unwrap();
    let repeats = std::env::var("RUSQSIEVE_BENCH_REPEATS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(1);
    let started = std::time::Instant::now();
    for _ in 0..repeats {
        let factors = rusqsieve::factor_with(
            n.clone(),
            FactorConfig::default().with_parallelism(Parallelism::threads(4).unwrap()),
        )
        .unwrap();
        assert!(factors.verify_product(&n));
    }
    eprintln!(
        "BENCH input_bits={} repeats={} elapsed={:.6}s",
        n.bit_len(),
        repeats,
        started.elapsed().as_secs_f64()
    );
}
