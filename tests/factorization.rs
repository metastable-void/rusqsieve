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
    for prime in [2u64, 3] {
        let input: Natural = Natural::from_u64(prime);
        let factors = factor(input.clone()).unwrap();
        assert_eq!(factors.total_len(), 1);
        assert!(factors.verify_product(&input));
    }
}

#[test]
fn widths_above_the_engine_capacity_use_the_real_siqs_path() {
    let input =
        Natural::<17>::from_decimal("18446744400127067027").expect("65-bit value fits 17 limbs");
    let factors = factor(input.clone()).unwrap();
    assert!(factors.verify_product(&input));
    assert_eq!(factors.distinct_len(), 2);
}

#[test]
fn perfect_power_and_recursive_multifactor_inputs() {
    let prime = Natural::<16>::from_decimal("170141183460469231731687303715884105727").unwrap();
    let square = prime.checked_mul(&prime).unwrap();
    let factors = factor(square.clone()).unwrap();
    assert!(factors.verify_product(&square));
    assert_eq!(factors.multiplicity(&prime).unwrap().get(), 2);

    let three_prime = Natural::<16>::from_decimal("10007000160112000630441").unwrap();
    // 10007 × 1000000007 × 1000000009
    let factors = factor(three_prime.clone()).unwrap();
    assert!(factors.verify_product(&three_prime));
    assert_eq!(factors.total_len(), 3);

    let odd = Natural::<16>::from_u64(1_000_003);
    let power_of_two_input = odd.clone() << 180;
    let factors = factor(power_of_two_input.clone()).unwrap();
    assert!(factors.verify_product(&power_of_two_input));
    assert_eq!(
        factors.multiplicity(&Natural::from_u64(2)).unwrap().get(),
        180
    );
    assert_eq!(factors.multiplicity(&odd).unwrap().get(), 1);
}

/// Width alone must not disqualify an input. Everything here is far wider than the quadratic
/// sieve's 400-bit range, and every one of them factors completely, because the factors are
/// reachable by trial division and Pollard-Brent and the sieve is never consulted.
///
/// The engine used to reject anything over 512 bits up front, on the caller's input width, which
/// refused work it could do trivially. The range limit now applies only to a composite that
/// actually reaches SIQS.
#[test]
fn wide_inputs_with_rho_reachable_factors_factor_completely() {
    // 506 bits, thirty-two distinct 16-bit primes. Sized so the whole ladder stays in Pollard-
    // Brent's reach: `rho_budget` shrinks with the cofactor, and factors much larger than this
    // outrun it partway down and fall through to the sieve. Verified to enter SIQS zero times.
    let expected: [&str; 32] = [
        "50119", "50231", "51239", "52121", "52223", "52837", "52967", "53527", "53549", "53591",
        "53917", "54331", "54347", "55001", "55291", "55763", "55787", "55949", "56359", "57587",
        "57689", "58897", "59183", "60427", "61561", "61717", "62653", "64591", "64747", "64877",
        "65011", "65371",
    ];
    let input = Natural::<16>::from_decimal(
        "1376775833903088655333039822117419079135686880866602737071940539886253922\
         79103037177076957972313736914320822855524800052427782940879523556486820043373181",
    )
    .unwrap();
    assert_eq!(input.bit_len(), 506);
    let factors = factor(input.clone()).unwrap();
    assert!(factors.verify_product(&input));
    assert_eq!(factors.total_len(), 32);
    assert_eq!(
        factors
            .expanded()
            .map(ToString::to_string)
            .collect::<Vec<_>>(),
        expected
    );

    // 511 bits, but the hard part is a 320-bit *prime*: the six small factors peel off and the
    // cofactor is settled by the primality test, so nothing is ever handed to the sieve. A
    // width-based gate would have rejected this outright.
    let mixed = Natural::<16>::from_decimal(
        "5095761221610806052403790675953969613924327896615930674429444570224327135\
         255024776878949868364526366809455906988948729433090635881939673406756856142258511",
    )
    .unwrap();
    assert_eq!(mixed.bit_len(), 511);
    let factors = factor(mixed.clone()).unwrap();
    assert!(factors.verify_product(&mixed));
    assert_eq!(factors.total_len(), 7);
    let big = Natural::<16>::from_decimal(
        "1822715629003050836886041093405870953882958356267076001214210352068467068\
         972779348077648099904357",
    )
    .unwrap();
    assert_eq!(big.bit_len(), 320);
    assert_eq!(factors.multiplicity(&big).unwrap().get(), 1);

    // A power-of-two shift is the degenerate case of the same rule: arbitrarily wide, trivially
    // factored, and nowhere near the sieve.
    let odd = Natural::<16>::from_u64(1_000_003);
    let wide = odd.clone() << 900;
    assert_eq!(wide.bit_len(), 920);
    let factors = factor(wide.clone()).unwrap();
    assert!(factors.verify_product(&wide));
    assert_eq!(
        factors.multiplicity(&Natural::from_u64(2)).unwrap().get(),
        900
    );
    assert_eq!(factors.multiplicity(&odd).unwrap().get(), 1);
}

/// Above the sieve's ceiling Pollard-Brent is not an opportunistic peel in front of SIQS — it is
/// the entire factoring attempt, because the sieve refuses the composite outright. Every input
/// here is a small prime times a wide prime, which is exactly the shape that used to be rejected
/// with `SiqsCompositeTooLarge` whenever the small prime outran a budget sized as a fraction of a
/// sieve run that never happens.
///
/// These three are the guaranteed tier: factors up to 32 bits, at the bottom, middle and top of the
/// supported width range. The budget above the ceiling is hundreds of times Brent's expected
/// `1.2·sqrt(p)` cost for a 32-bit factor at every one of those widths.
#[test]
fn wide_composites_are_split_by_rho_rather_than_refused() {
    // 401 bits — one over the sieve ceiling — carrying a 16-bit factor, which is above the 10^4
    // trial-division bound and so genuinely reaches rho.
    let input = Natural::<16>::from_decimal(
        "3902936505210420261021964138186048598645757055007853121964024023082987217643\
         286159953866070270752263666313381145517080159",
    )
    .unwrap();
    assert_eq!(input.bit_len(), 401);
    let factors = factor(input.clone()).unwrap();
    assert!(factors.verify_product(&input));
    assert_eq!(factors.total_len(), 2);
    assert_eq!(
        factors
            .multiplicity(&Natural::from_u64(56_857))
            .unwrap()
            .get(),
        1
    );

    // 512 bits with a 32-bit factor: the width a caller is most likely to hand a library capped at
    // 400 bits of sieve.
    let input = Natural::<16>::from_decimal(
        "9501012405705509564680437712617447440170980081112656222237073910419870316392\
         859702111963091481439276805995800801743430916377894473378632368751322056628119",
    )
    .unwrap();
    assert_eq!(input.bit_len(), 512);
    let factors = factor(input.clone()).unwrap();
    assert!(factors.verify_product(&input));
    assert_eq!(factors.total_len(), 2);
    assert_eq!(
        factors
            .multiplicity(&Natural::from_u64(3_667_435_003))
            .unwrap()
            .get(),
        1
    );

    // 1024 bits, the engine's capacity, with the same 32-bit factor class. Per-iteration cost is
    // roughly three times the 512-bit case, which is why the budget tiers by width: the guarantee
    // is stated in factor bits and has to hold at the expensive end too.
    let input = Natural::<16>::from_decimal(
        "1401022297307954297991679230211882827730660578524007095304227670281026951677\
         4811096642399661825372091824156874528352425797533201535316521537410076121003\
         4074640519308908102560817045658913006325713539052727661523689985681266686274\
         209226917233091828080256346478868845059891172761884504662361609922427332741249613",
    )
    .unwrap();
    assert_eq!(input.bit_len(), 1024);
    let factors = factor(input.clone()).unwrap();
    assert!(factors.verify_product(&input));
    assert_eq!(factors.total_len(), 2);
    assert_eq!(
        factors
            .multiplicity(&Natural::from_u64(3_479_286_313))
            .unwrap()
            .get(),
        1
    );
}

/// The tier the raised budget was actually bought for. A 48-bit factor costs about 20 M iterations,
/// three times the 6.29 M the sieve-derived budget allowed at every width above the ceiling, so
/// before this change the 512-bit input below was refused after 2.7 s of work that was nearly deep
/// enough; it now splits in 30 s. Unoptimized builds run rho at roughly a fifteenth of release
/// speed, hence the profile.
#[test]
#[ignore = "40- and 48-bit factors cost seconds to minutes of rho: cargo test --profile release-test"]
fn wide_composites_split_40_and_48_bit_factors_with_the_default_budget() {
    for (decimal, small) in [
        (
            "1015936249400950287360808203240867131305086909045547514124726594384293137\
             1677588336226503124482280315907519033170704907241245468573613367599493866551320501",
            971_259_191_779u64,
        ),
        (
            "1012682424890306230360672028802799323297238247563131601312546154540901930\
             5100244454696718048789411013895449112084116270765951273303169020877632325189728337",
            263_495_738_435_519,
        ),
    ] {
        let input = Natural::<16>::from_decimal(decimal).unwrap();
        assert_eq!(input.bit_len(), 512);
        let factors = factor(input.clone()).unwrap();
        assert!(factors.verify_product(&input));
        assert_eq!(factors.total_len(), 2);
        assert_eq!(
            factors
                .multiplicity(&Natural::from_u64(small))
                .unwrap()
                .get(),
            1
        );
    }
}

/// A wide composite made of many middling primes — the shape that has no balanced semiprime in it
/// anywhere, and that the ladder used to abandon partway down.
///
/// Rho peeled factors while the cofactor was above the sieve's ceiling and had the deep budget, then
/// the cofactor crossed 400 bits, the budget collapsed to a fraction of a sieve run, and the
/// remainder — still six or eight 50-bit primes — went to SIQS at the 369..=400 tier, which wanted
/// 206,403 relations at about two per second. That is weeks of sieving on a number whose every
/// factor rho finds in seconds. A split under rho now marks the cofactor as unbalanced and keeps the
/// deep budget down to `DEEP_RHO_MIN_BITS`, below which the sieve really is the faster tool: this
/// input peels five factors in rho and hands a 250-bit remainder to a sieve that returns it in
/// under three seconds.
///
/// Measured end to end on an x86-64 Xeon 8259CL with 96 workers: 65 s through the release CLI, 67 s
/// as this test under `--profile release-test`.
#[test]
#[ignore = "a minute of rho and one 250-bit sieve: cargo test --profile release-test"]
fn wide_products_of_many_middling_primes_do_not_stall_in_the_sieve() {
    let input = Natural::<16>::from_decimal(
        "6921330803157597027523689283187318825197185831714411528408744449716399940722\
         77047941800262750998157133588227123597503809544337586905787678816852754733",
    )
    .unwrap();
    assert_eq!(input.bit_len(), 498);
    let factors = factor(input.clone()).unwrap();
    assert!(factors.verify_product(&input));
    assert_eq!(factors.total_len(), 10);
    assert_eq!(
        factors
            .expanded()
            .map(ToString::to_string)
            .collect::<Vec<_>>(),
        [
            "846008445527897",
            "865990488546691",
            "871893255711629",
            "901642578578681",
            "923403409425967",
            "957330526028417",
            "1031228931434011",
            "1087834615977859",
            "1099924897280071",
            "1101716055838291",
        ]
    );
}

/// What the default budget deliberately does not pay for. Each factor bit doubles Brent's cost, so
/// 56 bits is minutes and 64 bits is tens of minutes to hours — a real search rather than a stage
/// in a ladder. It is reachable, and this pins the mechanism that reaches it: the same override the
/// CLI exposes as `RUSQSIEVE_RHO_ITERATIONS`.
///
/// Measured single-threaded on an x86-64 Xeon 8259CL, `--profile release-test`: the two cases
/// together took 1,452 s. Run with `--nocapture` for the per-case split.
#[test]
#[ignore = "tens of minutes of Pollard-Brent by design: cargo test --profile release-test"]
fn raised_budgets_reach_56_and_64_bit_factors_above_the_ceiling() {
    // 4 G iterations is about three times Brent's expected cost for a 64-bit factor, which covers
    // the spread around the expectation rather than only its centre.
    let deep = FactorConfig::default()
        .with_parallelism(Parallelism::threads(1).unwrap())
        .with_tuning_overrides(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(4_000_000_000),
            false,
        );
    for (decimal, bits, small) in [
        (
            "1252065390128543989678566233540615804226711829449772216415528509501425921\
             5604536975970501806596235777062580778663894798784169636410299447575447913433618819",
            512,
            67_876_963_450_791_707u64,
        ),
        (
            "3669910812722365764460471742391282485897850073502478531673082844008509335\
             074097393459666891883936817708908023289961021419",
            401,
            15_590_429_627_348_277_103,
        ),
    ] {
        let input = Natural::<16>::from_decimal(decimal).unwrap();
        assert_eq!(input.bit_len(), bits);
        let started = std::time::Instant::now();
        let factors =
            factor_with_progress(input.clone(), deep.clone(), |_| ProgressAction::Continue)
                .unwrap();
        eprintln!(
            "BENCH deep_rho input_bits={bits} factor_bits={} elapsed={:.1}s",
            64 - small.leading_zeros(),
            started.elapsed().as_secs_f64()
        );
        assert!(factors.verify_product(&input));
        assert_eq!(factors.total_len(), 2);
        assert_eq!(
            factors
                .multiplicity(&Natural::from_u64(small))
                .unwrap()
                .get(),
            1
        );
    }
}

/// The other half of the contract: a composite that genuinely needs the sieve and is too wide for
/// it is refused with an error naming the composite, not the input.
#[test]
#[ignore = "burns a full Pollard-Brent budget — about a minute above the ceiling — on a 416-bit semiprime before the sieve is asked"]
fn oversized_hard_composites_are_refused_by_the_sieve_not_the_input_width() {
    let hard = Natural::<16>::from_decimal(
        "1265686468695484903648964277331152191512075117088221193348532259113389592\
         93458582681292290648687613986290058064435220042047901",
    )
    .unwrap();
    assert_eq!(hard.bit_len(), 416);
    assert!(matches!(
        factor(hard),
        Err(FactorError::SiqsCompositeTooLarge(416))
    ));
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

/// Every entry of the supplied corpus, parsed. The arity varies from 1 to 10, so the whole line is
/// read and no entry may be assumed to be a semiprime.
fn corpus_entries() -> Vec<(&'static str, Vec<&'static str>)> {
    let corpus = include_str!("data/rusqsieve-factorization-corpus.txt");
    let mut entries = Vec::new();
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
        entries.push((fields[0], fields[1..].to_vec()));
    }
    assert_eq!(entries.len(), 309, "unexpected corpus entry count");
    entries
}

#[test]
fn browser_tuning_corpus_has_exact_balanced_products() {
    let corpus = include_str!("data/browser-balanced-corpus.txt");
    let mut counts = [0usize; 7];
    for (line_index, line) in corpus.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split_ascii_whitespace().collect::<Vec<_>>();
        assert_eq!(
            fields.len(),
            4,
            "bad browser corpus line {}",
            line_index + 1
        );
        let bits: usize = fields[0].parse().unwrap();
        let n = Natural::<16>::from_str(fields[1]).unwrap();
        let p = Natural::<16>::from_str(fields[2]).unwrap();
        let q = Natural::<16>::from_str(fields[3]).unwrap();
        assert_eq!(n.bit_len(), bits, "wrong width on line {}", line_index + 1);
        assert_eq!(p.bit_len(), bits / 2);
        assert_eq!(q.bit_len(), bits / 2);
        assert_eq!(p.checked_mul(&q).unwrap(), n);
        counts[match bits {
            216 => 0,
            224 => 1,
            232 => 2,
            240 => 3,
            256 => 4,
            272 => 5,
            288 => 6,
            _ => panic!("unexpected browser tier {bits}"),
        }] += 1;
    }
    assert_eq!(counts, [5, 5, 5, 5, 5, 5, 1]);
}

/// Where the default corpus test stops and the `#[ignore]`d one picks up. Every one of the 117
/// entries in the 65-85-bit `choose_a` dead zone — the regression this corpus exists for — is below
/// it, with margin. Above it the entries are genuine SIQS work: a 256-bit factorization is tens of
/// seconds in a release build and minutes in the default unoptimized test profile, which would make
/// `cargo test` too slow to be run.
const DEFAULT_TIER_BITS: usize = 128;

/// The corpus through [`DEFAULT_TIER_BITS`]. Note this exercises the ladder's cheap stages, not the
/// sieve: since the Pollard-Brent stage was added, everything in this band splits before SIQS is
/// reached. The dead zone is pinned against the sieve itself by `engine`'s
/// `siqs_builds_polynomials_across_the_dead_zone` and `siqs_alone_factors_the_dead_zone`.
#[test]
fn supplied_factorization_corpus_through_128_bits() {
    let mut tested = 0usize;
    let mut in_dead_zone = 0usize;
    for (n, factors) in corpus_entries() {
        let bits = Natural::<16>::from_str(n).unwrap().bit_len();
        if bits > DEFAULT_TIER_BITS {
            continue;
        }
        assert_factorization(n, &factors);
        tested += 1;
        if (65..=85).contains(&bits) {
            in_dead_zone += 1;
        }
    }
    assert_eq!(
        in_dead_zone, 117,
        "the 65-85-bit regression band is not fully covered"
    );
    assert!(tested >= 270, "only {tested} entries ran");
}

/// The remaining corpus entries, up to 256 bits. Run with
/// `cargo test --profile release-test --test factorization -- --ignored`.
#[test]
#[ignore = "minutes of sieving; run in CI on a schedule"]
fn supplied_factorization_corpus_above_128_bits() {
    let mut tested = 0usize;
    for (n, factors) in corpus_entries() {
        if Natural::<16>::from_str(n).unwrap().bit_len() <= DEFAULT_TIER_BITS {
            continue;
        }
        assert_factorization(n, &factors);
        tested += 1;
    }
    assert_eq!(tested, 29, "unexpected large corpus entry count");
}

#[test]
#[ignore = "manual interleaved performance measurement"]
fn profile_single_input() {
    // `--ignored` runs this alongside the slow corpus tier, so an absent input is "nothing to
    // measure" rather than a failure.
    let Ok(input) = std::env::var("RUSQSIEVE_BENCH_INPUT") else {
        eprintln!("BENCH skipped: set RUSQSIEVE_BENCH_INPUT to a decimal integer");
        return;
    };
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
