//! Lenstra's elliptic curve method, opt-in.
//!
//! ECM is the right tool for exactly the composites the rest of this crate is wrong for: ones with
//! a medium-size factor. Pollard-Brent costs `O(sqrt p)` in the smallest factor and runs out around
//! 2^50; the quadratic sieve charges by the size of `N` and refuses anything past 400 bits. ECM
//! costs subexponentially in `p` rather than `N`, so it reaches 25- to 35-digit factors of numbers
//! either of the others would decline.
//!
//! A balanced semiprime inside the sieve's range never pays for a curve unless the caller asks,
//! and that is a hard requirement rather than a default: it is the workload this crate's claim
//! rests on, and no curve can succeed on it. Everything else reaches ECM on its own, because
//! everything else has already given away that it is a different shape — the composite is wider
//! than the sieve accepts, or trial division peeled a factor, or Pollard–Brent split an ancestor.
//! See `FactorConfig::with_ecm` and `engine::factor_node`.
//!
//! # What is implemented
//!
//! Montgomery curves in `(X : Z)` coordinates with Suyama's `σ` parameterization, which forces the
//! group order divisible by 12. Stage 1 raises the point to `lcm(1..B1)` one prime power at a time
//! using PRAC addition chains — Montgomery's practical chains, about 1.55 point operations per bit
//! against the binary ladder's 2. Stage 2 is the standard continuation: a wheel of baby steps
//! coprime to `D`, giant steps of `[D]Q`, and one gcd for the whole stage rather than one per
//! prime.
//!
//! That is materially stronger than a textbook implementation — FLINT's, for comparison, uses a
//! binary ladder and takes a gcd after every prime — and short of GMP-ECM, whose stage 2 evaluates
//! a polynomial at many points at once and whose stage 1 has assembly modular arithmetic. The gap
//! is stage 2's asymptotics, not correctness.
use crate::natural::{LIMB_CAP, Limb, MontgomeryContext};
use crate::{Natural, PARTS};

/// A projective point `(X : Z)` on a Montgomery curve, in Montgomery residue form.
///
/// `Y` is never needed: doubling and differential addition determine `x`-coordinates alone, which
/// is the whole reason this curve shape is used for factoring.
#[derive(Clone, Copy)]
struct Point {
    x: [Limb; LIMB_CAP],
    z: [Limb; LIMB_CAP],
}

impl Point {
    const ZERO: Self = Self {
        x: [0 as Limb; LIMB_CAP],
        z: [0 as Limb; LIMB_CAP],
    };
}

/// Work bounds for one ECM run.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct EcmParams {
    /// Stage 1 bound: the run finds `p` when `p + 1 − a_E` is `B1`-smooth on some curve.
    pub(crate) b1: u64,
    /// Stage 2 bound, which additionally catches a single prime factor in `(B1, B2]`.
    pub(crate) b2: u64,
    /// Curves to try before giving up.
    pub(crate) curves: u32,
}

impl EcmParams {
    /// Bounds for a composite of `bits` bits.
    ///
    /// The schedule is deliberately asymmetric around the sieve's ceiling, because what ECM is
    /// competing against changes completely there.
    ///
    /// At or below it, SIQS is going to run and finish, so curves are optional extra work on top of
    /// a solution that already exists. They are budgeted to seconds: enough to catch an unbalanced
    /// input cheaply, never enough to delay the sieve meaningfully. This range is opt-in.
    ///
    /// Above it there is no sieve — the alternative to a curve is refusing the number. The budget
    /// there is minutes and the run is unconditional, because a bounded wait beats an answer of
    /// "too large" on a composite whose factor is 25 or 30 digits.
    ///
    /// Bounds follow the standard `B1`/curve-count levels — 2k for 15 digits, 11k for 20, 50k for
    /// 25, 250k for 30 — with `B2 = 100·B1`.
    pub(crate) fn for_composite(bits: usize) -> Self {
        let (b1, curves) = if bits <= 256 {
            // The sieve handles these in seconds; a handful of cheap curves only has to beat the
            // chance that the input is unbalanced.
            (2_000, 16)
        } else if bits <= crate::engine::MAX_SIQS_BITS {
            (11_000, 48)
        } else if bits <= 512 {
            // Past the ceiling: no fallback, so this is the level that reaches 25-digit factors.
            (50_000, 300)
        } else {
            (250_000, 500)
        };
        Self {
            b1,
            b2: b1.saturating_mul(100),
            curves,
        }
    }
}

/// The wheel modulus for stage 2's baby steps. 210 = 2·3·5·7 leaves 24 of every 210 residues
/// coprime to it, so the table holds 24 points per 210 integers covered.
const WHEEL: u64 = 210;

/// Searches for a nontrivial factor of `n`.
///
/// Returns `Ok(None)` when every curve is exhausted without one, which is a result rather than a
/// failure: the caller falls through to the sieve, or reports that the composite is out of range.
pub(crate) fn factor(
    n: &Natural,
    params: EcmParams,
    seed: u64,
    mut keep_going: impl FnMut() -> bool,
) -> Result<Option<Natural>, crate::engine::EngineError> {
    if n.is_even() {
        return Ok(Some(Natural::from_u64(2)));
    }
    let Some(montgomery) = MontgomeryContext::new(n) else {
        return Ok(None);
    };
    // Stage 1 walks every prime power up to B1; stage 2 needs the primes in (B1, B2].
    let stage1_primes = crate::smallfactor::sieve_primes(params.b1.min(u64::from(u32::MAX)) as u32);
    let stage2_primes = prime_bitmap(params.b1, params.b2);

    for curve in 0..params.curves {
        if !keep_going() {
            return Err(crate::engine::EngineError::Cancelled);
        }
        // Suyama's σ must avoid 0, 1, 2, 3 and 5, which give degenerate curves. Curves are drawn
        // from a counter rather than an entropy source so a run is reproducible.
        let sigma = 6 + seed
            .wrapping_add(u64::from(curve))
            .wrapping_mul(2_654_435_761)
            % 1_000_000;
        let (mut point, a24) = match build_curve(&montgomery, n, sigma) {
            CurveSetup::Ready(ready) => (ready.0, ready.1),
            // The inverse did not exist, which means the denominator shared a factor with `n`.
            CurveSetup::Factor(factor) => return Ok(Some(factor)),
            CurveSetup::Degenerate => continue,
        };

        if let Some(factor) = stage_one(
            &montgomery,
            n,
            &a24,
            &mut point,
            &stage1_primes,
            params.b1,
            &mut keep_going,
        )? {
            return Ok(Some(factor));
        }
        if let Some(factor) = stage_two(
            &montgomery,
            n,
            &a24,
            &point,
            &stage2_primes,
            params,
            &mut keep_going,
        )? {
            return Ok(Some(factor));
        }
    }
    Ok(None)
}

/// Outcome of building a curve from a `σ`.
enum CurveSetup {
    /// Boxed because a point plus a curve constant is several hundred bytes against a `Natural`
    /// handful, and this is returned once per curve rather than in any loop.
    Ready(Box<(Point, [Limb; LIMB_CAP])>),
    Factor(Natural),
    Degenerate,
}

/// Suyama's parameterization: `u = σ² − 5`, `v = 4σ`, the starting point is `(u³ : v³)`, and the
/// curve constant is `a24 = (A + 2)/4 = (v − u)³(3u + v) / (16u³v)`.
///
/// The single inversion per curve is the only division in the whole method, and when it fails it
/// fails usefully: a denominator sharing a factor with `n` *is* the factor.
fn build_curve(montgomery: &MontgomeryContext<PARTS>, n: &Natural, sigma: u64) -> CurveSetup {
    let mut scratch = Scratch::new();
    let sigma = montgomery.encode(
        &Natural::from_u64(sigma)
            .div_rem(n)
            .map_or(Natural::ZERO, |(_, r)| r),
    );
    let mut u = [0 as Limb; LIMB_CAP];
    let mut v = [0 as Limb; LIMB_CAP];
    let mut five = [0 as Limb; LIMB_CAP];
    montgomery.load(&sigma, &mut scratch.a);
    montgomery.load(&montgomery.encode(&Natural::from_u64(5)), &mut five);

    // v = 4σ
    montgomery.add_raw(&scratch.a, &scratch.a, &mut v);
    let doubled = v;
    montgomery.add_raw(&doubled, &doubled, &mut v);
    // u = σ² − 5
    montgomery.sqr_raw(&scratch.a, &mut scratch.b);
    montgomery.sub_raw(&scratch.b, &five, &mut u);

    // X0 = u³, Z0 = v³
    let mut point = Point::ZERO;
    montgomery.sqr_raw(&u, &mut scratch.b);
    montgomery.mul_raw(&scratch.b, &u, &mut point.x);
    montgomery.sqr_raw(&v, &mut scratch.b);
    montgomery.mul_raw(&scratch.b, &v, &mut point.z);

    // numerator = (v − u)³(3u + v)
    montgomery.sub_raw(&v, &u, &mut scratch.a);
    montgomery.sqr_raw(&scratch.a, &mut scratch.b);
    montgomery.mul_raw(&scratch.b, &scratch.a, &mut scratch.c);
    montgomery.add_raw(&u, &u, &mut scratch.a);
    montgomery.add_raw(&scratch.a, &u, &mut scratch.b);
    montgomery.add_raw(&scratch.b, &v, &mut scratch.a);
    let mut numerator = [0 as Limb; LIMB_CAP];
    montgomery.mul_raw(&scratch.c, &scratch.a, &mut numerator);

    // denominator = 16·u³·v
    montgomery.mul_raw(&point.x, &v, &mut scratch.a);
    for _ in 0..4 {
        let doubled = scratch.a;
        montgomery.add_raw(&doubled, &doubled, &mut scratch.a);
    }
    let denominator = montgomery.decode(&montgomery.store(&scratch.a));
    if denominator.is_zero() {
        return CurveSetup::Degenerate;
    }
    let inverse = match mod_inverse(&denominator, n) {
        Inversion::Inverse(inverse) => inverse,
        Inversion::Factor(factor) => return CurveSetup::Factor(factor),
    };

    let mut a24 = [0 as Limb; LIMB_CAP];
    montgomery.load(&montgomery.encode(&inverse), &mut scratch.a);
    montgomery.mul_raw(&numerator, &scratch.a, &mut a24);
    CurveSetup::Ready(Box::new((point, a24)))
}

/// Reusable temporaries, so no point operation allocates or clears a buffer it does not need.
struct Scratch {
    a: [Limb; LIMB_CAP],
    b: [Limb; LIMB_CAP],
    c: [Limb; LIMB_CAP],
}

impl Scratch {
    const fn new() -> Self {
        Self {
            a: [0 as Limb; LIMB_CAP],
            b: [0 as Limb; LIMB_CAP],
            c: [0 as Limb; LIMB_CAP],
        }
    }
}

/// `2P` from the sum and difference of `P`'s coordinates, which the caller usually has already.
fn double_from(
    montgomery: &MontgomeryContext<PARTS>,
    a24: &[Limb; LIMB_CAP],
    sum: &[Limb; LIMB_CAP],
    difference: &[Limb; LIMB_CAP],
) -> Point {
    let mut out = Point::ZERO;
    let mut squared_difference = [0 as Limb; LIMB_CAP];
    let mut squared_sum = [0 as Limb; LIMB_CAP];
    montgomery.sqr_raw(difference, &mut squared_difference);
    montgomery.sqr_raw(sum, &mut squared_sum);
    montgomery.mul_raw(&squared_difference, &squared_sum, &mut out.x);
    let mut spread = [0 as Limb; LIMB_CAP];
    montgomery.sub_raw(&squared_sum, &squared_difference, &mut spread);
    let mut scaled = [0 as Limb; LIMB_CAP];
    montgomery.mul_raw(&spread, a24, &mut scaled);
    let mut inner = [0 as Limb; LIMB_CAP];
    montgomery.add_raw(&scaled, &squared_difference, &mut inner);
    montgomery.mul_raw(&inner, &spread, &mut out.z);
    out
}

/// `2P`.
fn double(montgomery: &MontgomeryContext<PARTS>, a24: &[Limb; LIMB_CAP], p: &Point) -> Point {
    let mut sum = [0 as Limb; LIMB_CAP];
    let mut difference = [0 as Limb; LIMB_CAP];
    montgomery.add_raw(&p.x, &p.z, &mut sum);
    montgomery.sub_raw(&p.x, &p.z, &mut difference);
    double_from(montgomery, a24, &sum, &difference)
}

/// Differential addition: `P1 + P2` given `P1 − P2`.
fn add(montgomery: &MontgomeryContext<PARTS>, p1: &Point, p2: &Point, difference: &Point) -> Point {
    let mut sum1 = [0 as Limb; LIMB_CAP];
    let mut diff1 = [0 as Limb; LIMB_CAP];
    let mut sum2 = [0 as Limb; LIMB_CAP];
    let mut diff2 = [0 as Limb; LIMB_CAP];
    montgomery.add_raw(&p1.x, &p1.z, &mut sum1);
    montgomery.sub_raw(&p1.x, &p1.z, &mut diff1);
    montgomery.add_raw(&p2.x, &p2.z, &mut sum2);
    montgomery.sub_raw(&p2.x, &p2.z, &mut diff2);

    let mut cross1 = [0 as Limb; LIMB_CAP];
    let mut cross2 = [0 as Limb; LIMB_CAP];
    montgomery.mul_raw(&diff1, &sum2, &mut cross1);
    montgomery.mul_raw(&sum1, &diff2, &mut cross2);

    let mut total = [0 as Limb; LIMB_CAP];
    let mut gap = [0 as Limb; LIMB_CAP];
    montgomery.add_raw(&cross1, &cross2, &mut total);
    montgomery.sub_raw(&cross1, &cross2, &mut gap);

    let mut squared_total = [0 as Limb; LIMB_CAP];
    let mut squared_gap = [0 as Limb; LIMB_CAP];
    montgomery.sqr_raw(&total, &mut squared_total);
    montgomery.sqr_raw(&gap, &mut squared_gap);

    let mut out = Point::ZERO;
    montgomery.mul_raw(&squared_total, &difference.z, &mut out.x);
    montgomery.mul_raw(&squared_gap, &difference.x, &mut out.z);
    out
}

/// Montgomery's PRAC: `[k]P` through an addition chain chosen by a golden-ratio-like split.
///
/// The chain costs about 1.55 point operations per bit of `k` against the binary ladder's 2, which
/// is most of stage 1's cost. The rule set is Montgomery's; the ratios are the standard ones, and
/// the chain is run for whichever ratio yields the cheapest chain for this `k`.
fn prac(
    montgomery: &MontgomeryContext<PARTS>,
    a24: &[Limb; LIMB_CAP],
    point: &Point,
    k: u64,
    ratio: f64,
) -> Point {
    if k == 0 {
        return *point;
    }
    let shift = k.trailing_zeros();
    let mut k = k >> shift;

    let mut result = *point;
    if k > 1 {
        let r = (k as f64 * ratio + 0.5) as u64;
        let mut d = k - r;
        let mut e = 2 * r - k;

        // A is [1]P doubled to [2]P; B and C stay at [1]P. The invariant through the loop is that
        // A − B = C, which is what makes every step a differential addition.
        let mut a = double(montgomery, a24, point);
        let mut b = *point;
        let mut c = *point;

        while d != e {
            if d < e {
                core::mem::swap(&mut d, &mut e);
                core::mem::swap(&mut a, &mut b);
            }
            if d - e <= e / 4 && (d + e).is_multiple_of(3) {
                d = (2 * d - e) / 3;
                e = (e - d) / 2;
                let t = add(montgomery, &a, &b, &c);
                let t2 = add(montgomery, &t, &a, &b);
                b = add(montgomery, &b, &t, &a);
                a = t2;
            } else if d - e <= e / 4 && (d - e).is_multiple_of(6) {
                d = (d - e) / 2;
                b = add(montgomery, &a, &b, &c);
                a = double(montgomery, a24, &a);
            } else if d.div_ceil(4) <= e {
                d -= e;
                let t = add(montgomery, &b, &a, &c);
                c = core::mem::replace(&mut b, t);
            } else if (d + e).is_multiple_of(2) {
                d = (d - e) / 2;
                b = add(montgomery, &b, &a, &c);
                a = double(montgomery, a24, &a);
            } else if d.is_multiple_of(2) {
                d /= 2;
                c = add(montgomery, &c, &a, &b);
                a = double(montgomery, a24, &a);
            } else if d.is_multiple_of(3) {
                d = d / 3 - e;
                let t = double(montgomery, a24, &a);
                let t2 = add(montgomery, &a, &b, &c);
                a = add(montgomery, &t, &a, &a);
                let t = add(montgomery, &t, &t2, &c);
                c = b;
                b = t;
            } else if (d + e).is_multiple_of(3) {
                d = (d - 2 * e) / 3;
                let t = add(montgomery, &a, &b, &c);
                b = add(montgomery, &t, &a, &b);
                let t = double(montgomery, a24, &a);
                a = add(montgomery, &a, &t, &a);
            } else if (d - e).is_multiple_of(3) {
                d = (d - e) / 3;
                let t = add(montgomery, &a, &b, &c);
                c = add(montgomery, &c, &a, &b);
                b = t;
                let t = double(montgomery, a24, &a);
                a = add(montgomery, &a, &t, &a);
            } else {
                e /= 2;
                c = add(montgomery, &c, &b, &a);
                b = double(montgomery, a24, &b);
            }
        }
        result = add(montgomery, &a, &b, &c);
        k = 1;
    }
    let _ = k;

    for _ in 0..shift {
        result = double(montgomery, a24, &result);
    }
    result
}

/// The ratios Montgomery's chains are searched over; the first is the golden ratio's reciprocal.
const PRAC_RATIOS: [f64; 4] = [
    0.618_033_988_749_894_9,
    0.723_606_797_749_979,
    0.580_178_728_295_464_1,
    0.632_839_808_081_518,
];

/// Stage 1: raise the point to `lcm(1..B1)`, one prime power at a time.
///
/// A single gcd at the end covers the whole stage. Taking one per prime, as a textbook version
/// does, would cost more than the point arithmetic it guards.
fn stage_one(
    montgomery: &MontgomeryContext<PARTS>,
    n: &Natural,
    a24: &[Limb; LIMB_CAP],
    point: &mut Point,
    primes: &[u32],
    b1: u64,
    keep_going: &mut impl FnMut() -> bool,
) -> Result<Option<Natural>, crate::engine::EngineError> {
    for (index, &prime) in primes.iter().enumerate() {
        if index % 256 == 0 && !keep_going() {
            return Err(crate::engine::EngineError::Cancelled);
        }
        let prime = u64::from(prime);
        // Every power of this prime that stays under B1 has to be applied.
        let mut power = prime;
        while power <= b1 {
            let ratio = PRAC_RATIOS[(prime % PRAC_RATIOS.len() as u64) as usize];
            *point = prac(montgomery, a24, point, prime, ratio);
            power = match power.checked_mul(prime) {
                Some(next) => next,
                None => break,
            };
        }
    }
    Ok(check(montgomery, n, &point.z))
}

/// Stage 2, standard continuation.
///
/// Stage 1 succeeds when the group order is `B1`-smooth. Stage 2 additionally catches an order that
/// is `B1`-smooth except for one prime in `(B1, B2]`, which is a far more common shape, and it does
/// so for a fraction of stage 1's cost: baby steps `[j]Q` for `j` coprime to the wheel, giant steps
/// `[iD]Q`, and a product of `x`-coordinate cross differences that is gcd-ed once at the end.
fn stage_two(
    montgomery: &MontgomeryContext<PARTS>,
    n: &Natural,
    a24: &[Limb; LIMB_CAP],
    q: &Point,
    primes: &PrimeBitmap,
    params: EcmParams,
    keep_going: &mut impl FnMut() -> bool,
) -> Result<Option<Natural>, crate::engine::EngineError> {
    // Baby steps: [j]Q for every j < WHEEL/2 coprime to WHEEL.
    let residues: Vec<u64> = (1..WHEEL / 2).filter(|j| gcd_u64(*j, WHEEL) == 1).collect();
    let mut baby = Vec::with_capacity(residues.len());
    for &j in &residues {
        baby.push(prac(montgomery, a24, q, j, PRAC_RATIOS[0]));
    }

    let step = prac(montgomery, a24, q, WHEEL, PRAC_RATIOS[0]);
    let first = params.b1 / WHEEL;
    let mut giant = prac(montgomery, a24, q, first.max(1) * WHEEL, PRAC_RATIOS[0]);
    // The differential addition that advances a giant step needs the previous one.
    let mut previous = prac(
        montgomery,
        a24,
        q,
        first.max(1).saturating_sub(1) * WHEEL,
        PRAC_RATIOS[0],
    );

    let mut accumulator = [0 as Limb; LIMB_CAP];
    montgomery.load(&montgomery.one(), &mut accumulator);
    let mut cross = [0 as Limb; LIMB_CAP];
    let mut left = [0 as Limb; LIMB_CAP];
    let mut right = [0 as Limb; LIMB_CAP];

    let mut center = first.max(1) * WHEEL;
    while center <= params.b2 + WHEEL {
        if !keep_going() {
            return Err(crate::engine::EngineError::Cancelled);
        }
        for (index, &j) in residues.iter().enumerate() {
            // Only pairs that straddle a prime in range contribute; the rest would multiply in a
            // value that cannot reveal anything.
            let below = center.checked_sub(j).is_some_and(|v| primes.contains(v));
            let above = primes.contains(center + j);
            if !below && !above {
                continue;
            }
            // x(R)·z(S) − z(R)·x(S) vanishes modulo p exactly when the two points coincide there.
            montgomery.mul_raw(&giant.x, &baby[index].z, &mut left);
            montgomery.mul_raw(&giant.z, &baby[index].x, &mut right);
            montgomery.sub_raw(&left, &right, &mut cross);
            if !montgomery.is_zero_raw(&cross) {
                montgomery.mul_assign(&mut accumulator, &cross);
            }
        }
        let next = add(montgomery, &giant, &step, &previous);
        previous = giant;
        giant = next;
        center += WHEEL;
    }
    Ok(check(montgomery, n, &accumulator))
}

/// Turns an accumulated residue into a factor, if it shares one with `n`.
fn check(
    montgomery: &MontgomeryContext<PARTS>,
    n: &Natural,
    value: &[Limb; LIMB_CAP],
) -> Option<Natural> {
    let candidate = montgomery.store(value).gcd(n);
    if candidate.is_one() || candidate == *n {
        None
    } else {
        Some(candidate)
    }
}

/// Result of inverting modulo `n`, where failure is informative.
enum Inversion {
    Inverse(Natural),
    Factor(Natural),
}

/// `a^-1 mod n` by the extended Euclidean algorithm, with the cofactor kept reduced so no signed
/// big integer is needed. A non-unit gcd is a factor of `n`, which is the outcome ECM is looking
/// for anyway.
fn mod_inverse(a: &Natural, n: &Natural) -> Inversion {
    let Some((_, mut remainder)) = a.div_rem(n) else {
        return Inversion::Factor(n.clone());
    };
    if remainder.is_zero() {
        return Inversion::Factor(n.clone());
    }
    let mut previous_remainder = n.clone();
    let mut cofactor = Natural::ZERO;
    let mut previous_cofactor = Natural::ONE;
    core::mem::swap(&mut previous_remainder, &mut remainder);
    core::mem::swap(&mut previous_cofactor, &mut cofactor);

    while !remainder.is_zero() {
        let Some((quotient, next_remainder)) = previous_remainder.div_rem(&remainder) else {
            break;
        };
        previous_remainder = core::mem::replace(&mut remainder, next_remainder);
        let product = quotient.mul_mod(&cofactor, n);
        let next_cofactor = if product.is_zero() {
            previous_cofactor.clone()
        } else {
            previous_cofactor.add_mod(&n.wrapping_sub(&product), n)
        };
        previous_cofactor = core::mem::replace(&mut cofactor, next_cofactor);
    }

    if previous_remainder.is_one() {
        Inversion::Inverse(previous_cofactor)
    } else {
        Inversion::Factor(previous_remainder)
    }
}

/// Odd-only bitmap of the primes in `(low, high]`, for stage 2's membership tests.
struct PrimeBitmap {
    low: u64,
    bits: Vec<u64>,
    high: u64,
}

impl PrimeBitmap {
    fn contains(&self, value: u64) -> bool {
        if value <= self.low || value > self.high || value.is_multiple_of(2) {
            return false;
        }
        let index = ((value - self.low - 1) / 2) as usize;
        self.bits
            .get(index / 64)
            .is_some_and(|word| word >> (index % 64) & 1 == 1)
    }
}

/// Sieves the primes in `(low, high]` by striking multiples of every prime up to `sqrt(high)`.
fn prime_bitmap(low: u64, high: u64) -> PrimeBitmap {
    if high <= low {
        return PrimeBitmap {
            low,
            bits: Vec::new(),
            high: low,
        };
    }
    let span = ((high - low) / 2 + 1) as usize;
    let mut bits = vec![u64::MAX; span.div_ceil(64)];
    let root = (high as f64).sqrt() as u32 + 1;
    for &prime in crate::smallfactor::sieve_primes(root).iter() {
        let prime = u64::from(prime);
        if prime == 2 {
            continue;
        }
        // First odd multiple of `prime` strictly above `low`.
        let mut multiple = prime * prime;
        if multiple <= low {
            multiple = (low / prime + 1) * prime;
            if multiple.is_multiple_of(2) {
                multiple += prime;
            }
        }
        while multiple <= high {
            if multiple > low {
                let index = ((multiple - low - 1) / 2) as usize;
                if let Some(word) = bits.get_mut(index / 64) {
                    *word &= !(1u64 << (index % 64));
                }
            }
            multiple += 2 * prime;
        }
    }
    PrimeBitmap { low, bits, high }
}

/// Euclid on machine words, for the stage 2 wheel.
fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The property that matters: curves, chains, and both stages together actually recover a
    /// factor. Every piece of this file is silent when wrong — a broken addition chain simply never
    /// finds anything — so this is the gate.
    #[test]
    fn finds_a_medium_factor() {
        // 4,294,967,311 × 340,282,366,920,938,463,463,374,607,431,768,211,507. The small factor is
        // ten digits, which stage 1 reaches at this bound within a few curves.
        let n = Natural::from_decimal("1461501642435138422017761784666902131351499047677").unwrap();
        let params = EcmParams {
            b1: 2_000,
            b2: 200_000,
            curves: 64,
        };
        let found = factor(&n, params, 12345, || true).unwrap();
        let found = found.expect("ECM found no factor");
        assert!(!found.is_one() && found != n, "trivial factor {found}");
        assert!(
            n.div_rem(&found).unwrap().1.is_zero(),
            "{found} does not divide"
        );
    }
    /// ECM's actual territory: a 20-digit factor, which Pollard-Brent cannot reach at any budget
    /// this crate would spend and which the sieve would pay for by the size of `N` instead.
    #[test]
    #[ignore = "tens of curves at B1=50k: cargo test --profile release-test"]
    fn finds_a_twenty_digit_factor() {
        let n = Natural::from_decimal(
            "85927972517868228536373302480860846469815062148749486446868087553483220690888693",
        )
        .unwrap();
        let params = EcmParams {
            b1: 50_000,
            b2: 5_000_000,
            curves: 300,
        };
        let started = std::time::Instant::now();
        let found = factor(&n, params, 7, || true).unwrap().expect("no factor");
        eprintln!(
            "BENCH ecm 20-digit factor {found} in {:.2}s",
            started.elapsed().as_secs_f64()
        );
        assert_eq!(found.to_string(), "61218436624818344687");
    }
}
