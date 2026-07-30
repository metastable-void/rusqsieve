//! Shared machine-word number-theory helpers.

#[inline(always)]
pub(crate) fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

pub(crate) fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64 % modulus;
    base %= modulus;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = mul_mod(result, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    result
}

/// Deterministic Miller–Rabin. The seven-base Jaeschke/Sinclair set is proven
/// for every input below 2^64, so this is an exact primality test.
pub(crate) fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for prime in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == prime {
            return true;
        }
        if n.is_multiple_of(prime) {
            return false;
        }
    }
    let mut odd_part = n - 1;
    let shifts = odd_part.trailing_zeros();
    odd_part >>= shifts;
    'witness: for witness in [2u64, 325, 9375, 28178, 450775, 9780504, 1795265022] {
        let witness = witness % n;
        if witness == 0 {
            continue;
        }
        let mut value = pow_mod(witness, odd_part, n);
        if value == 1 || value == n - 1 {
            continue;
        }
        for _ in 1..shifts {
            value = mul_mod(value, value, n);
            if value == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

pub(crate) fn xorshift(mut value: u64) -> u64 {
    value ^= value << 13;
    value ^= value >> 7;
    value ^= value << 17;
    value
}

/// Shanks' square-forms factorization for inputs below 2^62.
///
/// Racing several small multipliers makes this substantially faster than rho
/// for the balanced word-size semiprimes left by double-large-prime sieving.
/// Failure is harmless: callers retain a bounded Pollard-Brent fallback.
pub(crate) fn squfof(n: u64) -> Option<u64> {
    const MULTIPLIERS: [u64; 16] = [
        1,
        3,
        5,
        7,
        11,
        3 * 5,
        3 * 7,
        3 * 11,
        5 * 7,
        5 * 11,
        7 * 11,
        3 * 5 * 7,
        3 * 5 * 11,
        3 * 7 * 11,
        5 * 7 * 11,
        3 * 5 * 7 * 11,
    ];
    const ONE_CYCLE: usize = 300;
    const MAX_CYCLES: usize = 100_000;

    if n < 4 {
        return None;
    }
    if n.is_multiple_of(2) {
        return Some(2);
    }

    let mut forms = Vec::with_capacity(MULTIPLIERS.len());
    for multiplier in MULTIPLIERS {
        let Some(scaled) = n.checked_mul(multiplier) else {
            break;
        };
        if scaled >= 1u64 << 62 {
            break;
        }
        let sqrt_n = isqrt(scaled);
        let q1 = scaled - sqrt_n * sqrt_n;
        if q1 == 0 {
            let factor = gcd(sqrt_n, n);
            return (factor > 1 && factor < n).then_some(factor);
        }
        forms.push(SqufofForm {
            multiplier,
            sqrt_n,
            cutoff: isqrt(2 * sqrt_n),
            q0: 1,
            p1: sqrt_n,
            q1,
            saved: Vec::with_capacity(50),
            failed: false,
        });
    }

    let mut iterations = 0usize;
    while iterations < MAX_CYCLES {
        let mut live = 0usize;
        for form in forms.iter_mut().rev() {
            if form.failed {
                continue;
            }
            live += 1;
            match squfof_cycle(n, form, ONE_CYCLE) {
                CycleResult::Factor(factor) => return Some(factor),
                CycleResult::Progress(done) => iterations += done,
                CycleResult::Failed(done) => {
                    iterations += done;
                    form.failed = true;
                }
            }
            if iterations >= MAX_CYCLES {
                break;
            }
        }
        if live == 0 {
            break;
        }
    }
    None
}

enum CycleResult {
    Factor(u64),
    Progress(usize),
    Failed(usize),
}

struct SqufofForm {
    multiplier: u64,
    sqrt_n: u64,
    cutoff: u64,
    q0: u64,
    p1: u64,
    q1: u64,
    saved: Vec<u64>,
    failed: bool,
}

fn squfof_cycle(n: u64, form: &mut SqufofForm, count: usize) -> CycleResult {
    let sqrt_n = form.sqrt_n;
    let multiplier = form.multiplier;
    let full_multiplier = 2 * multiplier;
    let coarse_cutoff = form.cutoff * full_multiplier;
    let mut q0 = form.q0;
    let mut p1 = form.p1;
    let mut q1 = form.q1;
    let mut p0 = 0u64;
    let mut square_root = 0u64;

    for iteration in 0..count {
        let quotient = 1 + (sqrt_n + p1 - q1) / q1;
        p0 = quotient * q1 - p1;
        q0 = if p1 >= p0 {
            q0 + (p1 - p0) * quotient
        } else {
            q0 - (p0 - p1) * quotient
        };
        if q1 < coarse_cutoff {
            let reduced = q1 / gcd(q1, full_multiplier);
            if reduced < form.cutoff {
                if form.saved.len() == 50 {
                    return CycleResult::Failed(iteration + 1);
                }
                form.saved.push(reduced);
            }
        }
        let zeros = q0.trailing_zeros();
        let odd = q0 >> zeros;
        if zeros & 1 == 0 && odd & 7 == 1 {
            square_root = isqrt(q0);
            if square_root * square_root == q0 && !form.saved.contains(&square_root) {
                if square_root == 1 {
                    return CycleResult::Failed(iteration + 1);
                }
                break;
            }
            square_root = 0;
        }

        let quotient = 1 + (sqrt_n + p0 - q0) / q0;
        p1 = quotient * q0 - p0;
        q1 = if p0 >= p1 {
            q1 + (p0 - p1) * quotient
        } else {
            q1 - (p1 - p0) * quotient
        };
    }

    if square_root == 0 {
        form.q0 = q0;
        form.p1 = p1;
        form.q1 = q1;
        return CycleResult::Progress(count);
    }

    q0 = square_root;
    p1 = p0 + square_root * ((sqrt_n - p0) / square_root);
    q1 = (n * multiplier - p1 * p1) / q0;
    loop {
        let quotient = 1 + (sqrt_n + p1 - q1) / q1;
        p0 = quotient * q1 - p1;
        q0 = if p1 >= p0 {
            q0 + (p1 - p0) * quotient
        } else {
            q0 - (p0 - p1) * quotient
        };
        if p0 == p1 {
            q0 = q1;
            break;
        }
        let quotient = 1 + (sqrt_n + p0 - q0) / q0;
        p1 = quotient * q0 - p0;
        q1 = if p0 >= p1 {
            q1 + (p0 - p1) * quotient
        } else {
            q1 - (p1 - p0) * quotient
        };
        if p0 == p1 {
            break;
        }
    }
    let factor = q0 / gcd(q0, full_multiplier);
    let factor = gcd(factor, n);
    if factor > 1 && factor < n {
        CycleResult::Factor(factor)
    } else {
        CycleResult::Failed(count)
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn isqrt(value: u64) -> u64 {
    if value < 2 {
        return value;
    }
    // Start at or above sqrt(value); starting one bit too low on values
    // whose floor-log2 is even makes Newton's decreasing convergence stop
    // immediately at an underestimate.
    let mut root = 1u64 << ((value.ilog2() + 2) / 2);
    loop {
        let next = (root + value / root) / 2;
        if next >= root {
            return root;
        }
        root = next;
    }
}

#[cfg(test)]
mod tests {
    use super::{isqrt, squfof};

    #[test]
    fn integer_square_root_is_exact_at_boundaries() {
        for value in [
            0,
            1,
            2,
            3,
            4,
            15,
            16,
            17,
            (1u64 << 40) + 123,
            u32::MAX as u64,
            (1u64 << 62) - 1,
        ] {
            let root = isqrt(value);
            assert!(root * root <= value);
            assert!((root + 1) * (root + 1) > value);
        }
    }

    #[test]
    fn squfof_splits_double_large_prime_scale_semiprimes() {
        for (left, right) in [
            (1_000_003u64, 1_000_033u64),
            (15_485_863, 15_485_867),
            (134_217_689, 134_217_757),
        ] {
            let n = left * right;
            let factor = squfof(n).expect("SQUFOF failed");
            assert!(factor == left || factor == right, "{n} -> {factor}");
        }
    }
}
