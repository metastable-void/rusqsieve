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
