//! Native fast path for integers that fit in a machine word.
//!
//! Deterministic Miller–Rabin primality and Pollard–Brent factorization using
//! only `u64`/`u128` arithmetic. This bypasses fixed-capacity `Natural` big-integer
//! arithmetic entirely for small cofactors, which dominates cost below the SIQS
//! range. Results are deterministic (fixed sequence seeds), as required by SPEC §20.
#![cfg(any(unix, windows))]

use std::sync::OnceLock;

/// Trial-division bound for the cached small-prime table.
pub const SMALL_PRIME_BOUND: u32 = 1 << 16;

/// Primes below [`SMALL_PRIME_BOUND`], computed once (Sieve of Eratosthenes).
pub fn small_primes() -> &'static [u32] {
    static PRIMES: OnceLock<Vec<u32>> = OnceLock::new();
    PRIMES.get_or_init(|| sieve_primes(SMALL_PRIME_BOUND))
}

/// Primes `<= limit` via the Sieve of Eratosthenes.
pub fn sieve_primes(limit: u32) -> Vec<u32> {
    let n = limit as usize;
    if n < 2 {
        return Vec::new();
    }
    let mut is_composite = vec![false; n + 1];
    let mut out = Vec::new();
    let mut i = 2usize;
    while i <= n {
        if !is_composite[i] {
            out.push(i as u32);
            let mut j = i * i;
            while j <= n {
                is_composite[j] = true;
                j += i;
            }
        }
        i += 1;
    }
    out
}

#[inline(always)]
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64 % m;
    a %= m;
    while e != 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}

#[inline]
fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Deterministic Miller–Rabin. The 7-base set is a proven witness set for all
/// n < 2^64 (Jaeschke / Sinclair), so this is an exact primality test here.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &p in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let mut d = n - 1;
    let mut s = 0u32;
    while d & 1 == 0 {
        d >>= 1;
        s += 1;
    }
    'witness: for &a in &[2u64, 325, 9375, 28178, 450775, 9780504, 1795265022] {
        let a = a % n;
        if a == 0 {
            continue;
        }
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// Brent's improvement to Pollard's rho with batched GCD. `n` must be an odd
/// composite; returns a nontrivial factor. Deterministic given `n`.
fn pollard_brent(n: u64) -> u64 {
    if n.is_multiple_of(2) {
        return 2;
    }
    let mut c = 1u64;
    loop {
        let f = |x: u64| ((x as u128 * x as u128 + c as u128) % n as u128) as u64;
        let mut y = 2u64;
        let mut r = 1u64;
        let mut q = 1u64;
        let mut g = 1u64;
        let mut x = 0u64;
        let mut ys = 0u64;
        while g == 1 {
            x = y;
            for _ in 0..r {
                y = f(y);
            }
            let mut k = 0u64;
            while k < r && g == 1 {
                ys = y;
                let lim = core::cmp::min(128, r - k);
                for _ in 0..lim {
                    y = f(y);
                    let diff = x.abs_diff(y);
                    if diff != 0 {
                        q = mulmod(q, diff, n);
                    }
                }
                g = gcd_u64(q, n);
                k += lim;
            }
            r *= 2;
        }
        if g == n {
            loop {
                ys = f(ys);
                g = gcd_u64(x.abs_diff(ys), n);
                if g != 1 {
                    break;
                }
            }
        }
        if g != n && g != 1 {
            return g;
        }
        c += 1;
    }
}

/// Fully factor `n` (which must fit in `u64`) into its prime factors, appended
/// to `out` (unsorted; caller sorts).
pub fn factor_u64(mut n: u64, out: &mut Vec<u64>) {
    if n <= 1 {
        return;
    }
    for &p in small_primes() {
        let p = p as u64;
        if p.saturating_mul(p) > n {
            break;
        }
        while n.is_multiple_of(p) {
            out.push(p);
            n /= p;
        }
    }
    if n > 1 {
        factor_rec(n, out);
    }
}

fn factor_rec(n: u64, out: &mut Vec<u64>) {
    if n == 1 {
        return;
    }
    if is_prime_u64(n) {
        out.push(n);
        return;
    }
    let d = pollard_brent(n);
    factor_rec(d, out);
    factor_rec(n / d, out);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primality_matches_known() {
        for p in [2u64, 3, 97, 104729, 1_000_003, (1u64 << 61) - 1] {
            assert!(is_prime_u64(p), "{p} should be prime");
        }
        for c in [1u64, 4, 91, 561, 1105, 1_000_004, (1u64 << 61) + 1] {
            assert!(!is_prime_u64(c), "{c} should be composite");
        }
    }

    #[test]
    fn factors_semiprimes_and_powers() {
        let cases: [u64; 6] = [
            360,
            1_000_003 * 1_000_033,
            3_366_523_423 * 3_390_825_731,
            2u64.pow(20) * 3u64.pow(5),
            18_446_744_073_709_551_557, // largest prime < 2^64
            15_485_863 * 15_485_867,
        ];
        for n in cases {
            let mut f = Vec::new();
            factor_u64(n, &mut f);
            let prod: u128 = f.iter().map(|&x| x as u128).product();
            assert_eq!(prod, n as u128, "product mismatch for {n}");
            for &p in &f {
                assert!(is_prime_u64(p), "non-prime factor {p} of {n}");
            }
        }
    }
}
