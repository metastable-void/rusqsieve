//! Strong probable-prime testing.
use crate::Natural;
use core::num::NonZero;

#[cfg(test)]
static LUCAS_TEST_CALLS: core::sync::atomic::AtomicUsize = core::sync::atomic::AtomicUsize::new(0);

#[derive(Clone, Debug)]
pub struct PrimalityConfig {
    pub rounds: NonZero<u32>,
    pub witnesses: WitnessPolicy,
}
#[derive(Clone, Debug)]
pub enum WitnessPolicy {
    FirstPrimes,
    Seeded { seed: [u8; 32] },
}
impl Default for PrimalityConfig {
    fn default() -> Self {
        Self {
            rounds: NonZero::new(16).unwrap(),
            witnesses: WitnessPolicy::FirstPrimes,
        }
    }
}

const SMALL: [u64; 32] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131,
];
pub fn is_probable_prime<const P: usize>(n: &Natural<P>, config: &PrimalityConfig) -> bool {
    if *n < Natural::from_u64(2) {
        return false;
    }
    for &p in &SMALL {
        let q = Natural::from_u64(p);
        if n == &q {
            return true;
        }
        if n.mod_u64(p) == 0 {
            return false;
        }
    }
    if n.is_even() {
        return false;
    }
    if n.bit_len() > 64 {
        if n.is_square() || !miller_rabin_witness(n, 2) || !strong_lucas_selfridge(n) {
            return false;
        }
    }
    let one = Natural::ONE;
    let nm1 = n.checked_sub(&one).unwrap();
    let s = nm1.trailing_zeros();
    let d = nm1.clone() >> s;
    let mut rng = seed_state(config, n);
    for round in 0..config.rounds.get() {
        let a = match config.witnesses {
            WitnessPolicy::FirstPrimes => Natural::from_u64(SMALL[round as usize % SMALL.len()]),
            WitnessPolicy::Seeded { .. } => {
                rng = xorshift(rng);
                Natural::from_u64(2u64.wrapping_add(rng))
            }
        };
        let a = a.div_rem(n).unwrap().1;
        if a.is_zero() {
            continue;
        }
        let mut x = a.pow_mod(&d, n);
        if x == one || x == nm1 {
            continue;
        }
        let mut composite = true;
        for _ in 1..s {
            x = x.mul_mod(&x, n);
            if x == nm1 {
                composite = false;
                break;
            }
            if x == one {
                return false;
            }
        }
        if composite {
            return false;
        }
    }
    true
}

fn miller_rabin_witness<const P: usize>(n: &Natural<P>, witness: u64) -> bool {
    let one = Natural::ONE;
    let nm1 = n.checked_sub(&one).unwrap();
    let s = nm1.trailing_zeros();
    let d = nm1.clone() >> s;
    let mut x = Natural::from_u64(witness).pow_mod(&d, n);
    if x == one || x == nm1 {
        return true;
    }
    for _ in 1..s {
        x = x.mul_mod(&x, n);
        if x == nm1 {
            return true;
        }
        if x == one {
            return false;
        }
    }
    false
}

fn strong_lucas_selfridge<const P: usize>(n: &Natural<P>) -> bool {
    #[cfg(test)]
    LUCAS_TEST_CALLS.fetch_add(1, core::sync::atomic::Ordering::Relaxed);

    let mut magnitude = 5i64;
    let mut positive = true;
    let d = loop {
        let candidate = if positive { magnitude } else { -magnitude };
        match jacobi_small_natural(candidate, n) {
            -1 => break candidate,
            0 => return false,
            _ => {}
        }
        magnitude += 2;
        positive = !positive;
    };
    let q = (1 - d) / 4;
    let n_plus_one = n.checked_add(&Natural::ONE).unwrap();
    let s = n_plus_one.trailing_zeros();
    let odd_part = n_plus_one >> s;
    let (u, mut v, mut qk) = lucas_sequence(n, d, q, &odd_part);
    if u.is_zero() || v.is_zero() {
        return true;
    }
    for _ in 1..s {
        v = sub_mod(&v.mul_mod(&v, n), &qk.add_mod(&qk, n), n);
        qk = qk.mul_mod(&qk, n);
        if v.is_zero() {
            return true;
        }
    }
    false
}

fn lucas_sequence<const P: usize>(
    n: &Natural<P>,
    d: i64,
    q: i64,
    k: &Natural<P>,
) -> (Natural<P>, Natural<P>, Natural<P>) {
    debug_assert!(!k.is_zero() && k.is_odd());
    let mut u = Natural::ONE;
    let mut v = Natural::ONE;
    let q_mod = signed_small_mod(q, n);
    let mut qk = q_mod.clone();
    let bits = k.bit_len();
    for bit in (0..bits.saturating_sub(1)).rev() {
        u = u.mul_mod(&v, n);
        v = sub_mod(&v.mul_mod(&v, n), &qk.add_mod(&qk, n), n);
        qk = qk.mul_mod(&qk, n);
        if natural_bit(k, bit) {
            let old_u = u;
            let old_v = v;
            u = half_mod(&old_u.add_mod(&old_v, n), n);
            let du = signed_mul_mod(d, &old_u, n);
            v = half_mod(&du.add_mod(&old_v, n), n);
            qk = qk.mul_mod(&q_mod, n);
        }
    }
    (u, v, qk)
}

fn jacobi_small_natural<const P: usize>(d: i64, n: &Natural<P>) -> i32 {
    debug_assert!(d != 0 && d & 1 != 0 && n.is_odd());
    let magnitude = d.unsigned_abs();
    let mut sign = 1;
    let n_mod_four = n.mod_u64(4);
    if d < 0 && n_mod_four == 3 {
        sign = -sign;
    }
    if magnitude & 3 == 3 && n_mod_four == 3 {
        sign = -sign;
    }
    sign * i32::from(crate::jacobi_u64(n.mod_u64(magnitude), magnitude))
}

fn natural_bit<const P: usize>(n: &Natural<P>, bit: usize) -> bool {
    n.as_parts()
        .get(bit / 64)
        .is_some_and(|limb| limb & (1u64 << (bit % 64)) != 0)
}

fn signed_small_mod<const P: usize>(value: i64, n: &Natural<P>) -> Natural<P> {
    let magnitude = Natural::from_u64(value.unsigned_abs())
        .div_rem(n)
        .unwrap()
        .1;
    if value < 0 && !magnitude.is_zero() {
        n.wrapping_sub(&magnitude)
    } else {
        magnitude
    }
}

fn signed_mul_mod<const P: usize>(value: i64, rhs: &Natural<P>, n: &Natural<P>) -> Natural<P> {
    let product = Natural::from_u64(value.unsigned_abs()).mul_mod(rhs, n);
    if value < 0 && !product.is_zero() {
        n.wrapping_sub(&product)
    } else {
        product
    }
}

fn sub_mod<const P: usize>(a: &Natural<P>, b: &Natural<P>, n: &Natural<P>) -> Natural<P> {
    if a >= b {
        a.wrapping_sub(b)
    } else {
        n.wrapping_sub(&b.wrapping_sub(a))
    }
}

fn half_mod<const P: usize>(value: &Natural<P>, n: &Natural<P>) -> Natural<P> {
    if value.is_even() {
        value.clone() >> 1
    } else {
        (value.clone() >> 1).add_mod(&((n.clone() >> 1).wrapping_add(&Natural::ONE)), n)
    }
}
fn seed_state<const P: usize>(c: &PrimalityConfig, n: &Natural<P>) -> u64 {
    let mut s = 0x9e3779b97f4a7c15;
    for &x in n.as_parts() {
        s ^= x;
        s = xorshift(s)
    }
    if let WitnessPolicy::Seeded { seed } = c.witnesses {
        for chunk in seed.chunks_exact(8) {
            s ^= u64::from_le_bytes(chunk.try_into().unwrap());
            s = xorshift(s)
        }
    }
    s
}
fn xorshift(mut x: u64) -> u64 {
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn known() {
        let c = PrimalityConfig::default();
        for p in [2, 3, 97, 104729] {
            assert!(is_probable_prime(&Natural::<2>::from_u64(p), &c))
        }
        for n in [0, 1, 4, 91, 561, 1105] {
            assert!(!is_probable_prime(&Natural::<2>::from_u64(n), &c))
        }
    }

    #[test]
    fn rejects_strong_pseudoprimes_and_exercises_lucas() {
        LUCAS_TEST_CALLS.store(0, core::sync::atomic::Ordering::Relaxed);
        let c = PrimalityConfig::default();
        for composite in ["318665857834031151167461", "3317044064679887385961981"] {
            let n = Natural::<16>::from_decimal(composite).unwrap();
            assert!(!is_probable_prime(&n, &c), "{composite} is composite");
        }
        let prime = Natural::<16>::from_decimal("170141183460469231731687303715884105727").unwrap();
        assert!(is_probable_prime(&prime, &c));
        assert!(
            LUCAS_TEST_CALLS.load(core::sync::atomic::Ordering::Relaxed) > 0,
            "the strong Lucas stage was not exercised"
        );
    }
}
