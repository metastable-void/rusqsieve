//! Strong probable-prime testing.
use crate::Natural;
use core::num::NonZero;

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
}
