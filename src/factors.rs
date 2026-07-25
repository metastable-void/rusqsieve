use crate::Natural;
use core::num::NonZero;
use std::collections::BTreeMap;

/// Sorted probable-prime factors and their multiplicities.
///
/// Values are produced by [`factor`](crate::factor). The iterator order is
/// ascending. Each entry stores a distinct factor and its nonzero
/// multiplicity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PrimeFactors<const P: usize = 16> {
    map: BTreeMap<Natural<P>, NonZero<usize>>,
}
impl<const P: usize> Default for PrimeFactors<P> {
    fn default() -> Self {
        Self::new()
    }
}
impl<const P: usize> PrimeFactors<P> {
    pub(crate) fn new() -> Self {
        Self {
            map: BTreeMap::new(),
        }
    }

    /// Iterates over distinct factors and their multiplicities in ascending
    /// order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (&Natural<P>, NonZero<usize>)> {
        self.map.iter().map(|(n, e)| (n, *e))
    }

    /// Returns the multiplicity of `factor`, or `None` when it is absent.
    pub fn multiplicity(&self, factor: &Natural<P>) -> Option<NonZero<usize>> {
        self.map.get(factor).copied()
    }

    /// Iterates over factors in ascending order, repeating each according to
    /// its multiplicity.
    pub fn expanded(&self) -> impl Iterator<Item = &Natural<P>> {
        self.map
            .iter()
            .flat_map(|(factor, exponent)| core::iter::repeat_n(factor, exponent.get()))
    }

    /// Returns the number of distinct prime factors.
    pub fn distinct_len(&self) -> usize {
        self.map.len()
    }

    /// Returns the number of prime factors including multiplicity.
    pub fn total_len(&self) -> usize {
        self.map.values().map(|exponent| exponent.get()).sum()
    }

    /// Returns `true` when this is the empty factorization of one.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Checks that multiplying all factors with their multiplicities
    /// reconstructs `original` without overflowing its capacity.
    pub fn verify_product(&self, original: &Natural<P>) -> bool {
        let mut out = Natural::ONE;
        for (p, e) in self.iter() {
            for _ in 0..e.get() {
                let Some(v) = out.checked_mul(p) else {
                    return false;
                };
                out = v
            }
        }
        &out == original
    }

    #[allow(dead_code)]
    pub(crate) fn get(&self, p: &Natural<P>) -> Option<NonZero<usize>> {
        self.map.get(p).copied()
    }
    pub(crate) fn insert_count(&mut self, p: Natural<P>, count: usize) {
        if count == 0 {
            return;
        }
        let old = self.map.get(&p).map_or(0, |x| x.get());
        self.map.insert(
            p,
            NonZero::new(old.checked_add(count).expect("factor exponent overflow")).unwrap(),
        );
    }
}
