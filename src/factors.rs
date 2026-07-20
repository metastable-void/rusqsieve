use crate::Natural;
use core::num::NonZero;
use std::collections::BTreeMap;

/// Sorted probable-prime factors and their multiplicities.
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
    pub fn new() -> Self {
        Self {
            map: BTreeMap::new(),
        }
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (&Natural<P>, NonZero<usize>)> {
        self.map.iter().map(|(n, e)| (n, *e))
    }
    pub fn get(&self, p: &Natural<P>) -> Option<NonZero<usize>> {
        self.map.get(p).copied()
    }
    pub fn len(&self) -> usize {
        self.map.len()
    }
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
    pub fn into_map(self) -> BTreeMap<Natural<P>, NonZero<usize>> {
        self.map
    }
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
