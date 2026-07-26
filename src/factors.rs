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

/// Iterator over distinct factors and their multiplicities.
pub struct PrimeFactorIter<'a, const P: usize> {
    inner: std::collections::btree_map::Iter<'a, Natural<P>, NonZero<usize>>,
}

impl<'a, const P: usize> Iterator for PrimeFactorIter<'a, P> {
    type Item = (&'a Natural<P>, NonZero<usize>);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|(factor, exponent)| (factor, *exponent))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<const P: usize> ExactSizeIterator for PrimeFactorIter<'_, P> {}

/// Iterator over factors repeated according to multiplicity.
pub struct ExpandedPrimeFactors<'a, const P: usize> {
    outer: std::collections::btree_map::Iter<'a, Natural<P>, NonZero<usize>>,
    current: Option<(&'a Natural<P>, usize)>,
}

impl<'a, const P: usize> Iterator for ExpandedPrimeFactors<'a, P> {
    type Item = &'a Natural<P>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some((factor, remaining)) = self.current
                && remaining != 0
            {
                self.current = Some((factor, remaining - 1));
                return Some(factor);
            }
            let (factor, exponent) = self.outer.next()?;
            self.current = Some((factor, exponent.get()));
        }
    }
}

impl<const P: usize> PrimeFactors<P> {
    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        Self {
            map: BTreeMap::new(),
        }
    }

    /// Iterates over distinct factors and their multiplicities in ascending
    /// order.
    #[must_use]
    pub fn iter(&self) -> PrimeFactorIter<'_, P> {
        PrimeFactorIter {
            inner: self.map.iter(),
        }
    }

    /// Returns the multiplicity of `factor`, or `None` when it is absent.
    #[must_use]
    pub fn multiplicity(&self, factor: &Natural<P>) -> Option<NonZero<usize>> {
        self.map.get(factor).copied()
    }

    /// Iterates over factors in ascending order, repeating each according to
    /// its multiplicity.
    #[must_use]
    pub fn expanded(&self) -> ExpandedPrimeFactors<'_, P> {
        ExpandedPrimeFactors {
            outer: self.map.iter(),
            current: None,
        }
    }

    /// Returns the number of distinct prime factors.
    #[must_use]
    pub fn distinct_len(&self) -> usize {
        self.map.len()
    }

    /// Returns the number of prime factors including multiplicity.
    #[must_use]
    pub fn total_len(&self) -> usize {
        self.map.values().map(|exponent| exponent.get()).sum()
    }

    /// Returns `true` when this is the empty factorization of one.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Checks that multiplying all factors with their multiplicities
    /// reconstructs `original` without overflowing its capacity.
    #[must_use]
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
    #[allow(dead_code)]
    pub(crate) fn insert_count(&mut self, p: Natural<P>, count: usize) {
        if count == 0 {
            return;
        }
        let old = self.map.get(&p).map_or(0, |x| x.get());
        let combined = old.saturating_add(count);
        debug_assert!(old.checked_add(count).is_some(), "factor exponent overflow");
        self.map.insert(p, NonZero::new(combined).unwrap());
    }
}

impl<'a, const P: usize> IntoIterator for &'a PrimeFactors<P> {
    type Item = (&'a Natural<P>, NonZero<usize>);
    type IntoIter = PrimeFactorIter<'a, P>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
