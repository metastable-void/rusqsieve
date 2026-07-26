use super::*;

impl EngineJobResult {
    /// Serialize this family's relations for transport to a coordinator (e.g. from a
    /// Web Worker back to the main thread). Format is little-endian:
    /// `family:u64, polynomials:u64, count:u32`, then per relation
    /// `root_len:u16, root[root_len], sign:u8,
    /// large:{tag:u8, 0/1/2 × u64}, powers_len:u32, [index:u32, exp:u16]…`.
    pub fn to_bytes(&self) -> Vec<u8> {
        let capacity = 20
            + self
                .inner
                .relations
                .iter()
                .map(|relation| {
                    2 + relation.root.bit_len().div_ceil(8)
                        + 2
                        + match relation.large {
                            LargePrime::None => 0,
                            LargePrime::One(_) => 8,
                            LargePrime::Two(_, _) => 16,
                        }
                        + 4
                        + relation.powers.len() * 6
                })
                .sum::<usize>();
        let mut v = Vec::with_capacity(capacity);
        v.extend_from_slice(&self.inner.family.to_le_bytes());
        v.extend_from_slice(&self.inner.polynomials.to_le_bytes());
        v.extend_from_slice(&(self.inner.relations.len() as u32).to_le_bytes());
        for r in &self.inner.relations {
            let mut root = [0u8; PARTS * 8];
            let root_len = r.root.write_le_bytes(&mut root).unwrap();
            v.extend_from_slice(&(root_len as u16).to_le_bytes());
            v.extend_from_slice(&root[..root_len]);
            v.push(r.sign as u8);
            match r.large {
                LargePrime::None => v.push(0),
                LargePrime::One(a) => {
                    v.push(1);
                    v.extend_from_slice(&a.to_le_bytes());
                }
                LargePrime::Two(a, b) => {
                    v.push(2);
                    v.extend_from_slice(&a.to_le_bytes());
                    v.extend_from_slice(&b.to_le_bytes());
                }
            }
            v.extend_from_slice(&(r.powers.len() as u32).to_le_bytes());
            for &(i, e) in &r.powers {
                v.extend_from_slice(&i.to_le_bytes());
                v.extend_from_slice(&e.to_le_bytes());
            }
        }
        v
    }
}

/// Inverse of [`EngineJobResult::to_bytes`].
pub(super) fn deserialize_family(b: &[u8]) -> Option<FamilyResult> {
    struct Cur<'a> {
        b: &'a [u8],
        o: usize,
    }
    impl Cur<'_> {
        fn take(&mut self, n: usize) -> Option<&[u8]> {
            let s = self.b.get(self.o..self.o + n)?;
            self.o += n;
            Some(s)
        }
        fn u8(&mut self) -> Option<u8> {
            Some(self.take(1)?[0])
        }
        fn u16(&mut self) -> Option<u16> {
            Some(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
        }
        fn u32(&mut self) -> Option<u32> {
            Some(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
        }
        fn u64(&mut self) -> Option<u64> {
            Some(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
        }
    }
    let mut c = Cur { b, o: 0 };
    let family = c.u64()?;
    let polynomials = c.u64()?;
    let count = c.u32()? as usize;
    let mut relations = Vec::with_capacity(count.min(1 << 20));
    for _ in 0..count {
        let root_len = c.u16()? as usize;
        if root_len > PARTS * 8 {
            return None;
        }
        let root = Natural::from_le_bytes(c.take(root_len)?).ok()?;
        let sign = c.u8()? != 0;
        let large = match c.u8()? {
            0 => LargePrime::None,
            1 => LargePrime::One(c.u64()?),
            2 => LargePrime::Two(c.u64()?, c.u64()?),
            _ => return None,
        };
        let plen = c.u32()? as usize;
        let mut powers = Vec::with_capacity(plen.min(1 << 16));
        for _ in 0..plen {
            let i = c.u32()?;
            let e = c.u16()?;
            powers.push((i, e));
        }
        relations.push(Relation {
            root,
            sign,
            powers,
            large,
        });
    }
    Some(FamilyResult {
        family,
        polynomials,
        relations,
        survivors: 0,
    })
}
