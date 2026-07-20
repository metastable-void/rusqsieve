//! Fixed-capacity unsigned integer arithmetic.

use core::cmp::Ordering;
use core::fmt;
use core::ops::*;
use core::str::FromStr;

/// A decimal parsing failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParseNaturalError {
    Empty,
    InvalidDigit { index: usize, byte: u8 },
    Overflow,
}

impl fmt::Display for ParseNaturalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty integer"),
            Self::InvalidDigit { index, byte } => {
                write!(f, "invalid decimal byte {byte:#x} at index {index}")
            }
            Self::Overflow => f.write_str("integer exceeds Natural capacity"),
        }
    }
}
impl std::error::Error for ParseNaturalError {}

/// Input bytes do not fit the selected capacity.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CapacityError;
impl fmt::Display for CapacityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("integer exceeds Natural capacity")
    }
}
impl std::error::Error for CapacityError {}

/// A serialization destination was too small.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BufferTooSmall {
    pub required: usize,
    pub available: usize,
}
impl fmt::Display for BufferTooSmall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "buffer has {} bytes; {} required",
            self.available, self.required
        )
    }
}
impl std::error::Error for BufferTooSmall {}

/// A fixed-capacity unsigned integer with little-endian 64-bit limbs.
#[repr(transparent)]
#[derive(Clone, Eq, PartialEq, Hash)]
pub struct Natural<const PARTS_64: usize = 16> {
    parts: [u64; PARTS_64],
}

impl<const P: usize> Natural<P> {
    pub const BITS: usize = P * 64;
    pub const ZERO: Self = Self { parts: [0; P] };
    pub const ONE: Self = Self::from_u64(1);
    pub const MAX: Self = Self {
        parts: [u64::MAX; P],
    };

    pub const fn from_u64(value: u64) -> Self {
        let mut parts = [0; P];
        if P != 0 {
            parts[0] = value;
        }
        Self { parts }
    }
    pub const fn as_parts(&self) -> &[u64; P] {
        &self.parts
    }
    /// The value as a `u64` if it fits (all limbs above the lowest are zero).
    pub fn to_u64(&self) -> Option<u64> {
        if P == 0 {
            return Some(0);
        }
        self.parts[1..].iter().all(|&x| x == 0).then_some(self.parts[0])
    }
    pub fn as_mut_parts(&mut self) -> &mut [u64; P] {
        &mut self.parts
    }
    pub fn is_zero(&self) -> bool {
        self.parts.iter().all(|&x| x == 0)
    }
    pub fn is_one(&self) -> bool {
        P != 0 && self.parts[0] == 1 && self.parts[1..].iter().all(|&x| x == 0)
    }
    pub fn is_even(&self) -> bool {
        P == 0 || self.parts[0] & 1 == 0
    }
    pub fn is_odd(&self) -> bool {
        !self.is_even()
    }
    pub fn bit_len(&self) -> usize {
        self.parts.iter().rposition(|&x| x != 0).map_or(0, |i| {
            i * 64 + (64 - self.parts[i].leading_zeros() as usize)
        })
    }
    pub fn trailing_zeros(&self) -> usize {
        self.parts
            .iter()
            .position(|&x| x != 0)
            .map_or(Self::BITS, |i| {
                i * 64 + self.parts[i].trailing_zeros() as usize
            })
    }
    pub fn bit(&self, index: usize) -> bool {
        index < Self::BITS && (self.parts[index / 64] >> (index % 64)) & 1 != 0
    }

    pub const fn from_decimal(value: &str) -> Result<Self, ParseNaturalError> {
        let bytes = value.as_bytes();
        if bytes.is_empty() {
            return Err(ParseNaturalError::Empty);
        }
        let mut out = Self::ZERO;
        let mut i = 0;
        while i < bytes.len() {
            let b = bytes[i];
            if b < b'0' || b > b'9' {
                return Err(ParseNaturalError::InvalidDigit { index: i, byte: b });
            }
            let mut carry = (b - b'0') as u128;
            let mut limb = 0;
            while limb < P {
                let v = out.parts[limb] as u128 * 10 + carry;
                out.parts[limb] = v as u64;
                carry = v >> 64;
                limb += 1;
            }
            if carry != 0 {
                return Err(ParseNaturalError::Overflow);
            }
            i += 1;
        }
        Ok(out)
    }

    pub fn from_be_bytes(bytes: &[u8]) -> Result<Self, CapacityError> {
        if bytes.len() > P * 8 && bytes[..bytes.len() - P * 8].iter().any(|&x| x != 0) {
            return Err(CapacityError);
        }
        let mut out = Self::ZERO;
        for (i, &b) in bytes.iter().rev().take(P * 8).enumerate() {
            out.parts[i / 8] |= (b as u64) << (8 * (i % 8));
        }
        Ok(out)
    }
    pub fn from_le_bytes(bytes: &[u8]) -> Result<Self, CapacityError> {
        if bytes.len() > P * 8 && bytes[P * 8..].iter().any(|&x| x != 0) {
            return Err(CapacityError);
        }
        let mut out = Self::ZERO;
        for (i, &b) in bytes.iter().take(P * 8).enumerate() {
            out.parts[i / 8] |= (b as u64) << (8 * (i % 8));
        }
        Ok(out)
    }
    fn byte_len(&self) -> usize {
        self.bit_len().div_ceil(8)
    }
    pub fn write_be_bytes(&self, out: &mut [u8]) -> Result<usize, BufferTooSmall> {
        let n = self.byte_len();
        if out.len() < n {
            return Err(BufferTooSmall {
                required: n,
                available: out.len(),
            });
        }
        for (i, slot) in out[..n].iter_mut().enumerate() {
            let j = n - 1 - i;
            *slot = (self.parts[j / 8] >> (8 * (j % 8))) as u8;
        }
        Ok(n)
    }
    pub fn write_le_bytes(&self, out: &mut [u8]) -> Result<usize, BufferTooSmall> {
        let n = self.byte_len();
        if out.len() < n {
            return Err(BufferTooSmall {
                required: n,
                available: out.len(),
            });
        }
        for (i, slot) in out[..n].iter_mut().enumerate() {
            *slot = (self.parts[i / 8] >> (8 * (i % 8))) as u8;
        }
        Ok(n)
    }

    pub fn overflowing_add(&self, rhs: &Self) -> (Self, bool) {
        let mut out = Self::ZERO;
        let mut carry = 0u128;
        for i in 0..P {
            let v = self.parts[i] as u128 + rhs.parts[i] as u128 + carry;
            out.parts[i] = v as u64;
            carry = v >> 64;
        }
        (out, carry != 0)
    }
    pub fn overflowing_sub(&self, rhs: &Self) -> (Self, bool) {
        let mut out = Self::ZERO;
        let mut borrow = false;
        for i in 0..P {
            let (a, b1) = self.parts[i].overflowing_sub(rhs.parts[i]);
            let (v, b2) = a.overflowing_sub(borrow as u64);
            out.parts[i] = v;
            borrow = b1 || b2;
        }
        (out, borrow)
    }
    pub fn widening_mul(&self, rhs: &Self) -> WideNatural<P> {
        let mut low = [0u64; P];
        let mut high = [0u64; P];
        // Only iterate over significant limbs so arithmetic cost scales with the
        // operands' actual magnitude, not the fixed capacity `P`.
        let alen = sig_len(&self.parts);
        let blen = sig_len(&rhs.parts);
        for i in 0..alen {
            let ai = self.parts[i] as u128;
            if ai == 0 {
                continue;
            }
            let mut carry = 0u128;
            for j in 0..blen {
                let k = i + j;
                let old = if k < P { low[k] } else { high[k - P] };
                let v = ai * rhs.parts[j] as u128 + old as u128 + carry;
                if k < P {
                    low[k] = v as u64;
                } else {
                    high[k - P] = v as u64;
                }
                carry = v >> 64;
            }
            let mut k = i + blen;
            while carry != 0 && k < 2 * P {
                let old = if k < P { low[k] } else { high[k - P] };
                let v = old as u128 + carry;
                if k < P {
                    low[k] = v as u64;
                } else {
                    high[k - P] = v as u64;
                }
                carry = v >> 64;
                k += 1;
            }
        }
        WideNatural { low, high }
    }
    pub fn widening_square(&self) -> WideNatural<P> {
        self.widening_mul(self)
    }
    pub fn overflowing_mul(&self, rhs: &Self) -> (Self, bool) {
        self.widening_mul(rhs).overflowing_narrow()
    }
    pub fn checked_add(&self, rhs: &Self) -> Option<Self> {
        let (v, o) = self.overflowing_add(rhs);
        (!o).then_some(v)
    }
    pub fn checked_sub(&self, rhs: &Self) -> Option<Self> {
        let (v, o) = self.overflowing_sub(rhs);
        (!o).then_some(v)
    }
    pub fn checked_mul(&self, rhs: &Self) -> Option<Self> {
        let (v, o) = self.overflowing_mul(rhs);
        (!o).then_some(v)
    }
    pub fn wrapping_add(&self, rhs: &Self) -> Self {
        self.overflowing_add(rhs).0
    }
    pub fn wrapping_sub(&self, rhs: &Self) -> Self {
        self.overflowing_sub(rhs).0
    }
    pub fn wrapping_mul(&self, rhs: &Self) -> Self {
        self.overflowing_mul(rhs).0
    }

    pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)> {
        let dtop = divisor.parts.iter().rposition(|&x| x != 0)?;
        // Single-limb divisor: use the dedicated fast path.
        if dtop == 0 {
            let (q, r) = self.div_rem_u64(divisor.parts[0]).unwrap();
            return Some((q, Self::from_u64(r)));
        }
        let (q, r) = knuth_divmod(&self.parts, &divisor.parts);
        let mut qn = Self::ZERO;
        for (slot, &x) in qn.parts.iter_mut().zip(q.iter()) {
            *slot = x;
        }
        let mut rn = Self::ZERO;
        for (slot, &x) in rn.parts.iter_mut().zip(r.iter()) {
            *slot = x;
        }
        Some((qn, rn))
    }
    pub fn div_rem_u64(&self, divisor: u64) -> Option<(Self, u64)> {
        if divisor == 0 {
            return None;
        }
        let mut q = Self::ZERO;
        let Some(top) = self.parts.iter().rposition(|&x| x != 0) else {
            return Some((q, 0));
        };
        let d = divisor as u128;
        let mut r = 0u128;
        for i in (0..=top).rev() {
            let v = (r << 64) | self.parts[i] as u128;
            let qi = (v / d) as u64;
            r = v - qi as u128 * d;
            q.parts[i] = qi;
        }
        Some((q, r as u64))
    }
    /// Binary (Stein) GCD: shifts and subtraction only, no division.
    pub fn gcd(&self, rhs: &Self) -> Self {
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }
        let mut a = self.clone();
        let mut b = rhs.clone();
        let shift = a.trailing_zeros().min(b.trailing_zeros());
        a >>= a.trailing_zeros();
        loop {
            b >>= b.trailing_zeros();
            if a > b {
                core::mem::swap(&mut a, &mut b);
            }
            b -= &a;
            if b.is_zero() {
                break;
            }
        }
        a << shift
    }
    pub fn extended_gcd(&self, rhs: &Self) -> ExtendedGcdResult<P> {
        ExtendedGcdResult { gcd: self.gcd(rhs) }
    }

    pub fn sqrt_rem(&self) -> (Self, Self) {
        if self.is_zero() {
            return (Self::ZERO, Self::ZERO);
        }
        let mut x = Self::ONE << self.bit_len().div_ceil(2);
        if x.is_zero() {
            x = Self::MAX;
        }
        loop {
            let q = self.div_rem(&x).unwrap().0;
            let sum = x.checked_add(&q).expect("square-root iterate fits");
            let next = sum >> 1usize;
            if next >= x {
                let sq = x.checked_mul(&x).unwrap();
                return (x, self.checked_sub(&sq).unwrap());
            }
            x = next;
        }
    }
    pub fn floor_sqrt(&self) -> Self {
        self.sqrt_rem().0
    }
    pub fn ceil_sqrt(&self) -> Self {
        let (s, r) = self.sqrt_rem();
        if r.is_zero() {
            s
        } else {
            s.checked_add(&Self::ONE).unwrap()
        }
    }
    pub fn is_square(&self) -> bool {
        self.sqrt_rem().1.is_zero()
    }
    pub fn checked_pow_u32(&self, mut exponent: u32) -> Option<Self> {
        let mut a = self.clone();
        let mut out = Self::ONE;
        while exponent != 0 {
            if exponent & 1 != 0 {
                out = out.checked_mul(&a)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                a = a.checked_mul(&a)?;
            }
        }
        Some(out)
    }
    pub fn perfect_power(&self) -> Option<(Self, u32)> {
        if self < &Self::from_u64(4) {
            return None;
        }
        let root = self.floor_sqrt();
        if root.checked_mul(&root).as_ref() == Some(self) {
            return Some((root, 2));
        }
        let bits = self.bit_len() as u32;
        for e in 3..=bits {
            if !is_prime_u32(e) {
                continue;
            }
            let mut lo = Self::ONE;
            let mut hi = Self::ONE << self.bit_len().div_ceil(e as usize);
            hi = hi.checked_add(&Self::ONE).unwrap_or(Self::MAX);
            while lo < hi {
                let sum = lo.checked_add(&hi)?;
                let mid = (sum + Self::ONE) >> 1usize;
                match mid.checked_pow_u32(e) {
                    Some(v) if v <= *self => lo = mid,
                    _ => hi = mid.wrapping_sub(&Self::ONE),
                }
            }
            if lo.checked_pow_u32(e).as_ref() == Some(self) {
                return Some((lo, e));
            }
        }
        None
    }

    pub(crate) fn add_mod(&self, rhs: &Self, m: &Self) -> Self {
        debug_assert!(self < m && rhs < m);
        let threshold = m.wrapping_sub(rhs);
        if self >= &threshold {
            self.wrapping_sub(&threshold)
        } else {
            self.wrapping_add(rhs)
        }
    }
    pub(crate) fn mul_mod(&self, rhs: &Self, m: &Self) -> Self {
        let a = self.div_rem(m).unwrap().1;
        let b = rhs.div_rem(m).unwrap().1;
        a.widening_mul(&b).rem_natural(m)
    }
    pub(crate) fn pow_mod(&self, e: &Self, m: &Self) -> Self {
        let mut a = self.div_rem(m).unwrap().1;
        let mut e = e.clone();
        let mut out = Self::ONE.div_rem(m).unwrap().1;
        while !e.is_zero() {
            if e.is_odd() {
                out = out.mul_mod(&a, m)
            }
            e >>= 1usize;
            if !e.is_zero() {
                a = a.mul_mod(&a, m)
            }
        }
        out
    }
    pub(crate) fn mod_u64(&self, m: u64) -> u64 {
        self.div_rem_u64(m).unwrap().1
    }
}

fn is_prime_u32(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 1
    }
    true
}

/// Number of significant (nonzero through) limbs in a little-endian slice.
#[inline]
fn sig_len(limbs: &[u64]) -> usize {
    limbs.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1)
}

/// Little-endian schoolbook long division (Knuth TAOCP Alg. D) on limb slices.
/// Returns `(quotient, remainder)` as trimmed little-endian limb vectors.
/// `den` must be nonzero. This is the shared normalized-long-division primitive
/// used by both `Natural::div_rem` and wide-product reduction (SPEC §6.11).
fn knuth_divmod(num: &[u64], den: &[u64]) -> (Vec<u64>, Vec<u64>) {
    const LOW: u128 = 0xffff_ffff_ffff_ffff;
    let num_len = sig_len(num);
    let n = sig_len(den);
    debug_assert!(n >= 1, "division by zero");
    if num_len < n {
        return (vec![0], num[..num_len.max(1)].to_vec());
    }
    // Single-limb divisor: straight long division.
    if n == 1 {
        let d = den[0] as u128;
        let mut q = vec![0u64; num_len];
        let mut r = 0u128;
        for i in (0..num_len).rev() {
            let cur = (r << 64) | num[i] as u128;
            let qi = (cur / d) as u64;
            r = cur - qi as u128 * d;
            q[i] = qi;
        }
        return (q, vec![r as u64]);
    }
    let m = num_len - n;
    let shift = den[n - 1].leading_zeros();
    // Normalize divisor and dividend so the divisor's top bit is set.
    let mut vn = vec![0u64; n];
    if shift == 0 {
        vn.copy_from_slice(&den[..n]);
    } else {
        for i in (1..n).rev() {
            vn[i] = (den[i] << shift) | (den[i - 1] >> (64 - shift));
        }
        vn[0] = den[0] << shift;
    }
    let mut un = vec![0u64; num_len + 1];
    if shift == 0 {
        un[..num_len].copy_from_slice(&num[..num_len]);
    } else {
        un[num_len] = num[num_len - 1] >> (64 - shift);
        for i in (1..num_len).rev() {
            un[i] = (num[i] << shift) | (num[i - 1] >> (64 - shift));
        }
        un[0] = num[0] << shift;
    }
    let mut q = vec![0u64; m + 1];
    let base = 1u128 << 64;
    for j in (0..=m).rev() {
        let top = ((un[j + n] as u128) << 64) | un[j + n - 1] as u128;
        let mut qhat = top / vn[n - 1] as u128;
        let mut rhat = top - qhat * vn[n - 1] as u128;
        while qhat >= base
            || qhat * vn[n - 2] as u128 > (rhat << 64) | un[j + n - 2] as u128
        {
            qhat -= 1;
            rhat += vn[n - 1] as u128;
            if rhat >= base {
                break;
            }
        }
        // Multiply and subtract qhat*divisor from the current window.
        let mut carry: u128 = 0;
        let mut borrow: i128 = 0;
        for i in 0..n {
            let p = qhat * vn[i] as u128 + carry;
            carry = p >> 64;
            let t = un[j + i] as i128 - borrow - (p & LOW) as i128;
            un[j + i] = t as u64;
            borrow = -((t >> 64) as i128);
        }
        let t = un[j + n] as i128 - carry as i128 - borrow;
        un[j + n] = t as u64;
        let mut qj = qhat as u64;
        if t < 0 {
            // qhat was one too large: add the divisor back.
            qj -= 1;
            let mut c: u128 = 0;
            for i in 0..n {
                let s = un[j + i] as u128 + vn[i] as u128 + c;
                un[j + i] = s as u64;
                c = s >> 64;
            }
            un[j + n] = (un[j + n] as u128 + c) as u64;
        }
        q[j] = qj;
    }
    // Denormalize the remainder.
    let mut rem = vec![0u64; n];
    if shift == 0 {
        rem.copy_from_slice(&un[..n]);
    } else {
        for i in 0..n - 1 {
            rem[i] = (un[i] >> shift) | (un[i + 1] << (64 - shift));
        }
        rem[n - 1] = un[n - 1] >> shift;
    }
    (q, rem)
}

impl<const P: usize> Ord for Natural<P> {
    fn cmp(&self, rhs: &Self) -> Ordering {
        for i in (0..P).rev() {
            match self.parts[i].cmp(&rhs.parts[i]) {
                Ordering::Equal => {}
                x => return x,
            }
        }
        Ordering::Equal
    }
}
impl<const P: usize> PartialOrd for Natural<P> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}
impl<const P: usize> FromStr for Natural<P> {
    type Err = ParseNaturalError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_decimal(s)
    }
}
impl<const P: usize> Default for Natural<P> {
    fn default() -> Self {
        Self::ZERO
    }
}
impl<const P: usize> From<u64> for Natural<P> {
    fn from(v: u64) -> Self {
        Self::from_u64(v)
    }
}

impl<const P: usize> fmt::Display for Natural<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        let mut n = self.clone();
        let mut chunks = Vec::new();
        while !n.is_zero() {
            let (q, r) = n.div_rem_u64(10_000_000_000_000_000_000).unwrap();
            chunks.push(r);
            n = q;
        }
        write!(f, "{}", chunks.pop().unwrap())?;
        while let Some(c) = chunks.pop() {
            write!(f, "{c:019}")?
        }
        Ok(())
    }
}
impl<const P: usize> fmt::LowerHex for Natural<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        let i = self.parts.iter().rposition(|&x| x != 0).unwrap();
        write!(f, "{:x}", self.parts[i])?;
        for x in self.parts[..i].iter().rev() {
            write!(f, "{x:016x}")?
        }
        Ok(())
    }
}
impl<const P: usize> fmt::UpperHex for Natural<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        let i = self.parts.iter().rposition(|&x| x != 0).unwrap();
        write!(f, "{:X}", self.parts[i])?;
        for x in self.parts[..i].iter().rev() {
            write!(f, "{x:016X}")?
        }
        Ok(())
    }
}
impl<const P: usize> fmt::Debug for Natural<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Natural<{P}>(0x{:x})", self)
    }
}

macro_rules! binop {
    ($trait:ident,$method:ident,$func:ident) => {
        impl<'a, 'b, const P: usize> $trait<&'b Natural<P>> for &'a Natural<P> {
            type Output = Natural<P>;
            fn $method(self, rhs: &'b Natural<P>) -> Self::Output {
                self.$func(rhs)
            }
        }
        impl<const P: usize> $trait for Natural<P> {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self {
                (&self).$method(&rhs)
            }
        }
    };
}
binop!(Add, add, wrapping_add);
binop!(Sub, sub, wrapping_sub);
binop!(Mul, mul, wrapping_mul);
impl<const P: usize> Add<&Natural<P>> for Natural<P> {
    type Output = Self;
    fn add(self, rhs: &Self) -> Self {
        (&self).add(rhs)
    }
}
impl<const P: usize> Sub<&Natural<P>> for Natural<P> {
    type Output = Self;
    fn sub(self, rhs: &Self) -> Self {
        (&self).sub(rhs)
    }
}
impl<const P: usize> Mul<&Natural<P>> for Natural<P> {
    type Output = Self;
    fn mul(self, rhs: &Self) -> Self {
        (&self).mul(rhs)
    }
}
impl<const P: usize> Add<Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn add(self, rhs: Natural<P>) -> Natural<P> {
        self.add(&rhs)
    }
}
impl<const P: usize> Sub<Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn sub(self, rhs: Natural<P>) -> Natural<P> {
        self.sub(&rhs)
    }
}
impl<const P: usize> Mul<Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn mul(self, rhs: Natural<P>) -> Natural<P> {
        self.mul(&rhs)
    }
}
impl<const P: usize> Div for Natural<P> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.div_rem(&rhs).expect("division by zero").0
    }
}
impl<const P: usize> Rem for Natural<P> {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self {
        self.div_rem(&rhs).expect("division by zero").1
    }
}
impl<const P: usize> Div<&Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn div(self, rhs: &Natural<P>) -> Natural<P> {
        self.div_rem(rhs).expect("division by zero").0
    }
}
impl<const P: usize> Rem<&Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn rem(self, rhs: &Natural<P>) -> Natural<P> {
        self.div_rem(rhs).expect("division by zero").1
    }
}
macro_rules! assign {($trait:ident,$method:ident,$op:tt)=>{impl<const P:usize>$trait<&Natural<P>> for Natural<P>{fn $method(&mut self,rhs:&Self){*self=&*self $op rhs;}}};}
assign!(AddAssign,add_assign,+);
assign!(SubAssign,sub_assign,-);
assign!(MulAssign,mul_assign,*);
impl<const P: usize> AddAssign for Natural<P> {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs
    }
}
impl<const P: usize> SubAssign for Natural<P> {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs
    }
}
impl<const P: usize> MulAssign for Natural<P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs
    }
}
impl<const P: usize> DivAssign<&Self> for Natural<P> {
    fn div_assign(&mut self, r: &Self) {
        *self = self.clone() / r.clone()
    }
}
impl<const P: usize> RemAssign<&Self> for Natural<P> {
    fn rem_assign(&mut self, r: &Self) {
        *self = self.clone() % r.clone()
    }
}
impl<const P: usize> DivAssign for Natural<P> {
    fn div_assign(&mut self, rhs: Self) {
        *self /= &rhs
    }
}
impl<const P: usize> RemAssign for Natural<P> {
    fn rem_assign(&mut self, rhs: Self) {
        *self %= &rhs
    }
}

impl<const P: usize> BitAnd for Natural<P> {
    type Output = Self;
    fn bitand(mut self, rhs: Self) -> Self {
        for i in 0..P {
            self.parts[i] &= rhs.parts[i]
        }
        self
    }
}
impl<const P: usize> BitOr for Natural<P> {
    type Output = Self;
    fn bitor(mut self, rhs: Self) -> Self {
        for i in 0..P {
            self.parts[i] |= rhs.parts[i]
        }
        self
    }
}
impl<const P: usize> BitXor for Natural<P> {
    type Output = Self;
    fn bitxor(mut self, rhs: Self) -> Self {
        for i in 0..P {
            self.parts[i] ^= rhs.parts[i]
        }
        self
    }
}
impl<const P: usize> BitAnd<&Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn bitand(self, r: &Natural<P>) -> Natural<P> {
        self.clone() & r.clone()
    }
}
impl<const P: usize> BitOr<&Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn bitor(self, r: &Natural<P>) -> Natural<P> {
        self.clone() | r.clone()
    }
}
impl<const P: usize> BitXor<&Natural<P>> for &Natural<P> {
    type Output = Natural<P>;
    fn bitxor(self, r: &Natural<P>) -> Natural<P> {
        self.clone() ^ r.clone()
    }
}
impl<const P: usize> BitAndAssign<&Natural<P>> for Natural<P> {
    fn bitand_assign(&mut self, r: &Natural<P>) {
        for i in 0..P {
            self.parts[i] &= r.parts[i]
        }
    }
}
impl<const P: usize> BitOrAssign<&Natural<P>> for Natural<P> {
    fn bitor_assign(&mut self, r: &Natural<P>) {
        for i in 0..P {
            self.parts[i] |= r.parts[i]
        }
    }
}
impl<const P: usize> BitXorAssign<&Natural<P>> for Natural<P> {
    fn bitxor_assign(&mut self, r: &Natural<P>) {
        for i in 0..P {
            self.parts[i] ^= r.parts[i]
        }
    }
}
impl<const P: usize> BitAndAssign for Natural<P> {
    fn bitand_assign(&mut self, r: Self) {
        *self &= &r
    }
}
impl<const P: usize> BitOrAssign for Natural<P> {
    fn bitor_assign(&mut self, r: Self) {
        *self |= &r
    }
}
impl<const P: usize> BitXorAssign for Natural<P> {
    fn bitxor_assign(&mut self, r: Self) {
        *self ^= &r
    }
}
impl<const P: usize> Not for Natural<P> {
    type Output = Self;
    fn not(mut self) -> Self {
        for x in &mut self.parts {
            *x = !*x
        }
        self
    }
}
impl<const P: usize> Shl<usize> for Natural<P> {
    type Output = Self;
    fn shl(self, s: usize) -> Self {
        if s >= Self::BITS {
            return Self::ZERO;
        }
        let mut out = Self::ZERO;
        let w = s / 64;
        let b = s % 64;
        for i in w..P {
            out.parts[i] = self.parts[i - w] << b;
            if b != 0 && i > w {
                out.parts[i] |= self.parts[i - w - 1] >> (64 - b)
            }
        }
        out
    }
}
impl<const P: usize> Shr<usize> for Natural<P> {
    type Output = Self;
    fn shr(self, s: usize) -> Self {
        if s >= Self::BITS {
            return Self::ZERO;
        }
        let mut out = Self::ZERO;
        let w = s / 64;
        let b = s % 64;
        for i in 0..P - w {
            out.parts[i] = self.parts[i + w] >> b;
            if b != 0 && i + w + 1 < P {
                out.parts[i] |= self.parts[i + w + 1] << (64 - b)
            }
        }
        out
    }
}
impl<const P: usize> ShlAssign<usize> for Natural<P> {
    fn shl_assign(&mut self, s: usize) {
        *self = self.clone() << s
    }
}
impl<const P: usize> ShrAssign<usize> for Natural<P> {
    fn shr_assign(&mut self, s: usize) {
        *self = self.clone() >> s
    }
}

/// The exact 2P-limb result of a multiplication.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct WideNatural<const P: usize> {
    low: [u64; P],
    high: [u64; P],
}
impl<const P: usize> WideNatural<P> {
    pub fn low(&self) -> Natural<P> {
        Natural { parts: self.low }
    }
    pub fn high(&self) -> Natural<P> {
        Natural { parts: self.high }
    }
    pub fn overflowing_narrow(self) -> (Natural<P>, bool) {
        (
            Natural { parts: self.low },
            self.high.iter().any(|&x| x != 0),
        )
    }
    /// Reduce this 2P-limb value modulo `m` (`m` must be nonzero).
    pub(crate) fn rem_natural(&self, m: &Natural<P>) -> Natural<P> {
        // Fast single-limb modulus.
        if m.parts[1..].iter().all(|&x| x == 0) {
            let d = m.parts[0] as u128;
            let mut r = 0u128;
            // Most-significant limb first: high limbs (descending) then low limbs.
            for limb in self.low.iter().chain(self.high.iter()).rev() {
                r = ((r << 64) | *limb as u128) % d;
            }
            return Natural::from_u64(r as u64);
        }
        let mut wide = Vec::with_capacity(2 * P);
        wide.extend_from_slice(&self.low);
        wide.extend_from_slice(&self.high);
        let (_, r) = knuth_divmod(&wide, &m.parts);
        let mut rn = Natural::ZERO;
        for (slot, &x) in rn.parts.iter_mut().zip(r.iter()) {
            *slot = x;
        }
        rn
    }
}

/// Extended-GCD result. Coefficients are intentionally private pending a public signed type.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExtendedGcdResult<const P: usize> {
    pub gcd: Natural<P>,
}

/// Invalid Montgomery modulus.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MontgomeryError {
    ZeroModulus,
    EvenModulus,
}
impl fmt::Display for MontgomeryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::ZeroModulus => "zero Montgomery modulus",
            Self::EvenModulus => "Montgomery modulus must be odd",
        })
    }
}
impl std::error::Error for MontgomeryError {}

/// Modular arithmetic context for an odd modulus.
pub struct Montgomery<const P: usize> {
    modulus: Natural<P>,
}
impl<const P: usize> Montgomery<P> {
    pub fn new(modulus: Natural<P>) -> Result<Self, MontgomeryError> {
        if modulus.is_zero() {
            Err(MontgomeryError::ZeroModulus)
        } else if modulus.is_even() {
            Err(MontgomeryError::EvenModulus)
        } else {
            Ok(Self { modulus })
        }
    }
    pub fn encode(&self, v: &Natural<P>) -> Natural<P> {
        v.div_rem(&self.modulus).unwrap().1
    }
    pub fn decode(&self, v: &Natural<P>) -> Natural<P> {
        self.encode(v)
    }
    pub fn mul(&self, a: &Natural<P>, b: &Natural<P>) -> Natural<P> {
        a.mul_mod(b, &self.modulus)
    }
    pub fn square(&self, a: &Natural<P>) -> Natural<P> {
        self.mul(a, a)
    }
    pub fn pow(&self, a: &Natural<P>, e: &Natural<P>) -> Natural<P> {
        a.pow_mod(e, &self.modulus)
    }
    pub fn inv(&self, v: &Natural<P>) -> Option<Natural<P>> {
        if v.gcd(&self.modulus) != Natural::ONE {
            return None;
        }
        // Extended Euclid with coefficients represented modulo the modulus.
        let mut old_r = self.modulus.clone();
        let mut r = v.div_rem(&self.modulus)?.1;
        let mut old_t = Natural::ZERO;
        let mut t = Natural::ONE;
        while !r.is_zero() {
            let (q, next_r) = old_r.div_rem(&r)?;
            old_r = r;
            r = next_r;
            let qt = q.mul_mod(&t, &self.modulus);
            let next_t = if old_t >= qt {
                old_t.wrapping_sub(&qt)
            } else {
                self.modulus.wrapping_sub(&qt.wrapping_sub(&old_t))
            };
            old_t = t;
            t = next_t;
        }
        Some(old_t)
    }
}

pub fn jacobi_u64(mut a: u64, mut n: u64) -> i8 {
    if n == 0 || n & 1 == 0 {
        return 0;
    }
    a %= n;
    let mut s = 1;
    while a != 0 {
        while a & 1 == 0 {
            a >>= 1;
            let r = n & 7;
            if r == 3 || r == 5 {
                s = -s
            }
        }
        core::mem::swap(&mut a, &mut n);
        if a & 3 == 3 && n & 3 == 3 {
            s = -s
        }
        a %= n
    }
    if n == 1 { s } else { 0 }
}
pub fn legendre_u32(n: u32, p: u32) -> i8 {
    if p == 2 {
        return if n & 1 == 0 { 0 } else { 1 };
    }
    jacobi_u64((n % p) as u64, p as u64)
}
pub fn tonelli_shanks_u32(n: u32, p: u32) -> Option<u32> {
    if p == 2 {
        return Some(n & 1);
    }
    let n = n % p;
    if n == 0 {
        return Some(0);
    }
    if legendre_u32(n, p) != 1 {
        return None;
    }
    if p & 3 == 3 {
        return Some(modpow_u64(n as u64, ((p + 1) / 4) as u64, p as u64) as u32);
    }
    let mut q = p - 1;
    let mut s = 0;
    while q & 1 == 0 {
        q >>= 1;
        s += 1
    }
    let mut z = 2;
    while legendre_u32(z, p) != -1 {
        z += 1
    }
    let mut c = modpow_u64(z as u64, q as u64, p as u64);
    let mut x = modpow_u64(n as u64, q.div_ceil(2) as u64, p as u64);
    let mut t = modpow_u64(n as u64, q as u64, p as u64);
    let mut m = s;
    while t != 1 {
        let mut i = 1;
        let mut tt = t * t % p as u64;
        while tt != 1 {
            tt = tt * tt % p as u64;
            i += 1;
            if i >= m {
                return None;
            }
        }
        let b = modpow_u64(c, 1u64 << (m - i - 1), p as u64);
        x = x * b % p as u64;
        t = t * b % p as u64 * b % p as u64;
        c = b * b % p as u64;
        m = i
    }
    Some(x as u32)
}
fn modpow_u64(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1;
    while e != 0 {
        if e & 1 != 0 {
            r = (r as u128 * a as u128 % m as u128) as u64
        }
        a = (a as u128 * a as u128 % m as u128) as u64;
        e >>= 1
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parse_arithmetic() {
        let n: Natural<2> = "340282366920938463463374607431768211455".parse().unwrap();
        assert_eq!(n, Natural::MAX);
        assert_eq!(n.to_string(), "340282366920938463463374607431768211455");
        assert_eq!(
            Natural::<1>::from_decimal("18446744073709551616"),
            Err(ParseNaturalError::Overflow)
        );
    }
    #[test]
    fn div_sqrt() {
        let n = Natural::<2>::from_u64(123456789);
        let d = Natural::from_u64(1234);
        let (q, r) = n.div_rem(&d).unwrap();
        assert_eq!(q.to_string(), "100046");
        assert_eq!(r.to_string(), "25");
        assert_eq!(n.floor_sqrt(), Natural::from_u64(11111));
    }
    #[test]
    fn modular() {
        let m = Natural::<2>::from_u64(97);
        assert_eq!(
            Natural::from_u64(88).mul_mod(&Natural::from_u64(77), &m),
            Natural::from_u64(83)
        );
        assert_eq!(tonelli_shanks_u32(10, 13), Some(7).or(Some(6)));
    }
    #[test]
    fn agrees_with_u128() {
        let mut state = 0x1234_5678_9abc_def0u128;
        for _ in 0..500 {
            state = state.wrapping_mul(0xda942042e4dd58b5).wrapping_add(1);
            let a = state;
            state = state.wrapping_mul(0xda942042e4dd58b5).wrapping_add(1);
            let b = state;
            let na = Natural::<2>::from_le_bytes(&a.to_le_bytes()).unwrap();
            let nb = Natural::<2>::from_le_bytes(&b.to_le_bytes()).unwrap();
            assert_eq!((&na + &nb).to_string(), a.wrapping_add(b).to_string());
            assert_eq!((&na * &nb).to_string(), a.wrapping_mul(b).to_string());
            if let Some(expected_q) = a.checked_div(b) {
                let (q, r) = na.div_rem(&nb).unwrap();
                assert_eq!(q.to_string(), expected_q.to_string());
                assert_eq!(r.to_string(), (a % b).to_string());
            }
        }
    }
    #[test]
    fn inverse_for_composite_modulus() {
        let m = Montgomery::<1>::new(Natural::from_u64(15)).unwrap();
        let inverse = m.inv(&Natural::from_u64(2)).unwrap();
        assert_eq!(inverse, Natural::from_u64(8));
    }
}

#[cfg(test)]
mod difftests {
    //! Randomized differential tests against `num-bigint` (dev-only oracle).
    //! Guards the fast arithmetic kernels (SPEC §19.3).
    use super::*;
    use num_bigint::BigUint;

    const P: usize = 16;

    fn to_big(n: &Natural<P>) -> BigUint {
        let mut bytes = Vec::with_capacity(P * 8);
        for &limb in n.as_parts() {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        BigUint::from_bytes_le(&bytes)
    }
    fn wide_to_big(w: &WideNatural<P>) -> BigUint {
        let mut bytes = Vec::with_capacity(2 * P * 8);
        for &limb in w.low.iter().chain(w.high.iter()) {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        BigUint::from_bytes_le(&bytes)
    }
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        /// Random Natural with `limbs` significant limbs (exercises the
        /// significant-limb-aware fast paths at every width).
        fn natural(&mut self, limbs: usize) -> Natural<P> {
            let mut n = Natural::<P>::ZERO;
            for i in 0..limbs.min(P) {
                n.as_mut_parts()[i] = self.next();
            }
            n
        }
    }

    #[test]
    fn diff_add_sub_mul() {
        let mut rng = Rng(0xdead_beef_1234_5678);
        let modulus = BigUint::from(1u8) << (64 * P);
        for _ in 0..3000 {
            let la = 1 + (rng.next() as usize % P);
            let lb = 1 + (rng.next() as usize % P);
            let a = rng.natural(la);
            let b = rng.natural(lb);
            let ba = to_big(&a);
            let bb = to_big(&b);
            assert_eq!(to_big(&a.wrapping_add(&b)), (&ba + &bb) % &modulus);
            assert_eq!(
                to_big(&a.wrapping_sub(&b)),
                (&ba + &modulus - (&bb % &modulus)) % &modulus
            );
            assert_eq!(to_big(&a.wrapping_mul(&b)), (&ba * &bb) % &modulus);
            assert_eq!(wide_to_big(&a.widening_mul(&b)), &ba * &bb);
        }
    }

    #[test]
    fn diff_div_rem() {
        let mut rng = Rng(0x1122_3344_5566_7788);
        for _ in 0..3000 {
            let la = 1 + (rng.next() as usize % P);
            let lb = 1 + (rng.next() as usize % P);
            let a = rng.natural(la);
            let mut b = rng.natural(lb);
            if b.is_zero() {
                b = Natural::ONE;
            }
            let (q, r) = a.div_rem(&b).unwrap();
            let (ba, bb) = (to_big(&a), to_big(&b));
            assert_eq!(to_big(&q), &ba / &bb, "quot a={ba} b={bb}");
            assert_eq!(to_big(&r), &ba % &bb, "rem a={ba} b={bb}");
        }
    }

    #[test]
    fn diff_div_rem_u64() {
        let mut rng = Rng(0x9988_7766_5544_3322);
        for _ in 0..3000 {
            let la = 1 + (rng.next() as usize % P);
            let a = rng.natural(la);
            let d = rng.next() | 1;
            let (q, r) = a.div_rem_u64(d).unwrap();
            let (ba, bd) = (to_big(&a), BigUint::from(d));
            assert_eq!(to_big(&q), &ba / &bd);
            assert_eq!(BigUint::from(r), &ba % &bd);
        }
    }

    fn big_gcd(mut a: BigUint, mut b: BigUint) -> BigUint {
        while b != BigUint::from(0u8) {
            let r = &a % &b;
            a = b;
            b = r;
        }
        a
    }
    #[test]
    fn diff_gcd() {
        let mut rng = Rng(0x0f0f_0f0f_f0f0_f0f0);
        for _ in 0..2000 {
            let la = 1 + (rng.next() as usize % P);
            let lb = 1 + (rng.next() as usize % P);
            let a = rng.natural(la);
            let b = rng.natural(lb);
            assert_eq!(to_big(&a.gcd(&b)), big_gcd(to_big(&a), to_big(&b)));
        }
    }

    #[test]
    fn diff_mul_mod_pow_mod() {
        let mut rng = Rng(0xabcd_1234_ef56_7890);
        for _ in 0..2000 {
            // Include single-limb moduli (exercise the fast rem path).
            let lm = 1 + (rng.next() as usize % P);
            let mut m = rng.natural(lm);
            if m < Natural::from_u64(2) {
                m = Natural::from_u64(3);
            }
            let a = rng.natural(P);
            let b = rng.natural(P);
            let bm = to_big(&m);
            let am = a.div_rem(&m).unwrap().1;
            let bmod = b.div_rem(&m).unwrap().1;
            assert_eq!(
                to_big(&am.mul_mod(&bmod, &m)),
                (to_big(&am) * to_big(&bmod)) % &bm
            );
            let e = rng.natural(2);
            assert_eq!(
                to_big(&a.pow_mod(&e, &m)),
                to_big(&am).modpow(&to_big(&e), &bm)
            );
        }
    }

    #[test]
    fn diff_montgomery() {
        let mut rng = Rng(0x5a5a_a5a5_1357_9bdf);
        for _ in 0..2000 {
            let lm = 1 + (rng.next() as usize % P);
            let mut m = rng.natural(lm);
            m.as_mut_parts()[0] |= 1; // odd
            if m < Natural::from_u64(3) {
                m = Natural::from_u64(3);
            }
            let mont = Montgomery::<P>::new(m.clone()).unwrap();
            let bm = to_big(&m);
            let a = rng.natural(P).div_rem(&m).unwrap().1;
            let b = rng.natural(P).div_rem(&m).unwrap().1;
            // encode/decode round-trip
            assert_eq!(to_big(&mont.decode(&mont.encode(&a))), to_big(&a));
            // multiplication in Montgomery domain equals modular product
            assert_eq!(
                to_big(&mont.mul(&a, &b)),
                (to_big(&a) * to_big(&b)) % &bm,
                "mont.mul m={bm}"
            );
            let e = rng.natural(2);
            assert_eq!(
                to_big(&mont.pow(&a, &e)),
                to_big(&a).modpow(&to_big(&e), &bm)
            );
        }
    }
}
