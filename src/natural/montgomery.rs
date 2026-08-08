//! Montgomery arithmetic for the Pollard-Brent stage.
//!
//! Rho spends essentially all of its time in two modular multiplications per iteration — `y² + c`
//! and `q · (x − y)` — so this module exists to make exactly those as cheap as they can be made in
//! portable Rust. Four things dominate, in the order they were measured to matter:
//!
//! 1. **Fixed per-call overhead.** The first implementation multiplied into a `2P`-word product,
//!    copied it into a 33-word scratch array that was zeroed on entry, reduced, and copied out —
//!    hundreds of bytes of memory traffic per multiply regardless of the modulus. At a 128-bit
//!    modulus (two limbs, eight word-multiplies of real work) that overhead was most of the cost.
//!    Everything here works in place over exactly the significant limbs.
//! 2. **Loop shape.** The limb count is fixed for the whole stage, so the inner loops are
//!    monomorphized over it and fully unrolled, with no `sig_len` scans and no per-word capacity
//!    branches.
//! 3. **Multiplication count.** Squaring is symmetric — `a[i]·a[j]` is computed once and doubled —
//!    which removes about a quarter of the word-multiplies from the squaring half of the loop.
//! 4. **Limb width per target.** x86-64 has a widening 64×64→128 multiply; wasm does not, so every
//!    `u128` product there is emulated out of 32-bit pieces. Using 32-bit limbs on wasm makes every
//!    product a single `i64.mul` on extended operands and removes 128-bit arithmetic entirely,
//!    measured at 1.7× on 512- to 1024-bit moduli. Both widths are generated from one macro and
//!    both are tested; [`Limb`] selects which one a target uses.
//!
//! The representation is the standard one: a residue is `a·R mod n` with `R = 2^(bits·limbs)`, and
//! values are kept canonical in `[0, n)`.
//!
//! ## Buffer invariant
//!
//! Every routine here reads and writes only the low `limbs` words of its operands and leaves the
//! rest untouched. Callers hold `[Limb; LIMB_CAP]` buffers whose upper words are zero from
//! construction and stay zero, so no routine has to clear them — which is the point, since the
//! capacity is 1024 bits and a typical modulus is a quarter of that or less.
use super::{Natural, sig_len};

/// Widest modulus this module handles, in 64-bit words. Matches the engine's `Natural` capacity;
/// `MontgomeryContext::new` refuses anything wider.
const MAX_WORDS: usize = 16;

/// Limb type for this target.
///
/// wasm32 has no widening 64×64 multiply, so 64-bit limbs there mean every product is emulated;
/// 32-bit limbs keep every product inside one `i64.mul`. Everything else has a real widening
/// multiply and wants the wider limb, which quarters the number of products.
#[cfg(target_arch = "wasm32")]
pub(crate) type Limb = u32;
#[cfg(not(target_arch = "wasm32"))]
pub(crate) type Limb = u64;

/// Buffer length, in [`Limb`]s, that holds the widest supported modulus.
pub(crate) const LIMB_CAP: usize = MAX_WORDS * (64 / Limb::BITS as usize);

/// Generates the limb arithmetic for one (limb, double-width) pair. Both instantiations are always
/// compiled so that the tests can check them against each other on any host.
macro_rules! limb_arithmetic {
    ($module:ident, $limb:ty, $wide:ty, $newton:literal, [$($width:literal),*]) => {
        pub(super) mod $module {
            /// Words a modulus of the maximum supported width occupies in this limb type.
            pub(super) const CAP: usize = super::MAX_WORDS * (64 / <$limb>::BITS as usize);
            /// CIOS needs two words of headroom above the modulus width.
            const CIOS_WORDS: usize = CAP + 2;
            /// Separated squaring needs a full double-width product plus a carry word.
            const SQUARE_WORDS: usize = 2 * CAP + 1;
            const BITS: u32 = <$limb>::BITS;

            /// `-n^-1 mod 2^BITS`, by Newton iteration: each round doubles the correct low bits,
            /// starting from 1, which is correct modulo 2.
            pub(super) fn negative_inverse(n0: $limb) -> $limb {
                let mut inverse: $limb = 1;
                for _ in 0..$newton {
                    inverse = inverse.wrapping_mul((2 as $limb).wrapping_sub(n0.wrapping_mul(inverse)));
                }
                debug_assert_eq!(n0.wrapping_mul(inverse.wrapping_neg()), <$limb>::MAX);
                inverse.wrapping_neg()
            }

            /// Subtracts `n` from `value` when `value >= n`, given whether the addition that
            /// produced it carried out of the top word.
            #[inline(always)]
            pub(super) fn conditional_subtract(value: &mut [$limb], n: &[$limb], k: usize, carry_out: $limb) {
                let mut difference = [0 as $limb; CAP];
                let mut borrow = 0 as $limb;
                for j in 0..k {
                    let (d, first) = value[j].overflowing_sub(n[j]);
                    let (d, second) = d.overflowing_sub(borrow);
                    difference[j] = d;
                    borrow = <$limb>::from(first) | <$limb>::from(second);
                }
                // Take the difference when the value did not underflow (so it was >= n), or when it
                // overflowed the top word in the first place. Branchless: the outcome is
                // data-dependent and would mispredict about half the time.
                let take = (carry_out | (1 ^ borrow)) & 1;
                let mask = take.wrapping_neg();
                for j in 0..k {
                    value[j] = (difference[j] & mask) | (value[j] & !mask);
                }
            }

            /// `out = (lhs + rhs) mod n`, for operands already in `[0, n)`.
            #[inline(always)]
            pub(super) fn add(lhs: &[$limb], rhs: &[$limb], n: &[$limb], k: usize, out: &mut [$limb]) {
                let mut carry = 0 as $limb;
                for j in 0..k {
                    let (sum, first) = lhs[j].overflowing_add(rhs[j]);
                    let (sum, second) = sum.overflowing_add(carry);
                    out[j] = sum;
                    carry = <$limb>::from(first) | <$limb>::from(second);
                }
                conditional_subtract(out, n, k, carry);
            }

            /// `out = (lhs − rhs) mod n`, for operands already in `[0, n)`.
            #[inline(always)]
            pub(super) fn sub(lhs: &[$limb], rhs: &[$limb], n: &[$limb], k: usize, out: &mut [$limb]) {
                let mut borrow = 0 as $limb;
                for j in 0..k {
                    let (difference, first) = lhs[j].overflowing_sub(rhs[j]);
                    let (difference, second) = difference.overflowing_sub(borrow);
                    out[j] = difference;
                    borrow = <$limb>::from(first) | <$limb>::from(second);
                }
                // Underflow means the true difference is negative; adding n brings it into range.
                let mask = borrow.wrapping_neg();
                let mut carry = 0 as $limb;
                for j in 0..k {
                    let (sum, first) = out[j].overflowing_add(n[j] & mask);
                    let (sum, second) = sum.overflowing_add(carry);
                    out[j] = sum;
                    carry = <$limb>::from(first) | <$limb>::from(second);
                }
            }

            /// Coarsely Integrated Operand Scanning: one pass interleaving multiplication with
            /// reduction, so no double-width product is ever materialized.
            #[inline(always)]
            pub(super) fn cios<const K: usize>(a: &[$limb], b: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                let mut t = [0 as $limb; CIOS_WORDS];
                for i in 0..K {
                    let bi = b[i] as $wide;
                    let mut carry = 0 as $wide;
                    for j in 0..K {
                        let value = t[j] as $wide + a[j] as $wide * bi + carry;
                        t[j] = value as $limb;
                        carry = value >> BITS;
                    }
                    let value = t[K] as $wide + carry;
                    t[K] = value as $limb;
                    t[K + 1] = (value >> BITS) as $limb;

                    // m makes the low word vanish, so the shift below is exact.
                    let m = t[0].wrapping_mul(n0inv) as $wide;
                    let value = t[0] as $wide + m * n[0] as $wide;
                    let mut carry = value >> BITS;
                    for j in 1..K {
                        let value = t[j] as $wide + m * n[j] as $wide + carry;
                        t[j - 1] = value as $limb;
                        carry = value >> BITS;
                    }
                    let value = t[K] as $wide + carry;
                    t[K - 1] = value as $limb;
                    t[K] = t[K + 1].wrapping_add((value >> BITS) as $limb);
                }
                out[..K].copy_from_slice(&t[..K]);
                conditional_subtract(out, n, K, t[K]);
            }

            /// The same reduction with the limb count only known at run time.
            pub(super) fn cios_dynamic(k: usize, a: &[$limb], b: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                let mut t = [0 as $limb; CIOS_WORDS];
                for i in 0..k {
                    let bi = b[i] as $wide;
                    let mut carry = 0 as $wide;
                    for j in 0..k {
                        let value = t[j] as $wide + a[j] as $wide * bi + carry;
                        t[j] = value as $limb;
                        carry = value >> BITS;
                    }
                    let value = t[k] as $wide + carry;
                    t[k] = value as $limb;
                    t[k + 1] = (value >> BITS) as $limb;

                    let m = t[0].wrapping_mul(n0inv) as $wide;
                    let value = t[0] as $wide + m * n[0] as $wide;
                    let mut carry = value >> BITS;
                    for j in 1..k {
                        let value = t[j] as $wide + m * n[j] as $wide + carry;
                        t[j - 1] = value as $limb;
                        carry = value >> BITS;
                    }
                    let value = t[k] as $wide + carry;
                    t[k - 1] = value as $limb;
                    t[k] = t[k + 1].wrapping_add((value >> BITS) as $limb);
                }
                out[..k].copy_from_slice(&t[..k]);
                conditional_subtract(out, n, k, t[k]);
            }

            /// Montgomery squaring: symmetric product, then reduction.
            ///
            /// `a[i]·a[j]` and `a[j]·a[i]` are the same word-multiply, so the off-diagonal half is
            /// computed once and doubled and only the `K` diagonal squares are added separately.
            /// That is `(K² + K)/2` word-multiplies against `K²` for a general product.
            #[inline(always)]
            pub(super) fn square<const K: usize>(a: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                let mut t = [0 as $limb; SQUARE_WORDS];
                for i in 0..K {
                    let ai = a[i] as $wide;
                    let mut carry = 0 as $wide;
                    for j in (i + 1)..K {
                        let value = t[i + j] as $wide + ai * a[j] as $wide + carry;
                        t[i + j] = value as $limb;
                        carry = value >> BITS;
                    }
                    t[i + K] = carry as $limb;
                }
                double_and_add_diagonal::<K>(a, &mut t);
                reduce::<K>(&mut t, n, n0inv, out);
            }

            /// Run-time-width squaring; same algorithm as [`square`].
            pub(super) fn square_dynamic(k: usize, a: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                let mut t = [0 as $limb; SQUARE_WORDS];
                for i in 0..k {
                    let ai = a[i] as $wide;
                    let mut carry = 0 as $wide;
                    for j in (i + 1)..k {
                        let value = t[i + j] as $wide + ai * a[j] as $wide + carry;
                        t[i + j] = value as $limb;
                        carry = value >> BITS;
                    }
                    t[i + k] = carry as $limb;
                }
                let mut carry = 0 as $limb;
                for slot in t[..2 * k].iter_mut() {
                    let high = *slot >> (BITS - 1);
                    *slot = (*slot << 1) | carry;
                    carry = high;
                }
                let mut carry = 0 as $wide;
                for i in 0..k {
                    let ai = a[i] as $wide;
                    let value = t[2 * i] as $wide + ai * ai + carry;
                    t[2 * i] = value as $limb;
                    let value = t[2 * i + 1] as $wide + (value >> BITS);
                    t[2 * i + 1] = value as $limb;
                    carry = value >> BITS;
                }
                t[2 * k] = carry as $limb;
                for i in 0..k {
                    let m = t[i].wrapping_mul(n0inv) as $wide;
                    let mut carry = 0 as $wide;
                    for j in 0..k {
                        let value = t[i + j] as $wide + m * n[j] as $wide + carry;
                        t[i + j] = value as $limb;
                        carry = value >> BITS;
                    }
                    let mut at = i + k;
                    let mut overflow = carry as $limb;
                    while overflow != 0 && at <= 2 * k {
                        let (value, next) = t[at].overflowing_add(overflow);
                        t[at] = value;
                        overflow = <$limb>::from(next);
                        at += 1;
                    }
                }
                out[..k].copy_from_slice(&t[k..2 * k]);
                conditional_subtract(out, n, k, t[2 * k]);
            }

            /// Doubles the off-diagonal half of a symmetric product and adds the diagonal squares.
            #[inline(always)]
            fn double_and_add_diagonal<const K: usize>(a: &[$limb], t: &mut [$limb]) {
                let mut carry = 0 as $limb;
                for slot in t[..2 * K].iter_mut() {
                    let high = *slot >> (BITS - 1);
                    *slot = (*slot << 1) | carry;
                    carry = high;
                }
                let mut carry = 0 as $wide;
                for i in 0..K {
                    let ai = a[i] as $wide;
                    let value = t[2 * i] as $wide + ai * ai + carry;
                    t[2 * i] = value as $limb;
                    let value = t[2 * i + 1] as $wide + (value >> BITS);
                    t[2 * i + 1] = value as $limb;
                    carry = value >> BITS;
                }
                t[2 * K] = carry as $limb;
            }

            /// Dispatches a multiply on the run-time limb count, specializing the widths this limb
            /// type actually sees; anything else falls to the run-time-width loop, which costs
            /// roughly a fifth.
            ///
            /// A 64-bit-limb modulus occupies 1..=16 limbs. A 32-bit-limb one occupies twice that
            /// and always an even number, and its table starts at ten rather than two because the
            /// only caller of the narrow backend is the browser's deep rho, which runs on
            /// composites from 257 bits up — specializing the four narrower widths added 7 KiB to
            /// the wasm artifact, a browser shipping gate, for arithmetic that path never performs.
            #[inline]
            pub(super) fn dispatch_mul(k: usize, a: &[$limb], b: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                dispatch_mul_inline(k, a, b, n, n0inv, out);
            }

            /// The same table, forced inline.
            ///
            /// A `#[target_feature]` wrapper only recompiles what is inlined into it, and the table
            /// above is far too large for the inliner to take on a hint. Without this the x86-64
            /// BMI2 wrapper compiled to a plain call into the baseline copy and emitted no `mulx`
            /// at all — the feature was enabled on a function containing nothing but a call.
            #[inline(always)]
            pub(super) fn dispatch_mul_inline(k: usize, a: &[$limb], b: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                match k {
                    $($width => cios::<$width>(a, b, n, n0inv, out),)*
                    _ => cios_dynamic(k, a, b, n, n0inv, out),
                }
            }

            /// The squaring counterpart of [`dispatch_mul`].
            #[inline]
            pub(super) fn dispatch_square(k: usize, a: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                dispatch_square_inline(k, a, n, n0inv, out);
            }

            /// The squaring counterpart of [`dispatch_mul_inline`].
            #[inline(always)]
            pub(super) fn dispatch_square_inline(k: usize, a: &[$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                match k {
                    $($width => square::<$width>(a, n, n0inv, out),)*
                    _ => square_dynamic(k, a, n, n0inv, out),
                }
            }

            /// Montgomery reduction of a double-width product held in `t`.
            #[inline(always)]
            fn reduce<const K: usize>(t: &mut [$limb], n: &[$limb], n0inv: $limb, out: &mut [$limb]) {
                for i in 0..K {
                    let m = t[i].wrapping_mul(n0inv) as $wide;
                    let mut carry = 0 as $wide;
                    for j in 0..K {
                        let value = t[i + j] as $wide + m * n[j] as $wide + carry;
                        t[i + j] = value as $limb;
                        carry = value >> BITS;
                    }
                    debug_assert_eq!(t[i], 0);
                    let mut at = i + K;
                    let mut overflow = carry as $limb;
                    while overflow != 0 && at <= 2 * K {
                        let (value, next) = t[at].overflowing_add(overflow);
                        t[at] = value;
                        overflow = <$limb>::from(next);
                        at += 1;
                    }
                }
                out[..K].copy_from_slice(&t[K..2 * K]);
                conditional_subtract(out, n, K, t[2 * K]);
            }
        }
    };
}

// u64 limbs for targets with a widening 64×64 multiply; u32 limbs for those without. Five Newton
// rounds reach 32 correct bits, six reach 64.
limb_arithmetic!(
    wide,
    u64,
    u128,
    6,
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
);
limb_arithmetic!(
    narrow,
    u32,
    u64,
    5,
    [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32]
);

/// The backend this target uses. Both are generated; the host test suite checks them against each
/// other so the wasm path is not the only thing that exercises the narrow one.
///
/// Measured through `qs_rho` on the shipped artifact under Node 24.15, same build settings and
/// moduli, only the limb type differing — 32-bit limbs against 64-bit: 320 bits 2.04 vs 2.10 M/s,
/// 400 bits 1.40 vs 1.08, 512 bits 0.99 vs 0.83, 768 bits 0.51 vs 0.37, 1024 bits 0.30 vs 0.22. The
/// browser's deep rho exists for composites the sieve refuses, which start at 257 bits and are
/// mostly wider than 400, so the crossover below 400 costs nothing that path cares about.
#[cfg(target_arch = "wasm32")]
use narrow as backend;
#[cfg(not(target_arch = "wasm32"))]
use wide as backend;

/// Montgomery context for one odd modulus.
///
/// Construction is the only place division-based arithmetic appears; everything after it is limb
/// multiplication and carry propagation.
#[cfg(any(unix, windows, target_arch = "wasm32", test))]
pub(crate) struct MontgomeryContext<const P: usize> {
    modulus: Natural<P>,
    /// The modulus again, in this target's limb type, which is what the arithmetic reads.
    modulus_limbs: [Limb; LIMB_CAP],
    negative_inverse: Limb,
    /// Significant length of the modulus in [`Limb`]s, not in 64-bit words.
    limbs: usize,
    r2: Natural<P>,
    one: Natural<P>,
}

#[cfg(any(unix, windows, target_arch = "wasm32", test))]
impl<const P: usize> MontgomeryContext<P> {
    /// Constructs a context for an odd modulus supported by the inline workspace. Shipped
    /// factorization uses `P == 16`; wider user-selected integer capacities retain the
    /// division-based arithmetic.
    pub(crate) fn new(modulus: &Natural<P>) -> Option<Self> {
        let words = sig_len(modulus.as_parts());
        if words == 0 || modulus.is_even() || P > MAX_WORDS {
            return None;
        }
        let per_word = 64 / Limb::BITS as usize;
        let limbs = words * per_word;
        let mut modulus_limbs = [0 as Limb; LIMB_CAP];
        split_words(modulus.as_parts(), &mut modulus_limbs);
        let negative_inverse = backend::negative_inverse(modulus_limbs[0]);

        // R² mod n is a one-time context cost. Repeated doubling avoids needing a 2P+1-bit
        // temporary merely to express 2^(2·bits·limbs).
        let mut r2 = Natural::ONE;
        if r2 >= *modulus {
            r2 = r2.div_rem(modulus)?.1;
        }
        for _ in 0..(2 * Limb::BITS as usize * limbs) {
            r2 = r2.add_mod(&r2, modulus);
        }

        let mut context = Self {
            modulus: modulus.clone(),
            modulus_limbs,
            negative_inverse,
            limbs,
            r2,
            one: Natural::ZERO,
        };
        context.one = context.encode(&Natural::ONE);
        Some(context)
    }

    pub(crate) fn encode(&self, value: &Natural<P>) -> Natural<P> {
        debug_assert!(value < &self.modulus);
        self.multiply(value, &self.r2)
    }

    #[cfg(test)]
    pub(crate) fn decode(&self, value: &Natural<P>) -> Natural<P> {
        self.multiply(value, &Natural::ONE)
    }

    pub(crate) fn one(&self) -> Natural<P> {
        self.one.clone()
    }

    pub(crate) fn multiply(&self, lhs: &Natural<P>, rhs: &Natural<P>) -> Natural<P> {
        let mut a = [0 as Limb; LIMB_CAP];
        let mut b = [0 as Limb; LIMB_CAP];
        let mut out = [0 as Limb; LIMB_CAP];
        split_words(lhs.as_parts(), &mut a);
        split_words(rhs.as_parts(), &mut b);
        self.mul_raw(&a, &b, &mut out);
        self.store(&out)
    }

    pub(crate) fn square(&self, value: &Natural<P>) -> Natural<P> {
        let mut a = [0 as Limb; LIMB_CAP];
        let mut out = [0 as Limb; LIMB_CAP];
        split_words(value.as_parts(), &mut a);
        self.sqr_raw(&a, &mut out);
        self.store(&out)
    }

    pub(crate) fn add(&self, lhs: &Natural<P>, rhs: &Natural<P>) -> Natural<P> {
        let mut a = [0 as Limb; LIMB_CAP];
        let mut b = [0 as Limb; LIMB_CAP];
        let mut out = [0 as Limb; LIMB_CAP];
        split_words(lhs.as_parts(), &mut a);
        split_words(rhs.as_parts(), &mut b);
        self.add_raw(&a, &b, &mut out);
        self.store(&out)
    }

    // ---------------------------------------------------------------------------------------
    // In-place primitives. These are what the rho loop calls; see the buffer invariant above.
    // ---------------------------------------------------------------------------------------

    /// `out = lhs · rhs · R^-1 mod n`.
    ///
    /// Dispatching on the limb count once per call buys a fully unrolled inner loop with no
    /// loop-carried bounds checks; see `dispatch_mul` for which widths are specialized.
    pub(crate) fn mul_raw(&self, lhs: &[Limb], rhs: &[Limb], out: &mut [Limb]) {
        backend::dispatch_mul(
            self.limbs,
            lhs,
            rhs,
            &self.modulus_limbs,
            self.negative_inverse,
            out,
        );
    }

    /// `out = value² · R^-1 mod n`.
    pub(crate) fn sqr_raw(&self, value: &[Limb], out: &mut [Limb]) {
        backend::dispatch_square(
            self.limbs,
            value,
            &self.modulus_limbs,
            self.negative_inverse,
            out,
        );
    }

    /// `out = (lhs + rhs) mod n`, for operands already in `[0, n)`.
    pub(crate) fn add_raw(&self, lhs: &[Limb], rhs: &[Limb], out: &mut [Limb]) {
        backend::add(lhs, rhs, &self.modulus_limbs, self.limbs, out);
    }

    /// `out = (lhs − rhs) mod n`, for operands already in `[0, n)`.
    ///
    /// The rho loop wants `|x − y|`, and this returns `x − y + n` where that would be negative.
    /// Those differ by a multiple of `n`, so every factor of `n` dividing one divides the other and
    /// the batched gcd is unaffected — which is worth a branch per iteration.
    pub(crate) fn sub_raw(&self, lhs: &[Limb], rhs: &[Limb], out: &mut [Limb]) {
        backend::sub(lhs, rhs, &self.modulus_limbs, self.limbs, out);
    }

    /// `value = value² + addend (mod n)`, the rho iteration itself.
    pub(crate) fn sqr_add_assign(&self, value: &mut [Limb], addend: &[Limb]) {
        let mut product = [0 as Limb; LIMB_CAP];
        self.sqr_raw(value, &mut product);
        self.add_raw(&product, addend, value);
    }

    /// `lhs = lhs · rhs · R^-1 mod n`.
    pub(crate) fn mul_assign(&self, lhs: &mut [Limb], rhs: &[Limb]) {
        let mut product = [0 as Limb; LIMB_CAP];
        self.mul_raw(lhs, rhs, &mut product);
        lhs[..self.limbs].copy_from_slice(&product[..self.limbs]);
    }

    /// Copies a value into a raw limb buffer, splitting words if this target uses narrow limbs.
    pub(crate) fn load(&self, value: &Natural<P>, out: &mut [Limb]) {
        split_words(value.as_parts(), out);
    }

    /// Rebuilds a `Natural` from a raw limb buffer, for the occasional gcd.
    pub(crate) fn store(&self, value: &[Limb]) -> Natural<P> {
        let mut out = Natural::ZERO;
        join_words(value, out.as_mut_parts());
        out
    }

    /// Whether the low `limbs` words are all zero.
    pub(crate) fn is_zero_raw(&self, value: &[Limb]) -> bool {
        value[..self.limbs].iter().all(|&word| word == 0)
    }
}

/// Splits 64-bit words into this target's limbs. A no-op copy when the limb is already 64 bits.
#[inline]
fn split_words(words: &[u64], out: &mut [Limb]) {
    let per_word = 64 / Limb::BITS as usize;
    for (index, &word) in words.iter().enumerate() {
        for part in 0..per_word {
            let slot = index * per_word + part;
            if slot < out.len() {
                out[slot] = (word >> (part as u32 * Limb::BITS)) as Limb;
            }
        }
    }
}

/// The inverse of [`split_words`].
// The widening cast below is a no-op where `Limb` is already `u64` and load-bearing where it is
// `u32`, so exactly one of the two targets always sees it as redundant.
#[allow(clippy::unnecessary_cast)]
#[inline]
fn join_words(limbs: &[Limb], out: &mut [u64]) {
    let per_word = 64 / Limb::BITS as usize;
    for (index, slot) in out.iter_mut().enumerate() {
        let mut word = 0u64;
        for part in 0..per_word {
            let at = index * per_word + part;
            if at < limbs.len() {
                word |= (limbs[at] as u64) << (part as u32 * Limb::BITS);
            }
        }
        *slot = word;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two limb widths must agree exactly. Only one of them ships per target, so without this
    /// the wasm path would be checked by nothing the host test suite runs.
    #[test]
    fn narrow_and_wide_limbs_agree() {
        let mut seed = 0x5deece66du64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            seed
        };
        for words in 1..=8usize {
            for _ in 0..64 {
                let mut n64 = [0u64; 16];
                for slot in n64[..words].iter_mut() {
                    *slot = next();
                }
                n64[0] |= 1;
                n64[words - 1] |= 1 << 63;
                let mut a64 = [0u64; 16];
                let mut b64 = [0u64; 16];
                for j in 0..words {
                    a64[j] = next() % n64[words - 1].max(2);
                    b64[j] = next() % n64[words - 1].max(2);
                }

                let n0inv_wide = wide::negative_inverse(n64[0]);
                let mut wide_out = [0u64; 16];
                wide::cios_dynamic(words, &a64, &b64, &n64, n0inv_wide, &mut wide_out);
                let mut wide_square = [0u64; 16];
                wide::square_dynamic(words, &a64, &n64, n0inv_wide, &mut wide_square);

                let mut n32 = [0u32; 32];
                let mut a32 = [0u32; 32];
                let mut b32 = [0u32; 32];
                for j in 0..words {
                    n32[2 * j] = n64[j] as u32;
                    n32[2 * j + 1] = (n64[j] >> 32) as u32;
                    a32[2 * j] = a64[j] as u32;
                    a32[2 * j + 1] = (a64[j] >> 32) as u32;
                    b32[2 * j] = b64[j] as u32;
                    b32[2 * j + 1] = (b64[j] >> 32) as u32;
                }
                let n0inv_narrow = narrow::negative_inverse(n32[0]);
                let mut narrow_out = [0u32; 32];
                narrow::cios_dynamic(words * 2, &a32, &b32, &n32, n0inv_narrow, &mut narrow_out);
                let mut narrow_square = [0u32; 32];
                narrow::square_dynamic(words * 2, &a32, &n32, n0inv_narrow, &mut narrow_square);

                // The two are the same value only after accounting for R: the narrow backend
                // reduces by 2^(32·2·words), which is the same R as the wide one at this width, so
                // the results must match word for word.
                for j in 0..words {
                    let joined =
                        u64::from(narrow_out[2 * j]) | (u64::from(narrow_out[2 * j + 1]) << 32);
                    assert_eq!(
                        joined, wide_out[j],
                        "multiply differs at {words} words, limb {j}"
                    );
                    let joined = u64::from(narrow_square[2 * j])
                        | (u64::from(narrow_square[2 * j + 1]) << 32);
                    assert_eq!(
                        joined, wide_square[j],
                        "square differs at {words} words, limb {j}"
                    );
                }
            }
        }
    }

    /// The rho inner loop in isolation, one Montgomery squaring and one Montgomery multiply per
    /// iteration, at each limb count that matters.
    ///
    /// This exists to be compared against the same loop written over GMP's mpn assembly, which is
    /// what both YAFU (through mpz) and FLINT (through mpn) run on. Keep the loop body identical to
    /// the reference when changing either.
    #[test]
    #[ignore = "manual Montgomery inner-loop measurement"]
    fn profile_montgomery_loop() {
        const ITERATIONS: u64 = 1_000_000;
        for words in [2usize, 3, 4, 5, 7, 8, 12, 16] {
            let mut seed = 0x9e37_79b9_7f4a_7c15u64;
            let mut modulus = Natural::<16>::ZERO;
            for slot in modulus.as_mut_parts()[..words].iter_mut() {
                seed = seed
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *slot = seed;
            }
            modulus.as_mut_parts()[0] |= 1;
            modulus.as_mut_parts()[words - 1] |= 1 << 63;

            let context = MontgomeryContext::new(&modulus).expect("odd modulus of supported width");
            let mut y = [0 as Limb; LIMB_CAP];
            let mut x = [0 as Limb; LIMB_CAP];
            let mut q = [0 as Limb; LIMB_CAP];
            let mut c = [0 as Limb; LIMB_CAP];
            let mut d = [0 as Limb; LIMB_CAP];
            context.load(&context.encode(&Natural::from_u64(2)), &mut y);
            context.load(&context.encode(&Natural::from_u64(3)), &mut x);
            context.load(&context.one(), &mut q);
            context.load(&context.encode(&Natural::from_u64(1)), &mut c);

            let started = std::time::Instant::now();
            for _ in 0..ITERATIONS {
                context.sqr_add_assign(&mut y, &c);
                context.sub_raw(&x, &y, &mut d);
                context.mul_assign(&mut q, &d);
            }
            let elapsed = started.elapsed().as_secs_f64();
            std::hint::black_box(&q);
            eprintln!(
                "BENCH mont_loop words={words:2} bits={:4} iterations={ITERATIONS} \
                 elapsed={elapsed:.3}s rate={:.0}/s",
                words * 64,
                ITERATIONS as f64 / elapsed
            );
        }
    }
}
