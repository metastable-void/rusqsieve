//! WebAssembly SIMD root advancement for the SIMD artifact.

use core::arch::wasm32::*;

#[inline]
pub(super) fn advance_add(root1: &mut [u32], root2: &mut [u32], delta: &[u32], primes: &[u32]) {
    let len = root1
        .len()
        .min(root2.len())
        .min(delta.len())
        .min(primes.len());
    // SAFETY: the SIMD artifact is built with simd128 enabled and the callee
    // operates only on the checked common prefix.
    unsafe { advance_simd(root1, root2, delta, primes, len, true) }
}

#[inline]
pub(super) fn advance_sub(root1: &mut [u32], root2: &mut [u32], delta: &[u32], primes: &[u32]) {
    let len = root1
        .len()
        .min(root2.len())
        .min(delta.len())
        .min(primes.len());
    unsafe { advance_simd(root1, root2, delta, primes, len, false) }
}

#[target_feature(enable = "simd128")]
unsafe fn advance_simd(
    root1: &mut [u32],
    root2: &mut [u32],
    delta: &[u32],
    primes: &[u32],
    len: usize,
    add: bool,
) {
    let invalid_value = u32x4_splat(u32::MAX);
    let mut index = 0usize;
    while index + 4 <= len {
        unsafe {
            let old1 = v128_load(root1.as_ptr().add(index).cast());
            let old2 = v128_load(root2.as_ptr().add(index).cast());
            let d = v128_load(delta.as_ptr().add(index).cast());
            let p = v128_load(primes.as_ptr().add(index).cast());
            let invalid = u32x4_eq(old1, invalid_value);
            let (next1, next2) = if add {
                let sum1 = u32x4_add(old1, d);
                let sum2 = u32x4_add(old2, d);
                (
                    u32x4_sub(sum1, v128_and(u32x4_ge(sum1, p), p)),
                    u32x4_sub(sum2, v128_and(u32x4_ge(sum2, p), p)),
                )
            } else {
                (
                    u32x4_add(u32x4_sub(old1, d), v128_and(u32x4_lt(old1, d), p)),
                    u32x4_add(u32x4_sub(old2, d), v128_and(u32x4_lt(old2, d), p)),
                )
            };
            let ordered = u32x4_lt(next1, next2);
            let minimum = v128_bitselect(next1, next2, ordered);
            let maximum = v128_bitselect(next2, next1, ordered);
            let out1 = v128_bitselect(old1, minimum, invalid);
            let out2 = v128_bitselect(old2, maximum, invalid);
            v128_store(root1.as_mut_ptr().add(index).cast(), out1);
            v128_store(root2.as_mut_ptr().add(index).cast(), out2);
        }
        index += 4;
    }
    for lane in index..len {
        if root1[lane] == u32::MAX {
            continue;
        }
        let p = primes[lane];
        let d = delta[lane];
        let (a, b) = if add {
            (
                super::add_mod_u32(root1[lane], d, p),
                super::add_mod_u32(root2[lane], d, p),
            )
        } else {
            (
                super::sub_mod_u32(root1[lane], d, p),
                super::sub_mod_u32(root2[lane], d, p),
            )
        };
        root1[lane] = a.min(b);
        root2[lane] = a.max(b);
    }
}
