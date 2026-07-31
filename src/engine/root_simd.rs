//! x86-64-baseline SIMD root advancement.
//!
//! The public-to-this-module wrappers are safe: slice lengths are checked and
//! x86-64 guarantees SSE2. Raw loads and stores stay confined here.

use core::arch::x86_64::*;

#[inline]
pub(super) fn advance_add(root1: &mut [u32], root2: &mut [u32], delta: &[u32], primes: &[u32]) {
    let len = root1
        .len()
        .min(root2.len())
        .min(delta.len())
        .min(primes.len());
    // SAFETY: x86-64 has SSE2 in its architecture baseline, and the callee
    // handles only the common in-bounds prefix established above.
    unsafe { advance_sse2(root1, root2, delta, primes, len, true) }
}

#[inline]
pub(super) fn advance_sub(root1: &mut [u32], root2: &mut [u32], delta: &[u32], primes: &[u32]) {
    let len = root1
        .len()
        .min(root2.len())
        .min(delta.len())
        .min(primes.len());
    // SAFETY: see `advance_add`.
    unsafe { advance_sse2(root1, root2, delta, primes, len, false) }
}

#[target_feature(enable = "sse2")]
unsafe fn advance_sse2(
    root1: &mut [u32],
    root2: &mut [u32],
    delta: &[u32],
    primes: &[u32],
    len: usize,
    add: bool,
) {
    let all_ones = _mm_set1_epi32(-1);
    let mut index = 0usize;
    while index + 4 <= len {
        // SAFETY: the loop condition proves all four-lane loads/stores are in
        // the checked common prefix. SSE2 unaligned operations accept the
        // allocation alignment supplied by Vec/Arc.
        unsafe {
            let old1 = _mm_loadu_si128(root1.as_ptr().add(index).cast());
            let old2 = _mm_loadu_si128(root2.as_ptr().add(index).cast());
            let d = _mm_loadu_si128(delta.as_ptr().add(index).cast());
            let p = _mm_loadu_si128(primes.as_ptr().add(index).cast());
            let invalid = _mm_cmpeq_epi32(old1, all_ones);

            let (next1, next2) = if add {
                let sum1 = _mm_add_epi32(old1, d);
                let sum2 = _mm_add_epi32(old2, d);
                let p_minus_one = _mm_sub_epi32(p, _mm_set1_epi32(1));
                let reduce1 = _mm_and_si128(_mm_cmpgt_epi32(sum1, p_minus_one), p);
                let reduce2 = _mm_and_si128(_mm_cmpgt_epi32(sum2, p_minus_one), p);
                (_mm_sub_epi32(sum1, reduce1), _mm_sub_epi32(sum2, reduce2))
            } else {
                let diff1 = _mm_sub_epi32(old1, d);
                let diff2 = _mm_sub_epi32(old2, d);
                let restore1 = _mm_and_si128(_mm_cmpgt_epi32(d, old1), p);
                let restore2 = _mm_and_si128(_mm_cmpgt_epi32(d, old2), p);
                (
                    _mm_add_epi32(diff1, restore1),
                    _mm_add_epi32(diff2, restore2),
                )
            };

            // SSE2 has no unsigned min/max, but every valid prime/root is far
            // below 2^31, so signed compare is identical here.
            let swap = _mm_cmpgt_epi32(next1, next2);
            let minimum = _mm_or_si128(_mm_and_si128(swap, next2), _mm_andnot_si128(swap, next1));
            let maximum = _mm_or_si128(_mm_and_si128(swap, next1), _mm_andnot_si128(swap, next2));
            let out1 = _mm_or_si128(
                _mm_and_si128(invalid, old1),
                _mm_andnot_si128(invalid, minimum),
            );
            let out2 = _mm_or_si128(
                _mm_and_si128(invalid, old2),
                _mm_andnot_si128(invalid, maximum),
            );
            _mm_storeu_si128(root1.as_mut_ptr().add(index).cast(), out1);
            _mm_storeu_si128(root2.as_mut_ptr().add(index).cast(), out2);
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
