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
    // SAFETY: feature detection guards AVX2; x86-64 has SSE2 in its
    // architecture baseline. Both callees handle only the checked prefix.
    unsafe {
        if is_x86_feature_detected!("avx2") {
            advance_avx2(root1, root2, delta, primes, len, true)
        } else {
            advance_sse2(root1, root2, delta, primes, len, true)
        }
    }
}

#[inline]
pub(super) fn advance_sub(root1: &mut [u32], root2: &mut [u32], delta: &[u32], primes: &[u32]) {
    let len = root1
        .len()
        .min(root2.len())
        .min(delta.len())
        .min(primes.len());
    // SAFETY: see `advance_add`.
    unsafe {
        if is_x86_feature_detected!("avx2") {
            advance_avx2(root1, root2, delta, primes, len, false)
        } else {
            advance_sse2(root1, root2, delta, primes, len, false)
        }
    }
}

#[target_feature(enable = "avx2")]
unsafe fn advance_avx2(
    root1: &mut [u32],
    root2: &mut [u32],
    delta: &[u32],
    primes: &[u32],
    len: usize,
    add: bool,
) {
    let all_ones = _mm256_set1_epi32(-1);
    let one = _mm256_set1_epi32(1);
    let mut index = 0usize;
    while index + 8 <= len {
        // SAFETY: the loop condition proves all eight-lane unaligned accesses
        // are inside the checked common prefix.
        unsafe {
            let old1 = _mm256_loadu_si256(root1.as_ptr().add(index).cast());
            let old2 = _mm256_loadu_si256(root2.as_ptr().add(index).cast());
            let d = _mm256_loadu_si256(delta.as_ptr().add(index).cast());
            let p = _mm256_loadu_si256(primes.as_ptr().add(index).cast());
            let invalid = _mm256_cmpeq_epi32(old1, all_ones);
            let (next1, next2) = if add {
                let sum1 = _mm256_add_epi32(old1, d);
                let sum2 = _mm256_add_epi32(old2, d);
                let p_minus_one = _mm256_sub_epi32(p, one);
                let reduce1 = _mm256_and_si256(_mm256_cmpgt_epi32(sum1, p_minus_one), p);
                let reduce2 = _mm256_and_si256(_mm256_cmpgt_epi32(sum2, p_minus_one), p);
                (
                    _mm256_sub_epi32(sum1, reduce1),
                    _mm256_sub_epi32(sum2, reduce2),
                )
            } else {
                let diff1 = _mm256_sub_epi32(old1, d);
                let diff2 = _mm256_sub_epi32(old2, d);
                let restore1 = _mm256_and_si256(_mm256_cmpgt_epi32(d, old1), p);
                let restore2 = _mm256_and_si256(_mm256_cmpgt_epi32(d, old2), p);
                (
                    _mm256_add_epi32(diff1, restore1),
                    _mm256_add_epi32(diff2, restore2),
                )
            };
            let swap = _mm256_cmpgt_epi32(next1, next2);
            let minimum = _mm256_or_si256(
                _mm256_and_si256(swap, next2),
                _mm256_andnot_si256(swap, next1),
            );
            let maximum = _mm256_or_si256(
                _mm256_and_si256(swap, next1),
                _mm256_andnot_si256(swap, next2),
            );
            let out1 = _mm256_or_si256(
                _mm256_and_si256(invalid, old1),
                _mm256_andnot_si256(invalid, minimum),
            );
            let out2 = _mm256_or_si256(
                _mm256_and_si256(invalid, old2),
                _mm256_andnot_si256(invalid, maximum),
            );
            _mm256_storeu_si256(root1.as_mut_ptr().add(index).cast(), out1);
            _mm256_storeu_si256(root2.as_mut_ptr().add(index).cast(), out2);
        }
        index += 8;
    }

    // Reuse the baseline kernel for the at-most-seven-lane tail. Passing
    // subslices keeps its bounds proof intact.
    unsafe {
        advance_sse2(
            &mut root1[index..],
            &mut root2[index..],
            &delta[index..],
            &primes[index..],
            len - index,
            add,
        )
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispatched_root_advancement_matches_scalar_modular_arithmetic() {
        let primes = [
            3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
        ];
        let delta: Vec<u32> = primes
            .iter()
            .enumerate()
            .map(|(i, &p)| (i as u32 * 7 + 2) % p)
            .collect();
        let original1: Vec<u32> = primes
            .iter()
            .enumerate()
            .map(|(i, &p)| if i == 8 { u32::MAX } else { (i as u32 * 5) % p })
            .collect();
        let original2: Vec<u32> = primes
            .iter()
            .enumerate()
            .map(|(i, &p)| if i == 8 { 0 } else { (i as u32 * 11 + 1) % p })
            .collect();

        for add in [true, false] {
            let mut actual1 = original1.clone();
            let mut actual2 = original2.clone();
            if add {
                advance_add(&mut actual1, &mut actual2, &delta, &primes);
            } else {
                advance_sub(&mut actual1, &mut actual2, &delta, &primes);
            }
            for i in 0..primes.len() {
                if original1[i] == u32::MAX {
                    assert_eq!((actual1[i], actual2[i]), (original1[i], original2[i]));
                    continue;
                }
                let first = if add {
                    super::super::add_mod_u32(original1[i], delta[i], primes[i])
                } else {
                    super::super::sub_mod_u32(original1[i], delta[i], primes[i])
                };
                let second = if add {
                    super::super::add_mod_u32(original2[i], delta[i], primes[i])
                } else {
                    super::super::sub_mod_u32(original2[i], delta[i], primes[i])
                };
                assert_eq!(actual1[i], first.min(second), "lane {i}, add={add}");
                assert_eq!(actual2[i], first.max(second), "lane {i}, add={add}");
            }
        }
    }
}
