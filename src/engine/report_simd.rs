//! x86-64 report scanning for the biased SIQS score stream.
//!
//! Every accepted score has its high bit set. SSE2/AVX2 movemasks therefore
//! reject 16/32 positions at once and enumerate only the rare survivors.

use core::arch::x86_64::*;

#[inline]
pub(super) fn collect_candidates(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    debug_assert!(threshold >= 128);
    candidates.clear();
    // SAFETY: feature detection guards AVX2 and x86-64 guarantees SSE2. Both
    // kernels use unaligned loads only within the supplied slice.
    unsafe {
        if is_x86_feature_detected!("avx2") {
            collect_avx2(scores, threshold, candidates)
        } else {
            collect_sse2(scores, threshold, candidates)
        }
    }
}

#[inline(always)]
fn append_mask(
    scores: &[u8],
    threshold: u8,
    start: usize,
    mut mask: u32,
    candidates: &mut Vec<u32>,
) {
    while mask != 0 {
        let lane = mask.trailing_zeros() as usize;
        let position = start + lane;
        if threshold == 128 || scores[position] >= threshold {
            candidates.push(position as u32);
        }
        mask &= mask - 1;
    }
}

#[target_feature(enable = "avx2")]
unsafe fn collect_avx2(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    let mut position = 0usize;
    while position + 32 <= scores.len() {
        // SAFETY: the loop condition proves the 32-byte load is in bounds;
        // loadu accepts the byte alignment of the score allocation.
        let values = unsafe { _mm256_loadu_si256(scores.as_ptr().add(position).cast()) };
        let mask = _mm256_movemask_epi8(values) as u32;
        if mask != 0 {
            append_mask(scores, threshold, position, mask, candidates);
        }
        position += 32;
    }
    collect_tail(scores, threshold, position, candidates);
}

#[target_feature(enable = "sse2")]
unsafe fn collect_sse2(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    let mut position = 0usize;
    while position + 16 <= scores.len() {
        // SAFETY: the loop condition proves the 16-byte load is in bounds.
        let values = unsafe { _mm_loadu_si128(scores.as_ptr().add(position).cast()) };
        let mask = _mm_movemask_epi8(values) as u32;
        if mask != 0 {
            append_mask(scores, threshold, position, mask, candidates);
        }
        position += 16;
    }
    collect_tail(scores, threshold, position, candidates);
}

fn collect_tail(scores: &[u8], threshold: u8, start: usize, candidates: &mut Vec<u32>) {
    for (offset, &score) in scores[start..].iter().enumerate() {
        if score >= threshold {
            candidates.push((start + offset) as u32);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispatched_scan_matches_scalar_at_boundaries_and_deep_thresholds() {
        for len in [0, 1, 15, 16, 17, 31, 32, 33, 127] {
            let scores: Vec<u8> = (0..len)
                .map(|index| (index as u8).wrapping_mul(37).wrapping_add(91))
                .collect();
            for threshold in [128, 137, 200, 255] {
                let expected: Vec<u32> = scores
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &score)| (score >= threshold).then_some(index as u32))
                    .collect();
                let mut actual = Vec::new();
                collect_candidates(&scores, threshold, &mut actual);
                assert_eq!(actual, expected, "len={len}, threshold={threshold}");
            }
        }
    }

    #[test]
    #[ignore = "manual SIMD-versus-scalar report-scan performance measurement"]
    fn profile_report_scan() {
        let mut scores = vec![12u8; 524_288];
        for position in (7_919..scores.len()).step_by(65_521) {
            scores[position] = 128;
        }
        let mut simd = Vec::new();
        let mut scalar = Vec::new();
        let started = std::time::Instant::now();
        for _ in 0..2_000 {
            collect_candidates(&scores, 128, &mut simd);
        }
        let simd_time = started.elapsed();
        let started = std::time::Instant::now();
        for _ in 0..2_000 {
            scalar.clear();
            for (position, &score) in scores.iter().enumerate() {
                if score >= 128 {
                    scalar.push(position as u32);
                }
            }
        }
        let scalar_time = started.elapsed();
        assert_eq!(simd, scalar);
        eprintln!("PROFILE report_scan simd={simd_time:?} scalar={scalar_time:?}");
    }
}
