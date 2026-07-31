//! WebAssembly SIMD report scanning for the biased SIQS score stream.
//!
//! `i8x16.bitmask` is WebAssembly's direct equivalent of the x86 movemask
//! used by the native kernel: it rejects sixteen non-survivors per branch and
//! leaves only set lanes for the scalar exact-threshold check.

use core::arch::wasm32::*;

#[inline]
pub(super) fn collect_candidates(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    debug_assert!(threshold >= 128);
    candidates.clear();
    // SAFETY: the SIMD build opts the entire artifact into simd128 and the
    // kernel bounds every unaligned load by the supplied slice length.
    unsafe { collect_simd(scores, threshold, candidates) }
}

#[target_feature(enable = "simd128")]
unsafe fn collect_simd(scores: &[u8], threshold: u8, candidates: &mut Vec<u32>) {
    let mut position = 0usize;
    while position + 16 <= scores.len() {
        // SAFETY: the loop condition proves this sixteen-byte load is within
        // the score slice. WebAssembly v128 loads permit unaligned addresses.
        let values = unsafe { v128_load(scores.as_ptr().add(position).cast()) };
        let mut mask = i8x16_bitmask(values) as u32;
        while mask != 0 {
            let lane = mask.trailing_zeros() as usize;
            let candidate = position + lane;
            if threshold == 128 || scores[candidate] >= threshold {
                candidates.push(candidate as u32);
            }
            mask &= mask - 1;
        }
        position += 16;
    }
    for (offset, &score) in scores[position..].iter().enumerate() {
        if score >= threshold {
            candidates.push((position + offset) as u32);
        }
    }
}
