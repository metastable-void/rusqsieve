use super::*;

pub(super) fn sieve_family(
    ctx: &Context,
    family: u64,
    scratch: &mut EngineScratch,
) -> FamilyResult {
    let empty = |family| FamilyResult {
        family,
        polynomials: 0,
        relations: Vec::new(),
        survivors: 0,
        timing: [0; 4],
    };
    let family_started = ctx.profile.then(std::time::Instant::now);
    let Some((a, aidx)) = choose_a(ctx, family) else {
        return empty(family);
    };
    let base = &ctx.base;
    let nfb = base.len();
    let s = aidx.len();
    let nvar = (s - 1).min(9); // number of sign bits varied per family
    let variants = 1u64 << nvar;

    // SIQS B-values: b = Σ ±Bⱼ, with Bⱼ ≡ ±sqrt(n) (mod qⱼ), 0 (mod other q).
    // Keep the true signed B instead of reducing it modulo A. This is the standard self-init
    // representation used by FLINT and makes every Gray-code root update one conditional
    // add/subtract, without a per-prime correction for modular-A wraps.
    let mut bvals: Vec<Natural> = Vec::with_capacity(s);
    for &i in &aidx {
        let q = base[i as usize].prime;
        let Some((ap, _)) = a.div_rem_u64(q as u64) else {
            return empty(family);
        };
        let Some(apinv) = inv_u32(ap.mod_u64(q as u64) as u32, q) else {
            return empty(family);
        };
        let mut coeff = (base[i as usize].sqrt_n as u64 * apinv as u64) % q as u64;
        coeff = coeff.min(q as u64 - coeff);
        bvals.push(ap.checked_mul(&Natural::from_u64(coeff)).unwrap());
    }
    let mut b = Natural::ZERO;
    for bj in &bvals {
        b = b.checked_add(bj).unwrap();
    }
    let mut bneg = false;
    let two_full: Vec<Natural> = bvals[..nvar].iter().map(|bj| bj.wrapping_add(bj)).collect();

    // Per-prime precompute for the initial polynomial: both roots and, for each
    // varying B-value, the O(1) root advance `2·Bⱼ·a⁻¹ mod p`.
    scratch.root1.clear();
    scratch.root1.resize(nfb, u32::MAX);
    scratch.root2.clear();
    scratch.root2.resize(nfb, 0);
    scratch.bainv.clear();
    scratch.bainv.resize(nvar * nfb, 0);
    for (idx, e) in base.iter().enumerate() {
        let p = e.prime;
        if p == 2 {
            continue;
        }
        let ap = a.mod_u64(p as u64) as u32;
        if ap == 0 {
            continue; // p | a: linear fallback per polynomial (root1 stays MAX)
        }
        let Some(ainvp) = inv_u32(ap, p) else {
            continue;
        };
        let mut bp = b.mod_u64(p as u64) as u32;
        if bneg && bp != 0 {
            bp = p - bp;
        }
        let xroot1 = mulmod_u32((e.sqrt_n + p - bp) % p, ainvp, p);
        let xroot2 = mulmod_u32(((p - e.sqrt_n) % p + p - bp) % p, ainvp, p);
        let r1 = add_mod_u32(xroot1, ctx.interval_mod_p[idx], p);
        let r2 = add_mod_u32(xroot2, ctx.interval_mod_p[idx], p);
        scratch.root1[idx] = r1.min(r2);
        scratch.root2[idx] = r1.max(r2);
        for (j, bj) in bvals.iter().take(nvar).enumerate() {
            let bjp = bj.mod_u64(p as u64) as u32;
            let two_bjp = (2 * bjp as u64 % p as u64) as u32;
            scratch.bainv[j * nfb + idx] = mulmod_u32(two_bjp, ainvp, p);
        }
    }

    // Sieve every polynomial in Gray-code order, advancing the roots in O(1) per
    // prime between consecutive polynomials instead of recomputing them.
    let mut relations = Vec::new();
    let mut survivors = 0u64;
    let mut timing = [0u64; 4];
    if let Some(started) = family_started {
        timing[0] = started.elapsed().as_nanos().min(u64::MAX as u128) as u64;
    }
    for v in 0..variants {
        survivors += sieve_one_poly(
            ctx,
            &a,
            &b,
            bneg,
            &aidx,
            &scratch.root1,
            &scratch.root2,
            &mut scratch.scores,
            &mut scratch.candidates,
            &mut relations,
            &mut timing,
        ) as u64;
        if v + 1 >= variants {
            break;
        }
        let j = (v + 1).trailing_zeros() as usize;
        let gray = v ^ (v >> 1);
        let flip_to_one = (gray >> j) & 1 == 0;
        let add_bainv = if flip_to_one {
            (b, bneg) = signed_add(&b, bneg, &two_full[j], true);
            true
        } else {
            (b, bneg) = signed_add(&b, bneg, &two_full[j], false);
            false
        };
        let off = j * nfb;
        let roots_started = ctx.profile.then(std::time::Instant::now);
        if add_bainv {
            for idx in 0..nfb {
                if scratch.root1[idx] == u32::MAX {
                    continue;
                }
                let p = base[idx].prime;
                let d = scratch.bainv[off + idx];
                let r1 = add_mod_u32(scratch.root1[idx], d, p);
                let r2 = add_mod_u32(scratch.root2[idx], d, p);
                scratch.root1[idx] = r1.min(r2);
                scratch.root2[idx] = r1.max(r2);
            }
        } else {
            for idx in 0..nfb {
                if scratch.root1[idx] == u32::MAX {
                    continue;
                }
                let p = base[idx].prime;
                let d = scratch.bainv[off + idx];
                let r1 = sub_mod_u32(scratch.root1[idx], d, p);
                let r2 = sub_mod_u32(scratch.root2[idx], d, p);
                scratch.root1[idx] = r1.min(r2);
                scratch.root2[idx] = r1.max(r2);
            }
        }
        if let Some(started) = roots_started {
            timing[3] += started.elapsed().as_nanos().min(u64::MAX as u128) as u64;
        }
    }
    FamilyResult {
        family,
        polynomials: variants,
        relations,
        survivors,
        timing,
    }
}

pub(super) fn signed_add(a: &Natural, aneg: bool, b: &Natural, bneg: bool) -> (Natural, bool) {
    if aneg == bneg {
        // A is selected near sqrt(kN)/interval and |B| remains below 5A, so
        // this intermediate retains roughly log2(interval)-3 bits of headroom
        // (more than 25 bits at every shipped tier).
        debug_assert!(
            a.checked_add(b).is_some(),
            "signed SIQS coefficient has tier-proved headroom"
        );
        let sum = a.wrapping_add(b);
        let neg = aneg && !sum.is_zero();
        (sum, neg)
    } else if a >= b {
        let diff = a.wrapping_sub(b);
        let neg = aneg && !diff.is_zero();
        (diff, neg)
    } else {
        let diff = b.wrapping_sub(a);
        let neg = bneg && !diff.is_zero();
        (diff, neg)
    }
}

pub(super) fn build_a_candidates(
    base: &[FactorBaseEntry],
    target_a: &Natural,
) -> (Arc<[usize]>, Arc<[usize]>, usize) {
    let target_bits = target_a.bit_len();
    let factor_count = target_bits.div_ceil(14).clamp(3, 10);
    let ideal_bits = target_bits.div_ceil(factor_count);
    let minimum_bits = ideal_bits.saturating_sub(1).max(2);
    let all: Vec<usize> = base
        .iter()
        .enumerate()
        .filter(|(_, e)| (32 - e.prime.leading_zeros()) as usize >= minimum_bits)
        .map(|(i, _)| i)
        .collect();
    if all.len() < factor_count {
        return (all.into(), Arc::from([]), factor_count);
    }
    let mut window = 1usize;
    let pool = loop {
        let candidates = all
            .iter()
            .copied()
            .filter(|&i| {
                let bits = (32 - base[i].prime.leading_zeros()) as usize;
                bits.abs_diff(ideal_bits) <= window
            })
            .collect::<Vec<_>>();
        if candidates.len() >= factor_count * 2 || window >= 31 {
            break candidates;
        }
        window += 1;
    };
    debug_assert!(!pool.is_empty(), "choose_a constraints must be satisfiable");
    (all.into(), pool.into(), factor_count)
}

pub(super) fn choose_a(ctx: &Context, family: u64) -> Option<(Natural, Vec<u32>)> {
    let all = &ctx.a_all;
    let pool = &ctx.a_pool;
    let factor_count = ctx.a_factor_count;
    if all.len() < factor_count || pool.len() < factor_count {
        return None;
    }
    let mut state = family ^ 0x9e3779b97f4a7c15;
    let mut best = None;
    for _ in 0..32 {
        let mut a = Natural::ONE;
        let mut idx = Vec::with_capacity(factor_count);
        while idx.len() + 1 < factor_count {
            state = crate::u64math::xorshift(state);
            let i = pool[state as usize % pool.len()];
            if idx.contains(&(i as u32)) {
                continue;
            }
            a = a.checked_mul(&Natural::from_u64(ctx.base[i].prime as u64))?;
            idx.push(i as u32)
        }
        let desired_u64 = ctx.target_a.div_rem(&a)?.0.to_u64()?;
        let last = all
            .iter()
            .copied()
            .filter(|&i| !idx.contains(&(i as u32)))
            .min_by_key(|&i| (ctx.base[i].prime as u64).abs_diff(desired_u64))?;
        a = a.checked_mul(&Natural::from_u64(ctx.base[last].prime as u64))?;
        idx.push(last as u32);
        let close = a
            .checked_mul(&Natural::from_u64(5))
            .zip(ctx.target_a.checked_mul(&Natural::from_u64(4)))
            .is_some_and(|(lhs, rhs)| lhs >= rhs)
            && ctx
                .target_a
                .checked_mul(&Natural::from_u64(5))
                .zip(a.checked_mul(&Natural::from_u64(4)))
                .is_some_and(|(lhs, rhs)| lhs >= rhs);
        if close {
            return Some((a, idx));
        }
        let distance = if a >= ctx.target_a {
            a.wrapping_sub(&ctx.target_a)
        } else {
            ctx.target_a.wrapping_sub(&a)
        };
        if best
            .as_ref()
            .is_none_or(|(prior, _, _): &(Natural, Natural, Vec<u32>)| distance < *prior)
        {
            best = Some((distance, a, idx));
        }
    }
    best.map(|(_, a, idx)| (a, idx))
}
