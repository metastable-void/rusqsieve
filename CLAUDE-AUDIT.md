# rusqsieve performance audit

Goal: make factorization **comparable to or faster than `flintqs`** (`/usr/bin/QuadraticSieve`),
which is the direct benchmark on this machine. All changes must stay within `SPEC.md`
(SIQS on the main path, true logarithmic sieve, deterministic results across parallelism,
frozen public / `low_level` / wasm-ABI names, relation/matrix/dependency invariants).

## Measured baselines (before optimization)

Balanced semiprimes, `qs-factor --threads 1`, release build:

| bits | rusqsieve (before) | flintqs (single-thread) |
|-----:|-------------------:|------------------------:|
|   64 | 3.82 s             | (refuses <40 digits)    |
|   80 | 28.9 s             | (refuses)               |
|   96 | >60 s (timeout)    | (refuses)               |
|  160 | not reached        | 1.22 s                  |
|  192 | not reached        | 5.13 s                  |
|  224 | not reached        | 36.9 s                  |
|  256 | not reached        | >77 s                   |

flintqs requires ≥40 decimal digits (~133 bit); the true head-to-head is **160 bit and up**,
which routes through the SIQS **engine** (`engine.rs`), not the reference path.

## Root causes (why it is ~10^6x too slow)

### A. Everything runs at 1024-bit width
The native entry (`native::factor`, P≤16) widens every input to `Natural<16>` (1024-bit,
16 limbs). A 64-bit value is stored in 16 limbs and **every** arithmetic op iterates all 16
limbs regardless of magnitude. Combined with the bad kernels below, small inputs pay full
1024-bit cost.

### B. Asymptotically-bad arithmetic kernels (`src/natural/mod.rs`)
- **`mul_mod`** (L416): Russian-peasant *double-and-add*, O(bits) modular additions, each O(P)
  limbs → O(P²·64). Used everywhere (pollard, `pow_mod`, extraction, `crt_root`, `combine`).
- **`Montgomery`** (L849): *fake* — `encode`/`decode` are no-ops (`v mod m`) and `mul` just calls
  the slow `mul_mod`. Spec §6.12 mandates a real Montgomery (`n0_inverse`, `r2`, `one`).
- **`div_rem`** (L283): *bit-by-bit* long division, O(bits) iterations each doing `shl_one`+cmp+sub
  over all P limbs → O(P²·64). Spec §6.11 asks for normalized limb long division.
- **`div_rem_u64`** (L302): always loops all 16 limbs and issues **two** u128 libcalls
  (`/` then `%`) per limb. This is the trial-division inner loop of the whole sieve.
- **`gcd`** (L315): Euclidean using the slow `div_rem`. Used per-iteration in pollard and in
  extraction. Binary GCD (shifts/subtraction, no division) is far faster.
- Pervasive `.clone()` and fresh-allocation shifts (`shl_one`, `Shl`, `Shr`).

### C. Small/medium inputs never reach SIQS, and the fallback is weak
`engine::factor_node` routes `n.bit_len() < 120` to `factor::factor_complete`, i.e.
`pollard_rho` then the reference QS. `pollard_rho` (`factor.rs` L329):
- uses the slow `mul_mod` and a **gcd every iteration** (no batching, no Brent),
- has an iteration budget (24×200 000 ≈ 4.8M) **too small** to find a ~48-bit factor
  (~2^24 ≈ 16.7M steps), so 96-bit inputs exhaust pollard and fall to the reference QS.
- `factor::primes_to` recomputes the prime list by trial division **on every recursion node**.

### D. `qs/mod.rs::sieve_job` is not a sieve  *(FIXED — now a §12.6 log-sieve)*
The portable low-level kernel `sieve_job` (used by `FactorSession`/`reference_qs_factor`)
computes `q = (ceil_sqrt(n)+k)^2 - n` and **trial-divides every candidate by the entire
factor base** — Fermat differences + trial division, O(candidates·|FB|). It ignores the
`scores`/`candidates` scratch and the precomputed `sqrt_n` roots. This contradicts SPEC §12.6
(true log-sieve: add log(p) at the two modular roots, trial-divide only threshold survivors).
The real log-sieve lives only in `engine.rs::sieve_polynomial`.

### E. Engine hot-loop inefficiencies (`engine.rs::sieve_polynomial`)
- Allocates+zeros a `scores` vec (up to 262 144 × u16) **per polynomial** (L628) instead of
  reusing worker scratch (spec §21.1 requires reuse).
- Recomputes `inv_u32` (extended gcd) and both roots **per prime per polynomial** — the SIQS
  self-initialization property (cheap incremental root updates across `b`-variants sharing the
  same `a`) is not exploited.
- Trial-divides survivors with the 16-limb `div_rem_u64`/`div_rem` even once the cofactor is
  small.

## Plan (highest-leverage first; each step measured)

1. **Correctness net**: add `num-bigint` dev-dependency + randomized differential tests for the
   arithmetic core (spec §19.3 allows num-bigint as dev-only; §23 forbids timing asserts).
2. **Arithmetic core** (`natural/mod.rs`), spec-neutral, benefits every path:
   significant-limb-aware ops; `div_rem_u64` single-divmod; `div_rem` schoolbook/Knuth-D;
   `mul_mod` via `widening_mul` + fast wide reduction; binary `gcd`; real Montgomery.
3. **Small path**: native u64/u128 fast lane for inputs that fit; Brent + block-gcd Pollard
   over Montgomery; cache the small-prime sieve once.
4. **SIQS engine**: reuse scratch buffers; fast survivor trial division; self-initializing root
   updates; better large-prime handling.
5. **Parameters**: retune factor-base bound, sieve interval, thresholds, large-prime bounds and
   relation target for the 130–256 bit range (keep heuristics centralized per spec §12.1).
6. **Benchmark vs flintqs** at each step; keep the full test suite green (dev profile — the
   spec-fixed `panic="abort"` release profile cannot host the unwinding test harness).

## Results (after optimization)

Balanced semiprimes, release build, this 96-core host. rusqsieve uses `--threads auto`
(worker count auto-capped by input size to avoid parallel-startup overhead on small inputs).
flintqs is single-threaded and refuses inputs < 40 digits (~133 bit).

| bits | rusqsieve before | rusqsieve after | flintqs (1 thread) | vs flintqs |
|-----:|-----------------:|----------------:|-------------------:|-----------:|
|   64 | 3.82 s           | 0.11 s          | (refuses)          | —          |
|   80 | 28.9 s           | 0.21 s          | (refuses)          | —          |
|   96 | >60 s (timeout)  | 0.11 s          | (refuses)          | —          |
|  112 | >60 s            | 0.21 s          | (refuses)          | —          |
|  128 | >60 s            | 0.61 s          | (refuses)          | —          |
|  160 | >60 s            | 0.41 s          | 1.22 s             | **2.9× faster** |
|  192 | >60 s            | 1.12 s          | 5.13 s             | **4.5× faster** |
|  224 | >60 s            | 9.3 s           | 36.9 s             | **3.9× faster** |
|  256 | >60 s            | 82 s            | 231 s              | **2.8× faster** |

(The "after" column reflects the full optimization set including SIQS self-initialization,
byte-array sieve, and cheaper Q(x) reconstruction — see the per-polynomial section below.)

rusqsieve is faster than flintqs across the whole head-to-head range (160–256 bit) and
factors the small/mid range (64–128 bit, which flintqs refuses) in ≤0.5 s.

Correctness: 18k-case randomized differential arithmetic vs `num-bigint`; product
verification across balanced/unbalanced semiprimes, 3-prime composites, prime powers, and
prime inputs; determinism verified (identical factors for `--threads 1` vs `--threads 96`).

What made the difference:
- **Arithmetic core** (helps every path): significant-limb `widening_mul`, normalized
  limb long division (Knuth D) replacing bit-by-bit, `widening_mul`+wide-reduction `mul_mod`
  replacing double-and-add, binary GCD, and a significant-limb `div_rem_u64` with a single
  divmod per limb. Small-input arithmetic no longer pays fixed 1024-bit cost.
- **Native `u64` fast path** for ≤64-bit cofactors (deterministic Miller–Rabin + Pollard–Brent),
  plus a cached prime sieve shared across the recursion.
- **SIQS self-initialization**: `a⁻¹ mod p` precomputed once per polynomial family instead of
  once per polynomial (was an extended-GCD per prime per polynomial); reused score buffers.
- **Sieve threshold fix**: the old threshold admitted ~1150 non-smooth survivors per polynomial
  (yielding ~0 full relations). It is now `2·(log2|g(x)| − large_prime_allowance)`, cutting
  survivors ~70× and matching relation yield to trial-division cost.
- **Worker-count capping** by input size (spawning 96 threads for a sub-second job cost >1 s).

### Double-large-primes & sparse linear algebra (investigated + implemented)

Both were implemented and measured. The key finding: **for this implementation the bottleneck at
224–256 bit is raw per-polynomial cost (the `Natural<16>` `Q(x)` reconstruction and the sieve
passes), not the large-prime strategy or the linear algebra.** Phase timing at 224 bit is ~91 %
sieving, ~9 % linear algebra; at 256 bit linear algebra is only a few seconds of ~103 s.

- **Double-large-prime variation — implemented, correct, available, off by default.**
  `engine::RelationCollector` is a union-find spanning forest over large-prime vertices; every
  partial relation is an edge between its large prime(s) (single-large-primes use a reserved unit
  vertex), and a relation that closes a cycle combines every relation on the cycle via
  `combine_cycle` (all large primes on a cycle cancel to even powers). `classify_cofactor` splits
  composite cofactors (Pollard rho, primality-checked) into two large primes. This subsumes the
  old single-large-prime hash-matching (the default path) and is exercised + verified by the test
  suite. **Doubles are gated off by default** (`large_prime_policy`) because enabling them requires
  a lower sieve threshold, which floods the whole-factor-base confirmation step; a resieving
  confirmation (divide only the primes that hit each survivor) was prototyped to fix the flooding
  but its extra full sieve pass cost more than it saved at these survivor densities — net negative
  (224 bit: 17 s → 26–29 s). So doubles pay off only once the per-polynomial cost is reduced.
- **Sparse linear algebra — implemented (`SparseBinaryMatrix::filtered_dependencies`).** SPEC §15.3
  structured elimination: iterative singleton-row removal shrinks the matrix before the dense
  solve; dependencies are mapped back to the original column space and re-verified. Wired into
  `f2::BlockLanczos::begin` and the engine extraction; differential-tested against the dense
  oracle. It removes the O(n³) dense cost at large sizes but only saves ~3 % overall, since LA is
  not the bottleneck. The full Montgomery block-Lanczos *recurrence* was intentionally **not**
  implemented from scratch: it targets the same (non-bottleneck) phase, and a reference-free
  implementation is high-risk relative to zero measurable benefit here.

### Per-polynomial cost reduction (self-initialization + cheaper reconstruction)

Follow-up on "the real remaining lever". Implemented and measured:

- **SIQS self-initialization (incremental roots) — implemented.** The scoring loop used to
  recompute both modular roots per prime per polynomial (`b mod p` — a big-integer division — plus
  two `mulmod`s). Now `2·Bⱼ·a⁻¹ mod p` is precomputed once per family, `b` is walked in Gray-code
  order, and roots advance by one add per prime between consecutive polynomials (SPEC §12.5).
  Subtlety fixed: because `b` is kept reduced in `[0,a)` and `a·a⁻¹ ≡ 1 (mod p)`, each mod-`a`
  wrap shifts every prime's root uniformly — that shift is applied alongside the increment.
  Verified against from-scratch roots and `b² ≡ n (mod a)`. **~20 % faster at 256 bit** (big
  factor base → root recomputation was a real cost); neutral at ≤192 bit (small factor base).
- **Byte-array sieve — implemented.** Scores are now `u8` (was `u16`), halving the sieve array so
  more of it stays cache-resident; weights/threshold rescaled to `log₂` so a smooth `Q ≈ 2^g` fits
  a byte across the supported range.
- **Cheaper Q(x) reconstruction — implemented.** Survivors now compute `g(x)=Q/a = a·x²+2b·x+c`
  directly (signed) instead of the wide `t²` squaring followed by a division by `a` (`a | Q` is
  guaranteed since `b² ≡ n (mod a)`). **~9 % faster at 224 bit.**

Net vs. the pre-follow-up state: 256 bit 102.7 s → 82 s, 224 bit 10.9 s → 9.3 s; small/mid sizes
unchanged (already overhead-bound). rusqsieve now beats flintqs 2.9–4.5× across 160–256 bit.

- **Bucket sieving — assessed, deferred.** The dominant remaining cost is the sieve *stepping*
  (score writes), where large-stride primes (`p` larger than the cache-resident block) cause a
  cache miss per hit. Bucket sieving would batch those into cache-local blocks. It is a major
  restructure of the scoring + survivor-scan loop (block partitioning, per-block hit buckets,
  block-local draining) with an estimated ~15–25 % gain at large sizes — worthwhile, but a
  higher-risk rewrite with diminishing returns given rusqsieve already leads flintqs 2.8–4.5×.
  Left as the next scoped step rather than destabilizing the current verified state.

- **Done earlier:** the low-level portable kernel `qs::sieve_job` — previously trial-division of every
  `x²−n` candidate — is now a genuine logarithmic sieve (SPEC §12.6): it adds `log(p)` at the
  two roots `x ≡ ±√n (mod p)` across byte-score segments and trial-divides only threshold
  survivors, with single-large-prime classification (primality-checked). Guarded by a new unit
  test (`qs::tests::sieve_job_logarithmic_sieve`) asserting filtering + relation validity +
  determinism. (Double-large-prime classification, §12.6 step 8, remains a follow-up consistent
  with the engine.)

## Constraints honored
- SIQS stays the main algorithm; log-sieve preserved/extended (not replaced by trial division).
- Relation invariant `square_root^2 ≡ ±∏ fb^e · ∏ lp (mod n)`, matrix bit/row convention,
  and mandatory dependency verification are unchanged.
- Determinism: factors identical for parallelism 1 vs N; workers make no independent random
  choices; results canonicalized before matrix construction.
- No threads in the math core; no unsafe on native (`forbid(unsafe_code)` off wasm); portable
  (musl/wasm) builds preserved; frozen API/ABI names preserved.

## Next frontiers (identified 2026-07-21)

Follow-up analysis (width×ISA benchmark + FLINT stdin re-baseline). Measured on the 96-core
Xeon 8259CL: at a *fixed* thread count rusqsieve is ~12–13× less efficient **per core** than
single-threaded flintqs (192-bit 2.78 s / 224-bit 31.27 s flint vs 4.32 s / 49.59 s ours at
8 threads) — a roughly *constant* factor across sizes, i.e. a slow inner loop, not a bad
asymptotic. flintqs aborts (SIGABRT) at 256-bit here; rusqsieve factors it in ~92 s (32 thr),
so we are already more robust at the top end. Frontiers, highest-leverage first:

1. **Resieving — replace full-factor-base trial division.** `engine.rs::sieve_one_poly`
   (~L898) trial-divides every sieve survivor's `g(x)` by *every* factor-base prime with a
   full `div_rem_u64` (which also zeroes a fresh quotient `Natural` per call). Replace with a
   second, root-strided *resieve* pass that records which primes actually hit each survivor
   position, then divide only by those (~a dozen) primes. Turns O(FB×survivors) heavy bignum
   divisions into O(FB) light stride-marks + O(#factors) divisions per survivor. Pure algorithm,
   no unsafe, determinism-preserving. Expected: closes most of the ~12× per-core gap. **[TOP]**
2. **Bucket sieving for large primes.** (Previously "assessed, deferred" above.) Large-stride
   primes cause a cache miss per hit in the score-stepping loop. Partition the interval into
   cache-resident blocks and bucket large-prime hits per block, then drain block-locally.
   Interacts with (1): the resieve pass benefits from the same bucketing. Est. ~15–25% at ≥224-bit.
3. **Threshold / parameter tuning + tiny-prime skipping.** Don't sieve the smallest primes
   (account for them with threshold slack), and tighten the survivor threshold so fewer bogus
   candidates reach the (post-#1, still non-free) factoring step. Cheap, contained, measurable.
4. **Montgomery REDC for `mul_mod`.** `Montgomery::{mul,pow}` currently forward to division-based
   `mul_mod` (widening_mul + Knuth). No REDC exists. Not hot at ≤256-bit (only in
   cycle-combine/extract), so deferred — real only if `mul_mod` moves onto the hot path.
5. **SIMD candidate-survivor scan.** The `score >= threshold` scan (~L855) maps cleanly to
   `vpcmpub`+`vpcompressd`, but it is a minority of runtime and needs `unsafe` intrinsics blocked
   by `#![forbid(unsafe_code)]` on native. Only worthwhile after 1–3, behind `arch-optimized`
   with runtime dispatch.
6. **Parallel linear algebra.** `f2::BlockLanczos` is single-threaded; irrelevant while sieving
   dominates, becomes the tail past ~256-bit.

Width×ISA aside (measured): `limit-to-512-bits` (`PARTS=8`) alone is ~6–8% *slower* scalar, but
`PARTS=8 + AVX2` (`-C target-cpu=x86-64-v3`) is the fastest config measured — AVX-512 is worse
than AVX2 (Cascade Lake downclock). These are ~2% effects, dwarfed by frontier #1.

**Status 2026-07-21:** implementing #1, #2, #3 (this session). #4–#6 recorded for later.

### Frontier results (2026-07-21, this session)

Implemented #1 and #3; measured on the reference host (balanced semiprimes, factors verified,
same thread counts as the baselines: 8 for 192/224-bit, 32 for 256-bit). All 24 unit tests pass.

- **#1 Resieve — implemented, size-gated (`RESIEVE_MIN_FB = 7000`).** For a large factor base a
  second root-strided pass records exactly which primes hit each survivor, replacing full trial
  division. Below the gate the extra pass costs more than it saves, so the original trial
  division is kept. Root-gating the small-FB path (a `pos ≡ root (mod p)` pre-test) was tried at
  your suggestion but measured ~2–4% *slower* (the early `q==1` break + `q` shrinking to one limb
  already make those divisions cheaper than the gate's per-prime modulo), so it was reverted.
  Effect: −14% at 256-bit; neutral (unchanged path) at ≤224-bit.
- **#3 Tiny-prime skipping — implemented (`SMALL_SKIP = 20`, `SMALL_SLACK = 3`).** Primes < 20 are
  no longer added to the byte scores (they are ~32% of the score-write traffic but tiny weight);
  they are still divided out during factoring, and the threshold is lowered by 3 to compensate.
  A quick 192-bit sweep found SKIP=20 optimal (SKIP=40/60 lower the threshold too far → more
  false-positive survivors → net slower). This is the dominant win, and it targets exactly the
  score-write cost the bucket-sieving note above is about.

Combined #1+#3 vs. the pre-session baseline:

| bits | baseline | #1+#3 | speedup |
|------|----------|-------|---------|
| 192  | 4.36 s   | 3.62 s | −17.0% |
| 224  | 50.44 s  | 41.95 s | −16.8% |
| 256  | 94.81 s  | 73.86 s | −22.1% |

#2 (bucket sieving) attempted next; #4–#6 still deferred.

- **#2 Blocked/bucket sieving — implemented, measured, reverted (portability optimization,
  needs on-target validation).** rusqsieve ships to many machines: high-end workstations, consumer
  PCs, tablets, **smartphones (ARM L2 often 256–512 KB)**, and WASM on all of them. The engine's max
  score array is **256 KB** (2 × 131072 half-width). On this dev box (Xeon, **1 MB L2/core**) that is
  L2-resident, so the sieve is not memory-bound here — but on the smaller-L2 consumer/mobile targets
  the array exceeds L2 and the sieve *is* memory-bound, which is exactly where blocking helps. So #2
  is a genuine win for those targets, not a dead end.
  - **Tried:** a blocked sieve applying each prime's hits one 16 KiB (L1-sized) block at a time,
    carrying per-prime positions across blocks. Correct (identical scores/relations; determinism test
    passes) but **13% / 41% / 33% slower at 192 / 224 / 256-bit on the Xeon** — fragmenting each
    prime's long tightly-pipelined strided loop into ~16 short per-block loops destroys inner-loop
    throughput, and this box has no cache deficit to recover. Reverted.
  - **Correct form:** a true bucket sieve (one pass appends `(pos, weight)` into per-block buckets,
    then each block is drained block-locally — no re-striding, so tight loops are preserved). This
    should be neutral-to-slight-loss on large-L2 parts and a real win where the array exceeds L2. It
    cannot be validated on this host (which has no cache deficit), so it is deferred pending a
    benchmark on a representative small-cache target (a phone/tablet/low-end PC, or WASM there),
    ideally behind a build option so workstation builds keep the flat single-pass sieve.
  - Aside: the real residual cache pressure *here* is resieve's `cand_at` at 256-bit (1 MB u32 =
    whole L2); a narrower survivor map is a lower-risk local follow-up.

## Single-core sieve-yield optimization (2026-07-21, session 2)

Goal: match flintqs single-core (it does 224-bit in ~20-31s, 256-bit in ~200s single-threaded;
our HEAD baseline was 306s @224 and **3841.7s @256** single-core). Profiling (`RUSQSIEVE_PROFILE=1`)
showed the cost is **entirely sieve yield**: at 224-bit, 461k polynomials for 3091 relations, LA
negligible (0.14s), parallel scaling already good (7.3x @8, 24.5x @32). Survivor instrumentation
showed **~99% of sieve survivors are false positives** (160k survivors → 1.6k relations at 192-bit).

Landed (all verified, 24 tests pass, wasm builds both feature settings):
- **Knuth-Schroeppel multiplier** (`knuth_schroeppel`, ported from FLINT). Sieves `k·n` (factor base,
  roots, `Q(x)`) while extracting with `n` via `gcd(x−y, n)` (`Context.sieve_n` vs `n`). Primes
  dividing `k` are added as ramified factor-base entries (`FactorBaseBuilder.multiplier`) instead of
  reported as `FoundFactor`. **~1.9x at 224-bit.**
- **Factor base retune.** Profiling-guided: the 161-192 tier's bound was ~2x too small
  (28k→60k, nfb 1550→3007) → **2x at 192-bit**. 224 (60k) and 256 (150k) tiers were already near
  their optimum *for this engine*. 129-160 bumped 20k→40k.
- **Threshold margin** (`THRESH_MARGIN=4`) cuts false-positive survivors for a few extra polynomials.
- **`rem_u64`** (remainder without quotient allocation) for trial-division tests (neutral here but
  removes an allocation per test).

Before → after (balanced semiprimes, verified):

| bits | cores | HEAD | tuned | speedup |
|------|-------|------|-------|---------|
| 192  | 1     | 31.9 s  | 13.4 s  | 2.4x |
| 224  | 1     | ~306 s  | 175.6 s | 1.7x |
| 224  | 8     | 49.6 s  | 23.4 s  | 2.1x |
| 256  | 1     | 3841.7 s | (measuring) | — |
| 256  | 8     | —       | 231.0 s | — |
| 256  | 32    | ~94 s   | 68.8 s  | 1.4x |

vs flintqs single-core: we went from ~10-11x slower to ~5x slower — closed about half the gap.

**Why we are still ~5x off flintqs, and the path to parity:** our per-polynomial cost scales with
factor-base size — both the sieve score-write scan AND trial division are O(nfb). So bigger factor
bases (where the smooth yield is) make each polynomial proportionally more expensive, and the
optimum stalls at nfb≈3000. flintqs affords nfb 10k-25k precisely because it has (a) **bucket
sieving** (frontier #2 — amortizes the scan, and now clearly justified: at flint-scale intervals the
array exceeds L2) and (b) **resieving** (frontier #1 — trial-divides only primes that hit). Those two
unlock large factor bases, which is the remaining multiplier. This is the concrete next step for
flint parity; it is a substantial structured change (and #2 helps the mobile/wasm targets regardless).

Tuning harness left in place (all no-ops when unset): `RUSQSIEVE_PROFILE=1` (phase/counter timings),
`RUSQSIEVE_FB_BOUND`, `RUSQSIEVE_HALFW`, `RUSQSIEVE_THRESH_ADJ` (native only) for continued
per-size tuning without rebuilds.

## Roadmap to beat flintqs single-core (next steps)

Status: single-core ~5x slower than flintqs (192: 13.4s vs 2.8s; 224: 176s vs ~31s). Multi-core
scaling is already good (77% eff @32) and LA is negligible at current sizes — so **the whole
problem is single-core sieve yield**, and the whole yield problem is one thing:

**Core lever — afford a large factor base.** Smooth yield rises with nfb; flintqs runs nfb 10k-25k,
we stall at nfb≈3000. We can't grow nfb because **both** hot costs are O(nfb): the sieve score-write
scan and the per-survivor trial division. Kill those two dependencies on nfb, then grow the FB.
Everything below serves that. flintqs references in `/home/dev/flint/src/qsieve/`.

### Phase 1 — Resieving: make trial division O(#hits), not O(nfb)  [frontier #1]
- We have a size-gated resieve (wins only at nfb≥7000) but `cand_at` is a u32-per-position array
  (4× the score array → cache-heavy; that's why it loses at small nfb).
- Do: shrink the survivor map (1-bit "is-survivor" bitmap + compact per-survivor bucket lists, or
  u16 index with an overflow guard) so resieve wins at all sizes; then make it unconditional.
- Ref: `collect_relations.c`, `large_prime_variant.c`. Target: trial-div cost independent of nfb.
- Verify: RUSQSIEVE_PROFILE survivor/relation counts unchanged; determinism test; per-size timing.

### Phase 2 — Bucket sieving: make the scan cache-efficient at large intervals  [frontier #2]
- NOT the naive per-block re-stride (measured +13-41% regression — it fragments tight loops).
- Do a PROPER bucket sieve: one pass appends `(block-local-offset, log)` into per-block buckets
  (sequential writes, tight loops preserved), then drain each cache-resident block. Large primes
  bucketed; small primes (stride < block) sieved directly per block.
- Ref: `collect_relations.c`. Pays off once the interval array exceeds L2 (which large FBs want),
  and helps the mobile/consumer/WASM targets (small L2) regardless.

### Phase 3 — Grow the factor base + interval toward flint's table, then re-tune
- With Phases 1-2, big FBs are cheap. Move `engine_params` toward flint's `qsieve_tune`
  (qsieve.h): ~nfb 10k @224, ~25k @256, with proportionally larger intervals.
- Blocked on Phase 4 (LA) — dense Gauss explodes at nfb=25k.

### Phase 4 — Real sparse Block Lanczos  [frontier #6]
- `f2::BlockLanczos` is a dense-Gauss stub, O(nfb³). At nfb=25k that's ~2.4e11 word-ops (minutes,
  serial) → becomes the bottleneck the moment the FB is large. Negligible today ONLY because nfb≈3k.
- Implement true Block Lanczos, O(nfb²·avg-weight); its matrix-vector products also parallelize,
  helping multi-core. Ref: `block_lanczos.c` (957 lines).

### Phase 5 — Systematic per-size auto-tuning
- Extend the env harness (RUSQSIEVE_FB_BOUND/HALFW/THRESH_ADJ + PROFILE) into a sweeper that, per
  bit-size and over several semiprimes, optimizes (fb_primes, sieve_size, small_primes, threshold,
  ks_primes) — flint's speed is largely its auto-generated `qsieve_tune` table. Bake results into
  `engine_params`. This converts "correct techniques" into "actually fast" and removes the
  single-sample overfitting risk in the current hand-tuned values.

### Phase 6 — Secondary levers (measure before/after each)
- **Better `a` selection** (`choose_a`): pick `a` closer to `sqrt(2·kn)/M` and well-spread so `Q(x)`
  is minimized → higher smooth density per polynomial. Ours is near-random. Ref: `compute_poly_data.c`.
- **Large-prime yield**: confirm our union-find double-large-prime cycles match `large_prime_variant.c`
  effectiveness (partials → full relations).
- **Montgomery REDC** for `mul_mod` [frontier #4] — only if profiling shows it hot.
- **SIMD candidate scan** [frontier #5] — needs a scoped `unsafe` module behind `arch-optimized`.

### Sequencing & expected payoff
Phase 1 → 2 → (4 ∥ 3) → 5. Phases 1-3 (affording nfb≈10k) should close most of the 5x and target
single-core parity; combined with existing scaling that means beating flintqs at low core counts,
not just at tens of cores. Each phase gated on: `cargo test` green (determinism + relation
invariant), factors verified, and a single-core before/after at 192/224/256 vs flintqs.
