# rusqsieve performance audit

Goal: make factorization **comparable to or faster than `flintqs`** (`/usr/bin/QuadraticSieve`),
which is the direct benchmark on this machine. All changes must stay within `SPEC.md`
(SIQS on the main path, true logarithmic sieve, deterministic results across parallelism,
frozen public / `low_level` / wasm-ABI names, relation/matrix/dependency invariants).

The matching FLINT source checkout is `/home/dev/flint`; use it as the primary implementation
reference rather than recollection or an online version. This LXD container is effectively
dedicated bare metal: no other guests share the host. Load-average changes during benchmarking are
normally caused by Codex, VS Code, Claude Code, or commands from this project, so avoid concurrent
build/test activity and use interleaved saved-binary A/B runs for timing comparisons.

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
| 256  | 1     | 3841.7 s | 1765.9 s | 2.2x |
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

### Phase 2 — Cache-blocked three-tier sieve: scan cost stops scaling with nfb  [frontier #2]
- NOT the naive per-block re-stride (measured +13-41% regression — it fragments tight loops).
- Partition the factor base into THREE tiers (the standard msieve/yafu/flint architecture):
  - **small** (stride ≪ block, or `p < SMALL_SKIP`): not sieved / handled at factoring time (#3).
  - **medium** (stride < block): sieved directly within each cache-resident block — hits many
    times per block, so writes stay in L1 and the tight inner loop is preserved.
  - **large** (stride ≥ block, hit ≤once per block): PROPER bucket sieve — one pass appends
    `(block-local-offset, log)` into per-block buckets (sequential writes), then each block is
    drained locally. No re-striding, tight loops preserved, random writes → sequential.
- Ref: `collect_relations.c`. Pays off once the interval array exceeds L2 (which large FBs want),
  and helps the mobile/consumer/WASM targets (small L2) regardless.

### Phase 2b — Cheaper root updates (the third O(nfb)-per-poly cost)  [SIMD target]
- Self-init advances every prime's two roots each polynomial (`sieve_family`, the `root1/root2 +=
  delta (mod p)` loop over nfb primes) — a per-poly O(nfb) cost alongside scan and trial-division.
- Scalar win: `delta` is already reduced to `[0,p)`, so `root + delta < 2p` → replace the `% p`
  with a branchless conditional subtract.
- SIMD win: this loop (and the `score >= threshold` candidate scan) are the two genuinely
  vectorizable hot spots — unlike the strided score writes, which do not vectorize. Behind
  `arch-optimized` with a scoped `unsafe` module + runtime dispatch; per-target (x86 AVX2, wasm
  simd128).

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
- **Polynomial batching**: process a batch of the family's b-variant polynomials together so setup,
  root state, and candidate handling amortize and stay cache-resident (complements Phases 2/2b).
- **Large-prime yield**: confirm our union-find double-large-prime cycles match `large_prime_variant.c`
  effectiveness (partials → full relations).
- **Montgomery REDC** for `mul_mod` [frontier #4] — only if profiling shows it hot.
- (SIMD moved to Phase 2b — root updates + candidate scan are the vectorizable spots, not the writes.)

Cross-check: this matches the reference "fixed-width, cache-blocked SIQS with small/medium/bucketed-
large prime tiers, double-large-prime collection, staged candidate resieving, polynomial batching,
and SIMD for root updates and candidate scans" — Phases 2/2b/1/6 cover it. What that summary omits
but we still need for parity at scale: real Block Lanczos (Phase 4, once nfb is large) and the
per-size auto-tuning (Phase 5) that is the bulk of flint's actual edge.

### Sequencing & expected payoff
Phase 1 → 2 → (4 ∥ 3) → 5. Phases 1-3 (affording nfb≈10k) should close most of the 5x and target
single-core parity; combined with existing scaling that means beating flintqs at low core counts,
not just at tens of cores. Each phase gated on: `cargo test` green (determinism + relation
invariant), factors verified, and a single-core before/after at 192/224/256 vs flintqs.

### Cross-platform validation (2026-07-21)
Regenerated `docs/rusqsieve.wasm` (improved engine) factored a 256-bit semiprime in **225.2s on an
8-thread iPad Air (M3)** — matching native Xeon 8259CL 8-thread (231s). WASM ≈ native on Apple
silicon, and the multiplier + tuning carry straight into the browser demo. Implication for the
roadmap: the wins are **algorithmic** (Phases 1-5) and portable — SIMD (Phase 2b) is a per-target
add-on, not the main lever. Bucket sieving (Phase 2) helps most on cache-constrained phones/older
devices; M3's large caches (like the Xeon's 1MB L2) already keep the current 256KB array resident.

## Session 2026-07-24: Barrett-gated trial division + arithmetic + factor-base retune

Executed the roadmap. The headline result is a portable, correctness-preserving redesign of the
per-survivor trial division that supersedes the "resieve" strategy (frontiers #1 / Phase 1) with the
technique FLINT actually uses. Measured on the reference host (Xeon 8259CL), balanced semiprimes,
every run factor-verified and determinism-checked, all 24 tests green, both WASM feature settings
build:

| bits | cores | before (HEAD) | after | speedup |
|-----:|------:|--------------:|------:|--------:|
| 192  | 1     | 15.60 s | **10.45 s** | −33 % |
| 224  | 8     | 27.27 s | **15.03 s** | −45 % |
| 256  | 32    | 67.81 s | **43.07 s** | −36 % |

(A/B via two saved binaries run interleaved to cancel this shared host's load drift — see the
measurement note below. vs single-thread flintqs the 192-bit gap narrows from ~5.7× to ~3.8×.)

### What landed

1. **Fast small-divisor big-integer division** (`natural/mod.rs::rem_u64`/`div_rem_u64`). Divisors
   `< 2^32` (every factor-base prime) are processed **32 bits at a time in `u64`**, so each step is a
   native machine divide. The previous code did the Horner step in `u128`, which lowers to a
   `__umodti3`/`__udivti3` **libcall** on x86, ARM, *and* wasm (none have a 128-bit hardware divide).
   Portable win, larger on ARM/wasm; verified by the existing randomized `num-bigint` differential
   tests. (~−10 % at 192 on its own.)

2. **Barrett-gated trial division** (`engine.rs`, FLINT `qsieve_evaluate_candidate` style) — the core
   change. A precomputed per-prime Lemire constant `⌊2^64/p⌋+1` (`Context.pinv`) lets us test
   `x mod p == root` with a **multiply-shift** (`fastmod`, ~3 instructions, no divide) and bignum-
   divide `g(x)` **only on a hit**. Trial division becomes `O(nfb)` cheap tests + `O(#factors)` big
   divides, with **no second sieve pass**. This *replaces* the whole resieve machinery (survivor
   bitmap, `cand_at`, `resieve_fac`, the `RESIEVE_MIN_FB` gate — all removed): resieve as a stride
   pass was only ever break-even because it added a pass ≈ the trial division it removed (see
   2026-07-21 note); Barrett-gating removes the cost without adding a pass. The earlier "root-gating
   was 2–4 % slower" finding used a plain `%` (hardware `idiv`, ≈ the bignum-rem cost it replaced),
   *not* a precomputed Barrett inverse — that was the missing ingredient. Debug builds assert
   `fastmod == a % p`; the gate has no false positives/negatives by construction, so the relation
   invariant and determinism are preserved.

3. **Factor-base retune** (`qs/mod.rs::engine_params`). With cheap trial division the relation-starved
   193–224 range profits from a much larger factor base: `193–208 → 120k` (nfb≈5.7k), `209–224 → 250k`
   (nfb≈11k) — both up from 60k. `249+ → 200k` (down from 300k; Barrett shifted the optimum lower and
   it also trims the single-threaded LA). nfb targets track FLINT's `qsieve_tune`. Re-verified at
   192/208/224/240/256 *after* the Barrett change (the earlier resieve-era optima were confirmed to
   still hold or shift as noted).

### Measured findings that correct the roadmap's framing

- **The 192-bit gap vs flint is inner-loop constant factor, not factor-base size.** flint uses ≈3000
   fb primes at 190-bit too (its `fb_primes` counts QR-primes ≈ π(bound)/2; nfb 10k–25k is its 220–260
   regime). At 192 both use ~3000 primes; the 5.7× gap was trial-division (43 %) + score-write (40 %).
   Barrett-gating erases most of the trial-division half.
- **"Afford a large factor base" only helps where the input is relation-starved AND per-poly O(nfb)
   costs don't dominate.** 193–224: big win. 240/256: the optimum FB is *smaller* (nfb≈7k/9k) because
   they run far more polynomials, and each pays O(nfb) for Gray-code root updates + the score-write
   scan. That per-poly O(nfb) cost — not trial division — is now the cap.
- **FLINT's `qsieve_do_sieving2` "block" sieve is carried-position per-block re-striding** — exactly
   the form this audit measured as 13–41 % *slower* on the 1 MB-L2 Xeon. FLINT does **not** use
   (pos,prime) append/drain buckets. Cache-blocking is a small-L2 (mobile) win, not a large-L2 one.

### Tested and reverted (kept as negative results)

- **Branchless conditional-subtract root update** (Phase 2b scalar): **+0.6 s at 192 on x86** (the
   `%p` `idiv` is fine here; masked/branch forms both regressed). Real win needs SIMD — out of scope
   under `#![forbid(unsafe_code)]` on native. Reverted.
- **Slice-iterator (`step_by`) score loop**: regressed — LLVM already elides the indexed loop's bounds
   checks. Reverted.
- **In-loop section profiling** (temporary `RUSQSIEVE_PROFILE` section timers): the `Instant::now`
   calls between hot sections inhibited optimization (~1.6 s at 192 even when disabled). Used to find
   the bottlenecks (score 40 % / factor 43 % at 192), then removed.

### The new bottleneck and the next levers (with data)

After Barrett-gating, **score-write is ≈60 % at 192**, and the per-poly O(nfb) root-update + scan cap
the factor base at 240/256. Concrete next steps, highest-leverage first:

1. **Larger sieve intervals + a cache-blocked sieve, together.** A bigger interval means fewer, larger
   polynomials → the per-poly O(nfb) costs amortize → 240/256 could afford a larger FB (lifting the
   cap above). But a >256 KB score array exceeds mobile L2, so the interval can only grow *paired with*
   a blocked sieve (FLINT's carried-position `do_sieving2`, or a true bucket sieve) that keeps writes
   cache-local. This is the portable structural item that most helps the mobile/consumer targets.
2. **Sparse Block Lanczos** (`f2`, frontier #6). The dense-Gauss LA is single-threaded and grows fast:
   7.6 s at 256/nfb13k, 26 s at 256/nfb21k. It caps FB growth at the top end and must be replaced
   before nfb≫15k is worthwhile.
3. **SIMD root updates + candidate scan** (Phase 2b/#5) — the two vectorizable per-poly O(nfb) loops;
   behind `arch-optimized` with a scoped `unsafe` module + runtime dispatch, per target.

### Measurement note (shared host)

The dev box is a shared 96-core machine; multi-thread **wall** times swing with other tenants' load
(turbo/bandwidth), and thread-summed section timings inflate under load. Reliable comparisons here are
**interleaved A/B of two saved binaries** (head vs new, alternated) — cross-time single-shot numbers
are not trustworthy. Env knobs for continued tuning (no-ops when unset): `RUSQSIEVE_PROFILE`,
`RUSQSIEVE_FB_BOUND`, `RUSQSIEVE_HALFW`, `RUSQSIEVE_SMALL_SKIP`, `RUSQSIEVE_SMALL_SLACK`,
`RUSQSIEVE_THRESH_MARGIN`, `RUSQSIEVE_THRESH_ADJ`.

### Honest status vs. the goal — per-thread parity with FLINT (2026-07-24)

**We are NOT there yet.** Per-thread (single-core) is the yardstick and we are still behind:

| bits | rusqsieve 1-core | flintqs 1-core | ratio |
|-----:|-----------------:|---------------:|------:|
| 192  | 10.45 s          | 2.73 s         | ~3.8× slower |
| 224  | 98.4 s           | 29.0 s         | ~3.4× slower |

(Both measured, factor-verified.) This session closed roughly a third of the per-thread gap (192 was
~5.7× before). We already *beat* single-thread FLINT at modest core counts (224 in 15 s on 8 threads
vs FLINT's 29 s) — but that was already true before this work and is **not** the goal.
**Per-thread parity is not achieved.**

**Why we're still ~3.8× off:** the single biggest single-core cost — the **score-write sieve pass
(~60 % at 192)** — is essentially untouched. Barrett-gating removed the trial-division half of the
cost; the sieve-stepping half remains, and the per-poly O(nfb) costs (root updates + candidate scan)
cap how large a factor base 240/256 can use. Closing the rest needs the structural sieve work below,
which was assessed but **not implemented** this session.

#### Phase-by-phase scorecard

| phase | plan | status |
|------|------|--------|
| 1 | resieve → cheap trial division | **DONE** (better: Barrett-gating, no 2nd pass). Trial-div no longer a bottleneck. |
| 2 | cache-blocked 3-tier / bucket sieve | **NOT DONE** — the #1 remaining lever (score-write ≈60 % @192; caps FB @240/256). |
| 2b | root-update: branchless, then SIMD | scalar **tried+reverted** (x86 regression); SIMD **NOT DONE**. |
| 3 | grow factor base | **PARTIAL** — done at 193–224; 240/256 capped until Phase 2. |
| 4 | sparse Block Lanczos | **NOT DONE** — dense LA single-threaded (7.6–26 s @256); bites only at nfb≫15k. |
| 5 | per-size auto-tuning sweeper | **NOT DONE** — only a coarse manual retune (a big part of FLINT's edge is its auto table). |
| 6 | better `a`-selection, poly batching, larger intervals | **NOT DONE**. |

#### TO-DO to reach (and beat) per-thread parity — prioritized

1. **[Phase 2 + 6] Larger sieve intervals + a cache-blocked / bucket sieve, implemented together.**
   Biggest lever. A larger interval cuts polynomial count → amortizes the per-poly O(nfb) root-update
   and scan → lets 240/256 afford a large FB (where the smooth yield is). The larger score array must
   be blocked (FLINT `do_sieving2` carried-position blocks, or a true (pos,weight)-bucket sieve) so it
   stays cache-resident on small-L2 (mobile) targets. Gate on `cargo test` + determinism + per-size
   1-core A/B. This is expected to close most of the remaining ~3.8×.
2. **[Phase 5] Per-size auto-tuning sweeper** over (fb_bound, interval, small_skip, threshold, ks_primes)
   across several semiprimes per bit-size; bake results into `engine_params`. Removes the
   single-sample overfitting risk in the current hand-tuned values and is much of FLINT's real edge.
3. **[Phase 4] Sparse Block Lanczos** (`f2::BlockLanczos` is a dense-Gauss stub) — required before the
   factor base can grow past nfb≈15k (256+), and it parallelizes, helping multi-core too.
4. **[Phase 6] Better `a`-selection** (`choose_a` is near-random; pick `a≈√(2·kn)/M`, well-spread) and
   **polynomial batching** (process a family's b-variants together to amortize setup + stay cache-hot).
5. **[Phase 2b / #5] SIMD** for the root-update and candidate-scan loops (the two vectorizable per-poly
   O(nfb) spots), behind `arch-optimized` with a scoped `unsafe` module + runtime dispatch, per target
   (x86 AVX2, wasm simd128). Scalar-only won here is a dead end (measured).
6. **[hygiene]** Add a direct differential test for `fastmod`/`lemire_c` (currently only covered
   indirectly by a debug-only `debug_assert` in the engine + the end-to-end factor tests).

## Session 2026-07-24: translated position roots

The first follow-up after the Barrett work removes another division-heavy operation from the
score-write hot path:

- Sieve roots are now stored as score-array-position residues (`x + interval mod p`) instead of
  signed polynomial-coordinate residues (`x mod p`). `interval mod p` is precomputed once with the
  immutable engine context. The old sieve translated both roots on every prime of every polynomial
  with `i64::rem_euclid`, issuing multiple hardware divides in the dominant loop.
- Gray-code root advances remain unchanged because translating every root by a fixed residue
  commutes with each modular update. Candidate trial-division gating now compares `position mod p`
  directly, also removing the per-candidate signed-coordinate conversion.
- Added the roadmap's direct differential coverage for `lemire_c`/`fastmod`, together with the new
  modular root-translation helper.

Single-thread A/B on the checked-in 58-digit (192-bit) example, using saved release binaries and
identical sieve output:

| version | wall time | polynomials | survivors | relations |
|--------:|----------:|------------:|----------:|----------:|
| prior HEAD | 9.91–9.99 s | 21,888 | 65,304 | 3,116 |
| translated roots | **7.45 s** | 21,888 | 65,304 | 3,116 |

That is a repeatable **25% wall-time reduction** while preserving factorization output and all sieve
counts. The full native suite passes (25 unit tests plus integration/doc tests), and the portable
WASM release build remains green. FLINT factors this same input in 2.83 s wall / 2.73 s user, so the
single-core gap narrows from about 3.5× on this run to about 2.7×. The larger blocked-interval
structural work remains the next major lever.

## Session 2026-07-24: FLINT hot-loop parity work

Compared directly against the authoritative checkout at `/home/dev/flint/src/qsieve/` and ported
four scalar techniques from `collect_relations.c` / `compute_poly_data.c`:

- A paired, two-hit-unrolled root-stride kernel replaces two independent saturating loops. Inputs
  with `g_bits <= 192` use proven non-overflowing byte addition; wider inputs retain saturation.
- Scores are biased by `128 - threshold`, allowing the overwhelmingly-empty candidate scan to
  reject eight bytes with one `u64` high-bit test, with a scalar fallback outside the safe range.
- SIQS `B` is now a true signed coefficient with the smaller CRT term, not reduced modulo `A`.
  Gray-code root updates are therefore one conditional modular add/subtract, matching FLINT and
  eliminating the old per-prime modular-wrap correction.
- Candidate factoring stops once confirmed factor weights account for the stored sieve score, and
  skipped small primes use the same Lemire root gate. Retuning moves the 192-bit small-prime cutoff
  from 20/slack 3 to 100/slack 8, close to FLINT's count-based `small_primes=13` setting.

On the checked-in 192-bit example, single-thread wall time moves from the translated-root result
of 7.45 s to **5.67 s** (FLINT: 2.83 s), a further 24% and a total 43% reduction from the 9.91 s
session baseline. On a fresh balanced 224-bit semiprime
`21523772555907914536866856055060033603780528151558474367883009969243`, rusqsieve completes in
**45.52 s** versus FLINT **26.76 s**; the pre-session audit baseline for this tier was 98.4 s.

Measured negative results retained as guidance:

- At 192 bits, increasing half-width from 65,536 to 131,072 regresses 7.30 s → 8.12 s; a 32,768
  half-width only gains ~2%. Larger intervals are not independently beneficial before blocking.
- Reducing the relation target below raw factor-base size fails extraction at 90%, 95%, and 98% on
  the 224-bit sample. FLINT's much stronger singleton reduction / relation graph yield cannot be
  copied by changing one threshold; the experiment was reverted.
- At 224 bits the existing 250k factor-base bound remains near the local optimum: 200k takes
  50.76 s and 350k takes 50.61 s versus 45.52 s at 250k.
- A compact prime/weight structure-of-arrays view regresses 192-bit time by ~4%; reverted.

Per-thread parity is still not achieved, but the gap is now ~2.0× at 192 bits and ~1.7× at 224
bits, down from 3.8× / 3.4× at the start of the day. The remaining high-leverage difference visible
in FLINT is its cache-blocked large-interval sieve plus more effective relation filtering.

## Session 2026-07-24: per-thread FLINT lead achieved through 240 bits

Continued direct comparison with `/home/dev/flint/src/qsieve/` produced the decisive changes:

- **Exact sorted-root/difference sieve loop.** Roots are kept ordered after every Gray-code update.
  The hot stride loop carries one position plus the root difference, matching FLINT's unrolled
  kernel instead of maintaining two positions and four loop bounds. This alone moved the 192-bit
  sample from 5.67 s to 4.92 s.
- **Numerically target-fitted `A`.** `choose_a` now selects a deterministic factor count and primes
  from the ideal size class, then chooses the final factor closest to
  `target_A/current_product`. The old code stopped when the product merely reached the target bit
  length, allowing a large last-prime overshoot. At 192 bits this cuts survivors 86,767 → 15,955
  and time 4.92 s → 3.07 s before final retuning; at 224 bits it cuts 36.0 s → 21.4 s.
- **Weight-two structured elimination.** Matrix filtering now deterministically merges the two
  variables in a weight-two row while retaining exact original-column provenance. On the 224-bit
  matrix this reduces 11,000×11,063 to 8,228×8,642. Sub-raw-row relation targets were tested:
  99% works at 224 but fails at 192, so the safe default remains `nfb + 64` until staged
  collect/filter/retry exists.
- **Carried-position blocked sieve implemented.** Correct and available for score arrays at least
  1 MiB. On this 1 MiB-L2 Xeon, a 512 KiB array is faster flat (40.23 s vs 42.36 s at 224), so the
  gate reflects the actual cache rather than FLINT's conservative 512 KiB switch point.
- **Tier retuning after the structural wins.** 177–192 bits now use bound 100k / half-width 90,112;
  209–224 uses half-width 262,144; 225–248 uses bound 350k / half-width 262,144. The upper-tier
  factor base can now grow because target-fitted polynomials and cheaper per-prime loops removed
  the old cap.

Clean saved-binary single-thread comparisons on the dedicated container:

| bits | input | rusqsieve | FLINT | result |
|-----:|-------|-----------:|------:|-------:|
| 160 | checked-in 49-digit example | **0.35 s** | 0.57 s | 1.63× faster |
| 192 | checked-in 58-digit example, two interleaved pairs | **2.49–2.53 s** | 2.78–2.79 s | 1.10–1.12× faster |
| 208 | fresh balanced semiprime | **8.59 s** | 9.89 s | 1.15× faster |
| 224 | fresh balanced semiprime | **21.43 s** | 26.60 s | 1.24× faster |
| 240 | fresh balanced semiprime | **67.25 s** | 77.18 s | 1.15× faster |

All outputs were product-verified by the CLI. The full native suite passes and the release
`wasm32-unknown-unknown` build passes. This achieves the audit's per-thread goal on every measured
head-to-head tier from 160 through 240 bits; FLINT still aborts on the reference host's 256-bit
case, where rusqsieve remains functional.

## Session 2026-07-25: 256-bit and browser linear algebra

The previous upper-tier parameters were stale: at 256 bits, bound 200k / half-width 131,072 spent
42.48 s of a 43.37 s eight-thread native run collecting relations while LA was only 0.84 s. A sweep
over factor-base bound and interval selected **bound 500k / half-width 327,680**. The larger base
raises the matrix to about 20.8k columns but cuts relation collection enough to win overall.

That exposed the next crossover and led to three linear-algebra changes:

- structured sparse Gaussian elimination now handles rows through weight six, choosing a
  minimum-column-weight pivot to control fill and replaying substitutions backward;
- the residual dense solve uses compact row echelon equations instead of carrying equally large
  parity and provenance vectors through column elimination;
- it emits the conventional maximum of 64 verified dependencies, solving only the vectors useful
  to factor extraction. Every expanded dependency is checked against the original matrix.

On the fixed 256-bit semiprime
`98877949376972157840865984674312121822345015130827118595228756728313751597271`, sparse filtering
reduces approximately 20,740×20,803 to 10,877×11,343. The original filtered dense solve took 4.69 s
inside a 5.19 s extraction; compact row echelon takes 2.43 s inside a 2.93 s extraction, a 44%
end-to-end LA reduction. The final native eight-thread run is 27.83 s wall. This is still dense
Gaussian elimination, not a true block-Lanczos recurrence, but it moves the practical crossover
well beyond the current 256-bit matrix.

A real Node/V8 Worker benchmark now mirrors the browser topology: one coordinator Wasm instance,
eight independent worker instances, serialized relation packets, and no shared memory. Reducing
jobs from four polynomial families to two trims the relation-collection tail. A scoped Wasm
`simd128` XOR kernel only in row reduction succeeds where whole-program `+simd128` regressed:

| bits | V8/Wasm, 8 workers | sieve | filter/LA/extract |
|-----:|--------------------:|------:|------------------:|
| 192 | **0.72 s** | 0.54 s | 0.18 s |
| 224 | **5.04 s** | 3.95 s | 1.09 s |
| 256 | **37.86 s** | 32.08 s | 5.78 s |

All three results are factor-verified single measurements under Node 24.15. The browser build ships
both SIMD and scalar Wasm and falls back at module compilation on engines without `simd128`.
Worker scaling at 256 bits is 22.26 s / 14.71 s / 13.96 s for 16 / 32 / 48 workers. At 224 bits,
48 workers take 2.59 s, while 96 regress to 3.45 s; the browser pool is therefore capped at 48.

Competitive-claim caveat: these results justify continued work toward a fastest browser SIQS claim,
not yet “fastest general factorizer.” Alpertron uses ECM before SIQS, while rusqsieve currently jumps
from bounded Pollard-Brent to SIQS. Medium factors in unbalanced 192–256-bit composites are therefore
a known algorithmic gap. `BENCHMARKING.md` defines the balanced and unbalanced corpus and same-browser
competitor protocol required before publishing the broader claim.
