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
|  128 | >60 s            | 0.51 s          | (refuses)          | —          |
|  160 | >60 s            | 0.31 s          | 1.22 s             | **3.9× faster** |
|  192 | >60 s            | 1.22 s          | 5.13 s             | **4.2× faster** |
|  224 | >60 s            | 10.7 s          | 36.9 s             | **3.4× faster** |
|  256 | >60 s            | 105.6 s         | 231.4 s            | **2.2× faster** |

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

### Further headroom (already beats flintqs, but per-core efficiency has room)
Single-threaded, rusqsieve is ~4–5× slower per core than flintqs; the win comes from using the
box's cores well. Two orthogonal, higher-effort optimizations would widen the lead, especially
at 224–256 bit (77 digits is the QS practical edge — the value sieved is ≈2^144, so smoothness
is rare):
- **Double-large-prime variation** (flintqs uses it): capture cofactors that are a product of
  two large primes and combine them via cycle-finding (union-find over the large-prime graph).
  Biggest yield lever at large sizes; higher risk (relation-combination + matrix correctness).
- **Real sparse linear algebra**: `f2::BlockLanczos` is currently a stub delegating to dense
  Gauss–Jordan (`dense_dependencies`, ~O(n³/64)); a genuine block-Lanczos recurrence (SPEC §15)
  would remove the dense-elimination cost that grows with factor-base size.
- **Done:** the low-level portable kernel `qs::sieve_job` — previously trial-division of every
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
