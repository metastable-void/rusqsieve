# rusqsieve v0.2.0 — Remediation Brief for a Coding Agent

You are working on `metastable-void/rusqsieve`, a self-initializing quadratic sieve (SIQS)
integer factorization crate for native Rust and WebAssembly. Published version is **0.2.0**.
Target workload: balanced RSA-style semiprimes of 192–256 bits, including in-browser
execution across Web Workers **without** `SharedArrayBuffer`.

This brief lists verified defects and the order in which to fix them. Every file:line
reference below was verified against the v0.2.0 source. If a reference does not match what
you find, **stop and report the discrepancy** rather than guessing at the intended target.

---

## 0. Ground rules

### 0.1 Semver discipline is mandatory

The crate is at `0.2.0`, so under Cargo's semver rules the **minor** field is the breaking
field. Classify every change before you make it:

| Change class | Allowed in | Notes |
|---|---|---|
| Internal-only (private items, function bodies, perf) | `0.2.1` patch | No change to any `pub` signature, no change to observable output other than fixing a documented bug |
| Fixing a wrong answer / a failure that should have succeeded | `0.2.1` patch | Bug fix, even though behaviour changes |
| Adding a `pub` item, adding an enum variant to a `#[non_exhaustive]` enum, adding a Cargo feature | `0.3.0` minor | Additive but still minor for a `0.x` crate; prefer to batch |
| Removing or renaming any `pub` item, removing an enum variant, removing a Cargo feature, changing a `pub` signature or a type's default generic parameter | `0.3.0` minor (breaking) | **Never in a patch release** |

Deliver the work as **two releases**, in this order:

- **`0.2.1`** — Phase 1 and Phase 2 only. Zero public API changes. This gets the correctness
  fix to users immediately.
- **`0.3.0`** — Phase 3 onward. All API removals and signature changes batched into one
  breaking release.

For anything you would remove in `0.3.0`, add `#[deprecated(since = "0.2.1", note = "...")]`
in the `0.2.1` release first where that is possible without changing a signature. Do not
deprecate items you are keeping.

Maintain `CHANGELOG.md`. Every entry states the version, the semver class, and for behaviour
changes what input now behaves differently.

### 0.2 Verification rules

- **"`cargo test` passes" is not a completion criterion.** The defects below all shipped in a
  tree where the full suite was green. Each task states its own acceptance criterion; meet
  that.
- Write the failing regression test **before** the fix, and confirm it fails for the stated
  reason.
- For every performance change, measure before and after on the same host, in the same
  session, with runs interleaved (A/B/A/B), and report both numbers. Do not claim an
  improvement you did not measure. If run-to-run variance exceeds the effect, say so and
  keep the change only if it is also a clarity or memory win.
- Never weaken a test to make it pass. Never delete a test without saying why in the
  changelog.

### 0.3 Baseline measurements

Reference figures from a 4-core Xeon @2.80GHz (L1d 32 KiB, L2 1 MiB/core, L3 33 MiB),
release build, 4 threads, `RUSQSIEVE_PROFILE=1`. Sieve wall times on this class of host vary
run to run by tens of percent — treat them as a starting point, re-measure on your own host,
and rely on interleaved A/B comparisons rather than absolute numbers.

```
192-bit: sieve+collect ≈ 0.83 s   nfb 4758    full 3844   partials 10961   cycles 979
224-bit: sieve+collect ≈ 6.17 s   nfb 10999   full 7682   partials 45134   cycles 3381
256-bit: sieve+collect ≈ 53.5 s   nfb 20844   full 11000  partials 261974  cycles 9803
256-bit: filter 20845x20908 -> 11163x11637, f2_filter 0.494 s, f2_dense 3.740 s,
         extract(LA) 5.182 s
```

Cost split of `sieve+collect` (CPU-seconds summed over threads): score writes 55–68%,
per-family setup 11–19%, trial division 7–12%, candidate scan 4–5%.

Fitting `T = S + P/w` to the published browser scaling points (8/16/32/48 workers →
37.86/22.26/14.71/13.96 s at 256-bit) gives **S ≈ 8.06 s serial, P ≈ 236 s·worker**. The
serial term is the linear algebra. Use this fit to sanity-check any claimed parallel-scaling
improvement.

### 0.4 Scope

Do not add new algorithms beyond those named here. Do not restructure modules beyond the
splits explicitly requested. Do not touch the benchmark corpus or `BENCHMARKING.md`
numbers except to append new measurements clearly labelled with your host. The factorization
corpus supplied with this brief is new to the repository; commit it as-is and do not prune
entries from it.

---

## Phase 1 — Correctness (`0.2.1`, blocking everything else)

### 1.1 The 65–85 bit dead zone

**Symptom.** Balanced semiprimes from 65 bits to roughly 85 bits fail deterministically:

```
$ echo 18446744400127067027 | qs-factor --progress never
qs-factor: factor engine failed: no nontrivial factor found     [exit 2]
# 18446744400127067027 = 4294967311 × 4294967357
```

64-bit and smaller inputs are fine (they take the `u64` Pollard-Brent fast path). 65–81 bits
always fail; 82–85 bits fail depending on the Knuth-Schroeppel multiplier; 86 bits and up are
fine. Reproduces at any thread count. The failure is not a wrong answer and not a hang — the
engine reports no factor found.

**33 of the 200 entries in the corpus used during the review fail** — that is the measured
reach of the bug, not a property of the corpus shipped with this brief. Most of those 33 are
not themselves in the dead zone: a larger `n` splits successfully once, and the *cofactor*
lands in the zone. Example, 151-bit input:

```
PROFILE nfb=2034 interval=65536 target=2098 k=5
PROFILE sieve+collect=0.043s polys=88 families=11 survivors=7655 relations=2098   <- ok
PROFILE nfb=226  interval=32768 target=290  k=2
PROFILE sieve+collect=1.660s polys=0 families=100000 survivors=0 relations=0      <- dead zone
qs-factor: factor engine failed: no nontrivial factor found
```

`rusqsieve-factorization-corpus.txt`, supplied alongside this brief, is a rebuilt and
independently verified replacement for that review corpus: **309 entries spanning 7–256
bits**, each line an `N` followed by its complete prime factorization (arity varies — see
the file's own header). **117 of them sit in the 65–85 bit dead-zone band, 87 of those
semiprimes that fail deterministically on v0.2.0.** It is *not* currently in the repository;
add it under `tests/data/`.

This affects the library API, not only the CLI: `rusqsieve::native::factor_impl` →
`engine::factor_node`.

**Root cause.** `choose_a` (`src/engine.rs:910-958`) can never build a candidate pool.
`src/engine.rs:914` draws candidates from factor-base primes `> 1000` — bit length ≥ 10 —
and then `src/engine.rs:923-932` keeps only primes whose bit length is within ±1 of
`ideal_bits`:

```rust
let target_bits = ctx.target_a.bit_len();                    // :921  == 18 for the min repro
let factor_count = target_bits.div_ceil(14).clamp(3, 10);    // :922  == 3
let ideal_bits = target_bits.div_ceil(factor_count);         // :923  == 6
// ... bits.abs_diff(ideal_bits) <= 1 ...                    // :930  needs 5..=7 bits
if pool.len() < factor_count * 2 { return None; }            // :933  pool is empty
```

`bits >= 10` and `bits <= 7` are unsatisfiable, so the pool is empty for every family,
`choose_a` returns `None` 100 000 times, and `polys` stays 0. Confirmed: `polys=0` persists
for every value of `RUSQSIEVE_HALFW` (32768 down to 128) and `RUSQSIEVE_FB_BOUND` (3000 down
to 200).

**Task.**

1. Add a regression test first. Table-drive it over at least one balanced semiprime at each
   of 65, 70, 75, 80, 85, and 90 bits, asserting the returned factors multiply back to `n`
   and each is prime. Include `18446744400127067027`. Confirm it fails before the fix.
2. Make the pool constraint satisfiable. Derive the lower prime bound from `ideal_bits`
   instead of hardcoding `1000` at `src/engine.rs:914`, and/or widen the `abs_diff <= 1`
   window when the pool comes up short. Whatever you choose, add an assertion or a
   `debug_assert` that the constructed constraint set is non-empty, so this class of bug
   cannot recur silently.
3. Add a corpus test driven by the supplied corpus file. Commit
   `rusqsieve-factorization-corpus.txt` into the repository at
   `tests/data/rusqsieve-factorization-corpus.txt` — its own header documents its provenance
   and the verification every entry passed. Each line is `N f1 f2 ... fk`; **the arity
   varies, so the parser must read the whole line and must not assume exactly two factors.**
   The test must factor every entry and verify each result (product equals `n`, every factor
   prime), with zero failures.

**Acceptance:** the size-sweep test and the corpus test both pass; `polys > 0` for
every size in the sweep under `RUSQSIEVE_PROFILE=1`.

### 1.2 A `choose_a` famine must not be reported as "no factor"

`choose_a` returning `None` for all 100 000 families (`src/engine.rs:662` sets
`max_families = 100_000`; `src/factor.rs:453` has a second, uncoordinated 100 000 cap) is a
parameter-selection failure, not a search failure. Today it burns ~1.6 s and reports
`NoFactor`, which is what made 1.1 hard to diagnose.

Surface it distinctly, and bail out early rather than retrying a deterministic failure
100 000 times. Because adding a public `FactorError` variant is a minor-version change, in
`0.2.1` implement this as an internal engine error plus a distinguishable log line, and add
the public variant in `0.3.0` (`FactorError` is `#[non_exhaustive]`, so that addition is
cheap).

**Acceptance:** with the 1.1 fix reverted locally, the min repro fails in well under 100 ms
with a message naming polynomial-coefficient selection.

### 1.3 No big-integer rho above 2^64

`factor_node` (`src/engine.rs:498-543`) is the complete dispatch ladder:

1. `:520-527` — `n < 2^64` → `smallfactor::factor_u64` (deterministic Miller-Rabin +
   Pollard-Brent, `src/smallfactor.rs:113,165`)
2. `:526` — `is_probable_prime`
3. `:530` — `perfect_power`
4. `:540` — `find_factor` → SIQS, unconditionally

There is **no trial division and no Pollard-rho for `n >= 2^64`** on this path. A generic
`pollard_rho` over `Natural<P>` already exists at `src/factor.rs:397` (reached from
`src/factor.rs:368-379`) but the engine never calls it. `src/engine.rs:1504` has
`pollard_u64`, used only for splitting double-large-prime cofactors that already fit in
`u64` (`src/engine.rs:1494`).

Add a bounded big-integer Pollard-Brent stage at `src/engine.rs:540`, before falling through
to SIQS. For 65–100 bit inputs rho beats SIQS outright, so this is a performance win as well
as defence in depth — it would have masked 1.1 entirely. Bound the iteration count and make
it cancellable.

**Acceptance:** inputs in 65–100 bits are factored by the rho stage (assert via profile
output that SIQS is not entered); no regression above 128 bits.

### 1.4 Cancellation and progress gaps that make 1.1-class bugs invisible

- `src/engine.rs:520-525` — for cofactors ≤ 64 bits, `smallfactor::factor_u64` runs with no
  cancellation check, and `src/smallfactor.rs:118`'s `pollard_brent` loop is unbounded and
  uninterruptible. Add a cancellation poll.
- `src/engine.rs:630` — `rx.lock().unwrap()` means one worker panic poisons the job mutex, so
  every other worker then panics with `PoisonError` instead of the real cause, and
  `src/engine.rs:705`'s `let _ = h.join()` discards the original payload. Use
  `lock().unwrap_or_else(|e| e.into_inner())` and propagate the first panic message.

### 1.5 Baillie-PSW as the final primality decision

This is a requirement from the project owner: add Baillie-PSW on top of a reasonable number of
Miller-Rabin rounds, in two places — the library's final primality decision, and the web demo's
RSA number generation — provided it does not slow things down materially.

**Why it matters.** `src/primality.rs`'s 16 default rounds draw their bases via
`SMALL[round % 32]`, i.e. bases 2..53. That is **deterministic, not probabilistic**, so a strong
pseudoprime to all of bases 2..53 is a guaranteed false positive — and such numbers are
constructible (Arnault's constructions; Albrecht, Massimo, Paterson and Somorovsky, "Prime and
Prejudice", 2018, which broke fixed-base implementations in named libraries). By contrast **no
Baillie-PSW pseudoprime is known to exist**; constructing one is an open problem.

**Scope, with one important exclusion.** Below 2^64 the 7-base Jaeschke/Sinclair witness set at
`src/smallfactor.rs:73` is **proven exact**, so BPSW adds nothing there. Branch on size and apply
BPSW for `n >= 2^64` only.

**Implementation notes.**

- BPSW is: trial division by small primes, then a base-2 Miller-Rabin, then a strong Lucas test
  with Selfridge Method A parameter selection.
- The strong Lucas test needs a **Jacobi symbol over `Natural<P>`**. The existing `jacobi_u64`
  (`src/natural/mod.rs:1152-1236`) is u64-only, so this is new code.
- Selfridge's `D` search does not terminate for perfect squares, so a square test must run first.
  **`is_square` (`src/natural/mod.rs:491`) is currently dead code and Phase 3.3 lists it for
  deletion — this gives it a user. Keep it.**
- Keep the Miller-Rabin rounds as an additional layer, but they are no longer what carries the
  guarantee.

**Cost.** BPSW is roughly one Miller-Rabin round plus one strong Lucas, and a strong Lucas costs
about 2–3× a Miller-Rabin round — so about 15–20% added to the primality path at 16 rounds.
Primality is **not** the sieve hot path: it runs at `factor_node` entry (`src/engine.rs:526`) and
on recovered factors, never per survivor. Measure the primality path before and after and report
both numbers. If it turns out to exceed a few percent of total runtime, reduce the Miller-Rabin
round count rather than dropping BPSW — BPSW is what carries the strength.

**Web demo.** Apply the same structure in `web/numtheory.js`'s RSA number generation. Candidates
are rejected quickly by trial division and the base-2 Miller-Rabin, so the Lucas test runs roughly
once per accepted prime; JS `BigInt` is adequate.

**Acceptance:** the extended primality tests from §4.3 (Carmichael numbers, strong pseudoprimes to
base 2) all reject; the primality-path timing is reported before and after. **And critically —
assert that the Lucas step is actually reached and exercised** (a counter, or a test hook).
Because no BPSW pseudoprime is known, there is no input that distinguishes a correct Lucas
implementation from one that silently never runs, so every test would also pass if the Lucas step
were dead. Given this crate's existing stubs, prove it executes.

---

## Phase 2 — Performance, internal only (`0.2.1`)

Everything in this phase is behind private items or function bodies. No public signature
changes. Do them in this order — the ordering matters, because items 2.1 and 2.2 are what
make 2.3 and 2.4 affordable.

### 2.1 Resieving (unblocks the rest)

`src/engine.rs:1394-1418` trial-divides each survivor by walking the **entire** factor base,
computing `fastmod(posu, p, pinv[idx])` (Lemire multiply-shift, `src/engine.rs:1663-1666`)
and comparing against both roots, with an early exit at `src/engine.rs:1395` once
`confirmed_score >= score_target`. There is no resieving and no bucket-recorded hits —
`grep -rniE "resiev|bucket" src/` returns a single comment at `src/engine.rs:1268`.

Measured marginal cost: **≈15 µs ≈ 43 000 cycles per survivor** at nfb=4758, about 9 cycles
per factor-base prime. The scan really is O(nfb): the early exit only fires once the
*largest* factor is found, and the base is ascending, so a typical survivor walks most of it.

Implement msieve-style resieving so the per-survivor cost becomes O(number of factors),
roughly 1 000 cycles — a ~40× reduction on that step. Directly this is only 7–12% of
runtime; its real value is that survivor cost currently caps 2.3 and 2.4.

Also in the same region:
- `src/engine.rs:1333` heap-allocates `powers` per survivor (≈340 000 allocations at
  256-bit) and `record` does a linear `find` over it.
- `src/engine.rs:1410-1412` divides twice per successful hit: `q.rem_u64(p)` followed by
  `q.div_rem_u64(p)`.

**Acceptance:** per-survivor trial-division cost drops by ≥10× measured; relation counts and
factorization results unchanged on the corpus.

### 2.2 Cheap arithmetic wins

- **`src/natural/mod.rs:568-570`** — `mod_u64` calls `div_rem_u64`, which materializes a full
  128-byte zeroed quotient (`src/natural/mod.rs:401`) that the caller discards. `rem_u64`
  exists directly above for exactly this purpose and its doc comment
  (`src/natural/mod.rs:367-369`) says so. **One-line fix.** Called from
  `src/engine.rs:803,810,821` roughly `(2 + nvar) · nfb` ≈ 167 000 times per family, times
  3 099 families at 256-bit.
- **`src/natural/mod.rs:548-552`** — `mul_mod` performs **three** Knuth divisions per
  multiply (`self.div_rem(m)`, `rhs.div_rem(m)`, then `rem_natural(m)`), although every
  caller in `extract` passes already-reduced operands. Drop the two redundant reductions;
  document the precondition and `debug_assert` it.
- **`src/natural/mod.rs:621,630,640,680`** — `knuth_divmod` allocates four `Vec`s per call.
  This is the workhorse of every multi-limb division, so those allocations should become
  stack buffers — but **do not write `[u64; 2*P+1]`**: an array length computed from a const
  generic parameter needs `generic_const_exprs`, which is unstable, and the crate is
  stable-only. Two stable ways to get the same result:
  **(a)** mirror what `WideNatural` (`src/natural/mod.rs:1032-1073`) already does and carry
  several `[u64; P]` arrays instead of one `2*P`-sized array — a standalone `P` is a legal
  array length; or
  **(b)** since `PARTS` is a crate constant, define a literal `const MAX_LIMBS: usize = 16;`
  and declare `[u64; 2 * MAX_LIMBS + 1]` — array lengths built from literal consts are fine
  on stable — then slice it down with `&mut buf[..2 * p + 1]` for the `p` actually in play.
  Option (b) costs a little stack when `P` is smaller than `MAX_LIMBS`, but it removes the
  heap traffic and keeps one buffer instead of several. This composes with Phase 3.4, which
  makes `PARTS` a fixed constant.
- **`src/natural/mod.rs:1063`** — `WideNatural::rem_natural` allocates a `Vec::with_capacity(2*P)`
  purely to concatenate `low`/`high` into a contiguous slice.
- **`src/natural/mod.rs:277-278,292`** — `widening_mul` zero-initializes `[u64; 16]` twice and
  returns 256 bytes by value, then `overflowing_narrow` copies 128 bytes out, even though
  `sig_len` (`:281-282`) correctly limits the multiply itself to a few limbs. Three such
  `checked_mul`s run per survivor (`src/engine.rs:1315,1322,1324`). On wasm, with no wide
  registers, this is pure `memory.copy`.
- **`src/engine.rs:166`** — `EngineJobResult::to_bytes` starts from `Vec::new()` and grows to
  ~4 KB per job by repeated realloc. The final size is computable up front; use
  `with_capacity`.
- **`src/factor.rs:388`** — `primes_to(limit)` rebuilds the prime list by trial division on
  every `factor_node` call, and `factor_node` recurses (`:385-386`). With the default
  `trial_division_limit = 10_000` (`src/factor.rs:97`) that is ~1 229 primes recomputed per
  node. `smallfactor::small_primes()` (`src/smallfactor.rs:15`) already has a cached
  Eratosthenes sieve; use it.
- **`src/qs/mod.rs:181-193`** — `prime_u32` builds the factor base by testing every odd number
  up to 500 000 with O(√p) trial division instead of a sieve of Eratosthenes.
  `fb_build ≈ 0.127 s` at 256-bit looks negligible, but **it is paid per Web Worker at
  browser startup**, against a 0.72 s total for the 192-bit case. Replace with a segmented
  Eratosthenes sieve.

### 2.3 Retune the sieve threshold

`small_slack = 8` (`src/engine.rs:970-973`) and `thresh_margin = 4` (`src/engine.rs:977-980`)
are hardcoded magic numbers, and the shipped default sits on the wrong side of the optimum.
Sweep at 224-bit, 4 threads:

| `THRESH_ADJ` | polys | survivors | wall |
|---:|---:|---:|---:|
| +8 | 39 168 | 22 088 | 8.64 s |
| +4 | 32 192 | 41 344 | 7.05 s |
| **0 (shipped)** | 25 984 | 72 540 | **6.00 s** |
| −4 | 21 824 | 122 740 | 5.30 s |
| **−8** | 19 200 | 206 265 | **5.16 s** |
| −12 | 17 664 | 354 551 | 5.77 s |

**≈14% available at 224-bit**, same direction at 192-bit. Note this interacts with 2.1: a
cheaper survivor makes the optimum deeper still, so retune *after* resieving lands.

A contributing cause is log resolution. The sieve log weight is `(32 - p.leading_zeros())` —
plain bit length, 1-bit resolution (`src/engine.rs:1061`, `src/engine.rs:1153`). That
over-estimates by up to 1 bit per prime, with ~15 primes per smooth value, which is precisely
why an 8-bit slack fudge is needed. `FactorBaseEntry::log_prime` already stores
`round(ln p · 8)` (`src/qs/mod.rs:154`) but **the engine never reads it** — only the
reference path does (`src/qs/mod.rs:562`). Switch the engine to a scaled log (2× or 4×) and
re-derive the slack.

Also: `small_slack` is a fixed constant where the true expectation is
`Σ log(p)/(p−1)` over the skipped set — ≈4.5 bits against the shipped 8. Derive it.

And note the overflow hazard at `src/engine.rs:1000-1002`: the non-saturating fast path's
soundness argument assumes every scored prime is at least 23, but `RUSQSIEVE_SMALL_SKIP` can
be set below that, silently turning the argument into wrapping overflow and false negatives.
Phase 3.1 removes that env var; until then, guard the invariant.

**Acceptance:** ≥10% wall-clock improvement at 224-bit and 256-bit, measured A/B interleaved;
the derived slack is computed from the skipped set rather than hardcoded.

### 2.4 Decouple the large-prime bound from the threshold slack

**`lp_allowance` is one number doing two incompatible jobs**: the sieve threshold slack
(`src/engine.rs:1260-1262`) and log2 of the large-prime acceptance bound
(`src/engine.rs:298`). Consequence: `single_limit = 1 << lp_allowance` makes the LP bound
**34 360× the factor-base bound** at 256-bit, against msieve's `large_prime_mult` of ~100–200.
Across the tier table (`src/qs/mod.rs:330-341`) the LP/FB ratio is incoherent — 22, 44, 105,
70, 42, 559, 268, 3 068, 34 360.

The cost: at 256-bit the crate retains **261 974 partials to produce 9 803 cycles — a 3.7%
yield**, roughly 50–80 MB of live relations per coordinator, nearly all never used. Partials
with a 30–34 bit large prime essentially never pair.

Split the two into separate parameters and set `large_prime_mult` in the 100–250 range.
Measured at 256-bit, 4 threads, changing only the LP bound:

| LP bound | wall |
|---|---:|
| shipped, 2^34 (34 360×) | 53.5 s |
| 256 × FB (2^27) | **50.9 s** |
| 64 × FB (2^25) | 52.8 s |

≈5% plus the memory.

While you are in this code: `Forest::reroot`/`Forest::path` (`src/engine.rs:1560-1582`)
**clone whole `Relation` values** — a 16-limb `Natural` plus a `Vec` — along every tree path
on every ingest. It is a re-rooting forest, not path-compressed union-find. Store indices and
clone nothing.

**Acceptance:** cycle yield at 256-bit improves substantially over 3.7%; peak partial-relation
memory drops by an order of magnitude; wall time does not regress.

### 2.5 Replace `pollard_u64`

`src/engine.rs:1504-1531` uses Floyd rather than Brent, a `u128 %` on every step, and a `gcd`
on **every iteration** with no batching. Measured ≈200 000 cycles per split. At a deep
threshold (`THRESH_ADJ=-20`, 192-bit) cofactor splitting accounted for ~14.5 s of a 19.2 s
run.

Replace with Brent's variant plus Montgomery multiplication and a batched GCD every ~128
steps, or with SQUFOF. Expect 20–40×. This is a prerequisite for 2.6.

**Acceptance:** ≥20× measured improvement in cycles per cofactor split; identical split
results on a fixed cofactor corpus.

### 2.6 Only now: enable double large primes

`LargePrime::Two`, `classify_cofactor`, and cycle combination all exist
(`src/engine.rs:77-91`, `:1487-1501`, `:1455-1483`) and are **unreachable in every shipped
configuration**. The gate at `src/engine.rs:297-302`:

```rust
let double_enabled = lp_allowance as u32 >= 2 * bound_bits + 2;
```

is satisfied by **no tier** in `src/qs/mod.rs:330-341`. Verified by patching
`double_enabled = true` and rebuilding: poly, survivor, and relation counts came out
**byte-identical** at 192-bit and 224-bit. A second, independent blocker: the threshold slack
is `lp_allowance + small_slack − thresh_margin` = 26 bits at 192-bit, while a genuine double
needs ≥ 2·log2(100 000) ≈ 34 bits of cofactor, so doubles cannot survive the threshold even
with the gate forced open.

**Do not enable DLP before 2.1 and 2.5 land.** Measured with a properly parameterized DLP
(msieve-style `LP = mult × FB_bound`, always on) at 192-bit, 1 thread, on the unfixed tree:

| config | polys | survivors | time |
|---|---:|---:|---:|
| shipped, adj=0 | 8 064 | 22 917 | 2.70 s |
| shipped, adj=−8 | 5 696 | 69 430 | 2.54 s |
| DLP on, adj=−10 | 5 152 | 86 975 | 2.78 s |
| DLP on, adj=−20 | 3 520 | 301 771 | **19.2 s** |

DLP cuts polynomials ~12% at matched thresholds but **loses on wall time** while survivors
and cofactor splits are expensive. After 2.1 and 2.5, re-measure and keep it only if it wins.

**Acceptance:** either a measured net win with the gate correctly parameterized, or a
documented decision to leave DLP disabled — in which case fix the docs per Phase 4 rather
than leaving dead code claiming a feature.

### 2.7 Per-family setup (11–19% of runtime)

`src/engine.rs:798-825` runs, for every factor-base prime: two `Natural::mod_u64`, one
`inv_u32`, then `nvar` further `mod_u64` plus two `mulmod_u32`. `inv_u32`
(`src/engine.rs:1624-1636`) is a textbook extended Euclid with a hardware 64-bit divide per
iteration — about 15 divides per prime. At 256-bit that is 20 739 primes × ~27 divides per
family × 2 444 families.

- Use the existing per-prime `pinv` Lemire constants for these reductions.
- Replace `inv_u32` with a binary (shift-subtract) xgcd to remove the divides.

### 2.8 Lift the `nvar` cap

`src/engine.rs:763`:

```rust
let nvar = (s - 1).min(6); // number of sign bits varied per family
```

SIQS uses all `2^(s−1)` B-values. At 256-bit `choose_a` picks `s = 8`
(`src/engine.rs:922`, `target_bits = 112`), so the cap **discards half the polynomials per A**
and halves the amortization of 2.7's setup cost. Confirmed: 156 416 polys / 2 444 families =
exactly 64.

Measured at 256-bit, interleaved:

| build | families | run 1 | run 2 |
|---|---:|---:|---:|
| shipped `min(6)` | 2 444 | 54.27 s | 54.36 s |
| `min(9)` → nvar = 7 | 1 232 | 53.71 s | 51.19 s |

≈3% for a one-character change (smaller than naive prediction because `bainv` grows to
~581 KB — check that growth on your target before raising the cap further).

### 2.9 `choose_a` quality

Separate from the dead-zone fix in 1.1, `src/engine.rs:910-958` is weak:

- `s−1` factors are drawn **at random** (xorshift on `family`) from primes within ±1 bit of
  `ideal_bits`, then one final prime is chosen to fit the residual. There is **no rejection or
  retry** when the resulting `A` lands far from `target_a`. With 7 draws each spanning 3 bit
  widths, `A` can miss the target by several bits, which directly costs smoothness.
  msieve/YAFU keep `A` within a few percent of optimal. Add rejection sampling with a bounded
  retry.
- `src/engine.rs:949` — `let desired_u64 = desired.as_parts()[0];` takes **only the low limb**.
  Safe at current parameters, silently wrong if `desired` ever exceeds 2^64. Add a checked
  conversion.
- No duplicate-`A` detection across families, so families can collide and re-sieve identical
  polynomials.
- `all` and `pool` (`src/engine.rs:911-931`) are rebuilt by scanning the whole factor base on
  **every family** — 2 444 × 20 739 at 256-bit. Hoist into `Context`.

### 2.10 De-duplicate `find_factor` and `prepare`

`src/engine.rs:576-613` (inside `find_factor`) is a near-line-for-line duplicate of
`src/engine.rs:260-291` (`prepare`): two copies of Knuth-Schroeppel selection, `sieve_n`,
factor-base build, `pinv`, `interval_mod_p`, `target_a`, `large_prime_policy`. **They already
differ**: `:273` computes `p.sieve_half_width % e.prime` while `:588` computes
`interval as u32 % e.prime`. One will get fixed and the other will not.

Make `find_factor` call `prepare`. Determine which of the two divergent expressions is
correct, state which in the changelog, and add a test that pins it.

### 2.11 Linear algebra — the actual serial bottleneck

Per §0.3, the serial term of the browser scaling fit is ≈8.06 s, essentially all linear
algebra, ≈58% of a 48-worker 256-bit run.

**2.11a — Move the coordinator off the browser main thread.** `web/index.js:95-104`'s
`finish()` calls `coord.qs_coord_extract(session)` directly on the main thread. Worse, the
preceding `report({phase: "linalg"})` sets `status.textContent` but there is **no
`await tick()`** — `web/index.js:167` defines exactly that helper and other paths in
`factorize` use it. So the browser never paints "Linear algebra…" before locking for ~8 s.
The page is unresponsive, showing a stale status, for well over half of a 256-bit run.
Run the coordinator in its own Worker. **Highest value per unit effort in this document, and
zero algorithmic risk.**

**2.11b — Replace the `BTreeSet` filtering structures.** `src/f2/mod.rs:265-278` allocates two
`Vec<BTreeSet<usize>>` that duplicate a matrix the struct already holds in **both** CSR and
CSC (`src/f2/mod.rs:39-42`). Of the 86 MB peak on the filtered path, the echelon basis is
~24 MB — the rest is these sets. Use flat CSR with in-place compaction, which is what msieve's
`filter_relations` and cado-nfs's `purge`/`merge` do.

**2.11c — Method of Four Russians for the dense residual.** `f2_dense = 3.740 s` at 256-bit.
M4RI-style Gray-code table blocking on `row_echelon_dependencies` (`src/f2/mod.rs:187-246`) is
a 4–8× constant-factor win, far cheaper than Block Lanczos, and the correct intermediate step.

**2.11d — Vectorize back-substitution.** `src/f2/mod.rs:231-234` is a scalar `fold` with
`count_ones()` and does not call `xor`, even though that loop is O(cols²/64) per dependency
and runs 64 times. `xor_wasm_simd` (`src/f2/mod.rs:422-437`) is already on the reduction hot
loop (`src/f2/mod.rs:208`, `:168-169`) but not this one.

**2.11e — Guard the dense fallback.** `src/f2/mod.rs:347-348`:

```rust
if alive_cols.len() == ncols || alive_cols.len() <= reduced_rows { return self.dense_dependencies(); }
```

silently drops to the O(n³), uncapped-dependency, 108.6 MB path on the **original** matrix.
The margin protecting it is `base_len + 64` (`src/engine.rs:461`), i.e. 64 surplus columns;
measured slack after filtering was 474. There is no dimension check, no memory probe, no
cancellation. **Wasm linear memory never shrinks**, so one fallback permanently inflates the
tab's heap for the session. `src/factor.rs:508` takes the raw `dense_dependencies()` path
unconditionally. Add an explicit guard that returns a resource-limit error instead.

Also: `src/f2/mod.rs:152,158` allocates two `vec![0u64; …]` per column **inside** the O(n³)
loop — 41 816 allocations at 256-bit — and `verify_dependency` (O(nnz)) runs inside the loop
once per dependency (`src/f2/mod.rs:162`).

**Do not implement Block Lanczos.** At the declared 256-bit target the matrix is
20 845 × 20 908 (not the 50 k some notes assume), where Lanczos's O(n·nnz) beats echelon's
O(n³/64) by only ~2–5×. The gap grows linearly in n, so Lanczos becomes correct only if the
target moves past ~288 bits. See Phase 3.3 for what to do with the existing stub.

### 2.12 Leave the sieve inner loop alone

Recorded here so it is not "optimized" by mistake.

- The flat kernel `sieve_root_pair` (`src/engine.rs:1004-1041`) sieves both roots in one pass
  via the root-difference stride, is 2× unrolled, and monomorphizes `const SATURATING` so the
  practical range uses `wrapping_add`. Measured ≈3.5–3.7 cycles per score write.
- **Do not add SIMD to the sieve.** The inner loop is a strided scatter-add
  (`scores[pos] += w; pos += p`); SIMD128 has no scatter and strides vary per prime. msieve
  and YAFU do not vectorize it either. `small_skip = 100` (`src/engine.rs:966-969`) already
  removes the small primes where a broadcast pattern would apply. The candidate scan
  (`src/engine.rs:1185-1221`) is the classic `pmovmskb` target but is only 4.1–5.2% of
  runtime, so vectorizing it buys under 2.5% overall. It is already word-at-a-time over `u64`
  with a high-bit test enabled by the `128 - threshold` bias at `src/engine.rs:1263-1264`.
- **Do not add bucket sieving at current interval sizes.** The existing negative result was
  independently reproduced: rebuilding the blocking at 32 KiB was 48% slower at 192-bit
  (4.01 s vs 2.70 s), and a from-scratch msieve-style bucket sieve (primes < 32 KiB carried
  per block, primes ≥ 32 KiB pre-sorted into per-block buckets) also lost — 3.29 s vs 2.60 s
  at 192-bit, 8.5–10.9 s vs 6.1 s at 224-bit. The score array is 176–640 KiB against 1 MiB
  private L2, so misses land in L2 (~14 cycles), not DRAM, and bucket bookkeeping costs more
  than it saves.
- **But fix the dead gate.** `score_polynomial_blocked` (`src/engine.rs:1101-1180`) uses
  `BLOCK = 256 KiB` (`:1114`) and is gated behind `BLOCK_GATE = 1 MiB` (`:1274`) on
  `scores.len() = 2·M`. The largest `sieve_half_width` in the table is 327 680, so
  `2M = 655 360 < 1 048 576`: **the blocked kernel is unreachable at every shipped tier**,
  reachable only via `RUSQSIEVE_HALFW`. Either delete it or lower the gate and re-measure —
  do not leave an unreachable second kernel in the tree.
- **Re-measure on a small-L2 target before concluding.** The current negative result is
  overfit to a large-L2 x86 host. Browsers on mobile have 256–512 KiB L2, often shared, and
  wasm adds bounds-checked linear memory. Record the measurement even if the conclusion does
  not change.

---

## Phase 3 — API, environment, and dead code (`0.3.0`, breaking)

### 3.1 Stop reading environment variables inside the library

The library reads eight environment variables that change **numerical behaviour**, none
documented in rustdoc:

```
src/engine.rs:456    RUSQSIEVE_REL_PERCENT    — changes the relation target (i.e. the success rate)
src/engine.rs:574    RUSQSIEVE_PROFILE        — eprintln! from a library
src/engine.rs:968    RUSQSIEVE_SMALL_SKIP
src/engine.rs:972    RUSQSIEVE_SMALL_SLACK
src/engine.rs:979    RUSQSIEVE_THRESH_MARGIN
src/engine.rs:1255   RUSQSIEVE_THRESH_ADJ
src/qs/mod.rs:344    RUSQSIEVE_FB_BOUND
src/qs/mod.rs:350    RUSQSIEVE_HALFW
src/f2/mod.rs:338    RUSQSIEVE_PROFILE
src/f2/mod.rs:392    RUSQSIEVE_PROFILE
```

Every one is cached in a process-global `OnceLock`: first call wins, frozen for the process
lifetime, so two `FactorConfig`s in one process cannot differ. `RUSQSIEVE_REL_PERCENT` is
clamped to `50..110`, and 50% of the factor-base size is below the dependency threshold — it
will produce `NoNontrivialFactor` on inputs that otherwise succeed. Separately,
`std::env::var` racing a host's `setenv` is undefined behaviour in glibc, and this is exposed
from a C-callable API on threads the crate spawns itself.

**Task.** Move every knob into `FactorConfig` (that is what it is for). Read the environment
**only** in `src/bin/qs-factor.rs`, mapping env → config there. Replace the `RUSQSIEVE_PROFILE`
`eprintln!` calls (`src/engine.rs:591,711,732`, `src/f2/mod.rs:339,393`) with `tracing` or
`log` behind an optional feature, or remove them.

Keep a documented way to reproduce the tuning sweeps you ran in Phase 2 — via the CLI or a
`#[doc(hidden)]` config constructor, not via ambient environment.

### 3.2 `FactorConfig` mostly does nothing on the default path

On the `P <= PARTS` path, `src/native.rs:92` calls `engine::factor(fast_input, workers, …)`,
which reads **none** of `FactorConfig`'s `primality`, `trial_division_limit`,
`small_factor_method`, `qs`, `limits`, or `seed` (`src/factor.rs:55-60`). For the default type,
`FactorConfig` is two effective knobs wearing eight fields.

Either make the engine honour the config (preferred, and 3.1 pushes you this way) or shrink
`FactorConfig` to what it actually controls. Do not leave fields that are silently ignored.

Related: `src/native.rs:61`'s `if P <= PARTS` silently switches algorithm. `Natural<17>` falls
through to `FactorSession` → `factor_complete` → `reference_qs_factor`
(`src/factor.rs:439`) — a completely different and dramatically slower sieve
(`src/qs/mod.rs:495`) with a different relation format, no type-level or documentation signal,
and **no test coverage** (no test uses `P > 16`). Document it loudly, test it, or gate it.

### 3.3 Delete or implement the dead and mis-named code

Roughly 15–20% of the crate is unreachable scaffolding. For each item: **implement it, or
delete it.** Do not leave a third state. All removals are breaking, so they belong in `0.3.0`;
`#[deprecated]` them in `0.2.1` where possible.

**Named for something it does not do — highest priority, because the names mislead readers
about the crate's performance characteristics:**

- `src/natural/mod.rs:1098-1150` — `Montgomery` holds only `modulus` (`:1099`). There is no
  R, no R², no `n′ = −m⁻¹ mod 2⁶⁴`, no REDC. `encode` (`:1111`) and `decode` (`:1114`) are
  both `v mod m`, i.e. the identity rather than Montgomery form; `mul` (`:1117`) is
  `a.mul_mod(b, m)`. **There is no Montgomery multiplication of any kind.** Implement CIOS, or
  delete the type and rename its call sites honestly. If you implement CIOS, its `s+2`-limb
  accumulator runs into the same stable-Rust limit as Phase 2.2's `knuth_divmod`:
  `[u64; P + 2]` does not compile, so size the buffer from a literal const and slice it.
  Note that
  `Montgomery::inv` (`:1126-1149`) hand-rolls Euclid with a `mul_mod` per iteration — an O(n)
  chain of O(n²) divisions where binary xgcd is O(n²) total.
- `src/natural/mod.rs:1439-1466` — `diff_montgomery` passes **because** encode/decode are the
  identity and `mul` is plain mulmod. It asserts agreement with `(a*b) % m`, so it structurally
  cannot detect the absence of Montgomery. Delete it or make it meaningful.
- `src/f2/mod.rs:475-520` — `BlockLanczos::begin` just calls `filtered_dependencies` (`:503`),
  `request` always returns `Complete` (`:508`), `submit_product` ignores its argument
  (`:510`). Zero callers outside `f2`. Per 2.11, do not implement Lanczos now — **delete** the
  type, along with `LanczosRequest`, `F2BlockVector` (`:457-473`), and `LinearAlgebraError`.
- `src/natural/mod.rs:457-459` — `extended_gcd` returns `ExtendedGcdResult` (`:1077-1079`) with
  a single `gcd` field, discarding the coefficients; the doc comment admits it. It is not an
  extended GCD. One reference (its own definition). Implement or delete both.
- `src/qs/mod.rs:539-560` — `sieve_job` sieves `x0 = ceil_sqrt(n) + seg_start`,
  `q0 = x0² − n`: plain single-polynomial QS with advancing segments, **no A/B polynomials**.
  `SPEC.md` §3 calls this module "reference SIQS types". Rename it. Its trial division is a
  full unguarded `div_rem_u64` loop over the whole factor base (`src/qs/mod.rs:617-630`), and
  `reference_qs_factor` discards every partial relation (`src/factor.rs:472`).

**Unreachable:**

- `src/work/mod.rs:110` — `execute_job` has one reference (its definition), so
  `KernelContexts`, `WorkerScratch`, `MatrixScratch`, `ArithmeticScratch`,
  `MatrixMultiplyJob`, `MatrixMultiplyResult`, and `KernelError` are all dead. Most of the
  155-line module is unreachable.
- `src/f2/mod.rs:107-130` — `mul_m_rows` / `mul_mt_columns`, the SpMV kernels Lanczos would
  need, are reachable only via `WorkJob::MatrixMultiply` (`src/work/mod.rs:127-152`), and
  **`WorkJob::MatrixMultiply` is never constructed anywhere in the crate.** The kernels were
  written before the recurrence that would call them.
- `src/f2/mod.rs:12-31` — `MatrixSolver` / `MatrixConfig { solver, dense_threshold,
  structured_elimination_limit }` are never read; only `MatrixConfig::default()` is
  constructed (`src/qs/mod.rs:58`).
- `src/factor.rs:223-308` — `FactorSession` is a stub state machine: `advance_local` (`:257`)
  **ignores its `LocalWorkBudget` and runs the entire factorization synchronously**,
  `take_jobs` (`:286`) returns `Vec::new()`, `submit` (`:289`) compares a generation and does
  nothing. Consequently `src/native.rs:168-181`'s progress loop executes exactly once —
  `Preprocessing`, then a multi-hour block, then `Complete` — and cancellation is impossible
  on that path. Same for the whole wasm `qs_session_*` API (`src/wasm.rs:175-184`).
- `src/wasm.rs` — seven stub exports in a **versioned ABI** (`ABI_VERSION = 1`, `:7`):
  `qs_session_export_context` (`:186`), `qs_session_take_jobs` (`:190`),
  `qs_session_submit` (`:194`), `qs_session_error` (`:227`),
  `qs_worker_context_import` (`:244`), `qs_worker_context_free` (`:248`),
  `qs_worker_execute` (`:250`). Removing exports from a versioned ABI requires bumping
  `ABI_VERSION`; do that.
- `src/qs/mod.rs:393` — `CombinedRelation`, one reference.
- `src/natural/mod.rs:491` — `is_square`, one reference today, but **keep it**: §1.5's
  Baillie-PSW needs a perfect-square test before Selfridge's `D` search, because that search does
  not terminate for squares. This is the one item in this list that is not a deletion candidate.
- Never-read config fields: `FactorLimits::{max_partial_relations, max_matrix_nonzeros,
  max_memory_bytes}` (`src/factor.rs:37-39`), `QsConfig::sieve_score_scale`
  (`src/qs/mod.rs:36`), `LargePrimeConfig::{double_product_limit, enable_double}`
  (`src/qs/mod.rs:20-21`).
- Never-constructed public variants: `FactorError::InvalidRelation` (`src/factor.rs:141`);
  `ResourceLimitKind::{PartialRelations, MatrixNonzeros, Memory, PollardRhoIterations}`
  (`src/factor.rs:116-124`) — only `Relations` (`:451`) and `PolynomialBatches` (`:486`) are
  ever produced; `ProgressPhase::{CombiningRelations, BuildingMatrix, FilteringMatrix,
  PrimalityTesting}` (`src/progress.rs:141-153`) — `src/native.rs:96-104` maps only five
  engine phases plus `Complete`. Note that for the default `Natural<16>`,
  `src/native.rs:132-135` collapses every engine error except `Cancelled` into
  `NoNontrivialFactor`, so the ten-variant public `FactorError` has **three** reachable
  variants on the default path.
- Five separate deterministic Miller-Rabin implementations
  (`src/smallfactor.rs:73`, `src/qs/mod.rs:436`, `src/engine.rs:1727`, plus trial-division
  tests at `src/qs/mod.rs:181` and `src/natural/mod.rs:573`), four `powmod`s
  (`src/smallfactor.rs:48`, `src/qs/mod.rs:449`, `src/engine.rs:1766`,
  `src/natural/mod.rs:1226`), and four `xorshift`s (`src/engine.rs:1667`,
  `src/primality.rs:96`, inlined at `src/factor.rs:409-411` and `src/f2/mod.rs:547-550`).
  Consolidate into one `src/u64math.rs`. Only `src/smallfactor.rs:71-72` documents the 7-base
  witness set's provenance (Jaeschke/Sinclair); carry that comment over.

### 3.4 Cargo features

- **Three features gate nothing.** `arch-optimized` (`Cargo.toml:62`), `reference-qs`
  (`:66`), `relation-checks` (`:67`) appear nowhere in `src/`, `Makefile`, or
  `build-release.sh` — `grep -rn 'cfg(feature' src/` returns five sites, all for
  `limit-to-512-bits` or `wasm-simd128`. `cargo add rusqsieve -F arch-optimized` succeeds and
  does nothing. `SPEC.md` §3 describes them as "internal development or portability
  switches", which overstates it. **Delete them, or wire them up.** If you implement
  `relation-checks`, it should `debug_assert` that each relation satisfies `x² ≡ y² mod N`
  and that each dependency's exponent vector sums to zero over F2.
- **`limit-to-512-bits` is a non-additive feature that mutates a public type's identity.**
  `src/natural/mod.rs:72-75` sets `PARTS` to 16 or 8, and `src/natural/mod.rs:84` uses it as
  the **default const-generic parameter**: `Natural<const PARTS_64: usize = PARTS>`. Cargo
  feature unification is global and additive, so any unrelated crate in the dependency graph
  enabling this silently changes what bare `Natural` means in *your* signatures, halves the C
  ABI's maximum input, and breaks downstream code that mixes `Natural` with `Natural<16>`.
  Make `PARTS` a fixed constant and produce the 512-bit build as a separate build
  (`Makefile:57-60` already does this for wasm), or at minimum remove the default type
  parameter so the width is always explicit.

  Note while you are here: PARTS = 16 means every `Natural` temporary is 128 bytes and every
  `WideNatural` is 256 bytes (`src/natural/mod.rs:1034-1037`) for values needing 32. A
  narrow, noisy experiment (n=3) suggested ~35% on the bignum portion of `extract` from
  halving the limb count. Measure properly on the reference host before acting.

### 3.5 C ABI safety

- **`panic = "abort"` makes the C ABI's documented error contract a lie.**
  `[profile.release]` applies whenever rusqsieve is the root package — `Makefile:31`,
  `build-release.sh`, `cargo install` — so **every shipped `cdylib`/`staticlib`/CLI artifact
  is `panic = "abort"`**. Therefore `src/capi.rs:107`'s
  `catch_unwind(AssertUnwindSafe(...))` is dead code in every artifact you ship, and
  `RUSQSIEVE_INTERNAL_ERROR` (`src/capi.rs:14`, `rusqsieve.h:23`) is unreachable: a panic
  anywhere in the engine becomes `SIGABRT` in the host process, not a `c_int`. **Pick one** —
  move `panic = "abort"` to a CLI-only profile, or delete the `catch_unwind` and document in
  `rusqsieve.h` that internal errors abort the process. Shipping both is the worst option.
- **Unbounded thread spawn from an attacker-controlled argument.**
  `src/capi.rs:129-135` caps only the auto path (`threads == 0` → `min(48)`); an explicit
  `threads` passes through unvalidated into `Parallelism::threads(workers)`.
  `src/engine.rs:555-561` caps by input size only up to 184 bits — above that the arm is
  `_ => threads` — and `src/engine.rs:619` does
  `for _ in 0..threads { std::thread::spawn(...) }`. So
  `rusqsieve_factor(n_200bit, SIZE_MAX, f)` attempts to spawn 2^64 threads.
  `rusqsieve.h:57-58` claims "the engine may still cap tiny inputs", which is actively
  misleading. Clamp on **both** paths.
- **Safe Rust functions that dereference caller-supplied pointers.** `src/lib.rs:3-6`
  disables `deny(unsafe_code)` for the entire wasm32 build. An `extern "C" fn` declared without
  `unsafe` is a *safe* function — that is true in every edition, not something edition 2024
  introduced (what 2024 changed is `unsafe extern` blocks and `unsafe_op_in_unsafe_fn` becoming
  warn-by-default). So `qs_dealloc`
  (`src/wasm.rs:116`) frees an arbitrary address and `qs_coord_submit`
  (`src/wasm.rs:331`), `qs_session_new`, and `qs_worker_prepare` dereference arbitrary
  addresses via `input()`. Make all of them `pub unsafe extern "C" fn`. Add `// SAFETY:`
  comments to the three unsafe blocks (`src/wasm.rs:85,112,121`), which currently have none,
  unlike `capi.rs`. Narrow the exemption from "all of wasm32" to
  `#[allow(unsafe_code)] mod wasm;`, matching `src/lib.rs:36-38`.
- **A fabricated `'static` lifetime.** `src/wasm.rs:79-86`'s
  `fn input(...) -> Option<&'static [u8]>` manufactures `'static` from a caller pointer with
  only a bounds-vs-memory-size check. `qs_coord_submit` then holds that `&'static [u8]` live
  across `r.borrow_mut()` on `COORDS` (`src/wasm.rs:332-352`) — if the caller points into a
  registry-owned buffer, that is aliasing UB. Bind the lifetime to a token, or copy.
- **Missing ABI hygiene.** `rusqsieve.h` has no version symbol (unlike wasm's
  `qs_abi_version`, `src/wasm.rs:101`), so a shared library cannot be checked for mismatch at
  load; no `rusqsieve_strerror(int)`, so C callers hard-code status values 0–5; and
  `rusqsieve.h:70` returns `int` rather than `enum rusqsieve_status`, costing callers
  `-Wswitch` exhaustiveness (returning the enum is ABI-identical on every supported
  platform). The C ABI also exposes neither progress nor cancellation although the Rust API
  has both (`src/native.rs:39`) — for an operation that can run for hours, a C caller has no
  way out short of killing the process.

### 3.6 Public API hardening before 1.0

- **`#[non_exhaustive]` is missing** on items that will need to grow. It is applied to six
  items (`src/progress.rs:12,106,134`, `src/factor.rs:111,129,178`) but not to:
  `Parallelism` (`src/factor.rs:12`, two variants; `Cores`/`Rayon`/`ThreadPool` are obvious
  additions); `ProgressTotal` (`src/progress.rs:95`, whose siblings `ProgressUnit` and
  `ProgressPhase` *are* non-exhaustive — the inconsistency is unjustifiable);
  `ParseNaturalError` (`src/natural/mod.rs:10`); and the **public struct-variant fields**
  `ParseNaturalError::InvalidDigit { index, byte }` (`:14-19`) and
  `BufferTooSmall { required, available }` (`:49-54`), which freeze the representation
  forever — make them accessor methods. `ProgressAction` (`src/factor.rs:171`) is correctly
  left exhaustive; leave it.
- **`#[must_use]` appears on two items in the whole crate** (`src/factor.rs:71,86`). The
  critical omission is `PrimeFactors::verify_product` (`src/factors.rs:62`) — the crate's
  correctness safety net, returning `bool`, silently discardable. Add it there, and on the
  pure queries: `Natural::checked_add/sub/mul` (`src/natural/mod.rs:322,327,332`), `gcd`
  (`:434`), `bit_len` (`:143`), `to_u64` (`:114`), `bit` (`:157`),
  `write_be_bytes`/`write_le_bytes` (`:224,241` — these return a byte count),
  `PrimeFactors::{iter, expanded, multiplicity, distinct_len, total_len, is_empty}`,
  `ProgressAmount::fraction` (`src/progress.rs:83`), `Parallelism::threads`
  (`src/factor.rs:21`).
- **Iterator and conversion conventions.** No `impl IntoIterator for &PrimeFactors<P>`, so
  `for f in &factors` does not compile even though `iter()` exists (`src/factors.rs:28`).
  `iter()` returns `impl Trait` (`:28,39`) so callers cannot name or store the type.
  `Default` is public (`src/factors.rs:14`) while `new()` and `insert_count` are
  `pub(crate)`, so users can construct an empty `PrimeFactors` they can never populate.
  `Natural` has `From<u64>` (`src/natural/mod.rs:719`) but exposes `to_u64` returning
  `Option` instead of `impl TryFrom<Natural<P>> for u64`.
- **No `serde` feature**, for a crate whose central use case is moving `Natural`s across a
  wasm `postMessage` boundary — currently hand-rolled at `src/engine.rs:165-256`. Consider an
  optional feature.
- **Wire-format padding.** `src/engine.rs:172-174` writes `root: PARTS × u64` = 128 bytes per
  relation regardless of `n` (read back at `:228`); at 256 bits the root needs 32, so ~33% of
  each ~286-byte record is zero padding, ~2.0 MB of 6.0 MB over a run. **Measure before
  prioritizing**: at ~1 550 jobs with ~150 ms of sieving each against ~1 ms round-trip,
  messaging overhead is ~0.7%. Fix it for tidiness, not speed. Changing the wire format
  requires an `ABI_VERSION` bump.
- **Inconsistent documented input bound.** `src/bin/qs-factor.rs:51-56` rejects inputs over
  512 bits with a clear message, while the library `factor()` and `rusqsieve_factor` accept
  up to 1024 bits (`Natural<16>`) with no guard and no documentation. State the supported
  range on the functions and enforce it consistently.

### 3.7 Module splits

`src/engine.rs` is 1 838 lines and `src/natural/mod.rs` is 1 467. Suggested line ranges — do
these as **pure moves with no behaviour change**, in separate commits from any logic change,
so the diffs stay reviewable:

`src/engine.rs`:

| lines | responsibility | destination |
|---|---|---|
| 160-256 | wire serialization (`to_bytes` / `deserialize_family`) | `engine/wire.rs` |
| 259-393 | context prep + scheduler-agnostic collector | `engine/session.rs` |
| 395-452 | `extract` — dependency → gcd factor recovery | `engine/extract.rs` |
| 465-748 | native OS-thread scheduler (all `#[cfg(any(unix,windows))]`) | `engine/native_scheduler.rs` |
| 750-892 | `sieve_family` — SIQS self-initialization / Gray code | `engine/siqs.rs` |
| 960-1221 | tuning knobs + sieve kernel | `engine/kernel.rs` |
| 1443-1623 | `combine_cycle`, `classify_cofactor`, `Forest`, `RelationCollector` | `engine/relations.rs` |
| 1624-1776 | `inv_u32`, `mulmod_u32`, `fastmod`, `lemire_c`, `xorshift`, `knuth_schroeppel`, `is_prime64`, `powmod64` | `src/u64math.rs` |

`src/natural/mod.rs`:

| lines | destination |
|---|---|
| 8-64 | `natural/error.rs` |
| 573-690 | `natural/divmod.rs` (`sig_len`, `knuth_divmod`) |
| 776-1030 | `natural/ops.rs` — 254 lines of operator boilerplate; the `binop!`/`assign!` macros should cover all four ref permutations, but six are hand-written at `:795-830` |
| 1032-1073 | `natural/wide.rs` |
| 1075-1150 | `natural/montgomery.rs` (only if 3.3 keeps it) |
| 1152-1236 | `jacobi_u64`, `legendre_u32`, `tonelli_shanks_u32`, `modpow_u64` → merge into `src/u64math.rs` |

Also split `sieve_one_poly` (`src/engine.rs:1224-1441`, 218 lines): the per-candidate body
(`:1302-1439`) is a standalone `evaluate_candidate`.

### 3.8 Seeded witness selection

**There is no OS entropy on `wasm32-unknown-unknown`.** `getrandom` needs JS glue, and this crate
has zero runtime dependencies and no `wasm-bindgen`; `Instant::now()` and `SystemTime::now()`
panic on that target (which is why the native scheduler at `src/engine.rs:465-748` is already
`#[cfg(any(unix, windows))]`). So Miller-Rabin witness bases **cannot** be nondeterministic here.
The only available construction is a **seeded deterministic PRNG**: an explicit seed-setting call,
with the seed defaulting to 0 when it is never called.

Check whether this is additive before designing new API: `WitnessPolicy::Seeded`
(`src/primality.rs:52-55`) and `FactorConfig::seed` (`src/factor.rs:55-60`) **already exist**.
Determine whether `WitnessPolicy` is `#[non_exhaustive]`. If the existing variants can express
this, prefer wiring them up — together with §3.2's finding that `seed` is ignored on the default
path — over adding public surface. Add API only if the existing surface genuinely cannot express
it.

**PRNG choice.**

- **Do not use Mersenne Twister.** Its internal state is recoverable from 624 outputs, its seed
  diffusion is weak, and its ~2.5 KB of state is a poor fit for wasm.
- SHA-256 in counter mode is sound but is ~200 lines of new code under the zero-dependency
  constraint.
- **Recommended: ChaCha8 in counter mode** — about 50 lines, smaller and faster than SHA-256, far
  stronger than MT. Use ChaCha20 if you want more margin.
- **Do not use the existing `xorshift`** (`src/engine.rs:1667` and the three duplicates) for
  witness selection. Consolidate those per §3.3 and keep them for polynomial family selection
  only.

Document the limit honestly: with the default seed of 0 the **default configuration is fully
deterministic and therefore targetable** by an adversary who chooses `N`. That is an acceptable,
documented tradeoff — but the docs must state plainly that in the default configuration the
strength comes from Baillie-PSW (§1.5), not from witness randomization, and must state what
setting a seed does and does not buy.

Fix the existing bug in the same code while you are there: `WitnessPolicy::Seeded` computes its
witness as `Natural::from_u64(2 + rng)` reduced mod `n`, which **can be 0**;
`src/primality.rs:59-60` then `continue`s, silently consuming a round, so a caller asking for 16
rounds can get fewer. Bases must be drawn uniformly from `[2, n-2]` — use `(rng mod (n-3)) + 2`.
This path is also completely untested today (§4.3).

**Sequencing: §1.5 does not depend on this subsection.** Do not let the seed API design block
shipping Baillie-PSW.

**Acceptance:** `WitnessPolicy::Seeded` produces a base in `[2, n-2]` on every round with no
silent skips, has tests, and the same seed reproduces the same base sequence on both a native
target and `wasm32-unknown-unknown`.

---

## Phase 4 — Tests, CI, and documentation truthfulness

### 4.1 CI does not exist

There is no `.github/workflows` in the published tree. (Caveat: the crates.io `include`
allowlist would strip `.github/` even if it exists upstream — **confirm against the git
repository before acting on this item.**) Every defect in this brief passed
`cargo clippy --all-targets --all-features -D warnings`, which is clean. That is the current
quality ceiling.

Add CI covering, at minimum:

- `cargo test` and `cargo test --release`
- **A feature matrix**: default, `--no-default-features`, each feature alone, and
  `limit-to-512-bits` (which changes `deserialize_family`'s `PARTS * 8` read at
  `src/engine.rs:228`, the wasm payload width at `src/wasm.rs:370`, and the C ABI's maximum
  input, and has never been built in CI)
- `cargo clippy -D warnings`, `cargo fmt --check`
- The wasm build, via both `Makefile` and `build-release.sh` (see 4.4)
- Compiling **and running** `tests/c_api_smoke.c`

### 4.2 Two build breakages CI would have caught

- **`cargo test --no-default-features` fails.** `Cargo.toml:85-86` registers
  `[[test]] name = "cli"` unconditionally, but `tests/cli.rs` uses
  `env!("CARGO_BIN_EXE_qs-factor")` and that binary has `required-features = ["cli"]`
  (`Cargo.toml:82`). Reproduce: `rm -f target/debug/qs-factor target/debug/deps/qs_factor-*`
  then `cargo test --locked --no-default-features --test cli` → 3 passed becomes 3 failed. It
  only appears to pass because a prior default-features build leaves a stale binary. Add
  `required-features = ["cli"]` to the `[[test]]` entry.
- **`cargo test --release` cannot link from a cold cache**:
  `rust-lld: error: undefined symbol: rusqsieve::smallfactor::factor_u64::…`, plus
  `warning: output filename collision at target/release/deps/librusqsieve.rlib`. Trigger is
  `crate-type = ["rlib", "cdylib", "staticlib"]` colliding on `librusqsieve.rlib`
  (cargo#6313) combined with `panic = "abort"` and `lto = "fat"` in `[profile.release]`.
  `cargo build --release --bin qs-factor` alone is fine. **Release testing currently only
  works when artifacts happen to be cached**, which is part of why the dead zone shipped.

### 4.3 Missing tests, specifically

Credit first: `src/natural/mod.rs:1296-1467`'s `difftests` module is a genuine differential
suite against `num-bigint` — add/sub/mul/`widening_mul`, `div_rem`, `div_rem_u64`, `gcd`,
`mul_mod`/`pow_mod`, at 2 000–3 000 iterations each. Extend it rather than replacing it.

**Extend the differential suite:**

- `perfect_power` against a `BigUint` nth-root oracle — currently zero coverage of
  `src/natural/mod.rs:508`, and it is on the hot path for `n = p^k`.
- `sqrt_rem` / `floor_sqrt` / `ceil_sqrt` against `BigUint::sqrt` — zero coverage.
- `from_decimal` ↔ `Display` round-trip against `BigUint::to_str_radix(10)`. `Display`
  (`src/natural/mod.rs:725`) uses a 10^19-chunked divide with `chunks.pop().unwrap()` and is
  untested against an oracle.
- `from_be_bytes` / `write_be_bytes` / `from_le_bytes` / `write_le_bytes` round-trips
  including the leading/trailing-zero-tolerance branches
  (`src/natural/mod.rs:196,209`) — zero coverage, and `src/native.rs:62-67,138-143` depends
  on them for every result.
- `Shl` / `Shr` at word boundaries (`s % 64 == 0`, `s >= BITS`) —
  `src/natural/mod.rs:985,1003`, zero coverage.

**Factorization correctness:**

- **The advertised 192–256 bit range has no test at all.** The largest `n` anywhere in the
  suite is 128 bits (`tests/factorization.rs:73-75`, `src/engine.rs:1799-1801,1812-1813,1826-1828`
  — all the same two 64-bit primes). Add one 192-bit and one 256-bit balanced semiprime with
  known factors, `#[ignore]`d if runtime demands, and run them in CI on a schedule.
- The 65–90 bit sweep and the corpus test from Phase 1.1.
- `n = p²` for a 100+ bit `p` (exercises `perfect_power` → `factor_node` recursion,
  `src/engine.rs:530-537`).
- `n = p·q·r` with three large primes — the recursive split at `src/engine.rs:543-544` is
  never exercised beyond trial division.
- `n = 2^k · m` for large `k` — `src/engine.rs:1346-1352` trailing-zero stripping.
- `n = 2` and `n = 3` through the public `factor()`; `tests/factorization.rs:31` covers only
  0 and 1.
- The `P > PARTS` reference-QS branch (`src/native.rs:161-181`), which has no test.

**Primality (`src/primality.rs`):**

- `WitnessPolicy::Seeded` (`src/primality.rs:52-55`) is **never tested**. Its witness is
  `Natural::from_u64(2 + rng)` reduced mod n; if that is 0 it `continue`s (`:59-60`), silently
  consuming a round without testing anything — so a caller asking for 16 rounds can get
  fewer, undocumented. Fix and test, per §3.8.
- Carmichael coverage is `[561, 1105]` (`src/primality.rs:112`), both under 2^11. Add 41041,
  62745, 162401, 825265, and strong pseudoprimes to base 2 (2047, 3277, 4033), plus a large
  Carmichael. These are tests of the **Baillie-PSW path** from §1.5, and every one of them must
  reject. They are what the 16 default rounds cannot be trusted to do on their own: those rounds
  use `SMALL[round % 32]`, bases 2..53 — **deterministic, not probabilistic** — so a strong
  pseudoprime to all of bases 2..53 is a guaranteed false positive.
- **A test that the strong Lucas step is actually reached**, per §1.5's acceptance criterion — a
  counter or a test hook, asserted non-zero. No known input distinguishes a correct Lucas
  implementation from one that never runs, so without this every other primality test above
  passes on a dead Lucas step.

**Fuzzing — none exists (no `fuzz/`, no `arbitrary` dependency):**

- `engine::deserialize_family` (`src/engine.rs:198`) is the highest-value target: it is fed
  **directly from cross-worker `postMessage` data** (`src/wasm.rs:331` → `:352`). It is
  capacity-bounded (`Vec::with_capacity(count.min(1 << 20))` at `:226`, `.min(1 << 16)` at
  `:237`), but the relation graph downstream then consumes attacker-chosen `LargePrime`
  values, and `Forest::path` (`src/engine.rs:1560-1565`) walks `self.parent` with `.unwrap()`
  on `self.edge[v]`.
- `Natural::from_decimal` (`src/natural/mod.rs:164`), which parses untrusted C strings via
  `src/capi.rs:118` and untrusted wasm input via `src/wasm.rs:154`.
- `rusqsieve_factor` with non-UTF-8 input, an embedded NUL, and a 10^6-digit string.

**C ABI:**

- Compile and run `tests/c_api_smoke.c` (see 4.4).
- Concurrency: two threads running two independent `rusqsieve_factors*`, per the contract at
  `rusqsieve.h:66-68`.
- `threads = SIZE_MAX` (see 3.5).

### 4.4 The C header is unverified, and there are two divergent wasm builds

- `tests/c_api_smoke.c` exists and ships (`Cargo.toml` `include = ["/tests/**"]`) but
  `Makefile:84-86`'s `test:` target runs only `cargo test` and `make wasm`. **No C compiler
  ever sees `rusqsieve.h`.** By inspection the header currently does match `src/capi.rs`
  (five functions, `size_t`↔`usize`, `int`↔`c_int`, status values 0–5), but it is
  hand-maintained with zero verification and will drift. Generate it with `cbindgen --verify`,
  or add a make rule that compiles and runs the smoke test.
- **`Makefile:60` and `build-release.sh:305` produce different artifacts under the same name.**
  `build-release.sh:305` passes `-C target-feature=+simd128` whole-program (so LLVM
  autovectorizes the plain XOR loop to the same code as the intrinsic);
  **`Makefile:60` passes only the cargo feature, no rustflag** — and `Makefile:60` is the
  `make docs` / GitHub Pages path. Reconcile them so "the SIMD artifact" means one thing.
- Consequently the safety comment at `src/f2/mod.rs:417-418` ("The feature explicitly opts the
  whole wasm artifact into the simd128 baseline") is **false for the Makefile path**. It
  happens to be sound anyway — a module containing v128 operations fails validation wholesale,
  and `web/index.js:37-43` falls back to the scalar artifact — but the stated reason is wrong.
  Fix the comment.

### 4.5 Documentation that claims things the code does not do

Fix each of these to match the code, or change the code to match the claim. Do not leave a
mismatch.

- `README.md` advertises "single/double-large-prime relation combination" and `SPEC.md` §7.5
  describes DLP. **DLP is unreachable at every bit size** (see 2.6).
- `src/qs/mod.rs:316-321` comments "≈7k at 240, ≈9k at 256", but the table at `:338-339` sets
  bounds of 350 000 / 500 000, giving ≈14.9k and a **measured 20 844**. The 208/224 figures
  (5.7k / 11k) do match, so the comment is stale for the ≥225 tiers — and it is the stated
  rationale for the tuning.
- `SPEC.md` §3 calls the three inert features "internal development or portability switches"
  (see 3.4) and `src/qs/mod.rs` "reference SIQS types" when it is single-polynomial QS
  (see 3.3).
- `rusqsieve.h:57-58` on thread capping (see 3.5).
- `README.md`'s statement that whole-program wasm SIMD is "not used because they regressed the
  measured sieve" contradicts `build-release.sh:305` (see 4.4).

### 4.6 Documentation quality

- **No `# Panics`, `# Safety`, `# Errors`, or `# Complexity` section exists anywhere**:
  `grep -rn '# Panics\|# Safety\|# Errors\|# Complexity' src/` returns zero hits.
  `#![deny(missing_docs)]` (`src/lib.rs:2`) only checks that a doc comment exists.
  The concrete reachable case: `src/natural/mod.rs:834,840,846,852` — the `Div`/`Rem` operator
  impls `.expect("division by zero")`. These are public and reachable on ordinary input
  (`a / Natural::ZERO`). Panicking there is correct and matches `u64`, but all four need a
  `# Panics` section. Every `pub unsafe extern` fn in `src/capi.rs` needs `# Safety`
  (`:39,50,65,88` have good prose and inline `// SAFETY:` justifications — credit — but no
  rustdoc heading stating the caller's obligations).
- **Unreachable-invariant panics need their invariants written down.**
  `src/natural/mod.rs:471,474,475` (`x ≤ 2^(BITS/2)` and `q ≤ n/x`, so `x+q ≤ 2^(BITS/2+1)`)
  and `src/engine.rs:896` (`b ≤ 5a`, `a ≈ √sieve_n / interval`, so SIQS intermediates sit
  ~25 bits below capacity) are all genuinely unreachable — but the reasoning appears nowhere
  in code, comments, or `SPEC.md`. Convert to `debug_assert!` with the bound written out.
  Under `panic = "abort"` (3.5) an error in that unwritten reasoning is a process abort.
  `src/factors.rs:86`'s `.expect("factor exponent overflow")` requires 2^64 factors —
  `debug_assert` plus `saturating_add`, or document the invariant.
- **`#![deny(missing_docs)]` covers only ~1 600 of 7 163 lines**, because `engine`, `qs`, `f2`,
  `work`, `primality`, `smallfactor`, `wasm`, and `capi` are all private modules
  (`src/lib.rs:12-44`). `src/work/mod.rs` has one doc line for 155 lines; `src/primality.rs`
  has one for 116.
- **No crate-level architecture doc.** `src/lib.rs:1` is
  `#![doc = include_str!("../README.md")]`, and the README is a usage document with no module
  map, no data-flow description, and no explanation of why the crate contains **two**
  quadratic sieves. `SPEC.md` (24 KB) is not in the rustdoc. Add a module map.
- **The non-constant-time caveat is item-invisible.** It appears once, at `README.md:31`. It is
  not on `factor`, `factor_with`, or `factor_with_progress`
  (`src/native.rs:22,27,39`), not on `Natural`, not on `is_probable_prime`, and not in
  `rusqsieve.h`. The crate is keyworded `rsa`; put the caveat on every entry point.
- **The C ABI has no rustdoc page at all** — `mod capi` is private (`src/lib.rs:38`), so
  `RusqsieveFactors` and the five exported functions never appear on docs.rs. The crate's
  headline feature is invisible to its documentation host.
- **Uncoordinated magic caps** that silently turn "would have succeeded" into
  `NoNontrivialFactor`: `src/engine.rs:662`'s `max_families = 100_000` and
  `src/factor.rs:453`'s `unwrap_or(100_000)`. Coordinate them and document them.

---

## Reporting

When you finish, produce a single report containing:

1. Every task from this brief, marked done / partial / skipped, with the reason for anything
   not done.
2. For each performance change: the before and after measurement, the host, and the
   methodology.
3. The semver classification of every change, and the resulting version numbers.
4. Any file:line reference in this brief that did not match the source you found.
5. Anything you found that is not in this brief and that you believe is a defect — with
   evidence, not inference.

Do not report a task complete on the strength of a green test suite alone.
