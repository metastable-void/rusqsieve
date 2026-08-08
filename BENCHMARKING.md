# Factorization benchmarks

The primary browser performance target is hard composite integers from 192
through 272 bits; native high-digit validation extends through the 364-bit
RSA-110 workload. Results are only comparable when they use the same input,
runtime/build, machine, worker count, warm-up policy, and factor verification.
Report both wall time and the complete prime factorization.

The established browser policy below 272 bits is a regression boundary.
Portable block Lanczos begins at 272 bits after a matched fixed-corpus run
reduced LA/extraction from 5.574 s to 1.721 s and wall time from 83.599 s to
80.072 s. A matched 256-bit control stayed within noise (27.146 s before,
26.986 s after).

For high-digit SIQS on native and Wasm, the initial automatic scale is:

| input bits | approximate digits | prime bound | half-width | LP multiplier | DLP |
|---:|---:|---:|---:|---:|:---:|
| 281–288 | 85–87 | 1,400,000 | 262,144 | 100 | no |
| 289–296 | 87–89 | 1,200,000 | 262,144 | 120 | yes |
| 297–304 | 90–92 | 1,500,000 | 262,144 | 120 | yes |
| 305–312 | 92–94 | 1,800,000 | 262,144 | 150 | yes |
| 313–319 | 95–97 | 2,250,000 | 262,144 | 150 | yes |
| 320 | 97 | 2,750,000 | 491,520 | 130 | yes |
| 321–333 | 97–100 | 3,250,000 | 524,288 | 145 | yes |
| 334–368 | 101–111 | 3,000,000 | 262,144 | 145 | yes |
| 369–400 | 112–120 | 6,000,000 | 524,288 | 145 | yes |

DLP products through 319 bits are capped at 12× or 16× the factor-base-bound
square. The exact 320-bit crossover uses an 873× cap (6.6021e15), a 100-bit
report cutoff, and 1,024-polynomial family packets. The 321–333 tier uses a
1,035× cap at B=3.25M; the 334–368 tier uses a measured 1,214× cap at B=3M.
Both give an approximately 1.093e16 product window and a 102-bit report cutoff
while retaining a 145B per-prime cap.

**The 369–400 tier is not performance-qualified.** RSA-110 at 364 bits remains
the highest tier with a competitive claim; 400 bits is the sieve's accepted
range limit (`engine::MAX_SIQS_BITS`), not a speed target, and 120-digit work
belongs to GNFS. The tier exists because the range previously inherited
RSA-110's parameters unchanged while its residues were ~2^17 larger, which
collapsed yield: a 384-bit balanced semiprime retained 276 relations from
4.9M polynomials against a 108,838 relation target.

Selection was measured on that 384-bit semiprime with 64 workers, 130–150 s per
configuration, reading relations and partials per second from
`RUSQSIEVE_PROFILE` checkpoints:

| prime bound | half-width | relations/s | partials/s | relation target |
|---:|---:|---:|---:|---:|
| 3,000,000 | 262,144 | 1.92 | 198 | 108,838 |
| 6,000,000 | 524,288 | 5.02 | 395 | 206,965 |
| 9,000,000 | 524,288 | 7.04 | 476 | 301,893 |
| 6,000,000 | 262,144 | 4.91 | 396 | 206,965 |

Cycle yield grows roughly as partials² / π(large-prime bound) and supplies most
of the late relations, so projected completion for the three fast rows lands
within about 20% of each other; the chosen row has the smallest matrix among
them. A large-prime multiplier sweep at 6M/524,288 over 145, 40, 12, and 4
moved projected completion by under 15%, so the multiplier stays at 145.

On the 96-thread Xeon 8259CL tuning host, fixed 289- and 304-bit balanced
semiprimes completed in 41.8 s and 103.5 s. The original RSA-100 result was
622.6 s (606.2 s collection, 14.9 s filtering/Lanczos/extraction). Complete
multiplier selection, Q/2 polynomials, 13-factor/2,048-variant A families,
flat report scratch, L2-sized resieve membership, and precomputed weights
first reduced the verified run to 424.9 s. Portable dense-prefix blocking,
the deeper DLP policy, allocation-free union-by-size partial combining, and
optional x86-64 SSE2 root advancement reduced the final-default run again to
355.4 s (317.2 s collection, 36.9 s filtering/Lanczos/extraction; 4.71M
polynomials). A later profiled run completed in 339.2 s (295.7 s collection
and 43.5 s filtering/Lanczos/extraction). The final 2026-07-31 release gate,
run without profiling, returned the exact RSA-100 factors in 320.1 s with a
1.48 GiB peak resident set. This does not match the clean portable YAFU
reference at 185.94 s; the remaining gap is relation-sieve throughput, not
block Lanczos.

The 364-bit RSA-110 validation used the challenge value
`35794234179725868774991807832568455403003778024228226193532908190484670252364677411513516111204504060317568667`
and 192 workers on a 384-logical-CPU host. The inherited RSA-100 implementation
returned the exact supplied factors in 392.50 s with 2,091,752 KiB peak RSS.
Runtime AVX2 root advancement and 4,096-variant RSA-110 families reduced that
to 369.70 s (−5.8%) and 2,086,912 KiB. Short profiling measured root-update CPU
per polynomial down about 21% and family-setup CPU down about 36%. Bounds of
2M, 4M, and 6M were sampled and all reduced early useful-column throughput, so
the qualified tier retains the 3M base rather than trading a larger matrix for
more smooth relations.

The subsequent matched 192-worker gate on a Xeon 6975P-C added movemask-based
report scanning and a row-major, eight-way block-Lanczos multiply. RSA-110
returned the supplied factors in 366.12 s (349.62 s collection, 15.72 s LA and
extraction) with 2,062,320 KiB peak RSS. The supplied portable YAFU 3.1.9 binary
on the same input, host, and worker count took 405.06 s and 9,412,624 KiB; its
333.77 s collection fed a 207k-row reduced matrix and its reported Lanczos time
was 210.83 s. Rusqsieve's reduced matrix was 96k rows.

At the lower crossover, a fixed balanced 320-bit semiprime
`1756278678942249845993617650632944122050045962599888588598430499675718186088530648734896935076713`
factored as
`1293851336743492807283211439665086943021713244083 ×
1357403767393126460672213357525616807202357674611`. The retained 192-worker
profiled run took 33.11 s (20.92 s collection, 11.69 s LA/extraction), compared
with portable YAFU at 33.60 s. An unprofiled 384-worker run took 30.50 s. The
older narrow-DLP default took 42.19 s; a bucket-sieve prototype regressed to
68.08 s and was removed.

The intervening 330-bit RSA-100 gate used
`1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139`.
The retained 3.25M/524,288 geometry, 102-bit report cutoff, 1.0932e16 DLP
window, and 2,048-pattern families returned its published factors in 54.01 s
(36.46 s collection, 16.81 s LA/extraction) with 1,995,804 KiB peak RSS. This
is an 86% reduction from the old 320.1 s release gate, but remains 6.5% behind
the matched portable YAFU result of 50.70 s (30.44 s collection and 19.16 s
reported Lanczos) on this one midpoint. Exact YAFU-scale 3.22M geometry took
57.35 s, 1,024-pattern families took 57.39 s, and deeper reporting took
55.78 s, so all three were rejected. The endpoint gates above and below this
input beat YAFU; no all-input parity claim is inferred from them.

The report-scan microbenchmark processed 2,000 524,288-byte streams in 19.90 ms
with AVX2 versus 270.02 ms scalar. On a synthetic 50k-row/5M-nonzero matrix,
row-major `B·x` took 22.6–24.7 ms per ten loops versus 47.6–48.0 ms for the old
column scatter; scoped eight-way multiplication reduced a paired forward/
transpose pass to about 9.6/9.6 ms. Node 26.5.1 end-to-end architecture checks
also returned the expected factors from both scalar Wasm and `simd128` Wasm.

The 281–288 tier removes the old drop from a 700k prime bound at 280 bits to
500k at 281. On an exact 288-bit balanced fixture, the old 500k/327,680 policy
needed 10,352 families and 420.223 s; 1.4M/262,144 with Lanczos needed 4,328
families and 258.512 s (−38.5%). This is a one-input boundary measurement, so
retain it as a strong correction to the discontinuity rather than a
multi-host optimum claim. The exact input and its two 144-bit factors are
recorded as the 288-bit row in `tests/data/browser-balanced-corpus.txt`. A
final real-Chromium run of the shipped SIMD artifacts took 23.702 s at 256
bits, 69.783 s at 272 bits, and 226.901 s at 288 bits; all three results were
factor-verified with eight Web Workers. Post-tuning checks of the same endpoint
fixtures took 67.433 s at 272 bits and 222.006 s at 288 bits. Wasm SIMD root
advancement plus the portable scorer subsequently reduced those verified
endpoints to 62.237 s and 185.830 s.

## 0.4 release-gate measurements

The final 0.4 source and artifacts were rechecked on 2026-07-31. These are
single-run release gates, not replacements for the multi-input tuning means:

| path | input | workers | wall | peak RSS |
|---|---:|---:|---:|---:|
| Node 24.15/V8 SIMD Wasm | 192 bits | 8 | 0.773 s | — |
| Node 24.15/V8 SIMD Wasm | 224 bits | 8 | 3.973 s | — |
| Node 24.15/V8 SIMD Wasm | 256 bits | 8 | 27.293 s | — |
| native release | 288 bits | 48 | 37.81 s | 297 MiB |
| native release | RSA-100 (330 bits) | 96 | 320.10 s | 1.48 GiB |

The shipped `docs/` frontend also completed the first 216-bit browser-corpus
case in headless Chromium with eight Web Workers and SIMD: 2.225 s total,
including 1.979 s through sieving and 0.246 s for linear algebra/extraction.
Every row returned factor-verified output. At this revision the scalar and
SIMD Wasm modules are 196,216 and 202,334 bytes.

## Reproducible rusqsieve harness

`tools/wasm-bench.mjs` runs the browser architecture under Node/V8: one coordinator Wasm instance,
independent worker instances, and serialized relation packets. It does not use shared memory.

```sh
make wasm
RUSQSIEVE_WASM=target/wasm-simd/wasm32-unknown-unknown/release/rusqsieve.wasm \
  node tools/wasm-bench.mjs DECIMAL 8 2
```

The fourth argument is polynomial families per worker job. Two is the measured default. The harness
rejects a result unless the returned divisor is nontrivial and divides the input.

The CI host's unpacked Playwright installation also needs its `/tmp` font
root:

```sh
FONTCONFIG_FILE=/tmp/rusqsieve-font-root/etc/fonts/fonts.conf \
FONTCONFIG_SYSROOT=/tmp/rusqsieve-font-root \
LD_LIBRARY_PATH=/tmp/rusqsieve-chromium-libs/root/usr/lib/x86_64-linux-gnu \
PLAYWRIGHT_MODULE=/tmp/rusqsieve-playwright/node_modules/playwright \
PLAYWRIGHT_BROWSERS_PATH=/tmp/rusqsieve-playwright-browsers \
  node tools/playwright-bench.mjs URL N P Q 8
```

Without the two fontconfig variables, Chromium's Skia font manager exits
during the first navigation.

For a real browser measurement, serve `docs/` and point the Playwright driver at
it. Playwright remains an external benchmark dependency rather than a crate or
website dependency:

```sh
make docs
make serve
# In another shell:
PLAYWRIGHT_MODULE=/path/to/node_modules/playwright \
PLAYWRIGHT_BROWSERS_PATH=/path/to/playwright-browsers \
  node tools/playwright-bench.mjs \
  http://127.0.0.1:8000/ DECIMAL EXPECTED_P EXPECTED_Q 8
```

The driver forces the requested `hardwareConcurrency`, exercises the real
coordinator and browser Worker pool, requires an exactly verified
factorization, and reports first-relation, sieve-end, and LA/extraction timing.
The fixed 216/224/232/240/256/272-bit tuning corpus is
`tests/data/browser-balanced-corpus.txt`.

## Fixed balanced-semiprime corpus

| bits | input | factors |
|---:|---|---|
| 192 | `5845354724375454473909137928398990449217655808523662886639` | `75335908545075305094962839541 × 77590551932854658187989536979` |
| 224 | `21523772555907914536866856055060033603780528151558474367883009969243` | `4146060183335910751156909939294247 × 5191379672301316010974896170794669` |
| 256 | `98877949376972157840865984674312121822345015130827118595228756728313751597271` | `303899915024639499827896288126367369941 × 325363530848487941099032348913090235131` |

The browser-target corpus contains five balanced, exact-width inputs at 216,
224, 232, 240, 256, and 272 bits. On the development server with headless
Chromium, SIMD, and eight workers, tier retuning changed the five-case means as
follows:

| bits | before | after | change |
|---:|---:|---:|---:|
| 216 | 5.097 s | 3.266 s | −35.9% |
| 224 | 6.417 s | 5.075 s | −20.9% |
| 232 | 11.258 s | 7.757 s | −31.1% |
| 240 | 14.733 s | 13.629 s | −7.5% |
| 256 | 38.334 s | 32.259 s | −15.8% |

Before the LA changes, the first 272-bit case measured 92.36 s; its phase split
was 84.80 s through sieving and 7.56 s for LA/extraction. Monotone echelon
truncation first reduced LA/extraction to 4.55 s. The retained eight-pivot M4RI
solver reduces it again to 2.79 s (−63.1% cumulatively), with a measured total
of 92.07 s. The nearby threshold, interval, and factor-base experiments at
that point predated M4RI and did not produce a repeatable end-to-end win; the
post-M4RI sweep below changes that balance. Data-layout and
dependency-back-substitution non-wins remain unretained.

The 224- and 256-bit “after” values include both retained LA optimizations.
After monotone truncation, M4RI further reduces five-case LA/extraction means
from 0.374 s to 0.317 s at 224 bits (−15.2%) and from 2.803 s to 1.769 s at
256 bits (−36.9%). A 3,200-column crossover keeps smaller residuals on scalar
echelon, and a 64 MiB estimated working-set cap prevents panel tables from
expanding unboundedly.

After M4RI landed, the factor-base/interval/threshold trade-off was swept again
in the same browser architecture. Five-case confirmation retained:

| bits | previous parameters | new parameters | previous mean | new mean | change |
|---:|---|---|---:|---:|---:|
| 224 | `(150k, 131072, 0)` | `(175k, 131072, 0)` | 5.176 s | 5.075 s | −2.0% |
| 256 | `(400k, 196608, −5)` | `(450k, 196608, −4)` | 32.917 s | 32.259 s | −2.0% |
| 272 | `(500k, 327680, −4)` | `(700k, 262144, −4)` | 105.110 s | 94.880 s | −9.7% |

The 216-, 232-, and 240-bit settings remain unchanged. In particular, the
post-M4RI 232-bit 200k baseline averaged 7.757 s while 250k averaged 7.774 s.
At 272 bits, 800k exceeded the bounded M4RI working set and fell back to scalar
elimination; 196,608 and 327,680 half-widths both lost to 262,144 at 700k.

On the development host with Node 24.15/V8 and eight workers, the scoped-SIMD build measured 0.72 s,
5.04 s, and 37.86 s respectively on 2026-07-25. These are factor-verified engineering measurements,
not cross-project records. Scaling the 256-bit case to 16, 32, and 48 workers measured 22.26 s,
14.71 s, and 13.96 s. The browser caps the pool at 48 because 96 workers regressed on the reference
host from startup, memory traffic, and relation overshoot.

## Pollard–Brent reach above the sieve ceiling

A composite wider than `engine::MAX_SIQS_BITS` (400) never reaches the sieve, so rho is the whole
factoring attempt there and its budget is a wall-clock decision rather than a fraction of a sieve
run. Iteration rates, release build, single-threaded, x86-64 Xeon 8259CL — reproduce with
`cargo test --profile release-test --lib -- --ignored profile_wide_rho_throughput --nocapture`:

| composite bits | iterations/s | budget | wall | smallest factor reached (`1.2·√p`) |
|---:|---:|---:|---:|---:|
| 400 | 6,881,194 | 6,291,456 (sieve-derived) | 0.9 s | 2^44.6 |
| 512 | 4,842,288 | 128,000,000 | 26 s | 2^53.0 |
| 768 | 2,203,329 | 72,000,000 | 33 s | 2^51.7 |
| 1024 | 1,338,216 | 48,000,000 | 36 s | 2^50.5 |

End to end through the release CLI, on inputs built as one small prime times one wide prime. The
old-budget column is the same binary run with `RUSQSIEVE_RHO_ITERATIONS=6291456`, which is exactly
what 0.4.2 spent at every width above the ceiling, so both columns are the same inputs on the same
host. The last row is the cost side of the decision: a hopeless composite is now refused only after
the budget is spent.

| input | smallest factor | 6.29M budget (0.4.2) | default (0.4.3) |
|---|---:|---:|---:|
| 401-bit | 16-bit | 0.11 s | 0.11 s |
| 512-bit | 32-bit | 0.11 s | 0.11 s |
| 1024-bit | 32-bit | 0.21 s | 0.21 s |
| 512-bit | 40-bit | 1.51 s | 1.51 s |
| 512-bit | 48-bit | refused, 2.71 s | 30.35 s |
| 512-bit balanced semiprime | — | refused, 2.71 s | refused, 61.99 s |

A cofactor that split under rho keeps the deep budget from 257 bits up, which is what lets a wide
product of middling primes finish at all. Measured through the release CLI with 96 workers:

| input | before | after |
|---|---:|---:|
| 498-bit, ten 50-bit primes | two factors peeled, then a 399-bit sieve wanting 206,403 relations at ~2/s | 65.3 s, all ten primes verified |

The run peels five factors in rho (498 → 448 → 399 → 349 → 299 bits) and hands the 250-bit
remainder to a sieve that builds its factor base in 0.03 s and collects in 2.9 s.

Nothing at or below the ceiling changed. The budget there is still the sieve-derived one, pinned
value-for-value by `budgets_at_and_below_the_ceiling_are_unchanged`, because rho finds nothing on a
balanced semiprime and every iteration it spends is overhead on the main workload. Confirmed against
a 0.4.2 binary built from the previous commit, interleaved with alternating order, median of three
after a warm-up, `--threads 8` (32 for the last two rows):

| input | 0.4.2 | 0.4.3 | delta |
|---|---:|---:|---:|
| 128-bit balanced | 0.032 s | 0.031 s | −3.1% |
| 192-bit balanced | 0.378 s | 0.371 s | −1.8% |
| 216-bit balanced | 1.498 s | 1.473 s | −1.7% |
| 224-bit balanced | 2.925 s | 2.892 s | −1.1% |
| 300-bit, 32-bit factor | 0.038 s | 0.038 s | +0.7% |
| 384-bit, 40-bit factor | 0.574 s | 0.565 s | −1.6% |
| 256-bit balanced | 6.766 s | 6.548 s | −3.2% |
| 272-bit balanced | 16.033 s | 15.720 s | −1.9% |

Every delta is host noise around an identical code path: below the ceiling the two binaries execute
the same budget through the same loop.

Deeper searches are the caller's decision, because each factor bit doubles Brent's cost:

```sh
RUSQSIEVE_RHO_ITERATIONS=4000000000 qs-factor < composite.txt
```

`raised_budgets_reach_56_and_64_bit_factors_above_the_ceiling` is the same search as a test: a
56-bit factor out of a 512-bit input plus a 64-bit factor out of a 401-bit input, single-threaded,
1,452 s for the pair on this host. Run it with `--nocapture` for the per-case split.

The browser runs the same loop in wasm across a pool of rho workers (`qs_rho`), each walking a
disjoint range of polynomial constants. Iteration rates measured under Node 24.15 on prime moduli:

| composite bits | native | wasm (scalar) | main-thread `BigInt` |
|---:|---:|---:|---:|
| 512 | 4,842,288/s | 1,075,000/s | 288,000/s |
| 1024 | 1,338,216/s | 315,000/s | 115,000/s |

Per-worker budget is 2^25 iterations up to 512 bits, 24M through 768, 16M above. `T` independent
walks collide in about `1.2·sqrt(p)/sqrt(T)` iterations, so eight workers reach a smallest factor of
roughly 2^52 at 512 bits and 2^50 at 1024 — parity with the native CLI — in 31 to 53 s of worker
time, with the main thread free throughout.

End to end on a 478-bit product of ten 48-bit primes, driving the real worker protocol on Node
worker threads with eight workers:

| ladder | result |
|---|---|
| `BigInt`, main thread, 2^23 budget | one factor in 12.1 s, then 33.2 s exhausted at 430 bits → sieve/refusal |
| `BigInt`, main thread, 2^26 budget | five factors, 140 s of blocked main thread |
| wasm pool, 2^25 per worker | five factors in **34.1 s** (5–8 s each), main thread free, 239-bit remainder to the sieve |

A runtime with no module or no `Worker` falls back to the main thread's sliced `BigInt` search,
which keeps the smaller 2^23/2^22 budget because it has to stay sliceable into ~50 ms macrotasks.

## Montgomery arithmetic and the rho inner loop

Rho is two modular multiplications per iteration, so `src/natural/montgomery.rs` is where its speed
lives. Reproduce the inner loop with
`cargo test --profile release-test --lib -- --ignored profile_montgomery_loop --nocapture`; it runs
one Montgomery squaring, one modular add, one modular subtract, and one Montgomery multiply per
iteration, which is exactly what the reference below runs.

The reference is the same loop over GMP 6.3.0's mpn assembly (`__gmpn_sqr`, `__gmpn_mul_n`,
`__gmpn_addmul_1`). That is the fair opponent for a "world-class" claim: YAFU's `montybrent` runs on
GMP through mpz, and FLINT's `flint_mpn_factor_pollard_brent_single` runs on mpn directly, so both
inherit this assembly. x86-64 Xeon 8259CL, one million iterations each:

| limbs | bits | rusqsieve | GMP mpn | ratio |
|---:|---:|---:|---:|---:|
| 2 | 128 | 19,337,101/s | 13,494,324/s | **1.43×** |
| 3 | 192 | 13,901,104/s | 9,579,733/s | **1.45×** |
| 4 | 256 | 10,276,357/s | 7,811,153/s | **1.32×** |
| 5 | 320 | 7,264,293/s | 5,513,344/s | **1.32×** |
| 7 | 448 | 4,073,624/s | 3,551,447/s | **1.15×** |
| 8 | 512 | 3,497,736/s | 3,407,755/s | **1.03×** |
| 12 | 768 | 1,515,986/s | 1,663,963/s | 0.91× |
| 16 | 1024 | 916,415/s | 1,088,155/s | 0.84× |

We are ahead through seven limbs, at parity at eight, and behind by 7–13% at twelve and sixteen,
where GMP's hand-written `mulx`/`adcx`/`adox` inner loops run two carry chains at once.

That technique was tried and rejected, and the sequence is worth recording because two plausible
measurements were both wrong:

1. A first probe using `_mulx_u64` and `_addcarry_u64` measured 1.15–1.19× on the multiply and
   looked like an ADX win. `objdump` showed **184 `mulx` and zero `adcx`/`adox`** — LLVM had lowered
   the carry intrinsics to plain `adc`, so the two-chain mechanism never ran. The gain was `mulx`
   alone, which frees the fixed `rdx:rax` operand pair.
2. Wiring `mulx` in properly needs care: a `#[target_feature(enable = "bmi2")]` wrapper only
   recompiles what is *inlined into it*, and the dispatch table is far too large for the inliner to
   take on a hint, so the first integration emitted no `mulx` at all. With an `#[inline(always)]`
   entry point the release binary contained 776.

Measured then against an otherwise identical build with selection forced off, six interleaved runs
each, medians:

| limbs | 2 | 3 | 4 | 5 | 7 | 8 | 12 | 16 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BMI2 vs portable | +2.2% | +8.1% | +8.5% | +4.1% | +5.5% | −0.8% | −1.0% | −4.0% |

It helps from 192 to 448 bits and hurts from 512 to 1024 — not a consistent win even on the machine
it was measured on, let alone a portable one, so it is not shipped. The dispatch also has to be an
inlined branch rather than a stored function pointer: selecting through a `fn` pointer blocked
inlining at the call site and cost 4–11%, more than the feature returned.

Two-chain ADX remains the real remaining gap at twelve and sixteen limbs, and reaching it means
inline assembly rather than intrinsics.

End to end against YAFU 3.1.9 (GMP 6.3.0, GMP-ECM 7.0.5), ten fixed 512-bit composites each with a
40-bit factor, `rho(N)` with `-rhomax 40000` against `qs-factor --threads 1`, both verified to
return the factor:

| | median | note |
|---|---:|---|
| rusqsieve | 0.442 s | 0.006 s of that is process startup |
| YAFU | 1.196 s | 0.227 s of that is process startup |

That is 2.2× on the rho work itself, and rusqsieve won 9 of the 10 inputs. Repeated at 1024 bits
with a 40-bit factor over eight inputs: 1.834 s against 3.309 s, winning 8 of 8. The walks differ
between implementations, so single inputs are luck; the medians are the claim.

Where the speed came from, measured as full-stage `profile_wide_rho_throughput` rates before and
after `src/natural/montgomery.rs` was rewritten:

| bits | before | after Montgomery | after gcd + batch | total |
|---:|---:|---:|---:|---:|
| 128 | 5,812,766/s | 17,776,441 | 30,272,029 | **5.21×** |
| 256 | 4,306,730/s | 9,210,792 | 13,367,638 | **3.10×** |
| 512 | 2,194,421/s | 3,729,081 | 4,842,288 | **2.21×** |
| 1024 | 730,330/s | 1,179,364 | 1,338,216 | **1.83×** |

Three changes account for it: the arithmetic works in place over the modulus's significant limbs
instead of copying through a zeroed 33-word scratch buffer (which is why the small widths gain
most), the inner loops are monomorphized over the limb count and unrolled, and squaring uses the
symmetric product. Correctness is pinned by `diff_montgomery_arithmetic`, which checks every limb
count from 1 to 16 against division-based modular arithmetic.

### Limb width on wasm

wasm has `i64.mul` but no widening 64×64 multiply, so 64-bit limbs there mean every product is
emulated. Measured through `qs_rho` on the shipped artifact under Node 24.15 — same build, same
moduli, only the limb type differing:

| bits | 320 | 400 | 512 | 640 | 768 | 896 | 1024 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 32-bit limbs | 2.004 M/s | 1.394 | 0.989 | 0.698 | 0.515 | 0.401 | 0.299 |
| 64-bit limbs | 2.097 M/s | 1.083 | 0.826 | 0.554 | 0.373 | 0.281 | 0.219 |

wasm therefore uses 32-bit limbs and every other target uses 64-bit ones. The crossover below 400
bits costs nothing the browser's deep rho cares about, since that path exists for composites the
sieve refuses. Both limb widths are generated from one macro and checked against each other by
`narrow_and_wide_limbs_agree`, so the host test suite covers the wasm arithmetic.

## Elliptic curve method

ECM's cost is governed by the size of the factor it is looking for, not the size of the input,
which is what makes it the right stage between Pollard–Brent and the sieve. Reproduce the reach
measurement with
`cargo test --profile release-test --lib -- --ignored ecm --nocapture`.

Single-threaded on an x86-64 Xeon 8259CL:

| factor | composite | bounds | result |
|---|---|---|---:|
| 10 digits (2^32) | 161-bit | `B1 = 2,000`, 64 curves | 0.11 s (debug build) |
| 20 digits (2^66) | 260-bit | `B1 = 50,000`, 300 curves | 8.9 s |

End to end through the release CLI, on a 466-bit composite that is a 20-digit prime times a 400-bit
prime — the shape that has no other stage:

| | before | after |
|---|---|---|
| 466-bit, 20-digit factor | `SiqsCompositeTooLarge` | 29.8 s |

Most of that 29.8 s is the deep rho that has to fail first; rho's budget at this width reaches about
2^53 and the factor is 2^66.

For scale against the neighbours: Pollard–Brent needs roughly `1.2·sqrt(p)` iterations, so a 2^66
factor is about 10^10 iterations — hours at the measured 4.8 M/s — and the sieve does not accept a
466-bit composite at all.

The browser runs the same stage through `qs_ecm` across a pool of workers. `tools/ecm-check.mjs`
splits a 161-bit composite in about 60 ms under Node 24.15.

### What enabling it costs

Curves run unasked only where a balanced semiprime cannot be: above the sieve's ceiling, or on a
composite trial division or rho has already shown to be unbalanced. Inside the sieve's range with no
such evidence the switch is the caller's, and the measured cost of flipping it is nil — a 256-bit
balanced semiprime, 32 threads, median of three:

| | default | `--enable-ecm` |
|---|---:|---:|
| 256-bit balanced semiprime | 6.72 s | 6.55 s |

The difference is noise: at that width the schedule is `B1 = 2,000` over 16 curves against a
6.7-second sieve run. The reason to keep it off by default is that no curve can succeed on a
balanced semiprime, not that the curves are expensive.

## The batched gcd

Pollard-Brent accumulates differences and takes one gcd per batch, so the batch size trades gcd cost
against work done past a collision before the next gcd notices it. The gcd is also the one part of
the loop that does not get cheaper when the Montgomery arithmetic does, so making the multiply 1.7×
to 3× faster raised the gcd's share of the stage and made this worth retuning.

`Natural::gcd` was a binary GCD written over whole `Natural`s, so every one of roughly a thousand
iterations touched all sixteen limbs whatever the operands' width. It now works over the significant
prefix, which only shrinks, and finishes in machine arithmetic once both operands fit in one word.
Reproduce with `cargo test --profile release-test --lib -- --ignored profile_gcd --nocapture`:

| bits | 128 | 256 | 512 | 768 | 1024 |
|---|---:|---:|---:|---:|---:|
| whole-`Natural` | 6,473 ns | 12,327 | 22,976 | 33,918 | 45,994 |
| prefix-aware | 1,398 ns | 3,552 | 8,440 | 14,719 | 22,870 |
| | **4.6×** | **3.5×** | **2.7×** | **2.3×** | **2.0×** |

The old cost is visible in the shape as well as the size: 128-bit gcds were only 7× cheaper than
1024-bit ones, where the algorithm says 64×, because fixed-width work dominated at the small end.

Full-stage rates at each batch size, before and after that fix:

| bits | 128 (old gcd) | 1024 (old gcd) | 128 (new gcd) | 512 (new gcd) |
|---|---:|---:|---:|---:|
| 128 | 20,558,118/s | 29,861,389 | 28,163,988 | 29,394,392 |
| 512 | 3,766,866/s | 4,780,142 | 4,502,406 | 4,895,467 |
| 1024 | 1,150,572/s | 1,331,569 | 1,243,044 | 1,331,044 |

The gcd was 16–45% of the stage at batch 128; the rewrite recovers most of that, leaving 5–9% for
the batch, which is now 512. End to end on the ten 512-bit composites with a 40-bit factor: 0.481 s
median at batch 128 against 0.443 s at 512, no failures either way. Batch sizes to 2048 were also
run against 28 short-cycle inputs (16- to 28-bit factors) with no failures, which is structural
rather than luck: `batch` is `min(r - k, B)`, so `B` only binds once `r` exceeds it and short
searches never reach that.

## Measurement failures worth not repeating

Every entry here produced a confident, wrong number first. They are recorded because the failure
modes are the reusable part.

**A benchmark that optimized itself away.** A microbenchmark of the Montgomery multiply reported
rates of 22–88 M/s *per million iterations* — physically impossible. `black_box` on the operands
fixed it. Any rate that implies less than a cycle of work is a dead loop, not a fast one.

**A layout probe that could not probe.** Testing whether a timing difference was code alignment, the
"perturbation" was a trailing comment. Comments do not reach codegen, so the experiment could only
ever return "no change". A probe has to emit instructions.

**A sweep that was not interleaved.** The first alignment comparison ran all of A then all of B and
reported the 224-bit case 14% *slower* under alignment. Interleaved with alternating order over five
runs, the same comparison was 3.1% faster. Drift on a shared host is larger than the effects being
measured; alternate the order or measure nothing.

**An edit that never applied.** A `sed` retuning the wasm specialization table matched nothing,
because `cargo fmt` had wrapped the macro invocation across lines. The build was byte-identical, and
the "no size change, no speed change" result was read as a finding rather than as a no-op. Confirm
the artifact changed before believing that it did not matter.

**A probe that measured the wrong thing.** 32-bit limbs looked 1.67× faster than 64-bit ones in
wasm. Both arms of that probe were run-time-width loops, while the shipped code specializes the
limb count; measured properly on the real artifact the win was 20–41% and negative below 400 bits.
A probe has to be built the way the thing it stands in for is built.

**An intrinsic that did not do what its name said.** `_addcarry_u64` alongside `_mulx_u64` measured
1.15–1.19× and was reported as an ADX result. Disassembly: 184 `mulx`, zero `adcx`/`adox`. The
mechanism named in the conclusion was absent from the binary.

**A feature flag enabled on nothing.** Wiring that same path in behind
`#[target_feature(enable = "bmi2")]` produced a binary with no `mulx` at all: the wrapper only
recompiles what is inlined into it, and the dispatch table was far past the inliner's threshold. The
fix was an `#[inline(always)]` entry point; the check that caught it was `objdump | grep -c mulx`.

**A dispatch that cost more than it dispatched.** Selecting the BMI2 routines through a stored `fn`
pointer blocked inlining at the call site and cost 4–11%, more than the feature returned, so the
"BMI2 build" was slower than the portable one for reasons unrelated to BMI2.

## Basic-block alignment in native builds

Timings of this crate's hot loops are sensitive to code alignment under `lto = "fat"` and
`codegen-units = 1`. Two inert lines added to the CLI's progress closure — in a match arm that never
executed, because progress was disabled — moved a 512-bit input with a 40-bit factor from 1.45 s to
1.57 s, while the same tree's library iteration rate, measured by `profile_wide_rho_throughput`
under `--profile release-test` (thin LTO, 16 codegen units), was unchanged at 2.18–2.80 M/s. Any A/B
of unrelated changes on this codebase can pick up several percent of that noise.

`make native` and `build-release.sh` therefore build with
`-C llvm-args=-align-all-nofallthru-blocks=5`. Medians of five interleaved runs against the same
tree built without it, x86-64 Xeon 8259CL:

| input | default | 32-byte alignment | delta |
|---|---:|---:|---:|
| 192-bit balanced | 0.210 s | 0.202 s | −4.0% |
| 216-bit balanced | 0.675 s | 0.658 s | −2.6% |
| 224-bit balanced | 1.163 s | 1.127 s | −3.1% |
| 232-bit balanced | 1.972 s | 1.931 s | −2.1% |
| 256-bit balanced | 6.089 s | 5.973 s | −1.9% |
| 272-bit balanced | 15.398 s | 15.049 s | −2.3% |
| 512-bit, 40-bit factor | 1.566 s | 1.431 s | −8.6% |
| 384-bit, 40-bit factor | 0.607 s | 0.569 s | −6.3% |

`qs-factor` grows from 773 KiB to 893 KiB. 64-byte alignment (`=6`) was measured at the same speed
for 38% more size and was rejected. A first, non-interleaved sweep of the same comparison reported
the 224-bit case 14% *slower* under alignment; interleaving with alternating order removed that
entirely, which is the reason the protocol above insists on it.

## What a “fastest general browser factorizer” claim requires

A balanced-semiprime corpus measures the SIQS worst case well, but it is not a representative
general-factorization corpus. Before making the broader claim, benchmark at every input tier:

- balanced two-prime composites;
- composites with 32-, 48-, 64-, 80-, and 96-bit smallest factors;
- three-or-more-prime composites;
- repeated prime powers;
- primes and near-prime perfect powers.

Run current Alpertron, the CrypTool Msieve WebAssembly port, and Yaffle's browser QS on the same
browser and hardware. Use at least five fixed inputs per cell, alternate implementation order, warm
each artifact once, publish medians and ranges, and retain raw results.

The current engine has trial division, Pollard-Brent rho, and SIQS, but no ECM. Independent
cross-factorizer testing on multiple mobile devices establishes it as the world's fastest browser
**SIQS / balanced-semiprime** implementation in this range. That specific result is distinct from a
"fastest general factorizer" claim: ECM changes the complexity for unbalanced composites with
medium-size factors. ECM must be delivered as an opt-in feature and separate general-purpose Wasm
artifact; the balanced-RSA build remains ECM-free, and its fixed-corpus runtime, artifact size, and
module startup are regression gates.
