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

DLP products through 319 bits are capped at 12× or 16× the factor-base-bound
square. The exact 320-bit crossover uses an 873× cap (6.6021e15), a 100-bit
report cutoff, and 1,024-polynomial family packets. The 321–333 tier uses a
1,035× cap at B=3.25M; the 334–368 tier uses a measured 1,214× cap at B=3M.
Both give an approximately 1.093e16 product window and a 102-bit report cutoff
while retaining a 145B per-prime cap.

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
