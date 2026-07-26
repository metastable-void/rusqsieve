# Browser factorization benchmarks

The performance target is hard composite integers from 192 through 272 bits. Results are only
comparable when they use the same input, browser build, machine, worker count, warm-up policy, and
factor verification. Report both wall time and the complete prime factorization.

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
| 224 | 6.417 s | 5.222 s | −18.6% |
| 232 | 11.258 s | 8.107 s | −28.0% |
| 240 | 14.733 s | 13.629 s | −7.5% |
| 256 | 38.334 s | 33.923 s | −11.5% |

Before the monotone-echelon change, the first 272-bit case measured 92.36 s;
its phase split was 84.80 s through sieving and 7.56 s for LA/extraction. The
new LA path measures 4.55 s on that same case (−39.8%), while sieve variance
left total wall time at 91.75 s. Nearby threshold, interval, factor-base,
data-layout, and dependency-back-substitution experiments did not produce a
repeatable 272-bit end-to-end win and were not retained.

The 224- and 256-bit “after” values include the later monotone-echelon
optimization. It reduces LA/extraction itself by 21.8% and 35.2% respectively;
the smaller matrix improves as well, so no size dispatch is required for that
change.

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

The current engine has trial division, Pollard-Brent rho, and SIQS, but no ECM. Consequently it is a
strong candidate for the fastest browser **SIQS / balanced-semiprime** implementation in this range;
it is not yet defensible to call it the fastest **general** implementation. An ECM stage is the main
algorithmic prerequisite because it changes the complexity for unbalanced composites with
medium-size factors. ECM must be delivered as an opt-in feature and separate general-purpose Wasm
artifact: the balanced-RSA build remains ECM-free, and its fixed-corpus runtime, artifact size, and
module startup are regression gates.
