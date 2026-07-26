# Changelog

## 0.3.0 — 2026-07-26

Breaking minor release for the API, environment, dead-code, safety, CI, and
documentation work deferred from 0.2.1.

### Configuration and supported path

- The library no longer reads `RUSQSIEVE_*` environment variables. All
  numerical tuning is owned by `FactorConfig`; only `qs-factor` maps the
  benchmark environment into a doc-hidden tuning constructor. Relation targets
  are clamped above the dependency threshold rather than accepting the old
  failure-inducing 50% setting.
- `FactorConfig` now contains only controls the optimized engine honors:
  parallelism, progress cadence, an optional deterministic witness seed, and
  internal benchmark tuning. The silently ignored reference-QS, resource-limit,
  primality, trial-division, and small-factor fields were removed.
- Every `Natural<P>` value through 512 significant bits uses the optimized SIQS
  engine, including `P > 16`. Wider values consistently return
  `FactorError::InputTooLarge`; the dramatically slower and uncancellable
  single-polynomial fallback was deleted.
- Seeded Miller–Rabin witnesses use ChaCha8 and full-width rejection sampling
  over `[2, n−2]`. No round can be silently consumed by a zero base. Equal
  seeds reproduce the same stream on native and Wasm; the default remains
  deterministic and Baillie–PSW remains the primary primality safeguard.

### Removed misleading and unreachable surfaces

- Deleted the public-in-private-module `Montgomery` façade and its
  non-Montgomery differential test. The old test only asserted ordinary
  `(a·b) mod m` semantics and could not detect the missing Montgomery domain.
- Deleted the fake `BlockLanczos` protocol, `F2BlockVector`, never-used sparse
  matrix multiply jobs/kernels, false identity `provenance` accessor,
  `MatrixSolver`/`MatrixConfig`, `extended_gcd` result without coefficients,
  always-zero work metrics, `CombinedRelation`, and the entire unreachable
  `work` module.
- Deleted the synchronous `FactorSession` state-machine stub and every
  `qs_session_*`/worker-import stub export. The active coordinator/worker Wasm
  protocol remains; the raw Wasm ABI is now version 2.
- Reduced `qs` to what the optimized engine actually consumes: factor-base
  construction and tier parameters. The segmented `x²−N` reference sieve and
  its SIQS-claiming names are gone.
- Consolidated deterministic `u64` Miller–Rabin, modular exponentiation, modular
  multiplication, and xorshift into `u64math`, retaining the documented
  Jaeschke/Sinclair witness provenance.
- `engine` wire serialization, dependency extraction, and SIQS
  self-initialization were split into focused submodules. The larger suggested
  engine-kernel/relation and `natural` pure-move splits remain deferred because
  combining them with this release's logic changes would make review harder.

### Public API and ABI hardening

- Removed inert Cargo features and the non-additive
  `limit-to-512-bits` feature. `Natural`'s default identity is fixed at
  `Natural<16>`. The remaining `cli`, `wasm-simd128`, and `fuzzing` features
  each gate real code.
- Added `#[non_exhaustive]` where public enums need room to grow; made parse and
  buffer error details accessor-based; added the requested `#[must_use]`
  coverage and `TryFrom<Natural<P>> for u64`.
- `PrimeFactors` is no longer publicly default-constructible but unpopulatable.
  Its distinct and expanded iterators are named public types, and
  `IntoIterator for &PrimeFactors` supports ordinary `for` loops.
- Release libraries use unwinding, so the C ABI's panic-to-internal-error
  contract is live. Explicit C and Rust worker counts are capped at 256.
- The C ABI is version 2 and adds `rusqsieve_abi_version`,
  `rusqsieve_strerror`, enum-typed statuses, progress/cancellation, and a
  cancellation status. The C module is visible in rustdoc with explicit
  safety contracts.
- Wasm exports that consume caller pointers are `unsafe extern "C"` internally;
  input is copied after checked linear-memory bounds rather than assigned a
  fabricated `'static` lifetime. Unsafe allowances are scoped to the raw Wasm
  boundary and the SIMD kernel.
- Worker relation packets encode only significant root bytes instead of
  padding every root to 128 bytes. This wire change is covered by the Wasm ABI
  bump; no speedup is claimed.

### Tests, CI, and documentation

- Added a deterministic five-case corpus at each browser product, boundary, and
  comparison target—216, 224, 232, 240, 256, and 272 bits—plus a real
  Chromium/Web Worker benchmark driver. Eight-worker SIMD A/B measurements
  retuned five SIQS tiers. Verified corpus means fell by 35.9% at 216 bits
  (5.097 s → 3.266 s), 16.8% at 224 bits (6.417 s → 5.340 s), 28.0% at 232 bits
  (11.258 s → 8.107 s), 7.5% at 240 bits (14.733 s → 13.629 s), and 8.0% at 256
  bits (38.334 s → 35.279 s). The retained `(factor-base bound, half-width,
  threshold adjustment)` settings are `(135k, 131072, 0)`,
  `(150k, 131072, 0)`, `(200k, 131072, −3)`, `(350k, 131072, −1)`, and
  `(400k, 196608, −5)` respectively.
- Rejected several measured non-wins instead of retaining speculative
  optimization code: a dense prime-only hot stream regressed the 240-bit corpus
  mean by 0.5%; 64-way bit-sliced dependency back-substitution left the 272-bit
  LA tail unchanged; and 272-bit interval/factor-base sweeps merely shifted
  time between sieving and LA without a reliable end-to-end gain.
- Added GitHub Actions for musl-native debug and release-profile tests, the
  feature matrix, fmt, clippy with warnings denied, both Wasm build paths, the
  C smoke executable, and scheduled slow 192–256-bit corpus coverage.
- Added an executable documentation audit: every `SPEC §` source citation must
  resolve to a real heading and source comments may not name absent Markdown
  files.
- Extended differential arithmetic tests to square roots, perfect powers,
  decimal formatting/parsing, endian serialization (including tolerated zero
  padding), and shifts at/above word and capacity boundaries.
- Added public-path tests for 2 and 3, a 127-bit prime square, recursive
  three-prime input, a large power of two, and `Natural<17>` routing through
  the optimized engine.
- Added hostile C tests for non-UTF-8, embedded NUL, one million digits,
  `SIZE_MAX` threads, independent concurrent result objects, ABI/status
  helpers, and progress cancellation. `make test` now compiles and runs the C
  smoke program.
- Added cargo-fuzz targets for worker-packet deserialization, decimal parsing,
  and the native C boundary.
- Reconciled the Makefile and release builder so both SIMD artifacts pass
  `-C target-feature=+simd128`. Linux release targets cover both GNU and musl.
- Updated README, SPEC, rustdoc, C header, feature descriptions, supported
  widths, module map, solver/provenance description, ABI exports, panic
  behavior, SIMD build claims, side-channel warnings, and release targets to
  describe the shipped code.

### Montgomery arithmetic

- Added a real fixed-limb Montgomery context for native big-integer
  Pollard–Brent. It computes `−N⁻¹ mod 2⁶⁴`, `R² mod N`, maintains encoded
  residues across the complete rho loop, and performs REDC with word
  multiplication and carry propagation. Randomized differential tests cover
  every significant width from one through sixteen limbs.
- At 256 bits, 100 000 modular squares measured 0.034880 s with
  division-based reduction and 0.011247 s with Montgomery reduction, a 3.10×
  kernel speedup. On the fixed unbalanced 224-bit rho input, 200 fresh CLI
  invocations measured 2.33 s before and 1.33–1.34 s after.
- The rho iteration budget was intentionally not raised. Interleaved fixed
  balanced 224-bit runs measured 5.55–5.58 s before and 5.55–5.62 s after,
  which is a wash; balanced RSA-like inputs do not pay for a deeper rho search.

### Deliberately not implemented

- ECM is not included. Real Montgomery reduction now exists, but ECM would
  still be unsuccessful overhead on balanced RSA-like semiprimes unless it is
  separately budgeted, measured on unbalanced inputs, and disabled by default
  for the balanced artifact. Pollard p−1 was likewise not inserted as
  unconditional overhead. The bounded Pollard–Brent-to-SIQS ladder remains.
- No optional serde dependency was added. The active worker wire is private,
  versioned, fuzzed, and substantially smaller; adding a general serialization
  contract was only a “consider” item and would expand the public format
  surface without serving the current raw-Wasm protocol.
- No SIQS speedup is claimed for 0.3.0. The measured Montgomery improvement is
  confined to native big-integer rho; balanced end-to-end time is unchanged.
  The extensive 0.2.1 interleaved measurements and rejected changes remain
  recorded below.

## 0.2.1 — 2026-07-25

Patch release: correctness fixes and internal-only performance/resource work.
Every change is behind a private item or a function body; no public Rust
signature, C ABI symbol, or Wasm export changed, and no Cargo feature was added
or removed.

### Correctness (behaviour changes)

- **Balanced semiprimes from 65 to 85 bits now factor.** They failed
  deterministically before, including `18446744400127067027 = 4294967311 ×
  4294967357`, and so did any larger input whose cofactor landed in that band.
  `choose_a` drew coefficient factors from factor-base primes above 1000 (10
  bits and up) while accepting only primes within one bit of `ideal_bits`, which
  is 6 there, so the candidate pool was empty for every polynomial family and no
  polynomial was ever built. The lower bound is now derived from `ideal_bits`,
  the acceptance window widens when the pool comes up short, and a
  `debug_assert` fails loudly if the constraint set is ever empty again.
- **A polynomial-coefficient famine is reported as a parameter-selection
  failure, not as "no factor found".** It is detected once in `prepare`, the
  single point both the native and the Wasm/session schedulers pass through, and
  names the factor-base size, the target `A` width, and the candidate count.
  With the fix above reverted locally the minimum reproducer fails in under a
  millisecond instead of sieving 100 000 families for ~1.6 s.
- **`EngineSession::take_jobs` no longer spins forever** when `choose_a`
  declines every family; it shares the native scheduler's 100 000-family bound,
  which was previously two uncoordinated constants. The unbounded loop also grew
  its buffered-result map without limit.
- **Baillie–PSW is the final primality decision above 2^64** — trial division,
  a base-2 Miller–Rabin, then a strong Lucas test with Selfridge Method A
  parameters — in both the library and the browser RSA generator. The 16
  Miller–Rabin rounds remain as an extra layer but no longer carry the
  guarantee: they draw bases from a fixed 32-entry table, so they are
  deterministic rather than probabilistic and a strong pseudoprime to all of
  bases 2..53 would be a guaranteed false positive. Below 2^64 the seven-base
  Jaeschke/Sinclair witness set is proven exact and is unchanged. The
  perfect-square test in front of Selfridge's `D` search is load-bearing, not an
  optimization: that search does not terminate for a square. Measured cost of
  the primality path itself: 0.390 s → 0.439 s per 800 calls above 2^64, **+12.6%**
  — but a factorization makes only a handful of such calls, at `factor_node`
  entry and on recovered factors, never per survivor, so this is 0.007% of a
  256-bit run and no round-count reduction is warranted.
- **Cofactors of 64 bits or less are cancellable.** `smallfactor::factor_u64`
  and its Pollard–Brent loop polled nothing and could not be interrupted.
- **One panicking worker no longer masks the cause.** The job mutex is recovered
  with `unwrap_or_else(|e| e.into_inner())` instead of poisoning every other
  worker with `PoisonError`, and the first panic payload is propagated instead of
  discarded by `let _ = h.join()`.
- **A bounded, cancellable Pollard–Brent stage runs before SIQS at every size.**
  See the note on its budget below.
- Oversized dense linear-algebra fallbacks return a resource-limit error instead
  of silently taking the O(n³) path on the unfiltered matrix, which on Wasm
  permanently inflated the tab's heap.
- **The browser demo could not factor anything.** Two separate faults, both from
  moving the coordinator into its own Worker:
  - `docs/`, which is committed and is what GitHub Pages serves, never received
    `coordinator.js`: the Makefile's published-asset list was hand-maintained and
    was not updated. The page 404s and `boot()` awaits a `ready` message that
    never arrives, so it hangs with the button disabled rather than degrading.
    The list is now derived from `web/` (excluding the local preview server), and
    `make docs-verify` fails if any same-directory reference in the published
    HTML or JS does not resolve — checked against both drift directions.
  - `EngineSession` dropped duplicate-`A` families in `take_jobs`, but the WASM
    coordinator numbers families itself and never calls `take_jobs`. Two families
    that pick the same `A` sieve identical polynomials, so ingesting both puts
    duplicate columns in the matrix and every dependency they form is trivial
    (`x ≡ ±y`) — extraction then reports "no factor" on an input that factors.
    A 110-bit semiprime the native path splits in 14 ms produced 3 duplicate
    families out of 56 and failed outright. The filter now runs where relations
    are ingested, so it protects any scheduler.
  - Added `tools/browser-arch-check.mjs`, which drives the real
    coordinator-Worker plus sieve-Worker protocol on node worker threads and
    asserts a known semiprime comes back correctly factored, and wired it plus
    `docs-verify` into `make test`. Neither fault was reachable by any existing
    test; the architecture was only ever exercised by hand in a browser.

### Tests

- Added the independently verified 309-entry, 7–256-bit factorization corpus at
  `tests/data/`. `cargo test` factors all 280 entries through 128 bits, including
  all 117 in the 65–85-bit regression band, in ~23 s; the 29 larger ones run from
  `supplied_factorization_corpus_above_128_bits` under
  `cargo test --profile release-test -- --ignored`, also ~23 s. All 309 entries
  are verified — product equals `n`, and every listed factor is itself accepted as
  a single prime.
- The dead zone is pinned against the sieve directly, by asserting that
  `prepare` yields a non-empty coefficient pool and that family 0 produces
  polynomials at 65, 70, 75, 80, 85 and 90 bits, and that SIQS alone recovers a
  factor there. An end-to-end `factor()` test would no longer catch a
  regression: the new Pollard–Brent stage splits that whole band first.
- Added a `#[profile.release-test]` profile. `cargo test --release` cannot link
  this crate from a cold cache — `crate-type = ["rlib", "cdylib", "staticlib"]`
  collides on `librusqsieve.rlib` and combines badly with `panic = "abort"` and
  fat LTO — so the slow tiers had no way to run optimized.
- Extended primality coverage to 12 Carmichael numbers and 8 base-2 strong
  pseudoprimes below 2^64, plus, above it, the two smallest strong pseudoprimes
  to all bases through 37 and 41, and a Carmichael number above 2^64. A test
  asserts the strong Lucas step actually executed: no Baillie–PSW pseudoprime is
  known, so every one of these tests would also pass if that step were dead code.
  The counter is thread-local, because `cargo test` runs test functions
  concurrently and a process-wide counter let one test satisfy another's
  assertion.
- Added a test pinning the sieve-root translation, which `find_factor`'s
  duplicate of `prepare` had already diverged on.

### Performance

Measured on a 48-core Xeon Platinum 8259CL @2.50GHz (L1d 32 KiB, L2 1 MiB/core,
L3 71.5 MiB shared), release build, 4 threads, interleaved A/B against the
published 0.2.0 binary, three repetitions, run-to-run spread under 2%. Balanced
corpus semiprimes at 192, 224 and 256 bits.

Sieve and linear algebra, seconds (0.2.0 → 0.2.1): 192-bit 0.70 → 0.68,
224-bit 5.07 → 4.87, 256-bit 37.6 → 36.7. End-to-end wall time is a wash,
because the new Pollard–Brent stage and Baillie–PSW spend 0.5–1.5% of it.
**This is well short of the 10% the remediation brief asked for, and several of
the changes it prescribed were measured to be losses and were not kept** — see
"Measured and rejected".

Kept:

- Sieve-threshold retune. The threshold's large-prime term is now `log2` of the
  large-prime acceptance bound rather than an independent per-tier constant:
  admitting a cofactor the relation collector will reject costs a full trial
  division for nothing. The small-prime give-back is derived from the skipped set
  (`Σ w(p)/(p−1)`, ≈4–6 bits) instead of a hardcoded 8. The remaining offset is
  a measured per-tier value — 0 at 192 bits, −2 at 224, −4 at 256 — replacing a
  single global constant; deeper thresholds trade survivors for polynomials and
  the optimum deepens with size. Worth ~5% at 224 bits over the shipped default.
- Candidate scan. Scores are biased by `128 − threshold` so the scan is one
  masked compare per eight positions. A byte-at-a-time scan had replaced this and
  cost 13.8% of sieve time at 224 bits against 2.5%; note that computing the bias
  inside the scan is not enough, since a runtime bias costs an add and an or per
  word and measured 2.7× slower than the masked test.
- Wrapping rather than saturating score writes, guarded by a bound proved from
  the smallest scored prime instead of assumed. Forcing saturating writes
  everywhere measured +1.5 cpu-s on 12.3 at 224 bits, about 7% of sieve time.
  This is the invariant `RUSQSIEVE_SMALL_SKIP` could previously break silently
  into wrapping overflow and false negatives.
- FLINT's `extra_bits < sieve[i]` stopping rule for trial division, which had
  been dropped. It ends the gated scan once the primes divided out account for
  the recorded score, cutting it to ~64% of the factor base on average.
- Cheap arithmetic: `mod_u64` no longer materializes a discarded 128-byte
  quotient; `mul_mod` drops two of three Knuth divisions and documents plus
  `debug_assert`s the reduced-operand precondition; `knuth_divmod` uses inline
  buffers instead of four heap `Vec`s per call; `WideNatural::rem_natural` and
  `overflowing_mul` no longer round-trip through 256-byte temporaries;
  `EngineJobResult::to_bytes` presizes; `primes_to` uses the cached sieve instead
  of rebuilding 1 229 primes per node; the factor base comes from a segmented
  sieve of Eratosthenes instead of O(√p) trial division per odd number
  (`fb_build` 0.066 s → 0.031 s at 256 bits, and it is paid per Web Worker at
  browser startup).
- Per-family setup: binary shift-subtract xgcd replaces the extended Euclid with
  a hardware divide per iteration.
- The `nvar` cap is `min(9)`, not `min(6)`, so a 256-bit family amortizes its
  setup over all 128 polynomials rather than half of them. Per-family setup fell
  from 19.0 to 10.2 cpu-s at 256 bits, the single largest measured win here.
  `bainv` grows to ~585 KiB at that size; check that against L2 on a target
  before raising the cap further.
- `choose_a` quality: rejection sampling with bounded retry keeps `A` within
  ±25% of target instead of accepting the first draw, `A` values are
  de-duplicated across families, the candidate pool is hoisted into the context
  instead of rescanning the factor base per family, and the residual conversion
  to `u64` is checked rather than truncating to the low limb.
- Large-prime acceptance is `256 × factor_base_bound`, replacing
  `1 << lp_allowance` — which reached 34 360× the factor-base bound at 256 bits,
  against msieve's 100–200. Retained partials at 256 bits fell from 263 024 to
  117 878 and cycle yield rose from 3.7% to 8.2%. The partial-relation forest
  stores relation indices, so re-rooting no longer clones a 16-limb `Natural` and
  a `Vec` along every tree path on every ingest.
- `pollard_u64` is Brent with Montgomery multiplication and a batched GCD every
  128 steps, replacing Floyd with a `u128 %` and a GCD on every iteration.
  **20.8× faster** on a representative double-large-prime cofactor (two 27-bit
  primes: 2.515 s → 0.121 s per 1 000 splits), 18.4× on a 48-bit one. Split
  results are identical on a fixed cofactor corpus.
- `find_factor` calls `prepare` instead of carrying a near-line-for-line copy of
  it. The pinned interval translation is `sieve_half_width % prime`; the
  duplicate's `interval as u32 % prime` was numerically equal only because
  `interval` is that same positive half width, and a test now pins it.
- Browser: relation collection and the GF(2) solve run in a dedicated
  coordinator Worker. Before, `coord.qs_coord_extract` ran on the main thread
  with no `await tick()` after reporting the phase, so the page froze on a stale
  status for the ~8 s serial term — over half of a 48-worker 256-bit run.
- The unreachable second sieve kernel is deleted. `score_polynomial_blocked` was
  gated at 1 MiB of scores while the largest shipped interval produces 640 KiB,
  so no tier could reach it, and forcing it on was slower.

Measured and rejected — each was implemented, measured, and removed:

- **msieve-style resieving** (brief §2.1). Replaying the sparse tail's root
  progressions over a candidate bitmap cut the gated scan from 2.40 to 1.49
  cpu-s at 224 bits and from 26.5 to 12.0 at 256, but the resieve pass itself
  cost 1.59 and 17.1 cpu-s — a net loss at every size and at every cutoff swept
  (interval/2 down to interval/32, and off). The reason is structural, not a
  tuning failure: the multiply-shift gate it competes against is only ~9 cycles
  per prime, and the resieve's cost falls per *polynomial* while its saving
  accrues per *survivor* — this engine sees about 5 survivors per polynomial.
  The brief's ≥10× per-survivor target is not reachable on that ratio.
- **Four Russians / M4RI for the dense residual** (brief §2.11c). Gray-code
  tables over 4-column blocks made `f2_dense` 75% *slower* at 256 bits (2.34 s →
  4.08 s). Pivots are installed incrementally, so every insertion invalidates its
  block's 16-row table and the rebuild cost dominates the XORs saved. A real
  M4RI needs the basis materialized densely first.
- **Twice-log2 sieve weights** (brief §2.3). Finer log resolution requires
  saturating score writes, since a smooth value's score then exceeds a byte, and
  that cost ~7% of sieve time while buying no measurable reduction in survivors.
  Weights are `ceil(log2 p)` with the resolution error absorbed by the tuned
  offset.
- **Raising the tiny-prime skip.** Swept 100/200/400/800/1600/3200 against the
  threshold offset at 224 bits; the flat optimum stays at 100. Each additional
  skipped prime widens the gap between a smooth value's score and its true log,
  so more smooth values fall under the threshold and polynomial count rises
  faster than score-write traffic falls.
- **Retuning the factor base and interval.** The table's own comment claimed the
  optimum above 224 bits was ≈9k primes at 256, against the 20 911 it builds.
  Re-measured at 256 bits (sieve + linear algebra, seconds): 150k → 79.6,
  250k → 50.8, 350k → 39.9, 500k → 35.9 (shipped), 700k → 37.2, 1M → 45.7;
  half-widths 327 680 → 35.9, 458 752 → 37.8, 655 360 → 43.8. The table is
  right and the comment was wrong; the comment is fixed.
- **Raising the number of `A` factors** to amortize per-family setup further.
  At 224 bits, 12-bit factors cut setup from 1.73 to 0.98 cpu-s and 11-bit to
  0.56, but wall time did not improve (4.48 / 4.47 / 4.62 s): `A` quality
  degrades and `bainv` grows.
- **Double large primes** (brief §2.6) remain disabled, now on our own
  measurement rather than an inherited claim. `LargePrime::Two`,
  `classify_cofactor` and cycle combination are all reachable — one constant is
  the whole switch — but at 192 bits, 4 threads (sieve+collect): off 0.605 s
  (5 632 polys, 1 945 cycles); on with the threshold matched to the single-prime
  bound 0.607 s (5 408 polys, 2 040 cycles); on with the threshold 4 bits deeper
  0.680 s; on with the threshold widened to `2 · log2(large-prime bound)`, which
  is what a genuine double actually needs, **did not finish in 300 s**. At a
  matched threshold it buys 4% fewer polynomials and the extra cofactor splits
  eat exactly that. `README.md` and `SPEC.md` no longer advertise it.
- **Per-column allocations in the dense O(n³) solver** are hoisted out of the
  loop (41 816 allocation pairs at 256 bits), moved out only on the two paths
  that consume them. No measurable wall-time change on the filtered path, where
  the dense solver runs on the reduced matrix; kept as a memory-traffic and
  clarity win on the unfiltered fallback.

Not reproduced from the brief's baseline figures: the ≈14% available at 224 bits
from threshold retuning alone. On this host, 0.2.0's shipped default measures
4.72 s and its own optimum 4.47 s — about 5%. The brief's reference host has a
quarter of the cores and half the L3 per thread, and the survivors-versus-
polynomials tradeoff is cache-sensitive.

Also worth recording, since it changes where the remaining time goes: the gated
trial-division scan is memory-bound, not compute-bound. It reads ~9.4 cycles per
factor-base prime while streaming 32 bytes per prime across four separate arrays
(`FactorBaseEntry`, `pinv`, `root1`, `root2`). Packing prime, roots and weight
into one 16-byte record would cut that to two streams and ~24 bytes and is the
largest identified remaining win at 256 bits, where the scan is ~11% of runtime.
It matters more on the small-L2 targets this crate aims at than on this host.

### Bounded Pollard–Brent stage, and a corrected premise

A bounded, cancellable Pollard–Brent stage over `Natural` now runs between the
perfect-power test and SIQS at every size. There was previously no trial division
and no rho above 2^64 on that path, even though a generic `pollard_rho` existed
unused elsewhere in the crate.

The brief expected this to beat SIQS outright from 65 to 100 bits. It does not,
and shipping it that way was a large regression. Measured on balanced corpus
semiprimes (release, single-threaded, seconds — rho with the brief's suggested
16 M-iteration budget, against SIQS alone): 70-bit 0.03 / 0.03, 80-bit **0.56 /
0.03**, 90-bit **2.16 / 0.01**. Rho costs `O(sqrt p)` in the smallest factor
while SIQS at those sizes is already trivial.

So the stage is budgeted rather than gated on bit length, at roughly 1% of the
estimated sieve cost for the input's tier — 1 024 iterations up to 128 bits,
~131 k at 224, ~328 k at 256, measured at 0.7–1.5% of total wall time. That
covers factors to about 2^26 at 192 bits and 2^34 at 256, above the 10^4
trial-division bound and below the 2^64 machine-word path. Where it pays is
unbalanced inputs, which SIQS is worst at: the corpus entry
`13695626177198106295200293487798368178679518660650179392786377544541`
(224 bits, smallest factor 11 624 449) now splits in 0.01 s against about 5 s of
sieving. The budget is the total across all polynomial constants tried; it was
per-constant at first, which made the stage cost 8× its nominal budget — 27 s of
a 64 s 256-bit run — while reporting the same number.

### Documentation and naming

- All eight `SPEC §x.y` citations in `src/` resolved to a wrong section or to
  nothing (`SPEC.md` has sections 1–17; they cited §19.3, §20, §21.1, §15.3,
  §12.5, §12.6, §6.11). All are repointed.
- The comment citing `CLAUDE-AUDIT.md`, a file absent from the crate, is gone
  with the blocked-sieve kernel it documented.
- `SPEC.md` §7.4 and `README.md` described byte scores, resieving and
  double-large-prime combination that the shipped engine does not do; both now
  match the code, including the threshold derivation.
- `qs::parameters`' comment claiming a smaller optimum above 224 bits than the
  table it documents is replaced with the measurements above.
- `smallfactor::pollard_brent` claimed `n` must be odd while handling even `n`.
- `sieve_root_pair` claimed an overflow bound from a hardcoded minimum prime
  that a tuning knob could invalidate; the bound is now computed and the comment
  describes what the code does.
- `legendre_u32` computes a Jacobi symbol, which equals the Legendre symbol only
  for an odd prime modulus. Every caller passes one; the precondition is now
  documented and `debug_assert`ed, and both symbols have direct tests against the
  Euler criterion — neither had any coverage.
- `tonelli_shanks_u32` short-circuits `p ≡ 3 (mod 4)` to `n^((p+1)/4)`, so the
  majority of calls never run the named algorithm. Documented, and its one
  assertion — `assert_eq!(tonelli_shanks_u32(10, 13), Some(7).or(Some(6)))`,
  which `Option::or` collapses to `== Some(7)` — is replaced by an explicit
  membership check plus exhaustive coverage of both branches.
- `filtered_dependencies` said "singleton-row elimination"; it eliminates rows of
  weight 1 through 6 with Markowitz-style pivot selection.
- The 64-dependency cap was justified by "what block solvers conventionally
  return"; this crate has no block solver, so the real reason is stated.
- `PrimalityConfig::rounds` accepts up to 2³²−1 while the witness table has 32
  entries, so rounds beyond 32 repeat witnesses for the cost of a full modexp.
  Documented.
- `ResourceLimitKind::PolynomialBatches`, `ProgressUnit::Polynomials` and
  `FactorConfig`'s doc all attributed SIQS to counters that only ever describe
  the reference `x² − N` path.

Addendum B.5 item 5 asked for the sweep to be extended to `web/*.js` and
`tools/*.mjs`, which it had not covered. Checked every algorithm-naming
identifier there — `pollardBrent`, `millerRabinWitness`, `strongLucasSelfridge`,
`jacobi`, `modPow`, `integerRoot`, `perfectPower`, `trialDivide`, `randomPrime`,
`rsaNumber`, `siqsParallel`, `MR_BASES`, `SMALL_PRIMES`. All implement what they
claim; `pollardBrent` in particular is genuine Brent (`r`-doubling epochs,
batched GCD over runs of ≤128, and the `g == n` backtrack replaying from `ys`),
matching the Rust side. Two comment defects found and fixed: `randomOdd`'s
derivation stated its own bound as `0.75·2^(bits-1)` when forcing the top two
bits gives `0.75·2^bits` (the conclusion it draws is right), and `randomPrime`
still described primality as "the same deterministic-then-strong Miller-Rabin as
isPrime" after `isPrime` became Baillie-PSW above 2^64.

One substantive browser defect came out of that sweep: `factorize` spent a 2^21
Pollard-Brent budget on anything below 84 bits, on the same premise the Rust side
had. Measured in node (BigInt, single-threaded), 2^21 against 2^15: an 80-bit
balanced semiprime 825 ms vs 44 ms, an 85-bit one 724 ms vs 44 ms, while the
unbalanced 127-bit case splits in 0.3 ms either way — and the sieve handles those
sizes in milliseconds. That was up to 825 ms of blocked main thread for work the
workers do faster. The bit-length special case is gone and the budget is a
uniform 2^15.

`Montgomery`, `BlockLanczos`, `extended_gcd`, `prepare_siqs`'s SIQS-claiming
names, `SparseBinaryMatrix::provenance`, `MatrixSolver` and the always-zero job
metrics are untouched here: resolving them removes or renames public items, so
they belong to 0.3.0.

### Attempted or deliberately not attempted, and still open

- **Per-prime `pinv` for the per-family reductions.** Not done. Lemire's fast
  remainder is exact only for a 32-bit dividend, and these reductions are a
  multi-limb `Natural` modulo a 32-bit prime. Precomputing `2^64 mod p` and
  `2^128 mod p` per prime to replace `rem_u64` was worked through and costs about
  the same number of hardware divides, while adding 8 bytes per prime to a loop
  that is already memory-bound. The binary xgcd part of that item did land.
- **`record`'s linear `find`** over the survivor's power list is unchanged. The
  list holds roughly 15 entries, so this is very likely noise; it was not
  measured and is not claimed either way.
- **The browser coordinator Worker (§2.11a) is unmeasured.** The change is
  clearly right — the main thread no longer blocks for the serial solve — but
  there is no browser in this environment, so no before/after number is claimed.
  Its supporting number, the ~8 s serial term, comes from the brief's published
  scaling fit, not from a measurement taken here.
- **No small-L2 re-measurement (§2.12).** The bucket-sieving and blocked-kernel
  negative results are re-confirmed only on a large-L3 host (71.5 MiB shared,
  four threads), which is precisely the configuration the brief warns those
  results are overfit to. The conclusion should be treated as provisional for
  mobile and Wasm targets until someone re-runs it on one.
- **§2.11b is a partial**: the `BTreeSet` filtering structures are now sorted
  `Vec`s, which removes the tree overhead, but not the flat CSR with in-place
  compaction the brief asked for. `f2_filter` at 256 bits moved 0.258 s → 0.225 s.
- **§2.11d is a partial**: back-substitution's parity dot product is unrolled
  four lanes wide but still scalar; it does not use `xor_wasm_simd`. Not measured
  in isolation.
- Acceptance criteria not met, stated plainly: §2.1's ≥10× per-survivor
  reduction (structurally unreachable here — see above), §2.3's ≥10% at 224 and
  256 bits (≈4% and ≈2% of sieve time), and §2.4's order-of-magnitude drop in
  partial-relation memory (2.2×, from 263 024 retained partials to 117 878).

## 0.2.0 — 2026-07-25

- Replaced the broad implementation-facing Rust API with a documented,
  high-level factorization surface.
- Made SIQS factor bases and relations, sparse-matrix kernels, scheduler
  sessions, worker packets, primality policy, and limb mutation crate-private.
- Encapsulated `FactorConfig`; supported tuning now uses builder methods for
  parallelism and progress cadence.
- Replaced `Parallelism::Exact` with the nonzero
  `Parallelism::Threads`/`Parallelism::threads` interface.
- Made progress snapshots read-only and added documented accessors.
- Clarified `PrimeFactors` cardinality with `distinct_len`, `total_len`, and
  `expanded`; removed direct map extraction.
- Added cooperative cancellation to the native SIQS engine.
- Enabled `deny(missing_docs)` for the complete public Rust API.
- Kept the decimal-string C ABI and opaque result ownership model unchanged.
- Replaced the aspirational implementation specification and historical
  performance audit with documentation of the shipped 0.2 architecture,
  verified results, and remaining work.
- Added complete crates.io metadata and an explicit publish allowlist containing
  the Rust/C/Wasm sources, production browser frontend, release tooling,
  licenses, and supporting documentation.
