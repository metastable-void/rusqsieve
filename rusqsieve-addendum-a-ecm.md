# Addendum A — ECM between Pollard-Brent rho and SIQS

**This addendum does not modify the main brief** (`rusqsieve-remediation-prompt.md`). It adds one
optional stage and states the conditions under which it is worth building. Read the main brief
first; §0.1's semver policy and §0.2's verification rules apply here unchanged. This work is a
**`0.3.0`** item: it adds a public tuning surface and changes which algorithm handles some inputs.

The question this answers: should the elliptic curve method sit between the big-integer rho stage
(main brief §1.3) and SIQS (`src/engine.rs:540`), and does it cost anything on the main target of
large balanced RSA-style semiprimes?

---

## A.0 Verdict

**Yes, worth building — but not as an unconditional pre-stage, and not before real Montgomery
multiplication exists.**

Three findings drive the design:

1. **ECM contributes exactly zero factors on the main target.** A balanced 256-bit semiprime has
   two ~128-bit factors, i.e. ~38–39 decimal digits each. Reaching a 35-digit level takes on the
   order of 1 800 curves at `B1 ≈ 10^6`; at roughly 3×10^6 modular multiplications per curve that
   is ~5×10^9 modmuls, which is **several times the cost of simply running SIQS on the same
   input** — and it still would not reach 38–39 digits. Record ECM factors sit in the low 80s of
   digits and were found with very large dedicated efforts. For a balanced RSA semiprime in the
   192–256 bit band, ECM is not a slower route to the answer; it is not a route to the answer.
   (Order-of-magnitude estimates, not measurements — see A.5 for what to measure.)
2. **Therefore any ECM budget spent on the main target is pure overhead.** It must be bounded, and
   bounded *relative to the estimated SIQS cost* rather than as a fixed absolute amount. A fixed
   budget is most damaging exactly where total runtime is smallest: a flat 1–2 s ECM stage would
   more than double the 192-bit browser case, which the project's own benchmarks put at 0.72 s.
3. **ECM is ~100% modular multiplication, and this crate's modular multiplication is currently
   very slow.** `mul_mod` (`src/natural/mod.rs:548-552`) performs **three** Knuth divisions per
   multiply, and `Montgomery` (`src/natural/mod.rs:1098-1150`) is a façade whose `encode`/`decode`
   are the identity and whose `mul` is plain `mul_mod` — there is no REDC. Building ECM on that
   arithmetic would produce a stage that is slow enough to be worthless at any budget.

So: **main brief §3.3's "implement CIOS or delete `Montgomery`" is a hard prerequisite.** If that
lands as a deletion rather than an implementation, do not build ECM — close this addendum and say
so.

---

## A.1 Where ECM actually pays

Not on the headline benchmark. On everything else the crate accepts:

- **Unbalanced `N`** — a moderate factor (roughly 20–35 digits) with a large cofactor. SIQS cost
  depends on the size of `N`, not on the size of its smallest factor, so a 200-bit `N = p·q` with
  a 25-digit `p` is enormously cheaper for ECM than for SIQS. ECM is the only stage in the ladder
  that exploits an unbalanced factorization.
- **`N` with three or more prime factors.** The supplied corpus
  (`rusqsieve-factorization-corpus.txt`) has 93 such entries out of 309, with arity up to 10.
- **Recursive cofactors.** `factor_node` (`src/engine.rs:498-543`) recurses after each split, and
  the ladder it recurses into today is: u64 fast path → primality → perfect power → SIQS. A
  90-bit cofactor carrying a 25-digit factor currently goes to SIQS. This is also the path where
  the main brief's §1.1 dead-zone failures cascaded, so it is worth making robust.
- **General-purpose use.** The public `factor()` accepts arbitrary `N` up to 1024 bits
  (`Natural<16>`). Users will not restrict themselves to balanced semiprimes, and the crate is
  categorized under `algorithms` and `mathematics` on crates.io.

## A.2 Do Pollard p−1 first — it is cheaper and it is nearly free

Before ECM, add **Pollard p−1**. One p−1 run at `B1 = 10^6` costs about the same as a single ECM
curve, and it finds any factor `p` whose `p − 1` is `B1`-smooth. On an RSA-style input `p − 1` is
not smooth by construction, so it will not help the main target — but at one-curve cost that is an
acceptable price for the inputs where it wins outright.

Optionally add Williams p+1 as well; it is the same shape at the same cost and catches a
different smoothness class. Lower priority than p−1.

Recommended full ladder after this addendum, replacing `src/engine.rs:498-543`:

1. small-prime trial division (cheap, currently absent from this path — see main brief §1.3)
2. primality (`is_probable_prime`, with Baillie-PSW per main brief §1.5)
3. perfect power
4. `n < 2^64` → `smallfactor::factor_u64` (deterministic, proven exact)
5. **Pollard p−1**, single run, small `B1`
6. **Pollard-Brent rho**, bounded (main brief §1.3)
7. **ECM**, budgeted per A.3
8. SIQS

## A.3 The budget rule — this is the part that protects the main target

**Do not give ECM a fixed absolute budget, and do not let it run to completion at any digit
level.** Budget it as a fraction of the *estimated* SIQS cost for the same input:

- Estimate the SIQS cost from the existing tier table (`src/qs/mod.rs:330-341`) — factor-base
  size and sieve interval are already size-indexed there, so a cost model calibrated against the
  §0.3 baseline measurements is straightforward.
- Cap the ECM stage at a **small single-digit percentage** of that estimate. Start at 2–3% and
  tune with measurements; make the cap a `FactorConfig` field (per main brief §3.1, **not** an
  environment variable) with a documented default.
- Choose `B1` and the curve count to fit the cap, working up the standard digit levels rather
  than starting deep. Approximate standard parameters: `t20` at `B1 ≈ 11×10^3` with ~75–90 curves;
  `t25` at `B1 ≈ 5×10^4` with ~200–300 curves; `t30` at `B1 ≈ 2.5×10^5` with ~700 curves; `t35` at
  `B1 ≈ 10^6` with ~1 800 curves. Verify these against a current GMP-ECM table rather than
  trusting them from this document.
- **Stop early and hand off to SIQS when the cap is reached.** ECM failing is the expected outcome
  on the main target; it must be a cheap expected outcome, not a long one.

Consequence, and the direct answer to the performance question: with the cap at a few percent, the
main target pays a few percent in the worst case and nothing at all once the cap is set to zero.
Make **zero a supported, documented setting** so anyone benchmarking balanced semiprimes can turn
the stage off outright.

## A.4 Implementation notes

- **Curve form and parameterization.** Montgomery curves with `(X : Z)` coordinates and Suyama's
  parameterization. This is the standard choice because it needs no modular inversion in the
  ladder — only multiplications, squarings, additions — which matters here since inversion is
  expensive on this crate's arithmetic.
- **Stage 1.** PRAC or a standard left-to-right ladder over the primes up to `B1`, with prime
  powers included.
- **Stage 2.** Standard continuation with a baby-step/giant-step table is sufficient and is much
  simpler than Brent-Suyama; measure before reaching for the latter. **Watch memory**: the table
  is the memory high-water mark of this stage, and wasm linear memory never shrinks (main brief
  §2.11e). Size it from a budget, not from `B2` alone.
- **Batch the GCDs.** Do not take a GCD per curve step. Accumulate the product of `Z` values
  modulo `N` and take one GCD per block (~64–128 steps), backing off to a per-step scan only when
  a block hits. This is the same batching the main brief §2.5 asks for in `pollard_u64`, so build
  it once and share it.
- **Stable Rust.** As in main brief §2.2 and §3.3, do not write array lengths computed from a
  const generic parameter — `[u64; 2*P+1]` and `[u64; P+2]` do not compile on stable
  (`generic_const_exprs` is unstable). Size buffers from a literal const and slice them, or carry
  several `[u64; P]` arrays as `WideNatural` (`src/natural/mod.rs:1032-1037`) already does.
- **No new dependencies.** The crate has zero runtime dependencies; keep it that way.
- **Determinism.** Curve selection needs a sigma sequence. Per main brief §3.8 there is no OS
  entropy on `wasm32-unknown-unknown`, so derive sigma from the seeded PRNG, defaulting to seed 0.
  Document that the default configuration therefore tries the same curves every time — which for
  ECM is a correctness-neutral, reproducibility-positive property, unlike the primality case.

## A.5 Parallelism — this is where ECM fits the architecture unusually well

ECM is embarrassingly parallel across curves with **zero inter-worker communication**: each worker
takes a disjoint range of sigma values and reports only success or failure. That is a better fit
for this crate's no-`SharedArrayBuffer` design than SIQS is, and it has **no serial phase at all** —
compare SIQS, where the linear algebra is ~58% of a 48-worker 256-bit run (main brief §2.11).

Two consequences worth stating in the docs:

- On inputs where ECM can win, parallel scaling is near-linear, unlike SIQS.
- The existing coordinator/worker protocol (`src/wasm.rs`, `web/worker.js`) can carry ECM with a
  much simpler message shape than relations — a sigma range in, a factor or nothing out.

Reuse the existing dispatch (`web/index.js:87-92`) rather than adding a second scheduler.

## A.6 Acceptance criteria

1. **The main target does not regress.** Measure 192, 224 and 256-bit balanced semiprimes with the
   ECM stage enabled at its default cap and with it disabled, A/B interleaved on the same host.
   Report both. A regression beyond the documented cap is a failure of the budget model, not an
   acceptable cost.
2. **Setting the ECM budget to zero produces byte-identical behaviour to a build without the
   stage.** Test this — it is the escape hatch benchmarkers need.
3. **ECM actually finds factors it should.** Build a targeted corpus of `N = p·q` where `p` is 20,
   25 and 30 digits and `q` is large, and assert the ECM stage — not SIQS — produces the split.
   Assert via an instrumented counter which stage returned the factor, not merely that the answer
   was right. The main brief's §1.5 acceptance criterion makes the same point for a different
   reason: a stage that silently never runs passes every correctness test.
4. **ECM beats SIQS on those inputs**, measured. If it does not, the budget or the parameters are
   wrong.
5. **Memory ceiling respected** in the browser at the 256-bit tier, with the stage-2 table sized
   from a budget.
6. **The unbalanced and ≥3-factor entries of the supplied corpus** are factored correctly, and the
   overall corpus wall time does not increase.

## A.7 When to abandon this addendum

Close it and report that ECM is not worth building if any of these hold:

- Main brief §3.3 resolves by **deleting** `Montgomery` rather than implementing CIOS. ECM on
  division-based modular arithmetic is not worth having.
- The measured cost of a `t20`-level stage cannot be brought under a few percent of SIQS time at
  the 192-bit tier — the smallest, most sensitive case.
- Criterion A.6.4 fails: ECM does not actually beat SIQS on deliberately unbalanced inputs.

Any of those means the ladder is better off going straight from rho to SIQS.
