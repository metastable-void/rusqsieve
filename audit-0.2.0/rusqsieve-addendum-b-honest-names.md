# Addendum B — Names that lie: standing rule and full inventory

**This addendum does not modify the main brief** (`rusqsieve-remediation-prompt.md`). It supersedes
the main brief's §3.3 subsection "Named for something it does not do" by **widening it into a rule
and completing the list**. Where this addendum and §3.3 disagree on scope, follow this addendum.
§0.1's semver policy and §0.2's verification rules apply unchanged.

Why this exists: the main brief named five such items, but that was a by-product of reading the
code, not the output of a deliberate sweep, and it gave you a list rather than a rule. A
subsequent systematic sweep enumerated ~390 declared items across all 17 source files, identified
**97 names that reference a named algorithm, a person's method, a named data structure, or a
specific mathematical object**, and read every one of their implementations. It found more.

---

## B.0 The standing rule

**Every identifier that names an algorithm, a person's method, a named data structure, or a
specific mathematical object must be verified against its implementation. Where they disagree,
either implement the named thing or rename/remove the identifier. There is no third state.**

This applies to types, traits, functions, methods, enum variants, modules, consts, and macros —
and to the doc comments and code comments attached to them, including attributions to FLINT,
msieve, YAFU, or a named author.

Apply the rule as a task in its own right, not only to the items listed below. When you add new
code, apply it to that too. Specifically:

- An identifier naming an algorithm the code does not run is a defect of the same class as a wrong
  answer, because it causes readers and future maintainers to reason about performance and
  correctness properties the code does not have.
- A doc comment that describes behaviour the function does not have is the same defect.
- An attribution (`SPEC §x.y`, "FLINT's ...", "per msieve") that does not resolve, or resolves to
  something else, is the same defect at lower severity.
- **A test that passes only because the named algorithm is absent is the most severe form**, because
  it converts the lie into a maintained invariant. See B.1.1.

For each item you resolve, record in `CHANGELOG.md` which of the two paths you took — implemented,
or renamed/removed — and its semver class. Removing or renaming a `pub` item is a `0.3.0` change.

---

## B.1 Outright lies — named for an algorithm they do not implement

Items 1–5 are already in main brief §3.3 and are restated here only where the sweep added detail.
Items 6–8 are **new** and are not in the main brief at all.

### B.1.1 `Montgomery<P>` — `src/natural/mod.rs:1098-1150` *(in §3.3; new detail below)*

Already covered, with one addition that changes how you must sequence the fix:

- `encode` (`:1111`) is `v.div_rem(&modulus).1`, `decode` (`:1114`) literally calls `encode`, and
  `mul` (`:1117`) is `mul_mod`, which routes through `knuth_divmod` (`:1052,1066`) — full long
  division, the exact cost Montgomery exists to remove.
- `MontgomeryError::EvenModulus` (`:1085`, enforced at `:1105`) is a restriction **nothing in the
  implementation needs** — plain `mul_mod` works fine for even moduli. It exists only to make the
  type resemble real Montgomery. Meanwhile `Montgomery::inv` is exercised with the *composite*
  modulus 15 at `:1289`.
- **The critical point: a correct Montgomery implementation would FAIL the existing test.**
  `diff_montgomery` (`:1439-1466`) asserts at `:1455-1459` that `mont.mul(a, b) == (a * b) % m`
  with `a` and `b` **not encoded**. In the Montgomery domain the product is `a·b·R⁻¹ mod m`, which
  is not `a·b mod m`. So this 2000-iteration differential test does not merely fail to detect the
  missing algorithm — **it pins the API to non-Montgomery semantics and will reject the fix.**
  Rewrite or delete that test as the *first* step of implementing CIOS, and say so in the
  changelog. If you leave it, you will conclude your correct implementation is broken.

### B.1.2 `BlockLanczos` and its protocol — `src/f2/mod.rs:475-520` *(in §3.3)*

Restated with the sweep's detail: `begin()` (`:501`) calls `filtered_dependencies()` and sets
`complete: true` in the same constructor (`:504`); `request()` (`:507`) returns `Complete`
unconditionally, so `LanczosRequest::{MultiplyM, MultiplyMt}` (`:480-481`) are never constructed
anywhere; `submit_product(&mut self, _: &[u64])` (`:510`) discards its parameter;
`LanczosProgress::Progressed` (`:487`) and `LinearAlgebraError::WrongProductLength` (`:491`) are
unreachable. No test names `BlockLanczos`, `LanczosRequest`, or `submit_product`.

`SPEC.md:385-387` and `README.md:214` both explicitly disclaim the recurrence. **The prose is
honest; the identifiers are not.** That is precisely the gap this addendum closes: a disclaimer in
the docs does not license a lying identifier.

### B.1.3 `F2BlockVector` — `src/f2/mod.rs:458-473` *(in §3.3's deletion list)*

`Box<[u64]>` with `new`/`as_slice`/`as_mut_slice`. No block width, no F2 operation. Never
constructed, never a parameter, never a return type — it exists only to make the Lanczos API look
populated.

### B.1.4 `extended_gcd` / `ExtendedGcdResult` — `src/natural/mod.rs:457-459`, `:1075-1079` *(in §3.3; new detail)*

The body is `ExtendedGcdResult { gcd: self.gcd(rhs) }` — the plain binary GCD, wrapped.

New: **the doc comment at `:1075` is itself false.** It says "Coefficients are intentionally
private pending a public signed type." There are no private coefficients; the struct has exactly
one field. A reader is told the Bézout data exists and is merely unexported. Fix the comment as
part of whichever path you take.

### B.1.5 `prepare_siqs` — `src/qs/mod.rs:245-288` — **NEW**

Named for *self-initializing* quadratic sieve setup. The defining SIQS artifacts are the
`a = ∏ qᵢ` coefficient, the `Bⱼ` values, and O(1) Gray-code root updates. `prepare_siqs` produces
**none** of them: it builds a factor base and computes an FNV hash of `n` as `context_id`. The only
polynomial state it emits is `PolynomialPlan { first_x_offset: 0 }` (`:284`), a one-field struct
hardcoded to 0 and never read.

Its consumer confirms it: `sieve_job` (`:495`) sieves `Q(x) = x² − N` at
`x = ceil_sqrt(N) + offset` (`:508`, `:543-551`) — **plain single-polynomial quadratic sieve**, not
even MPQS. The function's own doc at `:487-494` correctly calls these "segments", but the
surrounding names all call segments *polynomials*: `polynomial_count`, `first_polynomial`,
`polynomial_families: 1` (`:512`), `RelationSource.polynomial` (`:648`).

Real SIQS **does** exist in this crate — `sieve_family` (`src/engine.rs:750-892`) genuinely
implements `a`-selection, signed B-values, and Gray-code root advance. It simply does not go
through `prepare_siqs`; `engine::prepare` (`src/engine.rs:270`) calls `prepare_siqs` only to
obtain a factor base and discards the rest.

**Task:** rename this module's SIQS-claiming identifiers to say what they are (segmented
single-polynomial QS), or delete the path per main brief §3.2's decision on `P > PARTS`. Rename
`polynomial*` to `segment*` throughout `src/qs/mod.rs` where it refers to sieve segments. Note
`sieve_job_logarithmic_sieve` (`:669-706`) validates it as a working *log sieve*, which it is —
that test is orthogonal, not vacuous, and should be kept.

### B.1.6 `SparseBinaryMatrix::provenance` — `src/f2/mod.rs:43, 92, 104` — **NEW**

Named for provenance tracking: which original relation columns a combined column came from.
`SPEC.md:379` lists "exact provenance tracking back to original relation columns" as property 2 of
the 0.2 solver.

It is set to the identity permutation `(0..c)` at construction (`:92`) and **never mutated
anywhere**. `pub fn provenance()` (`:104`) always returns `[0, 1, 2, …, columns-1]`.

The contradiction is inside the same file: `row_echelon_dependencies`'s doc (`:186`) and
`filtered_dependencies`'s doc (`:282`) both state that dependencies are recovered *without*
maintaining provenance during elimination. So the code documents that it does not do the thing the
field is named for, while the field and `SPEC.md` both claim it does. No test.

To be precise about severity: this is **not** a correctness bug. `dense_dependencies` (`:146`)
carries its own provenance vectors internally and every dependency is re-verified against the
original matrix (`:162`, `:387`; `src/engine.rs:414`). The defect is the false public accessor and
the false SPEC claim.

**Task:** delete the field and the accessor and correct `SPEC.md:379`, or implement real
provenance tracking through the elimination. Deleting a `pub fn` is a `0.3.0` change.

### B.1.7 `MatrixSolver` + `MatrixConfig` — `src/f2/mod.rs:12-31` — **NEW as a naming defect**

Main brief §3.3 lists these under unreachable/never-read config. The sweep establishes something
sharper: **`MatrixSolver` is a solver-selection API that selects nothing.**

`MatrixConfig` is stored in `QsConfig::matrix` (`src/qs/mod.rs:38,58`) and never read. Per field:
`solver` — never constructed as `DenseGaussian` or `BlockLanczos`, never matched on;
`dense_threshold: 512` (`:20,27`) — never read; `structured_elimination_limit: 10_000` (`:21,28`) —
never read, while `filtered_dependencies` hardcodes `MAX_STRUCTURED_WEIGHT = 6` (`:284`) and
`row_echelon_dependencies(64)` (`:372`).

**A caller who sets `solver: MatrixSolver::BlockLanczos` silently gets the same code path as
`Auto`.** That is worse than a dead field: it is a public switch that appears to choose an
algorithm and does not. Delete the enum and the config, or wire them up.

### B.1.8 `MatrixJobMetrics::nonzeros_visited` — `src/work/mod.rs:62` — **NEW**

Hardcoded to `0` at `:153` despite the name. Also `SieveJobMetrics::double_large_prime_relations`
(`:51`) is permanently 0 because `qs::sieve_job` never captures doubles (it `continue`s at `:637`).
A metric that always reports zero reads as "this never happens" rather than "this is not measured".
Remove them or populate them.

---

## B.2 Stubs — real intended behaviour, no-op body

Main brief §3.3 covers `FactorSession`, the seven wasm exports, and `execute_job` under
"Unreachable". The sweep adds detail worth having, and the standing rule in B.0 applies to them as
naming defects too — a "session" that runs to completion in one call is misnamed, not merely
unreachable.

- **`FactorSession`** (`src/factor.rs:223-308`): `advance_local` (`:257`) ignores its
  `LocalWorkBudget` and calls `factor_complete` (`:267`); `LocalWorkBudget`'s five fields
  (`:193-199`) are never read. `take_jobs` (`:286`) returns `Ok(Vec::new())` always, which is why
  the entire `work::WorkJob` protocol is unreachable. `submit` (`:289`) reads only
  `header.generation` and discards the relations. `self.generation` is initialized to 0 (`:245`) and
  never incremented, so `IgnoredObsolete` is reachable only by a caller forging a nonzero
  generation. `SessionPhase` (`:179-191`) has 11 variants of which 3 are ever assigned;
  `AdvanceOutcome::{Progressed, NeedsWorkers}` (`:213-215`) are never returned.
- **The seven wasm exports** (`src/wasm.rs:186,190,194,227,244,248,250`): each ignores its
  arguments and returns an empty packet, `0`, or `-1`. Mitigating: `SPEC.md:491-517`'s "Active raw
  ABI" lists none of them, and `web/` and `tools/` never call them — so they are undisclosed dead
  exports rather than advertised lies. They still must go, and removing exports from a versioned
  ABI (`ABI_VERSION = 1`, `src/wasm.rs:7`) requires bumping it.
  Separately, `qs_session_advance_local` (`:175`) discards the budget pointer and length and
  substitutes `LocalWorkBudget::default()`, then runs the entire factorization synchronously — an
  "advance" that does not return until done.
- **`execute_job`** (`src/work/mod.rs:110-157`): the dispatch body is real and correct, but the
  function is never called. `WorkerScratch.matrix` / `.arithmetic` (`:80-84`) are never used — the
  matrix branch allocates a fresh `vec![0; end - start]` at `:136` instead of reusing
  `scratch.matrix.words`.

---

## B.3 Misleading but not false — fix the text, not the code

### B.3.1 Every `SPEC §x.y` citation in the source is stale — **NEW, systematic**

`SPEC.md` contains sections **1 through 17 only**. All eight citations in the source resolve to
the wrong section or to nothing:

| Citation | Site | Reality |
|---|---|---|
| `SPEC §21.1` | `src/engine.rs:120` | no §18+ exists |
| `SPEC §20` | `src/smallfactor.rs:6` | no §20 exists |
| `SPEC §19.3` | `src/natural/mod.rs:1299` | no §19 exists |
| `SPEC §15.3` | `src/f2/mod.rs:248` | §15 is "Safety and invariants"; the content is at §8 |
| `SPEC §12.5` | `src/engine.rs:132` | §12 is "CLI"; the content is at §7.3 |
| `SPEC §12.6` | `src/engine.rs:1269`, `src/qs/mod.rs:486` | §12 is "CLI"; the content is at §7.4 |
| `SPEC §6.11` | `src/natural/mod.rs:596` | §6 has no numbered subsections |

Fix all eight, and add a CI check (main brief §4.1) that every `SPEC §` reference in `src/`
resolves to a heading that exists in `SPEC.md`. Without the check they will drift again.

### B.3.2 A dangling reference to a file that does not exist — **NEW**

`src/engine.rs:1268-1273` cites `CLAUDE-AUDIT.md` for the blocked-sieving measurement. **That file
is not in the crate.** The comment also reads as if blocked sieving had been rejected, while nine
lines later (`:1275-1284`) the code does dispatch to `score_polynomial_blocked` above `BLOCK_GATE`.
Either commit the audit file or inline its conclusion, and reconcile the comment with main brief
§2.12's finding that the gate is unreachable at every shipped tier.

### B.3.3 `legendre_u32` computes a Jacobi symbol — `src/natural/mod.rs:1174-1179` — **NEW**

Named for the Legendre symbol; the body delegates to `jacobi_u64` (`:1178`), which is the Jacobi
symbol — equal to Legendre only when the modulus is an odd prime. It also special-cases `p == 2`
(`:1175`) returning `n & 1`, for which the Legendre symbol is undefined. Its only caller
(`src/qs/mod.rs:151`) does pass primes, so it is correct in practice. Rename it, or document the
precondition and `debug_assert` it. Neither `legendre_u32` nor `jacobi_u64` has a direct test —
add one (main brief §4.3).

### B.3.4 `tonelli_shanks_u32` short-circuits the named algorithm — `src/natural/mod.rs:1180-1225` — **NEW**

The main loop (`:1194-1223`) is genuine Tonelli-Shanks, but `:1191-1193` short-circuits
`p ≡ 3 (mod 4)` to `n^((p+1)/4)`, which is not Tonelli-Shanks. Roughly half the factor base is
such primes, so the majority of calls never execute the named algorithm. This is common practice
and not a defect — but document it, because the reader's cost model is wrong otherwise.

**And fix the test, which is a broken idiom.** `src/natural/mod.rs:1267`:

```rust
assert_eq!(tonelli_shanks_u32(10, 13), Some(7).or(Some(6)));
```

`Option::or` on a `Some` returns the receiver, so this means exactly `== Some(7)`. The author
evidently meant "7 or 6". It happens to pass, so it is not vacuous — but it is the **only**
assertion covering Tonelli-Shanks in the crate, and it will mislead the next editor. Write it as an
explicit membership check against the set of valid roots.

### B.3.5 `filtered_dependencies`'s doc understates itself — `src/f2/mod.rs:248-256` — **NEW**

The doc says "iterative singleton-row elimination". The code eliminates every row of weight
**1 through 6** (`MAX_STRUCTURED_WEIGHT = 6` at `:284,286,291,319`) with Markowitz-style pivot
selection (`:295-300`). `SPEC.md:376-377` describes this correctly; only the code doc is wrong. It
under-claims rather than over-claims, but it will mislead anyone reasoning about fill-in.

### B.3.6 A rationale appealing to a solver the crate does not contain — `src/f2/mod.rs:369-371` — **NEW**

"Block solvers conventionally return up to 64 independent QS dependencies" justifies the hardcoded
`64` by appeal to a class of solver this crate does not have. State the real reason.

### B.3.7 `PrimalityConfig::rounds` overstates — `src/primality.rs:7, 51` — **NEW**

`rounds: NonZero<u32>` accepts up to 2³²−1, but with the default `WitnessPolicy::FirstPrimes` the
witness is `SMALL[round % SMALL.len()]` (`:51`) over a 32-entry table (`:24`). **Rounds 33 and
beyond re-run identical witnesses**, adding exactly zero confidence at the cost of full modexps.
Cap it, or document the ceiling. This interacts with main brief §1.5 (Baillie-PSW) and §3.8
(seeded witnesses) — after those land, the fixed table is no longer what carries the strength.

### B.3.8 SIQS attribution on non-SIQS counters — **NEW**

`ResourceLimitKind::PolynomialBatches` (`src/factor.rs:121`) is documented "Maximum SIQS polynomial
batches", but bounds segments of the plain `x² − N` sieve (B.1.5). Same at
`src/qs/mod.rs:512`; `FactorConfig`'s doc (`src/factor.rs:46-47`) says releases "can improve SIQS
parameters" when the `QsConfig` it wraps drives only the non-SIQS reference path; and
`ProgressUnit::Polynomials` (`src/progress.rs:110`, "SIQS polynomials processed") is never emitted.

### B.3.9 `score_polynomial_blocked`'s comment mildly overstates — `src/engine.rs:1096-1099` — **NEW**

"No modular arithmetic or re-striding setup is repeated." Carried positions genuinely are preserved
across blocks (`pos1`/`pos2` at `:1163,1173`), so the substantive claim holds — but the per-block
loop does re-walk the full factor base and recompute `weight = 32 - p.leading_zeros()` for every
prime in every block (`:1147-1153`).

### B.3.10 A stale precondition — `src/smallfactor.rs:111` — **NEW**

Says "`n` must be an odd composite" while `:114` handles even `n`. Fix the comment.

---

## B.4 Explicitly checked and HONEST — do not "fix" these

The sweep specifically investigated two suspicions and **refuted both**. Do not change them.

- **`pollard_u64` (`src/engine.rs:1503-1531`) is honestly documented.** Its doc says "Pollard's rho
  (**Floyd**)" and the body is Floyd (`:1521-1522`). Nothing anywhere in `src/`, `SPEC.md`,
  `README.md`, `BENCHMARKING.md`, `CHANGELOG.md`, `web/`, or `tools/` claims Brent for it. The
  separate claim at `src/engine.rs:517-519` that inputs under 64 bits use "Pollard-Brent" **is
  true**: that path calls `smallfactor::factor_u64` (`:522`), and `src/smallfactor.rs:113-161` is
  **genuine Brent** — `r`-doubling epochs (`:128,145`), batched GCD over runs of ≤128 (`:134-142`),
  and the `g == n` backtrack replaying from `ys` (`:147-155`). Main brief §2.5's criticism of
  `pollard_u64` stands on **performance** grounds (Floyd is weaker than Brent, and it takes a
  `u128 %` and a GCD every iteration), not on naming. Fix it for speed; do not "correct" a name
  that is accurate.
- **`Forest` (`src/engine.rs:1533-1588`) is honestly named.** Its doc says "a spanning forest over
  large-prime vertices", which is what it is. Nothing in the crate claims union-find, disjoint-set,
  or path compression. `root()` (`:1554`) deliberately walks parents **without** compression, and it
  must: `path()` (`:1560`) reads the edge stored at each node on the way up in order to reconstruct
  the cycle's relations, so compressing the path would destroy the data. Main brief §2.4's
  criticism stands on the `Relation` **cloning** along tree paths, not on the absence of path
  compression. Do not add path compression.

Also verified honest and worth knowing, so you do not disturb them while working nearby:
`knuth_divmod` (`src/natural/mod.rs:597`) is real TAOCP Algorithm D including normalization, the
two-limb `qhat` estimate with correction, and multiply-subtract with add-back;
`Natural::gcd` (`:434`) is real binary/Stein; `lemire_c`/`fastmod` (`src/engine.rs:1656,1663`) is
genuine Lemire fast remainder and is **well tested and non-vacuous** (`:1783-1795` checks all
`p ∈ [2, 10000]`); `knuth_schroeppel` (`:1677`) is a faithful FLINT port with FLINT's 29-entry
multiplier table verbatim and an accurate attribution; `sieve_family` (`:750-892`) is **real
SIQS**; `smallfactor::sieve_primes` (`:15`) is a real sieve of Eratosthenes; `Registry<T>`
(`src/wasm.rs:14`) is a real generational-index slot map whose generation guard genuinely prevents
stale-handle aliasing; and all of `src/capi.rs` matches `rusqsieve.h`.

---

## B.5 Acceptance criteria

1. Every item in B.1, B.2 and B.3 is resolved by **implementing** or by **renaming/removing** —
   with the choice and its semver class recorded in `CHANGELOG.md`.
2. **`diff_montgomery` is rewritten or deleted before CIOS is implemented** (B.1.1). State in the
   changelog that the old test asserted non-Montgomery semantics.
3. A CI check exists that every `SPEC §` reference in `src/` resolves to a heading present in
   `SPEC.md` (B.3.1), and that no source comment references a file absent from the repository
   (B.3.2).
4. Nothing in B.4 was changed.
5. You report the result of applying the B.0 rule yourself — the identifiers you checked beyond
   this list, and anything new you found. The sweep behind this addendum covered `src/` only; it
   did **not** read the bodies of `web/*.js` or `tools/*.mjs`, where `web/numtheory.js:85`
   (`pollardBrent`) and the surrounding preprocessing path are unverified. Cover them.
