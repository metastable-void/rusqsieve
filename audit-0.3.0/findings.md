# Detailed findings

Severity meanings:

- **Critical**: release is broken or a central advertised capability is absent.
- **High**: realistic correctness, availability, or release-integrity risk.
- **Medium**: material quality, maintainability, or performance-evidence gap.
- **Low**: contained hardening or cleanup issue.

## Medium packaging finding

### A030-001 — The additional Wasm demo glue omits the coordinator Worker

**Status: remediated after the audit.** `build-release.sh` now includes
`coordinator.js`; a rebuilt archive was extracted and checked for all referenced
frontend files.

**Evidence**

- `web/index.js:46` constructs `./coordinator.js`.
- `build-release.sh:246-248` copies a hand-maintained list that contains
  `worker.js` but not `coordinator.js`.
- `rusqsieve-0.3.0-wasm32-unknown-unknown.tar.gz` contains no
  `web/coordinator.js`.
- `Makefile:24-32` explicitly describes the same historical omission and
  avoids it for `docs/` by deriving the asset list, but the release builder
  retained a second manual list.

**Scope and impact**

The primary artifacts in this archive are the scalar and SIMD raw Wasm modules;
both are present and passed the audit's build and execution checks. The
JavaScript files are additional/informational glue rather than the normative
artifact.

That bundled frontend is not self-contained: if deployed as supplied, its
module-Worker request for `coordinator.js` fails and the page never becomes
ready. This is therefore a defect in the optional demo/integration material,
not a defect in the raw Wasm modules.

**Remediation**

Choose and document one contract:

- if the glue is informational, label it explicitly and avoid describing the
  archive as a directly deployable frontend;
- if it is intended to be a self-contained demo, derive its asset list from
  `web/` just as the Makefile does.

For the latter, extract the archive, run the reference audit, serve `web/`, and
execute the browser architecture test against it.

## High

### A030-002 — The opening mobile-performance claim has no supporting evidence

`README.md:8` claims the fastest publicly documented browser SIQS timings "on
a consumer mobile device." The recorded measurements in `BENCHMARKING.md:92-96`
are Node 24.15/V8 on the development host; that host is a dual-socket,
96-thread Xeon Platinum 8259CL. `SPEC.md:728` and `CLAUDE-AUDIT.md:141-142,
177-178` say representative mobile cache measurements are still outstanding.
No phone model, browser version, thermal state, power mode, raw mobile results,
or competitor run is recorded.

This is both an evidence and release-trust problem. Replace it with the narrower
claim already used later in the README—"plausibly fastest-class" on the
documented host—or publish reproducible on-device and same-device competitor
measurements.

### A030-003 — Stale worker errors can fail a later factorization

**Status: remediated after the audit.** All run responses, including errors,
must now match the active generation.

`web/index.js:91` exempts coordinator errors from the generation check.
`web/index.js:133-138` likewise handles every sieve-worker error before checking
`data.gen`. A completed run leaves other jobs in flight. If one of those old
jobs reports an error after the user starts another run, the old error rejects
the new run.

Every per-run message, including errors, must carry and match the active
generation. Initialization errors should use a separate boot lifecycle. Add a
two-run regression test with delayed success/error messages from the first
generation.

### A030-004 — Family-budget exhaustion can leave the browser promise pending forever

**Status: remediated after the audit.** The scheduler tracks active jobs and
pending submissions and rejects once the bounded family domain is drained.

`web/index.js:84-89` silently returns from `dispatch` when the next family
exceeds `MAX_FAMILIES`. Once all workers take that path, no worker is active and
the SIQS promise is neither resolved nor rejected. The rejection at
`web/index.js:156-159` only runs when a worker returns a null payload, not when
normal results consume the budget.

The architecture test does not reproduce production here:
`tools/browser-arch-check.mjs:61-66` explicitly rejects on exhaustion.

Track active jobs and dispatched families centrally. Reject with a typed
`InsufficientRelations` error when the budget is exhausted and no jobs remain.
Add a tiny injected budget to make the test fast and deterministic.

### A030-005 — Scheduler limits and documented error semantics disagree

**Status: partially remediated after the audit.** JavaScript now uses the
engine's 100,000-family bound. The native exhaustion path still needs to return
`InsufficientRelations` explicitly.

- Rust defines `MAX_FAMILIES = 100_000` at `src/engine.rs:251`.
- The browser defines `MAX_FAMILIES = 2_000_000` at `web/index.js:12`.
- `EngineSession::submit_bytes` accepts caller-numbered families without
  enforcing the Rust budget.
- The comment at `src/engine.rs:246-250` says the native and session schedulers
  use one bound and surface `EngineError::InsufficientRelations`.
- The native loop exits at the limit (`src/engine.rs:694`) and immediately
  proceeds to linear algebra (`src/engine.rs:775-802`); it never constructs
  `EngineError::InsufficientRelations`.

The public `FactorError::InsufficientRelations` mapping therefore does not
represent the native exhaustion path, and browser/native behavior can diverge
by 20×.

Put the budget and exhaustion state in `EngineSession`, expose bounded job
allocation to Wasm, and explicitly return `InsufficientRelations` before
linear algebra when the relation target was not reached.

### A030-006 — Linear algebra is a destructive one-shot operation in the browser

**Status: remediated for the reference browser glue after the audit.** A
no-factor extraction retains the session and requests at least 64 additional
relations before retrying.

At the first relation target, `web/coordinator.js:45-53` runs extraction,
immediately frees the session, and throws if no factor is returned. A relation
target is a heuristic surplus, not a proof that one of the bounded dependencies
will yield a nontrivial GCD. In-flight relations are also discarded.

The native path has the same one-shot policy, but the browser is especially
well positioned to recover because more family results are already in flight.
On `NoFactor`, retain the session, raise the target by a bounded surplus,
ingest outstanding results, and resume dispatch. Distinguish invalid
dependencies, resource limits, and "need more relations" in the Wasm ABI.

### A030-007 — Worker initialization can hang boot indefinitely

**Status: remediated after the audit.** Boot handles Worker errors,
message errors, cancellation, and a 30-second timeout; failed runs reset the
entire Worker runtime.

The worker-ready promise at `web/index.js:55-65` only resolves on `ready`.
It does not reject on a worker `error` message, a Worker `error` event, or a
timeout. A failed Wasm instantiation therefore leaves `Promise.all` and the UI
pending forever.

Use a shared request helper that handles expected messages, protocol errors,
Worker errors/message errors, termination, and a timeout. Terminate the partly
created pool if boot fails.

### A030-008 — The browser accepts inputs outside the documented 512-bit range

**Status: remediated in the reference demo after the audit.** The UI rejects
inputs wider than 512 significant bits before creating a Wasm session.

Native entry points reject inputs above 512 bits in `src/native.rs:79-82`.
The Wasm path parses into the default 1024-bit `Natural` and calls
`engine::prepare` without a bit-length check (`src/wasm.rs:155-169, 200-209`).
The UI accepts arbitrary decimal text. A 513–1024-bit hard composite is thus
accepted into a parameter regime that is undocumented and impractical for
SIQS, potentially tying up dozens of workers until the broken JavaScript
budget path hangs.

Enforce one explicit browser limit before context creation and in the UI.
Return a typed range error rather than handle `0`.

## Additional medium findings

### A030-009 — CI does not test the shipped browser product

**Status: partially remediated after the audit.** CI now runs the Worker
architecture smoke test and focused glue failure checks. Extracted-archive and
real-browser CI remain outstanding.

The Wasm job compiles modules and runs `build-release.sh`, but does not inspect
or execute the resulting archive. It therefore passed while omitting a required
file. The Node architecture check is not run in CI, and no real-browser test is
present.

The CI file also omits several checks declared mandatory in `SPEC.md:680-689`:
strict rustdoc, all-target/all-feature check, and `cargo package`. The current
workflow has no musl job despite musl being a supported release target.
`actions/checkout@v4` and the Rust toolchain are not immutable.

Add:

- archive extraction plus reference/import audit;
- Node architecture and one small Playwright/Web Worker smoke test;
- scalar and SIMD runs, not build-only checks;
- strict rustdoc, `cargo check`, and `cargo package`;
- native musl compile/test or at least archive execution;
- pinned action SHAs and a pinned MSRV/current toolchain policy.

### A030-010 — Generated fuzz build output is committed

`fuzz/target` contains 162 tracked files totaling 83,880,588 bytes (80.0 MiB
from `git ls-tree`). The root `.gitignore` only ignores `/target/`, not nested
target directories. Running `cargo check --manifest-path fuzz/Cargo.toml`
rewrites and deletes tracked compiler artifacts.

Remove the directory from version control and ignore `**/target/` or
`fuzz/target/`. Keep only fuzz manifests, sources, seed corpora, and intentional
crash artifacts.

### A030-011 — The browser packet reader does not validate its envelope

**Status: remediated after the audit.** The glue validates QSV1 magic, version,
kind, exact payload length, memory bounds, and relation-batch framing, while
freeing owned handles on every path.

`web/abi.js:32-46` does not check `QSV1` magic, kind, version, nonzero
pointer/length consistency, or that `payloadLen <= len - 12`. `Uint8Array.slice`
silently truncates an oversized declared payload. Rust deserialization also
accepts trailing bytes and silently ignores malformed submitted families.

Validate the complete envelope and expected packet kind, require exact lengths,
and return typed protocol errors. Add mutation tests around every length and
tag field.

### A030-012 — Browser performance evidence is strong but not portable enough

The tier table is based on five fixed balanced inputs, which is good, but most
optimization decisions and public timings come from one large-cache Xeon.
There is no retained raw-result artifact, automated regression threshold,
cross-browser matrix, mobile ARM result, cold-start measure, memory peak, or
same-hardware competitor result.

Until those exist, "world-class browser readiness" and mobile superlatives are
not demonstrated. See `browser-siqs-readiness.md`.

### A030-013 — JavaScript preprocessing duplicates policy and blocks the UI

The page performs trial division, primality, perfect-power detection, and rho
synchronously on the main thread (`web/index.js:167-202`). On the fixed
192-bit case this audit measured about 91 ms before SIQS; the trial divider
tests every odd integer rather than primes (`web/numtheory.js:162-169`), and
perfect-power detection tests every exponent rather than only prime exponents
(`web/numtheory.js:152-159`).

This is modest beside a 256-bit sieve but material to the sub-second 192-bit
claim and to responsiveness. Move preprocessing into a Worker (ideally reuse
the Rust ladder to prevent policy drift), precompute primes, and test only prime
power exponents.

## Low

### A030-014 — Wasm handle generations eventually alias

`src/wasm.rs:10-69` uses a 16-bit generation and wraps it back to one. After
65,535 reuse cycles of a slot, an ancient stale handle can become valid again,
contradicting the absolute stale-handle language in `SPEC.md:544-546`.

Use a wider generation/slot split, retire wrapped slots, or narrow the contract
to bounded stale-handle detection.

### A030-015 — Release integrity metadata is minimal

Archives are reproducibly ordered on GNU tar and the release tag is signed,
which are positives. The repository does not provide a checked-in release
manifest with hashes, artifact sizes, toolchain versions, SBOM/provenance, or
an automated two-build reproducibility comparison.

Add SHA-256 manifests and signed provenance at minimum. For a "world-class"
release, generate SBOM and SLSA-style build attestations from pinned builders.

## Positive findings

- `#![deny(unsafe_code)]` is the default, with isolated C/Wasm exceptions and
  concrete safety comments.
- The C ABI exports exactly the eight intended `rusqsieve_*` symbols.
- Arithmetic has randomized `num-bigint` differential tests.
- Primality combines BPSW with fixed/seeded Miller–Rabin policy.
- SIQS relation collection is deterministic across out-of-order workers.
- Sparse filtering retains provenance and every dependency is re-verified
  against the original matrix.
- Factor outputs are sorted, probable-prime, and product-verifiable.
- The fixed browser corpus contains five exact balanced products at each of
  six tiers.
- Scalar/SIMD artifacts are small (156,823 and 160,402 bytes in this checkout).
- GNU and musl x86_64 release CLIs both executed successfully from their
  archives during the audit.
