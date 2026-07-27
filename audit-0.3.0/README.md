# rusqsieve 0.3.0 audit

Audit date: 2026-07-27  
Audited revision: `6ddc06f86c48d5e38870fb7bb2c358d0e59e32cb` (`v0.3.0`)  
Starting tree: clean; signed tag verified

## Verdict

**Go for the primary raw Wasm module and the native Rust/C/CLI core, subject to
the documented caveats. No-go only for treating the additional bundled
JavaScript frontend as a directly deployable application in its current
archive form.**

The core implementation is substantially better than the release verdict may
suggest. It has a real logarithmic SIQS sieve, deterministic family merging,
well-tested fixed-capacity arithmetic, bounded preprocessing, verified matrix
dependencies, factor-product verification, a narrow unsafe boundary, and an
impressively disciplined history of measured performance work. All routine
Rust quality gates passed during this audit.

The raw scalar and SIMD Wasm modules build and execute successfully. The
additional/informational JavaScript frontend in the Wasm archive nevertheless
has a packaging defect:

- `web/index.js:46` imports `./coordinator.js`;
- `build-release.sh:246-248` does not copy `coordinator.js`;
- the existing `rusqsieve-0.3.0-wasm32-unknown-unknown.tar.gz` confirms that
  the file is absent.

Deploying the bundled frontend without supplying the missing glue produces a
404 for the coordinator Worker and leaves startup waiting for a readiness
message. This blocks use of that frontend as a self-contained demo, but it does
not invalidate the primary raw Wasm binaries or consumers that provide their
own integration.

The browser orchestration also has untested failure states: stale errors can
poison a later run, family-budget exhaustion can hang forever, worker startup
errors can hang boot, and the coordinator destroys the session after a
one-shot linear-algebra attempt instead of collecting more relations.

Finally, the README's opening claim about timings "on a consumer mobile
device" is not supported by the repository's evidence. The documented
measurements are from a 96-thread Xeon under Node/V8 or headless Chromium, and
the project itself says representative mobile measurements remain future
work.

## Scorecard

| Area | Score | Assessment |
|---|---:|---|
| Native correctness and safety | 8.5/10 | Strong invariants, differential tests, narrow unsafe surface |
| SIQS implementation quality | 8/10 | Serious, measured engineering for balanced 192–272-bit inputs |
| Browser orchestration reliability | 4/10 | Happy path works; lifecycle and exhaustion paths are incomplete |
| Release engineering | 6/10 | Native and raw Wasm artifacts are healthy; bundled demo packaging is incomplete |
| Performance evidence | 6/10 | Good fixed-corpus work; weak portability and competitor evidence |
| Maintainability | 7/10 | Clear contracts, but large kernel modules and duplicated JS/native policy remain |
| Overall release readiness | 6.5/10 | Core artifacts are credible; bundled browser demo is not self-contained |

## Highest-priority actions

1. Clearly label JavaScript files as informational glue, or add
   `coordinator.js` and test the extracted demo end to end if the archive is
   intended to remain self-contained.
2. Remove or substantiate the mobile-device performance claim.
3. Make browser generations, timeouts, worker failures, and budget exhaustion
   terminal and testable.
4. Unify the Rust and JavaScript family budget and return
   `InsufficientRelations` when it is exhausted.
5. Keep the coordinator session alive after an unsuccessful extraction and
   resume sieving with a larger relation surplus.
6. Remove the 80.0 MiB of tracked `fuzz/target` build output and ignore nested
   Cargo target directories.
7. Add real-browser, packaged-artifact, musl, rustdoc, package, and failure-path
   gates to CI.

## Reports

- [Detailed findings](findings.md)
- [Browser SIQS readiness and performance roadmap](browser-siqs-readiness.md)
- [Verification record](verification.md)
