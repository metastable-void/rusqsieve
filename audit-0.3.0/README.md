# rusqsieve 0.3.0 audit

Audit date: 2026-07-27  
Audited revision: `6ddc06f86c48d5e38870fb7bb2c358d0e59e32cb` (`v0.3.0`)  
Starting tree: clean; signed tag verified

## Verdict

**Go for the primary raw Wasm module and the native Rust/C/CLI core, subject to
the documented caveats. The additional bundled JavaScript frontend's packaging
defect found by this audit has been remediated in the post-audit working tree.**

The core implementation is substantially better than the release verdict may
suggest. It has a real logarithmic SIQS sieve, deterministic family merging,
well-tested fixed-capacity arithmetic, bounded preprocessing, verified matrix
dependencies, factor-product verification, a narrow unsafe boundary, and an
impressively disciplined history of measured performance work. All routine
Rust quality gates passed during this audit.

The raw scalar and SIMD Wasm modules build and execute successfully. At the
audited revision, the additional/informational JavaScript frontend had a
packaging defect:

- `web/index.js:46` imports `./coordinator.js`;
- `build-release.sh:246-248` does not copy `coordinator.js`;
- the existing `rusqsieve-0.3.0-wasm32-unknown-unknown.tar.gz` confirms that
  the file is absent.

This has since been fixed by adding `coordinator.js` to `package_wasm`. The
rebuilt archive contains all referenced frontend files. The original defect
did not invalidate the primary raw Wasm binaries or consumers providing their
own integration.

The browser orchestration failure states found by the audit have since been
hardened: messages are generation-scoped, jobs and submissions are accounted,
timeouts terminate runs, the runtime resets after failure, extraction requests
more relations instead of destroying the session, and packet/input boundaries
are validated. Focused failure checks now run in CI.

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
| Browser orchestration reliability | 7/10 | Major lifecycle/exhaustion defects are remediated; broader browser coverage remains |
| Release engineering | 7/10 | Native/raw Wasm artifacts are healthy; bundled demo omission is remediated |
| Performance evidence | 6/10 | Good fixed-corpus work; weak portability and competitor evidence |
| Maintainability | 7/10 | Clear contracts, but large kernel modules and duplicated JS/native policy remain |
| Overall release readiness | 7.5/10 | Core artifacts and glue are credible; evidence/portability gaps remain |

## Highest-priority actions

1. Remove or substantiate the mobile-device performance claim.
2. Make the native engine return `InsufficientRelations` explicitly when its
   family budget is exhausted.
3. Remove the 80.0 MiB of tracked `fuzz/target` build output and ignore nested
   Cargo target directories.
4. Add real-browser, packaged-artifact, musl, rustdoc, and package
   gates to CI.

## Reports

- [Detailed findings](findings.md)
- [Browser SIQS readiness and performance roadmap](browser-siqs-readiness.md)
- [Verification record](verification.md)
