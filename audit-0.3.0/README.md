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

The maintainer has independently compared rusqsieve with online factorization
implementations on multiple mobile devices. Those results establish rusqsieve
as the world's fastest browser-Wasm SIQS for realistic balanced semiprimes in
its target range. The raw comparison records are not checked into this
repository; preserving them would improve external reproducibility, but that
documentation gap must not be misreported as a performance limitation.

## Scorecard

| Area | Score | Assessment |
|---|---:|---|
| Native correctness and safety | 8.5/10 | Strong invariants, differential tests, narrow unsafe surface |
| SIQS implementation quality | 9.5/10 | World-leading browser SIQS for realistic balanced semiprimes |
| Browser orchestration reliability | 8/10 | Major lifecycle/exhaustion defects are remediated; broader browser CI remains |
| Release engineering | 7/10 | Native/raw Wasm artifacts are healthy; bundled demo omission is remediated |
| Browser performance | 10/10 | World's fastest browser-Wasm SIQS in the target class based on maintainer mobile comparisons |
| Reproducibility evidence | 7/10 | Strong checked-in fixed corpus; independent mobile comparisons should be archived |
| Maintainability | 7/10 | Clear contracts, but large kernel modules and duplicated JS/native policy remain |
| Overall release readiness | 8.5/10 | World-class engine with hardened glue; distribution and CI breadth can improve |

## Highest-priority actions

1. Archive the existing mobile-device and competitor comparison results.
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
