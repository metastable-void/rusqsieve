# Verification record

## Context

```text
revision: 6ddc06f86c48d5e38870fb7bb2c358d0e59e32cb
tag:      v0.3.0 (good EdDSA signature)
host:     Linux x86_64, 96 logical CPUs
CPU:      dual-socket Intel Xeon Platinum 8259CL
Node:     v24.15.0
```

The working tree was clean before audit commands. The audit only adds the
`audit-0.3.0/` reports.

## Passing checks

| Check | Result |
|---|---|
| `cargo fmt --check` | pass |
| `cargo clippy --locked --all-targets --all-features -- -D warnings` | pass |
| `cargo test --locked` | pass |
| `cargo test --locked --no-default-features` | pass |
| `cargo test --locked --profile release-test` | pass |
| `cargo check --locked --all-targets --all-features` | pass |
| `cargo check --locked --manifest-path fuzz/Cargo.toml` | pass |
| strict `cargo doc --no-deps --all-features` | pass |
| documentation-reference audit | pass |
| `cargo package --locked --allow-dirty` | pass; 54 files, 509.6 KiB |
| `make c-api-smoke` | pass |
| `make wasm` | pass/up to date |
| `make docs` and `make docs-verify` | pass; nine published files |
| `node tools/browser-arch-check.mjs` | pass; coordinator + four workers |
| all eight existing 0.3.0 tarballs pass gzip/tar integrity | pass |
| x86_64 GNU archive CLI execution | pass |
| x86_64 musl static-PIE archive CLI execution | pass |

Default tests reported:

```text
unit:        42 passed, 3 ignored
CLI:          3 passed
integration: 11 passed, 2 ignored
doctest:      3 passed
```

The scheduled large factorization corpus was not run in full during this audit;
its 29 above-128-bit entries are intentionally ignored outside scheduled CI.

## Browser-shaped benchmark

Command:

```sh
RUSQSIEVE_WASM=target/wasm-simd/wasm32-unknown-unknown/release/rusqsieve.wasm \
  node tools/wasm-bench.mjs \
  5845354724375454473909137928398990449217655808523662886639 8 2
```

Result:

```text
5845354724375454473909137928398990449217655808523662886639
= 75335908545075305094962839541
× 77590551932854658187989536979

relations=4829/4822
families=216
sieve=0.749 s
finish=0.189 s
wall=0.938 s
```

The harness verifies that the returned divisor is nontrivial and divides the
input.

A headless-Chromium audit attempt was inconclusive in the execution
environment: the sandboxed launch was denied and the approved unsandboxed
retry closed the page during navigation. This is recorded as an environmental
limitation, not as a rusqsieve failure. The Node Worker architecture check and
Node/V8 benchmark both completed.

## Artifact evidence

Current Wasm sizes:

```text
scalar: 156,823 bytes
SIMD:   160,402 bytes
```

The primary raw Wasm artifacts are present and passed the build and
Node/Worker checks. At the audited revision, the archive's additional
JavaScript glue was incomplete as a self-contained demo:

```text
index.js imports: ./coordinator.js
archive contains: no web/coordinator.js
```

This did not invalidate the raw Wasm modules or custom integrations.

After remediation, `build-release.sh wasm32-unknown-unknown` produced an
archive containing:

```text
abi.js
coordinator.js
index.css
index.html
index.js
numtheory.js
rusqsieve-simd.wasm
rusqsieve.wasm
serve.mjs
worker.js
```

The extracted package passed a same-directory frontend-reference check.

Post-audit browser-glue verification also passed:

```text
node tools/browser-arch-check.mjs
node tools/browser-glue-failure-check.mjs
node --check web/abi.js web/coordinator.js web/index.js web/worker.js
```

The focused failure check covers invalid/stale packet handles, bad magic,
wrong packet kinds, length mismatches, allocator rollback, commands in the
wrong state, unknown commands, and malformed relation batches. The 192-bit
eight-worker SIMD benchmark remained factor-correct at 0.937 seconds in the
post-hardening run.

The GNU release CLI is a dynamically linked x86-64 PIE. The musl release CLI is
a static PIE. Both factored small smoke inputs correctly from freshly extracted
archives.

## Repository hygiene evidence

```text
tracked fuzz/target files: 162
tracked logical size:      83,880,588 bytes (80.0 MiB)
```

Running the fuzz check modifies these tracked compiler outputs, demonstrating
why generated target directories must not be versioned.

## Scope limitations

- No external competitor implementation was available for controlled
  same-host comparison.
- No mobile device was connected.
- No full 216–272-bit Playwright corpus was rerun.
- Cross-target archives other than x86_64 GNU/musl were inspected for tar
  integrity but not executed on native hardware.
- This was a source, test, artifact, and targeted performance audit; it was not
  a formal proof or side-channel review.
