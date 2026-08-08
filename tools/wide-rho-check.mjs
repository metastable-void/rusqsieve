// The frontend's half of the deep Pollard-Brent policy.
//
// A composite above the sieve's ceiling is never handed to the coordinator, so the BigInt rho on
// the main thread is the only thing between it and an outright refusal. The opening peel is sized
// for a sieve run that is about to happen; above the ceiling none is, and `index.js` spends a much
// deeper budget before giving up. This checks the three properties that policy depends on: the
// deep budget reaches factors the opening peel cannot, the sliced generator returns exactly what
// the unsliced function does, and it actually yields often enough to keep a page painting.
import assert from "node:assert/strict";

import { pollardBrent, pollardBrentSliced } from "../web/numtheory.js";

// Same inputs as the native `wide_composites_are_split_by_rho_rather_than_refused` test: one small
// prime times one wide prime, at 512 and 1024 bits.
const CASES = [
  {
    bits: 512,
    n: 9501012405705509564680437712617447440170980081112656222237073910419870316392859702111963091481439276805995800801743430916377894473378632368751322056628119n,
    factor: 3667435003n,
    budget: 8 << 20,
    // Brent's cost varies by several-fold around its `1.2·sqrt(p)` expectation, so whether the
    // opening peel happens to reach a given 32-bit factor is luck. This one it misses.
    openingPeelReaches: false,
  },
  {
    bits: 1024,
    n: 140102229730795429799167923021188282773066057852400709530422767028102695167748110966423996618253720918241568745283524257975332015353165215374100761210034074640519308908102560817045658913006325713539052727661523689985681266686274209226917233091828080256346478868845059891172761884504662361609922427332741249613n,
    factor: 3479286313n,
    budget: 4 << 20,
    // ...and this one it happens to hit, which is the whole argument for a guaranteed tier sized
    // at hundreds of times the expected cost rather than at it.
    openingPeelReaches: true,
  },
];

// The wrapper must keep behaving exactly as it did before it was expressed in terms of the
// generator: below the ceiling this is the cheap peel and its budget is deliberately small.
assert.equal(pollardBrent(10403n, 1 << 15), 101n);
assert.equal(pollardBrent(2n * 5717n, 1 << 15), 2n);

for (const { bits, n, factor, budget, openingPeelReaches } of CASES) {
  // The gap this closes: the opening peel reaches a 32-bit factor only by luck, and refusing the
  // composite on the strength of it is what made a findable factor look unfindable.
  assert.equal(
    pollardBrent(n, 1 << 15),
    openingPeelReaches ? factor : null,
    `${bits}-bit: the opening peel's reach is not what this check assumes`,
  );

  const started = performance.now();
  const run = pollardBrentSliced(n, budget);
  let slices = 0;
  let worstSlice = 0;
  let last = performance.now();
  let step = run.next();
  while (!step.done) {
    const now = performance.now();
    worstSlice = Math.max(worstSlice, now - last);
    last = now;
    slices += 1;
    step = run.next();
  }
  const elapsed = performance.now() - started;
  assert.equal(step.value, factor, `${bits}-bit: deep budget returned ${step.value}`);
  assert.ok(slices > 0, `${bits}-bit: the search never yielded to the event loop`);
  // A slice that runs long enough to drop frames defeats the point of slicing. The threshold is
  // loose because this runs on unknown CI hardware; the measured worst slice is around 50 ms.
  assert.ok(worstSlice < 500, `${bits}-bit: worst slice ${worstSlice.toFixed(0)}ms`);
  console.log(
    `ok: ${bits}-bit composite split by deep rho -> ${factor} ` +
      `(${elapsed.toFixed(0)}ms, ${slices} slices, worst ${worstSlice.toFixed(0)}ms)`,
  );
}
