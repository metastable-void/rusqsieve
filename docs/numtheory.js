// Arbitrary-precision number theory in BigInt for the coordinator: trial division,
// primality, perfect powers and bounded Pollard-rho. Hard composites are handed to
// the parallel wasm quadratic sieve; everything else is resolved here.

export function modPow(base, exp, mod) {
  base %= mod;
  let r = 1n;
  while (exp > 0n) {
    if (exp & 1n) r = (r * base) % mod;
    base = (base * base) % mod;
    exp >>= 1n;
  }
  return r;
}

export function gcd(a, b) {
  while (b) {
    [a, b] = [b, a % b];
  }
  return a < 0n ? -a : a;
}

// Deterministic Miller-Rabin (exact for n < 3.3·10^24 with these bases; strong
// probable-prime beyond).
const MR_BASES = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n];
const TWO64 = 1n << 64n;

function millerRabinWitness(n, a) {
  let d = n - 1n;
  let s = 0n;
  while ((d & 1n) === 0n) {
    d >>= 1n;
    s++;
  }
  let x = modPow(a % n, d, n);
  if (x === 1n || x === n - 1n) return true;
  for (let i = 1n; i < s; i++) {
    x = (x * x) % n;
    if (x === n - 1n) return true;
    if (x === 1n) return false;
  }
  return false;
}

function jacobi(a, n) {
  a %= n;
  if (a < 0n) a += n;
  let result = 1;
  while (a !== 0n) {
    while ((a & 1n) === 0n) {
      a >>= 1n;
      const n8 = n & 7n;
      if (n8 === 3n || n8 === 5n) result = -result;
    }
    [a, n] = [n, a];
    if ((a & 3n) === 3n && (n & 3n) === 3n) result = -result;
    a %= n;
  }
  return n === 1n ? result : 0;
}

function strongLucasSelfridge(n) {
  let magnitude = 5n;
  let positive = true;
  let D;
  for (;;) {
    D = positive ? magnitude : -magnitude;
    const symbol = jacobi(D, n);
    if (symbol === -1) break;
    if (symbol === 0) return false;
    magnitude += 2n;
    positive = !positive;
  }
  const Q = (1n - D) / 4n;
  let odd = n + 1n;
  let s = 0n;
  while ((odd & 1n) === 0n) {
    odd >>= 1n;
    s++;
  }
  const mod = (x) => ((x % n) + n) % n;
  const half = (x) => {
    x = mod(x);
    return (x & 1n) === 0n ? x >> 1n : (x + n) >> 1n;
  };
  let U = 1n;
  let V = 1n;
  let Qk = mod(Q);
  const bits = odd.toString(2);
  for (let i = 1; i < bits.length; i++) {
    U = mod(U * V);
    V = mod(V * V - 2n * Qk);
    Qk = mod(Qk * Qk);
    if (bits[i] === "1") {
      const oldU = U;
      const oldV = V;
      U = half(oldU + oldV);
      V = half(D * oldU + oldV);
      Qk = mod(Qk * Q);
    }
  }
  if (U === 0n || V === 0n) return true;
  for (let r = 1n; r < s; r++) {
    V = mod(V * V - 2n * Qk);
    Qk = mod(Qk * Qk);
    if (V === 0n) return true;
  }
  return false;
}

export function isPrime(n) {
  if (n < 2n) return false;
  for (const p of MR_BASES) {
    if (n === p) return true;
    if (n % p === 0n) return false;
  }
  if (n >= TWO64) {
    const root = integerRoot(n, 2);
    if (root * root === n) return false;
    if (!millerRabinWitness(n, 2n) || !strongLucasSelfridge(n)) return false;
  }
  let d = n - 1n;
  let s = 0n;
  while ((d & 1n) === 0n) {
    d >>= 1n;
    s++;
  }
  witness: for (const a of MR_BASES) {
    let x = modPow(a, d, n);
    if (x === 1n || x === n - 1n) continue;
    for (let i = 1n; i < s; i++) {
      x = (x * x) % n;
      if (x === n - 1n) continue witness;
    }
    return false;
  }
  return true;
}

export function integerRoot(n, k) {
  if (n < 2n) return n;
  let lo = 1n;
  let hi = 1n << (BigInt(n.toString(2).length) / BigInt(k) + 1n);
  while (lo < hi) {
    const mid = (lo + hi + 1n) >> 1n;
    if (mid ** BigInt(k) <= n) lo = mid;
    else hi = mid - 1n;
  }
  return lo;
}

// If n = base^k for some k >= 2, return {base, k}; else null.
export function perfectPower(n) {
  const bits = n.toString(2).length;
  for (let k = 2; k <= bits; k++) {
    const r = integerRoot(n, k);
    if (r ** BigInt(k) === n) return { base: r, k };
  }
  return null;
}

// Peel prime factors below `limit` (default 100000). Returns the remaining cofactor.
export function trialDivide(n, out, limit = 100000n) {
  for (let p = 2n; p <= limit && p * p <= n; p += p === 2n ? 1n : 2n) {
    while (n % p === 0n) {
      out.push(p);
      n /= p;
    }
  }
  return n;
}

// The polynomial constants Brent tries in turn, `y^2 + c`. Each is an independent walk over the
// same modulus, which is what makes them safe to hand to separate workers: `deepPollard` in
// index.js deals them out so that N workers run N different walks and the first collision wins.
export const RHO_CONSTANTS = Array.from({ length: 31 }, (_, index) => BigInt(index + 1));

// Bounded Pollard-Brent. Returns a nontrivial factor of the composite `n`, or null
// if the iteration budget was exhausted (hand the number to the sieve instead).
export function pollardBrent(n, budget = 1 << 20, constants = RHO_CONSTANTS) {
  const run = pollardBrentSliced(n, budget, Infinity, constants);
  let step = run.next();
  while (!step.done) step = run.next();
  return step.value;
}

// The same search, sliced. It yields the iteration count every `slice` steps so a caller on the
// main thread can await a macrotask between slices and keep the page painting; the generator's
// return value is the factor or null, exactly as above.
//
// This exists for composites the sieve refuses. Below its ceiling a budget worth slicing is wasted
// work — the sieve is about to run and is better at those inputs — but above it rho is the whole
// attempt, and a budget deep enough to matter is a budget too long to block on.
// The default slice is measured, not nominal: 2^14 inner iterations is about 50 ms of BigInt work
// on a 512-bit modulus, below the frame budget a user would notice, while leaving the macrotask
// overhead of yielding at a few percent of the run.
export function* pollardBrentSliced(n, budget, slice = 1 << 14, constants = RHO_CONSTANTS) {
  if (n % 2n === 0n) return 2n;
  let steps = 0;
  // Slicing is counted separately from the budget: Brent's cycle-advance loop is not part of
  // `steps` and never has been, but it is the loop that grows without bound as `r` doubles, so it
  // is exactly the one that must not be allowed to block. Counting it as budget instead would
  // quietly halve the reach of every existing caller.
  let sinceYield = 0;
  for (const c of constants) {
    let y = 2n;
    let r = 1n;
    let q = 1n;
    let g = 1n;
    let x = 0n;
    let ys = 0n;
    const f = (v) => (v * v + c) % n;
    while (g === 1n && steps < budget) {
      x = y;
      for (let i = 0n; i < r; i++) {
        y = f(y);
        if (++sinceYield >= slice) {
          sinceYield = 0;
          yield steps;
        }
      }
      let k = 0n;
      while (k < r && g === 1n && steps < budget) {
        ys = y;
        const lim = r - k < 128n ? r - k : 128n;
        for (let i = 0n; i < lim; i++) {
          y = f(y);
          const d = x > y ? x - y : y - x;
          if (d !== 0n) q = (q * d) % n;
          steps++;
        }
        g = gcd(q, n);
        k += lim;
        sinceYield += Number(lim);
        if (sinceYield >= slice) {
          sinceYield = 0;
          yield steps;
        }
      }
      r <<= 1n;
    }
    if (g === n) {
      do {
        ys = f(ys);
        g = gcd(x > ys ? x - ys : ys - x, n);
        if (++sinceYield >= slice) {
          sinceYield = 0;
          yield steps;
        }
      } while (g === 1n && steps++ < budget);
    }
    if (g > 1n && g < n) return g;
    if (steps >= budget) return null;
  }
  return null;
}

// Bit length of a positive BigInt.
export const bitLength = (n) => (n <= 0n ? 0 : n.toString(2).length);

// Cryptographically-random odd BigInt of exactly `bits` bits, with the top two
// bits forced high so that a product of two such numbers lands in the intended
// combined bit length. Setting bits `bits-1` and `bits-2` gives n >= 3·2^(bits-2)
// = 0.75·2^bits, so a product of one `pbits` and one `qbits` value is at least
// 0.5625·2^(pbits+qbits) > 2^(pbits+qbits-1) and strictly below 2^(pbits+qbits) —
// exactly pbits+qbits bits.
function randomOdd(bits) {
  const bytes = new Uint8Array((bits + 7) >> 3);
  crypto.getRandomValues(bytes);
  let n = 0n;
  for (const b of bytes) n = (n << 8n) | BigInt(b);
  n &= (1n << BigInt(bits)) - 1n; // trim to `bits` bits
  n |= (1n << BigInt(bits - 1)) | (1n << BigInt(bits - 2)); // top two bits
  return n | 1n; // odd
}

// Small odd primes for a quick composite pre-screen before Miller-Rabin.
const SMALL_PRIMES = [3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n, 41n, 43n, 47n];

// A random prime of exactly `bits` bits, found by scanning upward from a random
// odd seed. Primality is whatever `isPrime` does, which above 2^64 is Baillie-PSW
// — a base-2 Miller-Rabin and a strong Lucas test — on top of the fixed-base rounds.
export function randomPrime(bits) {
  for (;;) {
    let n = randomOdd(bits);
    for (let tries = 0; tries < 20000; tries++, n += 2n) {
      if (SMALL_PRIMES.some((p) => n % p === 0n)) continue;
      if (isPrime(n)) return n;
    }
    // Astronomically unlikely to fall through; reseed and try again.
  }
}

// An RSA-style modulus of exactly `bits` bits: the product of two distinct random
// primes of roughly half that size each.
export function rsaNumber(bits) {
  const pbits = bits >> 1;
  const qbits = bits - pbits;
  let p, q;
  do {
    p = randomPrime(pbits);
    q = randomPrime(qbits);
  } while (p === q);
  return p * q;
}

// [prime, ...] with multiplicity -> sorted [{ prime, exponent }].
export function groupFactors(primes) {
  primes.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  const grouped = [];
  for (const p of primes) {
    const last = grouped[grouped.length - 1];
    if (last && last.prime === p) last.exponent++;
    else grouped.push({ prime: p, exponent: 1 });
  }
  return grouped;
}
