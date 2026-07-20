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
export function isPrime(n) {
  if (n < 2n) return false;
  for (const p of MR_BASES) {
    if (n === p) return true;
    if (n % p === 0n) return false;
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

// Bounded Pollard-Brent. Returns a nontrivial factor of the composite `n`, or null
// if the iteration budget was exhausted (hand the number to the sieve instead).
export function pollardBrent(n, budget = 1 << 22) {
  if (n % 2n === 0n) return 2n;
  let steps = 0;
  for (let c = 1n; c < 32n; c++) {
    let y = 2n;
    let r = 1n;
    let q = 1n;
    let g = 1n;
    let x = 0n;
    let ys = 0n;
    const f = (v) => (v * v + c) % n;
    while (g === 1n && steps < budget) {
      x = y;
      for (let i = 0n; i < r; i++) y = f(y);
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
      }
      r <<= 1n;
    }
    if (g === n) {
      do {
        ys = f(ys);
        g = gcd(x > ys ? x - ys : ys - x, n);
      } while (g === 1n && steps++ < budget);
    }
    if (g > 1n && g < n) return g;
    if (steps >= budget) return null;
  }
  return null;
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
