//! Portable quadratic-sieve setup, relation generation, and verification.
use crate::progress::FactorBaseProgress;
use crate::{Natural, legendre_u32, tonelli_shanks_u32};
use core::fmt;

#[derive(Clone, Copy, Debug)]
pub enum AutoOr<T> {
    Auto,
    Value(T),
}
#[derive(Clone, Copy, Debug)]
pub enum MultiplierChoice {
    Auto,
    Value(u32),
}
#[derive(Clone, Debug)]
pub struct QsConfig {
    pub multiplier: MultiplierChoice,
    pub factor_base_bound: AutoOr<u32>,
}
impl Default for QsConfig {
    fn default() -> Self {
        Self {
            multiplier: MultiplierChoice::Auto,
            factor_base_bound: AutoOr::Auto,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FactorBaseEntry {
    pub prime: u32,
    pub log_prime: u8,
    pub sqrt_n: u32,
}
#[derive(Clone, Debug)]
pub struct FactorBase {
    entries: Box<[FactorBaseEntry]>,
}
impl FactorBase {
    pub fn entries(&self) -> &[FactorBaseEntry] {
        &self.entries
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
#[derive(Clone, Debug)]
pub enum FactorBaseError {
    InvalidBound,
    FoundFactor(u32),
    NotFinished,
}
impl fmt::Display for FactorBaseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "factor-base error: {self:?}")
    }
}
impl std::error::Error for FactorBaseError {}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FactorBaseBuildStatus {
    InProgress,
    Complete,
}
pub struct FactorBaseBuilder<const P: usize> {
    n: Natural<P>,
    bound: u32,
    candidate: u32,
    primes: Vec<u32>,
    prime_cursor: usize,
    entries: Vec<FactorBaseEntry>,
    tested: u64,
    nonresidue: u64,
    finished: bool,
    /// Knuth-Schroeppel multiplier `k` (1 if none). Primes dividing `k` divide the working
    /// modulus `k·n` but are not factors of `n`; they are added as ramified factor-base entries
    /// (`sqrt_n = 0`) rather than reported as `FoundFactor`.
    multiplier: u64,
}
impl<const P: usize> FactorBaseBuilder<P> {
    pub fn new(n: Natural<P>, bound: u32) -> Result<Self, FactorBaseError> {
        if bound < 2 {
            return Err(FactorBaseError::InvalidBound);
        }
        let primes = segmented_primes(bound);
        Ok(Self {
            n,
            bound,
            candidate: 2,
            primes,
            prime_cursor: 0,
            multiplier: 1,
            entries: Vec::new(),
            tested: 0,
            nonresidue: 0,
            finished: false,
        })
    }
    pub fn step(&mut self, budget: usize) -> Result<FactorBaseBuildStatus, FactorBaseError> {
        for _ in 0..budget {
            let Some(&p) = self.primes.get(self.prime_cursor) else {
                self.candidate = self.bound.saturating_add(1);
                self.finished = true;
                return Ok(FactorBaseBuildStatus::Complete);
            };
            self.prime_cursor += 1;
            self.candidate = self
                .primes
                .get(self.prime_cursor)
                .copied()
                .unwrap_or_else(|| self.bound.saturating_add(1));
            self.tested += 1;
            let r = self.n.mod_u64(p as u64) as u32;
            if r == 0 && self.n != Natural::from_u64(p as u64) {
                // `p | working`. If `p | k` it only divides the multiplier, not `n` — fall through
                // and add it as a ramified prime (`r == 0` ⇒ `sqrt_n = 0`). Otherwise it is a real
                // factor of `n`.
                if !self.multiplier.is_multiple_of(p as u64) {
                    return Err(FactorBaseError::FoundFactor(p));
                }
            }
            if p == 2 || legendre_u32(r, p) >= 0 {
                self.entries.push(FactorBaseEntry {
                    prime: p,
                    log_prime: ((p as f64).ln() * 8.0).round().min(255.0) as u8,
                    sqrt_n: tonelli_shanks_u32(r, p).unwrap_or(r & 1),
                })
            } else {
                self.nonresidue += 1
            }
        }
        Ok(FactorBaseBuildStatus::InProgress)
    }
    pub fn progress(&self) -> FactorBaseProgress {
        FactorBaseProgress {
            bound: self.bound,
            searched_through: self.candidate.min(self.bound),
            primes_tested: self.tested,
            primes_accepted: self.entries.len() as u64,
            nonresidue_primes: self.nonresidue,
        }
    }
    pub fn finish(self) -> Result<FactorBase, FactorBaseError> {
        if !self.finished {
            return Err(FactorBaseError::NotFinished);
        }
        Ok(FactorBase {
            entries: self.entries.into_boxed_slice(),
        })
    }
}
fn segmented_primes(limit: u32) -> Vec<u32> {
    if limit < 2 {
        return Vec::new();
    }
    let root = (limit as f64).sqrt() as usize + 1;
    let mut composite = vec![false; root + 1];
    let mut base = Vec::new();
    for value in 2..=root {
        if composite[value] {
            continue;
        }
        base.push(value as u32);
        if value <= root / value {
            for multiple in (value * value..=root).step_by(value) {
                composite[multiple] = true;
            }
        }
    }
    const SEGMENT: u32 = 32 * 1024;
    let mut primes = Vec::new();
    let mut low = 2u32;
    while low <= limit {
        let high = low.saturating_add(SEGMENT - 1).min(limit);
        let mut marked = vec![false; (high - low + 1) as usize];
        for &prime in &base {
            let square = prime.saturating_mul(prime);
            let first = square.max(low.div_ceil(prime).saturating_mul(prime));
            if first > high {
                continue;
            }
            for multiple in (first..=high).step_by(prime as usize) {
                marked[(multiple - low) as usize] = true;
            }
        }
        primes.extend(
            marked
                .iter()
                .enumerate()
                .filter(|(_, is_composite)| !**is_composite)
                .map(|(offset, _)| low + offset as u32),
        );
        low = high.saturating_add(1);
    }
    primes
}

#[derive(Clone, Debug)]
pub struct PreparedFactorBase<const P: usize> {
    pub(crate) factor_base: FactorBase,
}
impl<const P: usize> PreparedFactorBase<P> {
    pub fn factor_base(&self) -> &FactorBase {
        &self.factor_base
    }
}
#[derive(Clone, Debug)]
pub enum QsError {
    InputTooSmall,
    Capacity,
    FactorFound(u32),
    FactorBase(FactorBaseError),
}
impl fmt::Display for QsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "quadratic-sieve setup error: {self:?}")
    }
}
impl std::error::Error for QsError {}
pub fn prepare_factor_base<const P: usize>(
    n: &Natural<P>,
    config: &QsConfig,
) -> Result<PreparedFactorBase<P>, QsError> {
    if *n < Natural::from_u64(2) {
        return Err(QsError::InputTooSmall);
    }
    let multiplier = match config.multiplier {
        MultiplierChoice::Auto => 1,
        MultiplierChoice::Value(k) => k.max(1),
    };
    let working = n
        .checked_mul(&Natural::from_u64(multiplier as u64))
        .ok_or(QsError::Capacity)?;
    let bound = match config.factor_base_bound {
        AutoOr::Value(v) => v,
        AutoOr::Auto => parameters::factor_base_bound(n.bit_len()),
    };
    let mut b = FactorBaseBuilder::new(working.clone(), bound).map_err(QsError::FactorBase)?;
    b.multiplier = multiplier as u64;
    loop {
        match b.step(4096) {
            Ok(FactorBaseBuildStatus::Complete) => break,
            Ok(FactorBaseBuildStatus::InProgress) => {}
            Err(FactorBaseError::FoundFactor(p)) => return Err(QsError::FactorFound(p)),
            Err(e) => return Err(QsError::FactorBase(e)),
        }
    }
    let base = b.finish().map_err(QsError::FactorBase)?;
    Ok(PreparedFactorBase { factor_base: base })
}
pub mod parameters {
    pub fn factor_base_bound(bits: usize) -> u32 {
        match bits {
            0..=40 => 200,
            41..=64 => 2_000,
            65..=96 => 10_000,
            97..=128 => 30_000,
            _ => 100_000,
        }
    }
    pub fn sieve_half_width(bits: usize) -> u32 {
        match bits {
            0..=64 => 4096,
            65..=128 => 32768,
            _ => 131072,
        }
    }

    /// Tuned SIQS engine parameters per input bit-length. These were selected by
    /// benchmarking balanced semiprimes against `flintqs`.
    ///
    /// The 193–224 range grows the factor base far beyond the older 60k bound. With Barrett-gated
    /// trial division (`engine.rs`, FLINT-style) the per-survivor factoring is cheap at any nfb, so
    /// these relation-starved bit-lengths benefit from a larger factor base: it raises smooth
    /// density, needing far fewer polynomials, dropping total sieve work (measured −15% at 208,
    /// −33% at 192, −45% at 224 vs the pre-optimization baseline). These nfb targets (≈5.7k at 208,
    /// ≈11k at 224) track FLINT's `qsieve_tune` table.
    ///
    /// The original above-224 native tuning study re-measured factor-base bounds at 256-bit,
    /// 4 threads (sieve + linear algebra, seconds): 150k → 79.6, 250k →
    /// 50.8, 350k → 39.9, **500k → 35.9**, 700k → 37.2, 1M → 45.7. Shrinking the base makes the
    /// sieve relation-starved; growing it makes the single-threaded dense solve explode (LA alone:
    /// 2.6 s at 500k, 5.5 s at 700k, 12.3 s at 1M). Sieve half-widths were checked the same way
    /// (327 680 → 35.9, 458 752 → 37.8, 655 360 → 43.8). An earlier revision of this comment
    /// claimed the optimum beyond 224 bits was *smaller* — ≈7k at 240 and ≈9k at 256 — which
    /// contradicted the table it documented and does not hold on measurement.
    ///
    /// The 209–264 tiers were subsequently tuned in the actual browser architecture: five fixed,
    /// balanced semiprimes per target, Chromium/V8 SIMD, eight independent workers. The retained
    /// `(bound, half-width, threshold)` settings and verified corpus changes are:
    ///
    /// - 216 bits: `(250k, 262144, −2)` → `(135k, 131072, 0)`, 5.097 s → 3.266 s (−35.9%);
    /// - 224 bits: `(250k, 262144, −2)` → `(150k, 131072, 0)`, 6.417 s → 5.340 s (−16.8%);
    /// - 232 bits: `(350k, 262144, −3)` → `(200k, 131072, −3)`, 11.258 s → 8.107 s (−28.0%);
    /// - 240 bits: `(350k, 262144, −3)` → `(350k, 131072, −1)`, 14.733 s → 13.629 s (−7.5%);
    /// - 256 bits: `(500k, 327680, −4)` → `(400k, 196608, −5)`, 38.334 s → 35.279 s (−8.0%).
    ///
    /// Nearby sweeps bracketed each retained point: 216-bit bounds of 120k and 150k, 224-bit bounds
    /// of 100k and 200k, 232-bit bounds of 175k and 250k, 240-bit widths of 98,304 and 196,608, and
    /// 256-bit bounds of 300k and 450k all regressed. At 232 bits, −3 and −4 thresholds were equal
    /// over the full corpus (8.107 s versus 8.103 s); −3 is retained as the less permissive setting.
    ///
    /// The browser tiers were re-swept after eight-pivot M4RI made large residual matrices cheaper.
    /// Five-case Chromium means retained three changes:
    ///
    /// - 224 bits: `(150k, 131072, 0)` → `(175k, 131072, 0)`, 5.176 s → 5.075 s (−2.0%);
    /// - 256 bits: `(400k, 196608, −5)` → `(450k, 196608, −4)`, 32.917 s → 32.259 s (−2.0%);
    /// - 272 bits: `(500k, 327680, −4)` → `(700k, 262144, −4)`, 105.110 s → 94.880 s (−9.7%).
    ///
    /// The 216-bit 150k boundary regressed. At 232 bits, 250k averaged 7.774 s versus 7.757 s for
    /// 200k. At 240 bits, 400k and 450k anchor runs regressed. At 256 bits, 500k was a wash and the
    /// five-case 450k/−3 gain was only 0.4%; −4 won the confirmation corpus. At 272 bits, 600k,
    /// 700k, and 800k bracketed the bound; 800k exceeded M4RI's working-set guard and fell back to
    /// scalar elimination. Half-widths 196,608 and 327,680 both lost to 262,144 at 700k. The
    /// measured 272-bit settings cover the 265–280 tier; unmeasured larger widths retain the prior
    /// conservative parameters.
    ///
    /// `thresh_adj` is the measured sieve-threshold offset in bits, added to
    /// `log2|g(x)| − log2(large-prime bound) − small-prime slack`. Deeper thresholds trade more
    /// survivors for fewer polynomials, and the optimum deepens with input size because
    /// per-polynomial cost grows faster than per-survivor cost. Measured optima on a 48-core Xeon
    /// 8259CL at 4 threads: 0 at 192-bit, −2 at 224-bit, −4 at 256-bit. The browser-tier values
    /// above supersede those native or interpolated values where their bit ranges overlap.
    #[derive(Clone, Copy, Debug)]
    pub struct EngineParams {
        pub factor_base_bound: u32,
        pub sieve_half_width: u32,
        pub thresh_adj: i32,
        pub large_prime_mult: u32,
    }
    pub fn engine_params(bits: usize) -> EngineParams {
        let (factor_base_bound, sieve_half_width, thresh_adj) = match bits {
            0..=100 => (3_000, 32_768, 0),
            101..=128 => (6_000, 32_768, 0),
            129..=160 => (40_000, 65_536, 0),
            161..=176 => (60_000, 65_536, 0),
            177..=192 => (100_000, 90_112, 0),
            193..=208 => (120_000, 131_072, -1),
            209..=216 => (135_000, 131_072, 0),
            217..=224 => (175_000, 131_072, 0),
            225..=232 => (200_000, 131_072, -3),
            233..=248 => (350_000, 131_072, -1),
            249..=264 => (450_000, 196_608, -4),
            265..=280 => (700_000, 262_144, -4),
            _ => (500_000, 327_680, -4),
        };
        EngineParams {
            factor_base_bound,
            sieve_half_width,
            thresh_adj,
            large_prime_mult: 256,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::parameters::engine_params;

    #[test]
    fn browser_tier_parameters_match_the_confirmed_m4ri_sweep() {
        let expected = [
            (216, 135_000, 131_072, 0),
            (224, 175_000, 131_072, 0),
            (232, 200_000, 131_072, -3),
            (240, 350_000, 131_072, -1),
            (256, 450_000, 196_608, -4),
            (272, 700_000, 262_144, -4),
            (281, 500_000, 327_680, -4),
        ];
        for (bits, bound, half_width, threshold) in expected {
            let params = engine_params(bits);
            assert_eq!(params.factor_base_bound, bound, "{bits}-bit bound");
            assert_eq!(params.sieve_half_width, half_width, "{bits}-bit width");
            assert_eq!(params.thresh_adj, threshold, "{bits}-bit threshold");
        }
    }
}
