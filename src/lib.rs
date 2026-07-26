#![doc = include_str!("../README.md")]
#![deny(missing_docs)]
#![deny(unsafe_code)]
#![deny(unsafe_op_in_unsafe_fn)]

// Internal kernels are intentionally compiled in different combinations by
// native, Wasm coordinator, Wasm worker, and reference-engine builds.
#[allow(dead_code)]
mod f2;
#[allow(dead_code)]
mod natural;
#[allow(dead_code)]
mod progress;
#[allow(dead_code)]
mod qs;

#[allow(dead_code)]
mod engine;
#[allow(dead_code)]
mod factor;
mod factors;
#[allow(dead_code)]
mod primality;

#[cfg(any(unix, windows))]
mod native;

// Raw pointers are confined to the native C ABI boundary. The rest of the
// native crate remains under `deny(unsafe_code)`.
#[cfg(any(unix, windows))]
#[allow(unsafe_code)]
pub mod capi;

mod smallfactor;
mod u64math;

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
#[allow(unsafe_code)]
mod wasm;

pub use factor::{FactorConfig, FactorError, Parallelism, ProgressAction, ResourceLimitKind};
pub use factors::{ExpandedPrimeFactors, PrimeFactorIter, PrimeFactors};
pub use natural::{BufferTooSmall, CapacityError, InvalidDigit, Natural, ParseNaturalError};
pub use progress::{ProgressAmount, ProgressPhase, ProgressSnapshot, ProgressTotal, ProgressUnit};

/// Validate a serialized worker packet without executing it.
///
/// This hook exists only for the `fuzzing` feature and is not a stable wire
/// format API.
#[cfg(feature = "fuzzing")]
#[doc(hidden)]
pub fn fuzz_validate_worker_packet(bytes: &[u8]) -> bool {
    engine::validate_worker_packet(bytes)
}

pub(crate) use natural::{PARTS, jacobi_u64, legendre_u32, tonelli_shanks_u32};
#[cfg(any(unix, windows))]
pub(crate) use primality::{PrimalityConfig, WitnessPolicy, is_probable_prime};

#[cfg(any(unix, windows))]
pub use native::{factor, factor_with, factor_with_progress};

/// Constructs a fixed-capacity integer from a decimal literal at compile time.
///
/// The second argument is the number of 64-bit limbs. Invalid decimal text or
/// a value wider than the requested capacity causes a compile-time error.
///
/// ```
/// use rusqsieve::{Natural, natural};
///
/// const N: Natural<2> = natural!("340282366920938463463374607431768211455", 2);
/// assert_eq!(N, Natural::<2>::MAX);
/// ```
#[macro_export]
macro_rules! natural {
    ($value:literal, $parts:literal) => {{
        const VALUE: $crate::Natural<$parts> = match $crate::Natural::<$parts>::from_decimal($value)
        {
            Ok(value) => value,
            Err(_) => panic!("invalid or overflowing Natural literal"),
        };
        VALUE
    }};
}
