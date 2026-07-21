#![doc = include_str!("../README.md")]
#![cfg_attr(
    not(all(target_arch = "wasm32", target_os = "unknown")),
    forbid(unsafe_code)
)]

pub mod f2;
pub mod natural;
pub mod progress;
pub mod qs;

mod factor;
mod factors;
mod primality;
mod work;

pub mod engine;

#[cfg(any(unix, windows))]
mod native;

#[cfg(any(unix, windows))]
mod smallfactor;

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
mod wasm;

pub use factor::{
    AdvanceOutcome, FactorConfig, FactorError, FactorLimits, FactorSession, LocalWorkBudget,
    Parallelism, ProgressAction, ResourceLimitKind, SessionPhase, SmallFactorMethod, SubmitOutcome,
};
pub use factors::PrimeFactors;
pub use natural::{
    BufferTooSmall, CapacityError, ExtendedGcdResult, Montgomery, MontgomeryError, Natural, PARTS,
    ParseNaturalError, WideNatural, jacobi_u64, legendre_u32, tonelli_shanks_u32,
};
pub use primality::{PrimalityConfig, WitnessPolicy, is_probable_prime};
pub use progress::*;

#[cfg(any(unix, windows))]
pub use native::{factor, factor_with, factor_with_progress};

/// Stable low-level interfaces for custom native or WebAssembly schedulers.
pub mod low_level {
    pub use crate::engine::{
        EngineContext, EngineJob, EngineJobResult, EngineSession, execute as execute_engine_job,
        prepare as prepare_engine,
    };
    pub use crate::f2::{BlockLanczos, DependencySet, MatrixOperation, SparseBinaryMatrix};
    pub use crate::qs::{
        FactorBase, RawRelation, SieveContext, SieveScratch, prepare_siqs, sieve_job,
    };
    pub use crate::work::{
        JobHeader, KernelContexts, MatrixMultiplyJob, MatrixMultiplyResult, SieveJob, SieveResult,
        WorkJob, WorkResult, WorkerScratch, execute_job,
    };
}

/// Construct a fixed-capacity integer from a decimal literal at compile time.
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
