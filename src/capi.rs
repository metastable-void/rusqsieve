//! Minimal owning C ABI for the high-level native factorization interface.

use crate::{
    FactorConfig, FactorError, Natural, PARTS, Parallelism, ParseNaturalError, ProgressAction,
    ProgressPhase, ProgressSnapshot, ProgressTotal, ProgressUnit, factor_with_progress,
};
use core::ffi::{c_char, c_int, c_void};
use core::ptr;
use std::ffi::{CStr, CString};
use std::panic::{AssertUnwindSafe, catch_unwind};

/// Status returned by the native C ABI.
#[repr(C)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RusqsieveStatus {
    /// Operation completed successfully.
    Ok = 0,
    /// A null pointer, zero worker count, or another invalid argument was supplied.
    InvalidArgument = 1,
    /// The input was not an unsigned decimal integer.
    InvalidDecimal = 2,
    /// The input was zero, exceeded capacity, or required the quadratic sieve on a composite
    /// wider than its supported range. Input width alone does not produce this status.
    InputOutOfRange = 3,
    /// No complete factorization was found.
    FactorizationFailed = 4,
    /// A panic was caught at the ABI boundary.
    InternalError = 5,
    /// A progress callback requested cancellation.
    Cancelled = 6,
}

const AUTO_THREAD_CAP: usize = 48;
const EXPLICIT_THREAD_CAP: usize = 256;
const C_ABI_VERSION: u32 = 2;

/// Return the native C ABI version.
#[unsafe(no_mangle)]
pub extern "C" fn rusqsieve_abi_version() -> u32 {
    C_ABI_VERSION
}

/// Return a static description for a status code.
#[unsafe(no_mangle)]
pub extern "C" fn rusqsieve_strerror(status: c_int) -> *const c_char {
    match status {
        0 => c"success".as_ptr(),
        1 => c"invalid argument".as_ptr(),
        2 => c"invalid decimal integer".as_ptr(),
        3 => c"input out of supported range".as_ptr(),
        4 => c"factorization failed".as_ptr(),
        5 => c"internal error".as_ptr(),
        6 => c"factorization cancelled".as_ptr(),
        _ => c"unknown rusqsieve status".as_ptr(),
    }
}

/// Rust-owned implementation of the incomplete C type `rusqsieve_factors`.
///
/// C only ever receives a pointer created by `Box::into_raw`. Exactly one call
/// to `rusqsieve_factors_free` converts it back with `Box::from_raw`.
pub struct RusqsieveFactors {
    allocation: Option<Box<FactorAllocation>>,
}

struct FactorAllocation {
    // These own the character buffers referenced by `pointers`.
    _strings: Box<[CString]>,
    pointers: Box<[*const c_char]>,
}

/// Progress snapshot passed to a native C callback.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct RusqsieveProgress {
    /// Stable phase code documented in `rusqsieve.h`.
    pub phase: u32,
    /// Completed work in `unit`.
    pub completed: u64,
    /// Total work, or zero when `total_kind` is unknown.
    pub total: u64,
    /// Zero for unknown, one for exact, and two for estimated.
    pub total_kind: u32,
    /// Stable unit code documented in `rusqsieve.h`.
    pub unit: u32,
}

/// Native progress callback. Return zero to continue or nonzero to cancel.
pub type RusqsieveProgressCallback =
    unsafe extern "C" fn(*const RusqsieveProgress, *mut c_void) -> c_int;

/// Allocate a new empty factorization result.
#[unsafe(no_mangle)]
pub extern "C" fn rusqsieve_factors_new() -> *mut RusqsieveFactors {
    Box::into_raw(Box::new(RusqsieveFactors { allocation: None }))
}

/// Destroy a factorization result and all strings it owns.
///
/// # Safety
///
/// A non-null pointer must have been returned by
/// [`rusqsieve_factors_new`] and must be freed exactly once.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rusqsieve_factors_free(factors: *mut RusqsieveFactors) {
    if factors.is_null() {
        return;
    }
    // SAFETY: the ABI requires a live pointer returned by
    // `rusqsieve_factors_new`, transferred here exactly once.
    unsafe { drop(Box::from_raw(factors)) };
}

/// Return the number of prime factors, including multiplicity.
///
/// # Safety
///
/// A non-null pointer must designate a live result object and must not be
/// concurrently mutated or destroyed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rusqsieve_factors_len(factors: *const RusqsieveFactors) -> usize {
    if factors.is_null() {
        return 0;
    }
    // SAFETY: the caller supplies a live constructor-returned object and does
    // not concurrently mutate or destroy it.
    let factors = unsafe { &*factors };
    factors
        .allocation
        .as_ref()
        .map_or(0, |allocation| allocation.pointers.len())
}

/// Return one borrowed decimal factor, or null when `index` is out of bounds.
///
/// # Safety
///
/// A non-null pointer must designate a live result object and must not be
/// concurrently mutated or destroyed. The returned string is borrowed from
/// that object.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rusqsieve_factors_get(
    factors: *const RusqsieveFactors,
    index: usize,
) -> *const c_char {
    if factors.is_null() {
        return ptr::null();
    }
    // SAFETY: the caller supplies a live constructor-returned object and does
    // not concurrently mutate or destroy it.
    let factors = unsafe { &*factors };
    factors
        .allocation
        .as_ref()
        .and_then(|allocation| allocation.pointers.get(index))
        .copied()
        .unwrap_or(ptr::null())
}

/// Factor a positive decimal integer (`1` has an empty factorization).
///
/// `threads == 0` selects available parallelism capped at 48. A nonzero value
/// is capped at 256; the engine may cap smaller inputs further.
///
/// This function is not constant-time and must not process secret values.
///
/// # Safety
///
/// `n` must point to a readable NUL-terminated string. `factors` must point to
/// a live object returned by [`rusqsieve_factors_new`] and must be exclusively
/// accessible for the duration of the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rusqsieve_factor(
    n: *const c_char,
    threads: usize,
    factors: *mut RusqsieveFactors,
) -> RusqsieveStatus {
    if n.is_null() || factors.is_null() {
        return RusqsieveStatus::InvalidArgument;
    }

    // Copy before clearing the prior result. This permits a caller to use a
    // string returned by `rusqsieve_factors_get` as the next input on the same
    // object.
    // SAFETY: the caller supplies a readable NUL-terminated C string.
    let input = unsafe { CStr::from_ptr(n) }.to_bytes().to_vec();
    // SAFETY: the caller supplies a live constructor-returned object with no
    // concurrent readers, mutation, or destruction during this call.
    let factors = unsafe { &mut *factors };
    factors.allocation = None;

    let mut observer = continue_progress;
    match catch_unwind(AssertUnwindSafe(|| {
        factor_impl(&input, threads, &mut observer)
    })) {
        Ok(Ok(allocation)) => {
            factors.allocation = allocation;
            RusqsieveStatus::Ok
        }
        Ok(Err(status)) => status,
        Err(_) => RusqsieveStatus::InternalError,
    }
}

/// Factor while reporting progress to a native callback.
///
/// A null callback behaves like [`rusqsieve_factor`]. A nonzero callback
/// return requests cooperative cancellation.
///
/// # Safety
///
/// The pointer requirements of [`rusqsieve_factor`] apply. In addition,
/// `context` must remain valid according to the callback's own contract for
/// the duration of this call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rusqsieve_factor_with_progress(
    n: *const c_char,
    threads: usize,
    factors: *mut RusqsieveFactors,
    callback: Option<RusqsieveProgressCallback>,
    context: *mut c_void,
) -> RusqsieveStatus {
    if n.is_null() || factors.is_null() {
        return RusqsieveStatus::InvalidArgument;
    }
    // SAFETY: inherited from this function's caller contract.
    let input = unsafe { CStr::from_ptr(n) }.to_bytes().to_vec();
    // SAFETY: inherited from this function's caller contract.
    let factors = unsafe { &mut *factors };
    factors.allocation = None;
    let mut observer = |snapshot: &ProgressSnapshot| {
        let Some(callback) = callback else {
            return ProgressAction::Continue;
        };
        let progress = c_progress(snapshot);
        // SAFETY: the function pointer and context validity are caller
        // obligations documented above.
        if unsafe { callback(&progress, context) } == 0 {
            ProgressAction::Continue
        } else {
            ProgressAction::Cancel
        }
    };
    match catch_unwind(AssertUnwindSafe(|| {
        factor_impl(&input, threads, &mut observer)
    })) {
        Ok(Ok(allocation)) => {
            factors.allocation = allocation;
            RusqsieveStatus::Ok
        }
        Ok(Err(status)) => status,
        Err(_) => RusqsieveStatus::InternalError,
    }
}

fn factor_impl(
    input: &[u8],
    threads: usize,
    observer: &mut dyn FnMut(&ProgressSnapshot) -> ProgressAction,
) -> Result<Option<Box<FactorAllocation>>, RusqsieveStatus> {
    let text = core::str::from_utf8(input).map_err(|_| RusqsieveStatus::InvalidDecimal)?;
    let n = Natural::<PARTS>::from_decimal(text).map_err(|error| match error {
        ParseNaturalError::Overflow => RusqsieveStatus::InputOutOfRange,
        ParseNaturalError::Empty | ParseNaturalError::InvalidDigit(_) => {
            RusqsieveStatus::InvalidDecimal
        }
        #[allow(unreachable_patterns)]
        _ => RusqsieveStatus::InvalidDecimal,
    })?;
    if n.is_zero() {
        return Err(RusqsieveStatus::InputOutOfRange);
    }

    let workers = if threads == 0 {
        std::thread::available_parallelism()
            .map_or(1, usize::from)
            .min(AUTO_THREAD_CAP)
    } else {
        threads.min(EXPLICIT_THREAD_CAP)
    };
    let parallelism = Parallelism::threads(workers).ok_or(RusqsieveStatus::InvalidArgument)?;
    let config = FactorConfig::default().with_parallelism(parallelism);
    let result = factor_with_progress(n, config, observer).map_err(|error| match error {
        FactorError::Cancelled => RusqsieveStatus::Cancelled,
        FactorError::SiqsCompositeTooLarge(_)
        | FactorError::InputTooLarge
        | FactorError::CapacityExceeded => RusqsieveStatus::InputOutOfRange,
        _ => RusqsieveStatus::FactorizationFailed,
    })?;

    let mut strings = Vec::with_capacity(result.distinct_len());
    let mut multiplicities = Vec::with_capacity(result.distinct_len());
    for (prime, exponent) in result.iter() {
        let decimal =
            CString::new(prime.to_string()).map_err(|_| RusqsieveStatus::InternalError)?;
        strings.push(decimal);
        multiplicities.push(exponent.get());
    }
    if strings.is_empty() {
        return Ok(None);
    }

    let strings = strings.into_boxed_slice();
    let total = multiplicities.iter().sum();
    let mut pointers = Vec::with_capacity(total);
    for (decimal, count) in strings.iter().zip(multiplicities) {
        pointers.extend(core::iter::repeat_n(decimal.as_ptr(), count));
    }
    Ok(Some(Box::new(FactorAllocation {
        _strings: strings,
        pointers: pointers.into_boxed_slice(),
    })))
}

fn continue_progress(_: &ProgressSnapshot) -> ProgressAction {
    ProgressAction::Continue
}

fn c_progress(snapshot: &ProgressSnapshot) -> RusqsieveProgress {
    let amount = snapshot.amount();
    let (total, total_kind) = match amount.total() {
        ProgressTotal::Unknown => (0, 0),
        ProgressTotal::Exact(total) => (total, 1),
        ProgressTotal::Estimated(total) => (total, 2),
        #[allow(unreachable_patterns)]
        _ => (0, 0),
    };
    RusqsieveProgress {
        phase: match snapshot.phase() {
            ProgressPhase::Preprocessing => 0,
            ProgressPhase::BuildingFactorBase => 1,
            ProgressPhase::Sieving => 2,
            ProgressPhase::LinearAlgebra => 6,
            ProgressPhase::ExtractingFactor => 7,
            ProgressPhase::Complete => 9,
            #[allow(unreachable_patterns)]
            _ => u32::MAX,
        },
        completed: amount.completed(),
        total,
        total_kind,
        unit: match amount.unit() {
            ProgressUnit::Candidates => 0,
            ProgressUnit::Primes => 1,
            ProgressUnit::SievePositions => 3,
            ProgressUnit::Relations => 4,
            ProgressUnit::MatrixRows => 5,
            ProgressUnit::MatrixColumns => 6,
            ProgressUnit::MatrixNonzeros => 7,
            ProgressUnit::Iterations => 8,
            ProgressUnit::MatrixProducts => 9,
            ProgressUnit::Tasks => 10,
            #[allow(unreachable_patterns)]
            _ => u32::MAX,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factor_output_is_sorted_repeated_and_reusable() {
        let output = rusqsieve_factors_new();
        assert!(!output.is_null());
        let input = CString::new("360").unwrap();
        let status = unsafe { rusqsieve_factor(input.as_ptr(), 1, output) };
        assert_eq!(status, RusqsieveStatus::Ok);
        assert_eq!(unsafe { rusqsieve_factors_len(output) }, 6);
        let actual: Vec<&str> = (0..6)
            .map(|index| unsafe {
                CStr::from_ptr(rusqsieve_factors_get(output, index))
                    .to_str()
                    .unwrap()
            })
            .collect();
        assert_eq!(actual, ["2", "2", "2", "3", "3", "5"]);
        assert!(unsafe { rusqsieve_factors_get(output, 6) }.is_null());

        let one = CString::new("1").unwrap();
        assert_eq!(
            unsafe { rusqsieve_factor(one.as_ptr(), 0, output) },
            RusqsieveStatus::Ok
        );
        assert_eq!(unsafe { rusqsieve_factors_len(output) }, 0);
        assert!(unsafe { rusqsieve_factors_get(output, 0) }.is_null());
        unsafe { rusqsieve_factors_free(output) };
    }

    #[test]
    fn invalid_input_clears_prior_result() {
        let output = rusqsieve_factors_new();
        let valid = CString::new("15").unwrap();
        assert_eq!(
            unsafe { rusqsieve_factor(valid.as_ptr(), 1, output) },
            RusqsieveStatus::Ok
        );
        let invalid = CString::new("not-a-number").unwrap();
        assert_eq!(
            unsafe { rusqsieve_factor(invalid.as_ptr(), 1, output) },
            RusqsieveStatus::InvalidDecimal
        );
        assert_eq!(unsafe { rusqsieve_factors_len(output) }, 0);
        unsafe { rusqsieve_factors_free(output) };
        unsafe { rusqsieve_factors_free(ptr::null_mut()) };
    }

    #[test]
    fn hostile_inputs_and_thread_count_are_bounded() {
        let output = rusqsieve_factors_new();
        let non_utf8 = CString::new(vec![0xff]).unwrap();
        assert_eq!(
            unsafe { rusqsieve_factor(non_utf8.as_ptr(), usize::MAX, output) },
            RusqsieveStatus::InvalidDecimal
        );

        let million_digits = CString::new(vec![b'9'; 1_000_000]).unwrap();
        assert_eq!(
            unsafe { rusqsieve_factor(million_digits.as_ptr(), usize::MAX, output) },
            RusqsieveStatus::InputOutOfRange
        );

        let embedded_nul = b"15\0ignored\0";
        assert_eq!(
            unsafe { rusqsieve_factor(embedded_nul.as_ptr().cast(), usize::MAX, output) },
            RusqsieveStatus::Ok
        );
        assert_eq!(unsafe { rusqsieve_factors_len(output) }, 2);
        unsafe { rusqsieve_factors_free(output) };
    }

    #[test]
    fn independent_results_are_concurrent() {
        let handles = ["1000036000099", "1000070001221"].map(|input| {
            std::thread::spawn(move || {
                let output = rusqsieve_factors_new();
                let input = CString::new(input).unwrap();
                let status = unsafe { rusqsieve_factor(input.as_ptr(), 2, output) };
                let len = unsafe { rusqsieve_factors_len(output) };
                unsafe { rusqsieve_factors_free(output) };
                (status, len)
            })
        });
        for handle in handles {
            let (status, len) = handle.join().unwrap();
            assert_eq!(status, RusqsieveStatus::Ok);
            assert_eq!(len, 2);
        }
    }
}
