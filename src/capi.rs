//! Minimal owning C ABI for the high-level native factorization interface.

use crate::{FactorConfig, Natural, PARTS, Parallelism, ParseNaturalError, factor_with};
use core::ffi::{c_char, c_int};
use core::ptr;
use std::ffi::{CStr, CString};
use std::panic::{AssertUnwindSafe, catch_unwind};

const RUSQSIEVE_OK: c_int = 0;
const RUSQSIEVE_INVALID_ARGUMENT: c_int = 1;
const RUSQSIEVE_INVALID_DECIMAL: c_int = 2;
const RUSQSIEVE_INPUT_OUT_OF_RANGE: c_int = 3;
const RUSQSIEVE_FACTORIZATION_FAILED: c_int = 4;
const RUSQSIEVE_INTERNAL_ERROR: c_int = 5;
const AUTO_THREAD_CAP: usize = 48;

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

/// Allocate a new empty factorization result.
#[unsafe(no_mangle)]
pub extern "C" fn rusqsieve_factors_new() -> *mut RusqsieveFactors {
    Box::into_raw(Box::new(RusqsieveFactors { allocation: None }))
}

/// Destroy a factorization result and all strings it owns.
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
/// requests that exact worker count (the engine may still cap tiny inputs).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rusqsieve_factor(
    n: *const c_char,
    threads: usize,
    factors: *mut RusqsieveFactors,
) -> c_int {
    if n.is_null() || factors.is_null() {
        return RUSQSIEVE_INVALID_ARGUMENT;
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

    match catch_unwind(AssertUnwindSafe(|| factor_impl(&input, threads))) {
        Ok(Ok(allocation)) => {
            factors.allocation = allocation;
            RUSQSIEVE_OK
        }
        Ok(Err(status)) => status,
        Err(_) => RUSQSIEVE_INTERNAL_ERROR,
    }
}

fn factor_impl(input: &[u8], threads: usize) -> Result<Option<Box<FactorAllocation>>, c_int> {
    let text = core::str::from_utf8(input).map_err(|_| RUSQSIEVE_INVALID_DECIMAL)?;
    let n = Natural::<PARTS>::from_decimal(text).map_err(|error| match error {
        ParseNaturalError::Overflow => RUSQSIEVE_INPUT_OUT_OF_RANGE,
        ParseNaturalError::Empty | ParseNaturalError::InvalidDigit { .. } => {
            RUSQSIEVE_INVALID_DECIMAL
        }
    })?;
    if n.is_zero() {
        return Err(RUSQSIEVE_INPUT_OUT_OF_RANGE);
    }

    let workers = if threads == 0 {
        std::thread::available_parallelism()
            .map_or(1, usize::from)
            .min(AUTO_THREAD_CAP)
    } else {
        threads
    };
    let parallelism = Parallelism::threads(workers).ok_or(RUSQSIEVE_INVALID_ARGUMENT)?;
    let config = FactorConfig::default().with_parallelism(parallelism);
    let result = factor_with(n, config).map_err(|_| RUSQSIEVE_FACTORIZATION_FAILED)?;

    let mut strings = Vec::with_capacity(result.distinct_len());
    let mut multiplicities = Vec::with_capacity(result.distinct_len());
    for (prime, exponent) in result.iter() {
        let decimal = CString::new(prime.to_string()).map_err(|_| RUSQSIEVE_INTERNAL_ERROR)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factor_output_is_sorted_repeated_and_reusable() {
        let output = rusqsieve_factors_new();
        assert!(!output.is_null());
        let input = CString::new("360").unwrap();
        let status = unsafe { rusqsieve_factor(input.as_ptr(), 1, output) };
        assert_eq!(status, RUSQSIEVE_OK);
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
            RUSQSIEVE_OK
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
            RUSQSIEVE_OK
        );
        let invalid = CString::new("not-a-number").unwrap();
        assert_eq!(
            unsafe { rusqsieve_factor(invalid.as_ptr(), 1, output) },
            RUSQSIEVE_INVALID_DECIMAL
        );
        assert_eq!(unsafe { rusqsieve_factors_len(output) }, 0);
        unsafe { rusqsieve_factors_free(output) };
        unsafe { rusqsieve_factors_free(ptr::null_mut()) };
    }
}
