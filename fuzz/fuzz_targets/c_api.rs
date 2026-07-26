#![no_main]

use core::ffi::{c_char, c_void};
use libfuzzer_sys::fuzz_target;

unsafe extern "C" {
    fn rusqsieve_factors_new() -> *mut c_void;
    fn rusqsieve_factors_free(factors: *mut c_void);
    fn rusqsieve_factor(n: *const c_char, threads: usize, factors: *mut c_void) -> i32;
}

fuzz_target!(|bytes: &[u8]| {
    // Keep individual cases bounded while covering non-UTF-8 and embedded-NUL
    // C-string behavior. The dedicated regression suite covers million-digit
    // rejection without making every fuzz iteration allocate a megabyte.
    let mut input = bytes[..bytes.len().min(4096)].to_vec();
    input.push(0);
    unsafe {
        let factors = rusqsieve_factors_new();
        if !factors.is_null() {
            let _ = rusqsieve_factor(input.as_ptr().cast(), 1, factors);
            rusqsieve_factors_free(factors);
        }
    }
});
