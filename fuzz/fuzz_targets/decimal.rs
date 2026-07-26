#![no_main]

use libfuzzer_sys::fuzz_target;
use rusqsieve::Natural;

fuzz_target!(|bytes: &[u8]| {
    if let Ok(text) = core::str::from_utf8(bytes) {
        let _ = Natural::<16>::from_decimal(text);
    }
});
