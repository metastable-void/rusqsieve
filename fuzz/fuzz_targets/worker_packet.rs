#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|bytes: &[u8]| {
    let _ = rusqsieve::fuzz_validate_worker_packet(bytes);
});
