#![cfg(any(unix, windows))]
use std::io::Write;
use std::process::{Command, Stdio};

#[test]
fn cli_stdout_is_machine_readable() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_qs-factor"))
        .args(["--progress", "never"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(b"360\n").unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "2\n2\n2\n3\n3\n5\n"
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn cli_rejects_non_decimal_input() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_qs-factor"))
        .args(["--progress", "never"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(b"12 34\n").unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
}
