#[cfg(target_arch = "wasm32")]
compile_error!("the qs-factor CLI is unavailable on wasm32 targets");

use rusqsieve::{FactorConfig, Natural, ProgressAction, factor_with_progress};
use std::io::{self, IsTerminal, Read, Write};

#[derive(Clone, Copy, Eq, PartialEq)]
enum ProgressMode {
    Auto,
    Always,
    Never,
}
fn main() {
    if let Err(e) = run() {
        eprintln!("qs-factor: {e}");
        std::process::exit(2)
    }
}
fn run() -> Result<(), String> {
    let mut mode = ProgressMode::Auto;
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--progress" => {
                mode = match args.next().as_deref() {
                    Some("auto") => ProgressMode::Auto,
                    Some("always") => ProgressMode::Always,
                    Some("never") => ProgressMode::Never,
                    _ => return Err("--progress requires auto, always, or never".into()),
                }
            }
            "-h" | "--help" => {
                println!(
                    "Usage: qs-factor [--progress auto|always|never]\nReads one unsigned decimal integer from standard input."
                );
                return Ok(());
            }
            _ => return Err(format!("unknown argument: {a}")),
        }
    }
    let mut text = String::new();
    io::stdin()
        .read_to_string(&mut text)
        .map_err(|e| format!("cannot read stdin: {e}"))?;
    let input = text.trim();
    if input.is_empty() || input.bytes().any(|b| !b.is_ascii_digit()) {
        return Err("input must be exactly one unsigned decimal integer".into());
    }
    let n = Natural::<16>::from_decimal(input).map_err(|e| e.to_string())?;
    let show =
        mode == ProgressMode::Always || (mode == ProgressMode::Auto && io::stderr().is_terminal());
    let factors = factor_with_progress(n, FactorConfig::default(), |p| {
        if show {
            let fraction = p
                .amount
                .fraction()
                .map(|x| format!(" {:5.1}%", x * 100.0))
                .unwrap_or_default();
            eprint!("\r{:?}{fraction}", p.phase);
            let _ = io::stderr().flush();
        }
        ProgressAction::Continue
    })
    .map_err(|e| e.to_string())?;
    if show {
        eprintln!()
    }
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for (p, e) in factors.iter() {
        for _ in 0..e.get() {
            if let Err(err) = writeln!(out, "{p}") {
                if err.kind() == io::ErrorKind::BrokenPipe {
                    return Ok(());
                }
                return Err(format!("cannot write stdout: {err}"));
            }
        }
    }
    Ok(())
}
