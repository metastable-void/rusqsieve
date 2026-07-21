#[cfg(target_arch = "wasm32")]
compile_error!("the qs-factor CLI is unavailable on wasm32 targets");

use rusqsieve::{Natural, engine};
use std::io::{self, IsTerminal, Read, Write};
use std::time::Instant;

#[derive(Clone, Copy, Eq, PartialEq)]
enum ProgressMode {
    Auto,
    Always,
    Never,
}

struct Options {
    progress: ProgressMode,
    threads: Option<usize>,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("qs-factor: {error}");
        std::process::exit(2)
    }
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let options = parse_options()?;
    let mut text = String::new();
    io::stdin()
        .read_to_string(&mut text)
        .map_err(|error| format!("cannot read stdin: {error}"))?;
    let input = text.trim();
    if input.is_empty() || input.bytes().any(|byte| !byte.is_ascii_digit()) {
        return Err("input must be exactly one unsigned decimal integer".into());
    }

    // Parse with the crate's public type first so CLI capacity and syntax stay
    // consistent with the library API. The SIQS engine has a deliberately
    // smaller practical bound than Natural's storage capacity.
    let natural = Natural::<16>::from_decimal(input).map_err(|error| error.to_string())?;
    if natural.is_zero() {
        return Err("zero has no prime factorization".into());
    }
    if natural.is_one() {
        return Ok(());
    }
    if natural.bit_len() > 512 {
        return Err(format!(
            "input is {} bits; qs-factor's SIQS engine supports at most 512 bits",
            natural.bit_len()
        ));
    }

    let show_progress = options.progress == ProgressMode::Always
        || options.progress == ProgressMode::Auto && io::stderr().is_terminal();
    let thread_count = options.threads.unwrap_or_else(default_parallelism);

    if show_progress {
        eprintln!(
            "Factoring {}-bit input with {thread_count} worker{}",
            natural.bit_len(),
            if thread_count == 1 { "" } else { "s" }
        );
    }
    let factors = engine::factor(natural.clone(), thread_count, |snapshot| {
        if show_progress {
            match snapshot.phase {
                engine::EnginePhase::Preprocessing => eprint!("\rpreprocessing"),
                engine::EnginePhase::BuildingFactorBase => eprint!("\rbuilding factor base"),
                engine::EnginePhase::Sieving => eprint!(
                    "\rsieving: {}/~{} relations, {} polynomials, {} workers",
                    snapshot.relations, snapshot.target, snapshot.polynomials, snapshot.workers
                ),
                engine::EnginePhase::LinearAlgebra => {
                    eprint!("\rlinear algebra: {} relations", snapshot.relations)
                }
                engine::EnginePhase::Extracting => eprint!("\rextracting factors"),
            }
            let _ = io::stderr().flush();
        }
    })
    .map_err(|error| format!("factor engine failed: {error}"))?;
    if show_progress {
        eprintln!();
    }

    // The engine checks its own result. Verify it once more through Natural
    // before exposing machine-readable output, so incomplete results never
    // look like a successful factorization.
    let mut verified_product = Natural::<16>::ONE;
    let mut rendered = Vec::with_capacity(factors.len());
    for value in factors {
        let decimal = value.to_string();
        verified_product = verified_product
            .checked_mul(&value)
            .ok_or_else(|| "factor product exceeded input capacity".to_string())?;
        rendered.push(decimal);
    }
    if verified_product != natural {
        return Err("factor engine returned an incomplete factorization".into());
    }

    let stdout = io::stdout();
    let mut output = stdout.lock();
    for factor in rendered {
        if let Err(error) = writeln!(output, "{factor}") {
            if error.kind() == io::ErrorKind::BrokenPipe {
                return Ok(());
            }
            return Err(format!("cannot write stdout: {error}"));
        }
    }

    if show_progress {
        eprintln!("elapsed: {:.3} s", started.elapsed().as_secs_f64());
    }
    Ok(())
}

fn default_parallelism() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}

fn parse_options() -> Result<Options, String> {
    let mut progress = ProgressMode::Auto;
    let mut threads = None;
    let mut args = std::env::args().skip(1);
    while let Some(argument) = args.next() {
        match argument.as_str() {
            "--progress" => {
                progress = match args.next().as_deref() {
                    Some("auto") => ProgressMode::Auto,
                    Some("always") => ProgressMode::Always,
                    Some("never") => ProgressMode::Never,
                    _ => return Err("--progress requires auto, always, or never".into()),
                }
            }
            "--threads" => {
                let value = args
                    .next()
                    .ok_or_else(|| "--threads requires auto or a positive integer".to_string())?;
                threads = if value == "auto" {
                    None
                } else {
                    let count = value
                        .parse::<usize>()
                        .map_err(|_| "--threads requires auto or a positive integer".to_string())?;
                    if count == 0 {
                        return Err("--threads must be greater than zero".into());
                    }
                    Some(count)
                };
            }
            "-h" | "--help" => {
                println!(
                    "Usage: qs-factor [OPTIONS]\n\
                     Reads one unsigned decimal integer from standard input.\n\n\
                     Options:\n  \
                       --progress auto|always|never\n  \
                       --threads auto|N\n  \
                       -h, --help"
                );
                std::process::exit(0)
            }
            _ => return Err(format!("unknown argument: {argument}")),
        }
    }
    Ok(Options { progress, threads })
}
