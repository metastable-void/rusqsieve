use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

fn files_below(path: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(path).expect("read source directory") {
        let path = entry.expect("read directory entry").path();
        if path.is_dir() {
            files_below(&path, out);
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            out.push(path);
        }
    }
}

fn main() {
    let spec = fs::read_to_string("SPEC.md").expect("read SPEC.md");
    let headings: BTreeSet<String> = spec
        .lines()
        .filter_map(|line| {
            let line = line.trim_start_matches('#').trim_start();
            let section = line.split_whitespace().next()?;
            section
                .chars()
                .all(|character| character.is_ascii_digit() || character == '.')
                .then(|| section.trim_end_matches('.').to_owned())
        })
        .collect();

    let mut sources = Vec::new();
    files_below(Path::new("src"), &mut sources);
    let mut failures = Vec::new();
    for path in sources {
        let source = fs::read_to_string(&path).expect("read Rust source");
        for (line_number, line) in source.lines().enumerate() {
            let Some(comment) = line.find("//").map(|offset| &line[offset..]) else {
                continue;
            };
            let mut rest = comment;
            while let Some(offset) = rest.find("SPEC §") {
                let reference = rest[offset + "SPEC §".len()..]
                    .split(|character: char| {
                        !(character.is_ascii_digit() || character == '.')
                    })
                    .next()
                    .unwrap_or("")
                    .trim_end_matches('.');
                if reference.is_empty() || !headings.contains(reference) {
                    failures.push(format!(
                        "{}:{}: unresolved SPEC §{}",
                        path.display(),
                        line_number + 1,
                        reference
                    ));
                }
                rest = &rest[offset + "SPEC §".len()..];
            }

            for token in comment.split_whitespace() {
                let candidate = token.trim_matches(|character: char| {
                    matches!(
                        character,
                        '`' | '\'' | '"' | '(' | ')' | '[' | ']' | ',' | ':' | ';'
                    )
                });
                if candidate.ends_with(".md")
                    && !Path::new(candidate).exists()
                    && !path
                        .parent()
                        .is_some_and(|parent| parent.join(candidate).exists())
                {
                    failures.push(format!(
                        "{}:{}: referenced file does not exist: {}",
                        path.display(),
                        line_number + 1,
                        candidate
                    ));
                }
            }
        }
    }

    if !failures.is_empty() {
        for failure in failures {
            eprintln!("{failure}");
        }
        std::process::exit(1);
    }
}
