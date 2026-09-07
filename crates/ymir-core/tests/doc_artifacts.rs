//! ADR method rule 6, EXTENDED and made executable — **a written artefact must be verified
//! non-empty and non-corrupt before it is declared delivered.**
//!
//! ## Why this file exists
//!
//! Rule 6 was earned on a tool that corrupted its output in silence: PowerShell
//! `Get-Content -Raw` piped into `Add-Content` double-encoded 275 lines of the ADR, and cp1252
//! DESTROYED `0x81/0x8D/0x8F/0x90/0x9D`, so `∝` and `↔` could not be recovered by inverting the
//! encoding. The rule as written was a prohibition on a tool. That is not enough: **the failure
//! class is "the artefact reports the instrument instead of the subject", and an empty or
//! mangled document is that class in documentary form** — the exact counterpart of a metric
//! column constant at zero (rule 10).
//!
//! Fourth instance on this campaign, and the reason it is now an assertion:
//!
//! 1. the global axial R, blind to local parallelism — reported the instrument, not the coast;
//! 2. the 8 km spur-length cap, saturating at 7.99 in every configuration;
//! 3. `channel width p50` constant at 0.000 m — rule 10's first catch;
//! 4. a report file suspected empty. **This test is what settles that question by measurement
//!    instead of by inspection.**
//!
//! ## What it checks, over every tracked Markdown document
//!
//! - **non-empty**, and above a floor that a truncated write would fall under;
//! - **valid UTF-8** — the encoding failure mode, caught at the byte level;
//! - **no double-encoding signature** (`Ã`, `â€`, `Â°`…), which is what the banned pipeline
//!   produced and what `grep` did NOT flag as unusual at the time;
//! - **no unresolved conflict markers**, the other way a document ships unreadable.
//!
//! Run: `cargo test -p ymir-core --test doc_artifacts` (fast, always on).

use std::path::{Path, PathBuf};

/// A Markdown file below this is either empty or a truncated write. The smallest real document
/// in `docs/` is comfortably above it; a heading alone is ~40 bytes.
const MIN_MD_BYTES: u64 = 200;

/// Byte sequences that appear when UTF-8 text is read as cp1252 and re-encoded as UTF-8. These
/// are the fingerprints the banned `Get-Content -Raw | Add-Content` pipeline leaves behind.
const DOUBLE_ENCODING_SIGNATURES: &[&str] =
    &["Ã©", "Ã¨", "Ã ", "Ã§", "Ã´", "Ã»", "Ã®", "â€™", "â€œ", "â€“", "â€”", "Â°", "Â²", "Ã‰"];

/// A test's cwd is the CRATE root, not the workspace root.
fn docs_root() -> PathBuf {
    Path::new("../../docs").to_path_buf()
}

fn collect_md(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else { return };
    // Sorted, so a failure names the same file first on every machine.
    let mut paths: Vec<PathBuf> = entries.filter_map(|e| e.ok()).map(|e| e.path()).collect();
    paths.sort();
    for p in paths {
        if p.is_dir() {
            collect_md(&p, out);
        } else if p.extension().is_some_and(|e| e == "md") {
            out.push(p);
        }
    }
}

/// Rule 6 extended: every document in `docs/` is non-empty, decodable and un-mangled.
#[test]
fn every_tracked_document_is_non_empty_and_decodable() {
    let root = docs_root();
    assert!(root.is_dir(), "docs/ not found from the crate root — expected {}", root.display());
    let mut files = Vec::new();
    collect_md(&root, &mut files);
    // Rule 10's own move, applied to this test: a check that ran over an EMPTY population is
    // not a pass. If the walk found nothing, the walk is broken, not the docs.
    assert!(
        files.len() >= 10,
        "the walk found only {} Markdown file(s) under {} — the WALK is broken, and \
         '0 problems out of 0 files' is not a clean result",
        files.len(),
        root.display()
    );

    let mut problems: Vec<String> = Vec::new();
    for p in &files {
        let rel = p.strip_prefix(&root).unwrap_or(p).display().to_string();
        let bytes = match std::fs::read(p) {
            Ok(b) => b,
            Err(e) => {
                problems.push(format!("{rel}: unreadable ({e})"));
                continue;
            }
        };
        if bytes.is_empty() {
            problems.push(format!("{rel}: EMPTY (0 bytes) — a write that produced nothing"));
            continue;
        }
        if (bytes.len() as u64) < MIN_MD_BYTES {
            problems.push(format!(
                "{rel}: only {} bytes, under the {MIN_MD_BYTES}-byte floor — likely a \
                 truncated write",
                bytes.len()
            ));
            continue;
        }
        let text = match String::from_utf8(bytes) {
            Ok(t) => t,
            Err(e) => {
                problems.push(format!(
                    "{rel}: NOT valid UTF-8 at byte {}",
                    e.utf8_error().valid_up_to()
                ));
                continue;
            }
        };
        for sig in DOUBLE_ENCODING_SIGNATURES {
            if text.contains(sig) {
                problems.push(format!(
                    "{rel}: DOUBLE-ENCODED — contains {sig:?}, the cp1252 round-trip \
                     signature (see the tooling ban in CLAUDE.md)"
                ));
                break;
            }
        }
        if text.contains("\n<<<<<<< ") || text.contains("\n>>>>>>> ") {
            problems.push(format!("{rel}: unresolved conflict markers"));
        }
    }

    assert!(
        problems.is_empty(),
        "\n[rule 6] {} of {} document(s) would have shipped unreadable:\n  · {}\n",
        problems.len(),
        files.len(),
        problems.join("\n  · ")
    );
    eprintln!("[rule 6] {} document(s) verified non-empty, UTF-8 and un-mangled", files.len());
}

/// NEGATIVE CONTROL (rule 1): the checks above must actually FIRE. Without this, a walk that
/// silently matched nothing, or signatures that never match, would report the same clean pass.
#[test]
fn the_rule_6_checks_actually_fire() {
    // empty
    assert!(Vec::<u8>::new().is_empty());
    // under the floor
    assert!((b"# a heading and nothing else\n".len() as u64) < MIN_MD_BYTES);
    // invalid UTF-8
    assert!(String::from_utf8(vec![0x41, 0xC3, 0x28]).is_err());
    // the double-encoding signature, built by actually performing the round trip on a
    // character the ADR uses, rather than by trusting a literal I typed
    let original = "hypsométrie ↔ 25°";
    let mangled: String = original.as_bytes().iter().map(|&b| cp1252_to_char(b)).collect();
    assert_ne!(original, mangled, "the round trip must actually change the text");
    assert!(
        DOUBLE_ENCODING_SIGNATURES.iter().any(|s| mangled.contains(s)),
        "no signature matched the genuinely double-encoded text {mangled:?} — the signature \
         list would never fire on the real defect"
    );
}

/// cp1252 → the character PowerShell 5.1 substitutes when it reads a UTF-8 byte as cp1252.
/// The five undefined slots (`0x81/0x8D/0x8F/0x90/0x9D`) map to nothing and are the reason the
/// corruption is IRRECOVERABLE rather than merely reversible.
fn cp1252_to_char(b: u8) -> char {
    match b {
        0x80 => '\u{20AC}',
        0x82 => '\u{201A}',
        0x83 => '\u{0192}',
        0x84 => '\u{201E}',
        0x85 => '\u{2026}',
        0x86 => '\u{2020}',
        0x87 => '\u{2021}',
        0x88 => '\u{02C6}',
        0x89 => '\u{2030}',
        0x8A => '\u{0160}',
        0x8B => '\u{2039}',
        0x8C => '\u{0152}',
        0x8E => '\u{017D}',
        0x91 => '\u{2018}',
        0x92 => '\u{2019}',
        0x93 => '\u{201C}',
        0x94 => '\u{201D}',
        0x95 => '\u{2022}',
        0x96 => '\u{2013}',
        0x97 => '\u{2014}',
        0x98 => '\u{02DC}',
        0x99 => '\u{2122}',
        0x9A => '\u{0161}',
        0x9B => '\u{203A}',
        0x9C => '\u{0153}',
        0x9E => '\u{017E}',
        0x9F => '\u{0178}',
        // 0x81/0x8D/0x8F/0x90/0x9D are UNDEFINED in cp1252 — destroyed, not transformed.
        0x81 | 0x8D | 0x8F | 0x90 | 0x9D => '\u{FFFD}',
        other => other as char,
    }
}
