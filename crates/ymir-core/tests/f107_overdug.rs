//! ADR 0001 Finding 107 — the "canyon" gate becomes TOPOLOGICAL. **Author, 2026-09-20.**
//!
//! The class Findings 87–104 scored is the class of **OVER-DUG DEPRESSIONS** (Finding 87's
//! *"canyon noyé"*): a body ≥ 1 km² whose floor lies more than 50 m **below its own sill**, with
//! walls at p50 > 30°. The instrument scored something else — the **CUT**, `pre-incision −
//! delivered` — and a cut cannot tell an over-dug hollow from a deep valley the drainage runs
//! straight through. The day an uplift or a glaciation digs a real gorge, a cut-based gate kills
//! it; a topological one does not.
//!
//! ⇒ **The discriminant is `filled − raw`**, the priority flood (Findings 13/14) on the ERODED
//! field, already computed in the chain. A through-cut bed has `filled == raw` along its length
//! whatever its depth, so it scores 0 and is NOT in the class.
//!
//! Every table now carries **two columns**: **cut** (`pre-incision − delivered`) and **fill**
//! (`filled − raw`). `klass` is the over-dug count; `cut_class` is the number every historical
//! Finding printed.
//!
//! **The two controls the author set:**
//! 1. Finding 88's **sixteen** bodies — every one of them must read fill > 50 m. If the
//!    redefinition dropped them, it is not a sharpening, it is a different gate.
//! 2. **A river bed must read fill ≈ 0** — Finding 87-A's biggest `Spillway` (92.97 m³/s at
//!    Finding 87's render table, ADR L10863). This is the control the whole redefinition is for.
//!
//! ⚠️ **Rule 11 / 11b, recorded before the first measurement.** `over-dug`, `sur-creusé`,
//! `canyon noyé` (as `noy`), `filled − raw`, `raw − filled`, `own sill`, `propre col`,
//! `topological` in this sense: **NOTHING FOUND** — the vocabulary is new to the dossier.
//! `drowned` 60 hits, earliest L1306 (Finding 35) — a DIFFERENT sense, submerged by a lake, not
//! over-dug. `priority flood` 28 hits, earliest **L673 Finding 13**, read first: it is the
//! conditioning this gate now rests on, and **L14265 (Finding 96) already measured the
//! population** — *"pits (cells the priority flood must raise): χ 0 · delivered 1 240 677"*.
//! `92.97` 11 hits, earliest **L10863**; `115 m` 2 hits, **neither is the river** (L1715 is a
//! lake altitude, L1818 a table cell) ⇒ the control river is found by its DISCHARGE, not its drop.
//!
//! Run: cargo test -p ymir-core --release --test f107_overdug -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, f95_criteria, set_dump};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const DELIVERED_P50_M: f32 = 424.2;
/// Finding 106's window, all three seeds: the rows whose canyon column is now re-read.
const WINDOW: [f32; 3] = [0.021, 0.024, 0.027];

#[test]
#[ignore]
fn f107_overdug() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 107 . the over-dug class, and the river it must spare  =====");

    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);

    let run = |f: &ymir_core::grid::GridF32, label: &str, dump: bool| {
        let dr = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(f, &dr.flow.filled, &dr.lake_map, SEA, w, h);
        eprintln!("\n--------  {label}  --------");
        set_dump(dump);
        let cr = f95_criteria(f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        set_dump(false);
        let rate = |k: usize| 100.0 * k as f32 / cr.scanned.max(1) as f32;
        eprintln!(
            "   ⇒ **{label}**: over-dug **{} of {}** = **{:.1} %** · old cut class **{} of {}** = \
             **{:.1} %** · **Δ {:+}**",
            cr.klass,
            cr.scanned,
            rate(cr.klass),
            cr.cut_class,
            cr.scanned,
            rate(cr.cut_class),
            cr.klass as i64 - cr.cut_class as i64
        );
        (cr.klass, cr.cut_class, cr.scanned)
    };

    // ── CONTROLS 1 and 2: the delivered field, where Finding 88 found its sixteen ──
    let delivered = build_field_seed(Knobs::passes(2), PSEED);
    let (dug, cut, scanned) = run(&delivered, "DELIVERED (seed 1) — Finding 88's sixteen", true);
    eprintln!(
        "\n   ⛔ **CONTROL 1**: Finding 88 / 102 read **16 of 54 = 29.6 %** on this field. Here the \
         cut class is **{cut} of {scanned}** and the over-dug class is **{dug}** ⇒ **{}**",
        // ⚠️ `cut == dug` is NOT the control. The control is SET INCLUSION -- every cut-class
        // body must also be over-dug -- and the counts can differ while it holds, which is what
        // happens here (16 of 16 kept, 2 added). Equality of counts would even be the WEAKER
        // reading: two classes of 16 could be disjoint. The per-body lines carry the inclusion;
        // this line reports the direction so nobody reads a widening as a failure.
        match dug.cmp(&cut) {
            std::cmp::Ordering::Equal =>
                "the counts match -- check the per-body lines for INCLUSION",
            std::cmp::Ordering::Greater => {
                "a WIDENING -- the per-body lines say whether the sixteen are all kept"
            }
            std::cmp::Ordering::Less => "**a NARROWING -- some of the sixteen were dropped**",
        }
    );

    // ── the Finding 106 window, both columns ──────────────────────────────
    for c in WINDOW {
        let f = build_field_seed(
            Knobs { slope_floor_uk: Some(c), depression_floor: true, ..Knobs::passes(2) },
            PSEED,
        );
        // ⚠️ **Rule 16** (Finding 104): a global statistic that decides admissibility must be
        // re-read on the population where the defect lives. The over-dug rate moves this window
        // from three constants to one, on a change of **2 bodies -> 3 out of 29** -- one body
        // decides a gate, exactly the small-number regime Finding 101 flagged. The dump is ON so
        // the deciding body is NAMED, not inferred from a rate.
        run(&f, &format!("c = {c} (Finding 106's window)"), true);
    }

    eprintln!("\n==========  end Finding 107 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// ADR Finding 107-B — does the {0.027} collapse TRANSFER, or is it seed 1's?
///
/// The main test re-reads Finding 106's window on seed 1 and finds the gate `rate <= 0.25 x
/// delivered` closes on two of its three constants once the delivered denominator is the over-dug
/// rate (33.3 %, not 29.6 %). **That claim is worth nothing on one seed** -- Finding 102 spent a
/// whole round establishing that a window read off seed 1 is a calibration, not a rule, and
/// Finding 103 then found the intersection empty on the other two. So both halves are measured
/// here on seeds 2 and 3: the delivered over-dug rate, which SETS the gate, and the three
/// candidate rows it judges.
///
/// Run: cargo test -p ymir-core --release --test f107_overdug -- --ignored transfer --nocapture
#[test]
#[ignore]
fn f107_transfer() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let t0 = Instant::now();
    let s2 = splitmix64(PSEED);
    let s3 = splitmix64(s2);
    eprintln!(
        "
==========  Finding 107-B . does the collapse transfer?  =========="
    );

    for (si, seed) in [s2, s3].into_iter().enumerate() {
        eprintln!(
            "
========  SEED {} ({seed})  ========",
            si + 2
        );
        let pre = build_field_seed(Knobs::no_incision(), seed);
        let (w, h) = (pre.width, pre.height);
        let mut row = |f: &ymir_core::grid::GridF32, label: &str| -> (f32, f32) {
            let dr = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
            let bf = breach_monotone(f, &dr.flow.filled, &dr.lake_map, SEA, w, h);
            let cr = f95_criteria(f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
            let d = cr.scanned.max(1) as f32;
            let (dug, cut) = (100.0 * cr.klass as f32 / d, 100.0 * cr.cut_class as f32 / d);
            eprintln!(
                "   **{label}**: over-dug **{} of {}** = **{dug:.1} %** · cut **{} of {}** =                  **{cut:.1} %** · Δ **{:+}**",
                cr.klass,
                cr.scanned,
                cr.cut_class,
                cr.scanned,
                cr.klass as i64 - cr.cut_class as i64
            );
            (dug, cut)
        };

        // The delivered field SETS the gate -- it is the denominator, so it is measured, never cited.
        let delivered = build_field_seed(Knobs::passes(2), seed);
        let (d_dug, d_cut) = row(&delivered, "DELIVERED");
        eprintln!(
            "   ⇒ gate `0.25 x delivered`: under **cut {:.1} %**, under **over-dug {:.1} %**",
            0.25 * d_cut,
            0.25 * d_dug
        );

        for c in WINDOW {
            let f = build_field_seed(
                Knobs { slope_floor_uk: Some(c), depression_floor: true, ..Knobs::passes(2) },
                seed,
            );
            let (dug, cut) = row(&f, &format!("c = {c}"));
            eprintln!(
                "      ⇒ **c = {c}**: under cut **{}** · under over-dug **{}**",
                if cut <= 0.25 * d_cut + 0.001 { "PASS" } else { "**FAIL**" },
                if dug <= 0.25 * d_dug + 0.001 { "PASS" } else { "**FAIL**" }
            );
        }
    }
    eprintln!(
        "
==========  end Finding 107-B . {:.1} s  ==========
",
        t0.elapsed().as_secs_f64()
    );
}
