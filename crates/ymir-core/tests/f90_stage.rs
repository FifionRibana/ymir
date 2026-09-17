//! ADR 0001 Finding 90 — the bisection of Finding 89-C5's Δ(≥ 1 km) = +115, BY STAGE.
//!
//! Three stages, named everywhere: **ERODED** (the incision's output), **BREACHED**
//! (`breach_monotone`, which `hd_assembly` delivers) and **u16** (the export quantisation).
//! Finding 89-C5 compared a BREACHED delivered mask against an UNBREACHED pre-incision one; this
//! file measures every combination so the +115 can be attributed to a stage rather than argued.
//!
//! Rule 12: no walk here reads `flow.direction`; the fields are declared by their construction.
//!
//! Run: cargo test -p ymir-core --release --test f90_stage -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, land_u16, majority, pct, sorted, to_mask};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, Spur, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;

/// Spur counts at BOTH length scales, on one boolean mask at one resolution.
fn spurs_at(mask: &[bool], w: usize, factor: usize) -> (usize, usize) {
    let (m, ww) = majority(mask, w, factor);
    let ck = CELL_KM * factor as f32;
    let polys = marching_squares(&to_mask(&m, ww), 0.5);
    let (km, _) = coast_spurs(&polys, ck, MIN_SPUR_KM, NECK_KM);
    let (cells, _) = coast_spurs(&polys, ck, 2.0 * ck, NECK_KM);
    (km.len(), cells.len())
}

/// The ≥ 1 km spur LIST at 8192², for the A3/A4 attribution.
fn spur_list(mask: &[bool], w: usize) -> Vec<Spur> {
    let polys = marching_squares(&to_mask(mask, w), 0.5);
    coast_spurs(&polys, CELL_KM, MIN_SPUR_KM, NECK_KM).0
}

fn median(v: &[f32]) -> f32 {
    pct(&sorted(v.to_vec()), 0.50)
}

#[test]
#[ignore]
fn f90_stage() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!(
        "\n==========  Finding 90 · A (the stage) · B (the breach under the sea)  =========="
    );
    eprintln!(
        "   MIN_SPUR_KM = {MIN_SPUR_KM} km = **{:.1} cells at 8192²**, {:.1} at 4096², {:.1} at \\
         2048², {:.1} at 1024² (NECK_KM {NECK_KM})",
        MIN_SPUR_KM / CELL_KM,
        MIN_SPUR_KM / (CELL_KM * 2.0),
        MIN_SPUR_KM / (CELL_KM * 4.0),
        MIN_SPUR_KM / (CELL_KM * 8.0)
    );

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let mut off = on.clone();
    off.merged_union_relevel = None;

    let raw = build_field(Knobs::passes(2));
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let pre = build_field(Knobs::no_incision());
    let raw3 = build_field(Knobs::passes(3));
    let raw_nb = build_field(Knobs { base_level_off: true, ..Knobs::passes(2) });
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);

    // one helper, one convention: breach a field with its OWN drainage, as production does
    let breach_of = |g: &GridF32, cfg: &C1DrainageConfig| -> GridF32 {
        let d = c1_drainage_windowed(g, None, cfg, &ss, DOMAIN_KM);
        breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, w, h)
    };

    let br_raw = breach_of(&raw, &on);
    let br_pre = breach_of(&pre, &on);

    let lm = |g: &GridF32| -> Vec<bool> { (0..n).map(|k| g.data[k] > SEA).collect() };
    let stages: [(&str, Vec<bool>, Vec<bool>); 3] = [
        ("ERODED   (incision out, f32)", lm(&raw), lm(&pre)),
        ("BREACHED (hd_assembly, f32)", lm(&br_raw), lm(&br_pre)),
        ("u16      (the export mask)", land_u16(&br_raw, &ss), land_u16(&br_pre, &ss)),
    ];

    // ══ A0 — the six lines, stage-matched ═══════════════════════════════════
    eprintln!("\n── A0 · Δ(delivered − pre-incision) AT THE SAME STAGE, four resolutions ──");
    for (label, del, pr) in &stages {
        let mut k_line = String::new();
        let mut c_line = String::new();
        for f in [1usize, 2, 4, 8] {
            let (kd, cd) = spurs_at(del, w, f);
            let (kp, cp) = spurs_at(pr, w, f);
            k_line.push_str(&format!(" {:>5}/{:<5}Δ{:<+6}", kd, kp, kd as i64 - kp as i64));
            c_line.push_str(&format!(" {:>5}/{:<5}Δ{:<+6}", cd, cp, cd as i64 - cp as i64));
        }
        eprintln!("   {label}");
        eprintln!("      ≥ 1 km   (8192 / 4096 / 2048 / 1024):{k_line}");
        eprintln!("      ≥ 2 cell (8192 / 4096 / 2048 / 1024):{c_line}");
    }
    // the line Finding 89-C5 actually published: breached delivered vs UNBREACHED pre-incision
    eprintln!(
        "\n   **the Finding 89-C5 line, reproduced (MIXED stages: u16 of the BREACHED \\
               delivered against u16 of the UNBREACHED pre-incision)**"
    );
    {
        let du = land_u16(&br_raw, &ss);
        let pu = land_u16(&pre, &ss);
        let mut l1 = String::new();
        let mut l2 = String::new();
        for f in [1usize, 2, 4, 8] {
            let (kd, cd) = spurs_at(&du, w, f);
            let (kp, cp) = spurs_at(&pu, w, f);
            l1.push_str(&format!(" {:>5}/{:<5}Δ{:<+6}", kd, kp, kd as i64 - kp as i64));
            l2.push_str(&format!(" {:>5}/{:<5}Δ{:<+6}", cd, cp, cd as i64 - cp as i64));
        }
        eprintln!("      ≥ 1 km  :{l1}");
        eprintln!("      ≥ 2 cell:{l2}");
    }

    // ══ A2 — the two knobs that span Finding 80 → 83 → 86, one variable each ══
    eprintln!("\n── A2 · the same line under the two knobs that changed between F80 and now ──");
    let br_nb = breach_of(&raw_nb, &on);
    let br_raw_off = breach_of(&raw, &off);
    for (label, del) in [
        ("shipped (bound ON, relevel ON)", &br_raw),
        ("bound OFF (the pre-F83 world)", &br_nb),
        ("relevel OFF (the pre-F86 world)", &br_raw_off),
    ] {
        let du = land_u16(del, &ss);
        let pu = land_u16(&br_pre, &ss);
        let (kd, cd) = spurs_at(&du, w, 1);
        let (kp, cp) = spurs_at(&pu, w, 1);
        eprintln!(
            "      {label:<32} u16 8192²: ≥ 1 km **{kd} (Δ {:+})** · ≥ 2 cells {cd} (Δ {:+})",
            kd as i64 - kp as i64,
            cd as i64 - cp as i64
        );
    }

    // ══ B1 — does the breach put land under the sea? ════════════════════════
    eprintln!("\n── B1 · cells at or below sea level, ERODED against BREACHED ──");
    let wc_e = water_class(&raw, SEA);
    let wc_b = water_class(&br_raw, SEA);
    let count = |g: &GridF32| (0..n).filter(|&k| g.data[k] <= SEA).count();
    let drowned: Vec<usize> =
        (0..n).filter(|&k| raw.data[k] > SEA && br_raw.data[k] <= SEA).collect();
    let raised: usize = (0..n).filter(|&k| raw.data[k] <= SEA && br_raw.data[k] > SEA).count();
    let depth: Vec<f32> = drowned.iter().map(|&k| -m(&br_raw, k)).collect();
    let enclosed_e = (0..n).filter(|&k| wc_e[k] == 2).count();
    let enclosed_b = (0..n).filter(|&k| wc_b[k] == 2).count();
    let enclosed_new = (0..n).filter(|&k| wc_b[k] == 2 && wc_e[k] != 2).count();
    eprintln!(
        "      h ≤ sea: eroded {} → breached {} (**{:+}**) · **the breach DROWNS {} cells** \\
         ({:.2} km²) and raises {raised}",
        count(&raw),
        count(&br_raw),
        count(&br_raw) as i64 - count(&raw) as i64,
        drowned.len(),
        drowned.len() as f32 * cell_km2
    );
    if !drowned.is_empty() {
        let ds = sorted(depth.clone());
        eprintln!(
            "      their depth below sea: p50 **{:.2} m** · p90 {:.2} m · max **{:.2} m** \\
             (Finding 89-D's −8.67 sits at the {:.1}th centile)",
            pct(&ds, 0.50),
            pct(&ds, 0.90),
            ds.last().copied().unwrap_or(0.0),
            100.0 * ds.iter().filter(|&&d| d < 8.67).count() as f64 / ds.len() as f64
        );
    }
    eprintln!(
        "      `wc == 2` (ENCLOSED below-sea): eroded {enclosed_e} → breached {enclosed_b} · \\
         **{enclosed_new} cells the breach makes enclosed-below-sea = {:.1} % of the breached \\
         total**",
        100.0 * enclosed_new as f64 / enclosed_b.max(1) as f64
    );

    // ══ B3 — the same count at three passes ═════════════════════════════════
    let br_raw3 = breach_of(&raw3, &on);
    let drowned3: usize = (0..n).filter(|&k| raw3.data[k] > SEA && br_raw3.data[k] <= SEA).count();
    let wc_b3 = water_class(&br_raw3, SEA);
    let wc_e3 = water_class(&raw3, SEA);
    let enc_new3 = (0..n).filter(|&k| wc_b3[k] == 2 && wc_e3[k] != 2).count();
    eprintln!(
        "\n── B3 · at `iterations = 3`: the breach drowns **{drowned3}** cells (against \\
         {} at 2 passes, {:+.1} %) and makes **{enc_new3}** enclosed-below-sea (against \\
         {enclosed_new})",
        drowned.len(),
        100.0 * (drowned3 as f64 - drowned.len() as f64) / drowned.len().max(1) as f64
    );

    // ══ A3 — WHERE the added spurs are, and are they the drowned cells? ═════
    eprintln!("\n── A3 · the added ≥ 1 km spurs on the u16 mask, named ──");
    let du = land_u16(&br_raw, &ss);
    let pu = land_u16(&pre, &ss); // the AUTHORITY of Findings 80/83: unbreached pre-incision
    let sd = spur_list(&du, w);
    let sp = spur_list(&pu, w);
    // an added spur: no pre-incision spur whose neck midpoint is within 20 cells (1 km)
    let near = 20.0f32;
    let added: Vec<&Spur> = sd
        .iter()
        .filter(|a| {
            !sp.iter().any(|b| {
                let (dx, dy) = (a.mid.0 - b.mid.0, a.mid.1 - b.mid.1);
                (dx * dx + dy * dy).sqrt() < near
            })
        })
        .collect();
    let drowned_set: std::collections::HashSet<usize> = drowned.iter().copied().collect();
    let mut with_drowned = 0usize;
    let mut lens: Vec<f32> = Vec::new();
    for a in &added {
        lens.push(a.len_km / CELL_KM); // length in 8192² cells
        let (cx, cy) = (a.mid.0.round() as i32, a.mid.1.round() as i32);
        let mut found = false;
        'scan: for dy in -4i32..=4 {
            for dx in -4i32..=4 {
                let x = (cx + dx).rem_euclid(w as i32) as usize;
                let y = (cy + dy).rem_euclid(h as i32) as usize;
                if drowned_set.contains(&(y * w + x)) {
                    found = true;
                    break 'scan;
                }
            }
        }
        if found {
            with_drowned += 1;
        }
    }
    eprintln!(
        "      delivered {} spurs ≥ 1 km · pre-incision authority {} · **{} ADDED** (no \\
         counterpart within 1 km)",
        sd.len(),
        sp.len(),
        added.len()
    );
    eprintln!(
        "      of the added, **{with_drowned} ({:.1} %) have a breach-drowned cell within 4 \\
         cells of their neck midpoint**",
        100.0 * with_drowned as f64 / added.len().max(1) as f64
    );
    // ══ A4 — their length in cells, against the threshold at each grid ══════
    if !lens.is_empty() {
        let ls = sorted(lens.clone());
        eprintln!(
            "\n── A4 · the added spurs' length in 8192² cells: p10 **{:.1}** · p50 **{:.1}** · \\
             p90 **{:.1}** · max {:.1} — the ≥ 1 km threshold is 20.5 cells at 8192², 41 at \\
             4096², 82 at 2048², 164 at 1024² (in 8192² cells)",
            pct(&ls, 0.10),
            pct(&ls, 0.50),
            pct(&ls, 0.90),
            ls.last().copied().unwrap_or(0.0)
        );
        for (grid, thr) in [("4096", 41.0f32), ("2048", 82.0), ("1024", 164.0)] {
            eprintln!(
                "      at {grid}²: **{:.1} %** of the added spurs are shorter than the threshold \\
                 and cannot be counted there by construction",
                100.0 * ls.iter().filter(|&&l| l < thr).count() as f64 / ls.len() as f64
            );
        }
        eprintln!("      median added spur length {:.2} km", median(&lens) * CELL_KM);
    }
    eprintln!("\n==========  end Finding 90 · A · B  ==========\n");
}
