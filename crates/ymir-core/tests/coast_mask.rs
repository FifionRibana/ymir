//! ADR Finding 79 blocks B and C — the fringe defined on the MASK, and the pre-clamp depths.
//!
//! **B1, the instrument, and it needs no new code.** `marching_squares` on a BINARY field
//! (land = 1.0, sea = 0.0) at level 0.5 interpolates to the exact midpoint of every crossed
//! edge, so the traced path is the 8-connected land/sea boundary staircase and nothing else.
//! It cannot see the value of the water. The same `coast_spurs` / `coast_shape_thresholds`
//! then run on it, so the mask column and the isoline column are the same detector on two
//! definitions of "the coast" — which is what the open question to the author needs.
//!
//! Run: cargo test -p ymir-core --release --test coast_mask -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted, spectrum};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds};
use ymir_core::terrain::contour::marching_squares;

/// The land/sea mask as a field the tracer can read: 1.0 on land, 0.0 on sea. Tracing it at
/// 0.5 gives the boundary staircase — the coast a consumer that draws from cell tags would see.
fn mask_field(f: &GridF32) -> GridF32 {
    let mut m = GridF32::new(f.width, f.height, 0.0);
    for k in 0..f.data.len() {
        if f.data[k] > SEA {
            m.data[k] = 1.0;
        }
    }
    m
}

struct Six {
    km: usize,
    cell: usize,
    coast: f32,
    p90: f32,
    r: f32,
}

fn six(polys: &[Vec<(f32, f32)>]) -> Six {
    let a = coast_shape_thresholds(polys, CELL_KM, 1.0, NECK_KM);
    let b = coast_shape_thresholds(polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    Six { km: a.count, cell: b.count, coast: a.coast_km, p90: a.p90_len_km, r: b.local_axis_r }
}

fn line(name: &str, iso: &Six, msk: &Six, base: Option<(&Six, &Six)>) {
    let d = |v: usize, b: Option<usize>| match b {
        Some(b) => format!("{:+}", v as i64 - b as i64),
        None => "—".into(),
    };
    eprintln!(
        "   {name:<26} | ISO {:>6} {:>6} {:>7.0} {:>6.2} {:>6.3} | MASK {:>6} {:>6} {:>7.0} \
         {:>6.2} {:>6.3} | Δmask {:>8} {:>8}",
        iso.km,
        iso.cell,
        iso.coast,
        iso.p90,
        iso.r,
        msk.km,
        msk.cell,
        msk.coast,
        msk.p90,
        msk.r,
        d(msk.km, base.map(|b| b.1.km)),
        d(msk.cell, base.map(|b| b.1.cell))
    );
}

#[test]
#[ignore]
fn coast_mask() {
    let ss =
        ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!(
        "\n==========  Finding 79 B/C · the fringe on the MASK, and the pre-clamp field  =========="
    );
    eprintln!(
        "   ISO = marching squares on the height field at sea level (what an interpolating \
         consumer draws).\n   MASK = marching squares on the BINARY land/sea field at 0.5 (what \
         a tag-drawing consumer draws). Same detector, two definitions."
    );

    // ── B1 · the control: the mask must not move when only the water moves ────
    // Finding 78 measured the land/sea mask unchanged across the shelf sweep, so a mask-based
    // count MUST be identical at 20 m and at 1 m. If it is not, the instrument still reads the
    // sea and block B2 is void.
    let f20 = build_field(Knobs { shelf_min_depth_m: Some(20.0), ..Knobs::shipped() });
    let f01 = build_field(Knobs { shelf_min_depth_m: Some(1.0), ..Knobs::shipped() });
    let m20 = six(&marching_squares(&mask_field(&f20), 0.5));
    let m01 = six(&marching_squares(&mask_field(&f01), 0.5));
    let i20 = six(&marching_squares(&f20, SEA));
    let i01 = six(&marching_squares(&f01, SEA));
    eprintln!("\n── B1 · CONTROL — shelf 20 m against shelf 1 m ──");
    eprintln!(
        "   MASK   20 m: {} / {} / {:.0} km   |   1 m: {} / {} / {:.0} km   ⇒ {}",
        m20.km,
        m20.cell,
        m20.coast,
        m01.km,
        m01.cell,
        m01.coast,
        if m20.km == m01.km && m20.cell == m01.cell {
            "IDENTICAL — the instrument does not read the sea"
        } else {
            "DIFFERENT — the instrument still reads the sea, B2 is VOID"
        }
    );
    eprintln!(
        "   ISO    20 m: {} / {} / {:.0} km   |   1 m: {} / {} / {:.0} km   ⇒ the isoline moves \
         by {:+} / {:+} ({:+.1} % on the cell scale)",
        i20.km,
        i20.cell,
        i20.coast,
        i01.km,
        i01.cell,
        i01.coast,
        i01.km as i64 - i20.km as i64,
        i01.cell as i64 - i20.cell as i64,
        100.0 * (i01.cell as f64 - i20.cell as f64) / i20.cell as f64
    );
    assert_eq!(m20.km, m01.km, "the mask detector moved with the water: B2 would be void");

    // ── B2 · the re-read of Findings 75 to 78, mask beside isoline ────────────
    eprintln!("\n── B2 · Findings 75-78 re-read on both definitions ──");
    eprintln!(
        "   {:<26} | {:^34} | {:^34} | {:^17}",
        "variant", "ISOLINE  km cell coast p90 R", "MASK  km cell coast p90 R", "Δ MASK vs pre-inc"
    );
    let pre = build_field(Knobs::no_incision());
    let pi = six(&marching_squares(&pre, SEA));
    let pm = six(&marching_squares(&mask_field(&pre), 0.5));
    line("PRE-INCISION (authority)", &pi, &pm, None);
    let base = Some((&pi, &pm));
    line("SHIPPED (shelf 20)", &i20, &m20, base);
    line("PRODUCTION (shelf 1)", &i01, &m01, base);
    let mut extra: Vec<(&str, GridF32)> = Vec::new();
    for (nm, kn) in [
        (
            "diffusion 1.28 (shelf 20)",
            Knobs {
                shelf_min_depth_m: Some(20.0),
                diffusion: Some(1.28),
                diffusion_substeps: Some(8),
                ..Knobs::shipped()
            },
        ),
        (
            "diffusion 1.28 (shelf 1)",
            Knobs {
                shelf_min_depth_m: Some(1.0),
                diffusion: Some(1.28),
                diffusion_substeps: Some(8),
                ..Knobs::shipped()
            },
        ),
    ] {
        let f = build_field(kn);
        let i = six(&marching_squares(&f, SEA));
        let m = six(&marching_squares(&mask_field(&f), 0.5));
        line(nm, &i, &m, base);
        extra.push((nm, f));
    }

    // ── B3 · the spectrum, on the mask ────────────────────────────────────────
    eprintln!("\n── B3 · the spectrum on the MASK (Finding 76-C method) ──");
    eprintln!(
        "   {:<26} {:>9} {:>8} {:>12} {:>9} {:>9}",
        "variant", "segments", "gaps", "median gap", "lambda", "x white"
    );
    for (nm, polys) in [
        ("PRE-INCISION mask", marching_squares(&mask_field(&pre), 0.5)),
        ("SHIPPED mask (shelf 20)", marching_squares(&mask_field(&f20), 0.5)),
        ("SHIPPED isoline (shelf 20)", marching_squares(&f20, SEA)),
        ("PRODUCTION mask (shelf 1)", marching_squares(&mask_field(&f01), 0.5)),
    ] {
        match spectrum(&polys) {
            Some((s, g, med, p, r)) => eprintln!(
                "   {nm:<26} {s:>9} {g:>8} {med:>12.1} {p:>9} {r:>9.2}{}",
                if r >= 3.0 { "  ⇒ wavelength" } else { "" }
            ),
            None => eprintln!("   {nm:<26} {:>9} — rule 10: under 32 gaps, NOT a reading", "—"),
        }
    }

    // ── C · the PRE-CLAMP depths of the drowned cells ─────────────────────────
    let noclamp = build_field(Knobs { bathymetry_off: true, ..Knobs::shipped() });
    let n = noclamp.data.len();
    let drowned: Vec<usize> =
        (0..n).filter(|&k| pre.data[k] > SEA && noclamp.data[k] <= SEA).collect();
    let d = sorted(drowned.iter().map(|&k| -(noclamp.data[k] - SEA) * n2m).collect());
    eprintln!(
        "\n── C · depth of the DROWNED cells on the PRE-CLAMP field (population: pre-incision \
         land AND pre-clamp sea) ──\n   {} cells | p1 {:.4} p10 {:.4} MEDIAN **{:.4}** p90 \
         {:.4} p99 {:.4} max {:.2} m below sea",
        d.len(),
        pct(&d, 0.01),
        pct(&d, 0.10),
        pct(&d, 0.50),
        pct(&d, 0.90),
        pct(&d, 0.99),
        d.last().copied().unwrap_or(f32::NAN)
    );
    for t in [0.05f32, 0.2, 0.5, 1.0, 2.0, 5.0] {
        eprintln!(
            "     shallower than {t:>5} m : {:>7} cells ({:>5.1} %)",
            d.iter().filter(|&&x| x <= t).count(),
            100.0 * d.iter().filter(|&&x| x <= t).count() as f64 / d.len() as f64
        );
    }

    // ── C-bis · the NON-DEGENERATE D: the Δ(threshold) curve ──────────────────
    eprintln!("\n── C-bis · D, non-degenerate: lift drowned cells shallower than a threshold ──");
    eprintln!(
        "   {:<22} {:>9} {:>9} {:>9} {:>11} {:>11}",
        "threshold", "lifted", "ISO km", "ISO cell", "Δ mask km", "Δ mask cell"
    );
    for t in [0.0f32, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 1e9] {
        let mut g = noclamp.clone();
        let mut lifted = 0usize;
        for &k in &drowned {
            if -(noclamp.data[k] - SEA) * n2m <= t {
                g.data[k] = SEA + 0.01 / n2m;
                lifted += 1;
            }
        }
        let gi = six(&marching_squares(&g, SEA));
        let gm = six(&marching_squares(&mask_field(&g), 0.5));
        eprintln!(
            "   {:<22} {lifted:>9} {:>9} {:>9} {:>+11} {:>+11}",
            if t > 1e8 { "ALL (reference)".to_string() } else { format!("{t} m") },
            gi.km,
            gi.cell,
            gm.km as i64 - pm.km as i64,
            gm.cell as i64 - pm.cell as i64
        );
    }
    eprintln!("   (threshold 0 m is the negative control: it lifts nothing and must not move Δ.)");
}
