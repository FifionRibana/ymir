//! ADR 0001 Finding 83 blocks B0/B1 — **is the coast ORGANIC?**, and the instrument that decides.
//!
//! B0 fixes the criterion before anything is measured: an organic coast has NO wavelength and NO
//! visible generator. The instrument is the **Richardson curve** — coastline length against the
//! length of the ruler used to measure it — read on the CONSUMER mask (post-u16 land/sea,
//! `marching_squares` at 0.5), because that is the coast Living Landz draws.
//!
//! **Compass (divider), not box-counting, and the reason is not taste.** The divider returns a
//! LENGTH at each ruler, which is the quantity already on file (1 634 / 7 056 / 1 639 km,
//! Finding 80-B1) and therefore directly comparable; box-counting returns a dimension for the
//! point set and is contaminated by the one-pixel thickness of the staircase, which on a mask
//! contour is an artefact of the encoding rather than a property of the shore. The divider is
//! also Richardson's own instrument and the one Mandelbrot's D is defined from.
//!
//! **The instrument has a floor and it is measured, not assumed**: a digitised smooth curve is a
//! staircase, and a staircase can inflate the length at the smallest rulers. So the same pipeline
//! is run on a rasterised CIRCLE (true D = 1.000) and on a rasterised KOCH snowflake (true
//! D = log 4 / log 3 = 1.2619). Any D from the real coast is read against those two.
//!
//! Run: cargo test -p ymir-core --release --test f83_richardson -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, DOMAIN_KM, Knobs, RICHARDSON_MIN_POLY_KM, RICHARDSON_OFFSETS, RICHARDSON_RULERS, SEA,
    build_field, circle_polygon, fill_polygon, koch_polygon, land_u16, majority, richardson,
    richardson_fit, richardson_line, spectrum, to_mask,
};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds};
use ymir_core::terrain::contour::marching_squares;

// ─────────────────────────────────────────────────────────────────────────────
// The two synthetic controls
// ─────────────────────────────────────────────────────────────────────────────

/// D of a synthetic mask, on the whole-population curve.
fn synth_d(land: &[bool], w: usize) -> f64 {
    let polys = marching_squares(&to_mask(land, w), 0.5);
    let keep: Vec<&Vec<(f32, f32)>> = polys.iter().collect();
    richardson_fit(&richardson(&keep, CELL_KM, 10.0)).0
}

#[test]
#[ignore]
fn f83_richardson() {
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 83 · B0/B1 — the Richardson curve  ==========");
    eprintln!(
        "\nB0, THE CRITERION, written before the first number:\n   An ORGANIC coast has no \
         wavelength and no visible generator. It is met when, on the CONSUMER mask,\n   (1) \
         log L against log(ruler) is LINEAR over 100 m - 10 km (residual declared),\n   (2) the \
         slope gives D in [1.10 ; 1.30] (Great Britain ~1.25; a smooth coast 1.00),\n   (3) the \
         spectrum of excursion positions has no peak above 3x white,\n   (4) Delta(spurs >= 2 \
         cells) against pre-incision stays 0 for everything that is not an estuary,\n   (5) and \
         the author's eye, in Living Landz, with the counts beside it.\n   Instrument: DIVIDER \
         (compass), {RICHARDSON_OFFSETS} start offsets, {RICHARDSON_RULERS} log-spaced rulers, \
         polygons >= {RICHARDSON_MIN_POLY_KM} km."
    );

    // ── the instrument's own floor and its own ceiling ────────────────────────
    eprintln!("\n-- B1 control · the instrument on two shapes of KNOWN dimension (8192 grid) --");
    let w = 8192usize;
    let circ = fill_polygon(&circle_polygon(w, 3000.0), w);
    let kch = fill_polygon(&koch_polygon(w, 5400.0, 7), w);
    richardson_line("CIRCLE (true D = 1.000)", &circ, w, CELL_KM, 10.0);
    richardson_line("KOCH g7 (true D = 1.2619)", &kch, w, CELL_KM, 10.0);
    let dc = synth_d(&circ, w);
    let dk = synth_d(&kch, w);
    eprintln!(
        "   => instrument FLOOR {dc:.3} (bias {:+.3}) · on a known fractal {dk:.3} (bias \
         {:+.3}). Every D below is read against these.",
        dc - 1.0,
        dk - 1.2619
    );

    // ── the three worlds ──────────────────────────────────────────────────────
    let pre = build_field(Knobs::no_incision());
    let del = build_field(Knobs::pre83());
    let bnd = build_field(Knobs::shipped());

    for (gtag, factor, cell_km) in [
        ("u16 8192", 1usize, CELL_KM),
        ("terrain 2048", 4, CELL_KM * 4.0),
        ("ocean 1024", 8, CELL_KM * 8.0),
    ] {
        eprintln!("\n-- B1 · {gtag} (cell {cell_km:.4} km, rulers {:.3}-10 km) --", 2.0 * cell_km);
        for (tag, f) in
            [("PRE-INCISION", &pre), ("DELIVERED (pre-83)", &del), ("PRODUCTION (bounded)", &bnd)]
        {
            let l8 = land_u16(f, &ss);
            let (land, gw) = if factor == 1 { (l8, w) } else { majority(&l8, w, factor) };
            richardson_line(tag, &land, gw, cell_km, 10.0);
        }
    }

    // ── the other three legs of the criterion, on the 8192 u16 mask ──────────
    eprintln!("\n-- B1 · spurs and spectrum on the u16 8192 mask (criterion legs 3 and 4) --");
    eprintln!(
        "   {:<22} {:>10} {:>12} {:>10} {:>12} {:>12}",
        "variant", ">= 1 km", ">= 2 cells", "coast km", "spur peak", "x white"
    );
    let mut base = (0usize, 0usize);
    for (tag, f) in
        [("PRE-INCISION", &pre), ("DELIVERED (pre-83)", &del), ("PRODUCTION (bounded)", &bnd)]
    {
        let polys = marching_squares(&to_mask(&land_u16(f, &ss), w), 0.5);
        let a = coast_shape_thresholds(&polys, CELL_KM, 1.0, NECK_KM);
        let b = coast_shape_thresholds(&polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
        let sp = spectrum(&polys);
        if tag == "PRE-INCISION" {
            base = (a.count, b.count);
        }
        eprintln!(
            "   {tag:<22} {:>10} {:>12} {:>10.0} {:>12} {:>12}   d>=1km {:+} · **d>=2c {:+}**",
            a.count,
            b.count,
            a.coast_km,
            sp.map_or("n/a".into(), |s| format!("{}", s.3)),
            sp.map_or("n/a".into(), |s| format!("{:.2}", s.4)),
            a.count as i64 - base.0 as i64,
            b.count as i64 - base.1 as i64
        );
    }
    eprintln!(
        "\n   (domain {DOMAIN_KM} km, sea {SEA}, production seed; stage DELIVERED field, post-u16 \
         mask, bound ON in production)"
    );
}
