//! ADR 0001 Finding 84 block B1 — **the candidate "organic" quantities, calibrated on shapes of
//! known character BEFORE any of them is called a criterion.**
//!
//! Finding 83 wrote a criterion in advance ("D in [1.10 ; 1.30]") and the fringe passed it. The
//! lesson is in the round's own words: *"pas de seuils d'« organique » posés avant la mesure sur
//! les côtes réelles"*. So this bench does one thing only — it asks, of each candidate quantity,
//! whether it can tell a PERIODIC coast (a smooth lobe with regular teeth at the fur's measured
//! geometry) from an ISOTROPIC one (a radial fBm boundary rough at every scale, D = 1.25).
//! A quantity that cannot is struck out. **No threshold is written here and no criterion is
//! declared**; the thresholds have to be read on real coasts, and
//!
//! ⚠️ **there is no real coastline raster in this repository.** Every `coastline.geojson` under
//! `exports/` is Ymir's own output; `docs/refs/` holds papers. So block B2 does not run, by the
//! round's own instruction, and the two Ymir coasts below are printed as CONTEXT, not as a test.
//!
//! Run: cargo test -p ymir-core --release --test f84_organic -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, build_field, circle_polygon, curvature_stats, excursion_sizes, fbm_polygon,
    fill_polygon, koch_polygon, land_u16, normal_entropy, random_koch_polygon, richardson,
    richardson_fit, rotated_square_polygon, spiky_circle_mask, to_mask,
};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::terrain::coast_metrics::coast_shape_thresholds;
use ymir_core::terrain::contour::marching_squares;

/// Arc-length window the normal histogram is read over. 20 km is the cut `spectrum()` already
/// uses for "a continuous coast segment long enough to read", so the two instruments share a
/// population convention.
const WIN_KM: f32 = 20.0;
/// "Within three cells" — the round's own wording for a curvature inversion that is a tooth
/// rather than a headland.
const NEAR_CELLS: f32 = 3.0;
/// Curvature stencil, in contour samples either side. Reported at two values because the first
/// calibration showed the instrument's answer IS the stencil: at 2 samples a rasterised circle
/// reads 55.8 % of its (zero) inversions inside three cells.
const STENCILS: [usize; 2] = [2, 8];

struct Row {
    tag: &'static str,
    d: f64,
    r: f32,
    entropy: f64,
    windows: usize,
    r50: [f32; 2],
    inv_near: [f64; 2],
    inv_n: [usize; 2],
    e10: f32,
    e50: f32,
    e90: f32,
    e_cv: f64,
    e_n: usize,
}

fn measure(tag: &'static str, land: &[bool], w: usize, cell_km: f32) -> Row {
    let polys = marching_squares(&to_mask(land, w), 0.5);
    let keep: Vec<&Vec<(f32, f32)>> = polys.iter().collect();
    let d = richardson_fit(&richardson(&keep, cell_km, 10.0)).0;
    let shape = coast_shape_thresholds(&polys, cell_km, 2.0 * cell_km, cell_km);
    let (entropy, windows) = normal_entropy(&polys, cell_km, WIN_KM);
    let a = curvature_stats(&polys, NEAR_CELLS, STENCILS[0]);
    let b = curvature_stats(&polys, NEAR_CELLS, STENCILS[1]);
    let (e10, e50, e90, e_cv, e_n) = excursion_sizes(&polys, cell_km);
    Row {
        tag,
        d,
        r: shape.local_axis_r,
        entropy,
        windows,
        r50: [a.1, b.1],
        inv_near: [a.3, b.3],
        inv_n: [a.4, b.4],
        e10,
        e50,
        e90,
        e_cv,
        e_n,
    }
}

fn print_rows(rows: &[Row]) {
    eprintln!(
        "   {:<26} {:>7} {:>7} {:>9} {:>6} {:>9} {:>10} {:>9} {:>10}",
        "shape", "D", "R", "H(norm)", "wins", "Rc s=2", "inv<3c s=2", "Rc s=8", "inv<3c s=8"
    );
    for r in rows {
        eprintln!(
            "   {:<26} {:>7.3} {:>7.3} {:>9.4} {:>6} {:>9.2} {:>9.1} % {:>9.2} {:>8.1} %",
            r.tag, r.d, r.r, r.entropy, r.windows, r.r50[0], r.inv_near[0], r.r50[1], r.inv_near[1]
        );
    }
    eprintln!(
        "\n   {:<26} {:>9} {:>9} {:>9} {:>9} {:>9}",
        "shape", "exc p10", "p50", "p90", "CV", "count"
    );
    for r in rows {
        eprintln!(
            "   {:<26} {:>9.2} {:>9.2} {:>9.2} {:>9.3} {:>9}",
            r.tag, r.e10, r.e50, r.e90, r.e_cv, r.e_n
        );
    }
}

#[test]
#[ignore]
fn f84_organic() {
    let ss = SteinSteinParams::default();
    let w = 8192usize;
    eprintln!("\n==========  Finding 84 · B1 — calibrating the candidates  ==========");
    eprintln!(
        "\n   Quantities: D (Richardson, divider) · R (local axial concentration of excursions, \
         the Finding 52 instrument on the mask) · H(norm) = normalised entropy of the contour \
         NORMAL's axial angle over {WIN_KM} km windows, 18 bins · R_c = curvature radius in \
         cells, and the share of curvature SIGN INVERSIONS closer than {NEAR_CELLS} cells of arc \
         · excursion sizes in cells from the same `coast_spurs` walk the counts use."
    );

    // ── the three shapes of KNOWN character ───────────────────────────────────
    let circ = fill_polygon(&circle_polygon(w, 3000.0), w);
    let sq = fill_polygon(&rotated_square_polygon(w, 5000.0, 30.0), w);
    let kch = fill_polygon(&koch_polygon(w, 5400.0, 7), w);
    // ── the two negative controls the round demanded, SECOND attempt ──────────
    // PERIODIC: Finding 76's measured fur geometry, built on the GRID — 3-cell spikes one cell
    // wide every 14 cells. The polygon version (2-cell teeth, 7-cell necks) found ZERO
    // excursions and is kept below only as the measured negative it is.
    let per = spiky_circle_mask(w, 3000.0, 3.0, 14.0);
    // ISOTROPIC: a RANDOMISED Koch (apex side thrown at random) — same D = 1.2619, no
    // periodicity, no preferred axis, and rough down to 2.5 cells.
    let iso = fill_polygon(&random_koch_polygon(w, 5400.0, 7, 0xC0FF_EE01), w);
    // the two that FAILED as controls, kept so the failure is on the record
    let per_bad = fill_polygon(&common::periodic_polygon(w, 3000.0, 2.0, 14.0), w);
    let iso_bad = fill_polygon(&fbm_polygon(w, 3000.0, 1.25, 0.05, 512), w);

    eprintln!("\n── B1 · the calibration shapes ──");
    let rows = vec![
        measure("CIRCLE (D 1.000)", &circ, w, CELL_KM),
        measure("SQUARE rot 30 (D 1.000)", &sq, w, CELL_KM),
        measure("KOCH g7 (D 1.2619)", &kch, w, CELL_KM),
        measure("**PERIODIC** spikes 3c/14c", &per, w, CELL_KM),
        measure("**ISOTROPIC** random Koch", &iso, w, CELL_KM),
        measure("(void) periodic teeth 2c", &per_bad, w, CELL_KM),
        measure("(void) fBm lobe D=1.25", &iso_bad, w, CELL_KM),
    ];
    print_rows(&rows);

    // ── the verdict, per quantity, on the one question B1 asks ────────────────
    let p = &rows[3];
    let i = &rows[4];
    eprintln!("\n── B1 · does the quantity separate PERIODIC from ISOTROPIC? ──");
    let sep = |name: &str, a: f64, b: f64, unit: &str| {
        let ratio = if b.abs() > 1e-12 { a / b } else { f64::NAN };
        eprintln!(
            "   {name:<34} periodic {a:>9.4}{unit}  isotropic {b:>9.4}{unit}  ratio {ratio:>7.2}"
        );
    };
    sep("D (Richardson)", p.d, i.d, "");
    sep("R (local axial concentration)", p.r as f64, i.r as f64, "");
    sep("H(normal axis)", p.entropy, i.entropy, "");
    sep("curvature radius p50, s=8", p.r50[1] as f64, i.r50[1] as f64, "");
    sep("inversions < 3 cells, s=8", p.inv_near[1], i.inv_near[1], " %");
    sep("excursion size p50 (cells)", p.e50 as f64, i.e50 as f64, "");
    sep("excursion size CV", p.e_cv, i.e_cv, "");

    // ── CONTEXT ONLY: the two Ymir coasts, on the same instruments ────────────
    eprintln!(
        "\n── CONTEXT, not a test (no real coast on disk ⇒ no thresholds): the two Ymir coasts ──"
    );
    let del = build_field(Knobs::pre83());
    let bnd = build_field(Knobs::shipped());
    let rows2 = vec![
        measure("YMIR delivered (the fur)", &land_u16(&del, &ss), w, CELL_KM),
        measure("YMIR production (bounded)", &land_u16(&bnd, &ss), w, CELL_KM),
    ];
    print_rows(&rows2);
    eprintln!(
        "\n   ⚠️ **No criterion is declared from these rows.** Finding 83 wrote a window in \
         advance and the fringe passed it; the thresholds for an organic coast have to be read \
         off REAL coastlines, and there is none in this repository. B2 does not run."
    );
}
