//! ADR Finding 79 — the law re-read on the MASK, and the four panels the visual section needs.
//!
//! Run: cargo test -p ymir-core --release --test coast_mask_panels -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field};
use ymir_core::grid::GridF32;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds};
use ymir_core::terrain::contour::marching_squares;

const OX: usize = 4160;
const OY: usize = 1600;
const CROP: usize = 640;

fn mask_field(f: &GridF32) -> GridF32 {
    let mut m = GridF32::new(f.width, f.height, 0.0);
    for k in 0..f.data.len() {
        if f.data[k] > SEA {
            m.data[k] = 1.0;
        }
    }
    m
}

fn counts(polys: &[Vec<(f32, f32)>]) -> (usize, usize, f32, f32) {
    let a = coast_shape_thresholds(polys, CELL_KM, 1.0, NECK_KM);
    let b = coast_shape_thresholds(polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    (a.count, b.count, a.coast_km, b.local_axis_r)
}

/// Crop rasteriser: land grey, traced line white. `land` decides the grey, `polys` the line, so
/// a mask panel and an isoline panel differ ONLY in the line.
fn crop(f: &GridF32, polys: &[Vec<(f32, f32)>], path: &std::path::Path) -> GridF32 {
    let mut img = GridF32::new(CROP, CROP, 0.0);
    for y in 0..CROP {
        for x in 0..CROP {
            let (sx, sy) = (OX + x, OY + y);
            if sx < f.width && sy < f.height && f.data[sy * f.width + sx] > SEA {
                img.set(x, y, 0.55);
            }
        }
    }
    for pl in polys {
        for &(px, py) in pl.iter() {
            let (ix, iy) = (px.round() as isize - OX as isize, py.round() as isize - OY as isize);
            if ix >= 0 && iy >= 0 && (ix as usize) < CROP && (iy as usize) < CROP {
                img.set(ix as usize, iy as usize, 1.0);
            }
        }
    }
    img.save_png_u8(path).unwrap();
    img
}

#[test]
#[ignore]
fn coast_mask_panels() {
    let out = std::path::Path::new("../../exports/coastal_closure/panels");
    std::fs::create_dir_all(out).unwrap();
    eprintln!("\n==========  Finding 79 · the law on the MASK, and the panels  ==========");
    eprintln!(
        "   crop {OX},{OY} .. {},{} — the Finding 76 LAW_ON-densest window",
        OX + CROP,
        OY + CROP
    );

    // ── the law, on both definitions — the block-B2 line that was missing ─────
    eprintln!("\n── B2 (completed) · the channel-head law on both definitions ──");
    eprintln!(
        "   {:<28} {:>8} {:>8} {:>9} {:>7} | {:>8} {:>8} {:>9} {:>7}",
        "variant",
        "ISO km",
        "ISO cell",
        "ISO coast",
        "ISO R",
        "MSK km",
        "MSK cell",
        "MSK coast",
        "MSK R"
    );
    let mut fields: Vec<(&str, GridF32)> = Vec::new();
    for (nm, kn) in [
        ("PRE-INCISION", Knobs::no_incision()),
        ("SHIPPED (shelf 20)", Knobs { shelf_min_depth_m: Some(20.0), ..Knobs::shipped() }),
        ("PRODUCTION (shelf 1)", Knobs { shelf_min_depth_m: Some(1.0), ..Knobs::shipped() }),
        (
            "LAW ON (shelf 20)",
            Knobs { shelf_min_depth_m: Some(20.0), a_c_law: true, ..Knobs::shipped() },
        ),
        (
            "LAW ON (shelf 1)",
            Knobs { shelf_min_depth_m: Some(1.0), a_c_law: true, ..Knobs::shipped() },
        ),
    ] {
        let f = build_field(kn);
        let (ik, ic, icoast, ir) = counts(&marching_squares(&f, SEA));
        let (mk, mc, mcoast, mr) = counts(&marching_squares(&mask_field(&f), 0.5));
        eprintln!(
            "   {nm:<28} {ik:>8} {ic:>8} {icoast:>9.0} {ir:>7.3} | {mk:>8} {mc:>8} {mcoast:>9.0} {mr:>7.3}"
        );
        fields.push((nm, f));
    }

    // ── the panels: pre-incision / shelf 20 / shelf 1, each traced BOTH ways ──
    eprintln!("\n── the panels ──");
    let mut row_iso: Vec<GridF32> = Vec::new();
    let mut row_msk: Vec<GridF32> = Vec::new();
    for (nm, f) in fields.iter().take(3) {
        let tag = nm.split_whitespace().next().unwrap();
        let iso =
            crop(f, &marching_squares(f, SEA), &out.join(format!("F79_{tag}_ISO_{CROP}.png")));
        let msk = crop(
            f,
            &marching_squares(&mask_field(f), 0.5),
            &out.join(format!("F79_{tag}_MASK_{CROP}.png")),
        );
        row_iso.push(iso);
        row_msk.push(msk);
        eprintln!("   {nm:<28} -> F79_{tag}_ISO_{CROP}.png and F79_{tag}_MASK_{CROP}.png");
    }
    // one mosaic, two rows: isoline on top, mask below, same three columns
    let gap = 8usize;
    let (mw, mh) = (CROP * 3 + gap * 2, CROP * 2 + gap);
    let mut m = GridF32::new(mw, mh, 0.25);
    for (i, g) in row_iso.iter().enumerate() {
        for y in 0..CROP {
            for x in 0..CROP {
                m.set(i * (CROP + gap) + x, y, g.data[y * g.width + x]);
            }
        }
    }
    for (i, g) in row_msk.iter().enumerate() {
        for y in 0..CROP {
            for x in 0..CROP {
                m.set(i * (CROP + gap) + x, CROP + gap + y, g.data[y * g.width + x]);
            }
        }
    }
    let mp = out.join(format!("F79_mosaic_ISO_over_MASK_{CROP}.png"));
    m.save_png_u8(&mp).unwrap();
    eprintln!(
        "\n   MOSAIC 3x2 (columns: PRE-INCISION - SHELF 20 - SHELF 1; top row ISOLINE, bottom \
         row MASK) -> {}",
        mp.display()
    );
    assert!(std::fs::metadata(&mp).map(|f| f.len()).unwrap_or(0) > 4096);
}
