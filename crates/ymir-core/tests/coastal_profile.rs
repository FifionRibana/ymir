//! ADR Finding 78 — the visual section: the shoreline seen in CROSS-SECTION, and three panels.
//!
//! The cross-section is the new artefact: a median height profile over 200 transects, 20 cells
//! of sea to 20 cells of land, walked perpendicular to the coast. It is what a player walking
//! into the water would feel, and it is where a 20 m step is undeniable.
//!
//! Run: cargo test -p ymir-core --release --test coastal_profile -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, sorted};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::contour::marching_squares;

/// The crop chosen on LAW_ON's densest spur window in Finding 76 — the one where the author
/// saw the saw-tooth, and therefore the one where a remedy has to be judged. Declared, fixed.
const OX: usize = 4160;
const OY: usize = 1600;
const CROP: usize = 640;
const TRANSECTS: usize = 200;
const REACH: i32 = 20;

fn crop_png(f: &GridF32, polys: &[Vec<(f32, f32)>], out: &std::path::Path, name: &str) -> GridF32 {
    let mut img = GridF32::new(CROP, CROP, 0.0);
    for y in 0..CROP {
        for x in 0..CROP {
            let (sx, sy) = (OX + x, OY + y);
            if sx < f.width && sy < f.height && f.data[sy * f.width + sx] > SEA {
                img.set(x, y, 0.65);
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
    img.save_png_u8(&out.join(format!("{name}.png"))).unwrap();
    img
}

/// Median height profile across the shoreline. For each of `TRANSECTS` sampled coastline points,
/// the local outward normal is estimated from the polyline tangent, and the field is sampled at
/// `-REACH..=REACH` cells along it. Negative offsets are seaward.
fn profile(f: &GridF32, polys: &[Vec<(f32, f32)>], n2m: f32) -> Vec<f32> {
    let (w, h) = (f.width, f.height);
    let mut cols: Vec<Vec<f32>> = vec![Vec::new(); (2 * REACH + 1) as usize];
    // longest polylines first, sampled evenly, so the transects are not all on one islet
    let mut idx: Vec<usize> = (0..polys.len()).filter(|&i| polys[i].len() >= 32).collect();
    idx.sort_by_key(|&i| std::cmp::Reverse(polys[i].len()));
    let mut taken = 0usize;
    'outer: for &i in &idx {
        let pl = &polys[i];
        let step = (pl.len() / 16).max(8);
        let mut j = 4usize;
        while j + 4 < pl.len() {
            let (ax, ay) = pl[j - 4];
            let (bx, by) = pl[j + 4];
            let (tx, ty) = (bx - ax, by - ay);
            let len = (tx * tx + ty * ty).sqrt().max(1e-6);
            // outward normal: rotate the tangent, then orient it so it points at LOWER ground
            let (mut nx, mut ny) = (-ty / len, tx / len);
            let (px, py) = pl[j];
            let probe = |dx: f32, dy: f32| -> f32 {
                let (sx, sy) = ((px + dx).round() as i64, (py + dy).round() as i64);
                if sx < 0 || sy < 0 || sx as usize >= w || sy as usize >= h {
                    return f32::NAN;
                }
                f.data[sy as usize * w + sx as usize]
            };
            let (a, b) = (probe(nx * 3.0, ny * 3.0), probe(-nx * 3.0, -ny * 3.0));
            if a.is_nan() || b.is_nan() {
                j += step;
                continue;
            }
            if a > b {
                nx = -nx;
                ny = -ny;
            }
            let mut col = Vec::with_capacity((2 * REACH + 1) as usize);
            let mut ok = true;
            for s in -REACH..=REACH {
                let v = probe(nx * s as f32, ny * s as f32);
                if v.is_nan() {
                    ok = false;
                    break;
                }
                col.push((v - SEA) * n2m);
            }
            if ok {
                for (c, v) in col.into_iter().enumerate() {
                    cols[c].push(v);
                }
                taken += 1;
                if taken >= TRANSECTS {
                    break 'outer;
                }
            }
            j += step;
        }
    }
    assert!(taken >= TRANSECTS / 2, "rule 10: only {taken} transects, not a median");
    cols.into_iter()
        .map(|c| {
            let s = sorted(c);
            if s.is_empty() { f32::NAN } else { s[s.len() / 2] }
        })
        .collect()
}

#[test]
#[ignore]
fn coastal_profile() {
    let ss =
        ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let out = std::path::Path::new("../../exports/coastal_closure/panels");
    std::fs::create_dir_all(out).unwrap();
    eprintln!("\n==========  Finding 78 · the shoreline in cross-section  ==========");
    eprintln!(
        "   crop {OX},{OY} .. {},{} ({CROP}x{CROP} = {:.1} km) — LAW_ON's densest spur window \
         (Finding 76), where the saw-tooth was seen",
        OX + CROP,
        OY + CROP,
        CROP as f32 * CELL_KM
    );

    let variants: [(&str, Knobs); 3] = [
        ("PRE_INCISION", Knobs::no_incision()),
        ("SHIPPED", Knobs::shipped()),
        (
            "SHELF1_DIFF128",
            Knobs {
                shelf_min_depth_m: Some(1.0),
                diffusion: Some(1.28),
                diffusion_substeps: Some(8),
                ..Knobs::shipped()
            },
        ),
    ];
    let mut profiles: Vec<(&str, Vec<f32>)> = Vec::new();
    let mut crops: Vec<GridF32> = Vec::new();
    for (name, kn) in variants {
        let f = build_field(kn);
        let polys = marching_squares(&f, SEA);
        crops.push(crop_png(&f, &polys, out, &format!("F78_{name}_{CROP}")));
        profiles.push((name, profile(&f, &polys, n2m)));
        eprintln!("   built {name}");
    }
    // the shelf-1 field without the diffusion change, for the step comparison alone
    let f1 = build_field(Knobs { shelf_min_depth_m: Some(1.0), ..Knobs::shipped() });
    let p1 = marching_squares(&f1, SEA);
    profiles.push(("SHELF1_ONLY", profile(&f1, &p1, n2m)));

    eprintln!("\n── THE CROSS-SECTION · median over {TRANSECTS} transects, metres above sea ──");
    eprintln!("   offset in cells: negative = seaward, 0 = the traced coastline");
    eprint!("   {:<16}", "offset");
    for s in (-REACH..=REACH).step_by(2) {
        eprint!("{s:>8}");
    }
    eprintln!();
    for (name, p) in &profiles {
        eprint!("   {name:<16}");
        for (i, v) in p.iter().enumerate() {
            if (i as i32 - REACH) % 2 == 0 {
                eprint!("{v:>8.1}");
            }
        }
        eprintln!();
    }

    // the mosaic, REFERENCE first
    let gap = 8usize;
    let mut m = GridF32::new(CROP * 3 + gap * 2, CROP, 0.25);
    for (i, g) in crops.iter().enumerate() {
        for y in 0..CROP {
            for x in 0..CROP {
                m.set(i * (CROP + gap) + x, y, g.data[y * g.width + x]);
            }
        }
    }
    let mp = out.join(format!("F78_mosaic_PRE_SHIPPED_SHELF1DIFF_{CROP}.png"));
    m.save_png_u8(&mp).unwrap();
    eprintln!(
        "\n   MOSAIC (left to right: PRE-INCISION · SHIPPED · SHELF1+DIFF1.28) -> {}",
        mp.display()
    );
    assert!(std::fs::metadata(&mp).map(|f| f.len()).unwrap_or(0) > 4096);
}
