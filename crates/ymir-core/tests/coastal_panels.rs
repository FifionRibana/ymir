//! ADR Finding 75 — the PANELS, rendered FROM THE EXPORTS themselves.
//!
//! Rule 9 wants a panel, a reference for "absent", a declared crop, and absolute counts that
//! agree with the picture. Finding 75 delivered the counts and three `.ymir` containers; this
//! renders the pictures.
//!
//! **It reads the exported `height.u16` and `lake_mask.u32` back off disk and traces the
//! coastline from them**, rather than re-building the fields. That is the point: a panel derived
//! from the same bytes the author's tool opens cannot disagree with the export the way a
//! separately-rebuilt one can — and it costs three file reads instead of three 8192² builds.
//!
//! ⚠️ The `exports/coastal_fringes/law_*.png` panels of Finding 56 are STALE with respect to
//! VALIDATION_NOTE §1 as amended: they were rendered with `S_ref = 0.3319` (the 2048² value).
//! These are the 0.1128 target-grid panels.
//!
//! ⚠️ PATHS. A test's cwd is the CRATE root, not the workspace root, so the export path is
//! `../../exports/...`. The first run of `coastal_closure` wrote 3.3 GB into
//! `crates/ymir-core/exports/` and the delivered path was therefore wrong; both benches now
//! write and read the workspace `exports/` where every other export lives.
//!
//! Run: cargo test -p ymir-core --release --test coastal_panels -- --ignored --nocapture

use ymir_core::export::raw;
use ymir_core::grid::GridF32;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds, densest_window};
use ymir_core::terrain::contour::marching_squares;

const ROOT: &str = "../../exports/coastal_closure";
const SEED: &str = "seed10481999410520546993_8192";
const CROP: usize = 640;
const CELL_KM: f32 = 400.0 / 8192.0;
/// Sea in the DECODED metric field. The export's contract is `sea_level_m = 0.0`
/// (`hd.rs:1227`, "sea anchored to 0 m"), so the isoline is traced at 0 m, not at 0.5.
const SEA_M: f32 = 0.0;

/// Decode `height.u16` back to metres with the manifest's own range — the exact inverse of
/// `height::metric_height_u16`, so this is the field the consumer sees.
fn load_metric_height(dir: &std::path::Path) -> (GridF32, usize) {
    let manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).unwrap()).unwrap();
    let w = manifest["continent"]["grid"]["width"].as_u64().unwrap() as usize;
    let h = manifest["continent"]["grid"]["height"].as_u64().unwrap() as usize;
    let layer = manifest["layers"]
        .as_array()
        .and_then(|a| a.iter().find(|l| l["id"] == "height"))
        .expect("the export must declare a height layer");
    let (lo, hi) = (
        layer["min_m"].as_f64().expect("height layer must carry min_m") as f32,
        layer["max_m"].as_f64().expect("height layer must carry max_m") as f32,
    );
    let codes = raw::load_u16(&dir.join("height.u16"), w * h).unwrap();
    let mut g = GridF32::new(w, h, 0.0);
    for (k, &c) in codes.iter().enumerate() {
        g.data[k] = lo + (c as f32 / 65535.0) * (hi - lo);
    }
    eprintln!(
        "   {:<10} {w}×{h}, metric range {lo:.1} … {hi:.1} m",
        dir.file_name().unwrap().to_string_lossy()
    );
    (g, w)
}

/// ONE rasterisation, used by BOTH the per-variant panel and the mosaic, so the two cannot
/// drift. It reproduces `render_coast_crop`'s rule verbatim — land above the sea reads 0.65,
/// a coastline sample reads 1.0.
///
/// ⚠️ **NAMED DEBT.** `terrain::coast_metrics::render_coast_crop` was promoted to the library
/// (Finding 56, rule 9) precisely so one definition serves every sweep, and this duplicates it —
/// because that function WRITES a file and cannot RETURN the image, so a mosaic cannot be built
/// from it. The fix is to split it into `coast_crop_image(...) -> GridF32` plus a save; that is
/// a production change and this round's budget of two was spent.
fn crop_image(f: &GridF32, polys: &[Vec<(f32, f32)>], sea: f32, ox: usize, oy: usize) -> GridF32 {
    let mut img = GridF32::new(CROP, CROP, 0.0);
    for y in 0..CROP {
        for x in 0..CROP {
            let (sx, sy) = (ox + x, oy + y);
            if sx < f.width && sy < f.height && f.data[sy * f.width + sx] > sea {
                img.set(x, y, 0.65);
            }
        }
    }
    for pl in polys {
        for &(px, py) in pl.iter() {
            let (ix, iy) = (px.round() as isize - ox as isize, py.round() as isize - oy as isize);
            if ix >= 0 && iy >= 0 && (ix as usize) < CROP && (iy as usize) < CROP {
                img.set(ix as usize, iy as usize, 1.0);
            }
        }
    }
    img
}

#[test]
#[ignore]
fn coastal_panels() {
    let root = std::path::Path::new(ROOT);
    let out = root.join("panels");
    std::fs::create_dir_all(&out).unwrap();
    let variants = ["REFERENCE", "SHIPPED", "LAW_ON"];

    eprintln!("\n==========  Finding 75 · the panels, from the exported rasters  ==========");
    let mut fields = Vec::new();
    for v in variants {
        let dir = root.join(format!("{SEED}_{v}.ymir"));
        assert!(dir.join("height.u16").exists(), "missing export: {}", dir.display());
        let (f, w) = load_metric_height(&dir);
        let polys = marching_squares(&f, SEA_M);
        fields.push((v, f, polys, w));
    }

    // THE CROP, declared once and used on all three: the densest spur window of SHIPPED, which
    // is Finding 56's own convention. Choosing it on the law panel would flatter a remedy that
    // merely displaces the defect.
    let shipped = fields.iter().find(|e| e.0 == "SHIPPED").unwrap();
    let (ox, oy) = densest_window(&shipped.2, shipped.3, shipped.3, CROP);
    eprintln!(
        "\n   CROP (densest spur window of SHIPPED, used on ALL THREE): cells \
         **{ox},{oy} .. {},{}** — {CROP}×{CROP} = {:.1} × {:.1} km",
        ox + CROP,
        oy + CROP,
        CROP as f32 * CELL_KM,
        CROP as f32 * CELL_KM
    );

    let mut crops: Vec<GridF32> = Vec::new();
    eprintln!(
        "\n   {:<12} {:>14} {:>16} {:>12} {:>14} {:>16}",
        "panel", "≥1 km (whole)", "≥2 cells (whole)", "coast km", "≥1 km (CROP)", "≥2 cells (CROP)"
    );
    for (v, f, polys, _) in &fields {
        let whole_km = coast_shape_thresholds(polys, CELL_KM, 1.0, NECK_KM);
        let whole_cell = coast_shape_thresholds(polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
        // the crop's OWN counts, so the picture and a number are comparable at the same extent
        let clipped: Vec<Vec<(f32, f32)>> = polys
            .iter()
            .map(|pl| {
                pl.iter()
                    .filter(|(x, y)| {
                        *x >= ox as f32
                            && *x < (ox + CROP) as f32
                            && *y >= oy as f32
                            && *y < (oy + CROP) as f32
                    })
                    .copied()
                    .collect()
            })
            .filter(|pl: &Vec<(f32, f32)>| pl.len() >= 8)
            .collect();
        let crop_km = coast_shape_thresholds(&clipped, CELL_KM, 1.0, NECK_KM);
        let crop_cell = coast_shape_thresholds(&clipped, CELL_KM, 2.0 * CELL_KM, CELL_KM);
        let path = out.join(format!("panel_{v}_{CROP}.png"));
        let img = crop_image(f, polys, SEA_M, ox, oy);
        img.save_png_u8(&path).unwrap();
        crops.push(img);
        eprintln!(
            "   {v:<12} {:>14} {:>16} {:>12.0} {:>14} {:>16}   → {}",
            whole_km.count,
            whole_cell.count,
            whole_km.coast_km,
            crop_km.count,
            crop_cell.count,
            path.display()
        );
        // rule 6, mechanical: a panel that was not written is not a deliverable
        let bytes = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
        assert!(bytes > 1024, "panel {v} is {bytes} bytes on disk — not written");
    }

    // ── the two lake crops of block E, same treatment ─────────────────────────
    // 1000004: both halves fit in one small frame. 1000056: the main body plus the single
    // cell 11.44 km away that receives 380 m³/s — the question is whether it is visible at all.
    let dir = root.join(format!("{SEED}_SHIPPED.ymir"));
    let (f, w) = (&shipped.1, shipped.3);
    let lake = raw::load_u32(&dir.join("lake_mask.u32"), w * w).unwrap();
    for (name, cx, cy, crop) in
        [("lake_1000004", 4422usize, 2022usize, 64usize), ("lake_1000056", 4029, 5086, 1024)]
    {
        let (lox, loy) = (cx.saturating_sub(crop / 2), cy.saturating_sub(crop / 2));
        let mut img = GridF32::new(crop, crop, 0.0);
        let (mut n_lake, mut n_land) = (0usize, 0usize);
        for y in 0..crop {
            for x in 0..crop {
                let (sx, sy) = (lox + x, loy + y);
                if sx >= w || sy >= w {
                    continue;
                }
                let k = sy * w + sx;
                if lake[k] != 0 {
                    img.set(x, y, 1.0); // water body, whatever its id
                    n_lake += 1;
                } else if f.data[k] > SEA_M {
                    img.set(x, y, 0.55); // land
                    n_land += 1;
                }
            }
        }
        // A single cell is invisible at this scale, and that IS the question, so it is MARKED
        // rather than left to be hunted: a hollow crosshair around cell 4128,5223 — the 0.002 km²
        // half of id 1000056 that receives 380.5 m³/s (Finding 75 block E).
        if name == "lake_1000056" {
            let (mx, my) = (4128i64 - lox as i64, 5223i64 - loy as i64);
            for d in 4..12i64 {
                for (px, py) in [(mx + d, my), (mx - d, my), (mx, my + d), (mx, my - d)] {
                    if px >= 0 && py >= 0 && (px as usize) < crop && (py as usize) < crop {
                        img.set(px as usize, py as usize, 0.0);
                    }
                }
            }
            eprintln!(
                "     ↑ the isolated 1-cell half is at PANEL pixel {mx},{my}, marked with a                  crosshair (the cell itself is 1 px)"
            );
        }
        let path = out.join(format!("{name}_{crop}.png"));
        img.save_png_u8(&path).unwrap();
        eprintln!(
            "   {name:<14} crop {crop}×{crop} at cells {lox},{loy} ({:.2} × {:.2} km) | \
             {n_lake} lake cells, {n_land} land cells → {}",
            crop as f32 * CELL_KM,
            crop as f32 * CELL_KM,
            path.display()
        );
        assert!(n_lake > 0, "rule 10: {name}'s crop contains no lake cell at all");
        assert!(std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0) > 256);
    }

    // ── the MOSAIC (rule 9: judging three files means judging from recollection) ──
    // REFERENCE first, left to right, so "absent" anchors the eye before the defect.
    {
        let gap = 8usize;
        let mw = CROP * 3 + gap * 2;
        let mut m = GridF32::new(mw, CROP, 0.25);
        for (i, g) in crops.iter().enumerate() {
            let x0 = i * (CROP + gap);
            for y in 0..CROP {
                for x in 0..CROP {
                    m.set(x0 + x, y, g.data[y * g.width + x]);
                }
            }
        }
        let mp = out.join(format!("mosaic_REFERENCE_SHIPPED_LAW_{CROP}.png"));
        m.save_png_u8(&mp).unwrap();
        eprintln!(
            "
   MOSAIC (left→right: REFERENCE · SHIPPED · LAW_ON, same crop) → {}",
            mp.display()
        );
        assert!(std::fs::metadata(&mp).map(|f| f.len()).unwrap_or(0) > 4096);
    }
}
