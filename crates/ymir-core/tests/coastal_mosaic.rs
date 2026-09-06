//! ADR rule 9 support — assemble the coastline panels into ONE comparison mosaic.
//!
//! Judging a visual defect from ten separate files means holding ten images in memory and
//! comparing from recollection, which is how three wrong attributions survived. Side by side in a
//! fixed order, with the REFERENCE first, the differences are direct.
//!
//! Run: cargo test -p ymir-core --release --test coastal_mosaic -- --ignored --nocapture

use std::path::Path;

const TILE: u32 = 640;
const LABEL_H: u32 = 44;
const GAP: u32 = 6;

/// 5×7 glyphs for the digits — `image` has no text rendering and a numbered tile plus a legend
/// is enough to make a mosaic readable.
const DIGITS: [[u8; 7]; 10] = [
    [0b111, 0b101, 0b101, 0b101, 0b101, 0b101, 0b111], // 0
    [0b010, 0b110, 0b010, 0b010, 0b010, 0b010, 0b111], // 1
    [0b111, 0b001, 0b001, 0b111, 0b100, 0b100, 0b111], // 2
    [0b111, 0b001, 0b001, 0b111, 0b001, 0b001, 0b111], // 3
    [0b101, 0b101, 0b101, 0b111, 0b001, 0b001, 0b001], // 4
    [0b111, 0b100, 0b100, 0b111, 0b001, 0b001, 0b111], // 5
    [0b111, 0b100, 0b100, 0b111, 0b101, 0b101, 0b111], // 6
    [0b111, 0b001, 0b001, 0b010, 0b010, 0b010, 0b010], // 7
    [0b111, 0b101, 0b101, 0b111, 0b101, 0b101, 0b111], // 8
    [0b111, 0b101, 0b101, 0b111, 0b001, 0b001, 0b111], // 9
];

fn draw_digit(img: &mut image::GrayImage, d: usize, ox: u32, oy: u32, scale: u32, v: u8) {
    for (r, row) in DIGITS[d].iter().enumerate() {
        for c in 0..3u32 {
            if row & (1 << (2 - c)) != 0 {
                for sy in 0..scale {
                    for sx in 0..scale {
                        let (x, y) = (ox + c * scale + sx, oy + r as u32 * scale + sy);
                        if x < img.width() && y < img.height() {
                            img.put_pixel(x, y, image::Luma([v]));
                        }
                    }
                }
            }
        }
    }
}

fn mosaic(files: &[(&str, &str)], cols: u32, out: &Path, legend_path: &Path) {
    let dir = Path::new("../../exports/coastal_fringes");
    let rows = (files.len() as u32).div_ceil(cols);
    let cell_w = TILE + GAP;
    let cell_h = TILE + LABEL_H + GAP;
    let mut canvas = image::GrayImage::from_pixel(
        cols * cell_w + GAP,
        rows * cell_h + GAP,
        image::Luma([28]), // dark chrome so the black sea stays distinguishable from the gaps
    );
    let mut legend = String::from("MOSAIC LEGEND — tiles left→right, top→bottom\n\n");
    for (i, (file, label)) in files.iter().enumerate() {
        let (cx, cy) = (i as u32 % cols, i as u32 / cols);
        let (ox, oy) = (GAP + cx * cell_w, GAP + cy * cell_h);
        legend.push_str(&format!("{:>2}  {label}\n", i + 1));
        let path = dir.join(file);
        let img = match image::open(&path) {
            Ok(im) => im.to_luma8(),
            Err(e) => {
                eprintln!("  MISSING {}: {e}", path.display());
                continue;
            }
        };
        // tile number, drawn in the label strip
        let n = i + 1;
        let scale = 5;
        let mut dx = ox + 8;
        if n >= 10 {
            draw_digit(&mut canvas, n / 10, dx, oy + 6, scale, 255);
            dx += 4 * scale;
        }
        draw_digit(&mut canvas, n % 10, dx, oy + 6, scale, 255);
        for y in 0..img.height().min(TILE) {
            for x in 0..img.width().min(TILE) {
                canvas.put_pixel(ox + x, oy + LABEL_H + y, *img.get_pixel(x, y));
            }
        }
    }
    canvas.save(out).expect("save mosaic");
    std::fs::write(legend_path, &legend).expect("legend");
    eprintln!("  {} ({}×{})", out.display(), canvas.width(), canvas.height());
    eprint!("{legend}");
}

#[test]
#[ignore]
fn build_coastal_mosaics() {
    let dir = Path::new("../../exports/coastal_fringes");
    eprintln!("\n=====  COASTLINE MOSAICS  =====\nsea black · land grey · coastline white");

    // REFERENCE FIRST in both: it is the target, and it is what made the target legible.
    eprintln!("\nLEVERS (8192², fixed crop):");
    mosaic(
        &[
            (
                "lever_REFERENCE_coarse_8192.png",
                "REFERENCE coarse — TARGET: 1.2 spurs/100km, p90 23.3 km",
            ),
            ("lever_SHIPPED_8192.png", "SHIPPED — 29.8 spurs/100km, p90 4.6 km"),
            ("lever_MFD_off__D8__8192.png", "MFD off (D8) — 27.3, p90 5.4"),
            ("lever_MFD_p___1_1__dispersed__8192.png", "MFD p=1.1 dispersed — 30.5, p90 4.6"),
            ("lever_MFD_p___4_0__concentrated__8192.png", "MFD p=4.0 concentrated — 29.5, p90 4.7"),
            ("lever_A_c___0_01_km2__x0_1__8192.png", "A_c x0.1 — 32.2, p90 4.7"),
            (
                "lever_A_c___1_0_km2__x10__8192.png",
                "A_c x10 — 9.3, p90 5.3 (but erosion largely off)",
            ),
            (
                "lever_A_c___10_0_km2__x100__8192.png",
                "A_c x100 — 2.6, p90 16.5 (EROSION OFF: hypsometry = reference)",
            ),
            ("lever_diffusion_x10__0_8__8192.png", "diffusion x10 — 38.6, WORSE at 8192"),
            ("lever_diffusion_x100__8_0__8192.png", "diffusion x100 — destroys the coast"),
        ],
        5,
        &dir.join("mosaic_levers.png"),
        &dir.join("mosaic_levers_legend.txt"),
    );

    eprintln!("\nWARP SWEEP (8192², same crop):");
    mosaic(
        &[
            ("reference_coarse_8192.png", "REFERENCE coarse — TARGET"),
            ("warp_1.50_8192.png", "warp 1.50 SHIPPED — 29.8 spurs/100km"),
            ("warp_1.00_8192.png", "warp 1.00"),
            ("warp_0.50_8192.png", "warp 0.50"),
            ("warp_0.25_8192.png", "warp 0.25"),
            ("warp_0.00_8192.png", "warp 0.00 — DENSER comb, more parallel: worst panel"),
        ],
        3,
        &dir.join("mosaic_warp.png"),
        &dir.join("mosaic_warp_legend.txt"),
    );
}
