//! ADR Finding 80 — the visual section: the consumer masks WITH SHADER HALOS.
//!
//! First panel in this campaign that shows what a shader does with the contour. Halos are drawn
//! on the CONSUMER grids, as the author specified: sand on the terrain 2048 mask, foam on the
//! ocean 1024 mask, superposed. The crop is the Finding 76 LAW_ON-densest window mapped down to
//! 2048 (4160/4 .. 4800/4), upscaled x4 for legibility — declared, so nobody reads the
//! upscaling as detail.
//!
//! Run: cargo test -p ymir-core --release --test base_level_panels -- --ignored --nocapture

mod common;

use common::{Knobs, SEA, build_field};
use ymir_core::export::height::metric_height_u16;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;

#[allow(dead_code)] // ADR Finding 83: shipped, so no longer passed explicitly
const EPS: f32 = 0.5;
/// 8192 crop 4160,1600 .. 4800,2240 mapped to the 2048 grid.
const OX2: usize = 1040;
const OY2: usize = 400;
const CROP2: usize = 160;
const ZOOM: usize = 4;

fn land_u16(f: &GridF32, ss: &SteinSteinParams) -> Vec<bool> {
    let hl = metric_height_u16(f, ss);
    hl.codes
        .iter()
        .map(|&c| hl.min_m + (c as f32 / 65535.0) * (hl.max_m - hl.min_m) > 0.0)
        .collect()
}

fn majority(land: &[bool], w: usize, f: usize) -> (Vec<bool>, usize) {
    let nw = w / f;
    let mut out = vec![false; nw * nw];
    let half = f * f / 2;
    for y in 0..nw {
        for x in 0..nw {
            let mut c = 0usize;
            for dy in 0..f {
                for dx in 0..f {
                    if land[(y * f + dy) * w + x * f + dx] {
                        c += 1;
                    }
                }
            }
            out[y * nw + x] = c > half;
        }
    }
    (out, nw)
}

/// Chebyshev distance to the land/sea boundary, in pixels of this grid.
fn dist_to_edge(land: &[bool], w: usize) -> Vec<u16> {
    let n = w * w;
    let mut d = vec![u16::MAX; n];
    let mut q: Vec<u32> = Vec::new();
    for k in 0..n {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let mut edge = false;
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx >= 0
                    && ny >= 0
                    && (nx as usize) < w
                    && (ny as usize) < w
                    && land[ny as usize * w + nx as usize] != land[k]
                {
                    edge = true;
                }
            }
        }
        if edge {
            d[k] = 0;
            q.push(k as u32);
        }
    }
    let mut i = 0usize;
    while i < q.len() {
        let k = q[i] as usize;
        i += 1;
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= w {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if d[nk] == u16::MAX {
                    d[nk] = d[k] + 1;
                    q.push(nk as u32);
                }
            }
        }
    }
    d
}

fn blit(dst: &mut GridF32, src: &GridF32, ox: usize, oy: usize) {
    for y in 0..src.height {
        for x in 0..src.width {
            dst.set(ox + x, oy + y, src.data[y * src.width + x]);
        }
    }
}

#[test]
#[ignore]
fn base_level_panels() {
    let ss = SteinSteinParams::default();
    let out = std::path::Path::new("../../exports/coastal_closure/panels");
    std::fs::create_dir_all(out).unwrap();
    eprintln!("\n==========  Finding 80 · the consumer masks with shader halos  ==========");
    eprintln!(
        "   crop on the 2048 grid: {OX2},{OY2} .. {},{} ({CROP2} px = 31.2 km), upscaled x{ZOOM} \
         for legibility",
        OX2 + CROP2,
        OY2 + CROP2
    );

    let side = CROP2 * ZOOM;
    let gap = 8usize;
    let mut mosaic = GridF32::new(side * 3 + gap * 2, side * 3 + gap * 2, 0.25);

    for (col, (nm, kn)) in [
        ("PRE-INCISION", Knobs::no_incision()),
        ("DELIVERED (pre-83)", Knobs::pre83()),
        ("BASE LEVEL (production)", Knobs::shipped()),
    ]
    .into_iter()
    .enumerate()
    {
        let f = build_field(kn);
        let l8 = land_u16(&f, &ss);
        let (l2, w2) = majority(&l8, f.width, 4);
        let (l1, w1) = majority(&l8, f.width, 8);
        let d2 = dist_to_edge(&l2, w2);
        let d1 = dist_to_edge(&l1, w1);

        // row 0 — the bare terrain mask
        let mut bare = GridF32::new(side, side, 0.0);
        // row 1 — the same, with halos: sand (land, d < 1 px on 2048) and foam (sea, d < 1 px on 1024)
        let mut halo = GridF32::new(side, side, 0.0);
        // row 2 — the XOR of the two shorelines
        let mut xor = GridF32::new(side, side, 0.0);
        let (mut n_sand, mut n_foam, mut n_xor) = (0usize, 0usize, 0usize);
        for y in 0..CROP2 {
            for x in 0..CROP2 {
                let (sx, sy) = (OX2 + x, OY2 + y);
                let k2 = sy * w2 + sx;
                let k1 = (sy / 2) * w1 + sx / 2;
                let land = l2[k2];
                let b = if land { 0.55f32 } else { 0.0 };
                let h = if land && d2[k2] == 0 {
                    n_sand += 1;
                    1.0 // sand
                } else if !land && d1[k1] == 0 {
                    n_foam += 1;
                    0.85 // foam
                } else {
                    b
                };
                let disagree = l2[k2] != l1[k1];
                if disagree {
                    n_xor += 1;
                }
                let xv = if disagree { 1.0 } else { b };
                for zy in 0..ZOOM {
                    for zx in 0..ZOOM {
                        let (px, py) = (x * ZOOM + zx, y * ZOOM + zy);
                        bare.set(px, py, b);
                        halo.set(px, py, h);
                        xor.set(px, py, xv);
                    }
                }
            }
        }
        blit(&mut mosaic, &bare, col * (side + gap), 0);
        blit(&mut mosaic, &halo, col * (side + gap), side + gap);
        blit(&mut mosaic, &xor, col * (side + gap), 2 * (side + gap));
        eprintln!(
            "   {nm:<14} in-crop: sand {n_sand} px · foam {n_foam} px · shoreline disagreement \
             (XOR) {n_xor} px of {}",
            CROP2 * CROP2
        );
    }
    let mp = out.join("F80_mosaic_MASK_HALO_XOR.png");
    mosaic.save_png_u8(&mp).unwrap();
    eprintln!(
        "\n   MOSAIC 3x3 -> {}\n   columns: PRE-INCISION - DELIVERED - BASE LEVEL\n   rows: bare \
         2048 mask / with halos (white = sand, light grey = foam) / XOR of the 2048 and 1024 \
         shorelines",
        mp.display()
    );
    assert!(std::fs::metadata(&mp).map(|f| f.len()).unwrap_or(0) > 4096);
}
