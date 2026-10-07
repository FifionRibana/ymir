//! ADR Finding 152-B3 — the ROCK GRID: one u8 class per HD cell, fixed for a world (the author: « la roche n'est pas
//! une variable du monde »). The nine classes of `spec_geologie_v1.md` §1.
//!
//! - Coarse sources (craton, belt, rift) never draw a contact on the 64² grid: the coarse mask is sampled bilinearly
//!   onto a 1024² grid with the terrain's own mapping, smoothed by a Gaussian of σ = 8 cells there (≈ 3 km, PROXY),
//!   sampled bilinearly at the cell and cut at 0.5.
//! - Volcanic rock follows the eroded cone's relief, not a disc.
//! - Loose deposits and evaporites come from the HD fields (thresholds PROXY).
//!
//! ADR Finding 154-S — two layers: the SUBSTRATUM ([`build_substratum`]: craton, belt, basaltic or arc volcanic, rift
//! fill, basement), from the tectonics and the edifices read on a relief, available from the construction stage on;
//! and the SURFACE ([`apply_surface`]: water, evaporites, loose deposits), from the end of the chain. The rock grid is
//! the surface over the substratum ([`build_rocks`]).

use rayon::prelude::*;

use crate::grid::GridF32;
use crate::tectonics_c1::closures::lithology::upscale_k_to_hd;
use crate::tectonics_c1::closures::volcanism::{Edifice, VolcanoSetting};
use crate::terrain::convergence::gaussian_smooth;

/// One rock class.
#[derive(Clone, Copy, Debug)]
pub struct RockClass {
    pub id: u8,
    /// The key the rules file uses.
    pub key: &'static str,
    pub name_fr: &'static str,
    pub color: [u8; 3],
    /// For the falls at lake sills (the gorge, paused): 2 hard, 1 medium, 0 soft. The contrast is ANCHORED (Stock &
    /// Montgomery 1999: crystalline vs mudstone / volcaniclastic), the values PROXY.
    pub hardness: u8,
}

pub const NONE: u8 = 0;
pub const CRATON: u8 = 1;
pub const BELT: u8 = 2;
pub const BASALTIC_VOLCANIC: u8 = 3;
pub const ARC_VOLCANIC: u8 = 4;
pub const RIFT_FILL: u8 = 5;
pub const LOOSE_DEPOSITS: u8 = 6;
pub const EVAPORITES: u8 = 7;
pub const BASEMENT: u8 = 8;

/// The nine classes, indexed by id.
pub const ROCK_CLASSES: [RockClass; 9] = [
    RockClass { id: NONE, key: "none", name_fr: "Eau (mer, lac)", color: [40, 60, 90], hardness: 0 },
    RockClass { id: CRATON, key: "craton", name_fr: "Craton (granite, gneiss)", color: [0xC9, 0xA5, 0x5A], hardness: 2 },
    RockClass { id: BELT, key: "belt", name_fr: "Ceinture (métamorphique)", color: [0xB5, 0x48, 0x7F], hardness: 2 },
    RockClass { id: BASALTIC_VOLCANIC, key: "basaltic_volcanic", name_fr: "Volcanique basaltique", color: [0x4A, 0x4A, 0x52], hardness: 2 },
    RockClass { id: ARC_VOLCANIC, key: "arc_volcanic", name_fr: "Volcanique d'arc (andésite, tufs)", color: [0x9B, 0x59, 0xB6], hardness: 1 },
    RockClass { id: RIFT_FILL, key: "rift_fill", name_fr: "Remplissage de rift", color: [0x2E, 0xB8, 0xA8], hardness: 0 },
    RockClass { id: LOOSE_DEPOSITS, key: "loose_deposits", name_fr: "Dépôts meubles récents", color: [0xE2, 0xC8, 0x96], hardness: 0 },
    RockClass { id: EVAPORITES, key: "evaporites", name_fr: "Évaporites", color: [0xF2, 0xF0, 0xE6], hardness: 0 },
    RockClass { id: BASEMENT, key: "basement", name_fr: "Socle continental", color: [0x8C, 0x8C, 0x8C], hardness: 2 },
];

/// The id of a rock key (`None` if unknown).
#[must_use]
pub fn rock_id(key: &str) -> Option<u8> {
    ROCK_CLASSES.iter().find(|r| r.key == key).map(|r| r.id)
}

/// The smoothing grid's size and the Gaussian's σ there (PROXY: 0.5 coarse cell ≈ 3 km at 64² over 400 km).
pub const SMOOTH_GRID: usize = 1024;
pub const SMOOTH_SIGMA: f32 = 8.0;

/// ADR Finding 152 (amendment after run 1) — the contacts' terrain snapping: `a = clamp((z − z̄) / SNAP_SCALE_M, −1, 1)`
/// with z̄ the altitude smoothed at [`SNAP_SIGMA`] cells of the [`SMOOTH_GRID`] grid; hard classes take `+β·a`, the soft
/// rift fill `−β·a`, so a contact runs along the valleys (differential erosion, PROXY). A straight coarse edge stays
/// straight under any symmetric smoothing (run 1's G-blocks: 4.01 × the control); the terrain bends it.
pub const SNAP_BETA: f32 = 0.25;
pub const SNAP_SIGMA: f32 = 4.0;
pub const SNAP_SCALE_M: f32 = 100.0;

/// The terrain anomaly `a` on the [`SMOOTH_GRID`] grid (the 8×8-block mean of the altitude, minus its Gaussian).
#[must_use]
pub fn terrain_anomaly(z_m: &[f32], w: usize, h: usize) -> Vec<f32> {
    let n = SMOOTH_GRID;
    let (fx, fy) = ((w / n).max(1), (h / n).max(1));
    let zb: Vec<f32> = (0..n * n)
        .into_par_iter()
        .map(|q| {
            let (bx, by) = (q % n, q / n);
            let mut s = 0f32;
            for dy in 0..fy {
                for dx in 0..fx {
                    s += z_m[((by * fy + dy).min(h - 1)) * w + (bx * fx + dx).min(w - 1)];
                }
            }
            s / (fx * fy) as f32
        })
        .collect();
    let sm = gaussian_smooth(&zb, n, n, SNAP_SIGMA);
    zb.iter().zip(&sm).map(|(a, b)| ((a - b) / SNAP_SCALE_M).clamp(-1.0, 1.0)).collect()
}

/// A coarse mask (`nx·ny`) turned into a smooth field on the [`SMOOTH_GRID`] grid with the terrain's mapping.
#[must_use]
pub fn smooth_coarse(mask: &[bool], nx: usize, ny: usize, sample_origin: [f64; 2], sample_size: f64) -> Vec<f32> {
    let coarse = GridF32 { width: nx, height: ny, data: mask.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect() };
    let mid = upscale_k_to_hd(&coarse, SMOOTH_GRID, SMOOTH_GRID, sample_origin, sample_size);
    gaussian_smooth(&mid, SMOOTH_GRID, SMOOTH_GRID, SMOOTH_SIGMA)
}

/// A smooth field on the [`SMOOTH_GRID`] grid read at an HD cell (bilinear, clamped).
#[inline]
#[must_use]
pub fn at_smooth(f: &[f32], x: usize, y: usize, w: usize, h: usize) -> f32 {
    let n = SMOOTH_GRID;
    let u = (x as f32 * n as f32 / w as f32).min((n - 1) as f32);
    let v = (y as f32 * n as f32 / h as f32).min((n - 1) as f32);
    let (x0, y0) = (u.floor() as usize, v.floor() as usize);
    let (x1, y1) = ((x0 + 1).min(n - 1), (y0 + 1).min(n - 1));
    let (tx, ty) = (u - x0 as f32, v - y0 as f32);
    let a = f[y0 * n + x0] * (1.0 - tx) + f[y0 * n + x1] * tx;
    let b = f[y1 * n + x0] * (1.0 - tx) + f[y1 * n + x1] * tx;
    a * (1.0 - ty) + b * ty
}

/// The HD inputs of the rock grid.
pub struct RockInputs<'a> {
    pub w: usize,
    pub h: usize,
    /// Altitude (m) per HD cell, and the land test (above the sea and not a lake).
    pub z_m: &'a [f32],
    pub sea: &'a [bool],
    pub lake_map: &'a [u32],
    /// Lake ids that are endorheic.
    pub endorheic: &'a std::collections::HashSet<u32>,
    /// Flow accumulation (HD cells, REAL area) and the HD cell area (km²).
    pub acc: &'a [f32],
    pub cell_km: f32,
    pub precip_mm: &'a [f32],
    /// The coarse sources, `nx·ny`.
    pub nx: usize,
    pub ny: usize,
    pub craton: &'a [bool],
    pub belt: &'a [bool],
    pub rift: &'a [bool],
    pub sample_origin: [f64; 2],
    pub sample_size: f64,
    pub edifices: &'a [Edifice],
}

/// The per-cell slope (m/m), central differences on the metres.
#[must_use]
pub fn slope_field(z: &[f32], w: usize, h: usize, cell_m: f32) -> Vec<f32> {
    let mut s = vec![0f32; w * h];
    s.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        for (x, o) in row.iter_mut().enumerate() {
            let (xm, xp) = (x.saturating_sub(1), (x + 1).min(w - 1));
            let (ym, yp) = (y.saturating_sub(1), (y + 1).min(h - 1));
            let gx = (z[y * w + xp] - z[y * w + xm]) / (cell_m * (xp - xm).max(1) as f32);
            let gy = (z[yp * w + x] - z[ym * w + x]) / (cell_m * (yp - ym).max(1) as f32);
            *o = (gx * gx + gy * gy).sqrt();
        }
    });
    s
}

/// A box dilation of `m` by `r` cells (separable).
#[must_use]
pub fn dilate(m: &[bool], w: usize, h: usize, r: usize) -> Vec<bool> {
    let mut t = vec![false; w * h];
    t.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let src = &m[y * w..(y + 1) * w];
        let mut last: i64 = -1_000_000;
        for x in 0..w {
            if src[x] {
                last = x as i64;
            }
            row[x] = x as i64 - last <= r as i64;
        }
        let mut next: i64 = 1_000_000;
        for x in (0..w).rev() {
            if src[x] {
                next = x as i64;
            }
            row[x] = row[x] || next - x as i64 <= r as i64;
        }
    });
    let mut o = vec![false; w * h];
    o.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let (y0, y1) = (y.saturating_sub(r), (y + r).min(h - 1));
        for (x, v) in row.iter_mut().enumerate() {
            *v = (y0..=y1).any(|yy| t[yy * w + x]);
        }
    });
    o
}

/// The volcanic cells of the edifices: within 1.5 basal radii, above the base level (the median of the 1.3–1.6 radii
/// ring) by > 10 % of the edifice's height (PROXY). Returns the class (`ARC_VOLCANIC` / `BASALTIC_VOLCANIC`) or 0.
#[must_use]
pub fn volcanic_cells(edifices: &[Edifice], z_m: &[f32], w: usize, h: usize, sample_origin: [f64; 2], sample_size: f64, cell_km: f32) -> Vec<u8> {
    let mut out = vec![0u8; w * h];
    let (so, ss) = ([sample_origin[0] as f32, sample_origin[1] as f32], sample_size as f32);
    for e in edifices {
        let fx = (e.center_uv.0 - so[0]).rem_euclid(1.0) / ss;
        let fy = (e.center_uv.1 - so[1]).rem_euclid(1.0) / ss;
        if fx >= 1.0 || fy >= 1.0 {
            continue;
        }
        let (cx, cy) = (fx * w as f32, fy * h as f32);
        let rb = (e.basal_diameter_km * 0.5) / cell_km;
        if rb < 1.0 {
            continue;
        }
        let r_out = 1.6 * rb;
        let (i0, i1) = ((cx - r_out).floor().max(0.0) as usize, ((cx + r_out).ceil() as usize).min(w - 1));
        let (j0, j1) = ((cy - r_out).floor().max(0.0) as usize, ((cy + r_out).ceil() as usize).min(h - 1));
        let mut ring: Vec<f32> = Vec::new();
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d >= 1.3 * rb && d <= 1.6 * rb {
                    ring.push(z_m[j * w + i]);
                }
            }
        }
        if ring.is_empty() {
            continue;
        }
        ring.sort_by(f32::total_cmp);
        let base = ring[ring.len() / 2];
        let class = if e.setting == VolcanoSetting::Arc { ARC_VOLCANIC } else { BASALTIC_VOLCANIC };
        let thr = 0.1 * e.height_m;
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d <= 1.5 * rb && z_m[j * w + i] - base > thr {
                    out[j * w + i] = class;
                }
            }
        }
    }
    out
}

/// The loose-deposit and floodplain tests (PROXY): floodplain = slope < 0.02 and A ≥ 1 km²; lake margin = land within
/// 2 cells of a lake; slope foot = slope < 0.05 within 3 cells of a slope > 0.3.
pub const FLOODPLAIN_SLOPE: f32 = 0.02;
pub const FLOODPLAIN_AREA_KM2: f32 = 1.0;
pub const LAKE_MARGIN_CELLS: usize = 2;
pub const SLOPE_FOOT_MAX: f32 = 0.05;
pub const STEEP_SLOPE: f32 = 0.3;
pub const SLOPE_FOOT_CELLS: usize = 3;
/// Evaporites: an endorheic lake's bed, and its land within 1 km where the precipitation is below 400 mm (PROXY).
pub const EVAPORITE_MARGIN_KM: f32 = 1.0;
pub const EVAPORITE_PRECIP_MM: f32 = 400.0;

/// The floodplain mask (also the rules' `valley_floor`, PROXY: the construction's floor is not in `run_hd`).
#[must_use]
pub fn floodplain(slope: &[f32], acc: &[f32], cell_km: f32) -> Vec<bool> {
    let a_min = FLOODPLAIN_AREA_KM2 / (cell_km * cell_km);
    slope.par_iter().zip(acc.par_iter()).map(|(&s, &a)| s < FLOODPLAIN_SLOPE && a >= a_min).collect()
}

/// ADR Finding 154-S — the substratum's inputs: the coarse tectonic masks and the edifices, and the relief they are read
/// on (the contacts' terrain snapping, the cones' relief test). At the construction stage that relief is the
/// construction's input; at the end of the chain, the final field.
pub struct SubstratumInputs<'a> {
    pub w: usize,
    pub h: usize,
    pub z_m: &'a [f32],
    pub cell_km: f32,
    pub nx: usize,
    pub ny: usize,
    pub craton: &'a [bool],
    pub belt: &'a [bool],
    pub rift: &'a [bool],
    pub sample_origin: [f64; 2],
    pub sample_size: f64,
    pub edifices: &'a [Edifice],
}

/// ADR Finding 154-S — the SUBSTRATUM: one of craton, belt, basaltic or arc volcanic, rift fill, basement on every cell
/// (water included). The volcanic rock first, then the rift fill, the belt, the craton; the basement elsewhere.
#[must_use]
pub fn build_substratum(inp: &SubstratumInputs) -> Vec<u8> {
    let (w, h) = (inp.w, inp.h);
    let craton = smooth_coarse(inp.craton, inp.nx, inp.ny, inp.sample_origin, inp.sample_size);
    let belt = smooth_coarse(inp.belt, inp.nx, inp.ny, inp.sample_origin, inp.sample_size);
    let rift = smooth_coarse(inp.rift, inp.nx, inp.ny, inp.sample_origin, inp.sample_size);
    let anom = terrain_anomaly(inp.z_m, w, h);
    let volc = volcanic_cells(inp.edifices, inp.z_m, w, h, inp.sample_origin, inp.sample_size, inp.cell_km);
    let mut out = vec![BASEMENT; w * h];
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        for (x, o) in row.iter_mut().enumerate() {
            let c = y * w + x;
            *o = if volc[c] != 0 {
                volc[c]
            } else if at_smooth(&rift, x, y, w, h) - SNAP_BETA * at_smooth(&anom, x, y, w, h) > 0.5 {
                RIFT_FILL
            } else if at_smooth(&belt, x, y, w, h) + SNAP_BETA * at_smooth(&anom, x, y, w, h) > 0.5 {
                BELT
            } else if at_smooth(&craton, x, y, w, h) + SNAP_BETA * at_smooth(&anom, x, y, w, h) > 0.5 {
                CRATON
            } else {
                BASEMENT
            };
        }
    });
    out
}

/// ADR Finding 154-S — the surface's inputs, all from the end of the chain.
pub struct SurfaceInputs<'a> {
    pub w: usize,
    pub h: usize,
    pub sea: &'a [bool],
    pub lake_map: &'a [u32],
    pub endorheic: &'a std::collections::HashSet<u32>,
    pub acc: &'a [f32],
    pub cell_km: f32,
    pub precip_mm: &'a [f32],
}

/// ADR Finding 154-S — the SURFACE over the substratum: water (an endorheic lake's bed is evaporites), the evaporites
/// of the arid endorheic margins, the loose deposits (floodplain, lake margin, slope foot); the substratum elsewhere.
/// `slope` is [`slope_field`] of the final relief.
#[must_use]
pub fn apply_surface(substratum: &[u8], inp: &SurfaceInputs, slope: &[f32]) -> Vec<u8> {
    let (w, h) = (inp.w, inp.h);
    let n = w * h;
    let flood = floodplain(slope, inp.acc, inp.cell_km);
    let lake: Vec<bool> = inp.lake_map.iter().map(|&l| l != 0).collect();
    let lake_margin = dilate(&lake, w, h, LAKE_MARGIN_CELLS);
    let steep: Vec<bool> = slope.iter().map(|&s| s > STEEP_SLOPE).collect();
    let near_steep = dilate(&steep, w, h, SLOPE_FOOT_CELLS);
    let endo: Vec<bool> = inp.lake_map.iter().map(|&l| l != 0 && inp.endorheic.contains(&l)).collect();
    let endo_margin = dilate(&endo, w, h, (EVAPORITE_MARGIN_KM / inp.cell_km).round() as usize);
    let mut out = vec![NONE; n];
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        for (x, o) in row.iter_mut().enumerate() {
            let c = y * w + x;
            let l = inp.lake_map[c];
            if l != 0 {
                *o = if endo[c] { EVAPORITES } else { NONE };
                continue;
            }
            if inp.sea[c] {
                *o = NONE;
                continue;
            }
            *o = if endo_margin[c] && inp.precip_mm[c] < EVAPORITE_PRECIP_MM {
                EVAPORITES
            } else if flood[c] || lake_margin[c] || (near_steep[c] && slope[c] < SLOPE_FOOT_MAX) {
                LOOSE_DEPOSITS
            } else {
                substratum[c]
            };
        }
    });
    out
}

/// Build the rock grid: the surface over the substratum. `slope` is [`slope_field`] of `z_m`.
#[must_use]
pub fn build_rocks(inp: &RockInputs, slope: &[f32]) -> Vec<u8> {
    let sub = build_substratum(&SubstratumInputs {
        w: inp.w,
        h: inp.h,
        z_m: inp.z_m,
        cell_km: inp.cell_km,
        nx: inp.nx,
        ny: inp.ny,
        craton: inp.craton,
        belt: inp.belt,
        rift: inp.rift,
        sample_origin: inp.sample_origin,
        sample_size: inp.sample_size,
        edifices: inp.edifices,
    });
    apply_surface(
        &sub,
        &SurfaceInputs {
            w: inp.w,
            h: inp.h,
            sea: inp.sea,
            lake_map: inp.lake_map,
            endorheic: inp.endorheic,
            acc: inp.acc,
            cell_km: inp.cell_km,
            precip_mm: inp.precip_mm,
        },
        slope,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Finding 154-S, rule 13 — the rock grid is the surface over the substratum, cell for cell, and the substratum
    /// alone is a different grid where the surface acts. A 128² world: a slope to the sea, a lake, a craton on the west
    /// half, a rift band, no edifice. The per-cell order of the F152 rock grid is restated here (the old monolithic
    /// `build_rocks`), and both paths must give it. Negative control: a lake-margin land cell's substratum is not its
    /// final class (loose deposits), and the lake's cells are water over a hard substratum.
    #[test]
    fn the_rock_grid_is_the_surface_over_the_substratum() {
        let (w, h) = (128usize, 128usize);
        let n = w * h;
        let z_m: Vec<f32> = (0..n).map(|k| 400.0 - 3.0 * (k / w) as f32 + 20.0 * ((k % w) as f32 * 0.3).sin()).collect();
        let sea: Vec<bool> = z_m.iter().map(|&z| z <= 0.0).collect();
        let lake_cell = |k: usize| {
            let (x, y) = ((k % w) as f32, (k / w) as f32);
            ((x - 40.0).powi(2) + (y - 40.0).powi(2)).sqrt() < 6.0
        };
        let lake_map: Vec<u32> = (0..n).map(|k| if lake_cell(k) { 7 } else { 0 }).collect();
        let endorheic = std::collections::HashSet::new();
        let acc: Vec<f32> = (0..n).map(|k| ((k * 7919) % 3000) as f32).collect();
        let precip: Vec<f32> = vec![800.0; n];
        let (nx, ny) = (64usize, 64usize);
        let craton: Vec<bool> = (0..nx * ny).map(|c| c % nx < nx / 2).collect();
        let belt: Vec<bool> = vec![false; nx * ny];
        let rift: Vec<bool> = (0..nx * ny).map(|c| (40..46).contains(&(c % nx))).collect();
        let inp = RockInputs {
            w,
            h,
            z_m: &z_m,
            sea: &sea,
            lake_map: &lake_map,
            endorheic: &endorheic,
            acc: &acc,
            cell_km: 0.4,
            precip_mm: &precip,
            nx,
            ny,
            craton: &craton,
            belt: &belt,
            rift: &rift,
            sample_origin: [0.0, 0.0],
            sample_size: 1.0,
            edifices: &[],
        };
        let slope = slope_field(&z_m, w, h, 400.0);
        // F152's monolithic order, restated
        let cr = smooth_coarse(&craton, nx, ny, [0.0, 0.0], 1.0);
        let be = smooth_coarse(&belt, nx, ny, [0.0, 0.0], 1.0);
        let ri = smooth_coarse(&rift, nx, ny, [0.0, 0.0], 1.0);
        let an = terrain_anomaly(&z_m, w, h);
        let flood = floodplain(&slope, &acc, 0.4);
        let lake: Vec<bool> = lake_map.iter().map(|&l| l != 0).collect();
        let margin = dilate(&lake, w, h, LAKE_MARGIN_CELLS);
        let steep: Vec<bool> = slope.iter().map(|&s| s > STEEP_SLOPE).collect();
        let near_steep = dilate(&steep, w, h, SLOPE_FOOT_CELLS);
        let old: Vec<u8> = (0..n)
            .map(|c| {
                let (x, y) = (c % w, c / w);
                if lake_map[c] != 0 || sea[c] {
                    NONE
                } else if flood[c] || margin[c] || (near_steep[c] && slope[c] < SLOPE_FOOT_MAX) {
                    LOOSE_DEPOSITS
                } else if at_smooth(&ri, x, y, w, h) - SNAP_BETA * at_smooth(&an, x, y, w, h) > 0.5 {
                    RIFT_FILL
                } else if at_smooth(&be, x, y, w, h) + SNAP_BETA * at_smooth(&an, x, y, w, h) > 0.5 {
                    BELT
                } else if at_smooth(&cr, x, y, w, h) + SNAP_BETA * at_smooth(&an, x, y, w, h) > 0.5 {
                    CRATON
                } else {
                    BASEMENT
                }
            })
            .collect();
        let rocks = build_rocks(&inp, &slope);
        assert_eq!(rocks, old, "the split gives F152's grid");
        let sub = build_substratum(&SubstratumInputs { w, h, z_m: &z_m, cell_km: 0.4, nx, ny, craton: &craton, belt: &belt, rift: &rift, sample_origin: [0.0, 0.0], sample_size: 1.0, edifices: &[] });
        let surf = apply_surface(&sub, &SurfaceInputs { w, h, sea: &sea, lake_map: &lake_map, endorheic: &endorheic, acc: &acc, cell_km: 0.4, precip_mm: &precip }, &slope);
        assert_eq!(surf, rocks, "the surface over the substratum is the rock grid");
        for k in [RIFT_FILL, CRATON, BASEMENT] {
            assert!(sub.contains(&k), "the substratum has class {k}");
        }
        // the negative control: a lake-margin land cell and a lake cell
        let m = (0..n).find(|&c| margin[c] && !lake[c] && !sea[c]).expect("a margin cell");
        assert_eq!(rocks[m], LOOSE_DEPOSITS, "the margin is loose deposits at the surface");
        assert_ne!(sub[m], LOOSE_DEPOSITS, "its substratum is a rock");
        let l = (0..n).find(|&c| lake[c]).expect("a lake cell");
        assert_eq!(rocks[l], NONE);
        assert!(ROCK_CLASSES[sub[l] as usize].hardness > 0 || sub[l] == RIFT_FILL, "a lake has a substratum");
        assert!(sub.iter().all(|&k| ![NONE, LOOSE_DEPOSITS, EVAPORITES].contains(&k)), "the substratum holds no surface class");
    }
}
