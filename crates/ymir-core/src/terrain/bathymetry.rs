//! #submarine — bathymetric profile: re-map the ocean floor toward the real
//! plateau → continental-slope → abyssal-plain envelope, while PRESERVING the
//! existing relief as texture (no "onion" concentric rings).
//!
//! The diagnostic (`c1_closure_morphology::probe_submarine_relief`) measured the
//! ocean floor as a near-uniform mid-depth SLAB (~−2600 m, 95 % in one 1000-3000 m
//! band, deepest only ~3000 m, shelf ~1.5 %): Stein-Stein sets a coast→offshore
//! gradient by crustal age, but the young C1 ages cap it shallow, so the floor
//! never reaches the abyssal plain and there is almost no continental shelf.
//!
//! This RE-MAPS (does not replace) each ocean cell's depth toward the
//! distance-to-coast envelope — a shallow shelf near the coast, a steep slope,
//! then the deep abyssal plain — modulated by the cell's existing RELATIVE depth
//! deviation so the abyss keeps relief and the shelf stays smooth. Only cells
//! below sea level are touched and they stay strictly below it, so the coastline,
//! the land/sea mask, and the emergent-land hypsometry are invariant.

use std::collections::VecDeque;

use crate::grid::GridF32;

/// Plateau→slope→abyss bathymetric envelope, anchored on real passive-margin
/// morphology (shelf ~−130 m / ~70 km wide, break ~−200 m; steep slope to the
/// abyssal plain ~−4500 m). All depths are POSITIVE metres below sea level.
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct BathymetryProfile {
    /// Depth at the coastline — the shelf is shallow but submerged (stays ocean).
    pub shelf_min_depth_m: f32,
    /// Depth at the shelf break (the shelf edge; Earth ~130-200 m).
    pub shelf_break_depth_m: f32,
    /// Shelf width (km) — coast to shelf break (~70 km on Earth's passive margins).
    pub shelf_width_km: f32,
    /// Continental-slope width (km) — the steep drop from the shelf break to the
    /// abyssal plain (~80 km).
    pub slope_width_km: f32,
    /// Abyssal-plain depth (m) far offshore (~−4500 m).
    pub abyss_depth_m: f32,
    /// Relative texture amplitude. The cell's existing depth deviation from the
    /// ocean mean is carried as a MULTIPLICATIVE modulation `envelope·(1+texture·
    /// dev)` — the abyss keeps relief, the shelf stays smooth, and the result is
    /// never flat rings (the FBM/Stein-Stein structure survives). `0.0` = pure
    /// envelope (the "onion" to avoid); `~1.0` = full existing texture.
    pub texture: f32,
}

impl Default for BathymetryProfile {
    fn default() -> Self {
        Self {
            shelf_min_depth_m: 20.0,
            shelf_break_depth_m: 200.0,
            // #submarine — DELIBERATE GAMEPLAY CHOICE, not a calibration miss. At
            // 30 km the shelf lands at ~11.6 % of ocean area (measured, 6 seeds) —
            // WIDER than Earth's real ~7-8 %, kept on purpose for playable coastal
            // waters (coastal navigation / fishing in the Living Landz sandbox).
            // This is the generator's ONE gameplay-tuned departure from the real
            // anchor (everything else is anchored-not-tuned); it is named here so
            // the gap reads as intent, not drift. `shelf_width_km` is the knob:
            // ~20 km tightens the shelf toward the real ~7-8 %; 70 km gave ~21 %
            // (too wide for these regional maps' dense coastlines).
            shelf_width_km: 30.0,
            slope_width_km: 80.0,
            abyss_depth_m: 4500.0,
            texture: 1.0,
        }
    }
}

fn smoothstep(t: f32) -> f32 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Envelope depth (m, positive down) at distance-to-coast `d_km`: a gentle shelf
/// ramp (shelf_min → shelf_break), then a steep slope smoothstepping to the
/// abyssal plain, then flat abyss.
pub fn envelope_depth(d_km: f32, p: &BathymetryProfile) -> f32 {
    if d_km < p.shelf_width_km {
        let t = d_km / p.shelf_width_km.max(1e-3);
        p.shelf_min_depth_m + (p.shelf_break_depth_m - p.shelf_min_depth_m) * t
    } else {
        let t = (d_km - p.shelf_width_km) / p.slope_width_km.max(1e-3);
        p.shelf_break_depth_m + (p.abyss_depth_m - p.shelf_break_depth_m) * smoothstep(t)
    }
}

/// Multi-source BFS distance (in cells, 4-connectivity) from the coast over the
/// OCEAN cells. Coast = an ocean cell 4-adjacent to LAND only — the domain edge
/// is NOT a coast (the ocean continues beyond the bounded map, so edge water is
/// open ocean / abyss, not shelf). Non-ocean cells get `0` (unused); ocean cells
/// with no land in their connected component keep `u32::MAX` (→ abyss envelope).
fn dist_to_coast_cells(ocean: &[bool], w: usize, h: usize) -> Vec<u32> {
    let mut dist = vec![u32::MAX; w * h];
    let mut q: VecDeque<(usize, usize)> = VecDeque::new();
    let idx = |i: usize, j: usize| j * w + i;
    for j in 0..h {
        for i in 0..w {
            if !ocean[idx(i, j)] {
                continue;
            }
            let mut coast = false;
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni >= 0
                    && nj >= 0
                    && (ni as usize) < w
                    && (nj as usize) < h
                    && !ocean[idx(ni as usize, nj as usize)]
                {
                    coast = true;
                }
            }
            if coast {
                dist[idx(i, j)] = 0;
                q.push_back((i, j));
            }
        }
    }
    while let Some((i, j)) = q.pop_front() {
        let d = dist[idx(i, j)];
        for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
            let (ni, nj) = (i as i32 + di, j as i32 + dj);
            if ni >= 0 && nj >= 0 && (ni as usize) < w && (nj as usize) < h {
                let nk = idx(ni as usize, nj as usize);
                if ocean[nk] && dist[nk] == u32::MAX {
                    dist[nk] = d + 1;
                    q.push_back((ni as usize, nj as usize));
                }
            }
        }
    }
    dist
}

/// Re-map ocean-cell depths toward the plateau→slope→abyss envelope, preserving
/// each cell's existing RELATIVE deviation as texture. In place on the normalised
/// heightmap. `sea_norm` = sea-level normalised value (0.5); `depth_per_norm`
/// converts a norm step below sea to metres (`2·HALF·depth_scale_m`, e.g. 11300);
/// `km_per_cell` = horizontal scale.
///
/// INVARIANTS: only `norm <= sea_norm` (ocean) cells are touched, and they stay
/// strictly below `sea_norm` (depth ≥ `shelf_min_depth_m`), so the coastline /
/// land-sea mask and the emergent-land hypsometry are unchanged.
pub fn apply_bathymetry_profile(
    h: &mut GridF32,
    sea_norm: f32,
    depth_per_norm: f32,
    km_per_cell: f32,
    p: &BathymetryProfile,
) {
    let (w, ht) = (h.width, h.height);
    if depth_per_norm.abs() < 1e-6 {
        return;
    }
    let ocean: Vec<bool> = h.data.iter().map(|&n| n <= sea_norm).collect();
    let dist = dist_to_coast_cells(&ocean, w, ht);

    // Existing ocean-depth mean (m) → the relative-texture reference.
    let mut sum = 0.0f64;
    let mut cnt = 0u64;
    for k in 0..w * ht {
        if ocean[k] {
            sum += ((sea_norm - h.data[k]) * depth_per_norm) as f64;
            cnt += 1;
        }
    }
    if cnt == 0 {
        return;
    }
    let mean_depth = (sum / cnt as f64) as f32;
    let floor_depth = sea_norm * depth_per_norm; // norm 0 → the deepest (~5650 m)
    let min_depth = p.shelf_min_depth_m.max(1.0);

    for k in 0..w * ht {
        if !ocean[k] {
            continue;
        }
        let cur_depth = (sea_norm - h.data[k]) * depth_per_norm;
        let dev = if mean_depth.abs() > 1e-3 { (cur_depth - mean_depth) / mean_depth } else { 0.0 };
        let env = envelope_depth(dist[k] as f32 * km_per_cell, p);
        let new_depth = (env * (1.0 + p.texture * dev)).clamp(min_depth, floor_depth);
        h.data[k] = sea_norm - new_depth / depth_per_norm;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn envelope_is_monotone_plateau_to_abyss() {
        let p = BathymetryProfile::default();
        // Shallow at the coast, deep offshore, monotone non-decreasing.
        assert!((envelope_depth(0.0, &p) - p.shelf_min_depth_m).abs() < 1.0);
        assert!(envelope_depth(35.0, &p) > envelope_depth(0.0, &p)); // shelf ramp
        assert!(envelope_depth(150.0, &p) > envelope_depth(70.0, &p)); // slope
        assert!((envelope_depth(400.0, &p) - p.abyss_depth_m).abs() < 1.0); // abyss
        let mut prev = 0.0;
        for d in 0..400 {
            let e = envelope_depth(d as f32, &p);
            assert!(e + 1e-3 >= prev, "envelope must be non-decreasing at {d} km");
            prev = e;
        }
    }

    #[test]
    fn remap_keeps_coastline_and_deepens_offshore() {
        // A 1-D-ish map: left third land (norm 0.6), rest a shallow ocean slab
        // (norm 0.46 ≈ 450 m). After the re-map, no ocean cell emerges, all land
        // stays land, and the far-offshore floor is much deeper than near-coast.
        let (w, ht) = (200usize, 4usize);
        let mut h = GridF32::new(w, ht, 0.0);
        for j in 0..ht {
            for i in 0..w {
                h.set(i, j, if i < 60 { 0.60 } else { 0.46 });
            }
        }
        let land_before: usize = h.data.iter().filter(|&&n| n > 0.5).count();
        let p = BathymetryProfile::default();
        apply_bathymetry_profile(&mut h, 0.5, 11300.0, 2.0, &p); // 2 km/cell

        let land_after: usize = h.data.iter().filter(|&&n| n > 0.5).count();
        assert_eq!(land_before, land_after, "coastline / land-sea mask must not move");
        assert!(h.data.iter().all(|&n| n.is_finite()));
        // Near-coast (i=62) shallow shelf, far-offshore (i=199) deep abyss.
        let depth = |i: usize| (0.5 - h.get(i as i32, 0)) * 11300.0;
        assert!(depth(62) < 600.0, "near-coast should be shelf-shallow, got {}", depth(62));
        assert!(depth(199) > 2500.0, "far offshore should be abyssal, got {}", depth(199));
        assert!(depth(199) > depth(62) + 1500.0, "must deepen strongly coast→offshore");
        // Every ocean cell stays submerged.
        for i in 60..w {
            assert!(h.get(i as i32, 0) < 0.5, "ocean cell {i} must stay below sea");
        }
    }
}
