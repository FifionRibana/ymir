//! ADR Finding 160-V -- the rivers and the lakes of a cascade level (the Cascade window, the benches' images and the
//! depression instrument). Read-only on the field: the D8 routing of `compute_flow` on the priority-flood-filled
//! surface (sea 0 m), the Strahler order on the network above a drainage threshold, and the lakes as the filled
//! depressions (filled − z > 0.5 m, 8-connected, ≥ 2 cells) with their area and spill.

use crate::grid::GridF32;
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};

/// One lake: a filled depression.
#[derive(Clone, Debug)]
pub struct Lake {
    pub cells: usize,
    pub area_km2: f32,
    /// The spill level (m), the filled surface.
    pub spill_m: f32,
    /// The deepest cell.
    pub lowest: usize,
}

/// The hydrology of one field (m, sea at 0 m).
#[derive(Clone, Debug)]
pub struct Hydro {
    pub n: usize,
    pub cell_km: f32,
    /// D8 drainage area (km²).
    pub acc_km2: Vec<f32>,
    /// D8 direction (`DIR_NONE` on the sea and the outlets).
    pub dir: Vec<u8>,
    /// Strahler order on the network (acc ≥ the threshold), 0 off it.
    pub strahler: Vec<u8>,
    /// 0 = no lake, i = lake i − 1.
    pub lake_id: Vec<u32>,
    pub lakes: Vec<Lake>,
    /// The network threshold (km²).
    pub river_min_km2: f32,
}

/// The declared river threshold: max(1 km², 4 cells).
pub fn river_threshold_km2(cell_km: f32) -> f32 {
    (4.0 * cell_km * cell_km).max(1.0)
}

pub fn hydrology(z_m: &GridF32, cell_km: f32) -> Hydro {
    let n = z_m.width;
    let flow = compute_flow(z_m, &FlowConfig { sea_level: 0.0, ..Default::default() });
    let c2 = cell_km * cell_km;
    let acc_km2: Vec<f32> = flow.accumulation.data.iter().map(|a| a * c2).collect();
    let thr = river_threshold_km2(cell_km);
    let recv = |k: usize| -> Option<usize> {
        let d = flow.direction[k];
        if d == DIR_NONE || z_m.data[k] <= 0.0 {
            return None;
        }
        let (x, y) = ((k % n) as i32 + D8_DX[d as usize], (k / n) as i32 + D8_DY[d as usize]);
        (x >= 0 && y >= 0 && x < n as i32 && y < n as i32).then(|| y as usize * n + x as usize)
    };
    // Strahler, upstream first (D8 accumulation grows strictly downstream)
    let mut order: Vec<usize> = (0..n * n).filter(|&k| z_m.data[k] > 0.0 && acc_km2[k] >= thr).collect();
    order.sort_unstable_by(|&a, &b| acc_km2[a].partial_cmp(&acc_km2[b]).unwrap().then(a.cmp(&b)));
    let (mut best, mut count) = (vec![0u8; n * n], vec![0u8; n * n]);
    let mut strahler = vec![0u8; n * n];
    for &k in &order {
        let s = if best[k] == 0 { 1 } else if count[k] >= 2 { best[k].saturating_add(1) } else { best[k] };
        strahler[k] = s;
        if let Some(r) = recv(k) {
            if s > best[r] {
                best[r] = s;
                count[r] = 1;
            } else if s == best[r] {
                count[r] = count[r].saturating_add(1);
            }
        }
    }
    // lakes
    let wet: Vec<bool> = (0..n * n).map(|k| z_m.data[k] > 0.0 && flow.filled.data[k] - z_m.data[k] > 0.5).collect();
    let mut lake_id = vec![0u32; n * n];
    let mut lakes = Vec::new();
    let mut stack = Vec::new();
    for s in 0..n * n {
        if !wet[s] || lake_id[s] != 0 {
            continue;
        }
        let id = lakes.len() as u32 + 1;
        let mut cells = Vec::new();
        lake_id[s] = id;
        stack.push(s);
        while let Some(k) = stack.pop() {
            cells.push(k);
            let (x, y) = ((k % n) as i32, (k / n) as i32);
            for d in 0..8 {
                let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                if nx < 0 || ny < 0 || nx >= n as i32 || ny >= n as i32 {
                    continue;
                }
                let j = ny as usize * n + nx as usize;
                if wet[j] && lake_id[j] == 0 {
                    lake_id[j] = id;
                    stack.push(j);
                }
            }
        }
        if cells.len() < 2 {
            for &k in &cells {
                lake_id[k] = 0;
            }
            continue;
        }
        let lowest = *cells.iter().min_by(|&&a, &&b| z_m.data[a].partial_cmp(&z_m.data[b]).unwrap().then(a.cmp(&b))).unwrap();
        let spill = cells.iter().map(|&k| flow.filled.data[k]).fold(f32::MIN, f32::max);
        lakes.push(Lake { cells: cells.len(), area_km2: cells.len() as f32 * c2, spill_m: spill, lowest });
    }
    // (a dropped single cell's id is reused by the next component, so the ids stay dense)
    Hydro { n, cell_km, acc_km2, dir: flow.direction, strahler, lake_id, lakes, river_min_km2: thr }
}

/// Draw the rivers (Strahler ≥ `min_order`; trunks ≥ 100 km² in a darker blue) and the lakes' outline onto RGBA rows
/// laid out **north up** (row 0 = the field's last row).
pub fn draw(rgba: &mut [u8], hy: &Hydro, min_order: u8, rivers: bool, lakes: bool) {
    let n = hy.n;
    let put = |rgba: &mut [u8], k: usize, c: [u8; 3]| {
        let o = ((n - 1 - k / n) * n + k % n) * 4;
        rgba[o..o + 3].copy_from_slice(&c);
    };
    if lakes {
        for k in 0..n * n {
            if hy.lake_id[k] == 0 {
                continue;
            }
            let (x, y) = (k % n, k / n);
            let edge = x == 0 || y == 0 || x + 1 == n || y + 1 == n || [k - 1, k + 1, k - n, k + n].iter().any(|&j| hy.lake_id[j] != hy.lake_id[k]);
            put(rgba, k, if edge { [0x4F, 0xC3, 0xF7] } else { [0x9C, 0xD8, 0xF0] });
        }
    }
    if rivers {
        for k in 0..n * n {
            let s = hy.strahler[k];
            // a river is not drawn across a lake (the D8 routing over a filled flat is a construction, not a river)
            if s == 0 || s < min_order || (lakes && hy.lake_id[k] != 0) {
                continue;
            }
            put(rgba, k, if hy.acc_km2[k] >= 100.0 { [0x0D, 0x3B, 0x8C] } else { [0x2E, 0x6F, 0xD8] });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Finding 162-0 -- the cascade reads no SIGNIFIED quantity (F147: 831 junctions announced, 465 real, from
    /// mixing the ×7.5 scale's units): its sources never name the geographic scale ratio, and `hydrology`'s areas are
    /// cells × the real cell² exactly.
    #[test]
    fn the_cascade_reads_no_signified_quantity() {
        let forbidden = concat!("geo_scale", "_ratio");
        for (name, src) in [
            ("cascade.rs", include_str!("../cascade.rs")),
            ("amplify.rs", include_str!("amplify.rs")),
            ("hydro.rs", include_str!("hydro.rs")),
            ("measure.rs", include_str!("measure.rs")),
            ("predict.rs", include_str!("predict.rs")),
        ] {
            assert!(!src.contains(forbidden), "{name} names the signified scale");
        }
        let n = 32;
        let mut g = GridF32::new(n, n, -10.0);
        for y in 1..n {
            for x in 0..n {
                g.data[y * n + x] = 3.0 * y as f32 + 0.01 * x as f32;
            }
        }
        let cell = 0.390625;
        let hy = hydrology(&g, cell);
        let flow = crate::terrain::flow::compute_flow(&g, &crate::terrain::flow::FlowConfig { sea_level: 0.0, ..Default::default() });
        for k in 0..n * n {
            assert_eq!(hy.acc_km2[k], flow.accumulation.data[k] * cell * cell, "areas are real (cells × cell²) at {k}");
        }
    }

    /// Two valleys meeting and a closed basin: the confluence raises the Strahler order, the basin is a lake whose spill
    /// is its rim, and a field with no depression has no lake (negative control).
    #[test]
    fn strahler_and_lakes_on_a_synthetic_field() {
        let n = 64;
        let mut g = GridF32::new(n, n, 0.0);
        for y in 0..n {
            for x in 0..n {
                // a valley along x = 32 falling south (y = 0 is the sea row), its flanks steeper than its floor so
                // every row drains to the axis: two lateral channels meet on it
                let v = 2.0 * y as f32 + 10.0 * ((x as f32 - 32.0).abs());
                g.data[y * n + x] = if y == 0 { -10.0 } else { v };
            }
        }
        let flat = hydrology(&g, 1.0);
        assert!(flat.lakes.is_empty(), "no depression, no lake");
        assert!(flat.strahler.iter().any(|&s| s >= 2), "the confluence must reach order 2");
        // a closed basin around (20, 40), 30 m deep
        for y in 36..45 {
            for x in 16..25 {
                g.data[y * n + x] -= 30.0 + 5.0 * (4.0 - ((x as f32 - 20.0).abs().max((y as f32 - 40.0).abs())));
            }
        }
        let hy = hydrology(&g, 1.0);
        assert!(!hy.lakes.is_empty(), "the basin must be a lake");
        let l = &hy.lakes[0];
        assert!(l.cells >= 2 && l.spill_m > g.data[l.lowest]);
    }
}
