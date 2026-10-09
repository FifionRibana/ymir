//! ADR Finding 163-I -- **what the author sees, named** (F162's tight crop: straight trunks, a fishbone of parallel
//! tributaries, flat interfluves):
//!
//! 1. `a_dir` -- the share of the order-≥ 2 rivers' length whose 8-step chords lie within ±5° of 0°, 45°, 90° or
//!    135° (isotropic directions read 10 / 45 = 22.2 %);
//! 2. `sinuosity` -- the mean path length / chord of D8 walks of W = 2, 5 and 10 km down the network draining
//!    ≥ 10 km² (a window shorter than 4 cells is not read);
//! 3. `flat` -- the share of the land with a slope < 0.05 that lies farther than 3 cells (Chebyshev) from an
//!    order-≥ 2 river.
//!
//! The network is `super::hydro` (F162's: D8 on the filled surface, REAL km², the threshold max(1 km², 4 cells)); the
//! cells read are the land (z > 0 m) minus `excl` (the D5 cells) minus a 2-cell border. Declared in
//! `f163_declared.md` (I), with the controls below.

use super::hydro::{Hydro, hydrology};
use crate::grid::GridF32;
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
use rayon::prelude::*;

/// I1: a stream is cut into pieces of this many D8 steps …
pub const PIECE_STEPS: usize = 8;
/// … and a last piece shorter than this is dropped.
pub const PIECE_MIN_STEPS: usize = 4;
/// I1: a chord within this many degrees of a D8 direction is aligned.
pub const ALIGN_DEG: f32 = 5.0;
/// I1, reported: the trunks (the author's « grands cours d'eau »).
pub const TRUNK_KM2: f32 = 100.0;
/// I2: the network the windows start from.
pub const SINU_KM2: f32 = 10.0;
/// I2: the windows (km).
pub const WINDOWS_KM: [f32; 3] = [2.0, 5.0, 10.0];
/// I2: a window shorter than this many cells is not read.
pub const WINDOW_MIN_CELLS: f32 = 4.0;
/// I3: x (cells, Chebyshev) and y (m/m).
pub const FLAT_X_CELLS: u32 = 3;
pub const FLAT_Y: f32 = 0.05;

/// One reading of the three instruments.
#[derive(Clone, Debug, Default)]
pub struct Planform {
    /// I1: the aligned share of the order-≥ 2 pieces' chord length (0–1).
    pub a_dir: f32,
    /// I1 on the trunks only (pieces starting on ≥ 100 km²), reported.
    pub a_dir_trunks: f32,
    pub pieces: usize,
    pub trunk_pieces: usize,
    /// I2: S_W for [`WINDOWS_KM`], NaN where the window is under 4 cells.
    pub sinuosity: [f32; 3],
    pub windows: [usize; 3],
    /// I3: the share (0–1).
    pub flat: f32,
}

impl Planform {
    /// What R reads for I2 (F163-I2, amended): the excess S_W − 1.
    pub fn excess(&self, i: usize) -> f32 {
        self.sinuosity[i] - 1.0
    }
}

/// The three instruments on a field in metres (sea ≤ 0 m).
pub fn planform(z: &GridF32, cell_km: f32, excl: Option<&[bool]>) -> Planform {
    let hy = hydrology(z, cell_km);
    planform_of(z, &hy, excl)
}

/// [`planform`] with the hydrology already computed (`hy` from `hydrology(z, cell_km)`).
pub fn planform_of(z: &GridF32, hy: &Hydro, excl: Option<&[bool]>) -> Planform {
    let n = z.width;
    let cell_m = hy.cell_km * 1000.0;
    let read: Vec<bool> = (0..n * n)
        .map(|k| {
            let (x, y) = (k % n, k / n);
            x >= 2 && y >= 2 && x + 2 < n && y + 2 < n && z.data[k] > 0.0 && !excl.is_some_and(|m| m[k])
        })
        .collect();
    let recv = |k: usize| -> Option<usize> {
        let d = hy.dir[k];
        if d == DIR_NONE || z.data[k] <= 0.0 {
            return None;
        }
        let (x, y) = ((k % n) as i32 + D8_DX[d as usize], (k / n) as i32 + D8_DY[d as usize]);
        (x >= 0 && y >= 0 && x < n as i32 && y < n as i32).then(|| y as usize * n + x as usize)
    };
    let mut out = Planform::default();

    // ── I1: the directions ──
    let mut same_donor = vec![false; n * n];
    for k in 0..n * n {
        let s = hy.strahler[k];
        if s >= 2 {
            if let Some(r) = recv(k) {
                if hy.strahler[r] == s {
                    same_donor[r] = true;
                }
            }
        }
    }
    let (mut all, mut aligned, mut t_all, mut t_aligned) = (0f64, 0f64, 0f64, 0f64);
    for k in 0..n * n {
        let s = hy.strahler[k];
        if s < 2 || same_donor[k] {
            continue;
        }
        let mut path = vec![k];
        let mut cur = k;
        while let Some(r) = recv(cur) {
            if hy.strahler[r] != s {
                break;
            }
            path.push(r);
            cur = r;
        }
        let mut i = 0;
        while path.len() - 1 >= i + PIECE_MIN_STEPS {
            let j = (i + PIECE_STEPS).min(path.len() - 1);
            if path[i..=j].iter().all(|&c| read[c]) {
                let dx = (path[j] % n) as f32 - (path[i] % n) as f32;
                let dy = (path[j] / n) as f32 - (path[i] / n) as f32;
                let len = (dx * dx + dy * dy).sqrt() as f64;
                let phi = dy.abs().atan2(dx.abs()).to_degrees() % 45.0;
                let al = phi.min(45.0 - phi) <= ALIGN_DEG;
                all += len;
                out.pieces += 1;
                if al {
                    aligned += len;
                }
                if hy.acc_km2[path[i]] >= TRUNK_KM2 {
                    t_all += len;
                    out.trunk_pieces += 1;
                    if al {
                        t_aligned += len;
                    }
                }
            }
            i = j;
        }
    }
    out.a_dir = if all > 0.0 { (aligned / all) as f32 } else { f32::NAN };
    out.a_dir_trunks = if t_all > 0.0 { (t_aligned / t_all) as f32 } else { f32::NAN };

    // ── I2: the sinuosity ──
    let starts: Vec<usize> = (0..n * n).filter(|&k| read[k] && hy.acc_km2[k] >= SINU_KM2).collect();
    for (wi, &w_km) in WINDOWS_KM.iter().enumerate() {
        let w_cells = w_km / hy.cell_km;
        if w_cells < WINDOW_MIN_CELLS {
            out.sinuosity[wi] = f32::NAN;
            continue;
        }
        let s: Vec<f64> = starts
            .par_iter()
            .filter_map(|&k| {
                let (mut cur, mut len) = (k, 0f32);
                while len < w_cells {
                    let r = recv(cur).filter(|&r| read[r])?;
                    let diag = (r % n != cur % n) && (r / n != cur / n);
                    len += if diag { std::f32::consts::SQRT_2 } else { 1.0 };
                    cur = r;
                }
                let dx = (cur % n) as f32 - (k % n) as f32;
                let dy = (cur / n) as f32 - (k / n) as f32;
                let chord = (dx * dx + dy * dy).sqrt();
                (chord > 0.0).then(|| (len / chord) as f64)
            })
            .collect();
        out.windows[wi] = s.len();
        out.sinuosity[wi] = if s.is_empty() { f32::NAN } else { (s.iter().sum::<f64>() / s.len() as f64) as f32 };
    }

    // ── I3: the flat interfluves ──
    let mut dist = vec![u32::MAX; n * n];
    let mut front: Vec<usize> = (0..n * n).filter(|&k| hy.strahler[k] >= 2 && z.data[k] > 0.0).collect();
    for &k in &front {
        dist[k] = 0;
    }
    for d in 1..=FLAT_X_CELLS + 1 {
        let mut next = Vec::new();
        for &k in &front {
            let (x, y) = ((k % n) as i32, (k / n) as i32);
            for dd in 0..8 {
                let (nx, ny) = (x + D8_DX[dd], y + D8_DY[dd]);
                if nx < 0 || ny < 0 || nx >= n as i32 || ny >= n as i32 {
                    continue;
                }
                let j = ny as usize * n + nx as usize;
                if dist[j] == u32::MAX {
                    dist[j] = d;
                    next.push(j);
                }
            }
        }
        front = next;
    }
    let zl = |k: usize| z.data[k].max(0.0);
    let (mut cells, mut flat) = (0usize, 0usize);
    for k in 0..n * n {
        if !read[k] {
            continue;
        }
        let gx = (zl(k + 1) - zl(k - 1)) / (2.0 * cell_m);
        let gy = (zl(k + n) - zl(k - n)) / (2.0 * cell_m);
        cells += 1;
        if (gx * gx + gy * gy).sqrt() < FLAT_Y && dist[k] > FLAT_X_CELLS {
            flat += 1;
        }
    }
    out.flat = if cells > 0 { flat as f32 / cells as f32 } else { f32::NAN };
    out
}

// ───────────────────────────── the controls (F163-I, metres, sea ≤ 0 m) ─────────────────────────────

/// A deterministic hash in [0, 1).
pub fn hash01(k: u64, seed: u64) -> f32 {
    let mut x = k.wrapping_add(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15)).wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^= x >> 31;
    (x >> 40) as f32 / (1u64 << 24) as f32
}

/// I1's « high » control: **the D8 network of a smooth inclined plane**. It falls 2 % toward the west with its
/// gradient at 10° from the x axis, so D8 runs due west while the true gradient is at 10° (the D8 artefact); a white
/// hash noise of 0.2 × the per-cell drop sends a few steps diagonal so the lines merge. The sea is the 2 western
/// columns and the 2 southern rows.
pub fn control_plane(n: usize, cell_km: f32) -> GridF32 {
    let cell_m = cell_km * 1000.0;
    let (s, a) = (0.02f32, 10f32.to_radians());
    let mut g = GridF32::new(n, n, 0.0);
    for y in 0..n {
        for x in 0..n {
            let k = y * n + x;
            g.data[k] = if x < 2 || y < 2 {
                -10.0
            } else {
                10.0 + s * cell_m * (x as f32 * a.cos() + y as f32 * a.sin()) + 0.2 * s * cell_m * hash01(k as u64, 163)
            };
        }
    }
    g
}

/// I2's controls: parallel V valleys (5 % down the axis, 10 % across, 8 cells apart, so a flank line drains < 4 cells)
/// along a straight line at `angle_deg` to the grid, or, with `sine = Some((a, p))` (cells), along y = a·sin(2πx/p)
/// at 0° and 16 cells apart. The sea is where the axis coordinate is below 4 cells.
pub fn control_valleys(n: usize, cell_km: f32, angle_deg: f32, sine: Option<(f32, f32)>) -> GridF32 {
    let cell_m = cell_km * 1000.0;
    let a = angle_deg.to_radians();
    let spacing = if sine.is_some() { 16.0 } else { 8.0 };
    let mut g = GridF32::new(n, n, 0.0);
    for y in 0..n {
        for x in 0..n {
            let k = y * n + x;
            let (fx, fy) = (x as f32, y as f32);
            let (along, mut across) = match sine {
                Some(_) => (fx, fy),
                None => (fx * a.cos() + fy * a.sin(), -fx * a.sin() + fy * a.cos()),
            };
            if let Some((amp, p)) = sine {
                across -= amp * (std::f32::consts::TAU * fx / p).sin();
            }
            let t = across / spacing;
            let d = (t - t.round()).abs() * spacing;
            g.data[k] = if along < 4.0 {
                -10.0
            } else {
                10.0 + cell_m * (0.05 * along + 0.10 * d) + 0.01 * hash01(k as u64, 1631)
            };
        }
    }
    g
}

/// I2's sinuous reference: S_W walked on the CONTINUOUS curve y = a·sin(2πx/p) (cells), windows of `w_cells` of arc
/// length started every 1/16 cell over one period, the mean of arc / chord.
pub fn sine_reference(a: f32, p: f32, w_cells: f32) -> f32 {
    let curve = |x: f64| a as f64 * (std::f64::consts::TAU * x / p as f64).sin();
    let dx = 1.0 / 64.0;
    let mut sum = 0f64;
    let starts = (16.0 * p) as usize;
    for i in 0..starts {
        let x0 = i as f64 / 16.0;
        let (mut x, mut len) = (x0, 0f64);
        while len < w_cells as f64 {
            let (y0, y1) = (curve(x), curve(x + dx));
            len += (dx * dx + (y1 - y0) * (y1 - y0)).sqrt();
            x += dx;
        }
        let chord = ((x - x0).powi(2) + (curve(x) - curve(x0)).powi(2)).sqrt();
        sum += len / chord;
    }
    (sum / starts as f64) as f32
}

/// I3's « high » control: a terraced field. A plane falling 3 % toward the south, quantised into 100 m treads (each
/// tilted 1 % toward its nearest valley, steps for risers), with V valleys 60 m under the tread and 8 cells wide every
/// 64 cells carrying the order-≥ 2 rivers to the 2 southern sea rows.
pub fn control_terraces(n: usize, cell_km: f32) -> GridF32 {
    let cell_m = cell_km * 1000.0;
    let mut g = GridF32::new(n, n, 0.0);
    for y in 0..n {
        for x in 0..n {
            let k = y * n + x;
            if y < 2 {
                g.data[k] = -10.0;
                continue;
            }
            let plane = 0.03 * cell_m * y as f32;
            // each tread tilts 1 % toward its nearest valley (flat: < y), so its water runs there in parallel lines
            let to_valley = ((x % 64) as f32 - 32.0).abs();
            let tread = 100.0 * (plane / 100.0).floor() + 20.0 + 0.01 * cell_m * to_valley;
            // the V valley: 60 m under the local tread (so its floor steps down with the treads), 15 m per cell up its
            // flanks, 8 cells wide
            let valley = tread - 60.0 + 15.0 * to_valley;
            g.data[k] = tread.min(valley) + 0.001 * hash01(k as u64, 1632);
        }
    }
    g
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The three instruments on their controls: a straight valley at 0° reads S = 1, the sinuous one reads near its
    /// continuous curve, the D8 plane reads more aligned than isotropic, the terraced field reads mostly flat
    /// interfluves. The ×1.5 validations against Corsica are the bench's (`f163_cascade`).
    #[test]
    fn the_controls_read_as_built() {
        let (n, c) = (256, 0.195);
        let straight = planform(&control_valleys(n, c, 0.0, None), c, None);
        assert!((straight.sinuosity[1] - 1.0).abs() < 0.01, "straight S_5 {}", straight.sinuosity[1]);
        let sine = planform(&control_valleys(n, c, 0.0, Some((4.0, 32.0))), c, None);
        let want = sine_reference(4.0, 32.0, 5.0 / c);
        assert!(sine.sinuosity[1] > 1.0 + 0.5 * (want - 1.0), "sine S_5 {} against {want}", sine.sinuosity[1]);
        let plane = planform(&control_plane(n, c), c, None);
        assert!(plane.pieces > 0 && plane.a_dir > 0.3, "plane A_dir {} ({} pieces)", plane.a_dir, plane.pieces);
        let terr = planform(&control_terraces(n, c), c, None);
        assert!(terr.flat > 0.4, "terraces F {}", terr.flat);
        // a window under 4 cells is not read
        let coarse = planform(&control_valleys(64, 1.5625, 0.0, None), 1.5625, None);
        assert!(coarse.sinuosity[0].is_nan() && coarse.sinuosity[1].is_nan());
    }
}
