//! ADR Finding 162-P -- **the predictability instruments** (the author, F161: « trop prévisibles … un schéma se
//! répéter »). A regular valley spacing is physical (Perron 2009); what is measured is the DISPERSION:
//!
//! 1. `cv_lambda` -- the CV of the gaps between successive valley minima, along the rows and along the columns
//!    apart, averaged;
//! 2. `cv_len` -- the mean over Strahler orders 1–3 of the CV of the streams' lengths;
//! 3. `cv_rb`, `cv_rl` -- the CVs across basins (maximum order ≥ 3) of Horton's bifurcation and length ratios;
//! 4. `theta_std` -- the standard deviation of the confluence angles (each branch read 3 steps up);
//! 5. `peakiness` -- the residual standard deviation (dex) of the ¼-octave radial spectrum about a power law,
//!    between 2 cells and 12.5 km;
//! 6. `r2_gully` -- the mean over 5 km windows (≥ 5 gully heads) of the gullies' axial concentration |⟨e^(2iθ)⟩|.
//!
//! The network is `super::hydro` (D8 on the filled surface, REAL km², the threshold max(1 km², 4 cells)); the cells
//! measured are the land minus `excl` (the D5 cells) minus a 2-cell border. Declared in `f162_declared.md` (P).

use super::amplify::ocean_mask;
use super::hydro::hydrology;
use super::measure::radial_power;
use crate::grid::GridF32;
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

/// One reading of the six instruments (and the means they come with).
#[derive(Clone, Debug, Default)]
pub struct Predict {
    pub cv_lambda: f32,
    pub n_gaps: usize,
    pub cv_len: f32,
    pub cv_len_by_order: [f32; 3],
    pub rb_mean: f32,
    pub cv_rb: f32,
    pub rl_mean: f32,
    pub cv_rl: f32,
    pub n_basins: usize,
    pub theta_mean: f32,
    pub theta_std: f32,
    pub n_junctions: usize,
    pub peakiness: f32,
    pub r2_gully: f32,
    pub n_windows: usize,
}

impl Predict {
    /// The instruments as (name, value, `true` when a SMALLER value is more dispersed / less regular).
    pub fn values(&self) -> [(&'static str, f32, bool); 7] {
        [
            ("CV_λ", self.cv_lambda, false),
            ("CV_L", self.cv_len, false),
            ("CV_Rb", self.cv_rb, false),
            ("CV_Rl", self.cv_rl, false),
            ("σ_θ", self.theta_std, false),
            ("P_s", self.peakiness, true),
            ("R2_g", self.r2_gully, true),
        ]
    }
}

fn mean_std(v: &[f64]) -> (f64, f64) {
    if v.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    let s = (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / v.len() as f64).sqrt();
    (m, s)
}
fn cv(v: &[f64]) -> f32 {
    let (m, s) = mean_std(v);
    if m.abs() < 1e-12 { f32::NAN } else { (s / m) as f32 }
}

pub fn predictability(z: &GridF32, cell_km: f32, excl: Option<&[bool]>) -> Predict {
    let n = z.width;
    let ocean = ocean_mask(z);
    let tex: Vec<bool> = (0..n * n)
        .map(|k| {
            let (x, y) = (k % n, k / n);
            !ocean[k] && !excl.is_some_and(|m| m[k]) && x >= 2 && y >= 2 && x + 2 < n && y + 2 < n
        })
        .collect();
    let zm: Vec<f32> = z.data.iter().zip(&ocean).map(|(v, o)| if *o { 0.0 } else { *v }).collect();
    let mut p = Predict::default();
    // 1. the valley-spacing gaps, rows and columns apart (F162-P1, amended at its control: mixing the two directions
    //    reads a perfectly regular fishbone -- trunks every 64 cells across, gullies every 8 along -- as bimodal)
    let mut gaps_dir: [Vec<f64>; 2] = [Vec::new(), Vec::new()];
    for li in 0..2 * n {
        let gaps = &mut gaps_dir[li / n];
        let idx = |i: usize| if li < n { li * n + i } else { i * n + (li - n) };
        let mut i = 0;
        while i < n {
            if !tex[idx(i)] {
                i += 1;
                continue;
            }
            let s = i;
            while i < n && tex[idx(i)] {
                i += 1;
            }
            let run: Vec<f32> = (s..i).map(|t| zm[idx(t)]).collect();
            if run.len() < 16 {
                continue;
            }
            let mut last: Option<usize> = None;
            for t in 1..run.len() - 1 {
                if !(run[t] < run[t - 1] && run[t] <= run[t + 1]) {
                    continue;
                }
                let mut lm = run[t];
                for u in (0..t).rev() {
                    if run[u] < run[t] {
                        break;
                    }
                    lm = lm.max(run[u]);
                }
                let mut rm = run[t];
                for &v in &run[t + 1..] {
                    if v < run[t] {
                        break;
                    }
                    rm = rm.max(v);
                }
                if lm.min(rm) - run[t] >= 1.0 {
                    if let Some(l) = last {
                        gaps.push((t - l) as f64 * cell_km as f64);
                    }
                    last = Some(t);
                }
            }
        }
    }
    let cvs: Vec<f32> = gaps_dir.iter().map(|g| cv(g)).filter(|v| v.is_finite()).collect();
    p.cv_lambda = if cvs.is_empty() { f32::NAN } else { cvs.iter().sum::<f32>() / cvs.len() as f32 };
    p.n_gaps = gaps_dir[0].len() + gaps_dir[1].len();
    // the network
    let hy = hydrology(z, cell_km);
    let chan = |k: usize| hy.strahler[k] > 0 && !ocean[k];
    let recv = |k: usize| -> Option<usize> {
        let d = hy.dir[k];
        if d == DIR_NONE {
            return None;
        }
        let (x, y) = ((k % n) as i32 + D8_DX[d as usize], (k / n) as i32 + D8_DY[d as usize]);
        if x < 0 || y < 0 || x >= n as i32 || y >= n as i32 {
            return None;
        }
        let r = y as usize * n + x as usize;
        chan(r).then_some(r)
    };
    let step = |a: usize, b: usize| -> f64 {
        let diag = (a % n != b % n) && (a / n != b / n);
        cell_km as f64 * if diag { std::f64::consts::SQRT_2 } else { 1.0 }
    };
    let mut donors: Vec<Vec<usize>> = vec![Vec::new(); n * n];
    let mut has_same = vec![false; n * n];
    for k in 0..n * n {
        if !chan(k) {
            continue;
        }
        if let Some(r) = recv(k) {
            donors[r].push(k);
            if hy.strahler[r] == hy.strahler[k] {
                has_same[r] = true;
            }
        }
    }
    // the outlets (path-compressed)
    let mut outlet = vec![usize::MAX; n * n];
    for k0 in 0..n * n {
        if !chan(k0) || outlet[k0] != usize::MAX {
            continue;
        }
        let mut path = vec![k0];
        let mut cur = k0;
        let o = loop {
            match recv(cur) {
                Some(r) if outlet[r] != usize::MAX => break outlet[r],
                Some(r) => {
                    path.push(r);
                    cur = r;
                }
                None => break cur,
            }
        };
        for c in path {
            outlet[c] = o;
        }
    }
    // 2. the streams
    struct Stream {
        order: u8,
        len: f64,
        head: usize,
        end: usize,
    }
    let mut streams: Vec<Stream> = Vec::new();
    for k in 0..n * n {
        if !chan(k) || has_same[k] {
            continue;
        }
        let w = hy.strahler[k];
        let (mut cur, mut len) = (k, 0.0f64);
        while let Some(r) = recv(cur) {
            if hy.strahler[r] != w {
                len += step(cur, r);
                break;
            }
            len += step(cur, r);
            cur = r;
        }
        streams.push(Stream { order: w, len, head: k, end: cur });
    }
    for o in 1..=3u8 {
        let v: Vec<f64> = streams.iter().filter(|s| s.order == o && tex[s.head]).map(|s| s.len).collect();
        p.cv_len_by_order[o as usize - 1] = cv(&v);
    }
    let present: Vec<f32> = p.cv_len_by_order.iter().copied().filter(|v| v.is_finite()).collect();
    p.cv_len = if present.is_empty() { f32::NAN } else { present.iter().sum::<f32>() / present.len() as f32 };
    // 3. Horton per basin
    let mut basins: std::collections::BTreeMap<usize, Vec<(u8, f64)>> = std::collections::BTreeMap::new();
    for s in &streams {
        basins.entry(outlet[s.head]).or_default().push((s.order, s.len));
    }
    let (mut rbs, mut rls) = (Vec::new(), Vec::new());
    for v in basins.values() {
        let om = v.iter().map(|s| s.0).max().unwrap_or(0);
        if om < 3 {
            continue;
        }
        let (mut lb, mut ll, mut c) = (0f64, 0f64, 0usize);
        for w in 1..om {
            let (n1, n2) = (v.iter().filter(|s| s.0 == w).count(), v.iter().filter(|s| s.0 == w + 1).count());
            if n1 == 0 || n2 == 0 {
                continue;
            }
            let l1 = v.iter().filter(|s| s.0 == w).map(|s| s.1).sum::<f64>() / n1 as f64;
            let l2 = v.iter().filter(|s| s.0 == w + 1).map(|s| s.1).sum::<f64>() / n2 as f64;
            if l1 <= 0.0 || l2 <= 0.0 {
                continue;
            }
            lb += (n1 as f64 / n2 as f64).ln();
            ll += (l2 / l1).ln();
            c += 1;
        }
        if c > 0 {
            rbs.push((lb / c as f64).exp());
            rls.push((ll / c as f64).exp());
        }
    }
    p.n_basins = rbs.len();
    p.rb_mean = mean_std(&rbs).0 as f32;
    p.rl_mean = mean_std(&rls).0 as f32;
    p.cv_rb = cv(&rbs);
    p.cv_rl = cv(&rls);
    // 4. the confluence angles
    let up3 = |start: usize, j: usize| -> usize {
        let mut cur = start;
        for _ in 0..2 {
            match donors[cur].iter().max_by(|a, b| hy.acc_km2[**a].partial_cmp(&hy.acc_km2[**b]).unwrap().then(a.cmp(b))) {
                Some(&d) => cur = d,
                None => break,
            }
        }
        let _ = j;
        cur
    };
    let mut angles: Vec<f64> = Vec::new();
    for j in 0..n * n {
        if donors[j].len() < 2 || !tex[j] {
            continue;
        }
        let mut d = donors[j].clone();
        d.sort_by(|a, b| hy.acc_km2[*b].partial_cmp(&hy.acc_km2[*a]).unwrap().then(a.cmp(b)));
        let v = |s: usize| -> (f64, f64) {
            let u = up3(s, j);
            ((j % n) as f64 - (u % n) as f64, (j / n) as f64 - (u / n) as f64)
        };
        let (a, b) = (v(d[0]), v(d[1]));
        let (na, nb) = ((a.0 * a.0 + a.1 * a.1).sqrt(), (b.0 * b.0 + b.1 * b.1).sqrt());
        if na == 0.0 || nb == 0.0 {
            continue;
        }
        let c = ((a.0 * b.0 + a.1 * b.1) / (na * nb)).clamp(-1.0, 1.0);
        angles.push(c.acos().to_degrees());
    }
    let (am, asd) = mean_std(&angles);
    p.theta_mean = am as f32;
    p.theta_std = asd as f32;
    p.n_junctions = angles.len();
    // 5. the spectrum's peakiness (¼-octave bins, 2 cells … 12.5 km)
    let land_cells = tex.iter().filter(|&&t| t).count();
    let bins = radial_power(&zm, n, 4);
    let pts: Vec<(f64, f64)> = bins
        .iter()
        .filter(|(lam, pw)| *pw > 0.0 && *lam >= 2.0 && *lam * cell_km as f64 <= 12.5)
        .map(|(lam, pw)| ((lam * cell_km as f64).log10(), pw.log10()))
        .collect();
    let _ = land_cells;
    p.peakiness = if pts.len() >= 4 {
        let m = pts.len() as f64;
        let (sx, sy) = pts.iter().fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        let (mx, my) = (sx / m, sy / m);
        let sxx = pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum::<f64>();
        let sxy = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>();
        let b = sxy / sxx.max(1e-12);
        let res = pts.iter().map(|p| { let r = p.1 - (my + b * (p.0 - mx)); r * r }).sum::<f64>() / m;
        res.sqrt() as f32
    } else {
        f32::NAN
    };
    // 6. the gullies' parallelism (5 km windows)
    let win = ((5.0 / cell_km).round() as usize).max(4);
    let nw = n.div_ceil(win);
    let mut acc: Vec<(f64, f64, usize)> = vec![(0.0, 0.0, 0); nw * nw];
    for s in streams.iter().filter(|s| s.order == 1 && tex[s.head]) {
        let (dx, dy) = ((s.end % n) as f64 - (s.head % n) as f64, (s.end / n) as f64 - (s.head / n) as f64);
        if dx * dx + dy * dy < 9.0 {
            continue;
        }
        let t = dy.atan2(dx);
        let w = (s.head / n / win) * nw + (s.head % n) / win;
        acc[w].0 += (2.0 * t).cos();
        acc[w].1 += (2.0 * t).sin();
        acc[w].2 += 1;
    }
    let r2s: Vec<f64> = acc.iter().filter(|a| a.2 >= 5).map(|a| (a.0 * a.0 + a.1 * a.1).sqrt() / a.2 as f64).collect();
    p.r2_gully = mean_std(&r2s).0 as f32;
    p.n_windows = r2s.len();
    p
}

/// The « regular » control (F162-P): a periodic fishbone on `n` cells. Trunks every 64 cells run south to a sea row,
/// gullies perpendicular every 8 cells, 30 m deep; a 0.1 m hash noise breaks the ties. Metres, sea ≤ 0.
pub fn control_regular(n: usize, cell_km: f32) -> GridF32 {
    let cell_m = cell_km * 1000.0;
    let mut g = GridF32::new(n, n, 0.0);
    for y in 0..n {
        for x in 0..n {
            let k = y * n + x;
            if y < 2 {
                g.data[k] = -50.0;
                continue;
            }
            let xm = (x % 64) as f32 - 32.0; // the trunk at x ≡ 32 (mod 64)
            let gy = (y % 8) as f32 - 4.0; // a gully at y ≡ 4 (mod 8)
            let mut h = (k as u32).wrapping_mul(0x9E37_79B1);
            h ^= h >> 15;
            let noise = (h % 1000) as f32 / 1000.0 * 0.1;
            g.data[k] = 5.0 + 0.01 * y as f32 * cell_m + 0.05 * xm.abs() * cell_m - 30.0 * (1.0 - gy.abs() / 4.0) + noise;
        }
    }
    g
}

/// The « dispersed » control: a 9-octave fractal (persistence 0.6, 1 500 m) on a tilt toward a sea row.
pub fn control_fractal(n: usize, cell_km: f32, seed: u32) -> GridF32 {
    let src = crate::terrain::noise::SeededNoise::new(seed, 9);
    let mut g = GridF32::new(n, n, 0.0);
    for y in 0..n {
        for x in 0..n {
            let k = y * n + x;
            if y < 2 {
                g.data[k] = -50.0;
                continue;
            }
            let f = src.fbm(x as f64 / 96.0, y as f64 / 96.0, 9, 2.0, 0.6) as f32;
            g.data[k] = 50.0 + 1500.0 * (f + 0.5).max(0.0) + 0.002 * y as f32 * cell_km * 1000.0;
        }
    }
    g
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The instruments run on both controls, and the clearest separations hold: the fishbone's gullies are parallel
    /// (R2 ≈ 1) where the fractal's are not, and its valley gaps are regular (CV small). The full ×1.5 validation of
    /// every instrument is the bench's (`f162_cascade`), reported and declared there.
    #[test]
    fn the_fishbone_reads_regular_and_the_fractal_does_not() {
        let (n, c) = (256, 0.195);
        let reg = predictability(&control_regular(n, c), c, None);
        let dis = predictability(&control_fractal(n, c, 7), c, None);
        assert!(reg.r2_gully > 0.8 && reg.r2_gully > 1.5 * dis.r2_gully, "R2_g {} against {}", reg.r2_gully, dis.r2_gully);
        assert!(reg.cv_lambda < dis.cv_lambda, "CV_λ {} against {}", reg.cv_lambda, dis.cv_lambda);
        assert!(reg.n_junctions > 0 && dis.n_junctions > 0);
    }
}
