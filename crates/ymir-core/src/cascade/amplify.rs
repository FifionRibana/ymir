//! ADR Finding 160 -- **the amplification**: Schott et al. 2024 (« Terrain Amplification using Multi-scale Erosion »,
//! ACM TOG) transposed to Ymir's scales, after one physics level (`super::physics_level`).
//!
//! Each level is A_k = D_k ∘ T_k ∘ E_k ∘ U_k (§3): a bicubic ×2, then blocks of **bounded erosion** (§4.2), **noisy
//! talus** (§4.3) and **deposition** (§4.4), **with no uplift** (§4.2: « we ignore the uplift component responsible for
//! the emergence of new mountain ranges »). A **retargeting** (§5.1) puts the crests and peaks back on the physics level
//! once, at the end. Every number is declared in `docs/reports/relief_method/f160_amplification/f160_declared.md`
//! (D2–D6); the deviations from the article and from its published code are named where they sit.
//!
//! Heights are in METRES here (the physics level's normalised field is converted at the boundary). The ocean is held at
//! 0 m during a level and its bathymetry restored at the end; land never goes below +0.5 m (production's base level).
//! Every update is Jacobi (computed from the previous buffer), so the rayon split cannot change a bit.

use crate::grid::GridF32;
use crate::seed::WorldSeed;
use crate::terrain::noise::SeededNoise;
use rayon::prelude::*;
use std::time::Instant;

/// The eight neighbours, as the code's `next8`.
const D8: [(i32, i32); 8] = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)];

/// The variants of D6.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Variant {
    /// Schott transposed, constant hardness.
    S,
    /// With a fractal-noise hardness ρ (§4.2).
    SRho,
    /// Ymir's fallback: a small uplift per erosion iteration and a 2×2 retargeting on the previous level.
    SU,
}

impl Variant {
    pub fn label(self) -> &'static str {
        match self {
            Variant::S => "S",
            Variant::SRho => "S+ρ",
            Variant::SU => "S+U",
        }
    }
}

/// Iterations per process (D3).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Budget {
    pub erosion: usize,
    pub talus: usize,
    pub deposit: usize,
}

impl Budget {
    /// The moyen budget: the published code's preset per resolution (`PredefinedErosion`); 4 096² and 8 192² are the
    /// declared extrapolation of its trend (PROXY, used for the cost estimate only).
    pub fn moyen(n: usize) -> Self {
        let (e, t, d) = match n {
            0..=256 => (3000, 600, 2000),
            257..=512 => (1500, 1000, 700),
            513..=1024 => (700, 2000, 200),
            1025..=2048 => (400, 6000, 150),
            2049..=4096 => (300, 6000, 100),
            _ => (200, 6000, 100),
        };
        Self { erosion: e, talus: t, deposit: d }
    }
    /// Every count × `f` (faible 0.5, fort 2).
    pub fn scaled(self, f: f32) -> Self {
        let s = |v: usize| ((v as f32 * f).round() as usize).max(1);
        Self { erosion: s(self.erosion), talus: s(self.talus), deposit: s(self.deposit) }
    }
}

/// The declared parameters (D2, D4, D6).
#[derive(Clone, Debug, PartialEq)]
pub struct AmpConfig {
    pub domain_km: f32,
    /// The multiple-direction routing exponent p (§4.1).
    pub p_flow: f32,
    pub n_exp: f32,
    pub m_exp: f32,
    pub s_max: f32,
    /// a_max, in CELLS of the level.
    pub a_max: f32,
    /// d_L at 256² (m), halving per level.
    pub depth_256_m: f32,
    pub kc: f32,
    pub kd: f32,
    pub talus_tan: f32,
    pub talus_lo: f32,
    pub talus_hi: f32,
    /// m_L = talus_eps · cell (m per neighbour count per iteration).
    pub talus_eps: f32,
    pub talus_wavelength_km: f32,
    pub talus_octaves: usize,
    pub land_floor_m: f32,
    pub variant: Variant,
    /// Budget scale (faible 0.5, moyen 1, fort 2).
    pub budget_scale: f32,
    pub rho_mean: f32,
    pub rho_amp: f32,
    pub rho_wavelength_km: f32,
    pub rho_octaves: usize,
    /// S+U: the uplift over an erosion block, as a fraction of d_L.
    pub su_fraction: f32,
    pub retarget_a0: f32,
    pub retarget_iters: usize,
    pub seed: u64,
    /// ADR Finding 160, diagnostic A1 (NOT blind, added after the declared runs): a multiplier on k_L. 1 = the
    /// declaration. The declared transposition measured an amplification that carves 0.02–0.08 % of its bound; this
    /// asks whether the missing lever is the intensity, which neither the budgets nor added iterations can supply.
    pub k_boost: f32,
}

impl AmpConfig {
    /// The F160 declaration.
    pub fn declared(seed: u64, variant: Variant, budget_scale: f32) -> Self {
        Self {
            domain_km: 400.0,
            p_flow: 1.3,
            n_exp: 2.0,
            m_exp: 0.8,
            s_max: 1.0,
            a_max: 250.0,
            depth_256_m: 200.0,
            kc: 0.1,
            kd: 0.1,
            talus_tan: 0.6494,
            talus_lo: 0.8,
            talus_hi: 1.4,
            talus_eps: 0.002,
            talus_wavelength_km: 8.7,
            talus_octaves: 3,
            land_floor_m: 0.5,
            variant,
            budget_scale,
            rho_mean: 0.5,
            rho_amp: 0.4,
            rho_wavelength_km: 25.0,
            rho_octaves: 6,
            su_fraction: 0.1,
            retarget_a0: 2.0,
            retarget_iters: 500,
            seed,
            k_boost: 1.0,
        }
    }
    pub fn cell_km(&self, n: usize) -> f32 {
        self.domain_km / n as f32
    }
    /// d_L (m) at a level of `n` cells.
    pub fn depth_m(&self, n: usize) -> f32 {
        self.depth_256_m * 256.0 / n as f32
    }
    /// The level's budget (moyen × the scale).
    pub fn budget(&self, n: usize) -> Budget {
        Budget::moyen(n).scaled(self.budget_scale)
    }
}

/// One amplification level.
#[derive(Clone, Debug)]
pub struct AmpLevel {
    pub n: usize,
    pub cell_km: f32,
    /// The bicubic of the previous level (m), sea and the D5 lift included.
    pub upscaled: GridF32,
    pub ocean: Vec<bool>,
    /// D5: cells ≤ 0 m lifted to land, and their mean original depth (m).
    pub lifted: usize,
    pub lifted_mean_depth_m: f32,
    /// The working field (m), the ocean at 0 m.
    pub z: GridF32,
    /// The iterative drainage (cells).
    pub area: Vec<f32>,
    /// The erodibility factor: (1 − ρ) normalised to a mean of 1 for S+ρ, 1 otherwise.
    pub kf: Vec<f32>,
    /// The talus's critical slope per cell.
    pub s0: Vec<f32>,
    pub k: f32,
    pub s_ref: f32,
    pub depth_m: f32,
    /// S+U: the uplift per erosion iteration (m), `None` otherwise.
    pub uplift: Option<Vec<f32>>,
    /// Σ of each process's change (m), 0 on the ocean.
    pub sum_uplift: Vec<f32>,
    pub sum_erosion: Vec<f32>,
    pub sum_talus: Vec<f32>,
    pub sum_deposit: Vec<f32>,
    /// S+U's 2×2 retargeting change (m).
    pub sum_recalage: Vec<f32>,
    /// Iterations applied so far.
    pub done: Budget,
    /// Wall seconds: upscale, erosion, talus, deposition.
    pub secs: [f64; 4],
    /// Erosion iterations' cell-updates where s ≥ s_max, and all eroding cell-updates.
    pub smax_bind: (u64, u64),
}

impl AmpLevel {
    /// The level's result (m): the working field with the ocean's bathymetry restored.
    pub fn result(&self) -> GridF32 {
        let mut g = self.z.clone();
        for (k, v) in g.data.iter_mut().enumerate() {
            if self.ocean[k] {
                *v = self.upscaled.data[k];
            }
        }
        g
    }
    /// The field after a sub-step (m): 0 upscaled, 1 + erosion, 2 + talus, 3 + deposition (= the result).
    pub fn substep_field(&self, s: usize) -> GridF32 {
        if s >= 3 {
            return self.result();
        }
        let mut g = self.upscaled.clone();
        for k in 0..g.data.len() {
            if self.ocean[k] {
                continue;
            }
            let mut v = self.upscaled.data[k] + self.sum_uplift[k];
            if s >= 1 {
                v += self.sum_erosion[k];
            }
            if s >= 2 {
                v += self.sum_talus[k];
            }
            g.data[k] = v;
        }
        g
    }
    /// Retained bytes (the memory estimate): 4 grids + 6 sums + the factors + the mask.
    pub fn retained_bytes(&self) -> usize {
        let c = self.n * self.n;
        c * 4 * 13 + c
    }
}

/// The ocean: cells ≤ 0 m 4-connected to the map border. The rest of the land, below 0 m or not, is land (D5).
pub fn ocean_mask(z: &GridF32) -> Vec<bool> {
    let (w, h) = (z.width, z.height);
    let mut ocean = vec![false; w * h];
    let mut stack = Vec::new();
    for x in 0..w {
        for y in [0, h - 1] {
            let k = y * w + x;
            if z.data[k] <= 0.0 && !ocean[k] {
                ocean[k] = true;
                stack.push(k);
            }
        }
    }
    for y in 0..h {
        for x in [0, w - 1] {
            let k = y * w + x;
            if z.data[k] <= 0.0 && !ocean[k] {
                ocean[k] = true;
                stack.push(k);
            }
        }
    }
    while let Some(k) = stack.pop() {
        let (x, y) = (k % w, k / w);
        let nb = [
            (x > 0).then(|| k - 1),
            (x + 1 < w).then(|| k + 1),
            (y > 0).then(|| k - w),
            (y + 1 < h).then(|| k + w),
        ];
        for j in nb.into_iter().flatten() {
            if !ocean[j] && z.data[j] <= 0.0 {
                ocean[j] = true;
                stack.push(j);
            }
        }
    }
    ocean
}

/// D5: lift the land cells ≤ 0 m to `to_m`; (count, mean original depth).
pub fn lift_interior(z: &mut GridF32, ocean: &[bool], to_m: f32) -> (usize, f32) {
    let (mut c, mut s) = (0usize, 0f64);
    for (k, v) in z.data.iter_mut().enumerate() {
        if !ocean[k] && *v <= to_m {
            if *v <= 0.0 {
                c += 1;
                s += -*v as f64;
            }
            *v = to_m;
        }
    }
    (c, if c > 0 { (s / c as f64) as f32 } else { 0.0 })
}

/// A fractal noise in world coordinates (km), rescaled to [0, 1] over the grid (min–max).
fn noise01(n: usize, cell_km: f32, wavelength_km: f32, octaves: usize, seed: u64, phase: &str) -> Vec<f32> {
    let s = WorldSeed::new(seed).derive_seed(phase) as u32;
    let src = SeededNoise::new(s, octaves);
    let mut v: Vec<f32> = (0..n * n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = ((k % n) as f64 * cell_km as f64, (k / n) as f64 * cell_km as f64);
            src.fbm(x / wavelength_km as f64, y / wavelength_km as f64, octaves, 2.0, 0.5) as f32
        })
        .collect();
    let (lo, hi) = v.iter().fold((f32::MAX, f32::MIN), |(a, b), &x| (a.min(x), b.max(x)));
    let span = (hi - lo).max(1e-9);
    for x in &mut v {
        *x = (*x - lo) / span;
    }
    v
}

#[inline]
fn dist(d: usize, cell_m: f32) -> f32 {
    if D8[d].0 != 0 && D8[d].1 != 0 { cell_m * std::f32::consts::SQRT_2 } else { cell_m }
}

/// The sum of the downhill weights s^p of every cell (the denominator of eq. 1); 0 on a pit or the ocean.
fn norms(z: &[f32], ocean: &[bool], n: usize, cell_m: f32, p: f32) -> Vec<f32> {
    (0..n * n)
        .into_par_iter()
        .map(|k| {
            if ocean[k] {
                return 0.0;
            }
            let (x, y) = ((k % n) as i32, (k / n) as i32);
            let mut s = 0f32;
            for d in 0..8 {
                let (nx, ny) = (x + D8[d].0, y + D8[d].1);
                if nx < 0 || ny < 0 || nx >= n as i32 || ny >= n as i32 {
                    continue;
                }
                let j = ny as usize * n + nx as usize;
                let sl = (z[k] - z[j]) / dist(d, cell_m);
                if sl > 0.0 {
                    s += sl.powf(p);
                }
            }
            s
        })
        .collect()
}

/// Σ over the uphill land neighbours q of p of w(q, p) · src(q) (eqs. 1–2's sum).
#[inline]
#[allow(clippy::too_many_arguments)]
fn gather(z: &[f32], ocean: &[bool], norm: &[f32], src: &[f32], k: usize, n: usize, cell_m: f32, p: f32) -> f32 {
    let (x, y) = ((k % n) as i32, (k / n) as i32);
    let mut acc = 0f32;
    for d in 0..8 {
        let (nx, ny) = (x + D8[d].0, y + D8[d].1);
        if nx < 0 || ny < 0 || nx >= n as i32 || ny >= n as i32 {
            continue;
        }
        let q = ny as usize * n + nx as usize;
        if ocean[q] || norm[q] <= 0.0 {
            continue;
        }
        let sl = (z[q] - z[k]) / dist(d, cell_m);
        if sl > 0.0 {
            acc += sl.powf(p) / norm[q] * src[q];
        }
    }
    acc
}

/// The steepest downhill neighbour of k: (slope, receiver index), `None` on a pit.
#[inline]
fn steepest(z: &[f32], k: usize, n: usize, cell_m: f32) -> Option<(f32, usize)> {
    let (x, y) = ((k % n) as i32, (k / n) as i32);
    let mut best: Option<(f32, usize)> = None;
    for d in 0..8 {
        let (nx, ny) = (x + D8[d].0, y + D8[d].1);
        if nx < 0 || ny < 0 || nx >= n as i32 || ny >= n as i32 {
            continue;
        }
        let j = ny as usize * n + nx as usize;
        let sl = (z[k] - z[j]) / dist(d, cell_m);
        if sl > 0.0 && best.is_none_or(|b| sl > b.0) {
            best = Some((sl, j));
        }
    }
    best
}

/// One exact multiple-direction accumulation (cells sorted by height, the weights of eq. 1): the level's starting
/// drainage (D2's one deviation) and the retargeting's ridge detection.
pub fn exact_area(z: &[f32], ocean: &[bool], n: usize, cell_m: f32, p: f32) -> Vec<f32> {
    let norm = norms(z, ocean, n, cell_m, p);
    let mut order: Vec<usize> = (0..n * n).filter(|&k| !ocean[k]).collect();
    order.sort_unstable_by(|&a, &b| z[b].partial_cmp(&z[a]).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(&b)));
    let mut a = vec![1f32; n * n];
    for &k in &order {
        if norm[k] <= 0.0 {
            continue;
        }
        let (x, y) = ((k % n) as i32, (k / n) as i32);
        for d in 0..8 {
            let (nx, ny) = (x + D8[d].0, y + D8[d].1);
            if nx < 0 || ny < 0 || nx >= n as i32 || ny >= n as i32 {
                continue;
            }
            let j = ny as usize * n + nx as usize;
            let sl = (z[k] - z[j]) / dist(d, cell_m);
            if sl > 0.0 {
                a[j] += sl.powf(p) / norm[k] * a[k];
            }
        }
    }
    a
}

/// p90 of the land slope (steepest descent, m/m).
fn land_slope_p90(z: &[f32], ocean: &[bool], n: usize, cell_m: f32) -> f32 {
    let mut s: Vec<f32> =
        (0..n * n).filter(|&k| !ocean[k]).map(|k| steepest(z, k, n, cell_m).map_or(0.0, |b| b.0)).collect();
    if s.is_empty() {
        return 0.0;
    }
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    s[((s.len() - 1) as f32 * 0.9) as usize]
}

impl AmpLevel {
    /// U_k: the bicubic ×2 of `prev` (m), the ocean, the D5 lift, the level's k, the talus slopes and the hardness.
    /// `phys_up` is the physics level at this resolution (m), used by S+U's uplift.
    pub fn upscale(prev: &GridF32, cfg: &AmpConfig, phys_up: Option<&GridF32>, peak_m: f32) -> Self {
        let t = Instant::now();
        let mut up = super::upsample2(prev);
        let n = up.width;
        let ocean = ocean_mask(&up);
        let (lifted, lifted_mean_depth_m) = lift_interior(&mut up, &ocean, cfg.land_floor_m);
        let mut z = up.clone();
        for k in 0..n * n {
            if ocean[k] {
                z.data[k] = 0.0;
            }
        }
        let cell_km = cfg.cell_km(n);
        let cell_m = cell_km * 1000.0;
        let area = exact_area(&z.data, &ocean, n, cell_m, cfg.p_flow);
        let s_ref = land_slope_p90(&z.data, &ocean, n, cell_m).max(1e-4);
        let depth_m = cfg.depth_m(n);
        let n_e = Budget::moyen(n).erosion as f32;
        let k = cfg.k_boost * depth_m / (n_e * s_ref.min(cfg.s_max).powf(cfg.n_exp) * cfg.a_max.powf(cfg.m_exp));
        let kf = match cfg.variant {
            Variant::SRho => {
                let f = noise01(n, cell_km, cfg.rho_wavelength_km, cfg.rho_octaves, cfg.seed, "cascade_hardness");
                // ρ = mean + amp·(2f − 1) ∈ [0.1, 0.9]; k(1 − ρ) / (1 − mean) keeps S's mean erodibility
                f.iter().map(|&v| (1.0 - (cfg.rho_mean + cfg.rho_amp * (2.0 * v - 1.0))) / (1.0 - cfg.rho_mean)).collect()
            }
            _ => vec![1.0; n * n],
        };
        let nu = noise01(n, cell_km, cfg.talus_wavelength_km, cfg.talus_octaves, cfg.seed, "cascade_talus");
        let s0 = nu.iter().map(|&v| cfg.talus_tan * (cfg.talus_lo + (cfg.talus_hi - cfg.talus_lo) * v)).collect();
        let uplift = (cfg.variant == Variant::SU).then(|| {
            let per_iter = cfg.su_fraction * depth_m / n_e;
            (0..n * n)
                .map(|k| {
                    if ocean[k] {
                        0.0
                    } else {
                        per_iter * phys_up.map_or(1.0, |p| (p.data[k].max(0.0) / peak_m.max(1.0)).min(1.0))
                    }
                })
                .collect()
        });
        let zeros = vec![0f32; n * n];
        Self {
            n,
            cell_km,
            upscaled: up,
            ocean,
            lifted,
            lifted_mean_depth_m,
            z,
            area,
            kf,
            s0,
            k,
            s_ref,
            depth_m,
            uplift,
            sum_uplift: zeros.clone(),
            sum_erosion: zeros.clone(),
            sum_talus: zeros.clone(),
            sum_deposit: zeros.clone(),
            sum_recalage: zeros,
            done: Budget::default(),
            secs: [t.elapsed().as_secs_f64(), 0.0, 0.0, 0.0],
            smax_bind: (0, 0),
        }
    }

    /// E_k (§4.2): `iters` iterations of the bounded explicit stream power, each after one routing iteration (§4.1).
    pub fn erode(&mut self, iters: usize, cfg: &AmpConfig, cancel: &dyn Fn() -> bool) {
        let t = Instant::now();
        let n = self.n;
        let cell_m = self.cell_km * 1000.0;
        let (p, smax_n, amax_m) = (cfg.p_flow, cfg.s_max.powf(cfg.n_exp), cfg.a_max.powf(cfg.m_exp));
        for _ in 0..iters {
            if cancel() {
                break;
            }
            if let Some(u) = &self.uplift {
                for k in 0..n * n {
                    self.z.data[k] += u[k];
                    self.sum_uplift[k] += u[k];
                }
            }
            let z = &self.z.data;
            let norm = norms(z, &self.ocean, n, cell_m, p);
            let area: Vec<f32> = (0..n * n)
                .into_par_iter()
                .map(|k| 1.0 + gather(z, &self.ocean, &norm, &self.area, k, n, cell_m, p))
                .collect();
            let out: Vec<(f32, u8)> = (0..n * n)
                .into_par_iter()
                .map(|k| {
                    if self.ocean[k] {
                        return (z[k], 0);
                    }
                    let Some((s, r)) = steepest(z, k, n, cell_m) else { return (z[k], 0) };
                    let sn = s.powf(cfg.n_exp);
                    let e = sn.min(smax_n) * area[k].powf(cfg.m_exp).min(amax_m);
                    let h = (z[k] - self.k * self.kf[k] * e).max(z[r]).max(cfg.land_floor_m.min(z[k]));
                    (h, if sn >= smax_n { 2 } else { 1 })
                })
                .collect();
            for k in 0..n * n {
                let (h, flag) = out[k];
                self.sum_erosion[k] += h - self.z.data[k];
                self.z.data[k] = h;
                if flag > 0 {
                    self.smax_bind.1 += 1;
                    if flag == 2 {
                        self.smax_bind.0 += 1;
                    }
                }
            }
            self.area = area;
            self.done.erosion += 1;
        }
        self.secs[1] += t.elapsed().as_secs_f64();
    }

    /// T_k (§4.3, the code's form): h ← h + m · (α − β) against the noisy critical slope s₀(p).
    pub fn talus(&mut self, iters: usize, cfg: &AmpConfig, cancel: &dyn Fn() -> bool) {
        let t = Instant::now();
        let n = self.n;
        let cell_m = self.cell_km * 1000.0;
        let m = cfg.talus_eps * cell_m;
        for _ in 0..iters {
            if cancel() {
                break;
            }
            let z = &self.z.data;
            let out: Vec<f32> = (0..n * n)
                .into_par_iter()
                .map(|k| {
                    if self.ocean[k] {
                        return z[k];
                    }
                    let (x, y) = ((k % n) as i32, (k / n) as i32);
                    let (mut a, mut b) = (0i32, 0i32);
                    for d in 0..8 {
                        let (nx, ny) = (x + D8[d].0, y + D8[d].1);
                        if nx < 0 || ny < 0 || nx >= n as i32 || ny >= n as i32 {
                            continue;
                        }
                        let j = ny as usize * n + nx as usize;
                        let dd = dist(d, cell_m);
                        if (z[j] - z[k]) / dd > self.s0[k] {
                            a += 1;
                        }
                        if (z[k] - z[j]) / dd > self.s0[k] {
                            b += 1;
                        }
                    }
                    (z[k] + m * (a - b) as f32).max(cfg.land_floor_m.min(z[k]))
                })
                .collect();
            for k in 0..n * n {
                self.sum_talus[k] += out[k] - self.z.data[k];
            }
            self.z.data = out;
            self.done.talus += 1;
        }
        self.secs[2] += t.elapsed().as_secs_f64();
    }

    /// D_k (§4.4, the paper's form, in metres): sediment created from the (bounded) stream power, carried by the same
    /// weights, deposited where it exceeds the stream power. One routing iteration per deposition iteration.
    pub fn deposit(&mut self, iters: usize, cfg: &AmpConfig, cancel: &dyn Fn() -> bool) {
        let t = Instant::now();
        let n = self.n;
        let cell_m = self.cell_km * 1000.0;
        let (p, smax_n, amax_m) = (cfg.p_flow, cfg.s_max.powf(cfg.n_exp), cfg.a_max.powf(cfg.m_exp));
        let mut g = vec![0f32; n * n];
        for _ in 0..iters {
            if cancel() {
                break;
            }
            let z = &self.z.data;
            let norm = norms(z, &self.ocean, n, cell_m, p);
            let out: Vec<(f32, f32, f32)> = (0..n * n)
                .into_par_iter()
                .map(|k| {
                    let a = 1.0 + gather(z, &self.ocean, &norm, &self.area, k, n, cell_m, p);
                    if self.ocean[k] {
                        return (z[k], 0.0, a);
                    }
                    let e_m = steepest(z, k, n, cell_m).map_or(0.0, |(s, _)| {
                        self.k * self.kf[k] * s.powf(cfg.n_exp).min(smax_n) * a.powf(cfg.m_exp).min(amax_m)
                    });
                    let tr = gather(z, &self.ocean, &norm, &g, k, n, cell_m, p);
                    let phi = tr - e_m;
                    let d = if phi > 0.0 { tr.min(cfg.kd * phi) } else { 0.0 };
                    (z[k] + d, cfg.kc * e_m + tr - d, a)
                })
                .collect();
            for k in 0..n * n {
                let (h, gg, a) = out[k];
                self.sum_deposit[k] += h - self.z.data[k];
                self.z.data[k] = h;
                g[k] = gg;
                self.area[k] = a;
            }
            self.done.deposit += 1;
        }
        self.secs[3] += t.elapsed().as_secs_f64();
    }

    /// S+U's 2×2 retargeting: z ← z + U₂(prev − R(z)) on land, so the level restricted equals the previous level.
    pub fn recalage_2x2(&mut self, prev: &GridF32, cfg: &AmpConfig) {
        let r = restrict(&self.result());
        let mut d = prev.clone();
        for (v, b) in d.data.iter_mut().zip(&r.data) {
            *v -= b;
        }
        let du = super::upsample2(&d);
        for k in 0..self.n * self.n {
            if self.ocean[k] {
                continue;
            }
            let h = (self.z.data[k] + du.data[k]).max(cfg.land_floor_m);
            self.sum_recalage[k] += h - self.z.data[k];
            self.z.data[k] = h;
        }
    }
}

/// The centred full-weighting restriction (1/4, 1/2, 1/4) of a 2N field to N, periodic.
pub fn restrict(g: &GridF32) -> GridF32 {
    let (w2, n) = (g.width, g.width / 2);
    let wt = [0.25f32, 0.5, 0.25];
    let mut out = GridF32::new(n, n, 0.0);
    for j in 0..n {
        for i in 0..n {
            let mut s = 0f32;
            for (a, wy) in wt.iter().enumerate() {
                let y = (2 * j + w2 + a - 1) % w2;
                for (b, wx) in wt.iter().enumerate() {
                    let x = (2 * i + w2 + b - 1) % w2;
                    s += wy * wx * g.data[y * w2 + x];
                }
            }
            out.data[j * n + i] = s;
        }
    }
    out
}

/// §5.1: the error δ = h₀ − h_A at the ridge points (exact drainage < a₀ cells), diffused over `iters` Jacobi
/// iterations with the constraints held, then added to the land. Returns (the retargeted field, the constraint count).
pub fn retarget(h_a: &GridF32, h0: &GridF32, ocean: &[bool], cfg: &AmpConfig) -> (GridF32, usize) {
    let n = h_a.width;
    let cell_m = cfg.cell_km(n) * 1000.0;
    let mut za = h_a.data.clone();
    for k in 0..n * n {
        if ocean[k] {
            za[k] = 0.0;
        }
    }
    let a = exact_area(&za, ocean, n, cell_m, cfg.p_flow);
    let cons: Vec<bool> = (0..n * n).map(|k| !ocean[k] && a[k] < cfg.retarget_a0).collect();
    let e0: Vec<f32> = (0..n * n).map(|k| if cons[k] { h0.data[k] - h_a.data[k] } else { 0.0 }).collect();
    let mut e = e0.clone();
    for _ in 0..cfg.retarget_iters {
        e = (0..n * n)
            .into_par_iter()
            .map(|k| {
                if cons[k] {
                    return e0[k];
                }
                let (x, y) = (k % n, k / n);
                let (xm, xp, ym, yp) = ((x + n - 1) % n, (x + 1) % n, (y + n - 1) % n, (y + 1) % n);
                0.25 * (e[y * n + xm] + e[y * n + xp] + e[ym * n + x] + e[yp * n + x])
            })
            .collect();
    }
    let mut out = h_a.clone();
    for k in 0..n * n {
        if !ocean[k] {
            out.data[k] = (h_a.data[k] + e[k]).max(cfg.land_floor_m);
        }
    }
    (out, cons.iter().filter(|&&c| c).count())
}

/// A whole chain: the physics level (m) then the amplification levels.
#[derive(Clone, Debug)]
pub struct Chain {
    pub cfg: AmpConfig,
    /// The physics level's result (m), the D5 lift included.
    pub physics_m: GridF32,
    pub peak_m: f32,
    pub levels: Vec<AmpLevel>,
}

impl Chain {
    pub fn new(physics_m: GridF32, cfg: AmpConfig) -> Self {
        let ocean = ocean_mask(&physics_m);
        let peak_m = physics_m.data.iter().zip(&ocean).filter(|(_, o)| !**o).map(|(v, _)| *v).fold(f32::MIN, f32::max);
        Self { cfg, physics_m, peak_m, levels: Vec::new() }
    }
    /// The last level's result (m), or the physics level's.
    pub fn last_result(&self) -> GridF32 {
        self.levels.last().map_or_else(|| self.physics_m.clone(), |l| l.result())
    }
    /// The physics level bicubic-upscaled to `n` cells (m).
    pub fn physics_at(&self, n: usize) -> GridF32 {
        super::upsample_to(&self.physics_m, n)
    }
    /// U_k only: start the next level (the viz runs the processes one block at a time).
    pub fn start_level(&mut self) -> &mut AmpLevel {
        let prev = self.last_result();
        let n = prev.width * 2;
        let phys = (self.cfg.variant == Variant::SU).then(|| self.physics_at(n));
        let lvl = AmpLevel::upscale(&prev, &self.cfg, phys.as_ref(), self.peak_m);
        self.levels.push(lvl);
        self.levels.last_mut().unwrap()
    }
    /// A whole level at its declared budget (D3): U, E, T, D, then S+U's 2×2 retargeting.
    pub fn next_level(&mut self, cancel: &dyn Fn() -> bool) -> Option<&AmpLevel> {
        let prev = self.last_result();
        let cfg = self.cfg.clone();
        let lvl = self.start_level();
        let b = cfg.budget(lvl.n);
        lvl.erode(b.erosion, &cfg, cancel);
        lvl.talus(b.talus, &cfg, cancel);
        lvl.deposit(b.deposit, &cfg, cancel);
        if cfg.variant == Variant::SU {
            lvl.recalage_2x2(&prev, &cfg);
        }
        if cancel() {
            self.levels.pop();
            return None;
        }
        self.levels.last()
    }
    /// The final retargeting of the last level on the physics level (§5.1, D4).
    pub fn retargeted(&self) -> Option<(GridF32, usize)> {
        let l = self.levels.last()?;
        let h0 = self.physics_at(l.n);
        Some(retarget(&l.result(), &h0, &l.ocean, &self.cfg))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn island(n: usize) -> GridF32 {
        // a 2 km dome on an ocean at −500 m, with a closed interior basin below 0 m
        let mut g = GridF32::new(n, n, -500.0);
        let c = n as f32 / 2.0;
        for y in 0..n {
            for x in 0..n {
                let r2 = ((x as f32 - c).powi(2) + (y as f32 - c).powi(2)) / (n as f32 * 0.3).powi(2);
                let v = 2500.0 * (-r2).exp() - 500.0;
                let b = ((x as f32 - c * 1.3).powi(2) + (y as f32 - c).powi(2)) / (n as f32 * 0.04).powi(2);
                g.data[y * n + x] = v - 2200.0 * (-b).exp();
            }
        }
        g
    }

    /// ADR Finding 160 -- the amplification's invariants on a synthetic island: the erosion only lowers and never below
    /// its receiver or the land floor, the deposition only raises, the carving stays under the article's bound, the
    /// ocean is restored, the interior basin below 0 m is land (D5), the decomposition adds up, and the retargeting
    /// moves the ridge points onto the reference. Negative control: with k = 0 nothing is carved.
    #[test]
    fn the_amplification_is_bounded_and_adds_up() {
        let mut cfg = AmpConfig::declared(7, Variant::S, 0.05);
        cfg.domain_km = 100.0;
        let phys = island(64);
        let mut ch = Chain::new(phys.clone(), cfg.clone());
        let l = ch.next_level(&|| false).unwrap().clone();
        assert_eq!(l.n, 128);
        assert!(l.lifted > 0, "the interior basin below 0 m must be land (D5)");
        let b = cfg.budget(128);
        let bound = b.erosion as f32 * l.k * cfg.s_max.powf(cfg.n_exp) * cfg.a_max.powf(cfg.m_exp);
        let r = l.result();
        for k in 0..128 * 128 {
            if l.ocean[k] {
                assert_eq!(r.data[k], l.upscaled.data[k], "the ocean is restored");
                continue;
            }
            assert!(l.sum_erosion[k] <= 0.0 && -l.sum_erosion[k] <= bound * 1.0001, "erosion {} vs bound {bound}", l.sum_erosion[k]);
            assert!(l.sum_deposit[k] >= 0.0);
            assert!(r.data[k] >= cfg.land_floor_m - 1e-3, "land under the floor at {k}");
            let sum = l.upscaled.data[k] + l.sum_uplift[k] + l.sum_erosion[k] + l.sum_talus[k] + l.sum_deposit[k];
            assert!((sum - r.data[k]).abs() < 1e-2, "the decomposition does not add up at {k}");
        }
        assert!(l.sum_erosion.iter().any(|&v| v < -0.01), "the erosion must carve something");
        // the retargeting puts the ridge points on the reference
        let h0 = ch.physics_at(128);
        let (rt, nc) = ch.retargeted().unwrap();
        assert!(nc > 0);
        let a = exact_area(
            &l.result().data.iter().zip(&l.ocean).map(|(v, o)| if *o { 0.0 } else { *v }).collect::<Vec<_>>(),
            &l.ocean,
            128,
            cfg.cell_km(128) * 1000.0,
            cfg.p_flow,
        );
        for k in 0..128 * 128 {
            if !l.ocean[k] && a[k] < cfg.retarget_a0 && h0.data[k] > cfg.land_floor_m {
                assert!((rt.data[k] - h0.data[k]).abs() < 1e-2, "a ridge point is not retargeted at {k}");
            }
        }
        // negative control: k = 0 carves nothing
        let mut l0 = AmpLevel::upscale(&phys, &cfg, None, 2000.0);
        l0.k = 0.0;
        l0.erode(5, &cfg, &|| false);
        assert!(l0.sum_erosion.iter().all(|&v| v == 0.0));
    }

    /// The hardness keeps the mean erodibility (S+ρ against S) and varies.
    #[test]
    fn the_hardness_keeps_the_mean_erodibility() {
        let cfg = AmpConfig::declared(7, Variant::SRho, 1.0);
        let l = AmpLevel::upscale(&island(64), &cfg, None, 2000.0);
        let m = l.kf.iter().sum::<f32>() / l.kf.len() as f32;
        assert!((m - 1.0).abs() < 0.25, "mean factor {m}");
        let (lo, hi) = l.kf.iter().fold((f32::MAX, f32::MIN), |(a, b), &x| (a.min(x), b.max(x)));
        assert!(lo < 0.5 && hi > 1.5, "the hardness must vary: [{lo}, {hi}]");
    }
}
