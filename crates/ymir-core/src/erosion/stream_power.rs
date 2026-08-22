//! Routed stream-power incision — Braun & Willett (2013) implicit "FastScape".
//!
//! PROTOTYPE, off by default (no pipeline wiring). Carves valleys along the
//! drainage network that already exists, deterministically and globally, so the
//! channel hierarchy is imposed by construction rather than hoped for from
//! uncorrelated droplets (see docs/adr/0001, Finding 3).
//!
//! Detachment-limited stream power `E = K · A^m · S^n` solved with the implicit /
//! stack scheme: build a topological stack of the flow tree (receivers before
//! donors), then update each node from base level upward. For `n = 1` the per-node
//! update is a closed form and UNCONDITIONALLY STABLE (no CFL timestep limit); for
//! `n ≠ 1` a few Newton iterations per node. O(n) per drainage↔incision pass.
//!
//! Routing reuses [`crate::terrain::flow::compute_flow`] (depression fill, D8
//! receivers, flow accumulation) — no re-implemented pit filling.

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;
use crate::terrain::flow::{D8_DIST, D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};

/// Stream-power incision tunables (prototype).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct StreamPowerConfig {
    /// Erodibility `K` (lumped with the timestep — see `dt`). Calibrated to a
    /// target channel-incision depth, not to appearance (see the diagnostic).
    pub k: f32,
    /// Drainage-area exponent `m` (standard detachment-limited: 0.5).
    pub m: f32,
    /// Slope exponent `n` (standard: 1.0 — closed-form implicit update).
    pub n: f32,
    /// Timestep per drainage↔incision iteration (lumped into `K·dt`).
    pub dt: f32,
    /// Number of drainage↔incision iterations (recompute flow between each, so the
    /// network can reorganise as the terrain changes — the staleness handling).
    pub iterations: usize,
    /// Sea level (norm). Cells at or below it are fixed base level (no incision).
    pub sea_level: f32,
    /// Linear hillslope diffusion coefficient `D` (0 = off) — the classic
    /// stream-power partner that smooths interfluves and opens V/U valley walls.
    pub diffusion: f32,
    /// Explicit diffusion sub-steps per iteration (stability: keep `diffusion /
    /// diffusion_substeps ≤ 0.2`).
    pub diffusion_substeps: usize,
    /// **Critical drainage area `A_c`** (cells) — the LEM channel-head criterion.
    /// Cells with accumulation `< min_area_cells` are in the HILLSLOPE regime and
    /// receive NO fluvial incision (pure stream power over-carves steep headwaters
    /// because `E = K·A^m·S^n` has no lower bound — those cells are physically
    /// hillslopes, not channels). `0` (default) = fluvial everywhere (legacy).
    pub min_area_cells: f32,
    /// **Incision threshold `θ`** — `E = K·max(0, A^m·S^n − θ)`: no incision below a
    /// critical stream power, so low-energy cells (again, headwaters) do not carve.
    /// `0` (default) = no threshold (legacy).
    pub threshold: f32,
}

impl Default for StreamPowerConfig {
    fn default() -> Self {
        Self {
            k: 1.0,
            m: 0.5,
            n: 1.0,
            dt: 1.0,
            iterations: 4,
            sea_level: 0.5,
            diffusion: 0.0,
            diffusion_substeps: 4,
            min_area_cells: 0.0,
            threshold: 0.0,
        }
    }
}

/// Incise `height` in place-of-a-clone and return the eroded field. Runs
/// `cfg.iterations` drainage↔incision passes.
pub fn incise(height: &GridF32, cfg: &StreamPowerConfig) -> GridF32 {
    incise_with_progress(height, cfg, &mut |_, _| {})
}

/// [`incise`] with a per-iteration callback `(iteration_index, &current_field)` —
/// used by the diagnostic to measure whether the network reorganises between
/// iterations (staleness) and to profile runtime.
pub fn incise_with_progress(
    height: &GridF32,
    cfg: &StreamPowerConfig,
    progress: &mut dyn FnMut(usize, &GridF32),
) -> GridF32 {
    let (w, h) = (height.width, height.height);
    let n = w * h;
    let mut field = height.clone();
    let flow_cfg = FlowConfig { sea_level: cfg.sea_level, ..Default::default() };

    for iter in 0..cfg.iterations {
        // 1. Route on the current surface (depression fill + D8 + accumulation).
        let flow = compute_flow(&field, &flow_cfg);

        // 2. Receiver + distance per cell (base = points off-grid / no outlet /
        //    sub-sea). `receiver[k] == k` marks a fixed base node.
        let mut receiver = vec![0usize; n];
        let mut dist = vec![1.0f32; n];
        for k in 0..n {
            let d = flow.direction[k];
            let (x, y) = (k % w, k / w);
            if d == DIR_NONE || field.data[k] <= cfg.sea_level {
                receiver[k] = k;
                continue;
            }
            let nx = x as i32 + D8_DX[d as usize];
            let ny = y as i32 + D8_DY[d as usize];
            if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                receiver[k] = k; // routes off the map edge → outlet/base
            } else {
                receiver[k] = ny as usize * w + nx as usize;
                dist[k] = D8_DIST[d as usize];
            }
        }

        // 3. Topological stack (receivers before donors) via the donor CSR tree.
        let mut ndonor = vec![0u32; n];
        for k in 0..n {
            if receiver[k] != k {
                ndonor[receiver[k]] += 1;
            }
        }
        let mut offset = vec![0usize; n + 1];
        for k in 0..n {
            offset[k + 1] = offset[k] + ndonor[k] as usize;
        }
        let mut donors = vec![0usize; offset[n]];
        let mut fill = offset.clone();
        for k in 0..n {
            let r = receiver[k];
            if r != k {
                donors[fill[r]] = k;
                fill[r] += 1;
            }
        }
        let mut stack = Vec::with_capacity(n);
        let mut work = Vec::new();
        for k in 0..n {
            if receiver[k] == k {
                work.push(k);
                while let Some(node) = work.pop() {
                    stack.push(node);
                    for d in offset[node]..offset[node + 1] {
                        work.push(donors[d]);
                    }
                }
            }
        }

        // 4. Implicit incision in stack order (base → up). h_r is already updated.
        let kdt = cfg.k * cfg.dt;
        for &k in &stack {
            let r = receiver[k];
            if r == k || field.data[k] <= cfg.sea_level {
                continue; // base level, fixed
            }
            let area = flow.accumulation.data[k].max(1.0);
            if area < cfg.min_area_cells {
                continue; // A_c: hillslope regime — no fluvial incision (channel head)
            }
            let hr = field.data[r];
            let ho = field.data[k];
            let a = area.powf(cfg.m) / dist[k].powf(cfg.n);
            // Threshold gate: skip cells whose stream power A^m·S^n is below θ.
            let sp_now = a * (ho - hr).max(0.0).powf(cfg.n);
            if sp_now <= cfg.threshold {
                continue;
            }
            let f = kdt * a;
            let hn = if (cfg.n - 1.0).abs() < 1e-6 {
                // Closed form, n=1, with the −θ threshold (the +kdt·θ term reduces
                // incision): E = K·max(0, A·S − θ).
                (ho + f * hr + kdt * cfg.threshold) / (1.0 + f)
            } else {
                // Newton for general n with threshold: g(x) = x − ho + kdt·max(0,
                // a·(x−hr)^n − θ) = 0, x ≥ hr.
                let mut x = ho;
                for _ in 0..8 {
                    let s = (x - hr).max(0.0);
                    let e = (a * s.powf(cfg.n) - cfg.threshold).max(0.0);
                    let g = x - ho + kdt * e;
                    let de = if a * s.powf(cfg.n) > cfg.threshold {
                        kdt * a * cfg.n * s.powf(cfg.n - 1.0)
                    } else {
                        0.0
                    };
                    x -= g / (1.0 + de);
                    if x < hr {
                        x = hr;
                    }
                }
                x
            };
            // Never incise below the receiver (no reversal) or ABOVE the old height
            // (the threshold term must not deposit).
            field.data[k] = hn.clamp(hr, ho);
        }

        // 5. Hillslope diffusion (explicit, a few sub-steps), interleaved with the
        // incision each iteration so it holds interfluves WHILE channels incise (the
        // coupled LEM, not a post-pass smooth). REGIME SPLIT: when `min_area_cells`
        // is set, diffusion runs ONLY on hillslope cells (A < A_c) and stream power
        // ONLY on channels (A ≥ A_c) — the standard channel-head partition that stops
        // the fluvial law from over-carving physical hillslopes. `min_area_cells = 0`
        // → diffuse everywhere (legacy).
        if cfg.diffusion > 0.0 && cfg.diffusion_substeps > 0 {
            let dsub = cfg.diffusion / cfg.diffusion_substeps as f32;
            for _ in 0..cfg.diffusion_substeps {
                let src = field.data.clone();
                for y in 1..h - 1 {
                    for x in 1..w - 1 {
                        let k = y * w + x;
                        if src[k] <= cfg.sea_level {
                            continue;
                        }
                        if cfg.min_area_cells > 0.0
                            && flow.accumulation.data[k] >= cfg.min_area_cells
                        {
                            continue; // channel cell — stream power only, no diffusion
                        }
                        let lap = src[k - 1] + src[k + 1] + src[k - w] + src[k + w] - 4.0 * src[k];
                        field.data[k] = src[k] + dsub * lap;
                    }
                }
            }
        }

        progress(iter, &field);
    }
    field
}

#[cfg(test)]
mod tests {
    use super::*;

    /// On a simple inclined ramp with a notch, incision lowers channel cells and
    /// never pushes a cell below its receiver (monotone drainage preserved).
    #[test]
    fn incision_lowers_channels_without_reversal() {
        let (w, h) = (32usize, 32usize);
        let mut d = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                // ramp high at y=0 (0.9) to sea at y=h-1 (0.45), tiny notch mid-x.
                let t = y as f32 / (h as f32 - 1.0);
                let notch = if x == w / 2 { -0.01 } else { 0.0 };
                d[y * w + x] = 0.9 - 0.45 * t + notch;
            }
        }
        let f0 = GridF32::from_vec(w, h, d);
        let cfg = StreamPowerConfig { k: 2.0, iterations: 3, sea_level: 0.5, ..Default::default() };
        let f1 = incise(&f0, &cfg);
        // Some land cell was lowered (incision happened).
        let lowered = (0..w * h).filter(|&k| f1.data[k] < f0.data[k] - 1e-6).count();
        assert!(lowered > 0, "incision must lower some channel cells");
        // No NaN/Inf, and land stays ≥ sea (no spurious drowning of ridge tops).
        assert!(f1.data.iter().all(|v| v.is_finite()));
        assert!(f1.data[w / 2] > 0.5, "the high ridge top must stay land");
    }

    /// n = 1 closed form and the Newton path agree at n = 1.0000001.
    #[test]
    fn newton_matches_closed_form_at_n1() {
        let (w, h) = (16usize, 16usize);
        let mut d = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                d[y * w + x] = 0.9 - 0.4 * (y as f32 / (h as f32 - 1.0));
            }
        }
        let base = GridF32::from_vec(w, h, d);
        let a = incise(&base, &StreamPowerConfig { n: 1.0, iterations: 2, ..Default::default() });
        let b = incise(&base, &StreamPowerConfig { n: 1.0000001, iterations: 2, ..Default::default() });
        let maxdiff = (0..w * h).map(|k| (a.data[k] - b.data[k]).abs()).fold(0.0f32, f32::max);
        assert!(maxdiff < 1e-3, "closed form vs Newton at n≈1 diverged: {maxdiff}");
    }
}
