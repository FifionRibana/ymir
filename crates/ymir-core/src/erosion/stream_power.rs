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
use crate::terrain::flow::{
    D8_DIST, D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, mfd_accumulation,
};

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
    /// PHYSICAL cell size (km): makes the law dimensional — drainage area in km²
    /// (`A_cells · cell_km²`) and slope a true gradient (`Δh_m / dist_m`), so K is
    /// resolution-invariant instead of drifting ~1.5× per resolution doubling (the
    /// per-cell `Δnorm` slope measured steeper on finer cells). See ADR 0001 §D.
    pub cell_km: f32,
    /// Vertical scale (m) for the physical slope: `Δh_m = Δnorm · 2 · 1.13 ·
    /// depth_scale_m` (the `c1_altitude_norm_to_metres` slope).
    pub depth_scale_m: f32,

    // ─── CLOSURE (a): nonlinear hillslope diffusion with a critical slope ───────
    /// **Critical slope `S_c`** (dimensionless gradient = `tan(angle)`; `0` = OFF,
    /// use the legacy LINEAR Laplacian). Roering-type nonlinear flux `q = D·S /
    /// (1 − (S/S_c)²)`: the effective diffusivity DIVERGES as `S → S_c`, so no slope
    /// can exceed the critical angle — this is the closure that BOUNDS the maximum
    /// slope (arêtes/planed ridges), which nothing else in the pipeline does. Start
    /// at `tan(33°) ≈ 0.649`. Solved implicitly (see `hillslope_picard`) because the
    /// term is stiff near `S_c` and an explicit scheme blows up. See docs/adr/0001.
    pub critical_slope: f32,
    /// Outer lagged-diffusivity (Picard) iterations for the nonlinear-diffusion
    /// implicit solve. Each freezes the edge diffusivities from the current surface,
    /// then does a backward-Euler solve; a few (3) converge the nonlinearity.
    pub hillslope_picard: usize,
    /// Inner Jacobi sweeps per backward-Euler solve. The implicit operator is
    /// diagonally dominant (unconditionally stable for any step), so Jacobi converges;
    /// ~40 sweeps suffice at these grid sizes.
    pub hillslope_implicit_iters: usize,

    // ─── CLOSURE (b): channel lateral (bank) widening ──────────────────────────
    /// **Lateral bank erodibility `K_lat`** (`0` = OFF). After vertical incision,
    /// channel cells erode their two banks (perpendicular to flow) at a lateral
    /// stream power `K_lat·A_km²^m·S_lat` (S_lat = the physical bank gradient), so
    /// trunk valleys (high `A`) grow wide floors while headwaters (low `A`) stay
    /// narrow gorges — turning the 1-px/1000-m slit into a valley with a floor whose
    /// width grows downstream. Banks are never cut below the channel floor.
    pub lateral_erosion: f32,

    /// **Apply hillslope diffusion EVERYWHERE**, including channel cells (`false` =
    /// legacy regime split: diffusion skips `A ≥ A_c`). The regime split's hard channel
    /// exclusion is the non-physical element behind the Smith–Bretherton parallel-rilling
    /// comb (ADR 0001 Finding 8): a rill that captures `A_c` escapes the diffusion that
    /// should damp it, forever. The standard LEM `dz/dt = U − K·A^m·S^n + D∇²z` runs BOTH
    /// terms on every cell — diffusion is simply dominated by incision where `A` is large,
    /// so real channels/gorges survive while sub-threshold rills get damped. `true` sets a
    /// finite valley spacing from the D/K balance instead of imposing it via `A_c`.
    pub diffuse_channels: bool,

    // ─── CLOSURE (a′): TALUS / angle of repose (the C1-consistent alternative) ──
    /// **Talus repose slope `S_c`** (dimensionless gradient = `tan(angle)`; `0` = OFF).
    /// A CLOSURE, not a solver: sort land cells high→low and, in a single sweep, shed
    /// the excess of any drop steeper than `S_c` to the downhill neighbours (mass
    /// conserving). Guarantees `S ≤ S_c` by construction (up to residuals — see
    /// `talus_passes`), O(n log n), deterministic, no convergence to monitor — unlike
    /// the nonlinear diffusion's Gauss-Seidel. Produces STRAIGHT repose slopes (vs the
    /// diffusion's convex ones). Runs everywhere; it BACKFILLS NOTHING into hollows (it
    /// only bounds slope), so it should preserve headwater vallons the diffusion fills.
    /// See docs/adr/0001 Finding 10.
    pub talus_slope: f32,
    /// Talus sweeps. One sorted sweep can leave residual over-steep cells (a transfer can
    /// re-steepen a downstream slope); a few bounded passes clear them. If `k` is small and
    /// bounded it is still a closure; if unbounded it is a solver in disguise (measured).
    pub talus_passes: usize,
    /// Fraction of the worst per-cell excess moved per pass (`≤ 0.5` keeps it stable /
    /// non-reversing). Lower = gentler, more passes; higher = faster, risk of overshoot.
    pub talus_factor: f32,

    /// **MFD partition exponent `p`** for the incision drainage area (`None` = single-flow
    /// D8, legacy). When `Some(p)`, the area `A` in `E = K·A^m·S^n` comes from MULTIPLE-
    /// flow-direction accumulation (`slopeᵖ` split among lower neighbours) instead of D8.
    /// Dispersing `A` breaks the rill capture→incise→capture feedback, so the Smith–
    /// Bretherton comb never forms (ADR 0001 Finding 10 — attack the cause). `p → ∞` ≈ D8
    /// (comb returns); small `p` disperses (channels blur). The stack/receiver update and
    /// rivers/lakes stay D8 — MFD drives ONLY the incision area. Dispersion lowers peak
    /// `A`, so `K` usually needs raising to keep trunk incision.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mfd_exponent: Option<f32>,
}

/// Relief-v1 reference: physical critical drainage area (km²) for the channel head.
/// Stored in KM², converted to cells per resolution (a cells-based `A_c` is not
/// resolution-stable). **0.1 km²** — a realistic humid-temperate drainage density; the
/// earlier 7.6 km² put channel heads far too low, leaving the upper massif slopes
/// undissected (channels reached only ~7–9 % of peak vs ~94 % now). See docs/adr/0001.
pub const RELIEF_V1_A_C_KM2: f32 = 0.1;

impl StreamPowerConfig {
    /// The `ref/relief-streampower-v1` reference config (ADR 0001) — the first
    /// setting that produces legible relief: routed stream power on channels +
    /// hillslope diffusion, droplets OFF, uncoupled vertical scale. `cell_km2` is the
    /// run's cell area (`(domain_km / grid)²`), used to convert the physical `A_c`
    /// (7.6 km²) to cells at THIS resolution. Physical law: `cell_km` + `depth_scale_m`
    /// make the slope a true gradient and area km², so K is resolution-invariant.
    /// K=3.7e-4 (physical), m=0.5, n=1, iters=3, D=0.05, θ=0.
    pub fn relief_v1(cell_km2: f32, depth_scale_m: f32) -> Self {
        let cell_km2 = cell_km2.max(1e-9);
        let cell_km = cell_km2.sqrt();
        Self {
            k: RELIEF_V1_K,
            m: 0.5,
            n: 1.0,
            dt: 1.0,
            iterations: 2, // bounded incision → floors sit at a plausible fraction of
            // the local ridge (iters=3 planed them toward base level). See ADR §sculpt.
            sea_level: 0.5,
            diffusion: 0.05,
            diffusion_substeps: 4,
            min_area_cells: RELIEF_V1_A_C_KM2 / cell_km2,
            threshold: 0.0,
            cell_km,
            depth_scale_m,
            // Closures OFF in v1 (byte-identical to the reviewed sculpt).
            critical_slope: 0.0,
            hillslope_picard: 3,
            hillslope_implicit_iters: 40,
            lateral_erosion: 0.0,
            diffuse_channels: false,
            talus_slope: 0.0,
            talus_passes: 1,
            talus_factor: 0.5,
            mfd_exponent: None,
        }
    }

    /// `relief-v2` — v1 plus the two bounding closures (ADR 0001, Finding 7): nonlinear
    /// hillslope diffusion with a critical slope `S_c = tan(33°)` (bounds the maximum
    /// slope → arêtes replace the near-vertical FBM/slit faces) and channel lateral
    /// widening (`K_lat`, trunks get floors, headwaters stay gorges). Still OFF by
    /// default in the pipeline — driven by the viz checkbox / diagnostic pending the
    /// author's visual verdict. Same `A_c`, iters and K as v1.
    pub fn relief_v2(cell_km2: f32, depth_scale_m: f32) -> Self {
        // NB: the 40 GS sweeps are enough — scaling them with resolution (tested up to
        // 160 at 8192²) barely moved the residual steep share (33.7→32.9 %) at 4× the
        // cost, so GS convergence is NOT the 8192² lever. The residual upper-slope
        // striation at 8192² is FBM-TEMPLATE-driven: the finer FBM striae seed the
        // initial drainage, so the dense fine channel network (A_c≈42 cells) incises
        // ALONG the striae. The 2048² preview is clean. See ADR 0001 Finding 7.
        Self {
            critical_slope: RELIEF_V2_CRITICAL_SLOPE,
            lateral_erosion: RELIEF_V2_LATERAL,
            // The nonlinear closure needs more transport than the linear D=0.05 to plane
            // arêtes; 0.15 collapses the steep share while keeping drainage relief ~350 m
            // (D=0.3 over-planed it to ~240 m). See ADR 0001 Finding 7.
            diffusion: RELIEF_V2_DIFFUSION,
            ..Self::relief_v1(cell_km2, depth_scale_m)
        }
    }
}

/// `relief-v2` critical slope — `tan(33°)`, a mid-range angle of repose for fractured
/// rock. Above it the nonlinear flux diverges and the slope cannot steepen further.
pub const RELIEF_V2_CRITICAL_SLOPE: f32 = 0.6494; // tan(33°)
/// `relief-v2` lateral width coefficient (m per √km²): channel floor half-width =
/// `K_lat · A_km²^m`. Tuned so trunks (A~10⁴ km²) get ~0.4–0.8 km floors while
/// headwaters (A~1 km²) stay sub-cell gorges — the width variety the author asked for.
pub const RELIEF_V2_LATERAL: f32 = 4.0;
/// `relief-v2` hillslope diffusivity for the NONLINEAR closure (higher than v1's linear
/// 0.05 — the critical-slope denominator modulates it, and it must plane arêtes).
/// Dimensionless at the reference cell [`HILLSLOPE_REF_CELL_M`]; scaled ∝ 1/cell² inside.
pub const RELIEF_V2_DIFFUSION: f32 = 0.15;

/// Reference cell size (m) at which the dimensionless hillslope [`StreamPowerConfig::diffusion`]
/// is calibrated (2048² over a 400 km domain). The nonlinear implicit weight scales by
/// `(HILLSLOPE_REF_CELL_M / cell_m)²` so the closure planes the same metres at any grid.
pub const HILLSLOPE_REF_CELL_M: f32 = 400_000.0 / 2048.0;

/// Relief-v1 erodibility K in the PHYSICAL law (`E = K·A_km²^m·S_phys^n`). **1500** —
/// half the incision-reproducing K=3000, the "bounded incision" that keeps valley
/// floors at a plausible fraction of the local ridge (with iters=2) instead of
/// planing them to base level. See ADR 0001 §sculpt.
pub const RELIEF_V1_K: f32 = 1500.0;

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
            cell_km: 1.0,
            depth_scale_m: 5000.0,
            critical_slope: 0.0,
            hillslope_picard: 3,
            hillslope_implicit_iters: 40,
            lateral_erosion: 0.0,
            diffuse_channels: false,
            talus_slope: 0.0,
            talus_passes: 1,
            talus_factor: 0.5,
            mfd_exponent: None,
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
        // MFD area for the incision ONLY (D8 receiver/stack + rivers/lakes stay D8). When
        // set, the drainage area A is dispersed so the rilling feedback cannot run away.
        let mfd_acc = cfg
            .mfd_exponent
            .map(|p| mfd_accumulation(&flow.filled, &flow.direction, cfg.sea_level, p, w, h));
        let acc: &GridF32 = mfd_acc.as_ref().unwrap_or(&flow.accumulation);

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
        // PHYSICAL law: E = K · A_km²^m · S_phys^n, S_phys = Δh_m / dist_m. Heights
        // stay in NORM; K/θ carry the units. For n=1 the norm→m factor cancels in the
        // area term (f = kdt·A_km²^m / dist_m) and appears only on θ.
        let kdt = cfg.k * cfg.dt;
        let norm_to_m = 2.0 * 1.13 * cfg.depth_scale_m;
        let cell_m = cfg.cell_km * 1000.0;
        let cell_km2 = cfg.cell_km * cfg.cell_km;
        for &k in &stack {
            let r = receiver[k];
            if r == k || field.data[k] <= cfg.sea_level {
                continue; // base level, fixed
            }
            let area = acc.data[k].max(1.0);
            if area < cfg.min_area_cells {
                continue; // A_c: hillslope regime — no fluvial incision (channel head)
            }
            let hr = field.data[r];
            let ho = field.data[k];
            let am = (area * cell_km2).powf(cfg.m); // A_km²^m
            let dist_m = dist[k] * cell_m;
            // Threshold gate: physical stream power A_km²^m · S_phys^n vs θ.
            let s_now = (ho - hr).max(0.0) * norm_to_m / dist_m;
            if am * s_now.powf(cfg.n) <= cfg.threshold {
                continue;
            }
            let hn = if (cfg.n - 1.0).abs() < 1e-6 {
                // n=1: E_norm = K·A_km²^m·(h−hr)/dist_m; θ term = kdt·θ/norm_to_m.
                let f = kdt * am / dist_m;
                (ho + f * hr + kdt * cfg.threshold / norm_to_m) / (1.0 + f)
            } else {
                // Newton (general n), norm heights, physical slope.
                let c = norm_to_m / dist_m; // Δnorm → S_phys
                let mut x = ho;
                for _ in 0..8 {
                    let s = ((x - hr).max(0.0)) * c;
                    let e = (am * s.powf(cfg.n) - cfg.threshold).max(0.0);
                    let g = x - ho + kdt * e / norm_to_m;
                    let de = if am * s.powf(cfg.n) > cfg.threshold {
                        kdt * am * cfg.n * s.powf(cfg.n - 1.0) * c / norm_to_m
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
            // Never incise below the receiver or ABOVE the old height (no deposition).
            field.data[k] = hn.clamp(hr, ho);
        }

        // 4b. CLOSURE (b) — channel lateral (bank) widening, as HYDRAULIC GEOMETRY.
        // A channel's floor half-width scales with discharge: W ∝ Q^b ∝ A^b (the
        // standard width–area law; use b = m = 0.5). Each channel cell planes the banks
        // perpendicular to flow out to a PHYSICAL half-width `K_lat · A_km²^m` (metres),
        // so trunks (high A) grow wide floors and headwaters (low A) stay narrow gorges —
        // and the reach is in metres, so it is RESOLUTION-INVARIANT (a ±1-cell bank
        // erosion would widen 4× less at 4× finer cells → the slit returns at 8192²).
        // Banks are planed toward the channel floor, never below it. Order-independent
        // (delta buffer, deterministic min).
        if cfg.lateral_erosion > 0.0 {
            let mut newz = field.data.clone();
            for &k in &stack {
                let r = receiver[k];
                if r == k || field.data[k] <= cfg.sea_level {
                    continue;
                }
                let area = acc.data[k].max(1.0);
                if area < cfg.min_area_cells {
                    continue; // hillslope, no channel → no banks
                }
                let d = flow.direction[k];
                if d == DIR_NONE {
                    continue;
                }
                let (fdx, fdy) = (D8_DX[d as usize], D8_DY[d as usize]);
                let zc = field.data[k]; // channel floor
                let half_w_m = cfg.lateral_erosion * (area * cell_km2).powf(cfg.m);
                let half_cells = (half_w_m / cell_m).floor() as i32;
                if half_cells < 1 {
                    continue; // sub-cell floor (headwater gorge) — leave the 1-px channel
                }
                let (x, y) = (k % w, k / w);
                for (px, py) in [(-fdy, fdx), (fdy, -fdx)] {
                    for step in 1..=half_cells {
                        let (bx, by) = (x as i32 + px * step, y as i32 + py * step);
                        if bx < 0 || by < 0 || bx >= w as i32 || by >= h as i32 {
                            break;
                        }
                        let b = by as usize * w + bx as usize;
                        if field.data[b] <= cfg.sea_level || field.data[b] <= zc {
                            break; // reached sea, or ground already at/below the floor
                        }
                        if zc < newz[b] {
                            newz[b] = zc; // plane the bank down to the channel floor
                        }
                    }
                }
            }
            field.data = newz;
        }

        // 4c. CLOSURE (a′) — TALUS / angle of repose. Sort land cells high→low; in one
        // sweep, shed the excess of any drop steeper than `talus_slope` to the downhill
        // neighbours, mass-conserving (transfer, not carve). O(n log n), deterministic,
        // no convergence loop — the C1-consistent alternative to the nonlinear-diffusion
        // solver. `talus_passes` clears the residual re-steepening a single sweep leaves.
        // Applied EVERYWHERE (no channel exclusion; Finding 8). See ADR 0001 Finding 10.
        if cfg.talus_slope > 0.0 && cfg.talus_passes > 0 {
            let base_drop = cfg.talus_slope * cell_m / norm_to_m; // cardinal repose drop (norm)
            for _ in 0..cfg.talus_passes {
                // Sort land cells by height desc (deterministic tiebreak on index).
                let mut idx: Vec<usize> =
                    (0..n).filter(|&k| field.data[k] > cfg.sea_level).collect();
                idx.sort_unstable_by(|&a, &b| {
                    field.data[b]
                        .partial_cmp(&field.data[a])
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(a.cmp(&b))
                });
                for &k in &idx {
                    let (x, y) = ((k % w) as i32, (k / w) as i32);
                    let zk = field.data[k];
                    let (mut total, mut maxe) = (0.0f32, 0.0f32);
                    let mut exc = [0.0f32; 8];
                    for d in 0..8 {
                        let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                        if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                            continue;
                        }
                        let j = ny as usize * w + nx as usize;
                        let zj = field.data[j];
                        if zj >= zk {
                            continue; // downhill only
                        }
                        // max allowed drop to this neighbour = S_c · dist (norm units).
                        let e = (zk - zj) - base_drop * D8_DIST[d];
                        if e > 0.0 {
                            exc[d] = e;
                            total += e;
                            if e > maxe {
                                maxe = e;
                            }
                        }
                    }
                    if total <= 0.0 {
                        continue;
                    }
                    // Move a fraction of the WORST excess off k, split among the over-steep
                    // neighbours in proportion to their deficit (mass conserving).
                    let mv = cfg.talus_factor * maxe;
                    field.data[k] = zk - mv;
                    for d in 0..8 {
                        if exc[d] > 0.0 {
                            let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                            let j = ny as usize * w + nx as usize;
                            field.data[j] += mv * exc[d] / total;
                        }
                    }
                }
            }
        }

        // 5. Hillslope diffusion (explicit, a few sub-steps), interleaved with the
        // incision each iteration so it holds interfluves WHILE channels incise (the
        // coupled LEM, not a post-pass smooth). REGIME SPLIT: when `min_area_cells`
        // is set, diffusion runs ONLY on hillslope cells (A < A_c) and stream power
        // ONLY on channels (A ≥ A_c) — the standard channel-head partition that stops
        // the fluvial law from over-carving physical hillslopes. `min_area_cells = 0`
        // → diffuse everywhere (legacy).
        if cfg.diffusion > 0.0 && cfg.critical_slope > 0.0 {
            // CLOSURE (a) — NONLINEAR hillslope diffusion with critical slope S_c.
            // Flux q = D·S/(1−(S/S_c)²): the effective edge diffusivity diverges as
            // S→S_c, so slopes cannot exceed S_c (arêtes; the missing bound). The term
            // is STIFF near S_c → an explicit scheme blows up, so we solve it IMPLICITLY:
            // backward Euler (unconditionally stable for any step) with LAGGED-DIFFUSIVITY
            // Picard (re-freeze edge weights each outer pass) and a Gauss-Seidel inner
            // solve (the operator is diagonally dominant → GS converges, deterministic in
            // row-major order). The denominator is floored (edge slope capped at 0.999·S_c)
            // so weights stay finite; the self-arresting nature — a bank that drops below
            // S_c gets denom→1, weight→D, and stops — is what leaves relief instead of
            // planing to base level. Runs only on hillslope cells (A<A_c); channel/sea
            // cells are fixed Dirichlet values, so slit tops get pulled down into flanks.
            let sc = cfg.critical_slope;
            let g = norm_to_m / cell_m; // Δnorm(edge) → physical gradient
            // Resolution invariance: the implicit diffusion weight is κΔt/dx², so the
            // dimensionless `diffusion` (calibrated at the reference cell) is scaled by
            // (REF/cell)² — otherwise the same D smooths 16× fewer metres at 4× finer
            // cells and the closure fails to plane at 8192² (see ADR 0001 Finding 7).
            let dscale = (HILLSLOPE_REF_CELL_M / cell_m).powi(2);
            let w_base = cfg.diffusion * dscale;
            let z_old = field.data.clone(); // backward-Euler RHS (fixed over Picard)
            for _p in 0..cfg.hillslope_picard.max(1) {
                let coeff = field.data.clone(); // slopes lagged within this Picard pass
                for _j in 0..cfg.hillslope_implicit_iters.max(1) {
                    for y in 1..h - 1 {
                        for x in 1..w - 1 {
                            let k = y * w + x;
                            if field.data[k] <= cfg.sea_level
                                || (!cfg.diffuse_channels
                                    && cfg.min_area_cells > 0.0
                                    && acc.data[k] >= cfg.min_area_cells)
                            {
                                continue; // sea, or (regime split) channel = fixed boundary
                            }
                            let (mut sw, mut swz) = (0.0f32, 0.0f32);
                            for nb in [k - 1, k + 1, k - w, k + w] {
                                let s = ((coeff[k] - coeff[nb]).abs() * g / sc).min(0.999);
                                let w_e = w_base / (1.0 - s * s).max(0.02);
                                sw += w_e;
                                swz += w_e * field.data[nb]; // Gauss-Seidel: latest values
                            }
                            field.data[k] = (z_old[k] + swz) / (1.0 + sw);
                        }
                    }
                }
            }
        } else if cfg.diffusion > 0.0 && cfg.diffusion_substeps > 0 {
            let dsub = cfg.diffusion / cfg.diffusion_substeps as f32;
            for _ in 0..cfg.diffusion_substeps {
                let src = field.data.clone();
                for y in 1..h - 1 {
                    for x in 1..w - 1 {
                        let k = y * w + x;
                        if src[k] <= cfg.sea_level {
                            continue;
                        }
                        if !cfg.diffuse_channels
                            && cfg.min_area_cells > 0.0
                            && acc.data[k] >= cfg.min_area_cells
                        {
                            continue; // regime split: channel = stream power only
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

    /// CLOSURE (a): nonlinear diffusion with a critical slope bounds the maximum
    /// slope — a cliff far above `S_c` is planed toward `S_c` and never left steeper
    /// than the linear scheme would leave it; result stays finite.
    #[test]
    fn critical_slope_bounds_max_slope() {
        let (w, h) = (24usize, 24usize);
        // A plateau (0.9) dropping to a bench (0.55) at mid-x: a near-vertical cliff.
        let mut d = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                d[y * w + x] = if x < w / 2 { 0.9 } else { 0.55 };
            }
        }
        let f0 = GridF32::from_vec(w, h, d);
        // No incision (K=0), diffusion only, so we isolate the hillslope closure.
        let cell_km = 0.2f32;
        let base = StreamPowerConfig {
            k: 0.0,
            iterations: 3,
            sea_level: 0.5,
            diffusion: 0.3,
            cell_km,
            depth_scale_m: 5000.0,
            ..Default::default()
        };
        let max_slope = |f: &GridF32| -> f32 {
            let norm_to_m = 2.0 * 1.13 * base.depth_scale_m;
            let cell_m = cell_km * 1000.0;
            // Interior only — boundary rows/cols are not diffused (loop is 1..w-1).
            (0..w * h)
                .filter(|&k| {
                    let (x, y) = (k % w, k / w);
                    x > 1 && y > 0 && y < h - 1
                })
                .map(|k| (f.data[k] - f.data[k - 1]).abs() * norm_to_m / cell_m)
                .fold(0.0f32, f32::max)
        };
        let sc = 0.6494f32; // tan(33°)
        let nonlin = StreamPowerConfig { critical_slope: sc, ..base.clone() };
        let a = incise(&f0, &base); // linear
        let b = incise(&f0, &nonlin); // nonlinear, bounded
        assert!(b.data.iter().all(|v| v.is_finite()), "nonlinear diffusion must stay finite");
        // The nonlinear closure planes the cliff MORE than linear (it targets S_c),
        // so its residual max slope is lower.
        assert!(
            max_slope(&b) < max_slope(&a),
            "critical-slope diffusion should reduce the max slope below linear: {} vs {}",
            max_slope(&b),
            max_slope(&a)
        );
    }

    /// CLOSURE (b): lateral erosion lowers channel banks, widening the corridor —
    /// more cells adjacent to the channel are lowered than with vertical incision alone.
    #[test]
    fn lateral_erosion_widens_channel() {
        let (w, h) = (40usize, 40usize);
        let cx = (w / 2) as f32;
        // A V-shaped valley: high on the flanks (|x−cx|), sloping down to the sea at
        // y=h−1. Off-centre cells drain TOWARD the centre column, so the centre is the
        // high-accumulation channel and the flanks are hillslopes (low accumulation, no
        // vertical incision) — lateral erosion is then the ONLY thing that lowers banks.
        let mut d = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let t = y as f32 / (h as f32 - 1.0);
                let flank = (x as f32 - cx).abs() / cx;
                d[y * w + x] = 0.62 - 0.10 * t + 0.28 * flank;
            }
        }
        let f0 = GridF32::from_vec(w, h, d);
        let base = StreamPowerConfig {
            k: 2.0,
            iterations: 3,
            sea_level: 0.5,
            cell_km: 0.2,
            depth_scale_m: 5000.0,
            ..Default::default()
        };
        // lateral_erosion is now a half-width coefficient (m per √km²); use a large
        // value so the physical floor half-width exceeds the 200 m cell here.
        let widened = StreamPowerConfig { lateral_erosion: 200.0, ..base.clone() };
        let a = incise(&f0, &base);
        let b = incise(&f0, &widened);
        // Total lowering on the two columns adjacent to the centre channel: lateral
        // erosion pulls the banks down toward the channel floor, deeper than vertical
        // incision alone (the count saturates, so measure depth).
        let bank_incision = |f: &GridF32| -> f32 {
            (1..h - 1)
                .flat_map(|y| [(w / 2 - 1, y), (w / 2 + 1, y)])
                .map(|(x, y)| (f0.data[y * w + x] - f.data[y * w + x]).max(0.0))
                .sum()
        };
        assert!(b.data.iter().all(|v| v.is_finite()));
        assert!(
            bank_incision(&b) > bank_incision(&a) * 1.05,
            "lateral erosion should deepen the banks: {} vs {}",
            bank_incision(&b),
            bank_incision(&a)
        );
    }

    /// CLOSURE (a′): talus bounds the max slope (angle of repose) and CONSERVES MASS
    /// (it transfers, never carves). A tall central spike is reduced to ≤ S_c to its
    /// neighbours within a few bounded passes; the total height sum is unchanged.
    #[test]
    fn talus_bounds_slope_and_conserves_mass() {
        let (w, h) = (32usize, 32usize);
        let mut d = vec![0.7f32; w * h];
        d[16 * w + 16] = 0.97; // a tall spike on a plateau
        let f0 = GridF32::from_vec(w, h, d);
        let sc = 0.6494f32; // tan(33°)
        let cell_km = 0.2f32;
        let cfg = StreamPowerConfig {
            k: 0.0, // no incision — isolate the talus closure
            iterations: 1,
            sea_level: 0.5,
            cell_km,
            depth_scale_m: 5000.0,
            talus_slope: sc,
            talus_passes: 30,
            talus_factor: 0.5,
            mfd_exponent: None,
            ..Default::default()
        };
        let f1 = incise(&f0, &cfg);
        assert!(f1.data.iter().all(|v| v.is_finite()));
        // Mass conserved (transfers only; spike not on the boundary).
        let (s0, s1): (f64, f64) =
            (f0.data.iter().map(|&v| v as f64).sum(), f1.data.iter().map(|&v| v as f64).sum());
        assert!((s0 - s1).abs() < 1e-3, "talus must conserve mass: {s0} vs {s1}");
        // Max cardinal slope now ≤ S_c (within tolerance) — the spike is a repose cone.
        let norm_to_m = 2.0 * 1.13 * cfg.depth_scale_m;
        let cell_m = cell_km * 1000.0;
        let max_s = (0..w * h)
            .filter(|&k| k % w > 0)
            .map(|k| (f1.data[k] - f1.data[k - 1]).abs() * norm_to_m / cell_m)
            .fold(0.0f32, f32::max);
        assert!(max_s <= sc * 1.15, "talus should bound slope to ~S_c: {max_s} vs {sc}");
        assert!(f1.data[16 * w + 16] < f0.data[16 * w + 16], "the spike must be lowered");
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
        let b =
            incise(&base, &StreamPowerConfig { n: 1.0000001, iterations: 2, ..Default::default() });
        let maxdiff = (0..w * h).map(|k| (a.data[k] - b.data[k]).abs()).fold(0.0f32, f32::max);
        assert!(maxdiff < 1e-3, "closed form vs Newton at n≈1 diverged: {maxdiff}");
    }
}
