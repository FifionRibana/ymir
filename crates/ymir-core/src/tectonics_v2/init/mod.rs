//! Step 8.6 Phase 8a — S̃ field initialisation modes.
//!
//! Pre-Phase-8a, the harness initialised S̃ via a deterministic
//! sinusoidal perturbation around per-plate-type means
//! ([`InitMode::Checkerboard`] in this module). The pattern preserves
//! Steps 0–10 numerical baselines but produces a global sinusoidal
//! signature that pollutes visual review (see Step 8.6 Phase 7
//! follow-up: continents appear "circular" rather than Voronoï-shaped).
//!
//! This module formalises the legacy algorithm as one of four init
//! modes and adds three alternatives that follow TDD §4.2 more
//! closely (flat per-plate-type, no sinusoidal artefact):
//!
//! - [`InitMode::Uniform`] — flat per-plate-type with optional
//!   smoothstep blending at inter-plate boundaries.
//! - [`InitMode::Gaussian`] — peak at each plate's Voronoï centroid,
//!   decaying Gaussian with periodic-aware Euclidean distance.
//! - [`InitMode::Convolution`] — binary classification mask convolved
//!   with a periodic Gaussian kernel.
//!
//! Default: [`InitMode::Uniform`] with `boundary_smoothing_width = 1.0`.
//!
//! Steps 0–10 regression tests opt into [`InitMode::Checkerboard`]
//! explicitly via their config builders to preserve numerical
//! baselines (strategy γ).

use serde::{Deserialize, Serialize};

use super::boundaries::{PlateType, PlateTypeField};
use super::field::Field2D;
use super::voronoi::{compute_dist_to_inter_plate_boundary, PlateIdField};

pub mod radial_profile;
pub mod radial_profile_fbm;
pub use radial_profile::{
    ProfileShape, CONTINENTAL_VALUE_DEFAULT, OCEANIC_VALUE_DEFAULT, POW_EXPONENT_DEFAULT,
};
pub use radial_profile_fbm::{
    FBM_AMPLITUDE_DEFAULT, FBM_AMPLITUDE_OCEANIC_DEFAULT, FBM_LACUNARITY_DEFAULT,
    FBM_OCTAVES_DEFAULT, FBM_PERSISTENCE_DEFAULT, FBM_SCALE_DEFAULT,
    FBM_SEED_DEFAULT, FBM_SEED_OCEANIC_XOR_MAGIC, OCEANIC_CLAMP_MAX,
};

/// Per-plate-type reference S̃ values, dimensionless. `0.2` for
/// oceanic (≈ 7 km), `1.0` for continental (≈ 35 km). Shared by the
/// new init modes so the same scale family applies as in the legacy
/// `Checkerboard` mode.
pub const OCEANIC_S_DEFAULT: f64 = 0.2;
pub const CONTINENTAL_S_DEFAULT: f64 = 1.0;

/// S̃ initialisation mode, persisted on [`super::diagnostics::harness::BaselineConfig`].
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum InitMode {
    /// Legacy: deterministic sinusoidal perturbation around
    /// per-plate-type means. Preserves Steps 0–10 numerical baselines
    /// bit-for-bit. Falls back to a plate-agnostic variant when no
    /// plate data is provided (Steps 0–4 callers).
    Checkerboard,
    /// Flat per-plate-type with smoothstep transition at inter-plate
    /// boundaries. `boundary_smoothing_width` is the half-width of the
    /// transition zone, in cells.
    Uniform { boundary_smoothing_width: f64 },
    /// Gaussian decay from each plate's Voronoï seed (centroid).
    /// `sigma_*` measured in cells; per-type sigmas allow different
    /// continental and oceanic profile widths.
    Gaussian {
        sigma_continental: f64,
        sigma_oceanic: f64,
    },
    /// Convolution of binary classification mask with a periodic
    /// Gaussian kernel. `sigma` measured in cells.
    Convolution { sigma: f64 },
    /// Step 13 — radial profile per continental plate. Continental
    /// cells get `S̃ = oceanic_value + (continental_value -
    /// oceanic_value) · profile(d / L_plate)` where `d` is the
    /// Chebyshev BFS distance to the nearest inter-plate boundary
    /// and `L_plate` is the per-plate max distance. Oceanic cells
    /// get `S̃ = oceanic_value` uniform. See
    /// [`radial_profile`] module docstring for the algorithm and
    /// degenerate-case behaviour.
    RadialProfile {
        continental_value: f64,
        oceanic_value: f64,
        profile_shape: ProfileShape,
    },
    /// Step 13 — radial profile + isotropic FBM noise on
    /// continental cells, for intra-plate thickness heterogeneity
    /// (province texture).
    ///
    /// Step 13.5 — extends with optional FBM on **oceanic** cells
    /// (opt-in via `apply_fbm_to_oceanic`). When the flag is
    /// `false` (default), the variant is bit-identical to its
    /// Step 13 form: oceanic cells stay at `oceanic_value` uniform.
    /// When `true`, oceanic cells receive
    /// `clamp(oceanic_value + fbm_amplitude_oceanic ·
    /// fbm_oceanic.get(x, y), 0, OCEANIC_CLAMP_MAX)` with
    /// `OCEANIC_CLAMP_MAX = 0.49` strictly preventing
    /// threshold-crossing to continental classification (D7;
    /// volcanic islands are out of scope for Step 13.5).
    ///
    /// Continental output clamped to `[0, 1]`. See
    /// [`radial_profile_fbm`] module docstring.
    ///
    /// Explicit `rename` overrides serde's default `snake_case`
    /// expansion of `RadialProfileWithFBM` (which would split
    /// "FBM" into `f_b_m`) — keeps the on-disk JSON kind tag
    /// readable as `radial_profile_with_fbm`.
    #[serde(rename = "radial_profile_with_fbm")]
    RadialProfileWithFBM {
        continental_value: f64,
        oceanic_value: f64,
        profile_shape: ProfileShape,
        fbm_amplitude: f64,
        fbm_octaves: u8,
        fbm_persistence: f64,
        fbm_lacunarity: f64,
        fbm_scale: f64,
        fbm_seed: u64,
        // Step 13.5 — oceanic FBM extension. `#[serde(default)]`
        // on each new field so legacy preset JSON written before
        // this step (which lacks the new keys) deserialises with
        // safe defaults: `apply_fbm_to_oceanic = false` →
        // bit-identical Step 13 behaviour. The other defaults
        // are only consumed when the flag is flipped on.
        #[serde(default)]
        apply_fbm_to_oceanic: bool,
        #[serde(default = "default_fbm_amplitude_oceanic")]
        fbm_amplitude_oceanic: f64,
        #[serde(default)]
        fbm_scale_oceanic: Option<f64>,
        #[serde(default)]
        fbm_seed_oceanic: Option<u64>,
    },
}

/// `#[serde(default)]` helper for the
/// `RadialProfileWithFBM::fbm_amplitude_oceanic` field — bare
/// `Default::default()` on `f64` is `0.0`, which would silently
/// disable the oceanic FBM perturbation when the user later flips
/// the flag on. The constant from `radial_profile_fbm` keeps the
/// value in one place.
fn default_fbm_amplitude_oceanic() -> f64 {
    radial_profile_fbm::FBM_AMPLITUDE_OCEANIC_DEFAULT
}

impl Default for InitMode {
    fn default() -> Self {
        InitMode::Uniform { boundary_smoothing_width: 1.0 }
    }
}

/// Voronoï plate metadata accessible to non-`Checkerboard` modes.
/// All fields borrow from `BoundaryConfig::Enabled.geometry`. The
/// `seed_coords` slice is `None` when the geometry is built from a
/// static layout (no Voronoï tessellation): in that case
/// [`InitMode::Gaussian`] is unavailable and panics with a clear
/// message.
pub struct PlateInitData<'a> {
    pub plate_id: &'a PlateIdField,
    pub plate_type: &'a PlateTypeField,
    pub seed_coords: Option<&'a [(f64, f64)]>,
}

/// Inputs available to all init modes. `plate_data` is `None` for
/// boundary-disabled (Steps 0–4) runs.
pub struct InitContext<'a> {
    pub nx: usize,
    pub ny: usize,
    pub seed: u64,
    /// Sinusoidal-perturbation amplitude. Honoured by
    /// [`InitMode::Checkerboard`] only; the other modes ignore it
    /// (their S̃ patterns do not derive from a perturbation around a
    /// uniform mean).
    pub amplitude: f64,
    pub plate_data: Option<PlateInitData<'a>>,
}

/// Initialise the S̃ field according to `mode`.
///
/// Panics with a helpful message if the mode requires plate data
/// (Uniform / Gaussian / Convolution) but `ctx.plate_data` is `None`.
/// The intended pairing is `BoundaryConfig::Enabled` ↔ any mode,
/// `BoundaryConfig::Disabled` ↔ `InitMode::Checkerboard`.
pub fn init_s_field(mode: InitMode, ctx: &InitContext<'_>) -> Field2D {
    match mode {
        InitMode::Checkerboard => match &ctx.plate_data {
            Some(p) => checkerboard_plate_aware(
                ctx.nx,
                ctx.ny,
                ctx.seed,
                ctx.amplitude,
                p.plate_type,
            ),
            None => checkerboard_agnostic(ctx.nx, ctx.ny, ctx.seed, ctx.amplitude),
        },
        InitMode::Uniform { boundary_smoothing_width } => {
            let p = ctx.plate_data.as_ref().expect(
                "InitMode::Uniform requires plate data — pair with \
                 BoundaryConfig::Enabled, or use InitMode::Checkerboard \
                 for boundary-disabled runs",
            );
            uniform(ctx.nx, ctx.ny, p, boundary_smoothing_width)
        }
        InitMode::Gaussian { sigma_continental, sigma_oceanic } => {
            let p = ctx.plate_data.as_ref().expect(
                "InitMode::Gaussian requires plate data — pair with \
                 BoundaryConfig::Enabled",
            );
            let coords = p.seed_coords.expect(
                "InitMode::Gaussian requires Voronoï seed coordinates \
                 — only available for BoundaryConfig::enabled_voronoi_closed",
            );
            gaussian(ctx.nx, ctx.ny, p, coords, sigma_continental, sigma_oceanic)
        }
        InitMode::Convolution { sigma } => {
            let p = ctx.plate_data.as_ref().expect(
                "InitMode::Convolution requires plate data — pair with \
                 BoundaryConfig::Enabled",
            );
            convolution(ctx.nx, ctx.ny, p, sigma)
        }
        InitMode::RadialProfile {
            continental_value,
            oceanic_value,
            profile_shape,
        } => {
            let p = ctx.plate_data.as_ref().expect(
                "InitMode::RadialProfile requires plate data — pair with \
                 BoundaryConfig::Enabled",
            );
            radial_profile::build(
                ctx.nx,
                ctx.ny,
                p,
                continental_value,
                oceanic_value,
                profile_shape,
            )
        }
        InitMode::RadialProfileWithFBM {
            continental_value,
            oceanic_value,
            profile_shape,
            fbm_amplitude,
            fbm_octaves,
            fbm_persistence,
            fbm_lacunarity,
            fbm_scale,
            fbm_seed,
            apply_fbm_to_oceanic,
            fbm_amplitude_oceanic,
            fbm_scale_oceanic,
            fbm_seed_oceanic,
        } => {
            let p = ctx.plate_data.as_ref().expect(
                "InitMode::RadialProfileWithFBM requires plate data — pair with \
                 BoundaryConfig::Enabled",
            );
            radial_profile_fbm::build(
                ctx.nx,
                ctx.ny,
                p,
                continental_value,
                oceanic_value,
                profile_shape,
                fbm_amplitude,
                fbm_octaves,
                fbm_persistence,
                fbm_lacunarity,
                fbm_scale,
                fbm_seed,
                apply_fbm_to_oceanic,
                fbm_amplitude_oceanic,
                fbm_scale_oceanic,
                fbm_seed_oceanic,
            )
        }
    }
}

#[inline]
fn s_value_for(plate_type: PlateType) -> f64 {
    match plate_type {
        PlateType::Oceanic => OCEANIC_S_DEFAULT,
        PlateType::Continental => CONTINENTAL_S_DEFAULT,
    }
}

#[inline]
fn legacy_phase(seed: u64) -> f64 {
    use std::f64::consts::PI;
    ((seed.wrapping_mul(2654435761u64)) as f64) / (u64::MAX as f64) * 2.0 * PI
}

fn checkerboard_agnostic(nx: usize, ny: usize, seed: u64, amplitude: f64) -> Field2D {
    use std::f64::consts::PI;
    let phase = legacy_phase(seed);
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) / nx as f64;
            let y = (j as f64 + 0.5) / ny as f64;
            let bump = 1.0 + amplitude * ((2.0 * PI * x + phase).sin() * (2.0 * PI * y).cos());
            s.set(i, j, bump);
        }
    }
    s
}

fn checkerboard_plate_aware(
    nx: usize,
    ny: usize,
    seed: u64,
    amplitude: f64,
    plate_types: &PlateTypeField,
) -> Field2D {
    use std::f64::consts::PI;
    let phase = legacy_phase(seed);
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) / nx as f64;
            let y = (j as f64 + 0.5) / ny as f64;
            let bump_scale = (2.0 * PI * x + phase).sin() * (2.0 * PI * y).cos();
            let mean = s_value_for(plate_types.get(i, j));
            s.set(i, j, mean + amplitude * mean * bump_scale);
        }
    }
    s
}

/// Flat per-plate-type with smoothstep blending across inter-plate
/// boundaries. Outside the `boundary_smoothing_width` band the cell
/// holds its own plate's reference value exactly; inside the band the
/// value blends towards the midpoint between own and across-boundary
/// values, with cubic smoothstep weight `3t² − 2t³` of normalised
/// distance.
///
/// Step 13 Phase 1: the BFS distance computation is delegated to
/// [`super::voronoi::compute_dist_to_inter_plate_boundary`] (shared
/// with `plate_kinematic::field::build` and the upcoming
/// `init::radial_profile`). Bit-identical with the pre-refactor
/// implementation by construction — see that utility's module
/// docstring.
fn uniform(
    nx: usize,
    ny: usize,
    p: &PlateInitData<'_>,
    boundary_smoothing_width: f64,
) -> Field2D {
    let bfs = compute_dist_to_inter_plate_boundary(nx, ny, p.plate_id);

    // Per-cell own value + per-plate-id S̃ value lookup, in one scan.
    // PlateInitData carries a per-cell plate_type field but no
    // per-plate table; build the lookup inline so the utility's
    // `target_plate_id` index can be translated to the across-
    // boundary S̃ reference value. Bit-identical with storing the
    // value directly during BFS because per-plate properties are
    // constant within a plate (`plate_type[i,j] =
    // per_plate_type[plate_id[i,j]]`).
    let n = nx * ny;
    let mut own_value = vec![0.0_f64; n];
    let mut per_plate_value: Vec<Option<f64>> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let v = s_value_for(p.plate_type.get(i, j));
            own_value[j * nx + i] = v;
            let pid = p.plate_id.get(i, j) as usize;
            if pid >= per_plate_value.len() {
                per_plate_value.resize(pid + 1, None);
            }
            if per_plate_value[pid].is_none() {
                per_plate_value[pid] = Some(v);
            }
        }
    }

    let w = boundary_smoothing_width.max(1e-12);
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let idx = j * nx + i;
            let d = bfs.distance.get(i, j);
            let own = own_value[idx];
            let tpid = bfs.target_plate_id[idx];
            let other = if tpid == u16::MAX {
                // BFS never reached this cell (degenerate single-
                // plate-on-torus). dist=INFINITY → the `d >= w`
                // branch short-circuits to `own`; the value picked
                // here is unused. Defensive default = own.
                own
            } else {
                per_plate_value[tpid as usize].expect(
                    "BFS propagated a plate id without an entry in \
                     per_plate_value — should be impossible since \
                     target_plate_id is set from plate_id.get(ni, nj)",
                )
            };
            if d >= w {
                s.set(i, j, own);
            } else {
                let t = (d / w).clamp(0.0, 1.0);
                let st = t * t * (3.0 - 2.0 * t);
                let midpoint = 0.5 * (own + other);
                let blend = own * st + midpoint * (1.0 - st);
                s.set(i, j, blend);
            }
        }
    }
    s
}

/// Per-plate Gaussian peaked at the plate's Voronoï seed coordinate.
/// Each cell takes its own plate's `peak * exp(-d² / (2σ²))` where
/// `d` is the periodic minimum-image distance (in cells) from the
/// cell centre `(i + 0.5, j + 0.5)` to its plate's seed.
fn gaussian(
    nx: usize,
    ny: usize,
    p: &PlateInitData<'_>,
    seed_coords: &[(f64, f64)],
    sigma_continental: f64,
    sigma_oceanic: f64,
) -> Field2D {
    let nx_f = nx as f64;
    let ny_f = ny as f64;
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let pid = p.plate_id.get(i, j) as usize;
            let (sx, sy) = seed_coords[pid];
            let cx = i as f64 + 0.5;
            let cy = j as f64 + 0.5;
            let mut dx = (cx - sx).abs();
            let mut dy = (cy - sy).abs();
            if dx > 0.5 * nx_f {
                dx = nx_f - dx;
            }
            if dy > 0.5 * ny_f {
                dy = ny_f - dy;
            }
            let d2 = dx * dx + dy * dy;
            let pt = p.plate_type.get(i, j);
            let (peak, sigma) = match pt {
                PlateType::Oceanic => (OCEANIC_S_DEFAULT, sigma_oceanic),
                PlateType::Continental => (CONTINENTAL_S_DEFAULT, sigma_continental),
            };
            let g = (-d2 / (2.0 * sigma * sigma)).exp();
            s.set(i, j, peak * g);
        }
    }
    s
}

/// Periodic Gaussian blur applied to the binary classification mask.
/// Separable kernel truncated at 3σ, normalised to unit sum so the
/// output stays in `[OCEANIC_S_DEFAULT, CONTINENTAL_S_DEFAULT]`.
fn convolution(nx: usize, ny: usize, p: &PlateInitData<'_>, sigma: f64) -> Field2D {
    let n = nx * ny;
    let mut raw = vec![0.0_f64; n];
    for j in 0..ny {
        for i in 0..nx {
            raw[j * nx + i] = s_value_for(p.plate_type.get(i, j));
        }
    }

    let radius = ((3.0 * sigma).ceil() as i32).max(1);
    let kernel_size = (2 * radius + 1) as usize;
    let mut kernel = vec![0.0_f64; kernel_size];
    let two_sigma_sq = 2.0 * sigma * sigma;
    let mut sum = 0.0;
    for k in 0..kernel_size {
        let d = (k as i32 - radius) as f64;
        let v = (-d * d / two_sigma_sq).exp();
        kernel[k] = v;
        sum += v;
    }
    for k in &mut kernel {
        *k /= sum;
    }

    // Horizontal pass (periodic in i).
    let mut tmp = vec![0.0_f64; n];
    for j in 0..ny {
        for i in 0..nx {
            let mut acc = 0.0;
            for k in 0..kernel_size {
                let di = k as i32 - radius;
                let ni = ((i as i32 + di).rem_euclid(nx as i32)) as usize;
                acc += kernel[k] * raw[j * nx + ni];
            }
            tmp[j * nx + i] = acc;
        }
    }

    // Vertical pass (periodic in j).
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let mut acc = 0.0;
            for k in 0..kernel_size {
                let dj = k as i32 - radius;
                let nj = ((j as i32 + dj).rem_euclid(ny as i32)) as usize;
                acc += kernel[k] * tmp[nj * nx + i];
            }
            s.set(i, j, acc);
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

    fn build_test_plates(nx: usize, ny: usize, seed: u64) -> crate::tectonics_v2::voronoi::VoronoiPlates {
        let cfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
        generate_voronoi(nx, ny, &cfg, seed)
    }

    fn ctx_with_plates<'a>(
        nx: usize,
        ny: usize,
        seed: u64,
        amplitude: f64,
        plates: &'a crate::tectonics_v2::voronoi::VoronoiPlates,
    ) -> InitContext<'a> {
        InitContext {
            nx,
            ny,
            seed,
            amplitude,
            plate_data: Some(PlateInitData {
                plate_id: &plates.plate_id,
                plate_type: &plates.plate_type,
                seed_coords: Some(&plates.seed_coords),
            }),
        }
    }

    /// Bit-identical reproduction of the legacy plate-agnostic init
    /// (`init_thickness` from harness.rs). This is the contract that
    /// allows Steps 0–4 regression tests to opt into Checkerboard.
    #[test]
    fn checkerboard_agnostic_legacy_preserved() {
        use std::f64::consts::PI;
        let nx = 32;
        let ny = 32;
        let seed = 42;
        let amplitude = 0.2;

        let ctx = InitContext { nx, ny, seed, amplitude, plate_data: None };
        let s = init_s_field(InitMode::Checkerboard, &ctx);

        let phase = legacy_phase(seed);
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) / nx as f64;
                let y = (j as f64 + 0.5) / ny as f64;
                let expected =
                    1.0 + amplitude * ((2.0 * PI * x + phase).sin() * (2.0 * PI * y).cos());
                assert!(
                    (s.get(i, j) - expected).abs() < 1e-15,
                    "checkerboard_agnostic mismatch at ({},{}): got {} expected {}",
                    i,
                    j,
                    s.get(i, j),
                    expected
                );
            }
        }
    }

    /// Bit-identical reproduction of the legacy plate-aware init
    /// (`init_thickness_plate_aware` from harness.rs). Required for
    /// Steps 5–10 regression tests using Checkerboard.
    #[test]
    fn checkerboard_plate_aware_legacy_preserved() {
        use std::f64::consts::PI;
        let nx = 32;
        let ny = 32;
        let seed = 42;
        let amplitude = 0.2;
        let plates = build_test_plates(nx, ny, seed);

        let ctx = ctx_with_plates(nx, ny, seed, amplitude, &plates);
        let s = init_s_field(InitMode::Checkerboard, &ctx);

        let phase = legacy_phase(seed);
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) / nx as f64;
                let y = (j as f64 + 0.5) / ny as f64;
                let bump_scale = (2.0 * PI * x + phase).sin() * (2.0 * PI * y).cos();
                let mean = match plates.plate_type.get(i, j) {
                    PlateType::Oceanic => 0.2,
                    PlateType::Continental => 1.0,
                };
                let expected = mean + amplitude * mean * bump_scale;
                assert!(
                    (s.get(i, j) - expected).abs() < 1e-15,
                    "checkerboard_plate_aware mismatch at ({},{}): got {} expected {}",
                    i,
                    j,
                    s.get(i, j),
                    expected
                );
            }
        }
    }

    /// Uniform mode produces values strictly in `[0.2, 1.0]` (no
    /// sinusoidal overshoot) and at least some interior cells equal
    /// the per-plate reference exactly.
    #[test]
    fn uniform_no_sinusoidal_artefacts() {
        let nx = 64;
        let ny = 64;
        let plates = build_test_plates(nx, ny, 42);
        let ctx = ctx_with_plates(nx, ny, 42, 0.0, &plates);
        let s = init_s_field(
            InitMode::Uniform { boundary_smoothing_width: 1.0 },
            &ctx,
        );

        for &v in s.data() {
            assert!(
                (OCEANIC_S_DEFAULT - 1e-12..=CONTINENTAL_S_DEFAULT + 1e-12).contains(&v),
                "uniform out of [0.2, 1.0]: {}",
                v
            );
        }

        let interior = s
            .data()
            .iter()
            .filter(|&&v| {
                (v - OCEANIC_S_DEFAULT).abs() < 1e-12
                    || (v - CONTINENTAL_S_DEFAULT).abs() < 1e-12
            })
            .count();
        assert!(
            interior > 0,
            "uniform produced no purely-interior cells in 64² grid (8 plates)"
        );
    }

    /// Gaussian peaks at each plate's Voronoï seed and decays with
    /// distance. The cell containing the seed must hold a value
    /// close to the per-type peak.
    #[test]
    fn gaussian_decays_from_centroid() {
        let nx = 64;
        let ny = 64;
        let plates = build_test_plates(nx, ny, 42);
        let sigma = 5.0;
        let ctx = ctx_with_plates(nx, ny, 42, 0.0, &plates);
        let s = init_s_field(
            InitMode::Gaussian {
                sigma_continental: sigma,
                sigma_oceanic: sigma,
            },
            &ctx,
        );

        for plate_idx in 0..plates.num_plates {
            let pt = plates.per_plate_type[plate_idx];
            let peak = match pt {
                PlateType::Oceanic => OCEANIC_S_DEFAULT,
                PlateType::Continental => CONTINENTAL_S_DEFAULT,
            };
            let (sx, sy) = plates.seed_coords[plate_idx];
            let i_c = (sx.floor() as usize).min(nx - 1);
            let j_c = (sy.floor() as usize).min(ny - 1);
            // The seed cell must (a) belong to the plate (true by
            // construction of `generate_voronoi` for the seed cell)
            // and (b) hold a near-peak value with σ=5 cells, where
            // `g(0.5²+0.5²) = exp(-0.5/50) ≈ 0.99 · peak`.
            assert_eq!(plates.plate_id.get(i_c, j_c) as usize, plate_idx);
            let v = s.get(i_c, j_c);
            assert!(v <= peak + 1e-12);
            assert!(
                v >= 0.95 * peak,
                "plate {} seed cell value {} too far below peak {}",
                plate_idx,
                v,
                peak
            );

            // Decay: a cell sigma cells away from the seed (clamped to
            // grid + still inside the plate or not — distance test only)
            // should hold ≈ exp(-1/2) · peak ≈ 0.61 · peak. Pick a
            // direction that stays inside the field.
            let (sx_off, sy_off) = (sx + sigma, sy);
            let i_o = (sx_off.floor() as usize) % nx;
            let j_o = (sy_off.floor() as usize) % ny;
            // We only check decay magnitude for cells *of the same
            // plate*; if the offset cell is on another plate, its
            // sigma differs, skip it.
            if plates.plate_id.get(i_o, j_o) as usize == plate_idx {
                let v_off = s.get(i_o, j_o);
                let expected = peak * (-0.5_f64).exp(); // 0.6065 · peak
                assert!(
                    (v_off - expected).abs() < 0.1 * peak,
                    "plate {} 1σ offset value {} far from {} · peak ({})",
                    plate_idx,
                    v_off,
                    expected / peak,
                    expected
                );
            }
        }
    }

    /// Convolution produces a smooth field: max nearest-neighbour
    /// difference is bounded well below the unconvolved step (0.8).
    #[test]
    fn convolution_smooth_everywhere() {
        let nx = 64;
        let ny = 64;
        let plates = build_test_plates(nx, ny, 42);
        let ctx = ctx_with_plates(nx, ny, 42, 0.0, &plates);
        let s = init_s_field(InitMode::Convolution { sigma: 2.0 }, &ctx);

        let data = s.data();
        let mut max_grad = 0.0_f64;
        for j in 0..ny {
            for i in 0..nx {
                let here = data[j * nx + i];
                let right = data[j * nx + ((i + 1) % nx)];
                let down = data[((j + 1) % ny) * nx + i];
                max_grad = max_grad.max((here - right).abs()).max((here - down).abs());
            }
        }
        // Unconvolved step at a boundary is 0.8 (1.0 - 0.2). With
        // σ=2 cells the Gaussian envelope at the step's centre
        // produces a finite-difference ≤ ~0.16. Keep a generous
        // bound at 0.4.
        assert!(
            max_grad < 0.4,
            "convolution produced a too-steep gradient: {}",
            max_grad
        );
    }

    /// Determinism: same `(seed, nx, ny, mode)` → byte-identical S̃
    /// (relevant for the bit-determinism contract).
    #[test]
    fn determinism_same_seed_same_output() {
        let nx = 32;
        let ny = 32;
        let plates_a = build_test_plates(nx, ny, 42);
        let plates_b = build_test_plates(nx, ny, 42);

        for &mode in &[
            InitMode::Checkerboard,
            InitMode::Uniform { boundary_smoothing_width: 1.0 },
            InitMode::Gaussian { sigma_continental: 4.0, sigma_oceanic: 4.0 },
            InitMode::Convolution { sigma: 1.5 },
            InitMode::RadialProfile {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: ProfileShape::Smoothstep,
            },
            InitMode::RadialProfileWithFBM {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: ProfileShape::Smoothstep,
                fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
                fbm_octaves: FBM_OCTAVES_DEFAULT,
                fbm_persistence: FBM_PERSISTENCE_DEFAULT,
                fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
                fbm_scale: FBM_SCALE_DEFAULT,
                fbm_seed: FBM_SEED_DEFAULT,
                // Step 13.5 — disabled (default) so this
                // determinism check exercises the Step 13 path.
                apply_fbm_to_oceanic: false,
                fbm_amplitude_oceanic: FBM_AMPLITUDE_OCEANIC_DEFAULT,
                fbm_scale_oceanic: None,
                fbm_seed_oceanic: None,
            },
            // Step 13.5 — second variant with the oceanic FBM
            // path enabled, so the determinism contract also
            // covers the new code path.
            InitMode::RadialProfileWithFBM {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: ProfileShape::Smoothstep,
                fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
                fbm_octaves: FBM_OCTAVES_DEFAULT,
                fbm_persistence: FBM_PERSISTENCE_DEFAULT,
                fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
                fbm_scale: FBM_SCALE_DEFAULT,
                fbm_seed: FBM_SEED_DEFAULT,
                apply_fbm_to_oceanic: true,
                fbm_amplitude_oceanic: FBM_AMPLITUDE_OCEANIC_DEFAULT,
                fbm_scale_oceanic: None,
                fbm_seed_oceanic: None,
            },
        ] {
            let s_a = init_s_field(mode, &ctx_with_plates(nx, ny, 42, 0.2, &plates_a));
            let s_b = init_s_field(mode, &ctx_with_plates(nx, ny, 42, 0.2, &plates_b));
            assert_eq!(s_a.data(), s_b.data(), "non-deterministic for {:?}", mode);
        }
    }

    #[test]
    fn default_is_uniform() {
        match InitMode::default() {
            InitMode::Uniform { boundary_smoothing_width } => {
                assert!((boundary_smoothing_width - 1.0).abs() < 1e-12);
            }
            other => panic!("expected Uniform default, got {:?}", other),
        }
    }
}
