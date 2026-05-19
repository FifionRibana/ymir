//! Slab-mass source term `Q̃_sub_conv`.
//!
//! This is the "accumulation" counterpart of the Step 5 `Q_sub`
//! drain. Physically they both measure convergence at a subducting
//! cell, but they feed different pipelines:
//!
//! ```text
//!   Q_sub      = k_sub       · max(0, -div v)    drains S̃ (Steps 5-6)
//!   Q_sub_conv = k_slab_accum · max(0, -div v)   sources m̃ (Step 7)
//! ```
//!
//! The two rates are kept **independent** — see D3 of the Step 7
//! spec. Coupling them strictly (e.g. `k_slab_accum =
//! (1 − arc − coll_v − rift_v) · k_sub`) would force `Sp` out of
//! the §4.8 target band `[0.5, 3.0]` and entangle the slab
//! calibration with the recycling fractions. Baseline
//! `k_slab_accum = 1.0` absorbs the normalisation into `Sp`; no
//! conservation of `S̃` is violated because `m̃` is not tracked in
//! the mass budget.
//!
//! The source fires only on cells that are simultaneously
//! **oceanic** (`plate_type == Oceanic`) and **subducting**
//! (`boundary_flag ∈ {Subduction, OceanicSubduction}`). Other
//! cells receive zero. This reuses the exact same predicate as
//! `Q_sub` in [`crate::tectonics_v2::boundaries::source_sink`], so the
//! two rates act in lock-step on the flag field.

use crate::tectonics_v2::boundaries::boundary_flag::BoundaryFlagField;
use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::boundaries::source_sink::{convergent_component, div_v_cell};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

/// Knobs for [`compute_q_sub_conv`].
#[derive(Clone, Copy, Debug)]
pub struct AccumulationConfig {
    /// Slab-mass accumulation rate. Defaults to 1.0 (D3).
    pub k_slab_accum: f64,
}

impl Default for AccumulationConfig {
    fn default() -> Self {
        Self { k_slab_accum: super::K_SLAB_ACCUM_DEFAULT }
    }
}

/// Compute `Q̃_sub_conv(i, j)` in place.
///
/// The caller supplies a scratch `div_v` buffer to avoid per-step
/// allocations in the time loop. `out` is overwritten (not
/// accumulated) with:
///
/// ```text
///   Q_sub_conv[i,j] = k_slab_accum · max(0, -div v)   on (Oceanic ∧ subducting) cells
///                   = 0                                elsewhere
/// ```
pub fn compute_q_sub_conv(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    vx: &[f64],
    vy: &[f64],
    plate_type: &PlateTypeField,
    boundary_flag: &BoundaryFlagField,
    div_scratch: &mut Field2D,
    conv_scratch: &mut Field2D,
    out: &mut Field2D,
    config: &AccumulationConfig,
) {
    debug_assert_eq!(plate_type.nx(), nx);
    debug_assert_eq!(plate_type.ny(), ny);
    debug_assert_eq!(boundary_flag.nx(), nx);
    debug_assert_eq!(boundary_flag.ny(), ny);
    debug_assert_eq!(out.nx(), nx);
    debug_assert_eq!(out.ny(), ny);

    // Pass 1: div(v) in `div_scratch`.
    div_v_cell(nx, ny, dx, dy, idx_x, idx_y, vx, vy, div_scratch);
    // Pass 2: max(0, -div) in `conv_scratch`.
    convergent_component(div_scratch, conv_scratch);
    // Pass 3: gate by (Oceanic ∧ subducting), scale by k.
    let k = config.k_slab_accum;
    for j in 0..ny {
        for i in 0..nx {
            let pt = plate_type.get(i, j);
            let bf = boundary_flag.get(i, j);
            let fires = matches!(pt, PlateType::Oceanic) && bf.is_subduction();
            let value = if fires { k * conv_scratch.get(i, j) } else { 0.0 };
            out.set(i, j, value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::crate::tectonics_v2::boundaries::boundary_flag::BoundaryFlag;
    use super::*;

    fn make_env(nx: usize, ny: usize) -> (PeriodicIndex, PeriodicIndex, f64) {
        (PeriodicIndex::new(nx), PeriodicIndex::new(ny), 1.0 / nx as f64)
    }

    /// Oceanic non-boundary cell: no accumulation, even if div v < 0.
    #[test]
    fn non_boundary_cell_is_zero() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y, dx) = make_env(nx, ny);

        // Any divergence works; pick something obvious.
        let vx = vec![0.0; nx * ny];
        let mut vy = vec![0.0; nx * ny];
        for i in 0..nx * ny {
            vy[i] = -0.5;
        }

        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let boundary_flag = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);

        let mut div_scratch = Field2D::new(nx, ny);
        let mut conv_scratch = Field2D::new(nx, ny);
        let mut out = Field2D::new(nx, ny);
        compute_q_sub_conv(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &plate_type,
            &boundary_flag,
            &mut div_scratch,
            &mut conv_scratch,
            &mut out,
            &AccumulationConfig::default(),
        );
        for v in out.data().iter() {
            assert_eq!(*v, 0.0);
        }
    }

    /// Continental subduction flag (shouldn't happen physically,
    /// but we test that the gate keeps it out): zero output.
    #[test]
    fn continental_subducting_is_zero() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y, dx) = make_env(nx, ny);

        // Set up convergent vx field: vx[i+1] < vx[i] ⇒ div < 0.
        let mut vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                vx[j * nx + i] = -(i as f64) * 0.1;
            }
        }

        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        let boundary_flag = BoundaryFlagField::filled(nx, ny, BoundaryFlag::Subduction);

        let mut div_scratch = Field2D::new(nx, ny);
        let mut conv_scratch = Field2D::new(nx, ny);
        let mut out = Field2D::new(nx, ny);
        compute_q_sub_conv(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &plate_type,
            &boundary_flag,
            &mut div_scratch,
            &mut conv_scratch,
            &mut out,
            &AccumulationConfig::default(),
        );
        for v in out.data().iter() {
            assert_eq!(*v, 0.0);
        }
    }

    /// Oceanic + Subduction + zero flow: div(v) = 0 ⇒
    /// max(0, -div) = 0 ⇒ Q = 0. On a periodic tore, the global
    /// integral of div(v) is always 0, so a "uniformly
    /// divergent" field is impossible; we probe the convergent
    /// gate with the trivial zero-field here and rely on
    /// `linear_in_k_slab_accum` for the positive case.
    #[test]
    fn zero_flow_gives_zero_on_oceanic_subducting_cell() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y, dx) = make_env(nx, ny);
        let vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let boundary_flag = BoundaryFlagField::filled(nx, ny, BoundaryFlag::OceanicSubduction);

        let mut div_scratch = Field2D::new(nx, ny);
        let mut conv_scratch = Field2D::new(nx, ny);
        let mut out = Field2D::new(nx, ny);
        compute_q_sub_conv(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &plate_type,
            &boundary_flag,
            &mut div_scratch,
            &mut conv_scratch,
            &mut out,
            &AccumulationConfig::default(),
        );
        for v in out.data().iter() {
            assert_eq!(*v, 0.0);
        }
    }

    /// Oceanic + OceanicSubduction + convergent flow: Q > 0 and
    /// scales linearly with `k_slab_accum`.
    #[test]
    fn linear_in_k_slab_accum() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y, dx) = make_env(nx, ny);

        let mut vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                vx[j * nx + i] = -(i as f64) * 0.1;
            }
        }
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let boundary_flag = BoundaryFlagField::filled(nx, ny, BoundaryFlag::OceanicSubduction);

        let mut div_scratch = Field2D::new(nx, ny);
        let mut conv_scratch = Field2D::new(nx, ny);

        let mut out_1 = Field2D::new(nx, ny);
        compute_q_sub_conv(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &plate_type,
            &boundary_flag,
            &mut div_scratch,
            &mut conv_scratch,
            &mut out_1,
            &AccumulationConfig { k_slab_accum: 1.0 },
        );

        let mut out_3 = Field2D::new(nx, ny);
        compute_q_sub_conv(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &plate_type,
            &boundary_flag,
            &mut div_scratch,
            &mut conv_scratch,
            &mut out_3,
            &AccumulationConfig { k_slab_accum: 3.0 },
        );

        // Every cell should scale exactly 3× — and at least one
        // cell must be non-zero (the divergence is uniform negative).
        let mut any_nonzero = false;
        for k in 0..nx * ny {
            let a = out_1.data()[k];
            let b = out_3.data()[k];
            assert!((b - 3.0 * a).abs() < 1e-14);
            if a > 0.0 {
                any_nonzero = true;
            }
        }
        assert!(any_nonzero, "expected at least one convergent cell");
    }
}
