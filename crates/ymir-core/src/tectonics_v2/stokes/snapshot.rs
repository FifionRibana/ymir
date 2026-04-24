//! Linear-solve state snapshot for benchmark replay (Step 8.5a Phase 0).
//!
//! Captures the exact inputs to the preconditioned CG at Newton outer
//! iter 0 of a single physics timestep:
//! - the Picard block (`eta_center`, optional `drag_diag`)
//! - the optional Newton tangent context (four strain-rate fields + `c`)
//! - the Jacobi diagonal `(diag_vx, diag_vy)` pre-computed for the precond
//! - the gauge-fixed RHS `(rhs_vx, rhs_vy)` and initial guess `(x0_vx, x0_vy)`
//! - the CG hyperparameters (`tol`, `max_iter`, `diag_floor`)
//!
//! Benchmark replay rebuilds the matvec as
//! `apply_momentum(eta_center, drag_diag) + apply_tangent(ctx)` so a
//! Jacobi-CG run on the snapshot reproduces the exact iter count the
//! physics run saw at that Newton outer iter.
//!
//! Storage format: `bincode` with a `format_version` header; bump the
//! constant on any breaking schema change.

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::tectonics_v2::field::Field2D;

/// Current snapshot schema version. Phase 0 ships v1. Increment on any
/// breaking change (field removal, semantic change, layout reshape).
pub const SNAPSHOT_FORMAT_VERSION: u32 = 1;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LinearStokesSnapshot {
    pub format_version: u32,
    pub case_label: String,
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    /// Picard-block viscosity at cell centres, frozen for this solve.
    pub eta_center: Vec<f64>,
    /// Newton tangent scalar prefactor `c_cc = η'(ε̇_II_cc) / (ε̇_II_cc + ε̇_min)`
    /// at cell centres. `None` for Picard-only snapshots (step0 /
    /// Poisson synthetic cases).
    pub tangent_c_center: Option<Vec<f64>>,
    pub tangent_exx_center: Option<Vec<f64>>,
    pub tangent_eyy_center: Option<Vec<f64>>,
    pub tangent_exy_corner: Option<Vec<f64>>,
    /// Optional basal-drag diagonal (`Br · S̃^exp`, cell-centered).
    pub drag_diag: Option<Vec<f64>>,
    /// Preconditioner diagonal per velocity component (pre-computed).
    pub diag_vx: Vec<f64>,
    pub diag_vy: Vec<f64>,
    /// Gauge-fixed RHS packed per component.
    pub rhs_vx: Vec<f64>,
    pub rhs_vy: Vec<f64>,
    /// Initial guess packed per component (zeros at Newton outer iter 0).
    pub x0_vx: Vec<f64>,
    pub x0_vy: Vec<f64>,
    /// CG hyperparameters in force at capture time.
    pub tol: f64,
    pub max_iter: usize,
    pub diag_floor: f64,
}

impl LinearStokesSnapshot {
    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), String> {
        let p = path.as_ref();
        if let Some(parent) = p.parent() {
            std::fs::create_dir_all(parent)
                .map_err(|e| format!("create_dir_all {}: {}", parent.display(), e))?;
        }
        let file = File::create(p).map_err(|e| format!("create {}: {}", p.display(), e))?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, self).map_err(|e| format!("bincode serialize: {}", e))
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self, String> {
        let p = path.as_ref();
        let file = File::open(p).map_err(|e| format!("open {}: {}", p.display(), e))?;
        let reader = BufReader::new(file);
        let snap: Self =
            bincode::deserialize_from(reader).map_err(|e| format!("bincode deserialize: {}", e))?;
        if snap.format_version > SNAPSHOT_FORMAT_VERSION {
            return Err(format!(
                "snapshot {} was written with format v{} (this build supports up to v{})",
                p.display(),
                snap.format_version,
                SNAPSHOT_FORMAT_VERSION
            ));
        }
        Ok(snap)
    }

    pub fn n_cells(&self) -> usize {
        self.nx * self.ny
    }

    pub fn has_tangent(&self) -> bool {
        self.tangent_c_center.is_some()
    }
}

/// Materialise a `Field2D` from a raw row-major `Vec<f64>` plus grid dims.
pub fn field_from_vec(data: Vec<f64>, nx: usize, ny: usize) -> Field2D {
    assert_eq!(data.len(), nx * ny, "Field2D size mismatch");
    let mut f = Field2D::new(nx, ny);
    f.data_mut().copy_from_slice(&data);
    f
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_snapshot(nx: usize, ny: usize, with_tangent: bool) -> LinearStokesSnapshot {
        let n = nx * ny;
        LinearStokesSnapshot {
            format_version: SNAPSHOT_FORMAT_VERSION,
            case_label: "unit_test".into(),
            nx,
            ny,
            dx: 1.0 / nx as f64,
            dy: 1.0 / ny as f64,
            eta_center: (0..n).map(|k| 1.0 + 0.1 * (k as f64)).collect(),
            tangent_c_center: with_tangent.then(|| vec![0.5; n]),
            tangent_exx_center: with_tangent.then(|| vec![0.1; n]),
            tangent_eyy_center: with_tangent.then(|| vec![-0.1; n]),
            tangent_exy_corner: with_tangent.then(|| vec![0.05; n]),
            drag_diag: Some(vec![0.03; n]),
            diag_vx: vec![2.0; n],
            diag_vy: vec![2.0; n],
            rhs_vx: (0..n).map(|k| 0.1 * (k as f64).sin()).collect(),
            rhs_vy: (0..n).map(|k| 0.1 * (k as f64).cos()).collect(),
            x0_vx: vec![0.0; n],
            x0_vy: vec![0.0; n],
            tol: 1e-8,
            max_iter: 2000,
            diag_floor: 1e-20,
        }
    }

    #[test]
    fn roundtrip_with_tangent() {
        let tmp = std::env::temp_dir().join("ymir_snapshot_rt_tangent.bin");
        let snap = sample_snapshot(8, 8, true);
        snap.save(&tmp).unwrap();
        let back = LinearStokesSnapshot::load(&tmp).unwrap();
        assert_eq!(back.format_version, SNAPSHOT_FORMAT_VERSION);
        assert_eq!(back.nx, snap.nx);
        assert_eq!(back.eta_center, snap.eta_center);
        assert_eq!(back.tangent_c_center, snap.tangent_c_center);
        assert_eq!(back.rhs_vx, snap.rhs_vx);
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn roundtrip_picard_only() {
        let tmp = std::env::temp_dir().join("ymir_snapshot_rt_picard.bin");
        let snap = sample_snapshot(8, 8, false);
        snap.save(&tmp).unwrap();
        let back = LinearStokesSnapshot::load(&tmp).unwrap();
        assert!(!back.has_tangent());
        assert_eq!(back.eta_center, snap.eta_center);
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn field_from_vec_reconstructs() {
        let data: Vec<f64> = (0..16).map(|k| k as f64 * 0.25).collect();
        let f = field_from_vec(data.clone(), 4, 4);
        assert_eq!(f.data(), &data[..]);
    }
}
