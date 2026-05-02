//! Step 8.5a Phase 0 sanity test — the `step0_quiescent` "synthetic"
//! benchmark case is a deterministic physics capture with
//! `YieldingConfig::Disabled`, `BasalDragConfig::Disabled`,
//! `BoundaryConfig::Disabled`, `SlabPullConfig::Disabled`, and
//! `MantleConfig::Disabled`. This test pins two invariants:
//!
//!  1. Two independent captures of the same config are byte-identical
//!     on disk (D9 determinism).
//!  2. The captured snapshot reconstitutes exactly into the Newton
//!     iter-0 state and round-trips through `LinearStokesSnapshot::{save, load}`.
//!
//! If either invariant fails, the "synthetic = exact" claim in the
//! Phase 0 report no longer holds and the benchmark reference loses
//! its determinism anchor.

use std::path::PathBuf;

use ymir_core::tectonics_v2::basal_drag::BasalDragConfig;
use ymir_core::tectonics_v2::boundaries::BoundaryConfig;
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, HarnessCaptureSpec, NonlinearChoice,
};
use ymir_core::tectonics_v2::mantle::MantleConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::snapshot::LinearStokesSnapshot;

fn build_step0_config(path: PathBuf) -> BaselineConfig {
    let scales = Scales::default();
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: 42,
        grid_nx: 64,
        grid_ny: 64,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: 1,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_step0_synthetic_parity_scratch"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Disabled,
        basal_drag: BasalDragConfig::Disabled,
        boundary: BoundaryConfig::Disabled,
        boundary_layout_name: String::new(),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: Some(HarnessCaptureSpec {
            at_step: 0,
            path,
            case_label: "step0_quiescent".into(),
        }),
        linear_solver: Default::default(),
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        continuation: None,
    }
}

#[test]
fn two_captures_produce_byte_identical_snapshots() {
    let a = std::env::temp_dir().join("ymir_step0_capture_a.bin");
    let b = std::env::temp_dir().join("ymir_step0_capture_b.bin");
    std::fs::remove_file(&a).ok();
    std::fs::remove_file(&b).ok();
    let _ = run_baseline(&build_step0_config(a.clone()));
    let _ = run_baseline(&build_step0_config(b.clone()));
    let bytes_a = std::fs::read(&a).expect("capture a");
    let bytes_b = std::fs::read(&b).expect("capture b");
    assert_eq!(
        bytes_a.len(),
        bytes_b.len(),
        "two captures produced different file sizes ({} vs {})",
        bytes_a.len(),
        bytes_b.len()
    );
    assert_eq!(
        bytes_a, bytes_b,
        "two captures of identical config produced different bytes (D9 determinism violated)"
    );
    std::fs::remove_file(&a).ok();
    std::fs::remove_file(&b).ok();
}

#[test]
fn snapshot_roundtrip_preserves_every_byte() {
    // The "synthetic is exact" claim rests on two pieces:
    //   (a) The capture is deterministic (covered by the byte-
    //       identity test above).
    //   (b) Saving and re-loading the bincode blob yields a struct
    //       with identical field values — no lossy compression, no
    //       sentinel rounding.
    // Together they pin the invariant that the benchmark `Vec<f64>`
    // values equal the exact Newton iter-0 state the harness saw.
    //
    // Note on "uniform η at iter 0": the harness runs an Ar-
    // continuation ramp (`run_continuation`) BEFORE the main step
    // loop, which mutates `v` away from zero. The captured η is
    // therefore heterogeneous by construction — the "synthetic"
    // label means "this is the exact post-continuation initial
    // state of the Step 0 config", not "uniform viscosity".
    let p = std::env::temp_dir().join("ymir_step0_roundtrip.bin");
    std::fs::remove_file(&p).ok();
    let _ = run_baseline(&build_step0_config(p.clone()));
    let snap = LinearStokesSnapshot::load(&p).expect("load snapshot");
    assert_eq!(snap.nx, 64);
    assert_eq!(snap.ny, 64);
    assert_eq!(snap.dx, 1.0 / 64.0);
    assert_eq!(snap.eta_center.len(), 64 * 64);
    assert!(snap.has_tangent(), "step0 must record Newton tangent fields");

    // Gauge-fixed RHS (zero mean per component). If this fails, the
    // projection inside NewtonSolver is broken — not a Phase 0 bug.
    let mean_vx: f64 = snap.rhs_vx.iter().sum::<f64>() / snap.rhs_vx.len() as f64;
    let mean_vy: f64 = snap.rhs_vy.iter().sum::<f64>() / snap.rhs_vy.len() as f64;
    assert!(mean_vx.abs() < 1e-10, "rhs_vx mean = {:.3e}", mean_vx);
    assert!(mean_vy.abs() < 1e-10, "rhs_vy mean = {:.3e}", mean_vy);

    // Save-load round trip — every f64 survives bincode byte-for-byte.
    let p2 = std::env::temp_dir().join("ymir_step0_roundtrip_copy.bin");
    snap.save(&p2).unwrap();
    let snap2 = LinearStokesSnapshot::load(&p2).unwrap();
    assert_eq!(snap.eta_center, snap2.eta_center);
    assert_eq!(snap.tangent_c_center, snap2.tangent_c_center);
    assert_eq!(snap.tangent_exx_center, snap2.tangent_exx_center);
    assert_eq!(snap.tangent_eyy_center, snap2.tangent_eyy_center);
    assert_eq!(snap.tangent_exy_corner, snap2.tangent_exy_corner);
    assert_eq!(snap.drag_diag, snap2.drag_diag);
    assert_eq!(snap.diag_vx, snap2.diag_vx);
    assert_eq!(snap.diag_vy, snap2.diag_vy);
    assert_eq!(snap.rhs_vx, snap2.rhs_vx);
    assert_eq!(snap.rhs_vy, snap2.rhs_vy);
    assert_eq!(snap.format_version, snap2.format_version);
    assert_eq!(snap.case_label, snap2.case_label);
    std::fs::remove_file(&p).ok();
    std::fs::remove_file(&p2).ok();
}
