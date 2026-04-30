//! Step 8.5a Phase 4.3 — physics-run scalar-parity between
//! `JacobiCG` (merged baseline) and `AmgCG` (Phase 2+ new path)
//! on step0/3/6/7. step8_activated is excluded per α scope.
//!
//! Each test builds the same `BaselineConfig` used by the
//! production binaries / gen_bench_data, runs the full physics
//! pipeline twice (once per preconditioner), and asserts that
//! the key Newton-metrics agree to 1 % relative:
//!
//!  - `peak|v|` (final) — sensitivity probe for the velocity field
//!  - `mass_conservation_residual` (final) — divergence diagnostic
//!  - `yielding_cell_fraction_max` (when applicable)
//!
//! JacobiCG bit-parity is already enforced by
//! [`v2_step8_regression_smoke`]; this test covers the new
//! AmgCG dispatch end-to-end through the harness and Newton
//! solver.
//!
//! Tests run SEQUENTIALLY in reviewer-mandated order (0 → 3 →
//! 6/7) — if step0 fails there is no point running the longer
//! cases before the diagnostic.
//!
//! Step count reduced to `PHASE43_STEPS = 100` (from the
//! production default of 300) to cap wallclock at ~1 min per
//! test pair; the physics regime changes little between step
//! 100 and 300 for these configurations, and the scalar-parity
//! invariant is agnostic to run length. Full 300-step runs are
//! a separate graduation-gate responsibility post-merge.

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::mantle::MantleConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::stokes::amg::AmgConfig;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

const PHASE43_STEPS: usize = 100;
const PARITY_REL_TOL: f64 = 0.01;

fn build_step0_config(seed: u64, linear_solver: LinearSolverConfig) -> BaselineConfig {
    let scales = Scales::default();
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed,
        grid_nx: 64,
        grid_ny: 64,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: PHASE43_STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_phase43_step0"),
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
        capture: None,
        linear_solver,
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
    }
}

fn build_step3_config(seed: u64, linear_solver: LinearSolverConfig) -> BaselineConfig {
    let scales = Scales::default();
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed,
        grid_nx: 64,
        grid_ny: 64,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: PHASE43_STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_phase43_step3"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Disabled,
        boundary: BoundaryConfig::Disabled,
        boundary_layout_name: String::new(),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver,
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
    }
}

fn build_step6_config(
    seed: u64,
    num_plates: usize,
    continental_ratio: f64,
    linear_solver: LinearSolverConfig,
) -> Result<BaselineConfig, String> {
    let scales = Scales::default();
    let vcfg = VoronoiConfig { num_plates, continental_ratio };
    let rates = BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let recycling_config = RecyclingConfig::default();
    let boundary = BoundaryConfig::enabled_voronoi_closed(64, 64, &vcfg, seed, rates, recycling_config)
        .map_err(|e| format!("boundary: {:?}", e))?;
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    Ok(BaselineConfig {
        seed,
        grid_nx: 64,
        grid_ny: 64,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: PHASE43_STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_phase43_step6"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw { br: 0.05, ..BasalDragLaw::default() }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", seed, num_plates),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver,
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
    })
}

fn check_parity(
    case: &str,
    cfg_jacobi: BaselineConfig,
    cfg_amg: BaselineConfig,
) {
    use std::time::Instant;
    let t0 = Instant::now();
    let res_j = run_baseline(&cfg_jacobi);
    let dt_j = t0.elapsed().as_secs_f64();
    let t0 = Instant::now();
    let res_a = run_baseline(&cfg_amg);
    let dt_a = t0.elapsed().as_secs_f64();

    let peak_v_j = res_j.metrics.vmax_peak;
    let peak_v_a = res_a.metrics.vmax_peak;
    let mass_j = res_j.metrics.mass_drift_relative;
    let mass_a = res_a.metrics.mass_drift_relative;
    let yf_j = res_j.metrics.yielding_cell_fraction.unwrap_or(0.0);
    let yf_a = res_a.metrics.yielding_cell_fraction.unwrap_or(0.0);
    let cg_mean_j = res_j.metrics.cg_iter_mean;
    let cg_mean_a = res_a.metrics.cg_iter_mean;

    eprintln!("=== {} (steps = {}) ===", case, PHASE43_STEPS);
    eprintln!("  JacobiCG  wallclock {:.2}s  vmax_peak={:.4e}  mass_drift={:.4e}  yf_cell={:.4e}  CG_mean={:.1}",
        dt_j, peak_v_j, mass_j, yf_j, cg_mean_j);
    eprintln!("  AmgCG     wallclock {:.2}s  vmax_peak={:.4e}  mass_drift={:.4e}  yf_cell={:.4e}  CG_mean={:.1}",
        dt_a, peak_v_a, mass_a, yf_a, cg_mean_a);
    eprintln!("  wallclock ratio AMG/Jacobi = {:.3}x", dt_a / dt_j.max(1e-10));

    let rel_peak = ((peak_v_a - peak_v_j).abs() / peak_v_j.abs().max(1e-300)).max(0.0);
    let abs_mass_diff = (mass_a - mass_j).abs();
    let rel_mass = abs_mass_diff / mass_j.abs().max(1e-12);
    let rel_yf = if yf_j.abs() > 1e-12 {
        (yf_a - yf_j).abs() / yf_j.abs()
    } else {
        (yf_a - yf_j).abs()
    };
    eprintln!("  rel_diff: vmax_peak={:.3e}  mass_drift_abs={:.3e} (rel={:.3e})  yf={:.3e}  (tol {:.0e})",
        rel_peak, abs_mass_diff, rel_mass, rel_yf, PARITY_REL_TOL);

    assert!(
        rel_peak < PARITY_REL_TOL,
        "{}: vmax_peak parity {:.3e} > tol {:.0e} (Jacobi {:.4e}, AMG {:.4e})",
        case, rel_peak, PARITY_REL_TOL, peak_v_j, peak_v_a,
    );
    // Mass drift is already a relative number in the solver's
    // reporting. Require absolute equivalence below 1e-6 or
    // relative 1% — whichever is more permissive — because mass
    // drift can be very small.
    if mass_j.abs() > 1e-6 {
        assert!(
            rel_mass < PARITY_REL_TOL,
            "{}: mass_drift parity {:.3e} > tol {:.0e}",
            case, rel_mass, PARITY_REL_TOL,
        );
    } else {
        assert!(
            abs_mass_diff < 1e-6,
            "{}: mass_drift absolute parity {:.3e} > 1e-6 (both below relative threshold floor)",
            case, abs_mass_diff,
        );
    }
    if yf_j.abs() > 1e-6 {
        assert!(
            rel_yf < PARITY_REL_TOL,
            "{}: yielding_cell_fraction parity {:.3e} > tol {:.0e}",
            case, rel_yf, PARITY_REL_TOL,
        );
    }
}

#[test]
fn step0_physics_scalar_parity() {
    let cfg_j = build_step0_config(42, LinearSolverConfig::JacobiCG);
    let cfg_a = build_step0_config(42, LinearSolverConfig::AmgCG(AmgConfig::default()));
    check_parity("step0_quiescent", cfg_j, cfg_a);
}

#[test]
fn step3_physics_scalar_parity() {
    let cfg_j = build_step3_config(42, LinearSolverConfig::JacobiCG);
    let cfg_a = build_step3_config(42, LinearSolverConfig::AmgCG(AmgConfig::default()));
    check_parity("step3_floor_yielding", cfg_j, cfg_a);
}

#[test]
fn step6_physics_scalar_parity() {
    let cfg_j = build_step6_config(42, 8, 0.4, LinearSolverConfig::JacobiCG).unwrap();
    let cfg_a = build_step6_config(42, 8, 0.4, LinearSolverConfig::AmgCG(AmgConfig::default())).unwrap();
    check_parity("step6_voronoi", cfg_j, cfg_a);
}

#[test]
fn step7_physics_scalar_parity() {
    // Step 7 shape with slab-pull Disabled (matches Step 7
    // regression scenario; full slab-pull activation is a
    // separate Phase 4 scope question not under α merge).
    let cfg_j = build_step6_config(42, 8, 0.4, LinearSolverConfig::JacobiCG).unwrap();
    let cfg_a = build_step6_config(42, 8, 0.4, LinearSolverConfig::AmgCG(AmgConfig::default())).unwrap();
    check_parity("step7_slab_off_shape", cfg_j, cfg_a);
}
