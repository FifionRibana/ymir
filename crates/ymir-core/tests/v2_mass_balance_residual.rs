//! Step 5 — mass-balance residual on a non-balanced layout.
//!
//! Drives a short 64²·50-step run with `horizontal_oceanic_strip` +
//! GPE + yielding Enabled + basal drag Enabled + boundary Enabled
//! (uncalibrated `k_spread`, so the layout is deliberately
//! unbalanced and exercises both the physical flux ∫Q and the
//! clamp-induced artificial flux). The acceptance is that
//! `mass_balance_residual < 1%` even when the clamp has fired on a
//! non-negligible fraction of cells.

use std::path::PathBuf;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{
    horizontal_oceanic_strip, BoundaryRates,
};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;

fn run_step5_mini(k_sub: f64) -> (f64, f64) {
    let nx = 64;
    let ny = 64;
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let layout = horizontal_oceanic_strip(nx, ny);
    let layout_name = layout.name.to_string();
    let rates = BoundaryRates::baseline_uncalibrated()
        .with_k_sub(k_sub);
    let boundary = layout.into_config(rates);
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let cfg = BaselineConfig {
        seed: 42,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: 50,
        cfl_factor: 0.3,
        total_time_nondim: 1.0,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from("target/v2_mass_balance_residual_scratch"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: 0.05,
            ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: layout_name,
        slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
    };
    let r = run_baseline(&cfg);
    let na = r.metrics.newton.as_ref().expect("newton aggregate");
    let residual = na.mass_balance_residual.expect("mass_balance_residual populated");
    let clamp_mean = na.clamp_activation_fraction_mean.unwrap_or(0.0);
    (residual, clamp_mean)
}

#[test]
fn residual_below_one_percent_on_baseline_k_sub() {
    let (residual, clamp_mean) = run_step5_mini(0.5);
    println!("baseline: residual={:.3e} clamp_mean={:.3e}", residual, clamp_mean);
    assert!(
        residual < 0.01,
        "mass_balance_residual = {} ≥ 1% at baseline (should be tight)",
        residual,
    );
}

#[test]
fn residual_below_one_percent_with_high_k_sub_and_active_clamp() {
    // Raise k_sub to drive more cells into the clamp floor and
    // confirm the balance still holds: the artificial flux is
    // properly accounted in the residual.
    let (residual, clamp_mean) = run_step5_mini(1.5);
    println!("high k_sub: residual={:.3e} clamp_mean={:.3e}", residual, clamp_mean);
    // This test's whole point is that the residual doesn't blow up
    // even when clamp_mean grows. The 1% bound is deliberate.
    assert!(
        residual < 0.01,
        "mass_balance_residual = {} ≥ 1% with active clamp (clamp_mean={}) — flux accounting under-captures",
        residual, clamp_mean,
    );
}
