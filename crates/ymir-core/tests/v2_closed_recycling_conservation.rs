//! Step 6 — conservation of the Closed-mode recycling pipeline.
//!
//! Runs a 64² × 50-step baseline in Closed mode with `mantle_loss =
//! 0` and verifies that the Step 6 mass-conservation residual
//! `|Δmass_obs + mantle_loss + buffer_fill + pending − clamp_flux|
//! / initial_mass` lands below `1e-6`. This confirms the 5-way
//! balance across the pipeline (Q through the grid, buffer fill,
//! pending accumulators, clamp artificial flux, and mantle loss).

use std::path::PathBuf;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{horizontal_oceanic_strip, BoundaryRates};
use ymir_core::tectonics_v2::boundaries::boundary_flag::{BoundaryConfig, RecyclingModeInit};
use ymir_core::tectonics_v2::boundaries::crust_geometry::CrustGeometry;
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;

fn run_closed_mini(mantle_loss_fraction: f64) -> (f64, f64, f64, f64, f64) {
    let nx = 64;
    let ny = 64;
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let layout = horizontal_oceanic_strip(nx, ny);
    // Static layout carried into Closed mode: we use the Step 5
    // horizontal_oceanic_strip so subduction + rift cells exist.
    let geometry = CrustGeometry::from_static(
        layout.plate_types, layout.flags, layout.name,
    );
    let rates = BoundaryRates::baseline_uncalibrated().with_k_spread(0.05);
    // Closed mode config: spread carries whatever is left after
    // immediate + loss. Adjust `spread_fraction` to match.
    let arc = 0.15;
    let coll = 0.03;
    let rift = 0.02;
    let spread = 1.0 - arc - coll - rift - mantle_loss_fraction;
    let recycling_config = RecyclingConfig {
        arc_fraction: arc,
        coll_v_fraction: coll,
        rift_v_fraction: rift,
        spread_fraction: spread,
        mantle_loss_fraction,
        mantle_delay_steps: 20,
    };
    recycling_config.validate().unwrap();
    let boundary = BoundaryConfig::Enabled {
        geometry: std::sync::Arc::new(geometry),
        rates,
        recycling_mode: RecyclingModeInit::Closed(recycling_config),
    };
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
        output_dir: PathBuf::from("target/v2_closed_recycling_scratch"),
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
        boundary_layout_name: "horizontal_oceanic_strip_closed".into(),
        slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
    };
    let r = run_baseline(&cfg);
    let na = r.metrics.newton.as_ref().expect("newton aggregate");
    let residual = na.mass_conservation_residual.expect("mass_conservation_residual");
    let mantle_loss_int = na.mantle_loss_integral.unwrap_or(0.0);
    let m_sub_total = na.m_sub_total.unwrap_or(0.0);
    let buffer_fill = na.recycling_buffer_fill_final.unwrap_or(0.0);
    let pending = na.immediate_pending_final.unwrap_or(0.0);
    (residual, mantle_loss_int, m_sub_total, buffer_fill, pending)
}

#[test]
fn conservation_with_no_mantle_loss_is_below_1ppm() {
    let (residual, _, m_sub, buf, pend) = run_closed_mini(0.0);
    println!(
        "no-loss: residual={:.3e} m_sub_total={:.3e} buffer_fill={:.3e} pending={:.3e}",
        residual, m_sub, buf, pend,
    );
    assert!(
        residual < 1.0e-6,
        "mass_conservation_residual = {} ≥ 1e-6 with mantle_loss=0",
        residual,
    );
}

#[test]
fn recycling_with_5pct_loss_accumulates_expected_mass() {
    let (residual, mantle_loss_int, m_sub, _, _) = run_closed_mini(0.05);
    println!(
        "5pct-loss: residual={:.3e} mantle_loss_int={:.3e} m_sub_total={:.3e} ratio={:.3e}",
        residual, mantle_loss_int, m_sub, mantle_loss_int / m_sub.max(1e-30),
    );
    // Conservation residual still < 1e-6 with the 5-way balance.
    assert!(
        residual < 1.0e-6,
        "residual {} should stay below 1e-6 even with mantle loss",
        residual,
    );
    // Observed loss / total subducted ≈ 0.05 ± 1% (exact to f64
    // rounding at this precision — the ratio is `mantle_loss_fraction`
    // applied linearly to M_sub each step).
    if m_sub > 0.0 {
        let ratio = mantle_loss_int / m_sub;
        assert!(
            (ratio - 0.05).abs() < 1e-12,
            "observed loss ratio = {} (expected exactly 0.05)",
            ratio,
        );
    }
}
