//! Diagnostic-only test (not an acceptance): isolate whether the
//! Step 8 runaway peak|v| = O(10^34) observed in the baseline
//! smoke is mantle-alone or slab+mantle-combined.
//!
//! Runs a short (20-step) Step 6 physics setup (GPE + yielding +
//! basal drag + Voronoï + Closed recycling) with
//! `SlabPullConfig::Disabled` and `MantleConfig::Enabled` at
//! baseline parameters. Prints key metrics so we can read them
//! in the test output.

use std::path::PathBuf;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, ForceKind, NonlinearChoice, build_force, run_baseline,
};
use ymir_core::tectonics_v2::mantle::MantleConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

fn build_cfg(
    slab: SlabPullConfig,
    mantle: MantleConfig,
    steps: usize,
    total_time: f64,
) -> BaselineConfig {
    let nx = 64;
    let ny = 64;
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let rates = BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(nx, ny, &vcfg, 42, rates, RecyclingConfig::default()).unwrap();
    BaselineConfig {
        seed: 42, grid_nx: nx, grid_ny: ny,
        domain_lx: 1.0, domain_ly: 1.0, steps,
        cfl_factor: 0.3, total_time_nondim: total_time,
        preset, nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(), picard_cfg: Default::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from("target/mantle_diag_scratch"),
        force: build_force(ForceKind::Gpe, &scales, 10.0, 1.0),
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0, s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw { br: 0.05, ..BasalDragLaw::default() }),
        boundary, boundary_layout_name: "voronoi_seed42_n8".into(),
        slab_pull: slab,
        mantle,
    }
}

fn print_metrics(label: &str, r: &ymir_core::tectonics_v2::diagnostics::harness::BaselineResult) {
    let m = &r.metrics;
    let na = m.newton.as_ref().unwrap();
    eprintln!("--- {} ---", label);
    eprintln!("peak|v|:             {:.3e}", m.vmax_peak);
    eprintln!("peak|v_solved|:      {:.3e}", na.peak_v_solved_mantle_run.unwrap_or(0.0));
    eprintln!("peak|v_mantle|:      {:.3e}", na.peak_v_mantle_pattern.unwrap_or(0.0));
    eprintln!("alignment:           {:.3}", na.v_solved_to_v_mantle_alignment.unwrap_or(0.0));
    eprintln!("peak|f_mantle|:      {:.3e}", na.peak_f_mantle_run.unwrap_or(0.0));
    eprintln!("peak|f_slab|:        {:.3e}", na.peak_f_slab_run.unwrap_or(0.0));
    eprintln!("peak|f_GPE|:         {:.3e}", na.peak_f_gpe_run.unwrap_or(0.0));
    eprintln!("yielding_cell_frac:  {:.3e}", na.yielding_cell_fraction_max.unwrap_or(0.0));
    eprintln!("ε̇_II/floor:          {:.3e}", na.epsilon_ii_max_to_floor_ratio.unwrap_or(0.0));
    eprintln!("CG iters mean/max:   {:.1} / {}", m.cg_iter_mean, m.cg_iter_max);
    eprintln!("Newton outer mean:   {:.1}", na.outer_iters_mean());
    eprintln!("mass_cons_residual:  {:.3e}", na.mass_conservation_residual.unwrap_or(0.0));
    eprintln!("div_v_mantle_max:    {:.3e}", na.div_v_mantle_max.unwrap_or(0.0));
    eprintln!();
}

#[test]
#[ignore = "diagnostic-only, run with --include-ignored"]
fn mantle_only_20_steps_reports_metrics() {
    let cfg = build_cfg(
        SlabPullConfig::Disabled,
        MantleConfig::Enabled { mf: 1.0, coupling: 1.0, num_modes: 6, seed: 42, evolution_rate: 0.0 },
        20, 0.4,
    );
    let r = run_baseline(&cfg);
    print_metrics("mantle-only (slab Disabled)", &r);
}

/// Diagnostic: does slab+mantle produce runaway at short
/// timescales? Expected per the Step 8 smoke: peak|v| grows
/// unboundedly.
#[test]
#[ignore = "diagnostic-only, run with --include-ignored"]
fn mantle_plus_slab_runaway_characterisation() {
    for steps in [5, 10, 15, 20] {
        let cfg = build_cfg(
            SlabPullConfig::Enabled { sp: 1.5, tau_slab: 0.5, k_slab_accum: 1.0, epsilon: 1.0e-6 },
            MantleConfig::Enabled { mf: 1.0, coupling: 1.0, num_modes: 6, seed: 42, evolution_rate: 0.0 },
            steps, steps as f64 * 0.02,
        );
        let r = run_baseline(&cfg);
        print_metrics(&format!("mantle+slab @ {} steps", steps), &r);
    }
}
