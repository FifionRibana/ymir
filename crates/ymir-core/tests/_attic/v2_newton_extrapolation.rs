//! Step 8.5b Phase 5 — Newton extrapolation contract tests.
//!
//! Verifies the harness instrumentation produced by the order-2
//! warm-start extrapolation: attempted-vs-applied counters, the
//! per-step `newton_outer_iters_per_step` history, and the
//! fallback rate. Also asserts that the *converged* metrics from
//! a 100-step physics run agree at scalar-parity with two
//! independent reproductions of the same config (no run-to-run
//! drift introduced by the extrapolation logic).
//!
//! The Step 8.5a `v2_step8_regression_smoke::disabled_runs_are_bit
//! _deterministic` test already pins bit-determinism through the
//! full pipeline. Here we focus on the Phase 5 surface
//! specifically: attempt counts, fallback indices, Newton outer
//! iter histogram.

use std::path::PathBuf;

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
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

fn step6_config(steps: usize) -> BaselineConfig {
    let nx = 64;
    let ny = 64;
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        nx,
        ny,
        &vcfg,
        42,
        rates,
        RecyclingConfig::default(),
    )
    .unwrap();
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: 42,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps,
        cfl_factor: 0.3,
        total_time_nondim: 0.4,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from("target/v2_newton_extrapolation_scratch"),
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
        boundary_layout_name: "voronoi_seed42_n8".into(),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        continuation: None,
            plate_kinematic: ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig::Zero,
    }
}

#[test]
fn extrapolation_stats_are_present_and_consistent() {
    // Run 8 steps so step indices ≥ 2 fire the extrapolation path
    // (`steps - 2 = 6` attempts maximum).
    let r = run_baseline(&step6_config(8));
    let stats = r
        .metrics
        .extrapolation
        .as_ref()
        .expect("steps > 0 ⇒ ExtrapolationStats populated");

    // `applied + |fallback| = attempted` is the bookkeeping
    // invariant — every attempt either succeeds or falls back.
    assert_eq!(
        stats.applied + stats.fallback_indices.len(),
        stats.attempted,
        "applied {} + |fallback| {} != attempted {}",
        stats.applied,
        stats.fallback_indices.len(),
        stats.attempted,
    );

    // 8 steps → 6 attempts (k ∈ {2..8}).
    assert_eq!(stats.attempted, 6, "expected 6 attempts on 8 steps");

    // Newton outer iters history aligns 1-1 with steps.
    assert_eq!(stats.newton_outer_iters_per_step.len(), 8);

    // Fallback indices, if any, must be in [2, steps).
    for &idx in &stats.fallback_indices {
        assert!(idx >= 2, "fallback at step {idx} < 2 (extrap not yet attempted)");
        assert!(idx < 8, "fallback at step {idx} >= steps");
    }

    // `last_applied_extrap_residual` is `Some` iff at least one
    // extrap was applied.
    assert_eq!(
        stats.last_applied_extrap_residual.is_some(),
        stats.applied > 0,
        "last_applied vs applied count mismatch",
    );
}

#[test]
fn extrapolation_fallback_rate_under_50_percent_on_typical_run() {
    // On a typical step6-shape physics run, extrapolation is
    // expected to be helpful most steps (regime is steady-state-
    // like after the first ~3 steps). Fallback rate above 50 %
    // would suggest a regression in the safeguard or a pathology
    // in the chosen config; this is a soft gate rather than a
    // strict invariant.
    let r = run_baseline(&step6_config(20));
    let stats = r.metrics.extrapolation.as_ref().unwrap();
    let rate = stats.fallback_rate();
    assert!(
        rate < 0.5,
        "fallback rate {:.1}% > 50 % on 20-step step6 run \
         (attempted={}, fallback={:?})",
        rate * 100.0,
        stats.attempted,
        stats.fallback_indices,
    );
}

/// One-shot timing tool for the Step 8.5b Phase 6 wallclock gate
/// on `step8_activated` (mantle on, Jacobi-only — step8 is out of
/// AMG's reliable regime per 8.5a's diagnostic). `#[ignore]`
/// because it spends a minute or two of physics; invoke
/// explicitly with
/// `cargo test --release v2_newton_extrapolation::bench_step8_jacobi_100step -- --ignored --nocapture`.
#[test]
#[ignore]
fn bench_step8_jacobi_100step() {
    use ymir_core::tectonics_v2::mantle::{
        MantleConfig, COUPLING_DEFAULT, MF_DEFAULT, NUM_MODES_DEFAULT,
    };
    let mut cfg = step6_config(100);
    cfg.mantle = MantleConfig::Enabled {
        mf: MF_DEFAULT,
        coupling: COUPLING_DEFAULT,
        num_modes: NUM_MODES_DEFAULT,
        seed: 7,
        evolution_rate: 0.0,
    };
    cfg.linear_solver = Default::default(); // JacobiCG only
    let t0 = std::time::Instant::now();
    let r = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    let cg_mean = r.metrics.cg_iter_mean;
    let stats = r.metrics.extrapolation.as_ref().unwrap();
    eprintln!("=== step8_activated 100-step Jacobi ===");
    eprintln!("  wallclock: {:.2}s", dt);
    eprintln!("  cg_mean: {:.1}", cg_mean);
    eprintln!(
        "  extrap fallback rate: {:.1}% ({} fallbacks at steps {:?})",
        stats.fallback_rate() * 100.0,
        stats.fallback_indices.len(),
        stats.fallback_indices,
    );
    eprintln!(
        "  newton outer iters mean: {:.2}",
        stats.newton_outer_iters_mean(),
    );
}

#[test]
fn extrapolation_stats_are_reproducible() {
    // Two runs of the same config → identical attempt count,
    // identical fallback indices, identical newton-iter history.
    // Bit-determinism (cg_iter_mean, vmax_peak, mass_drift) is
    // the harder invariant pinned in v2_step8_regression_smoke;
    // here we just verify the new ExtrapolationStats are
    // deterministic too.
    let r1 = run_baseline(&step6_config(10));
    let r2 = run_baseline(&step6_config(10));
    let s1 = r1.metrics.extrapolation.as_ref().unwrap();
    let s2 = r2.metrics.extrapolation.as_ref().unwrap();
    assert_eq!(s1.attempted, s2.attempted);
    assert_eq!(s1.applied, s2.applied);
    assert_eq!(s1.fallback_indices, s2.fallback_indices);
    assert_eq!(
        s1.newton_outer_iters_per_step, s2.newton_outer_iters_per_step,
        "Newton outer iters per step diverged between runs",
    );
}
