//! Step 9 Phase 7 — physics baseline (Cr=0.3, K=5) + Cr sweep at
//! 64² × 100 steps. Marked `#[ignore]` because it takes minutes.
//!
//! Two ignored tests:
//!
//! - `step9_physics_baseline_64sq` — single 64² run with the
//!   default `CratonicConfig::Enabled(default)` parameters,
//!   prints the acceptance metrics in a one-shot summary line and
//!   writes a JSON-ish trace. Prints the comparison table vs Step 8
//!   baseline (extracted from `step8_physics_report.md` per the
//!   issue's "procedure for filling the comparison table").
//!
//! - `step9_cr_sweep_64sq` — runs the same 64² × 100 step shape
//!   for `Cr ∈ {0.1, 0.2, 0.3, 0.4, 0.5}` and reports the per-Cr
//!   metric series. Acceptance #9 requires
//!   `cratonic_cell_fraction` to be monotone non-decreasing in Cr.
//!
//! Run via:
//! ```text
//! cargo test --release -p ymir-core \
//!     --test v2_step9_physics_and_sweep \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! `--test-threads=1` because the two tests both compete for CPU
//! and a serial run gives a cleaner wallclock signal.

use std::path::PathBuf;
use std::time::Instant;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::BoundaryConfig;
use ymir_core::tectonics_v2::cratonic::{CratonicConfig, CratonicConfigEnabled};
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, BaselineResult, ForceKind, NonlinearChoice, build_force, run_baseline,
};
use ymir_core::tectonics_v2::mantle::{
    COUPLING_DEFAULT, MF_DEFAULT, MantleConfig, NUM_MODES_DEFAULT,
};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;
use ymir_core::tectonics_v2::boundaries::BoundaryRates;

const NX: usize = 64;
const NY: usize = 64;
const STEPS: usize = 100;
const SEED: u64 = 42;

/// Step 9 baseline — Step 7 shape (continental + oceanic Voronoï,
/// drag + yielding active, no slab, no mantle) plus the cratonic
/// configuration. We do NOT activate mantle for Step 9 baseline
/// because Step 8 found slab+mantle co-calibration is unresolved
/// (per `project_slab_mantle_cocalibration.md`); Step 9 ships on
/// the Step 7 baseline shape per the issue.
fn build_step9_config(cratonic: CratonicConfig, label: &str) -> BaselineConfig {
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let rates = BoundaryRates {
        k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0,
    };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        NX, NY, &vcfg, SEED, rates, RecyclingConfig::default(),
    )
    .expect("recycling config valid");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: SEED,
        grid_nx: NX,
        grid_ny: NY,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: vec![0.0, 1.0],
        output_dir: PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../docs/reports/step9_phase7_{}", label)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw {
            bi: 0.15, ..Default::default()
        }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: 0.05, ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n8_step9_{}", SEED, label),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic,
        capture: None,
        linear_solver: LinearSolverConfig::default(),
    }
}

fn print_summary(label: &str, dt: f64, r: &BaselineResult) {
    let m = &r.metrics;
    let na = m.newton.as_ref().expect("newton aggregate");
    println!();
    println!("=== Step 9 {} ===", label);
    println!("  wallclock           : {:.2} s ({:.1} ms/step)", dt, dt * 1000.0 / STEPS as f64);
    println!("  CG iters mean       : {:.1}", m.cg_iter_mean);
    println!("  Newton outer mean   : {:.2}", na.outer_iters_mean());
    println!("  peak|v|             : {:.3e}", m.vmax_peak);
    println!(
        "  yielding_cell_fraction_max : {:.4}",
        na.yielding_cell_fraction_max.unwrap_or(0.0)
    );
    if let Some(cr) = na.cr_diagnostic {
        println!("  --- Step 9 cratonic metrics ---");
        println!("  Cr (config)                 : {}", cr);
        println!("  K_viscous (config)          : {}", na.k_viscous_diagnostic.unwrap_or(0.0));
        println!(
            "  cratonic_cell_fraction      : {:.4} (expected {:.4} = Cr·cont_frac)",
            na.cratonic_cell_fraction.unwrap_or(0.0),
            cr * na.continental_cell_fraction.unwrap_or(0.0)
        );
        println!(
            "  continental_cell_fraction   : {:.4}",
            na.continental_cell_fraction.unwrap_or(0.0)
        );
        println!(
            "  peak_yielding_in_craton     : {:.4}  (acceptance #6 ≤ 0.01)",
            na.peak_yielding_in_craton.unwrap_or(0.0)
        );
        println!(
            "  peak_yielding_in_mobile     : {:.4}",
            na.peak_yielding_in_mobile_belt.unwrap_or(0.0)
        );
        println!(
            "  peak_eta_contrast_at_bdry   : {:.3}  (acceptance #3 ≤ K·1.05 = {:.2})",
            na.peak_eta_contrast_at_boundary.unwrap_or(1.0),
            na.k_viscous_diagnostic.unwrap_or(5.0) * 1.05
        );
    }
}

#[test]
#[ignore]
fn step9_physics_baseline_64sq() {
    let cratonic = CratonicConfig::Enabled(CratonicConfigEnabled::default());
    let cfg = build_step9_config(cratonic, "baseline");
    let t0 = Instant::now();
    let r = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    print_summary("Cr=0.3, K=5 baseline", dt, &r);
}

#[test]
#[ignore]
fn step9_baseline_disabled_reference_64sq() {
    // Companion run: same shape but `CratonicConfig::Disabled` —
    // anchors the comparison vs Step 8 (mobile-belt yielding
    // baseline) for acceptance #7.
    let cfg = build_step9_config(CratonicConfig::Disabled, "disabled_reference");
    let t0 = Instant::now();
    let r = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    print_summary("Cratonic Disabled (regression anchor)", dt, &r);
}

#[test]
#[ignore]
fn step9_cr_sweep_64sq() {
    let cr_values = [0.1, 0.2, 0.3, 0.4, 0.5];
    let mut points: Vec<(f64, f64, f64, f64, f64)> = Vec::with_capacity(cr_values.len());
    for &cr in &cr_values {
        let crcfg = CratonicConfigEnabled { cr, ..Default::default() };
        let cfg = build_step9_config(
            CratonicConfig::Enabled(crcfg),
            &format!("cr_{}", (cr * 10.0) as u32),
        );
        let t0 = Instant::now();
        let r = run_baseline(&cfg);
        let dt = t0.elapsed().as_secs_f64();
        let na = r.metrics.newton.as_ref().expect("newton aggregate");
        let crat_frac = na.cratonic_cell_fraction.unwrap_or(0.0);
        let yield_in_craton = na.peak_yielding_in_craton.unwrap_or(0.0);
        let yield_in_mobile = na.peak_yielding_in_mobile_belt.unwrap_or(0.0);
        let contrast = na.peak_eta_contrast_at_boundary.unwrap_or(1.0);
        points.push((cr, crat_frac, yield_in_craton, yield_in_mobile, contrast));
        eprintln!(
            "[cr_sweep] Cr={:.2} done in {:.1}s — crat_frac={:.4}, yield_craton={:.4}, yield_mobile={:.4}, contrast={:.2}",
            cr, dt, crat_frac, yield_in_craton, yield_in_mobile, contrast
        );
    }
    println!();
    println!("=== Step 9 Cr sweep at 64x64, 100 steps ===");
    println!(
        "{:>5} | {:>14} | {:>14} | {:>14} | {:>14}",
        "Cr", "crat_frac", "yield_craton", "yield_mobile", "eta_contrast"
    );
    for (cr, cf, yc, ym, ec) in &points {
        println!(
            "{:>5.2} | {:>14.4} | {:>14.4} | {:>14.4} | {:>14.3}",
            cr, cf, yc, ym, ec
        );
    }
    // Acceptance #9: cratonic_cell_fraction monotone non-decreasing in Cr.
    let mut prev = 0.0_f64;
    for (cr, cf, _, _, _) in &points {
        assert!(
            *cf >= prev - 1e-12,
            "cratonic_cell_fraction not monotone at Cr={}: {} < prev {}",
            cr, cf, prev
        );
        prev = *cf;
    }
    println!();
    println!("Monotonicity acceptance #9: PASS");
}

// ---------------------------------------------------------------------
// Section 2 — Immunity demonstration (Step 8 shape, 32²)
//
// Per the Step 9 acceptance #6 / #7 reformulation: the immunity
// demonstration runs on a *different* shape from the regression
// baseline, because Step 7 shape (no mantle, no slab) does not
// activate yielding (peak|v| ~ 3e-5 < activation threshold) so
// `peak_yielding_in_craton ≤ 0.01` would be vacuously satisfied.
//
// Step 8 shape (mantle on at MF_DEFAULT, slab off — matches the
// Step 8 baseline regime accepted in the milestone) drives the
// system into an active yielding regime. Two runs at 32² × 100
// steps:
//   - `step9_immunity_demo_step8_disabled_32sq`: anchor with
//     `CratonicConfig::Disabled`. Sanity precondition for the
//     immunity test: `yielding_cell_fraction_max > 0` here. If
//     this fails, the immunity test is invalid and a remontée is
//     required before drawing conclusions.
//   - `step9_immunity_demo_step8_enabled_32sq`: same shape, with
//     `CratonicConfig::Enabled(default)`. Acceptances #6
//     (`peak_yielding_in_craton ≤ 0.01`) and #7 (mobile-belt
//     yielding within 10 % of the Disabled value) are evaluated
//     on this run.
//
// 32² Step 8 baseline at Jacobi runs in ~91 s (per
// `step8_5b_baseline_32sq.md`); two runs ≈ 3 min total.

const IMMUNITY_NX: usize = 32;
const IMMUNITY_NY: usize = 32;

fn build_step9_immunity_demo_config(cratonic: CratonicConfig, label: &str) -> BaselineConfig {
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let rates = BoundaryRates {
        k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0,
    };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        IMMUNITY_NX, IMMUNITY_NY, &vcfg, SEED, rates, RecyclingConfig::default(),
    )
    .expect("recycling config valid");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: SEED,
        grid_nx: IMMUNITY_NX,
        grid_ny: IMMUNITY_NY,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: vec![0.0, 1.0],
        output_dir: PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../docs/reports/step9_phase7_immunity_{}", label)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw {
            bi: 0.15, ..Default::default()
        }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: 0.05, ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n8_step9_immunity_{}", SEED, label),
        slab_pull: SlabPullConfig::Disabled,
        // Step 8 shape — mantle enabled at the same defaults as
        // the milestone Step 8 baseline (MF_DEFAULT, COUPLING_DEFAULT,
        // NUM_MODES_DEFAULT, mantle seed 7, evolution_rate = 0).
        mantle: MantleConfig::Enabled {
            mf: MF_DEFAULT,
            coupling: COUPLING_DEFAULT,
            num_modes: NUM_MODES_DEFAULT,
            seed: 7,
            evolution_rate: 0.0,
        },
        cratonic,
        capture: None,
        linear_solver: LinearSolverConfig::default(),
    }
}

fn print_immunity_summary(label: &str, dt: f64, r: &BaselineResult) {
    let m = &r.metrics;
    let na = m.newton.as_ref().expect("newton aggregate");
    println!();
    println!("=== Step 9 immunity demo {} (32x32 Step 8 shape) ===", label);
    println!("  wallclock           : {:.2} s", dt);
    println!("  CG iters mean       : {:.1}", m.cg_iter_mean);
    println!("  Newton outer mean   : {:.2}", na.outer_iters_mean());
    println!("  peak|v|             : {:.3e}", m.vmax_peak);
    println!(
        "  yielding_cell_fraction_max : {:.4}",
        na.yielding_cell_fraction_max.unwrap_or(0.0)
    );
    if let Some(cr) = na.cr_diagnostic {
        println!("  --- Step 9 cratonic metrics ---");
        println!("  Cr (config)                 : {}", cr);
        println!("  K_viscous (config)          : {}", na.k_viscous_diagnostic.unwrap_or(0.0));
        println!(
            "  cratonic_cell_fraction      : {:.4}",
            na.cratonic_cell_fraction.unwrap_or(0.0)
        );
        println!(
            "  peak_yielding_in_craton     : {:.6}  (acceptance #6 ≤ 0.01)",
            na.peak_yielding_in_craton.unwrap_or(0.0)
        );
        println!(
            "  peak_yielding_in_mobile     : {:.4}",
            na.peak_yielding_in_mobile_belt.unwrap_or(0.0)
        );
        println!(
            "  peak_eta_contrast_at_bdry   : {:.3}",
            na.peak_eta_contrast_at_boundary.unwrap_or(1.0)
        );
    }
}

#[test]
#[ignore]
fn step9_immunity_demo_step8_disabled_32sq() {
    let cfg = build_step9_immunity_demo_config(CratonicConfig::Disabled, "disabled");
    let t0 = Instant::now();
    let r = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    print_immunity_summary("Disabled (anchor)", dt, &r);

    // Sanity precondition: the Step 8 shape MUST drive yielding
    // for the immunity test to be meaningful. If yielding is
    // zero here, the immunity demonstration is invalid and a
    // remontée is required before evaluating acceptance #6/#7.
    let yfrac = r
        .metrics
        .newton
        .as_ref()
        .expect("newton aggregate")
        .yielding_cell_fraction_max
        .unwrap_or(0.0);
    assert!(
        yfrac > 0.0,
        "SANITY PRECONDITION FAILED: yielding_cell_fraction_max = {} \
         on Step 8 shape Disabled run. Immunity demo at acceptance #6/#7 \
         is invalid in this regime — remontée required.",
        yfrac
    );
    eprintln!(
        "[immunity_demo] sanity precondition PASS: \
         yielding_cell_fraction_max = {:.4} on Disabled anchor",
        yfrac
    );
}

#[test]
#[ignore]
fn step9_immunity_demo_b_factor_sweep_32sq() {
    // B_factor sweep on Step 8 shape at 32²: characterises the
    // primary plastic-immunity mechanism. Per the §4.10 amendment:
    //   B_factor = 1  : reduces to "K viscous mult only" (the
    //                   pre-amendment behaviour) — expected fail #6,
    //                   validates that Bi elevation is the active
    //                   mechanism.
    //   B_factor ∈ {3, 5, 8, 10} : range bounded by the amendment.
    // The smallest B_factor that drops `peak_yielding_in_craton`
    // below 0.01 wins the "acceptance #6 PASS" annotation; if even
    // B_factor = 10 (top of the amended range) does not pass, a
    // diagnostic remontée is required (per discipline) before any
    // further extension.
    let b_values = [1.0_f64, 3.0, 5.0, 8.0, 10.0];
    let mut points: Vec<(f64, f64, f64, f64, f64, f64)> = Vec::with_capacity(b_values.len());
    for &b in &b_values {
        let crcfg = CratonicConfigEnabled { b_factor: b, ..Default::default() };
        let cfg = build_step9_immunity_demo_config(
            CratonicConfig::Enabled(crcfg),
            &format!("b_factor_{}", b as u32),
        );
        let t0 = Instant::now();
        let r = run_baseline(&cfg);
        let dt = t0.elapsed().as_secs_f64();
        let na = r.metrics.newton.as_ref().expect("newton aggregate");
        let yc = na.peak_yielding_in_craton.unwrap_or(0.0);
        let ym = na.peak_yielding_in_mobile_belt.unwrap_or(0.0);
        let yt = na.yielding_cell_fraction_max.unwrap_or(0.0);
        let cg_mean = r.metrics.cg_iter_mean;
        let peak_v = r.metrics.vmax_peak;
        points.push((b, yc, ym, yt, cg_mean, peak_v));
        eprintln!(
            "[b_sweep] B_factor={:.0} done in {:.1}s — yc={:.6}, ym={:.4}, yt={:.4}, cg_mean={:.0}, peak_v={:.2}",
            b, dt, yc, ym, yt, cg_mean, peak_v
        );
    }
    println!();
    println!("=== Step 9 B_factor sweep, 32x32 Step 8 shape, 100 steps ===");
    println!(
        "{:>8} | {:>12} | {:>12} | {:>12} | {:>10} | {:>10}",
        "B_factor",
        "peak_yc (#6)",
        "peak_ym",
        "y_total_max",
        "CG mean",
        "peak|v|"
    );
    for (b, yc, ym, yt, cg, pv) in &points {
        println!(
            "{:>8.0} | {:>12.6} | {:>12.4} | {:>12.4} | {:>10.0} | {:>10.3}",
            b, yc, ym, yt, cg, pv
        );
    }
    // Acceptance #6: B_factor=1 should fail (validates Bi
    // elevation IS the mechanism); the smallest B_factor reaching
    // ≤ 0.01 is reported. We do NOT panic the test if no value
    // passes — the table is the diagnostic, and a fail-everywhere
    // result triggers a remontée per discipline.
    let passing = points
        .iter()
        .find(|(_, yc, _, _, _, _)| *yc <= 0.01)
        .map(|(b, _, _, _, _, _)| *b);
    println!();
    match passing {
        Some(b) => println!(
            "Acceptance #6 (peak_yielding_in_craton ≤ 0.01) PASS at B_factor = {}",
            b
        ),
        None => println!(
            "Acceptance #6 NOT MET at any B_factor in [1, 10]. \
             Remontée required — do NOT extend B_factor beyond 10 \
             without architectural review."
        ),
    }
    // Sanity: B_factor=1 must fail (pre-amendment baseline). If
    // it passes, the test is degenerate (no yielding regime) and
    // results are not informative.
    let yc_b1 = points[0].1;
    assert!(
        yc_b1 > 0.01,
        "SANITY: B_factor=1 (pre-amendment) shows peak_yielding_in_craton = {} ≤ 0.01. \
         Either the regime is not yielding-active (sanity precondition fail in disabled \
         anchor) or the metric is broken. Remontée required.",
        yc_b1
    );
}

#[test]
#[ignore]
fn step9_immunity_demo_step8_enabled_32sq() {
    let cratonic = CratonicConfig::Enabled(CratonicConfigEnabled::default());
    let cfg = build_step9_immunity_demo_config(cratonic, "enabled");
    let t0 = Instant::now();
    let r = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    print_immunity_summary("Enabled (Cr=0.3, K=5)", dt, &r);

    let na = r
        .metrics
        .newton
        .as_ref()
        .expect("newton aggregate");
    let peak_yc = na.peak_yielding_in_craton.unwrap_or(0.0);
    let peak_ym = na.peak_yielding_in_mobile_belt.unwrap_or(0.0);
    let yfrac_total = na.yielding_cell_fraction_max.unwrap_or(0.0);

    // Sanity: enabled run must also yield (otherwise the cratonic
    // K mechanism has accidentally suppressed all yielding, which
    // would mask the immunity test as vacuously passing).
    assert!(
        yfrac_total > 0.0,
        "SANITY: enabled run shows zero yielding everywhere \
         (peak yielding total = {}). Cratonic K mechanism may have \
         globally suppressed yielding — remontée required.",
        yfrac_total
    );

    // Acceptance #6: peak_yielding_in_craton ≤ 0.01.
    // PER ISSUE DISCIPLINE: failure here triggers diagnostic
    // remontée, NOT silent threshold relaxation. Don't mask
    // a fail by widening the bound.
    assert!(
        peak_yc <= 0.01,
        "ACCEPTANCE #6 FAILED: peak_yielding_in_craton = {} > 0.01 \
         (cratons are yielding). Diagnostic remontée required — do NOT \
         tune Cr / K / smoothing_width to mask this.",
        peak_yc
    );
    eprintln!("[immunity_demo] ACCEPTANCE #6 PASS: peak_yielding_in_craton = {:.6}", peak_yc);
    eprintln!(
        "[immunity_demo] mobile-belt yielding peak = {:.4} (compare to Disabled run for #7)",
        peak_ym
    );
}

