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
use ymir_core::tectonics_v2::mantle::MantleConfig;
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
