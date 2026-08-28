//! Step 9 Phase 7 — smoothing-width calibration probe at 64².
//! Acceptance #3 requires `peak_eta_contrast_at_boundary ≤ K · 1.05`.
//! Default `smoothing_width = 0.05` gives `~5.5–5.9` (10 % over).
//! This probe runs a 64²/20-step truncated baseline at several
//! widths to pick the smallest value that meets the bound. Cr=0.3
//! (default), K=5 (default), seed=42.
//!
//! Output is a small table the reviewer can use to decide whether
//! to bump the default to 0.08, 0.10, or larger.
//!
//! Marked `#[ignore]`. Run via:
//! ```text
//! cargo test --release -p ymir-core --test v2_step9_smoothing_probe \
//!     -- --ignored --nocapture
//! ```

use std::path::PathBuf;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::{CratonicConfig, CratonicConfigEnabled};
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, ForceKind, NonlinearChoice, build_force, run_baseline,
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

#[test]
#[ignore]
fn step9_smoothing_width_calibration_probe() {
    let nx = 64usize;
    let ny = 64usize;
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };

    let widths = [0.05_f64, 0.08, 0.10, 0.12, 0.15, 0.20];
    println!();
    println!("Step 9 smoothing-width calibration probe");
    println!("  64x64, 20 steps, Cr=0.3, K=5, seed=42, target ≤ K·1.05 = 5.25");
    println!();
    println!("{:>14} | {:>14} | {:>14}", "smoothing_w", "eta_contrast", "<= 5.25 ?");
    for w in widths {
        let crcfg = CratonicConfigEnabled { smoothing_width: w, ..Default::default() };
        let boundary = BoundaryConfig::enabled_voronoi_closed(
            nx,
            ny,
            &vcfg,
            42,
            rates,
            RecyclingConfig::default(),
        )
        .expect("recycling config valid");
        let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
        let cfg = BaselineConfig {
            seed: 42,
            grid_nx: nx,
            grid_ny: ny,
            domain_lx: 1.0,
            domain_ly: 1.0,
            steps: 20,
            cfl_factor: 0.3,
            total_time_nondim: 1.2,
            preset: preset.clone(),
            nonlinear: NonlinearChoice::Newton,
            newton_cfg: NewtonConfig::default(),
            picard_cfg: PicardConfig::default(),
            heightmap_fractions: Vec::new(),
            output_dir: PathBuf::from(format!(
                "target/v2_step9_smoothing_probe_w{}",
                (w * 100.0) as u32
            )),
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
            boundary_layout_name: format!("voronoi_seed42_n8_w{}", (w * 100.0) as u32),
            slab_pull: SlabPullConfig::Disabled,
            mantle: MantleConfig::Disabled,
            cratonic: CratonicConfig::Enabled(crcfg),
            age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
            capture: None,
            linear_solver: LinearSolverConfig::default(),
            init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
            continuation: None,
            plate_kinematic: ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig::Zero,
        };
        let r = run_baseline(&cfg);
        let na = r.metrics.newton.as_ref().expect("newton aggregate");
        let contrast = na.peak_eta_contrast_at_boundary.unwrap_or(1.0);
        let pass = if contrast <= 5.25 { "PASS" } else { "FAIL" };
        println!("{:>14.3} | {:>14.3} | {:>14}", w, contrast, pass);
    }
}
