//! `gen_bench_data` — Step 8.5a Phase 0 capture driver.
//!
//! Runs short physics scenarios matching each benchmark case's regime
//! and serialises the Stokes linear-solve inputs at a representative
//! timestep via the harness capture hook (`HarnessCaptureSpec`). The
//! resulting files live in `bench_data/<case>.bin` and are consumed
//! by both the `amg_benchmark` criterion harness and the phase 0
//! Jacobi reference tables.
//!
//! Cases (see `step8_5a_issue.md` §Benchmark suite specification):
//!
//!   step0_quiescent        | 64² | step 0  (synthetic η uniform)
//!   step3_floor_yielding   | 64² | step 5  (yielding active, rest off)
//!   step6_voronoi          | 64² | step 50 (Voronoi closed, yielding+drag)
//!   step7_slab_off         | 64² | step 100 (Step 7 regression, slab off)
//!   step8_activated        | 64² | step 100 (Step 8 physics, mantle on)
//!   step8_activated_128    | 128² | step 50 (Step 8 physics, mantle on)
//!
//! Poisson synthetic cases (poisson_constant, poisson_contrast_100,
//! poisson_contrast_10000) are built on-demand inside the benchmark
//! harness and need no on-disk snapshot.
//!
//! ```bash
//! cargo run --release --bin gen_bench_data -- --all
//! cargo run --release --bin gen_bench_data -- --case step3_floor_yielding
//! ```

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, HarnessCaptureSpec, NonlinearChoice,
};
use ymir_core::tectonics_v2::mantle::{MantleConfig, COUPLING_DEFAULT, MF_DEFAULT};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;

const ALL_CASES: &[&str] = &[
    "step0_quiescent",
    "step3_floor_yielding",
    "step6_voronoi",
    "step7_slab_off",
    "step8_activated",
    "step8_activated_128",
];

struct Args {
    cases: Vec<String>,
    output_dir: PathBuf,
    seed: u64,
    num_plates: usize,
    continental_ratio: f64,
    mantle_seed: u64,
    mantle_num_modes: usize,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        cases: Vec::new(),
        output_dir: PathBuf::from("bench_data"),
        seed: 42,
        num_plates: 8,
        continental_ratio: 0.4,
        mantle_seed: 7,
        mantle_num_modes: 6,
    };
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < raw.len() {
        match raw[i].as_str() {
            "--all" => {
                a.cases = ALL_CASES.iter().map(|s| s.to_string()).collect();
            }
            "--case" => {
                i += 1;
                a.cases.push(raw[i].clone());
            }
            "--output-dir" => {
                i += 1;
                a.output_dir = PathBuf::from(&raw[i]);
            }
            "--seed" => {
                i += 1;
                a.seed = raw[i].parse().map_err(|e| format!("bad --seed: {e}"))?;
            }
            "--help" | "-h" => {
                eprintln!(
                    "Usage: gen_bench_data [--all | --case NAME] [--output-dir PATH] [--seed N]\n\
                     \n\
                     Available cases:"
                );
                for c in ALL_CASES {
                    eprintln!("  {c}");
                }
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    if a.cases.is_empty() {
        return Err("pass --all or --case <name>; see --help".into());
    }
    Ok(a)
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("argument error: {e}");
            return ExitCode::from(2);
        }
    };
    if let Err(e) = std::fs::create_dir_all(&args.output_dir) {
        eprintln!("failed to create output dir {}: {}", args.output_dir.display(), e);
        return ExitCode::FAILURE;
    }
    for case in &args.cases {
        if let Err(e) = dispatch(case, &args) {
            eprintln!("[gen_bench_data] case '{}' failed: {}", case, e);
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}

fn dispatch(case: &str, args: &Args) -> Result<(), String> {
    let path = args.output_dir.join(format!("{}.bin", case));
    println!("[gen_bench_data] capturing '{}' → {}", case, path.display());
    let start = Instant::now();
    match case {
        "step0_quiescent" => capture_step0_quiescent(args, &path)?,
        "step3_floor_yielding" => capture_step3_floor_yielding(args, &path)?,
        "step6_voronoi" => capture_step6_voronoi(args, &path)?,
        "step7_slab_off" => capture_step7_slab_off(args, &path)?,
        "step8_activated" => capture_step8_activated(args, 64, 100, &path)?,
        "step8_activated_128" => capture_step8_activated(args, 128, 50, &path)?,
        other => return Err(format!("unknown case '{}'", other)),
    };
    let meta = std::fs::metadata(&path)
        .map_err(|e| format!("stat {}: {}", path.display(), e))?;
    println!(
        "[gen_bench_data]  ✓ '{}' {:.2} MB in {:.1}s",
        case,
        meta.len() as f64 / (1024.0 * 1024.0),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

// --- Case builders -------------------------------------------------------

/// step0_quiescent: uniform η, no yielding, no drag/boundary/slab/mantle.
/// Captured at physics step 0 (first Newton solve of a fresh simulation).
/// Serves as the simplest Stokes benchmark; the report's sanity test
/// cross-checks that a synthetic reconstruction matches this capture
/// to 1e-10.
fn capture_step0_quiescent(args: &Args, path: &PathBuf) -> Result<(), String> {
    let scales = Scales::default();
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let cfg = BaselineConfig {
        seed: args.seed,
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
        output_dir: args.output_dir.join("scratch_step0"),
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
        capture: Some(HarnessCaptureSpec {
            at_step: 0,
            path: path.clone(),
            case_label: "step0_quiescent".into(),
        }),
        linear_solver: Default::default(),
    };
    let _ = run_baseline(&cfg);
    Ok(())
}

/// step3_floor_yielding: yielding Enabled (bi = 0.15), rest off.
/// Captured at step 5 (yielding is active early in the Step 3 regression).
fn capture_step3_floor_yielding(args: &Args, path: &PathBuf) -> Result<(), String> {
    let scales = Scales::default();
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let cfg = BaselineConfig {
        seed: args.seed,
        grid_nx: 64,
        grid_ny: 64,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: 6,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: args.output_dir.join("scratch_step3"),
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
        capture: Some(HarnessCaptureSpec {
            at_step: 5,
            path: path.clone(),
            case_label: "step3_floor_yielding".into(),
        }),
        linear_solver: Default::default(),
    };
    let _ = run_baseline(&cfg);
    Ok(())
}

/// step6_voronoi: Step 6 shape — Voronoi closed boundaries, yielding +
/// basal drag. Captured at step 50 (Voronoi pattern stabilised).
fn capture_step6_voronoi(args: &Args, path: &PathBuf) -> Result<(), String> {
    let cfg = build_step6_shape_config(args, 64, 64, 51, HarnessCaptureSpec {
        at_step: 50,
        path: path.clone(),
        case_label: "step6_voronoi".into(),
    })?;
    let _ = run_baseline(&cfg);
    Ok(())
}

/// step7_slab_off: Step 7 regression shape (slab Disabled) — identical
/// to Step 6 shape, running slightly longer. Captured at step 100.
fn capture_step7_slab_off(args: &Args, path: &PathBuf) -> Result<(), String> {
    let cfg = build_step6_shape_config(args, 64, 64, 101, HarnessCaptureSpec {
        at_step: 100,
        path: path.clone(),
        case_label: "step7_slab_off".into(),
    })?;
    let _ = run_baseline(&cfg);
    Ok(())
}

/// step8_activated: Step 8 physics — Voronoi closed + yielding + drag +
/// mantle Enabled (slab held Disabled per Step 8 co-calibration issue).
/// Captured at the specified step; supports both 64² (step 100) and
/// 128² (step 50) variants.
fn capture_step8_activated(
    args: &Args,
    grid_n: usize,
    capture_step: usize,
    path: &PathBuf,
) -> Result<(), String> {
    let scales = Scales::default();
    let vcfg = VoronoiConfig {
        num_plates: args.num_plates,
        continental_ratio: args.continental_ratio,
    };
    let rates = BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let recycling_config = RecyclingConfig::default();
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        grid_n, grid_n, &vcfg, args.seed, rates, recycling_config,
    )
    .map_err(|e| format!("boundary config invalid: {:?}", e))?;
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let case_label = if grid_n == 128 { "step8_activated_128" } else { "step8_activated" };
    let cfg = BaselineConfig {
        seed: args.seed,
        grid_nx: grid_n,
        grid_ny: grid_n,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: capture_step + 1,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: args.output_dir.join(format!("scratch_step8_{}", grid_n)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw { br: 0.05, ..BasalDragLaw::default() }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", args.seed, args.num_plates),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Enabled {
            mf: MF_DEFAULT,
            coupling: COUPLING_DEFAULT,
            num_modes: args.mantle_num_modes,
            seed: args.mantle_seed,
            evolution_rate: 0.0,
        },
        capture: Some(HarnessCaptureSpec {
            at_step: capture_step,
            path: path.clone(),
            case_label: case_label.into(),
        }),
        linear_solver: Default::default(),
    };
    let _ = run_baseline(&cfg);
    Ok(())
}

fn build_step6_shape_config(
    args: &Args,
    nx: usize,
    ny: usize,
    steps: usize,
    capture: HarnessCaptureSpec,
) -> Result<BaselineConfig, String> {
    let scales = Scales::default();
    let vcfg = VoronoiConfig {
        num_plates: args.num_plates,
        continental_ratio: args.continental_ratio,
    };
    let rates = BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let recycling_config = RecyclingConfig::default();
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        nx, ny, &vcfg, args.seed, rates, recycling_config,
    )
    .map_err(|e| format!("boundary config invalid: {:?}", e))?;
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    Ok(BaselineConfig {
        seed: args.seed,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: args.output_dir.join("scratch_step6_shape"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw { br: 0.05, ..BasalDragLaw::default() }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", args.seed, args.num_plates),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        capture: Some(capture),
        linear_solver: Default::default(),
    })
}
