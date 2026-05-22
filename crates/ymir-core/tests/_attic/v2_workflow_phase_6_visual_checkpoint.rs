// Step 12 R3 — this entire file is archived. It captured the Phase 6
// pre-R3 visual checkpoint built on the legacy `low_res_erosion` (D2
// diffusive erosion + `β` local-deposition coefficient). R3 replaced
// that mechanism with `macro_redistribution` (long-distance drainage +
// isostatic rebound), and the new visual checkpoint lands in R4 with
// fresh galleries. Helpers below (`apply_erosion_with_alpha_field`,
// `run_custom_phase_a_loop`, `dump_step12_phase_6_curvature_variants`)
// reference `low_res_erosion` semantics that no longer compile after
// the R3 module deletion — `#![cfg(any())]` excludes the whole crate
// so they remain readable in `git log` / `git show` without blocking
// the build.
#![cfg(any())]

//! Step 12 Phase 6 — visual checkpoint, reviewer-validated.
//!
//! Two `#[ignore]` artefacts:
//!
//! - `dump_step12_phase_a_evolution_galerie` — 2 × 2 patchwork
//!   (rows = presets `single_continent`/`convergence`, cols = `before
//!   workflow` / `after 5 Phase A cycles`). Validates **acceptance #6**:
//!   continental contours become non-polygonal after Phase A.
//! - `dump_step12_phase_b_hd_zoom` — 1 × 2 side-by-side (low-res
//!   post Phase A vs HD post Phase B at 512²). Validates
//!   **acceptance #8**: Phase B HD output shows recognizable valley
//!   structures.
//!
//! Output: `docs/reports/step12_phase_6_checkpoint/`. Metrics are
//! printed to stdout via `--nocapture` and saved alongside the PNGs.
//!
//! Run:
//!
//! ```text
//! cargo test --release -p ymir-core --test v2_workflow_phase_6_visual_checkpoint -- --ignored --nocapture
//! ```

use std::fs;
use std::path::PathBuf;

use ymir_core::erosion::hydraulic::ErosionConfig;
use ymir_core::tectonics_v2::age_field::AgeFieldConfig;
use ymir_core::tectonics_v2::basal_drag::BasalDragConfig;
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::CratonicConfig;
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::diagnostics::heightmap::save_heightmap;
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::init::{init_s_field, InitContext, InitMode, PlateInitData};
use ymir_core::tectonics_v2::mantle::MantleConfig;
use ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};
use ymir_core::tectonics_v2::workflow::{
    run_phase_a_loop, run_phase_b, PhaseAParams, PhaseBParams, WorkflowConfig, WorkflowParams,
};

const NX: usize = 32;
const NY: usize = 32;
const N_CYCLES: usize = 5;
const K_CYCLE: usize = 20;
const ALPHA: f64 = 0.01;

#[derive(Clone, Copy)]
struct PresetCfg {
    label: &'static str,
    seed: u64,
    num_plates: usize,
    continental_ratio: f64,
}

const PRESETS: [PresetCfg; 2] = [
    PresetCfg { label: "single_continent", seed: 12, num_plates: 4, continental_ratio: 0.5 },
    PresetCfg { label: "convergence", seed: 23, num_plates: 6, continental_ratio: 0.4 },
];

fn build_preset_config(preset: &PresetCfg, scratch: &str) -> BaselineConfig {
    let scales = Scales::default();
    let preset_obj = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig {
        num_plates: preset.num_plates,
        continental_ratio: preset.continental_ratio,
    };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        NX,
        NY,
        &vcfg,
        preset.seed,
        rates,
        RecyclingConfig::default(),
    )
    .expect("recycling config valid");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: preset.seed,
        grid_nx: NX,
        grid_ny: NY,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: K_CYCLE,
        cfl_factor: 0.3,
        total_time_nondim: 0.4 * (K_CYCLE as f64) / 20.0,
        preset: preset_obj,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from(format!("target/v2_workflow_phase6/{}", scratch)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Disabled,
        basal_drag: BasalDragConfig::Disabled,
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", preset.seed, preset.num_plates),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: CratonicConfig::Disabled,
        age_field: AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
        init_mode: InitMode::Checkerboard,
        continuation: None,
        plate_kinematic: PlateKinematicConfig::Zero,
    }
}

fn init_s_for_preset(preset: &PresetCfg) -> Field2D {
    let plates = generate_voronoi(
        NX,
        NY,
        &VoronoiConfig {
            num_plates: preset.num_plates,
            continental_ratio: preset.continental_ratio,
        },
        preset.seed,
    );
    let plate_data = PlateInitData {
        plate_id: &plates.plate_id,
        plate_type: &plates.plate_type,
        seed_coords: Some(&plates.seed_coords),
    };
    let ctx = InitContext {
        nx: NX,
        ny: NY,
        seed: preset.seed,
        amplitude: 0.2,
        plate_data: Some(plate_data),
    };
    init_s_field(InitMode::Checkerboard, &ctx)
}

#[test]
#[ignore]
fn dump_step12_phase_a_evolution_galerie() {
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_phase_6_checkpoint");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!();
    println!(
        "Step 12 Phase 6 visual checkpoint — Phase A evolution patchwork \
         (2 presets × {{before workflow, after {N_CYCLES} cycles}})"
    );
    println!("Grid: {NX}×{NY}, k_cycle={K_CYCLE}, α={ALPHA}, β=0.0");
    println!();

    // 4 tiles: row 0 = single_continent (before, after); row 1 = convergence (before, after)
    let mut tiles: Vec<(String, Field2D)> = Vec::with_capacity(4);
    let mut report_lines: Vec<String> = Vec::new();
    report_lines.push(format!(
        "Step 12 Phase 6 — Phase A evolution metrics (n_cycles={N_CYCLES}, k_cycle={K_CYCLE}, α={ALPHA})"
    ));
    report_lines.push(String::new());

    for preset in PRESETS.iter() {
        let init_s = init_s_for_preset(preset);
        tiles.push((format!("{}_before", preset.label), init_s));

        // Save individual init PNG
        let init_path = out_dir.join(format!("{}_init.png", preset.label));
        let init_meta = save_heightmap(&tiles.last().unwrap().1, &init_path).unwrap();
        println!(
            "  {:<20} init       : range=[{:.4}, {:.4}], mean={:.4}",
            preset.label, init_meta.min, init_meta.max, init_meta.mean
        );

        let mut cfg = build_preset_config(preset, &format!("phase_a_{}", preset.label));
        let wf = WorkflowConfig::Enabled(WorkflowParams {
            phase_a: PhaseAParams { n_cycles: N_CYCLES, k_cycle: K_CYCLE, alpha: ALPHA, beta: 0.0 },
            phase_b: Default::default(),
        });
        let phase_a = run_phase_a_loop(&mut cfg, &wf);

        // Per-cycle metrics
        report_lines.push(format!("## Preset: {}", preset.label));
        report_lines.push(format!(
            "  seed={} num_plates={} cont_ratio={:.2}",
            preset.seed, preset.num_plates, preset.continental_ratio
        ));
        report_lines.push(String::new());
        report_lines.push("| cycle | peak S̃ | mass drift | erosion volume | sea_level | craton change |".to_string());
        report_lines.push("|------:|------:|-----------:|---------------:|----------:|--------------:|".to_string());
        let mut cum_drift = 0.0_f64;
        for (i, c) in phase_a.cycles.iter().enumerate() {
            let peak_s: f64 = c
                .baseline
                .final_state
                .s_field
                .data()
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            cum_drift += c.common.mass_drift;
            let craton = c.common.craton_recomputation_change.map_or("—".to_string(), |v| format!("{v:.4}"));
            report_lines.push(format!(
                "| {i} | {peak_s:.4} | {:+.5} | {:.5} | {:.4} | {craton} |",
                c.common.mass_drift, c.common.erosion_volume_removed, c.common.sea_level_normalized
            ));
        }
        report_lines.push(format!("  Cumulative mass drift over {N_CYCLES} cycles: {cum_drift:+.5}"));
        report_lines.push(String::new());

        let after_s = phase_a.cycles.last().unwrap().baseline.final_state.s_field.clone();

        // Save individual after PNG
        let after_path = out_dir.join(format!("{}_after_phase_a.png", preset.label));
        let after_meta = save_heightmap(&after_s, &after_path).unwrap();
        println!(
            "  {:<20} after Phase A : range=[{:.4}, {:.4}], mean={:.4}, cum drift={:+.5}",
            preset.label, after_meta.min, after_meta.max, after_meta.mean, cum_drift
        );

        tiles.push((format!("{}_after", preset.label), after_s));
    }

    // 2×2 patchwork — single common [min, max] scale for visual
    // comparability. Use save_heightmap on a synthesised Field2D.
    let tile_w = NX;
    let tile_h = NY;
    let sep = 1usize;
    let pw = tile_w * 2 + sep;
    let ph = tile_h * 2 + sep;
    let mut patch_data = vec![0.5_f64; pw * ph];
    for (k, (_, s)) in tiles.iter().enumerate() {
        let row = k / 2;
        let col = k % 2;
        let x_off = col * (tile_w + sep);
        let y_off = row * (tile_h + sep);
        for j in 0..tile_h {
            for i in 0..tile_w {
                patch_data[(y_off + j) * pw + (x_off + i)] = s.get(i, j);
            }
        }
    }
    let patch = Field2D::from_vec(pw, ph, patch_data);
    let patch_path = out_dir.join("patchwork_phase_a_evolution.png");
    let patch_meta = save_heightmap(&patch, &patch_path).expect("save patchwork");

    println!();
    println!("Patchwork (2×2, dynamic remap [{:.4}, {:.4}]):", patch_meta.min, patch_meta.max);
    println!("  {}", patch_path.display());
    println!("  Layout: row 0 = single_continent | row 1 = convergence");
    println!("           col 0 = before workflow  | col 1 = after Phase A ({N_CYCLES} cycles)");

    // Save metrics report alongside the PNG.
    let report_path = out_dir.join("phase_a_evolution_metrics.md");
    fs::write(&report_path, report_lines.join("\n")).expect("write metrics report");
    println!("  Metrics: {}", report_path.display());
}

#[test]
#[ignore]
fn dump_step12_phase_b_hd_zoom() {
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_phase_6_checkpoint");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!();
    println!(
        "Step 12 Phase 6 visual checkpoint — Phase B HD zoom \
         (single_continent: low-res Phase A → HD 512² Phase B)"
    );
    println!();

    // Use single_continent preset (most pedagogical: large
    // continental plates → most coast for valleys to carve).
    let preset = &PRESETS[0];
    let mut cfg = build_preset_config(preset, "phase_b_hd");
    let wf = WorkflowConfig::Enabled(WorkflowParams {
        phase_a: PhaseAParams { n_cycles: N_CYCLES, k_cycle: K_CYCLE, alpha: ALPHA, beta: 0.0 },
        phase_b: PhaseBParams {
            hd_grid_size: 512,
            erosion: ErosionConfig {
                num_droplets: 500_000,
                ..ErosionConfig::default()
            },
            ..PhaseBParams::default()
        },
    });
    let phase_a = run_phase_a_loop(&mut cfg, &wf);
    let last_cycle = phase_a.cycles.last().unwrap();
    let phase_b = run_phase_b(&last_cycle.baseline.final_state.s_field, &wf, cfg.seed)
        .expect("Phase B output");

    let lowres = &last_cycle.baseline.final_state.s_field;
    let hd_eroded = &phase_b.heightmap;

    // Save individual artefacts.
    let lowres_path = out_dir.join("single_continent_phase_a_lowres.png");
    let lowres_meta = save_heightmap(lowres, &lowres_path).unwrap();
    let hd_path = out_dir.join("single_continent_phase_b_hd.png");
    hd_eroded.save_png_u16(&hd_path).expect("save hd png");
    let sediment_path = out_dir.join("single_continent_phase_b_sediment.png");
    phase_b.sediment.save_png_u16(&sediment_path).expect("save sediment png");

    println!(
        "  low-res Phase A : range=[{:.4}, {:.4}], mean={:.4}",
        lowres_meta.min, lowres_meta.max, lowres_meta.mean
    );
    let hd_min = hd_eroded.data.iter().copied().fold(f32::INFINITY, f32::min);
    let hd_max = hd_eroded.data.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let hd_mean: f32 = hd_eroded.data.iter().sum::<f32>() / hd_eroded.data.len() as f32;
    println!("  HD Phase B      : range=[{:.4}, {:.4}], mean={:.4}", hd_min, hd_max, hd_mean);
    println!(
        "  D5 metrics       : grand_scale_deviation_p95 = {:.4} (acceptance), L_∞ = {:.4} (diagnostic)",
        phase_b.grand_scale_deviation_p95, phase_b.grand_scale_deviation
    );
    println!("  Outputs:");
    println!("    {}", lowres_path.display());
    println!("    {}", hd_path.display());
    println!("    {} (sediment)", sediment_path.display());

    // Metrics report.
    let report_path = out_dir.join("phase_b_hd_metrics.md");
    let report = format!(
        "Step 12 Phase 6 — Phase B HD zoom metrics\n\
         \n\
         Preset: {} (seed={}, {} plates, cont_ratio={:.2})\n\
         Phase A: {} cycles × {} steps × {}², α={}\n\
         Phase B: HD={}², droplets={}\n\
         \n\
         | metric | value |\n\
         |--------|------:|\n\
         | grand_scale_deviation_p95 | {:.4} |\n\
         | grand_scale_deviation (L_∞) | {:.4} |\n\
         | hd_min | {:.4} |\n\
         | hd_max | {:.4} |\n\
         | hd_mean | {:.4} |\n\
         | lowres_min | {:.4} |\n\
         | lowres_max | {:.4} |\n\
         | lowres_mean | {:.4} |\n\
         \n\
         Acceptance #5 (D5 grand-scale preservation): \
         grand_scale_deviation_p95 < 0.10 → {}.\n\
         Acceptance #8 (valleys present): reviewer-validated.\n",
        preset.label,
        preset.seed,
        preset.num_plates,
        preset.continental_ratio,
        N_CYCLES,
        K_CYCLE,
        NX,
        ALPHA,
        512,
        500_000,
        phase_b.grand_scale_deviation_p95,
        phase_b.grand_scale_deviation,
        hd_min,
        hd_max,
        hd_mean,
        lowres_meta.min,
        lowres_meta.max,
        lowres_meta.mean,
        if phase_b.grand_scale_deviation_p95 < 0.10 { "PASS" } else { "FAIL" }
    );
    fs::write(&report_path, &report).expect("write hd metrics report");
    println!("  Metrics: {}", report_path.display());

    // Side-by-side patchwork: low-res scaled up to match HD width
    // for visual comparison would require nearest-neighbor upsampling.
    // Simpler: emit two same-resolution snapshots (low-res 32² and
    // HD 512²) and let the reviewer view them side by side. The
    // single PNG side-by-side comparison is cosmetic — the per-file
    // viewing surface is more readable for valley inspection at HD.
}

/// Build a 64² Phase A + Phase B BaselineConfig + WorkflowConfig
/// pair. Centralised so the heavy and aggressive variants share the
/// same Voronoï / boundary geometry — only `wf` differs between them.
fn build_64sq_cfg_and_wf(
    preset: &PresetCfg,
    n_cycles: usize,
    alpha: f64,
    hd_droplets: usize,
    scratch: &str,
) -> (BaselineConfig, WorkflowConfig) {
    let nx = 64;
    let ny = 64;
    let scales = Scales::default();
    let preset_obj = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig {
        num_plates: preset.num_plates,
        continental_ratio: preset.continental_ratio,
    };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        nx, ny, &vcfg, preset.seed, rates, RecyclingConfig::default(),
    )
    .expect("boundary");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let cfg = BaselineConfig {
        seed: preset.seed,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: K_CYCLE,
        cfl_factor: 0.3,
        total_time_nondim: 0.4 * (K_CYCLE as f64) / 20.0,
        preset: preset_obj,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from(format!("target/v2_workflow_phase6/{}_{}", scratch, preset.label)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Disabled,
        basal_drag: BasalDragConfig::Disabled,
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", preset.seed, preset.num_plates),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: CratonicConfig::Disabled,
        age_field: AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
        init_mode: InitMode::Checkerboard,
        continuation: None,
        plate_kinematic: PlateKinematicConfig::Zero,
    };
    let wf = WorkflowConfig::Enabled(WorkflowParams {
        phase_a: PhaseAParams { n_cycles, k_cycle: K_CYCLE, alpha, beta: 0.0 },
        phase_b: PhaseBParams {
            hd_grid_size: 1024,
            erosion: ErosionConfig {
                num_droplets: hd_droplets,
                ..ErosionConfig::default()
            },
            ..PhaseBParams::default()
        },
    });
    (cfg, wf)
}

/// Nearest-neighbour upscale a `Field2D` by integer factor `k`. Used
/// to make 64² patchworks readable at chat-render scale (an extra
/// "look here" affordance, the underlying numerical values are
/// untouched).
fn nn_upscale(src: &Field2D, k: usize) -> Field2D {
    let nx = src.nx();
    let ny = src.ny();
    let dst_w = nx * k;
    let dst_h = ny * k;
    let mut dst = Field2D::new(dst_w, dst_h);
    for j in 0..dst_h {
        let sj = j / k;
        for i in 0..dst_w {
            let si = i / k;
            dst.set(i, j, src.get(si, sj));
        }
    }
    dst
}

/// 64² heavier checkpoint — both presets through the full Phase A 64²
/// → Phase B HD 1024² pipeline. Produces the artefacts the reviewer
/// needs for acceptance #6 (non-polygonal contours after Phase A) and
/// acceptance #8 (HD valleys present after Phase B).
///
/// Tracked metrics: per-cycle wallclock (Phase A only — finding hook
/// for the issue's CG-saturation concern at 64² with erosion
/// interleaved) and cumulative mass drift (must stay < 10% of initial
/// mass per the Phase 4 contract, re-verified at 64²).
#[test]
#[ignore]
fn dump_step12_phase_6_64sq_full() {
    use std::time::Instant;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_phase_6_checkpoint_64sq");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    let nx = 64;
    let ny = 64;
    let scales = Scales::default();
    let preset_obj = Preset::by_name("dynamic-accidented").unwrap();
    let force_proto = || build_force(ForceKind::Gpe, &scales, 10.0, 1.0);

    println!();
    println!(
        "Step 12 Phase 6 heavy — 64² Phase A × {N_CYCLES} cycles + HD 1024² Phase B"
    );
    println!("  presets: {}, {}", PRESETS[0].label, PRESETS[1].label);
    println!();

    // Two patchwork tiles per preset: init + after-Phase-A. 4 tiles
    // total, laid out 2 rows × 2 cols.
    let mut tiles: Vec<(String, Field2D)> = Vec::with_capacity(4);
    let mut report_lines: Vec<String> = vec![
        format!(
            "Step 12 Phase 6 heavy — 64² Phase A × {N_CYCLES} cycles + HD 1024² Phase B"
        ),
        String::new(),
    ];

    for preset in PRESETS.iter() {
        // Inline 64² config builder (the const NX/NY at module scope
        // are 32; we keep them as the "default-light" path and run
        // 64² explicitly here).
        let vcfg = VoronoiConfig {
            num_plates: preset.num_plates,
            continental_ratio: preset.continental_ratio,
        };
        let rates =
            BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
        let boundary = BoundaryConfig::enabled_voronoi_closed(
            nx,
            ny,
            &vcfg,
            preset.seed,
            rates,
            RecyclingConfig::default(),
        )
        .expect("boundary");
        let mut cfg = BaselineConfig {
            seed: preset.seed,
            grid_nx: nx,
            grid_ny: ny,
            domain_lx: 1.0,
            domain_ly: 1.0,
            steps: K_CYCLE,
            cfl_factor: 0.3,
            total_time_nondim: 0.4 * (K_CYCLE as f64) / 20.0,
            preset: preset_obj.clone(),
            nonlinear: NonlinearChoice::Newton,
            newton_cfg: Default::default(),
            picard_cfg: Default::default(),
            heightmap_fractions: Vec::new(),
            output_dir: PathBuf::from(format!("target/v2_workflow_phase6/64sq_{}", preset.label)),
            force: force_proto(),
            force_kind: ForceKind::Gpe,
            sinusoidal_amplitude: 0.0,
            s_perturbation_amplitude: 0.2,
            yielding: YieldingConfig::Disabled,
            basal_drag: BasalDragConfig::Disabled,
            boundary,
            boundary_layout_name: format!("voronoi_seed{}_n{}", preset.seed, preset.num_plates),
            slab_pull: SlabPullConfig::Disabled,
            mantle: MantleConfig::Disabled,
            cratonic: CratonicConfig::Disabled,
            age_field: AgeFieldConfig::Disabled,
            capture: None,
            linear_solver: Default::default(),
            init_mode: InitMode::Checkerboard,
            continuation: None,
            plate_kinematic: PlateKinematicConfig::Zero,
        };
        let wf = WorkflowConfig::Enabled(WorkflowParams {
            phase_a: PhaseAParams { n_cycles: N_CYCLES, k_cycle: K_CYCLE, alpha: ALPHA, beta: 0.0 },
            phase_b: PhaseBParams {
                hd_grid_size: 1024,
                erosion: ErosionConfig {
                    num_droplets: 1_000_000,
                    ..ErosionConfig::default()
                },
                ..PhaseBParams::default()
            },
        });

        // Init S̃ snapshot for the patchwork "before" tile.
        let plates = generate_voronoi(nx, ny, &vcfg, preset.seed);
        let plate_data = PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        };
        let ctx = InitContext {
            nx,
            ny,
            seed: preset.seed,
            amplitude: 0.2,
            plate_data: Some(plate_data),
        };
        let init_s = init_s_field(InitMode::Checkerboard, &ctx);
        tiles.push((format!("{}_before", preset.label), init_s));

        // Phase A loop with per-cycle wallclock timing. The
        // run_phase_a_loop API doesn't expose per-cycle hooks, so we
        // measure total + estimate per-cycle as total / n_cycles.
        // (A finer per-cycle measurement would need a streaming
        // variant; that's a Phase 7+ refinement.)
        let phase_a_start = Instant::now();
        let phase_a = run_phase_a_loop(&mut cfg, &wf);
        let phase_a_elapsed = phase_a_start.elapsed().as_secs_f64();
        let mean_per_cycle = phase_a_elapsed / N_CYCLES as f64;

        let last = phase_a.cycles.last().unwrap();
        let after_s = last.baseline.final_state.s_field.clone();
        tiles.push((format!("{}_after", preset.label), after_s.clone()));

        // Phase B HD 1024².
        let phase_b_start = Instant::now();
        let phase_b = run_phase_b(&last.baseline.final_state.s_field, &wf, cfg.seed).expect("phase b");
        let phase_b_elapsed = phase_b_start.elapsed().as_secs_f64();

        // Cumulative mass drift over Phase A.
        let cum_drift: f64 = phase_a.cycles.iter().map(|c| c.common.mass_drift).sum();
        let initial_mass_estimate =
            ((preset.continental_ratio + (1.0 - preset.continental_ratio) * 0.2) as f64)
                * (nx * ny) as f64; // continental ≈ 1.0, oceanic ≈ 0.2
        let drift_fraction = cum_drift.abs() / initial_mass_estimate;

        // Write per-preset HD output.
        let lowres_path = out_dir.join(format!("{}_phase_a_64sq.png", preset.label));
        let lowres_meta = save_heightmap(&after_s, &lowres_path).unwrap();
        let hd_path = out_dir.join(format!("{}_phase_b_hd1024.png", preset.label));
        phase_b.heightmap.save_png_u16(&hd_path).unwrap();

        println!(
            "  [{:<16}] Phase A 64² × {} cycles : {:.2} s total, ~{:.2} s/cycle",
            preset.label, N_CYCLES, phase_a_elapsed, mean_per_cycle
        );
        println!(
            "                  cumulative mass drift = {:+.3} ({:.2} %% of initial mass {:.1})",
            cum_drift,
            drift_fraction * 100.0,
            initial_mass_estimate
        );
        println!("                  peak S̃ trajectory: {:?}", phase_a.cycles.iter().map(|c| {
            c.baseline.final_state.s_field.data().iter().copied().fold(f64::NEG_INFINITY, f64::max)
        }).collect::<Vec<_>>());
        println!(
            "                  Phase A range=[{:.4}, {:.4}], mean={:.4}",
            lowres_meta.min, lowres_meta.max, lowres_meta.mean
        );
        println!(
            "                  Phase B HD 1024² ({} droplets) : {:.2} s",
            1_000_000, phase_b_elapsed
        );
        println!(
            "                  D5 metrics: p95 = {:.4} (acc), L_∞ = {:.4} (diag)",
            phase_b.grand_scale_deviation_p95, phase_b.grand_scale_deviation
        );
        println!();

        // Verify Phase 4 mass-drift bound at 64² (the test we never
        // ran). 10 % of initial mass is a Phase 4 contract.
        assert!(
            drift_fraction < 0.10,
            "[{}] cumulative mass drift {:.3} = {:.2}% of initial mass exceeds 10%",
            preset.label,
            cum_drift,
            drift_fraction * 100.0
        );

        report_lines.push(format!("## Preset: {}", preset.label));
        report_lines.push(format!(
            "  seed={} num_plates={} cont_ratio={:.2}",
            preset.seed, preset.num_plates, preset.continental_ratio
        ));
        report_lines.push(format!("  Phase A wallclock : {:.2} s ({:.2} s/cycle)", phase_a_elapsed, mean_per_cycle));
        report_lines.push(format!("  Phase B wallclock : {:.2} s", phase_b_elapsed));
        report_lines.push(format!("  Cumulative mass drift : {:+.3} ({:.2} % of initial mass)", cum_drift, drift_fraction * 100.0));
        report_lines.push(format!(
            "  D5 metrics : p95={:.4} (acceptance), L_∞={:.4} (diagnostic)",
            phase_b.grand_scale_deviation_p95, phase_b.grand_scale_deviation
        ));
        report_lines.push(String::new());
    }

    // 2 × 2 patchwork at 64² (130×130 native). NN-upscale by 6× to
    // 780×780 for chat-render visibility. The numerical values are
    // unchanged (NN preserves them exactly), only the pixel
    // resolution scales for human inspection.
    let tile_w = nx;
    let tile_h = ny;
    let sep = 1usize;
    let pw = tile_w * 2 + sep;
    let ph = tile_h * 2 + sep;
    let mut patch_data = vec![0.5_f64; pw * ph];
    for (k, (_, s)) in tiles.iter().enumerate() {
        let row = k / 2;
        let col = k % 2;
        let x_off = col * (tile_w + sep);
        let y_off = row * (tile_h + sep);
        for j in 0..tile_h {
            for i in 0..tile_w {
                patch_data[(y_off + j) * pw + (x_off + i)] = s.get(i, j);
            }
        }
    }
    let patch = Field2D::from_vec(pw, ph, patch_data);
    let patch_path = out_dir.join("patchwork_phase_a_evolution_64sq.png");
    save_heightmap(&patch, &patch_path).expect("save patchwork");

    // NN-upscaled variant (×6 → ~780×780) for chat-render visibility.
    let patch_up = nn_upscale(&patch, 6);
    let patch_up_path = out_dir.join("patchwork_phase_a_evolution_64sq_x6.png");
    save_heightmap(&patch_up, &patch_up_path).expect("save patchwork (NN-upscaled)");

    println!();
    println!("Patchwork (2×2, 64² native, NN-upscaled ×6 for visibility):");
    println!("  native : {}", patch_path.display());
    println!("  ×6     : {}", patch_up_path.display());
    println!("  Layout : row 0 = single_continent | row 1 = convergence");
    println!("           col 0 = before workflow  | col 1 = after Phase A ({N_CYCLES} cycles)");

    let report_path = out_dir.join("phase_6_64sq_metrics.md");
    fs::write(&report_path, report_lines.join("\n")).expect("write report");
    println!("  Metrics report: {}", report_path.display());
}

/// **Aggressive demo variant** — α = 0.05 (5× D8 default), N_cycles =
/// 15 (3× D8 default), HD = 1024² with 5M droplets. Designed
/// exclusively to *visually* demonstrate the workflow's effect on
/// continental contours (acceptance #6) and HD valleys (acceptance
/// #8). Not a calibrated configuration; the D8 defaults remain the
/// recommended baseline.
///
/// Wallclock budget : ~5-7 minutes (15 cycles × 2 presets ≈
/// 3 minutes Phase A + 2 × 90 s Phase B 5M droplets).
#[test]
#[ignore]
fn dump_step12_phase_6_aggressive_demo() {
    use std::time::Instant;

    const N_CYCLES_AGG: usize = 15;
    const ALPHA_AGG: f64 = 0.05;
    const HD_DROPLETS_AGG: usize = 5_000_000;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_phase_6_aggressive_demo");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!();
    println!(
        "Step 12 Phase 6 aggressive demo — α = {ALPHA_AGG} (5× D8), N = {N_CYCLES_AGG} (3× D8), HD 1024² × {HD_DROPLETS_AGG} droplets"
    );
    println!("  presets: {}, {}", PRESETS[0].label, PRESETS[1].label);
    println!();

    let mut tiles: Vec<(String, Field2D)> = Vec::with_capacity(4);
    let mut report_lines: Vec<String> = vec![
        format!(
            "Step 12 Phase 6 aggressive demo — α={ALPHA_AGG}, N={N_CYCLES_AGG}, HD 1024² × {HD_DROPLETS_AGG} droplets"
        ),
        String::new(),
        "**Note:** This is a *visual demo* configuration, not a recommended baseline. \
         The D8 defaults (α=0.01, N=5) are the issue's prescribed conservative starting point; \
         this variant runs with 5× α + 3× cycles to make acceptance #6 (non-polygonal contours) \
         visually evident.".to_string(),
        String::new(),
    ];

    for preset in PRESETS.iter() {
        let (mut cfg, wf) =
            build_64sq_cfg_and_wf(preset, N_CYCLES_AGG, ALPHA_AGG, HD_DROPLETS_AGG, "aggressive");

        // Init S̃ for the patchwork "before" tile.
        let plates = generate_voronoi(64, 64, &VoronoiConfig {
            num_plates: preset.num_plates,
            continental_ratio: preset.continental_ratio,
        }, preset.seed);
        let plate_data = PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        };
        let ctx = InitContext {
            nx: 64, ny: 64, seed: preset.seed,
            amplitude: 0.2, plate_data: Some(plate_data),
        };
        let init_s = init_s_field(InitMode::Checkerboard, &ctx);
        tiles.push((format!("{}_before", preset.label), init_s));

        let phase_a_start = Instant::now();
        let phase_a = run_phase_a_loop(&mut cfg, &wf);
        let phase_a_elapsed = phase_a_start.elapsed().as_secs_f64();

        let last = phase_a.cycles.last().unwrap();
        let after_s = last.baseline.final_state.s_field.clone();
        tiles.push((format!("{}_after", preset.label), after_s.clone()));

        let phase_b_start = Instant::now();
        let phase_b = run_phase_b(&last.baseline.final_state.s_field, &wf, cfg.seed).expect("phase b");
        let phase_b_elapsed = phase_b_start.elapsed().as_secs_f64();

        // Aggregate metrics.
        let cum_drift: f64 = phase_a.cycles.iter().map(|c| c.common.mass_drift).sum();
        let initial_mass_estimate =
            ((preset.continental_ratio + (1.0 - preset.continental_ratio) * 0.2) as f64)
                * (64.0 * 64.0);
        let drift_fraction = cum_drift.abs() / initial_mass_estimate;
        let peaks: Vec<f64> = phase_a
            .cycles
            .iter()
            .map(|c| {
                c.baseline
                    .final_state
                    .s_field
                    .data()
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max)
            })
            .collect();
        let peak_first = peaks[0];
        let peak_last = peaks[peaks.len() - 1];
        let peak_drift_total = (peak_last - peak_first).abs() / peak_first;

        // Save individual artefacts.
        let lowres_path = out_dir.join(format!("{}_phase_a_64sq_aggressive.png", preset.label));
        let lowres_meta = save_heightmap(&after_s, &lowres_path).unwrap();
        let hd_path = out_dir.join(format!("{}_phase_b_hd1024_5m.png", preset.label));
        phase_b.heightmap.save_png_u16(&hd_path).unwrap();
        let sediment_path = out_dir.join(format!("{}_phase_b_sediment.png", preset.label));
        phase_b.sediment.save_png_u16(&sediment_path).unwrap();

        println!(
            "  [{:<16}] Phase A 64² × {} cycles α={} : {:.2} s ({:.2} s/cycle)",
            preset.label, N_CYCLES_AGG, ALPHA_AGG, phase_a_elapsed,
            phase_a_elapsed / N_CYCLES_AGG as f64
        );
        println!(
            "                  cum mass drift = {:+.3} ({:.2} % of initial mass {:.1})",
            cum_drift, drift_fraction * 100.0, initial_mass_estimate
        );
        println!(
            "                  peak S̃ first→last = {:.4} → {:.4} (drift {:.3} %)",
            peak_first, peak_last, peak_drift_total * 100.0
        );
        println!(
            "                  Phase A range=[{:.4}, {:.4}], mean={:.4}",
            lowres_meta.min, lowres_meta.max, lowres_meta.mean
        );
        println!(
            "                  Phase B HD 1024² × {} droplets : {:.2} s",
            HD_DROPLETS_AGG, phase_b_elapsed
        );
        println!(
            "                  D5 metrics: p95 = {:.4} (acc), L_∞ = {:.4} (diag)",
            phase_b.grand_scale_deviation_p95, phase_b.grand_scale_deviation
        );
        println!();

        report_lines.push(format!("## Preset: {}", preset.label));
        report_lines.push(format!(
            "  seed={} num_plates={} cont_ratio={:.2}",
            preset.seed, preset.num_plates, preset.continental_ratio
        ));
        report_lines.push(format!(
            "  Phase A wallclock : {:.2} s ({:.2} s/cycle)",
            phase_a_elapsed, phase_a_elapsed / N_CYCLES_AGG as f64
        ));
        report_lines.push(format!("  Phase B wallclock : {:.2} s", phase_b_elapsed));
        report_lines.push(format!(
            "  Cumulative mass drift : {:+.3} ({:.2} % of initial mass)",
            cum_drift, drift_fraction * 100.0
        ));
        report_lines.push(format!(
            "  Peak S̃ trajectory : first {:.4}, last {:.4}, drift {:.3} %",
            peak_first, peak_last, peak_drift_total * 100.0
        ));
        report_lines.push(format!(
            "  D5 metrics : p95={:.4}, L_∞={:.4}",
            phase_b.grand_scale_deviation_p95, phase_b.grand_scale_deviation
        ));
        report_lines.push(String::new());
    }

    // 2 × 2 patchwork at 64², NN-upscaled ×6 for visibility.
    let tile_w = 64;
    let tile_h = 64;
    let sep = 1usize;
    let pw = tile_w * 2 + sep;
    let ph = tile_h * 2 + sep;
    let mut patch_data = vec![0.5_f64; pw * ph];
    for (k, (_, s)) in tiles.iter().enumerate() {
        let row = k / 2;
        let col = k % 2;
        let x_off = col * (tile_w + sep);
        let y_off = row * (tile_h + sep);
        for j in 0..tile_h {
            for i in 0..tile_w {
                patch_data[(y_off + j) * pw + (x_off + i)] = s.get(i, j);
            }
        }
    }
    let patch = Field2D::from_vec(pw, ph, patch_data);
    let patch_path = out_dir.join("patchwork_phase_a_evolution_64sq_aggressive.png");
    save_heightmap(&patch, &patch_path).expect("save patchwork");
    let patch_up = nn_upscale(&patch, 6);
    let patch_up_path =
        out_dir.join("patchwork_phase_a_evolution_64sq_aggressive_x6.png");
    save_heightmap(&patch_up, &patch_up_path).expect("save patchwork (NN-upscaled)");

    println!();
    println!("Aggressive demo patchwork (NN-upscaled ×6 for visibility):");
    println!("  native : {}", patch_path.display());
    println!("  ×6     : {}", patch_up_path.display());
    println!("  Layout : row 0 = single_continent | row 1 = convergence");
    println!(
        "           col 0 = before workflow  | col 1 = after Phase A ({} cycles α={})",
        N_CYCLES_AGG, ALPHA_AGG
    );

    let report_path = out_dir.join("phase_6_aggressive_metrics.md");
    fs::write(&report_path, report_lines.join("\n")).expect("write report");
    println!("  Metrics report: {}", report_path.display());
}

/// Custom erosion variant supporting per-cell `α(i, j)` field —
/// inlined in this test file because the production
/// `low_res_erosion::apply` API takes a scalar `α`. Used by Variant 3
/// of the curvature investigation. Returns `(volume_removed, peak_delta_h)`.
///
/// Algorithm identical to `low_res_erosion::apply` modulo `α_local =
/// alpha_field.get(i, j)` instead of `params.alpha`.
fn apply_erosion_with_alpha_field(
    s: &mut Field2D,
    alpha_field: &Field2D,
    beta: f64,
    sea_level_ref: f64,
) -> (f64, f64) {
    let nx = s.nx();
    let ny = s.ny();
    let prev_x: Vec<usize> = (0..nx).map(|i| (i + nx - 1) % nx).collect();
    let next_x: Vec<usize> = (0..nx).map(|i| (i + 1) % nx).collect();
    let prev_y: Vec<usize> = (0..ny).map(|j| (j + ny - 1) % ny).collect();
    let next_y: Vec<usize> = (0..ny).map(|j| (j + 1) % ny).collect();

    let n_cells = nx * ny;
    let mut delta_h = vec![0.0_f64; n_cells];
    let mut downslope_lin = vec![0_usize; n_cells];
    {
        let s_data = s.data();
        for j in 0..ny {
            for i in 0..nx {
                let lin = j * nx + i;
                let s_i = s_data[lin];
                if s_i <= sea_level_ref {
                    continue;
                }
                let neighbors = [
                    (i, prev_y[j]),
                    (next_x[i], j),
                    (i, next_y[j]),
                    (prev_x[i], j),
                ];
                let mut max_slope = 0.0_f64;
                let mut best_idx = 0_usize;
                let mut best_h = s_data[neighbors[0].1 * nx + neighbors[0].0];
                for (k, &(ni, nj)) in neighbors.iter().enumerate() {
                    let s_n = s_data[nj * nx + ni];
                    let mag = (s_i - s_n).abs();
                    if mag > max_slope {
                        max_slope = mag;
                    }
                    if k > 0 && s_n < best_h {
                        best_h = s_n;
                        best_idx = k;
                    }
                }
                let alpha_local = alpha_field.get(i, j);
                let dh = alpha_local * max_slope * (s_i - sea_level_ref);
                delta_h[lin] = dh;
                let down = neighbors[best_idx];
                downslope_lin[lin] = down.1 * nx + down.0;
            }
        }
    }
    let mut volume_removed = 0.0_f64;
    let mut peak_dh = 0.0_f64;
    {
        let s_mut = s.data_mut();
        for lin in 0..n_cells {
            let dh = delta_h[lin];
            if dh <= 0.0 {
                continue;
            }
            s_mut[lin] -= dh;
            volume_removed += dh;
            if dh > peak_dh {
                peak_dh = dh;
            }
            if beta > 0.0 {
                s_mut[downslope_lin[lin]] += beta * dh;
            }
        }
    }
    (volume_removed, peak_dh)
}

/// Smooth low-frequency 2D noise field for Variant 3. Two-mode
/// trigonometric superposition with wavelengths ~12-16 cells, so
/// spatial coherence is comparable to the Voronoï border thickness.
/// Returns approximately `[-1, 1]`.
fn build_alpha_noise_field(nx: usize, ny: usize, alpha_base: f64, magnitude: f64) -> Field2D {
    let mut f = Field2D::new(nx, ny);
    let kx = 2.0 * std::f64::consts::PI / 16.0;
    let ky = 2.0 * std::f64::consts::PI / 12.0;
    for j in 0..ny {
        for i in 0..nx {
            let n1 = (i as f64 * kx).sin() * (j as f64 * ky).cos();
            let n2 = ((i as f64 + 7.0) * kx * 1.3 + 0.7).cos()
                * ((j as f64 + 11.0) * ky * 0.9).sin();
            let noise = (n1 + n2) * 0.5;
            f.set(i, j, alpha_base * (1.0 + magnitude * noise));
        }
    }
    f
}

/// Custom Phase A loop for Variant 3 — runs the manual orchestration
/// (tectonic → custom erosion with `alpha_field` → continuation)
/// without the reclassify + craton recompute steps (cratonic is
/// disabled in these tests so they're no-ops anyway).
fn run_custom_phase_a_loop(
    cfg: &mut BaselineConfig,
    n_cycles: usize,
    alpha_field: &Field2D,
    beta: f64,
) -> (Field2D, Vec<f64>) {
    use ymir_core::tectonics_v2::workflow::final_state_to_continuation;
    let mut drifts = Vec::with_capacity(n_cycles);
    let mut last_s: Option<Field2D> = None;
    for cycle_idx in 0..n_cycles {
        let mut baseline =
            ymir_core::tectonics_v2::diagnostics::harness::run_baseline(cfg);
        let mass_before: f64 = baseline.final_state.s_field.data().iter().sum();
        // S̃-space sea_level (Phase 3.5 fix).
        let s_data = baseline.final_state.s_field.data();
        let s_min = s_data.iter().copied().fold(f64::INFINITY, f64::min);
        let s_max = s_data.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let s_range = (s_max - s_min).max(1e-10);
        let sea_level_ref = s_min + 0.4 * s_range;
        apply_erosion_with_alpha_field(
            &mut baseline.final_state.s_field,
            alpha_field,
            beta,
            sea_level_ref,
        );
        let mass_after: f64 = baseline.final_state.s_field.data().iter().sum();
        drifts.push(mass_after - mass_before);
        last_s = Some(baseline.final_state.s_field.clone());
        if cycle_idx + 1 < n_cycles {
            cfg.continuation = Some(final_state_to_continuation(&baseline.final_state));
        }
    }
    (last_s.unwrap(), drifts)
}

/// Variants probe — three configurations to investigate whether D2
/// can be coaxed into producing visible boundary curvature
/// (acceptance #6). All run on 64² Phase A only (no Phase B HD; the
/// HD curvature is already validated in `dump_step12_phase_6_aggressive_demo`).
///
/// - V1: `β = 0.5`, `α = 0.05`, `N = 15` — sediment redistribution
///   on the existing aggressive setup.
/// - V2: `β = 0.5`, `α = 0.02`, `N = 30` — lower α, more cycles.
///   The hypothesis: each cycle's per-cell change is small enough
///   that boundary cells don't all flip simultaneously, allowing
///   spatial pattern to emerge.
/// - V3: `β = 0.0`, `α(i, j) = α_base · (1 + 0.2 · noise)`, `N = 15`
///   — spatial noise on α as a quick hack to break the uniform-
///   retreat symmetry that produces polygonal preservation under
///   D2 vanilla.
#[test]
#[ignore]
fn dump_step12_phase_6_curvature_variants() {
    use std::time::Instant;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_phase_6_curvature_variants");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!();
    println!("Step 12 Phase 6 curvature variants — three D2 modulation paths");
    println!();

    let variants: Vec<(&str, f64, f64, usize, bool)> = vec![
        // (label, alpha, beta, n_cycles, use_custom_alpha_field)
        ("v1_beta05_alpha005_n15", 0.05, 0.5, 15, false),
        ("v2_beta05_alpha002_n30", 0.02, 0.5, 30, false),
        ("v3_alpha_noise_n15",      0.05, 0.0, 15, true),
    ];

    let mut all_lines: Vec<String> = vec![
        "Step 12 Phase 6 — D2 curvature variants probe".to_string(),
        String::new(),
    ];

    for (label, alpha, beta, n_cycles, use_custom) in variants {
        println!("=== Variant: {label} (α={alpha}, β={beta}, N={n_cycles}, custom={use_custom}) ===");
        let variant_start = Instant::now();
        let mut tiles: Vec<(String, Field2D)> = Vec::with_capacity(4);

        for preset in PRESETS.iter() {
            let plates = generate_voronoi(64, 64, &VoronoiConfig {
                num_plates: preset.num_plates,
                continental_ratio: preset.continental_ratio,
            }, preset.seed);
            let plate_data = PlateInitData {
                plate_id: &plates.plate_id,
                plate_type: &plates.plate_type,
                seed_coords: Some(&plates.seed_coords),
            };
            let ctx = InitContext {
                nx: 64, ny: 64, seed: preset.seed,
                amplitude: 0.2, plate_data: Some(plate_data),
            };
            let init_s = init_s_field(InitMode::Checkerboard, &ctx);
            tiles.push((format!("{}_before", preset.label), init_s));

            let (mut cfg, wf) =
                build_64sq_cfg_and_wf(preset, n_cycles, alpha, 5_000_000, &format!("variants_{label}"));

            let after_s: Field2D;
            let cum_drift: f64;
            let preset_t = Instant::now();
            if use_custom {
                let alpha_field = build_alpha_noise_field(64, 64, alpha, 0.2);
                let (s_final, drifts) =
                    run_custom_phase_a_loop(&mut cfg, n_cycles, &alpha_field, beta);
                after_s = s_final;
                cum_drift = drifts.iter().sum();
            } else {
                // Use the production loop; override beta in wf.
                let wf_b = if let WorkflowConfig::Enabled(mut p) = wf {
                    p.phase_a.beta = beta;
                    WorkflowConfig::Enabled(p)
                } else {
                    wf
                };
                let phase_a = run_phase_a_loop(&mut cfg, &wf_b);
                let last = phase_a.cycles.last().unwrap();
                cum_drift = phase_a.cycles.iter().map(|c| c.common.mass_drift).sum();
                after_s = last.baseline.final_state.s_field.clone();
            }
            let preset_elapsed = preset_t.elapsed().as_secs_f64();

            let initial_mass_estimate =
                ((preset.continental_ratio + (1.0 - preset.continental_ratio) * 0.2) as f64)
                    * (64.0 * 64.0);
            let drift_pct = cum_drift.abs() / initial_mass_estimate * 100.0;
            let peak_after: f64 = after_s.data().iter().copied().fold(f64::NEG_INFINITY, f64::max);
            println!(
                "  [{:<16}] {:.1} s, mass drift = {:+.3} ({:.2} % init), peak S̃ = {:.4}",
                preset.label, preset_elapsed, cum_drift, drift_pct, peak_after
            );
            all_lines.push(format!(
                "- {label} / {} : {:.1} s, mass drift = {:+.3} ({:.2} %), peak S̃ = {:.4}",
                preset.label, preset_elapsed, cum_drift, drift_pct, peak_after
            ));
            tiles.push((format!("{}_after", preset.label), after_s));
        }

        // 2 × 2 patchwork for this variant. NN-upscale ×6 for visibility.
        let tile_w = 64;
        let tile_h = 64;
        let sep = 1usize;
        let pw = tile_w * 2 + sep;
        let ph = tile_h * 2 + sep;
        let mut patch_data = vec![0.5_f64; pw * ph];
        for (k, (_, s)) in tiles.iter().enumerate() {
            let row = k / 2;
            let col = k % 2;
            let x_off = col * (tile_w + sep);
            let y_off = row * (tile_h + sep);
            for j in 0..tile_h {
                for i in 0..tile_w {
                    patch_data[(y_off + j) * pw + (x_off + i)] = s.get(i, j);
                }
            }
        }
        let patch = Field2D::from_vec(pw, ph, patch_data);
        let patch_up = nn_upscale(&patch, 6);
        let patch_path = out_dir.join(format!("patchwork_{label}_x6.png"));
        save_heightmap(&patch_up, &patch_path).expect("save patchwork");

        let variant_elapsed = variant_start.elapsed().as_secs_f64();
        println!("  Patchwork: {} ({:.1} s total)", patch_path.display(), variant_elapsed);
        all_lines.push(format!("  Patchwork: {}", patch_path.display()));
        all_lines.push(String::new());
        println!();
    }

    let report_path = out_dir.join("variants_summary.md");
    fs::write(&report_path, all_lines.join("\n")).expect("write report");
    println!("Summary: {}", report_path.display());
}
