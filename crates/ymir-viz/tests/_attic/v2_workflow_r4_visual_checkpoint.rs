//! Step 12 R4 — visual checkpoint, reviewer-validated.
//!
//! Replaces the pre-R3 `v2_workflow_phase_6_visual_checkpoint`
//! (archived behind `#![cfg(any())]` in core/tests). Produces a
//! 2-preset × 6-state × 2-view gallery + a per-cycle metrics table,
//! plus a verdict against the R4.1–R4.5 acceptance criteria.
//!
//! Preset list (per R4 brief): `active_medley` + `single_continent`.
//! Mantle is forced **ON** on both presets — that's the regime that
//! exposed P_aplatissement on the pre-R3 codebase.
//!
//! Workflow params: Step 12 R3 defaults (`α = 0.01`,
//! `isostatic_rebound_ratio = 0.80`, `max_drainage_distance = 10`,
//! `n_cycles = 5`, `k_cycle = 20`) — 100 effective tectonic steps +
//! 5 macro-redistribution passes.
//!
//! Output: `docs/reports/step12_r4_visual_checkpoint/<preset>/`
//! ```
//!   cycle_0_s.png        cycle_0_altitude.png         (INIT, pre-step-1)
//!   cycle_1_s.png        cycle_1_altitude.png         (post-cycle 1)
//!   ...
//!   cycle_5_s.png        cycle_5_altitude.png         (post-cycle 5, final)
//!   metrics.md
//! ```
//!
//! Run with:
//! ```bash
//! cargo test --release -p ymir-viz --test v2_workflow_r4_visual_checkpoint \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! Estimated cost: 64² × 100 steps × mantle-on ≈ 4–8 min per preset
//! (CG-heavy Newton iterations at 64²) → 10–20 min total.

use std::path::{Path, PathBuf};

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::init::{init_s_field, InitContext, PlateInitData};
use ymir_core::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};
use ymir_core::tectonics_v2::workflow::{
    drainage::compute_drainage_targets, run_phase_a_loop_v2,
};

use ymir_viz::bridge::v2::{
    presets, V2FinalState, V2MantleSpec, V2PhaseAParams, V2PhaseBParams, V2WorkflowSpec,
};
use ymir_viz::visualization::v2_viz::{save_field_png, V2Field};

/// Number of cycles to run per preset.
const N_CYCLES: usize = 5;
/// Tectonic steps per cycle.
const K_CYCLE: usize = 20;

fn out_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r4_visual_checkpoint")
}

/// Build a minimal `V2FinalState` from an `init_s_field` result so the
/// V2 raster path can paint the pre-cycle-1 state. `vx/vy/strain_rate`
/// stay zero and the optional fields stay `None` — only S̃ + altitude
/// are needed for the gallery.
fn make_init_state(spec: &ymir_viz::bridge::v2::V2RunSpec) -> V2FinalState {
    let nx = spec.grid_nx;
    let ny = spec.grid_ny;
    let cfg = VoronoiConfig {
        num_plates: spec.num_plates,
        continental_ratio: spec.continental_ratio,
    };
    let plates = generate_voronoi(nx, ny, &cfg, spec.seed);
    let plate_data = PlateInitData {
        plate_id: &plates.plate_id,
        plate_type: &plates.plate_type,
        seed_coords: Some(&plates.seed_coords),
    };
    let init_ctx = InitContext {
        nx,
        ny,
        seed: spec.seed,
        amplitude: spec.s_perturbation_amplitude,
        plate_data: Some(plate_data),
    };
    let s = init_s_field(spec.init_mode.into_core(), &init_ctx);

    V2FinalState {
        nx,
        ny,
        dx: 1.0 / nx as f64,
        dy: 1.0 / ny as f64,
        s_field: s.data().to_vec(),
        vx: vec![0.0; nx * ny],
        vy: vec![0.0; nx * ny],
        strain_rate_invariant: vec![0.0; nx * ny],
        age_field: None,
        cratonic_factor: None,
        plate_id: Some(plates.plate_id.data().to_vec()),
        plate_type: None,
        boundary_flag: None,
    }
}

#[derive(Clone, Copy, Debug)]
struct CycleMetrics {
    cycle: usize,
    peak_s: f64,
    mean_s_continental: f64,
    fraction_above_0_8: f64,
    mass_total: f64,
    mass_continental: f64,
    mass_oceanic: f64,
    sea_level_ref: f64,
    max_path_length: u8,
    total_eroded: f64,
    mass_drift: f64,
}

fn compute_metrics(
    cycle: usize,
    s: &Field2D,
    sea_level_ref: f64,
    total_eroded: f64,
    mass_drift: f64,
) -> CycleMetrics {
    let data = s.data();
    let mut peak = f64::NEG_INFINITY;
    let mut sum_cont = 0.0_f64;
    let mut sum_total = 0.0_f64;
    let mut sum_ocean = 0.0_f64;
    let mut count_cont = 0_usize;
    let mut count_above_0_8 = 0_usize;
    for &v in data {
        if v > peak {
            peak = v;
        }
        sum_total += v;
        if v > sea_level_ref {
            sum_cont += v;
            count_cont += 1;
            if v > 0.8 {
                count_above_0_8 += 1;
            }
        } else {
            sum_ocean += v;
        }
    }
    let mean_cont = if count_cont > 0 {
        sum_cont / count_cont as f64
    } else {
        0.0
    };
    let fraction = if data.len() > 0 {
        count_above_0_8 as f64 / data.len() as f64
    } else {
        0.0
    };

    // Drainage diagnostic — peak path length on the post-cycle field.
    let drainage = compute_drainage_targets(s, sea_level_ref, 15);
    let max_path = *drainage.path_length.iter().max().unwrap_or(&0);

    CycleMetrics {
        cycle,
        peak_s: peak,
        mean_s_continental: mean_cont,
        fraction_above_0_8: fraction,
        mass_total: sum_total,
        mass_continental: sum_cont,
        mass_oceanic: sum_ocean,
        sea_level_ref,
        max_path_length: max_path,
        total_eroded,
        mass_drift,
    }
}

fn run_gallery_for_preset(preset_name: &str, out_dir: &Path) -> Vec<CycleMetrics> {
    use ymir_viz::bridge::v2::build_config;

    std::fs::create_dir_all(out_dir).expect("create output dir");
    println!("\n=== Preset: {preset_name} ===");

    // Load + amend the preset
    let mut spec = presets::load(preset_name)
        .unwrap_or_else(|e| panic!("preset '{preset_name}' load failed: {e}"));
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    // R4 brief: force mantle ON on both presets (the regime that
    // exposed the pre-R3 P_aplatissement). `single_continent` ships
    // with mantle off by default — override.
    if matches!(spec.mantle, V2MantleSpec::Off) {
        // V2MantleSpec::default() is `On { mf=1.0, coupling=1.0,
        // num_modes=6, seed=7, evolution_rate=0.0 }` — the Step 8
        // baseline ON profile. Same as active_medley's mantle block.
        spec.mantle = V2MantleSpec::default();
        println!("[r4] {preset_name}: forced mantle ON (default was Off)");
    }
    // Pin grid + steps to the brief: 64², 5 × 20 = 100 steps.
    spec.grid_nx = 64;
    spec.grid_ny = 64;
    spec.steps = N_CYCLES * K_CYCLE;
    spec.total_time_nondim = 6.0;

    // INIT capture (cycle 0).
    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("cycle_0_s.png"))
        .expect("save init s");
    save_field_png(&init_state, V2Field::Altitude, &out_dir.join("cycle_0_altitude.png"))
        .expect("save init alt");

    // Build configs and run the workflow loop.
    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    println!("[r4] {preset_name}: running 5 cycles × 20 steps (mantle ON, 64²)…");
    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r4] {preset_name}: completed in {elapsed:.1}s, {} cycles", output.cycles.len());

    // INIT metrics (use init field, no erosion stats).
    let init_field = Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let mut metrics = vec![compute_metrics(0, &init_field, init_sea, 0.0, 0.0)];

    // Per-cycle metrics + PNG.
    for (i, cycle) in output.cycles.iter().enumerate() {
        let cyc_idx = i + 1;
        let v2_state = V2FinalState::from_harness(&cycle.baseline.final_state);
        save_field_png(
            &v2_state,
            V2Field::SThickness,
            &out_dir.join(format!("cycle_{cyc_idx}_s.png")),
        )
        .expect("save s");
        save_field_png(
            &v2_state,
            V2Field::Altitude,
            &out_dir.join(format!("cycle_{cyc_idx}_altitude.png")),
        )
        .expect("save altitude");

        metrics.push(compute_metrics(
            cyc_idx,
            &cycle.baseline.final_state.s_field,
            cycle.common.sea_level_normalized,
            cycle.common.erosion_volume_removed,
            cycle.common.mass_drift,
        ));
    }

    metrics
}

fn write_metrics_table(preset_name: &str, metrics: &[CycleMetrics], out_dir: &Path) {
    let mut md = String::new();
    md.push_str(&format!("# Step 12 R4 metrics — {preset_name}\n\n"));
    md.push_str("64² grid, 5 cycles × 20 steps (mantle ON), Phase A defaults: α=0.01, ");
    md.push_str("isostatic_rebound_ratio=0.80, max_drainage_distance=10.\n\n");
    md.push_str("| cycle | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path | total_eroded | mass_drift |\n");
    md.push_str("|---|---|---|---|---|---|---|---|---|---|---|\n");
    for m in metrics {
        md.push_str(&format!(
            "| {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {} | {:.4e} | {:.3e} |\n",
            m.cycle,
            m.peak_s,
            m.mean_s_continental,
            m.fraction_above_0_8,
            m.mass_total,
            m.mass_continental,
            m.mass_oceanic,
            m.sea_level_ref,
            m.max_path_length,
            m.total_eroded,
            m.mass_drift,
        ));
    }
    md.push('\n');

    // Verdict R4.1–R4.5 (per-preset)
    md.push_str("## Verdict R4.1–R4.5\n\n");

    let init = &metrics[0];
    let final_ = metrics.last().unwrap();

    // R4.1 — continents émergés step 100 : peak S̃ final > sea_level + 0.1
    let r4_1 = final_.peak_s > final_.sea_level_ref + 0.1;
    md.push_str(&format!(
        "- **R4.1 — Continents émergés step 100** : peak S̃ final = {:.3} vs sea_level = {:.3} → {}\n",
        final_.peak_s,
        final_.sea_level_ref,
        if r4_1 { "**PASS**" } else { "**FAIL — peak S̃ insufficient above sea level**" }
    ));

    // R4.2 — cratons préservés : fraction (S̃ > 0.8) maintenu > 50 % de l'init
    let r4_2_ratio = if init.fraction_above_0_8 > 0.0 {
        final_.fraction_above_0_8 / init.fraction_above_0_8
    } else {
        0.0
    };
    let r4_2 = r4_2_ratio >= 0.5;
    md.push_str(&format!(
        "- **R4.2 — Cratons préservés (S̃ > 0.8 retained > 50 %)** : init = {:.3}, final = {:.3}, retention = {:.1} % → {}\n",
        init.fraction_above_0_8,
        final_.fraction_above_0_8,
        100.0 * r4_2_ratio,
        if r4_2 { "**PASS**" } else { "**FAIL — too few cratons retained**" }
    ));

    // R4.3 — bordures irrégulières : critère subjectif, signal proxy = max_path_length stable > 1
    // (if drainage stays at 1, no long-distance redistribution → coastal pattern matches init Voronoï)
    let final_max_path = final_.max_path_length;
    let r4_3 = final_max_path >= 2;
    md.push_str(&format!(
        "- **R4.3 — Bordures irrégulières (subjective)** : final max_path_length = {} (proxy ≥ 2 → drainage spans coast cells) → {}\n",
        final_max_path,
        if r4_3 { "**PASS (proxy)**" } else { "**REVIEW VISUAL**" }
    ));

    // R4.4 — conservation totale 5 cycles : drift cumulé < 1e-9 · mass_init
    let cumulative_drift: f64 = metrics
        .iter()
        .skip(1)
        .map(|m| m.mass_drift.abs())
        .sum();
    let budget = init.mass_total.abs() * 1e-9;
    let r4_4 = cumulative_drift < budget;
    md.push_str(&format!(
        "- **R4.4 — Conservation totale (drift < 1e-9 · mass)** : cumulative drift = {:.3e}, budget = {:.3e} → {}\n",
        cumulative_drift,
        budget,
        if r4_4 { "**PASS**" } else { "**FAIL**" }
    ));

    // R4.5 — drainage actif : max_path > 1 sur au moins un cycle post-init
    let max_path_over_cycles = metrics.iter().skip(1).map(|m| m.max_path_length).max().unwrap_or(0);
    let r4_5 = max_path_over_cycles >= 2;
    md.push_str(&format!(
        "- **R4.5 — Drainage actif (max_path > 1)** : max across cycles 1-5 = {} → {}\n",
        max_path_over_cycles,
        if r4_5 { "**PASS**" } else { "**FAIL — drainage limited to immediate neighbour**" }
    ));

    std::fs::write(out_dir.join("metrics.md"), &md).expect("write metrics.md");
    println!("\n{md}");
}

#[test]
#[ignore]
fn v2_workflow_r4_active_medley() {
    let out_dir = out_root().join("active_medley");
    let metrics = run_gallery_for_preset("active_medley", &out_dir);
    write_metrics_table("active_medley", &metrics, &out_dir);
}

#[test]
#[ignore]
fn v2_workflow_r4_single_continent() {
    let out_dir = out_root().join("single_continent");
    let metrics = run_gallery_for_preset("single_continent", &out_dir);
    write_metrics_table("single_continent", &metrics, &out_dir);
}

// ── R4b diagnostic — calibration probes ──────────────────────────────
//
// Step 12 R4 surfaced that the default macro_redistribution params
// (α=0.01, rebound=0.80, max_drainage_distance=10) over 5 × 64²
// active_medley cycles produce a diffuse blue-green patchwork rather
// than identifiable continents — the relief drowns under horizontal
// redistribution. R4b runs four diagnostic configurations to bracket
// which mechanism dominates the diffusion:
//
//   Test 1 — rebound_ratio 0.80 → 0.95: weaker per-cycle redistribution
//   Test 2 — max_drainage_distance 10 → 3: local-only drainage
//   Test 3 — α 0.01 → 0.001: 10× weaker erosion
//   Test 4 — mantle ON → Off: isolate macro_redistribution from
//            tectonic step (the mass-loss source diagnosis)
//
// All run on active_medley at **32²** for tractable runtime (~5-15 min
// each vs 141 min at 64²). Output: `step12_r4b_diagnostic/test_N/`.
//
// Each test captures INIT + final cycle 5 (S̃ + altitude) and metrics.
// The verdict per test reads: "continent reconnaissable" (the gallery
// shows a stable continental shape with relief) vs "diffus" (uniform
// blue-green patchwork like the failed R4 active_medley).

const R4B_GRID: usize = 32;
const R4B_N_CYCLES: usize = 5;
const R4B_K_CYCLE: usize = 20;

fn r4b_out_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r4b_diagnostic")
}

#[derive(Clone, Copy, Debug)]
struct R4bMetrics {
    init: CycleMetrics,
    final_: CycleMetrics,
    total_runtime_s: f64,
    cumulative_mass_drift: f64,
}

fn run_r4b_variant(
    label: &str,
    description: &str,
    out_dir: &Path,
    customize: impl FnOnce(&mut ymir_viz::bridge::v2::V2RunSpec),
) -> R4bMetrics {
    use ymir_viz::bridge::v2::build_config;

    std::fs::create_dir_all(out_dir).expect("create output dir");
    println!("\n=== R4b: {label} — {description} ===");

    // Always base on active_medley (the preset that exposed P_aplatissement),
    // grid 32² for tractable runtime.
    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = R4B_GRID;
    spec.grid_ny = R4B_GRID;
    spec.steps = R4B_N_CYCLES * R4B_K_CYCLE;
    spec.total_time_nondim = 6.0;
    // active_medley already has mantle ON by default; only Test 4
    // (`customize` closure) flips it Off.
    customize(&mut spec);

    // INIT capture
    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("init_s.png"))
        .expect("save init s");
    save_field_png(&init_state, V2Field::Altitude, &out_dir.join("init_altitude.png"))
        .expect("save init alt");

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r4b/{label}] completed in {elapsed:.1}s, {} cycles", output.cycles.len());

    // Final cycle capture
    let final_cycle = output.cycles.last().expect("at least one cycle");
    let v2_final = V2FinalState::from_harness(&final_cycle.baseline.final_state);
    save_field_png(&v2_final, V2Field::SThickness, &out_dir.join("final_s.png"))
        .expect("save final s");
    save_field_png(&v2_final, V2Field::Altitude, &out_dir.join("final_altitude.png"))
        .expect("save final alt");

    let init_field = Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let init_metrics = compute_metrics(0, &init_field, init_sea, 0.0, 0.0);
    let final_metrics = compute_metrics(
        R4B_N_CYCLES,
        &final_cycle.baseline.final_state.s_field,
        final_cycle.common.sea_level_normalized,
        final_cycle.common.erosion_volume_removed,
        final_cycle.common.mass_drift,
    );
    let cumulative_drift: f64 = output.cycles.iter().map(|c| c.common.mass_drift.abs()).sum();

    R4bMetrics {
        init: init_metrics,
        final_: final_metrics,
        total_runtime_s: elapsed,
        cumulative_mass_drift: cumulative_drift,
    }
}

fn write_r4b_report(label: &str, description: &str, metrics: &R4bMetrics, out_dir: &Path) {
    let mut md = String::new();
    md.push_str(&format!("# R4b {label} — {description}\n\n"));
    md.push_str(&format!(
        "32² grid, {} cycles × {} steps = {} steps total. Runtime: {:.1}s.\n\n",
        R4B_N_CYCLES,
        R4B_K_CYCLE,
        R4B_N_CYCLES * R4B_K_CYCLE,
        metrics.total_runtime_s
    ));
    md.push_str("|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path |\n");
    md.push_str("|---|---|---|---|---|---|---|---|---|\n");
    let stamp = |label: &str, m: &CycleMetrics| {
        format!(
            "| {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {} |\n",
            label,
            m.peak_s,
            m.mean_s_continental,
            m.fraction_above_0_8,
            m.mass_total,
            m.mass_continental,
            m.mass_oceanic,
            m.sea_level_ref,
            m.max_path_length,
        )
    };
    md.push_str(&stamp("INIT", &metrics.init));
    md.push_str(&stamp("cycle 5", &metrics.final_));
    md.push('\n');

    let mass_loss = metrics.init.mass_total - metrics.final_.mass_total;
    let mass_loss_pct = 100.0 * mass_loss / metrics.init.mass_total.abs().max(1e-12);
    md.push_str(&format!(
        "Mass loss: {:.3} ({:.1} %).  Cumulative macro_redistribution drift: {:.3e}.\n\n",
        mass_loss, mass_loss_pct, metrics.cumulative_mass_drift
    ));

    // Quantitative signal of "continent stable + relief".
    // frac > 0.8 retention at the cycle-5 absolute level (not ratio,
    // since the init can be 0 for this metric on RadialProfile).
    md.push_str("## Diagnostic\n\n");
    md.push_str(&format!(
        "- final `frac S̃>0.8` = **{:.3}** (proxy for cratonic patches: ≥ 0.05 → patches visible, < 0.02 → flattened)\n",
        metrics.final_.fraction_above_0_8
    ));
    md.push_str(&format!(
        "- mass loss {:.1} % over 5 cycles ({} non-conservative if > 5 %)\n",
        mass_loss_pct,
        if mass_loss_pct > 5.0 { "**>5%**" } else { "≤ 5%" }
    ));
    md.push_str(&format!(
        "- max drainage path (final): **{}** (≤ 3 → local, ≥ 6 → long-distance)\n",
        metrics.final_.max_path_length
    ));
    md.push_str("\n**Verdict subjectif** : voir le patchwork (`init_*.png` vs `final_*.png`) — la métrique scalaire ne tranche pas seule.\n");

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r4b_test_1_high_rebound() {
    let out_dir = r4b_out_root().join("test_1_high_rebound");
    let metrics = run_r4b_variant(
        "test_1_high_rebound",
        "rebound 0.95 (weaker redistribution)",
        &out_dir,
        |spec| {
            if let V2WorkflowSpec::On { phase_a, .. } = &mut spec.workflow {
                phase_a.isostatic_rebound_ratio = 0.95;
            }
        },
    );
    write_r4b_report(
        "test_1_high_rebound",
        "rebound 0.95 (weaker redistribution)",
        &metrics,
        &out_dir,
    );
}

#[test]
#[ignore]
fn r4b_test_2_short_drainage() {
    let out_dir = r4b_out_root().join("test_2_short_drainage");
    let metrics = run_r4b_variant(
        "test_2_short_drainage",
        "max_drainage_distance 3 (local drainage)",
        &out_dir,
        |spec| {
            if let V2WorkflowSpec::On { phase_a, .. } = &mut spec.workflow {
                phase_a.max_drainage_distance = 3;
            }
        },
    );
    write_r4b_report(
        "test_2_short_drainage",
        "max_drainage_distance 3 (local drainage)",
        &metrics,
        &out_dir,
    );
}

#[test]
#[ignore]
fn r4b_test_3_low_alpha() {
    let out_dir = r4b_out_root().join("test_3_low_alpha");
    let metrics = run_r4b_variant(
        "test_3_low_alpha",
        "alpha 0.001 (10x weaker erosion)",
        &out_dir,
        |spec| {
            if let V2WorkflowSpec::On { phase_a, .. } = &mut spec.workflow {
                phase_a.alpha = 0.001;
            }
        },
    );
    write_r4b_report(
        "test_3_low_alpha",
        "alpha 0.001 (10x weaker erosion)",
        &metrics,
        &out_dir,
    );
}

#[test]
#[ignore]
fn r4b_test_4_tectonic_off() {
    let out_dir = r4b_out_root().join("test_4_tectonic_off");
    let metrics = run_r4b_variant(
        "test_4_tectonic_off",
        "mantle OFF — isolate macro_redistribution from tectonic driver",
        &out_dir,
        |spec| {
            spec.mantle = V2MantleSpec::Off;
        },
    );
    write_r4b_report(
        "test_4_tectonic_off",
        "mantle OFF — isolate macro_redistribution from tectonic driver",
        &metrics,
        &out_dir,
    );
}

// ── R4b.5 — mantle rate sweep ────────────────────────────────────────
//
// Tests 1-3 (rebound, drainage, alpha overrides) all produced near-
// identical results: 14.7 % mass loss, frac S̃ > 0.8 ≈ 0.008, max_path
// = 1. Test 4 (mantle OFF) preserved continents perfectly (0.1 % mass
// loss, frac S̃ > 0.8 = 0.353 unchanged). The R4b verdict pinpoints
// the tectonic + macro_redistribution interaction as the destruction
// vector — specifically the mantle convection rate.
//
// R4b.5 sweeps `V2MantleSpec::On.mf` (mantle-flow magnitude) at
// {1.0, 0.5, 0.1} on the active_medley base with macro_redistribution
// defaults. Hypothesis: does an intermediate `mf` produce
// destruction-rate ≈ reconstruction-rate equilibrium?
//
// Output: `docs/reports/step12_r4b5_mantle_sweep/mf_<value>/`.

fn run_r4b5_mantle_sweep(mf_value: f64, label: &str, description: &str) {
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r4b5_mantle_sweep")
        .join(label);
    let metrics = run_r4b_variant(label, description, &out_dir, |spec| {
        // Preserve coupling / num_modes / seed / evolution_rate from
        // the active_medley default; override `mf` only.
        if let V2MantleSpec::On {
            coupling,
            num_modes,
            seed,
            evolution_rate,
            ..
        } = spec.mantle
        {
            spec.mantle = V2MantleSpec::On {
                mf: mf_value,
                coupling,
                num_modes,
                seed,
                evolution_rate,
            };
        }
    });
    write_r4b_report(label, description, &metrics, &out_dir);
}

#[test]
#[ignore]
fn r4b5_mf_1_0() {
    run_r4b5_mantle_sweep(
        1.0,
        "mf_1_0",
        "mantle flow magnitude mf=1.0 (default baseline, macro defaults)",
    );
}

#[test]
#[ignore]
fn r4b5_mf_0_5() {
    run_r4b5_mantle_sweep(
        0.5,
        "mf_0_5",
        "mantle flow magnitude mf=0.5 (half intensity)",
    );
}

#[test]
#[ignore]
fn r4b5_mf_0_1() {
    run_r4b5_mantle_sweep(
        0.1,
        "mf_0_1",
        "mantle flow magnitude mf=0.1 (very low intensity)",
    );
}

// ── R5b sweep — mf {0.5, 0.6, 0.7, 0.8, 0.9} × full D2+D1-ter ─────
//
// R4b.5 mf sweep was conducted BEFORE the D2 + D1-ter solver fixes;
// the verdict (mf=0.5 = sweet spot at 32²) may differ now that the
// solver converges properly across the regime. This sweep re-measures
// 5 mf values on 32² × 5 cycles × 20 steps × workflow ON × D2+D1-ter,
// capturing **peak |v| per cycle** explicitly (the metric missing
// from R4b.5 that hid the quasi-static failure of mf=0.5 at 64²).
//
// Output: `docs/reports/step12_r5b_mf_sweep_post_d1_ter/<label>/`.
//
// Multi-dim acceptance per axis (cf. user brief) :
//   1. Preservation     — frac S̃>0.8 retention > 50 %
//   2. Dynamics         — peak |v| > 0.1 nondim on ≥ 3 of 5 cycles
//   3. Stability        — mass loss cumul < 5 %
//   4. Convergence      — Stalled/Capped < 20 % of steps
//   5. Visual           — continent shape persists (out of scope of
//                         this scalar test — judged on cycle_5 PNG)

fn run_r5b_mf_sweep(mf_value: f64, label: &str) {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_mf_sweep_post_d1_ter")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!("\n=== R5b mf sweep: {label} (mf={mf_value}) — 32² × 5 cycles × D2+D1-ter ===");

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    if let V2MantleSpec::On {
        coupling,
        num_modes,
        seed,
        evolution_rate,
        ..
    } = spec.mantle
    {
        spec.mantle = V2MantleSpec::On {
            mf: mf_value,
            coupling,
            num_modes,
            seed,
            evolution_rate,
        };
    }
    spec.grid_nx = R4B_GRID;
    spec.grid_ny = R4B_GRID;
    spec.steps = R4B_N_CYCLES * R4B_K_CYCLE;
    spec.total_time_nondim = 6.0;

    // INIT capture
    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("cycle_0_s.png"))
        .expect("save init s");
    save_field_png(&init_state, V2Field::Altitude, &out_dir.join("cycle_0_altitude.png"))
        .expect("save init alt");

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r5b-sweep/{label}] completed in {elapsed:.1}s");

    let init_field = Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let mut metrics = vec![compute_metrics(0, &init_field, init_sea, 0.0, 0.0)];

    for (i, cycle) in output.cycles.iter().enumerate() {
        let cyc_idx = i + 1;
        let v2_state = V2FinalState::from_harness(&cycle.baseline.final_state);
        save_field_png(
            &v2_state,
            V2Field::SThickness,
            &out_dir.join(format!("cycle_{cyc_idx}_s.png")),
        )
        .expect("save s");
        save_field_png(
            &v2_state,
            V2Field::Altitude,
            &out_dir.join(format!("cycle_{cyc_idx}_altitude.png")),
        )
        .expect("save altitude");
        metrics.push(compute_metrics(
            cyc_idx,
            &cycle.baseline.final_state.s_field,
            cycle.common.sea_level_normalized,
            cycle.common.erosion_volume_removed,
            cycle.common.mass_drift,
        ));
    }

    let peak_v: Vec<f64> = output
        .cycles
        .iter()
        .map(|c| c.baseline.metrics.vmax_peak)
        .collect();
    let newton_outcomes: Vec<(usize, usize, usize, usize)> = output
        .cycles
        .iter()
        .map(|c| {
            let n = c.baseline.metrics.newton.as_ref();
            (
                n.map(|n| n.converged).unwrap_or(0),
                n.map(|n| n.stalled).unwrap_or(0),
                n.map(|n| n.diverged).unwrap_or(0),
                n.map(|n| n.capped).unwrap_or(0),
            )
        })
        .collect();

    let mut md = String::new();
    md.push_str(&format!("# R5b mf sweep — {label} (mf={mf_value})\n\n"));
    md.push_str(&format!(
        "32² active_medley, workflow ON (D2 + D1-ter), 5 cycles × 20 steps. Runtime: {:.1}s.\n\n",
        elapsed
    ));

    md.push_str("## Per-cycle solver health\n\n");
    md.push_str("| cycle | peak \\|v\\| | Newton C/S/D/Cap | peak S̃ | frac>0.8 | mass total | max_path |\n");
    md.push_str("|---|---|---|---|---|---|---|\n");
    for (i, cycle) in output.cycles.iter().enumerate() {
        let m = &metrics[i + 1];
        let outc = newton_outcomes[i];
        md.push_str(&format!(
            "| {} | {:.3e} | {}/{}/{}/{} | {:.3} | {:.3} | {:.2} | {} |\n",
            i + 1,
            peak_v[i],
            outc.0,
            outc.1,
            outc.2,
            outc.3,
            m.peak_s,
            m.fraction_above_0_8,
            m.mass_total,
            m.max_path_length,
        ));
    }
    md.push('\n');

    // Multi-dim verdict
    let init_frac = metrics[0].fraction_above_0_8;
    let final_frac = metrics.last().unwrap().fraction_above_0_8;
    let retention = if init_frac > 1e-9 { final_frac / init_frac } else { 0.0 };
    let n_cycles_dynamic = peak_v.iter().filter(|&&v| v > 0.1).count();
    let mass_init = metrics[0].mass_total;
    let mass_final = metrics.last().unwrap().mass_total;
    let mass_loss_pct = (mass_init - mass_final) / mass_init.abs().max(1e-12) * 100.0;
    let total_steps_all_cycles: usize = newton_outcomes
        .iter()
        .map(|(c, s, d, cap)| c + s + d + cap)
        .sum();
    let stalled_or_capped: usize = newton_outcomes
        .iter()
        .map(|(_c, s, _d, cap)| s + cap)
        .sum();
    let stall_cap_pct = if total_steps_all_cycles > 0 {
        100.0 * stalled_or_capped as f64 / total_steps_all_cycles as f64
    } else {
        0.0
    };

    md.push_str("## Multi-dim acceptance\n\n");
    md.push_str(&format!(
        "1. **Preservation** : frac>0.8 retention = {:.1} % (init {:.3} → final {:.3}) → {}\n",
        100.0 * retention,
        init_frac,
        final_frac,
        if retention > 0.5 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str(&format!(
        "2. **Dynamics** : peak |v| > 0.1 on {}/5 cycles → {}\n",
        n_cycles_dynamic,
        if n_cycles_dynamic >= 3 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str(&format!(
        "3. **Stability** : mass loss = {:.2} % → {}\n",
        mass_loss_pct,
        if mass_loss_pct.abs() < 5.0 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str(&format!(
        "4. **Convergence** : Stalled+Capped = {:.1} % of steps → {}\n",
        stall_cap_pct,
        if stall_cap_pct < 20.0 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str("5. **Visual** : inspect `cycle_5_altitude.png` vs `cycle_0_altitude.png` (judged externally)\n");

    std::fs::write(out_dir.join("report.md"), &md).expect("write report");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r5b_sweep_mf_0_5() {
    run_r5b_mf_sweep(0.5, "mf_0_5");
}

#[test]
#[ignore]
fn r5b_sweep_mf_0_6() {
    run_r5b_mf_sweep(0.6, "mf_0_6");
}

#[test]
#[ignore]
fn r5b_sweep_mf_0_7() {
    run_r5b_mf_sweep(0.7, "mf_0_7");
}

#[test]
#[ignore]
fn r5b_sweep_mf_0_8() {
    run_r5b_mf_sweep(0.8, "mf_0_8");
}

#[test]
#[ignore]
fn r5b_sweep_mf_0_9() {
    run_r5b_mf_sweep(0.9, "mf_0_9");
}

// ── R5b evolution_rate sweep — γ.1 investigation upstream ──────────
//
// R5b mf sweep verdict: aucun mf ∈ {0.5, 0.6, 0.7, 0.8, 0.9} ne
// passe les 3 critères (preservation, dynamics, stability)
// simultanément. Pattern mf=0.7/0.8 révélateur : cycle 1 explosion
// (peak |v| ~ 0.3-1) puis amortissement cycles 2-5. Le mantle field
// avec `evolution_rate = 0.0` (default `active_medley`) est figé,
// seulement scalé par mf — l'énergie n'est pas réinjectée
// cycle-à-cycle.
//
// (γ.1) Hypothesis : un `mantle.evolution_rate > 0` réinjecte
// temporellement l'énergie convective et maintient peak |v| > 0.1
// sur les 5 cycles. Sweet spot possible en combinant `mf` modéré
// (0.7 — assez pour produire la dynamique cycle 1) et evolution_rate
// non-nul (sustenance).
//
// Sweep : `evolution_rate ∈ {0.05, 0.10, 0.20}` avec `mf = 0.7` fixe.
//
// Output: `docs/reports/step12_r5b_evo_sweep_mf_0_7/<label>/`.

fn run_r5b_evo_sweep(mf_value: f64, evo_rate: f64, label: &str) {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_evo_sweep_mf_0_7")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!(
        "\n=== R5b evo sweep: {label} (mf={mf_value}, evolution_rate={evo_rate}) ==="
    );

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    if let V2MantleSpec::On { coupling, num_modes, seed, .. } = spec.mantle {
        spec.mantle = V2MantleSpec::On {
            mf: mf_value,
            coupling,
            num_modes,
            seed,
            evolution_rate: evo_rate,
        };
    }
    spec.grid_nx = R4B_GRID;
    spec.grid_ny = R4B_GRID;
    spec.steps = R4B_N_CYCLES * R4B_K_CYCLE;
    spec.total_time_nondim = 6.0;

    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("cycle_0_s.png"))
        .expect("save init s");
    save_field_png(&init_state, V2Field::Altitude, &out_dir.join("cycle_0_altitude.png"))
        .expect("save init alt");

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r5b-evo/{label}] completed in {elapsed:.1}s");

    let init_field =
        Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let mut metrics = vec![compute_metrics(0, &init_field, init_sea, 0.0, 0.0)];
    for (i, cycle) in output.cycles.iter().enumerate() {
        let cyc_idx = i + 1;
        let v2_state = V2FinalState::from_harness(&cycle.baseline.final_state);
        save_field_png(&v2_state, V2Field::SThickness, &out_dir.join(format!("cycle_{cyc_idx}_s.png"))).expect("save s");
        save_field_png(&v2_state, V2Field::Altitude, &out_dir.join(format!("cycle_{cyc_idx}_altitude.png"))).expect("save alt");
        metrics.push(compute_metrics(
            cyc_idx,
            &cycle.baseline.final_state.s_field,
            cycle.common.sea_level_normalized,
            cycle.common.erosion_volume_removed,
            cycle.common.mass_drift,
        ));
    }

    let peak_v: Vec<f64> = output.cycles.iter().map(|c| c.baseline.metrics.vmax_peak).collect();
    let newton_outcomes: Vec<(usize, usize, usize, usize)> = output
        .cycles
        .iter()
        .map(|c| {
            let n = c.baseline.metrics.newton.as_ref();
            (
                n.map(|n| n.converged).unwrap_or(0),
                n.map(|n| n.stalled).unwrap_or(0),
                n.map(|n| n.diverged).unwrap_or(0),
                n.map(|n| n.capped).unwrap_or(0),
            )
        })
        .collect();

    let mut md = String::new();
    md.push_str(&format!(
        "# R5b evolution_rate sweep — {label} (mf={mf_value}, evolution_rate={evo_rate})\n\n"
    ));
    md.push_str(&format!(
        "32² active_medley, workflow ON (D2 + D1-ter), 5 cycles × 20 steps. Runtime: {:.1}s.\n\n",
        elapsed
    ));
    md.push_str("## Per-cycle solver health\n\n");
    md.push_str("| cycle | peak \\|v\\| | Newton C/S/D/Cap | peak S̃ | frac>0.8 | mass total | max_path |\n");
    md.push_str("|---|---|---|---|---|---|---|\n");
    for (i, _cycle) in output.cycles.iter().enumerate() {
        let m = &metrics[i + 1];
        let o = newton_outcomes[i];
        md.push_str(&format!(
            "| {} | {:.3e} | {}/{}/{}/{} | {:.3} | {:.3} | {:.2} | {} |\n",
            i + 1, peak_v[i], o.0, o.1, o.2, o.3, m.peak_s, m.fraction_above_0_8, m.mass_total, m.max_path_length,
        ));
    }
    md.push('\n');

    let init_frac = metrics[0].fraction_above_0_8;
    let final_frac = metrics.last().unwrap().fraction_above_0_8;
    let retention = if init_frac > 1e-9 { final_frac / init_frac } else { 0.0 };
    let n_dynamic = peak_v.iter().filter(|&&v| v > 0.1).count();
    let mass_init = metrics[0].mass_total;
    let mass_final = metrics.last().unwrap().mass_total;
    let mass_loss_pct = (mass_init - mass_final) / mass_init.abs().max(1e-12) * 100.0;
    let total_steps: usize = newton_outcomes.iter().map(|(c, s, d, cap)| c + s + d + cap).sum();
    let stall_cap: usize = newton_outcomes.iter().map(|(_c, s, _d, cap)| s + cap).sum();
    let stall_pct = if total_steps > 0 { 100.0 * stall_cap as f64 / total_steps as f64 } else { 0.0 };

    md.push_str("## Multi-dim acceptance\n\n");
    md.push_str(&format!(
        "1. **Preservation** : frac>0.8 retention = {:.1} % → {}\n",
        100.0 * retention,
        if retention > 0.5 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str(&format!(
        "2. **Dynamics (sustained)** : peak |v| > 0.1 on {}/5 cycles → {}\n",
        n_dynamic,
        if n_dynamic >= 3 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str(&format!(
        "3. **Stability** : mass loss = {:.2} % → {}\n",
        mass_loss_pct,
        if mass_loss_pct.abs() < 5.0 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str(&format!(
        "4. **Convergence** : Stalled+Capped = {:.1} % → {}\n",
        stall_pct,
        if stall_pct < 20.0 { "**PASS**" } else { "**FAIL**" }
    ));
    md.push_str("5. **Visual** : inspect `cycle_5_altitude.png`\n");

    std::fs::write(out_dir.join("report.md"), &md).expect("write report");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r5b_evo_mf_0_7_evo_0_05() {
    run_r5b_evo_sweep(0.7, 0.05, "evo_0_05");
}

#[test]
#[ignore]
fn r5b_evo_mf_0_7_evo_0_10() {
    run_r5b_evo_sweep(0.7, 0.10, "evo_0_10");
}

#[test]
#[ignore]
fn r5b_evo_mf_0_7_evo_0_20() {
    run_r5b_evo_sweep(0.7, 0.20, "evo_0_20");
}

// ── R5.0 — dt sweep at mf=1.0 ────────────────────────────────────────
//
// R4b.5 showed that mf=1.0 runtime is 12× longer than mf=0.5 for a
// 2× factor on the mantle flow magnitude — signature of non-linear
// solver saturation (Newton iter and CG iter blowing up beyond the
// design regime). Before accepting "mf=1.0 destroys continents" as a
// physical finding, R5.0 sweeps dt at mf=1.0 to test whether a
// smaller dt keeps the solver in its convergent regime (Newton ~5
// iter, CG ~50 iter) and yields a different physical outcome.
//
// Same total simulated time (`total_time_nondim = 6.0`); only the
// step count and per-step dt vary:
//
//   baseline : 5 cycles × k_cycle=20  = 100 steps, dt = 0.06
//   dt/2     : 5 cycles × k_cycle=40  = 200 steps, dt = 0.03
//   dt/4     : 5 cycles × k_cycle=80  = 400 steps, dt = 0.015
//
// Output: `docs/reports/step12_r5_dt_sweep/<label>/`. Each report
// adds solver-health metrics (Newton iter mean, CG iter mean,
// kappa estimate) so the dt sweet spot — if any — is identifiable
// from numerical signatures, not just visual outcome.

#[derive(Clone, Copy, Debug)]
struct R5SolverMetrics {
    runtime_s: f64,
    total_steps: usize,
    cg_iter_mean: f64,
    cg_iter_max: usize,
    kappa_estimate: f64,
    vmax_peak: f64,
}

fn run_r5_dt_variant(label: &str, k_cycle: usize) -> (R4bMetrics, R5SolverMetrics) {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5_dt_sweep")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");
    println!(
        "\n=== R5.0: {label} — k_cycle={}, total_steps={}, dt= 6.0/{} ===",
        k_cycle,
        5 * k_cycle,
        5 * k_cycle
    );

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams {
            n_cycles: 5,
            k_cycle,
            ..V2PhaseAParams::default()
        },
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = R4B_GRID;
    spec.grid_ny = R4B_GRID;
    spec.steps = 5 * k_cycle;
    spec.total_time_nondim = 6.0; // fixed → dt = total_time / steps

    // INIT capture (mf=1.0 default mantle from active_medley)
    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("init_s.png"))
        .expect("save init s");
    save_field_png(&init_state, V2Field::Altitude, &out_dir.join("init_altitude.png"))
        .expect("save init alt");

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r5/{label}] completed in {elapsed:.1}s, {} cycles", output.cycles.len());

    // Final cycle PNG + metrics
    let final_cycle = output.cycles.last().expect("at least one cycle");
    let v2_final = V2FinalState::from_harness(&final_cycle.baseline.final_state);
    save_field_png(&v2_final, V2Field::SThickness, &out_dir.join("final_s.png"))
        .expect("save final s");
    save_field_png(&v2_final, V2Field::Altitude, &out_dir.join("final_altitude.png"))
        .expect("save final alt");

    let init_field = Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let init_metrics = compute_metrics(0, &init_field, init_sea, 0.0, 0.0);
    let final_metrics = compute_metrics(
        5,
        &final_cycle.baseline.final_state.s_field,
        final_cycle.common.sea_level_normalized,
        final_cycle.common.erosion_volume_removed,
        final_cycle.common.mass_drift,
    );
    let cumulative_drift: f64 = output.cycles.iter().map(|c| c.common.mass_drift.abs()).sum();

    let r4b = R4bMetrics {
        init: init_metrics,
        final_: final_metrics,
        total_runtime_s: elapsed,
        cumulative_mass_drift: cumulative_drift,
    };

    // Solver health — average over the 5 baseline runs.
    let n = output.cycles.len() as f64;
    let cg_iter_mean = output.cycles.iter().map(|c| c.baseline.metrics.cg_iter_mean).sum::<f64>() / n;
    let cg_iter_max = output
        .cycles
        .iter()
        .map(|c| c.baseline.metrics.cg_iter_max)
        .max()
        .unwrap_or(0);
    let kappa = output.cycles.iter().map(|c| c.baseline.metrics.kappa_estimate).sum::<f64>() / n;
    let vmax_peak = output
        .cycles
        .iter()
        .map(|c| c.baseline.metrics.vmax_peak)
        .fold(f64::NEG_INFINITY, f64::max);

    let solver = R5SolverMetrics {
        runtime_s: elapsed,
        total_steps: 5 * k_cycle,
        cg_iter_mean,
        cg_iter_max,
        kappa_estimate: kappa,
        vmax_peak,
    };

    write_r5_report(label, k_cycle, &r4b, &solver, &out_dir);

    (r4b, solver)
}

fn write_r5_report(
    label: &str,
    k_cycle: usize,
    m: &R4bMetrics,
    s: &R5SolverMetrics,
    out_dir: &Path,
) {
    let dt = 6.0 / (5.0 * k_cycle as f64);
    let mut md = String::new();
    md.push_str(&format!("# R5.0 {label}\n\n"));
    md.push_str(&format!(
        "32² active_medley (mf=1.0 default), workflow ON. k_cycle={}, total_steps={}, dt={:.4}.\n\n",
        k_cycle, s.total_steps, dt
    ));
    md.push_str("## Solver health (averaged over 5 cycles)\n\n");
    md.push_str(&format!(
        "- Runtime total: **{:.1}s**  (per step: {:.3}s)\n",
        s.runtime_s,
        s.runtime_s / s.total_steps as f64
    ));
    md.push_str(&format!(
        "- CG iter mean: **{:.1}**  (max: {})\n",
        s.cg_iter_mean, s.cg_iter_max
    ));
    md.push_str(&format!("- Kappa estimate: {:.2e}\n", s.kappa_estimate));
    md.push_str(&format!("- Peak |v|: {:.3e}\n\n", s.vmax_peak));

    md.push_str("## Physical state\n\n");
    md.push_str("|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | max_path |\n");
    md.push_str("|---|---|---|---|---|---|---|\n");
    let row = |label: &str, c: &CycleMetrics| {
        format!(
            "| {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {} |\n",
            label, c.peak_s, c.mean_s_continental, c.fraction_above_0_8, c.mass_total, c.mass_continental, c.max_path_length
        )
    };
    md.push_str(&row("INIT", &m.init));
    md.push_str(&row("cycle 5", &m.final_));

    let mass_loss = m.init.mass_total - m.final_.mass_total;
    let mass_loss_pct = 100.0 * mass_loss / m.init.mass_total.abs().max(1e-12);
    md.push_str(&format!(
        "\nMass loss: {:.3} ({:.1} %). Cumulative macro drift: {:.3e}.\n",
        mass_loss, mass_loss_pct, m.cumulative_mass_drift
    ));

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r5_0_dt_baseline() {
    // k_cycle=20 → 100 steps, dt=0.06 (matches R4b.5 mf_1_0 setup)
    run_r5_dt_variant("dt_baseline_k20", 20);
}

#[test]
#[ignore]
fn r5_0_dt_half() {
    // k_cycle=40 → 200 steps, dt=0.03 (dt/2 from baseline)
    run_r5_dt_variant("dt_half_k40", 40);
}

#[test]
#[ignore]
fn r5_0_dt_quarter() {
    // k_cycle=80 → 400 steps, dt=0.015 (dt/4 from baseline)
    run_r5_dt_variant("dt_quarter_k80", 80);
}

#[test]
#[ignore]
fn r5_0_dt_tenth() {
    // k_cycle=200 → 1000 steps, dt=0.006 (dt/10 from baseline).
    // Run conditionally: only useful if dt/4 hints at solver saturation
    // (Newton/CG metrics drop sharply between baseline and dt/4).
    run_r5_dt_variant("dt_tenth_k200", 200);
}

// ── R5.0.1 — short-run dt sweep to locate solver desaturation ─────────
//
// R5.0 found CG iter cap (2000) saturated on every dt ∈ {0.06, 0.03,
// 0.015} configuration — none in the convergent regime. R5.0.1 runs
// short 50-step probes at much smaller dt to LOCATE the desaturation
// threshold without paying the full simulation runtime.
//
// Config A: 50 steps × dt=0.006 (total_time_nondim=0.3, ≈9 Ma) — same
//           dt as R5.0 dt_tenth but only 1 cycle / 50 steps so the
//           run finishes in ~2-5 min.
// Config B: 50 steps × dt=0.0006 (total_time_nondim=0.03) — extreme
//           small dt to bracket the saturation if Config A still hits
//           the cap.
//
// These are NOT physical-state validation runs (total simulated time
// is too short for macro_redistribution to fully act). The point is
// to read `cg_iter_mean` and `cg_iter_max` and decide whether a
// convergent regime exists at all for mf=1.0.
//
// Output: `docs/reports/step12_r5_0_1_short_dt/<label>/`.

fn run_r5_0_1_short_dt(total_time: f64, label: &str) {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5_0_1_short_dt")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    let n_cycles: usize = 1;
    let k_cycle: usize = 50;
    let total_steps = n_cycles * k_cycle;
    let dt = total_time / total_steps as f64;
    println!(
        "\n=== R5.0.1: {label} — total_steps={}, total_time={:.4}, dt={:.6} ===",
        total_steps, total_time, dt
    );

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams {
            n_cycles,
            k_cycle,
            ..V2PhaseAParams::default()
        },
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = R4B_GRID;
    spec.grid_ny = R4B_GRID;
    spec.steps = total_steps;
    spec.total_time_nondim = total_time;

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();

    let cycle = output.cycles.first().expect("at least one cycle");
    let m = &cycle.baseline.metrics;
    println!(
        "[r5.0.1/{label}] {:.1}s, CG iter mean={:.1} max={}, kappa={:.2e}, peak|v|={:.3e}",
        elapsed, m.cg_iter_mean, m.cg_iter_max, m.kappa_estimate, m.vmax_peak
    );

    let mut md = String::new();
    md.push_str(&format!("# R5.0.1 {label}\n\n"));
    md.push_str(&format!(
        "32² active_medley (mf=1.0), workflow ON, n_cycles={}, k_cycle={}, total_steps={}, total_time={}, dt={:.6}.\n\n",
        n_cycles, k_cycle, total_steps, total_time, dt
    ));
    md.push_str("## Solver health\n\n");
    md.push_str(&format!("- Runtime: **{:.1}s** ({:.3}s/step)\n", elapsed, elapsed / total_steps as f64));
    md.push_str(&format!("- CG iter mean: **{:.1}**  (max: {})\n", m.cg_iter_mean, m.cg_iter_max));
    md.push_str(&format!("- Kappa estimate: {:.2e}\n", m.kappa_estimate));
    md.push_str(&format!("- Peak |v|: {:.3e}\n", m.vmax_peak));
    md.push_str(&format!("- Mass drift over 1 cycle: {:.3e}\n", cycle.common.mass_drift));
    md.push_str(&format!("- Mass drift cumulative (1 cycle = same): {:.3e}\n\n", cycle.common.mass_drift.abs()));

    let cap_2000 = m.cg_iter_max >= 2000;
    let nominal = m.cg_iter_mean < 500.0;
    let regime = if cap_2000 && !nominal {
        "**SATURATED** (CG max = 2000 cap, mean not converged)"
    } else if nominal {
        "**CONVERGENT** (CG mean < 500)"
    } else {
        "**TRANSITIONAL** (cap not hit but mean still high)"
    };
    md.push_str(&format!("\n## Regime\n\n{}\n", regime));

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r5_0_1_config_a_dt_006() {
    // 50 steps × dt=0.006 → total_time=0.3 (~ 9 Ma at Earth scale)
    run_r5_0_1_short_dt(0.3, "config_a_dt_0_006");
}

#[test]
#[ignore]
fn r5_0_1_config_b_dt_0006() {
    // 50 steps × dt=0.0006 → total_time=0.03 (~ 0.9 Ma)
    run_r5_0_1_short_dt(0.03, "config_b_dt_0_0006");
}

// ── R5b D1 — AMG vs JacobiCG benchmark on mf=1.0 regime ─────────────
//
// Objective (per D0 audit): measure the gain of the AMG preconditioner
// (2812 lines, opt-in, never benchmarked in workflow regime) versus
// the saturated JacobiCG default, on the same exact 10-step × 32² ×
// mf=1.0 setup. Workflow OFF — focus on the linear solver in
// isolation, not on macro_redistribution interaction.
//
// Three possible outcomes (per user decision tree):
//   1. AMG saturates too (CG > 1500 at 10 steps): historical note
//      "AMG ineffective at 32²" confirmed. D2 reintegration of #49
//      state/oscillation criteria becomes the primary fix.
//   2. AMG helps but still saturates (CG ≈ 500-1500): combine AMG +
//      D2. AMG opt-in for workflow regimes only.
//   3. AMG desaturates alone (CG < 500): surprise — historical note
//      based on Poisson tests, not coupled Stokes in this regime.
//      D2 still useful as defense-in-depth.
//
// Implementation note: `bridge::v2::build_config` has a known TODO
// (build_config.rs:108) mapping `V2LinearSolverSpec::Amg` to
// `LinearSolverConfig::default()` (= JacobiCG). To benchmark AMG
// we bypass the preset → build_config pipeline and construct
// BaselineConfig directly with the chosen LinearSolverConfig.

fn run_r5b_d1_solver(
    label: &str,
    linear_solver: ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig,
) {
    use std::path::PathBuf;
    use ymir_core::tectonics_v2::diagnostics::harness::{
        run_baseline, BaselineConfig,
    };
    use ymir_core::tectonics_v2::mantle::MantleConfig;
    use ymir_core::tectonics_v2::scales::Scales;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_d1_solver_bench")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    let scales = Scales::default();
    let mut cfg = BaselineConfig::dynamic_accidented_defaults(&scales);
    // Same seed (42) for both runs — `dynamic_accidented_defaults`
    // already pins it. Setting explicitly for documentation.
    cfg.seed = 42;
    cfg.grid_nx = 32;
    cfg.grid_ny = 32;
    cfg.steps = 10;
    cfg.total_time_nondim = 0.6; // dt = 0.06 (matches R5.0 baseline dt)
    cfg.heightmap_fractions = vec![];
    cfg.output_dir = std::path::PathBuf::from(format!(
        "target/step12_r5b_d1/{label}"
    ));
    // Mantle ON at default mf=1.0 — same regime that saturated
    cfg.mantle = MantleConfig::Enabled {
        mf: 1.0,
        coupling: 1.0,
        num_modes: 6,
        seed: 7,
        evolution_rate: 0.0,
    };
    // The linear solver under test
    cfg.linear_solver = linear_solver;

    println!(
        "\n=== R5b D1: {label} — 10 steps × 32² × mf=1.0, workflow OFF ==="
    );
    let t0 = std::time::Instant::now();
    let result = run_baseline(&cfg);
    let elapsed = t0.elapsed().as_secs_f64();
    let m = &result.metrics;

    // Per-step distributions (from NewtonAggregate when present).
    let (newton_per_step, cg_per_newton_step, converged, stalled, diverged, capped) =
        if let Some(newt) = m.newton.as_ref() {
            (
                newt.outer_iters.clone(),
                newt.cg_iters_per_newton_step.clone(),
                newt.converged,
                newt.stalled,
                newt.diverged,
                newt.capped,
            )
        } else {
            (vec![], vec![], 0, 0, 0, 0)
        };

    println!(
        "[r5b/{label}] {:.1}s ({:.2}s/step), CG mean={:.1} max={}, kappa={:.2e}, peak|v|={:.3e}, outcomes (C/S/D/Cap)={}/{}/{}/{}",
        elapsed, elapsed / 10.0, m.cg_iter_mean, m.cg_iter_max, m.kappa_estimate, m.vmax_peak,
        converged, stalled, diverged, capped
    );

    let cap_2000 = m.cg_iter_max >= 2000;
    let mean_high = m.cg_iter_mean > 1500.0;
    let mean_mid = m.cg_iter_mean > 500.0 && m.cg_iter_mean <= 1500.0;
    let mean_low = m.cg_iter_mean <= 500.0;
    let regime = match (cap_2000, mean_high, mean_mid, mean_low) {
        (_, true, _, _) => "**SATURATED** (CG mean > 1500)",
        (true, false, true, _) => "**SATURATED partial** (CG mean 500-1500, cap still hit)",
        (false, false, true, _) => "**TRANSITIONAL** (CG mean 500-1500, no cap)",
        (_, _, _, true) => "**CONVERGENT** (CG mean ≤ 500)",
        _ => "intermediate (read raw metrics)",
    };

    let mut md = String::new();
    md.push_str(&format!("# R5b D1 — {label}\n\n"));
    md.push_str("32² × 10 steps × mf=1.0 × workflow OFF, seed=42.\n\n");
    md.push_str("## Solver health (aggregate)\n\n");
    md.push_str(&format!("- Linear solver: `{label}`\n"));
    md.push_str(&format!(
        "- Total runtime: **{:.1}s** ({:.2}s/step)\n",
        elapsed,
        elapsed / 10.0
    ));
    md.push_str(&format!(
        "- CG iter mean: **{:.1}**  (max: {})\n",
        m.cg_iter_mean, m.cg_iter_max
    ));
    md.push_str(&format!("- Kappa estimate: {:.2e}\n", m.kappa_estimate));
    md.push_str(&format!("- Peak |v|: {:.3e}\n", m.vmax_peak));
    md.push_str(&format!(
        "- Outcomes (Converged / Stalled / Diverged / Capped): **{} / {} / {} / {}**\n\n",
        converged, stalled, diverged, capped
    ));

    md.push_str("## Per-step Newton outer iterations\n\n");
    if newton_per_step.is_empty() {
        md.push_str("_NewtonAggregate not populated (legacy Picard path)._\n\n");
    } else {
        md.push_str("| step | Newton outer iter |\n|---|---|\n");
        for (i, n) in newton_per_step.iter().enumerate() {
            md.push_str(&format!("| {} | {} |\n", i + 1, n));
        }
        md.push('\n');
    }

    md.push_str("## CG iterations per Newton inner solve\n\n");
    if cg_per_newton_step.is_empty() {
        md.push_str("_no inner-CG samples recorded_\n\n");
    } else {
        let total = cg_per_newton_step.len();
        let mean = cg_per_newton_step.iter().sum::<usize>() as f64 / total as f64;
        let min = *cg_per_newton_step.iter().min().unwrap_or(&0);
        let max = *cg_per_newton_step.iter().max().unwrap_or(&0);
        md.push_str(&format!(
            "- Total Newton inner solves: **{}**\n- CG iter min/mean/max: **{} / {:.1} / {}**\n\n",
            total, min, mean, max
        ));
        // First 20 samples explicit
        let take = total.min(20);
        md.push_str(&format!(
            "First {take} inner-CG iter counts: `{:?}`\n\n",
            &cg_per_newton_step[..take]
        ));
    }

    md.push_str("## Histogram (CG iter per inner solve)\n\n");
    md.push_str(&format!(
        "- bin edges (≤): `{:?}`\n- counts: `{:?}`\n\n",
        m.cg_iter_histogram.bin_edges, m.cg_iter_histogram.counts
    ));

    md.push_str("## Cost breakdown (AMG-specific, not instrumented)\n\n");
    md.push_str(
        "AMG setup cost (hierarchy build) and per-V-cycle time are **not** \
         exposed by `Metrics` in this iteration. Inferring from total \
         wallclock − Newton iter count × estimated CG iter time is approximate. \
         If decisive, instrument in D2.\n\n",
    );

    md.push_str(&format!("## Regime\n\n{regime}\n"));

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r5b_d1_jacobi_cg() {
    use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
    run_r5b_d1_solver("jacobi_cg", LinearSolverConfig::JacobiCG);
}

#[test]
#[ignore]
fn r5b_d1_amg_cg() {
    use ymir_core::tectonics_v2::stokes::amg::AmgConfig;
    use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
    run_r5b_d1_solver("amg_cg", LinearSolverConfig::AmgCG(AmgConfig::default()));
}

// ── R5b D1-bis — per-step instrumentation, workflow ON, 2 cycles ────
//
// D1 found that the solver converges fine in isolation (workflow OFF,
// 10 steps × mf=1.0 → CG mean = 342, all 15 Newton outcomes Converged).
// R5.0 / R5.0.1 measured saturation only on workflow ON. The bug is
// specific to the macro_redistribution → tectonic continuation
// interaction.
//
// D1-bis: instrument run_phase_a_loop_v2 on a short 2-cycle × 20-step
// configuration. Extract per-step CG iter and Newton outer iter,
// labelled by `(cycle_idx, step_idx_in_cycle)`. Goal: localise
// where in the cycle the CG saturates.
//
// Three possible patterns (per user decision tree):
//   (a) Transient — saturation concentrated on steps 1-3 post-macro,
//       decays within the cycle. Smoothing post-macro (D1-ter) should
//       fix.
//   (b) Persistent — saturation high throughout cycle. Smoothing
//       insufficient, deeper investigation needed.
//   (c) Cumulative — cycle 2 worse than cycle 1 across the board.
//       Continuation state itself is corrupted by macro between cycles.
//
// Output: `docs/reports/step12_r5b_d1_bis_per_step/`.

fn run_r5b_d1_bis() {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_d1_bis_per_step");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!(
        "\n=== R5b D1-bis: per-step CG/Newton on workflow ON, mf=1.0, 2 cycles × 20 steps ==="
    );

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams {
            n_cycles: 2,
            k_cycle: 20,
            ..V2PhaseAParams::default()
        },
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = 32;
    spec.grid_ny = 32;
    spec.steps = 40; // 2 × 20
    spec.total_time_nondim = 2.4; // dt = 0.06 per step

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r5b/d1-bis] completed in {elapsed:.1}s, {} cycles", output.cycles.len());

    // Build per-step table: for each cycle, walk outer_iters and slice
    // cg_iters_per_newton_step into per-step chunks.
    let mut md = String::new();
    md.push_str("# R5b D1-bis — per-step solver health, workflow ON\n\n");
    md.push_str(&format!(
        "32² active_medley (mf=1.0), workflow ON, 2 cycles × 20 steps, JacobiCG. Runtime: {:.1}s.\n\n",
        elapsed
    ));
    md.push_str("## Per-step CG + Newton iterations\n\n");
    md.push_str("| cycle | step in cycle | newton outer | inner solves | CG sum | CG mean | CG max | hit cap (2000) |\n");
    md.push_str("|---|---|---|---|---|---|---|---|\n");

    for (cycle_idx, cycle) in output.cycles.iter().enumerate() {
        let m = &cycle.baseline.metrics;
        let newt = match m.newton.as_ref() {
            Some(n) => n,
            None => continue,
        };
        // Walk outer_iters, slice cg_iters_per_newton_step accordingly.
        let mut cg_cursor: usize = 0;
        for (step_in_cycle, &outer) in newt.outer_iters.iter().enumerate() {
            let outer_n = outer as usize;
            let cg_slice = &newt.cg_iters_per_newton_step
                [cg_cursor..(cg_cursor + outer_n).min(newt.cg_iters_per_newton_step.len())];
            let cg_sum: usize = cg_slice.iter().sum();
            let cg_mean = if outer_n > 0 { cg_sum as f64 / outer_n as f64 } else { 0.0 };
            let cg_max = *cg_slice.iter().max().unwrap_or(&0);
            let hit_cap = cg_slice.iter().filter(|&&c| c >= 2000).count();
            md.push_str(&format!(
                "| {} | {} | {} | {} | {} | {:.0} | {} | {} |\n",
                cycle_idx + 1,
                step_in_cycle,
                outer,
                outer_n,
                cg_sum,
                cg_mean,
                cg_max,
                hit_cap,
            ));
            cg_cursor += outer_n;
        }
        md.push_str(&format!(
            "\n_Cycle {} aggregate: Converged={}, Stalled={}, Diverged={}, Capped={}, mass_drift={:.3e}_\n\n",
            cycle_idx + 1,
            newt.converged,
            newt.stalled,
            newt.diverged,
            newt.capped,
            cycle.common.mass_drift,
        ));
    }

    // Pattern detection heuristic: compare CG mean for step 0 (post-macro)
    // vs steps 5-10 vs steps 15-19 (end-of-cycle) for each cycle.
    md.push_str("## Pattern detection\n\n");
    for (cycle_idx, cycle) in output.cycles.iter().enumerate() {
        let m = &cycle.baseline.metrics;
        let newt = match m.newton.as_ref() {
            Some(n) => n,
            None => continue,
        };
        let mut cg_per_step: Vec<usize> = Vec::with_capacity(newt.outer_iters.len());
        let mut cg_cursor: usize = 0;
        for &outer in &newt.outer_iters {
            let outer_n = outer as usize;
            let cg_slice = &newt.cg_iters_per_newton_step
                [cg_cursor..(cg_cursor + outer_n).min(newt.cg_iters_per_newton_step.len())];
            cg_per_step.push(cg_slice.iter().sum());
            cg_cursor += outer_n;
        }
        // Quartile-ish aggregation
        let n = cg_per_step.len();
        if n < 4 {
            continue;
        }
        let early = &cg_per_step[..n / 4];
        let mid = &cg_per_step[n / 4..3 * n / 4];
        let late = &cg_per_step[3 * n / 4..];
        let avg = |v: &[usize]| -> f64 { v.iter().sum::<usize>() as f64 / v.len() as f64 };
        md.push_str(&format!(
            "- Cycle {} CG sum/step — early-{} (steps 0..{}): mean={:.0}; mid (steps {}..{}): mean={:.0}; late (steps {}..{}): mean={:.0}\n",
            cycle_idx + 1,
            early.len(),
            n / 4,
            avg(early),
            n / 4,
            3 * n / 4,
            avg(mid),
            3 * n / 4,
            n,
            avg(late),
        ));
    }

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r5b_d1_bis_per_step_workflow_on() {
    run_r5b_d1_bis();
}

/// Post-D2 comparison run: identical config to `r5b_d1_bis_per_step_workflow_on`
/// but exercised against the D2-patched solver. The output dir is
/// suffixed `_post_d2` so the pre-D2 baseline (committed in `c6d21a3`)
/// stays available side-by-side for the comparison table.
#[test]
#[ignore]
fn r5b_d1_bis_per_step_workflow_on_post_d2() {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_d1_bis_per_step_post_d2");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!(
        "\n=== R5b D1-bis POST-D2: per-step CG/Newton, workflow ON, 2 cycles ==="
    );

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams {
            n_cycles: 2,
            k_cycle: 20,
            ..V2PhaseAParams::default()
        },
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = 32;
    spec.grid_ny = 32;
    spec.steps = 40;
    spec.total_time_nondim = 2.4;

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r5b/d1-bis post-D2] completed in {elapsed:.1}s");

    let mut md = String::new();
    md.push_str("# R5b D1-bis POST-D2 — per-step solver health, workflow ON\n\n");
    md.push_str(&format!(
        "32² active_medley (mf=1.0), workflow ON, 2 cycles × 20 steps, JacobiCG + D2 criteria. Runtime: {:.1}s.\n\n",
        elapsed
    ));
    md.push_str("## Per-step CG + Newton iterations\n\n");
    md.push_str("| cycle | step in cycle | newton outer | inner solves | CG sum | CG mean | CG max | hit cap (2000) |\n");
    md.push_str("|---|---|---|---|---|---|---|---|\n");

    for (cycle_idx, cycle) in output.cycles.iter().enumerate() {
        let m = &cycle.baseline.metrics;
        let newt = match m.newton.as_ref() {
            Some(n) => n,
            None => continue,
        };
        let mut cg_cursor: usize = 0;
        for (step_in_cycle, &outer) in newt.outer_iters.iter().enumerate() {
            let outer_n = outer as usize;
            let cg_slice = &newt.cg_iters_per_newton_step
                [cg_cursor..(cg_cursor + outer_n).min(newt.cg_iters_per_newton_step.len())];
            let cg_sum: usize = cg_slice.iter().sum();
            let cg_mean = if outer_n > 0 { cg_sum as f64 / outer_n as f64 } else { 0.0 };
            let cg_max = *cg_slice.iter().max().unwrap_or(&0);
            let hit_cap = cg_slice.iter().filter(|&&c| c >= 2000).count();
            md.push_str(&format!(
                "| {} | {} | {} | {} | {} | {:.0} | {} | {} |\n",
                cycle_idx + 1,
                step_in_cycle,
                outer,
                outer_n,
                cg_sum,
                cg_mean,
                cg_max,
                hit_cap,
            ));
            cg_cursor += outer_n;
        }
        md.push_str(&format!(
            "\n_Cycle {} aggregate: Converged={}, Stalled={}, Diverged={}, Capped={}, mass_drift={:.3e}_\n\n",
            cycle_idx + 1,
            newt.converged,
            newt.stalled,
            newt.diverged,
            newt.capped,
            cycle.common.mass_drift,
        ));
    }

    // Same quartile aggregation
    md.push_str("## Pattern detection\n\n");
    for (cycle_idx, cycle) in output.cycles.iter().enumerate() {
        let m = &cycle.baseline.metrics;
        let newt = match m.newton.as_ref() {
            Some(n) => n,
            None => continue,
        };
        let mut cg_per_step: Vec<usize> = Vec::with_capacity(newt.outer_iters.len());
        let mut cg_cursor: usize = 0;
        for &outer in &newt.outer_iters {
            let outer_n = outer as usize;
            let cg_slice = &newt.cg_iters_per_newton_step
                [cg_cursor..(cg_cursor + outer_n).min(newt.cg_iters_per_newton_step.len())];
            cg_per_step.push(cg_slice.iter().sum());
            cg_cursor += outer_n;
        }
        let n = cg_per_step.len();
        if n < 4 { continue; }
        let early = &cg_per_step[..n / 4];
        let mid = &cg_per_step[n / 4..3 * n / 4];
        let late = &cg_per_step[3 * n / 4..];
        let avg = |v: &[usize]| -> f64 { v.iter().sum::<usize>() as f64 / v.len() as f64 };
        md.push_str(&format!(
            "- Cycle {} CG sum/step — early-{}: mean={:.0}; mid: mean={:.0}; late: mean={:.0}\n",
            cycle_idx + 1, early.len(), avg(early), avg(mid), avg(late),
        ));
    }

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

// ── R5b D2-bis — Oscillating physical evolution analysis ────────────
//
// D2 post-fix surfaced 26 Oscillating outcomes on cycle 2 of the
// active_medley 32² × 2 cycles × mf=1.0 workflow ON regime. D2-bis
// instruments step-to-step physical state (peek_s, peek_vx, peek_vy
// via the harness on_progress callback) to classify the oscillation:
//
//   A. Bénin — `‖Δv‖/‖v‖` small throughout, smooth evolution
//   B. Problème non-unique — `‖Δv‖/‖v‖` > 10 % on Oscillating
//   C. Init transient — first Oscillating step `‖Δv‖/‖v‖` >> mean
//
// Output: `docs/reports/step12_r5b_d2_bis_oscillating/report.md`.

#[derive(Clone, Debug)]
struct StepSnapshot {
    step_global: usize,
    cycle: usize,
    step_in_cycle: usize,
    peak_v: f64,
    s_mean: f64,
    s_max: f64,
    mass_total: f64,
    rel_dv: f64,
    rel_ds: f64,
}

#[test]
#[ignore]
fn r5b_d2_bis_oscillating_evolution() {
    use std::cell::RefCell;
    use ymir_core::tectonics_v2::workflow::{
        final_state_to_continuation_v2, run_phase_a_cycle_with_progress_v2,
    };
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_d2_bis_oscillating");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!("\n=== R5b D2-bis: per-step Oscillating physical evolution ===");

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams {
            n_cycles: 2,
            k_cycle: 20,
            ..V2PhaseAParams::default()
        },
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = 32;
    spec.grid_ny = 32;
    spec.steps = 40;
    spec.total_time_nondim = 2.4;

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let captured: RefCell<Vec<StepSnapshot>> = RefCell::new(Vec::new());
    let prev_state: RefCell<Option<(Vec<f64>, Vec<f64>, Vec<f64>)>> = RefCell::new(None);
    let cycle_counter: RefCell<usize> = RefCell::new(0);
    let step_in_cycle_counter: RefCell<usize> = RefCell::new(0);
    let step_global_counter: RefCell<usize> = RefCell::new(0);

    let mut record_step = |progress: &ymir_core::tectonics_v2::diagnostics::harness::StepProgress<'_>| -> bool {
        let peek_s = progress.peek_s.data().to_vec();
        let peek_vx = progress.peek_vx.to_vec();
        let peek_vy = progress.peek_vy.to_vec();
        let n = peek_s.len();
        let peak_v = (0..n)
            .map(|i| (peek_vx[i].powi(2) + peek_vy[i].powi(2)).sqrt())
            .fold(0.0_f64, f64::max);
        let s_mean = peek_s.iter().sum::<f64>() / n as f64;
        let s_max = peek_s.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mass_total: f64 = peek_s.iter().sum();

        let (rel_dv, rel_ds) = {
            let prev = prev_state.borrow();
            if let Some((pvx, pvy, ps)) = prev.as_ref() {
                let v_norm_sq: f64 = peek_vx.iter().map(|x| x * x).sum::<f64>()
                    + peek_vy.iter().map(|x| x * x).sum::<f64>();
                let dv_norm_sq: f64 = peek_vx
                    .iter()
                    .zip(pvx)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    + peek_vy
                        .iter()
                        .zip(pvy)
                        .map(|(a, b)| (a - b).powi(2))
                        .sum::<f64>();
                let s_norm_sq: f64 = peek_s.iter().map(|x| x * x).sum::<f64>();
                let ds_norm_sq: f64 = peek_s
                    .iter()
                    .zip(ps)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>();
                let rel_dv = (dv_norm_sq / v_norm_sq.max(1e-30)).sqrt();
                let rel_ds = (ds_norm_sq / s_norm_sq.max(1e-30)).sqrt();
                (rel_dv, rel_ds)
            } else {
                (0.0, 0.0)
            }
        };

        captured.borrow_mut().push(StepSnapshot {
            step_global: *step_global_counter.borrow(),
            cycle: *cycle_counter.borrow(),
            step_in_cycle: *step_in_cycle_counter.borrow(),
            peak_v,
            s_mean,
            s_max,
            mass_total,
            rel_dv,
            rel_ds,
        });

        *prev_state.borrow_mut() = Some((peek_vx, peek_vy, peek_s));
        *step_in_cycle_counter.borrow_mut() += 1;
        *step_global_counter.borrow_mut() += 1;
        true
    };

    let t0 = std::time::Instant::now();

    // Cycle 1
    *cycle_counter.borrow_mut() = 1;
    *step_in_cycle_counter.borrow_mut() = 0;
    let cycle_1 = run_phase_a_cycle_with_progress_v2(&cfg, &wf, &mut record_step);
    cfg.continuation =
        Some(final_state_to_continuation_v2(&cycle_1.baseline.final_state));

    // Cycle 2
    *cycle_counter.borrow_mut() = 2;
    *step_in_cycle_counter.borrow_mut() = 0;
    let cycle_2 = run_phase_a_cycle_with_progress_v2(&cfg, &wf, &mut record_step);

    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r5b/d2-bis] completed in {elapsed:.1}s");

    let newt_1 = cycle_1.baseline.metrics.newton.as_ref().expect("newton 1");
    let newt_2 = cycle_2.baseline.metrics.newton.as_ref().expect("newton 2");

    let mut md = String::new();
    md.push_str("# R5b D2-bis — per-step physical evolution + Newton outcomes\n\n");
    md.push_str(&format!(
        "Runtime: {:.1}s. 2 cycles, JacobiCG + D2 criteria, mf=1.0 active_medley 32².\n\n",
        elapsed
    ));
    md.push_str("## Cycle aggregates\n\n");
    md.push_str(&format!(
        "- Cycle 1: Converged={}, Stalled={} (Oscillating ⊂ Stalled), Diverged={}, Capped={}\n",
        newt_1.converged, newt_1.stalled, newt_1.diverged, newt_1.capped
    ));
    md.push_str(&format!(
        "- Cycle 2: Converged={}, Stalled={} (Oscillating ⊂ Stalled), Diverged={}, Capped={}\n\n",
        newt_2.converged, newt_2.stalled, newt_2.diverged, newt_2.capped
    ));

    md.push_str("## Per-step physical state evolution\n\n");
    md.push_str("| cycle | step | newton outer | ‖Δv‖/‖v‖ | ‖ΔS̃‖/‖S̃‖ | peak\\|v\\| | S̃ mean | S̃ max | mass |\n");
    md.push_str("|---|---|---|---|---|---|---|---|---|\n");

    let snaps = captured.borrow();
    let outer_iters_cycle = |cy: usize| -> &[u32] {
        if cy == 1 { &newt_1.outer_iters } else { &newt_2.outer_iters }
    };
    for snap in snaps.iter() {
        let outer = outer_iters_cycle(snap.cycle)
            .get(snap.step_in_cycle)
            .copied()
            .unwrap_or(0);
        md.push_str(&format!(
            "| {} | {} | {} | {:.3e} | {:.3e} | {:.3e} | {:.3} | {:.3} | {:.3} |\n",
            snap.cycle, snap.step_in_cycle, outer, snap.rel_dv, snap.rel_ds,
            snap.peak_v, snap.s_mean, snap.s_max, snap.mass_total,
        ));
    }
    md.push('\n');

    md.push_str("## Sub-case classification (cycle 2)\n\n");
    let cycle2_snaps: Vec<&StepSnapshot> = snaps.iter().filter(|s| s.cycle == 2).collect();
    if !cycle2_snaps.is_empty() {
        let dv_first = cycle2_snaps.first().map(|s| s.rel_dv).unwrap_or(0.0);
        let dv_max = cycle2_snaps.iter().map(|s| s.rel_dv).fold(0.0_f64, f64::max);
        let dv_mean = cycle2_snaps.iter().map(|s| s.rel_dv).sum::<f64>()
            / cycle2_snaps.len() as f64;
        let mass_init = cycle2_snaps.first().map(|s| s.mass_total).unwrap_or(0.0);
        let mass_final = cycle2_snaps.last().map(|s| s.mass_total).unwrap_or(0.0);
        let mass_change_pct = (mass_final - mass_init) / mass_init.max(1e-12) * 100.0;
        md.push_str(&format!(
            "- ‖Δv‖/‖v‖ : first step = {:.3e}, max = {:.3e}, mean = {:.3e}\n",
            dv_first, dv_max, dv_mean
        ));
        md.push_str(&format!(
            "- Mass evolution : init = {:.3}, final = {:.3} (Δ {:.3} %)\n",
            mass_init, mass_final, mass_change_pct
        ));

        let verdict = if dv_max < 1e-2 {
            "**Sous-cas A (bénin)** — `‖Δv‖/‖v‖ max < 1 %` ; oscillation autour d'un attracteur stable, état committé acceptable."
        } else if dv_first > 5.0 * dv_mean && dv_max > 0.05 {
            "**Sous-cas C (init transient)** — premier step Oscillating montre `‖Δv‖/‖v‖` >> moyenne ; init continuation mauvais. D1-ter peut aider mais cible la cause."
        } else if dv_max > 0.1 {
            "**Sous-cas B (problème non-unique)** — `‖Δv‖/‖v‖` > 10 % sur Oscillating ; état committé arbitraire ; D1-ter smoothing critique."
        } else {
            "**Intermédiaire** — pattern mixte, examen visuel nécessaire."
        };
        md.push_str(&format!("\n{verdict}\n"));
    }

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

// ── R5b D1-ter — init variant diagnostic for cycle 2 oscillation ────
//
// D2-bis classified cycle 2 oscillation as sub-case C (init transient,
// `‖Δv‖/‖v‖ max = 1143 %`). Before committing to a specific fix
// (reinit v=0, smoothing post-macro, or other), this diagnostic runs
// three init variants on the SAME post-cycle-1 state and measures
// Newton iter / CG iter / outcome on cycle 2 steps 0-4 (where the
// oscillation peaks).
//
// Variants:
//   A. warm_start : v = v_final_cycle_1 (current behaviour, baseline)
//   B. zero       : v = 0
//   C. smoothed   : v = 3x3 Gaussian-smoothed v_final_cycle_1 on
//                   periodic grid
//
// Output: `docs/reports/step12_r5b_d1_ter_init_variants/<label>/`.

fn apply_3x3_gauss_periodic(data: &mut [f64], nx: usize, ny: usize) {
    let src: Vec<f64> = data.to_vec();
    let k: [[f64; 3]; 3] = [
        [1.0, 2.0, 1.0],
        [2.0, 4.0, 2.0],
        [1.0, 2.0, 1.0],
    ];
    let norm = 16.0;
    for j in 0..ny {
        for i in 0..nx {
            let mut sum = 0.0;
            for dj in 0..3 {
                let jj = (j + ny + dj - 1) % ny;
                for di in 0..3 {
                    let ii = (i + nx + di - 1) % nx;
                    sum += k[dj][di] * src[jj * nx + ii];
                }
            }
            data[j * nx + i] = sum / norm;
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum InitVariant {
    WarmStart,
    Zero,
    Smoothed,
}

fn run_d1_ter_init_variant(variant: InitVariant, label: &str) {
    use std::cell::RefCell;
    use ymir_core::tectonics_v2::workflow::{
        final_state_to_continuation_v2, run_phase_a_cycle_with_progress_v2,
    };
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_d1_ter_init_variants")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!("\n=== R5b D1-ter init variant: {label} ===");

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams {
            n_cycles: 2,
            k_cycle: 20,
            ..V2PhaseAParams::default()
        },
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = 32;
    spec.grid_ny = 32;
    spec.steps = 40;
    spec.total_time_nondim = 2.4;

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    // Cycle 1 — identical for all variants (baseline)
    let no_op =
        |_: &ymir_core::tectonics_v2::diagnostics::harness::StepProgress<'_>| -> bool { true };
    let t0 = std::time::Instant::now();
    let cycle_1 = run_phase_a_cycle_with_progress_v2(&cfg, &wf, no_op);
    let cycle_1_elapsed = t0.elapsed().as_secs_f64();
    println!("[d1-ter/{label}] cycle 1 done in {cycle_1_elapsed:.1}s");

    let nx = cycle_1.baseline.final_state.s_field.nx();
    let ny = cycle_1.baseline.final_state.s_field.ny();
    let mut continuation =
        final_state_to_continuation_v2(&cycle_1.baseline.final_state);
    match variant {
        InitVariant::WarmStart => {}
        InitVariant::Zero => {
            for v in continuation.vx.iter_mut() { *v = 0.0; }
            for v in continuation.vy.iter_mut() { *v = 0.0; }
        }
        InitVariant::Smoothed => {
            apply_3x3_gauss_periodic(&mut continuation.vx, nx, ny);
            apply_3x3_gauss_periodic(&mut continuation.vy, nx, ny);
        }
    }
    cfg.continuation = Some(continuation);

    let captured: RefCell<Vec<StepSnapshot>> = RefCell::new(Vec::new());
    let prev_state: RefCell<Option<(Vec<f64>, Vec<f64>, Vec<f64>)>> =
        RefCell::new(None);
    let step_counter: RefCell<usize> = RefCell::new(0);

    let mut record_step =
        |progress: &ymir_core::tectonics_v2::diagnostics::harness::StepProgress<'_>| -> bool {
        let peek_s = progress.peek_s.data().to_vec();
        let peek_vx = progress.peek_vx.to_vec();
        let peek_vy = progress.peek_vy.to_vec();
        let n = peek_s.len();
        let peak_v = (0..n)
            .map(|i| (peek_vx[i].powi(2) + peek_vy[i].powi(2)).sqrt())
            .fold(0.0_f64, f64::max);
        let s_mean = peek_s.iter().sum::<f64>() / n as f64;
        let s_max = peek_s.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mass_total: f64 = peek_s.iter().sum();
        let (rel_dv, rel_ds) = {
            let prev = prev_state.borrow();
            if let Some((pvx, pvy, ps)) = prev.as_ref() {
                let v_norm_sq: f64 = peek_vx.iter().map(|x| x * x).sum::<f64>()
                    + peek_vy.iter().map(|x| x * x).sum::<f64>();
                let dv_norm_sq: f64 = peek_vx.iter().zip(pvx).map(|(a, b)| (a - b).powi(2)).sum::<f64>()
                    + peek_vy.iter().zip(pvy).map(|(a, b)| (a - b).powi(2)).sum::<f64>();
                let s_norm_sq: f64 = peek_s.iter().map(|x| x * x).sum::<f64>();
                let ds_norm_sq: f64 = peek_s.iter().zip(ps).map(|(a, b)| (a - b).powi(2)).sum::<f64>();
                (
                    (dv_norm_sq / v_norm_sq.max(1e-30)).sqrt(),
                    (ds_norm_sq / s_norm_sq.max(1e-30)).sqrt(),
                )
            } else { (0.0, 0.0) }
        };
        captured.borrow_mut().push(StepSnapshot {
            step_global: *step_counter.borrow(),
            cycle: 2,
            step_in_cycle: *step_counter.borrow(),
            peak_v, s_mean, s_max, mass_total, rel_dv, rel_ds,
        });
        *prev_state.borrow_mut() = Some((peek_vx, peek_vy, peek_s));
        *step_counter.borrow_mut() += 1;
        true
    };

    let t1 = std::time::Instant::now();
    let cycle_2 = run_phase_a_cycle_with_progress_v2(&cfg, &wf, &mut record_step);
    let cycle_2_elapsed = t1.elapsed().as_secs_f64();
    println!("[d1-ter/{label}] cycle 2 done in {cycle_2_elapsed:.1}s");

    let newt_2 = cycle_2.baseline.metrics.newton.as_ref().expect("newton 2");

    let mut md = String::new();
    md.push_str(&format!("# R5b D1-ter init variant — {label}\n\n"));
    md.push_str(&format!(
        "Cycle 1 runtime: {:.1}s. Cycle 2 runtime: {:.1}s.\n\n",
        cycle_1_elapsed, cycle_2_elapsed
    ));
    md.push_str(&format!(
        "Cycle 2 outcomes: Converged={}, Stalled={}, Diverged={}, Capped={}\n\n",
        newt_2.converged, newt_2.stalled, newt_2.diverged, newt_2.capped
    ));

    md.push_str("## Cycle 2 first 10 steps\n\n");
    md.push_str("| step | newton outer | ‖Δv‖/‖v‖ | ‖ΔS̃‖/‖S̃‖ | peak\\|v\\| | S̃ max | mass |\n");
    md.push_str("|---|---|---|---|---|---|---|\n");
    let snaps = captured.borrow();
    let take = snaps.len().min(10);
    for (idx, snap) in snaps.iter().take(take).enumerate() {
        let outer = newt_2.outer_iters.get(idx).copied().unwrap_or(0);
        md.push_str(&format!(
            "| {} | {} | {:.3e} | {:.3e} | {:.3e} | {:.3} | {:.3} |\n",
            idx, outer, snap.rel_dv, snap.rel_ds, snap.peak_v, snap.s_max, snap.mass_total,
        ));
    }
    md.push('\n');

    let first_5: Vec<&StepSnapshot> = snaps.iter().take(5).collect();
    if !first_5.is_empty() {
        let dv_max = first_5.iter().map(|s| s.rel_dv).fold(0.0_f64, f64::max);
        let dv_mean = first_5.iter().map(|s| s.rel_dv).sum::<f64>() / first_5.len() as f64;
        let newton_total: u32 = newt_2.outer_iters.iter().take(5).sum();
        let cg_total: usize = newt_2.cg_iters_per_newton_step.iter().take(newton_total as usize).sum();
        md.push_str("## Aggregate first 5 steps cycle 2\n\n");
        md.push_str(&format!("- ‖Δv‖/‖v‖ max : {:.3e}\n", dv_max));
        md.push_str(&format!("- ‖Δv‖/‖v‖ mean : {:.3e}\n", dv_mean));
        md.push_str(&format!("- Newton outer iter total : {}\n", newton_total));
        md.push_str(&format!("- CG iter total : {}\n", cg_total));
    }

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

#[test]
#[ignore]
fn r5b_d1_ter_init_a_warm_start() {
    run_d1_ter_init_variant(InitVariant::WarmStart, "a_warm_start");
}

#[test]
#[ignore]
fn r5b_d1_ter_init_b_zero() {
    run_d1_ter_init_variant(InitVariant::Zero, "b_zero");
}

#[test]
#[ignore]
fn r5b_d1_ter_init_c_smoothed() {
    run_d1_ter_init_variant(InitVariant::Smoothed, "c_smoothed");
}

/// Step 12 R5b — re-measure D1-bis pattern AFTER the D1-ter implement
/// (reinit v=0 post-macro inside `phase_a::run_phase_a_cycle_with_progress_v2`).
/// Output dir is distinct from the pre-D1-ter version so the
/// comparison report stays available.
/// Step 12 R5b validation — R4 active_medley 64² × 5 cycles after the
/// full D2 + D1-ter stack. This is the product-validation gate: if
/// continents are visible at cycle 5 with no Oscillating outcomes,
/// Step 12 is mergeable.
#[test]
#[ignore]
fn r5b_validation_r4_active_medley_64sq() {
    let out_dir = out_root().join("active_medley_post_d1_ter");
    let metrics = run_gallery_for_preset("active_medley", &out_dir);
    write_metrics_table("active_medley_post_d1_ter", &metrics, &out_dir);
}

/// Step 12 R5b validation — R4 active_medley 64² × 5 cycles with
/// mantle override mf=0.5 (R4b.5 finding) + the full D2 + D1-ter
/// solver stack. Decides whether the mf=0.5 calibration produces a
/// viable Living Landz product (continents preserved AND tectonic
/// dynamics visible: peak |v| non-trivial, mountain ranges at plate
/// boundaries, evolving coastlines). Output dir is suffixed
/// `_mf_0_5_post_d1_ter`.
#[test]
#[ignore]
fn r5b_validation_r4_active_medley_64sq_mf_0_5() {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = out_root().join("active_medley_mf_0_5_post_d1_ter");
    std::fs::create_dir_all(&out_dir).expect("create output dir");
    println!("\n=== Preset: active_medley (mf=0.5 override) ===");

    let mut spec = presets::load("active_medley")
        .unwrap_or_else(|e| panic!("preset 'active_medley' load failed: {e}"));
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    // R4b.5 sweet-spot override : mf 1.0 → 0.5 (preserving other
    // mantle fields).
    if let V2MantleSpec::On {
        coupling,
        num_modes,
        seed,
        evolution_rate,
        ..
    } = spec.mantle
    {
        spec.mantle = V2MantleSpec::On {
            mf: 0.5,
            coupling,
            num_modes,
            seed,
            evolution_rate,
        };
    }
    spec.grid_nx = 64;
    spec.grid_ny = 64;
    spec.steps = N_CYCLES * K_CYCLE;
    spec.total_time_nondim = 6.0;

    // INIT capture (mf=0.5, default mantle is already On in active_medley)
    let init_state = make_init_state(&spec);
    save_field_png(
        &init_state,
        V2Field::SThickness,
        &out_dir.join("cycle_0_s.png"),
    )
    .expect("save init s");
    save_field_png(
        &init_state,
        V2Field::Altitude,
        &out_dir.join("cycle_0_altitude.png"),
    )
    .expect("save init alt");

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    println!("[r5b validation] running 5 cycles × 20 steps (64², mf=0.5, D2+D1-ter active)…");
    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "[r5b validation] completed in {elapsed:.1}s, {} cycles",
        output.cycles.len()
    );

    let init_field =
        Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output
        .cycles
        .first()
        .map(|c| c.common.sea_level_normalized)
        .unwrap_or(0.5);
    let mut metrics = vec![compute_metrics(0, &init_field, init_sea, 0.0, 0.0)];

    for (i, cycle) in output.cycles.iter().enumerate() {
        let cyc_idx = i + 1;
        let v2_state = V2FinalState::from_harness(&cycle.baseline.final_state);
        save_field_png(
            &v2_state,
            V2Field::SThickness,
            &out_dir.join(format!("cycle_{cyc_idx}_s.png")),
        )
        .expect("save s");
        save_field_png(
            &v2_state,
            V2Field::Altitude,
            &out_dir.join(format!("cycle_{cyc_idx}_altitude.png")),
        )
        .expect("save altitude");

        metrics.push(compute_metrics(
            cyc_idx,
            &cycle.baseline.final_state.s_field,
            cycle.common.sea_level_normalized,
            cycle.common.erosion_volume_removed,
            cycle.common.mass_drift,
        ));
    }

    // Add peak|v| per cycle to the report (D2+D1-ter answered solver
    // convergence; this report needs to answer "is the tectonic
    // dynamic non-trivial?").
    let peak_v_per_cycle: Vec<f64> = output
        .cycles
        .iter()
        .map(|c| c.baseline.metrics.vmax_peak)
        .collect();

    write_metrics_table("active_medley_mf_0_5_post_d1_ter", &metrics, &out_dir);

    // Append peak|v| info to the metrics.md as a dynamics-check row.
    let dyn_extra = format!(
        "\n## Dynamics probe (peak |v| per cycle)\n\n{}\n\nIf all peak|v| > 0.5 (nondim) the run carries non-trivial tectonic dynamics. < 0.1 → quasi-static (bad).\n",
        peak_v_per_cycle
            .iter()
            .enumerate()
            .map(|(i, v)| format!("- Cycle {}: peak |v| = {:.3e}", i + 1, v))
            .collect::<Vec<_>>()
            .join("\n")
    );
    let path = out_dir.join("metrics.md");
    let existing =
        std::fs::read_to_string(&path).unwrap_or_default();
    std::fs::write(&path, format!("{existing}{dyn_extra}"))
        .expect("append dynamics probe");
    println!("Dynamics probe appended to {}", path.display());
}

#[test]
#[ignore]
fn r5b_d1_bis_per_step_workflow_on_post_d1_ter() {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r5b_d1_bis_per_step_post_d1_ter");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!(
        "\n=== R5b D1-bis POST-D1-ter: per-step CG/Newton, workflow ON, 2 cycles ==="
    );

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams {
            n_cycles: 2,
            k_cycle: 20,
            ..V2PhaseAParams::default()
        },
        phase_b: V2PhaseBParams::default(),
    };
    spec.grid_nx = 32;
    spec.grid_ny = 32;
    spec.steps = 40;
    spec.total_time_nondim = 2.4;

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r5b/d1-bis post-D1-ter] completed in {elapsed:.1}s");

    let mut md = String::new();
    md.push_str("# R5b D1-bis POST-D1-ter — per-step solver health\n\n");
    md.push_str(&format!(
        "32² active_medley (mf=1.0), workflow ON, 2 cycles × 20 steps, JacobiCG + D2 criteria + D1-ter reinit v=0. Runtime: {:.1}s.\n\n",
        elapsed
    ));
    md.push_str("## Per-step CG + Newton iterations\n\n");
    md.push_str("| cycle | step in cycle | newton outer | inner solves | CG sum | CG mean | CG max | hit cap (2000) |\n");
    md.push_str("|---|---|---|---|---|---|---|---|\n");

    for (cycle_idx, cycle) in output.cycles.iter().enumerate() {
        let m = &cycle.baseline.metrics;
        let newt = match m.newton.as_ref() {
            Some(n) => n,
            None => continue,
        };
        let mut cg_cursor: usize = 0;
        for (step_in_cycle, &outer) in newt.outer_iters.iter().enumerate() {
            let outer_n = outer as usize;
            let cg_slice = &newt.cg_iters_per_newton_step
                [cg_cursor..(cg_cursor + outer_n).min(newt.cg_iters_per_newton_step.len())];
            let cg_sum: usize = cg_slice.iter().sum();
            let cg_mean = if outer_n > 0 { cg_sum as f64 / outer_n as f64 } else { 0.0 };
            let cg_max = *cg_slice.iter().max().unwrap_or(&0);
            let hit_cap = cg_slice.iter().filter(|&&c| c >= 2000).count();
            md.push_str(&format!(
                "| {} | {} | {} | {} | {} | {:.0} | {} | {} |\n",
                cycle_idx + 1, step_in_cycle, outer, outer_n, cg_sum, cg_mean, cg_max, hit_cap,
            ));
            cg_cursor += outer_n;
        }
        md.push_str(&format!(
            "\n_Cycle {} aggregate: Converged={}, Stalled={}, Diverged={}, Capped={}, mass_drift={:.3e}_\n\n",
            cycle_idx + 1,
            newt.converged, newt.stalled, newt.diverged, newt.capped,
            cycle.common.mass_drift,
        ));
    }

    md.push_str("## Pattern detection\n\n");
    for (cycle_idx, cycle) in output.cycles.iter().enumerate() {
        let m = &cycle.baseline.metrics;
        let newt = match m.newton.as_ref() {
            Some(n) => n,
            None => continue,
        };
        let mut cg_per_step: Vec<usize> = Vec::with_capacity(newt.outer_iters.len());
        let mut cg_cursor: usize = 0;
        for &outer in &newt.outer_iters {
            let outer_n = outer as usize;
            let cg_slice = &newt.cg_iters_per_newton_step
                [cg_cursor..(cg_cursor + outer_n).min(newt.cg_iters_per_newton_step.len())];
            cg_per_step.push(cg_slice.iter().sum());
            cg_cursor += outer_n;
        }
        let n = cg_per_step.len();
        if n < 4 { continue; }
        let early = &cg_per_step[..n / 4];
        let mid = &cg_per_step[n / 4..3 * n / 4];
        let late = &cg_per_step[3 * n / 4..];
        let avg = |v: &[usize]| -> f64 { v.iter().sum::<usize>() as f64 / v.len() as f64 };
        md.push_str(&format!(
            "- Cycle {} CG sum/step — early-{}: mean={:.0}; mid: mean={:.0}; late: mean={:.0}\n",
            cycle_idx + 1, early.len(), avg(early), avg(mid), avg(late),
        ));
    }

    std::fs::write(out_dir.join("report.md"), &md).expect("write report.md");
    println!("\n{md}");
}

// ----------------------------------------------------------------------------
// Step 12 R6.3 — mantle.evolution_rate × mf product validation sweep.
//
// 12 configs: mf ∈ {0.5, 0.7, 0.8, 1.0} × evo ∈ {0.05, 0.10, 0.20}, at 32²,
// workflow ON (D2 + D1-ter active), 5 cycles × 20 steps each. Per-config
// outputs go under `docs/reports/step12_r6_3_sweep/<label>/`, plus a
// top-level `summary.md` with the comparison table and EVO.A/B/C/D verdict.
//
// Acceptance criteria R4.1–R4.6 (6 axes from R6 plan):
//   R4.1 — Continents émergés         (peak S̃_final > sea_level)
//   R4.2 — Cratons préservés          (frac>0.8 retention > 50 % init)
//   R4.3 — Bordures continentales évoluées (VISUAL, manual)
//   R4.4 — Conservation               (|mass loss| < 1 % per cycle)
//   R4.5 — Drainage actif             (max_path ≥ 5)
//   R4.6 — Dynamique soutenue         (peak|v| > 0.1 nondim on ≥ 3 cycles)
//
// A config passes 5/5 auto-criteria → EVO.A candidate (R4.3 visual still
// required). 5 auto-fail → EVO.D. Visual fail on R4.3 with auto-pass → EVO.C.
//
// Reference for evo=0: R5b mf sweep results in
// `docs/reports/step12_r5b_mf_sweep_post_d1_ter/` (mf ∈ {0.5, 0.7, 0.8})
// and R4 active_medley_post_d1_ter for mf=1.0.

#[allow(dead_code)]
#[derive(Clone, Debug)]
struct R6ConfigSummary {
    mf: f64,
    evo_rate: f64,
    label: String,
    runtime_s: f64,
    peak_v_per_cycle: Vec<f64>,
    peak_s_per_cycle: Vec<f64>,
    frac_above_0_8_per_cycle: Vec<f64>,
    mass_total_per_cycle: Vec<f64>,
    max_path_per_cycle: Vec<u8>,
    sea_level_per_cycle: Vec<f64>,
    newton_outcomes_per_cycle: Vec<(usize, usize, usize, usize)>,
    cg_mean_quartiles_per_cycle: Vec<(f64, f64, f64)>,
    init_frac_above_0_8: f64,
    init_peak_s: f64,
}

fn run_r6_3_config(mf: f64, evo_rate: f64, label: &str) -> R6ConfigSummary {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r6_3_sweep")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!("\n=== R6.3 sweep: {label} (mf={mf}, evolution_rate={evo_rate}) ===");

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    if let V2MantleSpec::On { coupling, num_modes, seed, .. } = spec.mantle {
        spec.mantle = V2MantleSpec::On {
            mf,
            coupling,
            num_modes,
            seed,
            evolution_rate: evo_rate,
        };
    }
    spec.grid_nx = R4B_GRID;
    spec.grid_ny = R4B_GRID;
    spec.steps = R4B_N_CYCLES * R4B_K_CYCLE;
    spec.total_time_nondim = 6.0;

    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("cycle_0_s.png"))
        .expect("save init s");
    save_field_png(&init_state, V2Field::Altitude, &out_dir.join("cycle_0_altitude.png"))
        .expect("save init alt");

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r6.3/{label}] completed in {elapsed:.1}s");

    let init_field =
        Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let init_metrics = compute_metrics(0, &init_field, init_sea, 0.0, 0.0);
    let mut metrics = vec![init_metrics];

    let mut peak_v_per_cycle: Vec<f64> = Vec::with_capacity(R4B_N_CYCLES);
    let mut newton_outcomes_per_cycle: Vec<(usize, usize, usize, usize)> =
        Vec::with_capacity(R4B_N_CYCLES);
    let mut cg_mean_quartiles_per_cycle: Vec<(f64, f64, f64)> = Vec::with_capacity(R4B_N_CYCLES);

    for (i, cycle) in output.cycles.iter().enumerate() {
        let cyc_idx = i + 1;
        let v2_state = V2FinalState::from_harness(&cycle.baseline.final_state);
        save_field_png(
            &v2_state,
            V2Field::SThickness,
            &out_dir.join(format!("cycle_{cyc_idx}_s.png")),
        )
        .expect("save s");
        save_field_png(
            &v2_state,
            V2Field::Altitude,
            &out_dir.join(format!("cycle_{cyc_idx}_altitude.png")),
        )
        .expect("save alt");
        metrics.push(compute_metrics(
            cyc_idx,
            &cycle.baseline.final_state.s_field,
            cycle.common.sea_level_normalized,
            cycle.common.erosion_volume_removed,
            cycle.common.mass_drift,
        ));
        peak_v_per_cycle.push(cycle.baseline.metrics.vmax_peak);
        if let Some(n) = cycle.baseline.metrics.newton.as_ref() {
            newton_outcomes_per_cycle.push((n.converged, n.stalled, n.diverged, n.capped));
            let mut cg_per_step: Vec<usize> = Vec::with_capacity(n.outer_iters.len());
            let mut cg_cursor: usize = 0;
            for &outer in &n.outer_iters {
                let outer_n = outer as usize;
                let cg_slice = &n.cg_iters_per_newton_step
                    [cg_cursor..(cg_cursor + outer_n).min(n.cg_iters_per_newton_step.len())];
                cg_per_step.push(cg_slice.iter().sum());
                cg_cursor += outer_n;
            }
            let k = cg_per_step.len();
            if k >= 4 {
                let early = &cg_per_step[..k / 4];
                let mid = &cg_per_step[k / 4..3 * k / 4];
                let late = &cg_per_step[3 * k / 4..];
                let avg = |v: &[usize]| -> f64 {
                    if v.is_empty() {
                        0.0
                    } else {
                        v.iter().sum::<usize>() as f64 / v.len() as f64
                    }
                };
                cg_mean_quartiles_per_cycle.push((avg(early), avg(mid), avg(late)));
            } else {
                cg_mean_quartiles_per_cycle.push((0.0, 0.0, 0.0));
            }
        } else {
            newton_outcomes_per_cycle.push((0, 0, 0, 0));
            cg_mean_quartiles_per_cycle.push((0.0, 0.0, 0.0));
        }
    }

    let mut md = String::new();
    md.push_str(&format!(
        "# R6.3 — {label} (mf={mf}, evolution_rate={evo_rate})\n\n"
    ));
    md.push_str(&format!(
        "32² active_medley, workflow ON (D2 + D1-ter), 5 cycles × 20 steps. Runtime: {:.1}s.\n\n",
        elapsed
    ));
    md.push_str("## Per-cycle solver health\n\n");
    md.push_str(
        "| cycle | peak \\|v\\| | Newton C/S/D/Cap | peak S̃ | frac>0.8 | mass | max_path | CG e/m/l |\n",
    );
    md.push_str("|---|---|---|---|---|---|---|---|\n");
    for (i, _cycle) in output.cycles.iter().enumerate() {
        let m = &metrics[i + 1];
        let o = newton_outcomes_per_cycle[i];
        let q = cg_mean_quartiles_per_cycle[i];
        md.push_str(&format!(
            "| {} | {:.3e} | {}/{}/{}/{} | {:.3} | {:.3} | {:.2} | {} | {:.0}/{:.0}/{:.0} |\n",
            i + 1,
            peak_v_per_cycle[i],
            o.0, o.1, o.2, o.3,
            m.peak_s,
            m.fraction_above_0_8,
            m.mass_total,
            m.max_path_length,
            q.0, q.1, q.2,
        ));
    }
    md.push('\n');

    let init_frac = metrics[0].fraction_above_0_8;
    let final_frac = metrics.last().unwrap().fraction_above_0_8;
    let retention = if init_frac > 1e-9 { final_frac / init_frac } else { 0.0 };
    let n_dynamic = peak_v_per_cycle.iter().filter(|&&v| v > 0.1).count();
    let mass_init = metrics[0].mass_total;
    let mass_final = metrics.last().unwrap().mass_total;
    let mass_loss_pct = (mass_init - mass_final) / mass_init.abs().max(1e-12) * 100.0;
    let mass_loss_per_cycle_pct = mass_loss_pct / R4B_N_CYCLES as f64;
    let final_peak_s = metrics.last().unwrap().peak_s;
    let final_sea = metrics.last().unwrap().sea_level_ref;
    let max_max_path = metrics
        .iter()
        .skip(1)
        .map(|m| m.max_path_length)
        .max()
        .unwrap_or(0);
    let stalled_total: usize = newton_outcomes_per_cycle.iter().map(|o| o.1 + o.3).sum();
    let total_steps: usize =
        newton_outcomes_per_cycle.iter().map(|o| o.0 + o.1 + o.2 + o.3).sum();
    let stall_pct = if total_steps > 0 {
        100.0 * stalled_total as f64 / total_steps as f64
    } else {
        0.0
    };

    let pass_r4_1 = final_peak_s > final_sea;
    let pass_r4_2 = retention > 0.5;
    let pass_r4_4 = mass_loss_per_cycle_pct.abs() < 1.0;
    let pass_r4_5 = max_max_path >= 5;
    let pass_r4_6 = n_dynamic >= 3;
    let auto_pass_count = [pass_r4_1, pass_r4_2, pass_r4_4, pass_r4_5, pass_r4_6]
        .iter()
        .filter(|&&p| p)
        .count();

    md.push_str("## Multi-dim acceptance (R4.1-R4.6)\n\n");
    md.push_str(&format!(
        "- **R4.1 - Continents émergés** : peak S̃_final = {:.3} vs sea_level = {:.3} → **{}**\n",
        final_peak_s, final_sea,
        if pass_r4_1 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "- **R4.2 - Cratons préservés** : retention = {:.1} % → **{}**\n",
        100.0 * retention,
        if pass_r4_2 { "PASS" } else { "FAIL" }
    ));
    md.push_str(
        "- **R4.3 - Bordures continentales évoluées** : VISUAL — inspect `cycle_5_altitude.png`\n",
    );
    md.push_str(&format!(
        "- **R4.4 - Conservation** : mass loss/cycle = {:.3} % → **{}**\n",
        mass_loss_per_cycle_pct,
        if pass_r4_4 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "- **R4.5 - Drainage actif** : max_path = {} (cycles 1-5) → **{}**\n",
        max_max_path,
        if pass_r4_5 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "- **R4.6 - Dynamique soutenue** : peak |v| > 0.1 on {}/5 cycles → **{}**\n",
        n_dynamic,
        if pass_r4_6 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "\nStall+Cap = {:.1} % over {} total Newton outcomes.\n",
        stall_pct, total_steps
    ));
    md.push_str(&format!(
        "\nAuto-criteria PASS count: **{} / 5** (R4.3 visual pending).\n",
        auto_pass_count
    ));

    std::fs::write(out_dir.join("report.md"), &md).expect("write per-config report");
    println!("{md}");

    R6ConfigSummary {
        mf,
        evo_rate,
        label: label.to_string(),
        runtime_s: elapsed,
        peak_v_per_cycle,
        peak_s_per_cycle: metrics.iter().skip(1).map(|m| m.peak_s).collect(),
        frac_above_0_8_per_cycle: metrics.iter().skip(1).map(|m| m.fraction_above_0_8).collect(),
        mass_total_per_cycle: metrics.iter().skip(1).map(|m| m.mass_total).collect(),
        max_path_per_cycle: metrics.iter().skip(1).map(|m| m.max_path_length).collect(),
        sea_level_per_cycle: metrics.iter().skip(1).map(|m| m.sea_level_ref).collect(),
        newton_outcomes_per_cycle,
        cg_mean_quartiles_per_cycle,
        init_frac_above_0_8: metrics[0].fraction_above_0_8,
        init_peak_s: metrics[0].peak_s,
    }
}

#[derive(Clone, Copy, Debug)]
struct AutoVerdict {
    r4_1: bool,
    r4_2: bool,
    r4_4: bool,
    r4_5: bool,
    r4_6: bool,
    auto_count: u8,
}

fn evaluate_auto(s: &R6ConfigSummary) -> AutoVerdict {
    let final_peak_s = s.peak_s_per_cycle.last().copied().unwrap_or(0.0);
    let final_sea = s.sea_level_per_cycle.last().copied().unwrap_or(0.5);
    let final_frac = s.frac_above_0_8_per_cycle.last().copied().unwrap_or(0.0);
    let retention = if s.init_frac_above_0_8 > 1e-9 {
        final_frac / s.init_frac_above_0_8
    } else {
        0.0
    };
    let n_dynamic = s.peak_v_per_cycle.iter().filter(|&&v| v > 0.1).count();
    let mass_init = s.mass_total_per_cycle.first().copied().unwrap_or(1.0);
    let mass_final = s.mass_total_per_cycle.last().copied().unwrap_or(mass_init);
    let mass_loss_per_cycle_pct =
        (mass_init - mass_final) / mass_init.abs().max(1e-12) * 100.0 / R4B_N_CYCLES as f64;
    let max_max_path = s.max_path_per_cycle.iter().copied().max().unwrap_or(0);

    let r4_1 = final_peak_s > final_sea;
    let r4_2 = retention > 0.5;
    let r4_4 = mass_loss_per_cycle_pct.abs() < 1.0;
    let r4_5 = max_max_path >= 5;
    let r4_6 = n_dynamic >= 3;
    let auto_count = [r4_1, r4_2, r4_4, r4_5, r4_6].iter().filter(|&&p| p).count() as u8;
    AutoVerdict { r4_1, r4_2, r4_4, r4_5, r4_6, auto_count }
}

fn write_r6_3_summary(summaries: &[R6ConfigSummary], out_root: &Path) {
    let mut md = String::new();
    md.push_str("# Step 12 R6.3 — mantle.evolution_rate × mf sweep\n\n");
    md.push_str("32² active_medley, workflow ON (D2 + D1-ter), 5 cycles × 20 steps. ");
    md.push_str("Phase A defaults (α=0.01, rebound=0.80, max_drainage=10). ");
    md.push_str(&format!(
        "{} configs total. Generated by `r6_3_sweep_all_configs`.\n\n",
        summaries.len()
    ));

    md.push_str("## Per-config auto verdict\n\n");
    md.push_str(
        "| mf | evo | runtime | R4.1 émergés | R4.2 cratons | R4.4 conserv. | R4.5 drainage | R4.6 dynamique | auto count |\n",
    );
    md.push_str("|---|---|---|---|---|---|---|---|---|\n");
    let mut verdicts: Vec<AutoVerdict> = Vec::with_capacity(summaries.len());
    for s in summaries {
        let v = evaluate_auto(s);
        verdicts.push(v);
        let p = |b: bool| if b { "PASS" } else { "FAIL" };
        md.push_str(&format!(
            "| {:.2} | {:.2} | {:.0}s | {} | {} | {} | {} | {} | **{} / 5** |\n",
            s.mf, s.evo_rate, s.runtime_s,
            p(v.r4_1), p(v.r4_2), p(v.r4_4), p(v.r4_5), p(v.r4_6),
            v.auto_count,
        ));
    }
    md.push('\n');

    md.push_str("## Per-config peak|v| series (cycles 1-5)\n\n");
    md.push_str("| mf | evo | c1 | c2 | c3 | c4 | c5 | n cycles > 0.1 |\n");
    md.push_str("|---|---|---|---|---|---|---|---|\n");
    for s in summaries {
        let n_dyn = s.peak_v_per_cycle.iter().filter(|&&v| v > 0.1).count();
        let vs: Vec<String> = s.peak_v_per_cycle.iter().map(|v| format!("{:.3e}", v)).collect();
        let row: String = vs.join(" | ");
        md.push_str(&format!(
            "| {:.2} | {:.2} | {} | **{}** |\n",
            s.mf, s.evo_rate, row, n_dyn,
        ));
    }
    md.push('\n');

    md.push_str("## Per-config preservation (frac>0.8 retention)\n\n");
    md.push_str("| mf | evo | init | c1 | c2 | c3 | c4 | c5 | retention |\n");
    md.push_str("|---|---|---|---|---|---|---|---|---|\n");
    for s in summaries {
        let final_frac = s.frac_above_0_8_per_cycle.last().copied().unwrap_or(0.0);
        let retention = if s.init_frac_above_0_8 > 1e-9 {
            final_frac / s.init_frac_above_0_8
        } else {
            0.0
        };
        let vs: Vec<String> = s
            .frac_above_0_8_per_cycle
            .iter()
            .map(|v| format!("{:.3}", v))
            .collect();
        let row: String = vs.join(" | ");
        md.push_str(&format!(
            "| {:.2} | {:.2} | {:.3} | {} | {:.1} % |\n",
            s.mf, s.evo_rate, s.init_frac_above_0_8, row, 100.0 * retention,
        ));
    }
    md.push('\n');

    md.push_str("## Reference — R5b evo=0 baseline (from earlier sweeps)\n\n");
    md.push_str(
        "These cells are the R6.3 anchor for `evo=0`. They are NOT re-run here \
         (evo=0 is bit-identical to the pre-R6 mantle path under the R6.2 wiring). \
         Source reports referenced below.\n\n",
    );
    md.push_str("| mf | evo | source | preservation | peak \\|v\\|_max | notes |\n");
    md.push_str("|---|---|---|---|---|---|\n");
    md.push_str(
        "| 0.50 | 0.00 | `step12_r5b_mf_sweep_post_d1_ter/mf_0_5/` | preserved | ~ 1.2e-3 | quasi-static (R5b finding) |\n",
    );
    md.push_str(
        "| 0.70 | 0.00 | `step12_r5b_mf_sweep_post_d1_ter/mf_0_7/` | preserved | ~ 1 cycle > 0.1 | borderline dynamic |\n",
    );
    md.push_str(
        "| 0.80 | 0.00 | `step12_r5b_mf_sweep_post_d1_ter/mf_0_8/` | transition | mixed | transition pathological (R5b) |\n",
    );
    md.push_str(
        "| 1.00 | 0.00 | `step12_r4_visual_checkpoint/active_medley_post_d1_ter/` | dissolves | high | continents dissolved (R4) |\n",
    );
    md.push('\n');

    let auto_passers: Vec<(&R6ConfigSummary, AutoVerdict)> = summaries
        .iter()
        .zip(verdicts.iter())
        .filter(|(_, v)| v.auto_count == 5)
        .map(|(s, v)| (s, *v))
        .collect();
    md.push_str("## EVO verdict (auto-criteria only — R4.3 visual still required)\n\n");
    match auto_passers.len() {
        0 => {
            md.push_str(
                "**Provisional EVO.D** — no config passes all 5 auto-criteria. \
                 Inspect `R4.3 visual` per config; if no chains observed either, \
                 confirm EVO.D and consider Phys.C pivot in Step 12.X. If R4.3 \
                 visual passes but auto-fail comes from a single criterion only, \
                 re-evaluate as marginal (no tuning).\n\n",
            );
        }
        1 => {
            let (s, _) = &auto_passers[0];
            md.push_str(&format!(
                "**Provisional EVO.A** — single config passes all 5 auto-criteria: \
                 `mf={:.2}, evo={:.2}`. Visual R4.3 inspection required on \
                 `step12_r6_3_sweep/{}/cycle_5_altitude.png`. If chains visible \
                 and migration coherent → 64² validation gate.\n\n",
                s.mf, s.evo_rate, s.label,
            ));
        }
        _ => {
            md.push_str(&format!(
                "**Provisional EVO.B** — {} configs pass all 5 auto-criteria. \
                 Priority for 64² validation (user rule): mf=0.7 + median evo. \
                 Visual R4.3 inspection on each. Candidates:\n\n",
                auto_passers.len(),
            ));
            for (s, _) in &auto_passers {
                md.push_str(&format!("- mf={:.2}, evo={:.2} ({})\n", s.mf, s.evo_rate, s.label));
            }
            md.push('\n');
        }
    }

    std::fs::write(out_root.join("summary.md"), &md).expect("write summary");
    println!("\n[r6.3] wrote summary.md");
}

/// R6.3 master sweep — 12 configs run sequentially. Estimated ~2h compute.
///
/// ```bash
/// cargo test --release -p ymir-viz --test v2_workflow_r4_visual_checkpoint \
///     r6_3_sweep_all_configs -- --ignored --nocapture --test-threads=1
/// ```
#[test]
#[ignore]
fn r6_3_sweep_all_configs() {
    let out_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r6_3_sweep");
    std::fs::create_dir_all(&out_root).expect("create r6_3 root");

    let configs: Vec<(f64, f64, &str)> = vec![
        (0.5, 0.05, "mf_0_5_evo_0_05"),
        (0.5, 0.10, "mf_0_5_evo_0_10"),
        (0.5, 0.20, "mf_0_5_evo_0_20"),
        (0.7, 0.05, "mf_0_7_evo_0_05"),
        (0.7, 0.10, "mf_0_7_evo_0_10"),
        (0.7, 0.20, "mf_0_7_evo_0_20"),
        (0.8, 0.05, "mf_0_8_evo_0_05"),
        (0.8, 0.10, "mf_0_8_evo_0_10"),
        (0.8, 0.20, "mf_0_8_evo_0_20"),
        (1.0, 0.05, "mf_1_0_evo_0_05"),
        (1.0, 0.10, "mf_1_0_evo_0_10"),
        (1.0, 0.20, "mf_1_0_evo_0_20"),
    ];

    let n = configs.len();
    let t_total = std::time::Instant::now();
    let mut summaries: Vec<R6ConfigSummary> = Vec::with_capacity(n);
    for (idx, (mf, evo, label)) in configs.iter().enumerate() {
        println!(
            "\n[r6.3] ===== {}/{} : mf={mf}, evo={evo} ({label}) =====",
            idx + 1,
            n,
        );
        let s = run_r6_3_config(*mf, *evo, label);
        summaries.push(s);
        // Write summary after every config so partial progress is visible
        // if the run is interrupted.
        write_r6_3_summary(&summaries, &out_root);
        println!(
            "[r6.3] cumulative elapsed: {:.1}s ({:.1} min)",
            t_total.elapsed().as_secs_f64(),
            t_total.elapsed().as_secs_f64() / 60.0,
        );
    }
    println!(
        "\n[r6.3] sweep done. {} configs in {:.1} min total.",
        n,
        t_total.elapsed().as_secs_f64() / 60.0,
    );
}

/// Single-config smoke test for the R6.3 scaffolding. Use before launching
/// the full sweep to catch wiring/compile bugs cheaply (~10 min).
///
/// ```bash
/// cargo test --release -p ymir-viz --test v2_workflow_r4_visual_checkpoint \
///     r6_3_smoke_single_config -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn r6_3_smoke_single_config() {
    let _ = run_r6_3_config(0.7, 0.10, "smoke_mf_0_7_evo_0_10");
}

/// R6.3 Option β — focused 2-config probe at mf=1.0 × evo ∈ {0.10, 0.20}.
///
/// Triggered by the smoke finding (mf=0.7, evo=0.10 → adiabatic equilibrium,
/// R4.6 dynamique FAIL). Hypothesis: at higher pattern amplitude (mf=1.0),
/// the system can no longer relax fast enough to track the drifting pattern
/// → sustained dynamics. Tests this with two evo points before declaring
/// EVO.D or escalating to a refined sweep.
///
/// Budget: ~50 min wall time (25 min × 2 configs from smoke runtime).
///
/// Decision tree:
/// - β.1 — sustained dynamics + mass loss < 5 % → EVO.A probable, recommend
///         refined mini-sweep mf ∈ {0.85..1.0} × evo ∈ {0.10, 0.20}
/// - β.2 — sustained dynamics but mass loss > 10 % → EVO.C, escalate to
///         Step 13.6 (stronger cratons) decision
/// - β.3 — peak |v| collapse at mf=1.0 too → EVO.D confirmed, pivot Phys.C
///         in Step 12.X
///
/// ```bash
/// cargo test --release -p ymir-viz --test v2_workflow_r4_visual_checkpoint \
///     r6_3_option_b_mf_1_0 -- --ignored --nocapture --test-threads=1
/// ```
#[test]
#[ignore]
fn r6_3_option_b_mf_1_0() {
    let out_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r6_3_sweep");
    std::fs::create_dir_all(&out_root).expect("create r6_3 root");

    let configs: Vec<(f64, f64, &str)> = vec![
        (1.0, 0.10, "mf_1_0_evo_0_10"),
        (1.0, 0.20, "mf_1_0_evo_0_20"),
    ];

    let t_total = std::time::Instant::now();
    let mut summaries: Vec<R6ConfigSummary> = Vec::with_capacity(configs.len());
    for (idx, (mf, evo, label)) in configs.iter().enumerate() {
        println!(
            "\n[r6.3-β] ===== {}/{} : mf={mf}, evo={evo} ({label}) =====",
            idx + 1,
            configs.len(),
        );
        let s = run_r6_3_config(*mf, *evo, label);
        summaries.push(s);
        write_r6_3_summary(&summaries, &out_root);
        println!(
            "[r6.3-β] cumulative elapsed: {:.1}s ({:.1} min)",
            t_total.elapsed().as_secs_f64(),
            t_total.elapsed().as_secs_f64() / 60.0,
        );
    }
    println!(
        "\n[r6.3-β] Option β done. {} configs in {:.1} min total.",
        configs.len(),
        t_total.elapsed().as_secs_f64() / 60.0,
    );
}

// ----------------------------------------------------------------------------
// Step 12 R6.3 Option C.4 — cratonic resistance amplification probe.
//
// Triggered by Option β verdict β.2/EVO.C: at mf=1.0, evo > 0 unlocks
// sustained dynamics but cratons retain only 31-35 % (vs 50 % target) and
// mass loss is ~2 %/cycle. Hypothesis: amplifying the existing cratonic
// resistance mechanisms (K viscous + B yield) lets cratons survive the
// mantle-pull dynamics that produce the wanted tectonic flow.
//
// Approach (minimum-surface): scale the existing `k_viscous` and `b_factor`
// fields of `V2CratonicSpec::On` by a `cratonic_amp` multiplier. No core
// library change; `amp = 1.0` is bit-identical to pre-C.4 behaviour.
//
//   amp=1   → K=5,    B=8    (defaults — bit-identical to Option β)
//   amp=2   → K=10,   B=16   (gate test, V1 vigilance)
//   amp=5   → K=25,   B=40
//   amp=10  → K=50,   B=80   (saturation-risk; CG cap probable)
//
// V1 — amp=2 first as a saturation gate. If CG mean > 1.5× the Option β
// baseline (~13000/cycle) OR retention stays below the β baseline, halt
// and remontée before running amp=5/10. Equivalent budget per config to
// Option β (~35 min).
//
// Output: `docs/reports/step12_r6_3_c4_craton_amp/<label>/` plus
// `summary.md` consolidating amp-by-amp comparison + verdict.

/// Run an R6.3 config with cratonic resistance scaled by `cratonic_amp`.
/// `cratonic_amp = 1.0` is bit-identical to `run_r6_3_config` at the same
/// (mf, evo).
fn run_r6_3_c4_config(
    mf: f64,
    evo_rate: f64,
    cratonic_amp: f64,
    label: &str,
) -> R6ConfigSummary {
    use ymir_viz::bridge::v2::build_config;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r6_3_c4_craton_amp")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!(
        "\n=== R6.3 C.4: {label} (mf={mf}, evo={evo_rate}, craton_amp={cratonic_amp}) ===",
    );

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    if let V2MantleSpec::On { coupling, num_modes, seed, .. } = spec.mantle {
        spec.mantle = V2MantleSpec::On {
            mf,
            coupling,
            num_modes,
            seed,
            evolution_rate: evo_rate,
        };
    }
    // Scale K and B by `cratonic_amp`. At amp=1.0 the spec is byte-equal
    // to the active_medley default cratonic block, preserving pre-C.4
    // regression for that case.
    if let ymir_viz::bridge::v2::V2CratonicSpec::On {
        cr,
        k_viscous,
        b_factor,
        smoothing_width,
        plate_area_min,
    } = spec.cratonic
    {
        spec.cratonic = ymir_viz::bridge::v2::V2CratonicSpec::On {
            cr,
            k_viscous: k_viscous * cratonic_amp,
            b_factor: b_factor * cratonic_amp,
            smoothing_width,
            plate_area_min,
        };
        println!(
            "[c.4/{label}] cratonic K {} → {:.1}, B {} → {:.1}",
            k_viscous,
            k_viscous * cratonic_amp,
            b_factor,
            b_factor * cratonic_amp,
        );
    }
    spec.grid_nx = R4B_GRID;
    spec.grid_ny = R4B_GRID;
    spec.steps = R4B_N_CYCLES * R4B_K_CYCLE;
    spec.total_time_nondim = 6.0;

    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("cycle_0_s.png"))
        .expect("save init s");
    save_field_png(&init_state, V2Field::Altitude, &out_dir.join("cycle_0_altitude.png"))
        .expect("save init alt");

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[c.4/{label}] completed in {elapsed:.1}s");

    let init_field =
        Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let init_metrics = compute_metrics(0, &init_field, init_sea, 0.0, 0.0);
    let mut metrics = vec![init_metrics];

    let mut peak_v_per_cycle: Vec<f64> = Vec::with_capacity(R4B_N_CYCLES);
    let mut newton_outcomes_per_cycle: Vec<(usize, usize, usize, usize)> =
        Vec::with_capacity(R4B_N_CYCLES);
    let mut cg_mean_quartiles_per_cycle: Vec<(f64, f64, f64)> = Vec::with_capacity(R4B_N_CYCLES);

    for (i, cycle) in output.cycles.iter().enumerate() {
        let cyc_idx = i + 1;
        let v2_state = V2FinalState::from_harness(&cycle.baseline.final_state);
        save_field_png(
            &v2_state,
            V2Field::SThickness,
            &out_dir.join(format!("cycle_{cyc_idx}_s.png")),
        )
        .expect("save s");
        save_field_png(
            &v2_state,
            V2Field::Altitude,
            &out_dir.join(format!("cycle_{cyc_idx}_altitude.png")),
        )
        .expect("save alt");
        metrics.push(compute_metrics(
            cyc_idx,
            &cycle.baseline.final_state.s_field,
            cycle.common.sea_level_normalized,
            cycle.common.erosion_volume_removed,
            cycle.common.mass_drift,
        ));
        peak_v_per_cycle.push(cycle.baseline.metrics.vmax_peak);
        if let Some(n) = cycle.baseline.metrics.newton.as_ref() {
            newton_outcomes_per_cycle.push((n.converged, n.stalled, n.diverged, n.capped));
            let mut cg_per_step: Vec<usize> = Vec::with_capacity(n.outer_iters.len());
            let mut cg_cursor: usize = 0;
            for &outer in &n.outer_iters {
                let outer_n = outer as usize;
                let cg_slice = &n.cg_iters_per_newton_step
                    [cg_cursor..(cg_cursor + outer_n).min(n.cg_iters_per_newton_step.len())];
                cg_per_step.push(cg_slice.iter().sum());
                cg_cursor += outer_n;
            }
            let k = cg_per_step.len();
            if k >= 4 {
                let early = &cg_per_step[..k / 4];
                let mid = &cg_per_step[k / 4..3 * k / 4];
                let late = &cg_per_step[3 * k / 4..];
                let avg = |v: &[usize]| -> f64 {
                    if v.is_empty() { 0.0 } else { v.iter().sum::<usize>() as f64 / v.len() as f64 }
                };
                cg_mean_quartiles_per_cycle.push((avg(early), avg(mid), avg(late)));
            } else {
                cg_mean_quartiles_per_cycle.push((0.0, 0.0, 0.0));
            }
        } else {
            newton_outcomes_per_cycle.push((0, 0, 0, 0));
            cg_mean_quartiles_per_cycle.push((0.0, 0.0, 0.0));
        }
    }

    let mut md = String::new();
    md.push_str(&format!(
        "# R6.3 C.4 — {label} (mf={mf}, evo={evo_rate}, craton_amp={cratonic_amp})\n\n"
    ));
    md.push_str(&format!(
        "32² active_medley, workflow ON (D2 + D1-ter), 5 cycles × 20 steps. \
         K_viscous={:.1}, B_factor={:.1} (amp={cratonic_amp}). Runtime: {:.1}s.\n\n",
        5.0 * cratonic_amp,
        8.0 * cratonic_amp,
        elapsed,
    ));
    md.push_str("## Per-cycle solver health\n\n");
    md.push_str(
        "| cycle | peak \\|v\\| | Newton C/S/D/Cap | peak S̃ | frac>0.8 | mass | max_path | CG e/m/l |\n",
    );
    md.push_str("|---|---|---|---|---|---|---|---|\n");
    for (i, _cycle) in output.cycles.iter().enumerate() {
        let m = &metrics[i + 1];
        let o = newton_outcomes_per_cycle[i];
        let q = cg_mean_quartiles_per_cycle[i];
        md.push_str(&format!(
            "| {} | {:.3e} | {}/{}/{}/{} | {:.3} | {:.3} | {:.2} | {} | {:.0}/{:.0}/{:.0} |\n",
            i + 1,
            peak_v_per_cycle[i],
            o.0, o.1, o.2, o.3,
            m.peak_s,
            m.fraction_above_0_8,
            m.mass_total,
            m.max_path_length,
            q.0, q.1, q.2,
        ));
    }
    md.push('\n');

    let init_frac = metrics[0].fraction_above_0_8;
    let final_frac = metrics.last().unwrap().fraction_above_0_8;
    let retention = if init_frac > 1e-9 { final_frac / init_frac } else { 0.0 };
    let n_dynamic = peak_v_per_cycle.iter().filter(|&&v| v > 0.1).count();
    let mass_init = metrics[0].mass_total;
    let mass_final = metrics.last().unwrap().mass_total;
    let mass_loss_pct = (mass_init - mass_final) / mass_init.abs().max(1e-12) * 100.0;
    let mass_loss_per_cycle_pct = mass_loss_pct / R4B_N_CYCLES as f64;
    let final_peak_s = metrics.last().unwrap().peak_s;
    let final_sea = metrics.last().unwrap().sea_level_ref;
    let max_max_path = metrics
        .iter()
        .skip(1)
        .map(|m| m.max_path_length)
        .max()
        .unwrap_or(0);
    let stalled_total: usize = newton_outcomes_per_cycle.iter().map(|o| o.1 + o.3).sum();
    let total_steps: usize =
        newton_outcomes_per_cycle.iter().map(|o| o.0 + o.1 + o.2 + o.3).sum();
    let stall_pct = if total_steps > 0 {
        100.0 * stalled_total as f64 / total_steps as f64
    } else {
        0.0
    };
    let cg_mean_run: f64 = cg_mean_quartiles_per_cycle
        .iter()
        .map(|(a, b, c)| (a + b + c) / 3.0)
        .sum::<f64>()
        / cg_mean_quartiles_per_cycle.len().max(1) as f64;

    // 5 critères C.4 (user spec) :
    //   1. cratons retention > 50 %
    //   2. mass loss < 5 % total (= 1 % / cycle)
    //   3. dynamique soutenue (peak |v| > 0.1 on ≥ 3 cycles)
    //   4. pas de régression solveur (CG iter mean < 1000)  ← user threshold
    //   5. visual (manual)
    let pass_c1 = retention > 0.5;
    let pass_c2 = mass_loss_pct.abs() < 5.0;
    let pass_c3 = n_dynamic >= 3;
    let pass_c4 = cg_mean_run < 1000.0;
    let auto_count = [pass_c1, pass_c2, pass_c3, pass_c4]
        .iter()
        .filter(|&&p| p)
        .count();

    md.push_str("## C.4 acceptance (5 criteria — last is visual)\n\n");
    md.push_str(&format!(
        "1. **Cratons retention > 50 %** : {:.1} % → **{}**\n",
        100.0 * retention,
        if pass_c1 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "2. **Mass loss < 5 % cumul** : {:.2} % ({:.3} %/cycle) → **{}**\n",
        mass_loss_pct,
        mass_loss_per_cycle_pct,
        if pass_c2 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "3. **Dynamique soutenue (≥ 3 cycles > 0.1)** : {}/5 → **{}**\n",
        n_dynamic,
        if pass_c3 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "4. **Pas de régression solveur (CG mean < 1000)** : {:.0} → **{}**\n",
        cg_mean_run,
        if pass_c4 { "PASS" } else { "FAIL" }
    ));
    md.push_str("5. **Visual** : inspect `cycle_5_altitude.png` (chains + bordures + cratons stables)\n");
    md.push_str(&format!(
        "\nAuxiliary: peak S̃ final = {:.3} vs sea {:.3}; max_path = {}; Stall+Cap = {:.1} % over {} outcomes.\n",
        final_peak_s, final_sea, max_max_path, stall_pct, total_steps,
    ));
    md.push_str(&format!(
        "\nAuto-criteria PASS count: **{} / 4** (visual #5 pending).\n",
        auto_count
    ));

    std::fs::write(out_dir.join("report.md"), &md).expect("write report");
    println!("{md}");

    R6ConfigSummary {
        mf,
        evo_rate,
        label: label.to_string(),
        runtime_s: elapsed,
        peak_v_per_cycle,
        peak_s_per_cycle: metrics.iter().skip(1).map(|m| m.peak_s).collect(),
        frac_above_0_8_per_cycle: metrics.iter().skip(1).map(|m| m.fraction_above_0_8).collect(),
        mass_total_per_cycle: metrics.iter().skip(1).map(|m| m.mass_total).collect(),
        max_path_per_cycle: metrics.iter().skip(1).map(|m| m.max_path_length).collect(),
        sea_level_per_cycle: metrics.iter().skip(1).map(|m| m.sea_level_ref).collect(),
        newton_outcomes_per_cycle,
        cg_mean_quartiles_per_cycle,
        init_frac_above_0_8: metrics[0].fraction_above_0_8,
        init_peak_s: metrics[0].peak_s,
    }
}

/// C.4 gate — amp=2 single run. Triggers V1 vigilance check before amp=5/10.
///
/// ```bash
/// cargo test --release -p ymir-viz --test v2_workflow_r4_visual_checkpoint \
///     r6_3_c4_amp_2 -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn r6_3_c4_amp_2() {
    let _ = run_r6_3_c4_config(1.0, 0.10, 2.0, "mf_1_0_evo_0_10_amp_2");
}

/// C.4 amp=5 — run only after amp=2 gate passes (CG mean < 1000 AND
/// retention improves vs Option β baseline).
#[test]
#[ignore]
fn r6_3_c4_amp_5() {
    let _ = run_r6_3_c4_config(1.0, 0.10, 5.0, "mf_1_0_evo_0_10_amp_5");
}

/// C.4 amp=10 — run only after amp=5 passes the gate.
#[test]
#[ignore]
fn r6_3_c4_amp_10() {
    let _ = run_r6_3_c4_config(1.0, 0.10, 10.0, "mf_1_0_evo_0_10_amp_10");
}

// ----------------------------------------------------------------------------
// Step 12 R7.A.1.2 step e — init-only visual preview.
//
// Builds the t=0 S̃ field for three init-mode variants under the
// active_medley plate seed (no simulation, no Stokes solve) and writes
// the corresponding S̃/altitude PNGs side-by-side. Cheap (< 1 s).
//
//   - active_medley (Uniform, default)            — reference shape of plates
//   - active_medley + RadialProfile (R7.A.1 Run A) — Step 13 baseline
//   - active_medley + Orogenic (R7.A.1 Run B)     — Step 12 R7.A.1 test
//
// Visual checks (manual, R7.A.1.2 remontée):
//   - Crest identifiable on the Orogenic PNG ?
//   - Oriented along plate's principal axis (not aligned with grid x) ?
//   - Degenerate plates (small / colinear) fall back to angle_rad=0 ?
//   - No artefacts from wrap (plates spanning the torus boundary)?
//
// Output: `docs/reports/step12_r7_a_orogenic_profile/init_preview/`.

#[test]
#[ignore]
fn r7_a_1_init_preview() {
    use ymir_viz::bridge::v2::{V2InitModeSpec, V2OrogenicOrientation};

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r7_a_orogenic_profile/init_preview");
    std::fs::create_dir_all(&out_dir).expect("create out dir");

    println!("\n=== R7.A.1 init preview ===");

    // Baseline preset — active_medley loads with V2InitModeSpec::Uniform
    // by default (no `init_mode` field in the JSON ⇒ Default trait).
    let base_spec = presets::load("active_medley").expect("load active_medley");

    // Variant 1 — Uniform (reference / plate shape only).
    {
        let mut spec = base_spec.clone();
        spec.grid_nx = R4B_GRID;
        spec.grid_ny = R4B_GRID;
        let state = make_init_state(&spec);
        save_field_png(&state, V2Field::SThickness, &out_dir.join("uniform_s.png"))
            .expect("save uniform s");
        save_field_png(&state, V2Field::Altitude, &out_dir.join("uniform_altitude.png"))
            .expect("save uniform altitude");
        println!("[r7.a.1] wrote uniform_s.png / uniform_altitude.png");
    }

    // Variant 2 — RadialProfile (Run A baseline for R7.A.1.3).
    {
        let mut spec = base_spec.clone();
        spec.grid_nx = R4B_GRID;
        spec.grid_ny = R4B_GRID;
        spec.init_mode = V2InitModeSpec::radial_profile_default();
        spec.s_perturbation_amplitude = 0.0;
        let state = make_init_state(&spec);
        save_field_png(&state, V2Field::SThickness, &out_dir.join("radial_s.png"))
            .expect("save radial s");
        save_field_png(&state, V2Field::Altitude, &out_dir.join("radial_altitude.png"))
            .expect("save radial altitude");
        println!("[r7.a.1] wrote radial_s.png / radial_altitude.png");
    }

    // Variant 3 — Orogenic (Run B test for R7.A.1.3).
    {
        let mut spec = base_spec;
        spec.grid_nx = R4B_GRID;
        spec.grid_ny = R4B_GRID;
        spec.init_mode = V2InitModeSpec::Orogenic {
            peak_value: 1.20,
            base_continental_value: 0.85,
            oceanic_value: 0.20,
            half_length_ratio: 0.40,
            width_sigma_ratio: 0.08,
            orientation: V2OrogenicOrientation::PlateMainAxisPca,
        };
        spec.s_perturbation_amplitude = 0.0;
        let state = make_init_state(&spec);
        save_field_png(&state, V2Field::SThickness, &out_dir.join("orogenic_s.png"))
            .expect("save orogenic s");
        save_field_png(&state, V2Field::Altitude, &out_dir.join("orogenic_altitude.png"))
            .expect("save orogenic altitude");
        // Sanity: max(S̃) over the orogenic field must exceed both the
        // base (0.85) and the radial peak (0.95) — confirms the ridge
        // amplitude is wired through and not clamped to the radial
        // ceiling.
        let max_s_oro: f64 = state.s_field.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("[r7.a.1] orogenic max(S̃) = {max_s_oro:.4} (peak target 1.20)");
        assert!(
            max_s_oro > 1.0,
            "orogenic max(S̃) = {max_s_oro} ≤ 1.0 — ridge did not form; \
             check periodic centroid / PCA / formula wiring",
        );
        // Also write a 64² version for the R7.A.1.3 sim grid.
        let mut spec64 = presets::load("active_medley").expect("reload medley");
        spec64.grid_nx = 64;
        spec64.grid_ny = 64;
        spec64.init_mode = spec.init_mode;
        spec64.s_perturbation_amplitude = 0.0;
        let state64 = make_init_state(&spec64);
        save_field_png(&state64, V2Field::SThickness, &out_dir.join("orogenic_s_64.png"))
            .expect("save 64 s");
        save_field_png(
            &state64,
            V2Field::Altitude,
            &out_dir.join("orogenic_altitude_64.png"),
        )
        .expect("save 64 altitude");
        let max_s_oro_64: f64 = state64
            .s_field
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        println!("[r7.a.1] orogenic 64² max(S̃) = {max_s_oro_64:.4}");
        println!("[r7.a.1] wrote orogenic_s.png / orogenic_altitude.png + 64 variants");
    }

    // Variant 4 — Orogenic 64² with width_sigma_ratio = 0.10
    // (Himalaya-like ~10 % of plate diameter, user proposal (I) for
    // R7.A.1.3). Compared to the σ=0.08 baseline above, σ=0.10 gives
    // ~25 % wider transverse profile and is the candidate for the
    // simulation run pre-validation visual.
    {
        let mut spec64w = presets::load("active_medley").expect("reload medley");
        spec64w.grid_nx = 64;
        spec64w.grid_ny = 64;
        spec64w.init_mode = V2InitModeSpec::Orogenic {
            peak_value: 1.20,
            base_continental_value: 0.85,
            oceanic_value: 0.20,
            half_length_ratio: 0.40,
            width_sigma_ratio: 0.10,
            orientation: V2OrogenicOrientation::PlateMainAxisPca,
        };
        spec64w.s_perturbation_amplitude = 0.0;
        let state64w = make_init_state(&spec64w);
        save_field_png(
            &state64w,
            V2Field::SThickness,
            &out_dir.join("orogenic_s_64_sigma_10.png"),
        )
        .expect("save 64 sigma=0.10 s");
        save_field_png(
            &state64w,
            V2Field::Altitude,
            &out_dir.join("orogenic_altitude_64_sigma_10.png"),
        )
        .expect("save 64 sigma=0.10 altitude");
        let max_s_64_w: f64 = state64w
            .s_field
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        println!(
            "[r7.a.1] orogenic 64² σ=0.10 max(S̃) = {max_s_64_w:.4} (peak target 1.20)",
        );
        println!("[r7.a.1] wrote orogenic_*_64_sigma_10.png");
    }

    println!("[r7.a.1] init preview done. Inspect {}", out_dir.display());
}

/// Render a fixed-global-scale "altitude" PNG that bypasses
/// `compute_isostasy`'s per-image `h_min`/`h_max` rescaling. Maps the
/// raw S̃ field through the isostasy piecewise (with **fixed**
/// `s_min`/`s_max` bounds shared across all renders) and the standard
/// hypsometric colormap. The point is to make composite ↔ radial ↔
/// orogenic PNGs directly comparable: the same S̃ value yields the
/// same color in every output (user's R7.A.2.3 finding).
///
/// `s_min` / `s_max` are S̃ bounds (not isostasy heights). Pick
/// `s_max ≥ orogenic_peak + margin` so the composite cap lands inside
/// the colormap's "mountain" band rather than saturating to white.
fn save_altitude_fixed_scale(
    state: &V2FinalState,
    path: &Path,
    s_min: f64,
    s_max: f64,
) -> std::io::Result<()> {
    use image::{ColorType, save_buffer};
    use ymir_core::tectonics::isostasy::IsostasyConfig;
    use ymir_viz::visualization::colormap::hypsometric_colormap;
    let cfg = IsostasyConfig::default();
    let buoyancy = 1.0 - (cfg.rho_crust / cfg.rho_mantle) as f64;
    let h_min = s_min * buoyancy;
    let h_max = s_max * buoyancy;
    let h_range = (h_max - h_min).max(1e-10);
    let h_sea = h_min + cfg.sea_level_fraction as f64 * h_range;
    let sea_norm = (cfg.max_depth_m / (cfg.max_depth_m + cfg.max_elevation_m)) as f64;
    let sea_clamped = sea_norm.clamp(1e-6, 1.0 - 1e-6);
    let nx = state.nx;
    let ny = state.ny;
    let mut rgba = vec![0u8; nx * ny * 4];
    for j in 0..ny {
        for i in 0..nx {
            let s = state.s_field[j * nx + i];
            let h_raw = (s * buoyancy).clamp(h_min, h_max);
            let normalized = if h_raw <= h_sea {
                let t = (h_raw - h_min) / (h_sea - h_min).max(1e-10);
                t * sea_norm
            } else {
                let t = (h_raw - h_sea) / (h_max - h_sea).max(1e-10);
                sea_norm + t * (1.0 - sea_norm)
            };
            let altitude = if normalized <= sea_clamped {
                0.5 * normalized / sea_clamped
            } else {
                0.5 + 0.5 * (normalized - sea_clamped) / (1.0 - sea_clamped)
            };
            let color = hypsometric_colormap(altitude.clamp(0.0, 1.0));
            // Match field_to_rgba: image row 0 maps to grid row
            // (ny - 1 - j). Without this Y-flip the *_altitude_fixed.png
            // renders are vertically inverted relative to the SThickness
            // and adaptive-altitude PNGs (user-reported R7.A.2.3
            // follow-up).
            let img_row = ny - 1 - j;
            let off = (img_row * nx + i) * 4;
            rgba[off..off + 4].copy_from_slice(&color);
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    save_buffer(path, &rgba, nx as u32, ny as u32, ColorType::Rgba8)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
    Ok(())
}

/// Numerical sanity stats for an init S̃ field. Printed alongside the
/// PNGs so the visual inspection can be cross-checked against
/// quantitative thresholds — protects against palette-adaptive
/// artefacts that would make e.g. composite "look less elevated than
/// radial" purely because its colour mapping is rescaled.
fn print_init_stats(label: &str, state: &V2FinalState) {
    let s = &state.s_field;
    let n = s.len();
    let max = s.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = s.iter().cloned().fold(f64::INFINITY, f64::min);
    let mean: f64 = s.iter().sum::<f64>() / n as f64;
    let count = |thr: f64| -> f64 {
        s.iter().filter(|&&v| v > thr).count() as f64 / n as f64
    };
    println!(
        "[r7.a.2/stats] {label:30} min={:.4} mean={:.4} max={:.4} | \
         f>0.85={:.3} f>0.95={:.3} f>1.00={:.3} f>1.10={:.3} f>1.15={:.3}",
        min,
        mean,
        max,
        count(0.85),
        count(0.95),
        count(1.00),
        count(1.10),
        count(1.15),
    );
}

/// R7.A.2.3 init-only visual preview — Composite mode vs the three
/// baselines under the same Voronoï layout at 64² (Living Landz
/// target resolution). The simulation run R7.A.2.4 is gated on this
/// preview passing the user's visual acceptance (a)+(b)+(c) criteria:
/// dome identifiable, ridge superposed distinct, piémont visible.
///
/// Output: `docs/reports/step12_r7_a_composite_profile/init_preview/`
///   - uniform_{s,altitude}.png      — reference plate shape
///   - radial_{s,altitude}.png       — Step 13 RadialProfile (Run A baseline)
///   - orogenic_{s,altitude}.png     — R7.A.1 Orogenic-seul (Run B comparator)
///   - composite_{s,altitude}.png    — R7.A.2 Composite (Run C candidate)
///   - *_altitude_fixed.png          — global-scale altitude (palette
///                                     bounds shared across runs;
///                                     adaptive per-image bounds is
///                                     the diagnostic gotcha that
///                                     prompted R7.A.2.3 follow-up).
///
/// SThickness ("*_s.png") already uses the fixed global palette
/// `[0, 2.5]` (see `field_to_rgba` in `v2_viz.rs`) — those PNGs are
/// comparable across runs as-is.
///
/// ```bash
/// cargo test --release -p ymir-viz --test v2_workflow_r4_visual_checkpoint \
///     r7_a_2_init_preview -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn r7_a_2_init_preview() {
    use ymir_viz::bridge::v2::V2InitModeSpec;

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r7_a_composite_profile/init_preview");
    std::fs::create_dir_all(&out_dir).expect("create out dir");

    println!("\n=== R7.A.2 init preview ===");

    let base_spec = presets::load("active_medley").expect("load active_medley");

    // Global altitude palette bounds (used by `save_altitude_fixed_scale`
    // below). `s_max = 1.5` keeps the composite cap (1.20) inside the
    // colormap's "mountain" band without saturating to white, leaving
    // headroom to spot anomalies above the cap if the wiring ever
    // misfires.
    let s_min_global: f64 = 0.0;
    let s_max_global: f64 = 1.5;
    println!(
        "[r7.a.2] altitude palette: SThickness adaptive=[0,2.5] (global, OK), \
         Altitude adaptive (BUGGY per-image isostasy), \
         altitude_fixed in this test = [{s_min_global}, {s_max_global}]",
    );

    // Variant 1 — Uniform reference (plate shape only).
    {
        let mut spec = base_spec.clone();
        spec.grid_nx = 64;
        spec.grid_ny = 64;
        let state = make_init_state(&spec);
        save_field_png(&state, V2Field::SThickness, &out_dir.join("uniform_s.png"))
            .expect("save uniform s");
        save_field_png(&state, V2Field::Altitude, &out_dir.join("uniform_altitude.png"))
            .expect("save uniform altitude");
        save_altitude_fixed_scale(
            &state,
            &out_dir.join("uniform_altitude_fixed.png"),
            s_min_global,
            s_max_global,
        )
        .expect("save uniform altitude fixed");
        print_init_stats("uniform", &state);
    }

    // Variant 2 — RadialProfile (Run A baseline).
    {
        let mut spec = base_spec.clone();
        spec.grid_nx = 64;
        spec.grid_ny = 64;
        spec.init_mode = V2InitModeSpec::radial_profile_default();
        spec.s_perturbation_amplitude = 0.0;
        let state = make_init_state(&spec);
        save_field_png(&state, V2Field::SThickness, &out_dir.join("radial_s.png"))
            .expect("save radial s");
        save_field_png(&state, V2Field::Altitude, &out_dir.join("radial_altitude.png"))
            .expect("save radial altitude");
        save_altitude_fixed_scale(
            &state,
            &out_dir.join("radial_altitude_fixed.png"),
            s_min_global,
            s_max_global,
        )
        .expect("save radial altitude fixed");
        print_init_stats("radial", &state);
    }

    // Variant 3 — Orogenic-seul σ=0.10 (Run B comparator from R7.A.1).
    {
        let mut spec = base_spec.clone();
        spec.grid_nx = 64;
        spec.grid_ny = 64;
        spec.init_mode = V2InitModeSpec::orogenic_default();
        if let V2InitModeSpec::Orogenic { width_sigma_ratio, .. } = &mut spec.init_mode {
            *width_sigma_ratio = 0.10;
        }
        spec.s_perturbation_amplitude = 0.0;
        let state = make_init_state(&spec);
        save_field_png(&state, V2Field::SThickness, &out_dir.join("orogenic_s.png"))
            .expect("save orogenic s");
        save_field_png(
            &state,
            V2Field::Altitude,
            &out_dir.join("orogenic_altitude.png"),
        )
        .expect("save orogenic altitude");
        save_altitude_fixed_scale(
            &state,
            &out_dir.join("orogenic_altitude_fixed.png"),
            s_min_global,
            s_max_global,
        )
        .expect("save orogenic altitude fixed");
        print_init_stats("orogenic σ=0.10", &state);
    }

    // Variant 4 — Composite (Run C candidate, R7.A.2 defaults).
    {
        let mut spec = base_spec;
        spec.grid_nx = 64;
        spec.grid_ny = 64;
        spec.init_mode = V2InitModeSpec::composite_default();
        spec.s_perturbation_amplitude = 0.0;
        let state = make_init_state(&spec);
        save_field_png(&state, V2Field::SThickness, &out_dir.join("composite_s.png"))
            .expect("save composite s");
        save_field_png(
            &state,
            V2Field::Altitude,
            &out_dir.join("composite_altitude.png"),
        )
        .expect("save composite altitude");
        save_altitude_fixed_scale(
            &state,
            &out_dir.join("composite_altitude_fixed.png"),
            s_min_global,
            s_max_global,
        )
        .expect("save composite altitude fixed");
        print_init_stats("composite", &state);
        let max_s: f64 = state
            .s_field
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        // Sanity probe: the composite max must exceed the orogenic-seul
        // max (it strictly adds dome to the ridge before capping). At
        // 64² σ=0.10 orogenic-seul reaches ~1.17; composite must clear
        // 1.18 at minimum (and likely hits ~1.20 thanks to the dome).
        assert!(
            max_s > 1.10,
            "composite max(S̃) = {max_s} ≤ 1.10 — dome + ridge wiring \
             check needed",
        );
    }

    println!(
        "[r7.a.2] init preview done. Inspect {}",
        out_dir.display(),
    );
    println!(
        "[r7.a.2] Compare *_altitude_fixed.png across the 4 variants — \
         shared palette [{s_min_global}, {s_max_global}]. The *_s.png \
         files are also globally comparable (fixed [0, 2.5] palette in \
         field_to_rgba). The *_altitude.png files (no `_fixed` suffix) \
         use the buggy per-image isostasy normalisation; do NOT compare \
         them across runs.",
    );
}

// ----------------------------------------------------------------------------
// Step 12 R7.A.2.4 — simulation A/B/C at 64² × 5 cycles × 20 steps under
// workflow ON × mf=1.0 × evo=0.10 × craton_amp=3 (the EVO.C sweet spot
// empirically located by R6.3 Option β + C.4). The init mode is the
// only variable across runs:
//
//   Run A — InitMode::RadialProfile (Step 13 baseline)
//   Run B — InitMode::Orogenic σ=0.10 (R7.A.1)
//   Run C — InitMode::Composite defaults (R7.A.2)
//
// Acceptance R4.1–R4.6 per cycle (the 6 multi-dim axes the user has
// repeatedly insisted on) + visual R4.3 (chains formed, ridges
// persistent, plate boundaries deformed). The verdict
// PASS / MARGINAL / FAIL / FAIL.deep is reported in
// `docs/reports/step12_r7_a_2_4_simulation/summary.md`.

/// Runtime palette bounds for the comparable altitude renders
/// (R7.A.2.3 finding — must be shared across all runs).
const R7_ALT_S_MIN: f64 = 0.0;
const R7_ALT_S_MAX: f64 = 1.5;

/// Cycles at which to capture PNGs + per-cycle metrics. Cycle 0 is
/// the init state; the loop's per-cycle output is `cycle ∈ 1..=5`.
const R7_CAPTURE_CYCLES: &[usize] = &[1, 3, 5];

fn r7_a_2_4_apply_craton_amp(spec: &mut ymir_viz::bridge::v2::V2RunSpec, amp: f64) {
    if let ymir_viz::bridge::v2::V2CratonicSpec::On {
        cr,
        k_viscous,
        b_factor,
        smoothing_width,
        plate_area_min,
    } = spec.cratonic
    {
        spec.cratonic = ymir_viz::bridge::v2::V2CratonicSpec::On {
            cr,
            k_viscous: k_viscous * amp,
            b_factor: b_factor * amp,
            smoothing_width,
            plate_area_min,
        };
    }
}

/// Run one of the three R7.A.2.4 configurations. Returns a summary
/// suitable for the consolidated comparison table.
fn run_r7_a_2_4_config(
    init_mode: ymir_viz::bridge::v2::V2InitModeSpec,
    label: &str,
) -> R6ConfigSummary {
    use ymir_viz::bridge::v2::{build_config, V2MantleSpec};

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r7_a_2_4_simulation")
        .join(label);
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!("\n=== R7.A.2.4: {label} ===");

    let mut spec = presets::load("active_medley").expect("load active_medley");
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    if let V2MantleSpec::On { coupling, num_modes, seed, .. } = spec.mantle {
        spec.mantle = V2MantleSpec::On {
            mf: 1.0,
            coupling,
            num_modes,
            seed,
            evolution_rate: 0.10,
        };
    }
    r7_a_2_4_apply_craton_amp(&mut spec, 3.0);
    spec.init_mode = init_mode;
    spec.s_perturbation_amplitude = 0.0;
    spec.grid_nx = 64;
    spec.grid_ny = 64;
    spec.steps = 5 * 20;
    spec.total_time_nondim = 6.0;

    // Init capture (cycle 0).
    let init_state = make_init_state(&spec);
    save_field_png(&init_state, V2Field::SThickness, &out_dir.join("cycle_0_s.png"))
        .expect("save init s");
    save_altitude_fixed_scale(
        &init_state,
        &out_dir.join("cycle_0_altitude_fixed.png"),
        R7_ALT_S_MIN,
        R7_ALT_S_MAX,
    )
    .expect("save init altitude_fixed");
    print_init_stats(&format!("{label}/init"), &init_state);

    let mut cfg = build_config::build(&spec);
    let wf = build_config::build_workflow(&spec.workflow);

    let t0 = std::time::Instant::now();
    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r7.a.2.4/{label}] simulation completed in {elapsed:.1}s");

    let init_field =
        Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.common.sea_level_normalized).unwrap_or(0.5);
    let init_metrics = compute_metrics(0, &init_field, init_sea, 0.0, 0.0);
    let mut metrics = vec![init_metrics];

    let mut peak_v_per_cycle: Vec<f64> = Vec::with_capacity(5);
    let mut newton_outcomes_per_cycle: Vec<(usize, usize, usize, usize)> =
        Vec::with_capacity(5);
    let mut cg_mean_quartiles_per_cycle: Vec<(f64, f64, f64)> = Vec::with_capacity(5);

    for (i, cycle) in output.cycles.iter().enumerate() {
        let cyc_idx = i + 1;
        if R7_CAPTURE_CYCLES.contains(&cyc_idx) {
            let v2_state = V2FinalState::from_harness(&cycle.baseline.final_state);
            save_field_png(
                &v2_state,
                V2Field::SThickness,
                &out_dir.join(format!("cycle_{cyc_idx}_s.png")),
            )
            .expect("save s");
            save_altitude_fixed_scale(
                &v2_state,
                &out_dir.join(format!("cycle_{cyc_idx}_altitude_fixed.png")),
                R7_ALT_S_MIN,
                R7_ALT_S_MAX,
            )
            .expect("save altitude_fixed");
        }
        metrics.push(compute_metrics(
            cyc_idx,
            &cycle.baseline.final_state.s_field,
            cycle.common.sea_level_normalized,
            cycle.common.erosion_volume_removed,
            cycle.common.mass_drift,
        ));
        peak_v_per_cycle.push(cycle.baseline.metrics.vmax_peak);
        if let Some(n) = cycle.baseline.metrics.newton.as_ref() {
            newton_outcomes_per_cycle.push((n.converged, n.stalled, n.diverged, n.capped));
            let mut cg_per_step: Vec<usize> = Vec::with_capacity(n.outer_iters.len());
            let mut cg_cursor: usize = 0;
            for &outer in &n.outer_iters {
                let outer_n = outer as usize;
                let cg_slice = &n.cg_iters_per_newton_step
                    [cg_cursor..(cg_cursor + outer_n).min(n.cg_iters_per_newton_step.len())];
                cg_per_step.push(cg_slice.iter().sum());
                cg_cursor += outer_n;
            }
            let k = cg_per_step.len();
            if k >= 4 {
                let early = &cg_per_step[..k / 4];
                let mid = &cg_per_step[k / 4..3 * k / 4];
                let late = &cg_per_step[3 * k / 4..];
                let avg = |v: &[usize]| -> f64 {
                    if v.is_empty() { 0.0 } else { v.iter().sum::<usize>() as f64 / v.len() as f64 }
                };
                cg_mean_quartiles_per_cycle.push((avg(early), avg(mid), avg(late)));
            } else {
                cg_mean_quartiles_per_cycle.push((0.0, 0.0, 0.0));
            }
        } else {
            newton_outcomes_per_cycle.push((0, 0, 0, 0));
            cg_mean_quartiles_per_cycle.push((0.0, 0.0, 0.0));
        }
    }

    // Per-config metrics.md with the per-cycle R4.1–R4.6 table.
    let mut md = String::new();
    md.push_str(&format!("# R7.A.2.4 — {label}\n\n"));
    md.push_str(&format!(
        "64² active_medley, workflow ON (D2 + D1-ter), mf=1.0, evo=0.10, \
         craton_amp=3. 5 cycles × 20 steps. Runtime: {:.1}s.\n\n",
        elapsed
    ));
    md.push_str("## Per-cycle solver + S̃ health\n\n");
    md.push_str(
        "| cycle | peak \\|v\\| | Newton C/S/D/Cap | peak S̃ | frac>0.85 | frac>0.95 | mass | max_path | CG e/m/l |\n",
    );
    md.push_str("|---|---|---|---|---|---|---|---|---|\n");
    for (i, _cycle) in output.cycles.iter().enumerate() {
        let m = &metrics[i + 1];
        let o = newton_outcomes_per_cycle[i];
        let q = cg_mean_quartiles_per_cycle[i];
        let frac_above_0_85 = m.mass_continental / (m.mass_continental.abs() + m.mass_oceanic.abs()).max(1e-12);
        md.push_str(&format!(
            "| {} | {:.3e} | {}/{}/{}/{} | {:.3} | {:.3} | {:.3} | {:.2} | {} | {:.0}/{:.0}/{:.0} |\n",
            i + 1,
            peak_v_per_cycle[i],
            o.0, o.1, o.2, o.3,
            m.peak_s,
            frac_above_0_85,
            m.fraction_above_0_8,
            m.mass_total,
            m.max_path_length,
            q.0, q.1, q.2,
        ));
    }
    md.push('\n');

    // 6-criteria verdict per run.
    let init_frac = metrics[0].fraction_above_0_8;
    let final_frac = metrics.last().unwrap().fraction_above_0_8;
    let retention = if init_frac > 1e-9 { final_frac / init_frac } else { 0.0 };
    let n_dynamic = peak_v_per_cycle.iter().filter(|&&v| v > 0.1).count();
    let mass_init = metrics[0].mass_total;
    let mass_final = metrics.last().unwrap().mass_total;
    let mass_loss_pct = (mass_init - mass_final) / mass_init.abs().max(1e-12) * 100.0;
    let mass_loss_per_cycle_pct = mass_loss_pct / 5.0;
    let final_peak_s = metrics.last().unwrap().peak_s;
    let final_sea = metrics.last().unwrap().sea_level_ref;
    let max_max_path = metrics
        .iter()
        .skip(1)
        .map(|m| m.max_path_length)
        .max()
        .unwrap_or(0);

    let r4_1 = final_peak_s > final_sea;
    let r4_2 = retention > 0.5;
    let r4_4 = mass_loss_per_cycle_pct.abs() < 1.0;
    let r4_5 = max_max_path >= 5;
    let r4_6 = n_dynamic >= 3;
    let auto_count = [r4_1, r4_2, r4_4, r4_5, r4_6].iter().filter(|&&p| p).count();

    md.push_str("## Multi-dim acceptance (R4.1–R4.6)\n\n");
    md.push_str(&format!(
        "- R4.1 Continents émergés: peak S̃_final = {:.3} > sea = {:.3} → **{}**\n",
        final_peak_s, final_sea,
        if r4_1 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "- R4.2 Cratons préservés: retention = {:.1} % → **{}**\n",
        100.0 * retention,
        if r4_2 { "PASS" } else { "FAIL" }
    ));
    md.push_str("- R4.3 Bordures + chaînes: VISUAL (inspect `cycle_5_altitude_fixed.png`)\n");
    md.push_str(&format!(
        "- R4.4 Conservation: mass loss/cycle = {:.3} % → **{}**\n",
        mass_loss_per_cycle_pct,
        if r4_4 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "- R4.5 Drainage actif: max_path = {} (cycles 1-5) → **{}**\n",
        max_max_path,
        if r4_5 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "- R4.6 Dynamique soutenue: peak |v| > 0.1 on {}/5 → **{}**\n",
        n_dynamic,
        if r4_6 { "PASS" } else { "FAIL" }
    ));
    md.push_str(&format!(
        "\nAuto count: **{} / 5** (R4.3 visual pending).\n",
        auto_count
    ));

    std::fs::write(out_dir.join("metrics.md"), &md).expect("write metrics");
    println!("{md}");

    R6ConfigSummary {
        mf: 1.0,
        evo_rate: 0.10,
        label: label.to_string(),
        runtime_s: elapsed,
        peak_v_per_cycle,
        peak_s_per_cycle: metrics.iter().skip(1).map(|m| m.peak_s).collect(),
        frac_above_0_8_per_cycle: metrics.iter().skip(1).map(|m| m.fraction_above_0_8).collect(),
        mass_total_per_cycle: metrics.iter().skip(1).map(|m| m.mass_total).collect(),
        max_path_per_cycle: metrics.iter().skip(1).map(|m| m.max_path_length).collect(),
        sea_level_per_cycle: metrics.iter().skip(1).map(|m| m.sea_level_ref).collect(),
        newton_outcomes_per_cycle,
        cg_mean_quartiles_per_cycle,
        init_frac_above_0_8: metrics[0].fraction_above_0_8,
        init_peak_s: metrics[0].peak_s,
    }
}

/// R7.A.2.4 master — runs A, B, C sequentially and writes
/// `summary.md` with the comparative R4.1–R4.6 table + EVO verdict.
///
/// ```bash
/// cargo test --release -p ymir-viz --test v2_workflow_r4_visual_checkpoint \
///     r7_a_2_4_simulation_abc -- --ignored --nocapture --test-threads=1
/// ```
#[test]
#[ignore]
fn r7_a_2_4_simulation_abc() {
    use ymir_viz::bridge::v2::V2InitModeSpec;

    let out_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r7_a_2_4_simulation");
    std::fs::create_dir_all(&out_root).expect("create out root");

    let configs: Vec<(V2InitModeSpec, &str)> = vec![
        (V2InitModeSpec::radial_profile_default(), "run_a_radial"),
        ({
            let mut m = V2InitModeSpec::orogenic_default();
            if let V2InitModeSpec::Orogenic { width_sigma_ratio, .. } = &mut m {
                *width_sigma_ratio = 0.10;
            }
            m
        }, "run_b_orogenic_sigma_10"),
        (V2InitModeSpec::composite_default(), "run_c_composite"),
    ];

    let t_total = std::time::Instant::now();
    let mut summaries: Vec<R6ConfigSummary> = Vec::with_capacity(configs.len());
    for (idx, (init_mode, label)) in configs.into_iter().enumerate() {
        println!(
            "\n[r7.a.2.4] ===== {}/3 : {label} =====",
            idx + 1,
        );
        let s = run_r7_a_2_4_config(init_mode, label);
        summaries.push(s);
        write_r6_3_summary(&summaries, &out_root); // partial-progress snapshot
        println!(
            "[r7.a.2.4] cumulative elapsed: {:.1} min",
            t_total.elapsed().as_secs_f64() / 60.0,
        );
    }
    println!(
        "\n[r7.a.2.4] simulation done. {} runs in {:.1} min total. Inspect {}",
        summaries.len(),
        t_total.elapsed().as_secs_f64() / 60.0,
        out_root.display(),
    );
}
