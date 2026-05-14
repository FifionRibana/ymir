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
    drainage::compute_drainage_targets, run_phase_a_loop,
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
    let output = run_phase_a_loop(&mut cfg, &wf);
    let elapsed = t0.elapsed().as_secs_f64();
    println!("[r4] {preset_name}: completed in {elapsed:.1}s, {} cycles", output.cycles.len());

    // INIT metrics (use init field, no erosion stats).
    let init_field = Field2D::from_vec(init_state.nx, init_state.ny, init_state.s_field.clone());
    let init_sea = output.cycles.first().map(|c| c.sea_level_normalized).unwrap_or(0.5);
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
            cycle.sea_level_normalized,
            cycle.erosion_volume_removed,
            cycle.mass_drift,
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
    let output = run_phase_a_loop(&mut cfg, &wf);
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
    let init_sea = output.cycles.first().map(|c| c.sea_level_normalized).unwrap_or(0.5);
    let init_metrics = compute_metrics(0, &init_field, init_sea, 0.0, 0.0);
    let final_metrics = compute_metrics(
        R4B_N_CYCLES,
        &final_cycle.baseline.final_state.s_field,
        final_cycle.sea_level_normalized,
        final_cycle.erosion_volume_removed,
        final_cycle.mass_drift,
    );
    let cumulative_drift: f64 = output.cycles.iter().map(|c| c.mass_drift.abs()).sum();

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
    let output = run_phase_a_loop(&mut cfg, &wf);
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
    let init_sea = output.cycles.first().map(|c| c.sea_level_normalized).unwrap_or(0.5);
    let init_metrics = compute_metrics(0, &init_field, init_sea, 0.0, 0.0);
    let final_metrics = compute_metrics(
        5,
        &final_cycle.baseline.final_state.s_field,
        final_cycle.sea_level_normalized,
        final_cycle.erosion_volume_removed,
        final_cycle.mass_drift,
    );
    let cumulative_drift: f64 = output.cycles.iter().map(|c| c.mass_drift.abs()).sum();

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
    let output = run_phase_a_loop(&mut cfg, &wf);
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
    md.push_str(&format!("- Mass drift over 1 cycle: {:.3e}\n", cycle.mass_drift));
    md.push_str(&format!("- Mass drift cumulative (1 cycle = same): {:.3e}\n\n", cycle.mass_drift.abs()));

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
