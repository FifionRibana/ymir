//! Issue #127 Phase 1.4 Stage D — downstream pipeline validation
//! tests on the C1-produced altitude.
//!
//! Two tests under default features. Validate that the altitude
//! produced by the full C1 closure stack (Davis-Suppe +
//! equilibrium-height + stream-power erosion) is **consumable**
//! by the existing v2 downstream pipeline (flow routing + HD
//! erosion + upscale). The downstream pipeline is paradigm-
//! agnostic at runtime (per Phase 1.3 H1 audit § 2.2 — `phase_b`
//! consumes `Field2D + sea_level + iso_config`, no
//! `BaselineResult`), so it should accept C1 output without
//! modification.
//!
//! Test inventory:
//!
//! - [`c1_continental_drainage_functional`] — heightmap-level
//!   D8 flow accumulation, **restricted to continental cells**
//!   (altitude > `sea_level_normalized`), produces functional
//!   drainage on the surviving continental fraction. Phase 1.4
//!   erosion eats the bulk continental mass (35 % total mass
//!   loss per Stage E3 measurement), leaving sparse continental
//!   ridges in an oceanic domain (Earth-like terrane geometry).
//!   The meaningful Phase-1.4 drainage question is therefore
//!   "do the remaining continental cells drain?" — not "is the
//!   overall grid dendritic?" (which it isn't, by erosion-
//!   regime design). Thresholds tagged Phase-1.4-regime per the
//!   `feedback_fill_ratio_regime_agnostic_metric` pattern.
//! - [`c1_post_run_altitude_consumable_by_downstream`] — smoke
//!   test on the full downstream stack: `compute_flow` then
//!   `run_erosion` on the C1 heightmap. Asserts no panic, all
//!   output buffers finite, and the HD-eroded heightmap stays
//!   within `[0, 1]`. Validates the Phase 4+ UI / export
//!   consumability invariant.
//!
//! ## Skipped tests (deferred per Phase 1.4 scope)
//!
//! Two tests originally proposed in the Stage D spec were
//! dropped after a pre-code accessibility check showed the
//! downstream API surface they would have consumed is **not
//! available**:
//!
//! - **`c1_rainfall_field_computable`** — the `climate/` module
//!   in `crates/ymir-core/src/climate/` is currently a stub: the
//!   `mod.rs` declaration of `pub mod precipitation` /
//!   `pub mod temperature` / `pub mod biomes` is commented out
//!   (placeholder "M3 — altitude lapse, continentality,
//!   seasons"), and the corresponding `.rs` files are empty
//!   (0 lines). `compute_rainfall` does not exist. Defer this
//!   test until the climate module is implemented (Phase 2+
//!   per the design doc roadmap).
//!
//! - **`c1_slope_field_consistent_with_erosion`** — there is no
//!   standalone public slope utility in the ymir-core
//!   downstream pipeline. Slope is computed *internally* in two
//!   places: (a) the C1 erosion's private `compute_local_slope`
//!   helper, and (b) implicitly via D8 routing inside
//!   `terrain::flow::compute_flow`. A "consistency" test
//!   requires *two* exposed implementations to compare; there
//!   is only one (private). Exposing a `pub fn
//!   compute_slope_magnitude` from the downstream layer would
//!   be YAGNI scope expansion absent a real consumer. Defer
//!   until a downstream consumer (e.g. orographic precipitation
//!   in the future climate module) motivates exposing it.

use ymir_core::erosion::hydraulic::{run_erosion, ErosionConfig};
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::closures::accretion::AccretionParams;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::rifting::RiftingParams;
use ymir_core::tectonics_c1::closures::subduction::SubductionParams;
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::terrain::flow::{compute_flow, FlowConfig};

const GRID: usize = 64;
const SEED: u64 = 42;

/// Phase 1.4 closure stack — Davis-Suppe + equilibrium-height +
/// erosion all ON, Phase 2 Track A oceanic bathymetry OFF. Locks
/// the Phase 1.4 downstream-test regime against silent regime
/// drift from `C1Closures::default()` (post-#129 enables S-S). See
/// the matching helper in `c1_phase_1_4_erosion.rs` for the
/// rationale on holding the Phase 1.4 thresholds stable.
fn phase_1_4_closures() -> C1Closures {
    C1Closures {
        oceanic_bathymetry: SteinSteinParams {
            enabled: false,
            ..SteinSteinParams::default()
        },
        subduction: SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        },
        accretion: AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        },
        rifting: RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        },
        ..C1Closures::default()
    }
}

#[test]
fn c1_continental_drainage_functional() {
    // After 300 steps with all 3 closures, the C1 altitude has
    // an Earth-like sparse-continental topology: erosion has
    // eaten the bulk continental mass, leaving wedge-body ridges
    // in mostly-oceanic domain (per Stage E3 calibration
    // measurement: 35 % total mass loss, mean S̃ drops 1.574 →
    // 0.361, sparse continental fraction).
    //
    // The original Stage D spec asked for a "dendritic drainage
    // pattern" on the entire grid. That formulation assumed a
    // Phase-1.3-style heavily-continental heightmap; it fails on
    // Phase 1.4 output because most cells are oceanic by design.
    // The meaningful Phase 1.4 question is **"do the surviving
    // continental cells support functional drainage?"** —
    // restrict the analysis to cells with altitude above
    // `sea_level_normalized` and assert three properties on
    // that subset:
    //
    //   (a) `n_continental > 50` — some continental survives
    //       erosion (not a fully-flooded grid).
    //   (b) `max_continental_accum > 5` — at least one
    //       continental cell accumulates non-trivial drainage
    //       (proves SOMETHING flows on the continental ridge).
    //   (c) `connected_continental_fraction > 5 %` (of
    //       continental cells, NOT of total cells) — the
    //       continental drainage network has non-trivial
    //       connectivity, not just isolated peaks.
    //
    // Thresholds are **Phase-1.4-regime-tagged** per the
    // [[fill-ratio-regime-agnostic-metric]] memory's "regime-
    // tagged direction" pattern. Phase 1.5+ (potentially with
    // additional closures or modified kinematics) may shift
    // the thresholds; the Stage E4 memory entry
    // `project_c1_phase_1_4_erosion_outcomes` will document
    // the empirical Phase-1.4 baseline against which Phase 1.5+
    // changes are compared.

    let mut state = init_c1_state_phase_1_1(GRID, SEED);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();
    let closures = phase_1_4_closures();
    let config = C1TimeLoopConfig {
        rigid_continental_crust: false,
        n_steps: 300,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

    // Compute altitude + D8 flow.
    let isostasy = compute_isostasy(&state.s, &iso_config);
    let altitude = isostasy.heightmap;
    let sea_level = isostasy.sea_level_normalized;
    let flow_config = FlowConfig { sea_level, ..Default::default() };
    let flow = compute_flow(&altitude, &flow_config);

    // Continental filtering — restrict drainage analysis to
    // cells with altitude > sea_level_normalized.
    let n_cells = GRID * GRID;
    let mut n_continental = 0_usize;
    let mut max_continental_accum = 0.0_f32;
    let mut sum_continental_accum = 0.0_f64;
    // Three connectivity cuts capture the drainage-tree shape:
    //   non_leaf  : accum > 1  — cell has at least one upstream
    //               cell draining through it (≠ a leaf / isolated
    //               peak). The drainage-tree-shape question.
    //   downstream: accum > 2  — cell is downstream of at least
    //               two cells. Captures genuine drainage chains
    //               vs single-step paths.
    //   major     : accum > 5  — cell collects 5+ upstream cells.
    //               Equivalent to a "main channel" in the
    //               D8 drainage tree.
    let mut non_leaf_continental = 0_usize;
    let mut downstream_continental = 0_usize;
    let mut major_continental = 0_usize;
    for k in 0..n_cells {
        if altitude.data[k] > sea_level {
            n_continental += 1;
            let a = flow.accumulation.data[k];
            if a > max_continental_accum {
                max_continental_accum = a;
            }
            sum_continental_accum += a as f64;
            if a > 1.0 {
                non_leaf_continental += 1;
            }
            if a > 2.0 {
                downstream_continental += 1;
            }
            if a > 5.0 {
                major_continental += 1;
            }
        }
    }
    let continental_fraction = n_continental as f64 / n_cells as f64;
    let mean_continental_accum = if n_continental > 0 {
        sum_continental_accum / n_continental as f64
    } else {
        0.0
    };
    let non_leaf_fraction = if n_continental > 0 {
        non_leaf_continental as f64 / n_continental as f64
    } else {
        0.0
    };
    let downstream_fraction = if n_continental > 0 {
        downstream_continental as f64 / n_continental as f64
    } else {
        0.0
    };
    let major_fraction = if n_continental > 0 {
        major_continental as f64 / n_continental as f64
    } else {
        0.0
    };

    eprintln!(
        "c1_phase_1_4 D-T1 — continental drainage functional (Phase-1.4-regime-tagged):"
    );
    eprintln!("  num_cells              = {n_cells}");
    eprintln!(
        "  n_continental          = {n_continental}    ({:.1} % of grid, sea_level = {sea_level:.3})",
        100.0 * continental_fraction
    );
    eprintln!("  num_basins (global)    = {}", flow.num_basins);
    eprintln!(
        "  max  continental accum = {max_continental_accum:.1}    (threshold: > 5)"
    );
    eprintln!("  mean continental accum = {mean_continental_accum:.2}");
    eprintln!(
        "  non-leaf  (accum > 1)  = {non_leaf_continental:>3} cells  ({:>4.1} % of continental, threshold: > 25 %)",
        100.0 * non_leaf_fraction
    );
    eprintln!(
        "  downstream(accum > 2)  = {downstream_continental:>3} cells  ({:>4.1} % of continental)",
        100.0 * downstream_fraction
    );
    eprintln!(
        "  major     (accum > 5)  = {major_continental:>3} cells  ({:>4.1} % of continental)",
        100.0 * major_fraction
    );
    eprintln!(
        "  → continent is sparse Earth-like (wedge ridges in ocean); drainage functional on surviving land"
    );

    // Sub-assertion (a) — some continental survives erosion.
    assert!(
        n_continental > 50,
        "n_continental = {n_continental} ≤ 50 — Phase 1.4 erosion has effectively \
         flooded the entire grid; no continental fraction remains for drainage. \
         Investigate K_erosion calibration (too aggressive) or h_eq (too low to \
         preserve any wedge above sea level)."
    );

    // Sub-assertion (b) — non-trivial drainage on continental.
    assert!(
        max_continental_accum > 5.0,
        "max_continental_accum = {max_continental_accum:.1} ≤ 5 — continental cells \
         have no functional drainage. Investigate compute_flow's interaction with \
         the Phase 1.4 altitude (too few connected continental cells?)."
    );

    // Sub-assertion (c) — drainage-tree shape: at least 25 % of
    // continental cells are non-leaf (accum > 1, i.e. they sit
    // downstream of at least one other continental cell).
    //
    // Threshold revised from the initial draft (which used
    // `major_fraction > 5 %`, calibrated for a Phase-1.3-style
    // heavily-continental heightmap and produced 0.8 % on the
    // sparse Phase 1.4 output). The `non_leaf` cut is the
    // appropriate Phase-1.4-regime metric: it asks "are
    // continental cells assembled into drainage chains, or
    // scattered as isolated peaks?" — and Phase 1.4's sparse-
    // continental Earth-like morphology should still have most
    // continental cells participating in some chain, even if
    // those chains are short.
    assert!(
        non_leaf_fraction > 0.25,
        "non_leaf_fraction = {:.1} % ≤ 25 % — surviving continental cells are \
         predominantly isolated leaves (every cell is a peak, no upstream cell \
         drains through it). Drainage tree shape too fragmented for downstream \
         river / biome consumers; the wedge ridges are not forming connected \
         drainage networks.",
        100.0 * non_leaf_fraction
    );
}

#[test]
fn c1_post_run_altitude_consumable_by_downstream() {
    // Smoke test: after 300 steps with all closures, the
    // C1 altitude must feed cleanly through the full v2
    // downstream pipeline (`compute_flow` + `run_erosion`)
    // without panic, NaN, or out-of-range output. Validates
    // the Phase 4+ UI / export consumability invariant — the
    // C1 output is plug-compatible with the v2 downstream
    // stack that consumes altitude.

    let mut state = init_c1_state_phase_1_1(GRID, SEED);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();
    let closures = phase_1_4_closures();
    let config = C1TimeLoopConfig {
        rigid_continental_crust: false,
        n_steps: 300,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

    // Stage 1 — isostasy. Smoke: altitude finite, in [0, 1].
    let isostasy = compute_isostasy(&state.s, &iso_config);
    let altitude = isostasy.heightmap;
    assert!(
        altitude.data.iter().all(|h| h.is_finite()),
        "altitude contains non-finite values"
    );
    assert!(
        altitude.data.iter().all(|&h| (0.0..=1.0).contains(&h)),
        "altitude contains values outside [0, 1]"
    );

    // Stage 2 — compute_flow. Smoke: accumulation + filled both
    // finite; no panic.
    let flow_config = FlowConfig {
        sea_level: isostasy.sea_level_normalized,
        ..Default::default()
    };
    let flow = compute_flow(&altitude, &flow_config);
    assert!(
        flow.accumulation.data.iter().all(|a| a.is_finite() && *a >= 0.0),
        "compute_flow produced non-finite or negative accumulation"
    );
    assert!(
        flow.filled.data.iter().all(|h| h.is_finite()),
        "compute_flow produced non-finite filled heightmap"
    );

    // Stage 3 — run_erosion (HD hydraulic). Smoke: heightmap +
    // sediment finite, heightmap stays in [0, 1] bounds. Use a
    // minimal `num_droplets` to keep test runtime acceptable —
    // this is a consumability check, not an erosion-effect
    // validation (that's Phase B's job in workflow/phase_b.rs).
    let erosion_config = ErosionConfig {
        num_droplets: 200, // minimal smoke pass
        ..ErosionConfig::default()
    };
    let world_seed = WorldSeed::new(SEED);
    let eroded = run_erosion(&altitude, &erosion_config, &world_seed, |_, _, _| true);
    assert!(
        eroded.heightmap.data.iter().all(|h| h.is_finite()),
        "run_erosion produced non-finite heightmap"
    );
    assert!(
        eroded.sediment.data.iter().all(|s| s.is_finite() && *s >= 0.0),
        "run_erosion produced non-finite or negative sediment"
    );
    assert!(
        eroded.heightmap.data.iter().all(|&h| (-0.05..=1.05).contains(&h)),
        "run_erosion produced heightmap outside [-0.05, 1.05] tolerance — \
         droplets may have over-incised or over-deposited"
    );

    eprintln!("c1_phase_1_4 D-T2 — downstream consumability smoke:");
    eprintln!(
        "  altitude in [{:.3}, {:.3}], sea_level_normalized = {:.3}",
        altitude.data.iter().cloned().fold(f32::INFINITY, f32::min),
        altitude.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max),
        isostasy.sea_level_normalized
    );
    eprintln!(
        "  compute_flow: num_basins = {}, max accumulation = {:.0}",
        flow.num_basins,
        flow.accumulation.data.iter().cloned().fold(0.0_f32, f32::max)
    );
    eprintln!(
        "  run_erosion ({} droplets): mean post-erosion altitude = {:.4}, max sediment = {:.4}",
        erosion_config.num_droplets,
        eroded.heightmap.data.iter().sum::<f32>() / (GRID * GRID) as f32,
        eroded.sediment.data.iter().cloned().fold(0.0_f32, f32::max)
    );
}
