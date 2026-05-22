//! Phase A — low-res loop orchestration.
//!
//! `run_phase_a_cycle` chains the 5-step single-cycle pipeline:
//!
//! 1. **Tectonic** — `run_baseline(cfg)`; the cfg may carry a
//!    [`ContinuationState`] for cycle-to-cycle warm-start (D3).
//! 2. **Isostasy** — `compute_isostasy(s_field)` to extract
//!    `sea_level_normalized` (Option 2 of Phase 0 finding E:
//!    adaptive threshold drives erosion + reclassification).
//! 3. **Macro mass redistribution** — `macro_redistribution::apply`
//!    in-place on `final_state.s_field` with the cycle's `α /
//!    isostatic_rebound_ratio / max_drainage_distance`. Step 12 R3
//!    replaced the legacy `low_res_erosion` (per-cycle diffusive
//!    erosion + local `β`-deposition) with a long-distance drainage
//!    + isostatic-rebound formulation — see
//!    [`super::macro_redistribution`] module docstring for the full
//!    arithmetic and conservation contract.
//! 4. **Reclassify** — per-cell `plate_type[i, j]` = `Continental`
//!    iff `s_field[i, j] > sea_level_normalized`; otherwise
//!    `Oceanic`.
//! 5. **Recompute cratonic factor** — clone the (static) `plate_id`
//!    field; rebuild `per_plate_type[p]` from the post-erosion
//!    continental fraction (D4: a plate retains its continental
//!    eligibility iff `frac >= plate_area_min`); call
//!    [`crate::tectonics_v2::cratonic::factor::build_cratonic_factor_field`]
//!    with this synthesised `VoronoiPlates` (Step 9 BFS, Manhattan
//!    4-conn, smoothstep — *not* `voronoi/distance.rs` which is
//!    Chebyshev 8-conn). The new factor field replaces
//!    `final_state.cratonic_factor`.
//!
//! Order is **strict**: erosion before reclassify, reclassify before
//! craton recompute. Otherwise the craton would reflect the pre-erosion
//! `S̃` and the cycle's effect would be invisible to Phase 4's
//! `craton_recomputation_change` metric.
//!
//! `WorkflowConfig::Disabled` short-circuits this entire pipeline:
//! the cycle is exactly `run_baseline(cfg)` with all extra scalars at
//! zero/`None`. The regression
//! `v2_workflow_disabled_regression::workflow_disabled_run_phase_a_cycle_is_bit_identical_to_run_baseline`
//! pins this contract byte-for-byte.

use super::{macro_redistribution, CycleOutput, CycleOutputCommon, PhaseAOutput, WorkflowConfig};
use crate::tectonics::isostasy::IsostasyConfig;
use crate::tectonics_v2::boundaries::{PlateType, PlateTypeField};
use crate::tectonics_v2::cratonic::factor::build_cratonic_factor_field;
use crate::tectonics_v2::cratonic::{CratonicConfig, CratonicConfigEnabled};
use crate::tectonics_v2::diagnostics::harness::{
    run_baseline_with_progress, BaselineConfig, ContinuationState, FinalState, StepProgress,
};
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::voronoi::{PlateIdField, VoronoiPlates};

/// Run a single Phase A cycle. Thin wrapper over
/// [`run_phase_a_cycle_with_progress`] with a no-op callback that
/// never aborts, preserving the bit-identical regression contract
/// (acceptance #15) byte-for-byte: `run_baseline_with_progress(cfg,
/// |_| true)` is itself a wrapper over `run_baseline` from Step 8.6
/// follow-up, so the call chain reduces to the same primitive.
///
/// `Disabled` → direct `run_baseline(cfg)` passthrough wrapped in a
/// [`CycleOutput`] with all extra scalars at zero/`None`.
///
/// `Enabled(params)` → 5-step pipeline (tectonic → isostasy → erosion →
/// reclassify → recompute craton). Returns the post-cycle state
/// suitable for cycle-to-cycle continuation via
/// [`final_state_to_continuation`].
pub fn run_phase_a_cycle(cfg: &BaselineConfig, wf: &WorkflowConfig) -> CycleOutput {
    run_phase_a_cycle_with_progress(cfg, wf, |_| true)
}

/// Streaming variant of [`run_phase_a_cycle`]. The callback fires
/// once per completed harness step inside the cycle's tectonic
/// sub-phase (steps 1 of the 5-step pipeline); returning `false`
/// requests a graceful abort of the harness step loop. Same
/// callback shape as
/// [`crate::tectonics_v2::diagnostics::harness::run_baseline_with_progress`].
///
/// Added in Step 12 follow-up so the v2 bridge can stream per-step
/// `V2Event::Progress` to the metrics dashboard during Phase A
/// (the dashboard previously froze between `WorkflowCycleCompleted`
/// events because `run_phase_a_cycle` invoked `run_baseline` —
/// the `|_| true` callback wrapper — with no streaming hook). The
/// post-tectonic substeps (isostasy, erosion, reclassify, craton
/// recompute) are not currently streamed; they're sub-second on
/// 64² mantle-on, so a single "cycle progress" tick is the
/// pragmatic granularity.
pub fn run_phase_a_cycle_with_progress<F>(
    cfg: &BaselineConfig,
    wf: &WorkflowConfig,
    on_progress: F,
) -> CycleOutput
where
    F: FnMut(&StepProgress<'_>) -> bool,
{
    match wf {
        WorkflowConfig::Disabled => {
            let baseline = run_baseline_with_progress(cfg, on_progress);
            CycleOutput {
                baseline,
                common: CycleOutputCommon::default(),
            }
        }
        WorkflowConfig::Enabled(params) => {
            // Step 1: Tectonic.
            let mut baseline = run_baseline_with_progress(cfg, on_progress);

            // Step 2: Adaptive sea-level threshold in S̃ space.
            //
            // Phase 3 originally piped
            // `compute_isostasy(s).sea_level_normalized` (a `f32` in
            // heightmap `[0, 1]` space), which on the default
            // IsostasyConfig (max_depth=500 m, max_elevation=4000 m)
            // resolves to `0.111` — well below S̃'s natural oceanic
            // floor of ≈ 0.2. As a result the continental/oceanic
            // mask was satisfied by *every* cell, the recompute D4
            // rule never fired, and the v2_workflow_cratonic_recompute_*
            // tests passed for the wrong reason.
            //
            // Phase 3.5 fix: compute the threshold *in S̃ units* via
            // the isostasy formula `h_sea = h_min + sea_level_fraction
            // · h_range`, applied to S̃ directly:
            //
            //     s_sea = s_min + sea_level_fraction · (s_max - s_min)
            //
            // At init (S̃ ∈ [≈0.2, ≈1.2]) this resolves to ≈ 0.6 —
            // close to the natural 0.5 midpoint between oceanic and
            // continental, and adaptive to the S̃ distribution as
            // cycles erode. Source: same formula as
            // `compute_isostasy::h_sea`, just in S̃ units rather
            // than buoyancy-scaled altitude.
            let iso_cfg = IsostasyConfig::default();
            let sea_level_ref = {
                let s_data = baseline.final_state.s_field.data();
                let (s_min, s_max) = s_data.iter().copied().fold(
                    (f64::INFINITY, f64::NEG_INFINITY),
                    |(lo, hi), v| (lo.min(v), hi.max(v)),
                );
                let s_range = (s_max - s_min).max(1e-10);
                s_min + (iso_cfg.sea_level_fraction as f64) * s_range
            };

            // Step 3: Capture original `per_plate_type` BEFORE
            // reclassification — needed for the craton recompute's
            // "plate was originally continental" gate (D4).
            let original_per_plate_type = baseline
                .final_state
                .plate_type
                .as_ref()
                .zip(baseline.final_state.plate_id.as_ref())
                .map(|(pt, pid)| extract_per_plate_type(pid, pt));

            // Step 4: Macro mass redistribution (in-place). Step 12 R3
            // — replaces the legacy `low_res_erosion::apply`. The new
            // call carries drainage + isostatic rebound, so mass drift
            // is now ~ IEEE-754 floor by construction (instead of the
            // legacy `-(1-β) · volume_removed` net change).
            let mass_before: f64 = baseline.final_state.s_field.data().iter().sum();
            let stats = macro_redistribution::apply(
                &mut baseline.final_state.s_field,
                &params.phase_a,
                sea_level_ref,
            );
            let mass_after: f64 = baseline.final_state.s_field.data().iter().sum();

            // Step 4b: R5b D1-ter — reset velocity to zero after
            // macro_redistribution.
            //
            // EMPIRICAL FINDING (counter-intuitive vs classical Stokes
            // wisdom): the warm-start `v = v_final_previous_cycle` is
            // not just sub-optimal post-macro, it is **actively
            // harmful**. `macro_redistribution::apply` shifts S̃ enough
            // that the GPE driver direction changes; `v_warm_start`
            // points in a direction now anti-useful for the next
            // tectonic step, and Newton oscillates trying to correct
            // it (sub-case C amplified in D2-bis classification).
            //
            // The 3-variant D1-ter benchmark (commit 4969de9) showed
            // `v = 0` gives:
            //   - cycle 2 Converged 45/45 vs 14/41 with warm-start
            //   - 0 Oscillating vs 26 with warm-start
            //   - CG iter total over first 5 cycle-2 steps: 28k vs 75k
            //   - ‖Δv‖/‖v‖ max: 1.00 vs 11.43
            // Variant C (Gaussian smoothing of v) is WORSE than
            // warm-start (90k CG iter, peak |v| explosion to 53) —
            // smoothing preserves the wrong direction.
            //
            // Gated by `WorkflowConfig::Enabled` (this match arm).
            // The Disabled branch in `Self::Disabled => …` is
            // untouched and the bit-identical regression
            // `v2_workflow_disabled_regression` continues to hold.
            //
            // Counter-intuitive — leave the comment block intact; a
            // future dev tempted to "re-enable warm-start because
            // it's faster on Stokes" would re-break the system. See
            // `docs/reports/step12_solver_audit.md` § F and
            // `docs/reports/step12_r5b_d1_ter_init_variants/` for the
            // full empirical record.
            for v in baseline.final_state.vx.iter_mut() {
                *v = 0.0;
            }
            for v in baseline.final_state.vy.iter_mut() {
                *v = 0.0;
            }

            // Step 5: Reclassify per-cell `plate_type` from new
            // `s_field` against `sea_level_ref`. The Voronoï
            // tessellation (`plate_id`) is static for the run; only
            // the per-cell type changes.
            if let Some(plate_type) = baseline.final_state.plate_type.as_mut() {
                let s = &baseline.final_state.s_field;
                for j in 0..plate_type.ny() {
                    for i in 0..plate_type.nx() {
                        let new_type = if s.get(i, j) > sea_level_ref {
                            PlateType::Continental
                        } else {
                            PlateType::Oceanic
                        };
                        plate_type.set(i, j, new_type);
                    }
                }
            }

            // Step 6: Recompute cratonic factor — only if the cfg has
            // CratonicConfig::Enabled. We synthesise an updated
            // `VoronoiPlates` whose `per_plate_type[p]` reflects the
            // post-erosion D4 retention rule, then call
            // `build_cratonic_factor_field` (Step 9 algorithm).
            let mut craton_change: Option<f64> = None;
            if let CratonicConfig::Enabled(crcfg) = cfg.cratonic {
                if let (Some(plate_id), Some(orig)) = (
                    &baseline.final_state.plate_id,
                    &original_per_plate_type,
                ) {
                    let new_factor = recompute_cratonic_factor_for_cycle(
                        plate_id,
                        orig,
                        &baseline.final_state.s_field,
                        sea_level_ref,
                        &crcfg,
                    );
                    if let Some(old_factor) = &baseline.final_state.cratonic_factor {
                        craton_change = Some(measure_craton_change(old_factor, &new_factor));
                    }
                    baseline.final_state.cratonic_factor = Some(new_factor);
                }
            }

            CycleOutput {
                baseline,
                common: CycleOutputCommon {
                    erosion_volume_removed: stats.total_eroded,
                    erosion_peak_delta_h: stats.peak_delta_h,
                    sea_level_normalized: sea_level_ref,
                    mass_drift: mass_after - mass_before,
                    craton_recomputation_change: craton_change,
                },
            }
        }
    }
}

/// Run the Phase A multi-cycle loop.
///
/// `Disabled` → exactly one cycle (single `run_baseline` passthrough).
/// The `&mut` requirement is preserved on this branch even though no
/// mutation actually fires, because the [`WorkflowConfig::Enabled`]
/// branch must mutate `cfg.continuation` between cycles to wire the
/// D3 warm-start contract.
///
/// `Enabled(params)` → loop `params.phase_a.n_cycles` cycles. After
/// each cycle (except the last) the loop sets
/// `cfg.continuation = Some(final_state_to_continuation(...))` so the
/// next cycle's `run_baseline` warm-starts from the prior cycle's
/// post-erosion state. The S̃ field, velocity, age and cratonic
/// factor all thread through (D3 contract pinned by
/// `v2_workflow_continuation_no_transient`).
///
/// `cfg.steps` is consumed as the number of tectonic steps per cycle.
/// The convention is to set `cfg.steps = params.phase_a.k_cycle`
/// before calling, but the loop does not enforce this — the two are
/// independently configurable.
pub fn run_phase_a_loop(cfg: &mut BaselineConfig, wf: &WorkflowConfig) -> PhaseAOutput {
    match wf {
        WorkflowConfig::Disabled => {
            let cycle = run_phase_a_cycle(cfg, wf);
            PhaseAOutput { cycles: vec![cycle] }
        }
        WorkflowConfig::Enabled(params) => {
            let n_cycles = params.phase_a.n_cycles.max(1);
            let mut cycles: Vec<CycleOutput> = Vec::with_capacity(n_cycles);
            for cycle_idx in 0..n_cycles {
                let cycle = run_phase_a_cycle(cfg, wf);
                // Set up the next cycle's warm-start *before* moving
                // `cycle` into the output vec. Skip the last cycle:
                // there is no next cycle to warm-start.
                if cycle_idx + 1 < n_cycles {
                    cfg.continuation =
                        Some(final_state_to_continuation(&cycle.baseline.final_state));
                }
                cycles.push(cycle);
            }
            PhaseAOutput { cycles }
        }
    }
}

/// Build a [`ContinuationState`] from a [`FinalState`].
///
/// The orchestrator (Phase 4) calls this at the end of cycle `N` to
/// build the input for cycle `N+1`'s `BaselineConfig.continuation`.
/// Step 8.6's `ContinuationState` carries everything `run_baseline`
/// needs to short-circuit re-init: `s, vx, vy, age, cratonic_factor`.
/// The Voronoï tessellation is implicitly preserved by the run's
/// `BoundaryConfig` (static for the run lifetime).
///
/// D3 contract: cycle `N+1` step 1 should produce a peak|v| within
/// 10 % of cycle `N` step `k_cycle` — pinned by the
/// `v2_workflow_continuation_no_transient` test.
pub fn final_state_to_continuation(fs: &FinalState) -> ContinuationState {
    ContinuationState {
        s: fs.s_field.clone(),
        vx: fs.vx.clone(),
        vy: fs.vy.clone(),
        age: fs.age_field.clone(),
        cratonic_factor: fs.cratonic_factor.clone(),
    }
}

/// Recover `per_plate_type[p]` from the per-cell `plate_type` field.
///
/// At `run_baseline` start the per-cell field is a deterministic
/// broadcast of `per_plate_type` (every cell of plate `p` has type
/// `per_plate_type[p]`). Sampling the first occurrence of each
/// `plate_id` is therefore a faithful recovery — the function returns
/// `Vec<PlateType>` of length `max(plate_id) + 1`.
///
/// A more defensive variant could sample by majority vote across cells
/// of each plate; not necessary at cycle 1 (broadcast is exact), and
/// for cycle 2+ the caller passes the *original* per_plate_type from a
/// prior call rather than re-sampling.
fn extract_per_plate_type(plate_id: &PlateIdField, plate_type: &PlateTypeField) -> Vec<PlateType> {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let max_id = plate_id.data().iter().copied().max().unwrap_or(0) as usize;
    let num_plates = max_id + 1;
    let mut per = vec![PlateType::Oceanic; num_plates];
    let mut seen = vec![false; num_plates];
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_id.get(i, j) as usize;
            if !seen[pid] {
                per[pid] = plate_type.get(i, j);
                seen[pid] = true;
            }
        }
    }
    per
}

/// Step 12 D4 cratonic recompute: synthesise an updated
/// [`VoronoiPlates`] in which each plate's `per_plate_type` reflects
/// the post-erosion D4 retention rule, then call
/// [`build_cratonic_factor_field`].
///
/// **D4 retention rule** (per `step12_issue.md`, Phase 3.5 disambiguation):
///
/// ```text
/// continental_fraction[p] = cells in plate p with S̃ > sea_level
///                         / total cells in plate p   (within-plate)
/// retained[p] = was_continental_at_init[p]
///             && continental_fraction[p] >= craton_retention_threshold
/// ```
///
/// The threshold is [`CratonicConfigEnabled::craton_retention_threshold`]
/// (Semantics 2, within-plate), **not** `plate_area_min` (Semantics 1,
/// fraction-of-domain). The two are independently configurable since
/// Step 12 Phase 3.5; the Phase 3 implementation overloaded
/// `plate_area_min` and could not isolate the D4 flip from the init
/// exclusion. See `step12_issue.md` D4 + the Phase 3.5 commit message
/// for the disambiguation rationale.
///
/// `was_continental_at_init` comes from `initial_per_plate_type` — the
/// caller is responsible for capturing this from the pre-erosion
/// `plate_type` field. After Phase 4 multi-cycle wiring, this becomes
/// the prior cycle's `per_plate_type`, threading the D4 retention rule
/// across the loop.
fn recompute_cratonic_factor_for_cycle(
    plate_id: &PlateIdField,
    initial_per_plate_type: &[PlateType],
    s: &Field2D,
    sea_level_reference: f64,
    cfg: &CratonicConfigEnabled,
) -> Field2D {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let num_plates = initial_per_plate_type.len();

    // Per-plate continental cell count under the new sea_level.
    let mut cont_count = vec![0u32; num_plates];
    let mut total_count = vec![0u32; num_plates];
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_id.get(i, j) as usize;
            total_count[pid] += 1;
            if s.get(i, j) > sea_level_reference {
                cont_count[pid] += 1;
            }
        }
    }

    // Apply D4 retention rule to derive the new per_plate_type. The
    // threshold is `craton_retention_threshold` (Semantics 2, within
    // plate), distinct from `plate_area_min` (Semantics 1, fraction
    // of domain) which the init-time `build_cratonic_factor_field`
    // uses internally.
    let mut updated_per_plate_type: Vec<PlateType> = Vec::with_capacity(num_plates);
    for p in 0..num_plates {
        let frac = if total_count[p] > 0 {
            cont_count[p] as f64 / total_count[p] as f64
        } else {
            0.0
        };
        let was_continental = matches!(initial_per_plate_type[p], PlateType::Continental);
        let new_type = if was_continental && frac >= cfg.craton_retention_threshold {
            PlateType::Continental
        } else {
            PlateType::Oceanic
        };
        updated_per_plate_type.push(new_type);
    }

    // Mirror per_plate_type into a per-cell PlateTypeField so the
    // synthesised `VoronoiPlates` is internally consistent (the BFS
    // in `build_cratonic_factor_field` reads `retained[plate_id]`
    // which is per-plate, but the test surface
    // `factor_zero_on_oceanic_cells` expects the per-cell field to
    // match for diagnostic coherence).
    let mut updated_plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_id.get(i, j) as usize;
            updated_plate_type.set(i, j, updated_per_plate_type[pid]);
        }
    }

    let plates = VoronoiPlates {
        num_plates,
        plate_id: plate_id.clone(),
        plate_type: updated_plate_type,
        per_plate_type: updated_per_plate_type,
        // `seed_coords` is unused by `build_cratonic_factor_field` —
        // empty Vec is safe.
        seed_coords: Vec::new(),
    };
    build_cratonic_factor_field(&plates, cfg)
}

/// Fraction of cells whose `cratonic_factor` changed by more than
/// `1e-9` between two snapshots. The threshold is well above the
/// rounding error of the smoothstep (a smooth function of distances)
/// while small enough to catch any genuine BFS reshuffle from the D4
/// retention rule kicking in or out for a plate.
fn measure_craton_change(old: &Field2D, new: &Field2D) -> f64 {
    let n = old.data().len();
    if n == 0 {
        return 0.0;
    }
    let mut changed = 0_usize;
    for k in 0..n {
        if (old.data()[k] - new.data()[k]).abs() > 1e-9 {
            changed += 1;
        }
    }
    changed as f64 / n as f64
}
