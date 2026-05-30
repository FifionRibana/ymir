//! Shared paradigm-agnostic post-tectonic-step pass for Phase A.
//!
//! Both the v2 path ([`super::phase_a_v2`]) and the C1 path
//! ([`super::phase_a_c1`], Commit 4 of Issue #125 H2) need the same
//! per-cycle pipeline AFTER the tectonic step has updated the
//! `S̃` field:
//!
//! 1. **Sea-level threshold** — compute the adaptive `sea_level_ref`
//!    in `S̃` space from the current min/max of `S̃` and
//!    `IsostasyConfig::sea_level_fraction` (the Phase 3.5 fix —
//!    the original heightmap-space sea_level produced wrong
//!    thresholds for the reclassification rule).
//! 2. **Macro mass redistribution** — Step 12 R3 erosion + drainage
//!    + isostatic rebound, mass-conserving. Mutates `S̃` in place;
//!    surfaces `erosion_volume_removed` + `erosion_peak_delta_h`
//!    in [`CycleOutputCommon`].
//! 3. **Reclassification** — per-cell `plate_type[i,j]` =
//!    `Continental` iff `S̃[i,j] > sea_level_ref`, else `Oceanic`.
//!    The Voronoi tessellation (`plate_id`) is static for the run;
//!    only the per-cell type changes.
//! 4. **Cratonic factor recompute** — D4 retention rule: a plate
//!    keeps its continental status iff (a) it was continental at
//!    init, AND (b) the post-erosion continental cell fraction
//!    inside the plate is `≥ craton_retention_threshold`. Rebuilds
//!    the `cratonic_factor: Field2D` via the Step 9 BFS smoothstep
//!    algorithm.
//!
//! All four steps use **paradigm-agnostic types**: `Field2D`,
//! `PlateIdField`, `PlateTypeField`, `PlateType`, `PhaseAParams`,
//! `IsostasyConfig`, `CratonicConfigEnabled`. The
//! `tectonics_v2::cratonic` namespace is data-only and is **not**
//! gated under `v2_legacy` (it has no Stokes / mantle / slab
//! dependency); see `docs/migrations/v2_to_c1_attic.md` § §4.8.
//!
//! ## Why this is a shared pass, not a trait
//!
//! Per the Phase 1.3 H1 audit `docs/migrations/
//! harness_paradigm_agnostic.md` (Option B chosen), the v2 and C1
//! tectonic *runners* have structurally different signatures (v2
//! takes 23-field `BaselineConfig` + returns `BaselineResult`; C1
//! mutates `&mut C1State` in place). But once the tectonic step
//! has updated the `S̃` field, **the post-step pipeline is
//! identical** — it operates on `S̃`, `plate_id`, `plate_type`,
//! and the per-plate continental classification. Extracting that
//! pipeline as a function takes care of W2 / W3 / R3 from the
//! audit doc: single implementation, two callers, no silent
//! divergence of the shared diagnostics
//! ([`CycleOutputCommon::mass_drift`],
//! [`CycleOutputCommon::erosion_volume_removed`], etc.).
//!
//! ## What stays paradigm-specific
//!
//! - v2's velocity-reset post-macro-redistribution (the
//!   counter-intuitive D1-ter `v = 0` empirical finding) stays in
//!   `phase_a_v2.rs` — C1 has no Stokes velocity field.
//! - v2's `BaselineResult::final_state` install of the new
//!   `cratonic_factor` stays in `phase_a_v2.rs`. The factor itself
//!   is paradigm-agnostic (`Field2D`), but where to store it is
//!   not (`baseline.final_state.cratonic_factor` for v2; C1 has
//!   no `cratonic_factor` field today — its `C1State.cratonic_mask`
//!   is a `BoolField`, a different data shape).
//!
//! Both paradigm-specific bits are surfaced via
//! [`PostTectonicOutput`]: the new `cratonic_factor: Option<Field2D>`
//! is returned, and the caller installs it (or ignores it) per its
//! state shape. C1 in Phase 1.3 has no cratonic-factor evolution,
//! so it discards. Phase 1.4+ may add one.

use crate::tectonics::isostasy::IsostasyConfig;
use crate::tectonics_v2::boundaries::{PlateType, PlateTypeField};
use crate::tectonics_v2::cratonic::factor::build_cratonic_factor_field;
use crate::tectonics_v2::cratonic::CratonicConfigEnabled;
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::voronoi::{PlateIdField, VoronoiPlates};

use super::{macro_redistribution, CycleOutputCommon, PhaseAParams};

/// Inputs to [`apply_post_tectonic`].
///
/// All types are paradigm-agnostic (used identically by the v2 and
/// C1 Phase A paths). The caller is responsible for capturing
/// `initial_per_plate_type` from the pre-reclassification
/// `plate_type` field via [`extract_per_plate_type`] **before**
/// invoking this function.
pub struct PostTectonicInput<'a> {
    /// `S̃` field, mutated in place by the macro-redistribution
    /// step.
    pub s_field: &'a mut Field2D,
    /// Voronoi plate-id field. Static for the run. `None` →
    /// reclassification and craton recompute skipped (no plate
    /// concept at this call site).
    pub plate_id: Option<&'a PlateIdField>,
    /// Per-cell plate-type field, mutated in place by the
    /// reclassification step. `None` → step skipped.
    pub plate_type: Option<&'a mut PlateTypeField>,
    /// Previous-cycle cratonic factor field (for the
    /// `craton_recomputation_change` diagnostic). `None` → no
    /// change measurement, only forward computation.
    pub previous_cratonic_factor: Option<&'a Field2D>,
    /// Per-plate type at run init (D4 "was continental at init"
    /// gate). Caller-captured via [`extract_per_plate_type`].
    /// `None` → craton recompute skipped.
    pub initial_per_plate_type: Option<&'a [PlateType]>,
    /// Phase A loop parameters (alpha, isostatic_rebound_ratio,
    /// max_drainage_distance) consumed by
    /// [`super::macro_redistribution::apply`].
    pub params: &'a PhaseAParams,
    /// Isostasy configuration. The function reads only
    /// `sea_level_fraction` from this — the Phase 3.5 S̃-space
    /// sea-level formula.
    pub iso_cfg: &'a IsostasyConfig,
    /// Cratonic configuration. `None` → craton recompute skipped.
    /// The `tectonics_v2::cratonic` namespace is data-only and
    /// not v2_legacy-gated (no Stokes coupling), so this is
    /// available to both paradigms.
    pub cratonic_cfg: Option<&'a CratonicConfigEnabled>,
}

/// Output of [`apply_post_tectonic`].
///
/// `common` is the paradigm-agnostic per-cycle diagnostics bundle.
/// `new_cratonic_factor` is the freshly-recomputed factor field; the
/// caller installs it where appropriate for its state shape
/// (v2: `baseline.final_state.cratonic_factor`; C1: discards in
/// Phase 1.3, may consume later).
pub struct PostTectonicOutput {
    pub common: CycleOutputCommon,
    pub new_cratonic_factor: Option<Field2D>,
}

/// Run the shared post-tectonic pass: sea-level → macro-redistribution
/// → reclassification → cratonic recompute. See module docstring
/// for the rationale and the v2/C1 boundary.
pub fn apply_post_tectonic(mut input: PostTectonicInput<'_>) -> PostTectonicOutput {
    // Step 1 — adaptive sea-level threshold in S̃ space.
    //
    // Phase 3.5 fix (see `phase_a_v2.rs` history for the pre-3.5
    // bug): compute the threshold *in S̃ units* via the isostasy
    // formula `h_sea = h_min + sea_level_fraction · h_range`,
    // applied to `S̃` directly:
    //
    //     s_sea = s_min + sea_level_fraction · (s_max - s_min)
    //
    // At init (S̃ ∈ [≈0.2, ≈1.2]) this resolves to ≈ 0.6 — close to
    // the natural 0.5 midpoint between oceanic and continental, and
    // adaptive to the S̃ distribution as cycles erode.
    let sea_level_ref = compute_sea_level_ref_s_space(input.s_field, input.iso_cfg);

    // Step 2 — macro mass redistribution. Step 12 R3 replaced the
    // legacy `low_res_erosion::apply` with this drainage + isostatic-
    // rebound formulation; mass drift is now ~ IEEE-754 floor by
    // construction.
    let mass_before: f64 = input.s_field.data().iter().sum();
    let stats = macro_redistribution::apply(input.s_field, input.params, sea_level_ref);
    let mass_after: f64 = input.s_field.data().iter().sum();

    // Step 3 — reclassify per-cell `plate_type` from the new
    // `s_field` against `sea_level_ref`.
    if let Some(plate_type) = input.plate_type.as_deref_mut() {
        reclassify_inplace(plate_type, input.s_field, sea_level_ref);
    }

    // Step 4 — recompute cratonic factor when the cfg + plate_id +
    // initial_per_plate_type are all present.
    let mut new_cratonic_factor: Option<Field2D> = None;
    let mut craton_change: Option<f64> = None;
    if let (Some(crcfg), Some(plate_id), Some(orig)) = (
        input.cratonic_cfg,
        input.plate_id,
        input.initial_per_plate_type,
    ) {
        let factor = recompute_cratonic_factor_for_cycle(
            plate_id,
            orig,
            input.s_field,
            sea_level_ref,
            crcfg,
        );
        if let Some(old_factor) = input.previous_cratonic_factor {
            craton_change = Some(measure_craton_change(old_factor, &factor));
        }
        new_cratonic_factor = Some(factor);
    }

    PostTectonicOutput {
        common: CycleOutputCommon {
            erosion_volume_removed: stats.total_eroded,
            erosion_peak_delta_h: stats.peak_delta_h,
            sea_level_normalized: sea_level_ref,
            mass_drift: mass_after - mass_before,
            craton_recomputation_change: craton_change,
        },
        new_cratonic_factor,
    }
}

/// Recover `per_plate_type[p]` from the per-cell `plate_type` field.
///
/// At `run_baseline` start (v2) or `init_c1_state_phase_1_1` (C1)
/// the per-cell field is a deterministic broadcast of `per_plate_
/// type` (every cell of plate `p` has type `per_plate_type[p]`).
/// Sampling the first occurrence of each `plate_id` is therefore a
/// faithful recovery — the function returns `Vec<PlateType>` of
/// length `max(plate_id) + 1`.
///
/// Callers should invoke this **before** [`apply_post_tectonic`] so
/// the captured per-plate type reflects the *pre-reclassification*
/// state, which is what the D4 "was continental at init" gate
/// needs.
pub fn extract_per_plate_type(
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
) -> Vec<PlateType> {
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

/// Phase 3.5 sea-level reference in `S̃` space.
///
/// Computes the adaptive sea-level threshold in **`S̃` units**
/// (not in heightmap `[0, 1]` units) via the isostasy formula
/// `h_sea = h_min + sea_level_fraction · h_range`, applied to
/// the current `S̃` field directly:
///
/// ```text
///     s_sea = s_min + sea_level_fraction · (s_max - s_min)
/// ```
///
/// At Phase 1.1 init (`S̃ ∈ [≈0.2, ≈1.0]`) this resolves to ≈ 0.5
/// — close to the natural midpoint between oceanic and continental
/// values, and **adaptive to the `S̃` distribution as cycles
/// erode**. This is the value that downstream consumers (drainage
/// classification, per-cell reclassification, cratonic recompute,
/// stream-power erosion in Phase 1.4) should use to discriminate
/// continental from oceanic cells inside C1's `S̃`-paradigm time
/// loop.
///
/// **Why this matters:** the alternative — using
/// `compute_isostasy(s).sea_level_normalized` (an `f32` in heightmap
/// `[0, 1]` space) — resolves to ≈ 0.111 for the default
/// `IsostasyConfig` (`max_depth=500 m, max_elevation=4000 m`),
/// well below `S̃`'s natural oceanic floor of ≈ 0.2. The Phase 3
/// (pre-3.5) workflow used the heightmap-space value and every
/// continental/oceanic mask was satisfied by every cell; the
/// reclassification rule never fired and the v2 cratonic-recompute
/// tests passed for the wrong reason.
///
/// **Single source of truth.** Called by
/// [`apply_post_tectonic`] (end-of-cycle in `phase_a_*`) AND by
/// the per-step Phase 1.4 stream-power erosion path (which needs
/// the same threshold to classify drainage targets). Keeping one
/// implementation prevents the two call sites from drifting under
/// future formula refinements.
///
/// The function is read-only on its arguments — pure computation,
/// thread-safe, no caching.
pub fn compute_sea_level_ref_s_space(s: &Field2D, iso_cfg: &IsostasyConfig) -> f64 {
    let s_data = s.data();
    let (s_min, s_max) = s_data.iter().copied().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(lo, hi), v| (lo.min(v), hi.max(v)),
    );
    let s_range = (s_max - s_min).max(1e-10);
    s_min + (iso_cfg.sea_level_fraction as f64) * s_range
}

/// Reclassify per-cell `plate_type` from the post-erosion `S̃`
/// field against `sea_level_ref`.
///
/// `pub(crate)` so the C1 viz facade
/// [`crate::tectonics_c1::reclassify::c1_reclassify_plate_type`]
/// can delegate to this function for snapshot-only reclassification
/// (Viz-0 Stage A bug fix, Issue #137). The C1 viz wraps this
/// behind a C1-facing public name to keep the workflow internals
/// from being called directly by viz-layer code.
pub(crate) fn reclassify_inplace(plate_type: &mut PlateTypeField, s: &Field2D, sea_level_ref: f64) {
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

/// Step 12 D4 cratonic recompute: synthesise an updated
/// [`VoronoiPlates`] in which each plate's `per_plate_type` reflects
/// the post-erosion D4 retention rule, then call
/// [`build_cratonic_factor_field`].
///
/// **D4 retention rule** (per `step12_issue.md`, Phase 3.5
/// disambiguation):
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

    // Apply D4 retention rule to derive the new per_plate_type.
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

    // Mirror `per_plate_type` into a per-cell PlateTypeField so the
    // synthesised `VoronoiPlates` is internally consistent.
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
/// while small enough to catch any genuine BFS reshuffle from the
/// D4 retention rule kicking in or out for a plate.
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
