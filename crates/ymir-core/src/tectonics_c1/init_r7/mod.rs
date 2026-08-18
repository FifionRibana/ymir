//! R7 generalised init for C1 Phase 2 Track B (Issue #131).
//!
//! Sibling of [`super::init`] (Phase 1.1 init, preserved verbatim).
//! Phase 2 Track B introduces three sub-components that compose on
//! top of the v2 Voronoï output:
//!
//! 1. **Boundary displacement** ([`boundary_displacement`]) — apply
//!    Perlin / Simplex noise displacement to the per-cell sampling
//!    position before re-querying the nearest Voronoï seed. Produces
//!    non-rectilinear plate boundaries while preserving the seed-
//!    based plate identity. Resolves the v1 / v2 visual failure mode
//!    of orogenic chains aligned along straight Voronoï edges
//!    (§6.1 design doc). Documented as the "boundary displacement"
//!    option of §6.1 — Lloyd relaxation and multi-scale overlay
//!    remain deferred.
//! 2. **Continental clustering** (Stage E2, separate file) — BFS
//!    cluster-based plate-type assignment producing a cadrable
//!    continental cluster.
//! 3. **Ridge-aligned age = 0** (Stage E3, separate file) — set
//!    `age = 0` on cells adjacent to divergent boundaries at init
//!    time. Resolves the Phase 2 Track A finding that flux-form
//!    advection of `age` produces ~1000× density pile-up at
//!    convergent boundaries (`feedback_age_advection_density_vs_lagrangian`).
//!
//! ## Determinism
//!
//! All R7 init sub-components are deterministic given
//! `(grid_size, params.seed)`. Stochastic elements:
//!
//! - **Boundary displacement**: Perlin / Simplex noise via two
//!   independent `Fbm<Perlin>` instances (one per displacement
//!   component) seeded from `params.seed`. Same seed → same noise
//!   field → same per-cell displacement.
//! - **Continental clustering**: ChaCha8Rng seeded from
//!   `cluster_params.seed` for seed-pick selection. BFS expansion
//!   itself is deterministic given the adjacency graph.
//!
//! No floating-point reproducibility caveats: `f64` arithmetic and
//! `noise::Fbm<Perlin>` are bit-deterministic on a given target
//! triple per the existing Phase 1.1 + Track A precedents.
//!
//! ## Composition with Phase 1.1 init
//!
//! [`super::init::init_c1_state_phase_1_1`] is **preserved
//! verbatim** as the Phase 1.1 regression baseline. Phase 2 Track B
//! ships a parallel entry point (Stage E4) that calls the same v2
//! Voronoï + S̃-init pipeline, then chains the three R7
//! sub-components in order:
//!
//! ```text
//!     generate_voronoi → R7 boundary displacement
//!                      → cluster-BFS plate_type override
//!                      → init_s_field (per overridden plate_type)
//!                      → ridge-aligned age init
//!                      → cratonic mask
//! ```
//!
//! Phase 1.x tests continue to call
//! `init_c1_state_phase_1_1(grid_size, seed)` directly; the Phase
//! 2 entry is a new function call.

pub mod age_init;
pub mod boundary_displacement;
pub mod clustering;
pub mod params;

pub use age_init::{AgeInitParams, init_age_field_ridge_aligned};
pub use boundary_displacement::apply_boundary_displacement;
pub use clustering::{
    ContinentalClusterParams, assign_continental_clusters, build_plate_adjacency,
};
pub use params::R7InitParams;

use crate::tectonics_v2::boundaries::plate_type::PlateTypeField;
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::init::{InitContext, InitMode, PlateInitData, init_s_field};
use crate::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

use super::boundary_classification::classify_boundaries;
use super::init::build_phase_1_1_cratonic_mask;
use super::kinematics::PlateKinematics;
use super::state::{BoolField, C1State};

/// Bundle of all Phase 2 R7 init tunables.
///
/// Composes the three sub-component parameter sets shipped in
/// Stages E1 / E2 / E3:
///
/// - [`r7`](Self::r7) — boundary displacement
///   ([`R7InitParams`])
/// - [`cluster`](Self::cluster) — continental clustering
///   ([`ContinentalClusterParams`])
/// - [`age`](Self::age) — ridge-aligned age = 0
///   ([`AgeInitParams`])
///
/// `Default::default()` enables all three sub-components with
/// the calibrated defaults documented in each sub-module. The
/// caller is responsible for setting sub-channel seeds
/// (`r7.seed`, `cluster.seed`) explicitly if a per-channel
/// override of the main `seed` argument is desired — the
/// dispatcher does **not** thread the main `seed` into the
/// sub-channel seeds, separating concerns.
#[derive(Clone, Copy, Debug)]
pub struct Phase2InitParams {
    /// R7 boundary displacement (sub-component 1).
    pub r7: R7InitParams,
    /// Continental clustering (sub-component 2).
    pub cluster: ContinentalClusterParams,
    /// Ridge-aligned age = 0 init (sub-component 3, Path 3.A).
    pub age: AgeInitParams,
    /// #155 A′ — cratonic crustal-thickness ratio. After the cratonic
    /// mask is built, continental cratonic cells get their initial S̃
    /// multiplied by this (Airy: thick craton crust → elevated). Default
    /// 1.25 (canonical C1; real crustal ratio ~40 km craton / ~32 km
    /// platform). Set to 1.0 to disable (byte-identical to pre-#155 init).
    /// Pairs with `ErosionParams::craton_resist` so the differential is
    /// not planed/inverted by erosion (the measured 2026 inversion).
    ///
    /// MAGNITUDE (measured, defaults A′-on): cratons render ~600-1100 m
    /// above non-craton land (conversion-dependent — the robust measure is
    /// the ~0.134 normalised difference / S̃ ratio ~1.5-1.8; the metres
    /// figure depends on the unpinned norm→m vertical scale). NOT
    /// geologically wrong — real cratons include high plateaus (~1000-
    /// 1500 m, e.g. southern African); the 300-500 m target was the MOST
    /// WORN shields. Lowering it = a MODEL refinement (worn-init / Jordan
    /// compositional isostasy / cumulative wear over eons — 300 steps ≠
    /// eons), NOT a knob (params stay anchored: 1.25 = crustal ratio,
    /// craton_resist mid-band). Worn-shield height = documented follow-up.
    pub craton_thickness_ratio: f64,
    /// Number of Voronoi plates (`VoronoiConfig::num_plates`). Default 8 (the
    /// TDD §3.4 `[5,15]` mid-point). A fragmentation lever for the "island
    /// continent" budget (M1 #190): more plates → smaller plates, and — paired
    /// with `cluster.seed_cluster_count > 1` — smaller separate landmasses at a
    /// constant land fraction. Folded into the tectonic cache key.
    pub num_plates: usize,
    /// #165 bimodal shield/platform — fraction of the CRATONIC AREA rendered as
    /// exposed HIGH SHIELD; the rest becomes LOW PLATFORM. The cratonic mask sets
    /// the AREA (old Precambrian lithosphere, ~50-70 % of continents — realistic),
    /// but the model rendered ALL of it as high shield, whereas Earth's exposed
    /// shield is only ~10-20 % of cratonic area; the rest (~80-90 %) is sediment-
    /// covered LOW platform (measured: high cratons were 58 % of emergent land at
    /// median 1814 m → the bulk of the missing-low-plain deficit, and lowering
    /// them uniformly hit a sea-level-percentile coupling — see
    /// `floor_thickness/VERDICT_fraction.md`). `Some(f)` keeps only ~`f` of the
    /// cratonic cells as shield (the high thick+resistant+dense treatment); the
    /// platform cells drop OUT of the mask and become normal low continental crust
    /// (already Earth-like). `None` → all cratonic cells are shield (byte-identical
    /// to pre-#165). Anchored on the exposed-shield fraction (~0.15), NOT tuned to
    /// a hypsometry target.
    pub craton_shield_fraction: Option<f64>,
}

impl Default for Phase2InitParams {
    fn default() -> Self {
        Self {
            r7: R7InitParams::default(),
            cluster: ContinentalClusterParams::default(),
            age: AgeInitParams::default(),
            num_plates: 8,
            craton_thickness_ratio: 1.25,
            // #165 bimodal shield/platform — ~15 % of the cratonic AREA stays
            // exposed HIGH shield, the rest becomes LOW platform (normal crust).
            // Anchored on Earth's exposed-shield share (~10-20 % of cratonic
            // area); VERIFIED via probe_craton_bimodal: lifts ALL-land<500m from
            // ~25 % to ~45 % (seed 42) with the platform/non-craton floor stable
            // (the sea-level-percentile coupling dissolved). None → all-shield
            // (byte-identical to pre-#165). Pairs with `craton_thickness_ratio`
            // (the high shields stay thick) + isostasy `craton_rho_crust`.
            craton_shield_fraction: Some(0.15),
        }
    }
}

/// Build a fresh [`C1State`] using the Phase 2 R7 init pipeline.
///
/// **New parallel entry — NOT a default-redirection of Phase
/// 1.1 init.** Phase 1.x tests continue to call
/// [`super::init::init_c1_state_phase_1_1`] directly; new Phase 2
/// code calls this function. The two functions co-exist by design
/// — see [`super`]'s module docstring on "Composition with Phase
/// 1.1 init".
///
/// ## Per-step pipeline (10 steps)
///
/// 1. `generate_voronoi(nx, ny, &VoronoiConfig::default(), seed)`
///    → `VoronoiPlates { plate_id, plate_type, per_plate_type,
///    seed_coords, num_plates }`. Phase 1.x baseline tessellation.
/// 2. [`apply_boundary_displacement`] on `plate_id` → curved
///    boundaries (sub-component 1).
/// 3. [`build_plate_adjacency`] on the displaced `plate_id` →
///    per-plate adjacency graph.
/// 4. [`assign_continental_clusters`] on `per_plate_type` via the
///    adjacency graph → cadrable cluster (sub-component 2).
/// 5. Broadcast the updated `per_plate_type` back into a
///    cell-level [`PlateTypeField`] via `plate_id` lookup.
/// 6. [`init_s_field`] with `InitMode::default()` using the
///    overridden `plate_type` — Phase 1.x baseline `S̃` per
///    plate-type-broadcasted field.
/// 7. Build kinematics via
///    [`PlateKinematics::preset_phase_1_1`] (Phase 1.x preset;
///    Track C/D will replace once those tracks land).
/// 8. [`classify_boundaries`] on `(plate_id, kinematics)` →
///    `BoundaryInfo` for stage 9.
/// 9. [`init_age_field_ridge_aligned`] with the trichotomy
///    decision tree (continental > divergent > oceanic) →
///    `age` field with ridge-aligned 0 at oceanic-divergent
///    cells (sub-component 3 / Path 3.A).
/// 10. `build_phase_1_1_cratonic_mask` — same rule as Phase
///     1.1, reused via the `pub(crate)`-promoted helper in
///     [`super::init`] (single-source-of-truth invariant).
///
/// ## Determinism
///
/// Given `(grid_size, seed, params)`, this function produces a
/// bit-identical `C1State`. Stochastic sub-channels (R7 FBM noise,
/// clustering seed pick) seed independently from
/// `params.r7.seed` and `params.cluster.seed`; the main `seed`
/// argument drives only the Voronoï tessellation. Caller-side
/// control of sub-channel seeds preserves user freedom to vary
/// noise / clustering texture independently of plate geometry.
///
/// ## Out of scope (Phase 2 Track B)
///
/// - **Constrained kinematics** (Track D): step 7 uses the Phase
///   1.1 preset; not every Voronoï layout is guaranteed to have
///   a divergent boundary, so the age=0 ridge cells may be empty
///   in some seeds (W7 architectural concern surfaced for
///   Stage V).
/// - **Boundary evolution** (Track C): the post-init kinematics
///   are static (Phase 1.1 contract). Track C will introduce
///   dynamic boundary updates.
/// #165 bimodal shield/platform — keep ~`fraction` of the cratonic cells as
/// exposed SHIELD (returned `true`); the rest become platform (`false`, dropped
/// from the high treatment). Selection is a deterministic, spatially COHERENT
/// value-noise mosaic (a coarse seeded lattice, bilinearly interpolated) so the
/// shields form contiguous patches (not salt-and-pepper) and the proportion is
/// hit by thresholding at the `fraction` upper quantile over the cratonic cells
/// — exact-ish regardless of seed / cratonic area. Distinct sub-seed so the
/// mosaic does not correlate with the tectonic / R7 channels.
fn select_shield_mask(
    craton: &BoolField,
    nx: usize,
    ny: usize,
    seed: u64,
    fraction: f64,
) -> BoolField {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    let frac = fraction.clamp(0.0, 1.0);
    // Coarse value-noise lattice (~1 node / 8 cells) → smooth coherent field.
    let lat = (nx / 8).max(2);
    let laty = (ny / 8).max(2);
    let mut rng = ChaCha8Rng::seed_from_u64(seed ^ 0x5165_B1_0DA1_C0DE);
    let nodes: Vec<f32> = (0..(lat + 1) * (laty + 1)).map(|_| rng.random::<f32>()).collect();
    let sample = |i: usize, j: usize| -> f32 {
        let gx = i as f32 * lat as f32 / nx as f32;
        let gy = j as f32 * laty as f32 / ny as f32;
        let x0 = gx.floor() as usize;
        let y0 = gy.floor() as usize;
        let tx = gx - x0 as f32;
        let ty = gy - y0 as f32;
        let at = |a: usize, b: usize| nodes[b * (lat + 1) + a];
        let v0 = at(x0, y0) + (at(x0 + 1, y0) - at(x0, y0)) * tx;
        let v1 = at(x0, y0 + 1) + (at(x0 + 1, y0 + 1) - at(x0, y0 + 1)) * tx;
        v0 + (v1 - v0) * ty
    };

    // Threshold = the (1 − frac) quantile of the noise over the cratonic cells.
    let mut vals: Vec<f32> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            if craton.get(i, j) {
                vals.push(sample(i, j));
            }
        }
    }
    if vals.is_empty() {
        return BoolField::filled(nx, ny, false);
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let k = (((1.0 - frac) * vals.len() as f64).floor() as usize).min(vals.len() - 1);
    let threshold = vals[k];

    let mut shield = BoolField::filled(nx, ny, false);
    for j in 0..ny {
        for i in 0..nx {
            if craton.get(i, j) && sample(i, j) >= threshold {
                shield.set(i, j, true);
            }
        }
    }
    shield
}

pub fn init_c1_state_phase_2_r7(grid_size: usize, seed: u64, params: &Phase2InitParams) -> C1State {
    let nx = grid_size;
    let ny = grid_size;

    // Step 1 — Voronoï tessellation. `num_plates` from the init params (default
    // 8); 30 % continental Bernoulli — the clustering step overrides
    // `per_plate_type` so the Bernoulli ratio is discarded; the Voronoï output is
    // still load-bearing for `plate_id`, `seed_coords`, and `num_plates`.
    let vor_config = VoronoiConfig { num_plates: params.num_plates, ..VoronoiConfig::default() };
    let voronoi = generate_voronoi(nx, ny, &vor_config, seed);
    let mut plate_id = voronoi.plate_id;
    let seed_coords = voronoi.seed_coords;
    let num_plates = voronoi.num_plates;

    // Step 2 — R7 boundary displacement. No-op if
    // `params.r7.enabled = false`.
    apply_boundary_displacement(&mut plate_id, &seed_coords, &params.r7);

    // Step 3 — adjacency from the (possibly displaced)
    // `plate_id`. Periodic `rem_euclid` wraparound is handled
    // inside `build_plate_adjacency`.
    let adjacency = build_plate_adjacency(&plate_id, num_plates);

    // Step 4 — cluster-based continental type override.
    let mut per_plate_type = voronoi.per_plate_type;
    assign_continental_clusters(&mut per_plate_type, &adjacency, &params.cluster);

    // Step 5 — broadcast the updated `per_plate_type` back into
    // cell-level `PlateTypeField`.
    let mut plate_type = PlateTypeField::filled(nx, ny, per_plate_type[0]);
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_id.get(i, j) as usize;
            plate_type.set(i, j, per_plate_type[pid]);
        }
    }

    // Step 6 — `S̃` field via v2's init dispatch using the
    // overridden plate_type. Same `InitMode::default()` as
    // Phase 1.1 (Uniform with smoothstep at boundaries).
    let init_ctx = InitContext {
        nx,
        ny,
        seed,
        amplitude: 0.0,
        plate_data: Some(PlateInitData {
            plate_id: &plate_id,
            plate_type: &plate_type,
            seed_coords: Some(&seed_coords),
        }),
    };
    let mut s = init_s_field(InitMode::default(), &init_ctx);

    // Step 7 — kinematics (Phase 1.1 preset; Track C/D
    // placeholder per Stage E4 W7 concern).
    let kinematics = PlateKinematics::preset_phase_1_1(num_plates);

    // Step 8 — boundary classification for ridge-aligned age
    // detection.
    let boundary_info = classify_boundaries(&plate_id, &kinematics);

    // Step 9 — age field with trichotomy (continental >
    // divergent > oceanic).
    let mut age = Field2D::new(nx, ny);
    init_age_field_ridge_aligned(&mut age, &plate_id, &plate_type, &boundary_info, &params.age);

    // Step 10 — cratonic mask (Phase 1.1 rule, reused via the
    // pub(crate)-promoted helper). This is the cratonic AREA.
    let cratonic_mask = build_phase_1_1_cratonic_mask(nx, ny, &plate_id, &plate_type, &seed_coords);

    // Step 10.5 (#165 bimodal) — narrow the cratonic AREA to the exposed SHIELD.
    // Only ~`craton_shield_fraction` of cratonic cells stay high (they keep the
    // thick+resist+dense treatment below + downstream); the rest become low
    // PLATFORM, dropping out of the mask → normal low continental crust. `None`
    // → all cratonic cells stay shield (byte-identical to pre-#165). See the
    // `craton_shield_fraction` docstring.
    let cratonic_mask = match params.craton_shield_fraction {
        Some(f) => select_shield_mask(&cratonic_mask, nx, ny, seed, f),
        None => cratonic_mask,
    };

    // Step 11 (#155 A′) — cratonic crustal-thickness differential. Thick
    // cratonic crust (Airy isostasy → elevated worn shields). Multiply the
    // initial continental cratonic S̃ by the ratio (mask is built only on
    // continental cells). ratio == 1.0 → no-op (byte-identical pre-#155).
    if params.craton_thickness_ratio != 1.0 {
        for j in 0..ny {
            for i in 0..nx {
                if cratonic_mask.get(i, j) {
                    s.set(i, j, s.get(i, j) * params.craton_thickness_ratio);
                }
            }
        }
    }

    C1State {
        s,
        age,
        plate_id,
        plate_type,
        cratonic_mask,
        num_plates,
        last_step_stats: Default::default(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Dispatcher smoke test — produces a valid `C1State` with
    /// the right dimensions and num_plates. Deeper assertions
    /// (R7 effects, clustering effects, ridge cells present)
    /// land in Stage V.
    #[test]
    fn phase_2_r7_init_dispatcher_smoke() {
        let state = init_c1_state_phase_2_r7(64, 42, &Phase2InitParams::default());
        assert_eq!(state.num_plates, 8);
        assert_eq!(state.s.nx(), 64);
        assert_eq!(state.s.ny(), 64);
        assert_eq!(state.age.nx(), 64);
        assert_eq!(state.age.ny(), 64);
        assert_eq!(state.plate_id.nx(), 64);
        assert_eq!(state.plate_id.ny(), 64);
        assert_eq!(state.plate_type.nx(), 64);
        assert_eq!(state.plate_type.ny(), 64);
        assert_eq!(state.cratonic_mask.nx(), 64);
        assert_eq!(state.cratonic_mask.ny(), 64);

        // All fields finite.
        for &v in state.s.data() {
            assert!(v.is_finite(), "non-finite S̃ value in dispatcher output");
        }
        for &v in state.age.data() {
            assert!(v.is_finite(), "non-finite age value in dispatcher output");
        }
    }

    /// Determinism — same `(grid_size, seed, params)` produces
    /// bit-identical `C1State`. Critical regression guard against
    /// any future refactor introducing non-determinism in the
    /// dispatcher's sub-channel coupling.
    #[test]
    fn phase_2_r7_init_dispatcher_deterministic() {
        let params = Phase2InitParams::default();
        let state_a = init_c1_state_phase_2_r7(64, 42, &params);
        let state_b = init_c1_state_phase_2_r7(64, 42, &params);

        for j in 0..64 {
            for i in 0..64 {
                assert_eq!(state_a.s.get(i, j), state_b.s.get(i, j), "S̃ mismatch at ({i}, {j})");
                assert_eq!(
                    state_a.age.get(i, j),
                    state_b.age.get(i, j),
                    "age mismatch at ({i}, {j})"
                );
                assert_eq!(
                    state_a.plate_id.get(i, j),
                    state_b.plate_id.get(i, j),
                    "plate_id mismatch at ({i}, {j})"
                );
            }
        }
    }
}
