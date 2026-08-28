//! Step 13 Phase 3 — radial-profile + FBM S̃ initialisation.
//!
//! Extends [`super::radial_profile`] with additive Fractional
//! Brownian Motion (FBM) noise on continental cells, producing
//! intra-plate thickness heterogeneity. Step 11 visual exploration
//! showed continents with quasi-uniform interior thickness; FBM
//! provides the "provinces" texture (older cratonic cores embedded
//! in younger surrounding belts) that real continents exhibit.
//!
//! ## Algorithm (issue D2)
//!
//! 1. Build the radial profile field via [`super::radial_profile::build`]
//!    (Phase 2). Continental cells: `S̃_radial[i, j] = oceanic +
//!    (continental − oceanic) · profile(d / L_plate)`. Oceanic cells:
//!    `S̃ = oceanic_value` flat.
//! 2. Sample isotropic FBM via `noise::Fbm<Perlin>` at every cell:
//!    `n[i, j] = fbm(x_norm / scale, y_norm / scale)` where
//!    `x_norm = (i + 0.5) / nx`, `y_norm = (j + 0.5) / ny` and
//!    `scale` is the largest-feature size in domain fractions
//!    (default `0.10 ⇒ feature size ≈ 1/10 of the domain ≈ 6
//!    cells on 64² grids`, smaller than typical plate `L_plate ≈
//!    10–15 cells` — calibration in the vigilance section
//!    below).
//! 3. Add the FBM perturbation **only on continental cells**:
//!    `S̃[i, j] = clamp(S̃_radial + amplitude · n[i, j], 0, 1)`.
//!    Oceanic cells stay at `oceanic_value` exactly (FBM never
//!    applied — issue D2).
//!
//! ## Vigilance
//!
//! - **`fbm_amplitude` ≤ continental − oceanic threshold**: the
//!   default amplitude `0.10` keeps the FBM perturbation well
//!   within `±0.10` of `S̃_radial`, far above the `0.5` continental
//!   threshold. Phase 5 UI clamps the slider to `[0.0, 0.40]`. With
//!   amplitude > 0.45, some continental cells may end up below
//!   threshold; the algorithm here does NOT clamp the input
//!   amplitude (anti-pattern D7 — no silent clamp), but does clamp
//!   the *output* to `[0, 1]` to keep S̃ in physical range.
//! - **Default calibration (Phase 6 amendment)**: the initial
//!   draft used `fbm_scale = 0.25` and `fbm_amplitude = 0.10`.
//!   Phase 6 acceptance probing on `single_continent` revealed
//!   two compounding issues:
//!
//!   1. `fbm_scale = 0.25` ⇒ wavelength `≈ 16 cells` on a 64²
//!      grid, ≥ typical plate `L_plate ≈ 10–15 cells`. The FBM
//!      does not actually oscillate intra-plate.
//!   2. `noise::Fbm<Perlin>` is auto-normalised and produces
//!      `σ ≈ 0.27 × amplitude`, not the `≈ amplitude / 1.5–2.0
//!      ≈ 0.05–0.07` the issue D2 estimate assumed.
//!
//!   Sweep `(scale × amplitude)` on `single_continent` (Phase 6
//!   probe `fbm_calibration_probe`):
//!
//!   | scale \ amp | 0.10  | 0.15  | 0.20  | 0.25  |
//!   |-------------|-------|-------|-------|-------|
//!   | 0.05        | 0.024 | 0.037 | 0.048 | 0.060 |
//!   | **0.10**    | 0.027 | 0.041 | **0.055** | 0.068 |
//!   | 0.20        | 0.018 | 0.027 | 0.036 | 0.045 |
//!   | 0.25        | 0.018 | 0.027 | 0.036 | 0.046 |
//!
//!   `scale = 0.10` is empirically optimal — `scale = 0.05`
//!   counter-productively introduces high-frequency Perlin grid
//!   artefacts that dilute the large-scale variance.
//!   `amplitude = 0.20` clears the acceptance lower bound
//!   `σ_fbm_isolated ≥ 0.040` with ~35 % margin across all
//!   continental plates of `single_continent`, while staying
//!   well within the Phase 5 UI clamp `[0.0, 0.40]` and well
//!   above the `0.5` continental threshold for interior cells.
//!   Users on larger grids (256²+) may raise `fbm_scale` to
//!   maintain a consistent physical wavelength relative to
//!   plate sizes — the defaults target the milestone's
//!   32²–64² validation grids.
//! - **`noise::Fbm<Perlin>` is not periodic on the torus**: any
//!   cell on a continental plate that wraps across the domain edge
//!   will see a noise discontinuity. Continental plates rarely
//!   wrap fully (they're surrounded by oceanic neighbours by
//!   construction in `Voronoi` mode), so this is not a visual
//!   issue in practice. Documented as expected variability.
//! - **`fbm_seed` is independent of the Voronoï seed**: lets the
//!   user explore noise textures without redrawing plate
//!   geometry.
//!
//! ## Determinism
//!
//! `noise::Fbm<Perlin>` is deterministic in its `u32` seed. The
//! `u64` `fbm_seed` is cast to `u32` (preserving entropy in the
//! low 32 bits — adequate for a noise seed); same `(plates,
//! continental/oceanic/profile, fbm_*) → byte-identical output.
//! Verified by [`tests::seed_reproducible`].

use noise::{Fbm, MultiFractal, NoiseFn, Perlin};

use super::super::boundaries::PlateType;
use super::super::field::Field2D;
use super::PlateInitData;
use super::radial_profile::{self, ProfileShape};

/// Default `fbm_amplitude` — the noise perturbation magnitude in
/// the same units as `S̃`. Calibrated empirically (Phase 6) so the
/// FBM contribution `σ_fbm_isolated ≥ 0.040` on the milestone's
/// validation grids. `noise::Fbm<Perlin>` is auto-normalised and
/// produces `σ ≈ 0.27 × amplitude` (not the `≈ amplitude / 1.5–2.0`
/// the issue D2 estimate assumed); `amplitude = 0.20` at
/// `scale = 0.10` gives `σ_fbm ≈ 0.05` with comfortable margin
/// across plate sizes. Stays well within the Phase 5 UI clamp
/// `[0.0, 0.40]` so a continental cell at the interior peak
/// (`S̃_radial ≈ 0.95`) cannot dip below the `0.5` continental
/// threshold via FBM.
pub const FBM_AMPLITUDE_DEFAULT: f64 = 0.20;

/// Default `fbm_octaves` — 4 layers of self-similar detail.
pub const FBM_OCTAVES_DEFAULT: u8 = 4;

/// Default `fbm_persistence` — amplitude ratio between successive
/// octaves (0.5 = halve).
pub const FBM_PERSISTENCE_DEFAULT: f64 = 0.5;

/// Default `fbm_lacunarity` — frequency ratio between successive
/// octaves (2.0 = double).
pub const FBM_LACUNARITY_DEFAULT: f64 = 2.0;

/// Default `fbm_scale` — largest feature size in fractions of the
/// domain. `0.10 ⇒ feature size ≈ 1/10 of the domain ≈ 6 cells on
/// 64² grids`, smaller than typical plate `L_plate ≈ 10–15 cells`,
/// so the FBM actually oscillates intra-plate. Calibration history
/// in the module docstring (Phase 6 amendment).
pub const FBM_SCALE_DEFAULT: f64 = 0.10;

/// Default `fbm_seed` — distinct from the Voronoï seed channel so
/// the user can vary the noise texture without redrawing the plate
/// geometry. Hex spelt as "FBA5EED" ≈ "FBA SEED".
pub const FBM_SEED_DEFAULT: u64 = 0x0FBA_5EED;

/// Step 13.5 — strict upper bound on `S̃` for oceanic cells under
/// FBM perturbation. Set to `0.49` (not `0.5`) to leave a defensive
/// margin for floating-point edge cases at the continental
/// classification threshold (`0.5`). Ensures oceanic cells cannot
/// cross to continental classification via FBM noise — volcanic
/// islands are out of scope for Step 13.5 (Step 13.6 if pursued).
pub const OCEANIC_CLAMP_MAX: f64 = 0.49;

/// Step 13.5 — default `fbm_amplitude_oceanic`. Calibrated
/// empirically by Phase 4 sweep (see `tests/v2_step13_5_acceptance.rs::
/// fbm_oceanic_calibration_probe`) on the `single_continent` 64²
/// preset, mid-target band of the issue's
/// `σ_fbm_oceanic_isolated ∈ [0.02, 0.08]`:
///
/// - `noise::Fbm<Perlin>` is auto-normalised, so `σ ≈ 0.27 ×
///   amplitude` (same finding as Step 13's continental side).
/// - At `amplitude = 0.15`, `σ_fbm_oceanic_isolated ≈ 0.040` —
///   sits in the middle of the target band with margin on both
///   sides.
/// - `max(S̃_oceanic) ≈ 0.31` — well under the strict
///   `OCEANIC_CLAMP_MAX = 0.49` threshold, no clipping at sane
///   amplitudes (the sweep showed clip-fraction = 0 % across the
///   full `amp × scale` grid up to `amp = 0.25`).
/// - σ is insensitive to `fbm_scale_oceanic` over `[0.05, 0.20]`
///   (variation < 5 %), so the `None` default for
///   `fbm_scale_oceanic` (= reuse continental `fbm_scale`) is
///   justified empirically rather than just by parsimony.
///
/// Visually: this amplitude produces a measurable bathymetric
/// texture without overwhelming the continental signature
/// (Phase 4 sanity-check patchwork). UI clamps the slider to
/// `[0.0, 0.40]`; the algorithm itself only enforces
/// `OCEANIC_CLAMP_MAX` on the perturbed value.
pub const FBM_AMPLITUDE_OCEANIC_DEFAULT: f64 = 0.15;

/// Step 13.5 — XOR magic constant for deriving `fbm_seed_oceanic`
/// from `fbm_seed` when the user does not supply an explicit
/// oceanic seed. `0xC0FFEE` — opaque enough to avoid trivial
/// correlation between the two seeds, memorable enough for the
/// reader to spot it in a diff. Documented in the seed
/// independence test (`oceanic_fbm_seed_default_derivation`).
pub const FBM_SEED_OCEANIC_XOR_MAGIC: u64 = 0xC0_FFEE;

/// Build the S̃ field for `InitMode::RadialProfileWithFBM`. See
/// module docstring for the algorithm.
///
/// Step 13.5 extends the function with optional FBM on **oceanic**
/// cells, opt-in via `apply_fbm_to_oceanic`. When `false` (default),
/// the function is **bit-identical** to its Step 13 form by
/// short-circuit — the oceanic FBM block is entirely skipped, no
/// second `Fbm<Perlin>` instance is constructed, and oceanic cells
/// retain `oceanic_value` from `radial_profile::build`. When `true`,
/// oceanic cells receive
/// `clamp(oceanic_value + fbm_amplitude_oceanic · fbm_oceanic.get(x, y),
/// 0, OCEANIC_CLAMP_MAX)` — with the strict upper bound `0.49`
/// preventing threshold-crossing to continental classification (D7).
///
/// `fbm_seed_oceanic = None` derives from `fbm_seed XOR
/// FBM_SEED_OCEANIC_XOR_MAGIC` (D3); `fbm_scale_oceanic = None`
/// reuses `fbm_scale` (D4). Either may be supplied explicitly.
///
/// # Panics
///
/// - `continental_value`, `oceanic_value`, or any `fbm_*` parameter
///   is not finite.
/// - `fbm_octaves == 0`.
/// - `fbm_scale <= 0`.
/// - `apply_fbm_to_oceanic = true` with non-finite
///   `fbm_amplitude_oceanic` or non-positive resolved `scale_oceanic`.
#[allow(clippy::too_many_arguments)]
pub fn build(
    nx: usize,
    ny: usize,
    p: &PlateInitData<'_>,
    continental_value: f64,
    oceanic_value: f64,
    profile_shape: ProfileShape,
    fbm_amplitude: f64,
    fbm_octaves: u8,
    fbm_persistence: f64,
    fbm_lacunarity: f64,
    fbm_scale: f64,
    fbm_seed: u64,
    // Step 13.5 — oceanic FBM extension (opt-in). Defaults at the
    // call site preserve Step 13 bit-identical behaviour by
    // short-circuiting the oceanic block when `apply_fbm_to_oceanic`
    // is false — no second FBM instance is constructed, no oceanic
    // cell is touched.
    apply_fbm_to_oceanic: bool,
    fbm_amplitude_oceanic: f64,
    fbm_scale_oceanic: Option<f64>,
    fbm_seed_oceanic: Option<u64>,
) -> Field2D {
    assert!(
        fbm_amplitude.is_finite()
            && fbm_persistence.is_finite()
            && fbm_lacunarity.is_finite()
            && fbm_scale.is_finite(),
        "RadialProfileWithFBM requires finite fbm_* parameters"
    );
    assert!(fbm_octaves >= 1, "fbm_octaves must be ≥ 1, got {fbm_octaves}");
    assert!(fbm_scale > 0.0, "fbm_scale must be > 0, got {fbm_scale}");

    // Phase 2 base — radial profile on continental cells, flat
    // oceanic_value on oceanic cells.
    let mut s = radial_profile::build(nx, ny, p, continental_value, oceanic_value, profile_shape);

    // FBM generator. `Fbm::<Perlin>::new(u32)` seeds the underlying
    // Perlin sources; subsequent `set_*` calls configure the multi-
    // fractal stack. `set_frequency(1.0)` is explicit for clarity —
    // we drive the spatial scale by dividing the sample coordinates
    // by `fbm_scale` instead.
    let fbm = Fbm::<Perlin>::new(fbm_seed as u32)
        .set_octaves(fbm_octaves as usize)
        .set_persistence(fbm_persistence)
        .set_lacunarity(fbm_lacunarity)
        .set_frequency(1.0);

    let nx_f = nx as f64;
    let ny_f = ny as f64;
    let inv_scale = 1.0 / fbm_scale;

    for j in 0..ny {
        for i in 0..nx {
            // FBM only on continental cells (D2). Oceanic cells
            // already hold oceanic_value from radial_profile::build.
            if matches!(p.plate_type.get(i, j), PlateType::Oceanic) {
                continue;
            }
            // Normalised cell-centre coordinates in [0, 1) divided
            // by `fbm_scale` so the largest noise feature has
            // wavelength ≈ `fbm_scale × domain_size` cells.
            let x = ((i as f64) + 0.5) / nx_f * inv_scale;
            let y = ((j as f64) + 0.5) / ny_f * inv_scale;
            let n = fbm.get([x, y]);
            let perturbed = s.get(i, j) + fbm_amplitude * n;
            // Clamp to physical range. With the default amplitude
            // (0.10) and `S̃_radial ∈ [0.20, 0.95]`, perturbed lies
            // in `[0.10, 1.05]` ≈ entirely inside [0, 1] modulo a
            // tiny upper-edge clip. Larger amplitudes may bite the
            // 0 floor at the boundary side — see vigilance note.
            s.set(i, j, perturbed.clamp(0.0, 1.0));
        }
    }

    // Step 13.5 — oceanic FBM extension. Strictly opt-in. When the
    // flag is false the entire block below is skipped — no second
    // `Fbm<Perlin>` instance, no oceanic cell modification — so the
    // function output is bit-identical to its Step 13 form.
    if apply_fbm_to_oceanic {
        // D3 — derive oceanic seed: explicit `Some(seed)` wins,
        // otherwise XOR the continental seed with the magic
        // constant for reasonable independence.
        let seed_oceanic = fbm_seed_oceanic.unwrap_or(fbm_seed ^ FBM_SEED_OCEANIC_XOR_MAGIC);
        // D4 — derive oceanic scale: explicit `Some(scale)` wins,
        // otherwise reuse the continental scale.
        let scale_oceanic = fbm_scale_oceanic.unwrap_or(fbm_scale);
        assert!(
            fbm_amplitude_oceanic.is_finite(),
            "fbm_amplitude_oceanic must be finite, got {fbm_amplitude_oceanic}"
        );
        assert!(
            scale_oceanic.is_finite() && scale_oceanic > 0.0,
            "resolved fbm_scale_oceanic must be a positive finite scalar, got {scale_oceanic}"
        );

        let fbm_oceanic = Fbm::<Perlin>::new(seed_oceanic as u32)
            .set_octaves(fbm_octaves as usize)
            .set_persistence(fbm_persistence)
            .set_lacunarity(fbm_lacunarity)
            .set_frequency(1.0);
        let inv_scale_oceanic = 1.0 / scale_oceanic;

        for j in 0..ny {
            for i in 0..nx {
                // Mirror of the continental loop: only oceanic
                // cells are touched here.
                if !matches!(p.plate_type.get(i, j), PlateType::Oceanic) {
                    continue;
                }
                let x = ((i as f64) + 0.5) / nx_f * inv_scale_oceanic;
                let y = ((j as f64) + 0.5) / ny_f * inv_scale_oceanic;
                let n = fbm_oceanic.get([x, y]);
                let perturbed = oceanic_value + fbm_amplitude_oceanic * n;
                // D7 — strict threshold protection. Upper bound
                // `OCEANIC_CLAMP_MAX = 0.49` keeps oceanic cells
                // strictly oceanic by classification regardless of
                // the FBM peak amplitude; volcanic islands (cells
                // crossing 0.5) are out of scope for Step 13.5.
                s.set(i, j, perturbed.clamp(0.0, OCEANIC_CLAMP_MAX));
            }
        }
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

    fn build_plates(
        nx: usize,
        ny: usize,
        seed: u64,
        num_plates: usize,
        continental_ratio: f64,
    ) -> crate::tectonics_v2::voronoi::VoronoiPlates {
        let cfg = VoronoiConfig { num_plates, continental_ratio };
        generate_voronoi(nx, ny, &cfg, seed)
    }

    fn make_init_data<'a>(
        plates: &'a crate::tectonics_v2::voronoi::VoronoiPlates,
    ) -> PlateInitData<'a> {
        PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        }
    }

    fn defaults() -> (f64, u8, f64, f64, f64, u64) {
        (
            FBM_AMPLITUDE_DEFAULT,
            FBM_OCTAVES_DEFAULT,
            FBM_PERSISTENCE_DEFAULT,
            FBM_LACUNARITY_DEFAULT,
            FBM_SCALE_DEFAULT,
            FBM_SEED_DEFAULT,
        )
    }

    /// Acceptance: oceanic cells are unaffected by the FBM noise.
    /// Run with two distinct FBM seeds; oceanic cells must be
    /// byte-identical between runs (= `oceanic_value`).
    #[test]
    fn continental_only() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, _) = defaults();

        let s_seed_1 = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            1,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );
        let s_seed_2 = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            999,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );

        let mut count_oceanic = 0;
        for j in 0..ny {
            for i in 0..nx {
                if matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    assert_eq!(
                        s_seed_1.get(i, j),
                        0.20,
                        "oceanic cell ({},{}) FBM seed 1 expected 0.20, got {}",
                        i,
                        j,
                        s_seed_1.get(i, j)
                    );
                    assert_eq!(
                        s_seed_2.get(i, j),
                        0.20,
                        "oceanic cell ({},{}) FBM seed 2 expected 0.20, got {}",
                        i,
                        j,
                        s_seed_2.get(i, j)
                    );
                    count_oceanic += 1;
                }
            }
        }
        assert!(count_oceanic > 0, "no oceanic cells found in 64² × 8 plates @ 40% continental");
    }

    /// Acceptance: variation introduced by FBM is bounded by
    /// `±amplitude` (modulo clamping to `[0, 1]`). We compare each
    /// cell to the FBM-free baseline produced by Phase 2 and assert
    /// `|S̃_fbm - S̃_radial| ≤ amplitude + ε` everywhere.
    #[test]
    fn amplitude_bounded() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (amplitude, o, persist, lac, scale, seed) = defaults();

        let s_radial = radial_profile::build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);
        let s_fbm = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            amplitude,
            o,
            persist,
            lac,
            scale,
            seed,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );

        for j in 0..ny {
            for i in 0..nx {
                let r = s_radial.get(i, j);
                let f = s_fbm.get(i, j);
                let delta = (f - r).abs();
                // Account for the [0, 1] clamp absorbing some
                // overshoot when r + amplitude·n would exceed the
                // bounds: in that case |f - r| < |amplitude·n|, so
                // the bound still holds.
                assert!(
                    delta <= amplitude + 1e-12,
                    "cell ({},{}) FBM perturbation {} exceeds amplitude {}",
                    i,
                    j,
                    delta,
                    amplitude
                );
            }
        }
    }

    /// Acceptance: same `fbm_seed` produces byte-identical output.
    #[test]
    fn seed_reproducible() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, seed) = defaults();

        let s_a = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );
        let s_b = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );
        assert_eq!(s_a.data(), s_b.data());
    }

    /// Acceptance #4: every cell value lies in `[0, 1]` (clamped).
    /// Stress with a large amplitude that would otherwise push
    /// values out of range.
    #[test]
    fn clamped() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);

        // Stress amplitude — would push values to ~[-0.05, 1.45]
        // without clamping (since radial range is [0.20, 0.95]).
        let s = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            0.50,
            FBM_OCTAVES_DEFAULT,
            FBM_PERSISTENCE_DEFAULT,
            FBM_LACUNARITY_DEFAULT,
            FBM_SCALE_DEFAULT,
            FBM_SEED_DEFAULT,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );

        for j in 0..ny {
            for i in 0..nx {
                let v = s.get(i, j);
                assert!(
                    v.is_finite() && (0.0..=1.0).contains(&v),
                    "cell ({},{}) value {} out of [0, 1]",
                    i,
                    j,
                    v
                );
            }
        }
    }

    /// Robustness: with `amplitude = 0.0`, the FBM mode reduces to
    /// the Phase 2 radial profile byte-for-byte.
    #[test]
    fn zero_amplitude_equals_radial() {
        let nx = 32;
        let ny = 32;
        let plates = build_plates(nx, ny, 42, 6, 0.4);
        let p = make_init_data(&plates);
        let s_radial = radial_profile::build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);
        let s_fbm = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            0.0,
            FBM_OCTAVES_DEFAULT,
            FBM_PERSISTENCE_DEFAULT,
            FBM_LACUNARITY_DEFAULT,
            FBM_SCALE_DEFAULT,
            FBM_SEED_DEFAULT,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );
        assert_eq!(s_radial.data(), s_fbm.data());
    }

    /// Step 13.5 acceptance #1 — `apply_fbm_to_oceanic = false`
    /// produces a field byte-identical to the equivalent Step 13
    /// build (no oceanic FBM contribution). The Step 13.5 contract
    /// rests on this short-circuit being structural, not numerical.
    #[test]
    fn oceanic_fbm_disabled_preserves_step13() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, seed) = defaults();

        // Step 13-equivalent build: oceanic FBM disabled, the
        // amplitude/scale/seed defaults are written but never read.
        let s_step13 = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );
        // Same shape, but with non-trivial oceanic params (any
        // non-default value the user might pass): the disabled
        // flag must short-circuit before they're read, so the
        // output stays byte-identical to the disabled-default
        // build above.
        let s_with_unused_params = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            false,        // <-- the gate
            0.42,         // bogus amplitude — must not be applied
            Some(0.07),   // bogus scale       — idem
            Some(0xDEAD), // bogus seed        — idem
        );
        assert_eq!(
            s_step13.data(),
            s_with_unused_params.data(),
            "apply_fbm_to_oceanic = false must short-circuit before \
             reading any of the oceanic FBM parameters"
        );
    }

    /// Step 13.5 acceptance #2 — with the flag enabled and any
    /// non-zero `fbm_amplitude_oceanic`, oceanic cells gain a
    /// measurable variance (var > 0). Confirms the FBM is
    /// actually sampled and added, not silently dropped.
    #[test]
    fn oceanic_fbm_enabled_varies() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, seed) = defaults();

        let s = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            true,
            0.10,
            None,
            None,
        );

        let mut sum = 0.0_f64;
        let mut count = 0usize;
        for j in 0..ny {
            for i in 0..nx {
                if matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    sum += s.get(i, j);
                    count += 1;
                }
            }
        }
        assert!(count > 0, "no oceanic cells found — preset assumption failed");
        let mean = sum / count as f64;
        let mut var = 0.0_f64;
        for j in 0..ny {
            for i in 0..nx {
                if matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    let d = s.get(i, j) - mean;
                    var += d * d;
                }
            }
        }
        var /= count as f64;
        assert!(
            var > 0.0,
            "oceanic FBM enabled but oceanic variance is zero (count={}, mean={:.6})",
            count,
            mean
        );
    }

    /// Step 13.5 acceptance #3 — across an amplitude sweep
    /// `{0.05, 0.10, 0.20, 0.30, 0.40}`, every oceanic cell
    /// stays at `S̃ ≤ OCEANIC_CLAMP_MAX = 0.49`. The strict
    /// upper bound prevents threshold-crossing to continental
    /// classification regardless of FBM peak amplitude (D7).
    #[test]
    fn oceanic_fbm_no_threshold_crossing() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, seed) = defaults();

        for amp_oceanic in [0.05, 0.10, 0.20, 0.30, 0.40] {
            let s = build(
                nx,
                ny,
                &p,
                0.95,
                0.20,
                ProfileShape::Smoothstep,
                a,
                o,
                persist,
                lac,
                scale,
                seed,
                true,
                amp_oceanic,
                None,
                None,
            );
            for j in 0..ny {
                for i in 0..nx {
                    if !matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                        continue;
                    }
                    let v = s.get(i, j);
                    assert!(
                        v <= OCEANIC_CLAMP_MAX + 1e-15,
                        "oceanic cell ({},{}) at amplitude={} crossed clamp: {} > {}",
                        i,
                        j,
                        amp_oceanic,
                        v,
                        OCEANIC_CLAMP_MAX
                    );
                    assert!(
                        v >= 0.0,
                        "oceanic cell ({},{}) at amplitude={} below zero: {}",
                        i,
                        j,
                        amp_oceanic,
                        v
                    );
                }
            }
        }
    }

    /// Step 13.5 acceptance #4 — different `fbm_seed_oceanic`
    /// values produce different oceanic fields, while the
    /// continental field stays byte-identical (continental seed
    /// unchanged). Verifies the two `Fbm<Perlin>` instances are
    /// independent: changing one does not perturb the other.
    #[test]
    fn oceanic_fbm_seed_independence() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, seed) = defaults();

        let s_a = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            true,
            0.10,
            None,
            Some(1),
        );
        let s_b = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            true,
            0.10,
            None,
            Some(999),
        );

        let mut oceanic_diff_count = 0usize;
        for j in 0..ny {
            for i in 0..nx {
                let pt = plates.plate_type.get(i, j);
                let va = s_a.get(i, j);
                let vb = s_b.get(i, j);
                if matches!(pt, PlateType::Oceanic) && va != vb {
                    oceanic_diff_count += 1;
                }
                if matches!(pt, PlateType::Continental) {
                    assert_eq!(
                        va, vb,
                        "continental cell ({},{}) differs between oceanic-seed \
                         variants — continental field must be insulated from \
                         oceanic seed: {} vs {}",
                        i, j, va, vb
                    );
                }
            }
        }
        assert!(
            oceanic_diff_count > 0,
            "no oceanic cells differ between distinct fbm_seed_oceanic values \
             — the oceanic FBM is not actually consuming the seed"
        );
    }

    /// Step 13.5 acceptance #5 — `fbm_seed_oceanic = None`
    /// derives the oceanic seed from `fbm_seed XOR
    /// FBM_SEED_OCEANIC_XOR_MAGIC`, so the explicit-derivation
    /// build and the default-derivation build produce identical
    /// fields.
    #[test]
    fn oceanic_fbm_seed_default_derivation() {
        let nx = 32;
        let ny = 32;
        let plates = build_plates(nx, ny, 42, 6, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, seed) = defaults();

        let derived = seed ^ FBM_SEED_OCEANIC_XOR_MAGIC;

        let s_default = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            true,
            0.10,
            None,
            None, // <-- derive from fbm_seed
        );
        let s_explicit = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            true,
            0.10,
            None,
            Some(derived), // <-- explicit derivation
        );
        assert_eq!(
            s_default.data(),
            s_explicit.data(),
            "fbm_seed_oceanic = None must derive from fbm_seed XOR 0x{:X}",
            FBM_SEED_OCEANIC_XOR_MAGIC
        );
    }

    /// Step 13.5 acceptance #6 — different `fbm_scale_oceanic`
    /// values produce different oceanic spectral content (the
    /// fields differ), confirming the scale parameter is
    /// honoured rather than silently shadowed by the continental
    /// scale. We use distinct fbm seeds to avoid coincidental
    /// agreement on a degenerate sample.
    #[test]
    fn oceanic_fbm_scale_independence() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, seed) = defaults();

        let s_short_wave = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            true,
            0.10,
            Some(0.05),
            Some(1),
        );
        let s_long_wave = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            seed,
            true,
            0.10,
            Some(0.20),
            Some(1),
        );

        let mut oceanic_diff_count = 0usize;
        for j in 0..ny {
            for i in 0..nx {
                if !matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    continue;
                }
                if s_short_wave.get(i, j) != s_long_wave.get(i, j) {
                    oceanic_diff_count += 1;
                }
            }
        }
        assert!(
            oceanic_diff_count > 0,
            "no oceanic cells differ between fbm_scale_oceanic = 0.05 and 0.20 — \
             the scale parameter is not being applied"
        );
    }

    /// Sanity: distinct seeds produce distinct fields (FBM is
    /// actually being sampled, not silently zeroed).
    #[test]
    fn different_seeds_differ() {
        let nx = 32;
        let ny = 32;
        let plates = build_plates(nx, ny, 42, 6, 0.4);
        let p = make_init_data(&plates);
        let (a, o, persist, lac, scale, _) = defaults();
        let s_1 = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            1,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );
        let s_2 = build(
            nx,
            ny,
            &p,
            0.95,
            0.20,
            ProfileShape::Smoothstep,
            a,
            o,
            persist,
            lac,
            scale,
            999,
            false,
            FBM_AMPLITUDE_OCEANIC_DEFAULT,
            None,
            None,
        );
        assert_ne!(s_1.data(), s_2.data());
    }
}
