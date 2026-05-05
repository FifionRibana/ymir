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

/// Build the S̃ field for `InitMode::RadialProfileWithFBM`. See
/// module docstring for the algorithm.
///
/// # Panics
///
/// - `continental_value`, `oceanic_value`, or any `fbm_*` parameter
///   is not finite.
/// - `fbm_octaves == 0`.
/// - `fbm_scale <= 0`.
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
    let mut s = radial_profile::build(
        nx,
        ny,
        p,
        continental_value,
        oceanic_value,
        profile_shape,
    );

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
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

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
        assert!(
            count_oceanic > 0,
            "no oceanic cells found in 64² × 8 plates @ 40% continental"
        );
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

        let s_radial =
            radial_profile::build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);
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
        let s_radial =
            radial_profile::build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);
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
        );
        assert_eq!(s_radial.data(), s_fbm.data());
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
        );
        assert_ne!(s_1.data(), s_2.data());
    }
}
