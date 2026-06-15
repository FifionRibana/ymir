//! Production C1 → anisotropic-FBM upscale wiring (Issue #147
//! FOLLOWUPS-#6 outcome).
//!
//! ## The robustness contract (enforced by the type, not just documented)
//!
//! The S̃ thickness field is mesh-NON-convergent in production
//! (downsample-correlation r~0.51 across 64²/128²/256²/512² — Issue
//! #147). The upscale is nonetheless ROBUST to this — 64²+upscale and
//! 256²+upscale produce the SAME world (structure r~0.90, measured in
//! `stage_upscale_robustness.md`) — **only because it reads the
//! laundered production ALTITUDE** (isostasy + Stein-Stein, which
//! converges at alt r~0.88), NOT the raw S̃ gradient (r~0.51). The
//! non-convergence is washed out by isostasy + Stein-Stein before the
//! FBM orientation ever sees it (same absorption that hides the
//! Davis-Suppe mass swing in production).
//!
//! Feeding raw S̃ to the upscale would silently break this: the FBM
//! would orient on a non-convergent gradient, 64²+FBM and 256²+FBM
//! would DIVERGE (invisible at a single resolution; only a 2-resolution
//! comparison reveals it), and it would reopen FOLLOWUPS #6 (the
//! deferred advection-scheme milestone).
//!
//! [`upscale_from_c1`] makes that illegal state **unrepresentable**: it
//! takes a [`C1State`] and builds the altitude INTERNALLY. There is no
//! S̃/heightmap parameter a caller could substitute. The contract is
//! enforced by the signature; the inline comment is belt-and-braces.
//!
//! The regression test
//! `c1_closure_morphology::upscale_from_c1_structure_converges`
//! asserts the precondition (64²/256² structure convergence) so a
//! future change that breaks it is caught, not just warned against.

use crate::erosion::hydraulic::run_erosion;
use crate::grid::GridF32;
use crate::seed::WorldSeed;
use crate::tectonics::isostasy::{compute_isostasy_c1, IsostasyConfig};
use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::Field2D;
use crate::terrain::upscale::{upscale_with_fbm, FbmUpscaleConfig, UpscaleResult};

use super::closures::oceanic_bathymetry::params::SteinSteinParams;
use super::closures::oceanic_bathymetry::source_term::apply_stein_stein_bathymetry;
use super::state::C1State;

/// The C1 PRODUCTION altitude (Architecture C): `compute_isostasy`
/// followed by the Stein-Stein bathymetry re-apply on oceanic cells.
///
/// **SINGLE SOURCE OF TRUTH** for the production altitude, shared by
/// the viz render (`bridge::c1` / `c1_viz::derive_altitude_field`, which
/// reconstructs these fields from a snapshot and calls this) AND the
/// upscale input ([`upscale_from_c1`]). Factoring it here (Issue #147
/// #6) prevents the upscaled HD terrain from silently diverging from
/// the rendered production altitude: the two agree BY CONSTRUCTION, not
/// by two copies of the same `isostasy + Stein-Stein` sequence drifting
/// apart. The regression test guards robustness; this shared function
/// guards render/upscale equality.
pub fn c1_production_altitude(
    s: &Field2D,
    age: &Field2D,
    plate_type: &PlateTypeField,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
) -> GridF32 {
    c1_production_altitude_inner(s, age, plate_type, iso, ss, None)
}

/// #155 B-Jordan — [`c1_production_altitude`] with spatial cratonic density.
/// `craton` (row-major `&[bool]`, the cratonic mask) routes those cells to
/// `compute_isostasy_craton` (denser cratonic crust → worn-shield altitude).
/// `iso.craton_rho_crust == None` → byte-identical to the plain function.
pub fn c1_production_altitude_craton(
    s: &Field2D,
    age: &Field2D,
    plate_type: &PlateTypeField,
    craton: &[bool],
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
) -> GridF32 {
    c1_production_altitude_inner(s, age, plate_type, iso, ss, Some(craton))
}

fn c1_production_altitude_inner(
    s: &Field2D,
    age: &Field2D,
    plate_type: &PlateTypeField,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    craton: Option<&[bool]>,
) -> GridF32 {
    // #155 land-ceiling repair — the C1 production ALWAYS supplies the
    // continental mask so the land ramp tops at the real terrestrial summit,
    // not the phantom oceanic advective spike (which Stein-Stein overwrites
    // below). Built from plate_type (row-major, length nx·ny). v2/export do
    // NOT go through this path → they keep `compute_isostasy` (no mask,
    // byte-identical). See `compute_isostasy_c1`.
    let continental: Vec<bool> = (0..plate_type.ny())
        .flat_map(|j| (0..plate_type.nx()).map(move |i| (i, j)))
        .map(|(i, j)| plate_type.get(i, j) == PlateType::Continental)
        .collect();
    let isostasy = compute_isostasy_c1(s, iso, craton, &continental);
    let mut altitude = isostasy.heightmap;

    // #155 Maillon 2 — sea-CENTRE the continental altitude to metres /
    // depth_scale (sea → 0), matching the Stein-Stein oceanic convention
    // (which already writes −depth/depth_scale_m, sea at 0). The isostasy
    // heightmap is [0,1] with continental sea at `sea_norm` (≈0.111); converting
    // it to a sea-centred metres/depth_scale field means the downstream FIXED
    // re-offset `(alt + 1.13)/2.26` (shared by the viz render AND
    // `upscale_from_c1`) maps sea → 0.5 EXACTLY — no more 0.5-vs-0.549 mismatch
    // between the two stages — and the whole field (land + ocean) shares ONE
    // linear norm→m scale (see `c1_altitude_norm_to_metres`). `peak_altitude_m`
    // (land summit) and `max_depth_m` (continental shelf) are the metre anchors
    // from the isostasy metadata; `depth_scale_m` is the shared unit. Oceanic
    // cells are overwritten by Stein-Stein just below, so converting them here
    // is harmless (they re-enter in the same metres/depth_scale convention).
    let sea_norm = isostasy.sea_level_normalized;
    let peak_m = isostasy.peak_altitude_m;
    let shelf_m = isostasy.max_depth_m;
    let scale = ss.depth_scale_m as f32;
    for v in altitude.data.iter_mut() {
        let n = *v;
        let metres = if n >= sea_norm {
            (n - sea_norm) / (1.0 - sea_norm).max(1e-10) * peak_m
        } else {
            -(sea_norm - n) / sea_norm.max(1e-10) * shelf_m
        };
        *v = metres / scale;
    }

    // #151 F3 fix — despike the age before Stein-Stein. The flux-form age
    // advection piles age up at convergent plate boundaries (sparse cells
    // with ~1000× the background age — the registered density-vs-Lagrangian
    // artefact). Stein-Stein turns those age spikes into the deepest cells,
    // rendering as dark dotted lines in the ocean (F3). A 3×3 median kills
    // the sparse spikes while preserving smooth age structure (verified: no
    // legitimate age gradient is flattened — the C1 model has none without
    // seafloor spreading). This is a RENDER-side band-aid: Stein-Stein is the
    // only consumer of age on the altitude path, so this fixes F3 fully for
    // both the viz render and the upscale. The internal advected age field is
    // unchanged (the root cause — age advection — is a deferred deep fix).
    let age_despiked = median_3x3(age);
    apply_stein_stein_bathymetry(&mut altitude, &age_despiked, plate_type, ss);
    altitude
}

/// 3×3 median filter (clamped edges). Removes sparse spikes (the age
/// pile-up cells) while preserving smooth structure.
fn median_3x3(field: &Field2D) -> Field2D {
    let nx = field.nx();
    let ny = field.ny();
    let mut out = field.clone();
    let mut nb: Vec<f64> = Vec::with_capacity(9);
    for j in 0..ny {
        for i in 0..nx {
            nb.clear();
            for dj in -1i32..=1 {
                for di in -1i32..=1 {
                    let ni = i as i32 + di;
                    let nj = j as i32 + dj;
                    if ni >= 0 && nj >= 0 && (ni as usize) < nx && (nj as usize) < ny {
                        nb.push(field.get(ni as usize, nj as usize));
                    }
                }
            }
            nb.sort_by(|a, b| a.partial_cmp(b).unwrap());
            out.set(i, j, nb[nb.len() / 2]);
        }
    }
    out
}

/// Fixed altitude→`[0,1]` normalisation half-range (sea level 0.0 maps
/// to 0.5). This MUST be a constant, NOT a per-call data range:
/// normalising 64² and 256² by their own min/max would make them
/// incomparable and break the measured cross-resolution robustness.
/// Value matches the production gallery palette half-range.
const ALTITUDE_NORM_HALF_RANGE: f32 = 1.13;

/// #155 Maillon 2 — **THE unified C1 vertical-scale contract.** Convert an HD
/// heightmap normalised value (sea = 0.5, the `upscale_from_c1` / viz-render
/// convention) to METRES. A SINGLE linear scale across land AND ocean — no
/// piecewise seam, no per-regime slope: the consumer reads metres without
/// knowing whether it is on land or sea.
///
/// ```text
///     metres = (norm − 0.5) · 2 · ALTITUDE_NORM_HALF_RANGE · depth_scale_m
/// ```
///
/// With the defaults (`1.13`, `5000`): `metres = (norm − 0.5) · 11300`. Sea
/// (0.5) → 0 m; `norm = 1.0` → +5650 m; `norm = 0.0` → −5650 m. Land occupies
/// `[0.5, ~0.85]` (the upper land bound is `max_elevation_m = 4000 m`, set by
/// isostasy); ocean occupies `[0, 0.5]` (down to the Stein-Stein asymptote
/// ~5651 m). The headroom `~0.85 → 1.0` is RESERVED for the separate
/// critical-wedge high-mountain chantier — the scale tells the truth, and the
/// gap it reveals is the diagnostic of the relief not yet produced.
///
/// This is the contract downstream consumers (rivers, biome, climate) read to
/// reason in defined metres. Coherent BY CONSTRUCTION with the sea-centred
/// production altitude (`c1_production_altitude`) — see its body.
pub fn c1_altitude_norm_to_metres(norm: f32, ss: &SteinSteinParams) -> f32 {
    (norm - 0.5) * 2.0 * ALTITUDE_NORM_HALF_RANGE * ss.depth_scale_m as f32
}

/// Inverse of [`c1_altitude_norm_to_metres`]: metres → normalised value.
pub fn c1_metres_to_altitude_norm(metres: f32, ss: &SteinSteinParams) -> f32 {
    metres / (2.0 * ALTITUDE_NORM_HALF_RANGE * ss.depth_scale_m as f32) + 0.5
}

/// #155 — **THE C1 HORIZONTAL scale (the coordinate contract's other half).**
/// The km side-length the unit (`1 × 1`) tectonic domain represents. This is the
/// pendant of the vertical `c1_altitude_norm_to_metres`: the metadata.json
/// coordinate contract (§9.3) lists meters-per-pixel ALONGSIDE the elevation
/// scale — both are ONE coordinate system, and a defined horizontal scale is the
/// shared prerequisite for every metric consumer (rivers' drainage-area
/// thresholds, biome zones, climate orographic distance, village spacing).
///
/// **Derivation (anchored, not arbitrary).** TDD §11 "Implicit physical scales"
/// states the C1 domain is regional (~1000-5000 km) with dx ≈ "~2 km at 512²".
/// Pinning the lower-end TDD anchor `2.0 km/cell × 512 = 1024 km` makes that
/// implicit statement explicit. dx then = 16 km @64² (tectonic), 0.5 km @2048²
/// (HD). Consistent with the §2.2 gameplay scale (a dense playable region —
/// cités, riverside villages, navigable valleys — i.e. "large region", not a
/// stretched whole continent). **Revisable product choice**: a "whole continent"
/// intent would set ~3000 km; change this ONE constant and every consumer
/// follows (nothing else encodes the horizontal scale).
pub const C1_DOMAIN_KM: f32 = 1024.0;

/// Kilometres per cell at the given grid resolution (the `1×1` domain is
/// [`C1_DOMAIN_KM`] on a side). Resolution-INdependent contract: the domain km
/// is fixed; km/cell = `C1_DOMAIN_KM / grid_size`. E.g. 2.0 km @512², 0.5 km
/// @2048².
pub fn c1_km_per_cell(grid_size: usize) -> f32 {
    C1_DOMAIN_KM / grid_size as f32
}

/// Cell area in km² at the given grid resolution — the unit for
/// resolution-independent drainage-area thresholds (flow accumulation × this =
/// upstream area in km², which does NOT change with grid_size for the same
/// physical catchment).
pub fn c1_cell_area_km2(grid_size: usize) -> f32 {
    let s = c1_km_per_cell(grid_size);
    s * s
}

/// Upscale a C1 simulation state to a detailed heightmap via the
/// anisotropic-FBM upscale, **through the robustness-preserving
/// altitude path** (Issue #147 FOLLOWUPS-#6 — see module docs).
///
/// Takes the [`C1State`] (not an S̃ field, not a heightmap) BY DESIGN:
/// the altitude is built internally so no caller can feed the
/// non-convergent raw S̃ to the FBM orientation.
///
/// Pipeline (identical to the measured-robust `stage_upscale_robustness`
/// path): `compute_isostasy(iso)` → `apply_stein_stein_bathymetry(ss)`
/// → fixed `[0,1]` normalisation (sea→0.5) → `upscale_with_fbm`.
pub fn upscale_from_c1(
    state: &C1State,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    seed: &WorldSeed,
    cfg: &FbmUpscaleConfig,
) -> UpscaleResult {
    // CONTRACT (Issue #147 #6): the upscale reads the laundered
    // ALTITUDE (isostasy + Stein-Stein, convergent r~0.88), NOT raw S̃
    // (non-convergent r~0.51). This is what makes the upscale robust to
    // S̃ mesh non-convergence. Built via the SHARED `c1_production_altitude`
    // (same source of truth as the viz render) — do NOT inline a raw-S̃
    // path here; it reopens FOLLOWUPS #6 and silently diverges.
    // #155 B-Jordan: route the cratonic mask so cratonic cells get the
    // dense-crust buoyancy (iso.craton_rho_crust). None → byte-identical.
    let mut coarse = c1_production_altitude_craton(
        &state.s, &state.age, &state.plate_type, state.cratonic_mask.data(), iso, ss,
    );

    // Fixed normalisation to [0,1] (sea 0.0 → 0.5), resolution-independent.
    let half = ALTITUDE_NORM_HALF_RANGE;
    for v in coarse.data.iter_mut() {
        *v = ((*v + half) / (2.0 * half)).clamp(0.0, 1.0);
    }
    let sea_level_normalized = 0.5_f32;

    let result = upscale_with_fbm(&coarse, sea_level_normalized, seed, cfg);

    // #155 méso — HD hydraulic erosion (the dendritic dissection that makes
    // the macro ridge read as credible eroded mountains). Applied AFTER the
    // FBM, ONLY when `cfg.erosion` is Some (the canonical C1 HD product config
    // `FbmUpscaleConfig::c1_hd_production` turns it on; default None →
    // byte-identical). Slope is RECOMPUTED post-erosion (it changed);
    // `sediment` is forwarded (the rivers/lakes hook, not consumed yet).
    let Some(ero) = &cfg.erosion else {
        return result;
    };
    let eroded = run_erosion(&result.heightmap, ero, seed, |_, _, _| true);
    let h = &eroded.heightmap;
    let mut slope = GridF32::new(h.width, h.height, 0.0);
    for j in 0..h.height {
        for i in 0..h.width {
            let (gx, gy) = h.gradient_at(i, j);
            slope.set(i, j, (gx * gx + gy * gy).sqrt());
        }
    }
    UpscaleResult {
        heightmap: eroded.heightmap,
        slope,
        sediment: Some(eroded.sediment),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_c1::init_r7::{init_c1_state_phase_2_r7, Phase2InitParams};

    /// Smoke: the contract path runs and produces a target-sized
    /// heightmap with finite values (no NaN/Inf), at a small target.
    #[test]
    fn upscale_from_c1_smoke() {
        let state = init_c1_state_phase_2_r7(64, 42, &Phase2InitParams::default());
        let cfg = FbmUpscaleConfig { target_size: 128, ..Default::default() };
        let out = upscale_from_c1(
            &state,
            &IsostasyConfig::c1_default(),
            &SteinSteinParams::default(),
            &WorldSeed::new(42),
            &cfg,
        );
        assert_eq!(out.heightmap.width, 128);
        assert_eq!(out.heightmap.height, 128);
        assert!(out.heightmap.data.iter().all(|v| v.is_finite()));
    }

    /// #155 horizontal contract — the TDD §11 anchor (2.0 km/cell @512²) +
    /// resolution-independence (domain km fixed; km/cell = domain/grid).
    #[test]
    fn horizontal_scale_contract() {
        assert_eq!(C1_DOMAIN_KM, 1024.0);
        assert!((c1_km_per_cell(512) - 2.0).abs() < 1e-6, "TDD §11 anchor: 2 km/cell @512²");
        assert!((c1_km_per_cell(64) - 16.0).abs() < 1e-6);
        assert!((c1_km_per_cell(2048) - 0.5).abs() < 1e-6);
        // cell area = (km/cell)²; a fixed physical catchment is resolution-invariant:
        // accumulation(cells) × cell_area_km2 = upstream area in km².
        assert!((c1_cell_area_km2(512) - 4.0).abs() < 1e-6);
        assert!((c1_cell_area_km2(2048) - 0.25).abs() < 1e-6);
    }

    /// #155 vertical contract — anchors + round-trip (Maillon 2).
    #[test]
    fn vertical_scale_contract() {
        let ss = SteinSteinParams::default();
        assert!((c1_altitude_norm_to_metres(0.5, &ss)).abs() < 1e-3, "sea = 0 m");
        assert!((c1_altitude_norm_to_metres(1.0, &ss) - 5650.0).abs() < 5.0);
        assert!((c1_altitude_norm_to_metres(0.0, &ss) + 5650.0).abs() < 5.0);
        let rt = c1_metres_to_altitude_norm(c1_altitude_norm_to_metres(0.73, &ss), &ss);
        assert!((rt - 0.73).abs() < 1e-5, "round-trip");
    }
}
