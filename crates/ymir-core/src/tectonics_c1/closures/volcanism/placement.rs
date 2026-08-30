//! C-2 placement — derive edifice sites from the C1 tectonic state (arcs,
//! hotspot chains, rifts). Pure and deterministic: a function of
//! `(C1State, PlateKinematics, seed, VolcanismConfig)` only, so it can be
//! recomputed identically wherever it is needed (HD injection AND lake typing)
//! without depending on a cache hit.
//!
//! All physical lengths stay in km and are converted to normalized coarse-torus
//! offsets at the point of use, via `torus_km` (the physical span of the full
//! coarse torus). Positions are emitted in normalized coarse coords `[0,1)`.

use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use super::{Edifice, VolcanismConfig, VolcanoKind, VolcanoSetting};
use crate::seed::WorldSeed;
use crate::tectonics_c1::boundary_classification::{
    BoundaryType, classify_boundaries, oc_override_seed_mask,
};
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_v2::boundaries::plate_type::PlateType;

/// **Trench→arc horizontal offset, PROVISIONAL (km).** Syracuse & Abers 2006 give
/// a slab DEPTH beneath the front (72–173 km, mean 105), NOT a horizontal
/// distance; our model has no slab-depth field, so this horizontal offset comes
/// from the secondary arc-trench-gap literature (100–300 km). It is therefore the
/// LEAST certain number in C-2 — kept here as a single named constant, in km,
/// converted at the point of use, so it is one line to change when a better source
/// or a slab field arrives. See ADR 0001 (C-2) for its uncertainty statement.
pub const TRENCH_ARC_OFFSET_KM_DEFAULT: f32 = 150.0;

/// Place all volcanic edifices derived from the tectonic state. `torus_km` is the
/// physical span of the FULL coarse torus (km) — km offsets divide by it to get
/// normalized coarse offsets. Deterministic for a fixed `(state, kin, seed, cfg)`.
#[must_use]
pub fn place_edifices(
    state: &C1State,
    kin: &PlateKinematics,
    seed: &WorldSeed,
    torus_km: f32,
    cfg: &VolcanismConfig,
) -> Vec<Edifice> {
    let mut out = Vec::new();
    if !cfg.enabled {
        return out;
    }
    let nx = state.plate_id.nx();
    let ny = state.plate_id.ny();
    let boundary = classify_boundaries(&state.plate_id, kin);
    let to_uv = |i: usize, j: usize| ((i as f32 + 0.5) / nx as f32, (j as f32 + 0.5) / ny as f32);
    // Minimum-separation thinning in normalized coords (periodic torus).
    let torus_dist2 = |a: (f32, f32), b: (f32, f32)| -> f32 {
        let mut dx = (a.0 - b.0).abs();
        let mut dy = (a.1 - b.1).abs();
        if dx > 0.5 {
            dx = 1.0 - dx;
        }
        if dy > 0.5 {
            dy = 1.0 - dy;
        }
        dx * dx + dy * dy
    };

    // ── Arcs: O-C convergent continental margin, offset INBOARD (away from the
    //    trench) by the trench-arc gap, thinned to `arc_spacing_km`. ──
    let arc_spacing_norm2 = (cfg.arc_spacing_km / torus_km).powi(2);
    let mut arc_kept: Vec<(f32, f32)> = Vec::new();
    for (_, pos) in arc_sites(state, &boundary, torus_km, cfg) {
        if arc_kept.iter().any(|&k| torus_dist2(k, pos) < arc_spacing_norm2) {
            continue;
        }
        arc_kept.push(pos);
        out.push(arc_edifice(pos, cfg));
    }

    // ── Rifts: Divergent ∧ Continental, thinned to `arc_spacing_km`. ──
    let rift_spacing_norm2 = arc_spacing_norm2;
    let mut rift_kept: Vec<(f32, f32)> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            if !matches!(boundary.boundary_type.get(i, j), BoundaryType::Divergent) {
                continue;
            }
            if !matches!(state.plate_type.get(i, j), PlateType::Continental) {
                continue;
            }
            let pos = to_uv(i, j);
            if rift_kept.iter().any(|&k| torus_dist2(k, pos) < rift_spacing_norm2) {
                continue;
            }
            rift_kept.push(pos);
            out.push(rift_edifice(pos, cfg));
        }
    }

    // ── Hotspot chains: pick intraplate seeds, extend the chain DOWNSTREAM in the
    //    −velocity direction (a point that sat over the hotspot earlier has since
    //    drifted −v·age), so age increases opposite to plate motion. ──
    let mut rng = ChaCha8Rng::seed_from_u64(seed.derive_seed("volcanism_hotspot"));
    let internal = |i: usize, j: usize| -> bool {
        matches!(boundary.boundary_type.get(i, j), BoundaryType::Internal)
    };
    let spacing_norm = cfg.hotspot_spacing_km / torus_km;
    let mut placed_hotspots = 0usize;
    let mut attempts = 0usize;
    while placed_hotspots < cfg.n_hotspots && attempts < cfg.n_hotspots * 64 {
        attempts += 1;
        let i = rng.gen_range(0..nx);
        let j = rng.gen_range(0..ny);
        if !internal(i, j) {
            continue; // intraplate only
        }
        let pid = state.plate_id.get(i, j);
        let (vx, vy) = kin.velocities[pid as usize];
        let vl = ((vx * vx + vy * vy).sqrt()) as f32;
        if vl < 1e-9 {
            continue; // a stationary plate builds no age-progressive chain
        }
        let vhat = (vx as f32 / vl, vy as f32 / vl);
        let (u0, v0) = to_uv(i, j);
        for k in 0..cfg.hotspot_chain_len {
            // Edifice k sits at hotspot − k·spacing·v̂ (downstream of motion).
            let pos = (
                (u0 - vhat.0 * spacing_norm * k as f32).rem_euclid(1.0),
                (v0 - vhat.1 * spacing_norm * k as f32).rem_euclid(1.0),
            );
            let age_frac = if cfg.hotspot_chain_len > 1 {
                k as f32 / (cfg.hotspot_chain_len - 1) as f32
            } else {
                0.0
            };
            out.push(hotspot_edifice(pos, age_frac, k == 0, cfg));
        }
        placed_hotspots += 1;
    }

    out
}

/// Arc sites: `(foot_uv, offset_uv)` for every O-C margin cell, BEFORE thinning.
/// The inboard direction is the SMOOTHED boundary normal — the distance-weighted
/// sum of directions to oceanic differing-plate cells within a radius-2 window,
/// negated (continent-ward). A radius-2 window approximates the true local normal
/// far better than the axis-aligned 4-neighbour sum (which can sit up to 45° off a
/// diagonal margin and shorten the perpendicular offset). Exposed so the placement
/// test can measure the applied offset directly. Offset magnitude is exactly
/// `trench_arc_offset_km` (unit normal × offset); its PERPENDICULAR component
/// depends on how well the normal is estimated — hence the smoothing.
#[must_use]
pub fn arc_sites(
    state: &C1State,
    boundary: &crate::tectonics_c1::boundary_classification::BoundaryInfo,
    torus_km: f32,
    cfg: &VolcanismConfig,
) -> Vec<((f32, f32), (f32, f32))> {
    let nx = state.plate_id.nx();
    let ny = state.plate_id.ny();
    let arc_mask = oc_override_seed_mask(boundary, &state.plate_id, &state.plate_type);
    let offset_norm = cfg.trench_arc_offset_km / torus_km;
    const R: i32 = 2;
    let mut sites = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            if !arc_mask.get(i, j) {
                continue;
            }
            let pid = state.plate_id.get(i, j);
            // Smoothed trench direction: distance-weighted sum of directions to
            // oceanic differing-plate cells within radius R.
            let mut trench = (0.0f32, 0.0f32);
            for dj in -R..=R {
                for di in -R..=R {
                    if di == 0 && dj == 0 {
                        continue;
                    }
                    let ni = ((i as i32 + di).rem_euclid(nx as i32)) as usize;
                    let nj = ((j as i32 + dj).rem_euclid(ny as i32)) as usize;
                    if state.plate_id.get(ni, nj) != pid
                        && matches!(state.plate_type.get(ni, nj), PlateType::Oceanic)
                    {
                        let d2 = (di * di + dj * dj) as f32;
                        let w = 1.0 / d2; // nearer cells weigh more
                        let l = d2.sqrt();
                        trench.0 += di as f32 / l * w;
                        trench.1 += dj as f32 / l * w;
                    }
                }
            }
            let tl = (trench.0 * trench.0 + trench.1 * trench.1).sqrt();
            if tl < 1e-6 {
                continue;
            }
            let inboard = (-trench.0 / tl, -trench.1 / tl);
            let foot = ((i as f32 + 0.5) / nx as f32, (j as f32 + 0.5) / ny as f32);
            let pos = (
                (foot.0 + inboard.0 * offset_norm).rem_euclid(1.0),
                (foot.1 + inboard.1 * offset_norm).rem_euclid(1.0),
            );
            sites.push((foot, pos));
        }
    }
    sites
}

/// Maximum crater depth/diameter ratio. Grosse's median crater depth (240 m)
/// applied INDEPENDENTLY of width (Wood: no depth/diameter correlation) is an
/// out-of-domain extrapolation on SMALL craters — a 0.37 km arc crater at 240 m
/// gives D/W ≈ 0.65, deeper than any observed crater lake. Real crater lakes
/// cluster at D/W ≈ 0.1–0.2 (Kawah Ijen 0.20, Poás 0.19; Wood 1980 shows D/W
/// decreasing with age). We cap depth at `CRATER_MAX_DW · diameter`, so Grosse's
/// 240 m is kept for craters ≳ 1 km (shields) and small craters get a realistic
/// depth. Same discipline as the `Wb ≥ 2 km` validity bound.
pub(crate) const CRATER_MAX_DW: f32 = 0.25;

/// Crater depth (m), capped so D/W ≤ [`CRATER_MAX_DW`] (see its doc). Also used by
/// the active-rim reconstruction to recover a crater's depth from its diameter.
pub(crate) fn crater_depth_capped(cfg: &VolcanismConfig, crater_diameter_km: f32) -> f32 {
    cfg.crater_depth_m.min(CRATER_MAX_DW * crater_diameter_km * 1000.0)
}

/// Arc stratocone: andesitic/viscous → steep. Wood constructional geometry.
fn arc_edifice(pos: (f32, f32), cfg: &VolcanismConfig) -> Edifice {
    let wb = cfg.arc_basal_km.max(2.0); // enforce Wood/Grosse validity (Wb ≥ 2 km)
    Edifice {
        center_uv: pos,
        basal_diameter_km: wb,
        height_m: Edifice::stratocone_height_m(wb),
        crater_diameter_km: Edifice::stratocone_crater_km(wb),
        crater_depth_m: crater_depth_capped(cfg, Edifice::stratocone_crater_km(wb)),
        kind: VolcanoKind::Stratocone,
        setting: VolcanoSetting::Arc,
        active: true, // an arc built by ongoing subduction is degassing
        age_frac: 0.0,
    }
}

/// Rift shield: basaltic/fluid → gentle. Grosse & Kervyn H/Wb.
fn rift_edifice(pos: (f32, f32), cfg: &VolcanismConfig) -> Edifice {
    let wb = cfg.shield_basal_km.max(2.0);
    Edifice {
        center_uv: pos,
        basal_diameter_km: wb,
        height_m: Edifice::shield_height_m(wb, cfg.shield_hwb),
        crater_diameter_km: (0.06 * wb + 0.05).min(wb * 0.2), // small central pit
        crater_depth_m: crater_depth_capped(cfg, (0.06 * wb + 0.05).min(wb * 0.2)),
        kind: VolcanoKind::Shield,
        setting: VolcanoSetting::Rift,
        active: true, // active spreading centre
        age_frac: 0.0,
    }
}

/// Hotspot shield, member `k` of a chain: only the youngest (`k = 0`, over the
/// plume) is active; the rest are extinct and increasingly subdued by age.
fn hotspot_edifice(
    pos: (f32, f32),
    age_frac: f32,
    youngest: bool,
    cfg: &VolcanismConfig,
) -> Edifice {
    let wb = cfg.shield_basal_km.max(2.0);
    Edifice {
        center_uv: pos,
        basal_diameter_km: wb,
        height_m: Edifice::shield_height_m(wb, cfg.shield_hwb),
        crater_diameter_km: (0.06 * wb + 0.05).min(wb * 0.2),
        crater_depth_m: crater_depth_capped(cfg, (0.06 * wb + 0.05).min(wb * 0.2)),
        kind: VolcanoKind::Shield,
        setting: VolcanoSetting::Hotspot,
        active: youngest,
        age_frac,
    }
}
