//! C-2 volcanism (closures roadmap §2) — built relief added at HD resolution,
//! AFTER the FBM and BEFORE erosion, from the tectonic state that already exists.
//!
//! Unlike the per-step S̃ source terms in the sibling `closures/` modules, this is
//! an HD-stage closure: a 10–20 km edifice is a single cell at the 64² tectonic
//! grid but 75–300 cells at HD, so it can only be resolved after the upscale.
//! Placing it before erosion lets the SAME erosion chain dissect the flanks and
//! breach the old craters (realistic ravined flanks), which is why a crater does
//! not need a depression-filter exemption: it is a genuine deep, wide hollow that
//! the normal lake detection captures and `breach_monotone` then holds flat.
//!
//! ## Anchoring (see ADR 0001, section C-2, for the full bibliography)
//!
//! Every geometric number carries a publication EXCEPT the extinct relief decay:
//! - **Stratocones** (arc, andesitic/viscous → steep): Wood 1978,
//!   `H = 0.122·Wb + 0.450` km and `Wcr = 0.027·Wb + 0.048` km, constructional
//!   geometry of 26 circum-Pacific cones. **Valid only for Wb ≥ 2 km** — below
//!   that they are cinder cones with different relations, so we never place a
//!   sub-2-km edifice with this law.
//! - **Shields** (hotspot/rift, basaltic/fluid → gentle): Grosse & Kervyn 2018,
//!   H/Wb in 0.01–0.1 (central 0.10), flank slopes 1–15°.
//! - **Crater depth**: Grosse et al. 2014, median 240 m, range 100–860 m, and
//!   INDEPENDENT of crater width (Wood: no significant depth/diameter relation) —
//!   so depth is drawn from the range, not computed from width.
//! - **Crater-lake acidity**: Varekamp et al. 2000, active degassing crater lakes
//!   have pH < 2; extinct ones are neutral freshwater — so `active` (not "is a
//!   crater") sets the lake chemistry.
//!
//! The ONE unanchored parameter is [`VolcanismConfig::extinct_relief_decay`]: an
//! explicit PROXY (marked as such) for the post-extinction erosion the single,
//! time-less erosion pass cannot date.

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;

pub mod placement;
pub use placement::{TRENCH_ARC_OFFSET_KM_DEFAULT, place_edifices};

/// Version of the volcanism placement/geometry code. ⚠️ BUMP on any change to
/// `place_edifices` or `apply_edifices` that moves an edifice or a crater, so a
/// volcanism-enabled eroded cache invalidates — WITHOUT touching the shared
/// `ALGO_UPSCALE_EROSION` (which would needlessly recompute every volcanism-OFF
/// terrain, byte-identical though it is). Only added to the key when volcanism is
/// enabled. v2: fixed the `apply_edifices` vertical-mirror + window-wrap bugs
/// (edifices were flipped and >half were dropped at a non-zero sample_origin).
/// v3: D/W-capped crater depth (arc craters were D/W 0.65) + active-rim
/// reconstruction (active craters re-close after erosion → hold lakes).
/// v4: active crater bowls are PROTECTED from the relief-v3 breach
/// (`breach_monotone_protected`) so they survive to the crater-lake stage — without
/// it the breach re-breached them and every crater came out dry (the export showed
/// 0 acidic lakes while the balance said one should fill).
pub const VOLCANISM_ALGO: u32 = 4;

/// Edifice kind, driven by the PHYSICAL setting (magma composition / viscosity),
/// not a style toggle: arc magmas are andesitic/viscous → steep stratocones;
/// hotspot and rift magmas are basaltic/fluid → gentle shields.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VolcanoKind {
    /// Steep composite cone (arc). Wood 1978 constructional law.
    Stratocone,
    /// Gentle basaltic shield (hotspot / rift). Grosse & Kervyn 2018.
    Shield,
}

/// Tectonic setting the edifice was derived from (for reporting / placement audit).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VolcanoSetting {
    /// Island arc along an oceanic-continental convergent margin.
    Arc,
    /// Intraplate hotspot (age-progressive chain along plate motion).
    Hotspot,
    /// Continental rift / divergent boundary.
    Rift,
}

/// One volcanic edifice. Position is in NORMALIZED coarse-torus coordinates
/// `[0,1)`; every size is PHYSICAL (km / metres), never cells — the conversion to
/// HD pixels happens at the point of use in [`apply_edifices`].
#[derive(Debug, Clone)]
pub struct Edifice {
    /// Centre in normalized coarse-torus coordinates `(u, v)`, each in `[0,1)`.
    pub center_uv: (f32, f32),
    /// Basal diameter Wb (km). Always `≥ 2` (Wood/Grosse validity domain).
    pub basal_diameter_km: f32,
    /// Constructional summit height above the base (m).
    pub height_m: f32,
    /// Summit crater diameter (km).
    pub crater_diameter_km: f32,
    /// Crater depth below the rim (m) — Grosse range, independent of width.
    pub crater_depth_m: f32,
    /// Stratocone (arc) or shield (hotspot/rift).
    pub kind: VolcanoKind,
    /// Tectonic setting it was derived from.
    pub setting: VolcanoSetting,
    /// Actively degassing → acidic crater lake (Varekamp: pH < 2). Extinct → the
    /// crater lake, if any, is neutral freshwater.
    pub active: bool,
    /// Age fraction along a hotspot chain: `0` = youngest (over the hotspot),
    /// `1` = oldest. Drives the extinct relief decay. `0` for arcs/rifts.
    pub age_frac: f32,
}

impl Edifice {
    /// Wood 1978 constructional stratocone height (m) from basal diameter (km).
    /// `H = 0.122·Wb + 0.450` km (n=17, r=0.95). Valid for Wb ≥ 2 km.
    #[must_use]
    pub fn stratocone_height_m(wb_km: f32) -> f32 {
        (0.122 * wb_km + 0.450) * 1000.0
    }

    /// Wood 1978 stratocone summit-crater diameter (km): `Wcr = 0.027·Wb + 0.048`.
    #[must_use]
    pub fn stratocone_crater_km(wb_km: f32) -> f32 {
        0.027 * wb_km + 0.048
    }

    /// Shield height (m) from basal diameter (km) and an H/Wb ratio in 0.01–0.1
    /// (Grosse & Kervyn 2018). Shields have no systematic summit crater at this
    /// scale — a small central pit is used ([`VolcanismConfig::crater_depth_m`]).
    #[must_use]
    pub fn shield_height_m(wb_km: f32, hwb: f32) -> f32 {
        hwb * wb_km * 1000.0
    }
}

/// Configuration for the C-2 volcanism closure. `enabled = false` (default) →
/// no edifice is placed, byte-identical to the pre-C-2 pipeline.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct VolcanismConfig {
    /// Master switch. Default `false` (OFF → byte-identical).
    pub enabled: bool,
    /// Physical span of the HD grid, in km (= the export `window_km`). The HD
    /// horizontal scale is `domain_km / target_size` km/cell — the SAME km/cell
    /// the export uses, so edifice sizes stay physically correct. Kept in km and
    /// converted to cells only at the point of use.
    pub domain_km: f32,
    /// Trench→arc horizontal offset (km). The arc-trench-gap literature gives
    /// 100–300 km; our model has no slab-depth field so we place by this
    /// horizontal offset (see ADR — Syracuse & Abers give a slab DEPTH, not this
    /// distance). Stored in km, converted at use.
    pub trench_arc_offset_km: f32,
    /// Along-arc spacing between arc edifices (km).
    pub arc_spacing_km: f32,
    /// Basal diameter assigned to arc stratocones (km, ≥ 2).
    pub arc_basal_km: f32,
    /// Number of intraplate hotspots.
    pub n_hotspots: usize,
    /// Edifices per hotspot chain (age-progressive along plate motion).
    pub hotspot_chain_len: usize,
    /// Along-chain spacing between hotspot edifices (km).
    pub hotspot_spacing_km: f32,
    /// Basal diameter assigned to hotspot / rift shields (km, ≥ 2).
    pub shield_basal_km: f32,
    /// Shield H/Wb ratio (Grosse & Kervyn 2018 range 0.01–0.1).
    pub shield_hwb: f32,
    /// Crater depth (m) — Grosse median 240, range 100–860, independent of width.
    pub crater_depth_m: f32,
    /// **PROXY, not from a publication.** Relief multiplier applied to a fully
    /// aged-out extinct edifice (`age_frac = 1`): stands in for the
    /// post-extinction erosion that the single, time-less erosion pass cannot
    /// date (a young active cone is pristine; an old extinct cone is subdued).
    /// `1.0` = no decay. This is an unanchored number in C-2 — see ADR.
    pub extinct_relief_decay: f32,
    /// **Active-rim reconstruction (PROXY, composed — no single named law).** The
    /// erosion pass breaches EVERY crater rim; the measurement showed active craters
    /// with sufficient inflow (margin > 1) yet no lake, because nothing maintains
    /// their rim. An ACTIVE volcano rebuilds its rim (construction, Wood 1978)
    /// faster than erosion tears it down (Wood 1980) — so after erosion the active
    /// crater rim is reconstructed to this fraction of the crater depth (a
    /// construction-vs-erosion RATE outcome, not a protection flag: `1.0` =
    /// construction fully wins → rim re-closed → the crater is ELIGIBLE and the
    /// water balance then decides; `0.0` = construction loses → breaches like an
    /// extinct one). Extinct edifices are never reconstructed. Default `1.0`.
    pub active_rim_rebuild: f32,
}

impl Default for VolcanismConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            domain_km: 400.0,
            // Provisional — see TRENCH_ARC_OFFSET_KM_DEFAULT for its provenance
            // and why it is the least-certain number in C-2.
            trench_arc_offset_km: TRENCH_ARC_OFFSET_KM_DEFAULT,
            // Real arc volcano along-strike spacing is ~50–80 km.
            arc_spacing_km: 70.0,
            arc_basal_km: 12.0, // Grosse composite WB median ~10–13 km
            n_hotspots: 2,
            hotspot_chain_len: 5,
            hotspot_spacing_km: 60.0,
            shield_basal_km: 20.0,
            shield_hwb: 0.08,           // central-ish in 0.01–0.1
            crater_depth_m: 240.0,      // Grosse median
            extinct_relief_decay: 0.35, // PROXY
            active_rim_rebuild: 1.0,    // PROXY (construction wins for active rims)
        }
    }
}

/// A placed crater, in HD pixel coordinates, for downstream lake classification.
/// Carried out of [`apply_edifices`] so the lake stage can type a lake by
/// intersection without re-deriving anything. `Serialize`/`Deserialize` so the
/// records can travel WITH the cached terrain (the `{ heightmap, craters }`
/// bundle) — a hit that returned the terrain without the craters would mistype
/// every crater lake silently.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CraterRecord {
    /// Crater centre in HD pixels (heightmap DATA space, row ↔ coarse row).
    pub center_px: (f32, f32),
    /// Crater radius in HD pixels (the rim) — a size proxy for the viz marker.
    pub radius_px: f32,
    /// Active degassing → the crater lake (if any) is acidic (pH < 2); the viz
    /// marks active volcanoes red, extinct grey. Kept minimal (no kind/setting) so
    /// the cached bundle format is unchanged — existing volcanism caches still load
    /// and show markers without a recompute.
    pub active: bool,
}

/// Result of applying the edifices to the HD heightmap.
pub struct VolcanismApplied {
    /// Craters placed, in HD pixel coordinates (for lake typing).
    pub craters: Vec<CraterRecord>,
    /// Number of edifices actually rendered into the window.
    pub rendered: usize,
}

/// Add the edifices to `heightmap` (norm units, sea = 0.5). Positions map from
/// normalized coarse-torus coords into the HD window `[sample_origin,
/// sample_origin+sample_size]`; sizes convert to cells via `km_per_hd_cell`.
///
/// `full_range_m` is the norm→metre span (e.g. `2·1.13·depth_scale_m`), so a
/// height in metres becomes `height_m / full_range_m` in norm. Returns the crater
/// records (HD pixels) for the lake stage. Land-only: an edifice on a sub-sea
/// cell still builds up (a volcanic island), but the profile is added to the
/// existing bed so a deep-ocean seamount that does not breach the surface simply
/// stays submerged.
#[must_use]
pub fn apply_edifices(
    heightmap: &mut GridF32,
    edifices: &[Edifice],
    sample_origin: [f64; 2],
    sample_size: f64,
    km_per_hd_cell: f32,
    full_range_m: f32,
    cfg: &VolcanismConfig,
) -> VolcanismApplied {
    let (w, h) = (heightmap.width, heightmap.height);
    let mut craters = Vec::new();
    let mut rendered = 0usize;
    let so = [sample_origin[0] as f32, sample_origin[1] as f32];
    let ss = sample_size as f32;

    for e in edifices {
        // Map normalized coarse-torus centre into HD window pixels. The coarse
        // field is sampled PERIODICALLY (like the FBM's sample_bilinear_periodic),
        // so the window `[origin, origin+size)` WRAPS the torus: the offset must be
        // taken mod 1 before dividing by the size, or every edifice at `u < origin`
        // (here > half the domain, since origin_v = 0.578) is wrongly dropped.
        let fx = (e.center_uv.0 - so[0]).rem_euclid(1.0) / ss;
        let fy = (e.center_uv.1 - so[1]).rem_euclid(1.0) / ss;
        if fx >= 1.0 || fy >= 1.0 {
            continue; // outside a sub-window (size < 1); full-domain (size 1) never skips
        }
        // heightmap.data is row-major with data row ↔ coarse row (the FBM writes
        // it with NO vertical flip — the north-up mirror is applied only when
        // rendering data → texture, not in the data array itself). lake_map shares
        // this layout, so the crater centre in data space is (fx·w, fy·h). Using
        // the render mirror here would place the cone at the vertically-flipped
        // position, off the real tectonic margin, and misalign it with lake_map.
        let cx = fx * w as f32;
        let cy = fy * h as f32;

        let rb_px = (e.basal_diameter_km * 0.5) / km_per_hd_cell;
        let rc_px = (e.crater_diameter_km * 0.5) / km_per_hd_cell;
        if rb_px < 1.0 {
            continue; // below one HD cell — nothing to render
        }
        // Extinct relief decay (PROXY): a young active cone is pristine, an old
        // extinct cone is subdued. Active edifices never decay.
        let relief_factor =
            if e.active { 1.0 } else { 1.0 - (1.0 - cfg.extinct_relief_decay) * e.age_frac };
        let height_norm = e.height_m / full_range_m * relief_factor;
        // Rim height above the base (the cone height at the crater rim radius).
        let rim_frac = (1.0 - (rc_px / rb_px)).max(0.0);
        let crater_depth_norm = (e.crater_depth_m / full_range_m) * relief_factor;

        let (imin, imax) =
            ((cx - rb_px).floor().max(0.0) as usize, ((cx + rb_px).ceil() as usize).min(w - 1));
        let (jmin, jmax) =
            ((cy - rb_px).floor().max(0.0) as usize, ((cy + rb_px).ceil() as usize).min(h - 1));
        for j in jmin..=jmax {
            for i in imin..=imax {
                let dx = i as f32 + 0.5 - cx;
                let dy = j as f32 + 0.5 - cy;
                let r = (dx * dx + dy * dy).sqrt();
                if r > rb_px {
                    continue;
                }
                // Linear cone: height falls from summit to 0 at the basal radius.
                // The linear profile reproduces the published AVERAGE flank slope
                // (H over Wb/2) that the H/Wb law encodes.
                let add = if r >= rc_px || rc_px < 1.0 {
                    height_norm * (1.0 - r / rb_px)
                } else {
                    // Inside the crater: a bowl from the rim height down to the
                    // floor at the centre.
                    let rim = height_norm * rim_frac;
                    let t = r / rc_px; // 0 at centre, 1 at rim
                    (rim - crater_depth_norm) + crater_depth_norm * t
                };
                let k = j * w + i;
                heightmap.data[k] = (heightmap.data[k] + add).clamp(0.0, 1.0);
            }
        }
        if rc_px >= 1.0 {
            craters.push(CraterRecord { center_px: (cx, cy), radius_px: rc_px, active: e.active });
        }
        rendered += 1;
    }
    VolcanismApplied { craters, rendered }
}

/// C-2 active-rim reconstruction — run AFTER the erosion pass, which breaches every
/// crater rim indiscriminately. An ACTIVE volcano rebuilds its rim faster than
/// erosion tears it down (construction, Wood 1978, vs destruction, Wood 1980), so
/// its crater re-closes and becomes ELIGIBLE for a lake (the water balance then
/// decides). EXTINCT craters are left breached. The rim ring is raised to
/// `floor + depth·active_rim_rebuild`, re-forming a closed bowl; the erosion-carved
/// FLANKS (outside the rim) are untouched, so an active edifice reads as dissected
/// flanks + an intact sharp rim (as real active craters do). `active_rim_rebuild`
/// is a construction/erosion RATE outcome (see [`VolcanismConfig`]), a labelled
/// PROXY. Physical km/m throughout. No-op when disabled or `active_rim_rebuild = 0`.
pub fn reconstruct_active_rims(
    heightmap: &mut GridF32,
    craters: &[CraterRecord],
    cfg: &VolcanismConfig,
    km_per_hd_cell: f32,
    full_range_m: f32,
) {
    if cfg.active_rim_rebuild <= 0.0 {
        return;
    }
    let (w, h) = (heightmap.width, heightmap.height);
    for c in craters.iter().filter(|c| c.active) {
        let (cx, cy, rc) = (c.center_px.0, c.center_px.1, c.radius_px.max(1.0));
        let diameter_km = 2.0 * rc * km_per_hd_cell;
        let depth_m = placement::crater_depth_capped(cfg, diameter_km);
        let depth_norm = depth_m / full_range_m * cfg.active_rim_rebuild;
        let (i0, i1) = (
            (cx - 1.2 * rc).floor().max(0.0) as usize,
            ((cx + 1.2 * rc).ceil() as usize).min(w - 1),
        );
        let (j0, j1) = (
            (cy - 1.2 * rc).floor().max(0.0) as usize,
            ((cy + 1.2 * rc).ceil() as usize).min(h - 1),
        );
        // Rim base = median elevation of the rim ring (the reconstructed summit the
        // cone rebuilds to), robust to the erosion notch. The bowl floor sits
        // `depth_norm` below it.
        let mut rim: Vec<f32> = Vec::new();
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if (0.9 * rc..=1.1 * rc).contains(&d) {
                    rim.push(heightmap.data[j * w + i]);
                }
            }
        }
        if rim.is_empty() {
            continue;
        }
        rim.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let rim_h = (rim[rim.len() / 2]).min(1.0); // median rim height
        let floor_h = (rim_h - depth_norm).max(0.0);
        // Re-stamp a CLEAN closed bowl (construction rebuilds the whole crater, not
        // just a thin ring — the incised interior otherwise keeps a drainage path
        // and never reads as a closed depression). Interior [0,rc]: a bowl from the
        // floor at the centre to the rim at rc (SET, overwriting the incision). Rim
        // ring [rc,1.1·rc]: raised to the rim height (max, keeping the outer flank).
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                let k = j * w + i;
                if d <= rc {
                    heightmap.data[k] = (floor_h + depth_norm * (d / rc)).clamp(0.0, 1.0);
                } else if d <= 1.1 * rc {
                    heightmap.data[k] = heightmap.data[k].max(rim_h);
                }
            }
        }
    }
}

/// C-2 lake typing: a crater is only a crater LAKE if it holds water. For each
/// crater, find the detected lake occupying it (the majority `lake_map` id within
/// the crater rim); if one exists, type it `CraterAcidic` (active degassing,
/// Varekamp pH < 2) or `CraterNeutral` (extinct → ordinary freshwater). A crater
/// with no lake (insufficient inflow / arid) stays dry relief — the accumulated
/// "a lake must have water" invariant applying here. Returns
/// `(craters_with_lake, dry_craters)`. Mutates only the matched lakes' type; the
/// footprint / level / outlet invariants they already satisfy are untouched.
pub fn classify_crater_lakes(
    lakes: &mut [crate::tectonics_c1::drainage::C1Lake],
    lake_map: &[u32],
    w: usize,
    h: usize,
    craters: &[CraterRecord],
) -> (usize, usize) {
    use crate::tectonics_c1::drainage::LakeType;
    use std::collections::HashMap;
    let (mut with_lake, mut dry) = (0usize, 0usize);
    for c in craters {
        let (cx, cy, r) = (c.center_px.0, c.center_px.1, c.radius_px.max(1.0));
        let imin = (cx - r).floor().max(0.0) as usize;
        let imax = ((cx + r).ceil() as usize).min(w.saturating_sub(1));
        let jmin = (cy - r).floor().max(0.0) as usize;
        let jmax = ((cy + r).ceil() as usize).min(h.saturating_sub(1));
        let mut counts: HashMap<u32, usize> = HashMap::new();
        for j in jmin..=jmax {
            for i in imin..=imax {
                let dx = i as f32 + 0.5 - cx;
                let dy = j as f32 + 0.5 - cy;
                if dx * dx + dy * dy > r * r {
                    continue;
                }
                let id = lake_map[j * w + i];
                if id != 0 {
                    *counts.entry(id).or_default() += 1;
                }
            }
        }
        match counts.into_iter().max_by_key(|&(_, n)| n) {
            Some((id, _)) => {
                if let Some(l) = lakes.iter_mut().find(|l| l.base.id == id) {
                    l.lake_type =
                        if c.active { LakeType::CraterAcidic } else { LakeType::CraterNeutral };
                    with_lake += 1;
                } else {
                    dry += 1;
                }
            }
            None => dry += 1,
        }
    }
    (with_lake, dry)
}

/// Minimum crater-LAKE area (km²) — the inventory floor for a crater lake. Crater
/// lakes are small by nature (Kawah Ijen 0.8 km², Pavin 0.4 km²), well below the
/// generic `lake_min_area_km2 = 5` that `detect_lakes` applies to filter FBM-noise
/// ponds. But too LOW a floor (an early 4-cell ≈ 0.01 km² version) admitted 0.06–
/// 0.11 km² PUDDLES in every wet climate — a majority of active craters holding,
/// the "acidic lakes everywhere" red flag. Set at 0.2 km²: below Pavin (0.4) and
/// small maars, above the marginal puddles, so a crater LAKE is a real lake and the
/// two thresholds stay coherent (5 km² excludes noise micro-lakes on generic
/// terrain; 0.2 km² admits real small crater lakes on the physical crater basin).
const CRATER_LAKE_MIN_AREA_KM2: f32 = 0.2;

/// C-2 crater-lake DETECTION + typing. The generic `detect_lakes` discards a
/// crater as too small (~1.2 km² < the 5 km² noise floor), so crater lakes never
/// reach the water balance. This is a dedicated pass (like `below_sea_basin_lakes`)
/// that, for each ACTIVE crater (reconstructed → a closed bowl), runs the SAME
/// balance the general lakes use — inflow (max catchment runoff) vs evaporation —
/// and fills it: exorheic full lake if `a_eq ≥ a_sill`, an endorheic pond of the
/// equilibrium area otherwise, DRY if the pond would be below the inventory floor.
/// Extinct craters are breached (not reconstructed) → never closed → no lake, so
/// only active craters can hold one, and it is acidic (Varekamp pH < 2). Adds the
/// lakes to `lakes` and marks `lake_map`. Returns `(held, dry)`. Physical units.
#[allow(clippy::too_many_arguments)]
pub fn detect_crater_lakes(
    eroded: &GridF32,
    filled: &GridF32,
    runoff: &[f32],
    temperature: &GridF32,
    craters: &[CraterRecord],
    cell_km2: f32,
    ss: &crate::tectonics_c1::closures::oceanic_bathymetry::SteinSteinParams,
    lakes: &mut Vec<crate::tectonics_c1::drainage::C1Lake>,
    lake_map: &mut [u32],
) -> (usize, usize) {
    use crate::lakes::detection::Lake;
    use crate::tectonics_c1::drainage::{C1Lake, LakeType, potential_evaporation_mm};
    use crate::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    let (w, h) = (eroded.width, eroded.height);
    let mut next_id: u32 = 2_000_001;
    let (mut held, mut dry) = (0usize, 0usize);
    // Crater-lake area floor in cells (see CRATER_LAKE_MIN_AREA_KM2).
    let min_cells = ((CRATER_LAKE_MIN_AREA_KM2 / cell_km2).ceil() as usize).max(1);
    for c in craters.iter().filter(|c| c.active) {
        let (cx, cy, rc) = (c.center_px.0, c.center_px.1, c.radius_px.max(1.0));
        let r = 1.15 * rc;
        let (i0, i1) = ((cx - r).floor().max(0.0) as usize, ((cx + r).ceil() as usize).min(w - 1));
        let (j0, j1) = ((cy - r).floor().max(0.0) as usize, ((cy + r).ceil() as usize).min(h - 1));
        // Flooded bowl cells (filled above the bed) + the fill sill.
        let mut basin: Vec<(usize, f32)> = Vec::new();
        let mut inflow = 0.0f32;
        let mut sill = 0.0f32;
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d > r {
                    continue;
                }
                let k = j * w + i;
                if filled.data[k] > eroded.data[k] + 1e-6 {
                    basin.push((k, eroded.data[k]));
                    inflow = inflow.max(runoff[k]);
                    sill = sill.max(filled.data[k]);
                }
            }
        }
        if basin.len() < min_cells {
            dry += 1; // not a closed bowl (breached) or too small
            continue;
        }
        basin.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let floor_norm = basin[0].1;
        let pe = potential_evaporation_mm(temperature.data[basin[0].0]).max(1.0);
        let a_eq_km2 = inflow / pe;
        let a_sill_km2 = basin.len() as f32 * cell_km2;
        let (n, level_norm) = if a_eq_km2 >= a_sill_km2 {
            (basin.len(), sill) // overflow → full crater lake at the sill
        } else {
            let n_eq = (a_eq_km2 / cell_km2).floor() as usize;
            (n_eq, basin.get(n_eq.saturating_sub(1)).map(|c| c.1).unwrap_or(floor_norm))
        };
        if n < min_cells {
            dry += 1; // inflow can't sustain a lake above the floor
            continue;
        }
        let id = next_id;
        next_id += 1;
        for &(k, _) in basin.iter().take(n) {
            lake_map[k] = id;
        }
        let level_m = c1_altitude_norm_to_metres(level_norm, ss);
        let floor_m = c1_altitude_norm_to_metres(floor_norm, ss);
        lakes.push(C1Lake {
            base: Lake {
                id,
                surface_elevation: level_norm,
                max_depth: level_norm - floor_norm,
                area: n,
                basin_id: 0,
                outlet: (cx as u32, cy as u32),
                shallow: false,
            },
            level_m,
            depth_m: level_m - floor_m,
            area_km2: n as f32 * cell_km2,
            lake_type: LakeType::CraterAcidic,
        });
        held += 1;
    }
    (held, dry)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sanity table: the constructional laws must reproduce named real volcanoes
    /// within the published scatter. Wood 1978 (stratocones) and Grosse & Kervyn
    /// 2018 (shields). This is the check that catches a units error (the kind that
    /// produced a 37 m channel for a 50 km² basin).
    #[test]
    fn morphometry_matches_named_volcanoes() {
        // Stratocones — Wood H = 0.122·Wb + 0.45 km.
        // Mayon: Wb ≈ 20 km, H ≈ 2.46 km. Fuji: Wb ≈ 40 km (broad), H ≈ 3 km.
        let mayon = Edifice::stratocone_height_m(20.0);
        assert!((mayon - 2890.0).abs() < 600.0, "Mayon ~2.4–2.9 km, got {mayon} m");
        // H/Wb stays in the Grosse composite band (0.01–0.30) for Wb 2–22 km.
        for wb in [2.0f32, 8.0, 15.0, 22.0] {
            let hwb = Edifice::stratocone_height_m(wb) / (wb * 1000.0);
            assert!((0.01..=0.35).contains(&hwb), "stratocone H/Wb {hwb} out of band at Wb {wb}");
        }
        // Crater width Wcr = 0.027·Wb + 0.048 km: for a composite the crater is a
        // small fraction of the base (~0.1 at Wb 2 km), NOT the wide cinder crater.
        let wcr_at2 = Edifice::stratocone_crater_km(2.0);
        assert!(
            (wcr_at2 - 0.10).abs() < 0.03,
            "composite crater at Wb 2 km ~0.1 km, got {wcr_at2}"
        );

        // Shields — H/Wb 0.01–0.1. Mauna Loa: Wb ≈ 120 km, subaerial H ≈ 4 km →
        // H/Wb ≈ 0.033, inside the band.
        let ml = Edifice::shield_height_m(120.0, 0.033);
        assert!((ml - 3960.0).abs() < 500.0, "Mauna Loa ~4 km, got {ml} m");
        for hwb in [0.01f32, 0.05, 0.1] {
            let h = Edifice::shield_height_m(20.0, hwb);
            assert!(h > 0.0 && (h / (20.0 * 1000.0) - hwb).abs() < 1e-6);
        }
    }

    /// Applying a single stratocone raises a summit with a central crater dip and
    /// leaves the far field untouched. Physical units: on a 400 km / 512-cell grid
    /// (0.78 km/cell) a 12 km cone spans ~15 cells.
    #[test]
    fn apply_single_stratocone_builds_cone_with_crater() {
        let n = 512usize;
        // HD-like horizontal scale so a 12 km cone (rb ~77 cells) and its
        // 0.37 km crater (rc ~2.4 cells) both resolve, as they do at 8192².
        let domain_km = 40.0f32;
        let km_per_cell = domain_km / n as f32;
        let full_range_m = 11302.0f32;
        let mut hm = GridF32::new(n, n, 0.5); // flat land at sea+0
        let wb = 12.0f32;
        let e = Edifice {
            center_uv: (0.5, 0.5),
            basal_diameter_km: wb,
            height_m: Edifice::stratocone_height_m(wb),
            crater_diameter_km: Edifice::stratocone_crater_km(wb),
            crater_depth_m: 240.0,
            kind: VolcanoKind::Stratocone,
            setting: VolcanoSetting::Arc,
            active: true,
            age_frac: 0.0,
        };
        let cfg = VolcanismConfig { enabled: true, domain_km, ..Default::default() };
        let out = apply_edifices(&mut hm, &[e], [0.0, 0.0], 1.0, km_per_cell, full_range_m, &cfg);
        assert_eq!(out.rendered, 1);
        assert_eq!(out.craters.len(), 1);
        // Centre pixel (north-up mirrored): summit region raised well above base.
        let cix = n / 2;
        let ciy = n / 2;
        let summit_ring = hm.get(cix as i32 + 8, ciy as i32); // on the flank, off the crater
        assert!(summit_ring > 0.55, "flank should be raised, got {summit_ring}");
        // Far field untouched (well beyond the basal radius).
        let far = hm.get(20, 20);
        assert!((far - 0.5).abs() < 1e-6, "far field must stay at base, got {far}");
        // Crater floor sits below the rim.
        let center = hm.get(cix as i32, ciy as i32);
        assert!(center < summit_ring + 0.05, "crater centre must dip vs flank, got {center}");
    }

    /// Crater-lake typing: a crater is a crater LAKE only if it holds water. A
    /// crater over a detected lake is typed by activity; a crater over dry terrain
    /// stays relief (counted dry), not a lake.
    #[test]
    fn crater_lake_typing_only_wet_craters() {
        use crate::lakes::detection::Lake;
        use crate::tectonics_c1::drainage::{C1Lake, LakeType};
        let (w, h) = (20usize, 20);
        let mut lake_map = vec![0u32; w * h];
        for j in 4..7 {
            for i in 4..7 {
                lake_map[j * w + i] = 7; // lake id 7 blob around (5,5)
            }
        }
        let mut lakes = vec![C1Lake {
            base: Lake {
                id: 7,
                surface_elevation: 0.5,
                max_depth: 0.01,
                area: 9,
                basin_id: 0,
                outlet: (0, 0),
                shallow: false,
            },
            level_m: 10.0,
            depth_m: 20.0,
            area_km2: 5.0,
            lake_type: LakeType::Exorheic,
        }];
        let cr = |cx, cy, active| CraterRecord { center_px: (cx, cy), radius_px: 2.0, active };
        let craters = vec![cr(5.0, 5.0, true), cr(15.0, 15.0, false)];
        let (with_lake, dry) = classify_crater_lakes(&mut lakes, &lake_map, w, h, &craters);
        assert_eq!(with_lake, 1, "the crater over the lake must be typed");
        assert_eq!(dry, 1, "the crater over dry ground must stay dry relief");
        assert_eq!(lakes[0].lake_type, LakeType::CraterAcidic, "active crater lake → acidic");
    }

    /// The extinct relief decay (PROXY) lowers an aged-out extinct edifice vs an
    /// identical active one.
    #[test]
    fn extinct_edifice_is_subdued() {
        let n = 256usize;
        let domain_km = 40.0f32; // HD-like scale (see the stratocone test)
        let km_per_cell = domain_km / n as f32;
        let mk = |active: bool, age: f32| {
            let mut hm = GridF32::new(n, n, 0.5);
            let e = Edifice {
                center_uv: (0.5, 0.5),
                basal_diameter_km: 20.0,
                height_m: Edifice::shield_height_m(20.0, 0.08),
                crater_diameter_km: 0.3,
                crater_depth_m: 240.0,
                kind: VolcanoKind::Shield,
                setting: VolcanoSetting::Hotspot,
                active,
                age_frac: age,
            };
            let cfg = VolcanismConfig { enabled: true, domain_km, ..Default::default() };
            let _ = apply_edifices(&mut hm, &[e], [0.0, 0.0], 1.0, km_per_cell, 11302.0, &cfg);
            hm.get((n / 2) as i32 + 6, (n / 2) as i32)
        };
        let active = mk(true, 0.0);
        let extinct_old = mk(false, 1.0);
        assert!(
            extinct_old < active,
            "extinct old edifice must be lower: {extinct_old} vs {active}"
        );
    }
}
