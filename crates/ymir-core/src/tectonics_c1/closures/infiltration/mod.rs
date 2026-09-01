//! H-1 infiltration — the first SUBSURFACE term. The water balance read
//! `runoff = max(0, precip − PE)`, i.e. the ENTIRE precipitation surplus became surface
//! runoff; Ymir had no subsurface water at all. This module derives a per-cell
//! PERMEABILITY from the tectonic closures and turns it into the fraction of the surplus
//! that INFILTRATES and never reaches a lake by surface flow.
//!
//! ## Permeability is NOT erodibility (the distinction that shapes this module)
//!
//! C-3's `build_coarse_k` is an ERODIBILITY. The two correlate only by accident: a basalt
//! is permeable and moderately erodible, a CLAY is highly erodible and NEARLY
//! IMPERMEABLE. So this module builds a SEPARATE class → PERMEABILITY mapping, even
//! though it reads the same lithology classes.
//!
//! ## Anchors — hydraulic conductivity by rock type (Heath 1983 / USGS), m/day
//!
//! | rock | K (m/day) |
//! |---|---|
//! | dense crystalline | 5×10⁻⁸ – 10⁻⁵ |
//! | FRACTURED crystalline | 10⁻³ – 10 |
//! | dense basalt | 10⁻⁶ – 10⁻³ |
//! | fractured basalt | 10⁻⁴ – 1 |
//! | sandstone | 5×10⁻⁵ – 20 |
//! | limestone (karstic) | 5×10⁻⁶ – 100 (10⁻¹ – 10³) |
//! | clay | 5×10⁻⁷ – 10⁻³ (nearly impermeable) |
//!
//! The load-bearing number: **fracturing raises K by 5–6 orders of magnitude** (dense
//! crystalline 10⁻⁸ → fractured 10). Fracture density is therefore the DOMINANT control on
//! rock-mass permeability; the lithological matrix is a minor parallel contribution.
//! Domain of validity: saturated bulk conductivity of the rock mass, shallow crust.
//!
//! ## Composition — DOUBLE POROSITY (Barenblatt 1960; Warren & Root 1963; Kazemi 1969)
//!
//! In a dual-porosity medium the FRACTURES are the source of permeability and
//! connectivity, and the MATRIX is the source of storage. Conductances in parallel add, so
//! `K_eff = K_matrix(class) + density · K_fractured` — NOT a maximum (a maximum would let
//! an erodible-but-impermeable class read as permeable). The exact functional form here is
//! a DERIVATION grounded on that concept, not a verbatim published equation — labelled.
//!
//! ## Slope is deliberately ABSENT, and that is a literature verdict
//!
//! The dedicated review (*Role of slope on infiltration: A review*, J. Hydrol.) reports
//! CONTRADICTORY effects — several studies find infiltration INCREASING with slope (thinner
//! surface seal, lower raindrop impact energy) — and field work finds "slope had minimal
//! predictive power", the real predictors being grass cover, leaf litter, soil organic
//! carbon and bulk density: SOIL properties Ymir does not model. The sources therefore do
//! NOT support a significant, consistently-signed slope term at this scale, so none is
//! included.
//!
//! ## Composition with Budyko (where an error would be invisible)
//!
//! Budyko bounds the SUPPLY: `AET = min(precip, PE)`, `runoff = precip − AET =
//! max(0, precip − PE)`. Infiltration does NOT touch AET; it SPLITS the post-Budyko runoff:
//! `surface_runoff = (precip − AET) · (1 − f_infil)`. The infiltrated part LEAVES the
//! surface balance and is never re-added (no double count — the trap that caught the
//! rain-credit vs catchment-runoff accounting once). First approximation: infiltrated water
//! is lost to deep groundwater, NOT returned as baseflow — which under-supplies slightly,
//! the desired direction. Refinement path: delayed baseflow return re-emerging downstream.

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;
use crate::tectonics_c1::closures::fracture::{FractureConfig, derive_coarse_density};
use crate::tectonics_c1::closures::lithology::{LithologyConfig, upscale_k_to_hd};
use crate::tectonics_c1::closures::volcanism::Edifice;
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_v2::boundaries::plate_type::PlateType;

/// Version of the infiltration derivation. ⚠️ BUMP on any change that moves the field.
pub const INFILTRATION_ALGO: u32 = 1;

// ── Heath (1983) anchors, m/day. Log-mid of the published range per class. PROXY: the
//    class→K assignment is anchored on the table but Ymir's classes are coarser than the
//    table's rock types, so these are representative values, not measurements.
/// Intact crystalline basement — dense crystalline, 5×10⁻⁸–10⁻⁵ → effectively impermeable.
pub const K_MATRIX_HARD: f32 = 1.0e-6;
/// Young rift crust (volcanic + clastic fill) — dense-basalt class, 10⁻⁶–10⁻³.
pub const K_MATRIX_RIFT: f32 = 1.0e-4;
/// Volcaniclastic (fragmental, poorly welded) — the permeable end of the basalt class.
pub const K_MATRIX_VOLCANIC: f32 = 1.0e-3;
/// Fully fractured rock mass — the UPPER end of Heath's fractured-crystalline range
/// (10⁻³–10 m/day). Reached at fracture density 1.
pub const K_FRACTURED: f32 = 10.0;
/// Negligible floor for the fracture term at density 0 (well below the intact matrix, so
/// unfractured rock contributes nothing). The fracture conductivity spans from here to
/// [`K_FRACTURED`] LOG-linearly with density — Heath's dense→fractured jump is 5–6 ORDERS,
/// so a linear interpolation would saturate at any appreciable density and erase the
/// contrast (measured: it did). PROXY: the floor value is a numerical convenience.
pub const K_FRAC_FLOOR: f32 = 1.0e-8;

/// Fracture conductivity (m/day) for a fracture density `d ∈ [0,1]`: log-linear from
/// [`K_FRAC_FLOOR`] (intact) to [`K_FRACTURED`] (fully fractured), i.e. `K = floor ·
/// (max/floor)^d`. Conductivity varies over ORDERS with fracturing, so the interpolation
/// must be geometric, not arithmetic. Labelled derivation (the SPAN is Heath's; the
/// functional form is ours).
#[inline]
#[must_use]
pub fn fracture_conductivity(d: f32) -> f32 {
    let d = d.clamp(0.0, 1.0);
    K_FRAC_FLOOR * (K_FRACTURED / K_FRAC_FLOOR).powf(d)
}

/// H-1 infiltration config. `enabled = false` (default) → no infiltration, the pre-H-1
/// water balance byte-identical.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct InfiltrationConfig {
    /// Master switch. Default `false` (OFF → byte-identical).
    pub enabled: bool,
    /// Upper bound of the infiltrated fraction, at saturating permeability. Anchored on the
    /// Baseflow-Index upper range (BFI 0.55–0.80, USGS Housatonic) — the groundwater share
    /// of streamflow in permeable-bedrock catchments. SWEPT; the success criterion is that
    /// the RESULTING per-class fractions stay inside the published range.
    pub f_cap: f32,
    /// Half-saturation conductivity (m/day) in `f = f_cap · K/(K + k_ref)`: the rainfall
    /// supply rate against which infiltration capacity saturates. PROXY: ~1000 mm/yr ≈
    /// 2.7×10⁻³ m/day. Below `k_ref` the ground cannot swallow the rain; far above it, it can.
    pub k_ref_m_per_day: f32,
}

impl Default for InfiltrationConfig {
    fn default() -> Self {
        Self { enabled: false, f_cap: 0.7, k_ref_m_per_day: 2.7e-3 }
    }
}

impl InfiltrationConfig {
    /// `true` when OFF — so a disabled config can be skipped from the cache key.
    #[must_use]
    pub fn is_disabled(&self) -> bool {
        !self.enabled
    }

    /// The infiltrated fraction for a bulk conductivity `k` (m/day). Saturating: the ground
    /// cannot infiltrate more than the rain supplies.
    #[inline]
    #[must_use]
    pub fn fraction_for_k(&self, k: f32) -> f32 {
        let kr = self.k_ref_m_per_day.max(1e-12);
        (self.f_cap * k / (k + kr)).clamp(0.0, 0.95)
    }
}

/// Coarse MATRIX permeability (m/day) per tectonic cell: intact crystalline basement
/// everywhere, young rift crust where the C-3 rift signal is (continental, `age <
/// rift_age_threshold`). Volcaniclastic is stamped at HD from the edifice footprints.
#[must_use]
pub fn derive_coarse_matrix_k(state: &C1State, litho: &LithologyConfig) -> GridF32 {
    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let mut k = GridF32::new(nx, ny, K_MATRIX_HARD);
    for j in 0..ny {
        for i in 0..nx {
            if matches!(state.plate_type.get(i, j), PlateType::Continental)
                && (state.age.get(i, j) as f32) < litho.rift_age_threshold
            {
                k.set(i, j, K_MATRIX_RIFT);
            }
        }
    }
    k
}

/// Build the per-cell INFILTRATED FRACTION at HD, registered to the terrain via
/// `(sample_origin, sample_size)`. Double porosity: `K_eff = K_matrix + density·K_FRACTURED`
/// (conductances in parallel), then `f = f_cap · K_eff/(K_eff + k_ref)`. Volcaniclastic
/// footprints raise the MATRIX term on the edifice basal discs. No slope term (see module
/// docs). Returns `f ∈ [0, 0.95]`.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn build_hd_infiltration(
    state: &C1State,
    kin: &PlateKinematics,
    litho: &LithologyConfig,
    fracture: &FractureConfig,
    cfg: &InfiltrationConfig,
    edifices: &[Edifice],
    dst_w: usize,
    dst_h: usize,
    sample_origin: [f64; 2],
    sample_size: f64,
    km_per_hd_cell: f32,
) -> Vec<f32> {
    // Matrix conductivity (lithology classes) and fracture density (C-3b), both upscaled
    // with the SAME mapping as the altitude so they register with the terrain.
    let coarse_matrix = derive_coarse_matrix_k(state, litho);
    let mut k_matrix = upscale_k_to_hd(&coarse_matrix, dst_w, dst_h, sample_origin, sample_size);
    let coarse_density = derive_coarse_density(state, kin, fracture, None);
    let density = upscale_k_to_hd(&coarse_density, dst_w, dst_h, sample_origin, sample_size);

    // Volcaniclastic footprints: raise the matrix term on the edifice basal discs (the same
    // discs the C-3 lithology and C-2 volcanism stamp).
    let (so, ss) = ([sample_origin[0] as f32, sample_origin[1] as f32], sample_size as f32);
    for e in edifices {
        let fx = (e.center_uv.0 - so[0]).rem_euclid(1.0) / ss;
        let fy = (e.center_uv.1 - so[1]).rem_euclid(1.0) / ss;
        if fx >= 1.0 || fy >= 1.0 {
            continue;
        }
        let (cx, cy) = (fx * dst_w as f32, fy * dst_h as f32);
        let rb = (e.basal_diameter_km * 0.5) / km_per_hd_cell;
        if rb < 1.0 {
            continue;
        }
        let (i0, i1) =
            ((cx - rb).floor().max(0.0) as usize, ((cx + rb).ceil() as usize).min(dst_w - 1));
        let (j0, j1) =
            ((cy - rb).floor().max(0.0) as usize, ((cy + rb).ceil() as usize).min(dst_h - 1));
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d <= rb {
                    let k = j * dst_w + i;
                    k_matrix[k] = k_matrix[k].max(K_MATRIX_VOLCANIC);
                }
            }
        }
    }

    // Double porosity + saturating fraction.
    k_matrix
        .iter()
        .zip(density.iter())
        .map(|(&km, &d)| cfg.fraction_for_k(km + fracture_conductivity(d)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fraction_saturates_and_orders_by_permeability() {
        let cfg = InfiltrationConfig { enabled: true, ..Default::default() };
        // Intact crystalline basement: K ≪ rainfall supply → essentially no infiltration.
        let f_hard = cfg.fraction_for_k(K_MATRIX_HARD);
        // Fully fractured rock: K ≫ supply → saturates near f_cap.
        let f_frac = cfg.fraction_for_k(K_FRACTURED);
        assert!(f_hard < 0.01, "intact crystalline barely infiltrates: {f_hard}");
        assert!(f_frac > 0.6, "fractured rock saturates near f_cap: {f_frac}");
        assert!(f_frac > f_hard, "permeability must order the fraction");
        // Double porosity is ADDITIVE, and the fracture term spans ORDERS with density, so
        // the contrast lands INSIDE the density range instead of saturating at once.
        let f_low = cfg.fraction_for_k(K_MATRIX_HARD + fracture_conductivity(0.15));
        let f_mid = cfg.fraction_for_k(K_MATRIX_HARD + fracture_conductivity(0.55));
        let f_high = cfg.fraction_for_k(K_MATRIX_HARD + fracture_conductivity(0.9));
        assert!(f_low < 0.1, "barely fractured craton stays tight: {f_low}");
        assert!(f_mid > f_low && f_high > f_mid, "monotone in density: {f_low}/{f_mid}/{f_high}");
        assert!(f_high > 0.6, "heavily fractured belt saturates: {f_high}");
    }

    #[test]
    fn disabled_config_is_flagged() {
        assert!(InfiltrationConfig::default().is_disabled());
    }
}
