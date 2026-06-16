//! #165 C1 climate — precipitation by a 1D CONSERVATIVE moisture transport
//! along the fixed prevailing wind. Complete (orographic soak + rain shadow)
//! yet light (one O(N) pass; straight streamlines, since the ~9° domain sits in
//! ONE wind belt → a single uniform zonal direction).
//!
//! Per streamline (a grid row, scanned from the upwind edge): carry a moisture
//! flux `M`. Over ocean, `M` charges toward the air's capacity at the SST
//! (Clausius-Clapeyron). Over ascending land (windward uplift along the wind),
//! precipitate ∝ `M · ascent` (Smith-Barstad form) and REMOVE it from `M`
//! (conservation — the model's law); descending/leeward land sees no uplift term
//! so the dried air casts a rain shadow downwind. If `M` exceeds the local
//! capacity (air cooled climbing), the excess falls out. Anchored forms
//! (Clausius-Clapeyron, Smith-Barstad); the two rate coefficients are calibrated
//! by the windward/leeward ratio (validated ~3-10×, real ranges) — the absolute
//! unit is relative moisture (a single mm/yr calibration is a documented scale).

use crate::grid::GridF32;

/// The C1 sea-level normalised value.
pub const SEA_LEVEL_NORM: f32 = 0.5;

/// Clausius-Clapeyron saturation vapour pressure (hPa) at `t` (°C) — the air's
/// moisture-carrying CAPACITY proxy (Magnus form).
pub fn e_sat(t_c: f32) -> f32 {
    6.112 * (17.62 * t_c / (243.12 + t_c)).exp()
}

/// Meridional precipitation belt factor (≈0.15–1.0): the real zonal-mean
/// precipitation SHAPE — wet ITCZ (~0°), dry subtropical highs (~25-30° → the
/// desert belts), moderate mid-latitude westerlies (~45-55°), dry poles. Anchors
/// the FRONTAL base's latitude dependence: the desert minimum is simply this at
/// ~30°, the temperate base is this at ~45° — one profile, the slider modulates
/// the base by belt through it. Anchored on the observed zonal-mean profile, not
/// a knob.
pub fn belt_factor(lat_deg: f32) -> f32 {
    let l = lat_deg.abs();
    let itcz = (-(l / 12.0).powi(2)).exp(); // equatorial / ITCZ peak
    let midlat = (-((l - 52.0) / 15.0).powi(2)).exp(); // mid-latitude westerly peak
    (0.15 + 0.85 * itcz + 0.5 * midlat).min(1.2)
}

/// Prevailing zonal wind direction by latitude belt: `+1` eastward (streamlines
/// W→E, westerlies); `−1` westward (trade / polar easterlies). The full belt
/// table is in place for the viz latitude slider, even though one belt is active
/// per ~9° domain. Default 45° → `+1` (westerlies).
pub fn wind_zonal_dir(lat_deg: f32) -> i32 {
    match lat_deg.abs() {
        l if l < 30.0 => -1, // trade easterlies
        l if l < 60.0 => 1,  // westerlies
        _ => -1,             // polar easterlies
    }
}

/// Precipitation tunables. Anchored FORMS (Clausius-Clapeyron, Smith-Barstad);
/// the rate coefficients are calibrated by the windward/leeward ratio.
#[derive(Clone, Copy, Debug)]
pub struct PrecipParams {
    /// Ocean charge rate per cell — `M` relaxes toward `e_sat(SST)` over the
    /// fetch (maritime air saturates over a long ocean run).
    pub k_evap: f32,
    /// Orographic precip efficiency: `precip = k_oro · M · ascent_slope`
    /// (Smith-Barstad uplift fallout; `ascent_slope` is the along-wind uphill
    /// gradient, dimensionless m/m). Calibrated to SPREAD the precip across the
    /// windward slope: too high dumps ~100 % in the first steep coastal cell
    /// (thin wet rim, dry interior); ~0.5 lets moisture penetrate inland.
    pub k_oro: f32,
    /// FRONTAL / synoptic base coefficient — `precip += k_frontal · belt_factor(lat)
    /// · e_sat(T_sea(lat))` on ALL land, INDEPENDENT of local slope. The
    /// physically-missing component the orographic-only transport omits: in each
    /// circulation belt large frontal systems deposit rain over BROAD areas (the
    /// westerlies wet a temperate interior like Europe), not just on relief. This
    /// is the LATITUDINAL base — it SUBSUMES the old convective floor (the desert
    /// minimum is just `belt_factor` at the subtropical ~30° highs; the temperate
    /// base is `belt_factor` at ~45° westerlies). Orography MODULATES it (the
    /// orographic term adds on windward; the lee keeps the frontal base, not the
    /// desert floor). A SEPARATE synoptic source (moisture from elsewhere in the
    /// belt) → additive, NOT part of the orographic conservation budget.
    pub k_frontal: f32,
}

impl Default for PrecipParams {
    fn default() -> Self {
        // #165 k_oro = 0.5. NB: a k_oro reduction (→0.2) was TRIED to lengthen the
        // orographic depletion (coastal dump), and the decisive re-measure FALSIFIED
        // it: the e-folding did NOT lengthen (it sharpened on cordillera seeds),
        // because lowering the uplift extraction feeds MORE moisture to the cold
        // high coastal peak where the Clausius-Clapeyron CAPACITY CAP (e_sat(−30°C)
        // ≈ 0.5) wrings it out regardless. The depletion is the cold-coastal-peak
        // cap (PHYSICAL — air crossing a 4.65 km cold range is wrung dry, Patagonia
        // behind the Andes), not k_oro. The interior behind such ranges is
        // legitimately orographic-dry and gets its rain from the frontal base. So
        // k_oro stays 0.5; the coastal dump + frontal-base interior is the correct
        // rain-shadow regime, not a bug.
        Self { k_evap: 0.20, k_oro: 0.5, k_frontal: 0.01 }
    }
}

/// Precipitation field (relative moisture units) by the conservative transport.
/// `altitude_m(n)` converts a normalised height to metres (the vertical
/// contract); `km_per_cell` gives the physical along-wind step (the horizontal
/// contract). `sst_of_row(j)` is the sea-surface temperature per row (for evap +
/// capacity), `temperature` the air T (for the capacity cap).
pub fn compute_precipitation(
    heightmap: &GridF32,
    temperature: &GridF32,
    lat_deg: f32,
    km_per_cell: f32,
    altitude_m: impl Fn(f32) -> f32,
    params: &PrecipParams,
) -> GridF32 {
    compute_precipitation_with_budget(heightmap, temperature, lat_deg, km_per_cell, altitude_m, params).0
}

/// As [`compute_precipitation`] but also returns the OROGRAPHIC moisture budget
/// `(evap_in, exit_out, orographic_precip_sum)`. The transport conservation law
/// is `evap_in == orographic_precip_sum + exit_out` (exactly — every unit
/// evaporated into the flux either precipitates orographically or exits). The
/// CONVECTIVE baseline (`k_conv`) is a separate additive source, included in the
/// returned grid but NOT in this budget. The validation the eye cannot do.
pub fn compute_precipitation_with_budget(
    heightmap: &GridF32,
    temperature: &GridF32,
    lat_deg: f32,
    km_per_cell: f32,
    altitude_m: impl Fn(f32) -> f32,
    params: &PrecipParams,
) -> (GridF32, f64, f64, f64) {
    let (w, h) = (heightmap.width, heightmap.height);
    let mut p = GridF32::new(w, h, 0.0);
    let dir = wind_zonal_dir(lat_deg);
    let dx_m = (km_per_cell * 1000.0).max(1.0);
    let mut evap_in = 0.0f64;
    let mut exit_out = 0.0f64;
    let mut oro_precip_sum = 0.0f64;
    // FRONTAL / synoptic base (per the domain's belt; ~9° span → one value).
    let frontal_base =
        params.k_frontal * belt_factor(lat_deg) * e_sat(super::temperature::sea_level_temperature(lat_deg));

    for j in 0..h {
        // Scan along the wind, upwind → downwind.
        let order: Vec<usize> = if dir > 0 { (0..w).collect() } else { (0..w).rev().collect() };
        let mut m = 0.0f32; // carried moisture flux
        let mut prev_alt = 0.0f32; // altitude of the previous (upwind) cell
        for &i in &order {
            let n = heightmap.get(i as i32, j as i32);
            let t = temperature.get(i as i32, j as i32);
            let cap = e_sat(t);
            let alt = if n <= SEA_LEVEL_NORM { 0.0 } else { altitude_m(n).max(0.0) };
            let mut precip = 0.0f32;
            if n <= SEA_LEVEL_NORM {
                // ocean: charge toward the SST capacity (e_sat at the sea-surface T).
                let add = params.k_evap * (cap - m).max(0.0);
                m += add;
                evap_in += add as f64;
            } else {
                // land: orographic precip on the along-wind ASCENT (windward).
                let ascent = ((alt - prev_alt) / dx_m).max(0.0); // m/m
                let oro = (params.k_oro * m * ascent).min(m);
                precip += oro;
                m -= oro;
            }
            // capacity cap: air can't carry more than e_sat(T) → excess falls out.
            if m > cap {
                precip += m - cap;
                m = cap;
            }
            // `precip` so far is the OROGRAPHIC (transport) component, conserved.
            oro_precip_sum += precip as f64;
            // FRONTAL/synoptic base (separate source, additive): on land, the
            // belt's broad frontal rain — INDEPENDENT of slope. NOT part of the
            // conserved orographic budget. Subsumes the old convective floor:
            // wets the temperate interior to a moderate base (not the desert
            // floor), with orography enhancing it on windward.
            let total = if n > SEA_LEVEL_NORM {
                precip + frontal_base
            } else {
                precip
            };
            if total > 0.0 {
                p.set(i, j, total);
            }
            prev_alt = alt;
        }
        exit_out += m as f64; // moisture leaving the downwind edge of this streamline
    }
    (p, evap_in, exit_out, oro_precip_sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e_sat_rises_with_temperature() {
        assert!(e_sat(30.0) > e_sat(0.0));
        assert!(e_sat(0.0) > e_sat(-30.0));
    }

    #[test]
    fn windward_wets_leeward_dries_and_conserves() {
        // A ridge in the middle of an ocean row; westerlies (W→E).
        let w = 40usize;
        let mut hm = GridF32::new(w, 1, 0.3); // ocean baseline
        // land hump rising to a peak at i=20, descending after.
        for i in 10..30usize {
            let d = (i as i32 - 20).abs() as f32;
            hm.set(i, 0, 0.55 + (10.0 - d) * 0.03);
        }
        let t = GridF32::new(w, 1, 10.0); // uniform mild air
        let p = compute_precipitation(&hm, &t, 45.0, 0.5, |n| (n - 0.5) * 11300.0, &PrecipParams::default());
        // windward (west of peak, i 10..20) should out-rain leeward (i 20..30).
        let windward: f32 = (10..20).map(|i| p.get(i as i32, 0)).sum();
        let leeward: f32 = (21..30).map(|i| p.get(i as i32, 0)).sum();
        assert!(windward > leeward, "windward {windward} must exceed leeward {leeward}");
        assert!(leeward < windward * 0.5, "leeward should be a rain shadow");
    }

    #[test]
    fn moisture_is_conserved() {
        // The model's law: evap_in == precip.sum() + exit_out (exactly).
        let w = 60usize;
        let mut hm = GridF32::new(w, 3, 0.3);
        for j in 0..3 {
            for i in 15..45usize {
                let d = (i as i32 - 30).abs() as f32;
                hm.set(i, j, 0.55 + (15.0 - d).max(0.0) * 0.02);
            }
        }
        let t = GridF32::new(w, 3, 8.0);
        let (_p, evap_in, exit_out, oro_sum) =
            compute_precipitation_with_budget(&hm, &t, 45.0, 0.5, |n| (n - 0.5) * 11300.0, &PrecipParams::default());
        // OROGRAPHIC conservation: evap_in == orographic_precip + exit (exact).
        // (The convective baseline is a separate additive source, not in this budget.)
        let residual = (evap_in - (oro_sum + exit_out)).abs();
        assert!(
            residual < 1e-3 * evap_in.max(1.0),
            "moisture not conserved: evap_in={evap_in:.4} oro_precip={oro_sum:.4} exit={exit_out:.4} residual={residual:.2e}"
        );
        assert!(oro_sum > 0.0, "some orographic precipitation expected");
    }
}
