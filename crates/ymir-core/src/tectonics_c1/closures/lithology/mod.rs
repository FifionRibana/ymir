//! C-3 lithological heterogeneity (closures roadmap §3) — a per-cell erodibility
//! MULTIPLIER for the stream-power incision, derived CAUSALLY from the tectonic
//! state, never from noise or geometry.
//!
//! ## The source nuance (Stock & Montgomery 1999) that shapes this
//!
//! K between granitoids and metasediments is NOT significant — the contrast is
//! HARD CLASS vs SOFT CLASS, not a continuum. So a continental basement treated as
//! uniformly HARD is physically correct (crystalline + metasedimentary are both
//! hard), and only the SOFT zones need marking. Those are minority BY NATURE:
//!   - **rift** (young, `age = 0`, already stamped by the rifting closure);
//!   - **volcaniclastic** (edifice footprints, from the C-2 placement).
//! Ymir's production erosion (relief-v3) is detachment-limited — no deposition — so
//! there is no causal sedimentary-basin signal (see ADR 0001, C-3). The soft class
//! is therefore rift + volcaniclastic; everything else is hard basement at the
//! reference K. This needs NO new advected field: both signals already exist.
//!
//! ## The multiplier scheme (hard = reference)
//!
//! `1.0` = hard basement = the relief-v3 reference (so the ~80 % hard bulk erodes
//! exactly as production — no global slowdown to disentangle from the contrast).
//! Soft rock is ABOVE 1.0 (erodes more). Stock & Montgomery measure 1-5 orders of
//! magnitude softer↔harder; the practical spread is swept (see the sweep bench),
//! not asserted. `m = 0.4, n = 1` (stable base level; NOT the Kauai exponents).

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;
use crate::tectonics_c1::closures::volcanism::Edifice;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_v2::boundaries::plate_type::PlateType;

/// Version of the lithology K derivation. ⚠️ BUMP on any change that moves the K
/// field, so a lithology-enabled eroded cache invalidates. Added to the eroded key
/// only when enabled.
pub const LITHOLOGY_ALGO: u32 = 1;

/// C-3 lithology config. `enabled = false` (default) → uniform K (no field), the
/// pre-C-3 pipeline byte-identical.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LithologyConfig {
    /// Master switch. Default `false` (OFF → byte-identical).
    pub enabled: bool,
    /// Erodibility multiplier for SOFT rift rock (young, `age = 0`) relative to the
    /// hard basement reference (1.0). Stock & Montgomery mudstone/granite spread is
    /// 1-5 orders; swept in the bench to find the visible-without-destruction value.
    pub soft_multiplier: f32,
    /// Erodibility multiplier for VOLCANICLASTIC rock (edifice footprints) —
    /// intermediate (Stock & Montgomery volcaniclastic sits between granite and
    /// mudstone).
    pub volcanic_multiplier: f32,
    /// A cell is "rift / young soft" when its (continental) age is below this
    /// (rift-spawned cells are stamped `age = 0`).
    pub rift_age_threshold: f32,
}

impl LithologyConfig {
    /// `true` when the closure is OFF — used by serde `skip_serializing_if` so a
    /// disabled config drops out of the eroded cache key (byte-identical to pre-C-3).
    #[must_use]
    pub fn is_disabled(&self) -> bool {
        !self.enabled
    }
}

impl Default for LithologyConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            soft_multiplier: 10.0,    // ~1 order, mid of the swept range
            volcanic_multiplier: 3.0, // intermediate
            rift_age_threshold: 1.0,
        }
    }
}

/// Build the COARSE (tectonic-grid) K multiplier field: `1.0` hard basement
/// everywhere, `soft_multiplier` on continental rift cells (`age < threshold`).
/// Volcaniclastic is stamped later at HD ([`stamp_volcanic_k`]) from the edifice
/// footprints. Causal: rift age is stamped by the rifting closure.
#[must_use]
pub fn build_coarse_k(state: &C1State, cfg: &LithologyConfig) -> GridF32 {
    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let mut k = GridF32::new(nx, ny, 1.0);
    for j in 0..ny {
        for i in 0..nx {
            if matches!(state.plate_type.get(i, j), PlateType::Continental)
                && (state.age.get(i, j) as f32) < cfg.rift_age_threshold
            {
                k.set(i, j, cfg.soft_multiplier);
            }
        }
    }
    k
}

/// Upscale the coarse K field to HD by sampling it at the SAME `(sx, sy)` mapping
/// `upscale_with_fbm` uses for the altitude (bilinear-periodic, window-aware), so
/// the K field is perfectly registered with the terrain. Returns a `dst_w·dst_h`
/// row-major multiplier field.
#[must_use]
pub fn upscale_k_to_hd(
    coarse_k: &GridF32,
    dst_w: usize,
    dst_h: usize,
    sample_origin: [f64; 2],
    sample_size: f64,
) -> Vec<f32> {
    let (src_w, src_h) = (coarse_k.width, coarse_k.height);
    let scale_x = sample_size * src_w as f64 / dst_w as f64;
    let scale_y = sample_size * src_h as f64 / dst_h as f64;
    let origin_x = sample_origin[0] * src_w as f64;
    let origin_y = sample_origin[1] * src_h as f64;
    let mut out = vec![1.0f32; dst_w * dst_h];
    for j in 0..dst_h {
        for i in 0..dst_w {
            let sx = origin_x + i as f64 * scale_x;
            let sy = origin_y + j as f64 * scale_y;
            out[j * dst_w + i] = coarse_k.sample_bilinear_periodic(sx as f32, sy as f32);
        }
    }
    out
}

/// Stamp VOLCANICLASTIC K on the edifice basal footprints (HD pixels), the same
/// discs [`crate::tectonics_c1::closures::volcanism::apply_edifices`] rasterises.
/// `set` (max), so a volcanic footprint over a rift cell keeps the higher soft K.
pub fn stamp_volcanic_k(
    k_field: &mut [f32],
    edifices: &[Edifice],
    sample_origin: [f64; 2],
    sample_size: f64,
    km_per_hd_cell: f32,
    w: usize,
    h: usize,
    cfg: &LithologyConfig,
) {
    let (so, ss) = ([sample_origin[0] as f32, sample_origin[1] as f32], sample_size as f32);
    for e in edifices {
        let fx = (e.center_uv.0 - so[0]).rem_euclid(1.0) / ss;
        let fy = (e.center_uv.1 - so[1]).rem_euclid(1.0) / ss;
        if fx >= 1.0 || fy >= 1.0 {
            continue;
        }
        let (cx, cy) = (fx * w as f32, fy * h as f32);
        let rb = (e.basal_diameter_km * 0.5) / km_per_hd_cell;
        if rb < 1.0 {
            continue;
        }
        let (i0, i1) =
            ((cx - rb).floor().max(0.0) as usize, ((cx + rb).ceil() as usize).min(w - 1));
        let (j0, j1) =
            ((cy - rb).floor().max(0.0) as usize, ((cy + rb).ceil() as usize).min(h - 1));
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d <= rb {
                    let k = j * w + i;
                    k_field[k] = k_field[k].max(cfg.volcanic_multiplier);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upscale_k_bilinear_registers_and_hard_is_one() {
        // A coarse field with a soft patch → HD sampling puts soft where the patch
        // is, hard (1.0) elsewhere, and interpolates smoothly at the edge.
        let mut ck = GridF32::new(8, 8, 1.0);
        ck.set(4, 4, 10.0);
        let hd = upscale_k_to_hd(&ck, 64, 64, [0.0, 0.0], 1.0);
        // The soft PEAK sits at grid point (4,4); with scale 8/64 = 0.125 cells/px it
        // lands exactly on HD pixel 32 (= 4.0 / 0.125), where bilinear returns the full
        // 10.0. (Neighbouring pixels interpolate down toward the hard 1.0 bulk.)
        assert!(hd[32 * 64 + 32] > 9.0, "soft patch peak upscales high: {}", hd[32 * 64 + 32]);
        // a far corner is hard.
        assert!((hd[0] - 1.0).abs() < 0.01, "hard basement stays 1.0: {}", hd[0]);
    }
}
