//! Metric `height` raster for the v1 `.ymir` container — the vertical pendant
//! of the vector coastline trace.
//!
//! Converts the normalized eroded field to METRES via the ONE vertical contract
//! ([`c1_altitude_norm_to_metres`]), anchored so `0 m` falls on the SAME
//! sea-level normalized constant the coastline is traced at
//! ([`vector::SEA_LEVEL_NORM`]) — so `0 m == the coastline` by construction,
//! even if that constant ever moves off the contract's built-in `0.5`. The
//! field is then quantised to `u16` linearly over its TRUE metric range
//! `[min_m, max_m]`:
//!
//! ```text
//!   encode:  code = round((m − min_m) / (max_m − min_m) · 65535)   (clamped)
//!   decode:  m    = min_m + (code / 65535) · (max_m − min_m)
//! ```
//!
//! No second vertical formula and no second sea-level constant: the anchor
//! subtracts the contract's own value at [`vector::SEA_LEVEL_NORM`].

use crate::export::vector::SEA_LEVEL_NORM;
use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;

/// Metric altitude (metres) of a normalized cell, re-anchored so that
/// [`SEA_LEVEL_NORM`] maps to exactly `0 m` (the coastline). Reuses the single
/// vertical contract and subtracts its value at sea level — no second formula.
/// (With the current constants `SEA_LEVEL_NORM == 0.5`, so the offset is `0`;
/// the subtraction keeps the two in lockstep if either ever changes.)
pub fn altitude_norm_to_metres_sea_anchored(norm: f32, ss: &SteinSteinParams) -> f32 {
    c1_altitude_norm_to_metres(norm, ss) - c1_altitude_norm_to_metres(SEA_LEVEL_NORM, ss)
}

/// Quantised metric `height` raster: the row-major `u16` codes plus the true
/// metric range `[min_m, max_m]` the codes map onto (see the module doc for the
/// exact linear encode/decode).
pub struct MetricHeight {
    /// Row-major `u16` codes, linear over `[min_m, max_m]`.
    pub codes: Vec<u16>,
    /// Lowest metric altitude in the field (the `code = 0` endpoint). Negative
    /// where the field dips below sea level; the honest `max_depth_m`.
    pub min_m: f32,
    /// Highest metric altitude in the field (the `code = 65535` endpoint); the
    /// honest `max_elevation_m`.
    pub max_m: f32,
}

impl MetricHeight {
    /// One quantisation step in metres: `(max_m − min_m) / 65535`.
    pub fn step_m(&self) -> f32 {
        (self.max_m - self.min_m) / u16::MAX as f32
    }

    /// Decode a `u16` code back to metres (inverse of the encode; see module doc).
    pub fn decode(&self, code: u16) -> f32 {
        self.min_m + (code as f32 / u16::MAX as f32) * (self.max_m - self.min_m)
    }
}

/// Build the quantised metric `height` raster for the normalized `eroded` field.
/// Deterministic: same field + params → same codes and range.
pub fn metric_height_u16(eroded_norm: &GridF32, ss: &SteinSteinParams) -> MetricHeight {
    let metres: Vec<f32> =
        eroded_norm.data.iter().map(|&n| altitude_norm_to_metres_sea_anchored(n, ss)).collect();

    let (mut min_m, mut max_m) = (f32::INFINITY, f32::NEG_INFINITY);
    for &m in &metres {
        min_m = min_m.min(m);
        max_m = max_m.max(m);
    }
    // Degenerate (perfectly flat) field: avoid a zero span / divide-by-zero.
    let span = (max_m - min_m).max(f32::EPSILON);

    let codes: Vec<u16> = metres
        .iter()
        .map(|&m| {
            (((m - min_m) / span) * u16::MAX as f32).round().clamp(0.0, u16::MAX as f32) as u16
        })
        .collect();

    MetricHeight { codes, min_m, max_m }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::terrain::contour::marching_squares;

    /// A central landmass on a sea background (norm 0.2), rising to 0.9 — the
    /// same shape the coastline test uses, so its metric coastline is a loop.
    fn island_field(w: usize, h: usize) -> GridF32 {
        let (cx, cy) = ((w as f32 - 1.0) * 0.5, (h as f32 - 1.0) * 0.5);
        let mut data = vec![0.2f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let d = ((x as f32 - cx).powi(2) + (y as f32 - cy).powi(2)).sqrt();
                data[y * w + x] = (0.2 + (1.0 - d / 12.0).max(0.0) * 0.6).min(0.9);
            }
        }
        GridF32::from_vec(w, h, data)
    }

    /// Bilinear sample of a row-major field at fractional cell coords — the same
    /// coordinate space marching-squares emits (`x`=col, `y`=row).
    fn sample_bilinear(data: &[f32], w: usize, h: usize, x: f32, y: f32) -> f32 {
        let xf = x.clamp(0.0, (w - 1) as f32);
        let yf = y.clamp(0.0, (h - 1) as f32);
        let (x0, y0) = (xf.floor() as usize, yf.floor() as usize);
        let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
        let (tx, ty) = (xf - x0 as f32, yf - y0 as f32);
        let v00 = data[y0 * w + x0];
        let v10 = data[y0 * w + x1];
        let v01 = data[y1 * w + x0];
        let v11 = data[y1 * w + x1];
        let a = v00 * (1.0 - tx) + v10 * tx;
        let b = v01 * (1.0 - tx) + v11 * tx;
        a * (1.0 - ty) + b * ty
    }

    /// Decoding a `u16` code recovers the contract metres within one quant step.
    #[test]
    fn decode_matches_contract_within_one_step() {
        let ss = SteinSteinParams::default();
        let (w, h) = (17usize, 13usize);
        // A smooth norm gradient spanning sea and land.
        let mut data = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                data[y * w + x] = 0.1 + 0.8 * (x as f32 / (w - 1) as f32);
            }
        }
        let field = GridF32::from_vec(w, h, data);
        let mh = metric_height_u16(&field, &ss);
        let step = mh.step_m();

        for (i, &n) in field.data.iter().enumerate() {
            let expected = altitude_norm_to_metres_sea_anchored(n, &ss);
            let decoded = mh.decode(mh.codes[i]);
            assert!(
                (decoded - expected).abs() <= step + 1e-3,
                "cell {i}: decoded {decoded} vs contract {expected} exceeds one step {step}"
            );
        }
    }

    /// THE regression: the metric height sampled along the traced coastline is
    /// ~0 m (within one quant step), i.e. the metric zero and the coastline
    /// isoline are the same surface by construction.
    #[test]
    fn coastline_samples_at_zero_metres() {
        let ss = SteinSteinParams::default();
        let (w, h) = (48usize, 48usize);
        let field = island_field(w, h);

        // Coastline traced exactly as the container does (same constant).
        let polylines = marching_squares(&field, SEA_LEVEL_NORM);
        assert!(!polylines.is_empty(), "island must have a coastline to sample");

        let mh = metric_height_u16(&field, &ss);
        // Decoded metric height field, so we sample the SHIPPED quantised values.
        let decoded: Vec<f32> = mh.codes.iter().map(|&c| mh.decode(c)).collect();
        let tol = mh.step_m();

        let mut n_pts = 0usize;
        for pl in &polylines {
            for &(x, y) in pl {
                let m = sample_bilinear(&decoded, w, h, x, y);
                assert!(
                    m.abs() <= tol + 1e-3,
                    "coastline point ({x},{y}) sampled {m} m, expected ~0 within {tol} m"
                );
                n_pts += 1;
            }
        }
        assert!(n_pts > 4, "coastline must have vertices to make the check meaningful");
    }

    /// The manifest extrema equal the field's true metric range.
    #[test]
    fn extrema_equal_field_range() {
        let ss = SteinSteinParams::default();
        let field = island_field(48, 48);
        let mh = metric_height_u16(&field, &ss);

        let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
        for &n in &field.data {
            let m = altitude_norm_to_metres_sea_anchored(n, &ss);
            lo = lo.min(m);
            hi = hi.max(m);
        }
        assert!((mh.min_m - lo).abs() < 1e-3, "max_depth_m must equal the field min");
        assert!((mh.max_m - hi).abs() < 1e-3, "max_elevation_m must equal the field max");
        // Sea background (0.2) is below sea level, land (0.9) above → signs honest.
        assert!(mh.min_m < 0.0 && mh.max_m > 0.0, "range must straddle sea level");
    }
}
