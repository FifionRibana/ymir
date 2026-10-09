//! ADR Finding 161 -- the measures the cascade shares with its benches and its window: the altitude spectrum by
//! octave of wavelength (the calibration of k_L reads its finest octave), and the D5 mask carried from level to level.
//!
//! **The spectrum** (F161-C3): the 2D periodogram of the altitude (radix-2 FFT, the sea at 0 m, no taper: every grid is
//! bordered by sea, so it is periodic at 0), summed by octave of wavelength — octave j holds the wavenumbers whose
//! wavelength is in [2^(j+1), 2^(j+2)) cells, octave 0 the 2–4 cell band, the « new » band of a level. By Parseval the
//! octaves' variances add up to the field's (DC excluded); each is reported as an RMS in metres **per land cell** (the
//! variance divided by the land share). After Perron, Kirchner & Dietrich 2008 (JGR, « spectral signatures »; not
//! fetched), without their detrending and windowing.
//!
//! F160's band-pass « bands » (differences of box blurs) are NOT this: their passbands overlap and their sinc lobes
//! move power between bands (a 3-cell sine read mostly in F160's band 1). F161 replaced them before any measurement.

use rayon::prelude::*;

/// In-place radix-2 complex FFT (forward), `re.len()` a power of two.
fn fft(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let ang = -2.0 * std::f64::consts::PI / len as f64;
        let (wr, wi) = (ang.cos(), ang.sin());
        for start in (0..n).step_by(len) {
            let (mut cr, mut ci) = (1.0f64, 0.0f64);
            for k in 0..len / 2 {
                let (a, b) = (start + k, start + k + len / 2);
                let (xr, xi) = (re[b] * cr - im[b] * ci, re[b] * ci + im[b] * cr);
                re[b] = re[a] - xr;
                im[b] = im[a] - xi;
                re[a] += xr;
                im[a] += xi;
                let t = cr * wr - ci * wi;
                ci = cr * wi + ci * wr;
                cr = t;
            }
        }
        len <<= 1;
    }
}

/// The spectrum by octave: element j is the RMS (m, per land cell) of the wavelengths [2^(j+1), 2^(j+2)) cells, for
/// every octave up to n cells. `zm` is the field (m) with the sea at 0 m; `land_cells` the normalising land count.
pub fn octave_rms(zm: &[f32], n: usize, land_cells: usize) -> Vec<f32> {
    assert!(n.is_power_of_two(), "the spectrum needs a power-of-two grid");
    // rows, then columns
    let mut rows: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .map(|y| {
            let mut re: Vec<f64> = zm[y * n..(y + 1) * n].iter().map(|&v| v as f64).collect();
            let mut im = vec![0f64; n];
            fft(&mut re, &mut im);
            (re, im)
        })
        .collect();
    let cols: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .map(|x| {
            let mut re: Vec<f64> = (0..n).map(|y| rows[y].0[x]).collect();
            let mut im: Vec<f64> = (0..n).map(|y| rows[y].1[x]).collect();
            fft(&mut re, &mut im);
            (re, im)
        })
        .collect();
    rows.clear();
    let octaves = (n.trailing_zeros() as usize).saturating_sub(1);
    let mut var = vec![0f64; octaves];
    let norm = (n as f64).powi(4);
    for (kx, (re, im)) in cols.iter().enumerate() {
        let fx = if kx < n / 2 { kx as f64 } else { kx as f64 - n as f64 };
        for ky in 0..n {
            let fy = if ky < n / 2 { ky as f64 } else { ky as f64 - n as f64 };
            let kk = (fx * fx + fy * fy).sqrt();
            if kk == 0.0 {
                continue;
            }
            let lam = n as f64 / kk; // wavelength in cells
            if lam < 2.0 {
                continue;
            }
            let j = (lam / 2.0).log2().floor() as usize;
            if j < octaves {
                var[j] += (re[ky] * re[ky] + im[ky] * im[ky]) / norm;
            }
        }
    }
    let share = land_cells.max(1) as f64 / (n * n) as f64;
    var.iter().map(|v| (v / share).sqrt() as f32).collect()
}

/// The wavelengths of octave j at a cell of `cell_km`: (shortest, longest) in km.
pub fn octave_km(j: usize, cell_km: f32) -> (f32, f32) {
    (cell_km * (2u32 << j) as f32, cell_km * (4u32 << j) as f32)
}

/// A boolean mask carried from an n-grid to its 2n bicubic successor: nearest neighbour (the fine pixel 2i+δ reads
/// the coarse pixel i).
pub fn upsample_mask2(m: &[bool], n: usize) -> Vec<bool> {
    let n2 = 2 * n;
    (0..n2 * n2).map(|k| m[(k / n2 / 2) * n + (k % n2) / 2]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A sine of wavelength λ cells lands in its octave and nowhere else (3 cells in octave 0, 12 in octave 2), the
    /// octaves add up to the variance (Parseval), and a flat field has no power (negative control).
    #[test]
    fn a_sine_lands_in_its_octave() {
        let n = 128;
        let sine = |cycles: f32| -> Vec<f32> {
            (0..n * n).map(|k| 100.0 * ((k % n) as f32 * std::f32::consts::TAU * cycles / n as f32).sin()).collect()
        };
        let s3 = octave_rms(&sine(43.0), n, n * n); // 128 / 43 ≈ 2.98 cells
        let s12 = octave_rms(&sine(11.0), n, n * n); // ≈ 11.6 cells
        let rms = 100.0 / 2f32.sqrt();
        assert!((s3[0] - rms).abs() < 0.5 && s3.iter().skip(1).all(|&v| v < 1e-3), "{s3:?}");
        assert!((s12[2] - rms).abs() < 0.5 && s12.iter().enumerate().all(|(j, &v)| j == 2 || v < 1e-3), "{s12:?}");
        let flat = octave_rms(&vec![500.0; n * n], n, n * n);
        assert!(flat.iter().all(|&v| v < 1e-3), "{flat:?}");
    }
}
