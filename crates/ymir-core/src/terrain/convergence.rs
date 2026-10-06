//! ADR Finding 148-M — flow convergence on a height field: the contour (planform) curvature
//! `C = div(∇z / |∇z|)` of the field smoothed at a declared scale.
//!
//! `C > 0` where the flow lines meet (a hollow, a valley), `C < 0` where they part (a nose, a ridge),
//! `C ≈ 0` on a planar slope. A channel head needs convergence (Montgomery & Dietrich's hollows; the
//! reviewer's citation, to check). At the cell scale the curvature is noise, hence the Gaussian.

use rayon::prelude::*;

/// Below this gradient (m/m) the flow direction is undefined: `C = 0`.
pub const FLAT_GRADIENT: f32 = 1e-4;

/// A separable Gaussian blur of standard deviation `sigma_cells` (radius ⌈3σ⌉), edges clamped.
/// `sigma_cells <= 0` returns the field unchanged.
pub fn gaussian_smooth(z: &[f32], w: usize, h: usize, sigma_cells: f32) -> Vec<f32> {
    if sigma_cells <= 0.0 {
        return z.to_vec();
    }
    let r = (3.0 * sigma_cells).ceil() as i64;
    let mut k: Vec<f32> = (-r..=r).map(|i| (-(i * i) as f32 / (2.0 * sigma_cells * sigma_cells)).exp()).collect();
    let s: f32 = k.iter().sum();
    k.iter_mut().for_each(|v| *v /= s);
    let mut tmp = vec![0f32; w * h];
    tmp.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        for (x, out) in row.iter_mut().enumerate() {
            let mut acc = 0f32;
            for (t, kv) in k.iter().enumerate() {
                let xx = (x as i64 + t as i64 - r).clamp(0, w as i64 - 1) as usize;
                acc += kv * z[y * w + xx];
            }
            *out = acc;
        }
    });
    let mut out = vec![0f32; w * h];
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        for (x, o) in row.iter_mut().enumerate() {
            let mut acc = 0f32;
            for (t, kv) in k.iter().enumerate() {
                let yy = (y as i64 + t as i64 - r).clamp(0, h as i64 - 1) as usize;
                acc += kv * tmp[yy * w + x];
            }
            *o = acc;
        }
    });
    out
}

/// The contour convergence `C = div(∇z / |∇z|)` (km⁻¹) of the field `z_m` (metres, cells of `cell_m`
/// metres) smoothed at `sigma_cells`. Central differences, edges clamped; `C = 0` where the
/// smoothed gradient is below [`FLAT_GRADIENT`].
pub fn contour_convergence(z_m: &[f32], w: usize, h: usize, sigma_cells: f32, cell_m: f32) -> Vec<f32> {
    let s = gaussian_smooth(z_m, w, h, sigma_cells);
    let at = |v: &[f32], x: i64, y: i64| v[(y.clamp(0, h as i64 - 1) as usize) * w + x.clamp(0, w as i64 - 1) as usize];
    // the unit gradient, interleaved (ux, uy)
    let mut u = vec![[0f32; 2]; w * h];
    u.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let y = y as i64;
        for (x, o) in row.iter_mut().enumerate() {
            let x = x as i64;
            let gx = (at(&s, x + 1, y) - at(&s, x - 1, y)) / (2.0 * cell_m);
            let gy = (at(&s, x, y + 1) - at(&s, x, y - 1)) / (2.0 * cell_m);
            let g = (gx * gx + gy * gy).sqrt();
            *o = if g < FLAT_GRADIENT { [0.0, 0.0] } else { [gx / g, gy / g] };
        }
    });
    drop(s);
    let flat: Vec<bool> = u.par_iter().map(|v| v[0] == 0.0 && v[1] == 0.0).collect();
    let ua = |x: i64, y: i64| u[(y.clamp(0, h as i64 - 1) as usize) * w + x.clamp(0, w as i64 - 1) as usize];
    let mut c = vec![0f32; w * h];
    c.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let yi = y as i64;
        for (x, o) in row.iter_mut().enumerate() {
            if flat[y * w + x] {
                continue;
            }
            let xi = x as i64;
            let d = (ua(xi + 1, yi)[0] - ua(xi - 1, yi)[0]) / (2.0 * cell_m) + (ua(xi, yi + 1)[1] - ua(xi, yi - 1)[1]) / (2.0 * cell_m);
            *o = d * 1000.0;
        }
    });
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cone pit converges (C > 0), a cone peak diverges (C < 0), and the negative control, a plane,
    /// reads C ≈ 0. A V-valley converges on its axis.
    #[test]
    fn convergence_signs_on_synthetic_fields() {
        let (w, h) = (64usize, 64usize);
        let cell_m = 50.0;
        let r = |x: usize, y: usize| (((x as f32 - 32.0).powi(2) + (y as f32 - 32.0).powi(2)).sqrt()) * cell_m;
        let pit: Vec<f32> = (0..w * h).map(|i| r(i % w, i / w)).collect();
        let peak: Vec<f32> = pit.iter().map(|v| -v).collect();
        let plane: Vec<f32> = (0..w * h).map(|i| 0.1 * (i % w) as f32 * cell_m + 0.05 * (i / w) as f32 * cell_m).collect();
        let valley: Vec<f32> = (0..w * h).map(|i| 0.3 * ((i % w) as f32 - 32.0).abs() * cell_m + 0.05 * (i / w) as f32 * cell_m).collect();
        let at = |c: &[f32], x: usize, y: usize| c[y * w + x];
        let cp = contour_convergence(&pit, w, h, 2.0, cell_m);
        let ck = contour_convergence(&peak, w, h, 2.0, cell_m);
        let cl = contour_convergence(&plane, w, h, 2.0, cell_m);
        let cv = contour_convergence(&valley, w, h, 2.0, cell_m);
        // 8 cells (400 m) from the centre: 1/r = 2.5 km⁻¹
        assert!((at(&cp, 40, 32) - 2.5).abs() < 0.3, "pit {}", at(&cp, 40, 32));
        assert!((at(&ck, 40, 32) + 2.5).abs() < 0.3, "peak {}", at(&ck, 40, 32));
        let max_plane = (8..56).flat_map(|y| (8..56).map(move |x| (x, y))).map(|(x, y)| at(&cl, x, y).abs()).fold(0f32, f32::max);
        assert!(max_plane < 1e-3, "negative control: a plane reads C ≈ 0 ({max_plane})");
        assert!(at(&cv, 32, 32) > 1.0, "the valley's axis converges ({})", at(&cv, 32, 32));
        assert!(at(&cv, 20, 32).abs() < 0.05, "its planar flank does not ({})", at(&cv, 20, 32));
    }

    /// The Gaussian keeps a constant field and the mean of a field.
    #[test]
    fn gaussian_keeps_constants() {
        let z = vec![3.0f32; 32 * 16];
        let s = gaussian_smooth(&z, 32, 16, 2.0);
        assert!(s.iter().all(|v| (v - 3.0).abs() < 1e-5));
    }
}
