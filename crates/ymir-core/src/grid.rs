//! Core 2D grid type used throughout the pipeline.
//!
//! [`GridF32`] wraps a flat `Vec<f32>` with width/height metadata and provides
//! indexed access, bilinear interpolation, gradient computation, and basic
//! statistics. Every phase of the pipeline reads and writes `GridF32` instances.

use std::path::Path;

/// Wrap coordinate to [0, size) range (toroidal).
#[inline]
fn wrap(x: i32, size: usize) -> usize {
    ((x % size as i32) + size as i32) as usize % size
}

/// A 2D grid of f32 values, stored in row-major order (y * width + x).
///
/// Convention: (0,0) is the top-left corner. X increases rightward, Y increases
/// downward. This matches image conventions and simplifies PNG I/O. The Bevy
/// visualization flips Y for rendering.
#[derive(Clone)]
pub struct GridF32 {
    pub width: usize,
    pub height: usize,
    pub data: Vec<f32>,
}

impl GridF32 {
    /// Create a new grid filled with a constant value.
    pub fn new(width: usize, height: usize, fill: f32) -> Self {
        Self { width, height, data: vec![fill; width * height] }
    }

    /// Create a grid from an existing data vector.
    ///
    /// # Panics
    /// Panics if `data.len() != width * height`.
    pub fn from_vec(width: usize, height: usize, data: Vec<f32>) -> Self {
        assert_eq!(
            data.len(),
            width * height,
            "GridF32::from_vec: data length {} does not match {}x{}",
            data.len(),
            width,
            height
        );
        Self { width, height, data }
    }

    /// Total number of cells.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the grid is empty (zero dimension).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    // ── Indexed access ───────────────────────────────────────────────────

    /// Convert (x, y) to a flat index. No bounds checking.
    #[inline]
    pub fn index(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    /// Get the value at (x, y). Returns 0.0 if out of bounds.
    #[inline]
    pub fn get(&self, x: i32, y: i32) -> f32 {
        if x < 0 || y < 0 || x >= self.width as i32 || y >= self.height as i32 {
            return 0.0;
        }
        self.data[y as usize * self.width + x as usize]
    }

    /// Set the value at (x, y). No-op if out of bounds.
    #[inline]
    pub fn set(&mut self, x: usize, y: usize, value: f32) {
        if x < self.width && y < self.height {
            self.data[y * self.width + x] = value;
        }
    }

    // ── Interpolation ────────────────────────────────────────────────────

    /// Bilinear interpolation at fractional coordinates (fx, fy).
    ///
    /// Coordinates are in pixel space: (0.0, 0.0) is the center of the top-left
    /// pixel, (width-1, height-1) is the center of the bottom-right pixel.
    /// Values outside the grid are clamped to the nearest edge pixel.
    pub fn sample_bilinear(&self, fx: f32, fy: f32) -> f32 {
        // Clamp to valid range
        let fx = fx.clamp(0.0, (self.width - 1) as f32);
        let fy = fy.clamp(0.0, (self.height - 1) as f32);

        let x0 = fx.floor() as usize;
        let y0 = fy.floor() as usize;
        let x1 = (x0 + 1).min(self.width - 1);
        let y1 = (y0 + 1).min(self.height - 1);

        // Fractional part within the cell
        let tx = fx - x0 as f32;
        let ty = fy - y0 as f32;

        // Four corner values
        let c00 = self.data[y0 * self.width + x0];
        let c10 = self.data[y0 * self.width + x1];
        let c01 = self.data[y1 * self.width + x0];
        let c11 = self.data[y1 * self.width + x1];

        // Interpolate horizontally, then vertically
        let top = c00 + (c10 - c00) * tx;
        let bot = c01 + (c11 - c01) * tx;
        top + (bot - top) * ty
    }

    /// Sample the grid at normalized UV coordinates (u, v) in [0, 1].
    ///
    /// (0,0) maps to the top-left corner, (1,1) to the bottom-right.
    /// Uses bilinear interpolation.
    pub fn sample_uv(&self, u: f32, v: f32) -> f32 {
        let fx = u * (self.width - 1) as f32;
        let fy = v * (self.height - 1) as f32;
        self.sample_bilinear(fx, fy)
    }

    // ── Gradient ─────────────────────────────────────────────────────────

    /// Compute the gradient (∂h/∂x, ∂h/∂y) at integer coordinates using
    /// central differences (Sobel-like, but unweighted for simplicity).
    ///
    /// Returns `(gx, gy)` where `gx` is the slope in the x direction and
    /// `gy` is the slope in the y direction. The magnitude `sqrt(gx²+gy²)`
    /// is the steepness; the direction `atan2(gy, gx)` points uphill.
    ///
    /// At grid edges, forward or backward differences are used instead.
    pub fn gradient_at(&self, x: usize, y: usize) -> (f32, f32) {
        let ix = x as i32;
        let iy = y as i32;

        // Central difference in x: (h(x+1) - h(x-1)) / 2
        let gx = (self.get(ix + 1, iy) - self.get(ix - 1, iy)) * 0.5;

        // Central difference in y: (h(y+1) - h(y-1)) / 2
        let gy = (self.get(ix, iy + 1) - self.get(ix, iy - 1)) * 0.5;

        (gx, gy)
    }

    /// Compute the gradient at fractional coordinates using bilinear samples.
    ///
    /// This is useful for erosion droplets that sit between grid cells.
    /// The step size `eps` controls the finite difference spacing (default ~1.0).
    pub fn gradient_at_f(&self, fx: f32, fy: f32, eps: f32) -> (f32, f32) {
        let gx =
            (self.sample_bilinear(fx + eps, fy) - self.sample_bilinear(fx - eps, fy)) / (2.0 * eps);
        let gy =
            (self.sample_bilinear(fx, fy + eps) - self.sample_bilinear(fx, fy - eps)) / (2.0 * eps);
        (gx, gy)
    }

    // ── Periodic (toroidal) methods ────────────────────────────────────

    /// Get value with toroidal wrapping.
    #[inline]
    pub fn get_periodic(&self, x: i32, y: i32) -> f32 {
        let wx = wrap(x, self.width);
        let wy = wrap(y, self.height);
        self.data[wy * self.width + wx]
    }

    /// Set value with toroidal wrapping.
    #[inline]
    pub fn set_periodic(&mut self, x: i32, y: i32, value: f32) {
        let wx = wrap(x, self.width);
        let wy = wrap(y, self.height);
        self.data[wy * self.width + wx] = value;
    }

    /// Bilinear interpolation with toroidal wrapping.
    pub fn sample_bilinear_periodic(&self, fx: f32, fy: f32) -> f32 {
        let w = self.width as f32;
        let h = self.height as f32;

        let fx = ((fx % w) + w) % w;
        let fy = ((fy % h) + h) % h;

        let x0 = fx.floor() as usize;
        let y0 = fy.floor() as usize;
        let x1 = (x0 + 1) % self.width;
        let y1 = (y0 + 1) % self.height;

        let tx = fx - x0 as f32;
        let ty = fy - y0 as f32;

        let c00 = self.data[y0 * self.width + x0];
        let c10 = self.data[y0 * self.width + x1];
        let c01 = self.data[y1 * self.width + x0];
        let c11 = self.data[y1 * self.width + x1];

        let top = c00 + (c10 - c00) * tx;
        let bot = c01 + (c11 - c01) * tx;
        top + (bot - top) * ty
    }

    /// Gradient at integer coordinates with periodic wrapping.
    pub fn gradient_at_periodic(&self, x: usize, y: usize) -> (f32, f32) {
        let w = self.width;
        let h = self.height;
        let xp = if x + 1 < w { x + 1 } else { 0 };
        let xm = if x > 0 { x - 1 } else { w - 1 };
        let yp = if y + 1 < h { y + 1 } else { 0 };
        let ym = if y > 0 { y - 1 } else { h - 1 };

        let gx = (self.data[y * w + xp] - self.data[y * w + xm]) * 0.5;
        let gy = (self.data[yp * w + x] - self.data[ym * w + x]) * 0.5;
        (gx, gy)
    }

    /// Gradient at fractional coordinates with periodic wrapping.
    pub fn gradient_at_f_periodic(&self, fx: f32, fy: f32, eps: f32) -> (f32, f32) {
        let gx = (self.sample_bilinear_periodic(fx + eps, fy)
            - self.sample_bilinear_periodic(fx - eps, fy))
            / (2.0 * eps);
        let gy = (self.sample_bilinear_periodic(fx, fy + eps)
            - self.sample_bilinear_periodic(fx, fy - eps))
            / (2.0 * eps);
        (gx, gy)
    }

    // ── Statistics ───────────────────────────────────────────────────────

    /// Minimum value in the grid.
    pub fn min(&self) -> f32 {
        self.data.iter().copied().reduce(f32::min).unwrap_or(0.0)
    }

    /// Maximum value in the grid.
    pub fn max(&self) -> f32 {
        self.data.iter().copied().reduce(f32::max).unwrap_or(0.0)
    }

    /// Mean value across all cells.
    pub fn mean(&self) -> f32 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.data.iter().sum::<f32>() / self.data.len() as f32
    }

    // ── Raw binary I/O ────────────────────────────────────────────────────

    /// Save the grid as raw little-endian f32 bytes. Lossless. Delegates to the
    /// shared codec (`export::raw`) so every `.raw` shares the §9 byte layout.
    pub fn save_raw(&self, path: &Path) -> Result<(), String> {
        crate::export::raw::save_f32(path, &self.data)
    }

    /// Load a grid from raw little-endian f32 bytes (via the shared codec).
    pub fn load_raw(path: &Path, width: usize, height: usize) -> Result<Self, String> {
        let data = crate::export::raw::load_f32(path, width * height)?;
        Ok(Self { width, height, data })
    }

    /// Load from raw f32 if available, otherwise fall back to 16-bit PNG.
    pub fn load_raw_or_png(
        raw_path: &Path,
        png_path: &Path,
        width: usize,
        height: usize,
    ) -> Result<Self, String> {
        if raw_path.exists() {
            Self::load_raw(raw_path, width, height)
        } else if png_path.exists() {
            Self::load_png(png_path)
        } else {
            Err(format!("Neither {} nor {} found", raw_path.display(), png_path.display()))
        }
    }

    // ── PNG I/O ──────────────────────────────────────────────────────────

    /// Load a grayscale PNG as a GridF32.
    ///
    /// Pixel values are normalized to [0, 1] regardless of the source bit
    /// depth (u8 → /255, u16 → /65535). For heightmaps, the caller scales
    /// by max_elevation after loading.
    pub fn load_png(path: &Path) -> Result<Self, String> {
        let img =
            image::open(path).map_err(|e| format!("Failed to open {}: {}", path.display(), e))?;

        let gray = img.to_luma16();
        let width = gray.width() as usize;
        let height = gray.height() as usize;

        let data: Vec<f32> = gray.pixels().map(|p| p.0[0] as f32 / 65535.0).collect();

        Ok(Self { width, height, data })
    }

    /// Save the grid as a 16-bit grayscale PNG.
    ///
    /// Values are clamped to [0, 1] and mapped to u16 [0, 65535].
    pub fn save_png_u16(&self, path: &Path) -> Result<(), String> {
        let mut img = image::ImageBuffer::<image::Luma<u16>, Vec<u16>>::new(
            self.width as u32,
            self.height as u32,
        );

        for y in 0..self.height {
            for x in 0..self.width {
                let val = self.data[y * self.width + x].clamp(0.0, 1.0);
                let u16_val = (val * 65535.0).round() as u16;
                img.put_pixel(x as u32, y as u32, image::Luma([u16_val]));
            }
        }

        img.save(path).map_err(|e| format!("Failed to save {}: {}", path.display(), e))
    }

    /// Save the grid as an 8-bit grayscale PNG.
    ///
    /// Values are clamped to [0, 1] and mapped to u8 [0, 255].
    /// Lower precision than u16, but compatible with more tools.
    pub fn save_png_u8(&self, path: &Path) -> Result<(), String> {
        let mut img = image::ImageBuffer::<image::Luma<u8>, Vec<u8>>::new(
            self.width as u32,
            self.height as u32,
        );

        for y in 0..self.height {
            for x in 0..self.width {
                let val = self.data[y * self.width + x].clamp(0.0, 1.0);
                let u8_val = (val * 255.0).round() as u8;
                img.put_pixel(x as u32, y as u32, image::Luma([u8_val]));
            }
        }

        img.save(path).map_err(|e| format!("Failed to save {}: {}", path.display(), e))
    }

    /// Separable Gaussian blur with toroidal (wrapping) boundary conditions.
    ///
    /// Sigma is in grid cells. Kernel radius = ceil(3 * sigma).
    /// Returns a new grid; the original is unchanged.
    pub fn gaussian_blur(&self, sigma: f32) -> Self {
        if sigma <= 0.0 {
            return self.clone();
        }

        let radius = (3.0 * sigma).ceil() as usize;
        let kernel_size = 2 * radius + 1;

        let mut kernel: Vec<f32> = (0..kernel_size)
            .map(|i| {
                let x = i as f32 - radius as f32;
                (-x * x / (2.0 * sigma * sigma)).exp()
            })
            .collect();
        let sum: f32 = kernel.iter().sum();
        for k in kernel.iter_mut() {
            *k /= sum;
        }

        // Horizontal pass
        let mut temp = Self::new(self.width, self.height, 0.0);
        for y in 0..self.height {
            for x in 0..self.width {
                let mut val = 0.0f32;
                for (ki, &w) in kernel.iter().enumerate() {
                    let sx = (x as i32 + ki as i32 - radius as i32).rem_euclid(self.width as i32)
                        as usize;
                    val += self.data[y * self.width + sx] * w;
                }
                temp.data[y * self.width + x] = val;
            }
        }

        // Vertical pass
        let mut result = Self::new(self.width, self.height, 0.0);
        for y in 0..self.height {
            for x in 0..self.width {
                let mut val = 0.0f32;
                for (ki, &w) in kernel.iter().enumerate() {
                    let sy = (y as i32 + ki as i32 - radius as i32).rem_euclid(self.height as i32)
                        as usize;
                    val += temp.data[sy * self.width + x] * w;
                }
                result.data[y * self.width + x] = val;
            }
        }

        result
    }
}

impl std::fmt::Debug for GridF32 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GridF32({}x{}, min={:.4}, max={:.4}, mean={:.4})",
            self.width,
            self.height,
            self.min(),
            self.max(),
            self.mean()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_and_access() {
        let grid = GridF32::new(4, 3, 0.5);
        assert_eq!(grid.len(), 12);
        assert_eq!(grid.get(0, 0), 0.5);
        assert_eq!(grid.get(3, 2), 0.5);
        // Out of bounds returns 0.0
        assert_eq!(grid.get(-1, 0), 0.0);
        assert_eq!(grid.get(4, 0), 0.0);
    }

    #[test]
    fn test_set() {
        let mut grid = GridF32::new(4, 3, 0.0);
        grid.set(2, 1, 1.0);
        assert_eq!(grid.get(2, 1), 1.0);
        assert_eq!(grid.get(0, 0), 0.0);
    }

    #[test]
    fn test_bilinear_center() {
        // A 2x2 grid with known values
        let grid = GridF32::from_vec(2, 2, vec![0.0, 1.0, 1.0, 0.0]);
        // Center of the grid: average of all four corners
        let center = grid.sample_bilinear(0.5, 0.5);
        assert!((center - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_bilinear_exact() {
        let grid = GridF32::from_vec(2, 2, vec![0.0, 1.0, 2.0, 3.0]);
        // Exact corner values
        assert!((grid.sample_bilinear(0.0, 0.0) - 0.0).abs() < 1e-6);
        assert!((grid.sample_bilinear(1.0, 0.0) - 1.0).abs() < 1e-6);
        assert!((grid.sample_bilinear(0.0, 1.0) - 2.0).abs() < 1e-6);
        assert!((grid.sample_bilinear(1.0, 1.0) - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_gradient_flat() {
        let grid = GridF32::new(5, 5, 1.0);
        let (gx, gy) = grid.gradient_at(2, 2);
        assert!((gx).abs() < 1e-6);
        assert!((gy).abs() < 1e-6);
    }

    #[test]
    fn test_gradient_slope_x() {
        // Linear ramp in x: 0, 1, 2, 3, 4
        let data: Vec<f32> = (0..25).map(|i| (i % 5) as f32).collect();
        let grid = GridF32::from_vec(5, 5, data);
        let (gx, gy) = grid.gradient_at(2, 2);
        assert!((gx - 1.0).abs() < 1e-6); // slope of 1 per pixel in x
        assert!((gy).abs() < 1e-6); // flat in y
    }

    #[test]
    fn test_statistics() {
        let grid = GridF32::from_vec(3, 1, vec![1.0, 2.0, 3.0]);
        assert!((grid.min() - 1.0).abs() < 1e-6);
        assert!((grid.max() - 3.0).abs() < 1e-6);
        assert!((grid.mean() - 2.0).abs() < 1e-6);
    }

    #[test]
    #[should_panic(expected = "data length")]
    fn test_from_vec_wrong_size() {
        GridF32::from_vec(3, 3, vec![0.0; 10]);
    }

    #[test]
    fn periodic_wrapping() {
        let mut g = GridF32::new(4, 4, 0.0);
        g.set(0, 0, 1.0);
        assert_eq!(g.get_periodic(-4, 0), 1.0);
        assert_eq!(g.get_periodic(4, 0), 1.0);
        assert_eq!(g.get_periodic(0, -4), 1.0);
    }

    #[test]
    fn periodic_bilinear_wraps() {
        let mut g = GridF32::new(8, 8, 0.0);
        g.set(7, 7, 1.0);
        g.set(0, 0, 1.0);
        let v = g.sample_bilinear_periodic(7.5, 7.5);
        assert!((v - 0.5).abs() < 0.01, "Periodic bilinear across boundary: {v}");
    }

    #[test]
    fn periodic_gradient_wraps() {
        let mut g = GridF32::new(8, 8, 0.0);
        g.set(0, 0, 1.0);
        g.set(1, 0, 0.0);
        g.set(7, 0, 0.5);
        let (gx, _) = g.gradient_at_periodic(0, 0);
        // gx = (val[1] - val[7]) / 2 = (0.0 - 0.5) / 2 = -0.25
        assert!((gx - (-0.25)).abs() < 1e-6, "Periodic gradient: {gx}");
    }
}
