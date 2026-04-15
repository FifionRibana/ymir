//! Centering: shift the continent to the middle of the toroidal grid.

use super::plates::Plate;
use super::solver::field::Field2D;

/// Compute the circular shift (dx, dy) needed to center the continental mass.
///
/// `threshold` separates continental cells (above) from oceanic (below).
/// Returns the shift in integer grid cells.
pub fn compute_centering_shift(field: &Field2D, threshold: f64) -> (i32, i32) {
    let n = field.n();
    let tau = std::f64::consts::TAU;

    let mut sum_sin_x = 0.0;
    let mut sum_cos_x = 0.0;
    let mut sum_sin_y = 0.0;
    let mut sum_cos_y = 0.0;
    let mut count = 0usize;

    for j in 0..n {
        for i in 0..n {
            if field.get(i, j) > threshold {
                let theta_x = tau * i as f64 / n as f64;
                let theta_y = tau * j as f64 / n as f64;
                sum_sin_x += theta_x.sin();
                sum_cos_x += theta_x.cos();
                sum_sin_y += theta_y.sin();
                sum_cos_y += theta_y.cos();
                count += 1;
            }
        }
    }

    if count == 0 {
        return (0, 0);
    }

    let mean_angle_x = sum_sin_x.atan2(sum_cos_x);
    let mean_angle_y = sum_sin_y.atan2(sum_cos_y);

    let center_x = mean_angle_x * n as f64 / tau;
    let center_y = mean_angle_y * n as f64 / tau;

    let dx = (n as f64 / 2.0 - center_x).round() as i32;
    let dy = (n as f64 / 2.0 - center_y).round() as i32;

    (dx, dy)
}

/// Circularly shift a Field2D by (dx, dy) cells.
pub fn shift_field(field: &Field2D, dx: i32, dy: i32) -> Field2D {
    let n = field.n();
    let mut result = Field2D::new(n);
    let ni = n as i32;

    for j in 0..n {
        for i in 0..n {
            let si = ((i as i32 - dx) % ni + ni) as usize % n;
            let sj = ((j as i32 - dy) % ni + ni) as usize % n;
            result.set(i, j, field.get(si, sj));
        }
    }
    result
}

/// Circularly shift a flat grid of plate IDs.
pub fn shift_ids(ids: &[usize], n: usize, dx: i32, dy: i32) -> Vec<usize> {
    let ni = n as i32;
    let mut result = vec![0usize; n * n];

    for j in 0..n {
        for i in 0..n {
            let si = ((i as i32 - dx) % ni + ni) as usize % n;
            let sj = ((j as i32 - dy) % ni + ni) as usize % n;
            result[j * n + i] = ids[sj * n + si];
        }
    }
    result
}

/// Shift plate seed positions by (dx, dy) with wrapping.
pub fn shift_plates(plates: &mut [Plate], n: usize, dx: i32, dy: i32) {
    let nf = n as f32;
    for plate in plates.iter_mut() {
        plate.seed_x = ((plate.seed_x + dx as f32) % nf + nf) % nf;
        plate.seed_y = ((plate.seed_y + dy as f32) % nf + nf) % nf;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrapping_continent_centers_correctly() {
        let n = 16;
        let mut field = Field2D::new(n);
        // Continental cells at columns 14-15 and 0-1 (straddling right edge)
        for j in 6..10 {
            for i in [14, 15, 0, 1] {
                field.set(i, j, 1.0);
            }
        }

        let (dx, _dy) = compute_centering_shift(&field, 0.4);
        let shifted = shift_field(&field, dx, 0);

        // After shift, the center of mass in x should be near n/2 = 8
        let mut sum = 0.0;
        let mut count = 0;
        for j in 0..n {
            for i in 0..n {
                if shifted.get(i, j) > 0.4 {
                    sum += i as f64;
                    count += 1;
                }
            }
        }
        let mean_x = sum / count as f64;
        assert!((mean_x - 8.0).abs() < 2.0, "Center of mass should be near 8, got {mean_x}");
    }

    #[test]
    fn already_centered_no_shift() {
        let n = 16;
        let mut field = Field2D::new(n);
        // Continental cells centered around (8, 8)
        for j in 6..10 {
            for i in 6..10 {
                field.set(i, j, 1.0);
            }
        }

        let (dx, dy) = compute_centering_shift(&field, 0.4);
        assert!(dx.abs() <= 1 && dy.abs() <= 1, "Should be near zero shift: ({dx}, {dy})");
    }

    #[test]
    fn shift_is_reversible() {
        let n = 8;
        let mut field = Field2D::new(n);
        for j in 0..n {
            for i in 0..n {
                field.set(i, j, (i * n + j) as f64);
            }
        }

        let shifted = shift_field(&field, 3, -2);
        let restored = shift_field(&shifted, -3, 2);

        for j in 0..n {
            for i in 0..n {
                assert!(
                    (field.get(i, j) - restored.get(i, j)).abs() < 1e-10,
                    "Mismatch at ({i}, {j})"
                );
            }
        }
    }
}
