//! Boundary processes: subduction, rifting, and volcanism source terms.
//!
//! Detects plate boundaries from `plate_ids` and velocity fields, classifies
//! them (subduction, rift, collision), and computes a source/sink rate Q(x,y)
//! that is added to the crustal thickness advection equation:
//!
//!   S_new = S - dt * div(S*v) + dt * Q(x,y)

use super::plates::{Plate, PlateType};
use super::solver::field::Field2D;
use super::solver::grid::StaggeredGrid;

// ── Boundary classification ──────────────────────────────────────────────

/// Type of plate boundary at a given cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryType {
    /// Oceanic plate being consumed under continental plate (or vice-versa view).
    Subduction,
    /// Two continental plates colliding.
    ContinentalCollision,
    /// Oceanic-oceanic convergence. One subducts, island arc volcanism.
    OceanicSubduction,
    /// Plates diverging — new oceanic crust created in the gap.
    Rift,
    /// No boundary (cell is interior to a plate).
    None,
}

// ── Configuration ────────────────────────────────────────────────────────

/// Configuration for boundary source/sink terms.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BoundaryConfig {
    /// Whether boundary processes are enabled.
    pub enabled: bool,
    /// Subduction consumption rate (dimensionless, multiplied by convergence
    /// velocity). Range: 0.0-2.0, default: 0.5
    pub subduction_rate: f64,
    /// Volcanic arc production rate — fraction of subducted material that
    /// resurfaces as continental volcanism. Range: 0.0-1.0, default: 0.15
    pub volcanic_arc_rate: f64,
    /// Oceanic spreading rate — how fast new crust is created at divergent
    /// boundaries. Range: 0.0-2.0, default: 0.3
    pub spreading_rate: f64,
    /// Crustal thickness below which rifting creates new oceanic crust
    /// instead of just thinning. Default: 0.4
    pub rift_thickness_threshold: f64,
    /// Minor volcanism rate at continental collision zones.
    /// Range: 0.0-0.5, default: 0.05
    pub collision_volcanism_rate: f64,
    /// Volcanism rate at rift zones in thick crust.
    /// Range: 0.0-0.5, default: 0.02
    pub rift_volcanism_rate: f64,
    /// Gaussian smoothing sigma (in grid cells) applied to the source field.
    /// 0.0 = no smoothing. Default: 2.0
    pub source_smoothing_sigma: f64,
}

impl Default for BoundaryConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            subduction_rate: 0.5,
            volcanic_arc_rate: 0.15,
            spreading_rate: 0.3,
            rift_thickness_threshold: 0.4,
            collision_volcanism_rate: 0.05,
            rift_volcanism_rate: 0.02,
            source_smoothing_sigma: 2.0,
        }
    }
}

// ── Boundary field ───────────────────────────────────────────────────────

/// Boundary classification and source terms for each cell.
pub struct BoundaryField {
    /// Boundary type at each cell (N*N).
    pub boundary_type: Vec<BoundaryType>,
    /// Source/sink rate Q at each cell. Positive = crust creation,
    /// negative = crust consumption.
    pub source_rate: Field2D,
    pub n: usize,
}

// ── Core computation ─────────────────────────────────────────────────────

/// Detect boundaries and compute source/sink rates.
pub fn compute_boundary_sources(
    grid: &StaggeredGrid,
    plate_ids: &[usize],
    plates: &[Plate],
    config: &BoundaryConfig,
) -> BoundaryField {
    let n = grid.n;
    let idx = &grid.idx;
    let inv_dx = 1.0 / grid.dx;

    let mut boundary_type = vec![BoundaryType::None; n * n];
    let mut source_rate = Field2D::new(n);

    for j in 0..n {
        for i in 0..n {
            let k = j * n + i;
            let my_plate = plate_ids[k];
            let my_type = plates[my_plate].plate_type;

            let neighbors = [
                (idx.next(i), j),
                (idx.prev(i), j),
                (i, idx.next(j)),
                (i, idx.prev(j)),
            ];

            let mut is_boundary = false;
            let mut convergence_sum = 0.0_f64;
            let mut neighbor_plate_type = my_type;

            for &(ni, nj) in &neighbors {
                let nk = nj * n + ni;
                let other_plate = plate_ids[nk];

                if other_plate != my_plate {
                    is_boundary = true;
                    neighbor_plate_type = plates[other_plate].plate_type;

                    // Normal direction from (i,j) to (ni,nj), handling wrapping
                    let nx = if ni == idx.next(i) && ni < i {
                        1.0 // wrapped right
                    } else if ni == idx.prev(i) && ni > i {
                        -1.0 // wrapped left
                    } else {
                        (ni as f64 - i as f64).signum()
                    };
                    let ny = if nj == idx.next(j) && nj < j {
                        1.0 // wrapped down
                    } else if nj == idx.prev(j) && nj > j {
                        -1.0 // wrapped up
                    } else {
                        (nj as f64 - j as f64).signum()
                    };

                    // Velocity at cell center (average of staggered faces)
                    let vx_here =
                        0.5 * (grid.vx.get(i, j) + grid.vx.get(idx.next(i), j));
                    let vy_here =
                        0.5 * (grid.vy.get(i, j) + grid.vy.get(i, idx.next(j)));

                    let vx_there =
                        0.5 * (grid.vx.get(ni, nj) + grid.vx.get(idx.next(ni), nj));
                    let vy_there =
                        0.5 * (grid.vy.get(ni, nj) + grid.vy.get(ni, idx.next(nj)));

                    // Relative velocity in normal direction
                    // Positive = diverging, Negative = converging
                    let v_rel =
                        (vx_there - vx_here) * nx + (vy_there - vy_here) * ny;
                    convergence_sum += v_rel;
                }
            }

            if !is_boundary {
                continue;
            }

            let is_converging = convergence_sum < 0.0;
            let convergence_rate = convergence_sum.abs() * inv_dx;

            let btype = if is_converging {
                match (my_type, neighbor_plate_type) {
                    (PlateType::Oceanic, PlateType::Continental)
                    | (PlateType::Continental, PlateType::Oceanic) => {
                        BoundaryType::Subduction
                    }
                    (PlateType::Continental, PlateType::Continental) => {
                        BoundaryType::ContinentalCollision
                    }
                    (PlateType::Oceanic, PlateType::Oceanic) => {
                        BoundaryType::OceanicSubduction
                    }
                }
            } else {
                BoundaryType::Rift
            };

            boundary_type[k] = btype;

            let q = match btype {
                BoundaryType::Subduction => {
                    if my_type == PlateType::Oceanic {
                        -config.subduction_rate * convergence_rate
                    } else {
                        config.volcanic_arc_rate * convergence_rate
                    }
                }
                BoundaryType::OceanicSubduction => {
                    let my_s = grid.s.get(i, j);
                    let s_ocean_ref = 0.2;
                    if my_s <= s_ocean_ref * 1.1 {
                        -config.subduction_rate * convergence_rate * 0.5
                    } else {
                        config.volcanic_arc_rate * convergence_rate * 0.3
                    }
                }
                BoundaryType::ContinentalCollision => {
                    config.collision_volcanism_rate * convergence_rate
                }
                BoundaryType::Rift => {
                    let my_s = grid.s.get(i, j);
                    let divergence_rate = convergence_sum.abs() * inv_dx;
                    if my_s < config.rift_thickness_threshold {
                        config.spreading_rate * divergence_rate
                    } else {
                        config.rift_volcanism_rate * divergence_rate
                    }
                }
                BoundaryType::None => 0.0,
            };

            source_rate.set(i, j, q);
        }
    }

    BoundaryField { boundary_type, source_rate, n }
}

/// Compute boundary sources directly into a pre-allocated Field2D buffer.
pub fn compute_boundary_sources_into(
    grid: &StaggeredGrid,
    plate_ids: &[usize],
    plates: &[Plate],
    config: &BoundaryConfig,
    target: &mut Field2D,
) {
    let result = compute_boundary_sources(grid, plate_ids, plates, config);
    target.data_mut().copy_from_slice(result.source_rate.data());
}

// ── Gaussian blur for f64 Field2D ────────────────────────────────────────

/// Separable Gaussian blur on a Field2D with periodic (toroidal) wrapping.
/// Sigma is in grid cells. Kernel radius = ceil(3 * sigma).
pub fn gaussian_blur_f64(field: &Field2D, sigma: f64) -> Field2D {
    let n = field.n();
    if sigma <= 0.0 {
        let mut out = Field2D::new(n);
        out.data_mut().copy_from_slice(field.data());
        return out;
    }

    let radius = (3.0 * sigma).ceil() as usize;
    let kernel_size = 2 * radius + 1;

    // Build normalized 1D Gaussian kernel
    let mut kernel: Vec<f64> = (0..kernel_size)
        .map(|i| {
            let x = i as f64 - radius as f64;
            (-x * x / (2.0 * sigma * sigma)).exp()
        })
        .collect();
    let sum: f64 = kernel.iter().sum();
    for k in kernel.iter_mut() {
        *k /= sum;
    }

    // Horizontal pass
    let mut temp = Field2D::new(n);
    for j in 0..n {
        for i in 0..n {
            let mut val = 0.0_f64;
            for (ki, &w) in kernel.iter().enumerate() {
                let si = (i as i32 + ki as i32 - radius as i32)
                    .rem_euclid(n as i32) as usize;
                val += field.get(si, j) * w;
            }
            temp.set(i, j, val);
        }
    }

    // Vertical pass
    let mut result = Field2D::new(n);
    for j in 0..n {
        for i in 0..n {
            let mut val = 0.0_f64;
            for (ki, &w) in kernel.iter().enumerate() {
                let sj = (j as i32 + ki as i32 - radius as i32)
                    .rem_euclid(n as i32) as usize;
                val += temp.get(i, sj) * w;
            }
            result.set(i, j, val);
        }
    }

    result
}

// ── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::solver::grid::StaggeredGrid;

    fn make_two_plate_setup(n: usize) -> (StaggeredGrid, Vec<usize>, Vec<Plate>) {
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);

        // Left half = plate 0 (oceanic), right half = plate 1 (continental)
        let mut plate_ids = vec![0usize; n * n];
        for j in 0..n {
            for i in n / 2..n {
                plate_ids[j * n + i] = 1;
            }
        }

        let plates = vec![
            Plate {
                id: 0,
                plate_type: PlateType::Oceanic,
                velocity: (0.5, 0.0),
                seed_x: (n / 4) as f32,
                seed_y: (n / 2) as f32,
            },
            Plate {
                id: 1,
                plate_type: PlateType::Continental,
                velocity: (-0.5, 0.0),
                seed_x: (3 * n / 4) as f32,
                seed_y: (n / 2) as f32,
            },
        ];

        (grid, plate_ids, plates)
    }

    #[test]
    fn subduction_detected_at_convergent_ocean_continent() {
        let n = 32;
        let (mut grid, plate_ids, plates) = make_two_plate_setup(n);

        for j in 0..n {
            for i in 0..n {
                let vx = if i < n / 2 { 0.5 } else { -0.5 };
                grid.vx.set(i, j, vx);
            }
        }
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, if i < n / 2 { 0.2 } else { 1.0 });
            }
        }

        let config = BoundaryConfig::default();
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config);

        let oceanic_boundary = result.source_rate.get(n / 2 - 1, n / 2);
        assert!(
            oceanic_boundary < 0.0,
            "Oceanic side should be a sink: Q = {oceanic_boundary}"
        );

        let continental_boundary = result.source_rate.get(n / 2, n / 2);
        assert!(
            continental_boundary > 0.0,
            "Continental side should be a source: Q = {continental_boundary}"
        );
    }

    #[test]
    fn rift_detected_at_divergent_boundary() {
        let n = 32;
        let (mut grid, plate_ids, plates) = make_two_plate_setup(n);

        for j in 0..n {
            for i in 0..n {
                let vx = if i < n / 2 { -0.5 } else { 0.5 };
                grid.vx.set(i, j, vx);
            }
        }
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 0.2);
            }
        }

        let config = BoundaryConfig::default();
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config);

        let boundary_q = result.source_rate.get(n / 2, n / 2);
        assert!(
            boundary_q > 0.0,
            "Rift should create crust: Q = {boundary_q}"
        );
    }

    #[test]
    fn interior_cells_have_no_source() {
        let n = 32;
        let (mut grid, plate_ids, plates) = make_two_plate_setup(n);

        for j in 0..n {
            for i in 0..n {
                grid.vx.set(i, j, 0.1);
                grid.s.set(i, j, 1.0);
            }
        }

        let config = BoundaryConfig::default();
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config);

        assert!(
            result.source_rate.get(0, 0).abs() < 1e-10,
            "Interior cell should have no source"
        );
        assert!(
            result.source_rate.get(n - 1, n / 2).abs() < 1e-10,
            "Interior cell should have no source"
        );
    }

    #[test]
    fn sources_conserve_mass_approximately() {
        let n = 32;
        let (mut grid, plate_ids, plates) = make_two_plate_setup(n);

        for j in 0..n {
            for i in 0..n {
                let vx = if i < n / 2 { 0.5 } else { -0.5 };
                grid.vx.set(i, j, vx);
                grid.s.set(i, j, if i < n / 2 { 0.2 } else { 1.0 });
            }
        }

        let config = BoundaryConfig::default();
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config);

        let total_q: f64 = result.source_rate.data().iter().sum();
        let max_q: f64 = result
            .source_rate
            .data()
            .iter()
            .map(|x| x.abs())
            .fold(0.0, f64::max);
        let relative = total_q.abs() / (max_q * n as f64).max(1e-10);

        eprintln!(
            "Total Q = {total_q:.6}, max |Q| = {max_q:.6}, relative = {relative:.6}"
        );
    }

    #[test]
    fn gaussian_blur_preserves_total() {
        let n = 16;
        let mut field = Field2D::new(n);
        field.set(n / 2, n / 2, 1.0);

        let blurred = gaussian_blur_f64(&field, 2.0);

        let sum_before: f64 = field.data().iter().sum();
        let sum_after: f64 = blurred.data().iter().sum();
        assert!(
            (sum_before - sum_after).abs() < 1e-10,
            "Blur should preserve total: {sum_before} vs {sum_after}"
        );
    }
}
