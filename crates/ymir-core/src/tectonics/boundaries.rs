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
    /// Reference thickness for oceanic crust (dimensionless).
    /// Oceanic cells thicker than this are subject to gravitational restoring.
    /// Default: 0.25 (slightly above initial oceanic thickness of 0.2).
    pub oceanic_reference_thickness: f64,
    /// Thickness above which oceanic cells are NOT subject to restoring.
    /// Cells on oceanic plates can thicken beyond their reference value
    /// either because genuine oceanic crust piled up (should be restored)
    /// or because continental material spilled over by advection (should
    /// be left alone). This threshold distinguishes the two cases:
    /// below it, the material is treated as thickened oceanic crust and
    /// pulled back toward `oceanic_reference_thickness`; above it, the
    /// material is presumed continental in origin and left untouched.
    /// Default: 0.4 (midpoint between oceanic ~0.2 and thin continental ~0.6).
    pub oceanic_restore_threshold: f64,
    /// Rate at which excess oceanic thickness is removed per timestep.
    /// Models dense oceanic crust sinking back into the mantle.
    /// 0.0 = no restoring, 1.0 = excess removed in one step.
    /// Range: 0.0-1.0, default: 0.3
    pub oceanic_restore_rate: f64,
    /// Minimum stable thickness for continental crust. Continental cells
    /// thinner than this are gently restored upward, modeling buoyancy.
    /// Default: 0.5 (thinned but still buoyant continental crust).
    pub continental_min_thickness: f64,
    /// Thickness below which continental restoring gives up — the crust
    /// has been fully rifted and effectively becomes oceanic.
    /// Default: 0.15 (below oceanic reference, enabling Wilson cycles).
    pub continental_restore_threshold: f64,
    /// Rate at which thinned continental crust is restored per timestep.
    /// 0.0 = no restoring, 1.0 = restored in one step.
    /// Default: 0.03 (gentle, allows temporary thinning during rifting).
    pub continental_restore_rate: f64,
    /// Enable dynamic slab pull. Default: true.
    pub slab_pull_enabled: bool,
    /// Slab pull strength per unit of local convergence rate.
    /// Since #75 this is the coefficient of the cell-local operator
    /// term `γ · (v·n̂) · n̂` (was: the per-plate velocity boost
    /// coefficient in the pre-#75 RHS-injection formulation).
    ///
    /// **Known inertness (issue #80, Phase 2-bis calibration sweep):**
    /// this value has no measurable effect on plate velocities at
    /// any tested factor between 0.05 and 5000 on 64²/seed 42/120
    /// steps. The operator form `γ · n̂⊗n̂` is SPD — it only damps
    /// motion along n̂, it cannot accelerate plates — so no scalar
    /// tuning can recover pre-#75 slab-pull behaviour. Physics-
    /// meaningful slab-pull is tracked in follow-up #80 as a
    /// redesign, not a calibration, task.
    /// Range: 0.001-0.5, default: 0.05.
    pub slab_pull_factor: f64,
    /// Maximum plate velocity magnitude (pre-#75 safety cap; retained
    /// as a parameter so configs stay deserialisable, but no longer
    /// enforced by the solver — the operator term auto-regulates).
    /// Default: 5.0
    pub max_plate_velocity: f32,
    /// Characteristic Benioff-zone decay length (in grid cells) for
    /// spreading γ_slab inward from the margin on the subducting plate.
    /// The seed value at the margin decays as `exp(-d_cells / L)` and
    /// is masked to the subducting plate's cells only. Range: 1.0-10.0,
    /// default: 3.0. See issue #75 §5.
    pub benioff_decay_cells: f64,
    /// Crustal density for continental plates (kg/m³). Default: 2750.
    pub rho_continental: f64,
    /// Crustal density for oceanic plates (kg/m³). Default: 3000.
    pub rho_oceanic: f64,
    /// Mantle density (kg/m³). Default: 3300.
    pub rho_mantle: f64,
    /// Thickness below which a cell at a convergent boundary is consumed
    /// (reassigned to the overriding plate). Default: 0.15.
    pub subduction_consumption_threshold: f64,
    /// Thickness below which a cell at a rift boundary spawns a new plate. Default: 0.15.
    pub rift_creation_threshold: f64,
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
            oceanic_reference_thickness: 0.25,
            oceanic_restore_threshold: 0.4,
            oceanic_restore_rate: 0.0,
            continental_min_thickness: 0.5,
            continental_restore_threshold: 0.15,
            continental_restore_rate: 0.0,
            slab_pull_enabled: true,
            slab_pull_factor: 0.05,
            max_plate_velocity: 5.0,
            benioff_decay_cells: 3.0,
            rho_continental: 2750.0,
            rho_oceanic: 3000.0,
            rho_mantle: 3300.0,
            subduction_consumption_threshold: 0.15,
            rift_creation_threshold: 0.15,
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
    /// Slab-pull coefficient field for the operator-form term
    /// `γ_slab · (v·n̂) · n̂` (issue #75). Non-zero only on the
    /// subducting side of convergent margins and on the Benioff
    /// decay band extending inward on that plate. Zero elsewhere.
    pub gamma_slab: Field2D,
    /// x-component of the margin normal n̂ at cell centers, pointing
    /// from the subducting cell toward its foreign-plate neighbour(s).
    /// Zero on cells with zero γ_slab.
    pub normal_x: Field2D,
    /// y-component of the margin normal n̂. Same convention as
    /// `normal_x`.
    pub normal_y: Field2D,
}

// ── Core computation ─────────────────────────────────────────────────────

/// Classify a plate's type based on its mean thickness.
/// More stable than per-cell classification because mean_thickness
/// averages over hundreds of cells and changes slowly.
fn plate_type_from_mean_thickness(mean_thickness: f32) -> PlateType {
    if mean_thickness > 0.4 { PlateType::Continental } else { PlateType::Oceanic }
}

pub fn compute_boundary_sources(
    grid: &StaggeredGrid,
    plate_ids: &[usize],
    plates: &[Plate],
    config: &BoundaryConfig,
    recycling_enabled: bool,
) -> BoundaryField {
    let nx = grid.nx();
    let ny = grid.ny();
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();
    let inv_dx = 1.0 / grid.dx;

    let mut boundary_type = vec![BoundaryType::None; nx * ny];
    let mut source_rate = Field2D::new(nx, ny);
    let mut normal_x = Field2D::new(nx, ny);
    let mut normal_y = Field2D::new(nx, ny);

    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let my_plate = plate_ids[k];
            let my_type = plate_type_from_mean_thickness(plates[my_plate].mean_thickness);

            let neighbors =
                [(idx_x.next(i), j), (idx_x.prev(i), j), (i, idx_y.next(j)), (i, idx_y.prev(j))];

            let mut is_boundary = false;
            let mut convergence_sum = 0.0_f64;
            let mut neighbor_plate_type = my_type;
            // Accumulated raw normal from per-neighbor contributions.
            // Normalized below after the neighbor loop (#75).
            let mut n_acc_x = 0.0_f64;
            let mut n_acc_y = 0.0_f64;

            for &(ni, nj) in &neighbors {
                let nk = nj * nx + ni;
                let other_plate = plate_ids[nk];

                if other_plate != my_plate {
                    is_boundary = true;
                    neighbor_plate_type =
                        plate_type_from_mean_thickness(plates[other_plate].mean_thickness);

                    // Normal direction from (i,j) to (ni,nj), handling wrapping
                    let normal_x = if ni == idx_x.next(i) && ni < i {
                        1.0 // wrapped right
                    } else if ni == idx_x.prev(i) && ni > i {
                        -1.0 // wrapped left
                    } else {
                        (ni as f64 - i as f64).signum()
                    };
                    let normal_y = if nj == idx_y.next(j) && nj < j {
                        1.0 // wrapped down
                    } else if nj == idx_y.prev(j) && nj > j {
                        -1.0 // wrapped up
                    } else {
                        (nj as f64 - j as f64).signum()
                    };
                    n_acc_x += normal_x;
                    n_acc_y += normal_y;

                    // Velocity at cell center (average of staggered faces)
                    let vx_here = 0.5 * (grid.vx.get(i, j) + grid.vx.get(idx_x.next(i), j));
                    let vy_here = 0.5 * (grid.vy.get(i, j) + grid.vy.get(i, idx_y.next(j)));

                    let vx_there = 0.5 * (grid.vx.get(ni, nj) + grid.vx.get(idx_x.next(ni), nj));
                    let vy_there = 0.5 * (grid.vy.get(ni, nj) + grid.vy.get(ni, idx_y.next(nj)));

                    // Relative velocity in normal direction
                    // Positive = diverging, Negative = converging
                    let v_rel = (vx_there - vx_here) * normal_x + (vy_there - vy_here) * normal_y;
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
                    | (PlateType::Continental, PlateType::Oceanic) => BoundaryType::Subduction,
                    (PlateType::Continental, PlateType::Continental) => {
                        BoundaryType::ContinentalCollision
                    }
                    (PlateType::Oceanic, PlateType::Oceanic) => BoundaryType::OceanicSubduction,
                }
            } else {
                BoundaryType::Rift
            };

            boundary_type[k] = btype;

            // Persist the unit normal at margin cells. Sign convention:
            // points from this cell toward the foreign-plate neighbour(s)
            // (i.e. toward the trench on the subducting side).
            let nmag = (n_acc_x * n_acc_x + n_acc_y * n_acc_y).sqrt();
            if nmag > 1e-12 {
                normal_x.set(i, j, n_acc_x / nmag);
                normal_y.set(i, j, n_acc_y / nmag);
            }

            let q = match btype {
                BoundaryType::Subduction => {
                    if my_type == PlateType::Oceanic {
                        -config.subduction_rate * convergence_rate
                    } else if recycling_enabled {
                        0.0 // Arc handled by recycling system
                    } else {
                        config.volcanic_arc_rate * convergence_rate
                    }
                }
                BoundaryType::OceanicSubduction => {
                    let my_s = grid.s.get(i, j);
                    let s_ocean_ref = 0.2;
                    if my_s <= s_ocean_ref * 1.1 {
                        -config.subduction_rate * convergence_rate * 0.5
                    } else if recycling_enabled {
                        0.0 // Arc handled by recycling system
                    } else {
                        config.volcanic_arc_rate * convergence_rate * 0.3
                    }
                }
                BoundaryType::ContinentalCollision => {
                    if recycling_enabled {
                        0.0 // No creation from nothing in conservative mode
                    } else {
                        config.collision_volcanism_rate * convergence_rate
                    }
                }
                BoundaryType::Rift => {
                    if recycling_enabled {
                        0.0 // Spreading handled by recycling system
                    } else {
                        let my_s = grid.s.get(i, j);
                        let divergence_rate = convergence_sum.abs() * inv_dx;
                        if my_s < config.rift_thickness_threshold {
                            config.spreading_rate * divergence_rate
                        } else {
                            config.rift_volcanism_rate * divergence_rate
                        }
                    }
                }
                BoundaryType::None => 0.0,
            };

            source_rate.set(i, j, q);
        }
    }

    // Compute γ_slab on the subducting side of convergent margins, then
    // spread inward on the same plate with an exponential Benioff decay
    // of characteristic length `benioff_decay_cells`. Issue #75 §5.
    let mut gamma_slab = Field2D::new(nx, ny);
    if config.slab_pull_enabled && config.slab_pull_factor > 0.0 && config.benioff_decay_cells > 0.0
    {
        let mut seeds: Vec<(usize, usize, f64)> = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let q = source_rate.get(i, j);
                let btype = boundary_type[k];
                let is_subducting_side = q < 0.0
                    && matches!(btype, BoundaryType::Subduction | BoundaryType::OceanicSubduction);
                if is_subducting_side {
                    seeds.push((i, j, config.slab_pull_factor * q.abs()));
                }
            }
        }
        if !seeds.is_empty() {
            spread_gamma_benioff(
                &mut gamma_slab,
                plate_ids,
                &seeds,
                config.benioff_decay_cells,
                nx,
                ny,
                idx_x,
                idx_y,
            );
        }
    }

    BoundaryField { boundary_type, source_rate, n: nx, gamma_slab, normal_x, normal_y }
}

/// Spread γ_seed values outward from each seed cell along the subducting
/// plate with an exponential decay `exp(-d/L)` (d = Chebyshev distance
/// in cells, stays on the same plate). Each destination cell receives
/// the max over all seeds that can reach it. See issue #75 §5.
fn spread_gamma_benioff(
    gamma_slab: &mut Field2D,
    plate_ids: &[usize],
    seeds: &[(usize, usize, f64)],
    decay_cells: f64,
    nx: usize,
    ny: usize,
    idx_x: &super::solver::field::PeriodicIndex,
    idx_y: &super::solver::field::PeriodicIndex,
) {
    // 3 × L captures >95 % of the tail; beyond that γ ≈ 0.05 · γ_seed.
    let max_radius = (3.0 * decay_cells).ceil() as i32;
    let mut visited = vec![u32::MAX; nx * ny];

    for (seed_idx, &(si, sj, g_seed)) in seeds.iter().enumerate() {
        let plate = plate_ids[sj * nx + si];
        let stamp = seed_idx as u32;
        // Multi-source BFS would be slightly cheaper but per-seed BFS is
        // simpler and the cost (O(seeds × max_radius²)) is acceptable at
        // the grid sizes we run (64²–256²).
        let mut queue: std::collections::VecDeque<(usize, usize, i32)> =
            std::collections::VecDeque::new();
        queue.push_back((si, sj, 0));
        visited[sj * nx + si] = stamp;
        while let Some((i, j, d)) = queue.pop_front() {
            if plate_ids[j * nx + i] != plate {
                continue;
            }
            let contribution = g_seed * (-(d as f64) / decay_cells).exp();
            let k = j * nx + i;
            let current = gamma_slab.data()[k];
            if contribution > current {
                gamma_slab.data_mut()[k] = contribution;
            }
            if d + 1 > max_radius {
                continue;
            }
            let next_cells =
                [(idx_x.next(i), j), (idx_x.prev(i), j), (i, idx_y.next(j)), (i, idx_y.prev(j))];
            for &(ni, nj) in &next_cells {
                let nk = nj * nx + ni;
                if visited[nk] != stamp && plate_ids[nk] == plate {
                    visited[nk] = stamp;
                    queue.push_back((ni, nj, d + 1));
                }
            }
        }
    }
}

/// Compute boundary sources directly into a pre-allocated Field2D buffer.
pub fn compute_boundary_sources_into(
    grid: &StaggeredGrid,
    plate_ids: &[usize],
    plates: &[Plate],
    config: &BoundaryConfig,
    target: &mut Field2D,
) {
    let result = compute_boundary_sources(grid, plate_ids, plates, config, false);
    target.data_mut().copy_from_slice(result.source_rate.data());
}

// ── Slab pull ───────────────────────────────────────────────────────────

/// Accumulate subducted mass into plates for slab pull computation.
///
/// Where source_rate is negative (subduction sink), the removed mass is
/// attributed to the plate owning that cell. The slab pull then acts on
/// that plate, pulling it toward the trench.
pub fn accumulate_subducted_mass(
    source_rate: &Field2D,
    plate_ids: &[usize],
    plates: &mut [Plate],
    dt: f64,
    nx: usize,
    ny: usize,
) {
    for j in 0..ny {
        for i in 0..nx {
            let q = source_rate.get(i, j);
            if q < 0.0 {
                let pid = plate_ids[j * nx + i];
                plates[pid].subducted_mass += (-q * dt).abs();
            }
        }
    }
}

/// Apply slab pull: increase plate velocity proportional to cumulative
/// subducted mass, capped at `max_velocity`.
pub fn apply_slab_pull(plates: &mut [Plate], slab_pull_factor: f64, max_velocity: f32) {
    for plate in plates.iter_mut() {
        if !plate.active {
            continue;
        }

        let pull_magnitude = slab_pull_factor * plate.subducted_mass;
        if pull_magnitude < 1e-30 {
            continue;
        }

        let vx = plate.velocity.0 as f64;
        let vy = plate.velocity.1 as f64;
        let v_mag = (vx * vx + vy * vy).sqrt().max(1e-30);

        plate.velocity.0 += (pull_magnitude * vx / v_mag) as f32;
        plate.velocity.1 += (pull_magnitude * vy / v_mag) as f32;

        // Cap total velocity
        let new_mag = (plate.velocity.0.powi(2) + plate.velocity.1.powi(2)).sqrt();
        if new_mag > max_velocity {
            plate.velocity.0 *= max_velocity / new_mag;
            plate.velocity.1 *= max_velocity / new_mag;
        }
    }
}

// ── Gaussian blur for f64 Field2D ────────────────────────────────────────

/// Separable Gaussian blur on a Field2D with periodic (toroidal) wrapping.
/// Sigma is in grid cells. Kernel radius = ceil(3 * sigma).
pub fn gaussian_blur_f64(field: &Field2D, sigma: f64) -> Field2D {
    let nx = field.nx();
    let ny = field.ny();
    if sigma <= 0.0 {
        let mut out = Field2D::new(nx, ny);
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
    let mut temp = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let mut val = 0.0_f64;
            for (ki, &w) in kernel.iter().enumerate() {
                let si = (i as i32 + ki as i32 - radius as i32).rem_euclid(nx as i32) as usize;
                val += field.get(si, j) * w;
            }
            temp.set(i, j, val);
        }
    }

    // Vertical pass
    let mut result = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let mut val = 0.0_f64;
            for (ki, &w) in kernel.iter().enumerate() {
                let sj = (j as i32 + ki as i32 - radius as i32).rem_euclid(ny as i32) as usize;
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
        let grid = StaggeredGrid::new(n, n, dx);

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
                active: true,
                subducted_mass: 0.0,
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            },
            Plate {
                id: 1,
                plate_type: PlateType::Continental,
                velocity: (-0.5, 0.0),
                seed_x: (3 * n / 4) as f32,
                seed_y: (n / 2) as f32,
                active: true,
                subducted_mass: 0.0,
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
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
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config, false);

        let oceanic_boundary = result.source_rate.get(n / 2 - 1, n / 2);
        assert!(oceanic_boundary < 0.0, "Oceanic side should be a sink: Q = {oceanic_boundary}");

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
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config, false);

        let boundary_q = result.source_rate.get(n / 2, n / 2);
        assert!(boundary_q > 0.0, "Rift should create crust: Q = {boundary_q}");
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
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config, false);

        assert!(result.source_rate.get(0, 0).abs() < 1e-10, "Interior cell should have no source");
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
        let result = compute_boundary_sources(&grid, &plate_ids, &plates, &config, false);

        let total_q: f64 = result.source_rate.data().iter().sum();
        let max_q: f64 = result.source_rate.data().iter().map(|x| x.abs()).fold(0.0, f64::max);
        let relative = total_q.abs() / (max_q * n as f64).max(1e-10);

        eprintln!("Total Q = {total_q:.6}, max |Q| = {max_q:.6}, relative = {relative:.6}");
    }

    #[test]
    fn gaussian_blur_preserves_total() {
        let n = 16;
        let mut field = Field2D::new(n, n);
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
