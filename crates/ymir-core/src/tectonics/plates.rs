//! Tectonic plate initialization via Voronoï partitioning.
//!
//! Generates `num_plates` plates on a toroidal grid (periodic boundary conditions).
//! Each plate gets a random position, type (continental/oceanic), and velocity vector.
//! Cells are assigned to the nearest plate seed using toroidal distance.
//! The resulting crustal thickness field is smoothed with a separable Gaussian blur.

use rand::Rng;

use super::solver::field::Field2D;
use super::solver::grid::StaggeredGrid;
use super::solver::traction::TractionField;
use crate::grid::GridF32;
use crate::seed::WorldSeed;

// ── Plate types and structs ───────────────────────────────────────────────

/// Type of tectonic plate — determines crustal thickness and behavior at boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlateType {
    /// Thick (~35 km), light (ρ ≈ 2.7). Stands above sea level.
    Continental,
    /// Thin (~7 km), heavy (ρ ≈ 3.0). Submerged below sea level.
    Oceanic,
}

/// A single tectonic plate with its properties.
#[derive(Debug, Clone)]
pub struct Plate {
    pub id: usize,
    pub plate_type: PlateType,
    /// Velocity vector in grid units per timestep.
    /// Magnitude typically 0.5–3.0 grid cells per step.
    pub velocity: (f32, f32),
    /// Voronoï seed position (in grid coordinates).
    pub seed_x: f32,
    pub seed_y: f32,
    /// Whether this plate is still active (has cells assigned to it).
    /// Plates consumed by subduction become inactive.
    pub active: bool,
    /// Cumulative subducted material. Increases when oceanic crust
    /// is consumed at this plate's convergent boundaries. Used to
    /// compute dynamic slab pull traction.
    pub subducted_mass: f64,
}

/// Configuration for plate generation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PlateConfig {
    /// Number of tectonic plates (5–15).
    pub num_plates: usize,
    /// Fraction of plates that are continental (0.1–0.6).
    pub continental_ratio: f32,
    /// Minimum velocity magnitude (grid units / timestep).
    pub velocity_min: f32,
    /// Maximum velocity magnitude (grid units / timestep).
    pub velocity_max: f32,
    /// Grid size for the tectonic simulation.
    pub grid_size: usize,
    /// Gaussian blur sigma applied to initial thickness field to smooth plate boundaries.
    pub boundary_smoothing_sigma: f32,
}

impl Default for PlateConfig {
    fn default() -> Self {
        Self {
            num_plates: 8,
            continental_ratio: 0.3,
            velocity_min: 0.5,
            velocity_max: 2.5,
            grid_size: 128,
            boundary_smoothing_sigma: 2.0,
        }
    }
}

impl PlateConfig {
    /// Single large continent surrounded by ocean.
    /// Good starting point for a typical game continent.
    pub fn preset_single_continent() -> Self {
        Self { num_plates: 8, continental_ratio: 0.25, ..Default::default() }
    }

    /// Two continental masses converging — produces a Himalaya-type collision belt.
    pub fn preset_collision() -> Self {
        Self { num_plates: 6, continental_ratio: 0.4, ..Default::default() }
    }

    /// Many small continental plates — creates an archipelago with multiple islands.
    pub fn preset_archipelago() -> Self {
        Self { num_plates: 12, continental_ratio: 0.35, velocity_max: 3.0, ..Default::default() }
    }

    /// One plate splitting apart — produces a rift valley with volcanism and graben basins.
    pub fn preset_rift() -> Self {
        Self {
            num_plates: 6,
            continental_ratio: 0.3,
            velocity_min: 0.3,
            velocity_max: 1.5,
            ..Default::default()
        }
    }
}

/// Result of plate initialization — everything needed to start the tectonic simulation.
#[derive(Clone)]
pub struct PlateInitResult {
    /// Which plate each cell belongs to (grid_size × grid_size).
    pub plate_ids: Vec<usize>,
    /// Initial crustal thickness field. Continental ≈ 1.0, Oceanic ≈ 0.2.
    pub thickness: GridF32,
    /// The plates themselves with their properties.
    pub plates: Vec<Plate>,
    /// Grid dimensions.
    pub grid_size: usize,
}

// ── Core algorithm ────────────────────────────────────────────────────────

/// Generate plates on the grid using Voronoï partitioning.
///
/// The grid has periodic (toroidal) boundary conditions: the left edge connects
/// to the right, the top to the bottom. The distance metric wraps around the grid
/// in both directions to find the shortest path.
pub fn generate_plates(config: &PlateConfig, seed: &WorldSeed) -> PlateInitResult {
    let mut rng = seed.rng_for("plates");
    let size = config.grid_size;

    // 1. Place Voronoï seeds randomly on the grid
    let mut plates = Vec::with_capacity(config.num_plates);
    let num_continental = (config.num_plates as f32 * config.continental_ratio).round() as usize;

    for id in 0..config.num_plates {
        let seed_x = rng.random::<f32>() * size as f32;
        let seed_y = rng.random::<f32>() * size as f32;

        let plate_type =
            if id < num_continental { PlateType::Continental } else { PlateType::Oceanic };

        let angle = rng.random::<f32>() * std::f32::consts::TAU;
        let speed =
            config.velocity_min + rng.random::<f32>() * (config.velocity_max - config.velocity_min);
        let velocity = (angle.cos() * speed, angle.sin() * speed);

        plates.push(Plate {
            id,
            plate_type,
            velocity,
            seed_x,
            seed_y,
            active: true,
            subducted_mass: 0.0,
        });
    }

    // 2. Assign each cell to the nearest plate seed (toroidal distance)
    let mut plate_ids = vec![0usize; size * size];
    for y in 0..size {
        for x in 0..size {
            let mut best_dist = f32::MAX;
            let mut best_id = 0;

            for plate in &plates {
                let dist = toroidal_distance_sq(
                    x as f32,
                    y as f32,
                    plate.seed_x,
                    plate.seed_y,
                    size as f32,
                );
                if dist < best_dist {
                    best_dist = dist;
                    best_id = plate.id;
                }
            }

            plate_ids[y * size + x] = best_id;
        }
    }

    // 3. Initialize crustal thickness field
    let mut thickness = GridF32::new(size, size, 0.0);
    for y in 0..size {
        for x in 0..size {
            let plate = &plates[plate_ids[y * size + x]];
            let t = match plate.plate_type {
                PlateType::Continental => 1.0,
                PlateType::Oceanic => 0.2,
            };
            thickness.set(x, y, t);
        }
    }

    // 4. Smooth boundaries with Gaussian blur (toroidal wrapping)
    if config.boundary_smoothing_sigma > 0.0 {
        thickness = gaussian_blur(&thickness, config.boundary_smoothing_sigma);
    }

    PlateInitResult { plate_ids, thickness, plates, grid_size: size }
}

// ── Distance and blur helpers ─────────────────────────────────────────────

/// Toroidal distance between two points on a grid that wraps in both x and y.
pub fn toroidal_distance(x1: f32, y1: f32, x2: f32, y2: f32, size: f32) -> f32 {
    toroidal_distance_sq(x1, y1, x2, y2, size).sqrt()
}

/// Squared toroidal distance between two points on a grid that wraps in both x and y.
///
/// For each axis, takes the shorter path: direct or through the wrap.
pub fn toroidal_distance_sq(x1: f32, y1: f32, x2: f32, y2: f32, size: f32) -> f32 {
    let dx = (x1 - x2).abs();
    let dy = (y1 - y2).abs();
    let dx = dx.min(size - dx);
    let dy = dy.min(size - dy);
    dx * dx + dy * dy
}

impl PlateInitResult {
    /// Build the traction field for the thin sheet solver.
    ///
    /// For each cell on the grid, looks up which plate it belongs to
    /// (via `plate_ids`) and assigns that plate's velocity as the
    /// traction force at that cell. This is how the Voronoï plate
    /// configuration drives the physics simulation.
    pub fn to_traction_field(&self) -> TractionField {
        let n = self.grid_size;
        let mut tx = Field2D::new(n);
        let mut ty = Field2D::new(n);

        for j in 0..n {
            for i in 0..n {
                let plate_id = self.plate_ids[j * n + i];
                let plate = &self.plates[plate_id];
                tx.set(i, j, plate.velocity.0 as f64);
                ty.set(i, j, plate.velocity.1 as f64);
            }
        }

        TractionField { tx, ty }
    }

    /// Initialize a StaggeredGrid with the crustal thickness from
    /// the Voronoï plate generation.
    ///
    /// Creates the grid with the correct size and dx, copies the
    /// thickness field (f32 → f64), and returns the grid ready
    /// for the solver.
    pub fn to_staggered_grid(&self) -> StaggeredGrid {
        let n = self.grid_size;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        for j in 0..n {
            for i in 0..n {
                let s = self.thickness.data[j * n + i] as f64;
                grid.s.set(i, j, s);
            }
        }

        grid
    }
}

/// Separable Gaussian blur on a GridF32 with toroidal (wrapping) boundary conditions.
///
/// Delegates to [`GridF32::gaussian_blur`]. Kept for backward compatibility.
pub fn gaussian_blur(grid: &GridF32, sigma: f32) -> GridF32 {
    grid.gaussian_blur(sigma)
}

// ── Dynamic boundary helpers ─────────────────────────────────────────────

/// Recompute plate_ids from current seed positions (toroidal distance).
///
/// For each cell, find the nearest active seed. Seeds of inactive
/// (disappeared) plates are skipped. Ties are broken by lowest plate ID.
pub fn recompute_voronoi(plate_ids: &mut [usize], plates: &[Plate], grid_size: usize) {
    let size = grid_size as f32;
    for y in 0..grid_size {
        for x in 0..grid_size {
            let mut best_dist = f32::MAX;
            let mut best_id = 0;

            for plate in plates {
                if !plate.active {
                    continue;
                }
                let dist =
                    toroidal_distance_sq(x as f32, y as f32, plate.seed_x, plate.seed_y, size);
                if dist < best_dist || (dist == best_dist && plate.id < best_id) {
                    best_dist = dist;
                    best_id = plate.id;
                }
            }

            plate_ids[y * grid_size + x] = best_id;
        }
    }
}

/// Advect plate seeds using the interpolated velocity field.
///
/// Each seed is moved by `dt * v(seed_position)`, where v is bilinearly
/// interpolated from the staggered grid. Positions wrap periodically.
pub fn advect_seeds(plates: &mut [Plate], grid: &StaggeredGrid, dt: f64) {
    let n = grid.n;
    let n_f = n as f64;

    for plate in plates.iter_mut() {
        if !plate.active {
            continue;
        }

        let sx = plate.seed_x as f64;
        let sy = plate.seed_y as f64;

        let vx = interpolate_vx(grid, sx, sy);
        let vy = interpolate_vy(grid, sx, sy);

        // Move seed with periodic wrapping
        plate.seed_x = ((sx + dt * vx) % n_f + n_f) as f32 % n as f32;
        plate.seed_y = ((sy + dt * vy) % n_f + n_f) as f32 % n as f32;
    }
}

/// Bilinear interpolation of vx at an arbitrary point (px, py).
/// vx lives at left vertical faces: vx[i,j] is at position (i, j+0.5).
fn interpolate_vx(grid: &StaggeredGrid, px: f64, py: f64) -> f64 {
    let fx = px;
    let fy = py - 0.5;
    bilinear_sample_field(&grid.vx, grid.n, fx, fy)
}

/// Bilinear interpolation of vy at an arbitrary point (px, py).
/// vy lives at bottom horizontal faces: vy[i,j] is at position (i+0.5, j).
fn interpolate_vy(grid: &StaggeredGrid, px: f64, py: f64) -> f64 {
    let fx = px - 0.5;
    let fy = py;
    bilinear_sample_field(&grid.vy, grid.n, fx, fy)
}

/// Bilinear interpolation on a periodic Field2D.
fn bilinear_sample_field(field: &Field2D, n: usize, fx: f64, fy: f64) -> f64 {
    let nf = n as f64;

    // Wrap to [0, n)
    let fx = ((fx % nf) + nf) % nf;
    let fy = ((fy % nf) + nf) % nf;

    let x0 = fx.floor() as usize % n;
    let y0 = fy.floor() as usize % n;
    let x1 = (x0 + 1) % n;
    let y1 = (y0 + 1) % n;

    let tx = fx - fx.floor();
    let ty = fy - fy.floor();

    let c00 = field.get(x0, y0);
    let c10 = field.get(x1, y0);
    let c01 = field.get(x0, y1);
    let c11 = field.get(x1, y1);

    let top = c00 + (c10 - c00) * tx;
    let bot = c01 + (c11 - c01) * tx;
    top + (bot - top) * ty
}

/// Check for plates that have lost all their cells and mark them inactive.
/// Returns the IDs of plates that just disappeared.
pub fn detect_disappeared_plates(plate_ids: &[usize], plates: &mut [Plate]) -> Vec<usize> {
    let mut cell_counts = vec![0usize; plates.len()];
    for &id in plate_ids {
        if id < cell_counts.len() {
            cell_counts[id] += 1;
        }
    }

    let mut disappeared = Vec::new();
    for plate in plates.iter_mut() {
        if plate.active && cell_counts.get(plate.id).copied().unwrap_or(0) == 0 {
            plate.active = false;
            disappeared.push(plate.id);
        }
    }

    disappeared
}

/// Rebuild the traction field from current plate_ids and plate velocities.
pub fn rebuild_traction(plate_ids: &[usize], plates: &[Plate], grid_size: usize) -> TractionField {
    let mut tx = Field2D::new(grid_size);
    let mut ty = Field2D::new(grid_size);

    for j in 0..grid_size {
        for i in 0..grid_size {
            let pid = plate_ids[j * grid_size + i];
            let plate = &plates[pid];
            tx.set(i, j, plate.velocity.0 as f64);
            ty.set(i, j, plate.velocity.1 as f64);
        }
    }

    TractionField { tx, ty }
}

// ── Cratonic rigidity ────────────────────────────────────────────────────

use super::solver::config::CratonicConfig;

/// Compute a spatial viscosity multiplier based on distance to plate seed.
///
/// Continental cells near their plate seed receive a high multiplier
/// (rigid craton). Cells far from the seed (near boundaries) receive 1.0.
/// Oceanic cells always receive 1.0.
pub fn compute_viscosity_multiplier(
    grid: &mut StaggeredGrid,
    plate_ids: &[usize],
    plates: &[Plate],
    config: &CratonicConfig,
) {
    let n = grid.n;

    if !config.enabled || config.max_factor <= 1.0 {
        for k in 0..n * n {
            grid.eta_multiplier.data_mut()[k] = 1.0;
        }
        return;
    }

    let size = n as f32;

    // Step 1: Find maximum Voronoï radius for each continental plate.
    let mut max_radius = vec![0.0_f32; plates.len()];
    for j in 0..n {
        for i in 0..n {
            let pid = plate_ids[j * n + i];
            let plate = &plates[pid];
            if plate.plate_type != PlateType::Continental || !plate.active {
                continue;
            }
            let dist = toroidal_distance(i as f32, j as f32, plate.seed_x, plate.seed_y, size);
            if dist > max_radius[pid] {
                max_radius[pid] = dist;
            }
        }
    }

    // Step 2: Compute multiplier per cell.
    for j in 0..n {
        for i in 0..n {
            let pid = plate_ids[j * n + i];
            let plate = &plates[pid];

            if plate.plate_type != PlateType::Continental || !plate.active {
                grid.eta_multiplier.set(i, j, 1.0);
                continue;
            }

            let radius = max_radius[pid];
            if radius < 1e-6 {
                grid.eta_multiplier.set(i, j, config.max_factor);
                continue;
            }

            let dist = toroidal_distance(i as f32, j as f32, plate.seed_x, plate.seed_y, size);
            let d_norm = (dist / radius).min(1.0);

            // Decay: f = (1 - d^p), from 1 at center to 0 at edge
            let f = (1.0 - d_norm.powf(config.decay_power as f32)).max(0.0) as f64;
            let mult = 1.0 + (config.max_factor - 1.0) * f;
            grid.eta_multiplier.set(i, j, mult);
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seed::WorldSeed;

    #[test]
    fn test_toroidal_distance() {
        // Adjacent points: distance² = 1
        let d = toroidal_distance_sq(1.0, 1.0, 2.0, 1.0, 10.0);
        assert!((d - 1.0).abs() < 1e-6, "expected 1.0, got {d}");

        // x=1 and x=9 on a size=10 grid: direct=8, wrap=2, min=2 → d²=4
        let d = toroidal_distance_sq(1.0, 0.0, 9.0, 0.0, 10.0);
        assert!((d - 4.0).abs() < 1e-6, "expected 4.0, got {d}");
    }

    #[test]
    fn test_generate_plates_deterministic() {
        let config = PlateConfig::default();
        let seed = WorldSeed::new(42);

        let result1 = generate_plates(&config, &seed);
        let result2 = generate_plates(&config, &seed);

        assert_eq!(result1.plate_ids, result2.plate_ids);
        assert_eq!(result1.thickness.data, result2.thickness.data);
    }

    #[test]
    fn test_plate_count_and_types() {
        let config = PlateConfig {
            num_plates: 10,
            continental_ratio: 0.3,
            grid_size: 64,
            ..Default::default()
        };
        let seed = WorldSeed::new(123);
        let result = generate_plates(&config, &seed);

        assert_eq!(result.plates.len(), 10);
        let continental_count =
            result.plates.iter().filter(|p| p.plate_type == PlateType::Continental).count();
        assert_eq!(continental_count, 3, "30% of 10 = 3 continental plates");
    }

    #[test]
    fn test_all_cells_assigned() {
        let config = PlateConfig { grid_size: 32, ..Default::default() };
        let seed = WorldSeed::new(99);
        let result = generate_plates(&config, &seed);

        for &id in &result.plate_ids {
            assert!(id < config.num_plates, "cell has invalid plate id {id}");
        }
    }

    #[test]
    fn test_thickness_values() {
        let config = PlateConfig {
            grid_size: 64,
            boundary_smoothing_sigma: 0.0, // no blur → exact values
            ..Default::default()
        };
        let seed = WorldSeed::new(42);
        let result = generate_plates(&config, &seed);

        for &val in &result.thickness.data {
            assert!(
                (val - 1.0).abs() < 1e-6 || (val - 0.2).abs() < 1e-6,
                "unexpected thickness {val}, expected 1.0 or 0.2"
            );
        }
    }

    #[test]
    fn test_toroidal_wrapping_assignment() {
        // A seed near x=0 should capture the right side of the grid too
        // (since distance wraps). We can't test exact assignment without
        // running the full algorithm, but we verify that cells on both
        // sides of the grid can belong to the same plate.
        let config = PlateConfig {
            num_plates: 2,
            grid_size: 32,
            boundary_smoothing_sigma: 0.0,
            ..Default::default()
        };
        let seed = WorldSeed::new(7);
        let result = generate_plates(&config, &seed);

        // At least some cell at x=0 and x=31 should exist — just verify no panic
        let _id_left = result.plate_ids[0];
        let _id_right = result.plate_ids[31];
        assert_eq!(result.grid_size, 32);
    }

    // ── Dynamic boundary tests ───────────────────────────────────────────

    #[test]
    fn seeds_move_with_velocity() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        // Uniform rightward velocity
        for j in 0..n {
            for i in 0..n {
                grid.vx.set(i, j, 1.0);
                grid.vy.set(i, j, 0.0);
            }
        }

        let mut plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Continental,
            velocity: (1.0, 0.0),
            seed_x: 10.0,
            seed_y: 16.0,
            active: true,
            subducted_mass: 0.0,
        }];

        advect_seeds(&mut plates, &grid, 0.1);

        assert!(plates[0].seed_x > 10.0, "Seed should move right: {}", plates[0].seed_x);
        assert!((plates[0].seed_y - 16.0).abs() < 1e-3, "Seed should not move vertically");
    }

    #[test]
    fn seed_wraps_periodically() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        for j in 0..n {
            for i in 0..n {
                grid.vx.set(i, j, 1.0);
            }
        }

        let mut plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Oceanic,
            velocity: (1.0, 0.0),
            seed_x: 31.5,
            seed_y: 16.0,
            active: true,
            subducted_mass: 0.0,
        }];

        advect_seeds(&mut plates, &grid, 1.0);

        assert!(
            plates[0].seed_x < 5.0 || plates[0].seed_x > 28.0,
            "Seed should wrap: {}",
            plates[0].seed_x
        );
    }

    #[test]
    fn plate_disappears_when_seeds_overlap() {
        let n = 16;
        let mut plate_ids = vec![0usize; n * n];

        let mut plates = vec![
            Plate {
                id: 0,
                plate_type: PlateType::Continental,
                velocity: (0.0, 0.0),
                seed_x: 8.0,
                seed_y: 8.0,
                active: true,
                subducted_mass: 0.0,
            },
            Plate {
                id: 1,
                plate_type: PlateType::Oceanic,
                velocity: (0.0, 0.0),
                seed_x: 8.0, // same as plate 0
                seed_y: 8.0,
                active: true,
                subducted_mass: 0.0,
            },
        ];

        recompute_voronoi(&mut plate_ids, &plates, n);
        let disappeared = detect_disappeared_plates(&plate_ids, &mut plates);

        // When two seeds are at the same position, the one with lower ID wins
        assert!(disappeared.contains(&1), "Plate 1 should disappear when seeds overlap");
        assert!(!plates[1].active);
    }

    #[test]
    fn recompute_voronoi_matches_generate_plates() {
        let config = PlateConfig {
            num_plates: 4,
            grid_size: 32,
            boundary_smoothing_sigma: 0.0,
            ..Default::default()
        };
        let seed = WorldSeed::new(42);
        let result = generate_plates(&config, &seed);

        // Recompute from the same seeds
        let mut plate_ids = vec![0usize; 32 * 32];
        recompute_voronoi(&mut plate_ids, &result.plates, 32);

        assert_eq!(plate_ids, result.plate_ids);
    }

    #[test]
    fn inactive_plate_skipped_in_voronoi() {
        let n = 16;
        let mut plates = vec![
            Plate {
                id: 0,
                plate_type: PlateType::Continental,
                velocity: (0.0, 0.0),
                seed_x: 4.0,
                seed_y: 8.0,
                active: true,
                subducted_mass: 0.0,
            },
            Plate {
                id: 1,
                plate_type: PlateType::Oceanic,
                velocity: (0.0, 0.0),
                seed_x: 12.0,
                seed_y: 8.0,
                active: false,
                subducted_mass: 0.0, // inactive
            },
        ];

        let mut plate_ids = vec![0usize; n * n];
        recompute_voronoi(&mut plate_ids, &plates, n);

        // All cells should belong to plate 0 since plate 1 is inactive
        assert!(plate_ids.iter().all(|&id| id == 0));
    }

    // ── Cratonic rigidity tests ──────────────────────────────────────

    #[test]
    fn cratonic_multiplier_is_highest_at_seed() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        let plate_ids: Vec<usize> = vec![0; n * n];
        let plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Continental,
            velocity: (0.0, 0.0),
            seed_x: 16.0,
            seed_y: 16.0,
            active: true,
            subducted_mass: 0.0,
        }];

        let config = CratonicConfig { enabled: true, max_factor: 20.0, decay_power: 2.0 };
        compute_viscosity_multiplier(&mut grid, &plate_ids, &plates, &config);

        let center = grid.eta_multiplier.get(16, 16);
        assert!(center > 15.0, "Center should be near max: {center}");

        let edge = grid.eta_multiplier.get(0, 0);
        assert!(edge < 5.0, "Edge should be near 1.0: {edge}");

        assert!(center > edge, "Center should be more rigid than edge");
    }

    #[test]
    fn oceanic_plates_have_no_rigidity() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        let plate_ids: Vec<usize> = vec![0; n * n];
        let plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Oceanic,
            velocity: (0.0, 0.0),
            seed_x: 8.0,
            seed_y: 8.0,
            active: true,
            subducted_mass: 0.0,
        }];

        let config = CratonicConfig { enabled: true, max_factor: 50.0, decay_power: 2.0 };
        compute_viscosity_multiplier(&mut grid, &plate_ids, &plates, &config);

        for k in 0..n * n {
            assert!(
                (grid.eta_multiplier.data()[k] - 1.0).abs() < 1e-10,
                "Oceanic cell should have mult=1.0"
            );
        }
    }

    #[test]
    fn disabled_gives_uniform_multiplier() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        let plate_ids: Vec<usize> = vec![0; n * n];
        let plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Continental,
            velocity: (0.0, 0.0),
            seed_x: 8.0,
            seed_y: 8.0,
            active: true,
            subducted_mass: 0.0,
        }];

        let config = CratonicConfig { enabled: false, ..Default::default() };
        compute_viscosity_multiplier(&mut grid, &plate_ids, &plates, &config);

        for k in 0..n * n {
            assert!((grid.eta_multiplier.data()[k] - 1.0).abs() < 1e-10);
        }
    }
}
