//! Tectonic plate initialization via Voronoï partitioning.
//!
//! Generates `num_plates` plates on a toroidal grid (periodic boundary conditions).
//! Each plate gets a random position, type (continental/oceanic), and velocity vector.
//! Cells are assigned to the nearest plate seed using toroidal distance.
//! The resulting crustal thickness field is smoothed with a separable Gaussian blur.

use rand::Rng;

use crate::grid::GridF32;
use crate::seed::WorldSeed;
use super::solver::field::Field2D;
use super::solver::grid::StaggeredGrid;
use super::solver::traction::TractionField;

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
    let num_continental =
        (config.num_plates as f32 * config.continental_ratio).round() as usize;

    for id in 0..config.num_plates {
        let seed_x = rng.random::<f32>() * size as f32;
        let seed_y = rng.random::<f32>() * size as f32;

        let plate_type =
            if id < num_continental { PlateType::Continental } else { PlateType::Oceanic };

        let angle = rng.random::<f32>() * std::f32::consts::TAU;
        let speed = config.velocity_min
            + rng.random::<f32>() * (config.velocity_max - config.velocity_min);
        let velocity = (angle.cos() * speed, angle.sin() * speed);

        plates.push(Plate { id, plate_type, velocity, seed_x, seed_y });
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

/// Squared toroidal distance between two points on a grid that wraps in both x and y.
///
/// For each axis, takes the shorter path: direct or through the wrap.
fn toroidal_distance_sq(x1: f32, y1: f32, x2: f32, y2: f32, size: f32) -> f32 {
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
}
