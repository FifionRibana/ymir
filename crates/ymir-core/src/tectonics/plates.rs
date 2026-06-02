//! Tectonic plate initialization via Voronoï partitioning.
//!
//! Generates `num_plates` plates on a toroidal grid (periodic boundary conditions).
//! Each plate gets a random position, type (continental/oceanic), and velocity vector.
//! Cells are assigned to the nearest plate seed using toroidal distance.
//! The resulting crustal thickness field is smoothed with a separable Gaussian blur.

use rand::Rng;
use tracing::info;

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
    // ── Runtime stats (updated by update_plate_stats) ──
    /// Number of cells belonging to this plate.
    pub cell_count: usize,
    /// Mean crustal thickness across all cells.
    pub mean_thickness: f32,
    /// Mean velocity across all cells.
    pub mean_velocity: (f32, f32),
    /// Centroid position (circular mean on torus).
    pub centroid_x: f32,
    pub centroid_y: f32,
}

impl Plate {
    /// Create a new plate with zero runtime stats.
    pub fn new(
        id: usize,
        plate_type: PlateType,
        velocity: (f32, f32),
        seed_x: f32,
        seed_y: f32,
    ) -> Self {
        Self {
            id,
            plate_type,
            velocity,
            seed_x,
            seed_y,
            active: true,
            subducted_mass: 0.0,
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: seed_x,
            centroid_y: seed_y,
        }
    }
}

/// Reduction factor applied to the raw continental/total plate ratio when
/// computing initial continental area. Without this factor, a 2/12 ratio
/// would give 16.7% continental area, which grows too large after mass
/// accumulates during simulation. 0.6 keeps initial continents compact
/// enough that the final area stays in a geologically reasonable range
/// after several hundred time steps.
pub const BASE_CONTINENTAL_RATIO: f64 = 0.6;

/// Configuration for plate generation.
///
/// `#[serde(default)]` (struct-level): metadata saved before a field
/// was added still deserializes — missing fields fall back to
/// [`PlateConfig::default`]. Guards against the #47-class legacy-
/// compat break (e.g. `continental_area_factor` added later without a
/// default broke `deserialize_legacy_metadata_without_upscale`).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct PlateConfig {
    /// Number of tectonic plates (5–15).
    pub num_plates: usize,
    /// Number of continental plates (0..=num_plates/2). The continental
    /// area fraction is derived as
    /// `(num_continental_plates / num_plates) * continental_area_factor`.
    pub num_continental_plates: usize,
    /// Land/sea ratio on continental plates: fraction of a continental
    /// plate's Voronoï region that becomes landmass (continental thickness).
    /// Remaining cells on the plate keep oceanic thickness (continental shelf).
    /// Range 0.1–1.0, default 0.6.
    pub continental_area_factor: f32,
    /// Minimum velocity magnitude (grid units / timestep).
    pub velocity_min: f32,
    /// Maximum velocity magnitude (grid units / timestep).
    pub velocity_max: f32,
    /// Grid width (x-axis cell count) for the tectonic simulation.
    pub grid_width: usize,
    /// Grid height (y-axis cell count) for the tectonic simulation.
    pub grid_height: usize,
    /// Gaussian blur sigma applied to initial thickness field to smooth plate boundaries.
    pub boundary_smoothing_sigma: f32,
}

impl Default for PlateConfig {
    fn default() -> Self {
        Self {
            num_plates: 8,
            num_continental_plates: 2,
            continental_area_factor: BASE_CONTINENTAL_RATIO as f32,
            velocity_min: 0.5,
            velocity_max: 2.5,
            grid_width: 128,
            grid_height: 128,
            boundary_smoothing_sigma: 2.0,
        }
    }
}

impl PlateConfig {
    /// Single large continent surrounded by ocean.
    /// Good starting point for a typical game continent.
    pub fn preset_single_continent() -> Self {
        Self { num_plates: 8, num_continental_plates: 1, ..Default::default() }
    }

    /// Two continental masses converging — produces a Himalaya-type collision belt.
    pub fn preset_collision() -> Self {
        Self { num_plates: 6, num_continental_plates: 2, ..Default::default() }
    }

    /// Many small continental plates — creates an archipelago with multiple islands.
    pub fn preset_archipelago() -> Self {
        Self { num_plates: 12, num_continental_plates: 3, velocity_max: 3.0, ..Default::default() }
    }

    /// One plate splitting apart — produces a rift valley with volcanism and graben basins.
    pub fn preset_rift() -> Self {
        Self {
            num_plates: 6,
            num_continental_plates: 2,
            velocity_min: 0.3,
            velocity_max: 1.5,
            ..Default::default()
        }
    }

    /// Set grid dimensions from a resolution (treated as width) and an aspect
    /// ratio. Height is derived as `round(resolution / aspect_ratio)`, clamped
    /// to at least 1. Examples:
    ///
    /// - `resolution=128, aspect=1.0`      → 128 × 128 (square)
    /// - `resolution=128, aspect=1.5`      → 128 × 85  (3:2)
    /// - `resolution=128, aspect=16.0/9.0` → 128 × 72  (16:9)
    ///
    /// The aspect ratio is consumed at construction time and not stored.
    pub fn with_resolution_aspect(mut self, resolution: usize, aspect_ratio: f32) -> Self {
        self.grid_width = resolution;
        self.grid_height =
            ((resolution as f32) / aspect_ratio).round().max(1.0) as usize;
        self
    }

    /// Set explicit grid width and height. Advanced mode — use this when the
    /// desired shape does not match a named aspect ratio.
    pub fn with_dimensions(mut self, width: usize, height: usize) -> Self {
        self.grid_width = width;
        self.grid_height = height;
        self
    }
}

/// Result of plate initialization — everything needed to start the tectonic simulation.
#[derive(Clone)]
pub struct PlateInitResult {
    /// Which plate each cell belongs to (grid_width × grid_height cells, row-major).
    pub plate_ids: Vec<usize>,
    /// Initial crustal thickness field. Continental ≈ 1.0, Oceanic ≈ 0.2.
    pub thickness: GridF32,
    /// The plates themselves with their properties.
    pub plates: Vec<Plate>,
    /// Grid width (x-axis cell count).
    pub grid_width: usize,
    /// Grid height (y-axis cell count).
    pub grid_height: usize,
}

// ── Core algorithm ────────────────────────────────────────────────────────

/// Generate plates on the grid using Voronoï partitioning.
///
/// The grid has periodic (toroidal) boundary conditions: the left edge connects
/// to the right, the top to the bottom. The distance metric wraps around the grid
/// in both directions to find the shortest path.
pub fn generate_plates(config: &PlateConfig, seed: &WorldSeed) -> PlateInitResult {
    let mut rng = seed.rng_for("plates");
    let nx = config.grid_width;
    let ny = config.grid_height;

    // 1. Place Voronoï seeds randomly on the grid.
    // Clamp the continental plate count defensively. UI enforces
    // num_continental_plates <= num_plates / 2, but we guard here in case
    // the config was constructed programmatically.
    let mut plates = Vec::with_capacity(config.num_plates);
    let num_continental = config.num_continental_plates.min(config.num_plates / 2);

    for id in 0..config.num_plates {
        let seed_x = rng.random::<f32>() * nx as f32;
        let seed_y = rng.random::<f32>() * ny as f32;

        let plate_type =
            if id < num_continental { PlateType::Continental } else { PlateType::Oceanic };

        let angle = rng.random::<f32>() * std::f32::consts::TAU;
        let speed =
            config.velocity_min + rng.random::<f32>() * (config.velocity_max - config.velocity_min);
        let velocity = (angle.cos() * speed, angle.sin() * speed);

        plates.push(Plate::new(id, plate_type, velocity, seed_x, seed_y));
    }

    // 2. Assign each cell to the nearest plate seed (toroidal distance)
    let mut plate_ids = vec![0usize; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let mut best_dist = f32::MAX;
            let mut best_id = 0;

            for plate in &plates {
                let dist = toroidal_distance_sq_2d(
                    x as f32,
                    y as f32,
                    plate.seed_x,
                    plate.seed_y,
                    nx as f32,
                    ny as f32,
                );
                if dist < best_dist {
                    best_dist = dist;
                    best_id = plate.id;
                }
            }

            plate_ids[y * nx + x] = best_id;
        }
    }

    // 3. Initialize crustal thickness field.
    // Derive the total continental cell fraction from the plate ratio.
    // This scales naturally: doubling num_plates while keeping
    // num_continental_plates fixed halves the continental area.
    let continental_fraction = if config.num_plates > 0 {
        (num_continental as f64 / config.num_plates as f64) * config.continental_area_factor as f64
    } else {
        0.0
    };
    let total_cells = nx * ny;
    let target_total = (continental_fraction * total_cells as f64) as usize;
    let cells_per_continent = if num_continental > 0 { target_total / num_continental } else { 0 };

    // Default every cell to oceanic thickness.
    let mut thickness = GridF32::new(nx, ny, 0.2);

    // For each continental plate, find all cells in its Voronoï region,
    // sort them by toroidal distance to the seed, and give continental
    // thickness only to the closest `cells_per_continent`. Remaining cells
    // in the same region stay on the plate but receive oceanic thickness,
    // effectively forming the continental shelf around the landmass.
    for plate_idx in 0..num_continental {
        let sx = plates[plate_idx].seed_x;
        let sy = plates[plate_idx].seed_y;

        let mut plate_cells: Vec<(usize, f32)> = Vec::new();
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                if plate_ids[k] == plate_idx {
                    let dist_sq = toroidal_distance_sq_2d(
                        x as f32,
                        y as f32,
                        sx,
                        sy,
                        nx as f32,
                        ny as f32,
                    );
                    plate_cells.push((k, dist_sq));
                }
            }
        }

        // Closest cells first — these form the compact continental core.
        plate_cells.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let n_core = cells_per_continent.min(plate_cells.len());
        for (idx, &(k, _)) in plate_cells.iter().enumerate() {
            if idx < n_core {
                let y = k / nx;
                let x = k % nx;
                thickness.set(x, y, 1.0);
            }
            // else: already 0.2 from the default fill — continental shelf.
        }
    }

    // 4. Smooth boundaries with Gaussian blur (toroidal wrapping)
    if config.boundary_smoothing_sigma > 0.0 {
        thickness = gaussian_blur(&thickness, config.boundary_smoothing_sigma);
    }

    PlateInitResult { plate_ids, thickness, plates, grid_width: nx, grid_height: ny }
}

// ── Distance and blur helpers ─────────────────────────────────────────────

/// Toroidal distance between two points on a rectangular grid that wraps
/// in both x and y. `period_x` is the x-axis wrap (grid_width), `period_y`
/// is the y-axis wrap (grid_height).
pub fn toroidal_distance_2d(
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    period_x: f32,
    period_y: f32,
) -> f32 {
    toroidal_distance_sq_2d(x1, y1, x2, y2, period_x, period_y).sqrt()
}

/// Squared toroidal distance between two points on a rectangular grid.
/// For each axis, takes the shorter path: direct or through the wrap.
pub fn toroidal_distance_sq_2d(
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    period_x: f32,
    period_y: f32,
) -> f32 {
    let dx = (x1 - x2).abs();
    let dy = (y1 - y2).abs();
    let dx = dx.min(period_x - dx);
    let dy = dy.min(period_y - dy);
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
        let nx = self.grid_width;
        let ny = self.grid_height;
        let mut tx = Field2D::new(nx, ny);
        let mut ty = Field2D::new(nx, ny);

        for j in 0..ny {
            for i in 0..nx {
                let plate_id = self.plate_ids[j * nx + i];
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
    /// Creates the grid with the correct dimensions and dx, copies the
    /// thickness field (f32 → f64), and returns the grid ready for the
    /// solver. Cell spacing `dx` is derived from the grid width so that
    /// the horizontal domain length stays at 1.0 regardless of aspect
    /// ratio.
    pub fn to_staggered_grid(&self) -> StaggeredGrid {
        let nx = self.grid_width;
        let ny = self.grid_height;
        let dx = 1.0 / nx as f64;
        let mut grid = StaggeredGrid::new(nx, ny, dx);

        for j in 0..ny {
            for i in 0..nx {
                let s = self.thickness.data[j * nx + i] as f64;
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
pub fn recompute_voronoi(
    plate_ids: &mut [usize],
    plates: &[Plate],
    grid_width: usize,
    grid_height: usize,
) {
    let w = grid_width as f32;
    let h = grid_height as f32;
    for y in 0..grid_height {
        for x in 0..grid_width {
            let mut best_dist = f32::MAX;
            let mut best_id = 0;

            for plate in plates {
                if !plate.active {
                    continue;
                }
                let dist = toroidal_distance_sq_2d(
                    x as f32,
                    y as f32,
                    plate.seed_x,
                    plate.seed_y,
                    w,
                    h,
                );
                if dist < best_dist || (dist == best_dist && plate.id < best_id) {
                    best_dist = dist;
                    best_id = plate.id;
                }
            }

            plate_ids[y * grid_width + x] = best_id;
        }
    }
}

/// Advect plate seeds using the interpolated velocity field.
///
/// Each seed is moved by `dt * v(seed_position)`, where v is bilinearly
/// interpolated from the staggered grid. Positions wrap periodically.
pub fn advect_seeds(plates: &mut [Plate], grid: &StaggeredGrid, dt: f64) {
    let nx = grid.nx();
    let ny = grid.ny();
    let nx_f = nx as f64;
    let ny_f = ny as f64;

    for plate in plates.iter_mut() {
        if !plate.active {
            continue;
        }

        let sx = plate.seed_x as f64;
        let sy = plate.seed_y as f64;

        let vx = interpolate_vx(grid, sx, sy);
        let vy = interpolate_vy(grid, sx, sy);

        // Move seed with periodic wrapping
        plate.seed_x = ((sx + dt * vx) % nx_f + nx_f) as f32 % nx as f32;
        plate.seed_y = ((sy + dt * vy) % ny_f + ny_f) as f32 % ny as f32;
    }
}

/// Advect plate IDs using accumulated sub-pixel displacement.
///
/// Each cell accumulates its forward displacement (dt × v) in `disp_x` and
/// `disp_y`. When the integer part of the accumulator is non-zero, the cell's
/// plate ID is replaced by tracing backward by that integer shift, and the
/// fractional remainder is kept in the accumulator for the next step.
///
/// This allows boundaries to move even when per-step displacement is a small
/// fraction of a pixel (e.g. 0.008 px/step → 1-pixel shift every ~125 steps).
pub fn advect_plate_ids(
    ids: &mut [usize],
    disp_x: &mut Field2D,
    disp_y: &mut Field2D,
    grid: &StaggeredGrid,
    dt: f64,
) {
    let nx = grid.nx();
    let ny = grid.ny();

    // Phase 1: accumulate forward displacement (dt × v) at each cell
    for j in 0..ny {
        for i in 0..nx {
            let vx = interpolate_vx(grid, i as f64, j as f64);
            let vy = interpolate_vy(grid, i as f64, j as f64);

            disp_x.set(i, j, disp_x.get(i, j) + dt * vx);
            disp_y.set(i, j, disp_y.get(i, j) + dt * vy);
        }
    }

    // Phase 2: shift IDs where the integer part of the accumulator is non-zero.
    // disp_x > 0 means material moved right, so the new ID at (i,j) comes from
    // the cell to the LEFT (i - shift_x). Work on a copy to avoid order artifacts.
    let old_ids = ids.to_vec();

    for j in 0..ny {
        for i in 0..nx {
            let dx = disp_x.get(i, j);
            let dy = disp_y.get(i, j);

            let shift_x = dx.round() as i64;
            let shift_y = dy.round() as i64;

            if shift_x == 0 && shift_y == 0 {
                continue;
            }

            let src_i = ((i as i64 - shift_x).rem_euclid(nx as i64)) as usize;
            let src_j = ((j as i64 - shift_y).rem_euclid(ny as i64)) as usize;

            ids[j * nx + i] = old_ids[src_j * nx + src_i];

            // Keep the fractional remainder
            disp_x.set(i, j, dx - shift_x as f64);
            disp_y.set(i, j, dy - shift_y as f64);
        }
    }
}

/// Morphological cleanup of plate IDs after advection.
///
/// At convergent boundaries, advection can create thick "mixed zones" where
/// plate IDs alternate in a checkerboard pattern. This function reassigns
/// cells that are isolated or in thin protrusions to the majority plate
/// among their 8-connected neighbors, keeping boundaries sharp at 1-cell width.
pub fn cleanup_plate_ids(ids: &mut [usize], nx: usize, ny: usize) {
    let old_ids = ids.to_vec();

    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let my_id = old_ids[k];

            // Count neighbors per plate ID; supports up to 256 plates
            let mut counts = [0u8; 256];
            let mut same_count = 0u8;

            for &(di, dj) in
                &[(-1isize, -1isize), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
            {
                let ni = ((i as isize + di).rem_euclid(nx as isize)) as usize;
                let nj = ((j as isize + dj).rem_euclid(ny as isize)) as usize;
                let nid = old_ids[nj * nx + ni];
                if nid < counts.len() {
                    counts[nid] += 1;
                }
                if nid == my_id {
                    same_count += 1;
                }
            }

            // Isolated cells (≤1 same) or thin protrusions (2 same) → reassign
            if same_count <= 2 {
                let mut best_id = my_id;
                let mut best_count = 0u8;
                for (pid, &c) in counts.iter().enumerate() {
                    if c > best_count {
                        best_count = c;
                        best_id = pid;
                    }
                }
                ids[k] = best_id;
            }
        }
    }
}

/// Compute per-plate runtime statistics from the cell data.
///
/// Updates `cell_count`, `mean_thickness`, `mean_velocity`, and
/// `centroid_x`/`centroid_y` (circular mean on the torus) for each plate.
pub fn update_plate_stats(ids: &[usize], plates: &mut [Plate], grid: &StaggeredGrid) {
    let nx = grid.nx();
    let ny = grid.ny();
    let tau = std::f64::consts::TAU;

    // Circular mean accumulators (local, not stored in Plate)
    let mut sin_x = vec![0.0_f64; plates.len()];
    let mut cos_x = vec![0.0_f64; plates.len()];
    let mut sin_y = vec![0.0_f64; plates.len()];
    let mut cos_y = vec![0.0_f64; plates.len()];

    for plate in plates.iter_mut() {
        plate.cell_count = 0;
        plate.mean_thickness = 0.0;
        plate.mean_velocity = (0.0, 0.0);
    }

    for j in 0..ny {
        for i in 0..nx {
            let pid = ids[j * nx + i];
            if pid >= plates.len() {
                continue;
            }
            let plate = &mut plates[pid];
            plate.cell_count += 1;
            plate.mean_thickness += grid.s.get(i, j) as f32;

            let vx = interpolate_vx(grid, i as f64, j as f64);
            let vy = interpolate_vy(grid, i as f64, j as f64);
            plate.mean_velocity.0 += vx as f32;
            plate.mean_velocity.1 += vy as f32;

            let theta_x = tau * i as f64 / nx as f64;
            let theta_y = tau * j as f64 / ny as f64;
            sin_x[pid] += theta_x.sin();
            cos_x[pid] += theta_x.cos();
            sin_y[pid] += theta_y.sin();
            cos_y[pid] += theta_y.cos();
        }
    }

    for (pid, plate) in plates.iter_mut().enumerate() {
        if plate.cell_count > 0 {
            let c = plate.cell_count as f32;
            plate.mean_thickness /= c;
            plate.mean_velocity.0 /= c;
            plate.mean_velocity.1 /= c;
            plate.centroid_x = (sin_x[pid].atan2(cos_x[pid]) * nx as f64 / tau) as f32;
            plate.centroid_y = (sin_y[pid].atan2(cos_y[pid]) * ny as f64 / tau) as f32;
            // Wrap to [0, nx) and [0, ny)
            plate.centroid_x = ((plate.centroid_x % nx as f32) + nx as f32) % nx as f32;
            plate.centroid_y = ((plate.centroid_y % ny as f32) + ny as f32) % ny as f32;
        }
    }
}

/// Bilinear interpolation of vx at an arbitrary point (px, py).
/// vx lives at left vertical faces: vx[i,j] is at position (i, j+0.5).
pub fn interpolate_vx(grid: &StaggeredGrid, px: f64, py: f64) -> f64 {
    let fx = px;
    let fy = py - 0.5;
    bilinear_sample_field(&grid.vx, grid.nx(), grid.ny(), fx, fy)
}

/// Bilinear interpolation of vy at an arbitrary point (px, py).
/// vy lives at bottom horizontal faces: vy[i,j] is at position (i+0.5, j).
pub fn interpolate_vy(grid: &StaggeredGrid, px: f64, py: f64) -> f64 {
    let fx = px - 0.5;
    let fy = py;
    bilinear_sample_field(&grid.vy, grid.nx(), grid.ny(), fx, fy)
}

/// Bilinear interpolation on a periodic Field2D.
fn bilinear_sample_field(field: &Field2D, nx: usize, ny: usize, fx: f64, fy: f64) -> f64 {
    let nxf = nx as f64;
    let nyf = ny as f64;

    // Wrap to [0, nx) × [0, ny)
    let fx = ((fx % nxf) + nxf) % nxf;
    let fy = ((fy % nyf) + nyf) % nyf;

    let x0 = fx.floor() as usize % nx;
    let y0 = fy.floor() as usize % ny;
    let x1 = (x0 + 1) % nx;
    let y1 = (y0 + 1) % ny;

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
/// At convergent boundaries, cells thin enough on the subducting side
/// are reassigned to the overriding plate's ID.
pub fn apply_subduction_consumption(
    ids: &mut [usize],
    grid: &StaggeredGrid,
    boundary_field: &crate::tectonics::boundaries::BoundaryField,
    threshold: f64,
) {
    use crate::tectonics::boundaries::BoundaryType;
    let nx = grid.nx();
    let ny = grid.ny();
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();

    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let btype = &boundary_field.boundary_type[k];

            let is_subducting =
                matches!(btype, BoundaryType::Subduction | BoundaryType::OceanicSubduction)
                    && grid.s.get(i, j) < threshold;

            if !is_subducting {
                continue;
            }

            let my_id = ids[k];
            let neighbors = [
                (idx_x.next(i), j),
                (idx_x.prev(i), j),
                (i, idx_y.next(j)),
                (i, idx_y.prev(j)),
            ];

            for &(ni, nj) in &neighbors {
                let nk = nj * nx + ni;
                if ids[nk] != my_id && grid.s.get(ni, nj) > grid.s.get(i, j) {
                    ids[k] = ids[nk];
                    break;
                }
            }
        }
    }
}

/// At divergent boundaries, very thin cells that are receiving material
/// from spreading become new plates. Only cells whose parent plate has
/// mean_thickness > 0.4 (continental being rifted apart) get new IDs.
/// Contiguous candidate cells are grouped into a single plate via BFS.
pub fn apply_rift_creation(
    ids: &mut [usize],
    plates: &mut Vec<Plate>,
    grid: &StaggeredGrid,
    boundary_field: &crate::tectonics::boundaries::BoundaryField,
    next_id: &mut usize,
    threshold: f64,
) {
    use crate::tectonics::boundaries::BoundaryType;
    let nx = grid.nx();
    let ny = grid.ny();

    // Pass 1: collect all cells that qualify for rift creation
    let mut candidates = vec![false; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            if boundary_field.boundary_type[k] != BoundaryType::Rift {
                continue;
            }
            if grid.s.get(i, j) > threshold {
                continue;
            }
            let parent_id = ids[k];
            if parent_id >= plates.len() || plates[parent_id].mean_thickness <= 0.4 {
                continue;
            }
            candidates[k] = true;
        }
    }

    // Pass 2: connected component labeling (4-connectivity, periodic)
    let mut visited = vec![false; nx * ny];
    let mut components: Vec<Vec<usize>> = Vec::new();

    for k in 0..nx * ny {
        if !candidates[k] || visited[k] {
            continue;
        }
        let mut group = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(k);
        visited[k] = true;

        while let Some(ck) = queue.pop_front() {
            group.push(ck);
            let ci = ck % nx;
            let cj = ck / nx;
            let neighbors = [
                ((ci + 1) % nx, cj),
                ((ci + nx - 1) % nx, cj),
                (ci, (cj + 1) % ny),
                (ci, (cj + ny - 1) % ny),
            ];
            for &(ni, nj) in &neighbors {
                let nk = nj * nx + ni;
                if candidates[nk] && !visited[nk] {
                    visited[nk] = true;
                    queue.push_back(nk);
                }
            }
        }
        components.push(group);
    }

    // Pass 3: one new plate per connected component
    for component in &components {
        let new_id = *next_id;
        *next_id += 1;

        let ci = component[0] % nx;
        let cj = component[0] / nx;
        plates.push(Plate::new(new_id, PlateType::Oceanic, (0.0, 0.0), ci as f32, cj as f32));

        for &k in component {
            ids[k] = new_id;
        }
    }
}

/// Detect plates that have been fragmented into disconnected components.
///
/// Connectivity uses a thickness threshold: cells with S below `connectivity_threshold`
/// do not propagate connectivity. This allows a thin rift zone to "break" a continent
/// into two pieces even though the thin cells still carry the same plate_id.
///
/// The largest component keeps the original ID, smaller ones get new IDs.
/// Thin cells in the rift zone are reassigned to the nearest fragment by adjacency.
pub fn detect_fragmentation(
    ids: &mut [usize],
    plates: &mut Vec<Plate>,
    next_id: &mut usize,
    nx: usize,
    ny: usize,
    grid: &StaggeredGrid,
    connectivity_threshold: f64,
) {
    let active_ids: Vec<usize> =
        plates.iter().filter(|p| p.active && p.cell_count > 0).map(|p| p.id).collect();

    for &pid in &active_ids {
        let cells: Vec<usize> = (0..nx * ny).filter(|&k| ids[k] == pid).collect();
        if cells.len() < 4 {
            continue;
        }
        // Only continental plates can fragment by rifting.
        if plates[pid].mean_thickness <= 0.5 {
            continue;
        }

        // BFS flood fill with thickness-aware connectivity
        let mut visited = vec![false; nx * ny];
        let mut components: Vec<Vec<usize>> = Vec::new();
        let mut thin_cells: Vec<usize> = Vec::new();

        for &start in &cells {
            if visited[start] {
                continue;
            }
            // Thin cells don't start components
            if grid.s.get(start % nx, start / nx) < connectivity_threshold {
                thin_cells.push(start);
                visited[start] = true;
                continue;
            }

            let mut component = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start);
            visited[start] = true;

            while let Some(k) = queue.pop_front() {
                component.push(k);
                let ci = k % nx;
                let cj = k / nx;
                let neighbors = [
                    ((ci + 1) % nx, cj),
                    ((ci + nx - 1) % nx, cj),
                    (ci, (cj + 1) % ny),
                    (ci, (cj + ny - 1) % ny),
                ];
                for &(ni, nj) in &neighbors {
                    let nk = nj * nx + ni;
                    if visited[nk] || ids[nk] != pid {
                        continue;
                    }
                    if grid.s.get(ni, nj) < connectivity_threshold {
                        thin_cells.push(nk);
                        visited[nk] = true;
                        continue;
                    }
                    visited[nk] = true;
                    queue.push_back(nk);
                }
            }

            if !component.is_empty() {
                components.push(component);
            }
        }

        if components.len() <= 1 {
            continue;
        }

        // Sort by size descending — largest keeps original ID
        components.sort_by(|a, b| b.len().cmp(&a.len()));

        for component in components.iter().skip(1) {
            let new_id = *next_id;
            *next_id += 1;

            let parent_vel = plates[pid].velocity;
            let ci = component[0] % nx;
            let cj = component[0] / nx;
            plates.push(Plate::new(
                new_id,
                plates[pid].plate_type,
                parent_vel,
                ci as f32,
                cj as f32,
            ));

            for &k in component {
                ids[k] = new_id;
            }

            info!(
                old_plate = pid,
                new_plate = new_id,
                cells = component.len(),
                "plate fragmented — continental breakup"
            );
        }

        // Reassign thin rift cells to nearest fragment by flood-fill from thick components
        let mut changed = true;
        while changed {
            changed = false;
            for &k in &thin_cells {
                if ids[k] != pid {
                    continue;
                }
                let ci = k % nx;
                let cj = k / nx;
                let neighbors = [
                    ((ci + 1) % nx, cj),
                    ((ci + nx - 1) % nx, cj),
                    (ci, (cj + 1) % ny),
                    (ci, (cj + ny - 1) % ny),
                ];
                for &(ni, nj) in &neighbors {
                    let nk = nj * nx + ni;
                    if ids[nk] != pid {
                        ids[k] = ids[nk];
                        changed = true;
                        break;
                    }
                }
            }
        }
    }
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

/// Build traction field with smooth interpolation from the displacement field.
///
/// Instead of assigning each cell the discrete velocity of its plate ID,
/// blend the traction based on the accumulated sub-pixel displacement.
/// As `disp` approaches ±0.5 (the shift threshold), the cell's traction
/// smoothly transitions toward the neighboring plate's velocity. The
/// transition is continuous — no discrete jump occurs when the ID shifts.
pub fn rebuild_traction_smooth(
    plate_ids: &[usize],
    plates: &[Plate],
    disp_x: &Field2D,
    disp_y: &Field2D,
    grid_width: usize,
    grid_height: usize,
) -> TractionField {
    let nx = grid_width;
    let ny = grid_height;
    let mut tx = Field2D::new(nx, ny);
    let mut ty = Field2D::new(nx, ny);

    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let pid = plate_ids[k];
            let plate = &plates[pid];

            let base_tx = plate.velocity.0 as f64;
            let base_ty = plate.velocity.1 as f64;

            let dx = disp_x.get(i, j);
            let dy = disp_y.get(i, j);

            // Weight: smoothly rises from 0 at disp=0 to 1 at |disp|=0.5
            let wx = smoothstep(dx.abs() * 2.0);
            let wy = smoothstep(dy.abs() * 2.0);

            let ni_x = if dx > 0.0 { (i + 1) % nx } else { (i + nx - 1) % nx };
            let ni_y = if dy > 0.0 { (j + 1) % ny } else { (j + ny - 1) % ny };

            let pid_nx = plate_ids[j * nx + ni_x];
            let pid_ny = plate_ids[ni_y * nx + i];

            let neighbor_tx_x = plates[pid_nx].velocity.0 as f64;
            let neighbor_ty_x = plates[pid_nx].velocity.1 as f64;
            let neighbor_tx_y = plates[pid_ny].velocity.0 as f64;
            let neighbor_ty_y = plates[pid_ny].velocity.1 as f64;

            let blended_tx = base_tx * (1.0 - wx) * (1.0 - wy)
                + neighbor_tx_x * wx * (1.0 - wy)
                + neighbor_tx_y * (1.0 - wx) * wy
                + neighbor_tx_x * wx * wy;
            let blended_ty = base_ty * (1.0 - wx) * (1.0 - wy)
                + neighbor_ty_x * wx * (1.0 - wy)
                + neighbor_ty_y * (1.0 - wx) * wy
                + neighbor_ty_y * wx * wy;

            tx.set(i, j, blended_tx);
            ty.set(i, j, blended_ty);
        }
    }

    TractionField { tx, ty }
}

/// Smooth hermite interpolation. Maps [0,1] → [0,1] with zero derivative at endpoints.
fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Rebuild the traction field from current plate_ids and plate velocities.
pub fn rebuild_traction(
    plate_ids: &[usize],
    plates: &[Plate],
    grid_width: usize,
    grid_height: usize,
) -> TractionField {
    let nx = grid_width;
    let ny = grid_height;
    let mut tx = Field2D::new(nx, ny);
    let mut ty = Field2D::new(nx, ny);

    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_ids[j * nx + i];
            let plate = &plates[pid];
            // Use plate.velocity (initial + slab pull), NOT mean_velocity (solved response).
            // The traction is the driving force from mantle convection, not the resulting flow.
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
    let nx = grid.nx();
    let ny = grid.ny();

    if !config.enabled || config.max_factor <= 1.0 {
        for k in 0..nx * ny {
            grid.eta_multiplier.data_mut()[k] = 1.0;
        }
        return;
    }

    let w = nx as f32;
    let h = ny as f32;

    // Step 1: Find maximum Voronoï radius for each continental plate.
    let mut max_radius = vec![0.0_f32; plates.len()];
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_ids[j * nx + i];
            let plate = &plates[pid];
            if plate.plate_type != PlateType::Continental || !plate.active {
                continue;
            }
            let dist = toroidal_distance_2d(
                i as f32,
                j as f32,
                plate.seed_x,
                plate.seed_y,
                w,
                h,
            );
            if dist > max_radius[pid] {
                max_radius[pid] = dist;
            }
        }
    }

    // Step 2: Compute multiplier per cell.
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_ids[j * nx + i];
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

            let dist = toroidal_distance_2d(
                i as f32,
                j as f32,
                plate.seed_x,
                plate.seed_y,
                w,
                h,
            );
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
        let d = toroidal_distance_sq_2d(1.0, 1.0, 2.0, 1.0, 10.0, 10.0);
        assert!((d - 1.0).abs() < 1e-6, "expected 1.0, got {d}");

        // x=1 and x=9 on a size=10 grid: direct=8, wrap=2, min=2 → d²=4
        let d = toroidal_distance_sq_2d(1.0, 0.0, 9.0, 0.0, 10.0, 10.0);
        assert!((d - 4.0).abs() < 1e-6, "expected 4.0, got {d}");
    }

    /// The 2D toroidal distance must use `period_x` for the x-difference
    /// and `period_y` for the y-difference. Swapping them goes unnoticed
    /// on square grids but produces the wrong distance — and so the
    /// wrong Voronoï assignment — on rectangular tori.
    #[test]
    fn test_toroidal_distance_rectangular() {
        // 10×4 grid. x-wrap period = 10, y-wrap period = 4.

        // x wraps: (9, 0) → (0, 0) is 1 cell via x-wrap (10 - 9 = 1).
        let d = toroidal_distance_2d(9.0, 0.0, 0.0, 0.0, 10.0, 4.0);
        assert!((d - 1.0).abs() < 1e-6, "x-wrap on 10×4: expected 1.0, got {d}");

        // y wraps: (0, 3) → (0, 0) is 1 cell via y-wrap (4 - 3 = 1).
        let d = toroidal_distance_2d(0.0, 3.0, 0.0, 0.0, 10.0, 4.0);
        assert!((d - 1.0).abs() < 1e-6, "y-wrap on 10×4: expected 1.0, got {d}");

        // No wrap needed in x (3 < 10/2): direct distance 3.
        let d = toroidal_distance_2d(3.0, 0.0, 0.0, 0.0, 10.0, 4.0);
        assert!((d - 3.0).abs() < 1e-6, "direct x on 10×4: expected 3.0, got {d}");

        // Axis-asymmetry: (5, 2) → (0, 0).
        //   - On a 10×4 torus: dx = min(5, 10-5) = 5, dy = min(2, 4-2) = 2 → √29
        //   - On a 4×10 torus: dx = min(5, 4-5) ... here 5 > 4/2 so dx = |4-5| = 1,
        //     dy = min(2, 10-2) = 2 → √5
        // A bug swapping period_x and period_y would return the same value
        // for both orientations; the correct implementation gives √29 vs √5.
        let d_wide = toroidal_distance_2d(5.0, 2.0, 0.0, 0.0, 10.0, 4.0);
        let d_tall = toroidal_distance_2d(5.0, 2.0, 0.0, 0.0, 4.0, 10.0);
        let expected_wide = (5.0_f32.powi(2) + 2.0_f32.powi(2)).sqrt();
        let expected_tall = (1.0_f32.powi(2) + 2.0_f32.powi(2)).sqrt();
        assert!(
            (d_wide - expected_wide).abs() < 1e-6,
            "10×4 (5,2)→(0,0): expected {expected_wide}, got {d_wide}"
        );
        assert!(
            (d_tall - expected_tall).abs() < 1e-6,
            "4×10 (5,2)→(0,0): expected {expected_tall}, got {d_tall}"
        );
        assert!(
            (d_wide - d_tall).abs() > 1e-3,
            "Distance should differ on 10×4 vs 4×10; a period swap would \
             make them identical. wide={d_wide}, tall={d_tall}"
        );
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

    /// Blur a delta function on a coprime rectangular grid and verify
    /// energy wraps correctly on both axes. Coprime dimensions ensure a
    /// horizontal-pass bug that uses `ny` for modulo instead of `nx`
    /// puts the wrapped energy at a different cell, which we can detect
    /// by sampling the expected wrap position.
    #[test]
    fn gaussian_blur_rectangular_wrap() {
        let nx = 13;
        let ny = 7;
        let mut field = GridF32::new(nx, ny, 0.0);
        // Delta at the right edge, middle row.
        field.set(nx - 1, 3, 1.0);

        let blurred = gaussian_blur(&field, 1.5);

        // Total mass is preserved (blur normalises the kernel).
        let total: f32 = blurred.data.iter().sum();
        assert!((total - 1.0).abs() < 1e-5, "blur mass drift: {total}");

        // The delta is at the RIGHT edge, so x-wrap carries energy to
        // column 0 in the SAME row. A mistake that wraps with `ny=7`
        // instead of `nx=13` would deposit energy at column (nx-1) % ny
        // = 12 % 7 = 5, not at column 0.
        let wrapped_x = blurred.data[3 * nx + 0]; // (0, 3) via x-wrap
        let bogus_wrap = blurred.data[3 * nx + 5]; // (5, 3) — wrap-with-ny bug
        assert!(
            wrapped_x > 1e-3,
            "x-wrap should deposit energy at (0, 3) but got {wrapped_x}"
        );
        assert!(
            wrapped_x > bogus_wrap * 5.0,
            "x-wrap landed at wrong cell: (0,3)={wrapped_x} vs (5,3)={bogus_wrap}. \
             Likely the blur used ny instead of nx for horizontal wrap."
        );

        // The delta is in the middle row (j=3), so the vertical pass
        // spreads energy into j=2 and j=4 in the same column.
        let up = blurred.data[2 * nx + (nx - 1)];
        let down = blurred.data[4 * nx + (nx - 1)];
        assert!(up > 1e-3 && down > 1e-3, "vertical spread missing: up={up}, down={down}");
    }

    /// Both PlateConfig input modes must agree on canonical storage.
    /// `with_resolution_aspect(w, r)` sets (w, round(w/r)); a manual
    /// `with_dimensions(w, h)` with the same h must produce an equivalent
    /// config (same grid_width, grid_height, and all other fields).
    #[test]
    fn plateconfig_input_modes_agree() {
        let r = 1.5_f32;
        let w = 128usize;
        let h = ((w as f32) / r).round() as usize;

        let a = PlateConfig::default().with_resolution_aspect(w, r);
        let b = PlateConfig::default().with_dimensions(w, h);

        assert_eq!(a.grid_width, b.grid_width);
        assert_eq!(a.grid_height, b.grid_height);
        assert_eq!(a.grid_width, 128);
        assert_eq!(a.grid_height, 85);

        // Non-dimensional fields also match (both came from default()).
        assert_eq!(a.num_plates, b.num_plates);
        assert_eq!(a.num_continental_plates, b.num_continental_plates);
    }

    #[test]
    fn test_plate_count_and_types() {
        let config = PlateConfig {
            num_plates: 10,
            num_continental_plates: 3,
            grid_width: 64,
            grid_height: 64,
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
        let config = PlateConfig { grid_width: 32, grid_height: 32, ..Default::default() };
        let seed = WorldSeed::new(99);
        let result = generate_plates(&config, &seed);

        for &id in &result.plate_ids {
            assert!(id < config.num_plates, "cell has invalid plate id {id}");
        }
    }

    #[test]
    fn test_thickness_values() {
        let config = PlateConfig {
            grid_width: 64,
            grid_height: 64,
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
            grid_width: 32,
            grid_height: 32,
            boundary_smoothing_sigma: 0.0,
            ..Default::default()
        };
        let seed = WorldSeed::new(7);
        let result = generate_plates(&config, &seed);

        // At least some cell at x=0 and x=31 should exist — just verify no panic
        let _id_left = result.plate_ids[0];
        let _id_right = result.plate_ids[31];
        assert_eq!(result.grid_width, 32);
        assert_eq!(result.grid_height, 32);
    }

    // ── Dynamic boundary tests ───────────────────────────────────────────

    #[test]
    fn seeds_move_with_velocity() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);

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
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: 0.0,
            centroid_y: 0.0,
        }];

        advect_seeds(&mut plates, &grid, 0.1);

        assert!(plates[0].seed_x > 10.0, "Seed should move right: {}", plates[0].seed_x);
        assert!((plates[0].seed_y - 16.0).abs() < 1e-3, "Seed should not move vertically");
    }

    #[test]
    fn seed_wraps_periodically() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);

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
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: 0.0,
            centroid_y: 0.0,
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
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            },
            Plate {
                id: 1,
                plate_type: PlateType::Oceanic,
                velocity: (0.0, 0.0),
                seed_x: 8.0, // same as plate 0
                seed_y: 8.0,
                active: true,
                subducted_mass: 0.0,
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            },
        ];

        recompute_voronoi(&mut plate_ids, &plates, n, n);
        let disappeared = detect_disappeared_plates(&plate_ids, &mut plates);

        // When two seeds are at the same position, the one with lower ID wins
        assert!(disappeared.contains(&1), "Plate 1 should disappear when seeds overlap");
        assert!(!plates[1].active);
    }

    #[test]
    fn recompute_voronoi_matches_generate_plates() {
        let config = PlateConfig {
            num_plates: 4,
            grid_width: 32,
            grid_height: 32,
            boundary_smoothing_sigma: 0.0,
            ..Default::default()
        };
        let seed = WorldSeed::new(42);
        let result = generate_plates(&config, &seed);

        // Recompute from the same seeds
        let mut plate_ids = vec![0usize; 32 * 32];
        recompute_voronoi(&mut plate_ids, &result.plates, 32, 32);

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
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            },
            Plate {
                id: 1,
                plate_type: PlateType::Oceanic,
                velocity: (0.0, 0.0),
                seed_x: 12.0,
                seed_y: 8.0,
                active: false,
                subducted_mass: 0.0, // inactive
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            },
        ];

        let mut plate_ids = vec![0usize; n * n];
        recompute_voronoi(&mut plate_ids, &plates, n, n);

        // All cells should belong to plate 0 since plate 1 is inactive
        assert!(plate_ids.iter().all(|&id| id == 0));
    }

    // ── Cratonic rigidity tests ──────────────────────────────────────

    #[test]
    fn cratonic_multiplier_is_highest_at_seed() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);

        let plate_ids: Vec<usize> = vec![0; n * n];
        let plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Continental,
            velocity: (0.0, 0.0),
            seed_x: 16.0,
            seed_y: 16.0,
            active: true,
            subducted_mass: 0.0,
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: 0.0,
            centroid_y: 0.0,
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
        let mut grid = StaggeredGrid::new(n, n, dx);

        let plate_ids: Vec<usize> = vec![0; n * n];
        let plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Oceanic,
            velocity: (0.0, 0.0),
            seed_x: 8.0,
            seed_y: 8.0,
            active: true,
            subducted_mass: 0.0,
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: 0.0,
            centroid_y: 0.0,
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
        let mut grid = StaggeredGrid::new(n, n, dx);

        let plate_ids: Vec<usize> = vec![0; n * n];
        let plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Continental,
            velocity: (0.0, 0.0),
            seed_x: 8.0,
            seed_y: 8.0,
            active: true,
            subducted_mass: 0.0,
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: 0.0,
            centroid_y: 0.0,
        }];

        let config = CratonicConfig { enabled: false, ..Default::default() };
        compute_viscosity_multiplier(&mut grid, &plate_ids, &plates, &config);

        for k in 0..n * n {
            assert!((grid.eta_multiplier.data()[k] - 1.0).abs() < 1e-10);
        }
    }
}
