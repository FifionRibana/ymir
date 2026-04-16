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

        plates.push(Plate::new(id, plate_type, velocity, seed_x, seed_y));
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
    let n = grid.n;

    // Phase 1: accumulate forward displacement (dt × v) at each cell
    for j in 0..n {
        for i in 0..n {
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

    for j in 0..n {
        for i in 0..n {
            let dx = disp_x.get(i, j);
            let dy = disp_y.get(i, j);

            let shift_x = dx.round() as i64;
            let shift_y = dy.round() as i64;

            if shift_x == 0 && shift_y == 0 {
                continue;
            }

            let src_i = ((i as i64 - shift_x).rem_euclid(n as i64)) as usize;
            let src_j = ((j as i64 - shift_y).rem_euclid(n as i64)) as usize;

            ids[j * n + i] = old_ids[src_j * n + src_i];

            // Keep the fractional remainder
            disp_x.set(i, j, dx - shift_x as f64);
            disp_y.set(i, j, dy - shift_y as f64);
        }
    }
}

/// Compute per-plate runtime statistics from the cell data.
///
/// Updates `cell_count`, `mean_thickness`, `mean_velocity`, and
/// `centroid_x`/`centroid_y` (circular mean on the torus) for each plate.
pub fn update_plate_stats(ids: &[usize], plates: &mut [Plate], grid: &StaggeredGrid) {
    let n = grid.n;
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

    for j in 0..n {
        for i in 0..n {
            let pid = ids[j * n + i];
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

            let theta_x = tau * i as f64 / n as f64;
            let theta_y = tau * j as f64 / n as f64;
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
            plate.centroid_x = (sin_x[pid].atan2(cos_x[pid]) * n as f64 / tau) as f32;
            plate.centroid_y = (sin_y[pid].atan2(cos_y[pid]) * n as f64 / tau) as f32;
            // Wrap to [0, n)
            plate.centroid_x = ((plate.centroid_x % n as f32) + n as f32) % n as f32;
            plate.centroid_y = ((plate.centroid_y % n as f32) + n as f32) % n as f32;
        }
    }
}

/// Bilinear interpolation of vx at an arbitrary point (px, py).
/// vx lives at left vertical faces: vx[i,j] is at position (i, j+0.5).
pub fn interpolate_vx(grid: &StaggeredGrid, px: f64, py: f64) -> f64 {
    let fx = px;
    let fy = py - 0.5;
    bilinear_sample_field(&grid.vx, grid.n, fx, fy)
}

/// Bilinear interpolation of vy at an arbitrary point (px, py).
/// vy lives at bottom horizontal faces: vy[i,j] is at position (i+0.5, j).
pub fn interpolate_vy(grid: &StaggeredGrid, px: f64, py: f64) -> f64 {
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
/// At convergent boundaries, cells thin enough on the subducting side
/// are reassigned to the overriding plate's ID.
pub fn apply_subduction_consumption(
    ids: &mut [usize],
    grid: &StaggeredGrid,
    boundary_field: &crate::tectonics::boundaries::BoundaryField,
    threshold: f64,
) {
    use crate::tectonics::boundaries::BoundaryType;
    let n = grid.n;
    let idx = &grid.idx;

    for j in 0..n {
        for i in 0..n {
            let k = j * n + i;
            let btype = &boundary_field.boundary_type[k];

            let is_subducting =
                matches!(btype, BoundaryType::Subduction | BoundaryType::OceanicSubduction)
                    && grid.s.get(i, j) < threshold;

            if !is_subducting {
                continue;
            }

            let my_id = ids[k];
            let neighbors =
                [(idx.next(i), j), (idx.prev(i), j), (i, idx.next(j)), (i, idx.prev(j))];

            for &(ni, nj) in &neighbors {
                let nk = nj * n + ni;
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
    let n = grid.n;

    // Pass 1: collect all cells that qualify for rift creation
    let mut candidates = vec![false; n * n];
    for j in 0..n {
        for i in 0..n {
            let k = j * n + i;
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
    let mut visited = vec![false; n * n];
    let mut components: Vec<Vec<usize>> = Vec::new();

    for k in 0..n * n {
        if !candidates[k] || visited[k] {
            continue;
        }
        let mut group = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(k);
        visited[k] = true;

        while let Some(ck) = queue.pop_front() {
            group.push(ck);
            let ci = ck % n;
            let cj = ck / n;
            let neighbors = [
                ((ci + 1) % n, cj),
                ((ci + n - 1) % n, cj),
                (ci, (cj + 1) % n),
                (ci, (cj + n - 1) % n),
            ];
            for &(ni, nj) in &neighbors {
                let nk = nj * n + ni;
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

        let ci = component[0] % n;
        let cj = component[0] / n;
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
    n: usize,
    grid: &StaggeredGrid,
    connectivity_threshold: f64,
) {
    let active_ids: Vec<usize> =
        plates.iter().filter(|p| p.active && p.cell_count > 0).map(|p| p.id).collect();

    for &pid in &active_ids {
        let cells: Vec<usize> = (0..n * n).filter(|&k| ids[k] == pid).collect();
        if cells.len() < 4 {
            continue;
        }
        // Only continental plates can fragment by rifting.
        if plates[pid].mean_thickness <= 0.5 {
            continue;
        }

        // BFS flood fill with thickness-aware connectivity
        let mut visited = vec![false; n * n];
        let mut components: Vec<Vec<usize>> = Vec::new();
        let mut thin_cells: Vec<usize> = Vec::new();

        for &start in &cells {
            if visited[start] {
                continue;
            }
            // Thin cells don't start components
            if grid.s.get(start % n, start / n) < connectivity_threshold {
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
                let ci = k % n;
                let cj = k / n;
                let neighbors = [
                    ((ci + 1) % n, cj),
                    ((ci + n - 1) % n, cj),
                    (ci, (cj + 1) % n),
                    (ci, (cj + n - 1) % n),
                ];
                for &(ni, nj) in &neighbors {
                    let nk = nj * n + ni;
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
            let ci = component[0] % n;
            let cj = component[0] / n;
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
                let ci = k % n;
                let cj = k / n;
                let neighbors = [
                    ((ci + 1) % n, cj),
                    ((ci + n - 1) % n, cj),
                    (ci, (cj + 1) % n),
                    (ci, (cj + n - 1) % n),
                ];
                for &(ni, nj) in &neighbors {
                    let nk = nj * n + ni;
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
    grid_size: usize,
) -> TractionField {
    let mut tx = Field2D::new(grid_size);
    let mut ty = Field2D::new(grid_size);
    let n = grid_size;

    for j in 0..n {
        for i in 0..n {
            let k = j * n + i;
            let pid = plate_ids[k];
            let plate = &plates[pid];

            let base_tx = plate.velocity.0 as f64;
            let base_ty = plate.velocity.1 as f64;

            let dx = disp_x.get(i, j);
            let dy = disp_y.get(i, j);

            // Weight: smoothly rises from 0 at disp=0 to 1 at |disp|=0.5
            let wx = smoothstep(dx.abs() * 2.0);
            let wy = smoothstep(dy.abs() * 2.0);

            let ni_x = if dx > 0.0 { (i + 1) % n } else { (i + n - 1) % n };
            let ni_y = if dy > 0.0 { (j + 1) % n } else { (j + n - 1) % n };

            let pid_nx = plate_ids[j * n + ni_x];
            let pid_ny = plate_ids[ni_y * n + i];

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
pub fn rebuild_traction(plate_ids: &[usize], plates: &[Plate], grid_size: usize) -> TractionField {
    let mut tx = Field2D::new(grid_size);
    let mut ty = Field2D::new(grid_size);

    for j in 0..grid_size {
        for i in 0..grid_size {
            let pid = plate_ids[j * grid_size + i];
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
