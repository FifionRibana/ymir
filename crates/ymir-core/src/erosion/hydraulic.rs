//! Beyer 2015 particle-based hydraulic erosion.
//!
//! Millions of simulated water droplets carve valleys, form drainage networks,
//! and deposit sediment on a heightmap.

use rand::Rng;
use serde::{Deserialize, Serialize};
use tracing::{debug, info, warn};

use crate::grid::GridF32;
use crate::seed::WorldSeed;

/// Configuration for hydraulic erosion.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ErosionConfig {
    pub num_droplets: usize,
    pub deposition_rate: f32,
    pub erosion_rate: f32,
    pub inertia: f32,
    pub gravity: f32,
    pub evaporation_rate: f32,
    pub max_lifetime: usize,
    pub min_slope: f32,
    pub erosion_radius: usize,
    pub coastal_deposition_range: usize,
    pub sea_level: f32,
    pub batch_size: usize,
    /// Reference grid size for parameter calibration. The algorithm scales
    /// step size automatically for larger or smaller grids. Default: 256.
    pub reference_size: usize,
}

impl Default for ErosionConfig {
    fn default() -> Self {
        Self {
            num_droplets: 5_000_000,
            deposition_rate: 0.35,
            erosion_rate: 0.4,
            inertia: 0.08,
            gravity: 6.0,
            evaporation_rate: 0.015,
            max_lifetime: 150,
            min_slope: 0.001,
            erosion_radius: 3,
            coastal_deposition_range: 12,
            sea_level: 0.1,
            batch_size: 50_000,
            reference_size: 256,
        }
    }
}

/// Result of the erosion simulation.
pub struct ErosionResult {
    pub heightmap: GridF32,
    pub sediment: GridF32,
    pub stats: ErosionStats,
}

/// Statistics collected during erosion.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ErosionStats {
    pub total_droplets: usize,
    pub total_eroded: f64,
    pub total_deposited: f64,
    pub avg_lifetime: f32,
}

/// Precomputed erosion/deposition brush — a weighted disk of pixel offsets.
struct ErosionBrush {
    offsets: Vec<(i32, i32)>,
    weights: Vec<f32>,
}

impl ErosionBrush {
    fn new(radius: usize) -> Self {
        let r = radius as i32;
        let rf = radius as f32;
        let mut offsets = Vec::new();
        let mut weights = Vec::new();

        for dy in -r..=r {
            for dx in -r..=r {
                let dist = ((dx * dx + dy * dy) as f32).sqrt();
                if dist <= rf {
                    let w = (rf - dist).max(0.0);
                    offsets.push((dx, dy));
                    weights.push(w);
                }
            }
        }

        let sum: f32 = weights.iter().sum();
        if sum > 0.0 {
            for w in &mut weights {
                *w /= sum;
            }
        }

        Self { offsets, weights }
    }
}

/// Run the hydraulic erosion simulation.
///
/// `progress_callback(completed, total, heightmap)` is called after each batch.
/// Return `false` from it to cancel early.
pub fn run_erosion(
    heightmap: &GridF32,
    config: &ErosionConfig,
    seed: &WorldSeed,
    mut progress_callback: impl FnMut(usize, usize, &GridF32) -> bool,
) -> ErosionResult {
    let w = heightmap.width;
    let h = heightmap.height;
    let mut hmap = heightmap.clone();
    let mut sediment_map = GridF32::new(w, h, 0.0);
    let brush = ErosionBrush::new(config.erosion_radius);

    let total = config.num_droplets;
    let num_batches = (total + config.batch_size - 1) / config.batch_size;

    let mut stats = ErosionStats::default();
    let mut total_lifetime: u64 = 0;

    let step_size = (w.max(h) as f32 / config.reference_size as f32).max(1.0);

    info!(
        "Erosion: {} droplets on {}x{} grid ({} batches, step_size={:.1})",
        total, w, h, num_batches, step_size
    );

    for batch in 0..num_batches {
        let batch_start = batch * config.batch_size;
        let batch_end = (batch_start + config.batch_size).min(total);
        let batch_count = batch_end - batch_start;

        let mut rng = seed.rng_for_indexed("erosion_batch", batch as u64);

        for _ in 0..batch_count {
            let lifetime = simulate_droplet(
                &mut hmap,
                &mut sediment_map,
                &brush,
                config,
                step_size,
                &mut rng,
                &mut stats,
            );
            total_lifetime += lifetime as u64;
        }

        stats.total_droplets = batch_end;

        if batch % 10 == 0 {
            debug!(
                "Erosion batch {}/{}: {:.0}% done",
                batch,
                num_batches,
                batch_end as f64 / total as f64 * 100.0
            );
        }

        if !progress_callback(batch_end, total, &hmap) {
            info!("Erosion cancelled at batch {}/{}", batch, num_batches);
            break;
        }
    }

    stats.avg_lifetime = if stats.total_droplets > 0 {
        total_lifetime as f32 / stats.total_droplets as f32
    } else {
        0.0
    };

    if stats.avg_lifetime < 5.0 && stats.total_droplets > 0 {
        warn!(
            "Very low average droplet lifetime ({:.1}). Parameters may need tuning.",
            stats.avg_lifetime
        );
    }

    info!(
        "Erosion complete: {} droplets, eroded={:.2}, deposited={:.2}, avg_life={:.1}",
        stats.total_droplets, stats.total_eroded, stats.total_deposited, stats.avg_lifetime
    );

    ErosionResult { heightmap: hmap, sediment: sediment_map, stats }
}

fn simulate_droplet(
    hmap: &mut GridF32,
    sediment_map: &mut GridF32,
    brush: &ErosionBrush,
    config: &ErosionConfig,
    step_size: f32,
    rng: &mut impl Rng,
    stats: &mut ErosionStats,
) -> usize {
    let w = hmap.width as f32;
    let h = hmap.height as f32;

    // Spawn only on land — rejection sampling
    let mut x: f32 = 0.0;
    let mut y: f32 = 0.0;
    let mut found_land = false;
    for _ in 0..100 {
        x = rng.random::<f32>() * (w - 1.0);
        y = rng.random::<f32>() * (h - 1.0);
        if hmap.sample_bilinear(x, y) > config.sea_level {
            found_land = true;
            break;
        }
    }
    if !found_land {
        return 0;
    }
    let mut dir_x: f32 = 0.0;
    let mut dir_y: f32 = 0.0;
    let mut speed: f32 = 1.0;
    let mut water: f32 = 1.0;
    let mut sediment: f32 = 0.0;
    let mut below_sea_steps: usize = 0;

    for step in 0..config.max_lifetime {
        // i. Gradient
        let (gx, gy) = hmap.gradient_at_f(x, y, step_size);

        // ii. Update direction with inertia
        dir_x = dir_x * config.inertia - gx * (1.0 - config.inertia);
        dir_y = dir_y * config.inertia - gy * (1.0 - config.inertia);
        let len = (dir_x * dir_x + dir_y * dir_y).sqrt();
        if len < 1e-10 {
            let angle: f32 = rng.random::<f32>() * std::f32::consts::TAU;
            dir_x = angle.cos();
            dir_y = angle.sin();
        } else {
            dir_x /= len;
            dir_y /= len;
        }

        // iii. Move
        let new_x = x + dir_x * step_size;
        let new_y = y + dir_y * step_size;
        if new_x < 0.0 || new_x >= w - 1.0 || new_y < 0.0 || new_y >= h - 1.0 {
            return step;
        }

        // iv. Height difference
        let h_old = hmap.sample_bilinear(x, y);
        let h_new = hmap.sample_bilinear(new_x, new_y);
        let delta_h = h_new - h_old;

        // v. Deposit or erode
        if delta_h > 0.0 {
            // Moving uphill — deposit
            let amount = sediment.min(delta_h);
            apply_brush(hmap, sediment_map, brush, x, y, amount);
            stats.total_deposited += amount as f64;
            sediment -= amount;
            if sediment < delta_h * 0.5 {
                return step;
            }
        } else {
            // Downhill or flat
            let capacity = (-delta_h).max(config.min_slope) * speed * water;
            if sediment > capacity {
                // Over capacity — deposit excess
                let amount = (sediment - capacity) * config.deposition_rate;
                apply_brush(hmap, sediment_map, brush, x, y, amount);
                stats.total_deposited += amount as f64;
                sediment -= amount;
            } else {
                // Under capacity — erode
                let amount = ((capacity - sediment) * config.erosion_rate).min(-delta_h);
                apply_brush(hmap, sediment_map, brush, x, y, -amount);
                stats.total_eroded += amount as f64;
                sediment += amount;
            }
        }

        // vii. Update velocity
        speed = (speed * speed - delta_h * config.gravity).max(0.0).sqrt();

        // viii. Evaporate
        water *= 1.0 - config.evaporation_rate;

        // ix. Coastal check
        if hmap.sample_bilinear(new_x, new_y) <= config.sea_level {
            below_sea_steps += 1;
            if below_sea_steps > config.coastal_deposition_range {
                if sediment > 0.0 {
                    apply_brush(hmap, sediment_map, brush, x, y, sediment);
                    stats.total_deposited += sediment as f64;
                }
                return step;
            }
        } else {
            below_sea_steps = 0;
        }

        // x. Advance
        x = new_x;
        y = new_y;
    }

    config.max_lifetime
}

/// Apply a brush at position (px, py). Positive amount = deposit, negative = erode.
fn apply_brush(
    hmap: &mut GridF32,
    sediment_map: &mut GridF32,
    brush: &ErosionBrush,
    px: f32,
    py: f32,
    amount: f32,
) {
    let ix = px.floor() as i32;
    let iy = py.floor() as i32;

    for (k, &(dx, dy)) in brush.offsets.iter().enumerate() {
        let bx = ix + dx;
        let by = iy + dy;
        let delta = amount * brush.weights[k];

        let old = hmap.get(bx, by);
        hmap.set(bx as usize, by as usize, old + delta);

        if delta > 0.0 {
            let old_sed = sediment_map.get(bx, by);
            sediment_map.set(bx as usize, by as usize, old_sed + delta);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brush_weights_sum_to_one() {
        let brush = ErosionBrush::new(3);
        let sum: f32 = brush.weights.iter().sum();
        assert!((sum - 1.0).abs() < 1e-5, "Brush weights should sum to 1.0, got {sum}");
        assert!(!brush.offsets.is_empty());
    }

    #[test]
    fn erosion_is_deterministic() {
        let mut hmap = GridF32::new(64, 64, 0.0);
        for j in 0..64 {
            for i in 0..64 {
                let cx = i as f32 - 32.0;
                let cy = j as f32 - 32.0;
                let dist = (cx * cx + cy * cy).sqrt();
                hmap.set(i, j, (1.0 - dist / 40.0).clamp(0.0, 1.0));
            }
        }

        let seed = WorldSeed::new(42);
        let config =
            ErosionConfig { num_droplets: 10_000, batch_size: 5_000, ..Default::default() };

        let r1 = run_erosion(&hmap, &config, &seed, |_, _, _| true);
        let r2 = run_erosion(&hmap, &config, &seed, |_, _, _| true);
        assert_eq!(r1.heightmap.data, r2.heightmap.data, "Erosion must be deterministic");
    }

    #[test]
    fn cone_erosion_carves_valleys() {
        let n = 128;
        let mut hmap = GridF32::new(n, n, 0.0);
        let center = n as f32 / 2.0;
        let max_dist = center * 0.8;
        for j in 0..n {
            for i in 0..n {
                let dx = i as f32 - center;
                let dy = j as f32 - center;
                let dist = (dx * dx + dy * dy).sqrt();
                hmap.set(i, j, (1.0 - dist / max_dist).clamp(0.0, 1.0));
            }
        }

        let original_mean = hmap.mean();
        let seed = WorldSeed::new(42);
        let config = ErosionConfig {
            num_droplets: 100_000,
            batch_size: 10_000,
            sea_level: 0.0,
            ..Default::default()
        };

        let result = run_erosion(&hmap, &config, &seed, |_, _, _| true);

        assert!(result.heightmap.mean() < original_mean, "Erosion should lower mean altitude");

        for &v in &result.heightmap.data {
            assert!(v >= -0.1 && v <= 1.1, "Height out of reasonable range: {v}");
        }

        assert!(result.stats.total_eroded > 0.0);
        assert!(result.stats.total_deposited > 0.0);
    }

    #[test]
    fn flat_terrain_minimal_change() {
        let hmap = GridF32::new(64, 64, 0.5);
        let seed = WorldSeed::new(42);
        let config = ErosionConfig {
            num_droplets: 10_000,
            batch_size: 5_000,
            sea_level: 0.0,
            ..Default::default()
        };

        let result = run_erosion(&hmap, &config, &seed, |_, _, _| true);
        let max_diff =
            result.heightmap.data.iter().map(|&v| (v - 0.5).abs()).fold(0.0f32, f32::max);

        assert!(max_diff < 0.05, "Flat terrain should change minimally, max_diff={max_diff}");
    }

    #[test]
    fn coastal_deposition_near_sea_level() {
        let n = 64;
        let mut hmap = GridF32::new(n, n, 0.0);
        for j in 0..n {
            for i in 0..n {
                let t = i as f32 / n as f32;
                hmap.set(i, j, 0.6 - t * 0.5);
            }
        }

        let seed = WorldSeed::new(42);
        let config = ErosionConfig {
            num_droplets: 50_000,
            batch_size: 10_000,
            sea_level: 0.3,
            coastal_deposition_range: 12,
            ..Default::default()
        };

        let result = run_erosion(&hmap, &config, &seed, |_, _, _| true);
        let max_sed = result.sediment.data.iter().cloned().fold(0.0f32, f32::max);
        assert!(max_sed > 0.0, "Should have coastal sediment deposits");
    }
}
