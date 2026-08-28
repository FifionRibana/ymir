//! Step 11 Phase 2 sanity check — visual ASCII dump of the velocity
//! field built for a 2-plate scenario with opposing velocities.
//!
//! Not a regression test (no asserts). Run via:
//!
//! ```text
//! cargo test --release -p ymir-core --lib \
//!     tectonics_v2::plate_kinematic::sanity_check::two_plates_opposing \
//!     -- --ignored --nocapture
//! ```
//!
//! This is the human-readable reviewer checkpoint requested by the
//! Step 11 implementation prompt — ensures the BFS + smoothstep
//! produces visibly converging velocity arrows around the inter-plate
//! boundary, no discontinuities. Once Phase 5 lands and the viz
//! panel can render arrows from `PlateKinematicConfig`, this dump
//! becomes redundant — the panel preview is the better artefact.

#[cfg(test)]
mod tests {
    use super::super::field::build;
    use crate::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

    /// 32×32 grid, 2 plates, opposing horizontal velocities. Dumps
    /// `vx` as ASCII (one char per cell, sign + 2-digit fractional)
    /// plus a `|v|` magnitude scan to confirm the field magnitude
    /// stays bounded by the assigned magnitude.
    #[test]
    #[ignore]
    fn two_plates_opposing() {
        let nx = 32;
        let ny = 32;
        let cfg = VoronoiConfig { num_plates: 2, continental_ratio: 0.5 };
        let plates = generate_voronoi(nx, ny, &cfg, 42);
        let velocities = vec![(0.5, 0.0), (-0.5, 0.0)];
        let width = 1.5;
        let (vx, vy) = build(nx, ny, &plates.plate_id, &velocities, width);

        println!("\n=== plate_id ===");
        for j in 0..ny {
            for i in 0..nx {
                print!("{}", plates.plate_id.get(i, j));
            }
            println!();
        }

        println!("\n=== vx (sign: '+' = +0.5, '-' = -0.5, '.' = 0, others = blend) ===");
        for j in 0..ny {
            for i in 0..nx {
                let v = vx[j * nx + i];
                let ch = if (v - 0.5).abs() < 1e-9 {
                    '+'
                } else if (v + 0.5).abs() < 1e-9 {
                    '-'
                } else if v.abs() < 1e-9 {
                    '.'
                } else if v > 0.25 {
                    '>'
                } else if v < -0.25 {
                    '<'
                } else if v > 0.0 {
                    'r'
                } else {
                    'l'
                };
                print!("{}", ch);
            }
            println!();
        }

        println!("\n=== vy ===");
        for j in 0..ny {
            for i in 0..nx {
                let v = vy[j * nx + i];
                let ch = if v.abs() < 1e-9 { '.' } else { '?' };
                print!("{}", ch);
            }
            println!();
        }

        let max_mag =
            vx.iter().zip(vy.iter()).map(|(&a, &b)| (a * a + b * b).sqrt()).fold(0.0_f64, f64::max);
        let mean_mag = {
            let total: f64 = vx.iter().zip(vy.iter()).map(|(&a, &b)| (a * a + b * b).sqrt()).sum();
            total / (nx * ny) as f64
        };
        println!("\nmax|v| = {:.4} (assigned magnitude = 0.5000)", max_mag);
        println!("mean|v| = {:.4}", mean_mag);

        // Boundary smoothness check — count cells whose vx differs
        // from both ±0.5 and 0 (i.e. they're inside the transition).
        let smoothing_cells: usize = vx
            .iter()
            .filter(|&&v| (v - 0.5).abs() > 1e-9 && (v + 0.5).abs() > 1e-9 && v.abs() > 1e-9)
            .count();
        let total_cells = nx * ny;
        println!(
            "smoothing band cells: {} / {} ({:.1}%)",
            smoothing_cells,
            total_cells,
            100.0 * smoothing_cells as f64 / total_cells as f64
        );
    }
}
