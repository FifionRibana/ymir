//! R7 boundary displacement — applies Perlin / Simplex noise
//! displacement to the per-cell sampling position before re-
//! querying the nearest Voronoï seed. Phase 2 Track B
//! sub-component 1 (Issue #131).
//!
//! See [`super`] for the module-level rationale.

use noise::{Fbm, MultiFractal, NoiseFn, Perlin};

use crate::tectonics_v2::voronoi::PlateIdField;

use super::params::R7InitParams;

/// Magic XOR constant for deriving the y-channel seed from the
/// x-channel seed. Borrowed in spirit from the FBM seed-derivation
/// pattern in `tectonics_v2::init::radial_profile_fbm`: any
/// non-zero magic suffices for reasonable channel independence;
/// the value below is arbitrary but stable.
const DY_CHANNEL_XOR_MAGIC: u64 = 0xDE13_4DDA_5EED_u64;

/// Pack a `u64` seed into a `u32` deterministically by XOR-folding
/// the high and low halves. Preserves entropy across the full
/// `u64` range while matching the `noise::Fbm<Perlin>::new(u32)`
/// constructor signature.
fn u64_to_u32_seed(seed: u64) -> u32 {
    ((seed >> 32) as u32) ^ (seed as u32)
}

/// Compute the squared Euclidean distance between two points on
/// the periodic domain `[0, nx) × [0, ny)` using the minimum-
/// image convention. Match-request-scope: reimplemented here
/// rather than promoted from the private `periodic_dist_sq` in
/// `tectonics_v2::voronoi` to avoid scope creep on a one-file
/// helper.
#[inline]
fn periodic_dist_sq(ax: f64, ay: f64, bx: f64, by: f64, nx: f64, ny: f64) -> f64 {
    let mut dx = (ax - bx).abs();
    let mut dy = (ay - by).abs();
    if dx > 0.5 * nx {
        dx = nx - dx;
    }
    if dy > 0.5 * ny {
        dy = ny - dy;
    }
    dx * dx + dy * dy
}

/// Apply Perlin / Simplex noise displacement to plate boundaries.
///
/// For each grid cell `(i, j)`:
///
/// 1. Sample two independent FBM noise channels at the cell's
///    normalised position to obtain a displacement vector
///    `(dx, dy)` (each component approximately in `[-1, +1]`).
/// 2. Scale by `params.amplitude` to obtain displacement in cell
///    units.
/// 3. Compute the displaced query position
///    `(i + 0.5 + amp · dx, j + 0.5 + amp · dy)` (cell-centre +
///    displacement).
/// 4. Find the nearest Voronoï seed to that displaced query
///    position under toroidal minimum-image distance.
/// 5. Reassign the cell's `plate_id` to that seed's id.
///
/// Cells whose displaced query position lands in the same Voronoï
/// region as their original cell-centre query are unchanged.
/// Cells near a boundary whose displacement carries them across
/// into a neighbouring region are reassigned — producing curved
/// boundaries.
///
/// `params.enabled = false` makes this a no-op (W4 closure-
/// isolation).
///
/// ## Determinism
///
/// Same `(plate_id_field, seed_coords, params)` → bit-identical
/// output `plate_id_field`. Two independent FBM channels are
/// derived from `params.seed` via a private `u64 → u32` XOR-fold
/// helper (x-channel uses `params.seed` directly; y-channel uses
/// `params.seed ^ DY_CHANNEL_XOR_MAGIC`).
///
/// ## Inputs
///
/// - `plate_id_field`: mutable `PlateIdField`, updated in place.
/// - `seed_coords`: per-plate Voronoï centre coordinates in cell
///   units, indexed by plate id. Typically
///   `VoronoiPlates::seed_coords` from
///   [`crate::tectonics_v2::voronoi::generate_voronoi`].
/// - `params`: R7 displacement tunables.
pub fn apply_boundary_displacement(
    plate_id_field: &mut PlateIdField,
    seed_coords: &[(f64, f64)],
    params: &R7InitParams,
) {
    if !params.enabled {
        return;
    }
    if seed_coords.is_empty() {
        return;
    }

    let nx = plate_id_field.nx();
    let ny = plate_id_field.ny();
    let nx_f = nx as f64;
    let ny_f = ny as f64;

    let seed_x = u64_to_u32_seed(params.seed);
    let seed_y = u64_to_u32_seed(params.seed ^ DY_CHANNEL_XOR_MAGIC);

    let fbm_x = Fbm::<Perlin>::new(seed_x)
        .set_octaves(params.octaves as usize)
        .set_persistence(params.persistence)
        .set_frequency(1.0);
    let fbm_y = Fbm::<Perlin>::new(seed_y)
        .set_octaves(params.octaves as usize)
        .set_persistence(params.persistence)
        .set_frequency(1.0);

    let freq = params.frequency;

    for j in 0..ny {
        for i in 0..nx {
            // Sample noise at the cell's normalised position
            // scaled by `frequency`. `(i + 0.5) / nx` puts the
            // sample at cell centre in `[0, 1)`; the `× frequency`
            // controls how many wavelengths fit in the domain.
            let nx_norm = (i as f64 + 0.5) / nx_f * freq;
            let ny_norm = (j as f64 + 0.5) / ny_f * freq;
            let dx_noise = fbm_x.get([nx_norm, ny_norm]);
            let dy_noise = fbm_y.get([nx_norm, ny_norm]);

            // Displaced query position in cell units.
            let qx = i as f64 + 0.5 + params.amplitude * dx_noise;
            let qy = j as f64 + 0.5 + params.amplitude * dy_noise;

            // Find nearest seed under periodic distance. Naive
            // O(N_seeds) — at 8 plates × 4096 cells = 32 768
            // distance computations, sub-millisecond at 64².
            let mut best_id: u16 = 0;
            let mut best_d2 = f64::INFINITY;
            for (sid, &(sx, sy)) in seed_coords.iter().enumerate() {
                let d2 = periodic_dist_sq(qx, qy, sx, sy, nx_f, ny_f);
                if d2 < best_d2 {
                    best_d2 = d2;
                    best_id = sid as u16;
                }
            }
            plate_id_field.set(i, j, best_id);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fixture: 2 plates split at column nx/2 via a clean Voronoï
    /// with seeds at the centres of each half. Produces a vertical
    /// boundary at `i = nx / 2`.
    fn two_plate_voronoi(nx: usize, ny: usize) -> (PlateIdField, Vec<(f64, f64)>) {
        let mut plate_id = PlateIdField::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                plate_id.set(i, j, if i < nx / 2 { 0 } else { 1 });
            }
        }
        let seed_coords = vec![
            (nx as f64 * 0.25, ny as f64 * 0.5), // plate 0 centre
            (nx as f64 * 0.75, ny as f64 * 0.5), // plate 1 centre
        ];
        (plate_id, seed_coords)
    }

    /// Test 1 — deterministic given seed.
    #[test]
    fn r7_boundary_displacement_deterministic() {
        let (mut p_a, seeds_a) = two_plate_voronoi(32, 32);
        let (mut p_b, seeds_b) = two_plate_voronoi(32, 32);
        let params = R7InitParams::default();

        apply_boundary_displacement(&mut p_a, &seeds_a, &params);
        apply_boundary_displacement(&mut p_b, &seeds_b, &params);

        for j in 0..32 {
            for i in 0..32 {
                assert_eq!(
                    p_a.get(i, j),
                    p_b.get(i, j),
                    "same (grid, params) must produce bit-identical output at ({i}, {j})"
                );
            }
        }
    }

    /// Test 2 — different seeds → different outputs.
    #[test]
    fn r7_boundary_displacement_different_seeds_differ() {
        let (mut p_a, seeds) = two_plate_voronoi(32, 32);
        let (mut p_b, _) = two_plate_voronoi(32, 32);

        let params_a = R7InitParams { seed: 42, ..R7InitParams::default() };
        let params_b = R7InitParams { seed: 1337, ..R7InitParams::default() };

        apply_boundary_displacement(&mut p_a, &seeds, &params_a);
        apply_boundary_displacement(&mut p_b, &seeds, &params_b);

        let mut diff_count = 0;
        for j in 0..32 {
            for i in 0..32 {
                if p_a.get(i, j) != p_b.get(i, j) {
                    diff_count += 1;
                }
            }
        }
        assert!(
            diff_count > 0,
            "different seeds (42 vs 1337) must produce different displacement; got 0 differences"
        );
    }

    /// Test 3 — `amplitude = 0` preserves Voronoï bit-identically.
    /// With zero amplitude the displacement vector is `(0, 0)`,
    /// the query position equals the cell centre, and the
    /// nearest-seed query returns the original Voronoï
    /// assignment.
    #[test]
    fn r7_boundary_displacement_amplitude_zero_preserves_voronoi() {
        let (mut p_displaced, seeds) = two_plate_voronoi(32, 32);
        let baseline = p_displaced.clone();
        let params = R7InitParams { amplitude: 0.0, ..R7InitParams::default() };

        apply_boundary_displacement(&mut p_displaced, &seeds, &params);

        for j in 0..32 {
            for i in 0..32 {
                assert_eq!(
                    p_displaced.get(i, j),
                    baseline.get(i, j),
                    "amplitude = 0 must preserve Voronoï at ({i}, {j})"
                );
            }
        }
    }

    /// Test 4 — `enabled = false` no-op (W4 closure-isolation).
    #[test]
    fn r7_boundary_displacement_disabled_no_op() {
        let (mut p_displaced, seeds) = two_plate_voronoi(32, 32);
        // Set up a non-uniform pre-state so a forgotten branch
        // can't pass on a default-zero plate_id.
        p_displaced.set(5, 5, 99);
        let baseline = p_displaced.clone();
        let params = R7InitParams {
            enabled: false,
            amplitude: 100.0, // pathological if it ran
            ..R7InitParams::default()
        };

        apply_boundary_displacement(&mut p_displaced, &seeds, &params);

        for j in 0..32 {
            for i in 0..32 {
                assert_eq!(
                    p_displaced.get(i, j),
                    baseline.get(i, j),
                    "`enabled = false` must not touch any cell at ({i}, {j})"
                );
            }
        }
    }

    /// Test 5 — boundaries deviate from Voronoï baseline, with
    /// healthy regime bounds. Counts reassigned cells (where
    /// post-displacement plate_id differs from pre-displacement)
    /// and asserts it stays in `(0, 20 %]` of total cells per the
    /// W7 surface threshold.
    #[test]
    fn r7_boundary_displacement_curves_boundaries() {
        let nx = 64;
        let ny = 64;
        let (mut p_displaced, seeds) = two_plate_voronoi(nx, ny);
        let baseline = p_displaced.clone();
        let params = R7InitParams::default();

        apply_boundary_displacement(&mut p_displaced, &seeds, &params);

        let total = nx * ny;
        let mut reassigned = 0;
        for j in 0..ny {
            for i in 0..nx {
                if p_displaced.get(i, j) != baseline.get(i, j) {
                    reassigned += 1;
                }
            }
        }
        let frac = 100.0 * reassigned as f64 / total as f64;
        eprintln!(
            "r7_boundary_displacement_curves_boundaries: reassigned = {reassigned} / {total} ({frac:.2} %)"
        );

        assert!(
            reassigned > 0,
            "displacement must reassign at least one boundary cell (got 0); \
             check amplitude / frequency defaults vs grid size"
        );
        let upper_bound = total / 5; // 20 %
        assert!(
            reassigned < upper_bound,
            "reassignment must stay under 20 % of cells (got {reassigned} / {total} = \
             {frac:.2} %); amplitude default may be too large for grid size {nx}"
        );
    }

    /// Test 6 — total cell count unchanged (no creation /
    /// destruction). Implicit because the function operates
    /// in-place on a fixed grid, but locked explicitly because
    /// any future "split a cell" extension would break this
    /// invariant silently.
    #[test]
    fn r7_boundary_displacement_preserves_total_cells() {
        let nx = 32;
        let ny = 32;
        let (mut p_displaced, seeds) = two_plate_voronoi(nx, ny);
        let pre_count = nx * ny;
        let params = R7InitParams::default();

        apply_boundary_displacement(&mut p_displaced, &seeds, &params);

        let post_count = p_displaced.data().len();
        assert_eq!(
            post_count, pre_count,
            "total cell count must be invariant under displacement; pre = {pre_count}, \
             post = {post_count}"
        );
    }
}
