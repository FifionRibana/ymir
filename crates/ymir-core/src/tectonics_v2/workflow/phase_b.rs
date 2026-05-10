//! Phase B — HD finalization (one-shot upscale + erosion).
//!
//! Takes the final-cycle output of Phase A (a low-res
//! [`BaselineResult`]), converts the `Field2D` S̃ field to `GridF32`
//! altitude via [`crate::tectonics::isostasy::compute_isostasy`], runs
//! [`crate::terrain::upscale::upscale_with_fbm`] to the configured
//! HD resolution, then [`crate::erosion::hydraulic::run_erosion`] on
//! the upscaled heightmap. Records the D5 grand-scale deviation
//! `‖HD_after - upscale(low_res)‖_∞`.
//!
//! Pipeline:
//!
//! ```text
//! Field2D S̃ (low-res, cellular)
//!     ↓ compute_isostasy
//! GridF32 altitude in [0, 1] + sea_level_normalized
//!     ↓ upscale_with_fbm (bicubic + anisotropic FBM)
//! GridF32 HD altitude (= "upscale(low_res)" baseline for the D5 metric)
//!     ↓ run_erosion (rain-drop, Beyer 2015)
//! GridF32 HD eroded altitude (= "HD_after")
//!     ↓ |HD_after - baseline|_∞
//! grand_scale_deviation
//! ```
//!
//! `Disabled` → `None` (Step 11 standalone contract).
//!
//! Phase B is one-shot: the user/orchestrator calls it once at the
//! end of the Phase A loop with the final cycle's `BaselineResult`.

use super::{PhaseBOutput, WorkflowConfig};
use crate::erosion::hydraulic::run_erosion;
use crate::grid::GridF32;
use crate::seed::WorldSeed;
use crate::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use crate::tectonics_v2::diagnostics::harness::BaselineResult;
use crate::terrain::upscale::upscale_with_fbm;

/// Run Phase B HD finalization on a Phase A output.
///
/// `Disabled` → returns `None`. The user is expected to consume the
/// low-res `BaselineResult` directly (Step 11 contract).
///
/// `Enabled(_)` → returns `Some(PhaseBOutput)`. The pipeline reads
/// `params.phase_b` (HD grid size, FBM config, erosion config, the
/// D5 tolerance is *not* enforced inside the pipeline — the caller
/// or test layer compares
/// `output.grand_scale_deviation < params.phase_b.grand_scale_tolerance`).
///
/// `seed` is shared between the upscale's FBM noise generators and
/// the rain-drop erosion's per-batch RNG. Using a single seed keeps
/// determinism predictable; the caller typically passes the host
/// `BaselineConfig.seed` so a workflow run is reproducible end-to-end
/// from a single user-supplied seed.
pub fn run_phase_b(
    input: &BaselineResult,
    wf: &WorkflowConfig,
    seed: u64,
) -> Option<PhaseBOutput> {
    match wf {
        WorkflowConfig::Disabled => None,
        WorkflowConfig::Enabled(params) => {
            let pb = &params.phase_b;
            // Step 1: S̃ → altitude. compute_isostasy returns a `[0, 1]`
            // GridF32 with `sea_level_normalized` known. The same call
            // we deliberately removed from the Phase A orchestrator
            // (Phase 3.5 finding 2: heightmap-space sea_level wrong
            // for S̃-space erosion) is correctly used here, because
            // upscale_with_fbm + run_erosion both operate in heightmap
            // [0, 1] space.
            let isostasy =
                compute_isostasy(&input.final_state.s_field, &IsostasyConfig::default());
            let coarse = isostasy.heightmap;
            let sea_level = isostasy.sea_level_normalized;

            // Step 2: Upscale with FBM injection. `pb.fbm.target_size`
            // is overridden by `pb.hd_grid_size` so the per-test
            // grid size is single-source-of-truth.
            let mut fbm_cfg = pb.fbm.clone();
            fbm_cfg.target_size = pb.hd_grid_size;
            let world_seed = WorldSeed::new(seed);
            let upscaled = upscale_with_fbm(&coarse, sea_level, &world_seed, &fbm_cfg);
            // Snapshot the post-upscale, pre-erosion HD heightmap as
            // the D5 grand-scale baseline reference.
            let baseline_hd = upscaled.heightmap.clone();
            let slope_hd = upscaled.slope;

            // Step 3: HD erosion. `sea_level` is overridden into the
            // erosion config so the rain-drop's coastal-deposition
            // and droplet-spawn logic sees the same threshold the
            // upscale used.
            let mut erosion_cfg = pb.erosion.clone();
            erosion_cfg.sea_level = sea_level;
            let erosion_result =
                run_erosion(&baseline_hd, &erosion_cfg, &world_seed, |_, _, _| true);

            // Step 4: D5 grand-scale deviation. Both metrics computed
            // in one pass via `compute_deviation_stats` —
            // `grand_scale_deviation` (L_∞) for diagnostic logging,
            // `grand_scale_deviation_p95` for the formal acceptance.
            // See PhaseBParams::grand_scale_tolerance docstring for
            // the Phase 5 reformulation rationale.
            let (grand_scale_deviation, grand_scale_deviation_p95) =
                compute_deviation_stats(&baseline_hd, &erosion_result.heightmap);

            Some(PhaseBOutput {
                heightmap: erosion_result.heightmap,
                sediment: erosion_result.sediment,
                slope: slope_hd,
                grand_scale_deviation,
                grand_scale_deviation_p95,
            })
        }
    }
}

/// Compute `(L_∞, p95)` of the per-cell `|a - b|` distribution in a
/// single pass over the slices. Both grids must have identical
/// shape — verified by the call site; here we just zip-iterate.
///
/// p95 is computed via partial sort (`select_nth_unstable`). The
/// allocation is `O(n)` for the deltas vec; not optimised for
/// memory because Phase B HD outputs are at most a few MB at 2048².
fn compute_deviation_stats(a: &GridF32, b: &GridF32) -> (f64, f64) {
    debug_assert_eq!(a.data.len(), b.data.len());
    let n = a.data.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mut deltas: Vec<f64> = a
        .data
        .iter()
        .zip(b.data.iter())
        .map(|(&x, &y)| (x - y).abs() as f64)
        .collect();
    let max = deltas.iter().copied().fold(0.0_f64, f64::max);
    // p95: index = floor(0.95 · (n - 1)). select_nth_unstable
    // partitions the slice in O(n) so we don't pay the full sort cost.
    let p95_idx = ((0.95 * (n - 1) as f64) as usize).min(n - 1);
    deltas
        .select_nth_unstable_by(p95_idx, |x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let p95 = deltas[p95_idx];
    (max, p95)
}
