//! Coarse tectonic DEBUG labels — a read-only derivation of the causal masks
//! (rift, subduction, collision, craton, lithology class) from a `C1State`, for the
//! viz "microscope" overlay so the author can SEE where each closure acts (where the
//! rifts are, which margins subduct, the hard craton, the lithology contrast).
//!
//! Not part of the production pipeline: it reads the same signals the closures read
//! (plate type, age, boundary classification, cratonic mask) and exposes them as
//! coarse per-cell masks the viz upscales to the HD window with the terrain's own
//! `(sample_origin, sample_size)` mapping — so every overlay is pixel-registered with
//! the terrain it explains.

use crate::tectonics_c1::boundary_classification::{
    BoundaryType, classify_boundaries, oc_override_seed_mask,
};
use crate::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use crate::tectonics_v2::boundaries::plate_type::PlateType;

/// Run the coarse tectonic pass and return the settled `(C1State, PlateKinematics)` —
/// the SAME derivation production erodes (`cached_product` miss path), so the labels
/// register with the terrain. Cheap (64² × ~300 steps, ~1 s) vs the HD erosion; the
/// viz calls it to feed the overlay. Deterministic in `seed`.
#[must_use]
pub fn run_coarse_tectonics(
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
) -> (C1State, PlateKinematics) {
    let mut state = init_c1_state_phase_2_r7(grid, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, run, closures, |_, _| {});
    (state, kin)
}

/// Per-cell lithology CLASS on the coarse grid (matches the C-3 K field's causal
/// source: hard basement, rift-soft, volcaniclastic). Volcaniclastic is stamped at HD
/// from edifice footprints, so at the COARSE grid only hard/rift are known — the viz
/// stamps volcanic separately from the edifice list (as the K field does).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LithoClass {
    /// Hard basement (crystalline + metasediments) — the erodibility reference.
    Hard,
    /// Rift-soft (young, `age ≈ 0`) — the C-3 soft class.
    RiftSoft,
}

/// Coarse tectonic labels for the overlay, all row-major `nx·ny` (the tectonic grid,
/// typically 64²). Bools are cheap; the viz nearest-samples them to the HD window.
#[derive(Clone, Debug)]
pub struct CoarseTectonicLabels {
    pub nx: usize,
    pub ny: usize,
    /// `true` = continental, `false` = oceanic (the binary C1 plate type).
    pub continental: Vec<bool>,
    /// Cratonic shield mask (the hard, old continental interior).
    pub craton: Vec<bool>,
    /// Rift-soft: continental cells with `age < 1` (the causal C-3 soft signal — the
    /// same cells `lithology::build_coarse_k` softens).
    pub rift: Vec<bool>,
    /// Subduction — overriding CONTINENTAL margin (O-C convergent, continental side;
    /// where arc volcanism seeds). Andes-type upper plate.
    pub subduction_upper: Vec<bool>,
    /// Subduction — subducting OCEANIC slab (O-C convergent, oceanic side).
    pub subduction_slab: Vec<bool>,
    /// Continental COLLISION (C-C convergent) — the accretion / orogenic-thickening
    /// proxy (accretion is a per-step event, not a static field; the standing C-C
    /// convergent boundary is its geographic footprint).
    pub collision: Vec<bool>,
    /// Divergent boundaries (spreading / oceanic rift axes).
    pub divergent: Vec<bool>,
    /// Advected crustal age, normalised to `[0,1]` by the field max (for a heat
    /// overlay: young = 0, old = 1). Degenerate on the current model (see ADR C-3),
    /// exposed anyway so the author can confirm that with their own eyes.
    pub age_norm: Vec<f32>,
}

impl CoarseTectonicLabels {
    #[inline]
    #[must_use]
    pub fn litho_class(&self, i: usize, j: usize) -> LithoClass {
        let k = j * self.nx + i;
        if self.rift[k] { LithoClass::RiftSoft } else { LithoClass::Hard }
    }
}

/// Derive the coarse tectonic labels from a settled `C1State` + its kinematics — the
/// same state production erodes. Read-only.
#[must_use]
pub fn derive_tectonic_labels(state: &C1State, kin: &PlateKinematics) -> CoarseTectonicLabels {
    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let n = nx * ny;
    let info = classify_boundaries(&state.plate_id, kin);
    // Overriding continental margin (O-C, continental side) = the subduction upper
    // plate — the same seed mask arc volcanism uses.
    let upper = oc_override_seed_mask(&info, &state.plate_id, &state.plate_type);

    let mut continental = vec![false; n];
    let mut craton = vec![false; n];
    let mut rift = vec![false; n];
    let mut subduction_upper = vec![false; n];
    let mut subduction_slab = vec![false; n];
    let mut collision = vec![false; n];
    let mut divergent = vec![false; n];
    let mut age_norm = vec![0.0f32; n];

    let mut age_max = 1e-6f32;
    for j in 0..ny {
        for i in 0..nx {
            age_max = age_max.max(state.age.get(i, j) as f32);
        }
    }

    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let cont = matches!(state.plate_type.get(i, j), PlateType::Continental);
            continental[k] = cont;
            craton[k] = state.cratonic_mask.get(i, j);
            let age = state.age.get(i, j) as f32;
            age_norm[k] = age / age_max;
            if cont && age < 1.0 {
                rift[k] = true;
            }
            match info.boundary_type.get(i, j) {
                BoundaryType::Divergent => divergent[k] = true,
                BoundaryType::Convergent => {
                    if upper.get(i, j) {
                        subduction_upper[k] = true;
                    } else {
                        // Convergent, not an O-C continental override → either the
                        // subducting oceanic slab (O-C oceanic side) or a C-C collision.
                        if cont {
                            collision[k] = true; // continental side of a C-C convergence
                        } else {
                            subduction_slab[k] = true; // oceanic side (subducting)
                        }
                    }
                }
                _ => {}
            }
        }
    }

    CoarseTectonicLabels {
        nx,
        ny,
        continental,
        craton,
        rift,
        subduction_upper,
        subduction_slab,
        collision,
        divergent,
        age_norm,
    }
}
