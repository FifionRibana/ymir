//! Initial condition for C1 Phase 1.1 — v2 init reused verbatim.
//!
//! ## What this module does
//!
//! Combines two v2 pieces preserved per §4.8 of the design doc:
//!
//! - [`crate::tectonics_v2::voronoi::generate_voronoi`] →
//!   `plate_id`, `plate_type`, `seed_coords` (and the implicit
//!   plate count).
//! - [`crate::tectonics_v2::init::init_s_field`] →
//!   the initial `S̃` field. Phase 1.1 uses [`InitMode::default()`]
//!   (`Uniform { boundary_smoothing_width: 1.0 }`) — sufficient
//!   transport-correctness signal without orogenic structure
//!   that closures haven't yet been built to interpret.
//!
//! The `age` field is hand-initialised here (continental = 7.0,
//! oceanic = 0.5, non-dim) to keep Phase 1.1 self-contained;
//! C1's age-field handling matures alongside Parsons-Sclater
//! bathymetry in Phase 2 (closure #3, §5.1 design doc).
//!
//! The cratonic mask is a Phase 1.1-local stand-in: continental
//! plates whose seed-coordinate `x` lies in the lower half of
//! the domain are tagged cratonic. This is **not** the future
//! C1 cratonic rule — it's a deterministic binary mask that
//! lets transport-rigidity be tested visually without depending
//! on the retired v2 cratonic-factor amplifier (§4.4 design doc).
//! The real C1 cratonic rule will be set when a closure needs it.
//!
//! ## Out of scope this phase
//!
//! - R7-generalised init (boundary perturbation, continental
//!   clustering): Phase 2.
//! - Cratonic rule based on age threshold or random-seed subset
//!   of continental plates: revisited when first closure needs it.

use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::init::{init_s_field, InitContext, InitMode, PlateInitData};
use crate::tectonics_v2::voronoi::{generate_voronoi, PlateIdField, VoronoiConfig};

use super::state::{BoolField, C1State};

/// Default continental age (non-dim) used by [`init_c1_state_phase_1_1`].
pub const CONTINENTAL_AGE_INIT: f64 = 7.0;
/// Default oceanic age (non-dim) used by [`init_c1_state_phase_1_1`].
pub const OCEANIC_AGE_INIT: f64 = 0.5;

/// Build a fresh [`C1State`] for Phase 1.1.
///
/// Wraps `generate_voronoi` (default 8 plates, 30% continental
/// ratio) and `init_s_field` ([`InitMode::default()`]). The age
/// field and cratonic mask are filled per the module docstring —
/// both stand-ins until later phases need different behaviour.
pub fn init_c1_state_phase_1_1(grid_size: usize, seed: u64) -> C1State {
    let nx = grid_size;
    let ny = grid_size;

    // Voronoï tessellation → plate_id, plate_type, seed coords.
    let vor_config = VoronoiConfig::default();
    let plates = generate_voronoi(nx, ny, &vor_config, seed);

    // Initial `S̃` field via v2's init dispatch. `Uniform { 1.0 }`
    // is the [`InitMode::default()`] — flat per-plate-type with
    // 1-cell smoothstep at inter-plate boundaries.
    let init_ctx = InitContext {
        nx,
        ny,
        seed,
        amplitude: 0.0,
        plate_data: Some(PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        }),
    };
    let s = init_s_field(InitMode::default(), &init_ctx);

    // Age: continental cells get CONTINENTAL_AGE_INIT, oceanic
    // get OCEANIC_AGE_INIT.
    let mut age = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let age_value = match plates.plate_type.get(i, j) {
                PlateType::Continental => CONTINENTAL_AGE_INIT,
                PlateType::Oceanic => OCEANIC_AGE_INIT,
            };
            age.set(i, j, age_value);
        }
    }

    let cratonic_mask = build_phase_1_1_cratonic_mask(
        nx,
        ny,
        &plates.plate_id,
        &plates.plate_type,
        &plates.seed_coords,
    );

    C1State {
        s,
        age,
        plate_id: plates.plate_id,
        plate_type: plates.plate_type,
        cratonic_mask,
        num_plates: plates.num_plates,
    }
}

/// Build the Phase 1.1 cratonic mask: continental cells whose
/// plate's seed-coordinate `x` is in the lower half of the
/// domain are marked `true`. Oceanic cells and continental
/// cells from "upper-half" plates are `false`.
fn build_phase_1_1_cratonic_mask(
    nx: usize,
    ny: usize,
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    seed_coords: &[(f64, f64)],
) -> BoolField {
    let x_threshold = nx as f64 / 2.0;
    let cratonic_plate: Vec<bool> =
        seed_coords.iter().map(|(x, _y)| *x < x_threshold).collect();

    let mut mask = BoolField::filled(nx, ny, false);
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_id.get(i, j) as usize;
            let is_continental = matches!(plate_type.get(i, j), PlateType::Continental);
            if is_continental && cratonic_plate.get(pid).copied().unwrap_or(false) {
                mask.set(i, j, true);
            }
        }
    }
    mask
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_produces_8_plate_default_layout() {
        let state = init_c1_state_phase_1_1(64, 42);
        assert_eq!(state.num_plates, 8);
        assert_eq!(state.s.nx(), 64);
        assert_eq!(state.s.ny(), 64);
        assert_eq!(state.age.nx(), 64);
        assert_eq!(state.age.ny(), 64);
        assert_eq!(state.cratonic_mask.nx(), 64);
        assert_eq!(state.cratonic_mask.ny(), 64);
    }

    #[test]
    fn age_distinguishes_continental_from_oceanic() {
        let state = init_c1_state_phase_1_1(64, 42);
        let mut saw_continental = false;
        let mut saw_oceanic = false;
        for j in 0..state.ny() {
            for i in 0..state.nx() {
                let a = state.age.get(i, j);
                match state.plate_type.get(i, j) {
                    PlateType::Continental => {
                        assert!((a - CONTINENTAL_AGE_INIT).abs() < 1e-12);
                        saw_continental = true;
                    }
                    PlateType::Oceanic => {
                        assert!((a - OCEANIC_AGE_INIT).abs() < 1e-12);
                        saw_oceanic = true;
                    }
                }
            }
        }
        assert!(saw_continental && saw_oceanic, "preset must produce both plate types");
    }

    #[test]
    fn cratonic_mask_is_subset_of_continental_cells() {
        let state = init_c1_state_phase_1_1(64, 42);
        for j in 0..state.ny() {
            for i in 0..state.nx() {
                if state.cratonic_mask.get(i, j) {
                    assert!(
                        matches!(state.plate_type.get(i, j), PlateType::Continental),
                        "cratonic cell at ({i},{j}) is not continental"
                    );
                }
            }
        }
    }
}
