//! `C1Snapshot` — thread-safe per-step state shipped from the
//! `bridge::c1` worker thread back to the Bevy UI.
//!
//! ## Raw-fields design (Issue #137 Q1.2)
//!
//! The snapshot carries **raw scalar fields** (`Vec<f64>` for `s`/`age`,
//! `Vec<u16>` for `plate_id`, `Vec<u8>` for `plate_type`/`cratonic_mask`)
//! — not pre-rendered RGBA. The Bevy UI thread renders fields on demand
//! from the cached snapshot, which enables view-switch during pause (W3
//! global watchpoint).
//!
//! ## Velocity field caveat (Issue #137 Viz-0-bis candidate)
//!
//! `plate_velocities` stores **init-time** per-plate velocities — NOT
//! the live mid-run kinematics. Track D's accretion mutates velocities
//! per merge event (mass-weighted average) and rifting splits push new
//! plate velocities; these mid-run mutations are NOT reflected in the
//! snapshot.
//!
//! Root cause: `run_with_closures` borrows `&mut PlateKinematics` for
//! the duration of the run; the `on_step` callback cannot also borrow
//! kinematics. The Stage E2 W7 surface investigated three workarounds
//! (signature change, field-on-state sibling, unsafe pointer) and
//! settled on init-time-only velocity capture for MVP. The worker
//! thread clones `kinematics.velocities` pre-run and the closure
//! captures the clone by reference (no borrow conflict because the
//! clone is a separate `Vec`).
//!
//! Acceptable because:
//! - Subduction (high-frequency) does NOT mutate kinematics.
//! - Rifting splits (rare, 0-3/seed per Track D Stage V) do change
//!   kinematics length; visualisation falls back to the init-time
//!   subset for newly-spawned plates.
//! - Accretion merges (moderate, 6-10/seed per Track D Stage V)
//!   change individual plate velocities; magnitude shifts are
//!   typically small (mass-weighted).
//!
//! Viz-0-bis follow-up: add `last_plate_velocities` sibling field on
//! `C1State` (mirroring `last_step_stats` Option B pattern), populated
//! by `run_with_closures` before each `on_step` call. Avoids signature
//! change while exposing live kinematics to the closure.
//!
//! ## PlateType encoding
//!
//! `plate_type[c] = 0` ↔ `PlateType::Oceanic`,
//! `plate_type[c] = 1` ↔ `PlateType::Continental` (mirrors
//! `V2FinalState.plate_type` convention).

use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::stats::C1StepStats;
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

#[derive(Clone, Debug)]
pub struct C1Snapshot {
    /// 0-based step index. `step + 1` is the number of completed
    /// steps (the C1 time loop calls `on_step(step, state)` AFTER
    /// the `step`-th iteration's mutations).
    pub step: usize,
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    /// Crust thickness `S̃`, row-major `nx × ny`.
    pub s: Vec<f64>,
    /// Cell age, row-major. Path 3.A init-only and Path 3.B
    /// event-driven `= 0` semantics inherited from Track B/D.
    pub age: Vec<f64>,
    /// Voronoï plate index per cell.
    pub plate_id: Vec<u16>,
    /// Plate type: 0 = Oceanic, 1 = Continental.
    pub plate_type: Vec<u8>,
    /// Cratonic mask: 0 = false, 1 = true.
    pub cratonic_mask: Vec<u8>,
    /// Init-time plate count (cached on `C1State.num_plates`).
    /// Live count requires `live_plate_count` (re-scanned).
    pub num_plates: usize,
    /// Distinct plate ids in `plate_id` at this step. Drops on
    /// accretion merges, grows on rifting splits.
    pub live_plate_count: usize,
    /// Per-plate `(vx, vy)` velocities. **Init-time only** (see
    /// module docstring caveat).
    pub plate_velocities: Vec<(f64, f64)>,
    /// Track D per-step diagnostic stats (Issue #137 Viz-D0
    /// Option B field on `C1State`).
    pub stats: C1StepStats,
}

impl C1Snapshot {
    /// Build a snapshot from the current `C1State` + the pre-run
    /// kinematics velocities clone.
    ///
    /// `step` is the 0-based step index just completed (matching
    /// `run_with_closures`'s `on_step(step, state)` callback
    /// convention). For the pre-run / cycle-0 snapshot the caller
    /// passes `step = 0` (the convention is "first event has
    /// step = 0; n_steps events total after the run completes").
    pub fn from_state(step: usize, state: &C1State, plate_velocities: &[(f64, f64)]) -> Self {
        let nx = state.nx();
        let ny = state.ny();
        let n_cells = nx * ny;

        let s = state.s.data().to_vec();
        let age = state.age.data().to_vec();
        let plate_id = state.plate_id.data().to_vec();

        let mut plate_type = Vec::with_capacity(n_cells);
        for &t in state.plate_type.data() {
            plate_type.push(match t {
                PlateType::Oceanic => 0,
                PlateType::Continental => 1,
            });
        }

        let mut cratonic_mask = Vec::with_capacity(n_cells);
        for &b in state.cratonic_mask.data() {
            cratonic_mask.push(if b { 1 } else { 0 });
        }

        // Live count via HashSet over plate_id.
        let mut seen = std::collections::HashSet::new();
        for &pid in &plate_id {
            seen.insert(pid);
        }
        let live_plate_count = seen.len();

        Self {
            step,
            nx,
            ny,
            // Phase 1.x convention: unit-domain, dx = dy = 1 / nx.
            dx: 1.0 / nx as f64,
            dy: 1.0 / ny as f64,
            s,
            age,
            plate_id,
            plate_type,
            cratonic_mask,
            num_plates: state.num_plates,
            live_plate_count,
            plate_velocities: plate_velocities.to_vec(),
            stats: state.last_step_stats.clone(),
        }
    }

    /// Expand `plate_velocities[plate_id[c]]` to per-cell `vx`/`vy`
    /// arrays. Helper for the velocity-overlay render path which
    /// reuses `visualization::overlay::draw_velocity_vectors`
    /// (expects per-cell vx/vy).
    ///
    /// Returns `(vx, vy)` row-major `nx × ny`. Both length
    /// `self.nx * self.ny`.
    pub fn expand_per_cell_velocities(&self) -> (Vec<f64>, Vec<f64>) {
        let n_cells = self.nx * self.ny;
        let mut vx = Vec::with_capacity(n_cells);
        let mut vy = Vec::with_capacity(n_cells);
        for &pid in &self.plate_id {
            // Defensive — rifting splits can push new plate ids
            // beyond the init-time `plate_velocities.len()`. New
            // plates fall back to (0, 0) — a known MVP limitation
            // (see module docstring).
            let (vxp, vyp) = self.plate_velocities.get(pid as usize).copied().unwrap_or((0.0, 0.0));
            vx.push(vxp);
            vy.push(vyp);
        }
        (vx, vy)
    }

    /// Convenience: total live plate count (re-scanned). Equivalent
    /// to `self.live_plate_count` — surface for symmetry with the
    /// init-time `num_plates`.
    pub fn live_plates(&self) -> usize {
        self.live_plate_count
    }

    /// Plate kinematics provided by the worker (helper for callers
    /// that want a `PlateKinematics` object — e.g., the Stage E4
    /// `compute_isostasy` driver reconstructs a temporary one for
    /// the per-frame altitude derivation).
    pub fn kinematics(&self) -> PlateKinematics {
        PlateKinematics { velocities: self.plate_velocities.clone() }
    }
}
