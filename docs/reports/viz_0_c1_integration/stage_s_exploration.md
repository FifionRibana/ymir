# Viz-0 C1 Integration — Stage S Exploration

Issue #137 — codebase exploration + bridge/v2 template alignment.

## C1 core API (Stage S grep results)

### Init + kinematics

- `crates/ymir-core/src/tectonics_c1/init_r7/mod.rs:181-185`:
  ```rust
  pub fn init_c1_state_phase_2_r7(
      grid_size: usize,
      seed: u64,
      params: &Phase2InitParams,
  ) -> C1State
  ```
- `crates/ymir-core/src/tectonics_c1/kinematics.rs:40`:
  `pub fn preset_phase_1_1(num_plates: usize) -> PlateKinematics`.

### Isostasy

- `crates/ymir-core/src/tectonics/isostasy.rs:15` — `pub struct IsostasyConfig` (8 fields: `rho_crust = 2750`, `rho_mantle = 3300`, `rho_water = 1025`, etc.). `impl Default` at line 36.
- `crates/ymir-core/src/tectonics/isostasy.rs:69` — `pub fn compute_isostasy(thickness: &Field2D, config: &IsostasyConfig) -> IsostasyResult`.

### Stein-Stein bathymetry (Architecture C live-altitude derivation)

- `crates/ymir-core/src/tectonics_c1/closures/oceanic_bathymetry/source_term.rs:126`:
  ```rust
  pub fn apply_stein_stein_bathymetry(
      altitude: &mut GridF32,
      age: &Field2D,
      plate_type: &PlateTypeField,
      params: &SteinSteinParams,
  )
  ```

### Time loop config + closure bundle

- `crates/ymir-core/src/tectonics_c1/time_loop.rs:87` — `pub struct C1TimeLoopConfig { n_steps, dx, dy, iso_config, drainage_max_distance }`.
- `crates/ymir-core/src/tectonics_c1/time_loop.rs:236` — `pub struct C1Closures` (7 fields including Track D trio).

## v2 bridge template (the cloning source)

`crates/ymir-viz/src/bridge/v2/` contains 9 files:
- `mod.rs` — module declaration + re-exports.
- `commands.rs` — `V2Command` enum.
- `events.rs` — `V2Event` enum + `V2FinalState` struct (per-cell `s_field`, `vx`, `vy`, `plate_id`, `plate_type` as `Vec<u8>` encoded 0/1).
- `snapshot.rs` — `V2RunSnapshot` + `V2ScalarMetrics`.
- `spec.rs` — `V2RunSpec` (huge config — C1 spec will be far smaller).
- `presets.rs` — JSON preset loader.
- `build_config.rs` — spec → `BaselineConfig`.
- `thread.rs` — `spawn_v2_thread` worker entry.
- `plugin.rs` — `V2BridgePlugin: Plugin` registers via `app.insert_resource(V2SolverBridge)` + `app.add_systems(Update, poll_v2_events)`.

Bridge module exposure: `crates/ymir-viz/src/bridge/mod.rs` = `pub mod v2;` — I'll add `pub mod c1;` alongside.

### Channel sizing (mirror)

`bounded::<V2Command>(4)` + `bounded::<V2Event>(256)` — Viz-0 plan says channel N=2 for events; v2 uses 256 because it ingests `Progress` per step at sub-second cadence. For C1 with per-step events, N=2 forces tight backpressure (= pause semantics per W1 Stage E2). Will pick N=2 unless Stage E2 testing surfaces issues.

### Cancel pattern

v2 uses `Arc<AtomicBool>` checked inside the harness step callback. C1's `run_with_closures` already takes a `FnMut(usize, &C1State)` callback — perfect injection point for cancel check. No core change needed for cancellation.

## Snapshot architecture (Q1.2 raw-fields)

`V2FinalState` (events.rs:104) carries per-cell `Vec<f64>` for `s_field`, `vx`, `vy`, `strain_rate_invariant`, optional `age_field`, `cratonic_factor`, `plate_id` (`Vec<u16>`), `plate_type` (`Vec<u8>` encoded `0 = Oceanic`, `1 = Continental`), `boundary_flag` (`Vec<u8>`).

For C1 snapshot (Q1.2 confirmed raw-fields):
- `s`, `age` — `Vec<f64>` (from `Field2D::data().to_vec()`).
- `plate_id` — `Vec<u16>` (from `PlateIdField::data().to_vec()`).
- `plate_type` — `Vec<u8>` (manual match: `PlateType::Oceanic → 0`, `Continental → 1`; the enum has no `#[repr(u8)]`).
- `cratonic_mask` — `Vec<u8>` (false → 0, true → 1; from `BoolField`).
- `vx`, `vy` per-cell — derived as `kinematics.velocities[plate_id[c]]` expanded to `Vec<f64>` (mirror v2 for velocity overlay reuse — `overlay::draw_velocity_vectors` takes per-cell vx/vy already, see W2 below).
- `stats: C1StepStats` — the Viz-D0 diagnostic field.
- `grid_size`, `dx`, `dy` — scalars.

### Velocity overlay (W2 Stage S answer)

`crates/ymir-viz/src/visualization/overlay.rs:56`:
```rust
pub fn draw_velocity_vectors(
    rgba: &mut [u8],
    nx: usize, ny: usize,
    vx: &[f64],
    vy: &[f64],
    plate_id: &[u16],
    arrow_scale_cells: f64,
)
```

Takes per-cell `vx`/`vy` + `plate_id`. **Decision**: C1Snapshot carries per-cell `vx`/`vy` expanded from `kinematics.velocities[plate_id[c]]` at snapshot creation. Cost ~32KB at 64² (acceptable). Reuses the v2 overlay function bit-identically. Track D accretion mutates `kinematics.velocities` per merge event — snapshot must be re-expanded each step (no caching).

Alternative (not chosen): carry `plate_id` + `kinematics: Vec<(f64, f64)>` and reconstruct vx/vy in the UI render system. Saves bandwidth but forks the overlay API. Reject — match-request-scope says reuse.

## PlateType encoding (W4 Stage S)

`crates/ymir-core/src/tectonics_v2/boundaries/plate_type.rs:16`:
```rust
pub enum PlateType {
    Oceanic,
    Continental,
}
```

NO `#[repr(u8)]`. Encoding for snapshot Vec<u8>:
```rust
fn plate_type_to_u8(t: PlateType) -> u8 {
    match t {
        PlateType::Oceanic => 0,
        PlateType::Continental => 1,
    }
}
```

Mirrors `V2FinalState` convention (`0 = Oceanic`, `1 = Continental`).

## W4 Stage E1 — `C1State { ... }` struct-literal call-sites

Grep `C1State \{` against `crates/ymir-core`:

| File:Line | Context |
|-----------|---------|
| `tectonics_c1/init.rs:99` | `init_c1_state_phase_1_1` constructor |
| `tectonics_c1/init_r7/mod.rs:264` | `init_c1_state_phase_2_r7` constructor |
| `tectonics_c1/time_loop.rs:767` | test helper `uniform_single_plate_state` |

**3 sites total.** All in `ymir-core`. None in `ymir-viz` (viz uses snapshots, not direct C1State construction). All 3 will need `last_step_stats: Default::default()` (or `..Default::default()` if I add `Default` impl).

No struct literals in test files or other tests — manageable mechanical change.

## Stage S W7 — questions surfacing for user before Stage E1

**Q-S.1 — C1StepStats location.** `tectonics_c1/state.rs` (alongside C1State) vs new `tectonics_c1/stats.rs` module file. Recommendation: **new `tectonics_c1/stats.rs`** for cleaner separation; state.rs already documents the dynamical/classification field grouping. Surface re-exports from `tectonics_c1/mod.rs`.

**Q-S.2 — C1StepStats Default vs Copy.** Track D's stats include `Vec<u16>` (`new_plate_ids_created` in `RiftingSplitStats`), preventing Copy. C1StepStats has 4 sub-stats inside it — would Clone+Debug+Default suffice (no Copy)? Recommendation: **Clone+Debug+Default**. The field on C1State is mutated in-place each step; no Copy needed.

**Q-S.3 — `C1Snapshot::from_state` allocations.** Each `Field2D::data().to_vec()` + `PlateIdField::data().to_vec()` + plate_type expansion is ~32-128KB at 64². N=2 channel = up to 4 snapshots buffered (current + 2 in flight + 1 receiving). 256KB-512KB transient memory. Acceptable.

**Q-S.4 — C1RunSpec configuration surface.** Minimal initial scope per match-request-scope:
- `grid_size: usize` (64 default)
- `seed: u64` (42 default)
- `n_steps: usize` (300 default)
- `init_params: Phase2InitParams` (default)
- `closures: C1Closures` (default, all enabled — including Track D)
- `kinematics_preset: Phase1_1` (only one preset for now)

Phase 2 R7 init is the only init mode for Viz-0. Track C kinematics presets deferred to a Viz-0-bis if needed.

**Q-S.5 — V2Field analog `C1Field`.** 5 modes per design exploration:
- `S` (crustal thickness, `[0, 3.0]` palette)
- `Age` (cell age, auto-normalized)
- `PlateId` (categorical, hash-modulo color)
- `PlateType` (2-color: cyan for Oceanic, beige for Continental)
- `Altitude` (Architecture C derived from `s`+`age`+`plate_type`, `[-1.13, +1.13]` symmetric bipolar)

Reuse `visualization/colormap.rs` for S / Age / Altitude. New 2-color palette for PlateType. PlateId hash-color: small new helper.

## Stage S W7 — answers I plan to commit if no user pushback

| Item | Answer |
|------|--------|
| C1StepStats location | new `tectonics_c1/stats.rs`, re-exported from `tectonics_c1/mod.rs` |
| C1StepStats derive | `Clone, Debug, Default` (NOT `Copy` — Vec<u16> in rifting splits) |
| Snapshot allocations | accept ~32-128KB per snapshot at 64²; N=2 channel ≤ 4 in-flight |
| Channel N | 2 (matches design Q1; ≤ 4 in-flight snapshots) |
| Velocity overlay | C1Snapshot carries expanded per-cell vx/vy (reuse v2 overlay bit-identically) |
| PlateType encoding | match-expr: Oceanic → 0, Continental → 1 |

## Architecture C live-altitude derivation (Stage E4 surface)

Per Stage E4 plan: per-frame altitude reconstruction from snapshot:
1. Reconstruct `Field2D` from `snapshot.s: Vec<f64>` (need `Field2D::from_vec(nx, ny, vec)` — confirm exists Stage E4; if not, build via `Field2D::new` + per-cell `set`).
2. `compute_isostasy(&s_field, &IsostasyConfig::default()) → IsostasyResult { heightmap, ... }`.
3. `let mut altitude = heightmap.clone();`
4. Reconstruct `Field2D` for age + `PlateTypeField` for plate_type (need to expand `Vec<u8>` back to enum via `match`).
5. `apply_stein_stein_bathymetry(&mut altitude, &age_field, &plate_type_field, &SteinSteinParams::default())`.
6. Apply altitude palette `[-1.13, +1.13]` symmetric → RGBA → Bevy Image update.

Match the Track A gallery code path verbatim (see `c1_phase_2_track_b_visual_gallery.rs::dump_snapshot`).

## Track D dynamic plate_id — cardinality-agnostic overlay handling (W6 global)

Track D mutates `plate_id` per-step (subduction reassignment, accretion merge, rifting split). The viz overlay must NOT assume static plate count:
- `draw_voronoi_boundaries` reads plate_id 4-connected, draws line where neighbour differs — naturally cardinality-agnostic. Reuse as-is.
- `draw_velocity_vectors` aggregates by plate_id key, draws one arrow per distinct id. Naturally cardinality-agnostic. Reuse as-is.
- PlateId field mode: hash-color by id value (small int mod palette). Cardinality-agnostic.
- New plate ids from rifting splits naturally render at their new positions.

## Effort estimate update

Stage S findings confirm 5-8 day estimate stands. No surprise hits.

## Out-of-scope deferrals (match-request-scope)

- Live closure toggle (run-locked per design Q3).
- Boundary-type overlay (fast-follow).
- Model 2 / multi-model selector (Viz-0 is C1 only).
- Workflow Phase A/B (C1 has no equivalent workflow yet).
- Track C kinematics presets UI.
- C1RunSpec JSON preset loader (Viz-0 uses programmatic defaults; presets land in Viz-0-bis if needed).
