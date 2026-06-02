# Issue #139 — Viz-0.5 Stage S exploration

Hover-to-inspect + workflow-mode pipeline + hypsometric lens + continuation-ready structure. Viz-only (no `tectonics_c1` changes; 9th bit-identical guard at Final). Branch `139-viz-05-hover-to-inspect-workflow-mode-pipeline-phase-3-validation-instruments` off `milestone/c1-lightweight-dynamic-tectonics` (Issue #137 merged via PR #138).

This document records the four Stage S confirmations the workflow requires before any code lands.

## 1. Cursor → cell transform

**Finding**: the world-space cursor already exists. [`camera.rs`](../../../crates/ymir-viz/src/camera.rs) defines `CursorWorldPos { pos: Option<Vec2> }` (a `Resource`), recomputed every frame in `update_cursor_world_pos` via `camera.viewport_to_world_2d(camera_transform, cursor)`. This is **world space, already camera-corrected** (pan + zoom folded in). The hover system reads this resource — no new cursor plumbing needed.

**Inversion math**: the C1 sprite ([`c1_plugin.rs`](../../../crates/ymir-viz/src/visualization/c1_plugin.rs)) is spawned at `Transform::from_xyz(0.0, 0.0, C1_SPRITE_Z=12.0)` — **centred at the world origin**. Its `custom_size = sprite_size(nx, ny)`:

```rust
let longer = nx.max(ny) as f32;
Vec2::new(BASE*nx/longer, BASE*ny/longer)   // BASE = C1_SPRITE_BASE_SIZE = 600.0
```

with a **nearest** sampler over the `nx × ny` texel image. Because the sprite is at the origin, `sprite_local = world` directly. Mapping world → cell:

```text
half = size / 2
u = (world.x + half.x) / size.x        // [0,1] left→right
v = (half.y - world.y) / size.y        // [0,1] top→bottom  (Bevy world +Y up; texel row 0 = top)
i = floor(u * nx)
j = floor(v * ny)
hover valid iff 0 ≤ u,v < 1  → else None (cursor outside sprite)
```

**W1 (cursor→cell at grid resize)**: `size` is recomputed from the *current* snapshot `nx/ny`, so the transform is resize-correct as long as the hover system reads `snapshot.nx/ny` (not a cached grid size). The `update_c1_texture` resize path already keys the sprite `custom_size` off `sprite_size(nx, ny)`, so the hover system must mirror that exact formula. Verification: cross-check a known cell (e.g. cursor at world origin → centre cell `(nx/2, ny/2)`).

**Y-axis caveat**: Bevy sprite texel row 0 renders at the **top**; world +Y points up. Hence `v = (half.y − world.y)/size.y` (NOT `(world.y + half.y)`). The render path writes `rgba[(j*nx+i)*4]` row-major with j=0 first → top row, consistent with this inversion. To verify, the W1 known-cell check must use an **off-centre, non-symmetric** cell to catch a Y-flip (a centre-cell check passes under either convention).

## 2. Split-borrow for the workflow worker

**Finding**: confirmed already compiling in [`phase_a_c1.rs:152`](../../../crates/ymir-core/src/tectonics_v2/workflow/phase_a_c1.rs). `apply_post_tectonic` takes:

```rust
pub struct PostTectonicInput<'a> {
    pub s_field: &'a mut Field2D,                 // &mut state.s
    pub plate_id: Option<&'a PlateIdField>,       // Some(&state.plate_id)
    pub plate_type: Option<&'a mut PlateTypeField>,// Some(&mut state.plate_type)
    ...
    pub params: &'a PhaseAParams,                 // ← NOTE: PhaseAParams, not WorkflowParams
    pub iso_cfg: &'a IsostasyConfig,
    pub cratonic_cfg: Option<&'a CratonicConfigEnabled>,
}
```

The struct-literal holds three **disjoint** borrows of `state` (`&mut s`, `&plate_id`, `&mut plate_type`) — the borrow checker accepts this under the splitting-borrow rule. The existing `run_phase_a_cycle_c1` does exactly this, so the worker can replicate it verbatim.

**W3 correction to the issue sketch**: `PostTectonicInput.params` is `&PhaseAParams`, NOT `&WorkflowParams`. The worker's `workflow_params` field must therefore be a **`PhaseAParams`** (the issue text says "workflow_params defaults to `PhaseAParams::default()`" — consistent; just be precise about the type). `WorkflowParams` wraps `{ phase_a: PhaseAParams, phase_b: PhaseBParams }`; we only need `phase_a` (Phase B is HD finalization, out of scope).

## 3. Calibration anchor (the A1-c guard)

**Finding**: [`PhaseAParams::default()`](../../../crates/ymir-core/src/tectonics_v2/workflow/mod.rs#L152) in `workflow/mod.rs`:

```rust
n_cycles = 5
k_cycle  = 20
alpha    = 0.01
isostatic_rebound_ratio = 0.80
max_drainage_distance   = 10
```

**Total tectonic steps = `n_cycles × k_cycle` = 100**, with **5** `macro_redistribution` passes (one per cycle). This is the calibrated cadence. A1-c failed by running `apply_post_tectonic` 6× over 300 steps (`steps_per_cycle=50`) — `macro_redistribution`'s `alpha=0.01` is calibrated for `k_cycle=20`; raising the per-pass step count (or pass frequency) without lowering `alpha` over-erodes. The workflow worker MUST use these defaults and compute total as `n_cycles × k_cycle`, **never** `n_steps`.

## 4. Altitude-field cache (reuse of Viz-0 E4 derivation)

**Finding**: [`render_altitude` in `c1_viz.rs`](../../../crates/ymir-viz/src/visualization/c1_viz.rs#L227) already performs the Architecture C derivation but **discards the raw altitude `GridF32`** (it only emits RGBA). Steps 1–3 of that function:

```rust
let s_field   = Field2D::from_vec(nx, ny, snapshot.s.clone());
let age_field = Field2D::from_vec(nx, ny, snapshot.age.clone());
let plate_type_field = /* decode snapshot.plate_type Vec<u8> → PlateTypeField */;
let iso = compute_isostasy(&s_field, &IsostasyConfig::default());
let mut altitude = iso.heightmap;                       // GridF32, non-dim
apply_stein_stein_bathymetry(&mut altitude, &age_field, &plate_type_field,
                             &SteinSteinParams::default());
// altitude.get(i as i32, j as i32) → f32 non-dim altitude per cell
```

**E1 refactor**: extract `derive_altitude_field(snapshot) -> GridF32` (steps 1–3) so both `render_altitude` (RGBA path) and the hover cache (raw value path) share one source of truth. The hover readout reads `altitude.get(i, j)` — the **non-dim** value (verification value, W3 global). The cache is keyed by `snapshot.step`; invalidated when a newer snapshot arrives (W2). Because the cache derives altitude independent of the active field view, **hover works in ALL views** (W4), not just Altitude.

## Extension points for later stages

- **`C1VizState`** ([`c1_plugin.rs`](../../../crates/ymir-viz/src/visualization/c1_plugin.rs)) currently: `texture_handle, field, pending_spec, show_voronoi_boundaries, show_velocity_vectors, arrow_scale, last_signature`. E1 adds an altitude cache (`Option<(usize, GridF32)>`) + hover readout state; E3 adds the pipeline toggle + cadence sliders surface. Hypsometric scale and `ActivePipeline` are proposed as **separate resources** (cleaner than overloading `C1VizState`).
- **Worker / command shape (E2)**: gallery `RunBaseline` must stay untouched (W4). Cleanest is a **new** `C1Command::RunWorkflow { spec, phase_a_params: PhaseAParams }` parallel entry (matches the `init-re-architecture-pattern` memory: ship re-architectured flow as a new parallel entry, preserve the old verbatim). The worker `Completed`-time state retention (E4) is a passive add — keep `Option<(C1State, PlateKinematics)>` after the run instead of dropping it.
- **`Started`/`total` counter (E2/E3)**: `poll_c1_events` currently sets `total = spec.n_steps`. Workflow per-step total is `n_cycles × k_cycle` (= 100 default), with `n_cycles` extra post-cycle snapshots. The counter must read the workflow total, not `n_steps`. Surfaced for E2 design (likely the `Started` event carries the resolved total, or the poll computes it from the pipeline + params).

## Open design choices to confirm before E1

1. **Hover readout placement** — fixed corner egui panel vs cursor-following tooltip. The sprite is a Bevy world entity (not an egui widget), so a cursor-following egui tooltip over the sprite is awkward (egui tooltips anchor to egui widgets; the pointer is *not* over an egui area when hovering the sprite). **Recommendation: fixed corner panel** (stable, no occlusion, trivial to implement). Confirm.
2. **E2 command shape** — new `RunWorkflow` command (recommended, preserves `RunBaseline` verbatim per init-re-architecture pattern) vs a `pipeline` field on a unified run command. Confirm.

## W7 surfaces resolved this stage

- Cursor→cell: existing `CursorWorldPos` (world space) + sprite-origin-centred inversion; Y-flip caveat noted; resize-correct via per-snapshot nx/ny.
- Split-borrow: confirmed compiling in `phase_a_c1.rs`; `params` is `&PhaseAParams` (issue-sketch type correction).
- Calibration anchor: `PhaseAParams::default()` = 5×20, alpha 0.01; total = n_cycles×k_cycle = 100.
- Altitude cache: factor `derive_altitude_field` out of `render_altitude`; key by step; available in all views.
