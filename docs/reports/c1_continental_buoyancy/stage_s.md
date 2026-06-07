# Issue #145 — Stage S audit (surface-before-implement)

Branch `145-c1-continental-buoyancy` off `origin/milestone/c1-lightweight-dynamic-tectonics`. Maps the v→advection wiring, the rigidity-mask nature, the Track-D/boundary interaction, and determinism — BEFORE the binary-vs-graded choice (point of design 1).

## S1 — Where velocity enters advection; where the fix attaches

**Single source.** `fill_velocity_field` ([time_loop.rs:161](../../../crates/ymir-core/src/tectonics_c1/time_loop.rs#L161)) builds per-cell `vx`/`vy` from `kinematics.velocities[plate_id[cell]]`. Called once for the whole run on the Phase-1.x / Track-A/B path (line 429), and **per step** on the Track-D path (line 505, because `plate_id` may have mutated). `step_upwind` ([advection.rs:52](../../../crates/ymir-core/src/tectonics_v2/advection.rs#L52)) is the only consumer of `vx`/`vy`; it computes **face** fluxes using the **neighbour cell's** velocity (lines 82–90).

**Two candidate hooks — NOT equivalent at boundaries:**

| hook | behaviour | mass |
|---|---|---|
| **(a) upstream v-masking** (zero/scale `vx,vy` where rigid, in/after `fill_velocity_field`) | rigid cell self-velocity 0, but `step_upwind` still fluxes across the rigid/oceanic face using the **oceanic neighbour's** velocity → ocean exchanges with the rigid margin (the boundary caveat) | **CONSERVED** (face fluxes still cancel; prototype Δmass −0.00%) |
| (b) in-advection skip (don't update marked cells in `step_upwind`) | perfect freeze, inbound flux discarded | **BREAKS** conservation (matches the craton-frozen +125% blow-up) |

→ **Use (a) upstream v-masking.** It is the proven prototype (92% craton area, Δmass −0.00%), a single hook, and preserves the conservation property `step_upwind` exists for. The boundary exchange is physical, routed to subduction (S3). Binary (`v=0`) and graded (`v*=f(S̃)`) BOTH attach here — same hook.

## S2 — Rigidity mask + evolution

- **Mark via `plate_type == Continental`** (broadest; the validated prototype basis → 92%). `cratonic_mask` (`BoolField`, init-time) is a strict **subset** — a "cores-only rigid" fallback, but it protects only cratons and leaves non-cratonic continental crust advecting.
- **Gallery / production:** `plate_type` is STATIC (no reclassify) → rigidity mask STATIC, built once alongside `fill_velocity_field`. Cheap, coherent.
- **Track D:** `plate_type` MUTATES per step (subduction promotes Oceanic→Continental on floor-trigger; accretion). `cratonic_mask` does NOT mutate (Track D never writes it; only erosion reads it). So if mask = `plate_type==Continental`, it must be **recomputed per step on the Track-D path** — which is free, alongside the existing per-step boundary/velocity rebuild (lines 481–506). Newly-promoted continental crust becomes rigid the next step → **point-of-design 3 satisfied by construction.**

→ **Rigidity = `plate_type==Continental`, recomputed where `plate_type` is rebuilt** (once in gallery, per-step in Track D). Auto-tracks accretion/rifting evolution.

## S3 — Track-D boundary routing (the caveat)

Subduction ([subduction/source_term.rs](../../../crates/ymir-core/src/tectonics_c1/closures/subduction/source_term.rs)) already reads `plate_type` (`&mut`, promotes Oceanic→Continental) + `kinematics` (relative velocity onto the outward normal) and **consumes oceanic crust converging on continental** — exactly the inflow that v-masking produces at a rigid margin.

**KEY INSIGHT:** the prototype's boundary accumulation (edge mean 0.60→1.99) was measured with **closures OFF** — no subduction to consume the inflow. In PRODUCTION (subduction enabled), advection runs first (step 1) pushing ocean toward the rigid margin, then subduction consumes/promotes it. **The caveat is largely self-resolving once subduction runs.** To confirm: measure the rigid fix WITH Track D enabled and check the edge accumulation drops. This is the most delicate integration — flag for careful re-validation.

## S4 — Determinism

`step_upwind` and `fill_velocity_field` are **serial nested loops — no rayon / par_iter / sort** (grep confirms). Deterministic by construction. The rigidity mask derives from `plate_type` / `cratonic_mask` (deterministic init + deterministic Track-D mutation). v-masking is a serial element-wise scale on `vx`/`vy` — **no new nondeterministic ordering.** C1 `Deterministic` invariant preserved. Risk: minimal.

## Point-of-design 1 — informed framing (binary vs graded)

Both insert at the SAME hook (S1a) and are EQUALLY deterministic (S4) — graded is **not** more architecturally risky, only a parameter choice (`f(S̃)` ramp vs hard 0).

- **Binary** (`v=0` on continental): proven 92%, simplest. Hard margin discontinuity = the caveat (largely handled by subduction, S3).
- **Graded** (`v *= f(S̃)`, e.g. `f = clamp((S̃−0.2)/(1.0−0.2), 0, 1)`): softer margin transition (ramp, not wall), more physical (a margin isn't a wall). Same hook, same determinism.

Decision input = this audit + the comparative measurement (RUN 1 config, closures off, seed 42, SPATIAL bar): does graded cut the boundary accumulation while keeping ≳90% craton area? (And does binary's caveat self-resolve with subduction on, per S3?)
