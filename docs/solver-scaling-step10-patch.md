# Step 10 — §4.11 amendments and clarifications

This document collects formal patches to `solver-scaling.md` §4.11
produced during Step 10 implementation. They will be folded into
`solver-scaling.md` itself when the Step 10 PR merges.

## Patch 1 — non-conservative advection scheme for `A`

§4.11 specifies the age-field equation in its **non-conservative**
(material-derivative) form:

```text
∂A/∂t + ṽ · ∇A = Γ
```

with the source `Γ` reset at boundary cells (ridge, arc,
collision) and `Γ = 1` (Lagrangian growth) elsewhere. The
"Transport scheme" paragraph then states "Same upwind scheme as
S advection (§4.6), with the same CFL limit", which the Step 10
issue (D1) re-states as "OR call the existing function on A as a
separate scalar after S̃ advection is complete. The latter is
simpler [...] and is the recommended path."

The literal recommendation conflates **flux-form** (conservative)
upwind with **non-conservative** upwind. The S̃ advection in
`tectonics_v2/advection::step_upwind` is conservative —
`∂_t S̃ + ∇·(S̃ ṽ) = 0` — because S̃ is a thickness (mass-like
density). Reusing that function on `A` would expand to

```text
∂A/∂t + ṽ · ∇A + A · ∇·ṽ = 0
```

introducing a spurious `A · ∇·ṽ` source whenever the velocity
has non-zero divergence. With mantle forcing (Step 8+) making
`∇·ṽ` regularly `O(1)` and `A` ranging up to `O(10)`, the
spurious term would add `O(10) · ∇·ṽ · dt` per cell per step —
dominant over the intended quiescent growth of `1 · dt` and
visibly wrong (`A` would track `S̃ · v` divergence rather than
geological history).

> **Patch text — non-conservative scheme.** The age-field
> advection uses a **separate** non-conservative first-order
> upwind scheme implemented in `tectonics_v2::age_field::advection
> ::step_age_advect`. For each cell `(i, j)`:
>
> ```text
> vx_c = ½(vx[i,j] + vx[ip,j])
> vy_c = ½(vy[i,j] + vy[i,jp])
> dA/dx = (A[i,j] - A[im,j]) / dx   if vx_c ≥ 0, else (A[ip,j] - A[i,j]) / dx
> dA/dy = (A[i,j] - A[i,jm]) / dy   if vy_c ≥ 0, else (A[i,jp] - A[i,j]) / dy
> A_next[i,j] = A[i,j] - dt · (vx_c · dA/dx + vy_c · dA/dy) + dt · 1.0
> ```
>
> The cell-centred velocity averages the staggered face velocities;
> the upwind one-sided differences pick the cell from which information
> is *coming*. The CFL bound is the same as the conservative scheme
> (`dt ≤ cfl_factor · min(dx, dy) / max|v|`); the existing
> `tectonics_v2::advection::cfl_dt` helper is reused.
>
> The uniform `+dt · 1.0` source applies the §4.11 quiescent
> growth `dA/dt = 1` (Lagrangian frame) to **every** cell each
> step. Boundary-event resets (ridge / arc / collision) overwrite
> specific cells *after* the advection step, so the quiescent
> growth applied at boundary cells is intentionally discarded by
> the event reset (the per-cell semantics is "either the cell
> experienced a boundary event (overwrite) or it did not (keep
> advected value with quiescent growth)").

## Patch 2 — arc-cell detection (continuation of D3 wording)

§4.11 / D3 lists "Volcanic arc resurfacing. Cell age reset to
`A = 0`" as an event type, but `BoundaryFlag` does not carry an
explicit `Arc` variant — the arc semantics is implicit in the
`Q_arc` source-term computation, which fires on continental cells
adjacent to subducting cells.

> **Patch text — arc-cell detection.** Arc cells are detected by
> mirroring the `Q_arc` computation in
> `boundaries::source_sink::compute_source_sink_terms` (pass 2):
> a cell `(i, j)` is an arc cell if
>
> 1. it is itself continental (per `plate_type[i, j] == Continental`,
>    or per `S̃[i, j] > 0.5` — the union to be robust to S̃-
>    driven plate-type drift over the run);
> 2. at least one of its 4 periodic neighbours has
>    `BoundaryFlag::Subduction` or `BoundaryFlag::OceanicSubduction`.
>
> When a cell is flagged both `Rift` and arc-eligible (rare edge
> case), the explicit `Rift` flag wins the attribution — both
> reset `A = 0` so the value is unaffected, but the
> diagnostic counter increments `ridge_resets` rather than
> `arc_resets`.

## Patch 3 — observed age-field statistics (Step 10 baseline)

§4.11 leaves the per-region age statistics at "soft check"
status. The Step 10 baseline (64² × 100 steps, Step 8 shape,
defaults) records:

> **Documented baseline behaviour.**
>
> - `age_at_continental_cells_mean_final ≈ 1.64`
> - `age_at_oceanic_cells_mean_final ≈ 0.79`
>
> Both means are **below** the initial values
> (`continental_age_init = 7.0`, `oceanic_age_init = 0.5`) because
> the Step 8 active regime drives substantial advection and
> boundary-event activity, transporting young (post-reset) ages
> into the bulk of the domain. The continental > oceanic
> ordering is preserved (`1.64 > 0.79`) — consistent with the
> §4.11 expectation that ridge resets fire frequently on oceanic
> cells while continental cells are reset only at arc / collision
> events (rarer).
>
> Event counts (run total, 4096 cells × 100 steps = 409,600
> cell-steps):
> - `ridge_resets_total ≈ 209,905` (~ 51 % cell-step rate)
> - `arc_resets_total ≈ 27,696`    (~ 6.8 %)
> - `collision_max_events_total ≈ 55,058` (~ 13.4 %)
> - `collision_max_age_mean ≈ 6.49` (close to
>   `continental_age_init = 7.0`, confirms the max-of-protolith
>   semantics is correctly picking up older neighbour ages)

The three numbers above are recorded for the baseline regime
only. Other cratonic / mantle / Mf settings can shift these
counts substantially without invalidating the design — they are
documentary, not acceptance bounds.
