# Solver scaling — Step 11 patch (§4.12 plate kinematic drift)

**Status:** addendum to `docs/solver-scaling.md` introduced by Step 11
(issue #108). Adds §4.12 documenting the plate kinematic drift
mechanism and its validity envelope. Patches the original "initial
velocity per plate" framing in the issue (which was inconsistent with
the quasi-static Stokes formulation).

---

## §4.12 — Plate kinematic drift

### Why the original "initial velocity" framing is wrong

Step 11 was originally specified as an "initial velocity per plate"
mechanism: the user assigns `(vx, vy)` per plate at simulation start,
and the system "moves accordingly" thereafter.

This framing is not consistent with the quasi-static Stokes solver
used in the milestone. At every macro time step, the harness solves

```text
∇ · (η(ε̇, S, …) · ∇v) = forcing(S, mantle, GPE, slab, …)
```

This equation has no time derivative of `v`: the velocity field is an
**instantaneous response** to the right-hand side at the current `S`.
It is not a state variable that carries history. Any "initial `v`" is
overwritten in full by the first solve and has zero effect on the
subsequent trajectory beyond serving as a Newton warm start (which
affects iteration count, not the converged solution).

Empirically: with `MantleConfig::Disabled` and a prescribed
`v(t = 0) = 0.5` on one plate via the original Phase-3 wiring,
`peak|v|` at step 1 was `3 · 10⁻⁵` — exactly the Step-7-baseline
quiescent régime. The prescribed magnitude was completely lost. The
issue's acceptance #5 ("`peak|v|` at step 1 ≈ assigned magnitude")
could not be met by any wiring of an "initial velocity" into a
quasi-static solve.

### The mechanism actually shipped: kinematic drift with deformation/transport split

Step 11 ships **plate kinematic drift**: a per-plate velocity field
`v_drift` is constructed once at init from the user's `(vx, vy)`
assignment with smoothstep blending across inter-plate boundaries
(see [`tectonics_v2::plate_kinematic::field::build`]), and is then
**added back to `vx, vy` only inside the advection scope of each
time-loop iteration**:

```text
loop iter (vx, vy = solver-only at iter start):
    slab pre-solve, force sample, drag-diag rebuild   // sees solver-only
    snapshot, extrapolation, Newton solve             // sees solver-only
    StrainRate::compute, eta rebuild, yielding metrics // sees solver-only
    add drift:        vx += v_drift_x;  vy += v_drift_y    // ← v_total
    cfl_dt, step_upwind (S advection),                 // sees v_total
    age advection, boundary source/sink, slab.advect   // sees v_total
    progress callback (peek_vx = v_total)              // user-visible
    strip drift:      vx -= v_drift_x;  vy -= v_drift_y    // ← solver-only
end loop
```

The drift exists only between the add-before-advection hook and
the strip-at-iter-end hook. Outside that scope, `vx, vy = v_solver`
so the deformation pipeline (Newton's η, the post-solve strain-rate
diagnostics, the yielding metrics) operates on a clean velocity
field. This split is the central design decision of the patch.

#### Why deformation must see `v_solver` only

`v_drift` is per-plate uniform — its gradient inside a plate is
identically zero. At inter-plate boundaries, the smoothstep
transition over `boundary_smoothing_width` cells produces a
gradient `Δv_drift / (width · dx)` that, on a 32² grid, can easily
exceed the yielding threshold even at `drift = 0.001`. If
`StrainRate::compute(vx, vy)` sees `v_total`, it picks up that
spurious smoothstep gradient and the rheology fires yielding
*artificially* — driving `η ↓ → v_solver ↑ → ε̇ ↑` runaway through
S̃ advection feedback. Empirically the runaway hits
`vmax_peak ≈ 10²⁸` over 20 steps regardless of how small the drift
magnitude is.

The fix is conceptual, not a tuning: `v_drift` is a per-plate
**rigid transport** (a change of reference frame), not a
deformation. Deformation is the local strain rate of the
solver-balanced velocity. The two must be evaluated separately:

| Pipeline | Velocity field used |
| --- | --- |
| Newton solve (η, residual, line-search) | v_solver |
| `StrainRate::compute` post-solve | v_solver |
| `eta_cc` rebuild for diagnostics | v_solver |
| Yielding metrics (`peak_yielding_*`, …) | v_solver |
| `peak_v_solved`, `vmax_peak` reporting | v_solver |
| `cfl_dt` (transport stability) | v_total |
| `step_upwind` (S̃ advection) | v_total |
| `step_age_advect` | v_total |
| Boundary source/sink Q computation | v_total |
| `slab.advect` (m_subducted transport) | v_total |
| `prev_iter_start_v` warm-start snapshot | v_solver |
| `peek_vx, peek_vy` in progress callback | v_total (user UI) |

### What this breaks and what survives

**Breaks:**

- **Strict momentum conservation.** Adding `v_drift` post-solve injects
  kinematic energy that is not balanced by a body force in the
  momentum equation. This is the price of the mechanism.
- **The interpretation of `peak|v|` as "Stokes solution norm".** It is
  now `peak|v_solver + v_drift|`. The two contributions need to be
  separated for any analysis that previously assumed `peak|v|` came
  from the solver alone (e.g., `peak_v_damping_ratio` in older
  reports; the new acceptance criteria #5/#6/#7 already account for
  this).

**Survives:**

- **Steps 0-10 regression bit-identity.** With
  `PlateKinematicConfig::Zero` (the default), the strip and add
  hooks branch on `is_zero()` and skip the per-cell loops; the harness
  takes the legacy zero-init path bit-for-bit. Acceptance #1, #11,
  #12 all hold.
- **Mantle, slab, age field, cratonic compatibility "by construction"
  for Zero.** None of those mechanisms see any change when drift is
  off.
- **Newton convergence and CG conditioning** when drift is on, *if*
  the drift stays inside the validity envelope below.

### Validity envelope

The mechanism is a small-perturbation forcing around the
quasi-static Stokes solution. Two régime boundaries determine where
it stays well-defined:

#### Cumulative displacement

The drift advects `S̃` at every step. The cumulative drift-driven
displacement after `N` steps is
`d_cum = drift · N · dt = drift · total_time_nondim`. To keep the
tessellation recognisable (plates have not periodically wrapped or
collapsed onto themselves) the displacement must stay sub-grid-length:

```text
drift · total_time_nondim  ≤  ~ 0.5      (grid widths)
```

For the milestone default `total_time_nondim = 6.0`, this gives
`|drift| ≤ ~0.083`. The Phase-4c `motion_without_mantle` test runs
at `drift = 0.5` over 30 steps (`total_time = 1.8`), giving
`d_cum = 0.9` grid widths — slightly over the envelope, observable
as significant `S̃` redistribution, but still bounded (no blowup
because that test runs with yielding *Disabled*, see next régime
boundary).

#### Strain-rate threshold (yielding ON)

When yielding is Enabled, the rheology drops `η_eff` sharply where
`ε̇_II > Bi · S / η_visc`. The drift creates a sustained shear at
inter-plate boundaries (the smoothstep transition in
`field::build`):

```text
ε̇_drift  ≈  Δv_drift / (boundary_smoothing_width · dx)
         =  21 · drift / boundary_smoothing_width  (on a 32² grid
                                                    with dx = 1/32)
```

To stay sub-yielding, the operator-level threshold demands
`ε̇_drift < Bi · S_typ / η_visc ≈ Bi`. Combined with the formula
above:

```text
drift / boundary_smoothing_width  <  Bi · dx
                                  ≈  Bi · domain_lx / nx
```

For the default `Bi = 0.15`, `width = 1.5 cells`, `nx = 32`:

```text
drift  <  0.15 · 1.5 / 21  ≈  0.011
```

**Note: this strain-rate threshold applied to the *original*
Phase-4b wiring** that fed `v_total` to `StrainRate::compute` and
to the rheology metrics. After the deformation/transport split
(this patch), `ε̇_II` measures `v_solver` only, so the drift
gradient at boundaries does not contribute to yielding metrics.
The test `with_cratonic` (yielding ON + cratonic ON + drift > 0)
is bounded after the fix: `vmax_peak = 3.6e-5` matches the Zero
baseline exactly, `peak_yielding_in_craton = 0.0`, variance of
`v_solver - 0` inside cratons is `9.8e-7`. The
`drift / width < Bi · dx` bound is preserved as a *user-facing*
note about how big a drift the user can dial in before the
*physical* dynamics (post-advection S̃ accumulation feeding the
next solve) gets the rheology near its threshold — but it is no
longer a hard validity boundary of the mechanism itself.

### Recommended user defaults

UI presets seeded with non-zero drift should clamp magnitudes to
`|drift| ≤ 0.5` (per the issue D1 range) and warn if
`drift · total_time_nondim > 1.0` (cumulative displacement exceeds
one grid width). The Phase-5 panel will surface a small
"cumulative shift" indicator next to each per-plate slider so the
user can see when they are approaching the validity boundary.

### Why we ship this rather than wait

A more rigorous alternative (Option C in the Phase-4 remontée)
would be to add a true momentum equation with a drag-anchored
target velocity per plate, recalibrate the rheology around it, and
re-baseline Steps 0-10 — a several-week effort, deferred. The
shipped mechanism is a controllable forcing knob orthogonal to
mantle convection, which is sufficient for the workflow scenarios
(convergence / divergence / shear / triple junction) Step 12 will
chain together. The validity envelope above is its scope of use.
