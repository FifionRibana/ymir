# C1 Phase 1.2 — Davis-Suppe orogenic closure outputs

Issue #123 visual + scalar outputs. Regenerate with:

```
cargo test -p ymir-core --test c1_phase_1_2_davis_suppe -- --nocapture
```

## Run configuration

- Grid: 64²
- Seed: 42
- Steps: 300 (forward Euler)
- Kinematics: hand-tuned 8-plate preset reused from Phase 1.1
  ([`crates/ymir-core/src/tectonics_c1/kinematics.rs::PlateKinematics::preset_phase_1_1`](../../../crates/ymir-core/src/tectonics_c1/kinematics.rs)).
- CFL: `Δt = 0.5 · dx / max|v| ≈ 0.69` non-dim/step (same as
  Phase 1.1).
- Davis-Suppe defaults (
  [`DavisSuppeParams::default`](../../../crates/ymir-core/src/tectonics_c1/closures/davis_suppe/source_term.rs)):
  - `coupling = 2.0`
  - `h_max = 2.5`
  - `L_taper = 4.0` cells
  - `L_decay = 6.0` cells
  - `max_distance = 30.0` cells (`= 5 · L_decay`)

## Snapshots

Each cycle dumps two PNGs:

| File | What it shows | Palette |
|---|---|---|
| `cycle_NNN_altitude.png` | Airy-isostasy heightmap of `S̃` through `tectonics::isostasy::compute_isostasy`. | **per-frame** auto-rescale (informational; the auto-rescale loses signal once boundary pile-up dominates the dynamic range) |
| `cycle_NNN_s.png` | Direct `S̃` heightmap with absolute palette `[0, 3.0]`. | **absolute**; saturating above 3.0 |

The fixed-palette `cycle_NNN_s.png` series is the
transport-correctness + closure-imprint check. **The visible
wedge formation** — green/brown patches in the upper-plate
interiors at the convergent edges — is the Davis-Suppe
signature; see `cycle_300_s.png` for the most readable end-state.

## Visual reading guide

- **cycle 000** — initial state: a single coherent Voronoï
  landmass (brown on dark blue ocean), 1-cell smoothstep at
  inter-plate boundaries (same start state as Phase 1.1).
- **cycle 050** — first wedges visible: bright green/brown
  patches appear on the upper-plate side of convergent
  boundaries, asymmetric vs the lower plate. Some saturation
  starts at convergence corners.
- **cycle 100** — wedges thicken; boundary corners begin to
  saturate the palette (white in `cycle_NNN_s.png`).
- **cycle 200 / 300** — steady-state-ish: wedge body on 5/8
  upper plates settles into the Davis-Suppe `h_critical(d)`
  fill profile (see below); boundary cells continue to
  accumulate via advection (no Phase 1.2 sink — see "On
  boundary cell brightness" section).

## Acceptance invariants (4 / 4 PASS)

The integration test
[`crates/ymir-core/tests/c1_phase_1_2_davis_suppe.rs`](../../../crates/ymir-core/tests/c1_phase_1_2_davis_suppe.rs)
asserts the four Phase 1.2 invariants below.

| Invariant | Acceptance | Empirical |
|---|---|---|
| (a) `wedge_p95 < 1.5 · h_max` | `< 3.75` | `0.376` ✓ |
| (b) `wedge_p99 > 1.0` | `> 1.0` | `5.83` ✓ |
| (c1) `fill_near (d∈0-5) > 0.5` | `> 0.5` | `0.778` ✓ |
| (c2) `mean(d∈0-5) / mean(d∈10-20) > 1.5` | `> 1.5` | `4.66` ✓ |
| (d) disabled-closure regression matches Phase 1.1 unbounded | `> 100` | `1079.7` ✓ |

Wedge-body filter: cells with `0 < wedge_distance < max_distance`
(intra-plate Dijkstra from upper-plate seeds, see Stage 3.1).
**Boundary cells** (Convergent type, `d = 0`) are excluded from
(a)/(b)/(c) by design — they are advection sinks at this stage;
Phase 1.4 erosion will balance them.

## Final wedge-body distribution

```
wedge body cells              = 2446 (59.7 %)
wedge S̃ min                  = 0.0006
wedge S̃ mean                 = 0.5598
wedge S̃ median               = 0.106
wedge S̃ p95                  = 0.376
wedge S̃ p99                  = 5.83
wedge S̃ max                  = 93.08
global_max (boundary pile-up) = 2297       (no Phase 1.2 bound)
```

## `h_critical(d)` profile imprint (per-bucket means)

The hallmark Phase 1.2 finding: Davis-Suppe imprints its
`h_critical(d)` physics on the **conditional** `S̃(d)`
distribution, but the **fill ratio** (mean / target) — not the
absolute mean — is what reveals the imprint in the
advection-dominated Phase 1.1 kinematics regime.

| Distance bucket | count | mean `S̃` | `h_crit` at mid | fill ratio |
|---|---|---|---|---|
| `d ∈ (0, 5]` | 687 | 0.904 | 1.162 | **0.778** |
| `d ∈ (5, 10]` | 764 | 0.752 | 2.117 | 0.355 |
| `d ∈ (10, 20]` | 843 | 0.194 | 2.441 | **0.079** |

→ Near-boundary cells reach 78 % of their h_crit target. Far
cells reach 8 %. The closure is alive everywhere but
spatially saturated only near boundaries.

## Advection-dominated regime finding

Phase 1.2 reveals that the Phase 1.1 hand-tuned kinematics
preset produces an **advection-dominated regime**: advection
rate `≈ 32 × source relaxation rate`. The Davis-Suppe closure
correctly imprints its `h_critical(d)` profile on the
conditional distribution `S̃(d)` (verified by invariants (c1)
and (c2)), but the **bulk wedge body drains continuously** as
advection moves cells through the wedge zone faster than the
source can fill them.

**Visual signature** (`cycle_300_s.png`):

- Sparse, asymmetric green/brown wedges on upper-plate sides
  of convergent boundaries — these are the cells where source
  saturation balances advective outflow.
- Bright "trap" cells along the boundary itself — these are
  the Convergent cells skipped by the closure (advection-only
  accumulators).
- Drained blue (ocean colour) covers the bulk of the
  upper-plate interior — source weak (envelope `≈ 0.1` at
  `d ≈ L_decay = 6` cells), advection drains faster than source
  fills.

This is **intrinsic to Phase 1.2 design** (Phase 1.1
kinematics + Davis-Suppe source only, no sink). Phase 1.3
equilibrium-height closure will face the same regime. Phase 1.4
erosion + isostasy will introduce a global mass sink, which
**may invert the spatial signature** (`mean(d∈10-20) >
mean(d∈0-5)`) once the source-vs-sink balance shifts. Phase 2's
constrained-kinematics sampling (§6.3 of the C1 design doc)
may also produce slower kinematics that naturally shift the
balance.

The acceptance criterion (c2) — `mean(d∈0-5) / mean(d∈10-20)
> 1.5` — encodes the advection-dominated direction. Future
phase tests should **re-evaluate the direction** rather than
copy this assertion verbatim.

## On boundary cell brightness

The Convergent boundary cells (`d = 0` by construction of the
Stage 3.1 intra-plate Dijkstra) are **not** bounded by
Davis-Suppe. The architectural skip is intentional:

> `h_critical(0) = h_max · (1 − exp(0)) = 0`. Applying the
> source term at `d = 0` would compute `driving = h_crit − h ≈
> −1` and *thin* the boundary instead of thickening the
> upper-plate interior — anti-geological.

Surfaced before Stage 4 implementation (Issue #123, Stage 3.1
finding); locked by
[`source_term_skips_boundary_cells`](../../../crates/ymir-core/src/tectonics_c1/closures/davis_suppe/source_term.rs)
unit test.

Boundary cells therefore continue to accumulate mass via
advection (no Phase 1.2 sink). The `global_max` of `2297`
(`≈ 2.1 ×` the Phase 1.1 baseline of `1080`) reflects this
plus the extra mass Davis-Suppe deposited on the upper-plate
side. **Phase 1.4 erosion** is the planned mass sink.

## Comparison Phase 1.1 vs Phase 1.2

| Metric | Phase 1.1 | Phase 1.2 |
|---|---|---|
| Closure | none (advection only) | Davis-Suppe orogenic |
| Final wedge-body p95 | n/a (no wedge concept) | 0.376 |
| Final wedge-body p99 | n/a | 5.83 |
| Final global max | 1080 (boundary advection) | 2297 (boundary + source) |
| Final mean `S̃` | 0.557 (constant — mass-conservation) | 1.574 (+1.0 from source) |
| Mass conservation | exact (drift 1.6 × 10⁻¹⁴) | source-broken (intentional, W3) |
| Wall time (300 steps) | 245 ms | 324 ms (+32 %) |
| Visible wedges | 0 | 5 plates (asymmetric near boundaries) |
| Silent plates | n/a | plates 0, 6, 7 (no upper-plate seeds) |
| Boundary cell control | none | none (skip-by-design, Phase 1.4 sink) |

The 32 % overhead in wall time comes from the once-outside-loop
boundary classification and intra-plate Dijkstra (`~ 50 ms`
total at 64²) plus the per-step Davis-Suppe sweep (`~ 100 ms`
over 300 steps).

## What this output is **not**

Still not a plausible continent. Wedges visible but boundary
saturation noisy; no equilibrium height (Phase 1.3); no
erosion sink (Phase 1.4); no upscale + climate downstream
(Phase 1.4). Phase 1.2 validates the orogenic-source-term
mechanism; Phase 1.3 will bound the asymptotic plateau by
gravitational collapse; Phase 1.4 will close the mass budget
with the sink and produce the first end-to-end heightmap.
