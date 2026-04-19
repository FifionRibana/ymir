# Slab-Pull Reformulation — Phase 1 Diagnostic

Issue #75 — *Reformulate slab-pull as an auto-regulated operator term instead of RHS forcing.*
Branch: `epic/75-slab-pull-operator/phase1-diagnostic`. Date: 2026-04-19.

This report maps the current slab-pull implementation so Phase 2 can proceed
without further code exploration. **Read §6 first**: one load-bearing claim in
the issue is contradicted by the code.

---

## 1. Slab-pull computation — current location and flow

Slab-pull is **not** a per-cell traction added to the RHS at subduction
margins. It is a **per-plate scalar boost of `Plate.velocity`**, which is then
copied verbatim to every cell owned by that plate during `rebuild_traction`.
`T_plates` is therefore a piecewise-constant velocity field (one value per
plate), averaged to faces in `compute_rhs`. It is not a spike along the
margin.

### Files involved

- [crates/ymir-core/src/tectonics/boundaries.rs](crates/ymir-core/src/tectonics/boundaries.rs) — margin classification, subducted-mass accumulator, `apply_slab_pull`.
- [crates/ymir-core/src/tectonics/plates.rs](crates/ymir-core/src/tectonics/plates.rs) — `rebuild_traction` / `rebuild_traction_smooth` write `TractionField` from per-plate velocities.
- [crates/ymir-core/src/tectonics/solver/traction.rs](crates/ymir-core/src/tectonics/solver/traction.rs) — `TractionField { tx, ty }` (two `Field2D`s).
- [crates/ymir-core/src/tectonics/solver/stokes.rs](crates/ymir-core/src/tectonics/solver/stokes.rs) — `compute_rhs` adds averaged `T_plates` to the face RHS.
- [crates/ymir-core/src/tectonics/solver/tectonics.rs](crates/ymir-core/src/tectonics/solver/tectonics.rs) — orchestration of the macro step.

### Where slab-pull enters

- Magnitude formula — [boundaries.rs:346-371](crates/ymir-core/src/tectonics/boundaries.rs#L346-L371):
  ```rust
  let pull_magnitude = slab_pull_factor * plate.subducted_mass;
  plate.velocity.0 += (pull_magnitude * vx / v_mag) as f32;
  plate.velocity.1 += (pull_magnitude * vy / v_mag) as f32;
  // then cap at max_plate_velocity
  ```
  Scalar per plate, in the direction of current velocity, capped at 5.0.
- Mass accumulator — [boundaries.rs:325-342](crates/ymir-core/src/tectonics/boundaries.rs#L325-L342). Reads `q < 0` cells in `source_rate`, attributes them to the owning plate's `subducted_mass` with `subducted_mass += |q| · dt`.
- RHS assembly — [stokes.rs:188-195](crates/ymir-core/src/tectonics/solver/stokes.rs#L188-L195): `row_vx[i] = dpdx + 0.5*(plates.tx.get(pi,j)+plates.tx.get(i,j))` (and ty equivalent).
- Traction rebuild from plate velocities — [plates.rs:1063-1086](crates/ymir-core/src/tectonics/plates.rs#L1063-L1086) (and smoothed variant at [plates.rs:999-1054](crates/ymir-core/src/tectonics/plates.rs#L999-L1054)): `tx.set(i,j, plate.velocity.0 as f64)`.

### Call chain (entry → slab-pull)

```
run_tectonics                        tectonics.rs:139
  └─ execute_tectonic_pass           tectonics.rs:624
       ├─ solve_velocity_* (Picard/Newton)   // reads T_plates built from LAST step's plate velocities
       │    └─ compute_rhs           stokes.rs:148
       ├─ compute_boundary_sources   boundaries.rs:159   // produces source_rate + boundary_type
       ├─ accumulate_subducted_mass  boundaries.rs:325   // updates plate.subducted_mass
       ├─ apply_slab_pull            boundaries.rs:346   // boosts plate.velocity
       └─ rebuild_traction[_smooth]  plates.rs:1063/999  // rewrites T_plates from plate.velocity
```

Entry: [tectonics.rs:760-774](crates/ymir-core/src/tectonics/solver/tectonics.rs#L760-L774) (slab-pull block), [tectonics.rs:813-827](crates/ymir-core/src/tectonics/solver/tectonics.rs#L813-L827) (traction rebuild).

### Inputs to slab-pull

- `plate.subducted_mass` — f64 scalar per plate, accumulated from **current-step** `workspace.source_rate` at `dt_rates` ([tectonics.rs:761-768](crates/ymir-core/src/tectonics/solver/tectonics.rs#L761-L768)).
- `plate.velocity` — (f32, f32), carries across steps; direction of previous velocity dictates the slab-pull direction.
- `workspace.source_rate` — produced by `compute_boundary_sources` from **current-step** `grid.vx/vy` (post-solve values from this step) and `plate_ids`.
- Config: `slab_pull_factor` (default 0.05), `max_plate_velocity` (default 5.0), `slab_pull_enabled` (default true) — [boundaries.rs:89-95,126-128](crates/ymir-core/src/tectonics/boundaries.rs#L89-L128).

**Crucial**: `apply_slab_pull` does not scale with local slab thickness, subduction angle, or age. It is a **global per-plate scalar** with one regulating mechanism (the velocity cap).

---

## 2. Subduction margin identification — current method

### Detection — [boundaries.rs:159-304](crates/ymir-core/src/tectonics/boundaries.rs#L159-L304)

For each cell, walk the four periodic neighbors ([L181-186](crates/ymir-core/src/tectonics/boundaries.rs#L181-L186)). If a neighbor has a different `plate_id`, mark the cell as a boundary and compute relative velocity in the normal direction:

```rust
let v_rel = (vx_there - vx_here) * normal_x + (vy_there - vy_here) * normal_y;
convergence_sum += v_rel;           // negative = converging
```
— [boundaries.rs:218-229](crates/ymir-core/src/tectonics/boundaries.rs#L218-L229).

Classification ([L240-251](crates/ymir-core/src/tectonics/boundaries.rs#L240-L251)) uses `plate_type_from_mean_thickness` (>0.4 → Continental, else Oceanic, [L155-157](crates/ymir-core/src/tectonics/boundaries.rs#L155-L157)):

- `is_converging && (Oceanic,Continental)` → `Subduction`.
- `is_converging && (Oceanic,Oceanic)` → `OceanicSubduction`.
- `is_converging && (Continental,Continental)` → `ContinentalCollision`.
- Else → `Rift` or `None`.

### Persistence & direction

- **Recomputed every macro step**; cached only for the step in `workspace.boundary_field: Option<BoundaryField>` ([tectonics.rs:743-751](crates/ymir-core/src/tectonics/solver/tectonics.rs#L743-L751)).
- Direction of subduction is **implicit**. The normal between cell and neighbor is axis-aligned (±1 on one axis, [boundaries.rs:201-215](crates/ymir-core/src/tectonics/boundaries.rs#L201-L215)); subducting side is inferred *per cell* by plate type — an oceanic cell next to continental is the sink (`q < 0` via `subduction_rate * convergence_rate`, [L256-263](crates/ymir-core/src/tectonics/boundaries.rs#L256-L263)); the continental side is the arc (`q > 0` or 0 under recycling).
- No explicit normal vector is stored. `convergence_sum` is a **scalar sum over up to four neighbors**, not a vector.

### Spatial resolution

Boundary is a **single-cell-thick layer** on *each* plate's side: the cell facing a foreign plate-id neighbor. A subduction trench is therefore two parallel cell strips (oceanic sink + continental arc). `source_rate` is Gaussian-smoothed post-detection with `source_smoothing_sigma = 2.0` ([tectonics.rs:752-758](crates/ymir-core/src/tectonics/solver/tectonics.rs#L752-L758)), which dilates the effective support.

### Can γ_slab be derived from existing data?

Yes, with two additions:

1. Persist a **vector normal** per margin cell. Currently `compute_boundary_sources` computes axis-aligned normals per neighbor and discards them. The aggregated normal (sum of per-neighbor normals, normalized) is enough in 2D.
2. Produce a **decayed scalar field** on the subducting side. Requires one extra pass (BFS or a distance-transform Gaussian kernel) over the margin-classified cells to spread γ_slab 3-5 cells inward.

Both can be appended to `BoundaryField` without touching `StaggeredGrid` or upstream code.

---

## 3. Mass recycling — current coupling to slab-pull

- File: [crates/ymir-core/src/tectonics/recycling.rs](crates/ymir-core/src/tectonics/recycling.rs) — only data structures (`RecyclingConfig`, `RecyclingBuffer`). All logic lives in [tectonics.rs:847-909](crates/ymir-core/src/tectonics/solver/tectonics.rs#L847-L909).
- Subducted mass flux — [tectonics.rs:848-859](crates/ymir-core/src/tectonics/solver/tectonics.rs#L848-L859):
  ```rust
  for j in 0..ny { for i in 0..nx {
      let q = workspace.source_rate.get(i, j);
      if q < 0.0 { total_subducted += (-q) * dt_rates; }
  }}
  ```
  Reads `source_rate`, which is computed from **velocity × convergence-rate** ([boundaries.rs:238, 256-274](crates/ymir-core/src/tectonics/boundaries.rs#L238)) — i.e. it already uses `v · n̂` implicitly via `convergence_sum`. It does **not** read `T_plates`.

**Consequence**: moving slab-pull out of the RHS does not touch recycling. Recycling consumes `source_rate`, which will remain correct as long as `compute_boundary_sources` still sees plausible velocities after the operator reformulation. The accumulator ([boundaries.rs:325-342](crates/ymir-core/src/tectonics/boundaries.rs#L325-L342)) can stay or be repurposed to drive γ_slab magnitude.

**Mass-conservation tests**:

- [boundaries.rs:560-580](crates/ymir-core/src/tectonics/boundaries.rs#L560-L580) `sources_conserve_mass_approximately` — prints the relative imbalance but does not assert a bound. Informational, not strict.
- [tests/rectangular_simulation.rs](crates/ymir-core/tests/rectangular_simulation.rs) — full-pipeline determinism tests; not targeted at slab-pull.

No test asserts mass conservation *under* slab-pull.

---

## 4. Tests that will become obsolete

Direct tests (in [tectonics.rs tests module](crates/ymir-core/src/tectonics/solver/tectonics.rs)):

| Test | Lines | Asserts | Fate |
|---|---|---|---|
| `slab_pull_increases_plate_velocity` | [1648-1676](crates/ymir-core/src/tectonics/solver/tectonics.rs#L1648-L1676) | `apply_slab_pull` boosts `plate.velocity` along its direction. | **Delete** — `plate.velocity` no longer carries slab-pull. |
| `slab_pull_capped` | [1678-1701](crates/ymir-core/src/tectonics/solver/tectonics.rs#L1678-L1701) | Velocity capped at `max_plate_velocity`. | **Delete** or rewrite as "γ_slab magnitude is bounded" once auto-regulation replaces the cap. |

Indirect tests (boundary classification, GPE/traction RHS assembly): none grep-positive on `T_plates` / `slab`. `compute_rhs` tests are unit-level on GPE and operator symmetry ([stokes.rs:447+](crates/ymir-core/src/tectonics/solver/stokes.rs#L447)) — they pass `TractionField::zero` or a uniform/convergent traction and remain valid because operator symmetry is orthogonal to slab-pull plumbing.

Boundary-detection tests at [boundaries.rs:483-596](crates/ymir-core/src/tectonics/boundaries.rs#L483-L596) (`subduction_detected_at_convergent_ocean_continent`, `rift_detected_at_divergent_boundary`, `interior_cells_have_no_source`, `sources_conserve_mass_approximately`, `gaussian_blur_preserves_total`) are untouched by the refactor.

Integration: `cfl_retry_succeeds_on_standard_configuration` ([tectonics.rs:1703-1739](crates/ymir-core/src/tectonics/solver/tectonics.rs#L1703-L1739)) uses `TractionField::two_plates_convergent`, so it bypasses the plate-velocity traction path; safe.

---

## 5. Proposed `γ_slab(x, y)` and `n̂(x, y)` design

### γ_slab(x, y) — computation path

Pseudocode, executed once per macro step right after `compute_boundary_sources`:

```
margin_cells = { (i,j) : boundary_type[i,j] ∈ {Subduction, OceanicSubduction}
                          AND source_rate[i,j] < 0 }            // subducting side only
for cell in margin_cells:
    γ_seed[cell] = slab_pull_factor * |source_rate[cell]|       // local, convergence-driven
spread γ_seed inward on the subducting plate with Benioff decay (see below)
```

Using `|source_rate|` (which is `subduction_rate × convergence_rate`) makes
γ_slab locally proportional to `|v · n̂|`. Combined with the operator term
`γ_slab · (v · n̂) · n̂`, this yields a damping force whose magnitude grows
with convergence — **auto-regulated** and SPD.

### n̂(x, y) — computation path

```
for cell (i,j) with boundary_type ∈ subduction types:
    sum_n = Σ over foreign-plate neighbors of (normal_to_neighbor)
    n̂[i,j] = normalize(sum_n)       // zero if interior
```

Normals are already computed on-the-fly in [boundaries.rs:201-215](crates/ymir-core/src/tectonics/boundaries.rs#L201-L215) — they just need to be persisted instead of discarded. In 2D with four axis-aligned neighbors, this produces one of eight canonical directions at a corner cell, or an axis-aligned direction on a straight margin.

### Benioff decay profile

**Recommendation: 1-D exponential on the subducting plate, decay perpendicular to the margin, characteristic width L = 3 cells.**

```
γ_slab(cell) = γ_seed[nearest_margin] * exp(-d_cells / L)    for cells on subducting plate
             = 0                                              elsewhere
```

- **Exponential** (not Gaussian) because it concentrates mass at the trench and fades monotonically — matches the physical picture that slab-pull is strongest where the slab first bends down, with diminishing horizontal projection deeper. A Gaussian would peak off-trench. A smooth-step would cut sharply at the tail and reintroduce a stencil spike.
- **L = 3 cells**: width 3 gives ~95 % of γ inside 9 cells, matching the issue's 3-5 cell band; narrow enough not to reach across a small plate.
- **Direction: perpendicular to the margin** (i.e., along `-n̂` relative to the margin cell), restricted to cells with `plate_id == subducting_plate_id`. No vertical dip angle since the model is 2D horizontal.

Implementation: iterative BFS up to 9 cells from each margin cell, or a separable Gaussian-blur kernel of `source_rate<0` followed by a mask restricting to subducting plate IDs. The Gaussian variant is cheaper and reuses `gaussian_blur_f64` ([boundaries.rs:377-428](crates/ymir-core/src/tectonics/boundaries.rs#L377-L428)); replace with exponential kernel or accept the Gaussian shape for simplicity in v1.

### Storage

**Recommendation: extend `BoundaryField`** rather than `StaggeredGrid` or a stand-alone struct.

```rust
pub struct BoundaryField {
    pub boundary_type: Vec<BoundaryType>,
    pub source_rate: Field2D,
    pub n: usize,
    // added:
    pub gamma_slab: Field2D,   // cell-centered, nonzero only near subduction margins
    pub normal_x: Field2D,     // cell-centered, unit on margin, 0 elsewhere
    pub normal_y: Field2D,
}
```

Rationale:

- `StaggeredGrid` is a physical-state container (thickness, velocity, plastic strain). γ_slab and n̂ are **derived margin metadata**; they belong with `boundary_type` and `source_rate`.
- `BoundaryField` already lives in `workspace.boundary_field`, recomputed every macro step — same cadence γ_slab needs.
- Keeps the operator term `γ_slab · (v·n̂) · n̂` in `apply_stokes` as a parameter passed alongside `eta`, rather than coupled to the grid.

### Recomputation cadence

Once per macro step, in the same block as `compute_boundary_sources` ([tectonics.rs:742-758](crates/ymir-core/src/tectonics/solver/tectonics.rs#L742-L758)), *before* the velocity solve of the **next** step. It does not need to be recomputed inside Newton iterations — γ_slab is frozen during the nonlinear solve, exactly as η is frozen in Picard and `T_plates` is today.

---

## 6. Open questions and risks

### Contradictions with the issue description

- **The issue states that slab-pull is "injected into the RHS via the precomputed `T_plates` field" creating "spiked RHS components along subduction margins".** The code does not match this description. `T_plates` is built from *per-plate scalar velocities* ([plates.rs:1080-1081](crates/ymir-core/src/tectonics/plates.rs#L1080-L1081)); slab-pull only boosts those scalars. The RHS contribution at a cell is `0.5·(v_plate(left) + v_plate(right))` — the average of two constant values, not a spike. Any RHS spike at a margin must come from elsewhere — most likely the **GPE gradient** `-g·(Φ(i,j) - Φ(pi,j))/dx` ([stokes.rs:188](crates/ymir-core/src/tectonics/solver/stokes.rs#L188)), which jumps across thin-oceanic / thick-continental margins. **Phase 2 must decide whether the refactor still pays off** given that the observed conditioning problem may not be slab-pull-driven.
- Phase 2 needs to validate this by instrumenting the linear-solve residual: if moving slab-pull into the operator does not improve BiCGSTAB convergence, the real culprit is GPE gradient magnitude at the trench, and a different remedy applies (GPE smoothing, bigger `source_smoothing_sigma`, etc.).

### Design ambiguities

- **Sign of `n̂`**: inward (toward trench) vs outward (away from trench) changes whether `γ · (v·n̂)·n̂` damps or excites subduction velocity. Recommendation: `n̂` points from subducting cell toward overriding cell (toward trench); combined with γ > 0 and the negative-sign operator convention, this damps convergent motion — which is the auto-regulated behavior the issue asks for. **Confirm in Phase 2.**
- **γ_slab magnitude tuning**: if seeded from `|source_rate|` the magnitude tracks `slab_pull_factor × subduction_rate`. The product must be re-tuned; existing `slab_pull_factor = 0.05` was calibrated against the scalar velocity-boost formulation and is not directly transferable.
- **`apply_slab_pull` replacement vs co-existence**: the current accumulator `plate.subducted_mass` and the per-plate velocity boost have downstream effects (plate velocity feeds `rebuild_traction`, which still feeds GPE-competitive traction in the RHS). Cleanly removing both is preferable — partial coexistence would double-count the pull.

### Coupling risks

- **Mantle convection** ([tectonics.rs:829-845](crates/ymir-core/src/tectonics/solver/tectonics.rs#L829-L845)) also modifies `plate_ctx.traction` by adding `coupling · mantle_flow`. Unaffected by the refactor, but will need co-existence (γ_slab term is additive in the operator; mantle flow stays additive in the RHS).
- **`rebuild_traction_smooth`** ([plates.rs:999-1054](crates/ymir-core/src/tectonics/plates.rs#L999-L1054)) blends velocities across plate boundaries based on accumulated sub-pixel displacement. Already reduces T_plates discontinuities. Interaction with γ_slab-based auto-regulation should be documented — possibly one renders the other redundant.
- **Plastic yielding + cratonic multiplier** ([picard.rs](crates/ymir-core/src/tectonics/solver/picard.rs)): η drops at yielded margins. Adding a `γ_slab · n̂n̂ᵀ` stencil contribution that scales only with γ (no η dependence) may dominate over the visco-plastic operator at yielded cells — could be a good thing (stabilizes the near-singular diagonal, reinforces #50's Jacobi floor), or a regression surface for continental-collision runs. No test currently covers this.

### Regression surface not covered by tests

- Wilson cycle / rift-to-subduction transitions — γ_slab would turn on and off as `boundary_type` switches.
- Plate fragmentation at step % 10 ([tectonics.rs:789-799](crates/ymir-core/src/tectonics/solver/tectonics.rs#L789-L799)) — margin topology changes discretely; γ_slab will flicker. Consider temporal smoothing.
- Plate consumption (`apply_subduction_consumption`) reassigns cell ownership. γ_slab recomputation after reassignment must be clean — ordering matters.
