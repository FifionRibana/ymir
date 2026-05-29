# C1 — Lightweight dynamic tectonics with empirical closures

**Status:** design document, not yet implemented. Successor architecture
to `tectonics_v2/` (the solver-reconstruction milestone). Issued at the
close of Step 12 once empirical evidence accumulated that no parameter
sweep within the v2 envelope unlocks the visual sweet spot Living Landz
needs.

**Scope:** Phase 1 (tectonics) only. Phases 2–6 (isostasy, upscale+FBM,
hydraulic erosion, climate, export) of the existing pipeline are
preserved and consume C1's output through the same interfaces v2 used.

**Companion documents:**
- `docs/solver-scaling.md` — physics reference for the v1/v2 thin-sheet
  formulation. Most of §1–§3 (scales, root-cause analysis) remains
  relevant to C1; §4 (per-mechanism formulations) is largely retired,
  replaced by the empirical closures listed in §5 of this document.
- `docs/solver-reconstruction-roadmap.md` — tracker for the v2
  step-by-step rebuild that this document closes.
- `docs/reports/step12_solver_audit.md` — the diagnostic that surfaced
  the structural mismatch and motivated this redesign.

---

## 1. Problem statement

### 1.1 What v2 delivered

The solver-reconstruction milestone (Steps 0–13.5) shipped a fully
non-dimensionalised thin-sheet solver with power-law rheology, plastic
yielding, basal drag, boundary sources/sinks, dynamic Voronoi
boundaries with conservative recycling, slab pull, mantle forcing,
cratonic immunity, geological age field, plate kinematic drift, and
continental + oceanic FBM heterogeneity. The non-dimensional core is
mathematically sound: discretisation MMS slopes pass at order 2,
Jacobian symmetry holds to 1e-14, mass conservation residuals sit at
machine noise.

### 1.2 What v2 did not deliver

Three pain points surfaced empirically and converge on the same
structural conclusion:

**P1 — Runtime is incompatible with the use case.** A 32² run takes
~35 minutes; 128² is hours; the 512² target is unreachable. The
Step 12 R5b D0 audit (`docs/reports/step12_solver_audit.md`) traced
this to CG cap saturation at the `mf = 1.0` nominal mantle regime
(`κ ~ 1.3-2.1 × 10⁴`), masked by every prior step's regression run
using mantle-off or weak-mantle baselines. The R5b D2 + D1-ter fixes
recovered ~85% runtime but the residual is still above target by
roughly two orders of magnitude.

**P2 — Visual output is inadequate.** Three independent parameter
sweeps (`mf` alone at R5b, `mf × evolution_rate` at R6.3, `cratonic_amp`
at C.4) reproduced the same binary pattern: preservation OR dynamics,
never both. The continental boundaries inherited from the rectilinear
Voronoi tessellation produce orogenic chains aligned along straight
edges — visually implausible at the Living Landz consumption scale.
Step 12 R3's reformulation of acceptance #6 ("border curvature is
out of Phase A scope") was not a calibration shortcut; it was the
correct framing of a structural limit of D2's averaging-by-construction.

**P3 — Class-of-problem mismatch.** The published thin-sheet
literature (England & McKenzie 1982, Houseman & England 1986, Flesch
et al. 2001) solves for **observed** geophysical phenomena (Tibet,
Andes, Basin & Range): the success criterion is reproducing measured
GPS velocities and observed topography within error bars. Ymir's
success criterion is different — produce plausible, legible
continental morphology at game-relevant scales. Thin-sheet captures
the former and emerges into the latter or does not. For Ymir it does
not, and no amount of solver engineering changes the class of
problem being solved.

### 1.3 Why this is a redesign, not a continuation

The three pain points are linked. P1 is a property of the equation
class (non-linear coupled Stokes with high viscosity contrast). P2 is
a property of the initial condition geometry (Voronoi rectilinear
boundaries inherited unchanged through v1 → v2). P3 says these are
not solvable by re-tuning either parameters or solver — the architecture
itself does not target what Ymir needs.

C1 is a different architecture, not a Step 14 of v2. The v2 module
tree (`tectonics_v2/stokes/`, `tectonics_v2/forcing/slab_pull.rs`,
`tectonics_v2/mantle/`, etc.) is retired. The v2 infrastructure that
is **not** Stokes-specific (advection, grid, init geometry, boundary
classification, age field, diagnostics framework, Bevy bridge,
hydraulic erosion HD, isostasy) is preserved and reused.

---

## 2. Use case constraints driving the redesign

These are explicit because they were the basis for choosing C1 over
the alternatives evaluated in §3.

### 2.1 Living Landz product scope

The generator produces continents for a medieval sandbox game. Output
is consumed at three scales: continental overview, regional (kingdom)
scale, local (hex cell ~30–50 m radius). The same heightmap drives
village placement, cliff/beach detection on coastlines, river network
navigability classification, biome assignment, and mineral resource
distribution. There is no separate authoring step — what the generator
produces is what the game consumes, after ×100 upscale.

### 2.2 Required morphological diversity

A single continent must contain:

- Coastal plains broad enough to host cités (1.5–2 km+ valley floors).
- Multiple mountain systems of **different ages**: young high ranges
  near active margins, older eroded massifs in the continental interior.
  Single-orogeny continents read as flat regardless of internal
  parameter variation.
- Variable coastal morphology: beaches, cliffs, fjord-like indentations
  in some regions, smooth coastlines in others. The variation must be
  emergent — a uniform FBM post-processing produces the characteristic
  "fake" coastline that Living Landz already rejects in the Azgaar
  output.
- A complete fluvial network linking interior to coast, which is
  already handled by the HD erosion stage downstream of C1.
- Distinguishable biomes emerging from orography × prevailing wind ×
  Whittaker classification (existing climate phase). Rain shadows
  require realistic orographic geometry.

### 2.3 Runtime envelope

Target: < 10 s per run at 512² on commodity desktop hardware. This is
not a soft target — it determines whether the tool is usable as a
design tool (interactive iteration on parameters) versus a batch
processor (set up, leave running, return). v1 was usable at 64² in
seconds but produced unusable visuals; v2 produces better-grounded
runs but is unusable at 32² in 35 minutes. C1's target sits in the
envelope v1 occupied, with the v2 lessons folded back in.

A secondary target: per-step display latency ~30–100 ms. This makes
the temporal evolution legible to a human watching the run, supports
mid-run abort if the trajectory is going wrong, and lets the user
form intuition about how parameters shape outcomes.

### 2.4 Domain shape

512² fully periodic torus, with user-controlled viewport scrolling
for visual centring. Continental tessellation must produce
geographically clustered continental plates so that **some** viewport
position shows a continent surrounded by ocean. Continental fraction
~29% (Earth-like) gives enough ocean for any reasonable continent to
be cadrable in a 512² window.

This design choice avoids the periodic-domain artifacts of v1/v2
(continent wraps across boundary) without paying the cost of a
finite-domain advection scheme. Seeds that fail to produce a cadrable
continent are rejected by the user.

---

## 3. Alternatives evaluated

Four architectures were considered. The evaluation was done at design
time before any implementation; the rationale is preserved here so
future work can re-evaluate the choice if assumptions change.

### 3.1 Option A — static with invented rules

Place mountain ranges, oceanic ridges, basins by geometric templates
modulated by procedural noise. No physical grounding, no temporal
dimension. Closest to what Azgaard and similar generators produce.

**Rejected because:** (a) cannot produce coherent chronological
diversity (young vs old mountain ranges on the same continent),
(b) parameter calibration drifts indefinitely with no physical anchor,
(c) "looks invented" is a real failure mode for users who can identify
imposed geometry by eye.

### 3.2 Option B — static with validated empirical closures

Same architecture as A but parameters anchored in geophysical
correlations (Davis-Suppe critical taper, Lallemand subduction-arc
relations, Parsons-Sclater oceanic bathymetry, etc.). No temporal
integration; the final morphology is computed in O(N) from the
plate kinematic configuration.

**Rejected as primary architecture because:** cannot produce
chronological diversity. Two mountain ranges of different apparent
age require either a chronological mechanism that produces them
sequentially (C1) or an explicit posing of "this range is old,
this one is young" (which reintroduces the painting metaphor B
explicitly rejected). The other limits of B are mitigatable by
adding closures, but the chronology limit is structural.

B is **partially adopted** as the closure source for C1's source
terms; see §5.

### 3.3 Option C1 — lightweight dynamic with validated closures (selected)

Advection + closed-form source terms for crustal thickness evolution,
under a prescribed plate kinematic field. No Stokes solver. No
non-linear rheology. No implicit coupling. Each closure (orogenic
profile, equilibrium height, oceanic bathymetry, macro erosion,
isostasy) is an independent additive term in the time-stepped
evolution of S̃. The temporal integration produces chronological
diversity emergently: mountain ranges formed at step 50 have
300-50 = 250 steps of weathering by the run's end; ranges formed
at step 250 have only 50 steps.

This is the system-code analogue (in the thermohydraulic-nuclear
sense) of the v2 thin-sheet solver: closures replace the resolution
of a coupled PDE, the cost per step drops from "Newton outer iters ×
CG inner iters × matvec" to "single pass advection + source
evaluation", and the morphological emergence still occurs because
the closures encode the validated physics of what they're closing.

### 3.4 Option C2 — repair v1

Reuse v1's thin-sheet linear solver, fix the three diagnosed defects
(mass leakage out of continental cells, no continental↔oceanic
reclassification, GPE-driven excessive spreading) without re-introducing
v2's full mechanism stack.

**Rejected because:** the diagnosed defects are precisely the
mechanisms v2 added (Step 6 conservative recycling, Step 9 cratonic
immunity, Step 3 non-linear yielding). Re-adding them lightly likely
re-introduces the saturation v2 demonstrated. C2 has the worst of
both worlds: v1's structural limits with v2's runtime regime.

C2 is **not foreclosed** — if C1 reveals an unexpected mechanism
limitation that C2 would address, the decision can be revisited.

### 3.5 Why C1 won

C1 dominates A and C2 on the dimensions that matter for Living
Landz, and dominates B on chronology while preserving B's
defensibility through shared empirical closures.

The key insight from the discussion: closures map naturally to
**source terms in a dynamic balance**, not to **geometric forms
posed statically**. This is the structural reason C1 is more
faithful to the system-code analogy than B. Dittus-Boelter is a
closure for a heat balance, not a temperature shape. The same logic
applies here: Davis-Suppe is a closure for an orogenic balance, not
a mountain shape — using it as a source term in C1 is closer to its
original semantic than using it as a static form in B.

---

## 4. Architecture

### 4.1 Pipeline overview

```text
[Plate init]            R7-generalised: clustered continental plates,
                        non-rectilinear boundary perturbation,
                        cratonic / non-cratonic classification per
                        cell.

[Kinematics]            Per-plate velocity assignment (translation
                        only; intra-plate strain emerges from
                        boundary forcing in §4.4).

[Time loop, N=300 steps]
    Advect S̃            Conservative upwind, periodic, O(N_cells).
    Advect age field    Same scheme.
    Apply closures      Sum of active terms (§5).
    Reclassify          Per-cell continental/oceanic from S̃ vs sea_level.
    Update boundaries   Plate fusion (accretion), consumption (subduction),
                        rifting if applicable.
    Optional: emit progress callback for UI streaming.

[Final state]           S̃, classification, age, plate id.

[Isostasy]              Existing v2 module — S̃ macro → altitude macro.

[Upscale + FBM]         Existing — macro heightmap → HD heightmap.

[HD erosion]            Existing — particle-based hydraulic erosion.

[Climate / biomes]      Existing — orography + wind + Whittaker.

[Export]                Existing.
```

The diagram makes explicit that **only the time-loop block is new**.
Everything upstream (plate init) is a generalisation of existing v2
work; everything downstream is preserved from v2 unchanged.

### 4.2 Domain and discretisation

512² fully periodic torus (target). Mise au point at 64²–128²
benefits from sub-second runtime. The discretisation is the same MAC
staggered grid v2 used for `S̃` at cell centres; velocity components
on cell faces. Periodicity is handled by `PeriodicIndex` from
`tectonics/solver/field.rs` (already re-exported via
`tectonics_v2::field`).

There is **no momentum equation**. Velocity is not a solved field;
it is a prescribed input derived from §4.4.

### 4.3 State

Per-cell:
- `S̃` (nondim crustal thickness)
- `plate_id` (which plate the cell belongs to; can change on accretion)
- `plate_type` (continental / oceanic / cratonic; derived from S̃ +
  history)
- `age` (geological age field, advected with S̃; reset on creation
  events at boundaries)

Per-plate:
- Translation velocity `(vx, vy)`
- Type classification (continental majority / oceanic majority,
  derived from cell histogram)
- Optional: euler pole for spherical-style rotation rather than pure
  translation (deferred — see §7)

### 4.4 Intra-plate velocity field (open question)

The simplest model is **constant velocity within each plate**: every
cell of plate `p` gets `v_p`. This produces velocity discontinuities
across plate boundaries, which is physically what creates the
convergence / divergence signal that drives source terms.

The discontinuity is a feature, not a defect — it concentrates the
deformation signal exactly where the closures need it. The
alternative (smoothed velocity across plate boundaries) blurs the
signal and may require parameter compensation.

**Decision deferred to prototype.** Start with constant-per-plate;
revisit if the source-term signal at boundaries is too sharp or too
noisy in practice. If smoothing is needed, a linear diffusion of
`v` near boundaries (small Laplacian, no rheology, no Newton) is
the natural extension and keeps the architecture in the C1 regime.

The cratonic case is special: cratonic cells **inherit the plate's
translation velocity exactly**, with no smoothing or modification.
The craton does not deform internally; it transports rigidly. This
is implemented as a per-cell mask, not as a viscosity contrast (the
v2 mechanism that drove κ up).

### 4.5 Boundary evolution

After each time step, plate boundaries are re-derived from `plate_id`
discontinuities. Three events can occur:

- **Subduction**: a cell with oceanic plate_id adjacent to a cell
  with continental plate_id, where the relative velocity is
  convergent. The oceanic cell's S̃ is consumed (transferred to the
  recycling budget); the continental cell's S̃ is incremented by
  the appropriate source term (orogenic + arc volcanism).
- **Accretion / continental collision**: two continental plates
  converging produce orogenic thickening on both sides of the
  boundary; the two plate_ids may eventually merge into one if the
  convergence is sustained.
- **Rifting**: divergent boundary inside a continental plate produces
  thinning + extension. If sustained, the plate splits into two
  plate_ids.

These mechanisms are inherited conceptually from v2 Step 6 but
implemented as **direct algorithmic updates to plate_id and S̃**
rather than as constraints inside a coupled solve. The Step 6
RecyclingBuffer pattern is preserved.

*Note: Phase 2 Track D (Issue #132) implements all three mechanisms
as parallel-entry closure modules in `tectonics_c1/closures/`:
`subduction/` (rate-based oceanic consumption + arc volcanism +
floor-triggered plate_id reassignment), `accretion/` (sustained-
convergence plate_id merge with mass-weighted velocity averaging,
no thickening — Davis-Suppe closure §5.1 already handles continental-
continental orogeny), and `rifting/` (dedicated thinning closure on
divergent continental boundaries + "chewing-gum cut" split mechanism
gated by both sustained-divergence time AND sub-threshold thinning).
Plate split events propagate `age = 0` along the new divergent
boundary via the Path 3.B event-driven extension of §6.5. See §5.2
for closure-table detail and §7.2 for delivery status.*

### 4.6 Time stepping

Forward Euler on the advection (`S̃_new = S̃ - Δt · ∇·(S̃·v) +
Δt · Σ closures`). The CFL condition is the only stability
constraint and is trivial to enforce because there is no implicit
coupling. `Δt = 0.5 · dx / max|v|` is conservative.

No Newton, no CG, no preconditioner. The cost per step is
dominated by the advection sweep + closure evaluations, all of
which are O(N_cells) with small constants.

### 4.7 Removed from v2

- `tectonics_v2/stokes/*` (all of it)
- `tectonics_v2/forcing/slab_pull.rs`
- `tectonics_v2/forcing/mantle_force.rs`
- `tectonics_v2/mantle/*`
- `tectonics_v2/rheology.rs` (the non-linear power-law closure)
- `tectonics_v2/basal_drag.rs` (in its v2 form; if a basal drag
  effect is needed in C1 it is reintroduced as a closure, not as
  an operator diagonal)
- All AMG, Picard, continuation, snapshot infrastructure

### 4.8 Preserved from v2

- `tectonics_v2/advection.rs` (the upwind scheme)
- `tectonics_v2/field.rs` and the `Field2D` / `PeriodicIndex` types
- `tectonics_v2/voronoi/*` (subject to R7 generalisation, see §6)
- `tectonics_v2/cratonic/*` (the BFS that builds the cratonic mask;
  the **factor** field is repurposed as a binary cratonic mask in
  C1, dropping the smoothstep amplification that v2 used)
- `tectonics_v2/age_field/*`
- `tectonics_v2/diagnostics/*` (instrumentation, reports, harness
  structure)
- `tectonics_v2/boundaries/{plate_type, boundary_flag, layouts}.rs`
- `tectonics_v2/recycling/*`
- `tectonics_v2/workflow/*` (the Step 12 Phase A / Phase B orchestrator;
  Phase A's tectonic sub-cycle becomes a C1 run, Phase B is unchanged)
- Everything in `crates/ymir-core/src/{erosion, terrain, climate,
  export}` (the downstream pipeline)
- The Bevy bridge plumbing (`crates/ymir-viz/src/bridge/`)

---

## 5. Closures

Each closure is an additive source term to `∂S̃/∂t`, evaluated
per cell from the current state. Closures are activatable
independently via UI toggle; see §6.3 for the rationale.

### 5.1 Minimum viable set (MVP)

These five form the initial implementation target. They are
sufficient to produce orogeny + ocean basins + chronological
weathering, which validates the architecture end to end. Other
closures are added incrementally after this MVP is visually
validated.

| Closure | Phenomenon | Reference |
|---|---|---|
| Orogenic profile | Davis-Suppe critical taper — asymmetric mountain profile from basal friction + internal friction | Davis, Suppe & Dahlen 1983, *JGR* 88(B2); review Dahlen 1990, *Annu. Rev. Earth Planet. Sci.* 18 |
| Equilibrium height | Gravitational collapse limit on plateau elevation; balance between compression and lateral spreading | Molnar & Lyon-Caen 1988, *GSA Special Paper* 218; England & Houseman 1989, *JGR* 94(B12) |
| Oceanic bathymetry | Thermal subsidence ~√age; depth = 2500 + 350·√age (m) for young crust | Parsons & Sclater 1977, *JGR* 82(5); revision Stein & Stein 1992, *Nature* 359 |
| Macro erosion | Stream-power law applied at the macro scale to age old ranges within the C1 time loop | Whipple & Tucker 1999, *JGR* 104(B8); Willett 1999 for orogen-erosion coupling |
| Isostasy | Airy local isostasy: altitude from S̃ and reference densities | Standard (Turcotte & Schubert *Geodynamics*); existing v2 `compute_isostasy` |

*Note: Phase 2 Track A (Issue #129) implements the oceanic-bathymetry
closure using Stein & Stein 1992 (the revision of Parsons-Sclater 1977
referenced above) under **Architecture C** — post-isostasy bathymetry
adjustment that replaces altitude on oceanic cells with the S-S depth
formula based on cell age, rather than entering as an additive source
term on `S̃`. Rationale and per-cell formula live in the
`tectonics_c1/closures/oceanic_bathymetry/` module docstring. If
Architecture C limitations surface during Stage D visual review,
fallback to additive-source Architecture A or hybrid Architecture B is
documented as a follow-up path.*

### 5.2 Second wave (post-MVP)

These are necessary for the morphological diversity Living Landz
requires but are not gating for architectural validation. Add
them once the MVP is visually validated.

| Closure | Phenomenon | Reference |
|---|---|---|
| Subduction arc | Volcanic arc position and amplitude from convergence velocity, slab age, dip | Lallemand, Heuret & Boutelier 2005, *G-cubed* 6(9); Syracuse & Abers 2006, *G-cubed* 7(5) |
| Rifting / passive margins | Continental thinning under divergent kinematics; passive margin morphology after rift maturity | McKenzie 1978, *EPSL* 40(1); Buck 1991, *JGR* 96(B12) |
| Foreland basin (flexure) | Lithospheric flexure under orogenic load creating foreland basin + forebulge | Beaumont 1981, *Geophys. J. R. Astron. Soc.* 65(2); standard treatment in Turcotte & Schubert |

These three matter for Living Landz specifically because the MVP
without them produces only one morphological class of mountain
(continent-continent collision) and one class of coast (compressive
margin). Andean arcs, Atlantic-style passive margins, and Po-Valley
foreland plains are all absent from the MVP output. Foreland basins
are particularly relevant because they produce fertile plains —
prime city-placement terrain.

*Note: Phase 2 Track D (Issue #132) implements the subduction and
rifting closures of the table above, plus a third **accretion**
mechanism. Track D is the first C1 work-track to mutate `plate_id`,
`plate_type`, and `kinematics` per-step — Phase 1.1–1.4 + Track A/B
all treated these as static-after-init. The three Track D mechanisms:*

- **Subduction (`closures/subduction/`)** — rate-based consumption
  `Δs = K_subduction · |v_convergence| · dt` on oceanic cells
  adjacent to convergent oceanic-continental boundaries. Consumed
  mass is redistributed as arc volcanism to the N nearest
  continental neighbours within `arc_distance` (local BFS). When
  the oceanic cell's S̃ drops below `plate_id_reassign_threshold`,
  the cell is reassigned to the adjacent continental plate (rest
  of S̃ contributes to arc). Foreland basin morphology (Beaumont
  flexure entry above) remains unimplemented in Track D — covered
  by §7.3 Phase 3.
- **Accretion (`closures/accretion/`)** — sustained-convergence
  merge: when two convergent plates remain convergent for at least
  `merge_time_threshold` steps (tracked by a per-pair
  `ConvergenceTracker`), the smaller-index plate absorbs the
  larger-index plate's cells. Post-merge velocity is the
  mass-weighted average of the two pre-merge velocities. **No
  thickening source** — the existing Davis-Suppe closure (§5.1)
  already produces orogenic morphology at convergent boundaries
  during the pre-merge phase; the merge itself only resolves the
  boundary topology.
- **Rifting (`closures/rifting/`)** — two-stage mechanism:
  (1) a thinning closure (negative `s_field` source) applied on
  divergent continental boundaries, mirroring the Davis-Suppe
  positive source on convergent ones; (2) a "chewing-gum cut"
  split mechanism gated by BOTH conditions: sustained-divergence
  time ≥ `split_time_threshold` AND boundary `S̃` < `split_thickness_threshold`.
  Either condition alone is insufficient — both must hold (per
  Phase 2 Track D Q3.2 hybrid-conditions decision). New plate_id
  allocated for the rifted-off cells; per-plate velocity inherits
  parent with a perpendicular offset. **Path 3.B event-driven
  age=0** propagation along the newly-spawned divergent boundary
  extends Track B's init-only Path 3.A to maintain the age-density
  pile-up mitigation across rift events.

*Mass conservation diagnostic (test-only): each subduction step
records `(consumed, arc_distributed)` deltas; the per-cycle
accumulator validates that `initial_total_mass - final_total_mass`
matches `(consumed - arc_distributed)` within `1e-6` tolerance.
Track D's other two mechanisms are exactly conservative by
construction (accretion mutates only `plate_id` / `kinematics`,
not S̃; rifting thinning is a closure source — its imprint shows
in the standard mass-budget envelope already validated for the
other closures).*

### 5.3 Optional enrichments

Lower priority; add only if specific morphologies are missing from
output.

| Closure | Phenomenon | Reference |
|---|---|---|
| Hotspot / island chains | Stationary mantle plume + moving plate produces linear volcanic island chain | Morgan 1971, *Nature* 230; Wilson 1963, *Can. J. Phys.* 41 |
| Terrane accretion | Oceanic-borne crustal fragments accreted to active margins, producing patchwork continental geology | Coney, Jones & Monger 1980, *Nature* 288 |
| Post-orogenic thermal relaxation | Long-timescale cooling-driven subsidence of old ranges, beyond simple stream-power erosion | Stüwe 2007 *Geodynamics of the Lithosphere*; no single canonical paper |

### 5.4 Inter-closure coherence

Closures interact indirectly through their shared target state
(`S̃`, `age`). Two practical consequences must be designed for:

**Calibration cross-validation.** Each closure was calibrated in its
own observational context; their joint behaviour is not guaranteed
to be quantitatively consistent. A global mass-balance diagnostic
(echo of v2 Step 6's conservation acceptance) is needed to catch
inter-closure drifts. Acceptance criterion: total continental mass
plus ocean recycling budget plus subducted accumulator stays
conserved to within 1e-6 per step.

**"Isolation" is partial.** Toggling closure X off does not measure
X's pure contribution — it measures the system without X relative
to the system with X. The other closures see a different state and
behave differently. This is the same ambiguity that caused the v2
Step 4 regression-target confusion. The UI documentation must make
this explicit.

---

## 6. Initial conditions

The init produces: a plate tessellation, per-plate kinematics, an
initial `S̃` field, an initial cratonic mask.

### 6.1 Plate tessellation — R7 generalised

The v2 Voronoi tessellation (`tectonics_v2/voronoi/`) produces
rectilinear plate boundaries. C1 needs non-rectilinear boundaries
to avoid the v1/v2 visual failure mode of orogenic chains aligned
along straight Voronoi edges.

Three plausible mechanisms, ordered by expected effort:

- **Boundary displacement.** Take the standard Voronoi, then perturb
  the boundary with a low-frequency noise field that pushes the
  boundary in/out by a controlled amplitude. Cheap, well-understood,
  predictable.
- **Multi-scale tessellation overlay.** Union of a coarse Voronoi
  (macro plate skeleton) with a fine Voronoi (boundary roughness).
  More expressive than displacement but harder to control.
- **Lloyd relaxation with stochastic seed perturbation.** Standard
  CVT relaxation with non-uniform target density (denser seeds
  where boundary detail is wanted). More principled but more code.

Start with boundary displacement. The Step 12.X follow-up issue
(noted in §4.14 of the v2 patch) covered this territory and can be
re-purposed.

*Note: Phase 2 Track B (Issue #131) implements the **boundary
displacement** option using Perlin/Simplex noise via the `noise`
crate. For each grid cell, a noise-sampled displacement vector
shifts the candidate sampling position before re-running the
nearest-Voronoï-seed query — produces curved boundaries while
preserving the seed-based plate identity. Defaults `amplitude =
grid_size / 8`, `frequency = 4.0`, `octaves = 3`, `persistence =
0.5`. Lloyd relaxation and multi-scale overlay remain deferred.*

### 6.2 Continental clustering

For the §2.4 viewport requirement, continental plates must be
**geographically clustered** on the torus — not uniformly scattered.
Two strategies:

- **Spatially biased seed sampling.** Sample plate seeds from a
  non-uniform distribution that concentrates them in one half of
  the torus.
- **Cluster-based type assignment.** Sample seeds uniformly, then
  classify continental/oceanic by graph connectivity: pick one seed,
  declare it continental, mark its neighbours continental too until
  the target continental fraction is reached. Remaining seeds are
  oceanic.

The second is simpler and produces more compact clusters. Start
there.

*Note: Phase 2 Track B (Issue #131) implements the **cluster-based
type assignment** option via a BFS expansion over the per-plate
adjacency graph derived from the (post-displacement) Voronoï
tessellation. Defaults `continental_fraction = 0.29` (Earth-like)
and `seed_cluster_count = 1` (single contiguous continent for the
§2.4 viewport-cadrable requirement). Spatially-biased seed sampling
remains deferred — empirically the BFS approach produces compact
clusters on the default 8-plate adjacency graph.*

### 6.3 Plate kinematics — open problem

Both v1 and v2, and now C1, share the same underlying problem: not
every seed produces an interesting cinematic configuration. A
random assignment of plate velocities can produce trivial flows
(all plates moving parallel — no convergence) or degenerate flows
(all plates moving toward a single point — single mass collapse).

Three plausible approaches:

- **Constrained random sampling.** Reject configurations that fail
  objective interest criteria (minimum convergent boundary length,
  at least one subduction zone, no parallel-velocity degeneracy).
  Coût: simple. Risk: criteria may be too strict and produce
  similar-looking maps.
- **Pre-filter scoring.** Generate K configurations cheaply, score
  each, present the top 3 to the user. Compatible with C1's
  sub-second 64² runtime.
- **Euler pole sampling.** Each plate's motion is a rotation around
  a randomly-placed pole on the (treated-as-spherical) domain. This
  is more physical than per-plate translation and tends to produce
  natural diversity of boundary types around a single configuration.
  Requires the kinematics to support rotation, not just translation.

Start with constrained random sampling. Revisit if the rejection
rate is too high.

### 6.4 Initial S̃ field

Per plate type: continental ≈ 1.0, oceanic ≈ 0.2, with smoothed
transitions across plate boundaries (Step 6 / Step 7 standard).
Cratonic cells start at the same continental value but are flagged
as immune. Add the R7-style boundary perturbation noise to break
the Voronoi rectilinearity in the S̃ field itself, not just in the
plate_id field.

### 6.5 Initial age field

The same §4.11 D2/D7 static identification from v2 Step 10 applies:
continental crust starts at some baseline age (the "geological past"
of the continent), oceanic crust ages from zero at oceanic ridges.
Reuse `age_field::init` from v2.

*Note: Phase 2 Track B (Issue #131) implements the **ridge-aligned
age = 0 initialisation** under **Path 3.A (init-only)** — at init
time, detect cells adjacent to divergent boundaries (reusing
[`tectonics_c1::boundary_classification::classify_boundaries`]
which already produces `BoundaryType::Divergent` from per-plate
kinematics), set `age = 0` on those cells, baseline elsewhere.
This addresses Phase 2 Track A's empirical finding that the
flux-form `∂_t·age + ∇·(age·v) = 0` advection produces ~1000×
density pile-up at convergent boundaries from initial uniform
oceanic `age = 0.5` (see `feedback_age_advection_density_vs_lagrangian`).
Path 3.A keeps the existing density-form advection unchanged;
escalation to Path 3.B (per-step ageing) or Path 3.C (Lagrangian
advection) is documented as fallback if Stage A reveals Path 3.A
does not preserve the Track A Spearman age-altitude correlation
sufficiently.*

---

## 7. Implementation plan

### 7.1 Phase 1 — prototype (target: 1–2 weeks)

Goal: a minimal C1 producing recognisable continental morphology on
one preset, at one grid size, validating the architectural choice.

- Strip v2's stokes module tree from the build (move to an `_attic`
  subdirectory rather than delete, so v2 reports remain runnable
  for regression).
- Add `tectonics_c1/` module skeleton with `time_loop.rs`,
  `closures/mod.rs`, `init/mod.rs`, `kinematics.rs`.
- Implement the time loop with constant-per-plate velocity, advection
  only (no closures yet).
- Visually verify that S̃ transports correctly: convergent boundaries
  produce visible thickening, divergent boundaries produce thinning.
  This is a sanity check on the advection alone; the result is not
  yet expected to be visually plausible.
- Add closure #1 (orogenic profile, Davis-Suppe). Re-run; verify
  orogenic morphology emerges at convergent boundaries.
- Add closure #2 (equilibrium height); verify the orogen plateaus
  at a physical height rather than growing unbounded.
- Add closure #4 (macro erosion) and #5 (isostasy + downstream
  pipeline); produce a heightmap and visually compare to v2 output.

Acceptance gate: Phase 1 succeeds if the output at 64² produces a
single continent that visually carries the chronological signature
(some older eroded regions, some younger sharper regions) without
having to "paint" it.

### 7.2 Phase 2 — boundary evolution + R7 init (1–2 weeks)

Phase 2 is split into four work tracks, each delivered as a separate
issue:

- **Track A — oceanic bathymetry (Stein & Stein 1992).**
  Status: ✓ **Complete (Issue #129, merged via PR #130).**
  Architecture C (post-isostasy bathymetry adjustment). MVP-table
  closure #3 is swapped from the original Parsons-Sclater 1977
  reference to its Stein-Stein 1992 successor (continuity with
  crossover at ~20 Ma). See §5.1 footnote.
- **Track B — R7 init (boundary displacement + continental
  clustering + ridge-aligned age, §6.1 / §6.2 / §6.5).**
  Status: ✓ **Complete (Issue #131, merged via PR #133).** Three
  sub-components: (1) Perlin/Simplex boundary displacement on
  Voronoï; (2) cluster-based BFS continental type assignment with
  cadrable-continent constraint; (3) Path 3.A ridge-aligned
  `age = 0` init at divergent boundaries (resolves Track A
  density-advection finding). Cadrable constraint deferred to
  Track B-bis (9 / 10 seeds wrap periodic boundary at 64²; see
  §6.2). Spearman ρ = -0.5233 IMPROVES on Track A baseline
  ρ = -0.476 (Δ -0.047) with 43 % age pile-up reduction.
- **Track D — boundary evolution: subduction, accretion, rifting
  (§4.5).** Status: ⏳ **In progress (Issue #132).** First C1
  work-track to mutate `plate_id`, `plate_type`, and `kinematics`
  per-step. Three closure modules in `tectonics_c1/closures/`:
  `subduction/`, `accretion/`, `rifting/`. Mass-conservation
  diagnostic test-only. Path 3.B event-driven age=0 propagation
  on rift-spawned divergent boundaries. *(Naming: the Issue #132
  rebrand swaps the historical Track C/D labels from earlier
  Track A + B project memory entries — Track D is now "boundary
  evolution" and Track C is now "kinematics sampling".)*
- **Track C — kinematics sampling (constrained random / Euler
  pole / scoring, §6.3).** Status: 📋 **Conditional (event-rarity
  escalation).** Phase 1.1 kinematics preset (8 plates, cardinal
  + diagonal velocities) is the default for Track D's acceptance
  runs. If Stage A's event-count diagnostic shows Track D events
  fire too rarely (< N per 300 steps systematically across seeds)
  under this default, Track C is escalated to produce more
  visibly active boundary evolution. Issue TBD.
- **Track B-bis — cadrable viewport offset for §2.4 compliance.**
  Status: 📋 **Pending.** Three remediation options documented in
  `c1_phase_2_track_b_acceptance.rs::acceptance_track_b2_continent_cadrable`
  docstring: constrained BFS seed selection, increased default
  plate count (8 → 12-16), spatially-biased seed sampling.
  Issue TBD.

Acceptance gate (cross-track): the same preset run multiple times
with different seeds produces visually distinct continents (not
just rotations of the same shape).

### 7.3 Phase 3 — second-wave closures (1–2 weeks)

- Subduction arc closure (Lallemand).
- Rift / passive margin closure (McKenzie-Buck).
- Foreland basin flexure (Beaumont).

Acceptance gate: the output contains all four major morphological
classes (continent-continent collision orogens, subduction arcs,
passive margins, foreland basins) on appropriate boundaries.

### 7.4 Phase 4 — UI + production (1 week)

- Per-closure toggle in the Bevy panel.
- Per-closure parameter sliders (one or two parameters each).
- Per-step display streaming (reuse v2's `StepProgress` callback).
- 512² target run profiled and tuned to < 10 s.
- Export pipeline end-to-end through HD erosion to PNG.

### 7.5 Phase 5 — optional enrichments

Hotspots, terrane accretion, post-orogenic thermal relaxation — only
if the Phase 4 output has visible gaps that these closures address.

### 7.6 Total estimated effort

5–7 weeks of focused work. This estimate has the standard caveats —
unforeseen calibration cross-coherence issues (§5.4) or closure
extraction difficulties from primary references can extend it by
50–100%. Compared to the v2 milestone (≈ 8 months of work across
Steps 0–13.5), this is roughly an order of magnitude less, reflecting
the much smaller algorithmic surface area.

---

## 8. Risks and open questions

### 8.1 Closure mathematical extraction

Several of the references in §5 publish their results in geophysical
notation with conventions tacit to the field. Davis-Suppe is the
classic example: the formulae involve dimensionless friction angles
and material properties that need careful interpretation. The risk
is that a naive transcription into C1 produces qualitatively wrong
results despite "implementing the formula".

**Mitigation:** for each closure, the implementation includes a
small unit test that reproduces a published figure from the source
paper (e.g., the Davis-Suppe Fig. 5 taper-angle prediction). If the
test passes, the implementation is at least consistent with the
paper's own examples.

### 8.2 Inter-closure cross-calibration

§5.4 anticipates this. The risk is the same as v2 Step 8's slab+mantle
co-calibration runaway: two closures individually validated produce
divergent joint behaviour. C1's regime is much more benign (no
implicit coupling, no condition number sensitivity), but quantitative
drift is still possible (e.g., orogenic source produces more mass
than oceanic recycling consumes).

**Mitigation:** global mass-balance diagnostic, enforced as a test.

### 8.3 Periodic-domain artifacts

The viewport scrolling strategy (§2.4) works if and only if the
continental cluster fits within a 512² window with ocean buffer
around it. A sufficiently large or elongated cluster may not fit.

**Mitigation:** the constrained kinematics sampler rejects
configurations whose initial continental footprint approaches the
torus extent. Recovery from bad seeds is by re-sampling, not by
algorithmic patching.

### 8.4 Kinematics seed quality

§6.3's open problem. Some fraction of seeds will produce visually
uninteresting maps regardless of closures.

**Mitigation:** filter at the kinematics-sampling stage. This pushes
the "interesting map" criterion into the input rather than the
output, which is more tractable but assumes the criterion can be
formalised.

### 8.5 Departure from physical fidelity

C1 trades physical fidelity for visual plausibility. The closures are
empirical fits, not solutions of a balance equation. The output is
defensible at the qualitative level (the closures encode the validated
physics of their source phenomena) but is not a quantitative
prediction of anything. This is the right trade for Living Landz; it
would be the wrong trade for a geophysics research tool.

**Mitigation:** documented as a property of the architecture, not a
bug to be fixed. C1 does not present itself as a physical predictor;
it presents itself as a procedural generator informed by physics.

### 8.6 R7 init not yet validated

Boundary displacement on Voronoi has not been tested. It may or may
not produce the morphological richness C1 needs. If it does not, the
second and third alternatives in §6.1 are available; if none of them
work, the architectural assumption "plate tessellation can be
non-rectilinear" needs revisiting.

**Mitigation:** R7 prototype is the first checkpoint in Phase 1. If
the boundary perturbation does not produce visible coastline
diversity at 64², the rest of C1 is built on shaky ground.

---

## 9. What v2 is good for, after C1 ships

v2 is not deleted. The non-dimensional solver remains in the source
tree (gated behind a feature flag or moved to an `_attic` namespace)
because it has two residual values:

- **Regression baseline.** The v2 reports (`docs/reports/step{0..13_5}_*.md`)
  document validated output at every step of the milestone. Re-running
  them periodically catches infrastructure regressions in the shared
  components (advection, init, isostasy, erosion HD).
- **Reference physics.** If a future question arises about whether
  a C1 closure deviates from "what the underlying thin-sheet physics
  would do", v2 is the available comparison point. The answer may be
  "C1 deviates, deliberately, because the closures are calibrated on
  observation rather than on the v2 formulation". But having the
  comparison available is useful for sanity-checking surprising C1
  output.

There is no plan to develop v2 further. Step 14 (oceanic FBM
heterogeneity) and the slab+mantle co-calibration follow-up issue
are both closed as won't-fix.

---

## 10. Decision record

The C1 architecture was selected over options A, B, C2 after a
structured discussion captured in the design log. The key drivers
of the decision were:

1. The chronological-diversity requirement (§2.2) rules out static
   architectures (A, B) unless chronology is explicitly painted.
2. The runtime envelope (§2.3) rules out v2-class implicit-coupled
   architectures.
3. The visual emergence requirement (§2.2, "no FBM-faked variation")
   rules out post-hoc geometric perturbation as a primary mechanism;
   variation must come from a system that produces it intrinsically.
4. The system-code analogy (§3.5) makes empirical closures the
   natural framing for source terms in a dynamic balance, rather
   than for static forms.

These four constraints are jointly satisfied by C1 and not by any
of the alternatives. The decision is not "C1 is the best possible
architecture" — it is "C1 is the architecture that satisfies the
known constraints, with the explicit understanding that some of
those constraints may shift as implementation reveals new
information". The plan §7 is structured to fail fast at the first
checkpoint (Phase 1 acceptance gate) so that a wrong-architecture
decision is detected at the cost of 1–2 weeks of work, not 6 months.

---

## 11. Implicit physical scales

C1 operates in non-dimensional units throughout — closures, tests,
visualisation, and reports all use non-dim values directly. The
table below makes the implicit physical scales explicit so that the
values shipped in code can be compared to literature without
chasing the conversion ad hoc.

The values in the "C1 default" column are read directly from the
shipped code (`DavisSuppeParams::default()`,
`EquilibriumHeightParams::default()`,
`PlateKinematics::preset_phase_1_1`,
`init_c1_state_phase_1_1`).

| Quantity                    | C1 default                         | Physical interpretation              |
|-----------------------------|------------------------------------|--------------------------------------|
| `S̃ = 1.0`                  | continental initial (`init_s_field`) | ~35 km — normal continental crust   |
| `S̃ = 0.2`                  | oceanic initial (`init_s_field`)   | ~7 km — normal oceanic crust         |
| `h_max = 2.5`               | Davis-Suppe wedge plateau (Phase 1.2) | ~87 km — orogen wedge ceiling       |
| `h_eq = 2.0`                | equilibrium-height collapse target (Phase 1.3) | ~70 km — observed Tibet plateau crustal thickness |
| `k_collapse = 2.0`          | quadratic-formula coefficient (Phase 1.3, post-E1.bis) | calibrated, not literature-derived (see §11.1) |
| Domain `1 × 1`              | non-dim grid extent                | ~1000–5000 km regional, 64²–512² resolution |
| `dx = 1 / N`                | cell size at `N × N`               | ~15–80 km at 64² (regional), ~2 km at 512² (sub-regional) |
| `\|v\| ≈ 0.005 – 0.011`     | `PlateKinematics::preset_phase_1_1` magnitudes (cardinal 0.01, diagonal `√2·0.008 ≈ 0.0113`) | ~5–10 cm/yr — typical convergent rate |
| `Δt ≈ 0.69`                 | CFL non-dim/step at 64² (`0.5·dx/max\|v\|`) | ~30–100 ka per step                |
| `300 steps`                 | typical Phase 1.1–1.3 run length    | ~10–30 Ma orogen evolution           |

Phase 2 Track A (Issue #129) introduces two additional scales,
documented separately because the values are tied to the
Stein-Stein 1992 closure rather than to the global time-loop /
grid setup:

- `SteinSteinParams::age_to_ma = 0.667` — `1 age step ~ 0.667 Ma`;
  default chosen so that the canonical `300 steps` run spans
  `~200 Ma`, matching the upper end of typical oceanic-plate
  lifetimes from ridge to subduction.
- `SteinSteinParams::depth_scale_m = 5000` — converts the S-S
  metric depth `(2600 m at ridge → 5651 m at asymptote)` to the
  non-dim altitude range `~ 0.52–1.13`, consistent with the
  Phase 1.4 isostatic altitude convention.

Both are consumed by `apply_stein_stein_bathymetry` (Architecture
C, see §5.1 footnote).

Phase 2 Track D (Issue #132) will introduce additional parameter
scales for the three new closures — `K_subduction` (rate
coefficient for oceanic consumption), `arc_efficiency` and
`arc_distance` (arc-volcanism distribution), `K_rift` (rate
coefficient for divergent continental thinning), and the time /
thickness thresholds for accretion-merge and rifting-split. The
defaults will land here as they are calibrated in Track D
Stages E1 / E2 / E3, alongside their physical interpretation.
Track D follows the calibration-via-visual-review discipline
(§11.1) under tier 2 (analytical first-pass + visual review, no
published universal coefficient) — same tier as the Phase 1.3
`k_collapse` and Phase 1.4 `K` (stream-power erosion)
calibrations.

The mapping is **not enforced in code** — there is no dimensional
analysis pass, no unit checking, no `Quantity` wrapper type.
Closures and tests operate on the bare `f64` non-dim values. The
table exists so that parameter calibration discussions can refer
to "Tibet-like" or "Andean-like" magnitudes without ambiguity, and
so that visual interpretation (cycle_NNN_s.png at palette
`[0, 3.0]`) carries an implicit "this is 0–~105 km of crustal
thickness".

### 11.1 Calibration via visual review, not dimensional derivation

`k_collapse = 2.0` (Phase 1.3) and the planned `K` for stream-power
erosion (Phase 1.4, §5.1 closure #4) are **calibrated for visual
balance**, not derived from the underlying physical scales above.
This is consistent with the design doc's overall position on
empirical closures (§3.5, §8.5): each closure is an additive
source term tuned so that its visual signature is legible against
the other active closures, not a quantitative prediction of any
specific geological context.

This matters for Phase 1.4 specifically because Lague 2014 *ESPL*
39(1), 38–61, the most recent compilation of stream-power
parameters, **explicitly declines** to publish a universal `K`
value, arguing that `K` aggregates lithology, climate, threshold
effects, and grain size, and that any single number is misleading.
C1 will ship a `K` calibrated for visual balance against the
Phase 1.2 + 1.3 closures, documented in the closure module
docstring with the analytical first-pass estimate plus the visual-
review iteration record.

### 11.2 Future formalisation (Phase 2+)

Phase 2 may need quantitative agreement with measured data (e.g.,
when comparing C1 output to the Lallemand subduction-arc relations
in closure §5.2 of this doc). At that point a thin
`tectonics_c1::scales` module could expose the table above as
compile-time constants and helper conversion functions
(`s_to_km`, `dt_to_ka`, etc.).

This is deliberately deferred: Phase 1.x ships without dimensional
infrastructure because the closures and acceptance tests do not
need it. Premature formalisation would add a typing surface that
the current code does not consume.

### 11.3 Cross-references

- §5.1 (MVP closures) — table entries `h_max` and `h_eq` originate
  from closures #1 (Davis-Suppe) and #2 (equilibrium height) of the
  MVP set.
- §8.5 (Departure from physical fidelity) — the rationale for
  empirical-calibration-over-dimensional-derivation lives there;
  §11.1 narrows it to the specific calibrated-parameter cases.
