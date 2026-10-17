# Finding 175 — declared before any measurement

Written after reading C1's subduction closure, Davis-Suppe's seeding and the time loop, before any code of F175 was written.
The rule follows `docs/reports/relief_method/f174_margin/f174_option5_note.md` with the reviewer's amendments (the F175 brief,
A1–A7).

## The frame

- **No solver.** The change is to two existing closures of C1 — the subduction closure (`apply_subduction_step`, extended by an
  ocean-ocean path) and Davis-Suppe's seeding (cut at ocean-ocean boundaries) — **behind a flag**:
  `OoSubductionParams`, carried by the continent profile (`cascade::margin::Margin::oo`), `None` by default.
- **The flag reaches C1 through new entry points** (`run_with_closures_terms_oo`, `run_c1_sources_mapped_oo`); `C1Closures`
  is **not** changed, because its `Debug` is part of the tectonic cache key (`cached_product.rs`, `tectonic_key`). The existing
  entry points pass `None`: C1 bit for bit.
- **ALGO_TECTONICS is bumped (1 → 2)**: the time loop's code changes even though the off path is bit for bit; the bump makes
  the cached tectonics recompute, so the guard 6/6 exercises the new code instead of reading the cache.
- The seaward margin (F174 B′, `Margin::seaward()`) is on in this finding; the sea stays out; K2, the age mapping, the torus and
  the 16 plates with their fallback are frozen.

## A1 — what is consumed

- **An ocean-ocean boundary cell**: a `Convergent` cell of `Oceanic` type with **no** 4-neighbour of another plate that is
  continental (those are the existing continent-ocean rule's). Its partner is the 4-neighbour of another plate, oceanic, with the
  largest positive v_rel · n̂ (the existing rule's convention).
- **Which side sinks**: the cell whose `state.age` is **larger** (the older, F174-D §1: Stern 2004 p. 276; Stern & Gerya 2018
  p. 26). **Ties** (equal ages) go to C1's existing velocity rule: the velocity-based upper-plate flag of `classify_boundaries`
  (read before A7 clears the O-O seeds): the cell sinks if it is not flagged upper and its partner is; if both or neither are, the
  lower plate id sinks (a declared deterministic tie-break).
- **The amount**: Δs = `consumption_rate` × v_rel · n̂ × dt, clamped at s — the existing continent-ocean rule's rate (0.5).
  **No new constant.**
- **A declared limit**: `age` is not a clock in C1 (it is advected, set to 0 at ridges and splits, 0.5 on initial ocean, never
  incremented): only the relative order at the boundary is read. Expected consequence, recorded before measuring: many ties
  (both sides at 0.5), decided by the velocity rule; the share decided by age is reported.

## A2 — the arcs

- The arc mass (Δs × `arc_efficiency`, 0.5) goes to the **overriding plate's** cells within `arc_distance` (3) of the consumed
  cell (BFS, 4-connected, depth ≥ 1): its continental cells if any are reached, **otherwise its oceanic cells** (an island arc).
  Nothing reached: lost, and counted (`SubductionStats`).
- **A declared consequence, before measuring**: the cascade reads C1's oceanic cells at s = 0.2 under Airy (F167,
  `Airy::oceanic_s`), so an arc built on oceanic cells has **no altitude**. Island arcs cannot emerge with A1–A7 alone; making them
  emerge needs a rule beyond A1–A7 (e.g. promoting an arc-loaded oceanic cell to continental, or a mapping change). If S3 fails
  for this reason, it is reported as a stop point, not implemented.

## A3 — the boundary migrates

- A consumed cell whose s falls under the existing `plate_id_reassign_threshold` (0.05) and that has **≥ 2** 4-neighbours of
  the overriding plate with positive convergence is **reassigned to the overriding plate**, its type staying oceanic.
- **The check asked by the brief**: the continent-ocean rule already reassigns a floor-triggered cell to the overriding plate,
  with the same threshold and the same ≥ 2-neighbour rule (`source_term.rs`, step 4 — « reassign `plate_id[c]` to the
  continental neighbour's plate id »). It also sets the type to continental, because its overriding plate is continental; the
  ocean-ocean path keeps the type oceanic, as the brief specifies (the overriding side is oceanic). **The plate reassignment is
  the existing mechanism**; only the type follows the overriding side. Not stopped.

## A4 — the trigger

- **Per unordered plate pair** (a, b), a cumulative convergence D_ab += mean(v_rel · n̂ > 0 over the pair's ocean-ocean contacts
  of this step) × dt, updated at each step's start from the step's boundary classification. The pair becomes subducting when
  **D_ab ≥ D_init**, and stays so.
- **D_init = 100 km (Stern & Gerya 2018 p. 44, after Hall et al. 2003) ÷ 7.5 = 13.3 km on the map = 2.13 cells of 64²**, i.e.
  2.13 × `config.dx` in C1's length units.
- **No rate threshold** (a declared limit): C1's time has no fixed physical duration, and a rate threshold would couple C1's
  outcome to the age slider.
- An **initiation** is the step a pair becomes subducting.

## A5 — the continent at the trench

No new rule: once a band (or rim) is consumed, the continental cells meet the other plate's ocean, and the **existing
continent-ocean rule** consumes that ocean beneath the continent (an Andean-type margin). A continental cell never sinks (#145).

## A6 — no spontaneous collapse

No initiation rule other than A4; a passive margin becomes active only through A1–A3.

## A7 — Davis-Suppe at ocean-ocean boundaries

With the flag, after `classify_boundaries` and `retarget_upper_plate_continental`, the upper-plate seed is **cleared on every
ocean-ocean convergent cell** (an oceanic convergent cell with no 4-neighbour of another plate that is continental), before the
wedge distance is computed. Davis-Suppe is unchanged at continent-ocean and continent-continent boundaries.

## A8 — the controls (permanent)

- The flag off: C1 bit for bit (s, age, plate_id, plate_type on a 300-step run) and F174 B′'s initial state and C1 bit for bit.
- D_init infinite: no ocean-ocean consumption (0 cells consumed, 0 reassigned); the run differs from F174 B′ by A7 only — the
  difference is reported (the bench).

## A9 — the viz

The Cascade window: a « Subduction océan-océan (C1) » box (with the seaward margin), and a map of the margins at the end (active,
passive, arcs) in the terms' frieze.

## B — the measures

**The grid**: the panel 8, 18, 19, 23, 29; four ages (B0, young, intermediate, old); the physics level at 256², the chain to
2 048², the torus, no sea, K2.
- **At a fixed continent** (F174 B′'s selected plates): Part A (no margin) | the seaward margin (B′) | the margin + option 5.
- **With the selection redone** (the same rules): Part A (its own selection) | B′ (its own) | the margin + option 5 (its own).

**The validation criteria (the note's, amended)**:
- **V1**: the share of initiations at the edge of a passive margin — an initiation whose pair's contact cells lie within 2 cells
  of a continental cell.
- **V2**: each initiation's distance to the nearest cell consumed (by either rule) at the previous step, against 200 km
  (1 500 km ÷ 7.5 = 32 cells of 64²); « none » when nothing was consumed before.
- **V3**: the number of initiations and their steps (no target).
- **V4**: every age-decided consumption has the older side sinking, i.e. the arcs on the younger side (counted; 0 violations
  expected by construction); the share decided by age vs by the velocity tie-break.
- **V5**: no one-cell finger (an oceanic cell with ≥ 3 of its 4 neighbours on another plate, at the end, against B′), and the arc
  mass conserved (consumed × efficiency = distributed + lost, both counted).
- **V6**: the active share of the continental perimeter at the end (a convergent boundary in the 3 × 3, F174's rule) and the
  **Andean margins**: 8-connected groups of ≥ 3 perimeter cells passive at t₀ and active at the end; counted per world.

**The other quantities**:
- Davis-Suppe's ocean-ocean seeds per step (0 expected with the flag) and Σ Δs_DS on the band (off Part A's continent).
- W₁₀ at young, global and **by segment type**: the land within 10 km of the coast near an active vs a passive perimeter cell
  (the 64² end map, nearest perimeter cell).
- Islands ≥ 1 km² off the main land: count, area, and the share lying where C1's subduction term deposited arc mass.
- Circumnavigation: r* of C1's end land (64²).
- The macro (classes, median, land; seed 8's mountain share against B′), the coast by rock type, the seaside slope, lakes,
  texture, images, cost — as in F174 C′.

**The success** (the brief's, on ≥ 4 of 5 seeds, read with the selection redone): 0 Davis-Suppe O-O seeds and no world above
80 % mountain at young; 10–60 % of the perimeter active at the end; V4 and V5 pass; circumnavigable. **The stops**: the
selection fails on ≥ 2 seeds; the arcs close the passage (circumnavigation lost on ≥ 2 seeds); any closure change beyond A1–A7.

### Amendments after the smoke run (seed 8, young only), before the full bench

- **The bench's success code** was brought back to the declaration: V5 is read against B′'s finger count (« against B′ »), not
  against 0; Davis-Suppe's O-O seeds are 0 by construction once the cut ran (the permanent test checks that a second cut finds
  none) — the smoke had compared two counts taken one step apart.
- **A counter added** (no rule changed): the age-decided consumptions whose two ages differ by < 0.01. The smoke decided every
  consumption by age (no exact tie): the advected ages are never exactly equal, so the order is often read from the advection's
  smearing; the counter reports how often.
- **The tracker's update** was moved before the continent-ocean subduction of the step (it reads the plates as they enter the
  subduction stage), so the bench's replica of it, fed with the previous step's end state, locates the same initiations (the
  pair counts agree: 30 / 30 and 24 / 24 on the smoke).

### Amendment after the measurement: the ALGO_TECTONICS bump reverted

- With the bump (1 → 2) the guard 6/6 recomputed every cached state (C1, then the HD upscale and erosion, whose keys chain from
  the tectonic key: 3 189 s instead of ~140 s) and found **all six fields and lakes different** (`checks_after_with_algo_bump.txt`).
- A direct check settles the cause: a hash of C1's production run (64², 300 steps, the default closures; seeds 42, 8 and the old
  témoin) is **identical bit for bit** between HEAD (`b080a3b`, a separate worktree) and the F175 tree
  (`tests/f175_c1_hash.rs`: c69acd85a0f7f484, 3565285d65e4c341, fd23c309498a9761 on both). F175 does not change C1's
  production output; the guard's references were matching **cached products that the current code no longer reproduces** — a
  stale cache older than F175, revealed by the forced recompute.
- So the bump is reverted (its only purpose was to make the guard exercise the new code; the hash check does that directly), and
  the stale cache is reported to the author, not fixed.
