# Finding 176 — declared before any measurement

Written after reading the accretion closure, the subduction closure and the sources, before any code of F176. Production keeps its
accretion rule; everything below is behind the continent profile's flags.

## Part 0 — the inventory

**The accretion closure** (`crates/ymir-core/src/tectonics_c1/closures/accretion/`):
- **The trigger** (`merge.rs`, `ConvergenceTracker::update`, l. 110–178): every step, for every pair of plates sharing a boundary,
  the 4-neighbour edges are classified by v_rel · n̂ (the same convention as `classify_boundaries`); a pair with more convergent
  than divergent edges increments its counter by 1, any other pair is reset to 0. **The boundary's type (continental or oceanic)
  is not read.**
- **The merge** (`merge.rs`, `apply_accretion_step`, l. 216 on): a pair whose counter reaches `merge_time_threshold` = **50**
  (`params.rs` l. 66) merges — the lower plate id wins, takes the mass-weighted mean velocity, and every cell of the loser joins
  it (`plate_type` kept per cell).
- **Why it was introduced** (`params.rs` l. 1–26, `merge.rs` l. 1–60): Issue #132, Track D Stage E2 — « sustained-convergence plate
  merge »; the threshold raised from 20 to 50 steps (« ≈ 33 Ma of sustained convergence — consistent with collision-zone
  timescales … Indian-Asian collision matured over ~50 Ma », and to suppress spurious merges). The rationale is a **collision**;
  the rule is applied to every convergent pair.
- **Its wiring** (`time_loop.rs`, step 5c/5d): `tracker.update(plate_id, kinematics)` after the subduction and the rift thinning,
  then `apply_accretion_step`.

**The continent-ocean rule takes over** (`closures/subduction/source_term.rs`, `apply_subduction_step`): an oceanic convergent cell
whose 4-neighbour of **another plate is continental** is consumed (Δs = rate × v_rel · n̂ × dt), the arc mass goes to the
continental cells within `arc_distance`, and a floor-triggered cell is promoted. So once option 5 has consumed a seaward band
(its distal cells reassigned to the overriding plate, F175-A3), the band's continental cells (the necking domain) touch the other
plate's ocean and the existing rule consumes that ocean beneath the continent — the Andean margin. F175 observed it on 2 seeds
(1 group each), the merge cutting most boundaries before.

**The terrestrial reference** (comparison only, never a target):
- **Bradley 2008**, Earth-Science Reviews 91, 1–26 (`docs/refs/Brad_08.pdf`): the present-day passive margins total **105 000 km**,
  « even longer than the spreading ridges (65 000 km) or the convergent plate boundaries (53 000 km) » (p. 2, §1; and p. 7, §3);
  their length-weighted mean age is ~104 Myr (p. 7); the 76 ancient margins whose lifespan could be measured average 181 Myr
  (abstract, p. 1).
- **Active continental margins**: no source in `docs/refs` splits the 53 000 km of convergent boundaries into ocean-continent and
  intra-oceanic. Syracuse & Abers 2006 (`syracuse2006.pdf`) compile 33 000 km of arcs in 45 segments but do not classify them, and
  classifying them by name would be my judgement. **So the comparison uses Bradley's bound**: if every convergent boundary bordered
  a continent, active margins would be at most 53 / (105 + 53) = **34 %** of the continental margins' length; the true share is
  lower (part of the 53 000 km is intra-oceanic). **The author is asked for a citable split** (e.g. Leat & Larter 2003, Geological
  Society Special Publication 219, or Stern 2002, Reviews of Geophysics 40, on intra-oceanic subduction systems — to check).

## Part A — the rule (behind the flag)

### A1. Merge only at collision

- With the flag (`Margin::merge_at_collision`, passed to C1 through `run_with_closures_terms_flags`), the tracker counts **only
  the edges whose two cells are continental** (`ConvergenceTracker::update_collision`): a pair is net-convergent when it has more
  convergent than divergent continent-continent edges, its counter then increments, otherwise it resets. The threshold (50) and
  the merge itself (the whole loser plate joins the winner) are unchanged.
- So a convergent **ocean-ocean or ocean-continent** boundary never merges: the subduction continues as long as the convergence
  lasts, carried by the existing continent-ocean rule and option 5.
- Declared reading of the brief's « le long des segments où les cellules de frontière sont continentales des deux côtés »: C1's
  merge is plate-wide, so the rule decides **when** a pair merges (on its continent-continent edges), not a partial merge.

### A2. Unchanged

Option 5 (F175's A1–A7, D_init 2.13 cells of 64²) and the seaward margin are on. Nothing else changes.

### A3. The controls (permanent)

- The flag off: C1 bit for bit (production: the existing entry points) and F175 bit for bit (`run_with_closures_terms_oo`).
- A world with no oceanic convergent boundary (every cell made continental) gives the same result with and without the flag.

### A4. The viz

A « Fusion seulement à la collision (C1) » box.

## Part B — the variant measured only: emerged arcs

- Under the cascade's Airy mapping, an oceanic cell that has received arc mass is read at **0.2 + its cumulative arc mass** (the
  positive part of C1's subduction term at that cell, s units, 1 ↔ 35 km) instead of 0.2 (`Margin::arcs_emerge`). It is C1's own
  mass budget; no constant is added.
- **The source asked by the brief** (the crustal thickness of arcs): Christensen & Mooney 1995 give **continental arcs 38.7 ±
  9.6 km** (Table 3); no source in `docs/refs` gives island-arc crust — the author is asked for one (e.g. Calvert 2011 or Suyehiro
  et al. 1996 on the Izu-Bonin arc — to check).
- Reported: the islands created (64² at C1, 256² at young) and the circumnavigation.

## Part C — the measures, by tiers (the budget rule)

**The seeds**: the panel 8, 18, 19, 23, 29; the screening adds the **10 smallest seeds of F172's catalogue outside the panel: 1, 2,
3, 4, 5, 6, 7, 9, 10, 11** (a fixed rule, no choice).

**The tiers** (announced before running; < 45 min wall; each world's lines written to `results/` as soon as ready):
1. **C1 only**, the panel at the fixed continent (F175's = F174 B′'s plates): F175 (`Margin::seaward_oo`) | F176 (F175 + merge at
   collision) | F176 + variant B.
2. **The screening**: the 15 seeds, the selection redone once under F176's rule (parallel); C1 on the selected continents.
3. **The physics level at 256²**: the panel, the fixed continent, young and old: F175 | F176; variant B at young.
4. **The cascade to 2 048²**: 1–2 F176 worlds with an Andean margin, young, for the images.

**The quantities**:
- **The margins**: the continental perimeter's active share (a convergent boundary in the 3 × 3, F174's rule) every step and at the
  end; **Andean margins** = 8-connected groups of ≥ 3 perimeter cells passive at t₀ and active at the end; their count and length
  (cells × 6.25 km) **per continent** — a continent = a 4-connected component of continental cells at the end with ≥ 20 cells;
  over the 15 seeds, the share of continents with ≥ 1 Andean margin and the distribution of the active share.
  - **A supplementary reading** (declared now): « active within 3 cells » — a perimeter cell with a cell consumed (either rule)
    within 3 cells (≈ 19 km on the map, ≈ 140 km real at ×7.5) at the last 10 steps; the 3 × 3 rule only sees a trench at ≤ 1 cell
    (≈ 47 km real), stricter than a real forearc (100–200 km).
- **The plates**: their count every step (16 → ?).
- **The terrestrial comparison** (Part 0), without target.
- **The rampart**: W₁₀ at young, global and active / passive (F175's `w10_by_type`).
- **The circumnavigation** and r*; **the selection**.
- **Option 5's validations V1–V6** (F175's definitions).
- **The macro**: classes, median, land.
- **Variant B**: islands, circumnavigation.
- **Images**: the margins at the end (active, passive, arcs) on the panel; a transect across an Andean margin; the whole maps
  F175 | F176.

**The success** (the brief's): ≥ 70 % of the continents with ≥ 1 Andean margin (over the 15 seeds); the median active share in
10–50 %; the controls pass. **The stops**: the selection fails on ≥ 3 of 15 or the circumnavigation is lost on ≥ 2; an active share
above 70 % on ≥ 3 seeds (the rampart back); any other closure change or solver. **No constant tuned to make Andean coasts appear.**
