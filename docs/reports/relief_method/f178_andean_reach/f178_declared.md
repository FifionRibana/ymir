# Finding 178 — declared before any measurement

**One closure modified — where Davis-Suppe's continent-ocean wedge starts (the seed's reach) — behind the continent profile's
flag; no solver; production unchanged.** The measured configuration is F176's profile (the seaward margin + option 5 + the merge
at collision) plus the new seed (`Margin::andean_gap`). Variant A2 (the arc at the arc-trench distance) is measured only.

## Part 0.1 — the guard

- **The exact verdicts at HEAD** are read from the per-state `Debug` line `f123_viz_guard` already prints (`{label}: field
  {:?} · lakes {:?}`), now kept by the checks' grep (`: field `, `panicked`, ` field: `, ` lakes: `).
- **F177's proposal applied**: `f123_viz_guard` prints « = banc » for `Match`, « ≠ (bench …, viz …) » for `Mismatch` and « non
  gardé (no reference) » for `NoReference`; the checks keep the panic line. This edit is in the test module of
  `crates/ymir-viz/src/ui/workspace.rs`, the author's uncommitted file; at commit only that hunk is staged (the author's own
  changes stay out).
- **Why F175 bumped `ALGO_TECTONICS`** (its `f175_declared.md:15`): « the time loop's code changes even though the off path is
  bit for bit; the bump makes the cached tectonics recompute, so the guard 6/6 exercises the new code instead of reading the
  cache ». It followed the memory rule « bump on any code change to a cached step » — a rule about invalidation, not about
  output. The output was proven unchanged by the C1 hash (`tests/f175_c1_hash.rs`). The bump changed the digest, so the guard
  read `NoReference`.
- **The re-registration rule** (for any state read `NoReference`): recompute fresh, compare the field and the lakes with the
  existing references' hashes regardless of the digest; bit-identical → re-register under the current digest (a
  re-registration, not a re-baseline); different → stop and report. At HEAD no digest moved (F177 reverted nothing that the
  guard keys on), so all six are expected `Match`, leaving nothing to re-register.
- **Budget**: the guard through the checks (~6 min, the cache warm); the 45 min two-tier fresh run only if a state reads
  `NoReference`.

## Part 0.2 — the sources (both in `docs/refs`)

- **Syracuse & Abers 2006** (G³ 7, Q05017, doi:10.1029/2005GC001045):
  - p. 1 (abstract): H, the slab depth beneath the arc, ranges from 72 to 173 km over 500 km arc segments, with a global
    average of 105 km; 839 volcanic centres, 33 000 km of arc;
  - **p. 10, §26**: « Horizontal distances from volcano to trench range from 85–90 km in Scotia and Vanuatu to 470 km in
    Alaska, with the middle 50% of arc-trench distances between 180 km and 275 km »; « 99% of H estimates lie between 58 and
    203 km and 50% between 85 and 119 km … H averages 105 km with a median of 100 km »;
  - **no median and no 90th percentile of the arc-trench distance are given in the text** (only Fig. 8's histogram, not
    read numerically; the per-volcano auxiliary dataset is not deposited).
- **England, Engdahl & Thatcher 2004** (GJI 156, 377–408, doi:10.1111/j.1365-246X.2003.02132.x):
  - p. 377 (summary): the depth to the top of the intermediate-depth seismicity beneath volcanoes ranges from 65 to 130 km;
  - p. 379, Table 1: per arc, that depth, the slab dip between 80 and 400 km (30–70°, also p. 382), the convergence rates,
    the age and the down-dip slab length;
  - **it gives no arc-trench distance.** A distance from depth and dip (D / tan δ) would assume a planar slab from the trench,
    which its own dips (measured below 80 km) do not support. So **England et al. is not used to set D_gap**; it confirms the
    depth range (65–130 km, consistent with Syracuse & Abers' 72–173 km segment means).

## D_gap (declared; one value tested for adoption; no tuning loop)

- **D_gap = the largest observed arc-trench distance**: 470 km (Alaska; Syracuse & Abers 2006 p. 10) ÷ 7.5 = 62.7 km world ÷
  6.25 km per cell = 10.03 → **10 cells of 64²** (a walk depth ≤ 10). The reasoning: beyond that distance no arc forms on the
  overriding continent.
- **The sensitivity** (the brief's « median, computed the same way »): **the source gives no median**, only the middle 50 %,
  180–275 km, which brackets it. So the sensitivity is run at both bounds: 180 km → 3.84 → **3 cells**; 275 km → 5.87 → **5
  cells** (C1 only; at 256² only if D_gap gives no cordillera).
- **A2's arc distance**: the arc on P at the middle 50 % of arc-trench distances, i.e. walk depths **4–5 cells** (180 and 275
  km → 3.84 and 5.87 cells, the integer depths inside).

## Part A — the closure

- **A1, the Andean seed** (`andean_seeds`, called in C1's time loop after A7's cut, only when `andean_gap` is set). A
  continental cell c of plate P becomes a continent-ocean seed of Davis-Suppe when:
  - **P overrides an active convergent boundary**:
    - **O-O**: an oceanic boundary cell b of P that is an ocean-ocean contact of a subducting pair (the D_init tracker as the
      step's Davis-Suppe stage sees it) and not its sinking side by option 5's rule (the older cell sinks; ties: the velocity
      flag read before A7, then the lower plate id; the partner of largest v_rel · n̂, as in `apply_oo_subduction_step`);
    - **C-O**: a continental cell of P touching another plate's oceanic cell on a convergent boundary is already a seed (the
      existing rule), unchanged;
  - **c is on P's continental front facing that trench**: a multi-source walk (4-connected, periodic) from the trench cells b
    through **P's cells only** stops at P's continental cells; the first continental cells reached are the front;
  - **c is within D_gap of the trench**: the walk depth (b at 0) ≤ D_gap.
  - The seeds join `upper_plate_mask` and the O-C seed mask, so the wedge is the O-C ridge (`h_max · exp(−d / l_taper)`).
    `wedge_distance_intra_plate_typed`, Davis-Suppe's formula, its reach, coupling and decay are unchanged.
  - **A declared consequence**: the Davis-Suppe step skips convergent cells (the boundary). A front cell is not a boundary
    cell, so it receives the source at d = 0 (the O-C target's peak, h_max). The existing O-C seed is a boundary cell and is
    skipped. Nothing in Davis-Suppe's step is changed.
  - No continental cell of P within D_gap → no seed. A7 stays, so a purely ocean-ocean boundary gets no dome.
- **A2, the variant (measured only)**: the O-O arc (`distribute_arc_mass_plate`) goes to P's cells at walk depths 4–5 from the
  consumed cell (continental first, else oceanic, as now), instead of within 3 cells; if none lie there, the existing rule
  applies. The C-O arc is untouched (its trench is at the coast).
- **A3, the controls (permanent)**: the flag off is F176 bit for bit (and production's hash, `f175_c1_hash`); `andean_gap = 0`
  (no continent within reach of any trench) and an all-continental world give the same result with and without the flag.
- **A4, the viz**: « Graine andine à distance arc-fosse (C1) ».

## Part B — the measures (no other closure modified)

The 15 seeds of F176's selection (the profiles cached by F177, `f177_accretion_close/results/sel_seed*.txt`). One C1 cell =
6.25 km world (46.9 km terrestrial).
- **Tier 1 (C1)**: F176 (baseline), F178 (D_gap 10), F178 + A2, and the sensitivity at 3 and 5 cells:
  - **the A1 share** with the two tie rules: « any » (F177's) and « majority » (strictly more active than passive boundary
    cells at the walk's first depth; a tie → passive);
  - **the uplift by origin by distance band** from the coast: every continental cell assigned to its nearest coast cell by a
    4-connected walk through continental cells, its depth d. The bands are 0–3 (d = 0, 1, 2), 3–6, 6–15 and 15–30 cells,
    active against passive coasts (the A1 class of the nearest coast cell). The origins (Σ Δs → Δh by the declared Airy, m):
    - **the C-O wedge**: Davis-Suppe on `oc_wedge` cells;
    - **collision**: the other Davis-Suppe;
    - **arcs**: every positive subduction increment;
    - **other**: the rest of Δs_src.
    This replaces F177's 30-cell band;
  - **the seeds whose C-O wedge is positive on the active coasts**: the mean over the 0–6 cell bands > 0 m;
  - **the trench continuity**: the active trench cells are the overriding-side O-O cells of subducting pairs (A1's rule) and
    the oceanic convergent cells touching another plate's continent. Segments are 8-connected; their length is the geodesic
    diameter × 6.25 km. Reported: the number of segments, the median, the 90th percentile and the maximum.
- **Tier 2 (256², young, the 15 seeds)**: F178, and F176 as the baseline:
  - **F177's tiers 2–4 measures**: B3's origin maps (Andean = the C-O wedge + the C-O arcs; collision; O-O arcs; other) and
    the formations with F177's complement definitions. Andean cordilleras are 8-connected Andean-origin mountain groups
    touching an A1-active coast, with their length along the coast. The other formations: the collisional highland
    (separated or not), the arc mountain, other. Condition (v) from those;
  - **W₁₀** at young: global, and on active and passive coasts (by the A1 map);
  - **the macro**: the class fractions, the median, the land; the circumnavigation (C1's end land, as the selection).
- **Tier 3 (2 048², young)**: the transect without a clamp (the range from the data) for seed 18 under F178, for the seed whose
  longest cordillera is new under F178, and for F176's seed 18 at the fixed continent (F176's image, unclipped).

## The budget (announced; < 45 min wall per tier; resumable per world)

| tier | content | estimate |
|---|---|---|
| 0.1 | the guard through the checks | ~6 min |
| 1 | C1: 15 seeds × 5 variants, parallel | ~10 min |
| 2 | 256² young: 15 seeds × F178 and F176 | ~18 min |
| 3 | 2 048²: three transects | ~5 min |

**Total ~40 min wall.**

## The rules (written before)

- **Stops**:
  - the rampart back: ≥ 3 seeds above 70 % active, or W₁₀ at young ≥ 0.43 on a panel seed;
  - the selection fails, or the circumnavigation is lost, on ≥ 2 seeds;
  - any other closure change, or a solver;
  - D_gap as declared gives no cordillera → **no other value tried**; report it with the sensitivity at 3 and 5 cells.
- **Adoption** (seaward margin + option 5 + merge at collision + the Andean seed in the continent profile) is proposed only if
  every condition holds:
  - (i) the A1 median is in 10–50 % with both rules;
  - (ii) ≥ 70 % of the continents have an active margin;
  - (iii) ≥ 5 of the 15 seeds have an Andean cordillera, and the separated highland is never the only formation;
  - (iv) no rampart;
  - (v) **the author validates B3's maps and the formations' plate.**
- **Queued, out of this finding**: the rift strip that becomes a one-cell plate (F177 B4, B6).
- Forbidden: as in F165–F177. No constant calibrated; D_gap is the source's maximum.
