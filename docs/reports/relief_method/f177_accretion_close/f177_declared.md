# Finding 177 — declared before any measurement (and before re-reading F176's data)

**This finding modifies no closure and adds no solver.** It measures, attributes and diagnoses. If a diagnosis shows a closure
defect, the round stops at the description (the mechanism in the code, the physical source of a possible correction); the author
writes the next finding. Deposition is F178. F176 is not committed (no feu vert in this brief): its code stays in the tree,
uncommitted, and F177 measures the configuration it built (`Margin::seaward_oo_collision()`).

## The author's criterion (the complement, verbatim; it goes into the ADR)

> « Oui c'est totalement accepté, si et seulement si ce n'est pas la seule formation qui persiste. Sur certaines seeds je dois
> récupérer des cordillères et autre formations. »

It concerns the highland separated from the continent by a strait (collision or back-arc basin).

## The unit convention (one, everywhere in F177)

**World km** (the cascade's world: a 400 km domain), with the terrestrial equivalent (× 7.5, the geographic scale ratio) in
parentheses:

| grid | cell, world km | terrestrial equivalent |
|---|---|---|
| C1, 64² | 6.25 km | 46.9 km |
| physics level, 256² | 1.5625 km | 11.7 km |
| 2 048² | 0.195 km | 1.46 km |

**F176's mixed units, corrected**: F176 wrote « a trench within ~47 km real » for the 3 × 3 rule — that is one C1 cell, **6.25
km world (46.9 km terrestrial)**; « within 3 cells (~140 km real) » is **18.75 km world (141 km terrestrial)**; the Andean
groups « 19–169 km on the map » were world km: **19–169 km world (141–1 268 km terrestrial)**.

## Part 0 — the guard (the « stale cache » queued since F175's commit)

**A reading of the code, before any run**:
- The guard looks its reference up **by digest**: `bench_guard::check_against` and `check_lakes_against`
  (`crates/ymir-core/src/tectonics_c1/bench_guard.rs`) return `NoReference` when no entry carries the run's digest, `Mismatch`
  only when an entry with the same digest has a different hash.
- `run_hd`'s digest is `eroded_key_full` (`crates/ymir-viz/src/bridge/c1/hd.rs:890`), chained onto `tectonic_key`, which folds
  `ALGO_TECTONICS` (`cached_product.rs:58–71`). **Bumping `ALGO_TECTONICS` changes every digest, so every state reads
  `NoReference`, whatever the field.**
- `f123_viz_guard` prints « = banc » for `Match` and **« ≠ / non gardé » for both `Mismatch` and `NoReference`**
  (`workspace.rs:7045`); the panic message carrying the variant was filtered out by `checks154.sh`'s grep. So F175's « all six
  fields and lakes different » (`checks_after_with_algo_bump.txt`) **is not established**: the output is what a bump gives with
  identical fields too.
- **The hypothesis to test**: there is no stale cache; F175 read the guard's `NoReference` as a mismatch.
- **The test (tier 0, ~9 min wall)**: the state « livré » recomputed **at unchanged digests with an empty cache** — the
  `ymir-viz` test binary run from an empty working directory (the cache is `.ymir_cache/` under the cwd), test
  `f177_fresh_guard` (test-only, in `hd.rs`; the author's `workspace.rs` untouched). Its field and lakes are compared with the
  references by the guard itself (`Match` / `Mismatch`, printed with the variant).
  - `Match` on both → no stale cache; the proposal: make `f123_viz_guard` print the variant (and not grep it away), and run
    the five other states fresh as a separate tier (~45 min) only if the author wants it.
  - `Mismatch` → find the first divergent stage (the tectonic state, then the eroded field, then the drainage and the lakes,
    comparing the fresh products with the cached ones) and the commit that introduced it (the cache entries' dates bound the
    range); **no reference is rewritten** without the author's feu vert.

## A1 — the active margin, without a parameter

The world at the end of C1 (64²); the measured configuration is F176's (`Margin::seaward_oo_collision()`).
- **A coast cell** of a continental plate P: a cell of continental type (`plate_type`) with ≥ 1 4-neighbour of oceanic type
  (F176's perimeter); P is its plate.
- **On a plate boundary**: the coast cell has a 4-neighbour of another plate → **F176's 3 × 3 rule** (active iff a convergent
  boundary cell lies in its 3 × 3).
- **Otherwise, the walk**: a breadth-first walk (4-connected, periodic) from the coast cell through **P's oceanic cells only**,
  until the first depth at which it reaches **a boundary cell of P** (an oceanic cell of P with a 4-neighbour of another plate).
  At that depth every boundary cell reached is classified by **the closures' own rules**:
  - **P overrides (active)**: the cell is an ocean-ocean contact of a subducting pair (`oo_contacts` with `OoTracker::is_active`,
    the tracker replayed by the observer as in F175) and is **not** the sinking side by option 5's rule (the older boundary cell
    sinks; ties: the velocity flag read before A7, then the lower plate id) — the partner sinks beneath P.
  - **P plunges (passive)**: the cell is convergent and touches a continental cell of another plate (the continent-ocean rule
    consumes it), or it is the sinking side of an O-O subducting pair.
  - **Passive**: divergent, transform, ambiguous, or convergent O-O on a pair whose cumulative convergence is still under
    D_init (no subduction yet).
  - **Ties at the first depth** (several boundary cells, different classes): active if at least one is active; the share of coast
    cells decided on a mixed tie is reported.
- **Enclosed**: the walk never reaches a boundary (an inner sea of P) → neither; excluded from the share and reported.
- **The distance coast → trench**: the walk's depth (cells × 6.25 km world) for the active cells; its distribution is reported.
- **The share** = active / (active + passive) over the coast cells of the continents (4-connected components of continental
  cells with ≥ 20 cells, F176's rule). **A continent has an active margin** if it carries ≥ 1 group of ≥ 3 8-connected
  A1-active coast cells (F176's group, without its t₀ condition); the F176 « Andean » reading (passive at t₀, active at the end)
  is reported beside it.

## A2 — the terrestrial reference (a comparison, never a target)

All four sources are in `docs/refs` (`leat2003.pdf`, `nishizawa2006.pdf`, `bird2003.pdf`, `Brad_08.pdf`).
- **Leat & Larter 2003** (Geol. Soc. London Spec. Publ. 219, 1–17), the exact wording: abstract (p. 1): intra-oceanic arcs
  « comprise some 40% of the subduction margins of the Earth »; introduction (p. 1): « The total length of margins where oceanic
  crust is being subducted (approximately the same as the total length of volcanic arcs) on Earth is about 43 500 km (von Huene
  & Scholl 1991) … Some 17 000 km, or nearly 40%, of the global length of volcanic arc is not situated on the continental
  margins, but is instead situated on oceanic crust »; the conclusions repeat « Nearly 40% … some 17 000 km of intra-oceanic
  arcs ». (17 000 / 43 500 = 39 %.)
- **Bradley 2008** (Earth-Sci. Rev. 91, 1–26), p. 2: « The present-day passive margins have an aggregate length of 105,000 km,
  even longer than the spreading ridges (65,000 km) or the convergent plate boundaries (53,000 km). » **No source and no
  definition is given for the 53,000 km** — whether collisions are included is not stated.
- **Bird 2003** (G³ 4(3), 1027, DOI 10.1029/2001GC000252), Table 3 (p. 46): subduction zones 51,310 km, continental convergent
  boundaries 23,003 km, oceanic convergent boundaries 17,449 km (all convergent: 91,762 km). Bradley's 53,000 km is 3 % from
  Bird's subduction zones and 42 % short of all three → **the reading (an inference, Bradley names no source): Bradley's total
  is the subduction zones, collisions excluded.**
- **The recalculation**: active continental margins = 60 % of the subduction length (Leat & Larter's ~40 % intra-oceanic):
  0.6 × 53,000 = 31,800 km → **31,800 / (31,800 + 105,000) = 23.2 %** of the continental margins (active + passive). With Leat &
  Larter's own totals (43,500 − 17,000 = 26,500 km): 26,500 / 131,500 = 20.2 %; with Bird's subduction zones: 0.6 × 51,310 =
  30,786 km → 22.7 %. **~20–23 %.** (The denominator omits transform margins and collision zones.)
- **Bird 2003's plates** (for B5): 52 plates, 14 large (abstract, p. 1); areas in Table 1 (p. 6); between 0.002 and 1 steradian
  the cumulative count ≈ 7 × (area in sr)^−1/3 (p. 43).
- **Variant B closed by the sources** (arithmetic on the declared Airy mapping, no measurement): Nishizawa et al. 2006 (Earth
  Planets Space 58, e33–e36), p. e33: the mid Izu-Ogasawara intra-oceanic arc crust is « about 20 km » (5 + 5 + 10 km), thinning
  to ~6 km within 30 km; Leat & Larter p. 2: intra-oceanic arcs are « typically mostly submerged below sea level, sometimes with
  only the tops of the largest volcanoes » emergent. Under F167's Airy (35 km per s, ρc 2 800, ρm 3 300, ρw 1 030, s = 1 →
  +0.45 km), 20 km of crust is s = 0.571 → −2.6 km; emergence needs s ≥ 0.915 (32 km). **A submerged arc is what the sources
  describe**: F176's « no arc island » is consistent with them, and variant B is closed (no emerged-arc rule is proposed).

## Part B — the measures (no closure modified)

- **B1, the share per A1**: the 15 screening seeds of F176 (C1 only, F176's selection redone: `select_continent_default_margin`
  with `Margin::seaward_oo_collision()`); per seed the share, the median over 15, the share of continents with an active margin;
  the A1 map on the panel (8, 18, 19, 23, 29; the selection's continents).
- **B2, the uplift at the coasts**: the tectonic source that drives the cascade, Δs_src = Σ over the run of every C1 term except
  the erosion (= s_src,end − s₀), turned into metres by the declared Airy mapping (Δh = e(s₀ + Δs) − e(s₀), continental cells).
  **The band**: the continental cells within C1's own Davis-Suppe wedge reach (`DavisSuppeParams::max_distance` scaled to 64²,
  an existing constant) of a coast cell, each assigned to its nearest coast cell (a 4-connected walk through continental cells)
  and typed by that cell's A1 class. Per seed: the mean Δh on the active band and on the passive band, and the share of the
  active band's Δs coming from the continent-ocean wedge (B3's attribution). **A defect is declared if the active band has no
  continent-ocean wedge uplift (mean ≤ 0) on a seed with active coasts**: the round stops and describes it.
- **B3, the attribution** (256², young, the panel at the fixed continent, F175 against F176; and the 15 screening seeds under
  F176's selection, young): the origin of each C1 cell's Σ Δs > 0:
  - **Andean (continent-ocean subduction)**: Davis-Suppe's increments on cells whose wedge comes from an O-C seed (C1's own
    typed wedge `oc_wedge`, recomputed by the observer from the state each step starts from) + the arc deposits whose nearest
    consumed cell of that step is a continent-ocean consumption (an oceanic cell touching a continental cell of another plate);
  - **collision (continent-continent)**: Davis-Suppe's increments on cells whose wedge comes from a non-O-C seed (C-C, including
    the band's continental cells; the O-O seeds are cut by option 5);
  - **ocean-ocean arcs**: the arc deposits whose nearest consumed cell is an O-O consumption (on the band's continental cells or
    on ocean);
  - **other**: none of the three positive.
  A mountain cell at 256² (the physics level's Kapos classifier, `MOUNTAIN`) takes its C1 cell's origin of largest Σ Δs > 0
  (other if none). **The gain** per origin = mountain cells (F176) − mountain cells (F175), in km² world. A map per seed: the
  relief shaded, the mountain cells coloured by origin.
- **The formations** (the complement, the 15 seeds, 256² young, no new parameter):
  - **Andean cordillera**: an 8-connected group of mountain cells of Andean origin that touches an A1-active coast (one of its
    cells lies in, or 8-touches, a 256² cell of an A1-active 64² coast cell); its length along the coast = the number of distinct
    A1-active coast cells it touches × 6.25 km world;
  - **collisional or separated highland**: the mountain cells of collision origin; « separated » when their 256² land component is
    not the largest one, the nearest B4 channel's class given;
  - **arc mountain on the continental band**: mountain cells of O-O-arc origin;
  - **other** (including the mountain inherited from the interior).
  A table per seed (area per formation, km² world; the cordilleras' lengths) and a plate of the 15 maps coloured by formation.
- **B4, the channels**: the sea at the end of C1 = the declared Airy mapping of s_src,end ≤ 0 (the cascade's coast at B0, F174);
  **a sea arm ≤ 2 cells wide** = the sea cells a 3 × 3 morphological opening removes (periodic); a channel = an 8-connected group
  of them; its length = its geodesic diameter (8-connected, two sweeps) × 6.25 km world; **one-cell-wide** = sea cells with land
  on both sides along one axis. The origin per cell, in this order:
  - (i) **rift**: the cell's cumulative rift-thinning Δs < 0 (the strip thinned by `apply_rifting_thinning`);
  - (ii) **ocean trapped by a merge**: oceanic at the end, not a boundary cell at the end, and an oceanic cell of either plate of a
    pair that merged during the run (merges detected by the observer: a plate id vanishing whole into another);
  - (iii) **the seaward band's remnant**: a cell whose plate or type differs between the profile's initial state with and without
    the seaward margin;
  - (iv) **other**.
  A channel takes its cells' majority class (ties: the order above). The number, the lengths, a map. **The straits of seeds 19
  (line 3 of F176's young montage), 18 and 29 first**: the channels that cross them are named.
- **B5, the plates**: the size distribution at the end (cells); the counts < 1 %, 1–5 % and > 5 % of the surface (< 41, 41–205,
  > 205 cells of 4 096); Bird 2003's Table 1 in the same classes (4π sr; 1 % = 0.126 sr); not a target.
- **B6, the rifting**: an inventory (file:line), its trigger, rate and source, and whether it was set assuming the merges
  compensate (the Track D report and the params' docs).
- **B7, the transect of seed 18**: a reading of the code first — `transect_png` clamps the altitude to [−5 500, 2 500] m
  (`f176_cascade.rs:1730`, `v.clamp(-5500.0, 2500.0)`), so a top at exactly 2.5 km is the plot's clip; the image is redone
  without a clamp (the vertical range from the data) and the true maximum along the transect is reported (and whether a run of
  equal values — a lake or a ceiling — exists).

## The budget (announced; < 45 min wall; resumable per world)

| tier | content | estimate |
|---|---|---|
| 0 | the guard: « livré » fresh at 8 192² (one `run_hd` from an empty cache) | ~9 min wall, ~1 h CPU |
| 1 | C1: the panel (fixed continent and selection) + the 10 screening seeds, the selection redone (parallel): A1, B1, B2, B4, B5 | ~2–3 min |
| 2 | 256² young: the panel × F175 / F176 at the fixed continent (B3) | ~3 min |
| 3 | 256² young: the 15 seeds under F176's selection (the formations) | ~5 min |
| 4 | 2 048² young, seed 18 (B7) | ~1.5 min |

**Total ~20 min wall.** Each world's lines are written to `results/` as soon as they are ready; a rerun skips them.

## The rules (written before)

- No closure modified, no solver. If B2, B4 or B6 shows a defect: stop, describe it (the mechanism in the code, the physical
  source of the correction); the author writes the next finding.
- **Adoption** (seaward margin + option 5 + merge at collision in the continent profile) is proposed only if: (i) the A1 median
  is in 10–50 %; (ii) ≥ 70 % of the continents have an active margin per A1; (iii) B4 classes no one-cell channel as an artefact
  without a physical explanation; (iv) **the author validates B3's and B4's maps**; (v) ≥ 5 of the 15 seeds have ≥ 1 Andean
  cordillera, and no seed has the collisional or separated highland as its only mountain formation. If (v) fails because the
  continent-ocean subduction builds no relief (B2): stop and report the mechanism. No constant is tuned to make cordilleras
  appear. **Nothing is adopted before the author's look.**
- Forbidden: as in F165–F176. No constant calibrated.
