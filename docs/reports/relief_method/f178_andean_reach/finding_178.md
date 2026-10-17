# Finding 178 — the Andean seed at the arc-trench distance: the wedge reaches every active coast, but it lifts the seaward band's submerged shelf and the circumnavigation is lost on 4 of 15 seeds — the declared stop fires after tier 1

**Status: built behind the continent profile's flag (one closure changed: where Davis-Suppe's continent-ocean wedge starts; no
solver), measured on C1; the circumnavigation stop fired, so tiers 2–3 (256², 2 048²) were not run; nothing adopted, nothing
committed without the feu vert.** Production unchanged: C1's production hash is F175's (`tests/f175_c1_hash.rs`); no cache key,
no `ALGO_*`.

- Declared before any measurement: `f178_declared.md` (D_gap from Syracuse & Abers 2006, p. 10). Predictions: `f178_predictions.md`.
- The bench: `crates/ymir-core/tests/f178_cascade.rs` (F177's helpers; the observer replays the new seeds). Results in
  `results/`, assembled in `f178_bench_output.txt`. Units: world km, terrestrial (× 7.5) in parentheses; a C1 cell is 6.25 km
  (46.9 km).

## Part 0.1 — the guard

- **The exact verdicts at HEAD (4593d8b)**: the six states read **`Match` on the field and on the lakes** (`checks_before.txt`).
  None reads `NoReference`, so there was nothing to re-register, and the 45 min fresh tier was not needed.
- **Applied**:
  - `f123_viz_guard` prints « = banc » for `Match`, « ≠ (banc …, viz …) » for `Mismatch` and « non gardé (aucune référence sous
    cette empreinte ; viz …) » for `NoReference` (the test module of `workspace.rs`, the author's file: at commit only this hunk
    would be staged);
  - the checks keep the per-state verdict line and the panic line.
- **Why F175 bumped `ALGO_TECTONICS`** (`f175_declared.md:15`): « the time loop's code changes even though the off path is bit
  for bit; the bump makes the cached tectonics recompute, so the guard 6/6 exercises the new code instead of reading the cache »
  — the « bump on any code change » rule, applied to a change whose output was unchanged. It moved the digest, hence F175's
  `NoReference`.

## Part 0.2 — the sources and D_gap (declared before measuring)

- **Syracuse & Abers 2006**:
  - p. 10, §26: « Horizontal distances from volcano to trench range from 85–90 km in Scotia and Vanuatu to 470 km in Alaska,
    with the middle 50% of arc-trench distances between 180 km and 275 km »;
  - H (the slab depth beneath the arc): 72–173 km for 500 km segments, an average of 105 km and a median of 100 km (p. 1, p. 10);
  - **no median of the arc-trench distance is given**.
- **England et al. 2004**: depths of 65–130 km (p. 377), and per-arc depth and dip in Table 1 (p. 379). **No arc-trench
  distance**, so it was not used for D_gap (D / tan δ would assume a planar slab from the trench).
- **D_gap = 470 km ÷ 7.5 = 62.7 km world = 10 cells** (one value, no tuning).
- The sensitivity at the middle-50 % bounds: 180 km → 3 cells, 275 km → 5 cells.
- A2's band: walk depths 4–5.

## Part A — built

- **A1**: `tectonics_c1/andean_seed.rs`.
  - `overriding_trench_cells` gives the overriding side of the active O-O trenches, by option 5's own rule;
  - `andean_seeds` gives P's continental cells first reached by a walk through P's cells from those trench cells, within D_gap.
    They join the upper-plate and O-C seed masks in C1's time loop after A7's cut, behind `C1Ext::andean_gap`
    (`Margin::andean_gap`, `Margin::seaward_oo_collision_andean()`, `ANDEAN_GAP_CELLS = 10`);
  - Davis-Suppe's formula is unchanged.
- **A2** (measured only): `OoSubductionParams::arc_band` places the O-O arc at walk depths 4–5 (`distribute_arc_mass_band`).
- **A3** (permanent, `the_andean_seed_controls`): the seed off is F176 bit for bit; a reach of 0 and an all-continental world
  are identical with and without the seed; at D_gap the seed acts. `f175_c1_hash` passes.
- **A4**: the viz box « Graine andine à distance arc-fosse (C1) ».

## Tier 1 — C1, the 15 seeds of F176's selection

| | F176 | **F178 (D_gap 10)** | F178 + A2 | reach 3 | reach 5 |
|---|---|---|---|---|---|
| A1 median, « any » | 27.4 % | 25.0 % | 25.0 % | 23.9 % | 24.9 % |
| A1 median, « majority » | 17.9 % | 18.3 % | 18.3 % | 17.9 % | 18.3 % |
| continents with an active margin | 15 / 15 | 15 / 15 | 15 / 15 | 15 / 15 | 15 / 15 |
| seeds with a positive C-O wedge on the active 0–6 cells | 12 | **15** | 15 | 15 | 15 |
| active ≥ 2 × passive on 0–6 cells | 4 | 4 | 4 | 15 | 11 |
| **circumnavigation lost** | 0 | **4 (seeds 23, 3, 4, 10)** | 4 | 1 | 3 |

- **The net tectonic uplift on the active 0–6 cells** (Δh_src): F176 0.4–4.7 km; **F178 12–22 km**. On the passive 0–6
  cells: F176 0.8–3.2 km; F178 7.6–15.4 km.
- **The land at the end of C1**: F176 50–57.5 %; **F178 55–70 %**; r* falls to 0 on 6 more seeds.
- **A2**: it changes the active coasts' uplift by −0.4 to +7.5 %, so T7 held.
- **The trench segments** (8-connected): medians 6–31 km world (1–5 cells), the same under F176 and F178.
- **No seed above 70 % active**; the A1 shares barely move, since the seed changes no kinematics.

## The stop, and its mechanism (measured; `f178_stop_mechanism`, seeds 23, 3, 4, 10, 18, 8)

**The circumnavigation rule fires (4 seeds ≥ 2).** Tiers 2–3 were not run. The measure explains it, in two parts.

1. **The « continental front » is the seaward band's submerged shelf, not the shoreline.**
   - The band's continental cells (the crust above 10 km, F174 B′) are where a walk from the trench first meets continental
     type. Their s at t₀ is 0.53–0.60 (19–21 km of crust), under the sea level (0.915).
   - F178 adds +233 to +622 cells of land over F176 (6–15 % of the map). Of that new land, **97–100 % lies in the band at
     t₀**, 74–88 % is of continental type, and 99–100 % received the C-O wedge.
   - So the shelf rises, the straits close, and the circumnavigation is lost.
   - The A1 coast (the type boundary) and the shoreline are different lines under the seaward margin. A1 and A1's seed both use
     the type boundary.
2. **The amplitude piles up in the cascade's source.** On the active 0–3 cells, mean Σ Δs over the run:

   | | Davis-Suppe | C1's erosion | C1's s at the end | the cascade's s_src |
   |---|---|---|---|---|
   | F176 | +0.13 to +1.09 | −0.47 to −1.36 | 0.19–0.41 | 0.67–1.69 |
   | F178 | **+4.4 to +5.8** | **−3.8 to −4.9** | 1.03–1.65 | **4.9–6.4** |

   - Davis-Suppe refills, at every step, what C1's erosion removes from the seeded cells: the target at d = 0 is h_max = 2.5
     (the declared consequence: a front cell is not a boundary cell), with the O-C rate boost × 6.
   - The cascade's source excludes C1's erosion (F166's « sources » mode, s_src = s − the erosion term), so it receives the
     whole pumped wedge: 4.9–6.4 in s, 170–220 km of crust.
   - F176's own O-C seeds show the same pumping, ten times weaker.

**What this says about the closure, for the author's next finding** (nothing changed here):
- The seed's reach is not at fault: the trenches are at 1–8 cells and D_gap reaches them.
- What is at fault, measured:
  - (a) under the seaward margin, the front the walk finds is the submerged shelf. If the wedge should start on the emerged
    continent, the front would be the first cells above sea level (or of thick crust), not the first of continental type. This
    is a choice for the author;
  - (b) the cascade's « sources » mode adds up Davis-Suppe increments that C1's erosion was balancing. Wherever a wedge seed sits
    on an eroding cell, the source runs away. This is also F176's case, weaker.

## Predictions

**The reviewer's**:

| | prediction | outcome |
|---|---|---|
| T1 | C-O wedge positive on ≥ 14 / 15 | held (15) |
| T2 | ≥ 8 / 15 seeds with an Andean cordillera | not measured (stop) |
| T3 | active ≥ 2 × passive on 0–6 cells on ≥ 10 / 15 | refuted (4) |
| T4 | W₁₀ young < 0.30 on the panel; no seed above 70 % active | W₁₀ not measured; no seed above 70 % (held) |
| T5 | condition (v) holds | not measured |
| T6 | majority median in 15–25 % | held (18.3 %) |
| T7 | A2 changes the active uplift by < 10 % | held (−0.4 to +7.5 %) |

**Mine**:
- held: 0.1 all `Match`; T1; T6 (17–24 %); T7; no seed above 70 %.
- refuted:
  - T3 « 9–12 seeds » (4);
  - « the land within ± 2 points » (+4 to +16);
  - « the circumnavigation holds on 15 / 15 » (lost on 4);
  - the trench segments « median 3–6 cells » (1–5).

## Open, for the author (nothing chosen)

- (a) The front: the first continental-type cells (the shelf, as declared) or the first emerged cells?
- (b) The cascade's source and C1's erosion: should the « sources » mode keep the part of Davis-Suppe that C1's erosion
  balanced? This affects F176's O-C seeds too.
- (c) The seed at d = 0 receives h_max, whereas a boundary seed is skipped. Should the new seed be skipped like a boundary
  cell?
- (d) Tiers 2–3 (~25 min) were not run; they are resumable.
- Queued (the author's order): the rift strip that becomes a plate.

**At commit (2026-10-18), decided**: the author's words: « Feu vert pour commiter le F178, drapeau coupé. Entrent : le banc, la graine derrière le drapeau, le contrôle, et le seul bloc de workspace.rs qui touche l'affichage de la garde. Entre aussi la correction de f175_cascade.rs (arc_band: None). config.rs et le reste de workspace.rs restent dehors. »
- **F178 is committed with its flag off**: `Margin::andean_gap = None` everywhere by default, so production, F176 and F177 are
  unchanged bit for bit (`the_andean_seed_controls`; production's C1 hash). `Margin::seaward_oo_collision_andean()` and the
  viz box « Graine andine à distance arc-fosse (C1) » are only used on request.
- Committed:
  - the bench (`f178_cascade.rs`);
  - the seed behind its flag (`andean_seed.rs`, `C1Ext`, `run_with_closures_terms_ext`, `run_c1_sources_mapped_ext`);
  - variant A2's band, measured only and off by default;
  - the permanent control, and the viz box;
  - the single hunk of `workspace.rs` that prints the guard's exact verdict;
  - the fix of `f175_cascade.rs` (`arc_band: None`).
  The rest of `workspace.rs` and `config.rs` stay out.
- **The reviewer's error, recorded (the author's; the same family as F174 B′)**: the reviewer defined the seed's front as « the
  first continental cells reached » without checking that, under the seaward margin, these are the submerged shelf's cells
  (s ≈ 0.55). Nor did the reviewer check how Davis-Suppe treats a seed at d = 0: a boundary seed is skipped, while a front cell
  receives h_max.
- **Still open** (nothing chosen): the front (continental type or emerged); the sources mode against C1's erosion (it touches
  F176 too); the seed at d = 0; tiers 2–3 on request.
