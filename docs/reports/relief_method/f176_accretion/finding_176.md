# Finding 176 — the accretion merges a pair only at collision, so that the subduction reaches the continents

**Status: built behind the continent profile's flag (one closure changed, as declared; no solver), measured by tiers; the success
holds on two of its three criteria — every continent of the 15 seeds has an Andean margin, the controls pass — but the median
active share at the end is 8.4 %, under the 10–50 % band; nothing committed without the feu vert.** Production unchanged: its
accretion rule is untouched, and C1's production hash is the F175 one (seeds 42, 8, the old témoin:
c69acd85a0f7f484, 3565285d65e4c341, fd23c309498a9761). No cache key touched; no ALGO_* bump.

- Declared before any measurement: `f176_declared.md` (the inventory, the rule, the sources). Predictions: `f176_predictions.md`.
- The bench: `crates/ymir-core/tests/f176_cascade.rs`, four resumable tiers; each world's lines in `results/`, assembled in
  `f176_bench_output.txt`. Images in `images/` (north up, checked).

## The budget (announced before running: ~22 min wall, ~50 min CPU)

| tier | content | wall |
|---|---|---|
| 1 | C1 only: the panel at the fixed continent (F175's plates) × F175 / F176 / F176 + B | 29 s |
| 2 | the screening: 15 seeds (the panel + 1, 2, 3, 4, 5, 6, 7, 9, 10, 11), the selection redone under F176 (parallel), C1 | 83 s |
| 3 | the physics level at 256²: the panel × F175 / F176 young and old, F176 + B young | 436 s |
| 4 | the cascade to 2 048²: seed 18 (the most Andean cells), young, F176 | 77 s |

Total **~10.5 min wall**.

## Part 0 — the inventory (details in `f176_declared.md`)

- **The accretion** (`closures/accretion/merge.rs`): `ConvergenceTracker::update` counts, per plate pair, convergent against
  divergent edges every step; a net-convergent step increments the counter, any other resets it; the pair merges at 50 (whole loser
  plate to the lower id, mass-weighted velocity). **The boundary's type is not read.** It was introduced in Issue #132 (Track D
  Stage E2) with a **collision** rationale (« ≈ 33 Ma … consistent with collision-zone timescales »).
- **The continent-ocean rule takes over** once a band is consumed (an oceanic convergent cell next to a continental cell of another
  plate is consumed beneath the continent) — observed here: the C-O consumption rises from 100–320 cells (F175) to 149–2 173 (F176)
  on the panel.
- **The terrestrial reference** (comparison, not a target): Bradley 2008 — passive margins 105 000 km, convergent plate boundaries
  53 000 km (p. 2), so active continental margins are **at most 34 %** of the continental margins' length (the convergent total
  includes intra-oceanic subduction). **No source in `docs/refs` splits the 53 000 km**; the author is asked for one (Leat & Larter
  2003 or Stern 2002 — to check).

## Part A — built

- **A1**: with `Margin::merge_at_collision`, `ConvergenceTracker::update_collision` counts **only continent-continent edges**; an
  O-O or O-C boundary never merges. The threshold (50) and the merge itself are unchanged. Through `run_with_closures_terms_flags`
  / `run_c1_sources_mapped_flags`; `Margin::seaward_oo_collision()` = the seaward margin + option 5 + the merge at collision.
- **A3, the controls** (permanent, `the_merge_at_collision_controls`): the flag off is F175 bit for bit (and production — the hash
  above); a world with every cell continental gives the same result with and without the flag; with the flag a real world keeps at
  least as many plates.
- **A4, the viz**: « Fusion seulement à la collision (C1) ».

## Tier 1 — C1, the panel at the fixed continent

| seed | plates at the end: F175 → F176 | active at the end: F175 → F176 | within 3 cells (F176) | Andean groups (F176) | C-O cells consumed | r* (F176) |
|---|---|---|---|---|---|---|
| 8 | 3 → **9** | 4.4 → **0.0 %** | 32.4 % | 0 | 128 → 149 | 4 |
| 18 | 2 → **14** | 2.3 → **15.3 %** | 38.1 % | 7 (194 km) | 135 → 2 173 | **0** |
| 19 | 2 → **13** | 0.0 → **8.5 %** | 30.5 % | 3 (81 km) | 100 → 1 234 | 4 |
| 23 | 2 → **10** | 0.0 → **5.5 %** | 26.9 % | 2 (56 km) | 100 → 200 | **0** |
| 29 | 3 → **10** | 0.0 → **10.2 %** | 38.4 % | 3 (112 km) | 320 → 629 | 3 |

- **The plates no longer collapse** (9–14 at the end against 2–3); the O-O consumption roughly doubles (22 000–34 000 cell-steps).
- **The coasts become active where the O-O consumption reaches them**: 4 of 5 continents with ≥ 1 Andean margin; seed 8 none (its
  bands are the overriding side — the trench stays at the band's outer edge, seen by the « within 3 cells » reading, 32 %).
- At the fixed continent the land grows on seeds 18 and 29 (+8 and +5 points: the arcs on the band's continental cells) and **r*
  falls to 0 on seeds 18 and 23** (still not wrapping); with the selection redone (tier 2) every r* is ≥ 1.
- **V1–V6**: V1 100 % of the initiations near a continent (21–41 per world); V2 100 % within 200 km; V4 by construction (the share
  decided on |Δage| < 0.01 rises to 30–52 %); V5 the arc mass conserved, but **one-cell fingers rise** (12–31 against 0–10 in F175);
  V6 above.

## Tier 2 — the screening (15 seeds, the selection redone under F176)

| seed | active at the end | within 3 cells | plates at the end | Andean groups (km) |
|---|---|---|---|---|
| 8 | 7.3 % | 30.9 % | 11 | 2 (44) |
| 18 | 15.8 % | 47.9 % | 14 | 5 (150) |
| 19 | 3.6 % | 36.0 % | 12 | 1 (44) |
| 23 | 2.5 % | 18.2 % | 10 | 1 (25) |
| 29 | 7.1 % | 19.3 % | 16 | 3 (88) |
| 1 | 10.5 % | 36.5 % | 11 | 3 (88) |
| 2 | 14.9 % | 37.1 % | 15 | 3 (100) |
| 3 | 8.9 % | 38.8 % | 13 | 3 (94) |
| 4 | 9.2 % | 24.5 % | 18 | 1 (31) |
| 5 | 9.1 % | 36.4 % | 12 | 3 (88) |
| 6 | 8.3 % | 29.6 % | 12 | 1 (19) |
| 7 | 8.4 % | 22.6 % | 24 | 5 (106) |
| 9 | 13.6 % | 34.3 % | 15 | 4 (169) |
| 10 | 7.0 % | 28.1 % | 13 | 3 (62) |
| 11 | 4.4 % | 22.6 % | 16 | 1 (25) |

- **The selection finds a continent on 15 / 15** (N 16 everywhere, land 50.0–57.5 %, r* 1–4); the circumnavigation holds on 15 /
  15 (r* 1–4 at the end).
- **Every continent (15 / 15) has ≥ 1 Andean margin** (1–5 groups, 19–169 km on the map).
- **The active share at the end** (3 × 3): median **8.4 %**, range 2.5–15.8 %; ≥ 10 % on 4 seeds; none above 70 %. The « within 3
  cells » reading (a trench within ~140 km real): median ~31 %, 18–48 % — the order of Bradley's ≤ 34 % bound.
- **The plates at the end**: 10–24 — above 16 on two seeds: without merges, C1's rifting splits are no longer compensated.

## Tier 3 — the physics level at 256² (the fixed continent)

| seed | mountain at young: F175 → F176 | W₁₀ young: F175 → F176 | W₁₀ active / passive (F176, young) |
|---|---|---|---|
| 8 | 71.5 → 53.2 % | 0.139 → 0.085 | — / 0.085 (no active coastal land) |
| 18 | 22.6 → 54.9 % | 0.044 → 0.275 | 0.416 / 0.230 (×1.8) |
| 19 | 40.3 → 36.4 % | 0.021 → 0.076 | 0.203 / 0.057 (×3.6) |
| 23 | 48.8 → 53.5 % | 0.039 → 0.207 | 0.000 / 0.215 (243 active cells, none steep) |
| 29 | 43.2 → 62.1 % | 0.105 → 0.238 | 0.201 / 0.245 (×0.8) |

- **No rampart back**: W₁₀ at young 0.08–0.28 (Part A's 0.43–0.67); rampart_W fires nowhere.
- **The Andean margins build mountains locally** (`maps256_young_F175_vs_F176.png`: white belts along some coasts of seeds 18, 23,
  29); the mountain share rises by 5–32 points there and falls 18 points on seed 8.
- **Active against passive**: ≥ 2 × on 1 seed (19); the active coastal land is small (0–1 779 cells), and part of the passive
  coast is steep from the band's own relief.
- **No island** at 256² from arcs (seed 19's islands are a land piece split by an arm, 6 000–8 000 km²).

## Tier 4 — the cascade to 2 048² (seed 18, young, F176)

The coast: dimension 0.978, length ratio 0.882 (smooth as before); no island; the transect (`transect_seed18.png`): the Andean
margin (red) carries a ~2.5 km coastal range then falls to −5 km; the passive margin (green) is low and falls more gently.

## Part B — emerged arcs (measured only)

**No effect**: islands 0 everywhere, the land within 0.03 % of F176. The arc load per oceanic cell is p99 0.22–0.48 and **max
0.26–0.59** in s units, under the **0.715** a cell needs to emerge under Airy (0.2 + arc ≥ 0.915); under 21 km of crust at most
(Christensen & Mooney's continental arcs: 38.7 ± 9.6 km, Table 3). The arc mass is spread over 324–628 oceanic cells
(`arc_distance` 3), too thin to emerge. A source on island-arc crust is still missing.

## The success and the stops

- **The stops**: the selection fails on 0 / 15; the circumnavigation is lost on 0 (with the selection redone); no seed above 70 %
  active → none fires. No other closure changed, no solver.
- **The success**: ≥ 70 % of the continents with an Andean margin ✓ (100 %); the controls ✓; **the median active share in 10–50 %
  ✗ (8.4 %)**. **Not met** (one criterion, narrowly): the three mechanisms do not enter the profile.

## Predictions

- **The reviewer's**: T1 held (9–14 plates on 5 / 5); T2 held (100 %); T3 refuted on the median (8.4 %), held on the ceiling (max
  15.8 %); T4 refuted (×2 on 1 seed); T5 held (15 / 15, circumnavigable); T6 refuted (no arc island).
- **Mine**: held: T1, T3's median in 8–25 % and the 70 % ceiling, T5, T6 refuted (the arc load too thin), the controls; refuted: T2
  « 40–70 % » (100 %), ≥ 10 % on 8–11 seeds (4), T4 on ≥ 3 (1), the mountain +0–10 points (−18 to +32), the « 1 criterion fails »
  being T2 (it was the median).

## Open, for the author (nothing chosen)

- (a) The success missed by the median active share (8.4 % against 10 %), while every continent has an Andean margin: whether the
  3 × 3 rule (a trench within ~47 km real) or the « within 3 cells » reading (~31 %, against Bradley's ≤ 34 % bound) is the right
  measure.
- (b) The active-margin source: a split of the 53 000 km of convergent boundaries into ocean-continent and intra-oceanic.
- (c) The plate count rising above 16 on two seeds (the rifting splits without merges).
- (d) The one-cell fingers rising (12–31) with the longer O-O consumption.
- (e) Emerged arcs need a concentrated arc load (the BFS spreads it over ~3 cells) — and a source on island-arc crust.
- Not in this brief, queued from F175's commit: the guard's stale downstream cache (F175 found the guard's references match a cache
  the current code no longer reproduces).

**At commit (2026-10-18), decided**: the author's words: « Feu vert pour commiter le F176 puis le F177 (config.rs et workspace.rs restent dehors ; le test ignoré de hd.rs peut entrer). »
- **F176 is committed behind its flag, off by default**: `Margin::merge_at_collision = false`, so production, F175 and every
  existing profile are unchanged bit for bit (the permanent control `the_merge_at_collision_controls`).
  `Margin::seaward_oo_collision()` is only used on request, with the viz box « Fusion seulement à la collision (C1) ».
- Variant B (`Margin::arcs_emerge`) is committed as measured: off by default, with no effect. F177 closed it on the sources (an
  island-arc crust of ~20 km is submerged).
- The measure of the active share, open at F176, was settled by F177's A1.
