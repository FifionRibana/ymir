# Finding 167 — the thickness → altitude mapping: physical isostasy and an absolute sea level (cascade only)

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ». **Production, C1 and the
guard are unchanged** (the guard 6 / 6, field and lakes, before and after: `checks_before.txt`, `checks_after.txt`).

- Declared before any measurement: `f167_declared.md`, with the three sources read and cited by page, and one
  amendment at its control (Hovius per basin).
- Predictions: `f167_predictions.md`.
- Bench: `f167_bench_output.txt` (825 s wall); the inventory: `f167_inventory_output.txt`.

**The verdict: the « mountains still above 37.2 % » stop fires on 4 / 4 seeds.** Decomposed, stopped, no new
mechanism.
- **The Airy mapping with an absolute sea level does what it was built for.**
  - The ocean–land coupling is gone: « sources » and « sources − EH » give the same land on three seeds (0 cells); on
    seed 9 they differ by 0.28 %, where the EH sink touched a few continental cells.
  - The land is kept (1–2 % lost against F164).
  - The peaks follow the tectonics (2.7–4.2 km) with T frozen at 5.62 Myr.
- **But the mountains RISE**: 53–78 % of the land, against F166's 28–64 % and Europe's 24.8 %. Plain + plateau fall to
  0–6 %.
  - **The decomposition**: 50–81 % of the mountain cells are Davis-Suppe zones (Σ Δs > 0.05), 5–10 % cratons, 0–2 %
    arcs, and 14–41 % the rest.
- **The reading**:
  - Airy's land slope is 5.30 km of altitude per unit S̃, against ~3.5 km in C1's ramp (which runs from its datum near
    S̃ ≈ 0.3 to 5.65 km at S̃ = 2).
  - So the same tectonic thickening makes ~50 % more relief. C1's Davis-Suppe thickening covers most of each
    continent, which is therefore classed mountain.
  - The middle of C1's ramp was too high (the reviewer's diagnosis holds: S̃ = 1 at 1.8 km), but lowering it while
    keeping the top raises the gradient, and the wedge's spread then dominates.

## Partie 0

1. **F166 committed** (`4ffec0d`), the author's edits left out.
   - Its ADR records:
     - the new order (F167 mapping, F168 sea-level rise, F169 heterogeneous K);
     - the mapping corrected in the cascade only;
     - the facets / walls (the talus threshold, a hypothesis) and the map-border rim queued;
     - the author's remark, verbatim.
2. **Checks**:
   - before: the guard 6 / 6, lib 641, viz 33;
   - after: the guard 6 / 6 identical (field and lakes = banc on all six states), lib 641 → 643 (the mapping's
     negative control, Hovius's synthetic range), viz 33.
   - **No ALGO_* bump**: nothing cached changes; the mapping lives in the cascade.

## The sources (read; the details and the page citations in `f167_declared.md`)

- **Whitehead & Clift 2009, p. 3**:
  - ~38 km average stable crust, ~35 km at sea level (seismic), ρc 2 800, ρm 3 300 kg/m³;
  - 42.5 km → ~1.5 km; 58 km → ~3.8 km;
  - p. 7: ρw 1 030.
- **Lachenbruch & Morgan 1990, p. 42, eq. (4)–(7)**: the ridge at E = −2.5 km as the reference column, with ρa 3.2,
  ρl 2.8, L 5.5 km, ρw 1.0.
- **Talling et al. 1997, p. 276**: Hovius's R = W/S between 1.91 and 2.23 over 11 orogens.

**The reviewer's citations, checked (errors recorded, the papers prevail)**:
1. **Whitehead's numbers are not one Airy line.** With slope 0.1515, 42.5 → 1.5 km puts sea level at 32.6 km, so 35 km
   → +0.36 km and 38 km → +0.82 km. « 35 km at sea level » is the seismic statement, not the Airy line, and « 38 km
   → ~+0.45 km » follows neither.
   - The paper cannot settle the anchor, so **the brief's central +0.45 km** was taken, with ±0.2 km measured.
2. **Lachenbruch & Morgan's −2.5 km ridge is a LITHOSPHERE column** over a 3.2 asthenosphere, not crust over a 3.3
   mantle. A crust-only Airy puts their 5.5 km ridge column at −5.84 km. It cannot anchor the continents' function and
   was reported as a check.
3. The paper is Whitehead **and Clift** 2009.
4. Hovius (via Talling): **confirmed**.

## A — the inventory (in the declaration)

- **C1's mapping**:
  - a land ramp from a field-percentile datum (h_sea = h_min + 0.4 · (p92 − h_min)) to 5 650 m at S̃ = 2;
  - a σ = 0.5 blur;
  - Stein-Stein bathymetry (by age) over the oceanic cells;
  - a normalised ceiling at 5 650 m.
  - On the témoin's final datum the ramp **crosses sea level at S̃ ≈ 0.3** (~10 km of crust; `f167_mapping.png`).
- **`age` is not a clock** (never incremented in time), so **no thermal subsidence** was added.
- **What reads h in C1**: only the erosion closure, for its slopes. Davis-Suppe, the EH sink, the subduction, the rift
  and the transport read s.
  - So the cascade-only correction leaves out only C1's own erosion, which the sources' drive already excludes. The
    drive is not approximated.
- **The inventory of s** (`f167_inventory_output.txt`): at the end of C1, **76–186 oceanic cells per seed have
  S̃ > 0.915** (up to 2.18). These are C1's documented « phantom oceanic advective spike », which C1 never maps.
  - **Main**: the oceanic-plate cells are read at C1's oceanic thickness S̃ = 0.2.
  - **Variant « own s »**: every cell by its own s.

## B — the mapping (`history::Airy`, `Mapping::Airy`, `run_c1_sources_mapped`)

- **e = 5.303 · S̃ − 4.853 km** above sea level (ρc 2 800, ρm 3 300, 35 km per S̃, S̃ = 1 → +0.45 km), water-loaded
  below it (× 3300/2270).
  - One continuous monotone function; the sea at S̃ = 0.915; **sea level 0 m absolute**.
  - S̃ = 0.2 → −5.51 km; S̃ = 1.25 (the cratons) → +1.78 km; **S̃ = 2 → +5.75 km** (the high check: inside 5–6.5 km;
    the normalised field clamps at 5.65).
- **The negative control** (permanent): `Mapping::C1` gives F166's record bit for bit.

## C — the measures

### 256² (Europe: plain 51.6 · plateau 4.7 · hill 18.9 · mountain 24.8)

| seed | run | plain | plateau | hill | mountain | peak | cratons p50 | lost / gained vs F164 | lakes |
|---|---|---|---|---|---|---|---|---|---|
| témoin | F166 | 1.0 | 0.0 | 37.7 | 61.3 | 2 704 | 832 m | 0.0 / 16.4 % | 2 / 88 km² |
| | **F167** | 0.0 | 0.0 | 25.6 | **74.4** | 2 670 | 1 165 m | 2.0 / 11.9 % | 3 / 51 |
| | own s | 0.0 | 0.1 | 21.8 | 78.0 | 4 022 | 1 331 m | 1.7 / **42.6** % | 7 / 1 570 |
| | anchor 0.25 / 0.65 | 0.0 / 0.0 | 0.1 / 0.0 | 27.0 / 23.9 | 72.9 / 76.1 | 2 643 / 2 916 | | 2.3 / 1.9 % lost | |
| 42 | F166 | 10.7 | 0.0 | 60.7 | 28.5 | 3 262 | 70 m | 4.5 / 20.2 | 4 / 208 |
| | **F167** | 0.9 | 0.4 | 45.3 | **53.3** | 3 004 | 367 m | 1.9 / 18.8 | 2 / 24 |
| | anchor 0.25 / 0.65 | 2.0 / 0.4 | 0.6 / 0.1 | 49.1 / 42.5 | 48.3 / 57.0 | 2 907 / 3 183 | | 2.4 / 1.6 | |
| 1 | F166 | 1.6 | 0.0 | 34.5 | 63.9 | 4 301 | 455 m | 0.0 / 17.1 | 5 / 876 |
| | **F167** | 0.0 | 0.0 | 22.4 | **77.6** | 4 195 | 921 m | 1.0 / 12.6 | 1 / 27 |
| | anchor 0.25 / 0.65 | 0.1 / 0.0 | 0.0 / 0.2 | 23.5 / 20.9 | 76.5 / 78.8 | 4 136 / 4 437 | | 1.3 / 0.8 | |
| 9 | F166 | 4.8 | 0.0 | 39.3 | 55.9 | 5 365 | 470 m | 3.8 / 12.9 | 17 / 4 026 |
| | **F167** | 0.3 | **5.6** | 23.1 | **70.9** | 4 182 | 858 m | 1.7 / 12.9 | 7 / 518 |
| | anchor 0.25 / 0.65 | 0.7 / 0.1 | 9.4 / 5.3 | 28.4 / 20.9 | 61.5 / 73.7 | 3 950 / 4 437 | | 1.9 / 1.4 | |

- **The decoupling test** (« sources » against « sources − EH »):
  - main: 0 / 0 / 0 / 59 cells differ (0 / 0 / 0 / 0.28 %);
  - « own s »: 4.3–14.5 %. The EH sink caps the phantom islands; without it they grow.
- **The phantom islands under « own s »**: 3 833–5 837 land cells on oceanic plates at 256² (against 408–713 on the
  main), +32 to +62 % of land. **The main's substitution is what keeps the ocean an ocean.**
- **The anchor ±0.2 km** moves the mountain share by 1–10 points (−9.4 at most, seed 9 at 0.25 km) and the land by
  < 1 %. It is not the lever.
- The interior at ≤ 0 m is ≤ 0.06 % on the main.

### 2 048² (the frozen cascade), against Corsica at 195 m

| seed | F167 whole | F167 mountain class | the slope guard (Corsica 0.318 / 0.550) |
|---|---|---|---|
| témoin | 15 / 19 | 15 / 19 (facets ×1.7) | 2 870 m · 0.324 / 0.500 |
| 42 | 15 / 19 | 14 / 19 (facets pass) | 3 079 m · 0.291 / 0.483 |
| 1 | 14 / 19 (facets, walls) | 16 / 19 (facets ×2.5, walls ×2.0) | 4 545 m · 0.385 / 0.522 |
| 9 | 17 / 19 | 15 / 19 (facets ×1.6) | 4 171 m · 0.305 / 0.498 |

- **The slopes stay Corsica's.**
- **The facets fail in the mountain class on 3 / 4 seeds** (seed 1 the walls too).
- **The coast is degraded against F166 on seed 42** (the length ratio, the block index).

### The regularity baseline (for F169), against Corsica

| | CV_λ (391 / 195 m) | CV_L | confluence θ mean, σ | R2_g (parallelism) | Hovius R (A ≥ 25 km²) |
|---|---|---|---|---|---|
| **Corsica** | 0.61 / 0.68 | 0.87 / 0.92 | 62–71°, **σ 27–30°** | 0.48 / 0.47 | **2.34** (34 basins; 2.33 at 10 km², 2.55 at 100 km²) |
| F164 | 0.55–0.66 / 0.59–0.64 | 0.82–0.98 / 0.82–0.97 | 60–69°, σ 22–24° | **0.57–0.63** / 0.45–0.49 | 1.47–1.95 |
| F166 | 0.60–0.68 / 0.61–0.69 | **0.70–1.01** / 0.72–1.06 | 62–71°, σ 22–24° | **0.53–0.65** / 0.46–0.49 | 1.91–2.42 |
| F167 | 0.60–0.65 / 0.61–0.64 | **0.73–0.77** / 0.76–0.78 | 64–72°, σ 22–25° | **0.50–0.56** / 0.45–0.46 | 2.14–2.53 |

**Where we are more regular than Corsica, and at which cell**:
- **The confluence angles at both cells**: σ_θ 22–25° against 27–30°, i.e. 15–25 % narrower. The means are Corsica's.
- **The gullies' parallelism at 391 m**: R2_g 0.50–0.65 against 0.48. At 195 m it is Corsica's.
- **The tributaries' lengths under the sources' drive** (F166, F167): CV_L 0.73–0.78 against 0.87–0.92, i.e. 15–20 %
  more uniform, at both cells. F164 was Corsica's.
- **Not more regular**:
  - the valley spacing (CV_λ ≈ Corsica's);
  - the basins' shape: Hovius's R is 2.1–2.5 for F167 against Corsica's 2.3–2.5, both near Hovius's 1.91–2.23.
    Corsica sits a little above the band; its fronts are mostly coastal.

## R

- **Macro success: 0 / 4** (the mountains 53–78 %, plain + plateau 0–6 %; the land lost 1–2 % passes).
- **« The mountains still above 37.2 % » (4 / 4): decomposed above, and stopped.** The land-lost stop does not fire.

## The reading (hypotheses to check)

- **The mountains are Davis-Suppe's**: 50–81 % of the mountain cells.
  - C1's wedge (`DavisSuppeParams`: l_taper 4 cells = 25 km, l_decay 6 cells = 37 km, `max_distance` 30 cells =
    187 km) thickens most of each small continent.
  - Under a physical mapping its Σ Δs becomes kilometres of relief everywhere.
  - In C1's ramp the gradient was ~2/3 of Airy's, which hid part of it.
- **Lowering the middle of the ramp (S̃ = 1 from 1.8 to 0.45 km) did not make plains.** The land never sits at S̃ ≈ 1:
  the sources' s_src is the initial crust plus the wedge, so the continents are thickened almost everywhere.
- **What C1's work (b) would have to give**: a wedge confined to the margins, collisions, and normal crust left normal.
  Only then can a physical mapping give plains. Not done here.

## Images (north up; Copernicus notice in `images/NOTICE.md`)

- `f167_seeds_256.png`, `f167_seeds_2048.png`, `f167_classes_seeds.png`: rows = the four seeds, columns F166 | F167.
  The ocean is darker under Airy (−5.5 km).
- `f167_mapping.png`: h(S̃) for S̃ 0–2.5. Red is C1's land ramp at the témoin's final datum (it crosses sea level near
  S̃ = 0.3); blue is Airy (sea at 0.915, the water-loaded slope below). Grey is the 5.65 km ceiling.
- `f167_coast_25km.png`: Corsica | F166 | F167.

## Predictions

**The reviewer's**:
- **R1** (identical within 0.1 % on 4 / 4): refuted by one seed (3 / 4 exact; seed 9 0.28 %).
- **R2**: refuted (0 / 4).
- **R3**: refuted (0 / 4).
- **R4**: held (1–2 % lost on 4 / 4).
- **R5** (peaks −15 % at most): held on 3 / 4 (−1 to −8 %); refuted on seed 9 (−22 %).
- **R6**:
  - the cratons > 1 km on 1 / 4 (the témoin 1 165 m; 367–921 m on the others);
  - the plateaus ≥ 2 % on 1 seed (9: 5.6 %);
  - refuted.
- **R7**: held (3 / 4).
- **Meta**: held.

**Mine**:
- **Held**:
  - R1 on the main except seed 9, and R1 refuted on « own s »;
  - R6's plateaus refuted;
  - R7;
  - the macro success not met (but for the mountains, not the land);
  - the regularity on σ_θ and R2_g (at 391 m);
  - Hovius within the band on ≤ 2 seeds (1 at A ≥ 25 km²);
  - the cost;
  - meta.
- **Refuted**:
  - **the land dropping sharply / R4 refuted** (1–2 % lost);
  - R2 holds;
  - R3 on ≥ 2 seeds;
  - R5's 20–40 % drop (−1 to −22 %);
  - the cratons above 1 km.

## Cost (CPU)

- The Airy record: 0.2–0.3 s per C1 run.
- The physics level: 9.2–15.4 s CPU, as F166's 11–13 s.
- A frozen chain: ~170–235 s CPU (~30 s wall).
- The bench took 825 s wall.

## Queued (from the brief)

- F168, the sea-level rise (now dependent on this result: no land sits near zero except the margins);
- F169, the heterogeneous K (the regularity baseline above is its reference);
- C1 (b): the wedge's extent, the collisions, the mapping in production;
- the facets and walls (the talus threshold?);
- the map-border rim;
- the coasts;
- the warp's depressions;
- the river width;
- the deposition; the alignment; the lakes; K per rock; the resources.

## At commit (2026-10-10)
- **The Airy mapping with the sea at 0 m absolute becomes the cascade's default** (the oceanic plates read at
  S̃ = 0.2). Production is unchanged.
- **The reviewer's errors, recorded** (above):
  - Whitehead & Clift's numbers are not one Airy line;
  - Lachenbruch & Morgan's ridge is a lithosphere column;
  - the paper is Whitehead **and Clift** 2009.
- **Seed 9's plateaus are craton blocks** (the reviewer's reading, a hypothesis to check): the 5.6 % of plateau sits on
  the 1.78 km craton shields, flat 64² blocks under Airy.
- **The queue**:
  - F168 = the land's size (C1 (b), part 1); F169 = the sea-level rise; F170 = the heterogeneous K;
  - then the rest of C1 (b): the wedge's reach if it is still at fault, the collisions, the switch in production;
  - the facets and walls, the map-border rim, the coasts, the warp's depressions, the river width, the deposition,
    the alignment, the lakes, K per rock, the resources.
- **The reviewer's hypothesis for F168** (to check): our continents are islands, 18–29 % of the map against ~58 % in
  the European window. With Corsica-like slopes a single 3 km range is ~85 km wide in mountain class, so on 30 000 km²
  of land one range is a quarter of it.

**The author's decisions for F168, verbatim**:
- « Tant qu'on peut assurer de la mer tout autour, pourquoi pas. »
- « On applique déjà un offset dans le viz pour centrer le continent et qu'on ait un chemin qui permet de faire le tour
  du continent par la mer. Ainsi, globalement on obtient toujours une ile à la fin. Par contre le % de terre vs mer peut
  être variable. Ce n'est pas pareil que de dire qu'on a aucune terre qui touche le bord de la carte. Toutes les seeds
  témoins sont dans le cas: devient une ile lorsqu'un offset est appliqué car on peut en faire le tour. »
