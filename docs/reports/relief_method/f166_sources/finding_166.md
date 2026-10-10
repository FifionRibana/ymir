# Finding 166 — the physics level driven by C1's tectonic sources only

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ». **Production is
unchanged**: the guard 6 / 6, field and lakes, was run before and after (`checks_before.txt`, `checks_after.txt`).
C1's outputs are bit-identical with the term observer (a permanent test).

- Declared before any measurement: `f166_declared.md`. Predictions: `f166_predictions.md`.
- Bench: `f166_bench_output.txt` (848 s wall), `f166_terms_output.txt`.

**The verdict in four sentences.**
- **Q1 holds: the subsidence is C1's own erosion.**
  - On the témoin's interior, C1's erosion closure makes 100 % of the negative Σ Δh (−1 205 m); Davis-Suppe adds
    +163 m.
  - The isostatic datum, a percentile of the whole field, lifts the land by +657 m.
  - The equilibrium-height sink never touches the land.
- **Driven by the tectonic sources only, the continent stops drowning.**
  - It loses 0–4.5 % of F164's land (F165's history lost 21–38 %).
  - Its interior at ≤ 0 m is < 1 % on three seeds.
  - The peaks vary with the tectonics: 2.7–5.4 km at 256², 2.6–5.1 km at 2 048², with T frozen.
- **But the macro fails on every seed**:
  - mountains 56–64 % on three seeds (28 % on seed 42), against Europe's 24.8 %;
  - plain + plateau 1–11 %, against Europe's 56.3 %;
  - **no plateau** (0.0 %).
- **The reading**:
  - C1 starts every continent as a uniform plateau at ~1.8 km (normal 35 km crust through C1's isostatic ramp);
  - without C1's erosion it stays high, and 5.6 Myr of the physics level's erosion dissects it into Kapos mountains;
  - no stop rule fires: Q1 holds, so not « attribution », and there is no macro success, so not « success without
    plateaus ». Reported.

## Partie 0

1. **F165 committed** (`ee841a8`), the author's seed and config edits left out.
   - Its ADR records:
     - the reviewer's decisions: the sources' drive allowed, collisions after F166, **the width test withdrawn** (a
       design error of the reviewer's), the warp's depressions measured only;
     - **the reviewer's errors**: the width test, and P6 (a late uplift is MORE drained);
     - the author's remarks and answers, verbatim.
2. **Checks**:
   - before: the guard 6 / 6 (field and lakes = banc on all six states), lib 639, viz 33;
   - after: the guard 6 / 6 identical (field and lakes = banc on all six states), lib 639 → 641, viz 33.
   - **No ALGO_* bump**: C1's loop gained a read-only observer, and its outputs are bit-identical (the test and the
     guard), so no cached product changes.
   - Two permanent tests are new:
     - `the_term_observer_leaves_c1_bit_identical`: C1 is bit-identical, and the terms' Δs add up to the run's change of
       s within 1e-9;
     - `the_sources_with_every_term_are_f165s_history`: with nothing excluded, the record is F165's bit for bit; with
       the erosion excluded, C1 itself is untouched.

## A — the terms of C1

### A1 — the inventory (in `f166_declared.md`, with file:line and the formulas)

| term | class |
|---|---|
| transport (upwind, continental v = 0) | tectonic |
| Davis-Suppe (orogenic thickening to the critical taper) | tectonic |
| **equilibrium height** (s > 2.0 ⇒ s ← max(2.0, s − 2·(s − 2)²·dt)) | **ambiguous**: a collapse of thick crust, but it removes mass and ends in a clamp; measured both ways |
| **erosion** (stream power, **down to a floor of s = 0.2, the oceanic thickness**) | **surface** |
| subduction (consumption, arc mass, type change) | tectonic |
| rift thinning | tectonic |
| accretion, rift split | tectonic (plate_id, age only) |
| **the isostatic datum** (the sea level a percentile of the whole field) | **other: a mapping, not a term on s** |

### A2 — F165's history decomposed (the témoin)

On the témoin's 640 land cells, the interior is only 20 cells: the land lies within 3 cells of the sea, a convergent
cell or a suture almost everywhere at 64².

| term | interior: mean Σ Δs (Δh) | margin: mean Σ Δs (Δh) |
|---|---|---|
| transport | +0.000 (0 m) | −0.006 (−19 m) |
| Davis-Suppe | +0.052 (+163 m) | **+0.905 (+2 995 m)** |
| equilibrium height | **0** | **0** |
| **erosion** | **−0.381 (−1 205 m)** | **−1.243 (−4 106 m)** |
| subduction | 0 | +0.004 (+12 m) |
| rift thinning | 0 | −0.030 (−99 m) |
| **datum drift** | **+657 m** | +710 m |
| observed Δh_iso | −564 m | −569 m |

- **Q1 holds**: 100 % of the interior's negative Σ Δh is the erosion closure.
- **The datum fell** (raw h_sea 0.087 → 0.050), lifting every land cell by ~700 m at constant s.
- On the whole land (`f166_terms_output.txt`):
  - Davis-Suppe +2 906 m on average and the erosion −4 015 m: **C1 piles crust at the margin and erodes it away,
    tens of km of crust at the extreme cells**;
  - the rift −95 m, the transport −18 m, the subduction +12 m.
- **The equilibrium-height sink is 0 on every land cell**: it acts only on oceanic pile-ups at the convergent
  boundaries.
  - Yet excluding it changes the land (the variant below), **through the datum**: removing the oceanic pile-ups moves
    the whole field's p92, hence the sea level.
  - This is a hidden coupling of C1's isostatic mapping: an oceanic process changes the continent's altitude
    everywhere.

### A3 — the initial state

- `init_c1_state_phase_2_r7`: continental s = 1.0 (35 km), oceanic 0.2, the cratons × 1.25, 15 % of them shields.
- **h_iso(t₀) on the land**:

  | seed | land mean | deciles 2–9 |
  |---|---|---|
  | témoin | 1 573 m | 1 522–1 832 m |
  | 42 | 1 561 m | 1 519–1 832 m |
  | 1 | 1 568 m | 1 541–1 832 m |
  | 9 | 1 626 m | 1 541–1 832 m |

  - Half of the land sits at exactly 1 832 m: **one uniform plateau, 35 km of crust rendered 1.8 km high by C1's
    land ramp** (sea datum → 5 650 m at S̃ = 2.0).

## B — the mechanism

- **The observer**: `time_loop::run_with_closures_terms`, a read-only hook after each term.
- **The record**: `history::run_c1_sources`, s_src = s − Σ(the excluded terms' Δs), h_src = C1's isostasy of s_src.
- **T, calibrated once on the témoin (sources): T = 5.62 Myr** (7 trials, 75.6 s CPU / 43.3 s wall):
  - 10⁵ yr → 5 762 m; 10⁶ → 5 529 m; 10⁷ → 839 m; 3.16 × 10⁶ → 4 029 m; then 5.62 × 10⁶ → 2 728 m.
  - **T is ×0.19–0.56 of C1's nominal 10–30 Ma**, against ×0.04–0.12 for F165's 1.16 Myr.
- **Frozen for the four seeds and the variant.** The témoin's F165 history recomputed matches F165 (peak 2 852 m,
  land 20 381 against 20 376 km²).

## C — the measures

### 256² (Europe: plain 51.6 · plateau 4.7 · hill 18.9 · mountain 24.8)

| seed | run | plain | plateau | hill | mountain | peak | land lost / gained vs F164 | interior ≤ 0 | lakes | mountain regions > 100 km² |
|---|---|---|---|---|---|---|---|---|---|---|
| témoin | F164 | 8.7 | 0.0 | 48.5 | 42.8 | 2 863 | — | 0 | 14 / 2 458 km² | 2 |
| | F165 | 1.1 | 0.1 | 39.9 | 58.8 | 2 852 | 23.2 / 2.9 % | 9.67 % | 13 / 354 | 3 |
| | **sources** | 1.0 | 0.0 | 37.6 | **61.3** | 2 728 | **0.0 / 16.4 %** | 0.03 % | 2 / 85 | 1 |
| | sources − EH | 2.6 | 0.0 | 44.1 | 53.3 | 2 349 | 0.0 / 16.1 % | 0.87 % | 1 / 88 | 1 |
| 42 | F164 | 30.6 | 0.0 | 57.7 | 11.6 | 1 815 | — | 0 | 11 / 3 835 | 6 |
| | F165 | 4.0 | 0.0 | 41.0 | 55.0 | 2 454 | 38.1 / 1.0 | 4.41 % | 11 / 591 | 6 |
| | **sources** | 10.5 | 0.0 | 61.1 | **28.3** | 3 261 | **4.5 / 20.2** | 5.64 % | 4 / 200 | 6 |
| | sources − EH | 12.6 | 0.0 | 63.0 | 24.4 | 3 281 | 14.3 / 16.4 | 5.71 % | 3 / 222 | 7 |
| 1 | F164 | 36.3 | 0.0 | 50.5 | 13.2 | 2 498 | — | 0 | 12 / 4 124 | 4 |
| | F165 | 2.6 | 0.0 | 36.9 | 60.6 | 3 475 | 22.6 / 0.6 | 10.60 % | 9 / 239 | 4 |
| | **sources** | 1.5 | 0.0 | 34.6 | **63.9** | 4 301 | **0.0 / 17.2** | 0.12 % | 5 / 879 | 2 |
| | sources − EH | 1.8 | 0.1 | 36.2 | 61.9 | 4 123 | 0.0 / 17.0 | 1.38 % | 5 / 745 | 4 |
| 9 | F164 | 38.9 | 0.0 | 44.0 | 17.1 | 2 545 | — | 0 | 12 / 1 755 | 5 |
| | F165 | 1.7 | 0.1 | 32.5 | 65.7 | 3 964 | 21.3 / 3.5 | 6.23 % | 39 / 2 361 | 6 |
| | **sources** | 4.8 | 0.0 | 39.2 | **56.0** | 5 365 | **3.8 / 12.9** | 0.32 % | 17 / 4 023 | 8 |
| | sources − EH | 9.3 | **6.1** | 40.6 | 44.0 | 5 403 | 21.0 / 10.7 | 4.04 % | 16 / 1 201 (14.8 % on plateaus) | 6 |

- **The sources keep the land**: 0–4.5 % lost, and they even gain 13–20 % over F164 (the steady level drowned none,
  but F164's land was its own).
- **Their mountains are as many as F165's on three seeds.** The initial 1.8 km plateau stays high and is dissected.
- **The variant « sources − EH » differs a lot, through the datum alone** (the EH sink touches no land):
  - seed 9 loses 21 % of the land and gains **6.1 % of plateau**, the only plateau of the round (14.8 % of its lake
    area on it);
  - the témoin's peak falls by 380 m.
- **The lakes do not follow one direction**: against F165, −76 % / −66 % / +268 % / +70 %.
- **The blocks**: the {Σ U > 0} boundary's 64² block index is 1.6 / 1.3 / 1.8 / 1.5 (F165's 1.7 / 1.2 / 1.6 / 1.4):
  mildly blocky, not the rectangle feared.
- **A plateau-edge profile**: there is no plateau of any extent in the sources' runs.
  - The profile drawn (`f166_plateau_profile.png`, seed 42) crosses two isolated plateau-class cells on a dissected
    slope, with no break of slope. Reported as « no plateau edge exists ».

### 2 048² (the frozen cascade), against Corsica at 195 m

| seed | sources, whole | sources, mountain class | the slope guard: peak, mountain-class slope p50 / p90 (Corsica 0.318 / 0.550) |
|---|---|---|---|
| témoin | 17 / 19 (the 6–12.5 km octave, A_dir) | 15 / 19 (**facets ×1.9, walls ×1.6**, octave, A_dir) | 2 914 m · 0.344 / 0.507 |
| 42 | 14 / 19 | 14 / 19 (**facets**) | 2 601 m · 0.326 / 0.497 |
| 1 | 15 / 19 (**facets, walls**) | 15 / 19 (**facets ×2.5, walls ×2.1**) | 4 925 m · 0.382 / 0.522 |
| 9 | 15 / 19 | 12 / 19 (**facets, walls**) | 5 098 m · 0.342 / 0.515 |

- **The slopes are Corsica's** (within ×1.2) even at 5 km peaks: the guard is quiet.
- **The facets and the 28° walls still fail in the mountain class on 4 / 4 seeds**, as in F165: **Q5 holds**. Their
  excess is lower than F165's (×1.6–2.5 against ×2.7–3.4).
- **F by class** (sources): plain 25–39 %, plateau 0–10 %, hill 7–12 %, mountain 0.7–1.8 %.
- **The coast**:
  - no measure is degraded on three seeds;
  - **seed 9 is degraded** on the alignment (59.9 % against F164's 47.2 %) and the block index (12.6 against 8.2).
    The kept coast is the 64² mask's, no longer softened by the drowning.
- **The warp's depressions surviving to 2 048²**:
  - of the new depressions at 512², 1 024² and 2 048²: 13–42 %, 17–40 % and 18–40 % still hold a lake in the final
    field (sources);
  - F164: 27–48 %, 13–23 %, 13–21 %.
- **Seed 9's rim**: 630 land cells on the map border, mean 1 476 m against the land's 638 m.
  - **The cause**: `erosion/stream_power.rs:651` makes a cell whose D8 receiver leaves the map its own receiver (a
    fixed base node, never incised), and `terrain/flow.rs:256` routes D8 across the periodic wrap. A border cell
    draining outward is such a node, and under a positive U it only rises (the history does not hold the sea).
  - **Diagnosis only**: no fix is trivial without touching the solver's boundary rule, which production uses.

## R

- **Macro success: 0 / 4 seeds.**
  - The mountains are within 16.5–37.2 % only on seed 42.
  - Plain + plateau reach 37.5 % on no seed (1.0–10.6 %).
  - The land lost is under 15 % on all four; no coast is degraded on three.
- **The « attribution » stop: no** (Q1 holds).
- **« Success without plateaus »: no**, there is no success. **No rule fires; reported.**

## The reading (hypotheses to check)

- **The plateau and plain deficit sits in C1's vertical mapping, not in the drive.**
  - C1 renders normal 35 km continental crust at ~1.8 km: its land ramp maps the datum → 5 650 m at S̃ = 2.0, so
    S̃ = 1 lands at about a third.
  - Earth's continental freeboard puts such crust at a few hundred metres.
  - With C1's erosion the land sank (F165); without it, it stays a 1.5–1.8 km block that the physics level dissects
    into mountain.
  - Either way, no low smooth ground is ever made.
  - This is C1's work (b), with the collisions and the wide thickening, and it echoes the earlier floor diagnosis
    (« le socle haut = cratons trop hauts »).
- **C1's isostatic datum couples the ocean to the land.** The sea level is a percentile of the whole field, so the
  equilibrium-height sink, which acts only on oceanic pile-ups, moves the continent's altitude by hundreds of metres.
- **The texture defect (facets, walls) is independent of the drive** (Q5), as the reviewer predicted.

## Images (north up; Copernicus notice in `images/NOTICE.md`)

- `f166_seeds_256.png`, `f166_seeds_2048.png`, `f166_classes_seeds.png`: rows = the four seeds, columns F164 | F165 |
  sources.
- `f166_terms.png`: the témoin's Σ Δh per C1 term on the land: transport, Davis-Suppe, equilibrium height, erosion,
  subduction, rift thinning, then the datum drift (one saturation at the pooled p98, ±8 831 m; red +, blue −).
- `f166_coast_25km.png`: Corsica | F164 | sources.
- `f166_plateau_profile.png`: no plateau edge exists; see above.

## Predictions

**The reviewer's**:
- **Q1**: held (100 %).
- **Q2**: held (land lost 0–4.5 % on 4 / 4; interior ≤ 0 under 1 % on 3 / 4).
- **Q3**: refuted (plateaus 0.0 %; mountains under 40 % on 1 / 4).
- **Q4**: held (5.62 Myr).
- **Q5**: held (4 / 4).
- **Q6**: refuted (+30 % on 2 / 4; −76 % and −66 % on the others).
- **Q7**: held (2 728–5 365 m, a 2.6 km spread).
- **Meta**: held.

**Mine**:
- **Held**:
  - Q1 (the erosion ≥ 80 %), the EH sink ≈ 0 on the interior, the tectonic terms ≥ 0 on the interior;
  - the seeds' initial state alike within ±15 %;
  - Q5, Q7;
  - the macro success ≤ 1 seed, no stop rule;
  - seed 9's rim cause;
  - the observer < 1 s (0.2 s);
  - meta.
- **Refuted**:
  - h_iso(t₀) 800–1 300 m (1 561–1 626 m mean, half at 1 832 m);
  - **Q4 refuted** (T = 5.62 Myr, not 1.5–4);
  - the land lost 10–25 % (0–4.5 %);
  - the mountains < 40 % on ≥ 3 seeds (1);
  - Q6 holds (2 / 4);
  - **the « − EH » variant within 5 %** (it moves the land by up to 21 % through the datum);
  - the block index ≥ 3 (1.3–1.8);
  - the warp's survivors ≤ 30 % (up to 42 %);
  - the physics level 5–12 s CPU (10.7–14.1 s).

## Cost (CPU seconds; the wall times beside, no suspension flagged)

- C1 with the observer: 0.2 s.
- The sources' physics level: **10.7–14.1 s CPU** (7–9 s wall), against F164's steady 3.6–13.4 s, so **≈ +0 to
  +10 s CPU per world** (the reference is 249.8 s).
- The T calibration: 75.6 s CPU (7 trials).
- A frozen chain to 2 048²: 160–235 s CPU (28–37 s wall).
- The bench took 848 s wall.

## Queued (from the brief, not done)

- F167, the sea-level rise;
- C1's (b) (the collisions, the wide thickening, the 64² sources, and now the vertical mapping and the datum's
  coupling);
- the coasts;
- the warp's depressions;
- the map-border rim;
- the river width (`bed_width_a`);
- the deposition; the alignment; the lakes; K per rock; the resources.

## At commit (2026-10-10)
- **A new order**: F167 = the thickness → altitude mapping; F168 = the sea-level rise; F169 = a heterogeneous K. The
  reason: the sea-level rise depends on the land near zero, which the mapping will change.
- **The mapping is corrected in the cascade only.** The production correction comes with C1's work (b), the
  collisions.
- **Queued**:
  - the facets and the walls (the reviewer's hypothesis, to check: the talus threshold);
  - the map-border rim.

**The author on F166, verbatim**: « Les terrains sont maintenant complètement modifiés par rapport à ce qu'on connait
mais semble beaucoup plus intéressant. Des plateaux apparaissent, des chaines, des vallées, des plaines. C'est de plus
en plus propre. Mon avis est qu'il faudra bien traiter la trop grosse régularité des motifs. un K bruité autour des
valeurs dépendant de la tectonique serait peut-être une solution? »
