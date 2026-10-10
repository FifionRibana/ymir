# Finding 168 — C1 (b), part 1: a continent that fills the map and can be sailed around

**Status: built, measured, gated. Production unchanged. Nothing committed without the feu vert.**

**Verdict**:
- The land target is reachable and monotone, but quantised: L2 = 62.5 % on the témoin; no fraction gives 55–60 %.
- **Two stops fire on 4 / 4 seeds**:
  - **« the mountains above 37.2 % »**: 75.0–79.5 % in L2. The mountain share does **not** fall as the land grows.
  - **« the circumnavigation lost »**: in L2 at 64² and 2 048², and in L1 on 3 / 4.
- An addition made after the measure locates the loss: **it is set at t₀ by the continental-plate mask, on 12 / 12
  worlds**. C1 neither makes nor breaks it.

**Files**:
- Declared before any measurement: `f168_declared.md`. Predictions: `f168_predictions.md` (written first).
- Bench: `crates/ymir-core/tests/f168_cascade.rs`, with three tests:
  - `f168_trials` → `bench_trials.txt`;
  - `f168_cascade` → `bench_output.txt`;
  - `f168_t0` → `bench_t0.txt`.
- Checks: `checks_before.txt`, `checks_after.txt`.

## Part 0

F167 committed (`147f5f3`) with its « At commit » block and the author's decisions for F168, verbatim.

## A — the inventory (declared before measuring; `f168_declared.md`)

**A1 — `continental_fraction` is a fraction of the plates** (`clustering.rs:58–78, 150–200`):
- round(8 f) continental plates out of 8, picked by BFS from one seed plate, so the land comes in whole Voronoi plates;
- the continents never move (rigid, F165-A1), so their outline is fixed at init;
- `craton_shield_fraction` (0.15) and `num_plates` (8) were not touched.

**A2 — Davis-Suppe** (`source_term.rs:53–84, 133–141`):
- l_taper 4, l_decay 6, max_distance 30 cells of the 64² reference grid;
- « Defaults selected for Phase 1.2 visual demonstration at 64² », i.e. a visual calibration;
- read only, nothing changed.

**A3 — the offset, the border and the test**:
- **The viz's offset** (`hd.rs:451–455`) centres the largest mass's circular extent. It only centres and guarantees
  nothing.
- **The benches up to F167 used the témoin's roll (6, 37) for every seed.** F168 applies the viz's rule to each world's
  own final 64² field.
- **The cascade is periodic, except the solver's off-map receivers** (`stream_power.rs:651`): land on the border is a
  fixed base, which is the F166 rim.
- **The circumnavigation test** (declared before measuring; `physio::circumnavigation`):
  - the largest periodic 8-connected mass must not wrap the torus (a BFS carrying unwrapped coordinates);
  - r* = the largest Chebyshev dilation of all the land after which its component still does not wrap, found by
    bisection. A sea loop at least 2r* + 1 cells wide then surrounds the mass.
  - Unit tests: a square (r* = 21, exact), a band, a diagonal, a square across both seams (r* = 26), the dilation ball.

## B — the variants

### B1, the trials on the témoin (`bench_trials.txt`; land = the 256² history physics level above 0 m, whole torus)

| f | continental plates | continental cells at t₀ | land 64² t₀ → end | **land 256²** | peak |
|---|---|---|---|---|---|
| 0.29 (L0) | 2 | 18.7 % | 14.8 → 19.1 % | **19.4 %** | 2 907 m |
| 0.375 (**L1**) | 3 | 35.9 % | 29.8 → 35.3 % | **35.7 %** | 2 190 m |
| 0.5 | 4 | 52.3 % | 45.2 → 49.7 % | **51.1 %** | 3 906 m |
| 0.625 (**L2**) | 5 | 64.5 % | 56.8 → 61.5 % | **62.5 %** | 6 294 m |
| 0.75 | 6 | 75.7 % | 69.4 → 72.0 % | **74.2 %** | 5 945 m |

- **Monotone.** Each plate adds 11–16 % of the map.
- **55–60 % is not reachable with 8 plates.** The declared rule (« inside 55–60 %, else the nearest ») takes f = 0.625
  (62.5 %), and 62.5 % lies inside 50–65 %, so the « unreachable » stop does not fire.
- **L1 = 0.375 (35.7 %) and L2 = 0.625 are frozen** in `cascade::profile`, as `CONTINENT_L1` / `CONTINENT_L2`.
- **On the other seeds** the same fractions give L1 42–45 % and L2 68–70 %.

### B3 — the negative control (permanent)

- `profile::tests::the_continent_profile_at_production_fraction_is_c1_bit_for_bit`: f = 0.29 gives C1's production run
  bit for bit (s, plate ids, plate types).
- In the bench, L0 under the canonical roll reproduces F167's témoin exactly: 2 670 m, 31 091 km², 74.4 % mountain.

### B4 — the viz

- The Cascade window adds a row « Profil C1 : production / continent L1 / continent L2 ».
- It applies to the history modes.
- A continent profile is framed by the viz's rule on its own final field. The production profile keeps the workspace's
  offset.

## C — the measures (`bench_output.txt`)

Europe for reference: plain 51.6 / plateau 4.7 / hill 18.9 / mountain 24.8 %.

### The land and the circumnavigation (64² end of C1 | 2 048²)

| seed | L0 land · circ. (r* → sea loop) | L1 land · circ. | L2 land · circ. |
|---|---|---|---|
| témoin | 19.4 % · yes, 156 km | 35.7 % · yes, 44 / 50 km | 62.5 % · **wraps** |
| 42 | 26.4 % · yes, 19 / 28 km | 45.2 % · **wraps** | 69.6 % · **wraps** |
| 1 | 28.5 % · yes, 19 km | 42.3 % · **wraps** | 68.4 % · **wraps** |
| 9 | 31.9 % · yes, 81 / 85 km | 43.5 % · **wraps** | 69.9 % · **wraps** |

- **One mass > 1 000 km² everywhere.** Where the mass wraps, it wraps at both sizes.
- **The addition after the measure** (`f168_t0`; it informs the proposal only):
  - the circumnavigation status is the same for the continental-plate mask at t₀, the land at t₀ and the land at the
    end, on 12 / 12 worlds;
  - where circumnavigable, the final r* equals the plates' r*.
  - **So the loss is decided at init, by which plates the BFS makes continental.**
- **Seed 1 L0 is circumnavigable (r* = 1) yet `land_topology` flags a band in x**: its mass occupies every column
  without closing a loop.
  - No rectangular window frames it free of the border: 229 border cells after the offset.
  - The author's « devient une île lorsqu'un offset est appliqué » holds on 4 / 4 in the sense of circumnavigation. It
    does not hold on 1 / 4 in the sense of « no land on the window border ».
  - Caveat: the bench's offset is computed on the cascade's Airy field, not on production's preview.

### The border rim

- Present wherever land touches the window after the offset: the border land is 1.15–2.1× the land's mean altitude
  (e.g. témoin L2: 1 920 m against 905 m).
- That is every L2, three L1, and seed 1 L0.
- Absent on the other L0 seeds and on the témoin's L1.

### The interior (share of the land > 20 km / > 40 km from the coast)

- **L0**: 46–63 % / 9–34 %.
- **L2**: 58–77 % / 24–53 %; the largest distance 87–124 km.

### Macro and the geometric hypothesis

| seed | mountain L0 → L1 → L2 | plain + plateau in L2 | slope L0 → L2 |
|---|---|---|---|
| témoin | 74.5 → 59.7 → **75.0 %** | 0.4 % | +0.08 pt / 10 000 km² |
| 42 | 52.3 → 74.7 → **79.5 %** | 1.6 % | +3.93 |
| 1 | 76.0 → 63.8 → **77.9 %** | 0.7 % | +0.30 |
| 9 | 75.7 → 73.9 → **78.8 %** | 9.6 % | +0.50 |

**The hypothesis recorded at F167's commit is refuted on 4 / 4**: « one range is a quarter of a small island, so a
larger continent has a smaller mountain share ». The mountain share does not fall with the land area.

### C1's tectonics (end of C1, 64²)

- **The plates at the end**: 2 on 11 / 12 worlds (1 on the témoin's L1).
- **C-C collision cells**: 0 in every L0. In L2 on 4 / 4: 32, 86, 58, 50. In L1 only on seed 42 (44).
- **The active-margin (upper-plate) cells**: 54–80.
- **Davis-Suppe Σ Δs > 0.05**: on 57–87 % of the land in L0, and 50–68 % in L2.
- **The cratons**: 57–220 cells; land p50 at 256² 273–1 129 m; block index 3.3–4.0, against 0.6–1.0 for the coasts.
  The 64² blocks show.
- **The Davis-Suppe profile** (`images/f168_ds_profile.png`):
  - in L0 it peaks at 1–2 cells from a convergent cell or suture (1.6–1.9) and falls under 0.05 by 6–9 cells on three
    seeds;
  - in L2 it peaks lower (0.59–0.97) but reaches further: still 0.2–0.5 at 8–10 cells on seeds 42 and 1.

### The L2 mountain class decomposed (256²)

| component | share |
|---|---|
| cratons | 6.5–8.7 % |
| Davis-Suppe zones ≤ 5 cells (31 km) from a convergent cell or suture | 37.9–58.7 % |
| Davis-Suppe zones 6–15 cells | 5.6–14.6 % |
| Davis-Suppe zones > 15 cells | 0 |
| arcs | 0–0.3 % |
| **the rest** | **28.6–39.3 %** |

### 2 048² against Corsica

- **The slope guard holds**: the mountain class's slope p50 is 0.27–0.37 against Corsica's 0.318.
- **The peaks**: L0 3.1–4.1 km; L2 4.3–6.6 km.
- **The facets fail in the mountain class in L2 on 3 / 4** (all but seed 9). The walls fail on 2 / 4.
- **A_dir fails everywhere** (alignment, queued).
- **The lakes at 256² in L2**: 83–6 541 km². The largest is seed 9's, on its plateau.
- **F by class (L2)**: plain 0–29 %, plateau 6–27 %, hill 7–10 %, mountain 1.2–2.4 %.
- **The coast: the « degraded > 10 % » clause cannot be read in L2.**
  - `physio::coast` treats outside the window as sea (checked in the code), so land on the window border adds straight
    edges on the 32-cell lattice.
  - The block index jumps to 7.4–11.0 and the aligned share to 0.45–0.55, wherever land touches the border. Seed 1 L0
    shows it too: 5.5.
  - **Where no land touches the border, the coast is F167's**: the témoin's L1 has block index 0.95, aligned 0.26,
    dimension 1.00; L0 has 0.61–0.96.
- **The regularity baseline in L2** (Corsica's values in parentheses):

  | measure | L2 | Corsica |
  |---|---|---|
  | σ_θ | 23.2–24.3° | 27.4–29.9° |
  | R2_g at 391 m | 0.62–0.71 | 0.48 |
  | CV_L | 0.88–1.11 | 0.87–0.92 |
  | Hovius R (A ≥ 25 km²) | 2.34–2.72 | 2.33–2.55 |

  - More regular than Corsica on the confluence angles and the gullies' parallelism at 391 m, as in F167.
  - No longer more regular on the tributaries' lengths.
- **The images** (north up, the viz's offset):
  - `f168_seeds_256.png`, with the r* sea loop in red;
  - `f168_seeds_2048.png`;
  - `f168_classes_seeds.png`;
  - `f168_coast_25km.png`, each world's own half-land window, against Corsica;
  - `f168_ds_profile.png`.
  - The Copernicus notice is in `images/NOTICE.md`.

### Cost (CPU seconds per world, C1 included)

- C1 0.2–0.3 s; the physics level 10–26 s; the frozen chain 168–427 s (two runs).
- **Per world: L0 180–252 s, L1 262–320 s, L2 370–448 s** (two identical runs; the CPU time varies by ≤ 10 %). The chain grows with the land.

## R — the stop rules

- **Macro success in L2: 0 / 4.**
- **The mountains above 37.2 % in L2 on 4 / 4**: decomposed above, stopped.
- **The circumnavigation lost on 4 / 4 in L2**: at both sizes, and already at t₀. Nothing is implemented. The
  proposals follow.
- **The land target**: reachable (62.5 % lies inside 50–65 %) and monotone. 55–60 % is not hit, because of the
  whole-plate quantum.

## The proposals for the circumnavigation (implemented: none)

1. **The most neutral guarantee: an initial condition on the continental plates' selection.**
   - Among the BFS clusters of round(8 f) plates, take the one whose union is circumnavigable with the largest r*. The
     selection stays deterministic from the seed, with ties broken by the current RNG order. If no cluster qualifies,
     report that the seed cannot host this fraction.
   - The physics is untouched, and the measured invariance (12 / 12, C1 preserves the status) carries it to the end.
   - It moves the continents' outline, not their behaviour.
2. **A geometric limit for the author.**
   - On the 400 km torus, a continent covering 57.5 % of the map leaves at most a 97 km strait if it were a perfect
     square, and 58 km for a disc. At 62.5 %: 84 km and 43 km.
   - So 55–60 % of land with sea all around requires a compact continent. The Voronoi clusters are elongated.
3. **Less neutral**: more plates (a finer quantum, also reaching 55–60 %). This changes the plate geometry and the
   tectonics.
4. **The window border**: being circumnavigable does not give « no land on the window border » (seed 1 L0). If the
   author wants the latter, it is a separate criterion: r* large enough for a straight frame.

## The reading (hypotheses to check)

- **The mountain share is not set by the continent's size.** Half the L2 mountain class is the Davis-Suppe wedge
  within 31 km of a convergent boundary or suture, and the C-C collisions add boundaries. A third is outside every
  tectonic zone (« the rest »).
  - **Hypothesis to check**: the physics level classifies dissected land above ~300 m as mountain whatever its source
    (Kapos: 300–1 000 m with R7 ≥ 300 m). That would point back at the wedge's reach and the uplift's extent, F169's
    alternative.
- **The témoin's L0 peak is 2 907 m under the viz's offset, against 2 670 m under the canonical roll**, although no land
  touches the border in either.
  - **Hypothesis to check**: the cascade's noise fields (roughness, warp) are tied to the window, not to the land, so the
    cascade is not translation-invariant.

## Predictions

| | held | refuted |
|---|---|---|
| **The reviewer's** | S1 (0.625, though 55–60 % is unreachable); S2's second clause (> 37.2 % on 4 / 4); S4 (4 / 4); S5 (50–68 %, 4 / 4); S6 (4 / 4); S7 (the facets on 3 / 4); the meta | **S2** (the mountains rise L0 → L2 on 4 / 4); **S3** (plain + plateau 0.4–9.6 %) |
| **Mine** | L1 = 0.375 (35.7 %); the monotone response; S1; « just outside 55–60 % »; S6; S4 (2 plates at the end); S5 within 40–75 %; S3 refuted; S7; no macro success, with the mountains' stop and the circumnavigation's; the rim where land touches the border; C1's cost; the meta | L2's 50–62 % (62.5 %, marginally); **S2 on both clauses** (I expected a fall of 10–25 points); the interior (L0 > 40 km is 9–34 %, not 0–5 %; L2 24–53 %, not 15–35 %); the physics level (10–26 s, not 10–15) and the chain's cost (L2 352–427 s, not ~200) |

## Open, for the author (nothing chosen)

- **(a) The circumnavigation**: the initial selection condition (proposal 1), knowing the geometric limit (proposal 2).
- **(b) The mountains**: with S2 refuted, F169 = the wedge's reach (the critical-wedge mass balance) before the
  sea-level rise. Or first check « the rest » (the hypothesis above).
- **(c) The coast measure**: whether it should exclude the window border.

## At commit (2026-10-10)
- **The reviewer's error, recorded**: F167's commit hypothesis (« one range is a quarter of a small island, so a
  larger continent has a smaller mountain share ») is refuted on 4 / 4. The mountain share does not fall with the land
  area.
  - The reviewer's new hypotheses, to check in F169:
    - the mountain share follows the internal sutures' length per unit area;
    - **everything is young**: the physics level spreads all of C1's history linearly over 5.62 Myr, so an early
      suture gets 5.6 Myr of erosion, whereas a range without uplift wears down in 1–10 Myr (×6 at most with isostasy,
      Baldwin, Whipple & Tucker 2003);
    - « the rest » (29–39 %) is the initial crust at +0.45 km, dissected in 5.6 Myr.
- **The plate quanta**: `continental_fraction` counts plates, so the land comes in steps of one plate. The land target
  is relaxed to 50–62 % according to the step.
- **The invariance bug, recorded as a defect** (no longer only a hypothesis):
  - the cascade is not translation-invariant: the témoin's L0 peak is 2 907 m under the viz's offset and 2 670 m under
    the benches' roll;
  - its noise fields are tied to the window;
  - the rim comes from the solver's off-map receivers.
- **The coast measure rule**: `physio::coast` counts the window border as coast. From F169 the coast is measured on the
  torus.
- **L1 / L2 stay profile constants** (`cascade::profile::CONTINENT_L1` = 0.375, `CONTINENT_L2` = 0.625).

**The author's answers, verbatim**:
1. On the plates' selection by an initial condition, with the land target relaxed to 50–62 % by step: « Oui. SI ça
   marche pas on fera plus de plaques. »
2. On the continent's age A as a world parameter: « OUi. Ca donnera un slider potentiellement. Et on part sur du vieux.
   A évaluer la différenec entre jeune, intermédiaire et vieux. »
3. On the order (F169 = the infrastructure, the chronology and the age; F170 = the sea-level rise; F171 = the
   heterogeneous K): « Oui »
