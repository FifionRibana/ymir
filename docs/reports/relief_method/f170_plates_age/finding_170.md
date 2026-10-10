# Finding 170 — more plates for a continent that can be sailed around, then the chronology and the continent's age

**Status: built, measured, gated. Production unchanged. Nothing committed without the feu vert.**

**Verdict**:
- **The selection works with 16 plates**: 3 of 4 seeds qualify (42, 1, 9; r* 1–3, 51–52 % land). 12 plates qualify none; the
  témoin's best is r* = 0.
- **The thickening is early**: 70–99 % of the land's Σ Δs_DS arrives in the first half of C1's steps.
- **The age overshoots the target from the other side**:
  - the mountain share falls from 52–73 % (B0) to 0.6–11.6 % already at 20 Myr, and to 0.3–6.2 % at 300 Myr;
  - the old continent is a plain at sea level (median 3–50 m, plains 76–97 %), not Europe (195 m);
  - between A = T and 20 Myr lies the mountain band 16.5–37.2 %, unmeasured.
- **The macro success is 0 / 3.** No stop rule fires as written (peak < 1.5 km on 2 / 3 at 256²). A is not tuned: it is the
  author's choice.

**Files**:
- Declared before any measurement: `f170_declared.md`. Predictions: `f170_predictions.md` (written first).
- The bench `crates/ymir-core/tests/f170_cascade.rs`:
  - `f170_cascade` → `bench_output.txt`;
  - `f170_cost` → `bench_cost.txt`, the torus / window cost run alone.
- `images/`.
- Checks: `checks_before.txt`, `checks_after.txt`.

## Part 0

- F169 committed (`c382299`) with its « At commit » block: the torus as default, the window-mode comparisons before F169, the
  reviewer's two-sided rim band, F168's window coast artefact, the selection impossible with 8 plates.
- **The torus is now the cascade's default**: `CascadeConfig::declared` and `AmpConfig::declared` carry `torus: Some([0, 0])`;
  the viz ticks « tore » by default.
  - The pin test asks explicitly for `torus: None` and is F168 bit for bit: the negative control.
  - The F159–F169 benches would now run on the torus if re-run.
- **The cost on an idle machine**, run alone and twice (the témoin's L2):

  | | window | torus | change |
  |---|---|---|---|
  | physics level | 23–26 s CPU | 29–31 s CPU | |
  | chain, CPU | 385–387 s | 520–522 s | **+35 %** |
  | chain, wall | 55 s | 62 s | **+11 %** |

  - The difference is mostly the 4D torus noise, computed in parallel at every level (a hypothesis, not profiled).
  - F169's « under load » numbers (466–554 s) were the torus mode's normal cost, not the load.

## A — the selection with more plates (`profile::select_continent_plates`)

### A1, A2 — the candidates and the criterion

- **The candidates**: N ∈ {12, 16}; k ∈ 5–9 (N = 12) and 6–12 (N = 16); two families × every start plate.
  - The BFS family (F169's).
  - **The compact family** (`clustering::assign_continental_compact`): the adjacent plate whose centroid is nearest, on the
    torus, to the start plate's.
- **The count**: 120 candidates per seed at N = 12, 224 at N = 16. They share 104–110 and 196–215 distinct C1 runs (memoised
  on the plate set).
- **The criterion**: 50–62 % land, circumnavigable with r* ≥ 1, the largest r*.

| N | seed | in 50–62 % | qualify r* ≥ 1 (BFS / compact) | r* = 0 | chosen |
|---|---|---|---|---|---|
| 12 | all | 30–38 | 0 / 0 | 0–2 | none on any seed |
| 16 | témoin | 64 | 0 / 0 | 7 | **none** (best r* 0) |
| 16 | 42 | 68 | 1 / 3 | 6 | compact, k 7, start 3: 51.3 %, **r* 2** |
| 16 | 1 | 61 | 1 / 0 | 6 | BFS, k 10, start 6: 50.4 %, **r* 1** |
| 16 | 9 | 66 | 0 / 2 | 4 | compact, k 7, start 15: 51.2 %, **r* 3** |

- **A3: N = 16** (3 seeds against 0).
- Every qualified continent sits at the low end of the range (50–52 %).
- The selection costs 27–30 s CPU per seed at N = 12, and 47–54 s at N = 16.
- **With r* = 0 the témoin would qualify** (7 candidates at 51–58 %).

### A4 — C1's tectonics (16 plates) against F168's L2

| seed | plates at the end | C-C cells | upper-plate cells | convergent + suture, km per 10 000 km² of land | F168 L2 | ratio |
|---|---|---|---|---|---|---|
| 42 | 6 | 0 | 59 | 473 | 263 | **1.80** |
| 1 | 6 | 58 | 10 | 576 | 241 | **2.39** |
| 9 | 4 | 8 | 16 | 585 | 226 | **2.59** |

- **The reviewer's hypothesis** (« the mountains follow the suture length per area »): the suture length doubles, yet B0's
  mountains are lower than F168 L2's (52–73 % against 75–80 %).
- Not supported in its direction, on 3 seeds.

### A5 — the negative controls (permanent)

- 8 plates without selection is F168 bit for bit (`the_continent_profile_at_production_fraction_is_c1_bit_for_bit`).
- `torus: None` is F168 bit for bit (the pin).
- At 12 and 16 plates, the BFS family from the RNG's first draw is C1's init bit for bit, and the compact family is connected
  with k plates (`more_plates_keep_the_rng_init_and_the_compact_family_is_connected`).

## B — the chronology (`images/f170_thickening_age.png`)

| seed | DS thirds of Σ Δs (steps 1–100 / 101–200 / 201–300) | the first half's share | the last 10 / 25 % | thickening age p10 / p50 / p90 (step) |
|---|---|---|---|---|
| 42 | 665 / 261 / 239 | **70.3 %** | 6.2 / 15.5 % | 22 / 45 / 179 |
| 1 | 818 / 257 / 7 | **99.0 %** | 0.2 / 0.5 % | 24 / 69 / 129 |
| 9 | 853 / 258 / 60 | **91.5 %** | 1.6 / 3.9 % | 24 / 48 / 125 |

- **The C-C collisions are early.**
  - Seed 42: 148 cells at step 10, 37 at 50, 3 at 150, 0 from 260.
  - Seed 1: 115 → 24 → 58, which persist.
  - Seed 9: 142 → 104 at 80 → 8 from 170.
  - The plates merge, from 16 to 4–6.
- **B2 — « the rest » in B0** (23–38 % of the mountain class):
  - altitude 553–756 m; R7 600–714 m; distance to the coast 42–51 km;
  - its share of Σ rebound ≈ its share of the land (15–29 %);
  - **its altitude at t₀: 446–476 m**, the initial crust's +0.45 km.
  - So the reviewer's reading holds: the rest is the initial plateau, dissected.

## C, D — the age (`images/f170_seeds_*.png`, `f170_classes_seeds.png`, `f170_hyps_<seed>.png`)

Europe: plain 51.6 / plateau 4.7 / hill 18.9 / mountain 24.8 %; deciles 40 / 81 / 117 / 153 / 195 / 262 / 383 / 580 / 889 m.

| seed | | mountain | plain + plateau | median land | peak (256² / 2 048²) | land |
|---|---|---|---|---|---|---|
| 42 | B0 | 72.8 % | 0.6 % | 525 m | 3 461 / 3 730 m | 53.1 % |
| | young 20 Myr | **11.6 %** | 43.0 % | 152 m | 3 334 / 2 890 m | 53.2 % |
| | intermediate 100 Myr | 7.1 % | 73.8 % | 63 m | 3 109 / 2 708 m | 53.2 % |
| | old 300 Myr | 6.2 % | 76.0 % | 50 m | 2 928 / 2 695 m | 53.0 % |
| 1 | B0 | 52.2 % | 1.9 % | 327 m | 1 673 / 1 988 m | 54.3 % |
| | young | **0.6 %** | 91.6 % | 46 m | 827 / 1 010 m | 53.9 % |
| | intermediate | 0.3 % | 96.3 % | 5 m | 647 / 823 m | 50.3 % |
| | old | 0.3 % | 96.6 % | 3 m | 592 / 761 m | 46.8 % |
| 9 | B0 | 59.9 % | 5.7 % | 455 m | 3 587 / 3 689 m | 52.8 % |
| | young | **8.1 %** | 53.3 % | 101 m | 1 768 / 2 104 m | 52.9 % |
| | intermediate | 2.6 % | 91.9 % | 17 m | 1 484 / 1 748 m | 51.5 % |
| | old | 2.4 % | 91.4 % | 9 m | 1 251 / 1 504 m | 47.5 % |

- **The mountain share falls with A on 3 / 3**, but almost entirely between T and 20 Myr.
- The old profile's surviving mountains are the late belts: 23–66 % of them were thickened in the last third of the steps,
  against 1–5 % in B0.
- Each decomposition lies within ≤ 31 km of a convergent cell or suture (94–100 %).
- **The plains on the old profile**: 4–46 m high, R7 5–50 m, 18–43 km from the coast: a peneplain at sea level. The land
  shrinks to 47–53 %.
- **Texture**:
  - in the mountain class the slope guard holds on the young profile (p50 0.24–0.32);
  - on the old, it fails on seed 1 (0.13);
  - on the whole world the slopes and octaves fall far under Corsica's from 20 Myr (slope p50 0.006–0.058 against 0.247);
  - the facets fail in the mountain class on 1 / 3 young (seed 42), and on 2 / 3 old.
- **The lakes**: B0 34–1 306 km²; old 10–476 km².
- **The coast on the torus, old profile**:
  - dimension 0.99 / 1.06 / 1.07; length ratio 0.95 / 1.10 / 1.18 (Corsica 1.12 / 1.37);
  - the lowland coasts grow more intricate, towards Corsica;
  - but the aligned share is 0.31–0.34 against F169's 0.27–0.31, which is « degraded > 10 % » on seeds 1 and 9.
- **The regularity on the old profile** (Corsica 391 m: CV_λ 0.61, CV_L 0.87, σ_θ 29.9°, R2_g 0.48):
  - CV_λ 0.62–1.13, CV_L 0.85–1.02, σ_θ 28.3–32.7°, R2_g 0.49–0.74;
  - the confluence angles are now as irregular as Corsica's (F167's baseline was 22–25°).
  - Hovius has too few mountain basins.
- **The cross-sections** (`f170_section_old_belt.png`, `f170_section_recent_belt.png`; seed 42; the old profile red, B0 grey):
  a belt thickened at step 7 against one thickened at step 272.
- **The plates** (`f170_plates.png`): the selected continents are compact clusters of 7–10 of 16 plates.
- **Cost, idle** (per world, C1 included): 465–504 s CPU at every age. The physics level is 25–27 s, unchanged by the
  diffusion sub-steps; the chain 440–478 s.

## R

| rule | result |
|---|---|
| the selection (≥ 3 seeds with r* ≥ 1) | passes with 16 plates (3 / 4) |
| the macro success, old | **0 / 3**: the mountains are too low (0.3–6.2 %); plain + plateau 91–97 % on 2 / 3 (above 84.5 %) |
| « the old profile erases the high mountains » (peak < 1.5 km on ≥ 3) | 2 / 3 at 256² (761 m, 1 251 m), 1 / 3 at 2 048² |
| « mountains > 37.2 % on the old » | 0 / 3 |

- With only 3 seeds selected, no rule's « ≥ 3 » fires: the macro success fails, and no stop is met as written.
- The chronology is reported. **A is not tuned.**

## The reading (hypotheses to check)

1. **The cascade's K wears a relief down in a few Myr** (K_SI 2 × 10⁻⁶ /yr, F159's declaration, an erodible-rock value).
   - Without uplift, the implicit incision grades the channels to base level within a few long steps, and the diffusion
     (√(D t) ≈ 15 km at 300 Myr) planes the interfluves.
   - So any A ≥ 20 Myr erases the early relief, which is most of the relief, since C1 thickens early (70–99 % in the first
     half).
2. **What would keep an old continent's relief is physical, not A**: resistant old crust (the cratons, F172's heterogeneous
   K / K per rock), or a late or continuing uplift.
3. **The rest is the initial crust**, confirmed by its t₀ altitude. On the old profile it is gone (2.8–14 % of a tiny mountain
   class).
4. **The suture density is not the mountains' measure here**: 1.8–2.6× more sutures, fewer mountains in B0.

## Predictions

### The reviewer's

| held | refuted |
|---|---|
| W1's first two clauses (16: 3 seeds; 12: 0) | W1's « compact qualifies more seeds » (2 against 2 seeds; more candidates, 5 against 2) |
| W2's first clause (×1.8–2.6) | W2's second clause (B0 mountains > 75 %: 0 / 3) |
| W3 (70–99 %, 3 / 3) | W6 (median 3–50 m, ×2 of 195 m on 0 / 3) |
| W4 (old mountains < 37.2 % on 3 / 3; plain + plateau > 37.5 % on 3 / 3) | W7 (young facets / walls fail on 1 / 3) |
| W5 (old peak < 2.5 km on 2 / 3 at 2 048²) | |
| W8 (decreasing on 3 / 3) | |
| the meta | |

### Mine

| held | refuted |
|---|---|
| 16 plates 3 seeds; the selection does not fail; the compact r* larger (2–3 against 1) | 12 plates 2–3 seeds (0) |
| W2's first clause | B0 mountains 65–80 % on ≥ 2 (one seed) |
| W3 | the old mountains 10–35 % (0.3–6.2 %); plain + plateau 45–80 % (76–97 %) |
| W4 | the peaks 1.5–2.5 km on ≥ 2 (2.7 / 0.76 / 1.5 km) |
| W8 | W6, W7 |
| the « erase » stop not firing | the macro success on ≥ 3 |
| the meta | the coast unchanged (it moved towards Corsica on the old profile); the costs (the torus chain +35 % CPU, not +15 %; the old physics level unchanged, not +20–60 %) |

## Open, for the author (nothing chosen)

- **(a) The age**:
  - at the declared values the old continent is a plain at sea level;
  - the mountain band lies between T and 20 Myr;
  - A is the author's choice. A choice of A inside that band would be a calibration, hence forbidden.
- **(b) Keep the relief physically instead**: resistant old crust (F172's K, K per rock, the cratons), or a continuing late
  uplift; to be measured with the old profile.
- **(c) The témoin** has no r* ≥ 1 continent at 16 plates (r* = 0 would qualify it); 20 or 24 plates are untried.
- **(d) 16 plates** as the continent profile's plate count.

## At commit (2026-10-11)
- **The torus is the cascade's default mode.** Its cost, measured idle and alone: **+35 % CPU, +11 % wall** on the frozen chain.
- **The selection works with 16 plates** (3 / 4 seeds), and with 12 plates on none. The témoin has only r* = 0 candidates.
- **The reviewer's error, recorded**: the hypothesis « the mountains follow the suture length per land area » is refuted. The
  sutures are ×1.8–2.6, and B0 has fewer mountains (52–73 % against F168 L2's 75–80 %).
- **The old profile is a peneplain at sea level** (median land 3–50 m, against Europe's 195 m).
  - **The reviewer's reading, a hypothesis to check**: with a single soft-rock K and n = 1, every un-uplifted relief vanishes
    in a few Myr. So the age A does not correspond to a geological age as long as K does not depend on the rock.
  - **A is not tuned between T and 20 Myr.**
- **Noted**: the old profile brings the coast closer to Corsica (dimension up to 1.07, Corsica 1.12) and the confluence angles
  too (σ_θ 28–33°, Corsica 27–30°).
- **The default plate count and the order of the next findings** depend on the author's answer. The F171 prompt will follow.
