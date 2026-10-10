# Finding 169 — a world that does not depend on its framing, a continent that can be sailed around, and the continent's age

**Status: built, partly measured, gated. Production unchanged. Nothing committed without the feu vert.**

**Verdict**:
- **The infrastructure passes.**
  - The torus cascade is translation-invariant **bit for bit**: in the permanent unit test, and at full scale on 4 / 4 seeds
    × 3 shifts, at 256² and at 2 048².
  - The map-border rim is gone. The coast is measured on the torus.
- **The plates' selection fails on 4 / 4 seeds**: no BFS cluster with 50–62 % land is circumnavigable. The stop fires: next,
  more plates. B, C and D were therefore **not measured**.
- The age mapping (C1, C2) and the viz controls (C4) are built and tested, but no age was run on a world.

**Files**:
- Declared before any measurement: `f169_declared.md`. Predictions: `f169_predictions.md` (written first).
- The bench `crates/ymir-core/tests/f169_cascade.rs`:
  - `f169_cascade` → `bench_output.txt` (it stopped at the selection);
  - `f169_infra` → `bench_infra.txt` (A1–A3 at full scale, infrastructure only).
- Checks: `checks_before.txt`, `checks_after.txt`.

## Part 0

F168 committed (`3271a44`). Its ADR carries the reviewer's error, the plate quanta, the invariance defect, the coast measure
rule and L1 / L2 as profile constants. Its « At commit » block makes the invariance defect and the reviewer's error explicit,
with the author's three answers verbatim.

## A1 — the fully periodic cascade (`CascadeConfig::torus`, `AmpConfig::torus`; `None` = F168 bit for bit)

### What was made periodic or anchored on the terrain

The inventory is in `f169_declared.md`. What changed:

- **The solver**: `incise_with_floor_torus` (production's `incise_with_floor` calls it with `None`).
  - No D8 receiver leaves the map; the base level is the sea alone.
  - The lateral erosion, the talus (`talus_sweep_torus`) and the linear diffusion (`linear_diffusion_torus`) are periodic.
  - The talus sweep's height ties, the random receiver's hash and the MFD accumulation's ties (`mfd_accumulation_framed`)
    read the terrain index.
  - `StreamPowerConfig` is not touched: it enters the eroded cache key.
- **The amplification**:
  - periodic neighbours in the routing weights, the steepest descent, the exact area, the talus and the deposition;
  - the ocean is the largest periodic sea component (`ocean_mask_t`);
  - the noises (π, the talus ν, the hardness, the warp, the physics roughness) are periodic on the torus (`fbm_torus`, a 4D
    Clifford-torus embedding of the same OpenSimplex sources) and sampled at the terrain's coordinates;
  - π's and the roughness's land mean / RMS are summed in terrain order.
- **The rebound's FFT** runs in the terrain frame.
- **Found by the test** (my prediction « not at the first run » held): the warp sampled u(x + d) in f32, so its fractional
  part depended on the cell index's magnitude (≈ 3 × 10⁻⁴ m). `warp_grid_t` takes x + ⌊d⌋ and d − ⌊d⌋ apart. The diagnostic
  `f169_debug_amp_stages` (ignored) is kept.

### The invariance test (declared tolerance: bit for bit)

- **The permanent unit test** `the_torus_cascade_is_translation_invariant_bit_for_bit`:
  - a synthetic island straddling both seams, with the roughness, the random receiver and the rebound on;
  - a 64 → 128² history physics level, then one 256² amplification level (warp, π, R4);
  - shifts (5, 0), (0, 23) and (37, 50): **0 differing cells**.
- **Full scale** (`f169_infra`; F168's L2 worlds, whose land touches the window; the history physics level 256², then the
  frozen chain to 2 048²): **0 differing cells, max |Δ| 0 m, on 4 / 4 seeds × 3 shifts, at 256² and 2 048²**.
- **The negative control**: `None` reproduces F168. All 645 earlier lib tests pass, the pins included (650 with F169's).

## A2 — the rim

The cascade computes in C1's frame, and the viz's offset is applied at the end. Measured on F168's L2 worlds with the torus
on: the border land's mean altitude over the land's mean, at 2 048².

| seed | F168 (window, rim) | F169 (torus), display window | C1's frame border |
|---|---|---|---|
| témoin | 1 805 m / 896 m = 2.01 | 1.21 | 1.21 |
| 42 | 1 483 / 1 148 = 1.29 | **0.61** | 0.63 |
| 1 | 1 049 / 776 = 1.35 | 0.84 | 0.84 |
| 9 | 1 360 / 929 = 1.46 | **0.72** | 0.72 |

- The bench's printed label quotes F168's 256² range (1.15–2.1); the table uses F168's 2 048² values.
- **No border is raised any more.** Invariance makes the border a place like any other: the output follows a shift exactly.
- **But the declared reading (0.8–1.25) holds on 2 / 4.** Seeds 42 and 9 read under 0.8: the window's edge crosses land lower
  than the mean (a wrapping continent cut through its lowlands). That is not a rim; still, the band as declared was
  two-sided, so it is reported as missed.

## A3 — the coast on the torus (`physio::coast_torus`)

| seed | torus: dimension / length ratio / aligned / block index | window (F168's measure) aligned / block |
|---|---|---|
| témoin | 0.99 / 0.92 / 0.31 / 0.69 | 0.49 / 7.6 |
| 42 | 0.98 / 0.91 / 0.31 / 1.16 | 0.55 / 10.4 |
| 1 | 1.00 / 0.94 / 0.28 / 0.88 | 0.55 / 11.0 |
| 9 | 0.99 / 0.90 / 0.27 / 0.79 | 0.56 / 11.0 |

- Corsica: 1.12 / 1.37 / 0.27 / 0.78.
- **F168's « degraded coast » was the window**: on the torus, the L2 coasts have Corsica's aligned share and block index.
  They remain smoother (dimension ≈ 1.0 against 1.12, length ratio ≈ 0.9 against 1.37), as in every round since F165.
- Unit test: a band counts 2n edges, no window border; a disc across the seam reads as the centred disc.

## A4 — the continental plates' selection (`profile::select_continent`)

- **The setup**: 16 candidates per seed (f = 0.625 and 0.5, each of the 8 start plates); C1 under the Airy sources' record;
  the end land at 64².
- **The entry point**: `init_c1_state_phase_2_r7_start`, with `start_plate_order` = C1's RNG draw order.
- **The negative control** (permanent): `the_start_plate_at_the_rng_draw_is_c1_bit_for_bit`.
  `ContinentalClusterParams` is unchanged (it enters the tectonic key).

| seed | candidates in 50–62 % | of them circumnavigable | the circumnavigable candidates (f, start: land, r*) |
|---|---|---|---|
| témoin | 5 | **0** | 0.5 / 5: 47.1 %, r* 0 · 0.5 / 0: 44.0 %, r* 0 |
| 42 | 10 | **0** | none of the 16 |
| 1 | 5 | **0** | 0.625 / 3: 46.9 %, r* 0 · 0.5 / 7 and 0.5 / 3: 31.2 %, r* 5 · 0.5 / 0: 40.5 %, r* 4 |
| 9 | 4 | **0** | 0.5 / 5: 29.2 %, r* 8 · 0.5 / 3: 35.5 %, r* 4 |

- **No candidate qualifies on any seed.** With 8 plates, every BFS cluster that reaches 50 % of the map wraps the torus.
- The circumnavigable ones stop at 47 %. Those near 47 % have r* = 0: a sea loop of a single cell, 6 km.
- **The selection costs 4–5 s CPU per seed** (16 C1 runs).

## B, C, D — not measured (the stop)

- The brief's rule: « si la sélection échoue sur au moins 2 graines, arrête-toi et rapporte ». The chronology, « the rest »,
  the ages and every physical measure were not run.
- **Built and tested, ready for the round with more plates**:
  - **C1 — the time translation**: `history::step_durations` and `HistoryRun::age_years`. One physics step per C1 step; the
    step's uplift is the interval's Δh × d_k / D_j.
    - The permanent test `the_age_at_t_is_linear_time_bit_for_bit` holds: A = T is linear time bit for bit, an older A is
      not, and r(300 Myr) = 1.019306.
  - **C2 — the time step**: the diffusion's sub-steps max(4, ⌈w / 0.2⌉). On the oldest step:

    | age | d_1 | w | sub-steps | the incision's K·A^m·dt/dx (A = 100 000 km², implicit, indicative) |
    |---|---|---|---|---|
    | 20 Myr | 0.163 Myr | 0.05 | 4 | 66 |
    | 100 Myr | 1.47 Myr | 0.47 | 4 | 593 |
    | 300 Myr | 5.70 Myr | 1.82 | 10 | 2 307 |

  - **C4 — the viz** (the Cascade window):
    - the slider « Âge du continent (Ma) », 5.62–300, logarithmic, with the presets linéaire / jeune / intermédiaire / vieux;
    - « tore (F169) »;
    - « sélection (F169) », which falls back to L2 when nothing qualifies;
    - the view « Âge d'épaississement », in the history's term views.
  - **The B bench** (`chrono`: Σ Δs per term and per step, the thickening age, the C-C collisions per snapshot) is written in
    `f169_cascade` and runs after a successful selection.

## Cost

- **The full-scale chains with the torus on** (L2 worlds): 466–554 s CPU, against 352–448 s in F168's L2.
- Part of it is the periodic ops. Part is machine load: a viz `cargo check` ran alongside.
- Not separated; to remeasure alone in the next round.

## R

- **Infrastructure**:
  - the invariance test passes, bit for bit, at both scales;
  - **the selection fails on 4 / 4 seeds → stop. Next: more plates** (the author: « SI ça marche pas on fera plus de
    plaques »).
- The macro stops (success, erased mountains, mountains > 37.2 %) were not reached.

## The reading and proposals (hypotheses to check; nothing chosen)

1. **More plates** (the author's next step).
   - With 12 or 16 plates the land quantum falls from ~12 % to ~6–8 % of the map, and the clusters are built of smaller,
     more numerous cells. They could reach 50–62 % with less elongation.
   - The selection and the B / C bench run unchanged with `num_plates` as the profile's second parameter.
   - **Hypothesis**: on a 400 km torus, F168's geometric limit (≤ 97 km of strait for a square at 57.5 %) leaves room for a
     compact cluster, not for an elongated one.
2. **To decide with the author**:
   - whether r* = 0 (a one-cell, 6 km strait at 64²) counts as « on peut en faire le tour ». The declared test says yes.
   - Three seeds have such a candidate at 44–47 %, just under the 50 % floor.

## Predictions

| | held | refuted | not measured (the stop) |
|---|---|---|---|
| **The reviewer's** | V1 (bit for bit on 4 / 4; no raised border) | **V2** (0 / 4) | V3–V7 |
| **Mine** | V1 bit for bit « not at the first run » (the warp's f32 sampling); the small r* (0 near 47 %); the meta | V2 (I said 3–4 seeds); the rim ratio within 0.8–1.25 on 4 / 4 (2 / 4: two read under 0.8) | the rest |

## Open, for the author

- (a) **More plates**: how many (12 or 16), the selection unchanged.
- (b) Whether r* = 0 at 44–47 % land would do in the meantime.
- (c) The torus as the cascade's default, and the viz's offset applied at the end only.

## At commit (2026-10-10)
- **The torus becomes the cascade's default mode**. The display offset is applied at the display and the export only.
  Production is unchanged.
- **Every comparison before F169 was made in window mode** (F159–F168: the solver's off-map receivers, the window-anchored
  noises, the window border as coast).
- **The reviewer's error, recorded**: the rim band 0.8–1.25 was two-sided. A border that cuts through lowlands reads under 0.8
  with no rim (seeds 42 and 9).
- **F168's « degraded coast » was a window artefact** (A3).
- **The selection is impossible with 8 plates** (A4: no 50–62 % cluster is circumnavigable on any seed).
- **The threshold becomes r* ≥ 1** (a sea ring of at least 3 cells of 64², ≈ 19 km); the r* = 0 candidates are reported too.
  12 and 16 plates are tried.
- **The cost**: the torus chains are to be remeasured on an idle machine.
- **The numbering**: F170 = the plates, the chronology and the age; F171 = the sea-level rise; F172 = the heterogeneous K.

**The author's decision, verbatim** (recalled): « SI ça marche pas on fera plus de plaques. »
