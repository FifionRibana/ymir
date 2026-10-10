# F170 — declared before any measurement (2026-10-10)

**More plates for a continent that can be sailed around, then the chronology and the continent's age.**

**Production is unchanged**: C1's default stays 8 plates. `num_plates` becomes the C1 « continent » profile's second parameter;
the guard 6 / 6 is not touched.

**Frozen**: T, K, Airy, the sources' drive, F164's cascade, Davis-Suppe, the torus. The ages A are design values.

**Not blind**: F169's readings (the 8-plate candidates) and the code. No 12- or 16-plate world has been run.

**Part 0**: F169 committed (`c382299`) with its « At commit » block: the torus as the cascade's default, the pre-F169
comparisons in window mode, the reviewer's two-sided rim band, the window artefact on the coast, the selection impossible with
8 plates.

**The torus as the cascade's default** (this round): `CascadeConfig::declared` and `AmpConfig::declared` carry
`torus: Some([0, 0])`, i.e. an input in C1's frame.
- **The negative control** (permanent): `torus: None` is F168 bit for bit. The F159–F168 pins move to explicit `torus: None`.
- **The viz**: « tore (F169) » ticked by default.

## A — the selection with more plates

### A1 — the candidates

For N ∈ {12, 16} plates (the same Voronoi generator; plate p's velocity is preset_phase_1_1's base[p mod 8], production's
preset unchanged), and for every number k of continental plates in k ∈ [⌈0.35 N⌉, ⌊0.75 N⌋] (N = 12: k = 5 … 9; N = 16:
k = 6 … 12; the 8-plate runs put the end land ≈ the continental share + 0–8 points, so this bracket covers 50–62 %):

- **The BFS family** (F169's): C1's BFS from each start plate, stopped at k plates.
- **The compact family** (new): from each start plate s, add one at a time the plate adjacent to the group whose centroid is
  nearest to s's centroid.
  - Centroids are circular means of the plate's cells per axis.
  - The distance is on the torus.
  - Ties go to the smaller plate id.
- **Count**: N × |k| × 2 candidates per seed, i.e. 120 at N = 12 and 224 at N = 16.
- **Shared runs**: candidates with the same plate set share one C1 run (memoised on the set).
- **Each candidate**: C1 under the Airy sources' record, its end land at 64² (h_src > 0, whole torus), and F168's
  circumnavigation test.

### A2 — the criterion

- The end land is within 50–62 %, circumnavigable with **r* ≥ 1** (a sea ring of ≥ 3 cells of 64², ≈ 19 km); the largest r*.
- **Ties**, in order:
  1. the start plate's rank in C1's RNG draw order (for that N);
  2. BFS before compact;
  3. the larger k.
- **Also reported**: the candidates that qualify with r* = 0, and per seed the best r* reached within 50–62 %.

### A3 — the choice of N (declared rule)

The N that qualifies the most seeds (r* ≥ 1). On a tie, 12 (fewer sutures).

### A4 — C1's tectonics for the chosen N

- The plates at the end.
- The C-C collision cells (F165's rule).
- The active margins (upper-plate cells).
- **The convergent + suture length per land area**: the number of 64² cells that are convergent at the end or were
  convergent at some step (`suture`), × 6.25 km, per 10 000 km² of end land.
  - Compared with F168's L2 (the same measure on the 8-plate L2 worlds, recomputed in this bench).
  - The reviewer's hypothesis: the mountains follow this length.

### A5 — the negative controls

- N = 8 without the selection gives F168 bit for bit: the existing test.
- `torus: None` gives F168 bit for bit: the pins.
- The 12- and 16-plate entry points with the start plate at C1's first RNG draw and the BFS family are C1's init at that N bit
  for bit (a permanent test).

**The cost on an idle machine** (Part 0's request): one torus chain and one window chain on the same world (the témoin's L2),
CPU s.

## B — the chronology (F169-B, on the chosen worlds)

- **B1**:
  - Σ Δs per term (Davis-Suppe, subduction, rift thinning, transport) per C1 step, over the end land;
  - the thickening age per land cell (the first step at which its cumulative Σ Δs_DS reaches 50 % of its total; none if the
    total is ≤ 0.05), and its map;
  - the shares received in the last 10 / 25 / 50 % of the steps;
  - the C-C collision cells at each snapshot.
- **B2 — « the rest »** (mountain class outside the cratons, the Davis-Suppe zones Σ Δs > 0.05 and the arcs, B0): altitude, R7,
  the distance to the coast, its share of Σ rebound, its altitude at t₀.

## C — the age (F169-C, unchanged)

- A ∈ {T (B0), 20, 100, 300 Myr}: d_k = d_N r^(N−k), d_N = T/300. One physics step per C1 step; the diffusion's sub-steps
  max(4, ⌈w / 0.2⌉).
- A = T is linear time bit for bit (the permanent test).

## D — the measures

The 4 seeds (témoin, 42, 1, 9) × B0 and the 3 ages; physics 256²; the frozen chain to 2 048²; torus on; the display offset at
the end.

- **The selection**: N, k, family, land, r*, start plate.
- **Macro**:
  - the classes against Europe;
  - the land deciles against Europe's (40 / 81 / 117 / 153 / 195 / 262 / 383 / 580 / 889 m);
  - the peaks;
  - the mountain class decomposed by thickening third (steps 1–100, 101–200, 201–300), the cratons, the arcs and the rest;
    and by the distance to a convergent cell or suture (≤ 31 km, 37–94 km, > 94 km);
  - **the mountain share against A**.
- **Plateaus and plains**: the share, the distance to the coast and to a convergent cell or suture, the altitude, R7.
- **The lakes.**
- **Texture**: Corsica's battery (whole and mountain class), F by class, the slope guard.
- **The coast on the torus** against Corsica.
  - « Degraded > 10 % » is measured against F169's L2 coasts on the torus: the same seed's F168 L2 world, torus on, linear
    time, recomputed here.
- **Regularity**: F167's baseline on the old profile.
- **Images** (north up, checked; the offset at the end):
  - the 4 seeds × B0 | young | intermediate | old, shaded and classes, at 256² and 2 048²;
  - the chosen plates' map;
  - the thickening-age map;
  - the hypsometric curves against Europe;
  - the old profile's cross-sections of an old belt and a recent belt.
- **The cost**: CPU s per world per age, idle machine.

## R — the stop rules (the brief's)

- **The selection fails** (< 3 seeds qualify with r* ≥ 1, for 12 and for 16): stop and report what r* = 0 would give and the
  best r* reached.
- **Macro success on the old profile**, on ≥ 3 of 4:
  - mountains 16.5–37.2 %;
  - plain + plateau 37.5–84.5 %;
  - circumnavigable with r* ≥ 1;
  - no coast degraded > 10 % against F169 (torus).
  - Then F171.
- **The old profile erases the high mountains** (peak < 1.5 km on ≥ 3): report B1. A is not tuned.
- **Mountains still > 37.2 % on the old profile**: decompose by thickening age and by distance to the sutures, and stop.
- **Forbidden**: as F165–F169.
