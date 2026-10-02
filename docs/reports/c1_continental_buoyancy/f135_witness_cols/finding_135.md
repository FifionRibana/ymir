# Finding 135 — the témoin already WAS the author's world; θ's fall is the construction's own based lakes, not run_hd's; ON keeps every lake but leaves a step below its col; the nearest sample lays 99 % of the construction

**Raw outputs, in this folder:**
- `f135_w1.txt`, `f135_w2.txt`, `f135_w3m.txt`, `f135_k.txt`, `f135_km2.txt`, `f135_c.txt`;
- the profiles `k2_profile_lake{1,2,11}.svg` / `.png`;
- the viz layers `v_*.png`;
- `checks_before.txt` and `checks_after.txt`.

The predictions were written before any grep or measurement (`f135_predictions.md`). The instruments were
declared before the measurements (`f135_declared.md`). Two additions came after the first results and are said
so: W2, and K3 / M2.

**Units**: 1 cell = 48.8 m (domain), 0.00238 km².

## Partie 0

**0.1** — **F134 committed** (`5080b92`).
- The meta line now quotes its own wording ("Disagreements: P-0.3, P-D, P-U. At least two are wrong."): it was
  about my three disagreements only, and one of them was wrong.
- Added to ADR F134:
  - O is closed: float noise through the normalisation, which includes the inland basins (an architecture
    fact);
  - T / U: an instrument defect;
  - D: the cut through the cols, invisible to the canyons instrument;
  - the author's look of 2026-10-02.

**0.2** — The legend of F134's Δz maps is in ADR F134 and in `f134_guard_lakes_ocean/maps_legend.md`.
- **Correction found in this round**: those maps are drawn in data rows, so they are **SOUTH up**, not north up
  as the legend first said.
- Lake 2's red basin is therefore **south**-east of the lake. The F134 report text is corrected, with the
  correction marked.

**0.3** — Before and after: guard 6 / 6 (field and lakes), hash `a8d2d538d692c2f0`, lib, check
(`checks_before.txt`, `checks_after.txt`).

## W — the témoin and the viz state

**W1** — every configuration difference (`f135_w1.txt`; JSON diff, recursive):
- **Upscale config: one difference**, `valley_construction.wall_sea_floor_m` (témoin 0.5, viz `null`). The eroded
  digest is `f57954eebde57b00` against the viz's `cf1d4539c7afe588`.
- **Drainage config: one difference**, `infiltration` (benches `null`, viz `{enabled, f_cap 0.70, k_ref 2.7e-3}`).
- **Assembly steps**, read from the code, already named in F134-0.3:
  - the plain `breach_monotone` against `run_hd`'s `breach_monotone_protected` with the active craters' mask;
  - no C-2 crater-lake pass;
  - no H-1 infiltration field.

**W2** (added) — **`wall_sea_floor_m = 0.5` is inert on this world.**
- The témoin's field and the viz state's are **bit-identical**: both `a8d2d538d692c2f0`, 0 cells differ.
- The construction alone (`carve(S1)`) also differs on 0 cells: the viz state's construction lays **no** land
  cell at or below the sea, so the clause never binds.
- ⇒ **No F131–F134 terrain number carries a field difference** (θ, R8, relief, coast, Δz, the residual, the
  cuts).
- The assembly differences carry the **lake counts and areas** of F131–F133 (+1 crater lake, +1.19 km²; F134
  already corrected this). They carry nothing else: canyons read 0 either way.

**W3** — F133's table on the realigned témoin (viz state + run_hd tail), under the guard:

| | OFF | F132 | extended |
|---|---|---|---|
| guard field / lakes | **Match / Match** | NoReference | NoReference (= F134's ON fingerprint) |
| canyons | 0 | 0 | 0 |
| spurs near a coastal wall (/km) | 0.0142 | 0.0142 | 0.0142 |
| R8 terrain | 0.0452 | 0.0395 | 0.0411 |
| R8 network (chord 8) | 0.5613 | 0.5757 | 0.5658 |
| lakes (km²) | 26 (3 492.8) | 34 (4 454.0) | 34 (4 446.2) |
| relief p50 | 515.4 m | 551.9 m | 546.6 m |
| Δz between flat resolutions (outside the present lakes) | 1 150 588 | 235 602 | 12 509 |

- Every terrain number equals F133's to the digit, as W2 requires.
- The lakes are +1 against F133's assembly (the crater lake).

**W4** — the mur ↔ mer clause is a gated candidate.
- The viz has a "Mur ↔ mer (F126-B)" checkbox (`wall_sea_floor_m = Some(0.5)`), off by default, so the guarded
  digest is unchanged.
- Nothing is promoted and no look is asked. On this world it is inert anyway (W2); it binds only on the worlds
  of F126 (B1 / B2 variants).

## M — trunk metrics without the lakes

**As asked**: the mask is the union of the three worlds' final lakes (`run_hd` tail, after the crater pass). The
instrument is F133's, on the OFF skeleton's 70 720 trunk links ≥ 10 km².

| | OFF | F132 | extended |
|---|---|---|---|
| θ (i) as built, all → without the lakes | 0.338 → **0.467** [0.462, 0.470] | 0.243 → **0.386** [0.377, 0.395] | 0.246 → **0.395** [0.387, 0.404] |
| θ (ii) + run_hd breach | 0.369 → **0.466** | 0.297 → **0.389** | 0.301 → **0.398** |
| θ (iii) the world | 0.351 → **0.453** [0.446, 0.460] | 0.245 → **0.356** [0.345, 0.367] | 0.243 → **0.363** [0.351, 0.374] |
| θ on the carved links | 0.493 → **0.494** | 0.436 → **0.433** | 0.439 → **0.435** |
| under-law share | 32.1 → **12.4 %** | 49.6 → **19.1 %** | 45.6 → **14.3 %** |

**With run_hd's lakes removed, the fall REMAINS.**
- θ (iii): OFF 0.453, ON 0.356 / 0.363, a difference of 0.09–0.10.
- The under-law share rises by 6.7 / 1.9 points.
- This, as measured, contradicts F134's "it is the lake beds".

**M2** (added, `f135_km2.txt`) — the mask that removes it is **the construction's own based set**, the input
lakes ∪ the closed depressions (the lakes the χ walk sees). On `carve(S1)`:

| mask | θ all trunk links (OFF / F132 / extended) | θ carved links | under-law share |
|---|---|---|---|
| none | 0.338 / 0.243 / 0.247 | 0.493 / 0.436 / 0.439 | 32.1 / 49.6 / 45.6 % |
| union of run_hd's final lakes | 0.467 / 0.386 / 0.395 | 0.494 / 0.433 / 0.435 | 12.4 / 19.1 / 14.3 % |
| **the based set** | **0.472 / 0.471 / 0.489** | **0.493 / 0.494 / 0.493** | **10.9 / 10.9 / 3.9 %** |
| both | 0.477 / 0.476 / 0.489 | 0.494 / 0.494 / 0.494 | 9.4 / 9.4 / 4.0 % |

- **Without the based set, θ is the same in all three states**: 0.493–0.494 on the carved links, within each
  other's CI.
- The F132 under-law share equals OFF's (10.9 %). The extended share is lower (3.9 %).
- ⇒ **F134's conclusion holds with the right mask; the law is not deformed.** The fall is carried by trunk links
  inside or touching the construction's present lakes.
  - 31 593 links touch the based set; only 26 418 of them also touch a final lake.
  - About 5 200 links lie in input lakes or depressions that are not final lakes (shrunk, absent, or merged):
    there ON's law is the lake level and the terrain is the bed.
- **The instrument defect is the mask**: the run_hd mask (the lakes that exist at the end) is not the set the
  construction bases on.

## K — the cuts through the cols (realigned témoin)

**K1** (`f135_k.txt`) — 35 items: 25 input lakes and 10 closed depressions.
- The col is the lowest cell of the item's outer ring (on S1, or on the breached S1 for depressions).
- Two cuts are measured: level − the lowest built cell of the ring, and level − the built floor at the col cell
  itself.

| | ring cut > 10 / > 100 / > 500 m | at the col cell > 10 / > 100 / > 500 m |
|---|---|---|
| OFF | 17 / 13 / 6 | 13 / 11 / 4 |
| extended | 18 / 12 / 5 | 8 / 5 / 1 |

- **The canyons instrument counts 0 of them** (W3: canyons 0 in every state). It looks at closed depressions,
  and a cut col leaves none.
- The largest OFF cuts (ring; at the col): input lake 2 (57.8 km², col 1 339.6 m) 1 228 / 900 m; lake 14 1 116 /
  1 077 m; lake 8 839 / 473 m; lake 10 785 / 364 m; lake 3 761 / 709 m.

**K2 — the profiles** (`k2_profile_lake{1,2,11}.svg`, PNG renders beside):
- **In OFF**, the construction lays one continuous floor from the sea **through both lake bowls to the
  sources**: a canyon about 1 km under the terrain along the whole path.
  - Lake 2's path: OFF 154 m where S1 and ON are at 1 330 m.
- **In ON extended**, each lake keeps its surface. Lakes 1 and 11 sit exactly at the S1 col of their own
  footprint (327.2 and 362.0 m); lake 2 at 1 292.0 m.
- **But the outlet valley downstream of each col is laid from the NEXT base downstream** (the sea, or the lower
  lake), with χ counted from that base. **It starts with a step of hundreds of metres right below the col**:
  - about 300 m under lake 1's col (km 40–41);
  - about 280 m under lake 11's (km 41);
  - about 900 m out of lake 2 (km 7–8.5 on lake 1's profile, down to lake 1's level + kχ).
- K1's ring cut in ON measures that step at the foot of the old col. The lowering of the lake itself is small:
  - input lake 1's spill 342.6 m → final 327.2 m (15 m);
  - input lake 2 1 339.6 m → final lake 2 1 292.0 m (48 m).

**K3 — lake 2's over-deepening by distance.**
- The first instrument (distance to the S1 col of the final lake 2) reached 27 cells only: OFF's drainage
  leaves through its own trench, 0.5 km from that col. It is redone to the lake's footprint (`f135_km2`).
- On the cells carved in **both** worlds, Δz = EXT − OFF is **1 183, 1 182, 1 171, 1 171, … 1 171, 1 200 m
  from 0 to 18 km**: uniform to about 1 %.
- It falls at 18–22 km (249, 118 m), where ON's network leaves lake 2's base.
- On all cells the p50 swings (134–1 056 m) because uncarved cells mix in.
- ⇒ **Uniform on the carved network.** This is the χ law's constant offset: floor_ON − floor_OFF =
  (b_ON − b_OFF) − k·χ_OFF(col), the same for every cell upstream.

**K4 — the extended cuts > 10 m** (ring metric), 18 items:
- **bases on a lake absent at the end**: 4 (input lakes 8, 12, 14, 21; 5 absent bases in all);
- **present, based lakes**: 14 (1, 2, 3, 4, 7, 9, 10, 11, 13, 15, 16, 18, 23, 24). They are the outlet steps of
  K2.

The 8 unbased final lakes are 1000001, 1000002, 1000006, 1000007, 1000010, 1000014, 1000015 and the crater lake
2000001.
- Their rings are cut by 0–1 m in ON, except 1000001: 500 cells at 3 121.6 m, ring cut 2 968 m in OFF and in ON
  alike. It sits on a ridge reached by a valley in both states: named, not explained.

## C — carve_diag, attribution without remedy

**Grep (rules 11 / 11c).**
- **"Nearest" in the construction**: the earliest ADR hit is **Finding 121** ("the NEAREST skeleton sample (a
  minimum over every cone flattened the long profile)", ADR:17966). It came into the code with `338e682`
  (Findings 113–122).
- **Finding 127** names it as the dam's cause ("`carve` lays a cell from its NEAREST sample", ADR:18565).
- **Finding 128-A**: the priority by area ("the LINE with the largest drained area lays the cell from its nearest
  covering sample", ADR:18612).
- **Finding 129**: "the construction lays its law exactly where it cuts".
- **`carve_diag`**: introduced by `595e7b6` (Finding 130); ADR:18704.
- **`line_samples`**: ADR:18945 (Finding 133); code `valley_construction.rs:864`.

**Measures.**
- The cone is rebuilt in the bench. It reproduces `out` **bit for bit on every carved cell** (0 mismatches in
  all four constructions).
- **Constructed cells of the realigned témoin: 2 697 596.**
  - **98.91 % are laid by the nearest sample** (the `who` propagation); **1.09 % by the cross-line minimum**.
  - On the floor: 323 093 / 7 791. On skeleton trunk cells: 46 096 / 2 032.
- **The residual of path dependence** (12 509 cells):
  - laid by the nearest sample in both A and B: **89.7 %**;
  - carved on one side only: 6.2 % + 3.5 %;
  - a cross-line minimum on either side: 0.6 %.
- **Finding 127's dams** (B2 → A_c realigned, each true col ± 1 cell, 36 cells): **15 nearest (42 %), 21
  cross-line minimum (58 %)**.
  - Body 3 and body 13 are laid by sub-km² lines (0.28–0.50 km²).
  - Body 7 and body 14 are laid by their trunk (69.7, 82.5 km²).
  - At the dams the cross-line minimum, rare elsewhere, carries the majority.

## V — the relief made visible (viz, the author's tool)

The workspace's canvas toolbar has a "🗻 Relief" menu when the Relief layer is shown: Hypsométrie / Ombrage /
Différence.
- **Ombrage**: a Lambertian hillshade on the conditioned field in metres (light azimuth 315°, altitude 45°),
  multiplied into the hypsometric colour.
- **Différence**:
  - "Mémoriser ce monde comme référence" stores the current world; the view shows current − reference in
    metres;
  - a blue–white–red signed scale, its saturation chosen (±10 / ±100 / ±1 000 m);
  - the legend shows the scale, "saturé au-delà de ± X m", and the counts of saturated and ≠ 0 cells.
- **The hash is unchanged**: the views read `hd.eroded` and never write it (asserted in `f135v_relief_capture`).
  The guard reads 6 / 6 after the change.
- **The capture** comes from `layer_color_image`, the very function the viz draws with, cropped on lake 2.
  **Orientation: north up**: the crop is taken after the buffer's row flip, as the viz shows it. F134's `d_lake*`
  maps are south up. It is **not a screen grab**; the author's own capture in the app is still to take.
  - `v_shade_lake2_OFF.png` against `v_shade_lake2_ON.png`: OFF's incised network reads at once (long shadows,
    deep valleys). ON is a smooth plateau around the lake.
  - `v_diff_lake2_1000m.png`: lake 2 and its upstream saturate red, 46 921 cells beyond ± 1 000 m (of 42.18 M
    ≠ 0).
  - `v_diff_lake2_100m.png`: the same at ± 100 m.

## Predictions

**Mine:**
- **P-W1**: held. One field difference plus three assembly differences (the drainage config's infiltration, the
  breach, the crater pass).
- **P-W2**: "the field difference carries every témoin number" is **refuted**: the field is bit-identical and
  carries none. The lake counts held.
- **P-W3**: lakes 26 / 34 / 34, canyons, R8, relief and the residual held. "Spurs up > 50 %" is **refuted**:
  identical.
- **P-M**: **refuted with the asked mask** (Δθ (iii) about 0.1). It holds with the based set (M2, added after).
- **P-K1**: OFF ≥ 8 cuts > 100 m and ≥ 3 > 500 m held; canyons 0 held. ON "≤ 5, only absent bases / unbased",
  and "present lakes' cols ≤ 1 m", are **refuted**: the outlet steps.
- **P-K3**: held, uniform at 1 %.
- **P-C**: "> 95 %" held (98.91 %); the residual ≈ 90 % held (89.7 %). F127's dams "all of them" is
  **refuted** (42 %).
- **P-V**: held.
- **Meta**: not held. Of my three disagreements, none is wrong: W1 held, W3's lakes held (26, not 25), C held.

**The reviewer's:**
- **W**: "the only difference" refuted; "lakes and canyons unchanged" half (canyons held, lakes +1 against F133);
  "the coast moves > 10 %" refuted (identical).
- **M**: with the asked mask, θ (iii) < 0.02 refuted (0.09–0.10); < 2 points refuted for F132 (6.7), held for
  extended (1.9). With the based set (M2): θ held; the share held for F132 (0.0), refuted for extended (−7.0).
- **K**: OFF "≥ 5 cols cut > 100 m, canyons count none" held (11 at the col, 13 on the ring; 0 counted). "ON:
  no present lake's col cut > 10 m" refuted (the outlet steps; 8 at the col cell).
- **K3**: held (uniform on the carved network over 18 km).
- **C** "10–40 %" refuted (98.91 %).
- **Meta**: held.

## Limitations, stated

1. **K's col** is the lowest ring cell of the INPUT item. Input and final lakes differ (input lake 2 is not final
   lake 2), so a K1 line and a K2 lake are not always the same body. K2 reads final lakes.
2. **The step below each col** is read on three profiles and on K1's at-col counts. Its height is not tabulated
   lake by lake.
3. **M2 and W2 were added** after the first results. M2's θ is on `carve(S1)` (the construction), not on the
   pipeline's stages.
4. **V's image** is the layer's own buffer. The interactive view was compiled and its hash checked, not
   screen-captured.
5. **1000001's 2 968 m ring cut** is named, not explained.

## State

**Uncommitted**, awaiting the go-ahead.
- `workspace.rs`: the V layers, the W4 toggle, `f135v_relief_capture`.
- `tests/common/mod.rs`: `bench_upscale_cfg`.
- `f126_coast.rs`: `f135_*`.
- The F134 corrections: the maps are south up (ADR, `maps_legend.md`, `finding_134.md`).
- This folder, and ADR Finding 135.
