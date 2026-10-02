# Finding 136 — the steps are real on the ON water's path (8 > 50 m, 6 > 200 m, 2 > 500 m), invisible to the canyons instrument; the lakes' lowering is their own bed, not the notch; the absent bases are emptied by the construction and the light pass; no hanging junction > 50 m

**Raw outputs, in this folder:**
- `grep_rule11.txt`;
- the benches `f136_k.txt` and `f136_k6b.txt`;
- the profiles `k5_profile_lake{1,2,11}.svg` / `.png`;
- `checks_before.txt` and `checks_after.txt`.

The predictions were written before any grep or measurement (`f136_predictions.md`). The instruments were
declared before the measurements (`f136_declared.md`). **One amendment** came after K5's profiles: K6b, said so.

**Units**: 1 cell = 48.8 m (domain).

**World**: the realigned témoin (the viz state + run_hd tail).
- OFF: field `a8d2d538d692c2f0`, lakes Match.
- ON extended: field `a620a9aa0d5882a0`, 34 final lakes.

## Partie 0

- **0.1** — F135 committed (`bfc3907`), with F134's corrections (maps south up). Added to ADR F135:
  - W, M (the based-set mask added after the results, not blind), C, K;
  - the author's words of 2026-10-02, recorded as given;
  - the reviewer's consequence: a gorge by default, a waterfall on a criterion and tagged; the present steps are
    an artefact, untagged.
- **0.2** — The orientation of F135's `v_*` images: **north up** (crops of `layer_color_image` after its row
  flip). F134's `d_lake*` maps are south up.
- **0.3** — Before: F135's after-checks, reused because no code changed since (declared in `checks_before.txt`).
  After: `checks_after.txt`. F136 changes no production code (benches only).

## G — rules 11 / 11c (`grep_rule11.txt`)

**1. What joins two base levels with a slope, or makes a break retreat?**
- **Nothing in the code joins two bases.** The construction lays each reach from its own base. The breach only
  makes a path monotone (it carves a slot; it never grades a step).
- What exists in the code:
  - **knickpoint retreat by the stream-power incision**: `stream_power.rs:1028` (celerity `K·A^m`); F115 measured
    the wave arriving in a fraction of a pass, held by gates, not by speed;
  - **"Sill incision (H-2)"** (ADR:2664, earliest hit): "the outlet carves its col, the lake retreats". F114
    computed it as exactly zero. It was never built as a construction.
- **The step downstream of the lakes was already MEASURED, in F114.** Its "family 1" (16 lakes) has "the drop …
  DOWNSTREAM (d@10 2.8–465 m, d@200 85–975 m, max S 0.20–1.40)", with F119's amendment of the col. F136 finds the
  same object in the constructed world.
- These terms have **no hits** in ADR and code: "rupture de pente", "marche", "rapids", "recul",
  "outlet incision", "knickzone". "Waterfall" and "chute" appear only in passing.

**2. A lithology or erodibility field?**
- **Yes. An erodibility field exists**, though not a map of rock types.
- The earliest ADR hit (ADR:259, "There is NO lithology / geology / erodibility field in core") is **out of
  date**. The code built C-3 since: `production_k_field` (`production_upscale.rs:662`), a per-cell multiplier.
  - Hard basement ×1, rift-soft ×10, volcaniclastic ×3, from the coarse classes (`lithology::build_coarse_k`).
  - Times the C-3b fracture density.
- **A volcanic mask exists**: the edifices' basal discs (`lithology::stamp_volcanic_k`, the viz's
  "volcaniclastic" overlay). It has 1 806 190 cells in this world.
- No "coulée" (lava flow) and no glacial field.

**3. Can the rivers' export carry a list of tagged points?**
- **Not today.** `rivers.json` (`export/hydro.rs:19`, `RiverSegmentView`) carries per segment its `points` and
  one per-point array, `profile_m`, parallel to `points`, plus per-segment scalars.
- There is no tag or feature field. A per-point tag would be a parallel array of the same kind (a format change,
  not built).

## K5 — the step on the ON water's own path

The path starts at `Lake::outlet`'s receiver (the true col, Finding 119) and runs down the ON drainage's D8. The
profiles are `k5_profile_lake{1,2,11}`: the ON path solid, K2's OFF-skeleton path dashed.

**The step exists on the path the ON water takes, for all three lakes.**

| lake | declared K6 (first steep zone) | amended K6b (to the first graded reach ≥ 2 km) | K2's reading (F135) |
|---|---|---|---|
| 1 | 136 m over 0.33 km | **290 m** at 1.94 km | ~300 m |
| 2 | 413 m over 0.61 km | **862 m** at 3.31 km | ~900 m |
| 11 | 276 m over 1.60 km | **260 m** at 1.55 km | ~280 m |

- **The declared instrument stops at short benches**: lake 2 falls 1 292 → 870 m, rests on a bench of under
  1 km, then falls to 400 m. That is why K6b was added after reading the profiles.
- With K6b the three are within 7 % of K2's readings, read off its plots.
- **OFF on the same cells has no step**: its floor is already low, and it climbs over the ridges the ON water
  crosses.

## K6 — the inventory (14 present lakes measured)

The present lakes are the 26 ON final lakes ≥ 1 km², not crater lakes. **12 of them (the below-sea-basin lakes
1000xxx) have no D8 receiver at their outlet**: they drain by spillways, outside the D8. They are not measured.

| lake | km² | level | declared step H / L | K6b drop | drop over 5 km | max S (2 km) ON / OFF |
|---|---|---|---|---|---|---|
| 1 | 132.8 | 327.2 | 136 / 0.33 km | 290 | 320 | 0.154 / 0.078 |
| 2 | 35.5 | 1 292.0 | 413 / 0.61 km | 862 | 899 | 0.261 / 0.188 |
| 3 | 159.7 | 1 314.5 | 0 | 0 | 0 | 0.349 / 0.193 |
| 4 | 45.8 | 168.2 | 119 / 0.39 km | 119 | 148 | 0.067 / 0.005 |
| 7 | 5.3 | 383.4 | 158 / 0.30 km | 317 | 352 | 0.166 / 0.101 |
| 8 | 12.3 | 520.5 | 58 / 0.22 km | 51 | 317 | 0.114 / 0.059 |
| 9 | 79.7 | 418.0 | 281 / 1.55 km | 266 | 296 | 0.142 / 0.059 |
| 10 | 249.0 | 974.5 | 819 / 3.29 km | 811 | 841 | 0.345 / 0.135 |
| 11 | 163.0 | 362.0 | 276 / 1.60 km | 260 | 286 | 0.138 / 0.079 |
| 12 | 10.2 | 118.9 | 88 / 2.95 km | 0 | 103 | 0.050 / 0.009 |
| 13, 15, 19, 20 | | | 0 | 0 / — / 0 / 24 | 0–64 | ≤ 0.086 |

- **Counts**:
  - declared K6, steps > 50 / > 200 / > 500 m: **9 / 4 / 1**;
  - amended K6b: **8 / 6 / 2**;
  - OFF on the same cells: 0 / 0 / 0.
- Lake 3's step (max S 0.35) lies beyond the declared 5 km start window: its outlet path is graded first, then
  falls.
- **The canyons instrument sees none** (canyons 0 in ON extended, F135-W3).

## K7 — the notch and the lowering (15–48 m): it is the lake's own bed

**For all 14 lakes, the true col (the final sill) is a cell of the input lake's OWN BED.** Its S1 z is already
under the input spill level L_in: the construction does not touch it (S1 = built at the col for all 14).
- The light pass lowers it by 0–3.2 m (lake 11: 365.2 → 362.0 m). The breach does nothing (the level is set on the eroded world).
- **The lowering**: lake 2 47.6 m, lake 7 28.7, lake 11 28.3, lake 10 27.0, lake 3 17.9, lake 1 15.4, lake 15
  14.8, lake 9 14.6; the others 2.8–8.0 m.
- **The notch**: the ring cut of the input body, made by the construction's downstream valley. It is **49–941 m
  deep, and lies 0.4–17 km from the true col**.
  - Lake 1: notch 329 m, 1.1 km away; lowering 15 m.
  - Lake 2: notch 941 m, 2.0 km away; lowering 48 m.
- ⇒ **The lowering is not the notch.** The notch opens a lower exit beyond the old col, and the lake then spills
  over the highest point of its own bed between the deep part and the notch. The lowering is the old col minus
  that bed sill.
- No stage lowers the sill itself (the cone rebuild is never called: no true col is constructed).

## L — the bases on absent lakes

| body | km² | L_in | under water: construction / eroded / final | emptied by | final floor slope p50 / < 0.1 % | median z − L_in |
|---|---|---|---|---|---|---|
| input lake 8 | 5.17 | 1 148.2 | 14 / 0 / 0 % | **the construction** | 0.410 / 0 % | −143.7 m |
| input lake 12 | 5.34 | 723.4 | 57 / 0 / 0 % | **the light pass / erosion** | 0.119 / 0 % | −56.2 m |
| input lake 14 | 8.20 | 1 257.4 | 6 / 0 / 0 % | **the construction** | 0.328 / 0 % | −151.1 m |
| input lake 21 | 6.09 | 163.7 | 82 / 0 / 20 % | **the light pass / erosion** | 0.050 / 0 % | −28.4 m |
| closed depression 3 | 261.73 | −5.6 (\*) | 0 / 0 / 1 % | (\*) | 0.050 / 0 % | +101.2 m |

- The four absent input lakes are small (5–8 km²) and their rims carry deep notches (46–1 126 m).
  - The construction empties lakes 8 and 14: its notch drains the bowl below its bed.
  - The light pass finishes lakes 12 and 21.
- **Their floors are not plains**: 0 % of each footprint is under 0.1 % slope, median slopes 0.05–0.41. The
  floor sits 28–151 m under the base the construction laid: the base is the vanished lake's level, above a
  drained bowl.
- (\*) **Closed depression 3 is an instrument limit.**
  - Its ring touches the inland below-sea basin it surrounds, so "L_in" (−5.6 m) is that basin's floor, not a
    col.
  - The fill instrument on `carve(S1)` treats z ≤ sea as ocean, so it can never read it as water.
  - It is listed, not judged.
- **Their share of F135-M's links**: 5 175 links are in the based set but touch no final lake (reproduced).
  - 1 744 (33.7 %) touch these 5 bodies.
  - **3 102 (59.9 %) touch the cells of OTHER input lakes**: the margins of the lakes that shrank (input
    footprint vs final footprint, e.g. lake 2 57.8 → 35.5 km²).
  - 329 (6.4 %) touch a closed depression only.

## K8 — contexts available without geology (counts only, nothing tagged)

- **Hanging junctions** (ON drainage, tributaries ≥ 10 km², height = z_tributary − z_junction): 452 junctions.
  - Counts > 10 / > 50 / > 200 m: **11 / 0 / 0**.
  - Of the 11 over 10 m: 2 dam-type (laid by a line no larger than the tributary), 4 clean breaks, 5 uncarved.
  - The construction keeps junctions graded (the cross-line minimum and the shared χ).
- **Mouths on a coastal wall**, above the sea (> 1 m) within 2 cells of a wall cell: **0**. The construction lays
  no wall below the sea and every mouth reaches it graded.
- **Volcanic terrain**: 2 183 trunk cells (≥ 10 km²) inside the mask, 69 520 outside.
  - Slope over 2 km inside p50 0.0085 / p90 0.035 / max 0.27; outside 0.0106 / 0.047 / 0.42.
  - **No steeper inside.** The erodibility multiplier acts on the incision, which the construction replaced.
- **K6's steps (> 50 m) by context**: lake outlet only 8, lake outlet + volcanic 1 (lake 7). None at a hanging
  junction or a coastal-wall mouth.

## Predictions

**Mine:**
- **P-G1** held: nothing joins two bases; H-2 unbuilt; F114 measured the drop.
- **P-G2** held: an erodibility field exists (C-3), and a volcanic mask.
- **P-G3** held: no tagged point list.
- **P-K5 refuted**:
  - with the declared instrument, lakes 1 and 2 fall outside ± 30 % (136 / 413 against 300 / 900);
  - lake 11's step is NOT under half (276 ≈ 280);
  - with K6b all three are within 7 %.
- **P-K6**: "6–12 > 200 m" and "≥ 3 > 500 m" refuted with the declared instrument (4, 1). With K6b, 6 held and
  2 > 500 m refuted. "Canyons see none" held; "most lakes ≥ 5 km² have a step > 50 m" held (9 of 14 measured).
- **P-K7** held: the lowering is the input bed, not the notch; the notch is hundreds of metres deep.
- **P-L**:
  - "≥ 3 of 5 emptied by the construction" held, but only 2 cleanly (depression 3 is an instrument limit);
  - "< 30 % flat" held (0 %);
  - "< 50 % of the links" held (33.7 %).
- **P-K8**:
  - "fewer hanging > 50 m than steps > 200 m" held (0 against 4–6);
  - "few dams" held;
  - "mouths ≤ 10, under 20 m" held (0);
  - "steeper in the volcanic mask" **refuted**.
- **Meta**: of my five disagreements (G2, K5, K7, L, K8), one is wrong (K5). **"At least two wrong" is not
  held.**

**The reviewer's:**
- **G**:
  - "nothing joins two bases; only H-2 (F114) speaks of outlet incision" held, with F115's knickpoint retreat
    beside it;
  - "no lithology field, but a volcanic mask" half: the mask exists, the field exists too;
  - "no tagged point list" held.
- **K5** ("the step on the ON path for the three, ± 20 % of OFF's heights"): **refuted with the declared
  instrument** (lakes 1, 2 at −55 %), **held with K6b**.
- **K6** ("≥ 8 > 200 m, canyons none"): ≥ 8 refuted (4 declared, 6 K6b); canyons none held.
- **K7** ("the construction of the downstream valley at the rim; lowering = notch ± 5 m"): **refuted**. The
  lowering is the bed, and the notch is 3–60 times the lowering.
- **L** ("> 80 % of the footprint flat; > 80 % of the links"): **refuted** (0 %; 33.7 %).
- **K8** ("more hanging > 50 m than outlet steps > 200 m"): **refuted** (0 against 4–6).
- **Meta** held.

## Limitations, stated

1. **K6b is an amendment.** It was declared after seeing the profiles that the declared step stops at benches.
   Both are reported.
2. **12 below-sea-basin lakes have no D8 outlet** and are not measured; their outflow is a spillway.
3. **Closed depression 3**: "L_in" and its fill read are wrong by construction of the instrument (named, not
   judged).
4. **K5's comparison with K2** uses K2's heights read off its plots (F135), not a computed K2 number.
5. **Hanging height is one D8 step** on the conditioned (monotone) field. A hanging valley whose drop is spread
   over several cells reads lower.

## State

**Uncommitted**, awaiting the go-ahead: `f126_coast.rs` (`f136_k`, `f136_k6b`), this folder, ADR Finding 136.
No production code changes, **no remedy, no tag**.
