# Finding 140 — the gorge's retreat built gated (off by default); its gate run STOPS at the first world: at ×1, p = 0.5, the HD assembly fires Finding 38's invariant, because the construction's outlet area drains large lakes at the default age and the breach then carves their residual pits below the sea

**Status: STOPPED and reported, by the author's choice.** Nothing is corrected after the failure; nothing is
promoted; nothing is committed.

**Raw outputs, in this folder:**
- `f140_g.txt`: the gate bench, which panicked;
- `f140_diag.txt`, `f140_diag2.txt`, `f140_diag3.txt`: the diagnosis;
- `checks_before.txt`, `checks_after.txt`.

The predictions were written before any code or measurement (`f140_predictions.md`). The instruments were declared
before the measurements (`f140_declared.md`). The specification is **`spec_gorge_age_v3.md`**.

## Partie 0

- **F139 committed** (`2857814`). Added to ADR F139:
  - C1, C2 and C4 were not blind: the round confirms tables more than it tests hypotheses;
  - lake 1000001 is queued as a named defect.
- **The author's decisions (2026-10-02), recorded as given**:
  - when room is lacking the gorge steepens;
  - the per-lake retreat: « Ça peut être un slider ou sélection (0, 0.25 ou 0.5) car là je ne sais pas. Par défaut
    j'irai sur 0.5. »;
  - the head falls as proposed;
  - not contested: map (a), m = 10, the below-sea lakes with the merged spill (except 1000001 and the 3 without a
    spillway segment), H_f's 10 m floor.
- **Spec v3** = v2 + these decisions + A_ref = 418.7 km² declared constant + the G-ring gate. It declares five
  deviations before building (§ 3).
- **Checks before** = F139's after-checks (no code had changed).

## B — what is built (gated, `gorge_retreat: None` by default)

**In `valley_construction.rs`:**
- `GorgeRetreat { r_world, p, m, a_ref_km2 }` (`v3(r, p)`: m = 10, A_ref = 418.7), `GorgeFall`, `GorgeFallKind`,
  `GorgeBody`.
- The gate acts only with `lake_base = InputLakesAndBasins` and `basin_base`.
- **The bodies** (`gorge_bodies`): the input lakes and the closed-depression components ≥ 1 km², each with:
  - L_in (the lake's surface, or the merged spill), L_bed (the minimax from the lowest cell to the outflow) and
    L_floor;
  - the outflow, the col and A (**the col's drained area**: deviation 1);
  - r_lake and L(r_lake);
  - φ (splitmix64 of the lowest cell).
- **The lake is a base at L(r_lake)**: the χ walk's two lake stops take it instead of L_in.
- **The invariant**: after the walk, `base ← max(base, z_gorge − k·χ)` on the outlet path. The samples' floor
  (`line_samples`, `zf = base + age_k·χ`) is then `max(law, z_gorge)` before the cones. The col is at L(r_lake).
- **The rim**: a clamp `v ≥ L(r_lake)` on the outer ring, in `carve_diag`.
- **The falls** (head at the lip, shortage at the base beyond the cap), tagged above H_f = max(1.5·S·100 m, 10 m).

**The plumbing**:
- `UpscaleResult.gorge_falls` → `ErodedProduct.gorge_falls` (a `gorge.json` sidecar, written only when
  non-empty);
- `HdResult.gorge_falls` (with the final lake at each outflow);
- `rivers.json`: a per-segment `features` list (`rivers_json_with_falls`; `rivers_json` is the no-fall call).

**The viz**:
- a "Recul de la gorge (âge)" checkbox and a p selector (0 / 0.25 / 0.5, default 0.5);
- r read from the age selector by map (a); the gorge turns the extended lake base on;
- cyan diamonds for the falls on the Relief layer, with their kind, height and lake on hover.

**Permanent tests** (lib 591 → 593; the key test is an extended existing test):
- `the_gorge_keeps_the_rim_at_the_lakes_level_and_retreats_with_age`. Negative controls first: without the gorge
  the outlet notches the rim by > 10 m; the gate alone is inert, bit for bit. At r = 0 the rim keeps ≥ L_in; at
  r = 2 the col is cut to the floor. The key is absent when off.
- The cache-key test: the gorge and each of r and p move the key; off does not.
- `rivers_json_carries_the_tagged_falls_only_when_there_are_any`: no falls gives the same bytes and no `features`
  key; one fall is exported.

## G — the gates: the run stops at its first world

`f140_g` started with GORGE ×1 (r_world 1, p 0.5). The run_hd tail's assembly panicked:

> ADR Finding 38/92-B: 2 enclosed below-sea component(s) carry no water body. (floor cell, cells):
> [(29094171, 121), (29378947, 1)]

The floor cells are (4379, 3551) and (2435, 3586). **The same panic would stop the viz** for the author's world
with the gate on at ×1, p = 0.5.

### The diagnosis (rule 12, stage by stage)

1. **The construction lays no new land cell below the sea** (`f140_diag`: 0 cells). My first hypothesis is
   refuted.
2. **The stage** (`f140_diag2`): (4379, 3551) is 319.7 m in ON at every stage, and the breach fills it to 325.8 m
   (ON's pre-breach lake 5). In GORGE it is 178.7 m after the construction, the light pass, the droplets and the
   bathymetry, **and −8.4 m after the protected breach**. (2435, 3586) is an inland below-sea cell (−0.3 → −1.0 m)
   in both; ON's lakes cover it, GORGE's do not.
3. **The mechanism** (`f140_diag3`, body 4, which holds the floor cell):
   - Its levels: L_in 325.8, L_bed 325.7 and L_floor 168.9 m.
   - **A = 1 571 km²**, the col's drained area in the construction, where the F139 bench measured 298 km² for the
     same lake. That gives **r_lake = 1.94 at ×1**, so L(r_lake) = 178.7 m, near the floor.
   - Only 0.1 km² of the input bowl lies under that level.
   - **The pre-breach drainage holds 0.6 km² as lakes over the 52.4 km² footprint** (ON: all 52.4 km², lake 5).
   - **The breach then carves the residual pits.** Each outlet ramp descends EPS per step towards its base, and
     2 821 land cells end ≤ sea (ON: 292), 458 of them around the floor cell.
4. **A second defect, found on the way**: in the eroded GORGE world, a cell of body 4's footprint lies at
   **91.1 m, 78 m under the input floor (168.9 m)**: something lays below a drained bowl's floor. That is what
   G-drained forbids. Its stage is not located.

### Where it comes from

> **[CORRECTED in F142, after F141.]** Point 3's "where the F139 bench measured 298 km² for the same lake" and the
> deviation-1 cause below are **false**: body 4 was paired with lake 4 by index. Lake 4 is body 3 (A 300.7 / 298.1 /
> 305.2 km²). Body 4 is half of the merged below-sea lake 1000011, with A 1 571.3 km² by the inflow too. F141 names
> the causes:
> - (α) the merged lake is split into bodies;
> - (β) a below-sea lake's body is its land fringe, so "drained" means ≈ the sea;
> - (γ) the breach's ramp is anchored on each pit's floor.
>
> The text below is kept as written.


**Deviation 1** (A = the col's drained area in the construction) makes r_lake reach ~2 at ×1 for the large outlets.
- With p = 0.5, F139-R predicted no empty lake at r_world = 1. The construction drains them.
- The input footprint of lake 4 likely intercepts a large river's path, which the final lake does not. The col's
  area then counts that river.
- **The breach's behaviour is unchanged.** It is the drained bowls (sub-threshold, not held) that it now carves
  below the sea.

### Not measured (the run stopped)

G-levels, G-rim, G-ring, G-slope100, G-drained, G-area, G-tag, the non-regression (canyons, coast, θ), rule 14,
rule 18 and the cost.

**G-inert HOLDS**: the after-checks' guard reads 6 / 6, field and lakes, with the gate off (`checks_after.txt`).
Lib 593, viz 30, check clean.

## Options, named for the reviewer (none applied, by the author's choice)

1. **A = the lake's inflow** (the area of the cells entering the body plus the body). This is closer to F138 /
   F139's instrument and corrects deviation 1. It still depends on the input footprint (point 3 above).
2. **Cap r_lake at r_world** (no lake ahead of its world). This removes the drain at ×1, but keeps the per-lake
   retreat only downwards.
3. **Hold a drained bowl's residual pits as lakes**, or keep L(r) ≥ the level whose area passes the lake floor
   (5 km²). Otherwise the breach carves them. This touches the breach / lake-detection pairing, F38's territory.
4. **Locate the cell 78 m under the drained floor** before any rerun (a construction cone, or the light pass).

## Predictions

Most are unjudged (the gates were not measured).
- Mine: G-inert held.
- The reviewer's: G-inert held; the others unjudged.
- **Meta**: unjudged.

## State

**Uncommitted**, awaiting the reviewer.
- `valley_construction.rs`, `production_upscale.rs`, `upscale.rs`, `cached_product.rs`, `export/hydro.rs`: the
  gated build and its tests.
- `hd.rs`, `workspace.rs`: the plumbing and the viz.
- `tests/common/mod.rs`: `World.gorge_falls`.
- `f126_coast.rs`: `f140_g`, `f140_diag`, `f140_diag2`, `f140_diag3`.
- This folder (with `spec_gorge_age_v3.md`), ADR Finding 140.
