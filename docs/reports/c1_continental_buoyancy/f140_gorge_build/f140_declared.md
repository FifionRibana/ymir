# F140 — the gates' instruments, declared BEFORE the measurements (2026-10-02)

## The worlds
- **The témoin** is the viz state C2 /10 col, at the canonical framing, with run_hd's tail. Field = the
  conditioned one (metres); the construction = `carve(S1)` with the world's ValleyConstruction.
- **OFF**: `new(k, Some(0.1))`. **ON**: + `lake_base: InputLakesAndBasins`.
- **GORGE(age, p)**: ON + `gorge_retreat: GorgeRetreat::v3(r(age), p)`, with k = F121_AGE_K × age and map (a)
  r(age):

  | age | ×0.7 | ×0.85 | ×1 | ×1.2 | ×1.4 |
  |---|---|---|---|---|---|
  | r | 0 | 0.5 | 1 | 1.5 | 2 |

- **Measured**: GORGE at the 5 ages with p = 0.5, and at ×1 with p = 0 and 0.25. OFF and ON at ×1.

## Per GORGE world
- **The bodies**: the construction's own diagnostics (`Skeleton::gorge_bodies`, `gorge_body_of`): L_in, L_bed,
  L_floor, A, r_lake, L(r_lake), φ, the slope, the falls.
- **"An input-lake body"**: a body whose cells are under the input pre-drainage's lake map (the 14 D8 lakes'
  family); otherwise a depression body.
- **G-levels**: for each body with an outflow:
  - the final lake = the run_hd-tail lake with the largest overlap of the body's cells;
  - |final level − L(r_lake)| ≤ 1 m;
  - a body with r_lake = 2 is "drained" if no final lake covers ≥ 50 % of it.
  - Counts, and the worst mismatches listed.
- **G-rim** (r_lake ≤ 1): cut = L(r_lake) − the lowest cell of the body's outer ring, on the construction and on the
  eroded world. Gate: no cut > 10 m.
  - The negative control is ON's cut against L_in (F135-K1's instrument, at ×1).
  - **Deviation from the brief's letter**: the cut is read against L(r_lake), not L_in, because the spec cuts the
    col down to L(r_lake) on purpose.
- **G-ring**: the ring cells the clamp bound are those where the construction's z = the clamp ± 0.05 m and S1 is
  above it.
  - For each, the slope to its outside neighbours (not in the body, not on the ring) on the construction.
  - Gate: none > 28°, except at the col (the tagged lip). Failures are located, not corrected.
- **G-slope100**: along each final lake's D8 path from its true col (the run_hd tail), on the conditioned field:
  - the max slope over 100 m over the first 10 km;
  - the first 2 cells are excluded when that body's head fall is tagged.
  - Gate: ≤ 28°. The negative control is ON at ×1.
- **G-drained** (r_lake = 2): the cells whose skeleton D8 path enters the body (its upstream). Gate: none on the
  construction below L_floor − 1 m.
- **G-area**: the total area of the final non-crater lakes, and of the input-lake bodies' final lakes, by age at
  p = 0.5. The latter is compared with F139-R's 1 308.5 / 1 202.7 / 521.3 / 103.0 / 19.0 km². Gate:
  non-increasing.
- **G-tag**:
  - the tagged falls (`World::gorge_falls`): counts and heights by kind and age;
  - the gorge's own cells (the construction path from each col, over the cells the gorge raised): no 2-cell drop
    over H_f = max(1.5·S·100 m, 10 m) beyond the lip.
- **Cost**: skeleton + carve, timed, GORGE against ON at the same age.

## At ×1, p = 0.5 (and OFF, ON)
- **Canyons**: `f95_criteria` + `over_dug_depression` (viz dcfg). **The coast**: spurs near a coastal wall.
- **θ on the carved links with the construction's mask**: F135-M2, the based set excluded. Gate: within OFF's CI
  0.493 [0.492, 0.495].
- **Rule 14**: the removal S1 − the eroded world on land (km³): total, inside the bodies, elsewhere.
- **Rule 18**: OFF / ON / GORGE side by side.

## G-inert
- **The guard** after the change: 6 / 6 field and lakes (`f123_viz_guard`).
- `rivers.json` without falls is byte-identical: a permanent test (`rivers_json_carries_the_tagged_falls_only_when_…`).
- The eroded cache writes no `gorge.json` when off.
- The gate alone (without the extended lake base) is inert: a permanent test.
