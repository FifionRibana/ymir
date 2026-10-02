# Specification v3 — the age sets how far the gorge has retreated (Finding 140, BUILT GATED, OFF by default)

**v3 = v2** (`f138_spec_v2/spec_gorge_age_v2.md`) **with the author's decisions of 2026-10-02 and the reviewer's two
corrections.** Everything in v2 not restated here stands. Provenance tags as in v2:
- **MEASURED**: this dossier;
- **ANCHORED**: a source;
- **PROXY**: a labelled stand-in;
- **DECISION**: the author's, or proposed to the author.

## 0. The decisions (recorded as given, ADR Finding 139)

- **When room is lacking, the gorge steepens** (no large fall). The shortage fall exists only beyond the 28° cap
  (none in the témoin, F139-C3).
- **The retreat per lake**: « Ça peut être un slider ou sélection (0, 0.25 ou 0.5) car là je ne sais pas. Par
  défaut j'irai sur 0.5. » So p ∈ {0, 0.25, 0.5}, **default 0.5** (DECISION).
- **The head falls: as proposed** (one lake in three without a fall).
- **Proposed by the reviewer, not contested by the author**:
  - the map (a): ×0.7 → r 0, ×1 → r 1, ×1.4 → r 2, piecewise linear (DECISION);
  - **m = 10** (DECISION);
  - the below-sea lakes included with L_in = the merged spill (F139-C5), except 1000001 and the 3 lakes without a
    spillway segment;
  - **H_f's floor of 10 m** (DECISION).

## 1. The parameters

- **The gate** `ValleyConstruction::gorge_retreat: Option<GorgeRetreat>`, `None` by default (skipped in serde, so
  the guarded digests are unchanged).
- **`GorgeRetreat { r_world, p, m, a_ref_km2 }`**:
  - r_world: from the age selector by map (a);
  - p ∈ {0, 0.25, 0.5};
  - m = 10;
  - **A_ref = 418.7 km², a DECLARED CONSTANT** (MEASURED, F139-R: the median outlet area of the témoin's 14 D8
    lakes). It is not the current world's median: a lake's retreat must not depend on the other lakes.
- **The retreat**: `r_lake = min(2, r_world · (A / A_ref)^p)`, with A = the lake's outlet area.
- **The level**: `L(r) = L_in − r·(L_in − L_bed)` for r ≤ 1, `L_bed − (r − 1)·(L_bed − L_floor)` beyond (DECISION,
  linear).
- **The slope**: `S = min(tan 28°, max(m · S_loi(A), S_req))`, with `S_loi(A) = k·A^−0.5`.
- **φ**: 0 with probability 1/3, otherwise U[0.1, 0.5], by splitmix64.
- **The tag**: `H_f = max(1.5 · S · 100 m, 10 m)`.

## 2. What the construction does (per present lake ≥ 1 km²)

**The present lakes are the construction's own**: the input lakes (the pre-drainage's lake map) and the
`basin_base` closed-depression components (the merged spill: F139-C5's L_in), those ≥ 1 km².
1. **L_in**: an input lake's surface, or a depression's spill.
2. **L_floor**: the lowest cell of the construction's input.
3. **L_bed**: the minimax height between the lowest cell and the outflow cell (the level the lake stands at once
   its col is cut).
4. **The outflow**: the body cell with the largest drained area whose receiver leaves the body; **the col** is
   that receiver; **A** = the col's drained area.
5. **The lake is a base at L(r_lake)**: the χ walk's lake stops use L(r_lake) instead of L_in.
   - Rivers in are based on the lowered lake.
   - At r_lake = 2 (drained) the upstream is based on L_floor, never on χ from the sea.
   - **This IS the drained-lake rule** (v2 § 2): "L_floor plus the gorge's floor at the bowl's exit" reduces to
     L_floor, because the gorge's floor at the exit is L(2) = L_floor.
6. **The gorge**: along the outlet path (from the col down the skeleton's D8, until the sea, another present lake,
   or the law below):
   - the head fall φ·D_g at the lip (between the col and the next cell);
   - then `z_gorge(s) = L' − S·s`, with `L' = L − φ·D_g` and S steepened as needed.
7. **The invariant** acts on the samples' floor, before the cones: `zf = max(base + k·χ, z_gorge(s))`.
   - Written as `base ← max(base, z_gorge − k·χ)` on the path cells, after the χ walk, so `line_samples`
     (`valley_construction.rs`, `let zf = base[c] + vc.age_k * chi[c]`), `traced_polylines` and `floor_m` read it
     unchanged.
   - At the col, `z_gorge = L(r_lake)`: the col is cut to the lake's level.
8. **The rim**: every cell of the body's outer ring keeps `≥ L(r_lake)`, a clamp in `carve_diag`'s final `v`.
9. **The falls**: computed by the construction (`Skeleton::gorge_falls`), carried in the eroded product (a
   `gorge.json` sidecar, written only when non-empty), exported in `rivers.json` per segment as
   `features: [{kind, point, x, y, height_m, lake_id}]`, tagged above H_f. Off: `rivers.json` is byte-identical.

## 3. Deviations from v2's letter, declared before building (each with its reason)

1. **A is the construction's drained area at the col** (the skeleton's D8), not the world's accumulation + inflow
   (F138 / F139's benches): the construction does not know the final world.
2. **L_bed is the construction's minimax on its input**, not "today's level" (which F136-K7 measured as the bed
   sill). It is the same object, computed upstream.
3. **φ is drawn from splitmix64 of the body's lowest cell index**, not of (seed, lake id). The construction has no
   seed and no final lake id. It stays deterministic per world.
4. **The 3 below-sea lakes without a spillway segment cannot be recognised in the construction.** Their depression
   components are treated like every other. 1000001 is not an input body and is untouched.
5. **A lake whose outflow is not on a trunk (A < 10 km²) gets no gorge sample.** Its col is not cut, so it keeps
   L_in at any r. G-levels counts these.

## 4. The gates (each with its control, rule 15)

- **G-inert**: the gate off → the field and lake hashes bit-identical; `rivers.json` byte-identical.
- **G-levels**: each present lake's final level = L(r_lake) ± 1 m. The control is the construction's own L(r_lake)
  (diagnostics), computed on the input.
- **G-rim** (r_lake ≤ 1): no cut > 10 m below L(r_lake) on the input ring. The negative control is today's ON
  (18 cuts > 10 m against L_in, F135-K1).
- **G-ring** (the v2 review's point 3): the ring cells the clamp raised, and their slope to the outside neighbours.
  - **Gate**: no slope > 28° from a ring cell to an outside neighbour, except at the tagged lip. If it fails,
    locate, do not correct.
- **G-slope100**: on every untagged descent, the max slope over 100 m ≤ 28°. The negative control is today's ON
  (2 steps > 40°, F138-T2).
- **G-drained**: no cell upstream of a drained lake under L(2) beyond the law from that base.
- **G-area**: the total lake area is non-increasing in age. The control is F139-R (p = 0.5: 1 308.5 → 1 202.7 →
  521.3 → 103.0 → 19.0 km²).
- **G-tag**: no gorge cell over H_f; each tagged fall on a ≤ 2-cell break. The counts and heights are given per
  age.
- **Non-regression**: canyons 0; the coast unchanged; θ on the carved links (the construction's mask) within OFF's
  CI (0.493 [0.492, 0.495]).
- **Rule 14** (the removal, km³, where) and **rule 18** (OFF / ON / gorge on one world).
- **Cost**: the time added per world, measured, against the < 2 s estimate.

## 4b. The viz

- A toggle "Recul de la gorge (âge)" and a selector p (0 / 0.25 / 0.5, default 0.5).
- r is read from the age selector by map (a).
- The tagged falls are markers on the Relief layers, with their height, type and lake on hover.
- The badge reads "non gardé" when the gate is on (no bench reference).
