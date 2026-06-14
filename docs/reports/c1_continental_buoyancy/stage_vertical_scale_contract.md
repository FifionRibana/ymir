# #155 — Vertical-scale contract + high-mountain ceiling (diagnostic)

Read-only diagnostic. Two **distinct** problems, deliberately not conflated:
1. **CONTRACT** — what is `1.0` in metres, what does `[0,1]` cover? (Definition by
   reading the model — no parameter change.)
2. **CEILING** — orogens render far below the 6000–8000 m (Himalaya) target. This is
   a MODEL limit, not a scale-definition gap; pinning the scale does not create relief.

Probe: `c1_closure_morphology::probe_vertical_scale_ceiling` (+ `probe_orogen_equilibrium`),
`#[ignore]`, seeds {42, 1988, 2026}, A′ defaults, `c1_hd_production`.

---

## Volet 1 — the vertical contract (reading, not a decision)

**There is NO single `norm→metres` mapping. The altitude is double-normalised with an
internally incoherent metre scale — this IS the ~×2 ambiguity** flagged in
`project_c1_vertical_scale_contract`.

The conversion chain (`c1_production_altitude` → `upscale_from_c1`):

1. **`compute_isostasy_inner`** produces an **asymmetric** `[0,1]` heightmap.
   - sea level sits at `sea_norm = max_depth_m/(max_depth_m+max_elevation_m)
     = 500/4500 ≈ **0.111**` (NOT 0.5).
   - LAND ramp `[h_sea, h_max] → [0.111, 1.0]`; the metre anchor is
     `IsostasyResult.peak_altitude_m = (h_max−h_sea)/(h_max−h_min)·max_elevation_m`
     — **field-relative** (depends on the seed's raw S̃ spread), NOT a fixed constant.
     So `1.0` ≠ a fixed number of metres on land.
   - OCEAN ramp `[h_min, h_sea] → [0, 0.111]`.

2. **`apply_stein_stein_bathymetry`** then **overwrites** oceanic cells with
   `−depth_m / depth_scale_m` (`depth_scale = 5000`) → a **second, FIXED** metre scale
   (`1 unit = 5000 m`), with **sea at 0** (not 0.111). The two halves now use different
   metres-per-unit AND different sea-level zeros — incoherent by construction (the S-S
   module doc explicitly notes it breaks the isostasy `[0,1]` contract; downstream
   consumers were gradient-only so it never mattered — until a metre scale is needed).

3. **`upscale_from_c1`** re-normalises `(v + 1.13)/2.26`, clamp `[0,1]`:
   - deepest S-S ocean `−1.13 → 0.0`; S-S sea `0 → 0.5`; isostasy-continental sea
     `0.111 → 0.549`; land top `1.0 → 0.943`.
   - **Two "sea levels" now disagree** (0.5 vs 0.549), and land tops at 0.943 not 1.0.

→ **Honest contract:** land metres = field-relative `peak_altitude_m`; ocean metres =
fixed 5000 m; spliced at incompatible zeros, then re-offset. Pinning a single scale is a
**model decision** (unify it), not a reading. Until then: **report relief in normalized
units / ratios, not metres.** The downstream guess `(v−0.5)·2·4000` (used in
`probe_craton_calibration`) is exactly *one arbitrary pick* of the ambiguous scale.

**Extra incoherence found by measurement:** `peak_altitude_m ≈ 3300 m` (stable across
seeds) is a **PHANTOM** — it is computed from the global raw `h_max`, which is the
**oceanic** advective pile-up cell (S̃ ≈ 2.15 at the rigid margin, #145). That cell is
then *overwritten by Stein-Stein into deep ocean*. So `peak_altitude_m` describes a cell
that becomes ocean — it does **not** state any actual land elevation. The real rendered
land max is ~half that (below).

---

## Volet 2 — the high-mountain gap, quantified + attributed (diagnostic, not a fix)

| seed | orogen margin S̃ p95 | cells ≥ 1.9 (EH region) | model `peak_altitude_m` (phantom) | **HD rendered land-max** (guess conv.) |
|------|--------------------:|------------------------:|----------------------------------:|---------------------------------------:|
| 42   | 1.302               | 0 / 1042                | 3333 m                            | **1746 m**                              |
| 1988 | 1.466               | 0 / 1160                | 3325 m                            | **2075 m**                              |
| 2026 | 1.351               | 0 / 1285                | 3274 m                            | **1864 m**                              |

**The gap:** rendered orogens top out at **~1750–2075 m** vs the **6000–8000 m** target
→ a factor of **~3–4×** short.

### Attribution (chain-of-innocence) — what binds, in order

1. **NOT the EH ceiling.** Orogens equilibrate at S̃ ~1.3–1.5; **0 of ~1100 margin cells
   reach 1.9**. The EH `h_eq = 2.0` headroom is **entirely unused**. → **Raising `h_eq`
   does nothing** (orogens never reach the current ceiling). *(This corrects the earlier
   "orogens EH-ceiling-bound at S̃~2.0" verdict — that S̃~2.0 figure was the OCEANIC
   margin spike, S-S-overwritten, not the continental orogen.)*

2. **The DS-deposition vs erosion balance is the dominant cap.** `probe_orogen_equilibrium`
   (erosion-OFF counterfactual): margin p95 rises only **1.30→1.59 / 1.47→1.73 /
   1.35→1.59**. So in-loop erosion brakes ~0.25, BUT even erosion-OFF orogens stall at
   **~1.6–1.7** — below the DS margin-peaked target (≈1.95 @ d=1) and far below `h_max=2.5`.
   The deposition rate + advective smearing of the margin band is the **primary** S̃ limiter.

3. **`max_elevation_m = 4000` is a hard codomain ceiling.** Even a hypothetically-maxed
   land cell cannot render above 4000 m. **6000–8000 m is outside the conversion codomain
   entirely** — independent of any S̃ gain.

### Ceiling chantier — scoped (for separate decision, NOT coded)

To reach the Himalaya, in causal order:
- **(a) PRIMARY — lift the orogen S̃ equilibrium from ~1.4 toward/past 2.0.** This is the
  **critical-wedge over-thickening** subject (more convergent-margin deposition / less
  advective smearing / less in-loop erosion at the wedge). *This is the real binding
  constraint — not the EH ceiling.* Watch [[feedback_recursive_tuning_signals_structural]].
- **(b) SECONDARY — raise the EH ceiling `h_eq`** at convergent margins, so the lifted S̃
  is not re-capped at 2.0. Moot until (a) lifts S̃ near 2.0.
- **(c) raise `max_elevation_m`** (and unify the scale per Volet 1) so the conversion admits
  >4000 m absolute. Note: raising `max_elevation_m` alone scales *everything* (plains too)
  → it adds absolute height but no prominence; it only matters once (a) widens the S̃ spread.

**Reframe of the deferred Phase-3 "EH/prominence" subject:** it is NOT "raise the EH
ceiling" — it is "**lift the orogen S̃ equilibrium (critical-wedge)**," with EH-raise +
`max_elevation_m`-raise + scale-unification as downstream follow-ons. Same bind as the
relative-prominence verdict, now measured in absolute height: the orogens are simply LOW
in S̃, not capped by EH.

---

---

## Maillon 2 — sea-level unification + the norm→m contract (FIX)

The diagnostic above left TWO repairs. Maillon 1 (land-ceiling, merged) fixed the
phantom peak. **Maillon 2 unifies the scale end-to-end and exposes the norm→m
function** — the prerequisite for rivers/biome/climate.

### W7 verdict (read before coding)
- **500 vs 5000 is NOT a conflict — two physical regimes.** `max_depth_m=500` =
  submerged **continental shelf** (~200-500 m, physical); Stein-Stein
  `depth_scale_m=5000` = **deep oceanic** lithosphere (~5000 m, physical). The S-S
  param doc states 5000 was chosen "consistent with the isostatic convention" — the
  intent was always a shared sea-zero. Unify on a common **sea=0, uniform
  metres/5000 nondim**; 500 stays as the shelf anchor for continental-submerged cells.
- **The 0.549 origin.** Isostasy outputs `[0,1]` with continental sea at
  `sea_norm≈0.111`; S-S writes oceanic sea-centred (sea=0). The fixed upscale
  `(v+1.13)/2.26` maps S-S-sea(0)→0.5 but isostasy-sea(0.111)→0.549. Two stages, two
  sea anchors.
- **The re-offset is NOT a bug.** `ALTITUDE_NORM_HALF_RANGE=1.13 ≈ 5651/5000` (S-S
  deepest ocean) + sea→0.5 is a deliberate **resolution-invariance** constant (#151).
  It stays fixed; the other stages align to it.

### The fix
`c1_production_altitude` now **sea-centres** the continental altitude to
metres/`depth_scale_m` (sea→0), matching the S-S oceanic convention. Both consumers
(viz render + `upscale_from_c1`) keep their fixed `(alt+1.13)/2.26`, which now maps
**sea→0.5 exactly** (0.549 gone) and the FBM coast logic (sea_level=0.5) finally
aligns with the real coastline. v2/export untouched (they use `compute_isostasy`
`[0,1]` + their own `upscale_with_fbm` sea_level — C1-specific path).

### The unified vertical scale (the contract)
`c1_altitude_norm_to_metres` (+ inverse `c1_metres_to_altitude_norm`):

```text
    metres = (norm − 0.5) · 2 · ALTITUDE_NORM_HALF_RANGE · depth_scale_m
           = (norm − 0.5) · 11300       [defaults 1.13, 5000]
```

| norm | metres | what |
|-----:|-------:|------|
| 1.0  | +5650  | scale ceiling (unreachable: `max_elevation_m=4000` caps land at norm ~0.85) |
| ~0.73| ~2300-2600 | **measured highest mountains** (seeds 42/1988/2026) |
| 0.5  | 0      | **sea level** (exact, all seeds) |
| 0.0  | −5650  | S-S ocean asymptote (~5651 m) |

**Single linear scale, no piecewise seam** — land AND ocean on one inverse-able
relation; rivers/biome/climate read metres without knowing the regime. Land occupies
`[0.5, ~0.73]` (≤0.85 maxed); ocean `[0, 0.5]`. The headroom `~0.73→1.0` is RESERVED
for the separate critical-wedge high-mountain chantier — the scale tells the truth,
and the gap it reveals is the diagnostic of the relief not yet produced.

### Acceptance (probe `probe_sea_unification_acceptance`)
1. Continental coast → norm **0.5000 exactly** (was 0.549), all 3 seeds. ✓
2. norm→m exposed + round-trips; anchors as table above. The craton/orogen heights
   (the Jordan caveat + the high-mountain gap) are now **defined metres**
   (2265-2596 m highest). ✓
3. Resolution invariance preserved BY CONSTRUCTION — sea-centring is on the 64²
   coarse, before upscale (the fixed 1.13/sea→0.5 untouched). ✓
4. Render re-judged (2-3 seeds): coast at 0.5 aligned (shelf band renders), land
   legible, no re-conversion artifact → improvement (coherent + aligned coast).
5. v2/export byte-identical (C1-only path); isostasy 13/13, lib 446/0, all C1 green,
   only pre-existing v2 `rectangular_simulation` red.

### Frontier (NOT done here — distinct chantiers)
- **Submarine relief**: unifying the scale ≠ creating ocean-floor structure. Oceans
  stay flat (IQR≈0) at the −5650 m scale; ridges/trenches = a separate chantier.
- **High mountains** (the norm `0.73→1.0` headroom): the critical-wedge orogen lift —
  a relief chantier, independent of the scale.

---

## Maillon 3 — the HORIZONTAL scale (coordinate contract's other half)

Surfaced while scoping drainage: km² river thresholds need a **km/cell**, and the
C1 product had **no pinned horizontal scale** — a ~5× ambiguity, the exact pendant
of the vertical gap:
- Legacy `GenerationConfig`: `continent_size_km=300`, `meters_per_pixel=40`.
- C1 TDD §11: domain `1×1` = "~1000-5000 km regional", dx ≈ "~2 km at 512²".

These disagree, and the C1 figure is itself a 5× range → km² thresholds would rest
on a disguised extent knob. Pinned as its own mini-maillon (the coordinate contract
is a shared foundation — biomes/climate/villages, not just rivers — like the
vertical, done as its own chantier, not slipped into the first consumer).

**The pin (anchored, revisable):** `C1_DOMAIN_KM = 1024` — the TDD §11 lower anchor
`2.0 km/cell × 512`, making its implicit "~2 km at 512²" explicit. dx = 16 km @64²,
0.5 km @2048². Consistent with §2.2's dense playable-region gameplay (not a stretched
whole continent). A "whole continent" product intent → set ~3000 km; **one constant,
every consumer follows** (nothing else encodes the horizontal scale).

Exposed: `c1_km_per_cell(grid)` = `C1_DOMAIN_KM/grid` (resolution-independent — domain
km fixed); `c1_cell_area_km2(grid)` — the unit for resolution-independent
drainage-area thresholds (`accumulation × cell_area_km2` = upstream km², invariant for
a fixed physical catchment). Pure additive functions → byte-identical (no existing path
changed). Unit test asserts the §11 anchor (2.0 km/cell @512²) + resolution-independence.

**The coordinate contract is now complete**: vertical (`c1_altitude_norm_to_metres`,
sea=0.5→0 m, ±5650 m) + horizontal (`C1_DOMAIN_KM=1024`, `c1_km_per_cell`). metadata.json
§9.3 (meters-per-pixel beside elevation) = these two together. Drainage (the first
consumer) can now anchor thresholds in km².

---

## Discipline notes
- Volet 1 = reading; the scale is what the model encodes (two-scale, phantom peak) — not
  tuned toward a target height. Pinning the scale does not create relief.
- Volet 2 = diagnostic; the critical-wedge lift is a model change, scoped separately.
- The metre figures use the downstream `(v−0.5)·2·4000` convention **only to quantify the
  gap** — they inherit the Volet-1 ambiguity and are not a pinned contract.
