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

## Discipline notes
- Volet 1 = reading; the scale is what the model encodes (two-scale, phantom peak) — not
  tuned toward a target height. Pinning the scale does not create relief.
- Volet 2 = diagnostic; the critical-wedge lift is a model change, scoped separately.
- The metre figures use the downstream `(v−0.5)·2·4000` convention **only to quantify the
  gap** — they inherit the Volet-1 ambiguity and are not a pinned contract.
