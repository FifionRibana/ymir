# F171 — declared before any measurement (2026-10-11)

**A rock that depends on depth: a soft cover on a hard basement (K1), then the heterogeneous K (K2).**

**Production is unchanged**: everything goes through the cascade and the « continent » profile; the guard 6 / 6 is not touched.

**Frozen**: T = 5.62 Myr, Airy, the sources' drive, F164's cascade, the torus, Davis-Suppe, the ages A.
- **K_couv** is the cascade's current K (F159's erodible-rock declaration).
- **ρ_K, d_c, σ_lnK and L** are declared here; their sensitivity is measured, never chosen.

**Not blind**: F170's readings and the code. No F171 variant has been run.

**Part 0**: F170 committed (`f21c5c2`). Its ADR carries the torus default and its cost, the 16-plate selection, the reviewer's
error on the sutures, the peneplain and the reviewer's reading, and the coast and confluence notes. The author's two answers
(« oui et oui ») are added under it.

## A — the plates by default

### A1 — the profile

- **The default**: 16 plates, F170's selection unchanged (the BFS and the compact families; k ∈ [⌈0.35 N⌉, ⌊0.75 N⌋]; 50–62 %
  land; circumnavigable with r* ≥ 1; the largest r*; ties as at F170).
- **The deterministic fallback**: if no candidate qualifies, 20 plates (k 7–15), then 24 (k 9–18).
- The retained N is reported per seed. `profile::select_continent_default`.

### A2 — viz

- A slider « Plaques » (8–24) in the Cascade window, used by « sélection ».
- At 16 (the default) it falls back to 20 then 24. At any other value it selects at that N only and shows whether a candidate
  qualifies.

### A3 — the negative control

8 plates without the selection is F168 bit for bit: the existing test.

## B — cover and basement (K1, the physics level)

### B1 — the law

- K(x, t) = K_couv · f(x, t), with f = 1 while E_cum(x, t) < d_c(x), and f = 1 / ρ_K beyond.
- **E_cum** = the cumulated lowering by the erosion alone, in metres: per step, max(0, −(Δz_incision + Δz_talus + Δz_diffusion)).
  - The uplift and the rebound are excluded.
  - It is cumulated in the physics level's frame, which is the rock's frame (the sources' drive does not move the land).
- **d_c = 0** on the cratonic shields at t₀ (C1's `cratonic_mask` after `craton_shield_fraction`, which is the shields only),
  upsampled by nearest cell.
- **K is evaluated at the start of each step.**
- **The step that crosses the contact**: the fraction-of-step rule, linear in K at n = 1.
  - A cell whose step lowering e exceeds its remaining cover r = d_c − E_cum keeps r at K_couv.
  - The excess (e − r) is divided by ρ_K: z ← z + (e − r)(1 − 1/ρ_K).
  - No smoothing. The donors see the uncorrected receiver for that step only.
  - At ρ_K = 1 the correction adds exactly 0.

### B2 — the values

| constant | central | sensitivity | source |
|---|---|---|---|
| ρ_K | **10** | 3, 30 | Stock & Montgomery 1999, p. 4983 (abstract) and Table 2, p. 4986. At m = 0.4: granite and metamorphic 10⁻⁷–10⁻⁶; volcaniclastic 10⁻⁵–10⁻⁴; mudstones 10⁻⁴–10⁻² (their units; only ratios are used). Hard over soft is 10 to 10⁵; 10 is the cautious end, as the brief asks. |
| d_c | **1.5 km** | 0.75, 3 km | **A design value**: the order of magnitude of platform covers. No source in `docs/refs`; not checked against a paper. |

The sensitivity runs on the témoin, old profile only.

### B3 — the amplification levels

They keep their Corsica-calibrated k, with no cover or basement (the limit of this variant). Only K2's multiplier reaches them.

### B4 — the negative control (permanent)

ρ_K = 1 gives F170 bit for bit (`physics_history_rock` with ρ_K = 1 against `physics_history`).

## C — the heterogeneous K (K2 = K1 × M)

### C1 — the field

- ln M(x) = σ_lnK · Σ_j a_j g_j(x).
  - g_j is a unit-variance band noise at wavelength λ_j = L / 2^j: `fbm_torus` with one octave, standardised over the grid.
  - a_j = 2^(−j), normalised so Σ a_j² = 1 (an fBm with H = 1).
  - So the mean of ln M is 0 and its standard deviation is σ_lnK.
- **L = 25 km**, j = 0 … 6 (λ = 25, 12.5, 6.25, 3.1, 1.6, 0.78, 0.39 km), down to 2 cells of 2 048². A design value; L is
  F162's hardness wavelength.
- **It is fixed to the rock.**
  - Octave j is born at the first level whose 2 cells ≤ λ_j (256²: j ≤ 3; 512²: j = 4; 1 024²: j = 5; 2 048²: j = 6), at the
    terrain's coordinates.
  - It is then carried to the finer levels by the same upsample and **the same warp W** as the terrain.

### C2 — where it applies

- The physics level's K (K1 × M).
- Every amplification level's k (k_level × M at the level), through the N1 solver's K field.

### C3 — the values

- **σ_lnK = 0.5**: within-lithology variability, a design value read from Stock & Montgomery 1999, Table 2, p. 4986.
  - The same lithology differs by ×1.9 (Kauai basalts, 3.8–7.3 × 10⁻⁶), ×2.1 (Tumbarumba against Tumut granitoids and
    metasediments) and ×10 (Wheeo against Lachlan) between catchments.
  - σ = 0.5 puts ±2σ at ×7.4.
- **Sensitivity**: 0.25 and 1.0, on the témoin, old profile.
- **No anisotropy.**

### C4 — the negative control (permanent)

σ_lnK = 0 gives K1 bit for bit.

## D — the measures

- **The grid**: the témoin (with the fallback), 42, 1, 9.
  - K0: F170's old profile re-run here as the coast reference and for the images; the other ages are read from F170, and run
    for the témoin.
  - K1 and K2 at B0, 20, 100 and 300 Myr. Physics 256², the frozen chain to 2 048², the torus.
- **The selection**: N, k, family, land, r*.
- **Macro**: the classes against Europe, the land deciles against Europe's, the peaks, the mountain share against A per variant.
- **The rock at the surface**, per age: the land shares of cover, basement and shield, and the altitude of the plains standing
  on the basement.
- **The lithological knickpoints (the detector, declared)**: on the 256² physics level, the Strahler ≥ 3 network
  (`hydro::hydrology`).
  - For each river cell whose D8 receiver lies on the other rock (cover ↔ basement, from the final E_cum), take the mean
    channel slope over the 3 cells upstream (S_u) and the 3 cells downstream (S_d) along the flow.
  - It is a knickpoint if the basement side is steeper by a factor ≥ 2.
  - Its height is the drop over the basement side's 3 cells.
  - Reported: the count and the median height.
- **Lakes.**
- **Texture**: Corsica's battery (whole and mountain class), F by class, the slope guard.
- **Regularity** for K1 and K2 at the same ages: σ_θ, CV_L, CV_λ, R2_g at 391 m, Hovius.
- **The coast on the torus** against Corsica. « Degraded » is measured against F170 (K0) on the same profile.
- **Images** (north up, checked, the offset at the end):
  - per seed, shaded and classes: K0 | K1 | K2 on the old profile, then K1 at the 4 ages;
  - the rock map;
  - a river profile crossing a contact;
  - 25 km crops: Corsica | K1 | K2.
- **The cost**: CPU s per world, idle.

## R — the stop rules (the brief's)

- **The fallback fails for the témoin**: report it and go on with the 3 others; the « ≥ 3 » rules become « ≥ 2 of 3 ».
- **Macro success (old, K1 or K2)**, on ≥ 3 of 4:
  - mountains 16.5–37.2 %;
  - plain + plateau 37.5–84.5 %;
  - median land 98–390 m;
  - circumnavigable with r* ≥ 1;
  - no coast degraded > 10 % against F170 on the same profile.
  - Then F172.
- **A persistent peneplain** (median < 50 m on ≥ 3, old, K1): report it. The next candidate is n > 1 or an incision threshold
  (Harel et al. 2016). Nothing is implemented.
- **Peaks off the scale** (> 7 km, or the slope guard failing on ≥ 3) in B0 or young: report it. T is not recalibrated.
- **Forbidden**: as F165–F170.
