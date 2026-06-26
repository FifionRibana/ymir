# #165 — C1 climate design: temperature + precipitation by moisture transport

Anchor-before-build (the create-analogue of measure-before-fix). The climate module
(`climate/{temperature,precipitation,biomes}.rs`) is empty (M3 stubs). This establishes
the physical anchors + the closure formulation BEFORE coding. Relief (#155) is finished
→ its orography (in metres via `c1_altitude_norm_to_metres`, slopes via `c1_km_per_cell`)
is the input. Climate DERIVES from relief + latitude → re-runnable when relief changes.

## 0. The latitude finding (reframes the wind question)

`C1_DOMAIN_KM = 1024` ÷ `111 km/°` ≈ **9.2° of latitude**. Wind belts are ~30° wide
(Hadley 0-30°, Ferrel 30-60°, polar 60-90°). So **the domain sits within a SINGLE belt**
— it does not span equator-to-pole. Therefore:
- **Latitude is a placement PARAMETER** (the continent's centre latitude on the globe),
  NOT a span across all belts.
- The chosen band fixes ONE prevailing wind direction (uniform across the domain) plus
  the local sea-level temperature + its mild ~7°C across-domain gradient.
- The DOMINANT temperature variation is ALTITUDINAL: a 4.65 km orogen × 6.5°C/km ≈ 30°C,
  vs ~7°C latitudinal across 9°. So orography drives the climate (as intended).

→ `c1_climate(relief, latitude_deg)`; `latitude_deg` is the one product-placement input
(the future viz slider controls it). **Decision surfaced below.**

## 1. Anchor table — every term on a real quantity, NO free knob

| term | physical anchor | value / formula |
|------|-----------------|-----------------|
| sea-level T(lat) | equator-pole gradient | ~+27°C @0° → ~−25°C @90°, ≈ a cosine: `T_sea = T_eq − (T_eq−T_pole)·sin²(lat)` |
| lapse rate | environmental adiabatic | **6.5 °C/km** · altitude_m |
| air moisture capacity | **Clausius-Clapeyron** | `e_sat(T) = 6.112·exp(17.62·T/(243.12+T))` hPa → saturation mixing ratio |
| prevailing wind | Hadley/Ferrel/polar | band(lat): trade **easterlies** 0-30°, **westerlies** 30-60°, polar **easterlies** 60-90° |
| ocean evaporation | bulk/SST | `E ∝ e_sat(SST)` (warm sea charges the air toward its capacity) |
| orographic precip | **Smith & Barstad 2004** | `P ∝ wind·∇h_↑` (uplift × moisture); fallout of condensed moisture |
| rain shadow | conservation | precipitated moisture REMOVED from the flux → dry leeward |

## 2. Temperature (computed first — precipitation needs it)

`T(i,j) = T_sea(lat(j)) − 6.5·altitude_km(i,j)`, where `altitude_km =
max(0, c1_altitude_norm_to_metres(h, ss))/1000` (sea/ocean at sea-level T). `lat(j)`
spans the ~9° band around `latitude_deg`. Output: `temperature: GridF32` (°C).

## 3. Precipitation (the conservative moisture-transport closure)

A **1D conservative transport along the fixed wind lines** — complete (rain shadow,
windward soak) yet light (one O(N) pass, straight lines since the wind is uniform in
the single belt). NOT isotropic orographic smear, NOT a 2D GCM.

**Algorithm** (wind direction `ŵ` from the belt; streamlines = rows/cols along `ŵ`):
```
for each streamline from the UPWIND edge, marching downwind cell by cell:
    M = moisture flux (carried)
    cap = saturation mixing ratio at the local T   # Clausius-Clapeyron
    if cell is OCEAN:  M += k_evap · e_sat(SST)               # pick up moisture
    ascent = max(0, ∇h·ŵ)            # along-wind UPHILL slope, m/km (the contracts)
    if ascent > 0:                                            # windward
        P = k_oro · M · ascent       # Smith-Barstad uplift fallout
        P = min(P, M)                # can't precipitate more than carried
        precip[cell] += P;  M -= P                            # CONSERVE
    if M > cap:  precip[cell] += (M-cap);  M = cap            # supersaturation falls out
    # descending (leeward): no uplift term → M unchanged → rain shadow downwind
```
Output: `precipitation: GridF32` (mm/yr after a unit calibration anchored to a wet-belt
reference, e.g. windward-tropical ~2000 mm/yr).

## 4. Orographic input (physical units — the contracts pay)

Along-wind slope `∇h·ŵ` computed on the metre heightmap: `Δaltitude_m / Δdistance_m`,
with `Δaltitude_m = c1_altitude_norm_to_metres(Δh)` and `Δdistance_m =
c1_km_per_cell(grid)·1000`. Windward (ascent>0) → precip; leeward (ascent<0) → shadow.
The finished #155 relief (high cordilleras) is what casts the shadows.

## 5. Re-calculability (the #165 founding principle)

`c1_climate(relief, latitude_deg) → (temperature, precipitation)` is a pure DERIVED
computation: no frozen field, no per-step state. Re-run it whenever the relief changes
(a future closure, or a different seed) → the climate follows. The viz latitude slider
controls the `latitude_deg` PARAMETER (and re-derives), it does not replace the field.
Same discipline as the coordinate contracts: one anchored computation, not a tuned output.

## 6. Scope of the maillon
`c1_climate(heightmap, ss, latitude_deg, &ClimateConfig) -> ClimateResult{temperature,
precipitation}`, in `climate/` (un-stub the modules). Gated: a new entry, no existing
path changed (byte-identical for everything else). Validation (planned): VISUAL —
leeward deserts (rain shadow behind the cordilleras), wet windward flanks, the
latitudinal T belt; + metric (precip conserved along streamlines; T lapse holds;
windward/leeward precip ratio). Biomes (Whittaker T×P) = the NEXT maillon, not this one.

## 7. Latitude decision — RESOLVED: default 45° (westerlies), product-anchored

Default `latitude_deg = 45` → **westerlies** (wind from the west). NOT arbitrary: Living
Landz §2.2 targets a dense TEMPERATE playable region (cités, villages, river valleys,
4 seasons) — 45° westerlies *is* that register, so the default = the product target; the
central case (temperate continent) is right without touching the slider.

**Belt-selection logic in place from the start** (cheap, one function) so the slider is
functional with no later rework:
```
wind_dir(lat) = | lat ∈ [0,30)   → trade EASTERLIES  (wind from E → streamlines E→W)
                | lat ∈ [30,60)  → WESTERLIES        (wind from W → streamlines W→E)
                | lat ∈ [60,90]  → polar EASTERLIES  (wind from E → streamlines E→W)
```
At any placement the ~9° domain sees ONE belt → ONE uniform prevailing direction →
parallel streamlines → the simplest conservative-transport case (and realistic). For the
45° default: moisture enters from the **western** ocean, precipitates on the windward
(west-facing) relief, and the flux dries **eastward** → eastern rain shadows. The slider
varies `latitude_deg` → the band → the direction flips; no structural change.

## 8. Next step
Implement `c1_climate(heightmap, ss, latitude_deg, &ClimateConfig)` per §1-§6 (un-stub
`climate/{temperature,precipitation}.rs`), gated; then the planned visual + metric
validation (§6). Biomes (Whittaker) follow as a separate maillon.
