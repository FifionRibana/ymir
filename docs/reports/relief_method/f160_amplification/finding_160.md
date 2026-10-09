# Finding 160 — the cascade, round 2: a physics level, then Schott 2024's amplification transposed

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ». **Production is unchanged**:
the guard 6 / 6, field and lakes, was run before and after (`checks_before.txt`, `checks_after.txt`).

**The verdict in two sentences.**
- The stop rule R fires at the **first amplification level of every chain** (256² from P128, 512² from P256), on the
  crest facets. These are the physics level's facets carried by the bicubic: they are budget-independent and already
  above the témoin at the physics level itself.
- The amplification as transposed is nearly inert:
  - it carves **0.02–0.16 % of its theoretical bound**;
  - it adds **2.5–9× less detail** than the témoin in the fine bands;
  - its valley spacing stays constant in km;
  - S+ρ changes nothing (ΔR8 ±1 %);
  - a diagnostic ×100 on k does give the témoin's texture (λ, R8), but planes the land (drift −35 %).

## Partie 0

1. **F159 committed** (`a7707fb`), with the author's verdict recorded in the ADR and in `finding_159.md`.
   - `2024-MultiScaleHydro-Author.pdf` stays out of the repo: an author version with no licence statement.
   - It is cited in `docs/refs/REFERENCES.md`, now noting it is the author's copy.
2. **What the article says**, read in the copy: every claim of the reviewer's checked, the section given.
   - **§3**: T_{k+1} = D_k ∘ T_k ∘ E_k ∘ U_k (T_k) — bicubic upsampling, then fluvial erosion, thermal stabilisation,
     deposition. **Held.**
   - **§4.2**: « we ignore the uplift component responsible for the emergence of new mountain ranges ».
     - The bounded form is h ← h − k ẽ, with ẽ = min(sⁿ, s_maxⁿ) min(aᵐ, a_maxᵐ). Table 4: n = 2, m = 0.8, k = 5e-4,
       s_max = 1, a_max = 250.
     - « Bounds guarantee that N iterations … remove less bedrock than N k s_maxⁿ a_maxᵐ ». **Held.**
   - **§4.1**: the iterative drainage (eqs. 1–2), p = 1.3 « to avoid sharp fluvial incision ». « Only one iteration
     of flow routing is required » per erosion iteration. **Held.**
   - **§4.2**: k(p) = k(1 − ρ(p)). A fractal-noise ρ « reduces the axis-aligned artifacts produced by the regular grid
     discretization ». **Held.**
   - **§4.3**: a noisy critical slope (Table 4: γ ∈ [0.8, 1.4]); « the prescribed value [is not] guaranteed ». **Held.**
   - **§4.4**: a separate deposition, c = k_c e, t = Σ w g, φ = t − e, d = min(t, k_d φ), k_c = k_d = 0.1. **Held.**
   - **§5.1**: the retargeting.
     - The points are a < a₀ = 2; δ = h₀ − h_A; E is diffused over 500 iterations with the constraints held;
       h_R = h_A + E. « One retargeting step at the end of the amplification is sufficient ». **Held.**
     - **Two slips in the article**: the diffusion is printed E_{i+1} = m E₀ + (1 − m)(E_i − ΔE_i), read as a sign
       slip; and §6.6 swaps the section numbers of the retargeting and the breach.
   - **§5.2**: the multi-scale breach. Partial breaches P_k with radius r_{k+1} = 0.5 r_k, the last at r = 1. **Held.**
   - **Table 1**: erosion iterations per level, decreasing at high resolution (256²: 1 200–7 400; 4 096²: 300–800).
     **Table 2**: GPU ms per iteration, 128²–8 192² (8 192²: erosion 19.2, thermal 3.5, deposition 41.5). **Held.**
   - **Fig. 18D**: an input « tectonic simulation » [Schott 2023]. **Held.**
   - **The scale, added**: the code's presets are 20 km square (`Box2(Null, 10 000)`), with 3 000–4 000 m of relief.
     So the article's values were set for **cells of 10–80 m**, against our 0.2–3 km.
3. **The article against its published code** (`github.com/H-Schott/MultiScaleErosion`, MIT; the shaders and
   `PredefinedErosion` read).
   - **The reviewer's points, all held**:
     - the code bounds each cell at its steepest receiver: `new_height = max(new_height, receiver_height)`;
     - it has no retargeting and no multi-scale breach;
     - its talus noise is [0.9, 1.4].
   - **Found besides**:
     - **the hardness buffer is declared but never read** by `erosion.glsl`: ρ is not applied in the published
       code;
     - **the clamp is not the paper's**: `clamp(pow(slope, 2), 0, 1)`, then the PRODUCT is clamped at `max_spe =
       10 000`, with no a_max = 250 on the area;
     - **the drainage counts the cell's diagonal length** (`stream = length(cellDiag)` + inflow), metres not cells:
       about 110 per cell at 78 m cells, so a^0.8 is about 40× the paper's at equal cell counts;
     - **the deposition differs**: stream^0.3, the sediment reset outside pits, `rain = 2.6 · cellArea · 1e-5`;
     - **the talus moves `5e-5 · cell²` m per neighbour count** per iteration;
     - **the preset sequence**: 256 E 3 000 / T 600 / D 2 000; 512 E 1 500 / T 1 000 / D 700; 1 024 E 700 / T 2 000 /
       D 200; 2 048 E 400 / T 6 000 / D 150.
   - **F158 corrected in place**: the clamped form, the bound and the code's differences.
4. **F159's diagnosis recorded** (hypotheses to check, the reviewer's):
   - each level replayed an orogeny to equilibrium;
   - without an area cap the big rivers re-carved at each level;
   - the grid dependence raised the relief by about √2 per level.
   - **F160 bears on it**: with no uplift and the area capped, the drift is −0.3 to −1.4 % per level (it was +227 /
     +379 m in F159), and the trunks hold. **The first two hypotheses are consistent with F160's measurements; the
     √2 is not tested.**
5. **Checks**:
   - before: the guard 6 / 6, lib 624, viz 33;
   - after: the guard 6 / 6 identical, lib 627 (+3), viz 33 (the window's worker test replaced).

## B — what was built (gated)

**`ymir_core::cascade::physics_level`** — D1:
- F159's level from the 64² bicubic, with U₀ calibrated on the peak (2 850 m);
- the D5 lift of the interior cells ≤ sea to sea + 1 m;
- returns its `PeakCalibration`.

**`ymir_core::cascade::amplify`** — D2–D6:
- `AmpConfig::declared`, `Budget::moyen` (the code's presets);
- `AmpLevel::{upscale, erode, talus, deposit, recalage_2x2, result, substep_field}`;
- `ocean_mask`, `lift_interior`, `exact_area`, `restrict`, `retarget`;
- `Chain::{new, start_level, next_level, retargeted}`.
- All Jacobi and rayon, so it is deterministic.
- **Diagnostic A1** (not blind, see below): `AmpConfig::k_boost`, 1 by default.

**`ymir_core::cascade::hydro`** — V:
- `hydrology`: D8 `compute_flow` at sea 0 m, the Strahler order on the network ≥ max(1 km², 4 cells), the lakes
  (filled − z > 0.5 m, 8-connected, ≥ 2 cells, with area, spill and lowest cell);
- `draw`: rivers ≥ an order, the trunks ≥ 100 km² darker, the lakes' outline. Rivers are not drawn across a lake.

**Tests**:
- `the_amplification_is_bounded_and_adds_up`:
  - the erosion stays under the article's bound and only lowers; the deposition only raises;
  - the floor and the ocean hold; D5 holds;
  - the decomposition adds up;
  - the retargeting puts the ridge points on h₀;
  - **negative control**: k = 0 carves nothing.
- `the_hardness_keeps_the_mean_erodibility`.
- `strahler_and_lakes_on_a_synthetic_field`, with a **negative control**: no depression, no lake.

**The viz Cascade window, rebuilt**:
- « Calculer la physique » at 128² or 256²;
- the variant (S / S+ρ / S+U), the budget (faible / moyen / fort) and « intensité k × » (A1);
- the next level's budget per process, editable; « Niveau suivant », « Tout (→ 1024²) », « Jusqu'à 2048² »,
  « Recalage final »;
- per amplification level: Agrandissement / Érosion / Talus / Dépôt with iterations and times, and the buttons
  « + érosion / + talus / + dépôt » on top of the level's state (the later levels are discarded, and the window says
  so);
- the views Ombrage / Hypsométrie / Différence (each sub-step's own Σ, saturating at p98);
- the layers Rivières (ordre ≥ k) and Lacs;
- the hover: altitude, drainage area, Strahler, lake (area, spill).
- Its worker test drives: the physics → one level → « + érosion » → the retargeting.

**The bench** `f160_cascade.rs`: `f160_cascade`, and the diagnostic `f160_diag`.

## M — the measurements

The full printout: `f160_bench_output.txt`, run 2. Run 1 is `f160_bench_output_run1.txt`.
- **The two runs' measurements of the 11 chains to 1 024² are identical line for line.**
- Run 2 differs by the ranking correction below, which moved the 2 048² runs to P256.

### The physics level (D1)

| | P128 | P256 | témoin (same level) |
|---|---|---|---|
| U₀ (trial peak) | 6.47e-4 m/yr per km of h_iso (441 m) | 4.42e-4 (645 m) | — |
| steps (trial / run) | 71 / 62, at equilibrium | 144 / 125, at equilibrium | — |
| peak / mean / p90 (m) | 2 864 / 447 / 1 175 | 2 873 / 392 / 1 026 | 3 099 / 618 / 1 399 (128²); 3 318 / 632 / 1 436 (256²) |
| relief p50, 12.5 km blocks (m) | 608 | 756 | 876; 1 090 |
| crest facets (%) | **8.49** | **13.11** | 4.32; 5.23 |
| R8 | 0.161 | 0.018 | 0.152; 0.015 |
| lakes | 7 / 1 270 km² | 14 / 2 581 km² | — |
| D5 lifted | 331 cells (mean 31 m deep) | 748 cells (39 m) | — |
| detail per band, finest → coarsest (RMS m) | 129 / 97 / 154 / 183 | 86 / 60 / 106 / 132 / 151 | 230 / 158 / 172 / 182; 175 / 126 / 158 / 173 / 184 |

- **The peak is calibrated** (2 864 / 2 873 m), but the mean land altitude is 30–40 % below the témoin's.
- **The physics level is already faceted beyond the témoin** (8.5 against 4.3 %; 13.1 against 5.2 %).

### The amplification: S, moyen, from P128 (the other chains below)

| | 256² | 512² | 1 024² | 2 048² (from P256) |
|---|---|---|---|---|
| carved / bound | 67 km³ / 0.02 % | 56 km³ / 0.04 % | 48 km³ / 0.07 % | 26 km³ / 0.13 % |
| deepest cell / per-cell bound (m) | 40 / 9 485 | 62 / 4 506 | 67 / 2 187 | 32 / 706 |
| s ≥ s_max (% of the erosion updates) | 0.000 | 0.000 | 0.000 | 0.000 |
| drift (restricted bias) | −4.2 m (−0.94 %) | −3.1 m (−0.68 %) | −1.9 m (−0.43 %) | −1.0 m (−0.25 %) |
| trunks → previous level, p90 (cells) | 1.12 | 0.71 | 0.71 | 1.00 |
| peak / after retargeting (m) | 2 861 / 2 864 | 2 860 / 2 864 | 2 869 / 2 872 | 2 905 / 2 907 |
| walls 28° (% / cells) · témoin | 0 / 0 · 0.248 | 0.017 / 8 · 0.770 | 0.004 / 8 · 2.51 | 0.009 / 69 · 5.55 |
| crest facets (%) · témoin | **16.3** · 5.2 | **11.9** · 9.4 | 3.9 · 20.0 | 3.8 · 27.7 |
| coherent windows (%) · témoin | 38.7 · 8.4 | 75.0 · 25.1 | 90.1 · 51.5 | 93.1 · 71.5 |
| terrain R8 · témoin | 0.154 · 0.015 | 0.124 · 0.075 | 0.101 · 0.040 | 0.006 · 0.042 |
| λ (cells) · témoin | 15.1 · 7.5 | 29.5 · 10.7 | 58.4 · 16.6 | 74.5 · 26.2 |
| finest band RMS (m) · témoin | 68 · 175 | 35 · 111 | 18 · 64 | 11 · 35 |
| second band RMS (m) · témoin | 36 · 126 | 12 · 72 | 4 · 35 | 1.9 · 16.5 |
| lakes (count / km²) · témoin | 28 / 3 755 · 85 / 2 312 | 38 / 3 814 · 231 / 1 840 | 41 / 3 817 · 365 / 1 548 | 89 / 3 392 · 316 / 992 |

- **The physics level's 5 largest depressions** keep their area and spill through every level (e.g. P256's
  467 → 562 → 573 km², spill 9 → 11 → 12 m; the 280 km² one at 322 m stays at 278). **Nothing drains or fills them.**
- **λ doubles in cells at each level**, so it is constant in km: about 23 km from P128, about 15 km from P256. **The
  amplification adds no new valley order.**

### Every chain, at a glance (the first amplification level, where R fires)

| chain | level | facets % (témoin) | drift % | trunk p90 | finest band (témoin) |
|---|---|---|---|---|---|
| P128 S faible / moyen / fort | 256² | 15.4 / 16.3 / 16.5 (5.2) | −0.69 / −0.94 / −1.41 | 1.12 | 68 / 68 / 68 (175) |
| P128 S+ρ faible / moyen / fort | 256² | 15.4 / 15.9 / 16.4 (5.2) | identical to S ±0.01 | 1.12 | 68 (175) |
| P128 S+U faible / moyen / fort | 256² | 15.6 / 15.4 / 16.2 (5.2) | −0.01 / −0.02 / −0.02 | 1.50 | 76 (175) |
| P256 S / S+ρ / S+U moyen | 512² | 16.2 / 16.3 / 17.0 (9.4) | −0.46 / −0.46 / 0.00 | 0.71 / 0.71 / 1.50 | 45 / 45 / 48 (111) |

- **The budgets barely move anything**: the finest band is identical across faible, moyen and fort; the deepest cell
  goes 40 → 79 m at 256².
- **S+U holds the coarse level exactly** (drift 0.00 % by construction), but its facets reach 22.7 % at 512² and its
  peak 3 030 m before retargeting.

### S against S+ρ (D7.5)

| level | ΔR8 | Δ coherent |
|---|---|---|
| 256² | −1.1 % | 0.0 |
| 512² | +1.1 % | +0.2 |
| 1 024² | +0.9 % | 0.0 |
| 2 048² | +0.1 % | 0.0 |

- The same holds from P256 at 512² and 1 024² (+0.3 %, +1.1 %).
- **The hardness noise changes nothing, because the erosion it modulates is negligible.**

### Cost (D7.7; alert only, by the author's addendum)

**Per iteration at 2 048²** (P256 S moyen):

| process | ms per iteration | iterations | time at 2 048² |
|---|---|---|---|
| erosion | 22.1 | 400 | 8.8 s |
| talus | 7.6 | 6 000 | 45.4 s |
| deposition | 25.2 | 150 | 3.8 s |

- Measured 512² → 2 048²: **80 s**. From P128, the 128² → 2 048² chain takes 84 s.
  - The bench's cost line says « 128² → 2048² » for the P256 chain: a label slip; it measured 512² → 2 048².
- **ESTIMATE to 8 192²**, per process: erosion 97 s, **talus 907 s**, deposition 50 s. **Total 1 134 s = 18.9 min**,
  above the 15 min target and under the 1 h alert.
  - **The talus is responsible**: its declared 6 000 iterations at 4 096² and 8 192² (PROXY).
  - **What would reduce it**: fewer talus iterations at the fine levels (it moves a fixed m_L = 0.002 cell per
    iteration), or a GPU (Table 2: thermal 3.5 ms at 8 192², against our 7.6 ms extrapolated × 16).
- **Memory retained**: 222 MB at 2 048² (13 grids).

### R — the stop rule

| chain | where R fires | criterion |
|---|---|---|
| every P128 chain | 256² | facets, 15.4–16.5 % > 5.23 % (≥ 528 crests) |
| P128 S+U fort | 256² | facets |
| every P256 chain | 512² | facets, 16.2–17.0 % > 9.44 % (≥ 1 744 crests) |

The other criteria never fire on the declared design:
- walls ≤ 20 cells or under the témoin;
- drift ≤ 1.4 %;
- peak after retargeting 2 864–2 907 m;
- trunks p90 ≤ 1.50.

**The best variant** (the reading: the highest level reached without R, the physics level counting):
- **P256 S moyen**, which reaches 256²; every P128 chain reaches 128². S+ρ ties with it.
- **Correction, stated**: run 1's ranking ignored the physics level, saw every chain at « 0 » and fell back to the
  first chain (P128 S faible). So its 2 048² runs came from P128. Corrected in run 2: the 2 048² runs and the images
  come from P256. The P128 2 048² runs are kept in `f160_bench_output_run1.txt`.

### Diagnostic A1 (NOT blind; it does not change R's verdict) — `f160_diag_output.txt`, `images/f160_diag_A.png`

Written after the declared runs, to ask whether the missing lever is the intensity. Neither the budgets ×2 nor added
iterations can supply it.

| P128 S moyen, at 1 024² | k ×1 | k ×10 | k ×100 | témoin |
|---|---|---|---|---|
| drift at 256² | −0.9 % | −5.1 % | **−35 %** | — |
| peak before retargeting | 2 869 | 2 828 | **907** | 3 492 |
| mean land | 447 | 401 | 163 | 644 |
| finest / second band (m) | 18 / 4 | 19 / 5.5 | 13 / 8 | 64 / 35 |
| λ (cells) | 58 | 50 | **16.8** | 16.6 |
| R8 / coherent (%) | 0.101 / 90 | 0.075 / 87 | 0.028 / 47 | 0.040 / 52 |
| facets (%) | 3.9 | 6.5 | 18.5 | 20.0 |

- **At ×100 the texture statistics reach the témoin's** (λ, R8, coherence, facets), but the explicit erosion
  **lowers the whole land**: the peak falls to 907 m and the relief p50 to 347 m.
- **On the crop**, ×100 shows straight parallel north–south gullies: the grid's comb, under a strong explicit
  erosion.
- **The reading** (hypothesis to check): with n = 2, our 0.2–3 km cells have gentle slopes (p90 0.11–0.15, against
  the article's ~0.5–1 at 10–80 m). The term s² then decides everything.
  - At the declared k, the channels carve tens of metres.
  - At a k large enough to texture, the hillslopes themselves erode, and the land goes down.
  - **The article's form, transposed at constant n, does not separate « carve the valleys » from « lower the
    land » at our scales.**

## What the measurements say (the causes are hypotheses to check)

1. **R's verdict is not the amplification's doing.**
   - The facets that fire it are budget- and variant-independent (15.4–16.5 % across 9 chains at 256²).
   - They are already in the physics level (8.5 % at 128², 13.1 % at 256²), carried by the bicubic.
   - Hypothesis: the steady state of the physics level, a stream power with a sub-grid diffusion, is a field of
     straight flanks and sharp crests. The images show the comb of parallel straight gullies on its flanks.
2. **The amplification as transposed adds almost nothing.**
   - It carves 0.02–0.16 % of its bound and never reaches s_max.
   - The two finest bands are 2.5–9× below the témoin's; the valley spacing is constant in km; the hypsometry
     curves of every level overlap (`f160_hypsometry.png`).
   - It keeps the large forms (drift < 1.5 %, trunks hold, the depressions keep their spill), which is the article's
     promise; but it does not add the detail.
3. **The k transposition is the weak point.**
   - k_L was anchored on the land's p90 slope and on a_max in cells (the paper's Table 4).
   - The saturated channels are much flatter than that p90, and the published code's area is about 40× larger (metres
     of diagonal). So the declared erosion is two to three orders weaker than the code's at equal cells.
   - A1 shows that raising k does not fix it either: the form erodes the hillslopes and lowers the land before it
     dissects.
4. **S+ρ does nothing here** (ΔR8 ±1 %): there is no erosion to modulate. Note: the published code never reads its
   hardness.
5. **The lakes**: 3 392 km² at 2 048² against the témoin's 992 km².
   - Most are the physics level's closed depressions: the coarse field's below-sea cells, lifted flat (D5), filled to
     their spill.
   - Some outlines are straight, the coarse grid's. Nothing in the chain drains or fills them.
6. **The cost is manageable**: 80 s to 2 048²; 18.9 min ESTIMATE to 8 192², the talus first.

## The predictions

**The reviewer's** (hypotheses, now checked):

| prediction | verdict |
|---|---|
| S: the peak never exceeds the physics level's before retargeting | **refuted, barely**: S at 1 024² 2 869–2 871 against 2 864 (+7 m); from P256 2 905 against 2 873 (+32 m) |
| the drift negative and under 10 % per level | **held** (−0.25 to −1.41 %) |
| the carved volume under the bound and decreasing with resolution | **held** (P128 S moyen 67 → 56 → 48 → 41 km³ to 2 048², run 1; 0.02–0.16 % of the bound) |
| S+ρ lowers the axis alignment (R8, comb) by at least a quarter | **refuted** (±1 %) |
| the facets stay under the témoin with the noisy talus | **refuted** at the first level (256² / 512²); held from 1 024² |
| the trunks under 1.5 cells (p90) | **held** (≤ 1.50, S+U at the limit) |
| on the CPU the talus dominates the fine levels; the extrapolation exceeds 15 min | **held** (talus 45 of 58 s at 2 048²; 18.9 min) |
| meta | **held** |

**Mine** (`f160_predictions.md`):

| prediction | verdict |
|---|---|
| P-P1: U₀ 0.6–0.8 × F159's | **refuted** (0.51× at 128², 0.35× at 256²) |
| P-P1: equilibrium in < 100 and < 150 steps | **held** (62, 125) |
| P-P2: P256's mean within 15 % of P128's | **held** (−12 %) |
| P-A1: the peak never above the physics level's | **refuted** (+7 / +32 m) |
| P-A1: within ±5 % after retargeting | **held** |
| P-A2: drift negative | **held** |
| P-A2: drift of −1 to −6 % | **refuted** (smaller: −0.25 to −0.94 % at moyen) |
| P-A3: 5–40 % of the bound | **refuted** (0.02–0.16 %) |
| P-A3: halving per level | **refuted** (−16 % per level) |
| P-A4: trunk p90 ≤ 1.0 everywhere | **refuted** (1.12 at 128 → 256) |
| P-A5: walls under 20 cells or the témoin | **held** |
| P-A5: facets under the témoin at 256² / 512² | **refuted** |
| P-A5: facets above the témoin at 1 024² | **refuted** |
| P-A6: the finest bands below the témoin by 1.5–4× | **held** for the finest band (2.6–3.9×); **refuted** for the second (6–9×) |
| P-A7: S+ρ −10 to −20 % | **refuted** (±1 %) |
| P-A8: faible / fort ×0.6 / ×1.6 on the finest band | **refuted** (×1.0) |
| P-R1: the best reaches 1 024² | **refuted** (R at the first level) |
| P-C1: the talus dominates; > 15 min | **held** |
| P-V1: the guard | **held**: field and lakes = banc on all 6 states, before and after; lib 624 → 627, viz 33 → 33 (the F159 worker test replaced by F160's) |
| meta | **held** |

## I — the images (`images/`, north up, rivers of order ≥ 2 and lakes drawn)

| image | content |
|---|---|
| `f160_{A_chaine,B_plateau,C_bassin_lac4}.png` | row 1, the best (P256 S moyen): physics 256², 512², 1 024², 2 048², retargeted; row 2, the témoin 128²–2 048². F159's crops at 512² origins (128, 320), (176, 192), (74, 160) |
| `f160_best_{256,512,1024,2048}.png`, `f160_best_2048_retargeted.png`, `f160_temoin_{128…2048}.png` | whole views per level |
| `f160_hypsometry.png` | the land's hypsometric curves, one per level (red physics, orange 256², green 512², blue 1 024², violet 2 048², black retargeted). The pink lines are the R band 2 295 / 3 450 m; the grey lines 1 / 2 / 3 km. **The curves overlap** |
| `f160_S_vs_Srho_A.png` | crop A per level, row 1 S, row 2 S+ρ (from P256) |
| `f160_diag_A.png`, `f160_diag_k{1,10,100}_1024.png` | A1: crop A, rows k ×1 / ×10 / ×100 / the témoin, columns 256²–1 024² |

**The questions for the author's look**, level by level, in the images or the viz (« Cascade » → « Calculer la
physique » → « Niveau suivant »):
1. Crêtes en facettes, murs, promontoires, « pinceau » ?
2. Stries, terrasses droites, blocs ? Moins avec la dureté bruitée ?
3. Les grands fleuves et leurs affluents ont-ils une forme naturelle, et restent-ils en place ?
4. Les hauteurs te semblent-elles justes pour une Corse ?
5. Assez de détail, ou faut-il plus d'itérations à ce niveau ? (The measurements say: not enough, and the iterations
   do not supply it; « intensité k × » in the viz lets you see what A1 measured.)

## Open, for the author's decision (nothing chosen)

- **(a) The physics level is the faceted part.**
  - Its steady state at 128² / 256² carries the facets and the comb.
  - Options: a hillslope regime resolved at that level (F159's open (a)), or a physics level at 64² only with the
    amplification doing more.
- **(b) The erosion form at our scales.**
  - n = 2 with gentle slopes either carves nothing or lowers everything (A1).
  - Options: n = 1 (m/n = 0.5, like our physics); the code's area units (metres) rather than the paper's Table 4; or
    normalising the slope per level.
- **(c) The lakes** from the coarse field's below-sea cells (D5): lift them to their spill only where the témoin has a
  basin, or let the amplification breach them (the multi-scale breach, not built).
- **(d) The cost**: the talus's iteration count at the fine levels.

## At commit
- **The author**: « Avec la physique à 256², on a des formes plus ramifiées dans les reliefs. » The author ran the
  window to 4 096² in 203 s (S, physics 128²) and in 417.8 s (S+U, physics 256²).
- **The stop rule's reference was badly chosen (the reviewer, recorded as given).**
  - The témoin is the production world the author rejected at F157.
  - It grows more faceted with the resolution itself: 5 % at 256², 28 % at 2 048².
  - So the rule fired on the physics level's facets before the amplification could act.
  - **It is replaced by a real reference, the relief of Corsica (F161).**
- **The track left untested**: the block retargeting was only tried with uplift (S+U), never alone.
- **F161 follows**: Corsica's DEM as the texture reference; an amplification that carves (n = 1, k calibrated per
  level on Corsica, the area capped) and a smooth retargeting that separates « creuser les vallées » from « abaisser
  le continent »; the physics level at 256² by default.
