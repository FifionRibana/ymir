# F160 — declared before any measurement (2026-10-09)

**The second round of the cascade.** It runs one physics level, then Schott 2024's amplification transposed: bicubic
upsampling, then bounded erosion, a noisy talus and deposition, **with no uplift**. A final retargeting follows; rivers
and lakes go into the Cascade window, and the chain goes to 1 024² (2 048² for the declared runs).

**Production is unchanged**: the chain is gated, and the guard 6 / 6 runs before and after.

Written after reading:
- Schott 2024, §3–§6, Tables 1, 2, 4;
- its published code (`github.com/H-Schott/MultiScaleErosion`, MIT): the shaders and the preset sequence in `main.cpp`;
- F159's measurements.

The author's cost addendum (2026-10-09) is integrated below, **before** any measurement. Nothing of the amplification
has been run.

## D1 — the physics level (P128 and P256, both measured)

**F159's level, unchanged in its processes**:
- uplift U = U₀ · max(0, h_iso) / 1 km;
- the existing implicit stream power (relief-v3 at the cell, K 2e-3, dt 1e5 yr);
- the talus, then the linear diffusion at D = max(D_phys, 0.16 K cell²);
- run to equilibrium (mean |Δz| < 1 % of mean U·dt for 5 steps); caps 400 (128²) and 300 (256²);
- on a fixed grid: P128 from the 64² bicubic ×2, P256 from the 64² bicubic ×4.

**U₀ is calibrated once per physics level on the peak**:
- **the target**: the highest land cell = **2 850 m**, the middle of the author's 2 700–3 000 m (DECISION of the
  author, « type Corse »);
- **the method**: a trial run at U₀ = 1e-4, then U₀ ← U₀ · 2 850 / peak_trial (n = 1 makes the steady state linear in U;
  the talus's activity is reported, as in F159).

**What it must provide**: the large forms the article takes as input (its Fig. 18D: an input « from a tectonic
simulation »).

## D2 — the amplification (S), per level up to 1 024², and 2 048² for the declared runs

**The order**: bicubic ×2 (F159's periodic Catmull-Rom, z and the sea) → **erosion** (N_E iterations) → **talus**
(N_T) → **deposition** (N_D). This is the article's A_k = D_k ∘ T_k ∘ E_k ∘ U_k. **No uplift.**

**The working field**:
- metres;
- the ocean is held at 0 m during a level (the base level) and its bathymetry restored at the end;
- the land floor is +0.5 m, production's base level.

**The routing** (§4.1):
- the iterative multiple-direction drainage with exponent p = 1.3 (eqs. 1–2);
- **a counted in cells** (each cell contributes 1);
- one routing iteration per erosion iteration and per deposition iteration;
- **one deviation**: at a level's start, a is initialised by one exact accumulation (cells sorted by height, the same
  weights) on the upscaled field. The code restarts its drainage at zero at each `Init`, so the first ~250 iterations
  would only rebuild it.

**Erosion** (§4.2), **explicit and Jacobi** (double buffer) like the article:

  h ← max(h − k_L · (1 − ρ) · ẽ, h_r, 0.5 m), ẽ = min(s², s_max²) · min(a, a_max)^0.8

- s is the steepest physical slope (m/m) and h_r that receiver's height. The receiver bound is **the code's, not the
  paper's**. It is kept because an explicit step can otherwise overshoot below its receiver and dig a pit.
- **Why not our implicit solver**:
  - the article's control is a bounded carving per level, which the clamps give; the closed-form n = 1 update cannot
    carry n = 2 and min(a, a_max);
  - the implicit solver re-floods (priority flood) at every step, O(n log n), while the iterative drainage is O(n) and
    parallel.
- **n = 2, m = 0.8** (Table 4).
- **s_max = 1** (45°, a physical angle; transposed as is). At our cells it binds rarely; the share of cells where it
  binds is reported.
- **a_max = 250 cells.**
  - The cap is the same number of cells at every level: 610 km² at 256², 153 at 512², 38 at 1 024², 9.5 at 2 048².
  - So each level differentiates only the channels of its own scale. The larger rivers, carved by the coarser levels,
    erode at the capped rate and are not re-deepened preferentially.
  - That is the cause the reviewer named for F159 (hypothesis to check).
- **k_L, from a declared carving depth**:
  - the target: d_L of a saturated channel (a ≥ a_max) at the reference slope s_ref,L, at the moyen budget;
  - **d_L = 200 m at 256², halving per level**: 100 at 512², 50 at 1 024², 25 at 2 048².
  - The halving is the self-affine bound H = 1 (a feature half the size gets half the depth). **PROXY.**
  - s_ref,L is the **p90 of the land slope** of the upscaled field at the level's start, measured at run time.
  - So k_L = d_L / (N_E^moyen · s_ref,L² · 250^0.8).
- **The article's bound per cell**: N_E · k_L · s_max² · a_max^0.8. Reported per level, with the volume it bounds.

**Talus** (§4.3, the code's form), Jacobi: h ← h + m_L · (α − β).
- α and β count the neighbours above and below that are steeper than s₀(p) (8 neighbours, physical distance).
- **s₀(p) = tan 33° · (0.8 + 0.6 · ν(p))**, ν ∈ [0, 1] a 2D fractal noise:
  - 3 octaves;
  - base wavelength 8.7 km, the code's 1 / 0.0023 m⁻¹ = 435 m as a fraction of its 20 km domain, carried over to our
    400 km;
  - source: `WorldSeed` "cascade_talus".
  - So s₀ ∈ [0.52, 0.91], i.e. 27.5°–42.3°.
- **What is transposed, and why**:
  - the RELATIVE range [0.8, 1.4] (the paper's) is applied to production's physical repose tan 33°;
  - the absolute [0.8, 1.4] (38.7°–54.5°) was set for the article's 10–80 m cells. At our 0.2–3 km cells, where a cell's
    slope averages much steeper metre-scale slopes, it would never act.
- **m_L = 0.002 · cell_L** (metres per neighbour count per iteration). The code moves 5e-5 · cell² m, i.e. 0.30 m at its
  78 m cells, which is 0.0039 · cell. Halved here for stability at our lower slopes. **PROXY.**

**Deposition** (§4.4), the paper's form in metres:
- e_m = k_L (1 − ρ) ẽ, the erosion rate in m;
- c = k_c e_m;
- t(p) = Σ_{q∈V⁺(p)} w(q, p) g(q);
- φ = t − e_m;
- d = min(t, k_d φ) if φ > 0, else 0;
- g ← c + t − d; h ← h + d;
- k_c = k_d = 0.1 (Table 4). g starts at 0 at each block (the paper's g₀ = 0). Sediment reaching the sea is lost.
- The code's form (exponent 0.3, sediment reset outside pits) is not used.

## D3 — the budgets

**The moyen budget = the code's preset per resolution** (`PredefinedErosion`):

| level | E | T | D | Table 1, erosion only (its five terrains) |
|---|---|---|---|---|
| 256² | 3 000 | 600 | 2 000 | 1 200–7 400 |
| 512² | 1 500 | 1 000 | 700 | 1 200–4 800 |
| 1 024² | 700 | 2 000 | 200 | 600–2 700 |
| 2 048² | 400 | 6 000 | 150 | 500–2 000 |

- **Faible** = ×0.5 and **fort** = ×2, for the three processes.
- k_L is fixed by the moyen budget, so the budgets scale the carving depth (and the talus and deposition work).
- The viz can add iterations per process and per level, « sur la base de l'existant ».

## D4 — the final retargeting (§5.1), once, at the end of each chain

- **Constraint points C**: land cells whose exact drainage (p = 1.3, in cells) is < **a₀ = 2** at the final level.
- **The error**: δ = h_phys↑ − h_A at C, with h_phys↑ the physics level bicubic-upscaled to the final level.
- **The diffusion**: E is diffused over 500 Jacobi iterations, E ← m·E₀ + (1 − m)·mean₄(E).
  - The paper prints « E − ΔE », read as a sign slip for a diffusion. Declared.
- **The result**: h_R = h_A + E on land, with the land floor at 0.5 m; the sea unchanged.

## D5 — the coarse field's below-sea cells

**At every level** (the physics level's start and each upscale), a cell ≤ 0 m not 4-connected to the map border is
**land**, a land depression, not open sea. It is lifted:
- to +1 m at the physics level (so the implicit solver's scalar sea test sees land; the depression stays closed and the
  priority flood routes it);
- to +0.5 m in the amplification.

The count and their mean original depth are reported.

## D6 — the variants

- **(S)**: Schott transposed, ρ = 0.
- **(S+ρ)**: ρ(p) = 0.5 + 0.4 · f(p), with f ∈ [−1, 1] a fractal noise:
  - 6 octaves, base wavelength 25 km down to ~0.4 km;
  - source: `WorldSeed` "cascade_hardness";
  - k is divided by (1 − 0.5), so the mean erodibility equals S's.
- **(S+U)**, the fallback, Ymir's own and absent from the article: S plus
  - a small uplift per erosion iteration, 10 % of the level's carving depth over the block (∝ max(0, h_phys↑) / peak);
  - a 2×2 retargeting at the end of each level: z ← z + U₂(z_prev − R(z)) on land (R the full-weighting restriction,
    U₂ the bicubic), so the level restricted equals the previous level.
- **The multi-scale breach is not built this round.** It is named for the lakes' resumption.

**The runs**:
- P128 × {S, S+ρ, S+U} × {faible, moyen, fort}, up to 1 024²;
- P256 × {S, S+ρ} × moyen, up to 1 024²;
- **2 048²** for {S, S+ρ, S+U} × moyen, from the physics level that reaches 1 024² best.

## D7 — the instruments (per level, per variant and budget; the témoin's block means at the same level)

**F159's**:
- walls 28° and the 33° band;
- crest cells and their −∇²z;
- crest facets, with the synthetic control;
- coherent windows (the comb stand-in);
- terrain R8 (windows of 8 cells, noise floor printed);
- λ;
- relief p50 (12.5 km blocks);
- mean and p90 land;
- trunks ≥ 100 km² (D8 `compute_flow`): n+1 → n, and against the témoin;
- the D6 restriction.

**Added**:
1. **Detail per scale band** (the spectrum's stand-in, no FFT): bands of [2^k, 2^(k+1)] cells, each the RMS over land of
   B_(2^k)(z) − B_(2^(k+1))(z), with B_r a periodic box blur of side r and B₁ = z. Reported in km, against the témoin
   at the same level.
2. **The drift**: the restricted bias and |Δ| (m), and the bias / level n's mean land altitude.
3. **The peak, the mean and the p90 of the land**, against the target, before and after the retargeting.
4. **The carved volume per level** (km³, Σ of the erosion over land × the cell area), against the bound
   N_E · k_L · s_max² · a_max^0.8 × the land area. The same for the deepest cell.
5. **The axis alignment with and without ρ**: R8 and the coherent windows, S against S+ρ at the same budget.
6. **The physics level's depressions**:
   - the 5 largest (filled − z > 0.5 m, 8-connected) at the physics level;
   - at each level, the depression containing the same point (the physics depression's lowest cell, mapped), with its
     area (km²) and spill (m);
   - the total count and area.
7. **The cost**:
   - wall ms per iteration and per process at each resolution;
   - the **per-process extrapolation to 8 192²**, scaled by the cells, with 4 096² / 8 192² budgets E 300 / 200,
     T 6 000 / 6 000, D 100 / 100 (the code's trend; PROXY), plus the total.
8. **The share of erosion cells** where s ≥ s_max (the slope bound binding).

## D8 — the reference

The production témoin (C2 /10 col, PSEED, CANONICAL_ORIGIN, 8 192²), brought to each level 128²–2 048² by centred
block means (F159's D9).

## R — the stop rule (the brief's, verbatim, with the author's cost addendum)

« Pour la meilleure variante, on s'arrête au niveau où l'une de ces conditions est vraie :
- murs à 28° ou facettes de crête au-delà de la référence au même niveau, **et** plus de 20 cellules ;
- dérive du relief d'un niveau à l'autre de plus de 10 % (biais restreint) ;
- point culminant hors de la cible de plus de 15 % après recalage ;
- troncs de plus de 100 km² : p90 de déplacement au-delà de 1,5 cellule ;
- coût extrapolé à 8192² au-delà de 1 h. »

**The author's addendum (2026-10-09)**: the cost is **no longer a stop condition** for the levels up to 2 048². They are
measured to the end, and their result judges the method. The 8 192² extrapolation is reported **per process** as an
alert, with the process responsible and what would reduce it.

**How it reads, declared**:
- **« la meilleure variante »**: the physics level × variant × budget that reaches the highest level without firing R.
  Ties go to the lowest terrain R8 at that level.
- **walls / facets**: the cascade's share above the témoin's at the same level **and** more than 20 cells (wall cells,
  or facet crest cells). The facets count only where their control passes.
- **the drift**: |restricted bias| > 10 % of level n's mean land altitude.
- **the peak after retargeting**: outside [0.85 × 2 700, 1.15 × 3 000] = [2 295, 3 450] m.
  - At an intermediate level it is that level retargeted on the physics level, for the evaluation only.
  - The chain itself retargets once, at its end (D4). This reading was added before any measurement.
- **the trunks**: the p90 of the n+1 → n distance, in cells of level n.

« Pas assez de détail » is not a stop: it is reported per band.

## V — the Cascade window

At each computed level:
- **rivers**: D8 accumulation ≥ max(1 km², 4 cells), Strahler order on that network, filter « ordre ≥ k » (default 2,
  the rivers LL filter the author likes), drawn on the shade. The trunks ≥ 100 km² are highlighted.
- **lakes**: filled − z > 0.5 m, 8-connected, ≥ 2 cells, with outline and area.
- **hover**: altitude, drainage area (km²), Strahler order, lake (area and spill) or not.
- **per-process budgets**, editable per level, and one button per process (érosion, talus, dépôt) that adds iterations
  to the selected level on top of its current state. The levels after it are discarded, which is said.
- **layers can be switched off.**
- **The production world is not touched**: the guard 6 / 6, and a worker test.

## I — the images

**For the best variant and for the témoin**:
- F159's three crops: A (128, 320), B (176, 192), C (74, 160), 512² origins;
- a whole view per level;
- north up, with rivers and lakes drawn.

**Also**:
- a hypsometry per level against the target;
- S against S+ρ side by side on crop A.
