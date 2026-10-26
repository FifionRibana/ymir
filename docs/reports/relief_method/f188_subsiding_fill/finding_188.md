# Finding 188 — the fill of the sinking zones: built behind its flag, with its controls, but the declared stop fires at A3 (between 80 and 20 kyr the plains of the fill's zone move 6.2 points, against F174's 1); B built and controlled, not measured; seed 18's « rectangular » lake is cut by the torus' edge

**Status: STOP at A3. Parts A and B built with their controls; Part 0.2 answered; Part 0.1, the tiers and the plates not run
(the rule: « arrête-toi et rapporte, sans élargir »). Nothing committed.**
- Declared: `f188_declared.md` (one amendment before any measure). Predictions: `f188_predictions.md`.
- Output: `f188_a3_output.txt`.

## The author's words (F187's plates), verbatim

« Seed 19: Ca donne un aspect wetland, des marécages ou des archpels cotiers. Seed 18: je n'ai pas vu ce lac. Seed 8: Oui
effectivement c'est un peu gros. Les haut fonds cotiers sont très bien c'est quelque chose qui me manquait. Pourquoi on avait
retiré le remblaiement des rivières? Autre point, on ne voit pas les rivières justement. »

- **The answer to his question**: the fluvial deposition went through K_eff. K_eff was removed because the cascade's K is an
  effective K fitted on natural rivers, whose slope already holds the deposition's effect (Yuan JGR eq. 25), so K_eff counted it
  twice (F187). Part A brings the fluvial deposition back only where the ground sinks.

## Part 0

- **0.1, the rivers on the plates**: built (`f188_part0`), not run because of the stop.
  - The threshold is the project's, `hydro::river_threshold_km2` (`hydro.rs:39`): max(1 km², 4 cells), declared at F162–F163.
    The viz's Cascade window draws with it (`ymir-viz/src/ui/cascade.rs:1478`).
  - The trunks are ≥ 100 km² (`planform::TRUNK_KM2`).
  - The straight reaches are read by `a_dir` (I1, F163).
- **0.2, seed 18's « rectangular » lake**:
  - the lake (1 761 cells on F187's map, `ll_lakes_full_seed18.png`) touches both the right edge (x = 255) and the left edge
    (x = 0);
  - **it is cut by the torus' edge** and continues on the left side;
  - its straight side is the image's frame; no quantity read by a parent cell is involved.

## Part A — the fill of the sinking zones (`Deposition::subsiding_fill`)

- **A1**: on a land cell whose non-erosive motion (the step's tectonic source + the load's flexure) is negative, the fluvial
  pass deposits after Yuan (JGR eq. 13–15), G = 1.55 (Guerit 2019, p. 855); elsewhere G = 0.
  - **Declared**: the load's flexure is the previous step's, since the step's own is computed after the pass.
  - **How**: the pass takes a per-cell G (`deposition::pass_routed_g`).
- **A2**: no iteration; the deposit enters the cover, and the existing controls count it.
- **A4's controls**:
  - **the flags off = F187 bit for bit**: seed 23's field and cover against F187's run saved before F188's code
    (`fam3_f_no_keff.bin`);
  - a world that never sinks deposits nothing (`the_f188_controls`);
  - the mass < 10⁻⁹ on a sinking zone with the load, the live marine step, A and B.
- **The control's first setup was wrong, and was corrected**: a step in the uplift (the land rising, the sea strip fixed) rings
  in the bicubic resampling and sinks the coast. The control now uses a uniform uplift.

## A3 — the convergence test: FAILS (the declared stop)

- **The setup**: seed 8, intermediate, F187 + A, the fluvial step bounded at 80 kyr against 20 kyr. The zone is the cells with a
  fluvial deposit at 20 kyr: 17 023 cells, 41 560 km².

| on the zone | 80 kyr | 20 kyr | the change | F174's criterion |
|---|---|---|---|---|
| density | 0.2721 | 0.2789 | 2.4 % | ≤ 5 % ✓ |
| λ | 20.65 km | 21.05 km | 1.9 % | ≤ 5 % ✓ |
| S₁₀ | 1.1502 | 1.1374 | 1.1 % | ≤ 5 % ✓ |
| plains / plateau / hills / mountain | 36.5 / 3.8 / 35.2 / 24.5 % | 42.7 / 1.5 / 33.3 / 22.5 % | **the largest move 6.2 points (the plains)** | < 1 point ✗ |

- The runs: 101 s at 80 kyr, 225 s at 20 kyr. Both mass controls ≤ 5 × 10⁻¹¹.
- **The explicit fluvial deposition is not converged in its sub-steps**, as F184 and F185 found for the whole pass. Restricting
  it to the sinking land does not change that.
- **An observation, not a correction**: the zone covers a quarter of the map. The load's flexure in Δz_v holds the torus'
  uniform part (the mean load's subsidence at wavenumber 0, F187), which sinks most of the land at once. So A acts far beyond the
  coastal lowlands it targets.
  - **Corrected (F188's closure, D2)**: this observation was wrong. The uniform mode is an uplift of +390.6 m (D2.1), not a
    subsidence. The broad zone comes from Δz_v being negative on some steps of cells that rise in total (D2.4).
- **My error in the bench, corrected before the report**: the first reading shifted the classes (PLAIN is 1, not 0). The
  recomputation from the saved fields gives the table above. The verdict is the same.

## Part B — the sediment's erodibility (`Deposition::soft_cover`)

- **Built**:
  - the rock law's factor is 1 (K_couv) on a cell with a cover at the step's start (`cascade.rs:684`);
  - the crossing step's correction counts the remaining cover as max(the rock's, H);
  - a per-cell record of the marine walk's stranding, `SedRecord::marine_stranded_m`, for B2.
- **Controls**:
  - the flags off = F187 bit for bit (above);
  - nothing deposited = no change;
  - the mass < 10⁻⁹.
- **Not measured** (the stop).

## The predictions

| | the claim | the result |
|---|---|---|
| Z3 | A3 holds | **refuted** (6.2 points) |
| W3 (mine) | A3 fails (F184–F185's non-convergence) | **held** |
| Z1, Z2, Z4, Z5, W1, W2, W4, W5, W6 | — | not measured (the stop) |
| W7 | the mass holds for A, B, A + B | held in the controls; not measured on the panel |

## For the author

1. **A as declared does not converge** between 80 and 20 kyr on its zone (+6.2 points of plains).
   - Two separate questions follow:
     - whether to bound the fluvial step at 20 kyr or below for the fill alone (a cost: 225 s against 101 s at the physics
       level, which F189's budget does not hold);
     - whether Δz_v should exclude the torus' uniform flexure, so the fill acts on the coastal lowlands only.
   - Both are the author's.
2. **B and Part 0.1 (the rivers on the plates) do not depend on A.** They are built and can run without A:
   - B on the panel, ≈ 16 min;
   - Part 0.1 on F187's plates, ≈ 16 min.
   - They wait for the go-ahead, since the rule stops the finding.

## The checks

- The flags off = F187 bit for bit on seed 23.
- `the_f188_controls`, `the_f187_controls`, `the_f186_controls` and `the_accommodation_controls` pass.
- The full checks (lib, viz, guard, production's hash) are run before any commit.

## After the stop: Part 0.1, B, A′ and its A3

**After the stop: Part 0.1, B on the panel, A′ and its A3** (on mains; outputs `f188_ref_output.txt`, `f188_part0_output.txt`,
`f188_tier1_b_output.txt`, `f188_a3p_output.txt`):

- **The flags off = F187 bit for bit**, rechecked after A′'s code (seed 23, against F187's saved run).
- **0.1, the rivers on the 2 048² plates** (F187's configuration; the project's threshold max(1 km², 4 cells), the trunks
  ≥ 100 km² in dark blue):
  - plates `images/p0_rivers2048_off_A_seed{8,18,19,23}.png` (without | with the deposition) and `images/p0_lakes_seed*.png`;
  - on the filled zones (the fine land under a covered 256² cell), without → with the deposition:

  | seed | the filled zone | S₁₀ | `a_dir` (the straight, D8-aligned share; isotropic 0.222) | the order-≥ 2 pieces |
  |---|---|---|---|---|
  | 8 | 17 860 km² | 1.294 → 1.189 | 0.437 → 0.571 | 3 163 → 2 754 |
  | 18 | 23 012 km² | 1.515 → 1.161 | 0.232 → 0.662 | 4 592 → 3 316 |
  | 19 | 6 524 km² | 1.508 → 1.259 | 0.269 → 0.394 | 3 119 → 539 |
  | 23 | 38 970 km² | 1.377 → 1.200 | 0.375 → 0.537 | 6 131 → 5 624 |

  - **Z5 held**: on the filled zones the rivers are straighter (S₁₀ −0.11 to −0.35) and bound to the D8 axes (`a_dir` up to
    0.66). The filled surfaces are the bilinear 256² surface (`bury_smooth`), and the D8 network follows its gradient.
- **B on the panel** (the world without the deposition / F187 / F187 + B):
  - the mass ≤ 5.4 × 10⁻¹¹; the circumnavigation lost on no seed; (v) 5 / 5; the physics level's wall median 97.2 s
    (F187 93.0 s);
  - **seed 8's east lake** (the lakes in the box x 150–239, y 50–159): 0 / 3 174 / **498 km²**;
  - the large lakes, F187 → B:
    - 8: 7 739 → 2 393 km²;
    - 18: 6 289 → 710 km²;
    - 19: 332 → 349 km²;
    - 23: 2 876 → 5 706 km²;
    - 29: 4 568 → 4 785 km²;
  - the unconnected depressions, F187 → B: 1 → 3, 0 → 4, 30 → 48, 5 → 10, 3 → 2;
  - the plains move −0.8 to +2.8 points; the median moves −44 to +7 m.
  - **Seed 23's piles**:
    - **the stranding is 83–99 % of their marine deposit in F187's run** (23007: 2 327 of 2 787 m; 26553: 2 080 of 2 093 m; …);
    - so the walk's stranding is the source (`deposition.rs:199-211`), and my F187 hypothesis (the rock law) was wrong;
    - **with B they are gone** (z −413 to −22 m, H ≤ 4 m, nothing stranded there);
    - **but the marine deposit above 100 m does not fall on seed 23** (4.44 → 4.38 % of the marine deposit): other piles
      remain. Z4 is refuted as written (by half), though the 7 cells it named are cleared.
- **A′ (`Deposition::accommodation_fill`)**:
  - its controls hold: a world that never sinks deposits nothing fluvial; the mass < 10⁻⁹ on a sinking zone; the flags off =
    F187 bit for bit;
  - **its A3 FAILS, the declared stop**: seed 8, F187 + A′, 80 against 20 kyr, on the zone with a fluvial deposit (43 765 km²):
    - density 0.2516 → 0.2783 (**9.6 %**);
    - λ 4.2 %; S₁₀ 0.9 %;
    - the classes 38.4 / 4.4 / 29.4 / 27.8 → 39.7 / 2.4 / 33.5 / 24.4 % (**4.1 points**, the hills);
    - the runs 130 and 286 s; the mass ≤ 4.3 × 10⁻¹¹.
  - **Not run, as the rule asks**: the panel with A′ and A′ + B, and the 2 048² plates with A′.
- **A remark on the test itself, for the author** (not run): A3 changes the step of the whole fluvial stage (the incision, the
  lakes, the routing), not only the deposit.
  - **F187 without any fluvial deposit has never been tested between 80 and 20 kyr.**
  - If it also fails, the criterion measures the splitting of the incision itself, and neither A nor A′ is what fails. That
    control costs ≈ 6 min.
- **Muto & Steel 1997** (A′'s principle) is not in `docs/refs`: the PDF is to be deposited.

## The rivers' plates, corrected

**The rivers' plates, corrected** (`f188_part0b_output.txt`; plates `images/p0b_rivers2048_off_A_seed{8,18,19,23}.png`, without |
with the deposition):
- **My error in the first overlay**: it drew the network on any land fine cell of a 2 × 2 pixel, and in the lakes. The corrected
  overlay draws only where the map's own sampled cell is land above the sea, and outside the lakes.
- **On the filled zones without the lakes**, without → with the deposition:

| seed | S₁₀ | `a_dir` (isotropic 0.222) | the order-≥ 2 pieces |
|---|---|---|---|
| 8 | 1.296 → 1.126 | 0.253 → 0.609 | 1 203 → 1 593 |
| 18 | 1.382 → 1.129 | 0.216 → 0.731 | 2 387 → 2 433 |
| 19 | 1.560 → 1.217 | 0.198 → 0.406 | 1 666 → 420 |
| 23 | 1.364 → 1.136 | 0.259 → 0.593 | 2 978 → 3 908 |

  - **Z5 holds more strongly once the lakes are out.** Without the deposition the filled zones' rivers are isotropic
    (`a_dir` 0.20–0.26); with it they run straight along the D8 axes (0.41–0.73).
- **What the plates still show (faithful to the fields)**:
  - **the networks « in the sea » of the world without the deposition** are F187's islet mosaic. The amplified basement pokes
    above 0 m on every other pixel of the shallow seas, and the network drains those islets.
  - **the solid patches and 45° lines with the deposition** are perfectly flat surfaces: basins filled to their spill at 256²
    and rendered as the plane 256² surface by `bury_smooth`. They are no longer lakes for the hydrology, and the flats'
    routing (Garbrecht–Martz) draws solid areas and diagonals there.
- **For the author**: the straight rivers on the filled plains (Z5) and the flat filled surfaces are F189's « rivers on the
  filled plains » question, which the reviewer conditioned on Z5.

## A3's control

**A3's control, the result** (seed 8, 80 against 20 kyr; `f188_control_output.txt`):

| 80 kyr relative to 20 kyr (80 − 20) / 20 | density | λ | S₁₀ | plains / plateau / hills / mountain | F174 |
|---|---|---|---|---|---|
| the control (F187, no fluvial deposit), A′'s zone | −0.8 % | −4.6 % | +0.0 % | +1.81 / +2.91 / −6.20 / +1.49 points | fails |
| the control, all the land | +0.7 % | +3.5 % | −0.5 % | +0.49 / +1.59 / −4.20 / +2.13 points | fails |
| A′, its zone | −9.6 % | −4.2 % | +0.9 % | −1.32 / +1.97 / −4.09 / +3.44 points | fails |
| **A′ minus the control** | **−8.8 %** | +0.5 % | +0.9 % | **−3.13** / −0.93 / **+2.11** / **+1.95** points | **fails** |

- **The control fails too**: F187 without any fluvial deposit is not converged between 80 and 20 kyr (the hills −6.2 points on
  A′'s zone, −4.2 on all the land).
  - So F174's criterion, at this change of step, measures the splitting of the incision itself.
  - F185's and F188's « non-convergence » verdicts mixed it with the deposit's.
- **A′ judged relatively, by the rule written before the computation, fails too**: its excess over the control is −8.8 % in
  density, −3.1 points of plains and +2.1 / +2.0 points of hills / mountain. **The stop stands.**
- Not run, as the rule asks: the panel with A′, and the 2 048² plates with A′.

**Queued (the reviewer's, not in this finding)**:
- **F189, the rivers on the filled plains**: the fine fill by differential compaction (F183's Athy law applied to each fine
  column) instead of a flat fill. No parameter.
- **The islet mosaic of the shallow seas** (the amplification's relief under the sea), already there without the deposition.
  It is an amplification question, to settle with the author (against the coasts' dissection).
- **Then the chain's budget.**

## The closure: D0–D3 and R1's prediction

**F188's closure: diagnostics without a new bench, except two short runs** (outputs `f188_d1_output.txt`, `f188_d2_d3_output.txt`)

- **The author's answers (the reviewer's relay)**:
  - **Q1, A and A′: closed.** Neither is adopted. Their code stays behind flags off until D2 is reported, then is removed in the
    commit that adopts B or not.
    - **The reason**: the step dependence comes from the construction itself. The flexure lags one step, and the min is
      applied step by step instead of the min of the cumulated amounts.
    - The fluvial aggradation stays the one the effective K holds (Yuan JGR eq. 25), plus the lakes' fill by volume.
  - **Q2, the physics level's non-convergence in step: accepted if it is normal.** The test R1 checks it.
  - **Q3, the shallow seas' islet mosaic: accepted for now** (the author sees wetlands). It stays at the bottom of the queue, as
    an amplification question, not a deposition one.
- **D0**: the tables' columns now read « 80 kyr relative to 20 kyr », i.e. (value at 80 − value at 20) / value at 20.
- **D1, the production's step and the Courant number** (read in the code and on saved fields):
  - **The production's fluvial step is the step's own duration**: one fluvial step per physics step
    (`fluvial_dt_max = ∞`, `fluvial_substeps = 1`, `cascade.rs:628-630`).
    - The durations are `history.rs:605` (`step_durations`): at the intermediate age, 300 steps from 18.7 to 1 466 kyr,
      median 167 kyr.
    - Bounded at 80 / 20 / 5 kyr: 200 / 295 / 300 steps split, 1 404 / 5 157 / 20 150 fluvial sub-steps (against 300).
  - **The times, corrected**: the control at 80 kyr takes 128 s, against 93 s for F187 unbounded. The « 101 s » was A's run,
    with its deposit.
  - **C = K_SI A^0.5 dt / dx** (K_SI 2 × 10⁻⁶ yr⁻¹, dx 1 562.5 m):

  | A | 80 kyr | 20 kyr | 5 kyr | the production's median step (167 kyr) |
  |---|---|---|---|---|
  | 100 km² | 1.02 | 0.26 | 0.06 | 2.14 |
  | 1 000 km² | 3.24 | 0.81 | 0.20 | 6.76 |
  | 10 000 km² | 10.24 | 2.56 | 0.64 | 21.37 |

  - **The land share with C > 1** (seed 8, the D8 area of the saved field): 15.3 % at the production's median step, 7.3 % at
    80 kyr, 1.9 % at 20 kyr.
  - **The control's drift is not where C > 1** (the mask at the production's median step):
    - on C > 1 (5 018 cells) the classes move ≤ 0.52 point (the most: the mountain, −0.52);
    - on C ≤ 1 (27 689 cells) the hills move −3.29 points.
    - So the step dependence is not concentrated on the cells where the implicit incision exceeds Courant 1.
- **D2, the uniform mode (k = 0) of F183 A1's flexure** (seed 8, intermediate, F187; one run; a record-only addition to
  `SedRecord`, `w0_erosion_n` / `w0_load_n` / `basement_incised_n`):
  - **D2.1**: at the run's end, the erosional part is +928.5 m and the load part (the deposits, the recycling) −537.9 m. Their
    sum, **+390.6 m, is a uniform uplift**; it equals the cumulated rebound's mean. The series over time are not recorded, so
    they were not recomputed.
  - **D2.2, the mean load's budget**:
    - the basement incised: ρ_c V = 4.90 × 10¹⁷ kg (181 573 km³);
    - the deposits on land: ρ_g S + ρ_w (H − S) = 7.8 × 10¹⁵ kg (S 2 243 km³);
    - the deposits at sea: (ρ_g − ρ_w) S = 2.76 × 10¹⁷ kg (S 163 433 km³).
    - The expected mode, (ρ_c V_inc − ΔL_land − ΔL_sea) / (ρ_m A), is +390.6 m against the recorded +390.6 m. They agree:
      there is no double counting.
    - **Where the non-zero mean comes from** (my estimates from these terms):
      - **the water**: a marine deposit loads (ρ_g − ρ_w) S, but the water it displaces is never put back. The sea level is held
        fixed instead of conserving the water volume. ρ_w S_sea / (ρ_m A) ≈ +319 m.
      - **the recycling**: the cover subducted leaves the load, ≈ +90 m (from the budget's subducted term, the norm → m
        factor ≈ 10.9 km inferred from the sea's deposit).
      - **the hillslope's export to the sea** (7 063 km³): it loads the sea once deposited, but its removal from the land is
        not in A1's unloading, which flexes the incision's basement only. As basement, it would add +36 m to w0, so its
        omission accounts for ≈ −36 m.
      - A residual of ≈ +18 m is not attributed.
  - **D2.3, seed 8's east lake** (1 300 cells, 3 174 km²):
    - its cumulated flexure is **+436 m, an uplift**: the uniform mode +391 m (90 %) + the non-uniform part +45 m (10 %);
    - its cumulated tectonic source is +315 m.
    - **The lake's cells do not sink**, at the cumulated level: the subsidence reading does not hold for this lake.
  - **D2.4, A′'s zone** (17 926 cells):
    - the cumulated Δz_v (the source + the flexure) is negative on 122 cells only; without the uniform mode, on 1 379;
    - so removing w0 would **add** cells, not remove them: the uniform mode is an uplift;
    - A′'s zone is much broader than the cumulated negative cells, because A′ reads Δz_v step by step and it is negative on
      some steps of cells that rise in total.
- **D3, B on seed 23** (two runs; the map `images/d3_lakes_diff_F187_B_seed23.png`: red new, blue lost, grey kept):
  - all the lakes go from 4 006 to 6 565 km²;
  - **the new lakes** (≥ 16 cells, not overlapping F187's) are only 3, together 164 cells (400 km²):
    - 181 km² at (240, 48), 14 km from the coast;
    - 176 km² at (200, 248), 66 km;
    - 44 km² at (241, 244), 12 km.
  - On them (B − F187): the surface −13 m, the cover +18 m, the cumulated flexure −2 m.
  - **The flexure track is ruled out for the new lakes** (−2 m).
  - Their basement stands ≈ 31 m lower (the surface −13 m with the cover +18 m), and their outlets' surfaces are lowered by
    10–14 m. That is consistent with the soft cover's faster erosion opening basins, but not demonstrated.
  - **Most of the +2 560 km² is the existing lakes growing**, which this diagnosis did not analyse.
- **R1's prediction (the reviewer's, written before the computation)**: a first-order scheme gives between 5 and 20 kyr about a
  quarter of the changes between 20 and 80 kyr.
  - On all the land the hills should therefore move about 1 point (against 4.2), and density, λ and S₁₀ less than 1 %.
  - **The reading**:
    - if the change between 5 and 20 kyr is between 1/8 and 1/2 of the one between 20 and 80 kyr, the scheme behaves as a
      first-order one. The limit is then documented as normal: an error of about 4 points of hills at the production's step,
      due to the implicit scheme at Courant > 1, accepted for the budget;
    - otherwise, stop and report without interpreting.
  - **Noted before the run**: D1 already finds the control's drift on the cells with C ≤ 1, not C > 1. The attribution « due to
    the implicit scheme at Courant > 1 » would then not be the right one, even if the ratio holds.

## R1

**R1, Richardson: the result** (the control, F187 without a fluvial deposit, seed 8, intermediate, the fluvial step at 5 kyr: 729 s,
estimated 923 s; the mass ≤ 4.5 × 10⁻¹¹; output `f188_r1_output.txt`):

| all the land | 80 kyr relative to 20 kyr | 20 kyr relative to 5 kyr | the ratio | in 1/8 … 1/2 |
|---|---|---|---|---|
| density | +0.65 % | −3.99 % | 6.14 | no |
| λ | +3.50 % | −5.23 % | 1.49 | no |
| S₁₀ | −0.47 % | +0.47 % | 0.99 | no |
| plains | +0.49 pt | −1.41 pt | 2.88 | no |
| plateau | +1.59 pt | +0.39 pt | 0.25 | yes |
| **hills** | **−4.20 pt** | **+0.90 pt** | **0.21** | **yes** |
| mountain | +2.13 pt | +0.12 pt | 0.05 | no (below) |

| A′'s zone | 80 kyr relative to 20 kyr | 20 kyr relative to 5 kyr | the ratio |
|---|---|---|---|
| density | −0.83 % | −4.75 % | 5.76 |
| λ | −4.64 % | −12.03 % | 2.59 |
| S₁₀ | +0.05 % | +0.93 % | 20.21 |
| plains | +1.81 pt | −1.42 pt | 0.78 |
| plateau | +2.91 pt | +0.42 pt | 0.15 |
| hills | −6.20 pt | +0.22 pt | 0.04 |
| mountain | +1.49 pt | +0.77 pt | 0.52 |

- **F174's criterion between 20 and 5 kyr fails** on both zones (λ −5.2 % on all the land, −12.0 % on A′'s zone).
- **The hills on C > 1 / C ≤ 1** (the mask at the production's median step):
  - C > 1: 31.09 → 30.85 → 29.43 % (80 / 20 / 5 kyr);
  - C ≤ 1: 24.70 → 27.99 → 27.74 %.
- **Against the reviewer's prediction** (the hills ≈ 1 point; density, λ, S₁₀ < 1 %):
  - the hills, +0.90 point on all the land, held;
  - S₁₀, +0.47 %, held;
  - density (−3.99 %) and λ (−5.23 %) did not.
- **The reading, by the rule written before**: the hills' ratio (0.21) is first-order-like, but density, λ and the plains change
  as much or more between 5 and 20 kyr as between 20 and 80 (ratios 1.5–6.1). The condition « between 1/8 and 1/2 » does not
  hold across the quantities. **So: stop and report, without interpreting.**
- **F188 stops here.** No F189, no panel, no commit. A and A′ stay behind their flags off, to be removed in the commit that adopts
  B or not.

## The closure, continued

**F188's closure, continued**:

- **The reviewer's errors (recorded)**:
  1. supposing a subsidence by the uniform mode without checking its sign. The hypothesis on seed 8's lakes is refuted by D2.3:
     the east lake's cells rise by +436 m of flexure.
  2. a badly designed R1:
     - a step at 80 kyr outside the asymptotic regime (C up to ~10);
     - and no noise floor: the mountain's ratio of 0.05 was counted as a failure, though it shows convergence.
- **E1, the production against the fine step**: **the production's run is not saved for this exact configuration** (seed 8,
  intermediate, F187 unbounded).
  - The saved seed 8 fields are the control at 80 / 20 / 5 kyr and A / A′'s runs. Not rerun (93 s), as asked: the reviewer
    decides.
- **E2, the operators of the physics level's step and of the fluvial sub-step** (a: ∝ dt; b: implicit; c: once per
  application, not dependent on dt):
  - **once per physics step**:
    - the uplift's first part, `cascade.rs:661` (a);
    - the land mask / the ocean by connectivity, `:665` (c);
    - the rock law's factor, `:689` (c, the contact read at the step's start);
    - the erosion rebound, `:856` (a in the step's erosion; off when the load is on);
    - **the talus, `:883` (c: 4 passes, factor 0.5, the repose slope, whatever dt)**;
    - **the linear diffusion, `:895`, with its weight D dt / cell² recomputed per step, `:634-635` (a)**;
    - the sea held at its pre-hillslope value, `:900` (c);
    - the marine step, `:958` (a, explicit sub-steps ⌈4 K_M dt / dx²⌉, `deposition.rs:160`; its walk c per sub-step);
    - the recycling, `:1023` (a, the step's share of the interval);
    - the compaction, `:1039` (c, a state function H = T⁻¹(S));
    - the load's flexure, `:1095` (c: the step's load change, once).
  - **once per fluvial sub-step** (ks of them):
    - the uplift's part, `cascade.rs:721` (a);
    - the land mask, `:725` (c);
    - in the incision (`stream_power.rs`):
      - the routing on the pit-filled surface (`compute_flow`, `stream_power.rs:671`; the priority flood, `flow.rs:211`) (c);
      - the flats' gradient (Garbrecht–Martz), `flow.rs:228` (c, routing only);
      - the MFD area on the filled surface, `stream_power.rs:684` (c);
      - the base nodes, `:696` (c);
      - **the base-level floor at sea + ε, `:774` (c, a projection near the coast)**;
      - the implicit incision with n = 1, `:876` (b; its clamp to [h_r, h_o] cannot fire at n = 1);
      - the rock contact's correction, `cascade.rs:740-743` (c);
    - the fluvial pass (`deposition.rs`): the routing `:621` (c), **the lakes' fill by volume up to the spill, `:686`, `:731` (c)**,
      the deposit `:778` (a in the flux; 0 in F187).
  - **Inactive at 256²**:
    - the A_c gate (`stream_power.rs:833`: 0.1 km² is below one cell);
    - the threshold θ = 0 (`:848`);
    - A1's depression floor (`:790`, off);
    - the lateral erosion (`:888`: a half-width K_lat A^0.5 ≥ one cell needs A > 1.5 × 10⁵ km²);
    - the talus and the diffusion inside the incision call (off by `CascadeConfig::solver`).
  - **Flagged, without concluding**, the operators that touch the small areas and depend on the number of sub-steps:
    - **the hillslope (the talus c and the diffusion a) runs once per physics step, not with the fluvial sub-steps**. With ks
      sub-steps, ks incisions follow each other between two hillslope passes. The coupling of the hillslopes (small areas,
      C ≤ 1) to the channels changes with ks;
    - the routing on the pit-filled surface and on the flats, redone each sub-step, sets the heads of the network;
    - the lakes' fill up to the spill, each sub-step;
    - the base-level floor near the coast.
- **E3's rule (the reviewer's, written before the computation)**:
  - **the noise floor**: a quantity whose change between 10 and 5 kyr is below 1 % (density, λ, S₁₀) or 0.5 point (the classes)
    is called converged and leaves the test;
  - **the others**: a first-order scheme gives r = change(10, 5) / change(20, 10) ≈ 1/2. The accepted range is 0.25–1, with the
    same sign for the two changes;
  - **the reading**:
    - if every remaining quantity passes, the non-convergence is normal (first order, the asymptotic regime reached under
      ~10–20 kyr). The production's distance to the limit is then estimated by Richardson from E1, and documented;
    - otherwise, stop and report, with E2's (c) operators that touch the failing quantities.
- **D2a, the hillslopes' unloading**:
  - **confirmed in the code**: A1's flexed load, `cascade.rs:1069`, holds −ρ_c e_b with e_b the incision's basement only
    (`eb_vec`, `cascade.rs:812`). The hillslope's lowering (the talus `:883` and the diffusion `:895`, summed at `:909`) enters
    only as the marine source `:919`.
    - The cover it removes does leave ΔL (through `solid64`).
    - Its basement part is not unloaded on land, but loads the sea once deposited.
  - **The local effect**, the missing rebound on the exporting cells, needs the production run's cumulated hillslope field
    (`sum_diffusion`). It is not saved (E1). It can be read on the same 93 s run as E1, if decided.
- **D2b, the closed world (information for the author's question)**:
  - **the marine deposits** (seed 8, F187): 163 433 km³ of solids. The bulk volume (with the pores) is not in the saved
    record. It lies between the solids alone and the solids at the surface porosity φ₀ = 0.63, 441 700 km³.
  - **the ocean's area at the end**: ≈ 49.8 % of the 160 000 km², ≈ 79 700 km².
  - **the sea level's rise if the water volume were conserved**:
    - at first order, 2.05 km (the solids alone) to 5.54 km (at φ₀);
    - with Airy's factor for the water's load on the floor, 1 / (1 − ρ_w / ρ_m) = 1.454, 2.98 to 8.06 km.
    - That is first order, with the ocean's area fixed: a rise floods the land and widens the ocean.
- **D2c, the residual**: it came from my attribution, not from a missing term. The exact split of the recorded w0 (+390.6 m) by the
  same masses is:
  - the water displaced at sea, ρ_w S_sea / (ρ_m A): **+318.8 m**;
  - the solids missing, ρ_g (V_inc − S_tot) / (ρ_m A), with V_inc − S_tot = 15 897 km³ (the recycled cover, less the
    hillslope-sourced deposits, whose basement part is not in V_inc): **+81.9 m**;
  - the density difference, (ρ_c − ρ_g) V_inc / (ρ_m A) with ρ_c 2 700 against ρ_g 2 720: **−6.9 m**;
  - the pore water on land, −ρ_w (H − S)_land / (ρ_m A): **−3.3 m**;
  - **the sum is +390.5 m**.
  - My first split had estimated the recycling (+90 m) through an inferred factor and left out the density and the pore water.
    It also subtracted the hillslope's −36 m, which is an omission of A1 (D2a), not a term of the recorded balance.

## E3 and D3b

**E3, the step at 10 kyr: the result** (the control, seed 8; 414 s; the mass ≤ 3.8 × 10⁻¹¹; output `f188_e3_output.txt`), by the rule
written before:

| | all the land: 20 vs 10 | 10 vs 5 | the verdict | A′'s zone: 20 vs 10 | 10 vs 5 | the verdict |
|---|---|---|---|---|---|---|
| density | −1.70 % | −2.33 % | r +1.37, fails | −1.54 % | −3.26 % | r +2.12, fails |
| λ | −1.62 % | −3.67 % | r +2.27, fails | −2.12 % | −10.12 % | r +4.76, fails |
| S₁₀ | −0.09 % | +0.56 % | converged | −0.47 % | +1.41 % | r −3.01, fails |
| plains | +0.15 pt | −1.56 pt | r −10.2, fails | +0.20 pt | −1.62 pt | r −8.0, fails |
| plateau | −0.37 pt | +0.76 pt | r −2.05, fails | −0.04 pt | +0.46 pt | converged |
| hills | −0.12 pt | +1.02 pt | r −8.6, fails | −0.86 pt | +1.08 pt | r −1.26, fails |
| mountain | +0.34 pt | −0.22 pt | converged | +0.70 pt | +0.07 pt | converged |

- **On the masks** (the production's median step):
  - **on C > 1**: the plains (r +0.96) and the hills (r +0.60) pass; λ, S₁₀ and the mountain are converged; density and the
    plateau fail.
  - **on C ≤ 1**: density passes (r +0.54); S₁₀ and the mountain are converged; λ, the plains, the plateau and the hills fail.
- **The reading, by the rule: a quantity fails, so stop and report.**
  - The change between 10 and 5 kyr is larger than between 20 and 10, often of the opposite sign: the behaviour is not first
    order.
  - Not interpreted.
- **E2's (c) operators that touch the failing quantities** (density, λ, the plains, the plateau, the hills):
  - **the network's heads and spacing (density, λ)**:
    - the routing on the pit-filled surface, redone each sub-step (`stream_power.rs:671`, `flow.rs:211`);
    - the flats' gradient (`flow.rs:228`).
  - **the classes (the hills, the plains, the plateau)**:
    - **the hillslope, not sub-stepped with the fluvial stage**: the talus `cascade.rs:883` (c) and the diffusion `:895` (a),
      once per physics step after ks incisions;
    - the lakes' fill up to the spill, each sub-step (`deposition.rs:686`, `:731`);
    - the base-level floor near the coast (`stream_power.rs:774`).
  - The number of applications grows with the sub-steps: 5 157 / ≈ 10 000 / 20 150 at 20 / 10 / 5 kyr.
- **D3b, the lakes of seed 23 that grow with B** (two runs, 185 s; output `f188_d3b_output.txt`):
  - **the two measures**: « all the lakes » = every cell of lakes182's mask (land, filled > z), 4 006 → 6 565 km². « The large
    lakes » = its 8-connected components of ≥ 64 cells, 2 876 → 5 706 km² (5 → 7 lakes).
  - **The five existing lakes that add the most** (F187 → B):

  | the lake (cell) | the area | the spill | the floor min / median | the lake fill received | at B's sill: the cover, the surface | the cumulated flexure | why |
  |---|---|---|---|---|---|---|---|
  | (142, 171) | 212 → 2 209 km² | 23 → 40 m | 12 → 11 / 22 → 27 m | 1.5 → 21.6 km³ | 1 312 → 4 530 m, 23 → 40 m | +988 → +988 m | **its spill rises** (the sill buried under 3.2 km more cover) |
  | (35, 39) | 291 → 2 026 km² | 55 → 49 m | 10 → 21 / 54 → 45 m | 15.4 → 128.7 km³ | 5 951 → 3 627 m, 33 → 49 m | −319 → −334 m | **it merges 11 of F187's lakes** |
  | (136, 117) | 564 → 735 km² | 42 → 40 m | 42 → 40 / 42 → 40 m | 56.4 → 74.9 km³ | 0 → 0 m, 42 → 40 m | +1 387 → +1 454 m | its floor lowers |
  | (73, 234) | 95 → 186 km² | 63 → 45 m | 38 → 28 / 57 → 44 m | 50.0 → 86.5 km³ | 1 375 → 904 m, 63 → 45 m | +452 → +492 m | its floor lowers |
  | (143, 25) | 98 → 146 km² | 28 → 28 m | 28 → 28 / 28 → 28 m | 0.6 → 0.6 km³ | 0 → 2 m, 28 → 28 m | +1 259 → +1 276 m | not explained by its spill or its floor |

  - The two largest gains (+1 997 and +1 736 km², 3 733 of the +2 830 km² of large lakes) come from a sill raised by cover and
    from a merge.
  - The flexure changes little (≤ 67 m).
- **F188 stops here**: no F189, no panel, no commit.

## The open window, and V's setup

**The model's hypothesis on the water (the author's answer to Q1, confirmed): the window is open for the water.**
- The sea level is fixed, and the water the deposits displace leaves the map.
- The +318.8 m of uniform uplift (D2c) are its assumed consequence.
- The closed world is excluded by D2b (a sea-level rise of 2 to 8 km).

**V, the step's verification on a synthetic world (the setup and the predictions, written before the computation)**:
- **The setup**:
  - 256², the torus, the production's physics level through the same code path (`physics_history_dep`, the history
    timing, the sea by connectivity);
  - a fixed sea strip (the rows 0–7 of the 64² input) as the base level;
  - a uniform uplift U = 0.5 mm/yr on the land. The coast's step in U is ramped over four 64² rows, so that the bicubic does
    not ring and sink the coast; this is declared.
  - a uniform rock and K (K_SI 2 × 10⁻⁶, m 0.5, n 1; no rock law);
  - C1, the flexure, the deposition and the recycling removed; the initial relief is the level's fixed roughness (20 m RMS,
    a fixed seed).
- **The steps (declared)**: the physics step fixed at 160 kyr, the fluvial step bounded at 160, 40 and 10 kyr (`fluvial_dt_max`
  on `Deposition::nothing()`, which deposits nothing).
  - That varies the fluvial sub-step at a fixed physics step, as E3 did. The hillslope stays once per physics step, as in
    production.
- **The duration**: 60 Myr. The steady state is first reached when max |∂h/∂t| < 1 % of U on the land; its time is reported.
- **V1, without the hillslope** (the talus off by the new bench knob `CascadeConfig::talus_off`, the diffusion off by α = 0 and
  D_phys = 0). **The reviewer's prediction**: the steady state does not depend on dt.
  - h equal within 10⁻⁴ (relative) cell by cell between the three steps;
  - density, λ, S₁₀ and the classes identical;
  - h = (U/K) χ within 1 % on the river cells, χ = ∫ A^(−m) dx from the base level, m = 0.5, A₀ = 1 m² (Perron & Royden 2013,
    ESPL 38, 570–576).
  - **Two declared adaptations**:
    - the area is the incision's own, the MFD area (p = 6) along the D8 receivers;
    - the U of the ramp enters the integral, h = ε + ∫ U(x) / (K A^m) dx, with ε = 0.5 m the base-level floor.
- **If V1 fails**: V1 again at 160 and 10 kyr, cutting one (c) operator at a time, in this order: the base-level floor at sea
  + ε, the lakes' fill, the land mask per sub-step, the flats' gradient.
- **V2, with the production's hillslope** (the talus and the diffusion, once per physics step). **The reviewer's prediction**:
  - at most a first-order dependence, from the operators' splitting (Lie);
  - on the quantities above E3's noise floor, r = change(40, 10) / change(160, 40) ≈ 1/4, between 0.1 and 0.5, with the same
    sign.
- **If V2 fails while V1 passes**: V2 again with the hillslope split at the fluvial sub-steps' pace.

## E1 and V

**E1, the production against the fine step** (seed 8, F187 unbounded, 89 s, saved as `results/prod.bin`; output
`f188_e1_output.txt`). The raw readings, with no Richardson extrapolation (E3 excludes it):

| all the land | density | λ (km) | S₁₀ | plains / plateau / hills / mountain (%) |
|---|---|---|---|---|
| the production (median step 167 kyr) | 0.2537 | 21.34 | 1.1291 | 23.27 / 6.37 / 25.78 / 44.58 |
| 20 kyr | 0.2336 | 19.13 | 1.1399 | 23.72 / 3.34 / 28.83 / 44.11 |
| 10 kyr | 0.2376 | 19.44 | 1.1409 | 23.57 / 3.71 / 28.95 / 43.77 |
| 5 kyr | 0.2433 | 20.18 | 1.1346 | 25.13 / 2.95 / 27.92 / 44.00 |

- A′'s zone, C > 1 and C ≤ 1 are in the output.
- **Against 5 kyr, the production stands at**: density +4.3 %, λ +5.8 %, the plateau +3.4 points, the hills −2.1 points on all the
  land. On C > 1, the plateau +3.1 and the hills −1.4 points.
- **D2a, the missing rebound** (the hillslope's basement unloading not flexed; the cover's share counted as basement), at the
  run's end, on the 17 516 cells the hillslope lowers:
  - max 66 m, median 29 m;
  - on the 980 coastal ones, 22 % above 50 m;
  - its mean over the map is +29.9 m (the k = 0 part; D2's estimate from the export alone was +36 m).

**V, the result** (output `f188_v_output.txt`; 1 189 s, against the 9 min estimated):
- **V1, without the hillslope** (fluvial steps 160 / 40 / 10 kyr, physics step 160 kyr):
  - **the steady state (max |∂h/∂t| < 1 % U on the land) is not reached within 60 Myr at any step.** The land's median is
    864 / 859 / 859 m, its max 2 159 / 2 120 / 2 128 m.
  - **h cell by cell differs between the steps**: the relative difference has a median of 0.19–0.23 and a max of 30–37. The
    prediction (10⁻⁴) fails.
  - **density, λ, S₁₀ and the classes are under E3's noise floor** between 40 and 10 kyr. That is the prediction « identical »,
    up to the floor.
  - **h = ε + ∫ U / (K A^m) dx holds within 1 % at every step**: median 0.17 / 0.12 / 0.10 %, p90 0.87 / 0.58 / 0.55 % on
    ≈ 13 800 river cells.
  - **The reading**: each step's rivers sit on the steady law for its own network, but the networks differ between the steps,
    so h differs cell by cell. The criterion max |∂h/∂t| < 1 % U on every land cell is never met.
  - This is not interpreted further: V1's premise (the steady state reached) does not hold.
- **V2, with the production's hillslope**:
  - the steady state is not reached either;
  - **an anomaly**: the land's max is 4 449 / 12 677 / 3 030 m (12.7 km at 40 kyr), while the median is 1 327 / 1 218 / 1 129 m;
  - the χ law fails (median 26–36 %), as expected with a hillslope;
  - the ratios r = change(40, 10) / change(160, 40): density −0.89, λ +0.55, the hills +2.08, the mountain +1.42. All are out of
    0.1–0.5, and density changes sign. S₁₀, the plains and the plateau are under the floor.
  - **The prediction (first order at most) fails.**
- **What was not run, and why**: V1's failure branch (V1 again cutting one (c) operator at a time).
  - The reviewer's estimate was 2–3 min per pair; the measured runs give ≈ 7.6 min per pair (27 s + 428 s), so ≈ 30 min for
    the four.
  - Three of the four operators have no knob yet (the base-level floor, the land mask per sub-step, the flats' gradient). The
    fourth, the lakes' fill, is already off in V (`Deposition::nothing()`).
  - The branch's premise is not met: V1's steady state is not reached.
  - **Stopped and reported for the reviewer's decision.**

## P

**P, the plates of B for the author** (1 181 s; 2 048², F187 | F187 + B, the rivers by the corrected overlay, the lakes in cyan):
- `images/pb_rivers_lakes2048_F187_B_seed{8,18,19,23}.png`;
- on seed 23 the two lakes of D3b are circled in red: the raised sill at 256² (142, 171) and the merge at (35, 39).
- What the plates still show, as at F188's Part 0:
  - the islet mosaic of the shallow seas;
  - the flat filled surfaces with their solid patches and 45° lines;
  - many small lakes (cyan specks).
- **F188 stops here.** No F189, no panel, no commit.

## The final closure

**F188's final closure: the river layer was flipped vertically** (noticed by the author on the plates P).

- **My error (recorded)**: the river layer of F188's plates (Part 0, p0b, P) is flipped vertically. It had been « corrected » at
  p0b without checking that it sat on the relief.
- **The reviewer's errors (recorded)**:
  1. reading the plates without checking the layer's registration. The networks at sea were taken for the islet mosaic, and a
     false reading of « smooth zones without rivers » was proposed;
  2. V1: predicting a steady state identical within 10⁻⁴ cell by cell. The network's topology depends on the path taken, and
     the criterion max |∂h/∂t| < 1 % U is not reachable. The relevant criterion was the χ law per network, with the statistics,
     and V1 passes it;
  3. costs underestimated by a factor 2 to 3 (V and the ablation pairs).
- **The reviewer's reading (recorded)**: V1 shows the incision coherent in χ at every step, with statistics under the noise
  floor. E3's step dependence comes from elsewhere than the incision's scheme. The ablation branch is not run.
- **F1, the flip**:
  1. **Where the axis flips**: the bench's map renderers put the grid's row 0 at the image's bottom.
     - `alt_png` (`f188_cascade.rs:1659`, `img.put_pixel(x, m − 1 − y)`), `alt_png_plain` (`:1715`), `lakes_png` and
       `thick_png` (`255 − y`).
     - The layers I wrote for F187–F188 put row y at the image's row y: `rivers_png188`, `rivers_lakes_png188`, and the red
       circles of P.
  2. **The viz's Cascade window does not have the defect**: `cascade::hydro::draw` writes the pixel ((n − 1 − y) n + x)
     (`hydro.rs:121`), the shading's north-up convention (`ymir-viz/src/ui/cascade.rs:1579`, `:1609`).
  3. **What used the flipped layers**:
     - **the misregistered plates**:
       - F188's Part 0 (`p0_rivers2048_*`), p0b (`p0b_rivers2048_*`) and P (`pb_rivers_lakes2048_*`, the circles included);
       - F187's `r2_seed18_recycling_boundaries.png` (the boundaries drawn unflipped over a flipped map).
     - **Flipped but internally consistent** (no relief under them): F187's `ll_deflection_seed*.png` and
       `ll_deflection_rel_seed*.png`, and F188's `d3_lakes_diff_F187_B_seed23.png`. They are upside down against the other maps,
       but not misregistered.
     - **The measures were computed on the fields, not on the layer, so none is recomputed**:
       - Z5's S₁₀, `a_dir` and the order-≥ 2 pieces (`planform` on the 2 048² field);
       - F187's islet counts (on the fields);
       - the seed 18 shelf edge against the recycling boundaries (in the grid's coordinates on both sides).
  4. **The correction and the registration check** (`calage188`, run before every plate; a failing panel stops the plate):
     - the layers write the grid's row y at the image's row 1 023 − y;
     - the check requires ≥ 99 % of the drawn river pixels on the base map's land, and every trunk mouth (≥ 100 km², its
       receiver the sea) within one pixel of a sea pixel.

## F2, F3, G1

**F2, the islet mosaic on the relief** (no layer; seed 8, F187 at 2 048²): **12 land components below 1 km² in the shallow seas**
(−200 … 0 m), 4.3 km² in all; 15 (5.8 km²) with no depth condition.
- The mosaic read on the plates was mostly the flipped river layer drawn over the sea.
- F187's « islets » (11 000–48 000 fine land cells under a 256² sea cell) are coastal fringe cells, not isolated islands.

**F3, the plates of B redone** (the corrected layer; each panel checked first; outputs `f188_f3_output.txt`):
- **The check's first version was too broad**: it counted as a mouth a river ending in an interior basin below the sea, which the
  map draws as land (seed 8: 31 / 38). Corrected: a mouth's receiver is the connected ocean, as the map draws it.
- **The registration check**: 100.00 % of the drawn river pixels on land on every panel. The mouths within a pixel of the sea:
  - seed 8: 30 / 30 (F187) and 31 / 31 (F187 + B);
  - seed 18: 25 / 25 and 33 / 33;
  - **seed 19 F187: 47 / 48, so the check fails and the plate stops there.**
- **The plates written**: `images/f3_rivers_lakes2048_F187_B_seed8.png` and `_seed18.png`. Seeds 19 and 23 are not written.
  - The one mouth of seed 19 is not located. A mouth into a sea channel narrower than the map's pixel (2 fine cells) can stand
    more than a pixel from a drawn sea pixel.
- What the registered plates show: the rivers on the land only; the flat filled surfaces with their solid patches and 45° lines.

**G1, V2's 12.7 km peak** (outputs `f188_g1_output.txt`):
- **G1.1**: 65 cells end above 5 × the land's median (1 218 m), in clusters (e.g. x 235–237, y 56–58, up to 10 599 m).
- **The domain's land max per step is a succession of spikes, not one growth**: 25 259 m at 10 Myr, 3 070 at 20, 5 335 at 30,
  2 836 at 40, 4 588 at 50, 12 677 at 60.
  - The first exceeds 3 km at step 25 (4.0 Myr), on cell 16300 (x 172, y 63): the max jumps from 1 922 to 10 736 m in one step,
    then decays over the next steps (5 474, 4 944, …).
- **G1.2, the cell in that step**:
  - at step 24, h 300 m, the receiver (0, −1), the area 2.4 km² (one cell);
  - at step 25, h 10 736 m, its 8 neighbours 4.4–6.1 km. The step's incision −120 m, talus −92 874 m, diffusion −1 260 m.
  - The talus shed 93 km in that step, so before the hillslope the cell stood ~100 km high.
  - That rise is in none of the recorded operators (the incision, the talus, the diffusion).
- **G1.3, why the cell rises without eroding**: the rise is a **deposit of the fluvial pass**.
  - In `deposition::pass_inner`, the flux reaching a cell whose receiver is itself, neither the sea nor a lake, is stranded on
    that cell (`deposition.rs:724`, `p.stranded[k] += tot`).
  - The stranded part enters z with the pass's deposit (`cascade.rs:810` `fo.pass.deposit(k)` = fluvial + lake + stranded; z =
    `fo.z`, `:835`).
  - **V ran `Deposition::nothing()`, whose lakes are off**: every closed depression the hillslope opens strands its whole
    catchment's eroded flux on one cell, in one sub-step. The spike follows, and the talus spreads it.
  - **My error (recorded)**: that setup (the lakes off in V) puts a stranding that production does not have into V2, and to a
    lesser degree into V1. **V2's verdict and its spikes are artefacts of the setup.**
- **G1.3, whether it can act in production**:
  - **the production cascade** (no deposition) runs no pass, so no.
  - **F187's configuration** (the lakes on): a closed depression is a lake filled by volume up to its spill. The stranding remains
    only on a cell whose receiver is itself outside the sea and the lakes (a D8 dead end on a flat). F186–F188's budgets give it
    0.2–2 % of the eroded matter.
- **G1.4, on the production run** (`results/prod.bin`, seed 8, F187): **2 land cells stand more than 500 m above all their
  neighbours**, (x 160, y 29) at 1 496 m and (x 216, y 123) at 553 m.
  - A hypothesis, not checked: the same stranding on a D8 dead end.
- **F188 stops here.** B's adoption waits for the author's look at F3's plates (seeds 8 and 18 written; 19 stopped by the check;
  23 not done). The F189 prompt (unfolding the torus) follows this report.

## The last closure

**F188's last closure**:
- **CC's errors (recorded)**:
  - attributing at p0b the networks « at sea » to the islet mosaic (the layer was flipped);
  - V's setup with the lakes off;
  - the registration check's first version, too broad.
- **The shallow seas' islet mosaic leaves the queue**: F2 finds 12 of them, 4.3 km². The author's Q3 (the mosaic accepted as
  wetlands) was answered on misregistered plates; the question has no object any more.
- **The author's decision on B (verbatim)**: « ok B adopté ».
  - B enters the adopted deposition profile, with F187's configuration without K_eff.
  - A's and A′'s code and flags are removed.
  - The commit comes after F3 bis, V2 bis and G1 bis, and after the full checks.
- **V2 bis's prediction (the reviewer's, written before the computation)**:
  - **the setup**: V2's (the production's hillslope, once per physics step), but **with the lakes filled by volume as in F187,
    with no fluvial deposit** (`Deposition { lakes: true, g: 0, marine: None, ..nothing() }`); fluvial steps 160, 40 and 10 kyr,
    the physics step 160 kyr, 60 Myr.
  - **the prediction**:
    - no land cell above 5 × the land's median, at any step;
    - on the quantities above E3's noise floor, r = change(40, 10) / change(160, 40) between 0.1 and 0.5, with the same sign
      (Lie splitting, first order).
  - **the reading**:
    - if V2 bis passes, the fluvial core and the hillslope are coherent. E3's step dependence comes from the coupling with
      what V lacks (C1's forcing applied once per step, the moving coast, the flexure, the marine step). It is documented as an
      accepted limit (the author's answer to Q2): the production stands from the 5 kyr step at +4 % density, +6 % λ, 2–3 points
      on the classes.
    - if V2 bis fails, stop and report.

## V2 bis, G1 bis, H

**V2 bis: the result** (output `f188_v2bis_output.txt`; 551 s):
- **The spikes stay with the lakes filled by volume.** The land's max over its median, the worst over the steps: **48.2** at 160 kyr
  (step 105), 22.3 at 40 kyr (step 23), 7.7 at 10 kyr (step 26). Every value is above the predicted 5.
- **The ratios r = change(40, 10) / change(160, 40)**: density −0.66, λ +2.47, S₁₀ +1.43, the plateau −0.95, the hills −19.2. All
  are out of 0.1–0.5. The plains and the mountain are under the noise floor.
- **The prediction fails: stop and report.** The hillslope's splitting at the sub-steps' pace is to be decided.

**G1 bis: the stranding in F187's configuration** (seed 8, one run, 73 s; output `f188_g1bis_output.txt`; with G = 0 the record's
`fluvial_m` is the stranded thickness alone):
- **the stranded total over the run: 3 487 km³.**
- **92 cells got more than 10 m stranded**, up to 577 522 m cumulated on (x 169, y 108), 264 773 m on (122, 116), 92 915 m on
  (160, 113)…
  - At the run's end these cells stand at 86–2 372 m, near their neighbours: the talus and the incision have spread the
    stranded mass.
- **The two spikes**:
  - **(160, 29)**, 1 496 m against 725 m for its highest neighbour: 3 368 m stranded there, the cover 1 501 m (lake fill 0 m,
    marine 14 m). **Its height comes from the fluvial stranding** (`deposition.rs:724`).
  - **(216, 123)**, 553 m against 30 m: 0 m of fluvial stranding, but 2 727 m of marine deposit, the cover 656 m. **Its height
    comes from the marine walk's stranding** (`deposition.rs:199-211`), not from `p.stranded`.
- **Why these cells were D8 dead ends outside the lakes**: not determined.
  - At the run's end every listed cell has a D8 receiver and is not in a lake.
  - The stranding happened on the steps when they were dead ends: in `pass_inner`, the cells with `receiver == k` outside the
    sea and the lakes, i.e. `DIR_NONE` from `compute_flow`.
  - The run does not record those steps; only an instrumented run can say whether it is the flats' routing, the lake mask or
    another cause.
- **Not corrected.** The candidate correction (any dead end outside the sea a one-cell lake filled by volume up to its spill) goes
  to a separate finding.

**H, B's adoption: not done.**
- The author's condition: « Si l'un de ces trois points demande un changement de code au-delà de l'affichage, ne commite pas :
  rapporte. »
- V2 bis fails (the stop) and G1 bis shows the stranding at work in F187's configuration: both call for a code decision (the
  hillslope's splitting, the dead ends as lakes).
- So no commit. A and A′ stay behind their flags off, and B behind its flag.

## F3 bis

**F3 bis: the result** (676 s; output `f188_f3bis_output.txt`):
- **The amended check passes on every panel**: 100.00 % of the river pixels on land. The mouths pass:
  - seed 19: 48 / 48 (F187) and 53 / 53 (F187 + B);
  - seed 23: 40 / 40 and 41 / 41.
- **The registered plates**:
  - `images/f3_rivers_lakes2048_F187_B_seed19.png` and `_seed23.png`;
  - with the D3b circles at their registered places on seed 23 (the raised sill in the upper half, the merge at the bottom
    left against the mountain);
  - and `_seed8.png` and `_seed18.png` from F3.
- **The dark-green steps of seed 8's east coast with B** (the right panel, x ≈ 1 850–1 950 px, y ≈ 250–480 px): **they are the
  256² cells `bury_smooth` flattens**.
  - In that region 627 of the 799 covered 256² cells have their cover thicker than their fine basement's relief, and their fine
    surface is replaced by the bilinear 256² surface (map `images/f3bis_seed8_B_eastcoast_smooth_edge_flatmask.png`: bury_smooth,
    bury_edge, the mask of the flattened cells).
  - These patches carry no fine relief, hence their uniform colour and their sparse network. Their edges follow the 256² grid,
    because the switch is decided per parent cell (`buried[p]`, `amplify.rs`, `Transmit::bury_with`): **a read through the
    parent cell, against the author's rule**.
  - On the crop, 6 516 of the 88 740 land fine cells differ by more than 1 m between `bury_smooth` and `bury_edge`.
  - Not corrected: it is a code change, for a separate finding.
- **F188 stops here, nothing committed.** B's adoption (« ok B adopté ») waits, as the author asked, since V2 bis and G1 bis call
  for code decisions:
  - the hillslope's splitting;
  - the dead ends as one-cell lakes;
  - and `bury_smooth`'s read through the parent cell.

## The step commit (B adopted)

**The step commit: B adopted** (the reviewer's go-ahead after F188's last closure):
- **The reviewer's error (recorded)**: the condition « si l'un de ces trois points demande un changement de code, ne commite
  pas » was badly worded. It aimed at a fix made during those points, not at existing defects they would bring to light. CC's
  cautious reading (no commit) was defensible.
- **CC's error (recorded)**: G1's conclusion, « V2's verdict and its spikes are artefacts of the setup », was premature. V2 bis
  shows that the spikes stay with the lakes on (48.2 / 22.3 / 7.7 × the land's median); cutting the lakes only made them worse.
- **The adoption** (the author: « ok B adopté »):
  - `Deposition::continent_default()` is the adopted deposition profile: F187's configuration without K_eff (the load, the
    compaction, the recycling, the live marine step, G = 0, `fluvial_keff: None`) with B (`soft_cover`) and the deposition
    handed to the amplifications (`transmit`).
  - It goes with `Sea { connected: true }` at the physics level, C1's consumption for the recycling, and in the chain
    `AmpConfig::sea_strict` with `Transmit::bury_smooth`.
  - The bench's `a188(false, true)` with `transmit` equals it (`f188_profile_is_b`).
  - The deposition stays off by default: production does not run it, and its hash does not move.
  - **Not wired in the viz**: the box « Dépôt (niveau de physique) » still runs F182's declaration (`Deposition::declared`).
    The viz passes no C1 consumption and runs no burial, so the profile cannot run there as adopted. That wiring is queued.
- **A and A′ removed**: the fields `subsiding_fill` and `accommodation_fill`, `pass_routed_g` and `pass_routed_cap`, and the
  per-cell G and cap fields with the last load flexure that fed them (`cascade.rs`). Their controls in `the_f188_controls` go
  too; B's control stays, plus the adopted profile's mass control. Their benches (`f188_a3`, `f188_a3p`, `f188_tier1_ap`) are
  removed, and their outputs stay in the report.
- **The known defects, carried into the adopted profile**:
  - the fluvial and marine stranding on dead ends (G1 bis: 3 487 km³ on seed 8 in F187's configuration, 92 cells above 10 m);
  - `bury_smooth`'s read through the parent cell (F3 bis: the 8 × 8 flattened blocks of seed 8's east coast);
  - the step dependence (E3; V2 bis fails);
  - the hillslopes' unloading (D2a).
- **The order agreed with the author**: F189 the stranding, F190 the torus' unfolding, F191 the hillslopes' unloading, F192
  the fill by compaction with `bury_smooth`'s parent read removed, then the budget.
