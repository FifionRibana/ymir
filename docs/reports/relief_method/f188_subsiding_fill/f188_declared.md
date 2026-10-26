# Finding 188 — declared before the measures

Starting point: F187's configuration, not committed. It is F186's deposition without K_eff (`a187`), with:
- the per-cell marine step, the load, the compaction, the interpolated recycling;
- A (the sea by connectivity);
- B with `bury_smooth`, and the chain aligned (`sea_strict`).

## The author's words (to quote in the ADR)

« Seed 19: Ca donne un aspect wetland, des marécages ou des archpels cotiers. Seed 18: je n'ai pas vu ce lac. Seed 8: Oui
effectivement c'est un peu gros. Les haut fonds cotiers sont très bien c'est quelque chose qui me manquait. Pourquoi on avait
retiré le remblaiement des rivières? Autre point, on ne voit pas les rivières justement. »

**The answer to his question**: the fluvial deposition went through K_eff. K_eff was removed because the cascade's K is an
effective K, fitted on natural rivers, whose slope already holds the deposition's effect (Yuan JGR eq. 25), so K_eff counted the
deposition twice (F187). Part A brings the fluvial deposition back only where the ground sinks.

## Part 0 — see before correcting

- **0.1, the rivers on the plates.**
  - **The threshold**: the project's river threshold, `cascade::hydro::river_threshold_km2` (`hydro.rs:39`): max(1 km², 4
    cells), declared at F162–F163.
    - It is the one the viz's Cascade window draws (`ymir-viz/src/ui/cascade.rs:1478`: « rivières : aire ≥ … km², troncs ≥
      100 km² en bleu foncé »).
    - The plates use the same: the network (≥ the threshold) in light blue, the trunks (≥ 100 km², `planform::TRUNK_KM2`) in
      dark blue.
  - **The filled zones**: the fine land cells whose 256² parent carries a cover (H > 0).
  - **The instruments, on those zones only** (`planform` with the complement as its exclusion):
    - S₁₀, the sinuosity of 10 km D8 walks on the network ≥ 10 km² (I2, F163);
    - **the share of straight reaches**: `a_dir` (I1, F163), the share of the order-≥ 2 rivers' 8-step chords within ±5° of
      a D8 direction (isotropic directions read 22.2 %). It is the project's instrument for straight, axis-bound reaches.
  - Both with and without the deposition; seeds 8, 18, 19, 23 at 2 048².
- **0.2, seed 18's « rectangular » lake (already read on F187's lake map, `ll_lakes_full_seed18.png`)**:
  - the lake (1 761 cells) touches both the right edge (x = 255) and the left edge (x = 0) of the map;
  - **it is cut by the torus' edge** and continues on the left side;
  - its straight side is the image's frame, not a quantity read by its parent cell.

## Part A — the fill of the sinking zones (`Deposition::subsiding_fill`)

- **A1, the split, without a parameter**. On each land cell and step, the non-erosive vertical motion is Δz_v = the tectonic
  source of the step (the uplift du) + the load's flexure.
  - **The load's flexure is the previous step's (declared)**. The step's own is computed at the step's end from the step's
    deposits (`cascade.rs`, F183 A1), after the fluvial pass, so the pass reads the last one computed. The first step reads 0.
  - **Δz_v ≥ 0**: G = 0 on the cell (the effective K holds the deposition's effect on the slope).
  - **Δz_v < 0**: the fluvial pass deposits after Yuan (JGR eq. 13–15), D = min((G/Ñ) Q_in, Q_in), G = 1.55 (Guerit et al.
    2019, p. 855). What is not deposited goes downstream.
  - **How**: the pass takes a per-cell G (G on the sinking land cells, 0 elsewhere). With the flag off it takes the scalar G
    as before (0 in `a187`).
- **A2**: no iteration. The deposit enters the cover (the compaction, the load); the mass is counted by the existing controls.
- **A3, the convergence test**:
  - seed 8, intermediate, the fluvial step bounded at 80 kyr against 20 kyr (`fluvial_dt_max`, F185);
  - F174's criterion restricted to the sinking zone: the cells that sank at least once in the 20 kyr run;
  - density, λ, S₁₀ within 5 %, the classes' largest move < 1 point.
  - **If it fails: stop and report**, without widening.
- **A4's controls**:
  - the flag off = F187 bit for bit;
  - no cell sinking (a world with U ≥ 0 and no load) = no change;
  - the mass within 10⁻⁹.

## Part B — the sediment's erodibility (`Deposition::soft_cover`)

- **B1**: the deposited sediment (the cover H) erodes with K2's soft cover's K (K_couv). That is ρ_K = 10 times the basement's,
  Stock & Montgomery 1999 (F171).
  - **Where**: the rock law's factor (`cascade.rs:684`) is 1 (K_couv) on a cell with H > 0 at the step's start.
  - **The crossing step's correction** (`cascade.rs`, the `incise` closure) counts the remaining cover as max(the rock's
    remaining cover, H): the lowering beyond it is divided by ρ_K.
- **B2, the target**: the ~1 km marine piles (seed 23, the north coast).
  - The stranding's part (`deposition.rs:199-211`) is measured: a per-cell record of the walk's stranded thickness,
    `SedRecord::marine_stranded_m` (new).
- **B3's controls**:
  - the flag off = F187 bit for bit;
  - H = 0 everywhere = no change;
  - the mass within 10⁻⁹.

## Part C — the measures (on mains; the budget announced before each tier; < 45 min each)

- **Tier 1**, the panel at 256², intermediate: F187; + A; + B; + A + B.
  - The quantities:
    - the mass, the circumnavigation against the same world without the deposition, (v), W₁₀;
    - the plains, the mountain, the median, the shelves;
    - the lakes: large, area, and seed 8's east lake (the largest lake of F187's run, its area);
    - the depressions left; the marine deposit above 100 m;
    - the cost.
  - Split in two sub-tiers (F187 + A, then + B and + A + B) to stay under 45 min.
- **Tier 2**, the chain to 2 048², seeds 8, 18, 19 and 23, the configuration retained: the maps without | with, the rivers
  overlaid, and the lake map.

## The rules (the author's, written before)

- **The stops**: the mass control fails; a solver or an iteration to convergence; A3 fails; the circumnavigation lost to the
  deposition on ≥ 2 seeds.
- **No adoption in this finding** without the author's look at the plates with the rivers.

---
**Amendment (before any measure)**: A3's zone is the cells with a fluvial deposit in the 20 kyr run, not « the cells that sank
at least once ».
- The run does not export the set of sinking cells.
- G ≠ 0 only on sinking cells, so the deposit's cells are the zone where the fill acted.
