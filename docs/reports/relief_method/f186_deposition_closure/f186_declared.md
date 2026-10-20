# Finding 186 — declared (Part 0 measured first, as the brief asks; A and B declared before building)

**The budget rule (the author's, verbatim)**: « Oui ok. ET, il faut une extrapolation à 8192². Ca ne peut pas durer des heures.
Dans la version avant multigrid j'étais à 5 minutes max. Ce serait bien de ne pas dépasser cette cible. » The continent profile's
complete chain, C1 to 8 192², must take at most 5 minutes of wall time per world.

## 0.1 — the chain's budget (seed 23, intermediate, F181's default, no deposition; `f186_budget`)

**Caveat: the machine ran on battery (discharging, 75 %) during this bench.** Its CPU was throttled. The same physics level took
66.1 s CPU here against 26.4 s at F184 and F185 on the same seed and age, so every figure below is about 2.5 × the machine's
mains speed. Re-measure on mains power before any decision.

| stage | what it does (file:line) | time steps? | wall (battery) | wall ÷ 2.5 (estimate) |
|---|---|---|---|---|
| the selection | `select_continent_default_margin` (`cascade/profile.rs:192`): every candidate's C1 run to its end (224 candidates, 203 runs, 16 plates) | C1's 300 steps per run | 271 s | ~110 s |
| C1 | `run_c1_sources_mapped_ext` (`cascade/history.rs:359`) on the chosen continent | 300 steps | 3.0 s | ~1 s |
| the physics level, 256² | `physics_history_sea` → `run_loop` (`cascade.rs`): uplift, implicit incision, rebound, hillslope | 300 steps | 62.3 s | ~25 s (F185 measured 19–21 s) |
| amplification 512² | `Chain::next_level` (`cascade/amplify.rs:1398`): upscale and warp, N1 incision (`:1071`), talus (`:1160`), deposit (`:1203`) | iterations, no physical time (6 / 1 000 / 700) | 16.5 s | ~7 s |
| amplification 1 024² | the same (6 / 500 / 200) | iterations | 27.1 s | ~11 s |
| amplification 2 048² | the same (6 / 1 500 / 150) | iterations | 116.8 s (818 s CPU) | ~47 s |
| amplification 4 096² | **does not exist** for the continent profile (the benches stop at 2 048²). Extrapolated: the 2 048² cost per cell-iteration × cells × `Budget::moyen`'s 4 096² row (a declared PROXY) | iterations | ~429 s | ~170 s |
| amplification 8 192² | **does not exist**; extrapolated the same way | iterations | ~1 714 s | ~690 s |

- **The total**: to 2 048², 497 s measured (≈ 200 s on mains). To 8 192², ≈ 2 640 s with the extrapolation (≈ 1 060 s on mains),
  **against the 300 s target: about 3.5 × over on mains**.
- **The two most expensive stages**:
  - the 8 192² amplification (extrapolated, ~65 % of the chain);
  - then the 4 096² amplification and the selection, about equal (~170 s and ~110 s on mains).
- **The extrapolation's law**: cost ∝ cells × the budget's iterations, operator by operator, from the 2 048² level. It assumes
  the per-cell-iteration cost stays the 2 048² one (the memory traffic may make it larger), and it takes `Budget::moyen`'s 4 096²
  and 8 192² rows, which the code calls a PROXY.

## 0.2 — the deposition elements at 256²

- **One run of F184's live combination** (seed 23, intermediate, k = 1, on battery): 189.4 s wall, against 62.3 s without the
  deposition. Of it, the marine step took 109.1 s wall (65 691 sub-steps). So the marine step is about 3/4 of the deposition's
  added cost; the rest is the fluvial pass, the lakes, the load (one FFT per step), the compaction and the recycling.
- **From F184's mains-power costs** (CPU per world, intermediate): without the deposition 21–25 s; F183's combination (the
  fluvial pass, the lakes, the fixed marine step, the load, the compaction, the recycling) 42–44 s; with the live marine domain
  56–61 s. So the live marine step adds ~15 s over the fixed one, and the deposition as a whole adds ~35 s.
- The run did not time the load, the compaction and the recycling separately. They run once per step on the grid (one FFT,
  table reads, one multiplication), against the marine step's 65 000 grid sweeps.

## A — the fluvial closure at 256² (`Deposition::fluvial_keff`)

- **A1**: the explicit fluvial pass's deposition is off (the pass runs with G = 0: it routes, fills the lakes and delivers to the
  sea). The implicit incision uses, cell by cell, **K_eff = K / (1 + G (Ñ − 1)/Ñ)**.
  - G = 1.55 (Guerit et al. 2019, p. 855), p̃ = 1.
  - **Ñ is the incision's own area in cells**: the MFD accumulation (`relief_v3`'s MFD exponent, the profile's `mfd_exponent`)
    that `incise_with_floor_torus` uses for A^m. It is read inside the incision's loop, so there is no extra routing.
  - The source: Yuan JGR eq. 24–25 (the steady state: K_f A^m S^n = U (1 + G/p̃), the discrete factor 1 + G (Ñ − 1)/Ñ from
    eq. 13–14), eq. 28 (the profile × (1 + G/p̃)), eq. 31 (the response time × (1 + G/p̃)). Dividing K by the factor reproduces
    both the steady slope and the response time.
  - No iteration, no extra pass.
- **A2**: the eroded matter goes downstream with no deposition on the hillslopes or in the channels. It deposits only in the
  lakes (by volume to the spill) and at sea (the live domain, the delta's progradation, the marine diffusion). The load, the
  compaction and the recycling stay as at F183.
- **A3, declared as lost**: the aggradation of low, non-uplifted plains outside the lakes and the deltas (Yuan's deposition
  there is a transient that K_eff does not carry). It is judged on the plates.
- **A4, the controls (permanent)**:
  - the flag off is F185 bit for bit;
  - G = 0 is the incision without deposition;
  - F183's amended mass control at 10⁻⁹;
  - **the analytic test** (F185's, the version without the hillslope): with K_eff the slope ratio equals the discrete factor
    within 5 %.
- **A5, the budget**: the physics level with the deposition at most 60 s wall per world at the intermediate age, measured on mains
  power. The old age is reported.

## B — the transmission to the amplifications (`Deposition::transmit`)

- **B1, the inventory** (`cascade/amplify.rs`):
  - each amplification level receives the previous level's result (m), the D5 mask and ln K (`Chain::start_level`, `:1368`);
  - it upscales and warps them (`AmpLevel::upscale_warped`, `:826`), then runs its own operators, iterations without physical
    time: the N1 incision (`erode_n1`, `:1071`: `incise_with_floor_torus` at dt = 1), the talus (`:1160`) and a
    transport-capacity deposit (`:1203`).
- **B2, the burial closure**:
  - the physics level hands over the cover thickness H, the basement z − H, the lakes' levels (the routing's filled surface) and
    the sea level;
  - **the chain runs on the basement**: `Chain::new` takes z − H, so the relief each level adds applies to the basement, and the
    levels' own operators work on it as today;
  - **each level's output** (not its working field) is max(the basement amplified, the cover surface interpolated) on the cells
    where the interpolated H > 0. The cover surface is the physics level's z, bicubic to the level. The sediment fills the fine
    hollows up to the surface computed at 256²: no parameter, and the mean volume is kept where the added relief is centred;
  - **the lakes**: inside a coarse lake (`filled > z` at 256², the nearest parent) **that carries a cover** (H > 0), the fine
    cells are raised to the coarse cell's flat lake bed. A coarse lake with no sediment is left as it is: its water covers the
    fine hollows. This amendment, made before any measurement, is what keeps the control « H = 0 changes nothing » true.
  - **the sea**: unchanged.
  - **Combined with the levels' own erosion**: the closure is applied to each level's output and is not fed back into the next
    level's working field. The next level amplifies the basement, so the cover is not eroded by the amplifications' N1 and
    deposit operators (a declared simplification).
- **B3, the controls**:
  - the flag off is today's amplifications bit for bit;
  - H = 0 everywhere changes nothing;
  - the volume of the fine cover against the coarse cover is reported.

## The fallback (the author's, declared before tier 1; used only if the 60 s budget fails because of the marine diffusion)

- **The marine diffusion runs on a 128² grid** of 2 × 2 blocks of 256² cells.
- **The restriction**: a block's source is the sum of its 4 cells' sources; its room under the sea is the sum of their rooms.
- **F184 A1's live-domain rules apply to the block**: the domain, the seaward walk, and the inflow capped at the room.
- **The prolongation**: a block's deposit is shared among its 4 cells in proportion to their room under the sea (h_sea − h). The
  mass is kept, and no cell goes above the sea.
- **The stability limit** is computed at 128² (Δx doubled, so a quarter of the sub-steps).
- **The controls**: the mass at 10⁻⁹; and at 256² the variant without blocks is F186 bit for bit.
- **The justification**: the diffusion length √(K_M d) exceeds 2.7 km from the shortest step (18.7 kyr, K_M = 400 m²/yr), against
  a 128² cell of 3.1 km. Reported, with the time saved.

## Tier 2, as the author adjusted it

- The 4 096² and 8 192² levels do not exist for the continent profile. The chain runs to 2 048², with the deposition and the
  transmission; 4 096² and 8 192² are extrapolated as in Part 0.
- The plates are at 2 048². The memory is checked before the run.
