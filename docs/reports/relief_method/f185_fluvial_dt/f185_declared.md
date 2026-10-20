# Finding 185 — declared before any measurement

The reference is F184's configuration with the live marine domain: F183's combination (load, compaction, recycling) and
`marine_live`, on F181's default profile, at 256². The new flag is `Deposition::fluvial_dt_max` (Δt_f, years). Δt_f = ∞ is F184
at k = 1 bit for bit.

## 0.1 — the equilibrium of Yuan's equation (Yuan et al. 2019, JGR, `docs/refs/yuan2019_2.pdf`)

- **The article gives it**: §3.1, p. 9, eq. 24 (the evolution equation is eq. 11, p. 6):

  S = [(1 + G/p̃) U / K_f]^(1/n) A^(−m/n) p̃^(−m/n),  i.e.  **K_f p̃^m A^m S^n = U (1 + G/p̃)**.

  - The steepness index becomes k′_s = (1 + G/p̃)^(1/n) (U/K_f)^(1/n) p̃^(−m/n) (eq. 25): « the effect of sediment deposition is
    to increase the steepness index by a factor of 2 » for G/p̃ = 1 and n = 1.
  - The steady profile is (1 + G/p̃)^(1/n) times higher (eq. 28, p. 10).
  - The response time is (1 + G/p̃)^(1/n) longer (eq. 31).
  - The origin, in Yuan's words (p. 6): « At steady state, the catchment area A in the deposition term vanishes, and the
    deposition term is equal to GU/p̃ ».
- **The reviewer's hypothesis** K A^m S^n = U (1 + G/p̃) is **confirmed**, up to notation: the factor p̃^m sits with K_f in
  Yuan's form, and p̃ = 1 here.
- **The cascade's exponents**: m = 0.5, n = 1 (`StreamPowerConfig::relief_v1`, inherited by `relief_v3`). With G = 1.55 and
  p̃ = 1 the predicted ratio of the slopes with and without the deposition, at the same drainage area, is
  **(1 + G)^(1/n) = 2.55**.
- **The scheme's own discrete equilibrium** (derived here, from eq. 13–14):
  - at a steady state every cell's net change is zero, so the net erosion of each upstream cell is U Δt;
  - the deposit at cell i is then (G/p̃)/Ñ_i × Σ_{j ∈ ups(i)} U Δt = (G/p̃) U Δt (Ñ_i − 1)/Ñ_i, since ups(i) holds the Ñ_i − 1
    cells strictly upstream (eq. 14, Ñ_i counting i);
  - the cell's balance U − E_i + D_i = 0 gives **E_i = U (1 + (G/p̃)(Ñ_i − 1)/Ñ_i)**, so the slope ratio at the cell is
    (1 + (G/p̃)(Ñ_i − 1)/Ñ_i)^(1/n);
  - with G = 1.55: 2.395 at Ñ = 10, 2.535 at Ñ = 100, and the continuous 2.55 as Ñ → ∞. It is within 2 % of eq. 24 from Ñ ≈ 31,
    and 6 % below it at Ñ = 10;
  - our pass deposits min((G/Ñ) Q_in, Q_in) with Q_in the upstream flux (eq. 13–15), so the same discrete factor holds.
- **The judged reference** (0.2): eq. 24's factor 2.55, **and** the discrete factor per area bin, both reported. The brief's
  5 % is judged on the discrete factor at A ≥ 10 cells (the scheme's own equilibrium); the continuous one is reported, and is
  expected within 5 % only from Ñ ≈ 20.

## 0.2 — the analytic test (permanent; `cascade::history::f169_tests::the_yuan_equilibrium`, `#[ignore]`d for its cost)

- **The basin**:
  - a synthetic history on the torus at 128² (cell 3.125 km on the 400 km domain): a land plateau rising gently from a sea
    strip (rows < 8 of the 64² snapshots), with the cascade's roughness (`roughness_m` = 20 m) to seed a network;
  - a uniform uplift U = 10⁻⁴ m/yr on the land (the sea strip held fixed), linear time, run until the change per step falls
    under 1 % of U dt (the cascade's `eq_ratio`) or to its end.
- **The real cascade**: the implicit incision (`relief_v3`, the MFD area for the incision), the hillslope (the talus, the linear
  diffusion), the deposition (the fluvial pass with the lakes, the live marine step) in fluvial sub-steps bounded by Δt_f. No
  rebound, no rock law, no load, no compaction, no recycling.
- **Two runs**: G = 1.55 against G = 0 (the same scheme, the same lakes and marine step).
- **The measure**:
  - on the final field, every land cell with a D8 receiver on land;
  - its slope S (m/m) and its D8 drainage area Ñ (cells);
  - the area bins Ñ ∈ [10, 20), [20, 40), [40, 80), … ;
  - per bin, the ratio of the median steepness S Ñ^(m/n) with and without the deposition.
  - The ratio is reported per bin against the discrete factor (1 + G (Ñ̄ − 1)/Ñ̄) (Ñ̄ the bin's median area) and against 2.55.
- **The rule**:
  - the ratio over all cells with Ñ ≥ 10 must equal the discrete factor at their median Ñ within 5 %, at the Δt_f retained in
    A;
  - it is measured also at 4 × Δt_f, to show that the gap falls as the step falls.
- **The departures I expect, labelled hypotheses**:
  - the cascade's hillslope runs on the channels (`diffuse_channels`), and it has a lateral erosion (`relief_v2`) and a base
    level floor near the sea; each adds terms Yuan's eq. 24 neglects (« we neglect the hillslope processes », p. 8);
  - if the test misses by more than 5 %, these are the first suspects, measured by switching them off in the test before
    anything is concluded.

## A — the fluvial sub-step bounded by a duration (`Deposition::fluvial_dt_max` = Δt_f)

- A step of duration d_i runs k_i = ⌈d_i / Δt_f⌉ equal fluvial sub-steps: F184's A2, with k_i in place of a uniform k. The rest of
  the step runs once (the rebound or the load, the hillslope, the marine step, the recycling, the compaction).
- **Δt_f = ∞** (the default) uses F184's `fluvial_substeps`, so Δt_f = ∞ with k = 1 is F184 bit for bit (a control).
- **The tested values**, in order: 80 kyr against 20 kyr, then 20 against 5 kyr if needed.
  - **Retained**: the largest Δt_f whose comparison with Δt_f / 4 holds on every panel world (5 seeds × intermediate and old).
  - The criterion is F174's (`f174_declared.md` A3 and Q8 / R8; `tests/f174bc_cascade.rs:2394`): the density, λ and S₁₀ within
    5 % relative, and the Kapos classes' largest move < 1 point.
  - If 20 against 5 kyr fails: the stop, with the cost and the extrapolation of the limit.
- **The live marine domain** (F184's A1) stays on.

## B — the measures

- **Tier 1** (the panel, intermediate and old):
  - the variants: without the deposition, Δt_f = 80, 20 (, 5) kyr;
  - the quantities:
    - the criterion;
    - against the same world without the deposition: the median altitude, the plains, the mountain, the large lakes;
    - the sea; the circumnavigation; the mass controls; the cost;
  - **the limit's extrapolation**: the sequence of the gaps (80 → 20 → 5), extrapolated geometrically, compared with 0.1's
    factor.
- **E1 / E2 as reported measures, split by zone**: the change of the plains and of the mountain, between the cells in net
  deposition (cover thickness > 0 at the end) and in net erosion.
- **Tier 2** (on Δt_f retained) and **tier 3** (2 048², seed 23, intermediate): as F182 and F183.

## The rules (the brief's)

- The circumnavigation is lost only against the same world without the deposition.
- **The stops**:
  - the mass control fails;
  - an iteration to convergence or a solver;
  - 20 against 5 kyr fails;
  - the analytic test fails at Δt_f retained;
  - the deposition loses the circumnavigation on ≥ 2 seeds.
- **The adoption** needs:
  - the mass control, the circumnavigation and (v);
  - the analytic test;
  - E3 held or its departure explained;
  - the author's look at the plates.
