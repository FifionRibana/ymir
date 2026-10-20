# Finding 184 — the marine domain re-read at each sub-step works; the fluvial stage does not converge in sub-steps (k = 4 → 16 fails the project's criterion on 10 / 10 panel worlds): the declared stop, after tier 1

**Status: A1 (the live marine domain) and A2 (k fluvial sub-steps) are built behind their flags and their controls pass.
Tier 1 ran on the panel at the intermediate and old ages. The brief's stop « k = 16 ne converge pas » fired, so tiers 2 and 3 did
not run.** Nothing is committed. The reference is F183's combination (load, compaction, recycling) on F181's default profile.

Declared before any measurement: `f184_declared.md` (A1, A2, the convergence criterion and its source). Predictions:
`f184_predictions.md`.

## What is built (`Deposition::marine_live`, `Deposition::fluvial_substeps`; both off / 1 by default)

- **A1** (`deposition::marine_step_live`):
  - at every marine sub-step the domain is the cells with h < h_sea;
  - a source part landing at or above the sea walks down the steepest descent to the sea, filling each sea cell up to its room
    and walking on with the rest. At a closed depression it is stranded (reported);
  - every inflow into a domain cell is capped at its room (h_sea − h), and the capped part stays in the donor;
  - no iteration to a tolerance.
- **A2**: the step's uplift is split in k parts. Each sub-step adds its part, re-reads the land mask and runs the fluvial stage
  (incision, pass, lakes, the cover's update) at dt / k.
  - Once per step: the rebound or the load's flexure, the hillslope, the marine step, the recycling, the compaction, E_cum.
  - A river's outflow on a cell a later sub-step raised above the sea goes to the live marine walk, or stays in place
    (stranded) without A1. Without this, the first k = 4 control lost 14 % of the mass. Found by the control, fixed before any
    tier.
- **Controls** (permanent):
  - `the_live_marine_step_fills_to_the_sea_and_walks`: mass exact; a sea cell fills only to the sea; a land source walks away;
    a pit strands.
  - `the_f184_controls`: k = 1 with the domain fixed is F183's field; the live domain alone and k = 4 alone each keep both
    mass controls < 10⁻⁹; the live domain at k = 4 runs 4 sub-steps per step.
  - **F183 bit for bit across the trees**: the same synthetic history of F183's combination at HEAD and at 022d3cc hashes
    alike (674074df3c49c6cd).
- **The cap's diagnostic**: `Pass::capped`, the fluvial deposit made where G / Ñ ≥ 1.

## The convergence criterion (the project's, F174)

- The source: « statistical convergence (density, λ, sinuosity, classes) » (ADR, F174 at the commit of Parts B′, C′, D;
  `f174_declared.md` A3 and the Q8 / R8 rule; `tests/f174bc_cascade.rs:2394`).
- Applied between k and 4k on the same world: the drainage density, λ and S₁₀ within 5 % relative, and the Kapos classes'
  largest move < 1 point.

## Tier 1 — the panel (5 seeds) at the intermediate and old ages

### Against the same world without the deposition (medians over the panel)

| age | variant | sea kept | circumnavigation lost | marine Σ > 100 m | Δ median altitude | Δ plains | Δ mountain | Δ large lakes' area | E1 / E2 / E3 |
|---|---|---|---|---|---|---|---|---|---|
| 100 Myr | F183 (k = 1) | 5 / 5 | 0 | 2.23 % (max 14.7) | +358 m | −30.9 | +20.4 | −9 % | 0 / 0 / 3 |
| 100 Myr | + A1 | 5 / 5 | 0 | 0.39 % (max 8.1) | +351 m | −29.6 | +23.1 | +4 % | 0 / 0 / 2 |
| 100 Myr | + A1, k = 4 | 5 / 5 | 0 | 0.24 % (max 2.5) | +84 m | −10.2 | +9.0 | −22 % | 0 / 0 / 3 |
| 100 Myr | + A1, k = 16 | 5 / 5 | 0 | 0.23 % (max 2.6) | +40 m | −5.2 | +7.1 | −53 % | 1 / 0 / 4 |
| 300 Myr | F183 (k = 1) | 2 / 5 | 4 | (sea gone) | +5 120 m | −54.8 | +77.1 | −100 % | 0 / 0 / 3 |
| 300 Myr | + A1 | 2 / 5 | 4 | (sea gone) | +3 401 m | −54.8 | +70.5 | −100 % | 0 / 0 / 3 |
| 300 Myr | + A1, k = 4 | 5 / 5 | 0 | 0.18 % (max 3.3) | +106 m | −16.0 | +11.0 | −73 % | 0 / 0 / 2 |
| 300 Myr | + A1, k = 16 | 5 / 5 | 0 | 0.22 % (max 1.8) | +70 m | −10.1 | +7.7 | −52 % | 0 / 0 / 2 |

- **The mass controls** hold on every world: flux and stock ≤ 8.9 × 10⁻¹¹.
- **A1's stranded share** of the sea's solid: median 0.4–1.1 % (10.7 % at the old age with A1 alone, where the sea vanishes).
- **The fluvial deposit under the cap**: **0.0 % on every world**. With G = 1.55, a cell that receives a flux has Ñ ≥ 2, so
  G / Ñ ≤ 0.775 < 1: the cap can never fire. Only the sensitivity's G = 3.1 reaches it, at Ñ = 2 or 3.
- **Lakes born of a deposit**: median 17–26 at the intermediate age, 10–12 at the old age with k ≥ 4.
- **The thickest cover**: about 6 km at k ≥ 4 (the Ogooué: ~2 km in 5.3 Myr); 16–22 km at the old age with k = 1.

### The convergence

| age | k → 4k | density (worst) | λ (worst) | S₁₀ (worst) | the classes' largest move | holds | fluvial volume (median) | median altitude (median) |
|---|---|---|---|---|---|---|---|---|
| 100 Myr | 1 → 4 | 6.5 % | 8.3 % | 1.1 % | 12.2–23.9 points | 0 / 5 | −17.5 % | −58.3 % |
| 100 Myr | 4 → 16 | 2.4 % | 10.2 % | 2.0 % | 4.2–6.8 points | 0 / 5 | −4.5 % | −29.5 % |
| 300 Myr | 1 → 4 | 32.7 % | 172 % | 2.6 % | 14.2–79.9 points | 0 / 5 | −70.6 % | −93.6 % |
| 300 Myr | 4 → 16 | 5.6 % | 6.3 % | 0.5 % | 3.9–10.4 points | 0 / 5 | −14.5 % | −20.6 % |

- **What fails is the classes.** The drainage density and the sinuosity are within the 5 % at 4 → 16 on most worlds; λ fails on
  seed 23 (10.2 % at the intermediate age) and 4 → 16 at the old age (6.3 %). The Kapos classes move 4–10 points.
- **The direction is monotone, towards the world without the deposition.** Each refinement lowers the land and returns plains:
  - Δ median altitude: +351 → +84 → +40 m at the intermediate age, +3 401 → +106 → +70 m at the old age;
  - the plains: −29.6 → −10.2 → −5.2 points.
  - The gaps shrink by a factor of ~3–6 per ×4 of k. Extrapolated (a hypothesis, not measured), k ≈ 64 might bring the classes'
    move near 1 point.
- **The mechanism** (the reviewer's D2 hypothesis, corrected):
  - the explicit splitting is not converged, but not through the cap (which never fires);
  - within a sub-step the pass deposits G / Ñ of the incoming flux on every cell, down to the hillslopes' second cell, and
    nothing erodes that deposit before the next sub-step;
  - the longer the sub-step, the more a deposit made in one sub-step outlives the incision it would meet in a coupled
    scheme. So the land is raised (k = 1: +351 m) and the plains lost.
  - The fluvial thickness map shows deposits over almost the whole land, lighter along the channels
    (`fluvial_k1_k4_k16_age2.png`).
- **The old age's collapse at F183** was mostly D2, not D1: with A1 alone (k = 1) the old profile still loses the sea on 3 / 5
  worlds; with k ≥ 4 it keeps it on 5 / 5.

### The cost per world (CPU s, the panel)

| age | without | F183 | + A1 | k = 4 | k = 16 |
|---|---|---|---|---|---|
| 100 Myr | 21–25 | 42–44 | 56–61 | 87–99 | 206–374 |
| 300 Myr | 22–26 | 60–96 | 145–169 | 183–566 | 377–586 |

- One intermediate world's wall time was a machine suspension: 9 692 s wall for 280 s CPU.
- **k = 64** (needed if the trend holds) would cost about 4 × k = 16: ~15–40 CPU min per world.

## The predictions

- **The reviewer's**:
  - **S1** (A1: < 1 % above +100 m on every world and age): refuted on the maxima (8.1 % at k = 1; 2.5–3.3 % at k = 4–16),
    held on the medians (0.2–0.4 %).
  - **S2** (A1 alone keeps the old sea on ≥ 4 / 5): refuted (2 / 5).
  - **S3** (> 10 % between k = 1 and 16): held. The median altitude falls by > 60 % at the intermediate age and > 90 % at the
    old age; the fluvial volume by ~21 % and ~75 %.
  - **S4** (E1–E3 on the converged configuration): no configuration converged. At k = 16: E1 1 / 5 and 0 / 5, E2 0 / 5,
    E3 4 / 5 and 2 / 5.
  - **S5**: not measured (the stop).
- **Mine**:
  - T1: refuted on the maxima.
  - T2: refuted (A1 alone does not save the old age).
  - T3: refuted (the volume moves > 10 %; the cap never fires).
  - T4: refuted (4 → 16 fails).
  - T5: E1 failing, held; E3 on ≥ 3 / 5, held at the intermediate age only.
  - T6: not measured.
  - T7: held (controls ≤ 8.9 × 10⁻¹¹; k = 16 costs 4–7 × k = 1 with A1).

## The stop and the author's questions

**The stop** (written before): « k = 16 ne converge pas : arrête-toi et rapporte le coût ». It fired at both ages, 0 / 5 each.
Tiers 2 and 3 were not run.

**What is established**:
- **A1 is right and cheap**: the marine deposit above 100 m falls to a median 0.2–0.4 %. It does not by itself save the old age.
- **The fluvial splitting converges only slowly in k.** The deposition's effect on the relief shrinks from +351 m to +40 m
  (intermediate) as k goes 1 → 16, and is not yet converged.

**For the author**:
1. **The sub-step by duration rather than by count.** The history's steps are geometric (18.7 kyr to 1.47 Myr at the
   intermediate age, to 5.7 Myr at the old age), so a uniform k refines the many short steps needlessly.
   - A fluvial sub-step bounded by a duration Δt_f (k_i = ⌈d_i / Δt_f⌉), tested against Δt_f / 4, would refine where it fails.
   - Its cost is the history's age over Δt_f sub-steps; at 100 Myr, Δt_f = 20 kyr gives ~5 000, about k = 16's.
   - A choice of Δt_f needs no source: it is a numerical resolution judged by the criterion.
2. **Or k = 64**, at ~15–40 CPU min per world (before tier 2's 25 worlds).
3. **The physics, beyond the numerics**: even at k = 16 the deposition raises the land (+40–70 m) and turns plains into hills
   (E1, E2 fail). Whether E1–E2 hold once converged is not known yet.

**Not changed, waiting for the decision**: the flags stay off by default, production unchanged; `tests/f184_cascade.rs` resumes
per world.

## The checks (after; the two flags off by default)

- `cargo check --workspace`: clean.
- 669 lib tests pass (F183's 667 plus the two F184 controls); 33 viz tests pass.
- The guard reads Match on all six states (field and lakes).
- Production's C1 hash is unchanged (`checks_f184.txt`).

**At commit** -- the author's words (verbatim): « Feu vert pour commiter le F184 une fois tes vérifications passées (A1 et A2
derrière leurs drapeaux, coupés par défaut ; les contrôles ; le rapport avec l'arrêt). Garde target_f184/. Ajoute à l'ADR deux
erreurs du relecteur : 1. D2 : le mécanisme supposé (le plafond G/Ñ ≥ 1) était faux, puisque le plafond ne se déclenche jamais.
La cause est le dépôt de chaque sous-pas jusque sur les versants, repris seulement au sous-pas suivant (ton diagnostic). 2. E1 et
E2 (« le dépôt aplanit ») ont été écrites sans vérification contre l'équation de Yuan. À l'équilibre, le dépôt raidit les rivières
des zones soulevées : hypothèse du relecteur, K A^m S^n = U (1 + G/p̃), à vérifier au F185. E1 et E2 deviennent des mesures
rapportées, plus des conditions. »
- **The reviewer's errors, recorded (the author's)**:
  1. **D2**: the supposed mechanism (the cap G / Ñ ≥ 1) was wrong: the cap never fires. The cause is each sub-step's deposit,
     reaching down to the hillslopes, which is taken up again only at the next sub-step.
  2. **E1 and E2** (« the deposition flattens ») were written without checking them against Yuan's equation.
     - At equilibrium the deposition steepens the rivers of the uplifted areas. The reviewer's hypothesis, **to verify at
       F185**: K A^m S^n = U (1 + G / p̃).
     - **E1 and E2 become reported measures, no longer conditions.**
- **Committed as is**:
  - A1 (`marine_live`) and A2 (`fluvial_substeps`) behind their flags, off by default;
  - the controls;
  - the report with the stop.
  - `target_f184/` is kept.
- Production is unchanged: the C1 hash, the guard 6/6 Match, 669 lib and 33 viz tests.
