# Finding 188 — the predictions, written before any measure

## The reviewer's (hypotheses; at least one is false)

- **Z1**: with A, the area of seed 8's large east lake falls by at least half.
- **Z2**: with A, seed 19's marsh look stays (the depressions and the islets remain).
- **Z3**: the convergence test A3 holds between 80 and 20 kyr.
- **Z4**: with B, the marine piles above 100 m fall by at least half on seed 23.
- **Z5**: without correcting the rivers, the sinuosity falls on the filled zones (straighter rivers).

## Mine

- **W1 (Z1)**: holds. The east lake sits on land that sinks under the load (F187: its outlet 3 km from the coast). A deposits
  G/Ñ of the flux there, enough to fill a closed depression, which by construction receives its catchment's whole flux.
  Expected: −50 to −90 %.
- **W2 (Z2)**: holds partly. The 256² depressions fill (A deposits on sinking land), so their count falls below 10. The 2 048²
  islets stay, being the basement's amplification (F187).
- **W3 (Z3)**: fails. F184 and F185 found the explicit fluvial deposition not converged in sub-steps (F185: 1.1–3.2 points on
  the classes between 80 and 20 kyr), and A is that same explicit pass, restricted. Expected: the classes move > 1 point on
  the sinking zone.
- **W4 (Z4)**: holds only if the piles are uplifted marine fill eroding at the basement's K. If they are stranded by the walk,
  B does nothing to their origin, but still erodes them faster. Expected: the stranding's share < 20 % and a fall ≥ 50 %.
- **W5 (Z5)**: holds on the filled zones at 2 048². The filled surfaces are bilinear (`bury_smooth`), so the D8 network on
  them follows the interpolated gradient: S₁₀ falls by 0.02–0.05 and `a_dir` rises.
- **W6 (the cost)**: A adds < 5 s per world at the physics level (one more pass reading a per-cell G). B adds nothing measurable.
  The 20 kyr run of A3 costs ~4 × the 80 kyr run.
- **W7 (the mass)**: holds for A, B and A + B (< 10⁻¹⁰).
