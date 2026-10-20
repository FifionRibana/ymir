# Finding 185 — the fluvial sub-step bounded by a duration: Part 0 accepted (Yuan's eq. 24 reproduced exactly by the deposition scheme); step 2 of the reduced plan: 80 against 20 kyr fails F174's criterion on 4 / 4 seeds at the intermediate age, on the classes only (1.1–3.2 points; density, λ and S₁₀ within 5 %)

**Status: Part 0 measured and accepted; step 2 measured; step 3 (20 against 5 kyr on seed 23) waits for the author's go-ahead;
nothing else runs; nothing committed.** Declared: `f185_declared.md`. Predictions: `f185_predictions.md`.

## The author's reduced plan (after stopping the first tier-1 bench)

- **The reviewer's error, recorded (the author's)**: the prompt asked for 80 / 20 / 5 kyr on 10 worlds at both ages without
  costing it, far beyond the 45-minute wall-time rule per tier.
- **The plan**:
  1. Part 0 alone, reported;
  2. if the analytic test holds: 80 against 20 kyr on seeds 18 and 23 at the intermediate age (seeds 8 and 19 added, already
     computed);
  3. only if step 2 fails, and with the author's agreement: 20 against 5 kyr on seed 23;
  4. nothing else (no old age, no tier 2, no 2 048²).

## Part 0 (accepted)

- **0.1**: Yuan et al. 2019 JGR §3.1, p. 9, eq. 24: K_f p̃^m A^m S^n = U (1 + G/p̃).
  - With m = 0.5, n = 1 and G = 1.55, the slope ratio is 2.55.
  - The scheme's discrete form is 1 + G (Ñ − 1)/Ñ (eq. 13–14: the deposit uses the Ñ − 1 cells strictly upstream).
- **0.2** (`cascade::history::f169_tests::the_yuan_equilibrium`, permanent, `#[ignore]`d):

  | run | ratio | the discrete factor | gap |
  |---|---|---|---|
  | without the hillslope diffusion, 20 kyr | 2.527 | 2.527 | +0.0 % (asserted < 5 %) |
  | the full cascade, 20 kyr | 2.26 | 2.53 | −10.6 % (printed) |
  | the full cascade, 80 kyr | 2.24 | 2.53 | −11.5 % (printed; the gap falls with the step, asserted) |

  - The test exports at the base level (Yuan JGR p. 7). Its first design filled its sea strip and never equilibrated.
  - **The reviewer's error, recorded (the author's)**: the 5 % threshold was set on the full cascade, mixing the scheme's
    validation with the model's other processes (the hillslope diffusion, which Yuan neglects, p. 8).
- **The author's question, by reading the code**: the material the hillslope moves onto a channel cell enters Q_s when the next
  step's incision takes it up.
  - The hillslope runs after the fluvial stage: `cascade.rs:825`, `:837`.
  - The next incision's e = before − after: `deposition.rs:684–686`, routed through `cascade.rs:727`; Q_in + e at
    `deposition.rs:473`, `:483`, `:537`.
  - It is booked as basement, not cover: the hillslope's gains add no cover (`cascade.rs:854–858`), so the re-incision is e_b
    (`:748`, `:755`), loaded at ρ_c in A1 (`:1006`), with no porosity under A2.

## Step 2 — 80 against 20 kyr, the intermediate age, seeds 8, 18, 19, 23

| seed | density | λ | S₁₀ | the classes' largest move | median altitude | converged |
|---|---|---|---|---|---|---|
| 8 | +0.3 % | +4.0 % | +0.5 % | 3.17 points | −59 m | no |
| 18 | −1.8 % | +3.9 % | +0.8 % | 2.32 points | −30 m | no |
| 19 | −1.1 % | +2.0 % | +0.9 % | 1.65 points | −12 m | no |
| 23 | −3.3 % | −1.8 % | +0.0 % | 1.08 points | −6 m | no |

- **F174's criterion holds on 0 / 4**, failing on the classes only: density, λ and S₁₀ are within 5 % on all four (worst 3.3 %,
  4.0 %, 0.9 %). The classes move 1.08–3.17 points against the 1-point bar.
- **Against F184's uniform k**: 4 → 16 moved the classes 4.2–6.8 points at the intermediate age. 80 → 20 kyr moves them 1.1–3.2.
- **Against the same world without the deposition** (medians of the 4):

  | Δt_f | Δ median altitude | Δ plains | Δ mountain | Δ large lakes' area | E3 |
  |---|---|---|---|---|---|
  | 80 kyr | +50 m | −7.0 points | +5.7 points | −20 % | 2 / 4 |
  | 20 kyr | +32 m | −4.7 points | +5.4 points | +3 % | 1 / 4 |

  - The sea is kept on 4 / 4; no circumnavigation is lost (seed 23 is not circumnavigable without the deposition either).
  - The mass controls hold ≤ 4.8 × 10⁻¹¹.
- **E1 / E2 by zone** (20 kyr, medians): the zone in net deposition is 31 % of the land, with Δ plains / mountain −4.0 / +4.4
  points; the zone in net erosion has −0.2 / +5.7 points.
  - The mountain gains in both zones; the plains fall mostly in the deposition zone.
  - Not the « filled basins vs steepened rivers » split the brief expected: the deposition zone loses plains too.
- **The limit**: two points only (80, 20); Aitken's extrapolation needs 5 kyr.
- **The cost** (CPU per world): 80 kyr 117–140 s, 20 kyr 251–265 s. Seed 23 ran here in 425 s wall.

## What ran before the stop (reported, not concluded)

The first tier-1 bench had computed, at the intermediate age, seeds 8, 18, 19 (without, 80, 20 kyr) and seed 23 (without)
before the author stopped it. Those worlds are the ones step 2 reused.

## For the author

- **Step 3** (20 against 5 kyr on seed 23 alone): ~5 kyr costs ~4 × 20 kyr, ~17 CPU min, under the 45-minute rule. It waits for
  the go-ahead.
- Seed 23's classes moved the least (1.08 points), so seed 23 alone may pass where seeds 8 and 18 would not.
