# Finding 185 — predictions, written before any measurement

The reference is F184's live configuration (F183's combination with the live marine domain) on F181's default profile, at 256²,
the panel (8, 18, 19, 23, 29) at the intermediate and old ages; the circumnavigation against the same world without the
deposition.

## The reviewer's (hypotheses; at least one is false)

- **V1**: 0.1 confirms the form K A^m S^n = U (1 + G/p̃), up to notation.
- **V2**: 20 against 5 kyr holds on the 10 panel worlds; 80 against 20 kyr fails on ≥ 1.
- **V3**: the analytic test holds within 5 % at Δt_f retained.
- **V4**: at convergence the median altitude gains +10 to +60 m on the world without the deposition (intermediate), and the
  plains lose < 5 points.
- **V5**: E3 holds on ≥ 4 of the 5 panel worlds at both ages.
- **V6**: the deposition loses the circumnavigation on none of the 15 seeds at the intermediate age.

## Mine

- **W1 (0.1)**: confirmed (eq. 24, p. 9; V1 held). The scheme's discrete factor is (1 + G (Ñ − 1)/Ñ).
- **W2 (the analytic test)**: it misses the 5 % at first, because the cascade's hillslope diffusion runs on the channels and the
  lateral erosion adds an erosion term Yuan neglects. Expected order: the measured ratio is 2.0–2.4 against the discrete 2.4–2.5,
  and it approaches the factor when the hillslope is switched off in the test. The gap at 4 × Δt_f is larger than at Δt_f.
- **W3 (the convergence)**:
  - 80 against 20 kyr fails on ≥ 3 of the 10 worlds (the classes);
  - 20 against 5 kyr holds on ≥ 8 of the 10 but fails on at least one old world (seed 23 or 29, where k = 16 moved the classes
    by 10 points);
  - so V2 is refuted and the stop risks firing.
  - The steps of F184's k = 16 lasted 92 kyr at most at the intermediate age (1.47 Myr / 16) and 356 kyr at the old age, so
    20 kyr refines the old age far more than k = 16 did.
- **W4 (at convergence)**: the median altitude gains +20 to +45 m (intermediate), the plains lose 3–6 points (V4 at its edge).
  Split by zone, the plains are gained in the zones of net deposition (the filled basins) and lost in the zones of net erosion
  (the steepened rivers raise the hill classes).
- **W5 (E3)**: holds on ≥ 4 / 5 at the intermediate age, and on 3 / 5 at the old age (V5 refuted at the old age).
- **W6 (the cost)**: 20 kyr costs ~5 min CPU per intermediate world and ~15 min per old world; 5 kyr four times that.
