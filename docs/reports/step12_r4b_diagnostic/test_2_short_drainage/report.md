# R4b test_2_short_drainage — max_drainage_distance 3 (local drainage)

32² grid, 5 cycles × 20 steps = 100 steps total. Runtime: 2502.5s.

|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path |
|---|---|---|---|---|---|---|---|---|
| INIT | 1.0000 | 1.0000 | 0.3525 | 567.6000 | 361.0000 | 206.6000 | 0.8685 | 1 |
| cycle 5 | 2.5443 | 1.3934 | 0.0078 | 484.2138 | 11.1472 | 473.0666 | 1.0502 | 1 |

Mass loss: 83.386 (14.7 %).  Cumulative macro_redistribution drift: 6.821e-13.

## Diagnostic

- final `frac S̃>0.8` = **0.008** (proxy for cratonic patches: ≥ 0.05 → patches visible, < 0.02 → flattened)
- mass loss 14.7 % over 5 cycles (**>5%** non-conservative if > 5 %)
- max drainage path (final): **1** (≤ 3 → local, ≥ 6 → long-distance)

**Verdict subjectif** : voir le patchwork (`init_*.png` vs `final_*.png`) — la métrique scalaire ne tranche pas seule.
