# R4b test_1_high_rebound — rebound 0.95 (weaker redistribution)

32² grid, 5 cycles × 20 steps = 100 steps total. Runtime: 2502.5s.

|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path |
|---|---|---|---|---|---|---|---|---|
| INIT | 1.0000 | 1.0000 | 0.3525 | 567.6000 | 361.0000 | 206.6000 | 0.8685 | 1 |
| cycle 5 | 2.5494 | 1.3936 | 0.0078 | 484.3090 | 11.1488 | 473.1602 | 1.0504 | 1 |

Mass loss: 83.291 (14.7 %).  Cumulative macro_redistribution drift: 1.023e-12.

## Diagnostic

- final `frac S̃>0.8` = **0.008** (proxy for cratonic patches: ≥ 0.05 → patches visible, < 0.02 → flattened)
- mass loss 14.7 % over 5 cycles (**>5%** non-conservative if > 5 %)
- max drainage path (final): **1** (≤ 3 → local, ≥ 6 → long-distance)

**Verdict subjectif** : voir le patchwork (`init_*.png` vs `final_*.png`) — la métrique scalaire ne tranche pas seule.
