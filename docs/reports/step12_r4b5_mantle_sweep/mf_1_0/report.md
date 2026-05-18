# R4b mf_1_0 — mantle flow magnitude mf=1.0 (default baseline, macro defaults)

32² grid, 5 cycles × 20 steps = 100 steps total. Runtime: 2421.9s.

|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path |
|---|---|---|---|---|---|---|---|---|
| INIT | 1.0000 | 1.0000 | 0.3525 | 567.6000 | 361.0000 | 206.6000 | 0.8685 | 1 |
| cycle 5 | 1.9314 | 0.9920 | 0.0469 | 498.1198 | 47.6158 | 450.5041 | 0.8133 | 1 |

Mass loss: 69.480 (12.2 %).  Cumulative macro_redistribution drift: 9.095e-13.

## Diagnostic

- final `frac S̃>0.8` = **0.047** (proxy for cratonic patches: ≥ 0.05 → patches visible, < 0.02 → flattened)
- mass loss 12.2 % over 5 cycles (**>5%** non-conservative if > 5 %)
- max drainage path (final): **1** (≤ 3 → local, ≥ 6 → long-distance)

**Verdict subjectif** : voir le patchwork (`init_*.png` vs `final_*.png`) — la métrique scalaire ne tranche pas seule.
