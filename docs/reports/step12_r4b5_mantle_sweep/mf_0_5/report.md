# R4b mf_0_5 — mantle flow magnitude mf=0.5 (half intensity)

32² grid, 5 cycles × 20 steps = 100 steps total. Runtime: 206.0s.

|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path |
|---|---|---|---|---|---|---|---|---|
| INIT | 1.0000 | 0.8645 | 0.3525 | 567.6000 | 472.0000 | 95.6000 | 0.5306 | 2 |
| cycle 5 | 1.1372 | 0.8748 | 0.3516 | 566.9117 | 448.7637 | 118.1480 | 0.5641 | 8 |

Mass loss: 0.688 (0.1 %).  Cumulative macro_redistribution drift: 2.046e-12.

## Diagnostic

- final `frac S̃>0.8` = **0.352** (proxy for cratonic patches: ≥ 0.05 → patches visible, < 0.02 → flattened)
- mass loss 0.1 % over 5 cycles (≤ 5% non-conservative if > 5 %)
- max drainage path (final): **8** (≤ 3 → local, ≥ 6 → long-distance)

**Verdict subjectif** : voir le patchwork (`init_*.png` vs `final_*.png`) — la métrique scalaire ne tranche pas seule.
