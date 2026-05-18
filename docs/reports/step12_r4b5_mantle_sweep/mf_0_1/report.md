# R4b mf_0_1 — mantle flow magnitude mf=0.1 (very low intensity)

32² grid, 5 cycles × 20 steps = 100 steps total. Runtime: 76.7s.

|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path |
|---|---|---|---|---|---|---|---|---|
| INIT | 1.0000 | 0.8645 | 0.3525 | 567.6000 | 472.0000 | 95.6000 | 0.5199 | 2 |
| cycle 5 | 1.0024 | 0.8630 | 0.3525 | 567.4829 | 471.2046 | 96.2783 | 0.5198 | 12 |

Mass loss: 0.117 (0.0 %).  Cumulative macro_redistribution drift: 3.638e-12.

## Diagnostic

- final `frac S̃>0.8` = **0.353** (proxy for cratonic patches: ≥ 0.05 → patches visible, < 0.02 → flattened)
- mass loss 0.0 % over 5 cycles (≤ 5% non-conservative if > 5 %)
- max drainage path (final): **12** (≤ 3 → local, ≥ 6 → long-distance)

**Verdict subjectif** : voir le patchwork (`init_*.png` vs `final_*.png`) — la métrique scalaire ne tranche pas seule.
