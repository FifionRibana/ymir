# R4b test_4_tectonic_off — mantle OFF — isolate macro_redistribution from tectonic driver

32² grid, 5 cycles × 20 steps = 100 steps total. Runtime: 47.9s.

|  | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path |
|---|---|---|---|---|---|---|---|---|
| INIT | 1.0000 | 0.8645 | 0.3525 | 567.6000 | 472.0000 | 95.6000 | 0.5198 | 2 |
| cycle 5 | 0.9988 | 0.8624 | 0.3525 | 566.8936 | 470.8870 | 96.0066 | 0.5188 | 13 |

Mass loss: 0.706 (0.1 %).  Cumulative macro_redistribution drift: 1.364e-12.

## Diagnostic

- final `frac S̃>0.8` = **0.353** (proxy for cratonic patches: ≥ 0.05 → patches visible, < 0.02 → flattened)
- mass loss 0.1 % over 5 cycles (≤ 5% non-conservative if > 5 %)
- max drainage path (final): **13** (≤ 3 → local, ≥ 6 → long-distance)

**Verdict subjectif** : voir le patchwork (`init_*.png` vs `final_*.png`) — la métrique scalaire ne tranche pas seule.
