# Step 12 R4 metrics — active_medley

64² grid, 5 cycles × 20 steps (mantle ON), Phase A defaults: α=0.01, isostatic_rebound_ratio=0.80, max_drainage_distance=10.

| cycle | peak S̃ | mean S̃ cont. | frac S̃>0.8 | mass total | mass cont. | mass ocean | sea_lvl | max_path | total_eroded | mass_drift |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1.0000 | 0.0000 | 0.0000 | 2280.4000 | 0.0000 | 2280.4000 | 1.0539 | 0 | 0.0000e0 | 0.000e0 |
| 1 | 2.5550 | 1.2651 | 0.0369 | 2256.6031 | 191.0322 | 2065.5710 | 1.0539 | 3 | 1.5017e-1 | -2.728e-12 |
| 2 | 2.2048 | 1.0864 | 0.1001 | 2213.6982 | 445.4286 | 1768.2697 | 0.9126 | 8 | 2.3906e-1 | 4.547e-13 |
| 3 | 2.4278 | 1.2298 | 0.0337 | 2180.9274 | 169.7185 | 2011.2089 | 1.0027 | 4 | 1.5586e-1 | 1.364e-12 |
| 4 | 2.1698 | 1.0809 | 0.0408 | 2135.8303 | 180.5081 | 1955.3222 | 0.8994 | 4 | 1.4025e-1 | -4.547e-13 |
| 5 | 2.1272 | 1.0567 | 0.0383 | 2089.2211 | 165.9000 | 1923.3210 | 0.8820 | 5 | 1.2361e-1 | -9.095e-13 |

## Verdict R4.1–R4.5

- **R4.1 — Continents émergés step 100** : peak S̃ final = 2.127 vs sea_level = 0.882 → **PASS**
- **R4.2 — Cratons préservés (S̃ > 0.8 retained > 50 %)** : init = 0.000, final = 0.038, retention = 0.0 % → **FAIL — too few cratons retained**
- **R4.3 — Bordures irrégulières (subjective)** : final max_path_length = 5 (proxy ≥ 2 → drainage spans coast cells) → **PASS (proxy)**
- **R4.4 — Conservation totale (drift < 1e-9 · mass)** : cumulative drift = 5.912e-12, budget = 2.280e-6 → **PASS**
- **R4.5 — Drainage actif (max_path > 1)** : max across cycles 1-5 = 8 → **PASS**
