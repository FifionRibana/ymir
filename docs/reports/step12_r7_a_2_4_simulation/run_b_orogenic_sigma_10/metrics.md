# R7.A.2.4 — run_b_orogenic_sigma_10

64² active_medley, workflow ON (D2 + D1-ter), mf=1.0, evo=0.10, craton_amp=3. 5 cycles × 20 steps. Runtime: 24127.7s.

## Per-cycle solver + S̃ health

| cycle | peak \|v\| | Newton C/S/D/Cap | peak S̃ | frac>0.85 | frac>0.95 | mass | max_path | CG e/m/l |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.534e0 | 36/66/0/3 | 1.967 | 0.482 | 0.266 | 1987.03 | 15 | 7958/9413/8667 |
| 2 | 8.315e-1 | 101/1/0/0 | 1.797 | 0.442 | 0.217 | 1987.93 | 15 | 15661/14350/13807 |
| 3 | 1.185e-1 | 101/1/0/0 | 2.016 | 0.294 | 0.156 | 1972.63 | 12 | 14873/16467/15231 |
| 4 | 9.797e-2 | 105/0/0/0 | 1.974 | 0.311 | 0.164 | 1960.09 | 12 | 14299/16792/15778 |
| 5 | 4.417e-2 | 105/0/0/0 | 2.079 | 0.226 | 0.115 | 1951.40 | 8 | 13207/15925/14879 |

## Multi-dim acceptance (R4.1–R4.6)

- R4.1 Continents émergés: peak S̃_final = 2.079 > sea = 0.862 → **PASS**
- R4.2 Cratons préservés: retention = 25.7 % → **FAIL**
- R4.3 Bordures + chaînes: VISUAL (inspect `cycle_5_altitude_fixed.png`)
- R4.4 Conservation: mass loss/cycle = 0.587 % → **PASS**
- R4.5 Drainage actif: max_path = 15 (cycles 1-5) → **PASS**
- R4.6 Dynamique soutenue: peak |v| > 0.1 on 3/5 → **PASS**

Auto count: **4 / 5** (R4.3 visual pending).
