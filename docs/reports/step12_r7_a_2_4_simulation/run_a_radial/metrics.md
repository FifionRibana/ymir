# R7.A.2.4 — run_a_radial

64² active_medley, workflow ON (D2 + D1-ter), mf=1.0, evo=0.10, craton_amp=3. 5 cycles × 20 steps. Runtime: 1742.2s.

## Per-cycle solver + S̃ health

| cycle | peak \|v\| | Newton C/S/D/Cap | peak S̃ | frac>0.85 | frac>0.95 | mass | max_path | CG e/m/l |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.849e-3 | 105/0/0/0 | 0.950 | 0.377 | 0.054 | 1277.93 | 9 | 7723/8124/8885 |
| 2 | 1.702e-3 | 105/0/0/0 | 0.950 | 0.380 | 0.051 | 1277.49 | 9 | 7931/8104/8866 |
| 3 | 1.577e-3 | 105/0/0/0 | 0.950 | 0.380 | 0.051 | 1277.05 | 9 | 7494/7947/8847 |
| 4 | 1.459e-3 | 105/0/0/0 | 0.950 | 0.380 | 0.049 | 1276.61 | 9 | 7588/8070/8496 |
| 5 | 1.362e-3 | 105/0/0/0 | 0.949 | 0.379 | 0.047 | 1276.16 | 9 | 7848/8117/8361 |

## Multi-dim acceptance (R4.1–R4.6)

- R4.1 Continents émergés: peak S̃_final = 0.949 > sea = 0.488 → **PASS**
- R4.2 Cratons préservés: retention = 78.4 % → **PASS**
- R4.3 Bordures + chaînes: VISUAL (inspect `cycle_5_altitude_fixed.png`)
- R4.4 Conservation: mass loss/cycle = 0.034 % → **PASS**
- R4.5 Drainage actif: max_path = 9 (cycles 1-5) → **PASS**
- R4.6 Dynamique soutenue: peak |v| > 0.1 on 0/5 → **FAIL**

Auto count: **4 / 5** (R4.3 visual pending).
