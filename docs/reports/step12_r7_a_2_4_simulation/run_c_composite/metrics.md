# R7.A.2.4 — run_c_composite

64² active_medley, workflow ON (D2 + D1-ter), mf=1.0, evo=0.10, craton_amp=3. 5 cycles × 20 steps. Runtime: 1925.2s.

## Per-cycle solver + S̃ health

| cycle | peak \|v\| | Newton C/S/D/Cap | peak S̃ | frac>0.85 | frac>0.95 | mass | max_path | CG e/m/l |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.901e-3 | 105/0/0/0 | 1.184 | 0.297 | 0.054 | 1282.09 | 7 | 7710/8263/8929 |
| 2 | 1.746e-3 | 105/0/0/0 | 1.169 | 0.316 | 0.052 | 1281.65 | 8 | 7595/8034/8911 |
| 3 | 1.613e-3 | 105/0/0/0 | 1.155 | 0.324 | 0.052 | 1281.20 | 8 | 7571/8076/8719 |
| 4 | 1.491e-3 | 105/0/0/0 | 1.142 | 0.331 | 0.050 | 1280.75 | 8 | 7831/8097/8538 |
| 5 | 1.387e-3 | 105/0/0/0 | 1.129 | 0.335 | 0.048 | 1280.29 | 8 | 7562/8101/8544 |

## Multi-dim acceptance (R4.1–R4.6)

- R4.1 Continents émergés: peak S̃_final = 1.129 > sea = 0.560 → **PASS**
- R4.2 Cratons préservés: retention = 80.5 % → **PASS**
- R4.3 Bordures + chaînes: VISUAL (inspect `cycle_5_altitude_fixed.png`)
- R4.4 Conservation: mass loss/cycle = 0.035 % → **PASS**
- R4.5 Drainage actif: max_path = 8 (cycles 1-5) → **PASS**
- R4.6 Dynamique soutenue: peak |v| > 0.1 on 0/5 → **FAIL**

Auto count: **4 / 5** (R4.3 visual pending).
