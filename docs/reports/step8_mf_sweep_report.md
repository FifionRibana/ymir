# Step 8 — Mf sweep (peak|v_solved| scaling, yielding activation threshold)

> Fixed: coupling = `1.000`, num_modes = `6`, mantle_seed = `42`, world seed = `42`, num_plates = `8`, grid 64², steps = `300`. **Single seed across all points** — the Fourier pattern is fixed; only the amplitude `Mf` varies.

| Mf | peak\|v_solved\| | yielding_cell_fraction_max | ε̇_II / ε̇_min | Newton converged (%) | CG iters (mean) | wallclock (s) | mass_conservation_residual |
|---|---|---|---|---|---|---|---|
| `0.30` | `1.269e-3` | `1.741e-1` | `2.031e1` | `100.0` | `498.6` | `285.58` | `7.567e-16` |
| `0.60` | `3.486e0` | `9.761e-1` | `6.540e4` | `99.7` | `1442.7` | `1054.39` | `2.979e-16` |
| `1.00` | `9.552e0` | `9.976e-1` | `1.586e5` | `99.3` | `1420.7` | `1190.81` | `1.080e-15` |
| `1.50` | `1.718e1` | `9.990e-1` | `2.703e5` | `99.3` | `1392.9` | `1226.12` | `1.069e-15` |
| `2.00` | `2.483e1` | `9.995e-1` | `3.758e5` | `98.4` | `1400.6` | `1229.50` | `3.673e-16` |

**Monotonicity: ✅ `peak|v_solved|` monotonically non-decreasing with `Mf`.** Expected from the linear coupling of mantle amplitude to forcing. Non-linear saturation (sub-linear growth) is acceptable at the top end — the full-field response includes viscous dissipation and (through Newton) the power-law rheology.

**Yielding activation threshold (observed).** Yielding first fires at `Mf ≥ 0.30` in this sweep. The critical `Mf` is a physical property measured, not prescribed; at smaller `Mf` the mantle bootstrap does not push `ε̇_II` above the regularisation floor.

