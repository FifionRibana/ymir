# Step 7 — Sp sweep (peak|v| monotonicity)

> Fixed: τ_slab = `0.500`, k_slab_accum = `1.000`, ε = `1.0e-6`, seed = `42`, num_plates = `8`, grid 64², steps = `300`.

| Sp | peak\|v\| | m_subducted_max | yielding_cell_fraction_max | Newton converged (%) | CG iters (mean) | wallclock (s) |
|---|---|---|---|---|---|---|
| `0.50` | `3.602e-5` | `1.007e-4` | `0.000e0` | `100.0` | `129.8` | `36.21` |
| `1.00` | `3.602e-5` | `1.007e-4` | `0.000e0` | `100.0` | `130.1` | `34.38` |
| `1.50` | `3.602e-5` | `1.007e-4` | `0.000e0` | `100.0` | `130.4` | `35.84` |
| `2.00` | `3.602e-5` | `1.007e-4` | `0.000e0` | `100.0` | `130.7` | `34.93` |
| `3.00` | `3.602e-5` | `1.008e-4` | `0.000e0` | `100.0` | `131.3` | `35.06` |

**Interpretation — flat across the Sp band (bootstrap failure regime).**

`peak|v|` is identical across `Sp ∈ [0.5, 3.0]` to f64 precision. This is the signature of the closed-loop gain `G = Sp · k_slab_accum · τ_slab / (η · L)` sitting `≪ 1` everywhere in the §4.8 target band, with the floor-dominated `η_newton ≈ 100`. The quiescent fixed point is linearly stable; the system remains at the Step 6 baseline regardless of `Sp`. `peak|f_slab|` does scale linearly with `Sp` (visible in the physics report), but the Stokes inversion `v = f · L²/η` damps it by `1/η ≈ 0.01`, so no measurable `peak|v|` response. Monotonicity is trivially satisfied (zero difference).

This is consistent with the amplifier-vs-initiator revision documented in `step7_physics_report.md §Yielding checkpoint`. A non-flat sweep is expected once Step 8 (mantle forcing) imposes `v_mantle` externally and breaks the floor-dominated regime; slab-pull will then amplify visibly.

