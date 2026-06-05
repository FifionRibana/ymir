# Issue #145 — Stage 5a: re-validation matrix (what changes under rigid)

Method: flipped **all** test configs `rigid_continental_crust: false → true` (22 configs / 16 files — enumerate, not sample), ran the full `ymir-core` suite `--no-fail-fast`, enumerated every change, reverted. This maps exactly what the production flip (5d) would affect.

## Result — only 3 tests change, ALL advection-imprint-dependent

| test | phase | why it changes |
|---|---|---|
| `davis_suppe_wedge_body_invariants` | 1.2 | wedge imprint depends on continental crust advecting into/through the wedge; rigid freezes it |
| `davis_suppe_imprint_preserved_with_equilibrium` | 1.3 | same — wedge-body imprint formed differently under rigid |
| `erosion_preserves_davis_suppe_imprint_partially` (`wedge_p95`) | 1.4 | already `#[ignore]`'d (step 1); wedge_p95 lift was floor-injection-dependent |

All three are the **Davis-Suppe wedge imprint** — which *legitimately* differs when continental crust is rigid (the orogenic wedge is no longer fed by advected continental material). This is the expected, correct consequence of the fix, not a regression.

## PASS under rigid (the safe surface — exhaustive)

- **Track D acceptance** — `boundary_events_fire_correctly` (subduction > 0, accretion > 0), `acceptance_phase_2_gate_seed_diversity`. Events still fire under rigid.
- **9th bit-identical** (`ninth_bit_identical_preservation_phase_2_r7`) — PASSES. This is a *determinism* contract (run-twice-identical, NOT legacy-byte-match), so it holds under rigid → **determinism preserved** (partial 5c confirmation).
- **Mass-conservation budget** (`c1_phase_2_track_d_mass_conservation`, 1e-6) — PASSES. The erosion clean-removal fix is tracked correctly by the budget; rigid advection still conserves (face-flux cancellation).
- **Bathymetry** (`bathymetry_modulated_by_age_after_300_steps`, `disabled_matches_phase_1_4`, `downstream_pipeline_accepts_phase_2_altitude`) — PASS. Confirms the production-render finding: Stein-Stein (age-derived depth) is unaffected by the S̃ swing.
- **Track B acceptance**, **init_r7**, **isostasy**, **boundary_evolution**, **Phase 1.1 advection** — all PASS.
- (`rectangular_simulation_smoke_test` — pre-existing v2-Stokes, unrelated.)

## Verdict

**The flip is SAFE.** The entire re-validation surface is the **3 Davis-Suppe wedge-imprint tests**, all of which change for the same legitimate reason (rigid continental crust no longer advects into the wedge). Track D, determinism, mass budget, and bathymetry are all robust to the flip.

## Next

- **5b**: re-baseline the 3 imprint tests WITH understanding — measure the new wedge imprint under rigid, confirm it is physically sensible (wedge still forms via orogenic source on rigid crust), document the new baselines (not blind-rebaseline).
- **5c**: confirm determinism explicitly on the rigid production path (run ×2 byte-identical). (9th bit-identical already passes under rigid — strong partial confirmation.)
- **5d**: flip production (`RunBaseline`, `phase_a_c1`) to `true`; decide flip-default vs remove-flag.
