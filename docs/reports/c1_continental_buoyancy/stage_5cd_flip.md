# Issue #145 — Stage 5c (determinism) + 5d (production flip)

## 5c — determinism on the rigid production path

`run_with_closures` rigid, run ×2, byte-identical on S̃ + age + plate_type (seeds 42, 2026). PASS. Combined with the 9th-bit-identical contract passing under rigid (5a), the rigid production path is deterministic — the C1 invariant holds (the rigidity mask + no-flux + ≥2-promotion are all serial, deterministic).

## 5d — production flip

The buoyancy fix is now LIVE in production: `rigid_continental_crust: true` at the three production sites —
- `crates/ymir-viz/src/bridge/c1/thread.rs` RunBaseline (gallery) + workflow worker,
- `crates/ymir-core/src/tectonics_v2/workflow/phase_a_c1.rs` (workflow wrapper).

**Flag kept (not removed)**, documented as the transitional regression-guard: `false` = the legacy (continents-collapse) behaviour, retained as the byte-identical A/B reference until the flag is removed (rigidity unconditional). Most tests keep `false` (legacy regression guard); the 3 imprint tests + production use `true`.

Suite green (only the pre-existing `rectangular_simulation` v2-Stokes failure); core + viz build. Production render = the validated clean state (final-visual + nb2 post-finger-fix: credible compact continental masses multi-seed, no false-land finger).

## #145 implementation complete

Cause proven (continental buoyancy) → fix (rigid crust) → erosion clean-removal (Step 1) → DS skip (justified) → equilibrium confirmed regulator → re-validation matrix (5a) → no-flux boundary + ≥2-promotion finger fix → bounded-oscillation verdict → imprint re-baseline (5b) → determinism (5c) → production flip (5d). All measured, not assumed.

## Registered follow-ups (NOT #145)

1. **Rigid-boundary refinement** — the bounded mesh oscillation ("curtain") at the sharp continental/oceanic contrast (upwind on a net contrast, capped by equilibrium; cosmetic, not false land). Possibly soften via a local boundary buoyancy transition.
2. **Oceanic ridge crust-creation** (seafloor spreading) — the raw-field deep oceans; fidelity (invisible in production via Stein-Stein).
3. **Conservative erosion deposition** (sediment routing → plains) — richness, builds on the clean base.
4. **Init: seeds with zero cratonic cells** (R7 clustering stochastic) — note.
