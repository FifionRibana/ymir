# Tracking issue — `rectangular_simulation_smoke_test` fails (Picard non-convergence)

Status: OPEN, do not fix yet. Filed from the erosion/terrace campaign (see
[ADR 0001](../adr/0001-erosion-coastal-sink-and-terraces.md)).

## Symptom
`cargo test -p ymir-core --test rectangular_simulation` fails, both release and
debug, on the **clean tree** (verified with the erosion campaign's changes stashed):

```
run_tectonics failed on rectangular grid: Some(NonlinearSolverDidNotConverge { step: 3 })
```

The failure is in the **tectonics Picard nonlinear solver**
(`crates/ymir-core/tests/rectangular_simulation.rs`, seed 42, 128×85 grid, 50 steps).
It is independent of the erosion/sink work — that test imports only
`ymir_core::tectonics::*`, no erosion symbol.

## Why it matters (three points)
1. **Permanent noise in the suite.** A permanently-red test eventually masks a real
   regression — nobody notices one more red among the known red. It should be made
   green or removed, not left indefinitely.
2. **Most upstream stage of the whole chain.** The Picard solver sits in the same
   tectonic layer that produces the equilibrium-height closures — the very thing the
   terrace chantier (ADR 0001, Finding 4) is about to touch. A solver failing to
   converge could plausibly be related to the flat equilibrium plateaux (a
   degenerate/near-singular state). **Re-check this during the terrace job.**
3. **Rectangular vs square.** The test name and setup use a RECTANGULAR (128×85)
   domain, while the whole production pipeline assumes a SQUARE grid. This may be a
   dead path exercising a configuration nobody ships — in which case the fix is to
   **REMOVE the test, not repair the solver.** **Check this hypothesis FIRST**; it
   may close the issue cheaply.

## Do first
Confirm whether any production path uses a non-square tectonic grid. If not, remove
`rectangular_simulation.rs` (and note it in the commit). If yes, the non-convergence
is a real solver bug and needs the tectonics owner.

## To file on GitHub (gh not installed in this environment)
```
gh issue create --repo FifionRibana/ymir \
  --title "rectangular_simulation_smoke_test: Picard NonlinearSolverDidNotConverge (step 3)" \
  --body-file docs/issues/picard-nonconvergence-rectangular-smoke.md
```
