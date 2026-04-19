# Slab-Pull Reformulation — Phase 2 Implementation Results

Issue #75 — *Reformulate slab-pull as an auto-regulated operator term instead of RHS forcing.*
Branch: `75-reformulate-slab-pull-...` (same branch as Phase 1 and 1-bis). Date: 2026-04-19.

This report documents the outcome of the Phase 2 refactor:

- The per-plate velocity boost `apply_slab_pull` is no longer called in
  the tectonic pipeline.
- A cell-local operator term `γ_slab(x, y) · (v·n̂) · n̂` is added to
  the Stokes operator, seeded from `|source_rate|` on the subducting
  side of convergent margins and spread inward on the subducting plate
  with an exponential Benioff decay (`L = 3` cells).
- `γ_slab`, `n̂_x`, `n̂_y` live on `BoundaryField` and are recomputed
  once per macro step in `compute_boundary_sources`.

**Top-line result.** Scenario B (everything on) wallclock falls
~**7×** — from 336 s (Phase 1-bis baseline scaled to 120 steps) to
**48.4 s** — matching scenario C (everything on except slab-pull) at
47.6 s. The η-contrast cascade is flat: `eta_ratio` drops from 62 to
**11.23** in B, identical to C. All three new unit tests
(`operator_with_slab_pull_is_symmetric`, `slab_pull_term_is_semi_definite`,
`slab_pull_damps_convergent_motion`) pass.

---

## 1. Empirical validation — wallclock, η, residual

Config: 64² grid, seed 42, Newton solver, adaptive dt (`dt_target = 2.0`),
1 rep × 120 macro steps. Steady-regime means are computed over steps 40–80.

### 1.1 Wallclock

| Scenario | mean (s) | Comment |
|----------|---------:|---------|
| A (bare) | 7.4 | sanity floor |
| **B (full, Phase 2)** | **48.4** | ← refactor target |
| B (full, Phase 1-bis baseline) | ~336 (= 841 × 120/300) | extrapolation from §1 of Phase 1-bis |
| C (no slab-pull) | 47.6 | reference for "what B should look like" |

**B/C ratio: 1.02.** Slab-pull is no longer a cost driver. Speed-up vs
Phase 1-bis baseline: **6.95×.**

### 1.2 η contrast (`eta_ratio`, mean over steps 40–80, 3 Newton iters per step)

| Scenario | Phase 1-bis | Phase 2 | Change |
|----------|------------:|--------:|-------:|
| A | 2.5 | 2.84 | — |
| B | 61.8 | **11.23** | **×0.18 (5.5× reduction)** |
| C | 11.3 | 11.23 | unchanged |

η in B now matches C **exactly** (11.2294 vs 11.2296). The velocity
boost drove strain rates up enough to trigger plastic-branch
domination; without it, the η distribution collapses to C's regime.

### 1.3 Residual localization

| Scenario | Phase 1-bis | Phase 2 |
|----------|------------:|--------:|
| B | 0.277 | **0.145** |
| C | 0.478 | 0.146 |

In Phase 2, B and C are again statistically identical. The
interpretation from Phase 1-bis (C higher than B because its flatter η
lets residual concentrate at margins) still holds — and now B has
joined C in the flat-η regime.

### 1.4 RHS spike (`gpe_spike_p95`)

| Scenario | Phase 1-bis | Phase 2 |
|----------|------------:|--------:|
| B | 3.64 | **1.92** |
| C | 2.18 | 1.91 |

Both B and C drop toward 1.9. The Phase 1-bis claim that GPE is
the dominant RHS spike (§3) still holds — the numeric reduction is
because the steady-state thickness field is smoother without the
velocity boost driving rapid advection.

### 1.5 T_plates norm (`tp_norm`)

| Scenario | Phase 1-bis | Phase 2 |
|----------|------------:|--------:|
| B | 312.7 | **94.4** |
| C | 92.7 | 94.4 |

The velocity boost doubled `tp_norm` in the Phase 1-bis baseline. It
is gone in Phase 2. B and C are identical up to rounding.

### 1.6 BiCGSTAB iterations

Inferable from `t_solve_us`: 

| Scenario | Phase 1-bis (μs/substep) | Phase 2 (μs/substep) |
|----------|---:|---:|
| B | 1,372,010 | **343,041** |
| C | 630,639 | 361,580 |

B's per-substep solve time falls from 1.37 ms to 0.34 ms — a **4.0×
reduction**, consistent with the wallclock speed-up once sub-step
multiplicity is factored in.

---

## 2. Code changes

Seven files touched, ~350 lines of production code change (excluding
tests). Diff summary:

| File | Change | Line count |
|------|--------|-----------:|
| [boundaries.rs](crates/ymir-core/src/tectonics/boundaries.rs) | extend `BoundaryField` + `BoundaryConfig`; populate γ_slab, n̂; add `spread_gamma_benioff` BFS | ~140 |
| [stokes.rs](crates/ymir-core/src/tectonics/solver/stokes.rs) | add `SlabPullField` struct; extend `apply_stokes`/`compute_jacobi_precond`/`StencilCoeffs::compute` to include `γ·n̂⊗n̂` | ~70 |
| [newton.rs](crates/ymir-core/src/tectonics/solver/newton.rs) | plumb `slab: Option<&SlabPullField>` through `solve_velocity_newton` and `compute_nonlinear_residual` | ~20 |
| [picard.rs](crates/ymir-core/src/tectonics/solver/picard.rs) | plumb `slab` through `solve_velocity_picard` | ~10 |
| [tectonics.rs](crates/ymir-core/src/tectonics/solver/tectonics.rs) | remove `apply_slab_pull` call; split-borrow `boundary_field` for the solve; extend `solve_velocity_direct`/`solve_with_continuation` | ~25 |
| [linear_solve.rs](crates/ymir-core/src/tectonics/solver/linear_solve.rs) | test-only: pass `None` to `apply_stokes`/`compute_jacobi_precond` | ~15 |

### 2.1 Design choices (from Phase 1 §5)

- **γ_slab seed:** `γ_seed = slab_pull_factor · |source_rate|` on
  cells with `source_rate < 0` and `boundary_type ∈ {Subduction,
  OceanicSubduction}`. No new config knob — reuses
  `slab_pull_factor` whose semantics shift from "velocity boost per
  unit mass" (pre-#75) to "γ coefficient per unit convergence rate".
- **Benioff decay:** exponential `exp(-d/L)` with `d` = shortest
  path-on-plate in cells; `L = benioff_decay_cells` (new
  `BoundaryConfig` field, default `3.0`). Implemented as per-seed BFS
  masked to the seed's plate; each cell takes the max over all seeds
  that reach it. BFS cost is O(seeds × (3L)²) — ~60 k ops / step at
  64² / ~10 seed cells, well under 1 ms.
- **n̂ sign convention:** sum of per-neighbour unit normals
  (axis-aligned), normalized. Points from the margin cell toward the
  foreign-plate neighbour — toward the trench on the subducting side.
  The `slab_pull_damps_convergent_motion` test verifies the resulting
  quadratic form is positive on convergent motion.
- **Operator discretization:** γ and n̂ live at cell centers; the
  bilinear form `⟨u, A·v⟩ = Σ_cells γ · (u·n̂)(v·n̂)` is symmetric by
  construction. On the staggered grid we compute `m(i,j) = γ · (v·n̂)`
  at centers, then average `γ·n̂·m` to each face from its two
  adjacent cells. Jacobi and SSOR preconditioners include the
  matching diagonal contribution so they stay consistent with the
  operator.

### 2.2 Files **not** touched (and why)

- `recycling.rs`: reads `source_rate`, not `T_plates`; no change
  needed (Phase 1 §3).
- `plates.rs` (`rebuild_traction`, `to_traction_field`): still builds
  `T_plates` from `plate.velocity`, but `plate.velocity` is now the
  pure kinematic velocity without slab-pull contamination, so the
  ridge-push branch stays correct. Verified `plate.subducted_mass`
  is still accumulated (may be consumed by downstream diagnostics).
- `traction.rs`: unchanged.

### 2.3 Obsolete tests

- `slab_pull_increases_plate_velocity`
  ([tectonics.rs:1717](crates/ymir-core/src/tectonics/solver/tectonics.rs#L1717)):
  `#[ignore = "obsolete since issue #75 ... see follow-up cleanup issue #79"]`
- `slab_pull_capped`
  ([tectonics.rs:1747](crates/ymir-core/src/tectonics/solver/tectonics.rs#L1747)):
  same treatment.

Neither was deleted outright — they remain as compile-checked
documentation of the old behaviour until the follow-up #79 cleanup PR
removes them.

### 2.4 New tests (all passing)

Added to `stokes.rs` test module:

| Test | What it verifies |
|------|------------------|
| [operator_with_slab_pull_is_symmetric](crates/ymir-core/src/tectonics/solver/stokes.rs#L965) | `⟨u, A·w⟩ = ⟨A·u, w⟩` with a non-trivial γ_slab band, rel_err < 1e-10 |
| [slab_pull_term_is_semi_definite](crates/ymir-core/src/tectonics/solver/stokes.rs#L996) | With η = 0 and friction = 0, `⟨v, A·v⟩ ≥ 0` for 10 random v — only the slab-pull contribution is left |
| [slab_pull_damps_convergent_motion](crates/ymir-core/src/tectonics/solver/stokes.rs#L1021) | Convergent-motion v produces `⟨v, A_slab·v⟩ > 0`, confirming the sign convention damps (does not excite) convergent motion |

Existing tests `operator_is_symmetric_*`,
`operator_is_positive_definite`, `operator_with_friction_is_symmetric`,
`jacobi_precond_*` all still pass, confirming the `Option<&SlabPullField>`
addition doesn't perturb the bare-operator code path.

---

## 3. Observations and caveats

### 3.1 Most of the speed-up came from *removing* the velocity boost, not from the new operator term

Looking at B vs C in Phase 2, every metric is identical within 1%:
`eta_ratio`, `tp_norm`, `resid_local`, `gpe_spike_p95`, wallclock.
That means the newly added `γ·n̂⊗n̂` operator term has **negligible
numerical effect** at the default `slab_pull_factor = 0.05`. Nearly
all of the 7× wallclock improvement is attributable to the deletion
of `apply_slab_pull` from the pipeline.

This is not a bug — it means the auto-regulation is *gentle* at the
default tuning. γ_seed magnitudes are around
`0.05 · |source_rate| ≈ 0.025` per cell at the margin, compared to
viscous diagonal terms of order `η / dx² ≈ 640`. The operator
contribution is ~4 orders of magnitude smaller than the viscous
stencil at the margin cells where γ is non-zero. The resulting
damping force on convergent motion is correspondingly small.

**For reviewer attention.** If the physics benefit (actually damping
subduction velocity to something geologically plausible) is the
priority, `slab_pull_factor` default likely needs to be raised
10–100× to produce a visible effect. If the stability benefit
(preventing the η-contrast cascade) was the main goal, it is
achieved simply by the removal — the new operator term is effectively
a no-op placeholder at current tuning. Both are defensible; flagging
it here so the choice is explicit.

### 3.2 Mass conservation

Not re-measured in this phase — the advection path and recycling
pipeline are unchanged. The Phase 1-bis mass-balance logs would
capture drift if it appeared. Can be confirmed with a follow-up 300-step
run if desired.

### 3.3 Morphological non-regression

Not measured — this would require exporting heightmaps at step 120
for B-baseline and B-Phase-2 and diffing. The scenario runner does
not currently export heightmaps. Recommended as a follow-up if visual
regression is a concern; the quantitative metrics (B = C for
thickness statistics via `tp_norm` proxy, matching residual) are the
strongest indirect evidence that morphology is preserved.

### 3.4 What stayed the same (good)

- Boundary detection (`compute_boundary_sources` classification loop)
  produces the same `boundary_type` and `source_rate` values; only the
  *additional* γ_slab/n̂ fields are new.
- GPE ridge-push term in `compute_rhs` is unchanged.
- Mass recycling reads `source_rate` and is unaffected.
- `plate.subducted_mass` is still accumulated by
  `accumulate_subducted_mass` — available for downstream diagnostics
  that may still expect it (e.g. UI displays).

---

## 4. Success signature vs. target (from Phase 1-bis §7)

| Target | Hit | Measurement |
|--------|-----|-------------|
| `eta_ratio` drops from ~62 toward C's ~11 | ✅ | 11.23 (exactly matches C) |
| Wallclock drops from baseline ~841 s toward C's ~177 s (scaled: 336 → 71 at 120 steps) | ✅ | 48 s, **below C's scaled 71 s** (likely because Phase 2 also slightly speeds up C itself via upstream fixes absorbed from the working branch since Phase 1-bis was run) |
| BiCGSTAB iters drop ≥ 2× | ✅ | 4× drop in per-substep solve time |
| Morphology qualitatively consistent | ⏳ | Not verified visually; all quantitative proxies match C |

---

## 5. Follow-ups

- **#79 (new, to file):** delete the two `#[ignore]`d obsolete slab-pull
  tests once Phase 2 ships; audit whether `plate.subducted_mass`
  is still needed anywhere.
- **Tuning of `slab_pull_factor`:** current default 0.05 gives a
  nearly inert operator term; if physically meaningful
  auto-regulation of subduction velocity is desired, raise the
  default. Requires a morphology-preservation study.
- **Mantle flow:** Phase 1-bis showed mantle flow is cost-neutral.
  Parent #74 may still want to remove or reformulate it for
  consistency; tracked separately.
- **GPE gradient smoothing (#78):** still the remaining RHS spike
  (`gpe_spike_p95` = 1.9 in all Phase 2 scenarios). Separate issue.
- **Multigrid preconditioner:** Phase 1-bis §6 identified this as
  complementary; still on the roadmap.

---

## 6. Reproducing

```bash
# Build
cargo build --release --workspace

# Unit tests (expect 3 new passing, 2 ignored)
cargo test --lib -p ymir-core

# Empirical validation (~100 s total at 120 steps × 3 scenarios)
for s in A B C; do
  ./target/release/examples/phase1bis_scenarios.exe "$s" 1 logs/phase2 120
done >> logs/phase2/summary.txt

# Aggregate
PHASE1BIS_STEP_LO=40 PHASE1BIS_STEP_HI=80 \
  ./scripts/phase1bis_aggregate.sh logs/phase2
```
