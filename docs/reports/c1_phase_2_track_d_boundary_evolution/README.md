# C1 Phase 2 Track D — Boundary Evolution

**Issue:** [#132](https://github.com/FifionRibana/ymir/issues/132)
**Branch:** `132-c1-phase-2-track-d-boundary-evolution-subduction-accretion-rifting`
**Status:** ✓ **Complete, pending PR merge.**

Track D is the third Phase 2 work track (after Track A oceanic
bathymetry and Track B R7 init). It delivers the **boundary
evolution closures** (subduction, accretion, rifting) that mutate
`plate_id`, `plate_type`, and `kinematics.velocities` per-step —
the first C1 work-track to lift the static-classification
optimisation in `run_with_closures`.

## Acceptance summary

| Stage | SHA | Description |
|-------|-----|-------------|
| S      | `e4f4005` | Track D doc prerequisites + codebase exploration (5 architectural concerns surfaced) |
| E1     | `74d3b9e` | Subduction closure module + 6 unit tests |
| E2     | `95dd4db` | Accretion mechanism + `ConvergenceTracker` + 6 unit tests |
| E3     | `81d9e1f` | Rifting closure + `DivergenceTracker` + split + 8 unit tests |
| E4     | `7271515` | Track D integration + mass-balance diagnostic + Phase 1.x regression preserved |
| V      | `66badaa` | Track D validation 14 tests + Track C SUFFICIENT verdict |
| A      | `07f3e5e` | Track D acceptance 3 tests + visual gallery `#[ignore]`'d |
| D      | `3e41df9` | Downstream regression sweep + visual review + downstream cross-reference |
| Final  | *(this commit)* | README + memory entries persisted |

**Total: 9 commits.** Effort estimate 16-21 days realised at the upper end.

### Cumulative test inventory (Track D contribution)

- **Lib closure tests (20):** 6 subduction + 6 accretion + 8 rifting unit tests in `tectonics_c1::closures::{subduction, accretion, rifting}::`.
- **Integration mass-balance (1):** `c1_phase_2_track_d_mass_conservation::mass_conservation_holds_per_step_100_run` — drift `6.253e-12` vs tolerance `1e-6`.
- **Stage V validation (14):** `c1_phase_2_boundary_evolution.rs` — subduction integration ×4, accretion integration ×3, rifting integration ×4, integration ×2, multi-seed Track C escalation ×1.
- **Stage A acceptance (3):** `c1_phase_2_track_d_acceptance.rs` — Q-V.1 Option B event-firing, 9th bit-identical decomposition preservation at Phase 2 R7 scope, Phase 2 milestone gate proxy.
- **Stage A visual gallery (2 `#[ignore]`'d):** main 5-cycle seed-42 gallery (10 PNGs) + multi-seed diversity at cycle 300 (6 PNGs).

**Total Track D active tests: 38.** All PASS. Cumulative ymir-core test count after Track D ≈ 88.

### Pre-existing untouched regressions

- `export::deserialize_legacy_metadata_without_upscale` — Issue #21, predates Phase 1.x.
- `rectangular_simulation_smoke_test` — v1 Stokes `NonlinearSolverDidNotConverge`. Verified pre-Track-D via `git stash` (fails identically without my changes).

Both unchanged by Track D scope.

## Sub-components delivered

### 1. Subduction closure (`tectonics_c1/closures/subduction/`)

Rate-based oceanic consumption + arc volcanism distribution + floor-triggered `plate_id` reassignment (Q1.2-Q1.4):

- `apply_subduction_step` — per-cell loop on `BoundaryType::Convergent` + `plate_type[c] == Oceanic`. Picks the continental neighbour with the largest positive `v_rel · n̂`. `Δs = consumption_rate × convergence × dt`, clamped at `s_before`.
- `distribute_arc_mass` — BFS from the consuming cell up to `arc_distance` cells, equal-share distribution on continental cells reached. Arc-mass-lost case (`continental_cells.is_empty()`) returns `0.0` (graceful — architectural concern documented inline).
- Floor reassignment when `s_after < plate_id_reassign_threshold` (default `0.05` = 1/4 oceanic baseline). Mutates `plate_id` + `plate_type` to the continental neighbour. **First C1 closure to mutate `plate_id` and `plate_type`.**

**Default parameters (analytical first-pass, Stage E1 W7):**
`consumption_rate = 0.5` (K_subduction), `arc_efficiency = 0.5`, `arc_distance = 3`, `plate_id_reassign_threshold = 0.05`. 5th C1 first-shot calibration success.

**Reference:** Lallemand, Heuret & Boutelier 2005, *G-cubed* 6(9), Q09006.

### 2. Accretion mechanism (`tectonics_c1/closures/accretion/`)

Sustained-convergence plate-id merge (Q2.4, Q3-revised):

- `ConvergenceTracker` — per-pair counter keyed by canonical `(a, b)` with `a < b`. `update()` re-scans grid, increments net-convergent pairs, resets non-convergent to 0.
- `apply_accretion_step` — when `count >= merge_time_threshold`, merge lower-index winner absorbs higher-index loser. Mass-weighted average velocity `v_new = (m_a v_a + m_b v_b) / (m_a + m_b)`. **First C1 closure to mutate `kinematics.velocities`.**
- No `S̃` thickening — Davis-Suppe handles orogenic morphology during the pre-merge phase.

**Q3 revision (20 → 50)**: original design exploration set `merge_time_threshold = 20`. Stage E1 W7 analytical refined to **50 steps (~33 Ma)** for spurious-merge suppression. Documented in `params.rs` + memory entry `feedback_recursive_tuning_signals_structural`.

**Reference:** Coney, Jones & Monger 1980, *Nature* 288, 329-333.

### 3. Rifting closure + split (`tectonics_c1/closures/rifting/`)

Two-stage mechanism (Q3.2 chewing-gum cut, Q3.4 perpendicular offset):

- `apply_rifting_thinning` — closure portion. Negative `S̃` source on continental cells classified as `BoundaryType::Divergent`. Mirror of Davis-Suppe orogenic source on convergent boundaries.
- `DivergenceTracker` — symmetric mirror of `ConvergenceTracker` (sign flipped on the verdict line).
- `apply_rifting_split` — event mechanism. Fires when BOTH conditions hold simultaneously:
  1. `divergence_count >= split_time_threshold` (sustained extension).
  2. `min(S̃) at rift_strip < split_thickness_threshold` (sub-threshold thinning).

Partitioning: pair.0's rift strip (continental cells of plate `a` 4-adjacent to plate `b`) becomes a new plate. `kinematics.velocities` extended with parent + right-hand-rule perpendicular offset. Path 3.B `age = 0` on all reassigned cells.

**Default parameters (Stage E2 W7):**
`thinning_rate = 1.0` (K_rift — 6th C1 first-shot), `split_time_threshold = 75` (~50 Ma sustained), `split_thickness_threshold = 0.7` (McKenzie β = 1.4 stretching cap), `split_velocity_offset = 0.005`, `plate_id_cap = 256`.

**References:** McKenzie 1978, *EPSL* 40(1), 25-32; Buck 1991, *JGR* 96(B12), 20161-20178.

### 4. Mass conservation diagnostic

`c1_phase_2_track_d_mass_conservation.rs::mass_conservation_holds_per_step_100_run`:

```text
Track-D-only mass-conservation invariant:
    delta = consumed - arc_distributed + thinning_removed

100-step run drift  = 6.253e-12  vs tolerance 1e-6  (~6 orders of magnitude tighter)
300-step run drift  = 1.606e-12  vs tolerance 1e-6  (Stage V mass_conservation_holds_300_steps_full_stack)
```

Machine-precision algorithmic mass balance. Per-cycle accumulator validates the per-step pipeline's mass budget against design doc §5.4.

### 5. Path 3.B event-driven `age = 0`

Track B's Path 3.A applied ridge-aligned `age = 0` at **init time only**. Path 3.B extends this to **mid-simulation rifting events**: every cell reassigned to a newly-spawned plate via `apply_rifting_split` has its `age` reset to `0`. Preserves Track B's Spearman age-altitude correlation under Track D mutation; without Path 3.B, rifted-off cells would carry their parent plate's advected age (implausible "ancient rift floors").

## Phase 1.4 / Track A / Track B / Track D comparison

| Metric | Phase 1.4 | Track A | Track B | Track D |
|--------|-----------|---------|---------|---------|
| Issue | #127 | #129 | #131 | #132 |
| `plate_id` mutation | static | static | static (R7 init only) | **dynamic per-step** |
| `plate_type` mutation | static | static | static | **dynamic** (subduction reassign) |
| `kinematics.velocities` | static | static | static | **dynamic** (accretion + rifting) |
| Per-step cost (64²) | 367 µs | 455 µs | 465 µs | **838.7 µs** |
| Cost increment vs predecessor | +271 µs | +88 µs | +10 µs | **+374 µs** |
| Cost source increment | erosion + isostasy | S-S anchor | R7 init dispatcher | Track D recompute + 3 closures |
| First-shot calibrations | `K = 0.001` | 5-pt anchor | R7 amplitude/freq | `K_subduction`, `K_rift`, thresholds |
| `S̃` global_max (seed 42, cycle 300) | ≈ 2.18 | ≈ 2.18 (S-S re-applied) | ≈ 2.18 | clamped per `apply_subduction_step` |
| Spearman ρ (cycle 300) | n/a | -0.476 | **-0.5233** | (Track D mutates plate_id; not directly comparable) |
| Plates remaining (seed 42) | 8 | 8 | 8 | **2** (6 merges) |
| Continental cells (seed 42, cycle 300) | varies | varies | similar to Phase 1.4 | **1123 / 4096** |
| Mass-balance invariant tolerance | 1e-6 (W-T floor caveat) | 1e-6 | 1e-6 | **6.25e-12** (Track-D-only path) |
| Visual signature | wedges + erosion | + bipolar altitude | + R7 curved boundaries | + Pangaea-like supercontinent |

## Architectural findings (8)

### 1. Per-step cost 838.7 µs at 64² × 300 — Phase 3 optimisation candidate

Track D forces per-step `classify_boundaries` + `wedge_distance_intra_plate` recompute when ANY Track D closure is enabled, lifting the static-classification optimisation from Phase 1.2. Cost: ~200 µs/step. Combined with Track D closure work (~170 µs) the total exceeds the 800 µs Phase 2 W4 budget by 4.8 %.

**Decision:** NOT blocking acceptance; Phase 3+ optimisation issue drafted (conditional skip when no Track D event fired previous step).

### 2. Plate-count collapse 8 → 1-2 by cycle 300 — Pangaea narrative ACCEPT

Default `merge_time_threshold = 50` produces 6-10 accretion merges per 300-step run, collapsing the 8 initial plates to 1-2 surviving supercontinents. Stage D visual review confirmed gradual progression (cycle 0 → 100 → 200 → 300), consistent with ~200 Ma of geological evolution (real Pangaea formed over ~150 Ma).

**Decision:** Option A — ACCEPT the Pangaea narrative. Track D-bis recalibration issue drafted as optional follow-up if Phase 3+ narrative wants multi-continent worlds.

### 3. Asymmetric divergence: plate_id ≠ geographic similarity

Multi-seed Test 3 measured pairwise plate_id divergence: 70.5 % / 70.5 % / 12.3 % (mean 51.1 %). Visual inspection of `seed_diversity/` PNGs showed seeds 1337 and 2026 produce **visually DISTINCT** continents despite the 12.3 % numeric divergence.

**Architectural insight:** plate_id divergence is a **structural** metric (canonical-low-wins merge rule produces dominant plate id `0` in both seeds after collapse), NOT a **geographic-similarity** metric. For Phase 2+ milestone gate proxies, plate_id divergence is a conservative lower bound; future gates should add `plate_type` mask comparison and altitude distribution diversity for direct geographic-diversity signal.

Captured in memory entry `feedback_fill_ratio_regime_agnostic_metric` as a regime-specific vs regime-agnostic metric distinction.

### 4. Rifting splits rare by design (0-3 / seed at 300 steps)

The chewing-gum-cut split mechanism requires BOTH sustained-divergence (≥ 75 steps) AND sub-threshold thinning (`S̃ < 0.7`). Stage V multi-seed scan: seeds 42/1337/9999 = 0 splits, seed 100 = 1 split, seed 2026 = 3 splits.

**Decision:** Rare-by-design (Atlantic-style rifting takes ~150 Ma in nature). NOT a bug. Path 3.B age=0 mechanism still verified by seed 2026's 2 visible rift basins in the gallery.

### 5. No palette-clip artifacts

Altitude `[-1.13, +1.13]` symmetric palette + `S̃ [0, 3.0]` (Q-V.3 Option A — Track A/B continuity). `print_stats` clip marker (`[clip-low]` / `[clip-high]`) NOT triggered at any snapshot in the gallery. Cross-phase visual comparability preserved.

### 6. Mass-balance machine precision

Algorithmic mass balance for the Track-D-only path: drift `6.253e-12` (100 steps) and `1.606e-12` (300 steps) — ~6 orders of magnitude tighter than the design doc §5.4 tolerance of `1e-6`. Verifies the invariant:

```text
Σ S̃_initial − Σ S̃_final = total_consumed − arc_distributed + thinning_removed
```

### 7. 9th bit-identical decomposition preserved

`c1_phase_a_decomposes_into_closures_then_post_tectonic` PASS EXACT through Track D scope addition. Also verified at Phase 2 R7 + Track D disabled scope via Stage A Test 2 (`max |S̃_wrapper − S̃_decomposition| = 0`). Stage E0 invariant preserved.

### 8. Per-cell `plate_type` consistency through Track D mutations

Initial Stage E4 W3 concern dismissed: per-cell plate_type is preserved through accretion + rifting (continental → continental, oceanic → oceanic). Subduction's floor-trigger reassignment re-types cells inline at the (i, j) being processed, so per-cell consistency holds at all times. No "recompute plate_type from plate_id" step needed.

## Phase 2 milestone closeout

| Track | Status | Notes |
|-------|--------|-------|
| A — oceanic bathymetry (Stein-Stein 1992) | ✓ **merged** PR #130 | Architecture C, ρ = -0.476 baseline |
| B — R7 init (boundary displacement + clustering + Path 3.A) | ✓ **merged** PR #133 | ρ = -0.5233 IMPROVES on Track A by 0.047 |
| D — boundary evolution (subduction + accretion + rifting + Path 3.B) | ✓ **this PR** | Pangaea narrative accepted; 5/5 seeds with Track D events (Track C SUFFICIENT) |
| B-bis — cadrable viewport offset | 📋 pending | 9/10 seeds wrap; user-discretion follow-up (Track B Issue #131 deferral) |
| C — kinematics sampling | 📋 not urgent | Stage V multi-seed scan resolved Track C escalation; not on Phase 2 critical path |
| D-bis — plate-count recalibration + Path 3.B subduction extension | 📋 optional | Pangaea narrative accepted; defer unless Phase 3+ narrative needs more diversity |
| Phase 3 — per-step cost optimisation | 📋 separate issue | 838.7 µs > 800 µs budget; conditional-recompute skip |

**Design-doc §7.2 acceptance gate**: *"the same preset run multiple times with different seeds produces visually distinct continents (not just rotations of the same shape)."* — **ACHIEVED.** Mean pairwise plate_id divergence 51.1 % across 3 seeds (Stage A Test 3); visually distinct continental geometries confirmed in seed_diversity gallery (Stage D Q3 REFUTE attractor hypothesis).

## Methodological consolidation

Track D applied (and validated) the following discipline patterns consolidated through Phase 1.4 / Track A / Track B:

- **Surface-before-implement (W7)** at every stage transition — 3 spec extensions surfaced + 1 fixture-pathology surfaced + 1 architectural-finding refutation.
- **Match-request-scope** — Phase 1.x + Track A/B regression preservation via Option B per-test-fixture disablement; Track B-bis NOT bundled.
- **Init-re-architecture-pattern** — `apply_subduction_step` / `apply_accretion_step` / `apply_rifting_split` shipped as new parallel entries; existing closures untouched.
- **Calibration-via-visual-review tier 2** — analytical first-pass for all 6 Track D defaults; visual review at Stage A confirmed Pangaea narrative; Q3 revision (20 → 50) documented as empirical refinement.
- **Recursive-tuning-signals-structural-limit** — Stage E1 W7 analytical refined Q3 (20 → 50) BEFORE shipping; avoided the iteration-creep pattern.
- **Visual-validation-supersedes-scalar-metrics** — Stage D plate_id divergence numeric metric refuted by visual inspection (insight added to fill-ratio-regime-agnostic memory).

## Cross-references

- Issue [#132](https://github.com/FifionRibana/ymir/issues/132)
- Design doc [`c1_lightweight_dynamic_tectonics.md`](../../c1_lightweight_dynamic_tectonics.md) §4.5 (boundary evolution), §5.2 (Track D footnote), §7.2 (Phase 2 progress), §11 (parameter scales).
- Memory entries (off-repo): `project_c1_phase_2_track_d_outcomes` (created Stage Final), `feedback_age_advection_density_vs_lagrangian` (Path 3.B section added), `feedback_calibration_via_visual_review` (Track D first-shot section added), `feedback_recursive_tuning_signals_structural` (Q3 revision case added), `feedback_fill_ratio_regime_agnostic_metric` (plate_id ≠ geographic insight added).
- Track A README: [`docs/reports/c1_phase_2_oceanic_bathymetry/`](../c1_phase_2_oceanic_bathymetry/) (Architecture C foundation).
- Track B README: [`docs/reports/c1_phase_2_track_b_init_r7/`](../c1_phase_2_track_b_init_r7/) (Path 3.A + R7 init foundation).

PNGs at this directory: `cycle_NNN_*.png` (5-cycle main gallery, seed 42) + `seed_diversity/seed_NNNNN_*.png` (3-seed diversity at cycle 300). NOT committed per Phase 1.x + Track A/B convention; regenerate via `cargo test --release -p ymir-core --test c1_phase_2_track_d_visual_gallery -- --ignored --nocapture`.
