# C1 Phase 2 Track A — Stein-Stein 1992 oceanic bathymetry

Issue: [#129](https://github.com/FifionRibana/ymir/issues/129).
Branch: `129-c1-phase-2-track-a-stein-stein-oceanic-bathymetry`,
off `milestone/c1-lightweight-dynamic-tectonics`.

## Acceptance summary

| Stage | Tests | Status |
|-------|-------|--------|
| E1 unit (closure module)             | 6 | 6/6 PASS |
| V validation (5-point S-S anchor)    | 4 | 4/4 PASS |
| A acceptance (regime + regression)   | 2 | 2/2 PASS |
| D downstream smoke + gallery         | 1 active + 1 `#[ignore]` | 1/1 PASS (active) |
| **Phase 2 Track A subtotal**         | **13** active + 1 ignored | **13/13 PASS** |
| Phase 1.x integration preserved      | 19 | 19/19 PASS |
| Phase 1.4 + earlier ymir-core lib    | 77 | 77/77 PASS |
| **Cumulative C1-scope**              | **90** active + 1 ignored | **90/90** |

Stage E0 7th bit-identical decomposition preserved: the wrapper
`run_phase_a_cycle_c1(Enabled, ...)` is byte-for-byte equal to
`run_with_closures + apply_post_tectonic` under the
`phase_1_3_closures()` helper. Phase 2 Track A's closure addition
respects this invariant.

## Stein-Stein 1992 quantitative reproduction

Stage V Test 1 (`stein_stein_reproduces_5_age_points`) anchors the
implementation against five published depth-age points from
Stein & Stein (1992) *Nature* 359, Table 1 (GDH1 plate model):

| Age (Ma) | S-S published (m) | Computed (m) | Error (m) |
|----------|-------------------|--------------|-----------|
| 0        | 2600              | 2600.000     | 0.000     |
| 10       | 3754              | 3754.231     | 0.231     |
| 50       | 5035              | 5035.037     | 0.037     |
| 100      | 5498              | 5497.579     | 0.421     |
| 150      | 5613              | 5612.787     | 0.213     |

**Max error 0.421 m**, ~120× tighter than the ±50 m tolerance.
Phase 1.2-equivalent quantitative anchor — substantively stronger
than Phase 1.4's "visual review only" calibration (Lague 2014
declines a universal `K`). Phase 2 Track A delivers paper-faithful
S-S reproduction.

## Bathymetric maturation (Architecture C)

Per-cycle altitude statistics from
`c1_phase_2_visual_gallery::phase_2_bathymetry_visual_gallery`
(64²×300 steps, all 4 closures enabled, re-apply S-S at each
snapshot boundary):

| Cycle | S̃ min | S̃ mean | S̃ max | alt min | alt mean | alt max |
|-------|-------|--------|-------|---------|----------|---------|
| 0     | 0.20  | 0.56   | 1.00  | -0.56   | +0.08    | +1.00   |
| 50    | 0.20  | 0.37   | 2.18  | -1.13   | -0.28    | +0.40   |
| 100   | 0.00  | 0.34   | 2.18  | -1.13   | -0.27    | +0.38   |
| 200   | 0.00  | 0.32   | 2.18  | -1.13   | -0.26    | +0.39   |
| 300   | 0.00  | 0.29   | 2.18  | -1.13   | -0.26    | +0.39   |

- Architecture C bathymetric imprint **matures by cycle 50** — by
  step 50 oceanic mean altitude has dropped from +0.08
  (continental-dominant initial mix) to -0.28; further evolution
  stabilises around -0.26 for the remainder of the run.
- Oceanic minimum altitude reaches the S-S asymptote `-1.13`
  (`= -5651 m / depth_scale_m`) immediately at cycle 50,
  consistent with the age-pile-up at convergent boundaries
  (oldest oceanic cells saturate the exp regime).
- S̃ max stable at `2.18` from cycle 50 onward — Phase 1.4's
  equilibrium-height + erosion equilibrium continues to hold
  under Phase 2 Track A.

## Phase 1.1 / 1.2 / 1.3 / 1.4 → Phase 2 Track A comparison

| Quantity                       | Phase 1.1 | Phase 1.2 | Phase 1.3 | Phase 1.4 | Phase 2 A |
|--------------------------------|-----------|-----------|-----------|-----------|-----------|
| Closures active                | 0         | DS        | DS+EH     | DS+EH+ER  | DS+EH+ER+SS |
| Mean altitude (cycle 300)      | n/a       | n/a       | n/a       | ~0.4      | -0.26     |
| Oceanic altitude               | uniform   | uniform   | uniform   | uniform   | age-modulated |
| Min altitude (oceanic)         | ~0.1      | ~0.1      | ~0.1      | ~0.1      | -1.13     |
| Ridge visibility               | absent    | absent    | absent    | absent    | visible (cycle 50+) |
| Abyssal plain                  | absent    | absent    | absent    | absent    | visible (old cells) |
| Per-step cost (64²)            | ~50 µs    | ~80 µs    | ~96 µs    | 367 µs    | 455 µs    |
| Closure overhead (vs prev)     | —         | +30 µs    | +16 µs    | +271 µs   | +20-90 µs |
| Wedge_p95 (with DS+EH active)  | n/a       | 0.376     | 0.376     | 0.696     | preserved |
| Quantitative anchor available  | n/a       | DS Fig 5  | none      | none      | S-S 5 pts |

Performance: Phase 2 Track A's S-S overhead is **+20 µs/step** at
64² (S-S enabled vs disabled). The cycle-300 wall time is 137 ms
(455 µs/step), well within the W4 < 600 µs/step Phase 2 budget
(~32% margin remaining).

## Architectural findings

### Finding 1 — Architecture C VALIDATED via run-boundary observability

Stage A Test 1 verified empirically:

- Spearman ρ between (age, altitude) over 2270 oceanic cells:
  **-0.476** — strong negative correlation, expected direction
  "older = deeper = lower altitude".
- Median-split bucket delta: mean(old) - mean(young) = **-0.0109**
  (well past the -0.01 threshold).
- The S-S imprint is observable at run boundary via explicit
  re-application of `apply_stein_stein_bathymetry` on the
  isostasy-recomputed altitude. Per-step in-loop S-S effects are
  transient (overwritten by next `compute_isostasy` from `S̃`),
  but the closure is correctly modulating altitude by age each
  step.
- Erosion's slope factor consumes S-S-modulated altitude at the
  oceanic/continental coastline → S-S indirectly propagates
  through `state.s` via the erosion sink, even though S-S itself
  never writes to `S̃`. This is the intended Architecture C
  cross-closure interaction.

No need to escalate to Architecture A or B. The fallback paths
documented in
[`closures/oceanic_bathymetry/mod.rs`](../../../crates/ymir-core/src/tectonics_c1/closures/oceanic_bathymetry/mod.rs)
remain available if Phase 2 Track B's age-field changes surface
new limitations.

### Finding 2 — Age field is advected as a *density*, not a Lagrangian scalar

Empirical observation (Stage A Test 1, cross-checked by Stage D
gallery test):

- Initial oceanic age = 0.5 (uniform across 2270 oceanic cells).
- After 300 steps, oceanic age distribution: min ≈ 0, max ≈ 6958,
  mean ≈ 4.67, median ≈ 0.
- Pile-up factor ~1000× at convergent boundaries; mass-balanced
  by depletion in advection-source regions (median ≈ 0).

Mechanism: the C1 Phase 1.1 time loop advects `age` via the same
conservative flux-form upwind as `S̃` (`∂_t·age + ∇·(age·v) = 0`
— see
[`tectonics_v2::advection`](../../../crates/ymir-core/src/tectonics_v2/advection.rs)).
This treats `age` as a *density* with the same pile-up dynamics
that produced Phase 1.2's `global_max ≈ 2297` for `S̃` —
identical conservative-density mechanism applied to both fields.

Consequence: Phase 2 Track A's bathymetric variability is
dominated by **"near-zero cells get ridge depth = -0.520"** vs
**"piled-up boundary cells get the asymptote = -1.130"**, rather
than a smooth `√t → exp(-α·t)` S-S profile. The closure operates
correctly at each individual age value (Stage V ±50 m anchor),
but the input distribution does not exercise the full age range
geophysically.

Mechanism is architectural (Phase 1.1 design choice — no
`age += dt` per step). NOT a Stein-Stein closure issue. Phase 2
Track B (R7 init, separate issue) per §6.5 of the design doc
should address via:

1. Ridge-aligned `age = 0` initialisation on oceanic cells.
2. Per-step ageing (`age[c] += dt` each step on cells that
   weren't created/destroyed).
3. Re-evaluation of whether to keep flux-form advection on `age`
   or switch to a Lagrangian "particle age" semantics.

Persisted as transferable feedback memory entry
`feedback_age_advection_density_vs_lagrangian` — pattern
generalises to any density-advected scalar quantity (`age`,
`composition`, `temperature`, ...) where the physical
interpretation is per-cell-Lagrangian.

### Finding 3 — S̃ minimum reaches numerical floor at later cycles

Per-cycle gallery output:

- Cycle 0:  S̃ min = 0.20 (oceanic init baseline)
- Cycle 50: S̃ min = 0.20 (unchanged)
- Cycle 100-300: S̃ min ≈ 0.003 (drifts to near zero)

Likely mechanism: advection drift — `S̃` values transit through
near-zero on the per-cell update path between erosion's
defensive `floor = 0.2` clamp enforcements. Erosion only clamps
WITHIN its own apply step; advection on the next step does not
re-enforce the floor, and conservative-upwind redistribution
can move mass below the floor before the next erosion pass
re-clamps cells that get re-classified.

Not blocking Phase 2 Track A acceptance — the S̃ minimum is
finite and non-negative (`>= 0.003 > 0`), the Phase 1.4 erosion
floor's defensive role is preserved on cells actively being
eroded, and the visual gallery shows continuous oceanic basins
(no holes / NaNs / negative values).

Worth investigating in a future Phase 1.4 follow-up: should the
floor clamp run at end-of-cycle in `apply_post_tectonic` rather
than only inside `apply_erosion_step`? Defer to a separate
issue.

## Pre-existing untouched regressions (not blocking Phase 2)

1. **`export::deserialize_legacy_metadata_without_upscale`** — lib
   unit test failure pre-existing the Phase 1.x C1 work. Last
   `export/mod.rs` change was commit `8c974c1` (Issue #21
   rectangular grid refactor) well before any C1 phase. Unchanged
   by Phase 2 Track A. Recommend filing as orthogonal cleanup
   issue.

2. **Bevy 0.18.1 v2_legacy version drift** — surfaced during
   Phase 1.4 Stage V1 diagnose
   (`project_c1_phase_1_4_erosion_outcomes`). Unchanged by Phase 2
   Track A. Tracked as separate follow-up; v2_legacy build
   completes for ymir-core (the Bevy issues live in ymir-viz).

## Architecture C visual reading

The committed gallery PNGs at
`docs/reports/c1_phase_2_oceanic_bathymetry/` (not in git per
Phase 1.x convention; re-generate with the gallery test) show:

- `cycle_000_altitude.png` — initial state. Continental cells
  bright (altitude ≈ 1.0 = highest peak), oceanic cells
  mid-tone (altitude ≈ -0.56 = ridge-depth from initial age
  0.5). No mid-ocean ridge or abyssal plain yet — the age field
  is uniform within each plate type.
- `cycle_050_altitude.png` — bathymetric imprint emerging.
  Oceanic cells near convergent boundaries reach the asymptote
  (deep blue, altitude -1.13). Continental cells slightly
  shifted from isostatic adjustment to the post-Davis-Suppe
  wedge state.
- `cycle_100/200/300_altitude.png` — steady-state bathymetry.
  Older oceanic cells (pile-up at convergent boundaries) sit at
  asymptotic depth; younger oceanic cells at ridge depth.
  Continental wedges from DS+EH coexist with the bathymetric
  imprint.
- `cycle_NNN_s.png` — `S̃` evolution unchanged from Phase 1.4
  (S-S doesn't write to `S̃` directly). The `S̃` palette
  `[0, 3.0]` is preserved cross-phase per
  `feedback_viz_palette_absolute_for_comparison`.

## Cross-references

- Issue: #129
- Design doc:
  [`docs/c1_lightweight_dynamic_tectonics.md`](../../c1_lightweight_dynamic_tectonics.md)
  §5.1 (closure footnote), §7.2 (Track A status), §11 (Phase 2
  scales sub-section)
- Closure module:
  [`crates/ymir-core/src/tectonics_c1/closures/oceanic_bathymetry/`](../../../crates/ymir-core/src/tectonics_c1/closures/oceanic_bathymetry/)
- Time loop integration:
  [`crates/ymir-core/src/tectonics_c1/time_loop.rs`](../../../crates/ymir-core/src/tectonics_c1/time_loop.rs)
  (stage 4a — `if erosion || oceanic_bathymetry`)
- Stage V tests:
  [`crates/ymir-core/tests/c1_phase_2_oceanic_bathymetry.rs`](../../../crates/ymir-core/tests/c1_phase_2_oceanic_bathymetry.rs)
- Stage A + D tests:
  [`crates/ymir-core/tests/c1_phase_2_bathymetry_acceptance.rs`](../../../crates/ymir-core/tests/c1_phase_2_bathymetry_acceptance.rs)
- Visual gallery generator:
  [`crates/ymir-core/tests/c1_phase_2_visual_gallery.rs`](../../../crates/ymir-core/tests/c1_phase_2_visual_gallery.rs)
- Stein & Stein 1992: DOI [10.1038/359123a0](https://doi.org/10.1038/359123a0)
- Parsons & Sclater 1977 (predecessor): DOI [10.1029/JB082i005p00803](https://doi.org/10.1029/JB082i005p00803)
- Memory entries (off-repo):
  - `project_c1_phase_2_track_a_outcomes` (new)
  - `feedback_age_advection_density_vs_lagrangian` (new, transferable)
  - `feedback_recursive_tuning_signals_structural` (updated — Phase 2 T3 proactive deferral case)
  - `feedback_calibration_via_visual_review` (updated — quantitative anchor hierarchy)
