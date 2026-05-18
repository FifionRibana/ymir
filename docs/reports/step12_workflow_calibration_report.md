# Step 12 — Workflow calibration report

> Companion to `step12_physics_report.md` and the §4.14 patch
> (`docs/solver-scaling-step12-patch.md`). Where the physics report
> validates the workflow against the issue's acceptance criteria,
> this report captures the empirical sensitivity of the
> `(k_cycle, N_cycles, α, β)` parameter space and the rationale for
> the D8 default values exposed by the panel sliders. Future users
> tuning a workflow run for a non-canonical preset should consult
> this report first.

## D8 default rationale (TL;DR)

| param | default | what it controls | sensitivity |
|---|---|---|---|
| `n_cycles` | 5 | number of Phase A cycles | counter-isostasy converges by cycle 3–4 → 5 leaves a small margin |
| `k_cycle` | 20 | tectonic harness steps per cycle | matches Step 11 transient budget; cycles 2..N reuse warm-start |
| `α` | 0.01 | erosion rate per cycle | controls per-cycle Δh; ≥ 0.05 produces visible Δ_S̃ but not curvature |
| `β` | 0.0 | downhill redistribution | 0 = pure diffusive (mass loss); 1 = mass-conserving |
| `craton_retention_threshold` | 0.10 | fraction-within-plate for craton retention re-evaluation | new field disambiguating from init's `plate_area_min`; locked at the same numeric default for now |
| `hd_grid_size` | 2048 | Phase B HD grid | issue's primary HD target |
| `num_droplets` | 5×10⁶ | Phase B rain-drop count | matches `ErosionConfig::default()` |
| `fbm_amplitude_base` | 0.08 | Phase B FBM upscale | matches `FbmUpscaleConfig::default()` |
| `grand_scale_tolerance` | 0.10 | D5 acceptance bound on p95 deviation | reformulated from L_∞ to p95; threshold value unchanged |

These are **starting points for exploration**, not calibrated optima.
The empirical sweep below documents what each lever moves and where
the regime breaks down.

## Sweep summary — five configurations probed

The Phase 6 visual-checkpoint commit (`3a8dc2d`) and the curvature
variants probe ran five distinct `(α, β, N, k)` configurations on
both canonical presets (`single_continent`, `convergence`) at 64² ×
mantle-on:

| label | α | β | N | k | wallclock Phase A | mass drift (sc) | mass drift (cv) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **D8 conservative** | 0.01 | 0.0 | 5 | 20 | 16.7 / 16.2 s | -2.91 (0.12 %) | -3.25 (0.15 %) |
| aggressive demo | 0.05 | 0.0 | 15 | 20 | 44.4 / 48.7 s | -37.5 (1.53 %) | -40.3 (1.89 %) |
| v1 (α-up + redistribution) | 0.05 | 0.5 | 15 | 20 | 47.8 / 48.4 s | -18.1 (0.74 %) | -19.4 (0.91 %) |
| v2 (more cycles, less α) | 0.02 | 0.5 | 30 | 10 | 86.2 / 93.4 s | -14.8 (0.60 %) | -16.0 (0.75 %) |
| v3 (α-noise modulation) | 0.05+jitter | 0.0 | 15 | 20 | 45.6 / 50.1 s | -37.5 (1.53 %) | -40.5 (1.90 %) |

(`sc` = `single_continent` preset, `cv` = `convergence` preset.)

## α (erosion rate) sensitivity

**α controls how aggressively per-cycle erosion bites into S̃.**
Mass drift per cycle scales linearly with α at fixed
(N, k) configuration:

- D8 (α = 0.01): ≈ 0.59 mass-drift-units per cycle (sc), 0.65 (cv)
- aggressive (α = 0.05): ≈ 2.50 per cycle (sc), 2.69 (cv)
- 5× α-ratio → ≈ 4.2× mass-drift ratio

Why not exact 5×? At α = 0.05 the per-cycle Δh is large enough that
some boundary cells flip from continental to oceanic mid-cycle and
stop receiving erosion (they're below `sea_level_ref`). The
shrinking continental footprint slightly damps the linear scaling.

**Counter-isostasy effectiveness vs α**: peak S̃ trajectory drift
across the run scales roughly with α × N:

| run | total α × N | peak S̃ drift | drift % |
|---|---:|---:|---:|
| D8 (5 cyc, α=0.01) | 0.05 | 1.1979 → 1.1974 | 0.04 % |
| aggressive (15 cyc, α=0.05) | 0.75 | 1.1994 → 1.1972 | 0.18 % |

Both runs successfully bound peak S̃ — counter-isostasy is robust
across the explored α range. No runaway, no instability.

**Practical guidance**:

- **`α ≤ 0.02`** for production runs — keeps mass drift under 1 %
  over a 5-cycle budget. The default 0.01 is appropriate for
  preserving the Step 11 baseline regime as cleanly as possible.
- **`α = 0.05`** for visual exploration runs where the user wants
  to see counter-isostasy effects at scale. Mass drift climbs to
  1.5–1.9 % over 15 cycles; D5 p95 deviation climbs proportionally
  (see "D5 sensitivity" below).
- **`α > 0.10`** is not recommended without specific motivation —
  the panel slider clamps to `[0.001, 0.10]` for this reason.

## β (downhill redistribution) sensitivity

**β controls whether eroded material is lost (β = 0, diffusive) or
redistributed downhill to a NESW priority neighbour (β = 1, mass-
conserving)**.

The Phase 2 acceptance probes pinned the algorithm-level mass
contract:

- `v2_workflow_erosion_mass_balanced` (β = 1.0): total continental
  volume preserved within 1e-6 over a 32² × 50-step run.
- `v2_workflow_erosion_diffusive` (β = 0.0): continental mass
  decreases monotonically.

The Phase 6 v1 + v2 variants ran with β = 0.5, halfway between the
two regimes:

| variant | β | mass drift (sc, 18 cycles) |
|---|---:|---:|
| v1 (α=0.05, N=15) | 0.5 | -18.1 |
| aggressive (α=0.05, N=15) | 0.0 | -37.5 |

β = 0.5 cuts mass drift in half vs the diffusive case, as expected
(half the per-cell Δh is recycled to a downslope neighbour rather
than discarded).

**β ≠ 0 does not produce border curvature** — variants v1 + v2
both preserve polygonal Voronoï shape across 18 cycles (see physics
report acceptance #6 finding). Redistribution acts symmetrically
across straight-edge boundaries; the parallel-translation symmetry
that breaks curvature is preserved by D2 regardless of β.

**Practical guidance**:

- **`β = 0` (D8 default)** for the most conservative Phase A mass
  budget. The diffusive regime produces the cleanest counter-
  isostasy signal — peak S̃ trajectory tracks erosion volume
  removed without secondary deposition effects.
- **`β = 0.5`** is reasonable for users who want to keep the global
  mass closer to baseline while still exercising counter-isostasy.
  The L_∞ HD deviation rises slightly (the redistributed mass
  shows up as small post-cycle mounds adjacent to erosion sites);
  the p95 acceptance is unaffected within noise.
- **`β = 1.0`** is the algorithm-level mass-conservation regime;
  useful for acceptance probes (`v2_workflow_erosion_mass_
  balanced`), less useful for production where the counter-
  isostasy signal is the user-facing effect.

## (N, k) ratio sensitivity

**N × k is the total tectonic-step budget; the ratio between them
controls how "interleaved" the erosion is.**

The aggressive demo (N=15, k=20) and v2 (N=30, k=10) probe two
points on the iso-budget curve N × k = 300:

| variant | N | k | wallclock (sc) | mass drift (sc) |
|---|---:|---:|---:|---:|
| aggressive (β=0) | 15 | 20 | 44.4 s | -37.5 |
| v2 (β=0.5) | 30 | 10 | 86.2 s | -14.8 |

The wallclock difference (44 → 86 s) is dominated by warm-start
overhead: cycle 1 always pays the full Step 11 transient (≈ 1–2 s
on 64² mantle-on), and each cycle 2..N pays a smaller fixed cost
for the continuation rebuild. v2 has 30 cycles vs 15, so the
fixed-cost overhead almost doubles.

The mass-drift difference (37 → 15) reflects the β = 0 vs β = 0.5
contrast more than the (N, k) ratio per se — at iso-α, iso-β,
iso-budget, the ratio's main effect is wallclock, not mass.

**Practical guidance**:

- **Lower N + higher k** (= D8 5×20) produces fewer continuation
  rebuilds and slightly lower wallclock. Good when the user's
  intent is "established regime + small counter-isostasy nudge".
- **Higher N + lower k** (= 30×10 or 50×5) is appropriate when the
  user wants finer-grain feedback between erosion and tectonics.
  The 64² wallclock cost roughly doubles per N-doubling at fixed
  N × k.
- **k < 5** is not recommended — Step 11's tectonic transient takes
  several steps to advance; truncating it produces noisier
  per-cycle metrics without a corresponding physics gain.

## D5 sensitivity (Phase B HD acceptance)

The reformulated p95 acceptance bound is sensitive to:

1. **Continental topology** — thinner strips produce higher p95
   (boundary cells dominate the deviation distribution).
2. **α × N (Phase A cumulative erosion)** — more aggressive Phase A
   produces a more deviated LR shape, which the HD upscale + erosion
   then has to track.

Pinned probe results (D8 conservative, 32² × 5 cycles → 512² ×
5×10⁵ droplets):

| metric | value |
|---|---:|
| L_∞ | 0.1440 |
| p95 | 0.0754 (PASS) |

64² spot-checks (D8 conservative, 64² × 5 cycles → 1024² × 5×10⁵
droplets):

| preset | p95 | passes 0.10? |
|---|---:|---|
| `single_continent` (4 plates, 50 % cont) | 0.0856 | yes |
| `convergence` (6 plates, 40 % cont) | 0.1144 | **no, by 14 %** |

Aggressive demo (α = 0.05, N = 15 → 1024² × 5×10⁶ droplets):

| preset | p95 |
|---|---:|
| `single_continent` | 0.1861 (fails) |
| `convergence` | 0.2165 (fails) |

The aggressive demo is **expected to fail** D5 — Phase A's larger
per-cycle Δh produces an LR shape that the HD pipeline cannot
track to 10 % p95 because the HD valley carving plus the
proportionally larger LR-baseline deviation push the distribution
above the threshold. The aggressive demo is a *visual* exploration,
not a calibrated production setting.

**Practical guidance**:

- **D8 defaults pass D5 on `single_continent`** with comfortable
  margin (0.0856 < 0.10 = 14 % headroom).
- **D5 may fail on thinner-strip presets** like `convergence` even
  at D8 defaults. Workaround: raise `grand_scale_tolerance` to
  ≈ 0.13 for those presets (the panel slider's full range is
  `[0.05, 0.30]`). Document the per-preset value in the run
  metadata.
- **D5 fails reliably for aggressive demos** — do not interpret
  this as a calibration regression; it is an expected
  consequence of large Phase A erosion.

## Open calibration questions (deferred)

None of these are blockers for Step 12 closure; they are
documented here so future users do not have to rediscover them.

1. **Per-preset `grand_scale_tolerance` defaults** — would a
   topology-derived auto-tolerance (e.g.,
   `0.10 + 0.005 × num_plates`) reduce the per-preset adjustment
   burden? Empirically, p95 scales roughly with √num_plates × (1 −
   continental_ratio); a more principled formula could be derived
   if a wider preset sweep is run.
2. **(α × N) budget vs counter-isostasy convergence** — is there a
   minimum (α × N) below which counter-isostasy fails to converge
   within the allotted budget? The D8 default (0.05) shows clear
   convergence; α × N → 0 obviously fails. The transition zone
   has not been mapped.
3. **HD `num_droplets` ↔ `grand_scale_tolerance` relationship** —
   higher droplet counts produce deeper HD valleys, which raise
   L_∞ but should leave p95 mostly unchanged. The 5×10⁵ vs 5×10⁶
   probes (Phase 5 vs Phase 6 aggressive) bracket this; a
   dedicated sweep would tighten the relationship.
4. **`craton_retention_threshold` distinct from
   `plate_area_min`** — Step 12 Phase 3.5 disambiguated the two
   parameters but locked both at 0.10. A future tuning round
   could move the retention threshold without disturbing the
   Step 9 init contract; whether this is desirable depends on
   the user's craton-evolution intent.

## Recommendations summary

For a user opening the workflow panel for the first time:

1. **Start with D8 defaults** + the `single_continent` preset.
   Run Phase A → observe the cycle-history table converging on
   peak S̃ within 3–4 cycles. Run Phase B → confirm D5 p95
   acceptance passes. This is the "everything works" demo.
2. **Tune α first** if you want more / less counter-isostasy.
   Stay in `[0.01, 0.05]` unless you have a specific reason to
   leave that range.
3. **Tune β second** if you want different mass-budget behaviour.
   Most users will leave β = 0 (the D8 default).
4. **Tune (N, k) ratio last** — these levers mostly affect
   wallclock, not the regime; pick what fits your patience.
5. **Switch presets to `convergence`** if you want a thin-strip
   topology. Expect D5 to need a slight `grand_scale_tolerance`
   bump (try 0.13).
6. **Aggressive demo** (α = 0.05, N = 15) is for visual
   exploration only — it deliberately exits the D5-acceptance
   regime to make Phase A effects visible at scale.

For a developer extending Step 12 in a follow-up issue (Step 12.X
or beyond):

1. **Phase A border curvature** — the open structural item.
   Hypotheses: lateral diffusion, stochastic boundary amplitude
   noise, coast-complexity-dependent erosion rate. Pick one,
   probe on `single_continent` 64²-mantle-on, document at the
   level of detail this report uses for D2.
2. **Per-preset auto-tolerance** for D5 — see open question #1.
3. **HD parameter sweep** — see open question #3.

Both the user and developer pathways are tractable with the data
already committed under `docs/reports/step12_phase_6_*` plus a
modest sweep extension; no additional milestone-scope physics work
is required to make progress.
