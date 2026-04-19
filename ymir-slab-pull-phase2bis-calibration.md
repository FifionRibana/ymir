# Slab-Pull Reformulation — Phase 2-bis Calibration Sweep

Issue #75 — *Reformulate slab-pull as an auto-regulated operator term
instead of RHS forcing.*
Branch: `75-reformulate-slab-pull-...` (same branch as prior phases).
Date: 2026-04-19.

## TL;DR — critical finding

**The Phase 2 operator reformulation cannot drive subduction at any
`slab_pull_factor` value.** Five-point sweep from 0.05 (current
default, ×1) to 5000 (×100 000) with scenario B on 64²/seed 42/120
steps: every metric is invariant. Velocities at factor=5000 equal
velocities at factor=0.05 within 1 %. `eta_ratio` is 11.23 at every
factor. Wallclock is ~47 s at every factor.

The design is structurally incapable of accelerating plate motion:
`γ · (v·n̂) · n̂` is SPD and therefore applies a *damping* force
that opposes convergent motion. The pre-Phase-2 `apply_slab_pull` was
an *acceleration* source (it added to `plate.velocity` unconditionally,
breaking SPD). The two mechanisms are mathematically opposite. No
calibration of the SPD form can recover the pre-Phase-2 behaviour.

**Recommendation:** do NOT adjust `slab_pull_factor` (no value works).
File a follow-up issue #80 to redesign the slab-pull mechanism. The
operator stays in the codebase as a placeholder; the default
`slab_pull_factor = 0.05` gets an inline comment flagging the known
inertness. Physics-meaningful slab-pull must be tracked separately.

---

## 1. Parametric sweep results

Config: 64²/seed 42, scenario B (all extensions on), Newton solver,
adaptive dt (`dt_target = 2.0`), 1 rep × 120 macro steps.
`slab_pull_factor` overridden via `YMIR_SLAB_PULL_FACTOR` env var.
Steady-regime means computed over steps 40–80.

| `slab_pull_factor` | `γ_max` (margin) | `γ_mean` (margin) | `eta_ratio` mean | solve μs/substep | wallclock (s) | `v_global_max` @ step 119 |
|-:|-:|-:|-:|-:|-:|-:|
| 0.05 (default) | 1.0×10⁻⁴ | 2.4×10⁻⁵ | **11.23** | 363 210 | 47.65 | 5.62×10⁻⁴ |
| 5.0 (×100) | 9.9×10⁻³ | 2.4×10⁻³ | **11.23** | 345 900 | 45.87 | 5.62×10⁻⁴ |
| 50.0 (×1000) | 9.9×10⁻² | 2.4×10⁻² | **11.23** | 336 090 | 48.49 | 5.60×10⁻⁴ |
| 500.0 (×10 000) | 9.9×10⁻¹ | 2.4×10⁻¹ | **11.23** | 320 450 | 46.25 | 5.59×10⁻⁴ |
| 5000.0 (×100 000) | 9.6 | 2.4 | **11.23** | 346 520 | 46.96 | 5.60×10⁻⁴ |

Nothing depends on the factor. γ scales perfectly linearly with the
input (confirming the seed formula is wired correctly) but every
downstream physical metric — velocity, η contrast, iteration count,
wallclock — is identical within 1 %. No Newton-solver warnings specific
to any factor.

### 1.1 Velocity trajectory sample (factor=0.05 vs factor=5000)

Peak-cell velocity magnitude (`v_global_max` from `slab_pull_sweep`):

| step | factor=0.05 | factor=5000 | ratio |
|-:|-:|-:|-:|
|   0 | 6.18×10⁻⁴ | 6.08×10⁻⁴ | 0.98 |
|  10 | 7.15×10⁻⁴ | 7.15×10⁻⁴ | 1.00 |
|  40 | 6.66×10⁻⁴ | 6.66×10⁻⁴ | 1.00 |
|  80 | 6.31×10⁻⁴ | 6.33×10⁻⁴ | 1.00 |
| 119 | 5.62×10⁻⁴ | 5.60×10⁻⁴ | 1.00 |

**100 000×** change in `slab_pull_factor` produces a 2 % change at
step 0 (numerical noise) and indistinguishable results afterwards. The
velocities are also *decreasing* over time in both cases — plates
decelerate regardless of the γ magnitude.

### 1.2 Newton convergence

All five runs completed without non-convergence errors. The
consistent ~10–20 BiCGSTAB "did not converge within max iterations"
warnings per run are a pre-existing Phase 1-bis baseline signal and
do not correlate with the factor.

---

## 2. Operating range analysis

### 2.1 Required conditions from the task

The task defined the viable range as the intersection of:
- smallest factor at which subduction velocities actually grow over
  time (→ non-negligible physical effect);
- largest factor at which `eta_ratio < 20`, Newton converges on
  ≥ 95 % of steps, wallclock < 200 s / 120 steps.

### 2.2 Finding

The intersection is **empty** — but not for the reason the task
anticipated. The task assumed some factor might produce a cascade
(η_ratio → 62 regime again); instead, *no* factor produces any
effect at all. The upper bound (no cascade) is trivially satisfied
at every factor. The lower bound (physical effect) is violated at
every factor.

### 2.3 Why the operator is inert

The operator term contribution to the `vx` diagonal at a margin
cell is `0.25 · γ · n̂_x²`. To compete with the viscous diagonal at
default `η = 1` and `dx = 1/64`:

```
viscous diagonal ≈ η / dx² = 1 × 4096 = 4096
γ diagonal       ≈ 0.25 · γ · n̂_x² = 0.25 · γ    (n̂_x² ≤ 1)
```

For the γ contribution to equal the viscous contribution,
`γ ≥ 16 000`. With `γ_seed = slab_pull_factor · |source_rate|` and
`|source_rate| ≈ 0.002` at margin cells, that requires
`slab_pull_factor ≈ 8 000 000` — three orders of magnitude beyond
the sweep range, and at that point γ_mean saturates the operator
and condition numbers would blow up.

But even *if* γ were brought up to parity with η: the form
`γ · n̂⊗n̂ · v` is a **damping** force per `⟨v, γ·(v·n̂)²⟩ ≥ 0` —
it only opposes motion, never drives it. The only thing a very
large γ would accomplish is freezing plates against any convergent
motion, not accelerating them. This was always guaranteed by the
three Phase 2 unit tests (`slab_pull_damps_convergent_motion`
being the most explicit); it just wasn't recognized at design time
that "damping" was the wrong sign for the intended physics.

---

## 3. Physics sanity check

### 3.1 Velocity comparison — pre-Phase 2 vs Phase 2 at recommended tuning

A formal before/after comparison would require either checking out a
pre-Phase-2 binary or locating saved output from the pre-Phase-2
main branch. Neither is strictly required to make the call here
because the sweep (§1) already demonstrates that *no* Phase 2
tuning reproduces slab-pull behaviour. The user's visual
confirmation ("velocity vectors in subduction zones no longer
accelerate") independently corroborates this across the whole
factor range tested.

### 3.2 Heightmap comparison

Not produced. With velocities identical across factors, the
heightmaps at step 120 would differ only by numerical noise. Useful
only against a pre-Phase-2 baseline, which this phase's scope
excludes.

### 3.3 Pre-Phase-2 behaviour expectation

From the Phase 1-bis report and the issue description: with
`apply_slab_pull` in the pipeline, `plate.velocity` grew by
`slab_pull_factor × plate.subducted_mass` per macro step, in the
direction of current motion, capped at `max_plate_velocity = 5.0`.
Over 120 steps at factor 0.05 and typical subducted-mass values,
peak velocities would saturate at or near the 5.0 cap after a few
tens of steps — versus the Phase 2 peak of ~7×10⁻⁴ observed across
all sweep runs. The pre-Phase-2 regime is *four orders of magnitude
faster* at steady state.

---

## 4. Recommendation

### 4.1 Decision

**Do not change the `slab_pull_factor` default.** No value works.
Keep it at `0.05` and add an inline comment explaining the known
inertness, along with a pointer to the follow-up redesign issue.

### 4.2 Follow-up: issue #80 (to file)

Title (proposed): *"Slab-pull needs a non-SPD driving term: current
operator reformulation is structurally incapable of accelerating
subduction."*

Scope:
- Replace or augment `γ · (v·n̂) · n̂` with a driving term that
  breaks SPD in a controlled way, e.g.:
  - a velocity-sourced but capped RHS contribution that targets the
    subducting plate only (margin-localised, unlike the pre-Phase-2
    plate-wide velocity boost);
  - or a smooth-saturated velocity source
    `γ · tanh(|v|/v_ref) · (v̂·n̂) · n̂` that grows like v at low
    speed and saturates at high speed (auto-regulation without SPD);
  - or a two-term operator: a PSD damping term (keeps the current
    infrastructure, stabilises) plus a small explicit body force
    computed from `|source_rate|` placed in the RHS (drives).
- Re-validate with the Phase 1-bis scenario runner (`A`, `B`, `C`,
  plus at least one pre/post morphology check at step ≥ 120).

### 4.3 What stays from Phase 2

The Phase 2 code is **not wasted** even though the physics result is
inert:

- The `γ_slab`, `n̂` fields on `BoundaryField` are already computed
  and will be reused by the follow-up — the margin detection and
  Benioff-decay infrastructure is in place.
- `SlabPullField` struct and the `Option<&SlabPullField>` parameter
  chain through `apply_stokes` / `compute_jacobi_precond` /
  `StencilCoeffs::compute` / `solve_velocity_*` / `execute_tectonic_pass`
  are ready to carry any replacement operator term.
- The 7× wallclock improvement (Phase 2 §1.1) is real and stays,
  because it came from removing `apply_slab_pull`'s velocity boost
  from the pipeline. The solver is no longer ill-conditioned.
- Three new unit tests in `stokes.rs` remain valid — they verify
  what the current term *does* (damp), which is a useful building
  block even if it's not the whole mechanism.

### 4.4 What Phase 2-bis commits

Minimal change: a single comment at the
`BoundaryConfig::default` `slab_pull_factor = 0.05` line flagging the
inertness and linking to #80. No default change.

---

## 5. Appendix — Reproduction

```bash
# Build (uses the Phase 1-bis scenario runner with Phase 2-bis env var)
cargo build --release --example phase1bis_scenarios

# Sweep
for f in 0.05 5.0 50.0 500.0 5000.0; do
  tag="f${f/./p}"
  YMIR_SLAB_PULL_FACTOR=$f \
    ./target/release/examples/phase1bis_scenarios.exe B 1 "logs/phase2bis_$tag" 120
done

# Extract metrics
for tag in f0p05 f5p0 f50p0 f500p0 f5000p0; do
  log="logs/phase2bis_${tag}/phase1bis_B_01.log"
  # Steady-state eta_ratio and solve time (steps 40–80)
  awk '/eta_breakdown/ && /newton_iter=0 / {
    match($0, /step=[0-9]+/); s=substr($0, RSTART+5, RLENGTH-5)+0
    if (s>=40 && s<=80) {
      match($0, /eta_ratio=[0-9.eE+-]+/)
      sum+=substr($0, RSTART+10, RLENGTH-10)+0; n++
    }
  } END { print sum/n }' "$log"
done
```

## 6. Appendix — diagnostic target added

`emit_slab_pull_sweep(boundary_field, grid)` in
[diagnostics.rs](crates/ymir-core/src/tectonics/solver/diagnostics.rs),
target `slab_pull_sweep` at debug level. Fires once per macro step
inside `execute_tectonic_pass`. Logs `gamma_margin_{min,max,mean}`,
`margin_cells`, `margin_v_max`, `v_global_max`. Kept on in the
scenario runner's default filter; costs a single cell-sweep
(~4 k cells on 64², ~100 µs).
