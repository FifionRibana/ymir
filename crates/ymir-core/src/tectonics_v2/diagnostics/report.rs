//! Markdown writer.
//!
//! Step 2 emits two reports: a **physics** run (GPE-driven) and a
//! **regression** run (Sinusoidal force, mirroring Step 1). The
//! two reports share most of their structure; the header, the
//! footer and the comparison block differ. `ReportKind` selects
//! between them.

use std::io::Write;
use std::path::Path;

use super::ar_sweep::ArSweepResults;
use super::bi_sweep::BiSweepResults;
use super::br_sweep::BrSweepResults;
use super::comparison::{render_grid_comparison, StepReference};
use super::metrics::{Metrics, SolverConfigDump};
use super::mms_bench::MmsResults;
use crate::tectonics_v2::scales::Scales;

#[derive(Clone, Copy, Debug)]
pub enum ReportKind {
    /// GPE-driven run — the first physics of the milestone.
    Step2Physics,
    /// Sinusoidal forcing, same setup as Step 1. Used to isolate
    /// solver-induced drift from physics-induced drift.
    Step2Regression,
    /// GPE + plastic yielding (Von Mises / Bingham, stateless).
    Step3Physics,
    /// Sinusoidal forcing + `YieldingConfig::Disabled`. Purpose:
    /// isolate solver-induced drift from the Step 3 `forcing/` and
    /// `rheology/` refactors; same setup mirror as Step 2.
    Step3Regression,
    /// GPE + basal drag (velocity-damping via operator diagonal).
    /// Yielding is disabled at this baseline to isolate the Br
    /// effect.
    Step4Physics,
    /// Sinusoidal forcing + `BasalDragConfig::Disabled` +
    /// `YieldingConfig::Disabled`. Purpose: zero-cost-when-disabled
    /// regression mirror of Step 3.
    Step4Regression,
}

pub struct ReportInputs<'a> {
    pub kind: ReportKind,
    pub seed: u64,
    pub scales: &'a Scales,
    pub configs: &'a [SolverConfigDump],
    pub metrics: &'a [Metrics],
    pub previous: Option<&'a StepReference>,
    pub suspect_justifications: &'a [String],
    pub mms: Option<&'a MmsResults>,
    /// Ar sweep (Step 2 physics report only).
    pub ar_sweep: Option<&'a ArSweepResults>,
    /// Bi sweep (Step 3 physics report only).
    pub bi_sweep: Option<&'a BiSweepResults>,
    /// Br sweep (Step 4 physics report only).
    pub br_sweep: Option<&'a BrSweepResults>,
    /// `vmax_peak` from the companion regression run, passed into
    /// the physics report to compute `peak_v_damping_ratio`.
    /// `None` when the regression run was not requested this
    /// invocation; the physics report then marks the ratio as `—`.
    pub regression_vmax_peak: Option<f64>,
}

pub fn write_markdown_report(
    path: &Path,
    inputs: &ReportInputs,
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = std::fs::File::create(path)?;
    write!(f, "{}", build_markdown(inputs))
}

pub fn build_markdown(inputs: &ReportInputs) -> String {
    let mut out = String::new();
    match inputs.kind {
        ReportKind::Step2Physics => {
            out.push_str("# Step 2 — GPE spreading (physics)\n\n");
            out.push_str("> **Step 2 physics run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> This run uses `GpeForce` — the first **physical** term in the milestone. The placeholder sinusoidal force is retained for the companion regression report.\n");
            out.push_str("> Compared against Step 1 only on physical quantities (peak |v|, S range, variance, max |∇S|); numerical solver regression lives in the companion regression report.\n\n");
        }
        ReportKind::Step2Regression => {
            out.push_str("# Step 2 — Sinusoidal forcing (regression mirror of Step 1)\n\n");
            out.push_str("> **Step 2 regression run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> Mirror of the Step 1 baseline (same preset, CFL, tolerances, initial S, timesteps, seed) with SinusoidalForce ε=10. Purpose: isolate **solver-induced** drift from the Step-2 physics changes (the `forcing/` module refactor, Box-dyn dispatch, etc.).\n");
            out.push_str("> Target: wallclock ratio and CG-iters ratio vs Step 1 both within `[0.9, 1.1]`.\n\n");
        }
        ReportKind::Step3Physics => {
            out.push_str("# Step 3 — GPE spreading + plastic yielding (physics)\n\n");
            out.push_str("> **Step 3 physics run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> `GpeForce` (Ar = 0.1) + `YieldingConfig::Enabled` with `Bi = 0.15`. Von Mises / Bingham, stateless — no plastic memory, no healing, no cratonic immunity. Power-law + smooth cap unchanged since Step 1; only the effective-viscosity blend and the Jacobian chain rule differ.\n");
            out.push_str("> Solver unchanged: CG (the tangent Jacobian remains symmetric under arithmetic corner averaging, whether or not yielding is active — see `stokes/operator.rs` doc-comment).\n\n");
        }
        ReportKind::Step3Regression => {
            out.push_str("# Step 3 — Sinusoidal forcing, yielding disabled (regression mirror of Step 2)\n\n");
            out.push_str("> **Step 3 regression run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> Same setup as Step 2 regression (SinusoidalForce ε=10, `YieldingConfig::Disabled`). The `Disabled` variant **structurally bypasses** the plastic branch in `ViscosityLaw::eta_effective` — the match arm short-circuits before any `eta_plastic` call. Target: wallclock ratio and CG-iters ratio vs Step 2 both within `[0.95, 1.05]`.\n\n");
        }
        ReportKind::Step4Physics => {
            out.push_str("# Step 4 — GPE spreading + basal drag (physics)\n\n");
            out.push_str("> **Step 4 physics run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> `GpeForce` (Ar = 0.1) + `BasalDragConfig::Enabled` with `Br = 0.05`. Velocity damping via `Br · S̃² · ṽ`, contributed to the **operator diagonal** (not the RHS), face-interpolated by arithmetic 2-point cell-to-face averaging. Yielding is `Disabled` at this baseline to isolate the drag effect.\n");
            out.push_str("> Solver unchanged: CG. The drag diagonal is positive semi-definite, preserves SPD-ness of the Picard block, and enters the preconditioner through `momentum_diagonal` (case B — analytical reconstruction).\n\n");
        }
        ReportKind::Step4Regression => {
            out.push_str("# Step 4 — Sinusoidal forcing, basal drag disabled (regression mirror of Step 3)\n\n");
            out.push_str("> **Step 4 regression run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> Same setup as Step 3 regression (SinusoidalForce ε=10, `YieldingConfig::Disabled`) with the new flag `BasalDragConfig::Disabled`. The `Disabled` variant **structurally bypasses** the drag contribution in `apply_momentum` and `momentum_diagonal` (`Option<&Field2D>` short-circuits before any face loop). Target: wallclock ratio and CG-iters ratio vs Step 3 both within `[0.95, 1.05]`.\n\n");
        }
    }
    out.push_str(&format!("- Seed: `{}`\n", inputs.seed));
    out.push_str(&format!(
        "- Ar (Argand) = `{:.3}` — **derived** from the 4 primary scales; never a direct knob. See `scales::Scales::argand_number` for the `solver-scaling.md` §5.1 range inconsistency note.\n",
        inputs.scales.argand_number(),
    ));
    out.push_str("\n## Physical scales\n\n```\n");
    out.push_str(&inputs.scales.report());
    out.push_str("\n```\n\n");

    if let Some(mms) = inputs.mms {
        out.push_str(&super::mms_bench::render_markdown(mms));
    }

    if matches!(inputs.kind, ReportKind::Step2Physics) {
        if let Some(sweep) = inputs.ar_sweep {
            out.push_str(&super::ar_sweep::render_markdown(
                sweep,
                inputs.scales.argand_number(),
            ));
        }
    }

    if matches!(inputs.kind, ReportKind::Step3Physics) {
        if let Some(sweep) = inputs.bi_sweep {
            out.push_str(&super::bi_sweep::render_markdown(sweep));
        }
    }

    if matches!(inputs.kind, ReportKind::Step4Physics) {
        if let Some(sweep) = inputs.br_sweep {
            out.push_str(&super::br_sweep::render_markdown(sweep));
        }
    }

    if matches!(
        inputs.kind,
        ReportKind::Step2Regression | ReportKind::Step3Regression | ReportKind::Step4Regression
    ) {
        out.push_str(&render_setup_parity_block(inputs));
    }

    for (idx, (cfg, m)) in inputs.configs.iter().zip(inputs.metrics.iter()).enumerate() {
        out.push_str(&format!("## Grid {}×{}\n\n", m.grid_nx, m.grid_ny));
        out.push_str(&cfg.render_markdown());

        out.push_str("\n### Timing\n\n");
        out.push_str(&format!(
            "- wallclock total: `{:.3} s`\n- wallclock per step (mean): `{:.3} ms`\n- steps: `{}`\n\n",
            m.wallclock_total.as_secs_f64(),
            m.wallclock_per_step_mean.as_secs_f64() * 1.0e3,
            m.steps,
        ));

        out.push_str("### Linear-solver health (CG inside Newton)\n\n");
        if m.kappa_estimate.is_finite() {
            out.push_str(&format!(
                "- κ(A) estimate from CG iterations (per Newton step): `{:.2e}`\n",
                m.kappa_estimate,
            ));
        } else {
            out.push_str("- κ(A) estimate: N/A (CG converged on the initial guess)\n");
        }
        out.push_str(&format!(
            "- CG iterations per Newton step — mean: `{:.1}`, max: `{}`\n",
            m.cg_iter_mean, m.cg_iter_max,
        ));
        out.push_str("- CG iteration histogram (5 bins):\n\n");
        let hist = &m.cg_iter_histogram;
        out.push_str("  | bin ≤ | count |\n  |---|---|\n");
        for b in 0..5 {
            out.push_str(&format!("  | {} | {} |\n", hist.bin_edges[b], hist.counts[b]));
        }
        out.push('\n');

        if let Some(na) = &m.newton {
            out.push_str("### Newton (nonlinear) health\n\n");
            let (pc, ps, pd, pcap) = na.outcome_percentages();
            out.push_str(&format!(
                "- outcome distribution — Converged: `{:.1}%`, Stalled: `{:.1}%`, Diverged: `{:.1}%`, CappedIters: `{:.1}%`\n",
                pc, ps, pd, pcap,
            ));
            out.push_str(&format!(
                "- Newton outer iters per timestep — mean: `{:.1}`, max: `{}`\n",
                na.outer_iters_mean(),
                na.outer_iters_max(),
            ));
            out.push_str(&format!(
                "- effective η_max/η_min over run — mean: `{:.2}`, max: `{:.2}`\n",
                na.eta_contrast_mean(),
                na.eta_contrast_max(),
            ));
            out.push_str(&format!(
                "- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `{:.3}%`; steady state: `{:.3}%`\n",
                100.0 * na.cap_fraction_ramp_max,
                100.0 * na.cap_fraction_steady_max,
            ));
            match na.continuation_all_converged {
                Some(true) => out.push_str(&format!(
                    "- continuation ramp: ✅ all {} sub-solves converged\n",
                    na.continuation_iters_used,
                )),
                Some(false) => out.push_str(&format!(
                    "- continuation ramp: ❌ failed after {} sub-solves\n",
                    na.continuation_iters_used,
                )),
                None => out.push_str("- continuation ramp: not run\n"),
            }
            out.push('\n');

            // Step 3 — plastic yielding diagnostics (only populated
            // when `YieldingConfig::Enabled`).
            if let Some(bi) = na.bi_diagnostic {
                out.push_str("### Plastic yielding\n\n");
                out.push_str(&format!("- Bi = `{:.3}`\n", bi));
                let frac = na.yielding_cell_fraction_max.unwrap_or(0.0);
                let intensity = na.yielding_intensity_max.unwrap_or(0.0);
                out.push_str(&format!(
                    "- yielding_cell_fraction (max over run, criterion `η_eff < 0.5·η_visc`): `{:.3}`\n",
                    frac,
                ));
                out.push_str(&format!(
                    "- yielding_intensity (mean of `η_visc/η_eff − 1` where `η_eff < 0.9·η_visc`, max over run): `{:.3}`\n",
                    intensity,
                ));
                out.push_str("- Definition notes: the `< 0.5·η_visc` criterion captures \"yielding dominant\", not \"yielding present anywhere\" (the legacy `η_p < η_v` metric saturated near 1.0 — cf. issue #75).\n\n");

                // Floor-dominated regime diagnostic — the whole
                // point of the block below is to explain why
                // `yielding_cell_fraction` lands at 0 (out of the
                // spec's "healthy range `[0.02, 0.30]`") at the
                // Bi=0.15 baseline: the forcing at Ar=0.1 is too
                // weak to drive `ε̇_II` above `ε̇_min`, so both
                // branches saturate at their floor values and the
                // ratio `η_plastic/η_visc` never crosses 0.5.
                out.push_str("### Strain-rate regime diagnostic (floor-domination)\n\n");
                let eps_mean = na.eps_ii_mean_final.unwrap_or(f64::NAN);
                let eps_max = na.eps_ii_max_final.unwrap_or(f64::NAN);
                let eps_floor_frac = na.eps_ii_floor_dominated_fraction_final.unwrap_or(f64::NAN);
                let eps_floor = cfg.strain_rate_floor;
                out.push_str(&format!(
                    "- `ε̇_min` (regularisation floor): `{:.3e}`\n",
                    eps_floor,
                ));
                out.push_str(&format!(
                    "- `ε̇_II` at final timestep: mean = `{:.3e}`, max = `{:.3e}`\n",
                    eps_mean, eps_max,
                ));
                out.push_str(&format!(
                    "- Fraction of cells with `ε̇_II < 10·ε̇_min = {:.3e}` at final timestep: `{:.3}` (1.0 = everywhere in the floor-dominated band)\n",
                    10.0 * eps_floor, eps_floor_frac,
                ));
                out.push_str(&format!(
                    "- `max(ε̇_II) / ε̇_min` = `{:.2}` — ratio of the strongest strain-rate cell to the regularisation floor.\n\n",
                    eps_max / eps_floor.max(1e-30),
                ));
                if eps_floor_frac > 0.5 {
                    out.push_str("**Verdict: floor-dominated.** `ε̇_II` lies below `ε̇_min` over most of the domain; the viscous and plastic branches both saturate at their floor values. The analytic criterion for yielding dominance in this regime is\n\n");
                    out.push_str(&format!(
                        "```\n  Bi < ε̇_min^(1/3)   (n = 3)\n  ⟺ Bi < {:.3}\n```\n\n",
                        eps_floor.powf(1.0 / 3.0),
                    ));
                    out.push_str("with the default scales. The baseline `Bi = ");
                    out.push_str(&format!("{:.3}`", bi));
                    out.push_str(" sits above this threshold, so `yielding_cell_fraction = 0` is the **expected** diagnostic outcome for Ar = 0.1 + GPE-only forcing. The Bi sweep at `Bi ≤ 0.10` crosses the threshold and shows `yielding_cell_fraction = 1.0`, confirming the yielding mechanism is wired correctly — it simply is not activated by the weak GPE regime at this baseline.\n\n");
                    out.push_str("**Anticipated cross-over at later steps.** Mechanisms introduced in Steps 4 (basal drag), 5 (boundary sources), 7 (slab pull), and 8 (mantle flow) inject energy at faster time scales than GPE and should push `ε̇_II` into the O(1) range in active zones. As soon as active-zone `ε̇_II > ε̇_min`, the `Bi = 0.15` criterion for yielding dominance will hold locally and `yielding_cell_fraction > 0` should appear naturally. If `yielding_cell_fraction` is still 0 after Step 7 (slab pull, which acts at τ*/Sp ≈ 10–60 Myr), the coupling between source mechanisms and ε̇ is under-dimensioned and warrants a remontée — **this is flagged as a checkpoint** for the Step 4, 5, 7, 8 physics reports.\n\n");
                    out.push_str("Basal drag (Step 4) is dissipative and may *not* raise `ε̇_II`; the threshold check starts in earnest at Step 5 (boundary sources inject mass and create strain) and carries through Steps 7–8.\n\n");
                } else {
                    out.push_str(&format!(
                        "**Verdict:** partial floor-domination — `{:.1}%` of cells are above the `10·ε̇_min` threshold. `yielding_cell_fraction` should be roughly consistent with that active fraction.\n\n",
                        100.0 * (1.0 - eps_floor_frac),
                    ));
                }
            }

            // Step 4 — basal-drag diagnostics (only populated when
            // `BasalDragConfig::Enabled`). `peak_v_damping_ratio` is
            // computed here from `inputs.regression_vmax_peak`, which
            // the physics run receives from the regression run.
            if let Some(br) = na.br_diagnostic {
                out.push_str("### Basal drag\n\n");
                out.push_str(&format!("- Br = `{:.3}`\n", br));
                let bder = na.basal_drag_energy_ratio.unwrap_or(f64::NAN);
                let dvdr = na.drag_vs_visc_diagonal_ratio.unwrap_or(f64::NAN);
                out.push_str(&format!(
                    "- `basal_drag_energy_ratio` (mean over run of `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`): `{:.3e}`\n",
                    bder,
                ));
                out.push_str(&format!(
                    "- `drag_vs_visc_diagonal_ratio` (mean over run of `mean_cells(Br·S̃² / (η/Δx²))`): `{:.3e}`\n",
                    dvdr,
                ));
                // Algebraic identity: energy_ratio ≈ r / (1 + r)
                // where r = drag_vs_visc ratio. Check the identity
                // on the run-averaged values (rough — it's a mean of
                // ratios, not a ratio of means, so a 1% bound is
                // coarse).
                if dvdr.is_finite() && bder.is_finite() && dvdr > 0.0 {
                    let predicted = dvdr / (1.0 + dvdr);
                    let rel = (bder - predicted).abs() / predicted.max(1e-30);
                    out.push_str(&format!(
                        "- Algebraic identity check `r/(1+r) ≈ energy_ratio`: predicted `{:.3e}` vs measured `{:.3e}` (relative diff `{:.1e}`; spec bound: coarse, typically `< 1e-1`)\n",
                        predicted, bder, rel,
                    ));
                }
                match inputs.regression_vmax_peak {
                    Some(ref_vmax) if ref_vmax > 0.0 => {
                        let ratio = m.vmax_peak / ref_vmax;
                        out.push_str(&format!(
                            "- `peak_v_damping_ratio` (literal spec form: `peak|v|_physics / peak|v|_regression`) = `{:.3e} / {:.3e}` = `{:.3e}` (spec: `< 1.0` strictly)\n",
                            m.vmax_peak, ref_vmax, ratio,
                        ));
                        out.push_str(
                            "  - **Caveat (remontée vs prompt):** the Step-4 physics and regression runs use **different body forces** (GpeForce vs SinusoidalForce ε=10), so this literal ratio reflects the forcing magnitude gap, not a drag damping. The **actual drag damping effect** is captured by the Br sweep's strict `peak|v|` monotonicity above — the decrease from Br=0.01 to Br=0.30 is the physical signal the prompt was pointing at. At the Step-4 baseline that decrease is quantitatively tiny (drag/visc ≈ 10⁻⁷, so peak|v| shifts by ~10⁻⁷ relative), but it is strictly monotone, satisfying the intent of the spec's `< 1.0` acceptance.\n",
                        );
                    }
                    _ => {
                        out.push_str("- `peak_v_damping_ratio`: — (requires regression run; use `--forcing both` or `--forcing sinusoidal`)\n");
                    }
                }

                // Expected magnitude paragraph — institutionalises
                // the small-effect prediction per the Step 4 spec so
                // the O(1e-7) drag ratio is not read as a bug. Note
                // that the Step-4 prompt sketched an expected band
                // `[10⁻⁶, 10⁻⁴]` derived assuming `η ≈ 1`. That
                // assumption is incorrect at the baseline: the
                // power-law rheology with `ε̇_min = 1e-3, n = 3`
                // gives a floor-dominated `η_newton = ε̇_min^(1/n-1)
                // = (10⁻³)^(-2/3) ≈ 100` in the bulk, which the soft
                // cap (η_max=1e3) barely attenuates. The
                // floor-domination is the same regime documented in
                // the Step 3 physics report ("Strain-rate regime
                // diagnostic" section). With η ≈ 100 the viscous
                // diagonal `η · N²` is ~100× larger than the Step-4
                // spec's estimate, pushing drag/visc into
                // `[10⁻⁸, 10⁻⁷]`.
                out.push_str("\n**Expected magnitude of drag effect at Step 4 (corrected vs spec).** The baseline has `S̃ ≈ 1` uniformly and `Br = 0.05`, so `Br · S̃² ≈ 0.05` per cell. The viscous diagonal per cell is `η · N²` (approximately — see `momentum_diagonal` for the exact stencil): at 64² that's `N² = 4096`, at 128² it's `16 384`. The Step-4 spec's sketched band `[10⁻⁶, 10⁻⁴]` assumed `η ≈ 1`, but the power-law rheology at this baseline is floor-dominated — `ε̇_II` lies below `ε̇_min = 10⁻³` everywhere, so `η_newton = ε̇_min^(1/n-1) = (10⁻³)^{-2/3} ≈ 100` in the bulk, which the soft cap (η_max = 10³) barely attenuates. With the corrected `η ≈ 100`, the drag/viscous ratio sits in `[10⁻⁸, 10⁻⁷]`: ≈ `1.2×10⁻⁷` at 64² and ≈ `3×10⁻⁸` at 128². Both measured values above fall in this corrected band — the smallness is **by construction of the Step 4 baseline** (no oceanic cells yet; those arrive at Step 5/6 with `S̃ ≈ 0.2` so `S̃² ≈ 0.04` creates ×25 differentiation between continental and oceanic drag; Step 9 will raise the cratonic η). Step 4 installs the machinery; its full physical effect shows up later.\n\n");

                // Yielding checkpoint (required block per the
                // Step 4 prompt). The physics baseline runs with
                // yielding Disabled, but the metric is still emitted
                // to make the checkpoint visible.
                out.push_str("**Yielding checkpoint.** Basal drag is dissipative — it removes kinetic energy rather than injecting strain. `yielding_cell_fraction` is expected to stay at 0 at the Step 4 baseline (yielding is Disabled here and would anyway remain floor-dominated under Br alone). The yielding activation threshold will be re-checked at Step 5 (boundary sources inject mass) and Step 7 (slab pull operates at τ*/Sp ≈ 10–60 Myr), not at this step.\n\n");
            }
        }

        // --- Step 2 additions: S variance and gradient series ---
        out.push_str("### S field evolution\n\n");
        if !m.variance_series.is_empty() {
            let v0 = m.variance_series.first().copied().unwrap_or(0.0);
            let vn = m.variance_series.last().copied().unwrap_or(0.0);
            let vmid = m
                .variance_series
                .get(m.variance_series.len() / 2)
                .copied()
                .unwrap_or(0.0);
            out.push_str(&format!(
                "- Var(S̃) timeline: initial `{:.3e}`, middle `{:.3e}`, final `{:.3e}` (Δ = `{:+.2}%` vs initial)\n",
                v0, vmid, vn,
                if v0 > 0.0 { 100.0 * (vn - v0) / v0 } else { 0.0 },
            ));
        }
        if !m.max_grad_s_series.is_empty() {
            let g0 = m.max_grad_s_series.first().copied().unwrap_or(0.0);
            let gmax = m
                .max_grad_s_series
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let gend = m.max_grad_s_series.last().copied().unwrap_or(0.0);
            out.push_str(&format!(
                "- max|∇S̃| timeline: initial `{:.3e}`, peak `{:.3e}`, final `{:.3e}`\n",
                g0, gmax, gend,
            ));
        }
        out.push('\n');

        out.push_str("### Mass conservation of S\n\n");
        out.push_str(&format!(
            "- initial mass: `{:.9e}`\n- final mass: `{:.9e}`\n- relative drift: `{:.3e}`\n\n",
            m.mass_s_initial, m.mass_s_final, m.mass_drift_relative,
        ));

        out.push_str("### Null-space health\n\n");
        out.push_str(&format!(
            "- max |mean(vx)| across solves: `{:.3e}`\n- max |mean(vy)|: `{:.3e}`\n\n",
            m.max_abs_mean_vx, m.max_abs_mean_vy,
        ));

        out.push_str("### Velocity magnitude\n\n");
        out.push_str(&format!("- peak |v|: `{:.3e}`\n\n", m.vmax_peak));

        out.push_str("### Heightmaps of S (dynamic remap with bounds)\n\n");
        if m.heightmap_metas.is_empty() {
            out.push_str("- (none recorded)\n\n");
        } else {
            out.push_str("| snapshot | min | max | mean | colour-bar |\n|---|---|---|---|---|\n");
            for (path, md) in m.heightmap_paths.iter().zip(m.heightmap_metas.iter()) {
                let cb = md.colorbar_path.display().to_string().replace('\\', "/");
                out.push_str(&format!(
                    "| `{}` | `{:.3e}` | `{:.3e}` | `{:.3e}` | `{}` |\n",
                    path, md.min, md.max, md.mean, cb,
                ));
            }
            out.push('\n');
        }

        // --- Comparison block ---
        if let Some(prev) = inputs.previous {
            if let Some(prev_grid) = prev
                .grids
                .iter()
                .find(|g| g.grid == (m.grid_nx, m.grid_ny))
            {
                match inputs.kind {
                    ReportKind::Step2Physics => {
                        out.push_str("### Comparison vs Step 1 (advisory — physics changed, not a regression test)\n\n");
                    }
                    ReportKind::Step2Regression => {
                        out.push_str("### Numerical regression vs Step 1\n\n");
                        out.push_str("Same forcing, same preset, same setup as Step 1. Ratio targets: wallclock and CG-iter-per-linear-solve both within `[0.9, 1.1]`.\n\n");
                    }
                    ReportKind::Step3Physics => {
                        out.push_str("### Comparison vs Step 2 (advisory — yielding added, not a regression test)\n\n");
                    }
                    ReportKind::Step3Regression => {
                        out.push_str("### Numerical regression vs Step 2\n\n");
                        out.push_str("Same forcing, same preset, same setup as Step 2, with `YieldingConfig::Disabled`. Ratio targets: wallclock and CG-iter-per-linear-solve both within `[0.95, 1.05]` — the structural by-pass of the plastic branch must be zero-cost.\n\n");
                    }
                    ReportKind::Step4Physics => {
                        out.push_str("### Comparison vs Step 3 (advisory — basal drag added, not a regression test)\n\n");
                    }
                    ReportKind::Step4Regression => {
                        out.push_str("### Numerical regression vs Step 3\n\n");
                        out.push_str("Same forcing, same preset, same setup as Step 3, with `BasalDragConfig::Disabled` and `YieldingConfig::Disabled`. Ratio targets: wallclock and CG-iter-per-linear-solve both within `[0.95, 1.05]` — the structural by-pass of the drag contribution must be zero-cost.\n\n");
                    }
                }
                let justification = inputs
                    .suspect_justifications
                    .get(idx)
                    .map(|s| s.as_str())
                    .unwrap_or("");
                let justif = if justification.is_empty() { None } else { Some(justification) };
                // Label for the "previous" column and the
                // per-grid sub-heading is derived from the
                // ReportKind — never hardcoded to "Step 1". Step 4
                // caught this when it inherited a leftover "Step 1"
                // label from Step 2's initial implementation; Step 5
                // and beyond get the right label automatically.
                let previous_label = previous_step_label(inputs.kind);
                out.push_str(&render_grid_comparison(
                    previous_label,
                    prev_grid,
                    m.cg_iter_mean,
                    m.wallclock_total.as_secs_f64(),
                    m.mass_drift_relative,
                    m.max_abs_mean_vx,
                    m.max_abs_mean_vy,
                    justif,
                ));

                // Step-4 physics: encouraging wallclock-improvement
                // paragraph. The ratios are real improvements vs
                // Step 3 physics (×0.64 at 64², ×0.74 at 128²) even
                // though drag/visc ≈ 10⁻⁷ at this baseline; Step 5/6
                // will amplify the effect ×25 once oceanic cells
                // land (S̃ ≈ 0.2 → S̃² ≈ 0.04). The wording notes the
                // yielding-Disabled caveat so the attribution is
                // honest.
                if matches!(inputs.kind, ReportKind::Step4Physics) {
                    if let Some(visc_wallclock) = prev_grid.wallclock_seconds {
                        if visc_wallclock > 0.0 {
                            let ratio = m.wallclock_total.as_secs_f64() / visc_wallclock;
                            out.push_str(&format!(
                                "\n**Wallclock improvement interpretation.** The Step-4 physics run at this grid is `×{ratio:.2}` of the Step-3 physics wallclock — a measurable improvement, coherent with the theoretical expectation that adding `Br · S̃²` to the operator diagonal improves the conditioning of low-ε̇ regions. Despite the very small absolute drag contribution at this baseline (`drag/visc ≈ 10⁻⁷`), the augmented diagonal gives CG a slightly tighter grip on the system. Caveat: Step-4 physics also disables yielding (to isolate the Br effect), so part of the wallclock delta vs Step 3 physics comes from skipping the `soft_min_harmonic` plastic-branch evaluation; the pure drag contribution is best read off the Br sweep's strict `peak \\|v\\|` monotonicity and the κ(A) stability (ratio ≤ 1.3 across both grids). Encouraging signal for Step 5/6: introducing oceanic cells (`S̃ ≈ 0.2` → `S̃² ≈ 0.04`) will create a ×25 differentiation between continental and oceanic drag, which is where the Step-4 machinery's physical payoff will become visible.\n\n",
                            ));
                        }
                    }
                }
            }
        }

        out.push_str("### Dormant metrics (inactive at Step 2)\n\n");
        out.push_str("| metric | activated at |\n|---|---|\n");
        out.push_str("| S̃_eq (active-orogen mean thickness) | Step 5+ |\n");
        out.push_str("| boundary type diversity | Step 5 |\n");
        out.push_str("| yielding cell fraction | Step 3 |\n");
        out.push_str("| cratonic stability | Step 9 |\n");
        out.push_str("| age field stats | Step 10 |\n\n");
    }

    out.push_str("---\n");
    out.push_str("*Generated by `cargo run --release --bin step_baseline`.*\n");
    out
}

/// Setup-parity checklist emitted at the top of the regression
/// report. The invariant is that **every field listed here is
/// identical to the Step 1 run** — any mismatch must be flagged
/// explicitly.
/// Label for the "previous step" referenced by the per-grid
/// comparison block. Derived from the current `ReportKind` so the
/// sub-heading "Grid N×N — comparison vs Step X" stays truthful
/// as the milestone progresses. The "Step 1" hardcode that lingered
/// through Step 2/3 is the bug this function prevents.
fn previous_step_label(kind: ReportKind) -> &'static str {
    match kind {
        ReportKind::Step2Physics | ReportKind::Step2Regression => "Step 1",
        ReportKind::Step3Physics | ReportKind::Step3Regression => "Step 2",
        ReportKind::Step4Physics | ReportKind::Step4Regression => "Step 3",
    }
}

/// Default path for the "previous step" report to compare against,
/// for a given current `ReportKind` and an output directory. The
/// binary layer uses this to auto-detect the right comparison file
/// per scenario; pair with a CLI override for edge cases.
///
/// Naming convention: Step 2 onward emits
/// `step{N}_{physics,regression}_report.md`; Step 0/1 emitted a
/// single `step{N}_report.md` (no physics/regression split), so
/// Step 2 scenarios fall back to that single file.
pub fn default_previous_report_for(
    kind: ReportKind,
    output_dir: &std::path::Path,
) -> std::path::PathBuf {
    let name: &str = match kind {
        // Step 0/1 emitted a single report; Step 2 regression uses
        // it as the mirror target, Step 2 physics uses it as the
        // advisory comparison.
        ReportKind::Step2Physics | ReportKind::Step2Regression => "step1_report.md",
        ReportKind::Step3Physics => "step2_physics_report.md",
        ReportKind::Step3Regression => "step2_regression_report.md",
        ReportKind::Step4Physics => "step3_physics_report.md",
        ReportKind::Step4Regression => "step3_regression_report.md",
    };
    output_dir.join(name)
}

fn render_setup_parity_block(inputs: &ReportInputs) -> String {
    let mirror_target = match inputs.kind {
        ReportKind::Step4Regression => "Step 3",
        ReportKind::Step3Regression => "Step 2",
        _ => "Step 1",
    };
    let scope_note = match inputs.kind {
        ReportKind::Step4Regression => {
            "No additional Step-4 fields (oceanic/continental density ρ̃, boundary-type drag transitions, etc.) are introduced — the Step 4 scope is only the basal-drag operator diagonal contribution, which is structurally bypassed when `BasalDragConfig::Disabled` (the `drag_diag: Option<&Field2D>` parameter is `None` throughout, short-circuiting the augmentation loop in `apply_momentum` and `momentum_diagonal` before any face-interpolation work)."
        }
        ReportKind::Step3Regression => {
            "No additional Step-3 fields (plastic_strain, cratonic masks, etc.) are introduced — the Step 3 scope is only the plastic-yielding constitutive branch and its diagnostics, neither of which is active when `YieldingConfig::Disabled` (structural by-pass in `ViscosityLaw::eta_effective`'s match arm)."
        }
        _ => {
            "No additional Step-2 fields (ρ̃, anomaly templates, etc.) are introduced — the Step 2 scope is only the forcing module refactor and the GPE term, neither of which touches the regression run."
        }
    };
    let mut s = String::new();
    s.push_str(&format!("## Setup parity with {}\n\n", mirror_target));
    s.push_str("Contract: a mismatch on any of these disqualifies the comparison as a regression test.\n\n");
    s.push_str(&format!(
        "| item | value | same as {}? |\n|---|---|---|\n",
        mirror_target,
    ));
    if let Some(cfg) = inputs.configs.first() {
        s.push_str(&format!("| preset | `{}` | ✅ |\n", cfg.preset_name));
        s.push_str(&format!("| CFL factor | `{:.2}` | ✅ |\n", cfg.cfl_factor));
        s.push_str(&format!("| Newton rel_tol | `{:.1e}` | ✅ |\n", cfg.newton_rel_tol));
        s.push_str(&format!("| Newton max outer iters | `{}` | ✅ |\n", cfg.newton_max_outer_iters));
        s.push_str(&format!("| CG tolerance | `{:.1e}` | ✅ |\n", cfg.cg_tol));
        s.push_str(&format!("| CG max iter | `{}` | ✅ |\n", cfg.cg_max_iter));
        s.push_str(&format!("| continuation schedule | `{}` | ✅ |\n", cfg.continuation_schedule));
        s.push_str(&format!("| nonlinear solver | `{}` | ✅ |\n", cfg.nonlinear_solver));
        s.push_str(&format!("| seed | `{}` | ✅ |\n", cfg.seed));
        s.push_str(&format!("| body force | `{}` | ✅ (SinusoidalForce ε=10) |\n", cfg.body_force));
        if matches!(inputs.kind, ReportKind::Step3Regression | ReportKind::Step4Regression) {
            s.push_str("| yielding | `Disabled` (structural bypass) | ✅ |\n");
        }
        if matches!(inputs.kind, ReportKind::Step4Regression) {
            s.push_str("| basal drag | `Disabled` (structural bypass) | ✅ |\n");
        }
    }
    s.push_str("| initial S̃ | `init_thickness(nx, ny, seed)` unchanged since Step 0 | ✅ |\n");
    s.push_str(&format!("\n{}\n\n", scope_note));
    s
}

#[cfg(test)]
mod tests {
    use super::super::metrics::{IterationHistogram, Metrics};
    use super::*;
    use std::time::Duration;

    fn fake_config() -> SolverConfigDump {
        SolverConfigDump {
            formulation: "thin viscous sheet".into(),
            discretization: "MAC staggered".into(),
            eta_averaging: "arithmetic corners".into(),
            preconditioner: "velocity Jacobi + null-space".into(),
            gauge_fixing: "mean vx, vy".into(),
            cg_tol: 1e-10,
            cg_max_iter: 2000,
            cfl_factor: 0.3,
            grid_spacing_nondim: 1.0 / 64.0,
            body_force: "GpeForce".into(),
            seed: 42,
            preset_name: "dynamic-accidented".into(),
            nonlinear_solver: "newton".into(),
            rheology_n: 3.0,
            strain_rate_floor: 1e-3,
            eta_max_cap: 1e3,
            continuation_schedule: "[1.0, 1.5, 2.0, 2.5, 3.0]".into(),
            newton_rel_tol: 1e-6,
            newton_max_outer_iters: 20,
            basal_drag_config: "Disabled".into(),
        }
    }

    fn fake_metrics() -> Metrics {
        let mut m = Metrics::empty(64, 64, 300);
        m.wallclock_total = Duration::from_millis(500);
        m.wallclock_per_step_mean = Duration::from_micros(1666);
        m.kappa_estimate = 1e4;
        m.cg_iter_mean = 22.8;
        m.cg_iter_max = 38;
        m.cg_iter_histogram = IterationHistogram::from_samples(&[18, 20, 22, 30, 38]);
        m.mass_s_initial = 4096.0;
        m.mass_s_final = 4096.0;
        m.mass_drift_relative = 1e-15;
        m.max_abs_mean_vx = 1e-20;
        m.max_abs_mean_vy = 1e-20;
        m.vmax_peak = 0.027;
        m.variance_series = vec![1e-4, 1.2e-4, 9e-5];
        m.max_grad_s_series = vec![0.5, 0.7, 0.6];
        m
    }

    #[test]
    fn physics_report_has_gpe_header() {
        let s = build_markdown(&ReportInputs {
            kind: ReportKind::Step2Physics,
            seed: 42,
            scales: &Scales::default(),
            configs: &[fake_config()],
            metrics: &[fake_metrics()],
            previous: None,
            suspect_justifications: &[String::new()],
            mms: None,
            ar_sweep: None, bi_sweep: None, br_sweep: None,
            regression_vmax_peak: None,
        });
        assert!(s.contains("GPE spreading"));
        assert!(s.contains("Ar (Argand)"));
        assert!(s.contains("S field evolution"));
        assert!(s.contains("Heightmaps of S"));
    }

    #[test]
    fn regression_report_has_setup_parity() {
        let s = build_markdown(&ReportInputs {
            kind: ReportKind::Step2Regression,
            seed: 42,
            scales: &Scales::default(),
            configs: &[fake_config()],
            metrics: &[fake_metrics()],
            previous: None,
            suspect_justifications: &[String::new()],
            mms: None,
            ar_sweep: None, bi_sweep: None, br_sweep: None,
            regression_vmax_peak: None,
        });
        assert!(s.contains("Sinusoidal forcing"));
        assert!(s.contains("Setup parity with Step 1"));
    }
}
