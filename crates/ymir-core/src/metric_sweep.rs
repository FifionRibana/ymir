//! ADR method rule 10 — **a metric block must prove it measured something.**
//!
//! ## The defect this exists to make impossible
//!
//! Twice now a bench has reported a column of zeros and it was read as a clean result. The
//! first H-1 bench measured lakes on the POST-breach drainage — production discards those, and
//! a breached field has no depressions by construction — so it printed `lakes 0` in every
//! configuration and the lake invariants trivially "passed": zero violations out of zero
//! lakes. The channel-head invariant suite (ADR Finding 56) reproduced the identical mistake,
//! on a suite whose entire purpose was to catch what the coastal metrics could not see.
//!
//! The rule was already written down after the first occurrence. Writing it down did not stop
//! the second occurrence. **A recorded rule does not protect if it is not executable** — so
//! rule 10 is this type rather than a paragraph.
//!
//! ## Scope
//!
//! **Any bench reporting more than one configuration.** With a single configuration there is
//! nothing to vary against and the guard cannot apply; with two or more, every reported
//! quantity must differ across at least one pair, and a quantity that is constant at zero
//! everywhere is a measurement failure, not a result.
//!
//! ## What it asserts, and what it deliberately does not
//!
//! - [`Sweep::push`] — **constant at zero across the whole sweep → FAIL.** Nothing was
//!   measured. This is the error above, and it is unconditional: there is no legitimate
//!   reading of a column that is zero in every configuration of a sweep designed to move it.
//!   **Constant at a non-zero value → FAIL too**, because it is either the same defect or a
//!   rule-7 invariant that was not declared as one.
//! - [`Sweep::invariant`] — a rule-7 control, where flatness is the *result* (land area under
//!   a change that must not move it). The guard is inverted: it fails if the column MOVED.
//! - [`Sweep::checked`] — **an invariant check, as a violation count WITH its population.**
//!   This is the pair that catches the real shape of the defect: a violation counter flat at
//!   zero is a good result, so the thing that must be guarded is the DENOMINATOR. `0
//!   violations out of 0 lakes` reads exactly like `0 out of 4 000` and is not a pass; an
//!   empty population is an unrun check.
//! - [`Sweep::pinned`] — the narrow escape hatch: the bench states in writing that a quantity
//!   is not expected to respond (a thresholded classification whose thresholds the model
//!   cannot reach yet). The reason is mandatory and reprinted, and **a pinned column at zero
//!   still fails** — otherwise the hatch would swallow the very case the rule exists for.
//! - A quantity printed in only SOME configurations is reported: it cannot be compared, and
//!   that is how "Strahler S5 59 → 0" hid as an absent column rather than a collapse.
//! - It does **not** check magnitudes, signs or plausibility. It answers one question only:
//!   *did this column respond to the sweep at all?*
//!
//! ```no_run
//! use ymir_core::metric_sweep::Sweep;
//! let mut sw = Sweep::new("channel-head law");
//! for cfg in ["off", "on"] {
//!     sw.config(cfg);
//!     sw.push("lakes", 12.0);
//!     sw.checked("footprint above level", 0, 12); // 0 violations over 12 lakes
//!     sw.invariant("land km²", 23_916.0); // must NOT move — rule 7
//! }
//! sw.verify(); // panics naming every column that measured nothing
//! ```

use std::collections::BTreeMap;

/// What a column promises, and therefore what failing looks like for it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Guard {
    /// [`Sweep::push`] — the sweep is expected to MOVE this. Constant across every
    /// configuration is a failure; constant at zero is the flagship failure.
    MustVary,
    /// [`Sweep::invariant`] — a rule-7 control. Flatness is the RESULT, so the guard is
    /// inverted: this fails if it moved.
    MustNotVary,
    /// The population half of [`Sweep::checked`] — must be non-empty in every configuration.
    /// It may legitimately be constant (the same lakes under both settings).
    MustBeNonZero,
    /// The violation half of [`Sweep::checked`] — zero everywhere is the DESIRED result and
    /// carries no promise of its own. Its meaning comes entirely from its population column.
    Free,
    /// [`Sweep::pinned`] — the bench states IN WRITING that this column is not expected to
    /// respond, and why. It may be flat or it may move; what it may NOT be is zero
    /// everywhere. That case stays unwaivable, because a column pinned at zero measured
    /// nothing whatever the excuse.
    Pinned,
}

/// One reported quantity, accumulated across the sweep's configurations.
#[derive(Clone, Debug)]
struct Column {
    /// `(configuration label, value)`, in report order.
    values: Vec<(String, f64)>,
    guard: Guard,
    /// [`Guard::Pinned`] only: the written reason this column is not expected to respond.
    why: String,
}

/// A sweep's metric ledger. Record every quantity you print, then [`verify`](Sweep::verify).
#[derive(Clone, Debug)]
pub struct Sweep {
    what: String,
    current: String,
    /// `BTreeMap` so the failure report is in a deterministic, readable order.
    cols: BTreeMap<String, Column>,
    order: Vec<String>,
}

impl Sweep {
    /// A ledger for the sweep of `what` (used in the panic message).
    pub fn new(what: impl Into<String>) -> Self {
        Sweep {
            what: what.into(),
            current: String::from("<unlabelled>"),
            cols: BTreeMap::new(),
            order: Vec::new(),
        }
    }

    /// Open a configuration. Every [`push`](Sweep::push) until the next `config` is attributed
    /// to `label`.
    pub fn config(&mut self, label: impl Into<String>) {
        self.current = label.into();
        if !self.order.contains(&self.current) {
            self.order.push(self.current.clone());
        }
    }

    /// Record a reported quantity that the sweep is expected to MOVE.
    pub fn push(&mut self, name: &str, value: impl Into<f64>) {
        self.record(name, value.into(), Guard::MustVary, "");
    }

    /// Record a quantity the change must NOT move — a rule-7 control. The guard is inverted:
    /// this column fails if it varies, not if it is flat.
    pub fn invariant(&mut self, name: &str, value: impl Into<f64>) {
        self.record(name, value.into(), Guard::MustNotVary, "");
    }

    /// Record an INVARIANT CHECK: `violations` found over a population of `population`.
    ///
    /// **This is the pair that catches the actual defect.** A violation counter constant at
    /// zero is a perfectly good result — an invariant that holds — so guarding it for
    /// variation would be wrong. What was wrong in both failures was the DENOMINATOR:
    /// `0 violations out of 0 lakes` passes vacuously and reads exactly like `0 out of
    /// 4 000`. So the violations are recorded [`Guard::Free`] and the population is required
    /// NON-EMPTY in every configuration. An empty population is never a clean invariant; it
    /// is an unrun check.
    pub fn checked(&mut self, name: &str, violations: u64, population: u64) {
        self.record(name, violations as f64, Guard::Free, "");
        self.record(&format!("{name} [population]"), population as f64, Guard::MustBeNonZero, "");
    }

    /// Record a quantity the bench does NOT claim will respond, stating `why` in writing.
    ///
    /// The escape hatch, deliberately narrow. It exists so that the honest response to a
    /// genuinely pinned quantity is to DOCUMENT it rather than to delete the metric — which is
    /// what a guard with no escape actually produces. Two things keep it from being a hole:
    /// the reason is mandatory and is reprinted in the report, and **a column pinned at zero
    /// still fails**. The motivating case is the navigability class count at 8192²: the
    /// thresholds are in m³/s and the discharge sits orders of magnitude below them (Finding
    /// 47), so the classification cannot move until the hypsometry converges — that is a
    /// result about the model, not a measurement that failed.
    pub fn pinned(&mut self, name: &str, value: impl Into<f64>, why: &str) {
        debug_assert!(!why.is_empty(), "`pinned` requires a written reason");
        self.record(name, value.into(), Guard::Pinned, why);
    }

    fn record(&mut self, name: &str, value: f64, guard: Guard, why: &str) {
        let cfg = self.current.clone();
        let c = self.cols.entry(name.to_string()).or_insert_with(|| Column {
            values: Vec::new(),
            guard,
            why: why.to_string(),
        });
        debug_assert_eq!(c.guard, guard, "column `{name}` registered under two different guards");
        c.values.push((cfg, value));
    }

    /// The columns that measured nothing, as `(name, reason)`. Empty means the block is sound.
    ///
    /// A column seen in only ONE configuration is reported too: rule 10's scope is a
    /// multi-configuration bench, so a quantity printed for one setting and not the others
    /// cannot be compared and is a hole in the report.
    pub fn failures(&self) -> Vec<(String, String)> {
        let n_cfg = self.order.len();
        if n_cfg < 2 {
            return Vec::new(); // out of scope: nothing to vary against
        }
        let mut out = Vec::new();
        for (name, col) in &self.cols {
            if col.values.len() < n_cfg {
                out.push((
                    name.clone(),
                    format!(
                        "reported in {} of {n_cfg} configurations — not comparable",
                        col.values.len()
                    ),
                ));
                continue;
            }
            let first = col.values[0].1;
            // Exact inequality is the right test: these are measured aggregates, and a change
            // too small to alter an f64 did not alter the reported number either.
            let varies = col.values.iter().any(|(_, v)| *v != first);
            let spread = || -> String {
                col.values.iter().map(|(c, v)| format!("{c}={v}")).collect::<Vec<_>>().join(", ")
            };
            match col.guard {
                Guard::Free => {}
                Guard::MustBeNonZero => {
                    let empty: Vec<&str> = col
                        .values
                        .iter()
                        .filter(|(_, v)| *v == 0.0)
                        .map(|(c, _)| c.as_str())
                        .collect();
                    if !empty.is_empty() {
                        out.push((
                            name.clone(),
                            format!(
                                "EMPTY POPULATION in [{}] — the invariant it carries was never \
                                 actually checked there (0 violations out of 0 is not a pass)",
                                empty.join(", ")
                            ),
                        ));
                    }
                }
                Guard::MustNotVary => {
                    if varies {
                        out.push((
                            name.clone(),
                            format!("declared INVARIANT but moved: {}", spread()),
                        ));
                    } else if first == 0.0 {
                        // An invariant flat at ZERO is indistinguishable from an unwired
                        // measurement, so it does not get the invariant exemption either.
                        out.push((
                            name.clone(),
                            String::from(
                                "declared INVARIANT and constant AT ZERO — an invariant that \
                                 holds at nothing is not evidence that it holds",
                            ),
                        ));
                    }
                }
                Guard::Pinned => {
                    if !varies && first == 0.0 {
                        out.push((
                            name.clone(),
                            format!(
                                "PINNED AT ZERO — the reason given (\"{}\") does not cover a \
                                 column that is zero everywhere; that is still an unwired \
                                 measurement",
                                col.why
                            ),
                        ));
                    }
                }
                Guard::MustVary => {
                    if !varies {
                        out.push((
                            name.clone(),
                            if first == 0.0 {
                                String::from(
                                    "CONSTANT AT ZERO in every configuration — this column \
                                     measured NOTHING (rule 10's flagship case: the \
                                     post-breach lake population)",
                                )
                            } else {
                                format!(
                                    "constant at {first} in every configuration — either it is \
                                     a rule-7 invariant (register it with `Sweep::invariant`) \
                                     or the measurement is not wired to the change"
                                )
                            },
                        ));
                    }
                }
            }
        }
        out
    }

    /// Panic naming every column that measured nothing. Call at the END of a metric block —
    /// after the numbers are printed, so the report is on screen when the guard fires.
    ///
    /// # Panics
    /// If any reported quantity is constant across the sweep (or any declared invariant moved).
    pub fn verify(&self) {
        let f = self.failures();
        if f.is_empty() {
            eprintln!(
                "[rule 10] {}: {} column(s) over {} configuration(s), all responded",
                self.what,
                self.cols.len(),
                self.order.len()
            );
            for (name, col) in self.cols.iter().filter(|(_, c)| c.guard == Guard::Pinned) {
                eprintln!("[rule 10]   PINNED  {name}: {}", col.why);
            }
            return;
        }
        let mut msg = format!(
            "\n[rule 10] {} — {} REPORTED QUANTITY/QUANTITIES MEASURED NOTHING:\n",
            self.what,
            f.len()
        );
        for (name, why) in &f {
            msg.push_str(&format!("  · {name}: {why}\n"));
        }
        msg.push_str(
            "A column that does not respond to the sweep is not a result. Fix the \
             measurement before reading the numbers above.\n",
        );
        panic!("{msg}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The exact defect: a column of zeros across the sweep must FAIL, not pass as clean.
    #[test]
    fn a_column_constant_at_zero_fails() {
        let mut sw = Sweep::new("t");
        sw.config("off");
        sw.push("lakes", 0.0);
        sw.config("on");
        sw.push("lakes", 0.0);
        let f = sw.failures();
        assert_eq!(f.len(), 1);
        assert!(f[0].1.contains("CONSTANT AT ZERO"), "{}", f[0].1);
    }

    /// NEGATIVE CONTROL (rule 1): the same shape, one value moved, must pass — otherwise the
    /// guard above would be firing on something other than the constancy.
    #[test]
    fn a_column_that_varies_passes() {
        let mut sw = Sweep::new("t");
        sw.config("off");
        sw.push("lakes", 0.0);
        sw.config("on");
        sw.push("lakes", 3.0);
        assert!(sw.failures().is_empty());
        sw.verify();
    }

    /// The defect in its ACTUAL shape: zero violations over an EMPTY population. The
    /// violation column is legitimately flat at zero, so only the population guard can catch
    /// it — and it must, in both configurations.
    #[test]
    fn zero_violations_over_an_empty_population_fails() {
        let mut sw = Sweep::new("t");
        for cfg in ["off", "on"] {
            sw.config(cfg);
            sw.checked("footprint above level", 0, 0);
        }
        let f = sw.failures();
        assert_eq!(f.len(), 1, "{f:?}");
        assert!(f[0].0.contains("[population]"), "{}", f[0].0);
        assert!(f[0].1.contains("EMPTY POPULATION"), "{}", f[0].1);
    }

    /// NEGATIVE CONTROL (rule 1): the SAME zero violation count, now over a real population,
    /// must pass — so the failure above is the empty population and nothing else. Note the
    /// population is CONSTANT here (4 000 both ways) and that is fine: a population needs to
    /// be non-empty, not to move.
    #[test]
    fn zero_violations_over_a_real_population_passes() {
        let mut sw = Sweep::new("t");
        for cfg in ["off", "on"] {
            sw.config(cfg);
            sw.checked("footprint above level", 0, 4_000);
        }
        assert!(sw.failures().is_empty(), "{:?}", sw.failures());
    }

    /// An empty population in ONE configuration only is still a failure — that
    /// configuration's check did not run.
    #[test]
    fn an_empty_population_in_one_configuration_fails() {
        let mut sw = Sweep::new("t");
        sw.config("off");
        sw.checked("orphan mouths", 0, 4_000);
        sw.config("on");
        sw.checked("orphan mouths", 0, 0);
        let f = sw.failures();
        assert_eq!(f.len(), 1, "{f:?}");
        assert!(f[0].1.contains("[on]"), "{}", f[0].1);
    }

    #[test]
    fn constant_non_zero_fails_unless_declared_invariant() {
        let mut sw = Sweep::new("t");
        sw.config("off");
        sw.push("n", 7.0);
        sw.config("on");
        sw.push("n", 7.0);
        assert_eq!(sw.failures().len(), 1);

        let mut ok = Sweep::new("t");
        ok.config("off");
        ok.invariant("n", 7.0);
        ok.config("on");
        ok.invariant("n", 7.0);
        assert!(ok.failures().is_empty());
    }

    #[test]
    fn a_declared_invariant_that_moves_fails() {
        let mut sw = Sweep::new("t");
        sw.config("off");
        sw.invariant("land", 100.0);
        sw.config("on");
        sw.invariant("land", 101.0);
        let f = sw.failures();
        assert_eq!(f.len(), 1);
        assert!(f[0].1.contains("declared INVARIANT but moved"), "{}", f[0].1);
    }

    #[test]
    fn a_column_missing_from_a_configuration_is_reported() {
        let mut sw = Sweep::new("t");
        sw.config("off");
        sw.push("s5", 117.0);
        sw.config("on");
        // `on` lost the S5 order entirely — printing nothing here is exactly how the
        // two-sided law's "Strahler S5 59 → 0" hid as an ABSENT column.
        let f = sw.failures();
        assert_eq!(f.len(), 1);
        assert!(f[0].1.contains("not comparable"), "{}", f[0].1);
    }

    /// The escape hatch does NOT cover a column of zeros — that case stays unwaivable.
    #[test]
    fn a_pinned_column_at_zero_still_fails() {
        let mut sw = Sweep::new("t");
        for cfg in ["off", "on"] {
            sw.config(cfg);
            sw.pinned("barge reaches", 0.0, "the discharge is orders below the threshold");
        }
        let f = sw.failures();
        assert_eq!(f.len(), 1, "{f:?}");
        assert!(f[0].1.contains("PINNED AT ZERO"), "{}", f[0].1);
    }

    /// NEGATIVE CONTROL (rule 1): the same pinned column at a NON-zero constant passes, so the
    /// failure above is the zero and not the pinning.
    #[test]
    fn a_pinned_column_at_a_non_zero_constant_passes() {
        let mut sw = Sweep::new("t");
        for cfg in ["off", "on"] {
            sw.config(cfg);
            sw.pinned("boat reaches", 2.0, "thresholds in m3/s, discharge orders below");
        }
        assert!(sw.failures().is_empty(), "{:?}", sw.failures());
    }

    #[test]
    fn an_invariant_flat_at_zero_fails() {
        let mut sw = Sweep::new("t");
        for cfg in ["off", "on"] {
            sw.config(cfg);
            sw.invariant("below-sea basins", 0.0);
        }
        let f = sw.failures();
        assert_eq!(f.len(), 1, "{f:?}");
        assert!(f[0].1.contains("constant AT ZERO"), "{}", f[0].1);
    }

    /// Out of scope with one configuration: nothing to vary against, so no verdict.
    #[test]
    fn a_single_configuration_is_out_of_scope() {
        let mut sw = Sweep::new("t");
        sw.config("only");
        sw.push("lakes", 0.0);
        assert!(sw.failures().is_empty());
    }
}
