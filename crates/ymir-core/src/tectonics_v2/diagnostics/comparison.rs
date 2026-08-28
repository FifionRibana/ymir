//! Parse a previous step's markdown report and emit a comparison
//! section against the current step's metrics.
//!
//! Not a general-purpose markdown parser — we only extract the
//! handful of figures we want to diff. The parser scans the file for
//! grid-qualified blocks ("`## Grid 64×64`") and then for known
//! leading substrings within those blocks.
//!
//! The classification of CG-iter ratio follows the Step-1 spec
//! palier rules: `< 1.3` idéal, `< 3` acceptable, `> 3` suspect,
//! `> 10` fail. Only the category is reported numerically; text in
//! `suspect` tier gets a note flagging it for investigation.

use std::fs;
use std::path::Path;

#[derive(Clone, Debug, Default)]
pub struct StepReference {
    /// One entry per grid size appearing in the previous report, in
    /// declaration order.
    pub grids: Vec<GridReference>,
}

#[derive(Clone, Debug, Default)]
pub struct GridReference {
    pub grid: (usize, usize),
    pub wallclock_seconds: Option<f64>,
    /// CG-iter-per-*linear-solve* mean. At Step 0 this is the mean
    /// CG iterations over the 300 sheet solves.
    pub cg_iters_mean: Option<f64>,
    pub mass_drift_relative: Option<f64>,
    pub max_abs_mean_vx: Option<f64>,
    pub max_abs_mean_vy: Option<f64>,
}

/// Parse a Step-0 (or later) report and return per-grid references.
pub fn parse_step_report(path: &Path) -> Result<StepReference, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("read {:?}: {}", path, e))?;
    let mut out = StepReference::default();
    let mut current: Option<GridReference> = None;

    for line in text.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("## Grid ") {
            if let Some(g) = current.take() {
                out.grids.push(g);
            }
            if let Some((a, b)) = rest.split_once('×') {
                let nx: usize = a.trim().parse().unwrap_or(0);
                let ny: usize = b.trim().parse().unwrap_or(0);
                current = Some(GridReference { grid: (nx, ny), ..Default::default() });
            }
        } else if let Some(rest) = line.strip_prefix("- wallclock total: `") {
            let val = extract_backtick_f64(rest, " s");
            if let Some(g) = current.as_mut() {
                g.wallclock_seconds = val;
            }
        } else if let Some(rest) = line.strip_prefix("- CG iterations per sheet solve — mean: `")
        {
            // Step 0 format.
            let val = extract_backtick_f64_leading(rest);
            if let Some(g) = current.as_mut() {
                g.cg_iters_mean = val;
            }
        } else if let Some(rest) = line.strip_prefix("- CG iterations per Newton step — mean: `")
        {
            // Step 1+ format (Newton wraps every linear solve).
            let val = extract_backtick_f64_leading(rest);
            if let Some(g) = current.as_mut() {
                g.cg_iters_mean = val;
            }
        } else if let Some(rest) = line.strip_prefix("- relative drift: `") {
            let val = extract_backtick_f64(rest, "`");
            if let Some(g) = current.as_mut() {
                g.mass_drift_relative = val;
            }
        } else if let Some(rest) = line.strip_prefix("- max |mean(vx)| across solves: `") {
            let val = extract_backtick_f64(rest, "`");
            if let Some(g) = current.as_mut() {
                g.max_abs_mean_vx = val;
            }
        } else if let Some(rest) = line.strip_prefix("- max |mean(vy)|: `") {
            let val = extract_backtick_f64(rest, "`");
            if let Some(g) = current.as_mut() {
                g.max_abs_mean_vy = val;
            }
        }
    }
    if let Some(g) = current.take() {
        out.grids.push(g);
    }
    Ok(out)
}

fn extract_backtick_f64(s: &str, terminator: &str) -> Option<f64> {
    let end = s.find(terminator)?;
    let raw = &s[..end];
    raw.trim().parse::<f64>().ok()
}

/// Extract the first backtick-enclosed number in `s`.
fn extract_backtick_f64_leading(s: &str) -> Option<f64> {
    let end = s.find('`')?;
    s[..end].trim().parse::<f64>().ok()
}

/// CG-iter ratio classifier per the Step 1 spec.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CgRatioTier {
    Ideal,      // < 1.3
    Acceptable, // < 3
    Suspect,    // 3..=10
    Fail,       // > 10
    Unknown,    // missing reference
}

impl CgRatioTier {
    pub fn classify(ratio: f64) -> Self {
        if !ratio.is_finite() {
            return CgRatioTier::Unknown;
        }
        if ratio < 1.3 {
            CgRatioTier::Ideal
        } else if ratio < 3.0 {
            CgRatioTier::Acceptable
        } else if ratio <= 10.0 {
            CgRatioTier::Suspect
        } else {
            CgRatioTier::Fail
        }
    }
    pub fn label(&self) -> &'static str {
        match self {
            CgRatioTier::Ideal => "idéal",
            CgRatioTier::Acceptable => "acceptable",
            CgRatioTier::Suspect => "suspect",
            CgRatioTier::Fail => "fail",
            CgRatioTier::Unknown => "N/A",
        }
    }
}

/// Comparison section snippet for one grid, suitable for embedding in
/// the step-N markdown report.
pub fn render_grid_comparison(
    prev_label: &str,
    prev: &GridReference,
    cur_cg_iters_mean: f64,
    cur_wallclock_seconds: f64,
    cur_mass_drift: f64,
    cur_mean_vx: f64,
    cur_mean_vy: f64,
    justification_if_suspect: Option<&str>,
) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "#### Grid {}×{} — comparison vs {}\n\n",
        prev.grid.0, prev.grid.1, prev_label,
    ));
    out.push_str("| metric | previous | current | ratio / note |\n");
    out.push_str("|---|---|---|---|\n");

    if let Some(w0) = prev.wallclock_seconds {
        let ratio = cur_wallclock_seconds / w0.max(1e-9);
        let note = if ratio > 20.0 { " (>20×, flag)" } else { "" };
        out.push_str(&format!(
            "| wallclock (s) | {:.3} | {:.3} | ×{:.2}{} |\n",
            w0, cur_wallclock_seconds, ratio, note,
        ));
    }
    if let Some(c0) = prev.cg_iters_mean {
        // If the previous baseline had a solver-trivial scenario
        // (~0 CG iterations), a naive ratio blows up and is not
        // physically meaningful. Declare the comparison not-
        // applicable in that case and report the current absolute
        // count instead.
        if c0 < 0.5 {
            out.push_str(&format!(
                "| CG iters / linear solve (mean) | {:.1} (solver-trivial) | {:.1} | N/A — no denominator; report absolute |\n",
                c0, cur_cg_iters_mean,
            ));
        } else {
            let ratio = cur_cg_iters_mean / c0;
            let tier = CgRatioTier::classify(ratio);
            let note = match tier {
                CgRatioTier::Suspect => match justification_if_suspect {
                    Some(j) => format!(" — suspect; justification: {}", j),
                    None => " — suspect (no justification on file)".into(),
                },
                _ => String::new(),
            };
            out.push_str(&format!(
                "| CG iters / linear solve (mean) | {:.1} | {:.1} | ×{:.2} [{}]{} |\n",
                c0,
                cur_cg_iters_mean,
                ratio,
                tier.label(),
                note,
            ));
        }
    }
    if let Some(m0) = prev.mass_drift_relative {
        out.push_str(&format!(
            "| S mass drift (relative) | {:.3e} | {:.3e} | gate 1e-10 |\n",
            m0, cur_mass_drift,
        ));
    }
    if let Some(v0) = prev.max_abs_mean_vx {
        out.push_str(&format!(
            "| max \\|mean(vx)\\| | {:.3e} | {:.3e} | bruit machine |\n",
            v0, cur_mean_vx,
        ));
    }
    if let Some(v0) = prev.max_abs_mean_vy {
        out.push_str(&format!(
            "| max \\|mean(vy)\\| | {:.3e} | {:.3e} | bruit machine |\n",
            v0, cur_mean_vy,
        ));
    }
    out.push('\n');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cg_ratio_tier_boundaries() {
        assert_eq!(CgRatioTier::classify(1.0), CgRatioTier::Ideal);
        assert_eq!(CgRatioTier::classify(1.29), CgRatioTier::Ideal);
        assert_eq!(CgRatioTier::classify(1.5), CgRatioTier::Acceptable);
        assert_eq!(CgRatioTier::classify(2.99), CgRatioTier::Acceptable);
        assert_eq!(CgRatioTier::classify(5.0), CgRatioTier::Suspect);
        assert_eq!(CgRatioTier::classify(10.0), CgRatioTier::Suspect);
        assert_eq!(CgRatioTier::classify(15.0), CgRatioTier::Fail);
        assert_eq!(CgRatioTier::classify(f64::NAN), CgRatioTier::Unknown);
    }

    #[test]
    fn parse_finds_grid_and_wallclock() {
        let fake = "# Report\n\n## Grid 64×64\n\n### Timing\n\n- wallclock total: `0.152 s`\n- wallclock per step (mean): `0.505 ms`\n- steps: `300`\n\n";
        let path = std::env::temp_dir().join("v2_comparison_test.md");
        std::fs::write(&path, fake).unwrap();
        let r = parse_step_report(&path).unwrap();
        assert_eq!(r.grids.len(), 1);
        assert_eq!(r.grids[0].grid, (64, 64));
        assert!((r.grids[0].wallclock_seconds.unwrap() - 0.152).abs() < 1e-9);
    }
}
