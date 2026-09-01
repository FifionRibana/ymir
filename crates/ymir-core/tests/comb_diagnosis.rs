//! COMB artefact diagnostic (read-only, the export IS the verdict). The river overlay
//! shows bundles of strictly parallel, regularly-spaced, rigorously AXIAL secondary
//! tributaries — a rake grafted on valley flanks. It predates C-3b. This reads the
//! SHIPPED export (`rivers.json` + `height.u16` + `flow_accumulation.f32`) and answers,
//! with measurement not argument:
//!   - THE FIRST DISCRIMINANT: is the comb in the TERRAIN or only in the POLYLINES?
//!     (cross-section under comb bundles, ⊥ the segments — grooved → incision artefact;
//!     smooth → routing/extraction artefact);
//!   - is it a GRID signature or a FLOW one? (axial concentration in the grid frame vs
//!     the local-gradient frame);
//!   - are the comb segments secondary? (Strahler + drainage-area distribution);
//!   - the flat-ground hypothesis (slope under comb cells vs normal channel cells).
//!
//! Run: cargo test -p ymir-core --test comb_diagnosis -- --ignored --nocapture

use std::path::Path;

use serde::Deserialize;

const DIR: &str = "../../exports/seed10481999410520546993_8192.ymir";
const N: usize = 8192;
const HMIN: f32 = -5535.70458984375;
const HMAX: f32 = 4145.34228515625;
const KM_PER_CELL: f32 = 400.0 / N as f32;

#[derive(Deserialize)]
struct Rivers {
    segments: Vec<Seg>,
}
#[derive(Deserialize)]
struct Seg {
    points: Vec<[i64; 2]>,
    strahler_order: u32,
    drainage_km2: f32,
}

fn load_height() -> Vec<f32> {
    let bytes = std::fs::read(Path::new(DIR).join("height.u16")).unwrap();
    let span = HMAX - HMIN;
    bytes
        .chunks_exact(2)
        .map(|b| HMIN + (u16::from_le_bytes([b[0], b[1]]) as f32 / 65535.0) * span)
        .collect()
}

#[inline]
fn at(h: &[f32], x: i64, y: i64) -> f32 {
    let xi = x.rem_euclid(N as i64) as usize;
    let yi = y.rem_euclid(N as i64) as usize;
    h[yi * N + xi]
}

/// Slope magnitude (m/m) at a cell via central differences (physical).
fn slope_at(h: &[f32], x: i64, y: i64) -> f32 {
    let cell_m = KM_PER_CELL * 1000.0;
    let gx = (at(h, x + 1, y) - at(h, x - 1, y)) / (2.0 * cell_m);
    let gy = (at(h, x, y + 1) - at(h, x, y - 1)) / (2.0 * cell_m);
    (gx * gx + gy * gy).sqrt()
}

/// Axial concentration R = |mean e^{i2θ}| of a set of directions (1 = one axis, 0 =
/// uniform), plus a 12-bin axial histogram (bin = ⌊(θ mod π)/(π/12)⌋).
fn axial_concentration(dirs: &[(f32, f32)]) -> (f32, [u32; 12]) {
    let (mut sc, mut ss) = (0.0f64, 0.0f64);
    let mut hist = [0u32; 12];
    for &(dx, dy) in dirs {
        let m = (dx * dx + dy * dy).sqrt();
        if m < 1e-9 {
            continue;
        }
        let th = dy.atan2(dx);
        sc += (2.0 * th as f64).cos();
        ss += (2.0 * th as f64).sin();
        let a = (th.rem_euclid(std::f32::consts::PI)) / (std::f32::consts::PI / 12.0);
        hist[(a as usize).min(11)] += 1;
    }
    let n = dirs.len().max(1) as f64;
    (((sc / n).powi(2) + (ss / n).powi(2)).sqrt() as f32, hist)
}

/// Per-segment axiality: the axial concentration of its consecutive step directions
/// (a straight line → 1, a meander → lower).
fn seg_axiality(seg: &Seg) -> f32 {
    let dirs: Vec<(f32, f32)> = seg
        .points
        .windows(2)
        .map(|w| ((w[1][0] - w[0][0]) as f32, (w[1][1] - w[0][1]) as f32))
        .collect();
    if dirs.len() < 3 {
        return 1.0;
    }
    axial_concentration(&dirs).0
}

#[test]
#[ignore]
fn comb_diagnosis() {
    if !Path::new(DIR).exists() {
        eprintln!("export missing — skip");
        return;
    }
    let rivers: Rivers =
        serde_json::from_slice(&std::fs::read(Path::new(DIR).join("rivers.json")).unwrap())
            .unwrap();
    let h = load_height();
    let segs = &rivers.segments;
    eprintln!("\n=== COMB DIAGNOSIS — export {DIR} ===");
    eprintln!("segments: {}", segs.len());

    // Strahler distribution.
    let mut ord = std::collections::BTreeMap::new();
    for s in segs {
        *ord.entry(s.strahler_order).or_insert(0u32) += 1;
    }
    eprintln!("Strahler orders: {ord:?}");

    // COMB candidates: order 1, short/low-area, HIGH axiality (straight rake teeth).
    // Report the axiality distribution first so the threshold is grounded.
    let mut ax: Vec<f32> = segs.iter().map(seg_axiality).collect();
    ax.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |p: f32| ax[((ax.len() as f32 * p) as usize).min(ax.len() - 1)];
    eprintln!(
        "segment axiality (R) p10/p50/p90/p99: {:.2}/{:.2}/{:.2}/{:.2}  (1.0 = perfectly straight)",
        q(0.1),
        q(0.5),
        q(0.9),
        q(0.99)
    );

    let comb: Vec<&Seg> = segs
        .iter()
        .filter(|s| s.strahler_order == 1 && seg_axiality(s) > 0.9 && s.points.len() >= 5)
        .collect();
    let normal: Vec<&Seg> = segs.iter().filter(|s| s.strahler_order >= 3).collect();
    eprintln!(
        "comb candidates (order 1, axiality>0.9, ≥5 pts): {} ({:.1}% of segments)",
        comb.len(),
        100.0 * comb.len() as f32 / segs.len() as f32
    );

    // ── DISCRIMINANT: cross-section ⊥ the segment, under comb vs trunk. Sample the mean
    //    elevation profile relative to the channel cell across ±12 cells. Groove → both
    //    sides rise; flat → ~0.
    let cross = |set: &[&Seg]| -> [f32; 25] {
        let mut acc = [0.0f64; 25];
        let mut cnt = 0u64;
        for s in set {
            // midpoint + local direction.
            let m = s.points.len() / 2;
            if m == 0 || m + 1 >= s.points.len() {
                continue;
            }
            let (px, py) = (s.points[m][0], s.points[m][1]);
            let dx = (s.points[m + 1][0] - s.points[m - 1][0]) as f32;
            let dy = (s.points[m + 1][1] - s.points[m - 1][1]) as f32;
            let dm = (dx * dx + dy * dy).sqrt().max(1e-6);
            let (perpx, perpy) = (-dy / dm, dx / dm); // ⊥ unit
            let h0 = at(&h, px, py);
            for (bi, t) in (-12i64..=12).enumerate() {
                let sx = (px as f32 + perpx * t as f32).round() as i64;
                let sy = (py as f32 + perpy * t as f32).round() as i64;
                acc[bi] += (at(&h, sx, sy) - h0) as f64;
            }
            cnt += 1;
        }
        let mut out = [0.0f32; 25];
        for i in 0..25 {
            out[i] = (acc[i] / cnt.max(1) as f64) as f32 * 1.0; // already metres
        }
        out
    };
    let ccomb = cross(&comb);
    let ctrunk = cross(&normal);
    eprintln!(
        "\ncross-section ⊥ segment, mean Δelevation vs channel cell (m), offsets −12..+12 cells (×{:.0} m):",
        KM_PER_CELL * 1000.0
    );
    let fmt = |c: &[f32; 25]| {
        [c[0], c[4], c[8], c[10], c[11], c[12], c[13], c[14], c[16], c[20], c[24]]
            .iter()
            .map(|v| format!("{v:>5.0}"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    eprintln!("  offset(cells):   -12   -8   -4   -2   -1    0   +1   +2   +4   +8  +12");
    eprintln!("  COMB  (o1):    {}", fmt(&ccomb));
    eprintln!("  TRUNK (o≥3):   {}", fmt(&ctrunk));
    let groove = |c: &[f32; 25]| ((c[8] + c[16]) * 0.5 - c[12]).max((c[4] + c[20]) * 0.5 - c[12]);
    eprintln!(
        "  groove depth (mean flank − channel): COMB {:.0} m | TRUNK {:.0} m",
        groove(&ccomb),
        groove(&ctrunk)
    );

    // ── SLOPE under comb cells vs normal channel cells (flat-ground hypothesis).
    let slope_stats = |set: &[&Seg]| -> (f32, f32) {
        let mut v = Vec::new();
        for s in set {
            for p in &s.points {
                v.push(slope_at(&h, p[0], p[1]).atan().to_degrees());
            }
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = v.get(v.len() / 2).copied().unwrap_or(0.0);
        let flat = v.iter().filter(|&&s| s < 1.0).count() as f32 / v.len().max(1) as f32;
        (med, 100.0 * flat)
    };
    let (cm, cf) = slope_stats(&comb);
    let (nm, nf) = slope_stats(&normal);
    eprintln!(
        "\nslope under cells: COMB median {cm:.2}° ({cf:.0}% < 1°) | TRUNK median {nm:.2}° ({nf:.0}% < 1°)"
    );

    // ── GRID frame vs GRADIENT frame: step directions of comb segments, raw (grid) and
    //    relative to the local downslope (−∇h). Grid-concentrated but gradient-spread =
    //    grid signature (steps set by tie-break, not slope).
    let mut grid_dirs = Vec::new();
    let mut grad_dirs = Vec::new();
    for s in &comb {
        for w in s.points.windows(2) {
            let (dx, dy) = ((w[1][0] - w[0][0]) as f32, (w[1][1] - w[0][1]) as f32);
            grid_dirs.push((dx, dy));
            // local downslope at the from-cell.
            let (x, y) = (w[0][0], w[0][1]);
            let cell_m = KM_PER_CELL * 1000.0;
            let gx = (at(&h, x + 1, y) - at(&h, x - 1, y)) / (2.0 * cell_m);
            let gy = (at(&h, x, y + 1) - at(&h, x, y - 1)) / (2.0 * cell_m);
            let gm = (gx * gx + gy * gy).sqrt();
            if gm < 1e-9 {
                continue;
            }
            // rotate step into the frame where downslope (−∇h) is the x-axis.
            let (ddx, ddy) = (-gx / gm, -gy / gm);
            let rx = dx * ddx + dy * ddy;
            let ry = -dx * ddy + dy * ddx;
            grad_dirs.push((rx, ry));
        }
    }
    let (rg, hg) = axial_concentration(&grid_dirs);
    let (rgr, hgr) = axial_concentration(&grad_dirs);
    eprintln!("\naxial concentration R of comb step directions:");
    eprintln!("  GRID frame     R={rg:.2}  hist(12 bins over 0..π): {hg:?}");
    eprintln!("  GRADIENT frame R={rgr:.2}  hist: {hgr:?}");
    eprintln!(
        "(GRID R high + GRADIENT R low → the steps align on grid axes, not on the slope → grid/tie-break\n signature. GROOVE depth COMB≈0 with TRUNK≫0 → the comb is in the POLYLINES, not the terrain.)"
    );
}
