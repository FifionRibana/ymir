//! ADR 0001 Finding 116-B/D — **separating the two operators Finding 41 measured as one.**
//!
//! ⛔ **Rule 11 found block B's table already written, at Finding 41 (L1802)**, with more depth
//! than the round asks for. *"A closed depression = an 8-connected component the priority-flood
//! raises by > 0.1 m."*
//!
//! | stage | pits | ≤ 2-cell | median depth |
//! |---|---|---|---|
//! | coarse post-isostasy | **16** | 11 | 115 m |
//! | **after FBM upscale** | **90 682** | 70 335 | 3.2 m |
//! | after relief-v3 incision (2 iter) | 75 060 | 51 968 | 2.1 m |
//!
//! …plus the maturity curve **78 464 / 75 060 / 64 936 / 53 379 / 47 777** at 1/2/4/8/16 passes,
//! *"decelerating to a plateau near ~45 k — still ~3000× the tectonic 16"*, and the verdict
//! *"maturity alone does NOT deliver France"*. Every erosion sub-process **reduces** the count;
//! none creates any.
//!
//! ## What is genuinely new here
//!
//! ⛔ **`16 → 90 682` crosses TWO operators in one step**: the bicubic interpolation from the
//! coarse grid to 8192², **and** the FBM. Finding 41 does not separate them, so "the FBM is the
//! SOLE creator" is a claim about a composite. **`amplitude_base = 0` separates them**, and that
//! is this bench's first column.
//!
//! ⚠️ Thresholds: Finding 41's **0.1 m** is carried so the numbers are comparable, and the round's
//! 5 m and 20 m are added. **Independent `if` flags, never `else if`** — the defect I reported
//! against my own Finding 115 classification, which the round adopted.
//!
//! Run: cargo test -p ymir-core --release --test f116_fbm -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, land_u16, majority, pct, sorted, to_mask,
};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const S_EQ: f32 = 0.024;
/// Finding 41's own threshold, carried so this table can be read against its table.
const T41: f32 = 0.1;

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

/// 8-connected components of `filled − raw > t`, with their areas and depths. Finding 41's
/// definition, at three thresholds.
fn pits(
    fill: &GridF32,
    raw: &GridF32,
    n2m: f32,
    t: f32,
    cell_km2: f32,
) -> (usize, usize, f32, f32) {
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let deep: Vec<bool> = (0..n).map(|k| (fill.data[k] - raw.data[k]) * n2m > t).collect();
    let mut seen = vec![false; n];
    let (mut count, mut tiny) = (0usize, 0usize);
    let mut depths = Vec::new();
    let mut areas = Vec::new();
    for s in 0..n {
        if !deep[s] || seen[s] {
            continue;
        }
        let mut q = vec![s];
        seen[s] = true;
        let mut cells = 0usize;
        let mut dmax = 0f32;
        while let Some(k) = q.pop() {
            cells += 1;
            dmax = dmax.max((fill.data[k] - raw.data[k]) * n2m);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nx = (x + dx).rem_euclid(w as i32) as usize;
                    let ny = (y + dy).rem_euclid(h as i32) as usize;
                    let nk = ny * w + nx;
                    if deep[nk] && !seen[nk] {
                        seen[nk] = true;
                        q.push(nk);
                    }
                }
            }
        }
        count += 1;
        if cells <= 2 {
            tiny += 1;
        }
        depths.push(dmax);
        areas.push(cells as f32 * cell_km2);
    }
    let sd = sorted(depths);
    let sa = sorted(areas);
    (
        count,
        tiny,
        if sd.is_empty() { 0.0 } else { pct(&sd, 0.50) },
        if sa.is_empty() { 0.0 } else { sa.iter().sum::<f32>() },
    )
}

#[test]
#[ignore]
fn f116_fbm() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!(
        "\n==========  Finding 116-B/D . the FBM sweep, and the operator Finding 41 hid  ==="
    );
    eprintln!(
        "   Finding 41, for comparison (its own > {T41} m threshold): coarse post-isostasy **16** \
         pits · after FBM upscale **90 682** · after incision **75 060**"
    );

    // ── B1 · the PRE-INCISION field, FBM amplitude swept, incision OFF ──────
    eprintln!(
        "\n   ── B1 · `Knobs::no_incision()` (tectonics + bicubic + FBM, NO stream power) ──\n   \
         {:>6} {:>10} {:>9} {:>9} {:>10} {:>9} {:>9} {:>10} {:>8}",
        "amp", ">0.1m n", "≤2cell", "p50 m", ">5m n", "p50 m", ">20m n", "area km²", "spurs"
    );
    for amp in [0.08f64, 0.04, 0.02, 0.01, 0.0] {
        let f = build_field_seed(Knobs { fbm_amp: Some(amp), ..Knobs::no_incision() }, PSEED);
        let dr = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
        let (c1, t1, d1, _) = pits(&dr.flow.filled, &f, n2m, T41, cell_km2);
        let (c5, _, d5, a5) = pits(&dr.flow.filled, &f, n2m, 5.0, cell_km2);
        let (c20, _, _, _) = pits(&dr.flow.filled, &f, n2m, 20.0, cell_km2);
        eprintln!(
            "   {amp:>6.2} {c1:>10} {t1:>9} {d1:>9.2} {c5:>10} {d5:>9.2} {c20:>9} {a5:>10.1} \
             {:>8}",
            spurs(&land_u16(&f, &ss), f.width)
        );
    }

    // ── B2 / D · the DELIVERED field at each amplitude ──────────────────────
    eprintln!(
        "\n   ── B2/D · the DELIVERED field (closure ON, s_eq {S_EQ}) at each amplitude ──\n   \
         {:>6} {:>10} {:>9} {:>9} {:>10} {:>9} {:>9} {:>8}",
        "amp", ">0.1m n", ">5m n", ">20m n", "relief p50", "R8", "σ local", "spurs"
    );
    for amp in [0.08f64, 0.04, 0.02, 0.01, 0.0] {
        let f = build_field_seed(
            Knobs { fbm_amp: Some(amp), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
            PSEED,
        );
        let (w, h) = (f.width, f.height);
        let dr = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
        let (c1, _, _, _) = pits(&dr.flow.filled, &f, n2m, T41, cell_km2);
        let (c5, _, _, _) = pits(&dr.flow.filled, &f, n2m, 5.0, cell_km2);
        let (c20, _, _, _) = pits(&dr.flow.filled, &f, n2m, 20.0, cell_km2);
        let land: Vec<bool> = f.data.iter().map(|&v| v > SEA).collect();
        let relief: Vec<f32> = (0..w * h)
            .filter(|&k| land[k])
            .map(|k| c1_altitude_norm_to_metres(f.data[k], &ss))
            .collect();
        // local σ over 3×3, land only — the texture column the author cares about
        let mut sig = Vec::new();
        for y in 1..h - 1 {
            for x in (1..w - 1).step_by(7) {
                let k = y * w + x;
                if !land[k] {
                    continue;
                }
                let (mut s, mut s2) = (0f64, 0f64);
                for dy in 0..3 {
                    for dx in 0..3 {
                        let v = (f.data[(y + dy - 1) * w + x + dx - 1] * n2m) as f64;
                        s += v;
                        s2 += v * v;
                    }
                }
                sig.push(((s2 / 9.0 - (s / 9.0) * (s / 9.0)).max(0.0)).sqrt() as f32);
            }
        }
        eprintln!(
            "   {amp:>6.2} {c1:>10} {c5:>9} {c20:>9} {:>10.1} {:>9.4} {:>9.3} {:>8}",
            pct(&sorted(relief), 0.50),
            aniso(&f, &land, 16).r8,
            pct(&sorted(sig), 0.50),
            spurs(&land_u16(&f, &ss), w)
        );
    }
    eprintln!(
        "\n==========  end Finding 116-B/D . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
