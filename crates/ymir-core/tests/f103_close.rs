//! ADR 0001 Finding 103 — the two builds that close the re-read. **No production change.**
//!
//! Finding 103's block B re-reads Findings 100 and 102 under declared gates and finds the
//! intersection empty, with **one gate standing between ×0.6 and a common factor**
//! (`erosion ≥ 70 %`). Two columns are missing from seed 1 at ×0.5 and ×0.6 — the `wc == 2`
//! COMPONENT count and `D_L`, because Finding 100 pre-dates both. These are the two builds the
//! round allows, spent where they decide.
//!
//! ⚠️ Everything else in block B needed no build: Finding 102 swept seeds 2 and 3 at 0.3–0.6 and
//! Finding 100 swept seed 1 at 0.3–0.7, so the round's assumed gaps (seed 2 at ×0.6, seed 3 at
//! ×0.5) already exist.
//!
//! Run: cargo test -p ymir-core --release --test f103_close -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, f95_criteria, pct, sorted};
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const DELIVERED_P50_M: f32 = 424.2;
/// Seed 1's own Flint intercept (Finding 97).
const KS: f32 = 0.0451;

fn components(f: &GridF32, w: usize, h: usize) -> usize {
    let n = w * h;
    let wc = water_class(f, SEA);
    let mut seen = vec![false; n];
    let mut count = 0usize;
    for s in 0..n {
        if wc[s] != 2 || seen[s] {
            continue;
        }
        count += 1;
        let mut stack = vec![s];
        seen[s] = true;
        while let Some(k) = stack.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for d in 0..8 {
                let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if wc[nk] == 2 && !seen[nk] {
                    seen[nk] = true;
                    stack.push(nk);
                }
            }
        }
    }
    count
}

/// `D_L = P / (2√(πA))` on a raster. ⚠️ Finding 102 calibrated this: a disk reads **1.284**, not
/// 1.000 (the 4/π staircase bias of a perimeter counted on cell edges), while an axis-aligned
/// rectangle is exact. Comparable between fields measured identically, never with textbook values.
fn lake_shapes(lake_map: &[u32], w: usize, h: usize, cell_km2: f32) -> (f32, f32, usize) {
    let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..w * h {
        if lake_map[k] != 0 {
            by.entry(lake_map[k]).or_default().push(k);
        }
    }
    let mut v: Vec<f32> = Vec::new();
    for cells in by.values() {
        if cells.len() as f32 * cell_km2 < 1.0 {
            continue;
        }
        let inside: HashSet<usize> = cells.iter().copied().collect();
        let mut edges = 0usize;
        for &k in cells {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for (dx, dy) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    edges += 1;
                } else if !inside.contains(&(ny as usize * w + nx as usize)) {
                    edges += 1;
                }
            }
        }
        let p = edges as f32 * CELL_M;
        let a = cells.len() as f32 * CELL_M * CELL_M;
        v.push(p / (2.0 * (std::f32::consts::PI * a).sqrt()));
    }
    let n = v.len();
    let hi = v.iter().filter(|&&x| x > 5.0).count() as f32;
    (pct(&sorted(v), 0.50), 100.0 * hi / n.max(1) as f32, n)
}

#[test]
#[ignore]
fn f103_close() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 103 . the two builds that close the re-read  ==========");
    eprintln!(
        "   seed 1 only. Reference from Findings 100/102: delivered relief 424.2 m, canyons \
         29.6 %, R8 0.0924, components **20**, catchment 9779, D_L p50 4.49 (>5: **32 %**)."
    );

    let t0 = Instant::now();
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let delivered = build_field_seed(Knobs::passes(2), PSEED);
    let intensity = |f: &GridF32| -> f64 {
        let c: Vec<usize> = (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] > SEA).collect();
        c.iter().map(|&k| ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64).sum::<f64>()
            / c.len().max(1) as f64
    };
    let i_del = intensity(&delivered);
    eprintln!("   prerequisites in {:.1} s", t0.elapsed().as_secs_f64());

    for m in [0.5f32, 0.6] {
        let t = Instant::now();
        let f = build_field_seed(
            Knobs { slope_floor_uk: Some(KS * m), depression_floor: true, ..Knobs::passes(2) },
            PSEED,
        );
        let cost = t.elapsed().as_secs_f64();
        let d = c1_drainage_windowed(&f, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(&f, &d.flow.filled, &d.lake_map, SEA, w, h);
        let c = f95_criteria(&f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
        let (dl, dlhi, dln) = lake_shapes(&d.lake_map, w, h, cell_km2);
        eprintln!(
            "   **seed 1 x{m}**: relief **{:.1} m ({:+.1} % of delivered)** . canyons **{:.1} %** \
             . erosion **{:.1} %** . R8 **{:.4}** ({:+.1} %) . **components {}** . catchment \
             **{:.0}** ({:+.1} %) . **D_L p50 {dl:.2}, share > 5 = {dlhi:.0} %** (n={dln}) . \
             {cost:.1} s",
            c.p50,
            100.0 * (c.p50 / 424.2 - 1.0),
            100.0 * c.klass as f32 / c.scanned.max(1) as f32,
            100.0 * intensity(&f) / i_del,
            aniso(&f, &land, 16).r8,
            100.0 * (aniso(&f, &land, 16).r8 / 0.0924 - 1.0),
            components(&f, w, h),
            c.max_catchment_km2,
            100.0 * (c.max_catchment_km2 / 9779.0 - 1.0)
        );
    }
    eprintln!("\n==========  end . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
