//! ADR 0001 Finding 102 — does the ×0.4 window TRANSFER, or is it a rule? **No production change.**
//!
//! Every measurement of this chantier has been made on one seed, and Finding 92 asked for a second
//! one without getting it. ⚠️ **The seeds here are not chosen, they are DERIVED**: `splitmix64`
//! applied once and twice to `PSEED`, printed below, so there is no way to have picked them for
//! their looks.
//!
//! ⚠️ **Which experiment this is, stated before it runs.** `k_s` is already per-seed by
//! construction — Finding 97 defined it as the median of `S·√A` over the DELIVERED channel cells.
//! So "apply ×0.4 to a new seed" has two readings: **(a)** the absolute floor `0.0451 × 0.4`, and
//! **(b)** the rule `k_s(seed) × 0.4`. **(b) is the closure as built and is what is swept**; each
//! seed's `k_s` is printed so (a) can be read off it.
//!
//! ⚠️ `D_L` is **NEW to the dossier**: case-sensitive with word boundaries it has **0 hits**; the 9
//! case-insensitive hits are all collisions with `channel_head_law` and `detected_lake`, both of
//! which contain the substring `d_l`. It is calibrated on a disk and a 10:1 rectangle here BEFORE
//! any field is read, because a rasterised perimeter counted on cell edges cannot return the
//! textbook 1.00 for a disk.
//!
//! Run: cargo test -p ymir-core --release --test f102_seeds -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Crit, Knobs, PSEED, SEA, aniso, build_field_seed, f95_criteria, land_u16, majority,
    pct, sorted, to_mask,
};
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const DELIVERED_P50_M: f32 = 424.2;
/// Finding 99: on seed 1 the oracle sits **+15.1 %** above the delivered field. ⚠️ PROXY — the
/// oracle does not exist for a new seed (1 h 43 each, refused), so the relief target is read as
/// this ratio against each seed's OWN delivered field and is labelled a proxy everywhere.
const ORACLE_OVER_DELIVERED: f32 = 1.151;

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn spurs_at(mask: &[bool], w: usize, min_km: f32) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, min_km, NECK_KM).0.len()
}

/// `wc == 2` components identified by BOTTOM CELL — the integrity column Finding 101 established
/// after `spillways.len()` turned out to count outflows rather than basins.
fn components(f: &GridF32, w: usize, h: usize) -> Vec<usize> {
    let n = w * h;
    let wc = water_class(f, SEA);
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for s in 0..n {
        if wc[s] != 2 || seen[s] {
            continue;
        }
        let (mut stack, mut cells) = (vec![s], Vec::new());
        seen[s] = true;
        while let Some(k) = stack.pop() {
            cells.push(k);
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
        out.push(
            *cells.iter().min_by(|&&a, &&b| f.data[a].partial_cmp(&f.data[b]).unwrap()).unwrap(),
        );
    }
    out
}

/// `D_L = P / (2·√(π·A))`, the shoreline development index, on a rasterised body: `P` counted as
/// 4-connected boundary edges × cell length. **1.00 for a perfect disk in the continuum** — and NOT
/// 1.00 on a raster, which is why this is calibrated before it is used.
fn d_l(cells: &[usize], inside: &HashSet<usize>, w: usize, h: usize) -> f32 {
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
    p / (2.0 * (std::f32::consts::PI * a).sqrt())
}

/// The median `D_L` over bodies ≥ 1 km², and the share above 5 (drowned valley networks).
fn lake_shapes(lake_map: &[u32], w: usize, h: usize, cell_km2: f32) -> (f32, f32, f32, usize) {
    let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..w * h {
        if lake_map[k] != 0 {
            by.entry(lake_map[k]).or_default().push(k);
        }
    }
    let mut v: Vec<f32> = Vec::new();
    for (_, cells) in by.iter() {
        if cells.len() as f32 * cell_km2 < 1.0 {
            continue;
        }
        let inside: HashSet<usize> = cells.iter().copied().collect();
        v.push(d_l(cells, &inside, w, h));
    }
    let n = v.len();
    let hi = v.iter().filter(|&&x| x > 5.0).count() as f32;
    let s = sorted(v);
    (pct(&s, 0.50), pct(&s, 0.90), 100.0 * hi / n.max(1) as f32, n)
}

/// The median channel slope and the Flint intercept `k_s = median(S·√A)` over `A ≥ A_c`.
fn channel_stats(
    f: &GridF32,
    on: &C1DrainageConfig,
    ss: &SteinSteinParams,
    cell_km2: f32,
) -> (f32, f32) {
    let (w, h) = (f.width, f.height);
    let d = c1_drainage_windowed(f, None, on, ss, DOMAIN_KM);
    let n2m = c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);
    let (mut s_all, mut ks) = (Vec::new(), Vec::new());
    for k in (0..w * h).step_by(7) {
        if f.data[k] <= SEA {
            continue;
        }
        let a = d.flow.accumulation.data[k] * cell_km2;
        if a < RELIEF_V1_A_C_KM2 {
            continue;
        }
        let dir = d.flow.direction[k];
        if dir == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[dir as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[dir as usize]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[dir as usize] != 0 && D8_DY[dir as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let s = (f.data[k] - f.data[ny * w + nx]).max(0.0) * n2m / dx_m;
        if s > 0.0 {
            s_all.push(s);
            ks.push(s * a.sqrt());
        }
    }
    (pct(&sorted(s_all), 0.50), pct(&sorted(ks), 0.50))
}

#[test]
#[ignore]
fn f102_seeds() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let s2 = splitmix64(PSEED);
    let s3 = splitmix64(s2);
    eprintln!("\n==========  Finding 102 . does the window transfer?  ==========");
    eprintln!(
        "   seeds, DERIVED not chosen (splitmix64 of PSEED, once and twice):\n     seed 1 \
         **{PSEED}**\n     seed 2 **{s2}**\n     seed 3 **{s3}**"
    );

    // ── E first: the shape instrument, calibrated before anything is read ───
    for (label, mk) in [("disk r=64", 0u8), ("rect 10:1 (40x400)", 1u8)] {
        let (mut cells, mut inside) = (Vec::new(), HashSet::new());
        let (w, h) = (1024usize, 1024usize);
        for y in 0..h {
            for x in 0..w {
                let ok = if mk == 0 {
                    let (dx, dy) = (x as f32 - 512.0, y as f32 - 512.0);
                    dx * dx + dy * dy <= 64.0 * 64.0
                } else {
                    (300..700).contains(&y) && (492..532).contains(&x)
                };
                if ok {
                    cells.push(y * w + x);
                    inside.insert(y * w + x);
                }
            }
        }
        eprintln!(
            "   E calibration . {label:<20} {} cells => **D_L = {:.3}** (continuum: disk 1.000, \
             10:1 rectangle 1.963)",
            cells.len(),
            d_l(&cells, &inside, w, h)
        );
    }

    // ── A/B: per seed ───────────────────────────────────────────────────────
    let t0 = Instant::now();
    let mut summary: Vec<(u64, f32, f32, f32, Vec<(f32, f32, i64, f32, f32, f32, usize)>)> =
        Vec::new();
    for (si, seed) in [PSEED, s2, s3].into_iter().enumerate() {
        eprintln!("\n========  SEED {} ({seed})  ========", si + 1);
        let pre = build_field_seed(Knobs::no_incision(), seed);
        let (w, h) = (pre.width, pre.height);
        let n = w * h;
        let delivered = build_field_seed(Knobs::passes(2), seed);
        let auth = spurs_at(&land_u16(&pre, &ss), w, 1.0);
        let intensity = |f: &GridF32| -> f64 {
            let c: Vec<usize> = (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] > SEA).collect();
            c.iter().map(|&k| ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64).sum::<f64>()
                / c.len().max(1) as f64
        };
        let i_del = intensity(&delivered);
        let (s_pre, ks_pre) = channel_stats(&pre, &on, &ss, cell_km2);
        let (s_del, ks_del) = channel_stats(&delivered, &on, &ss, cell_km2);
        let ddel = c1_drainage_windowed(&delivered, None, &on, &ss, DOMAIN_KM);
        let bdel = breach_monotone(&delivered, &ddel.flow.filled, &ddel.lake_map, SEA, w, h);
        let cdel =
            f95_criteria(&delivered, &bdel, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        let ldel: Vec<bool> = (0..n).map(|k| delivered.data[k] > SEA).collect();
        let (dl50, dl90, dlhi, dln) = lake_shapes(&ddel.lake_map, w, h, cell_km2);
        let target = cdel.p50 * ORACLE_OVER_DELIVERED;
        eprintln!(
            "   CONTROL delivered: relief **{:.1} m** . coast >=1km **{:+}** . >=2cells **{}** . \
             canyons **{}/{} ({:.1} %)** . lakes {:.2} % . R8 **{:.4}** . components **{}** . \
             catchment {:.0} . D_L p50 **{:.2}** p90 {:.2} (>5: {:.0} %, n={dln}) . cut/cell \
             {i_del:.1} m",
            cdel.p50,
            spurs_at(&land_u16(&bdel, &ss), w, 1.0) as i64 - auth as i64,
            spurs_at(&land_u16(&bdel, &ss), w, 2.0 * CELL_KM),
            cdel.klass,
            cdel.scanned,
            100.0 * cdel.klass as f32 / cdel.scanned.max(1) as f32,
            cdel.lake_pct,
            aniso(&delivered, &ldel, 16).r8,
            components(&delivered, w, h).len(),
            cdel.max_catchment_km2,
            dl50,
            dl90,
            dlhi
        );
        eprintln!(
            "   channels: median slope pre **{s_pre:.4}** delivered **{s_del:.4}** . k_s pre \
             {ks_pre:.4} **delivered {ks_del:.4}** (seed 1 at Finding 97: 0.0451) . relief PROXY \
             target = {:.1} m (delivered x{ORACLE_OVER_DELIVERED})",
            target
        );

        let base = components(&delivered, w, h).len() as f32;
        let facs: Vec<f32> = if si == 0 { vec![0.4] } else { vec![0.3, 0.4, 0.5, 0.6] };
        let mut rows = Vec::new();
        for m in facs {
            let t = Instant::now();
            let f = build_field_seed(
                Knobs {
                    slope_floor_uk: Some(ks_del * m),
                    depression_floor: true,
                    ..Knobs::passes(2)
                },
                seed,
            );
            let cost = t.elapsed().as_secs_f64();
            let d = c1_drainage_windowed(&f, None, &on, &ss, DOMAIN_KM);
            let bf = breach_monotone(&f, &d.flow.filled, &d.lake_map, SEA, w, h);
            let c = f95_criteria(&f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
            let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
            let r8 = aniso(&f, &land, 16).r8;
            let inten = intensity(&f);
            let sp = spurs_at(&land_u16(&bf, &ss), w, 1.0) as i64 - auth as i64;
            let sp2 = spurs_at(&land_u16(&bf, &ss), w, 2.0 * CELL_KM);
            let comp = components(&f, w, h).len();
            let (dl, _, dlh, _) = lake_shapes(&d.lake_map, w, h, cell_km2);
            let cut: Vec<f32> = (0..n)
                .filter(|&k| pre.data[k] > SEA && f.data[k] > SEA)
                .map(|k| ((pre.data[k] - f.data[k]) * n2m).max(0.0))
                .collect();
            let med = pct(&sorted(cut), 0.50);
            let rel_pct = 100.0 * (c.p50 - target) / target;
            eprintln!(
                "   **x{m}** (k_s {:.4}): relief **{:.1} m ({:+.1} % vs PROXY target)** . coast \
                 **{sp:+}** . >=2cells **{sp2}** . canyons **{:.1} %** . lakes {:.2} % . erosion \
                 **{:.1} %** . R8 **{r8:.4}** . components **{comp}** ({:+.0} %) . catchment \
                 {:.0} . D_L {dl:.2} (>5: {dlh:.0} %) . median cut {med:.1} m (mean/median \
                 {:.1}) . {cost:.1} s",
                ks_del * m,
                c.p50,
                rel_pct,
                100.0 * c.klass as f32 / c.scanned.max(1) as f32,
                c.lake_pct,
                100.0 * inten / i_del,
                100.0 * (comp as f32 / base - 1.0),
                c.max_catchment_km2,
                inten as f32 / med.max(0.1)
            );
            rows.push((m, rel_pct, sp, 100.0 * inten as f32 / i_del as f32, r8, dl, comp));
        }
        summary.push((seed, s_del, ks_del, target, rows));
    }

    // ── B: the window per seed ──────────────────────────────────────────────
    eprintln!(
        "\n-- B . the window per seed: relief within +-5 % of the PROXY target, coast <= +5, rule 14 >= 70 %, R8 < 0.06 --"
    );
    for (seed, s_del, ks_del, _t, rows) in &summary {
        let win: Vec<f32> = rows
            .iter()
            .filter(|(_, rel, sp, ero, r8, _, _)| {
                rel.abs() <= 5.0 && *sp <= 5 && *ero >= 70.0 && *r8 < 0.06
            })
            .map(|r| r.0)
            .collect();
        eprintln!(
            "   seed {seed}: median channel slope {s_del:.4} . k_s {ks_del:.4} . **window = {:?}**",
            win
        );
    }
    eprintln!("\n==========  end Finding 102 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
