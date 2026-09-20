//! ADR 0001 Finding 106 — is there an ABSOLUTE `S_eq = c·A^(−1/2)` that passes the author's gates
//! on three seeds, with no calibration pass? **No promotion in this bench.**
//!
//! ⚠️ **The point-3 control is ill-posed as the round writes it and is run in a repaired form.**
//! `c = 0.6 × 0.045095 = 0.027057` cannot be bit-identical to the ×0.6 row, because **0.045095 is
//! a six-digit PRINTING** of a full-precision `f32` and the bench built its floor from the
//! unrounded value. Both floors are printed with their bit patterns and the field difference is
//! measured: a micron-scale discrepancy is a decimal rounding, not a different closure.
//!
//! ⚠️ **The delivered rows are CITED, not re-measured** (Findings 102/104): canyons 29.6 / 29.4 /
//! 33.3 %, components 20 / 43 / 15, largest catchment 9 779 / 735 / 5 215 km², D_L share above 5
//! 32 / 42 / 24 %. Re-running their criteria chains would cost 9 minutes to reproduce numbers the
//! dossier already carries at the digit.
//!
//! Run: cargo test -p ymir-core --release --test f106_constant -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, f95_criteria, land_u16, majority, pct,
    sorted, to_mask,
};
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, flint_intercept};
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const DELIVERED_P50_M: f32 = 424.2;
const CS: [f32; 4] = [0.018, 0.021, 0.024, 0.027];
/// Findings 102 / 104, cited rather than re-measured: canyon rate, components, largest catchment,
/// D_L share above 5, and the delivered relief, per seed.
const DEL: [(f32, usize, f32, f32, f32); 3] = [
    (29.6, 20, 9779.0, 32.0, 424.2),
    (29.4, 43, 735.0, 42.0, 246.9),
    (33.3, 15, 5215.0, 24.0, 363.0),
];

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn hash(f: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for v in &f.data {
        for b in v.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
    }
    h
}

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

fn components(f: &GridF32, w: usize, h: usize) -> usize {
    let n = w * h;
    let wc = water_class(f, SEA);
    let mut seen = vec![false; n];
    let mut c = 0usize;
    for s in 0..n {
        if wc[s] != 2 || seen[s] {
            continue;
        }
        c += 1;
        let mut st = vec![s];
        seen[s] = true;
        while let Some(k) = st.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for d in 0..8 {
                let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if wc[nk] == 2 && !seen[nk] {
                    seen[nk] = true;
                    st.push(nk);
                }
            }
        }
    }
    c
}

/// Share of lakes ≥ 1 km² with `D_L > 5`. ⚠️ Finding 102 calibrated the raster bias: a disk reads
/// 1.284, not 1.000. Comparable between fields measured identically, never with textbook values.
fn dl_share(lake_map: &[u32], w: usize, h: usize, cell_km2: f32) -> f32 {
    let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..w * h {
        if lake_map[k] != 0 {
            by.entry(lake_map[k]).or_default().push(k);
        }
    }
    let (mut n, mut hi) = (0usize, 0usize);
    for cells in by.values() {
        if cells.len() as f32 * cell_km2 < 1.0 {
            continue;
        }
        let inside: HashSet<usize> = cells.iter().copied().collect();
        let mut e = 0usize;
        for &k in cells {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for (dx, dy) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    e += 1;
                } else if !inside.contains(&(ny as usize * w + nx as usize)) {
                    e += 1;
                }
            }
        }
        let p = e as f32 * CELL_M;
        let a = cells.len() as f32 * CELL_M * CELL_M;
        n += 1;
        if p / (2.0 * (std::f32::consts::PI * a).sqrt()) > 5.0 {
            hi += 1;
        }
    }
    100.0 * hi as f32 / n.max(1) as f32
}

#[test]
#[ignore]
fn f106_constant() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let s2 = splitmix64(PSEED);
    let s3 = splitmix64(s2);
    eprintln!("\n==========  Finding 106 . the absolute constant  ==========");
    let t0 = Instant::now();

    for (si, seed) in [PSEED, s2, s3].into_iter().enumerate() {
        let (dc, dcomp, dcatch, ddl, drelief) = DEL[si];
        let pre = build_field_seed(Knobs::no_incision(), seed);
        let (w, h) = (pre.width, pre.height);
        let n = w * h;
        let delivered = build_field_seed(Knobs::passes(2), seed);
        let auth = spurs(&land_u16(&pre, &ss), w);
        let ks = flint_intercept(&delivered, &ss, SEA, cell_km2, DOMAIN_KM);
        let intensity = |f: &GridF32| -> f64 {
            let c: Vec<usize> = (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] > SEA).collect();
            c.iter().map(|&k| ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64).sum::<f64>()
                / c.len().max(1) as f64
        };
        let i_del = intensity(&delivered);
        eprintln!(
            "\n========  SEED {} ({seed})  ========\n   k_s = **{ks:.8}** (bits {:08x}) . delivered \
             CITED from Findings 102/104: canyons {dc} %, components {dcomp}, catchment {dcatch}, \
             D_L>5 {ddl} %, relief {drelief} m . coast authority {auth}",
            si + 1,
            ks.to_bits()
        );

        // ── the point-3 control, repaired, on seed 1 only ───────────────────
        if si == 0 {
            let prod = ks * 0.6f32;
            let dec = 0.027057f32;
            eprintln!(
                "   ⚠️ CONTROL: k_s×0.6 = **{prod:.9}** (bits {:08x}) . the round's decimal \
                 0.027057 = **{dec:.9}** (bits {:08x}) ⇒ **{}**",
                prod.to_bits(),
                dec.to_bits(),
                if prod.to_bits() == dec.to_bits() {
                    "SAME f32 -- the control is well posed after all"
                } else {
                    "**DIFFERENT f32** -- the decimal is a 6-digit rounding, as predicted"
                }
            );
            let a = build_field_seed(
                Knobs { slope_floor_uk: Some(prod), depression_floor: true, ..Knobs::passes(2) },
                seed,
            );
            let b = build_field_seed(
                Knobs { slope_floor_uk: Some(dec), depression_floor: true, ..Knobs::passes(2) },
                seed,
            );
            let d: Vec<f32> = (0..n)
                .filter(|&k| a.data[k] != b.data[k])
                .map(|k| ((a.data[k] - b.data[k]) * n2m).abs())
                .collect();
            let s = sorted(d.clone());
            eprintln!(
                "   ⚠️ CONTROL fields: hashes {} . cells differing **{}** . |delta| p50 **{:.6} m** \
                 max **{:.4} m** ⇒ the control passes **{}**",
                if hash(&a) == hash(&b) { "**MATCH**" } else { "differ" },
                d.len(),
                if s.is_empty() { 0.0 } else { pct(&s, 0.50) },
                s.last().copied().unwrap_or(0.0),
                if s.last().copied().unwrap_or(0.0) < 0.01 {
                    "IN SUBSTANCE (sub-centimetre)"
                } else {
                    "**NOT AT ALL -- stop**"
                }
            );
        }

        // ── A: the four constants ───────────────────────────────────────────
        for c in CS {
            let t = Instant::now();
            let f = build_field_seed(
                Knobs { slope_floor_uk: Some(c), depression_floor: true, ..Knobs::passes(2) },
                seed,
            );
            let cost = t.elapsed().as_secs_f64();
            let dr = c1_drainage_windowed(&f, None, &on, &ss, DOMAIN_KM);
            let bf = breach_monotone(&f, &dr.flow.filled, &dr.lake_map, SEA, w, h);
            let cr = f95_criteria(&f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
            let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
            let rate = 100.0 * cr.klass as f32 / cr.scanned.max(1) as f32;
            let comp = components(&f, w, h);
            let dl = dl_share(&dr.lake_map, w, h, cell_km2);
            let sp = spurs(&land_u16(&bf, &ss), w) as i64 - auth as i64;
            let ero = 100.0 * intensity(&f) / i_del;
            let eff = c / ks;
            let pass = sp <= 5
                && rate <= 0.25 * dc + 0.001
                && comp == dcomp
                && cr.max_catchment_km2 >= dcatch * 0.98
                && dl <= 0.1;
            eprintln!(
                "   **c = {c}** (effective ×{eff:.3}{}): coast **{sp:+}** . canyons **{:.1} %** \
                 (gate ≤ {:.1}) . components **{comp}** . catchment **{:.0}** ({:+.1} %) . D_L>5 \
                 **{dl:.0} %** . relief {:.1} m ({:+.1} % of delivered) . R8 {:.4} . erosion \
                 **{ero:.1} %**{} . {cost:.0} s ⇒ **{}**",
                if !(0.3..=0.7).contains(&eff) { ", ⚠️ OUTSIDE the swept range" } else { "" },
                rate,
                0.25 * dc,
                cr.max_catchment_km2,
                100.0 * (cr.max_catchment_km2 / dcatch - 1.0),
                cr.p50,
                100.0 * (cr.p50 / drelief - 1.0),
                aniso(&f, &land, 16).r8,
                if ero < 43.8 { " ⛔ UNDER the B1 alert" } else { "" },
                if pass { "PASS" } else { "**FAIL**" }
            );
        }
    }
    eprintln!("\n==========  end Finding 106 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
