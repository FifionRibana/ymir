//! ADR 0001 Finding 104 — channel-cut retention, with the negative control the previous gate
//! never had. **No production change.**
//!
//! ⚠️ **The round asked for "no construction, on the existing fields". No field in this campaign
//! is persisted** — `build_field` re-runs the tectonics, the upscale and the incision every time,
//! and `cache.rs` covers the HD pipeline, not this path. So every row here is rebuilt, and the
//! departure is declared rather than hidden.
//!
//! ⚠️ **The population is fixed on the DELIVERED field's network.** Classifying each candidate by
//! its own flow would make the channel sets differ between fields, so the ratio would mix a
//! population change with a cut change — the unpaired-comparison trap of Findings 63, 64, 81, 88.
//!
//! The gate under test replaces `erosion ≥ 70 %` (Finding 103: the only gate that closed the
//! intersection, and the one Finding 101 showed hides an 81 % fall in the median cell). Its
//! **negative control** is Finding 99's B1 area cap, which must FAIL it.
//!
//! Run: cargo test -p ymir-core --release --test f104_retention -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, build_field_with_floor, pct, sorted};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone, propagation_order};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
/// Finding 88-D2's two published points, as used at Finding 97 B for the area cap.
const CAP_C: f32 = 308.9;

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The Flint intercept `k_s = median(S·√A)` over `A ≥ A_c` on a field's own channels.
fn k_s(f: &GridF32, on: &C1DrainageConfig, ss: &SteinSteinParams, cell_km2: f32) -> f32 {
    let (w, h) = (f.width, f.height);
    let d = c1_drainage_windowed(f, None, on, ss, DOMAIN_KM);
    let n2m = c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);
    let mut v = Vec::new();
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
            v.push(s * a.sqrt());
        }
    }
    pct(&sorted(v), 0.50)
}

#[test]
#[ignore]
fn f104_retention() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let s2 = splitmix64(PSEED);
    let s3 = splitmix64(s2);
    eprintln!(
        "\n==========  Finding 104 . channel-cut retention, with its negative control  =========="
    );
    eprintln!(
        "   ⚠️ every field is REBUILT: none is persisted. Population FIXED on the delivered \
         field's network (A >= A_c = {RELIEF_V1_A_C_KM2} km2)."
    );

    let t0 = Instant::now();
    for (si, seed) in [PSEED, s2, s3].into_iter().enumerate() {
        let pre = build_field_seed(Knobs::no_incision(), seed);
        let (w, h) = (pre.width, pre.height);
        let n = w * h;
        let delivered = build_field_seed(Knobs::passes(2), seed);
        let ks = k_s(&delivered, &on, &ss, cell_km2);

        // the FIXED classification, from the delivered field's own flow
        let d_del = c1_drainage_windowed(&delivered, None, &on, &ss, DOMAIN_KM);
        let is_chan: Vec<bool> = (0..n)
            .map(|k| {
                delivered.data[k] > SEA
                    && d_del.flow.accumulation.data[k] * cell_km2 >= RELIEF_V1_A_C_KM2
            })
            .collect();
        let land: Vec<bool> =
            (0..n).map(|k| pre.data[k] > SEA && delivered.data[k] > SEA).collect();
        let nc = is_chan.iter().filter(|&&b| b).count();
        let nh = (0..n).filter(|&k| land[k] && !is_chan[k]).count();

        let split = |f: &GridF32| -> (f64, f64) {
            let (mut c, mut hs) = (0f64, 0f64);
            for k in 0..n {
                if !land[k] || f.data[k] <= SEA {
                    continue;
                }
                let cut = ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64;
                if is_chan[k] {
                    c += cut;
                } else {
                    hs += cut;
                }
            }
            (c, hs)
        };
        let (c_del, h_del) = split(&delivered);
        eprintln!(
            "\n========  SEED {} ({seed})  ========\n   channels **{nc}** cells, hillslopes \
             **{nh}** . delivered cut: channels **{:.3e} m** ({:.1} m/cell), hillslopes \
             **{:.3e} m** ({:.1} m/cell) . k_s {ks:.4}",
            si + 1,
            c_del,
            c_del / nc.max(1) as f64,
            h_del,
            h_del / nh.max(1) as f64
        );
        eprintln!(
            "   **delivered** (positive control): channels **100.0 %** . hillslopes **100.0 %** . total 100.0 %"
        );

        // ── the negative control: Finding 99's B1 area cap, seed 1 only ─────
        if si == 0 {
            let d_pre = c1_drainage_windowed(&pre, None, &on, &ss, DOMAIN_KM);
            let br = breach_monotone(&pre, &d_pre.flow.filled, &d_pre.lake_map, SEA, w, h);
            let d_br = c1_drainage_windowed(&br, None, &on, &ss, DOMAIN_KM);
            let (_o, un) = propagation_order(&d_br.flow.direction, |k| br.data[k] > SEA, w, h);
            assert_eq!(un, 0);
            let q = (151.6f32 / CAP_C).ln() / (0.01f32).ln();
            let floor: Vec<f32> = (0..n)
                .map(|k| {
                    if pre.data[k] > SEA {
                        let a = (d_br.flow.accumulation.data[k] * cell_km2).max(cell_km2);
                        pre.data[k] - (CAP_C * a.powf(q)) / n2m
                    } else {
                        f32::NEG_INFINITY
                    }
                })
                .collect();
            let b1 = build_field_with_floor(Knobs::passes(2), Some(std::sync::Arc::new(floor)));
            let (c, hs) = split(&b1);
            eprintln!(
                "   ⛔ **B1 area cap** (NEGATIVE CONTROL, Finding 99: 43.7 % total): channels \
                 **{:.1} %** . hillslopes **{:.1} %** . total **{:.1} %** ⇒ gate at 80 % says \
                 **{}**",
                100.0 * c / c_del,
                100.0 * hs / h_del,
                100.0 * (c + hs) / (c_del + h_del),
                if 100.0 * c / c_del < 80.0 {
                    "FAIL (control fires)"
                } else {
                    "**PASS -- the gate does not detect it**"
                }
            );
        }

        for m in [0.4f32, 0.5, 0.6] {
            let f = build_field_seed(
                Knobs { slope_floor_uk: Some(ks * m), depression_floor: true, ..Knobs::passes(2) },
                seed,
            );
            let (c, hs) = split(&f);
            eprintln!(
                "   **x{m}**: channels **{:.1} %** . hillslopes **{:.1} %** . total **{:.1} %** \
                 ⇒ gate at 80 % says **{}**",
                100.0 * c / c_del,
                100.0 * hs / h_del,
                100.0 * (c + hs) / (c_del + h_del),
                if 100.0 * c / c_del >= 80.0 { "PASS" } else { "FAIL" }
            );
        }
    }
    eprintln!("\n==========  end Finding 104 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
