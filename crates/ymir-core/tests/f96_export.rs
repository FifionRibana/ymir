//! ADR 0001 Finding 96 — the visual validation of the best candidate. **No production change.**
//!
//! Finding 96's table C leaves ONE candidate standing on relief, texture, the coastal criterion
//! and cost: **B3, the χ-profile floor**. It fails the canyon class, and the reason is named in
//! the finding — χ is built on the BREACHED pre-incision field, so inside the footprints the
//! breach carved through, χ is low and the floor permits the cut it should stop.
//!
//! This bench exists to put that candidate in front of an eye rather than a table. It builds the
//! delivered field and the B3 field, exports both plus their signed difference, and costs
//! **~3 minutes** instead of the 1 h 43 the Finding 95 oracle field would cost — which the author
//! has already refused once.
//!
//! Run: cargo test -p ymir-core --release --test f96_export -- --ignored --nocapture

mod common;

use common::{Knobs, SEA, build_field, build_field_with_floor};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone, propagation_order};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
/// A test's working directory is the CRATE root, not the repo root -- the first run of this bench
/// wrote its PNGs to `crates/ymir-core/docs/...`. Anchor on `CARGO_MANIFEST_DIR`.
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f96_closure");
/// Finding 95's oracle, as a number: the paired relief p50 the χ calibration targets.
const ORACLE_P50_M: f32 = 488.1;

/// Land-only grey ramp: sea flat black, land stretched over its own p1..p99 so the valleys are
/// legible instead of being crushed by the volcanic maxima.
fn land_ramp(f: &GridF32, lo: f32, hi: f32) -> GridF32 {
    let mut g = GridF32::new(f.width, f.height, 0.0);
    for k in 0..f.width * f.height {
        g.data[k] = if f.data[k] <= SEA {
            0.0
        } else {
            0.08 + 0.92 * ((f.data[k] - lo) / (hi - lo)).clamp(0.0, 1.0)
        };
    }
    g
}

fn pctile(v: &mut Vec<f32>, p: f64) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if v.is_empty() { 0.0 } else { v[(((v.len() - 1) as f64) * p) as usize] }
}

fn crop(f: &GridF32, x0: usize, y0: usize, side: usize) -> GridF32 {
    let mut g = GridF32::new(side, side, 0.0);
    for y in 0..side {
        for x in 0..side {
            g.data[y * side + x] = f.data[(y0 + y) * f.width + (x0 + x)];
        }
    }
    g
}

fn save(g: &GridF32, name: &str) {
    let dir = Path::new(OUT);
    std::fs::create_dir_all(dir).expect("create the report directory");
    let p = dir.join(name);
    g.save_png_u8(&p).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!("   wrote {}", p.display());
}

#[test]
#[ignore]
fn f96_export() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = (DOMAIN_KM / 8192.0) * (DOMAIN_KM / 8192.0);
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 96 · the visual validation of B3  ==========");

    // ── the χ floor, rebuilt exactly as `f96_closure` builds it ─────────────
    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let d_pre = c1_drainage_windowed(&pre, None, &on, &ss, DOMAIN_KM);
    let br = breach_monotone(&pre, &d_pre.flow.filled, &d_pre.lake_map, SEA, w, h);
    let d_br = c1_drainage_windowed(&br, None, &on, &ss, DOMAIN_KM);
    let land_pre: Vec<bool> = (0..n).map(|k| br.data[k] > SEA).collect();
    let (order, unreached) = propagation_order(&d_br.flow.direction, |k| land_pre[k], w, h);
    assert_eq!(unreached, 0, "the χ traversal must reach every land cell");

    let a_c = RELIEF_V1_A_C_KM2;
    let mut chi = vec![0.0f32; n];
    for &ku in order.iter().rev() {
        let k = ku as usize;
        if !land_pre[k] {
            continue;
        }
        let d = d_br.flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let r = ny * w + nx;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let a_km2 = (d_br.flow.accumulation.data[k] * cell_km2).max(a_c);
        chi[k] = if land_pre[r] { chi[r] } else { 0.0 } + (1.0 / a_km2).sqrt() * dx_m;
    }
    // the single free parameter, calibrated to Finding 95's oracle median
    let mut cs: Vec<f32> = (0..n).filter(|&k| land_pre[k]).map(|k| chi[k]).collect();
    let c_uk = ORACLE_P50_M / pctile(&mut cs, 0.50).max(1e-6);
    let z_chi: Vec<f32> = (0..n)
        .map(|k| if land_pre[k] { SEA + (c_uk * chi[k]) / n2m } else { br.data[k] })
        .collect();
    eprintln!("   χ floor rebuilt · (U/K)^(1/n) = {c_uk:.4} · {:.1} s", t0.elapsed().as_secs_f64());

    // ── the two fields ──────────────────────────────────────────────────────
    let floor = std::sync::Arc::new(z_chi);
    let delivered = build_field(Knobs::passes(2));
    let b3 = build_field_with_floor(Knobs::passes(2), Some(floor.clone()));
    // the union: A1's exclusion of the closed depressions ON TOP of B3's floor -- the candidate
    // Finding 96's table C names, because chi cannot see the population A1 excludes.
    let union =
        build_field_with_floor(Knobs { depression_floor: true, ..Knobs::passes(2) }, Some(floor));

    // one shared ramp for both, so the images are comparable rather than each self-normalised
    let mut land_m: Vec<f32> =
        (0..n).filter(|&k| delivered.data[k] > SEA).map(|k| delivered.data[k]).collect();
    let hi = pctile(&mut land_m, 0.99);
    let lo = pctile(&mut land_m, 0.01);
    eprintln!(
        "   shared ramp: {:.0} m .. {:.0} m of altitude",
        c1_altitude_norm_to_metres(lo, &ss),
        c1_altitude_norm_to_metres(hi, &ss)
    );
    save(&land_ramp(&delivered, lo, hi), "delivered_8192.png");
    save(&land_ramp(&b3, lo, hi), "b3_chi_floor_8192.png");
    save(&land_ramp(&union, lo, hi), "a1_b3_union_8192.png");

    // the signed difference, centred at mid-grey. The floor can only STOP incision, so the
    // difference is expected positive almost everywhere -- but not provably so at the cell level:
    // the diffusion and the talus run AFTER the incision and redistribute mass, so a cell can end
    // up lower than in the unbounded run. The ramp is therefore two-sided on purpose.
    let mut dif = GridF32::new(w, h, 0.5);
    let mut dv: Vec<f32> = Vec::new();
    for k in 0..n {
        if delivered.data[k] > SEA || b3.data[k] > SEA {
            dv.push((b3.data[k] - delivered.data[k]) * n2m);
        }
    }
    let span = pctile(&mut dv, 0.99).max(1.0);
    for k in 0..n {
        let d = (b3.data[k] - delivered.data[k]) * n2m;
        dif.data[k] = (0.5 + 0.5 * (d / span)).clamp(0.0, 1.0);
    }
    save(&dif, "difference_b3_minus_delivered_8192.png");
    eprintln!(
        "   difference: p50 {:.1} m · p99 {:.1} m · max {:.1} m (mid-grey = 0, white = +{span:.0} m)",
        pctile(&mut dv, 0.50),
        pctile(&mut dv, 0.99),
        dv.last().copied().unwrap_or(0.0)
    );

    // ── and a 1024² crop where the floor bit hardest, for both fields ───────
    let mut best = (0usize, 0usize, -1.0f32);
    for by in (0..h - 1024).step_by(1024) {
        for bx in (0..w - 1024).step_by(1024) {
            let mut s = 0.0f32;
            for y in (by..by + 1024).step_by(8) {
                for x in (bx..bx + 1024).step_by(8) {
                    s += (b3.data[y * w + x] - delivered.data[y * w + x]).max(0.0);
                }
            }
            if s > best.2 {
                best = (bx, by, s);
            }
        }
    }
    eprintln!("   the crop the floor changed most: ({}, {}) 1024²", best.0, best.1);
    save(&crop(&land_ramp(&delivered, lo, hi), best.0, best.1, 1024), "crop_delivered_1024.png");
    save(&crop(&land_ramp(&b3, lo, hi), best.0, best.1, 1024), "crop_b3_chi_floor_1024.png");
    save(&crop(&land_ramp(&union, lo, hi), best.0, best.1, 1024), "crop_a1_b3_union_1024.png");
    eprintln!("\n==========  end · total {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
