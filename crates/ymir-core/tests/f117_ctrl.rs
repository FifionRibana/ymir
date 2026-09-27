mod common;
use common::{Knobs, PSEED, build_field_seed};
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
fn pits(fill: &GridF32, raw: &GridF32, n2m: f32, t: f32) -> usize {
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let deep: Vec<bool> = (0..n).map(|k| (fill.data[k] - raw.data[k]) * n2m > t).collect();
    let mut seen = vec![false; n];
    let mut count = 0usize;
    for s in 0..n {
        if !deep[s] || seen[s] {
            continue;
        }
        let mut q = vec![s];
        seen[s] = true;
        while let Some(k) = q.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nk = (y + dy).rem_euclid(h as i32) as usize * w
                        + (x + dx).rem_euclid(w as i32) as usize;
                    if deep[nk] && !seen[nk] {
                        seen[nk] = true;
                        q.push(nk);
                    }
                }
            }
        }
        count += 1;
    }
    count
}
/// ADR Finding 117 — rule 18 control for B2: does the INCISION add pits without the age closure?
#[test]
#[ignore]
fn f117_incision_off_closure() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    for (label, knobs) in [
        ("closure OFF, full K", Knobs::passes(2)),
        (
            "closure OFF, uniform K",
            Knobs {
                lithology_off: true,
                fracture_off: true,
                bathymetry_off: true,
                ..Knobs::passes(2)
            },
        ),
    ] {
        let f = build_field_seed(knobs, PSEED);
        let dr = c1_drainage_windowed(&f, None, &dcfg, &ss, 400.0);
        eprintln!(
            "   RULE-18 · {label}: pits > 0.1 m **{}** · > 20 m **{}**   (S3 pre-incision: 6 895 / 49 · S7 closure ON: 10 867 / 1 383)",
            pits(&dr.flow.filled, &f, n2m, 0.1),
            pits(&dr.flow.filled, &f, n2m, 20.0)
        );
    }
}
