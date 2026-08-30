//! C-2 measurement bench (closures roadmap §2). Reconstruction WITH vs WITHOUT
//! volcanism at the production config, both resolutions — a RELATIVE before/after
//! (not the shipped product; the author validates the export visually). Reports:
//!   - closed depressions per stage, crater contribution as the WITH−WITHOUT delta
//!     (must be a small increment over the C-1 baseline, not a flood);
//!   - hydrological displacement (drainage density, lake/pit count, local relief) —
//!     intended around the cones (radial divergence, basin capture); the finding is
//!     any UNEXPLAINED shift;
//!   - young-vs-old hotspot edifice: flank slope + crater-bowl survival — the test
//!     of whether the relief-decay PROXY does measurable work.
//!
//! Run: cargo test -p ymir-core --test c2_volcanism_bench --release -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{
    Edifice, VolcanismConfig, VolcanoSetting, apply_edifices, place_edifices,
};

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const FULL_RANGE_M: f32 = 11302.0;

fn count_closed_depressions(field: &GridF32) -> usize {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let flow =
        compute_flow(field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let thr = 0.1f32 / FULL_RANGE_M;
    let mut seen = vec![false; n];
    let mut count = 0usize;
    let mut stack = Vec::new();
    for s in 0..n {
        if seen[s] || flow.filled.data[s] - field.data[s] <= thr {
            continue;
        }
        count += 1;
        seen[s] = true;
        stack.push(s);
        while let Some(k) = stack.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if !seen[nk] && flow.filled.data[nk] - field.data[nk] > thr {
                            seen[nk] = true;
                            stack.push(nk);
                        }
                    }
                }
            }
        }
    }
    count
}

/// Drainage density proxy: fraction of LAND cells whose flow accumulation exceeds
/// the 99th percentile (channel cells), ‰.
fn drainage_density_permille(field: &GridF32) -> f32 {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let flow =
        compute_flow(field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let mut acc = flow.accumulation.data.clone();
    acc.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let thr = acc[(acc.len() as f64 * 0.99) as usize];
    let (mut land, mut chan) = (0u64, 0u64);
    for k in 0..field.data.len() {
        if field.data[k] > 0.5 {
            land += 1;
            if flow.accumulation.data[k] >= thr {
                chan += 1;
            }
        }
    }
    1000.0 * chan as f32 / land.max(1) as f32
}

fn local_relief_median_m(field: &GridF32) -> f32 {
    let (w, h) = (field.width, field.height);
    let r = 5i32;
    let mut rel = Vec::new();
    for y in (r as usize..h - r as usize).step_by(4) {
        for x in (r as usize..w - r as usize).step_by(4) {
            if field.data[y * w + x] <= 0.5 {
                continue;
            }
            let (mut mn, mut mx) = (f32::MAX, f32::MIN);
            for dy in -r..=r {
                for dx in -r..=r {
                    let v = field.data[(y as i32 + dy) as usize * w + (x as i32 + dx) as usize];
                    mn = mn.min(v);
                    mx = mx.max(v);
                }
            }
            rel.push((mx - mn) * FULL_RANGE_M);
        }
    }
    rel.sort_by(|a, b| a.partial_cmp(b).unwrap());
    rel.get(rel.len() / 2).copied().unwrap_or(0.0)
}

/// Mean flank slope (degrees) of an edifice: sampled on the annulus between the
/// crater rim and the basal radius, on the eroded field, in physical units.
fn edifice_flank_slope_deg(
    field: &GridF32,
    cx: f32,
    cy: f32,
    rb_px: f32,
    rc_px: f32,
    km_per_cell: f32,
) -> f32 {
    let (w, h) = (field.width, field.height);
    let cell_m = km_per_cell * 1000.0;
    let (mut sum, mut n) = (0.0f64, 0u64);
    let (i0, i1) = ((cx - rb_px).max(1.0) as usize, ((cx + rb_px) as usize).min(w - 2));
    let (j0, j1) = ((cy - rb_px).max(1.0) as usize, ((cy + rb_px) as usize).min(h - 2));
    for j in j0..=j1 {
        for i in i0..=i1 {
            let r = (((i as f32 - cx).powi(2)) + ((j as f32 - cy).powi(2))).sqrt();
            if r < rc_px || r > rb_px {
                continue;
            }
            let gx = (field.data[j * w + i + 1] - field.data[j * w + i - 1]) * 0.5 * FULL_RANGE_M
                / cell_m;
            let gy =
                (field.data[(j + 1) * w + i] - field.data[(j - 1) * w + i]) * 0.5 * FULL_RANGE_M
                    / cell_m;
            sum += (gx * gx + gy * gy).sqrt().atan().to_degrees() as f64;
            n += 1;
        }
    }
    if n == 0 { 0.0 } else { (sum / n as f64) as f32 }
}

/// Crater-bowl depth (m) at the edifice centre: rim (max in a ring at rc) minus
/// the floor (min within rc). > 0 ⇒ the crater still holds a bowl (not breached).
fn crater_bowl_depth_m(field: &GridF32, cx: f32, cy: f32, rc_px: f32) -> f32 {
    let (w, h) = (field.width, field.height);
    let (mut floor, mut rim) = (f32::MAX, f32::MIN);
    let r = rc_px.max(1.0);
    let (i0, i1) = ((cx - r - 1.0).max(0.0) as usize, ((cx + r + 1.0) as usize).min(w - 1));
    let (j0, j1) = ((cy - r - 1.0).max(0.0) as usize, ((cy + r + 1.0) as usize).min(h - 1));
    for j in j0..=j1 {
        for i in i0..=i1 {
            let d = (((i as f32 - cx).powi(2)) + ((j as f32 - cy).powi(2))).sqrt();
            let v = field.data[j * w + i];
            if d <= r * 0.4 {
                floor = floor.min(v);
            }
            if (r * 0.8..=r * 1.2).contains(&d) {
                rim = rim.max(v);
            }
        }
    }
    if floor == f32::MAX || rim == f32::MIN { 0.0 } else { (rim - floor) * FULL_RANGE_M }
}

#[test]
#[ignore]
fn c2_volcanism_bench() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use ymir_core::tectonics_c1::kinematics::PlateKinematics;
    use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
    use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use ymir_core::terrain::upscale::{FbmUpscaleConfig, upscale_with_fbm};

    let ss = SteinSteinParams::default();
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});
    let coarse = c1_coarse_normalized_altitude(&state, &run.iso_config, &ss, None);
    let seed = WorldSeed::new(PSEED);

    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let (arc, hot, rift) = (
        edifices.iter().filter(|e| e.setting == VolcanoSetting::Arc).count(),
        edifices.iter().filter(|e| e.setting == VolcanoSetting::Hotspot).count(),
        edifices.iter().filter(|e| e.setting == VolcanoSetting::Rift).count(),
    );
    eprintln!("\n=== C-2 VOLCANISM BENCH (production seed) ===");
    eprintln!("placement: {arc} arc + {hot} hotspot + {rift} rift = {} edifices", edifices.len());
    eprintln!(
        "(arc volcanism is barely visible on this seed — look at the hotspot chains and rifts)"
    );

    for &target in &[2048usize, 8192usize] {
        let km_per_cell = DOMAIN_KM / target as f32;
        let mut cfg = FbmUpscaleConfig::c1_hd_production(target);
        cfg.amplitude_base = 0.04;
        cfg.sample_origin = [0.09375, 0.578125];
        cfg.sample_size = 1.0;
        cfg.stream_power = None;
        cfg.erosion = None;
        let mut sp =
            StreamPowerConfig::relief_v3(km_per_cell * km_per_cell, ss.depth_scale_m as f32);
        sp.mfd_exponent = Some(2.0);
        sp.iterations = 2;

        // WITHOUT volcanism.
        let fbm0 = upscale_with_fbm(&coarse, 0.5, &seed, &cfg).heightmap;
        let er0 = incise(&fbm0, &sp);
        // WITH volcanism (inject after FBM, before erosion — the production order).
        let mut fbm1 = upscale_with_fbm(&coarse, 0.5, &seed, &cfg).heightmap;
        let applied = apply_edifices(
            &mut fbm1,
            &edifices,
            cfg.sample_origin,
            cfg.sample_size,
            km_per_cell,
            FULL_RANGE_M,
            &volc,
        );
        let er1 = incise(&fbm1, &sp);

        eprintln!("\n--- {target}² ---  ({} crater records resolved)", applied.craters.len());
        eprintln!(
            "closed depressions  post-FBM: {:>6} → {:>6} (Δ {:+})  | post-relief: {:>6} → {:>6} (Δ {:+})",
            count_closed_depressions(&fbm0),
            count_closed_depressions(&fbm1),
            count_closed_depressions(&fbm1) as i64 - count_closed_depressions(&fbm0) as i64,
            count_closed_depressions(&er0),
            count_closed_depressions(&er1),
            count_closed_depressions(&er1) as i64 - count_closed_depressions(&er0) as i64,
        );
        eprintln!(
            "drainage density ‰: {:.1} → {:.1}  | local relief 11² m: {:.0} → {:.0}",
            drainage_density_permille(&er0),
            drainage_density_permille(&er1),
            local_relief_median_m(&er0),
            local_relief_median_m(&er1),
        );

        // Young vs old hotspot edifice — the relief-decay PROXY test. Take the first
        // contiguous hotspot chain; compare member 0 (young, active) to the last.
        if target == 8192 {
            young_vs_old(&edifices, &er1, cfg.sample_origin, cfg.sample_size, km_per_cell);
        }
    }
}

fn young_vs_old(edifices: &[Edifice], eroded: &GridF32, so: [f64; 2], ss: f64, km_per_cell: f32) {
    let (w, h) = (eroded.width, eroded.height);
    let hot: Vec<&Edifice> =
        edifices.iter().filter(|e| e.setting == VolcanoSetting::Hotspot).collect();
    if hot.len() < 2 {
        eprintln!("young-vs-old: no hotspot chain to compare");
        return;
    }
    let young = hot[0];
    let old = *hot.last().unwrap();
    let to_px = |e: &Edifice| -> Option<(f32, f32, f32, f32)> {
        let fx = (e.center_uv.0 - so[0] as f32).rem_euclid(1.0) / ss as f32;
        let fy = (e.center_uv.1 - so[1] as f32).rem_euclid(1.0) / ss as f32;
        if fx >= 1.0 || fy >= 1.0 {
            return None;
        }
        let cx = fx * w as f32;
        let cy = fy * h as f32;
        let rb = (e.basal_diameter_km * 0.5) / km_per_cell;
        let rc = (e.crater_diameter_km * 0.5) / km_per_cell;
        Some((cx, cy, rb, rc))
    };
    eprintln!("\n--- young vs old hotspot edifice (relief-decay PROXY test, 8192²) ---");
    for (label, e) in [("young(active)", young), ("old(extinct)", old)] {
        match to_px(e) {
            Some((cx, cy, rb, rc)) => {
                let slope = edifice_flank_slope_deg(eroded, cx, cy, rb, rc, km_per_cell);
                let bowl = crater_bowl_depth_m(eroded, cx, cy, rc);
                eprintln!(
                    "  {label:<13} age {:.2} active={} | flank slope {:.1}° | crater bowl {:.0} m | Wb {:.0} km H {:.0} m",
                    e.age_frac, e.active, slope, bowl, e.basal_diameter_km, e.height_m
                );
            }
            None => eprintln!("  {label:<13} outside the render window"),
        }
    }
    eprintln!(
        "  VERDICT: the decay proxy works iff the old edifice is measurably gentler and/or its crater\n\
         \x20 bowl is shallower (breached) than the young one; identical numbers = the proxy does no work."
    );
}
