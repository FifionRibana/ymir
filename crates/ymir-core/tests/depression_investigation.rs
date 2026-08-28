//! #190 INVESTIGATION (read-only) — closed-depression population: transient vs steady state.
//! Reproduces the EXACT production eroded field (coarse tectonics grid 64 -> upscale_with_fbm at
//! amplitude 0.04, offset [0.09375,0.578125] -> relief-v3 incise) and counts closed depressions after
//! each stage / iteration. Answers whether erosion maturity (more incision iterations) integrates the
//! drainage (transient -> a maturity knob works) or a process regenerates hollows every pass (steady).
//!
//! Run: cargo test -p ymir-core --test depression_investigation --release -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;

struct DepStats {
    count: usize,
    single_cell: usize,
    le2cell: usize,
    depths_m: Vec<f32>,
    areas: Vec<u32>,
    vols_norm: Vec<f64>,
    floors_norm: Vec<f32>,
}

fn count_depressions(field: &GridF32) -> DepStats {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let flow =
        compute_flow(field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let raise: Vec<f32> = (0..n).map(|k| flow.filled.data[k] - field.data[k]).collect();
    let thr = 0.1f32 / 11302.0; // 0.1 m in norm — below this is float/quantisation residue
    let mut seen = vec![false; n];
    let mut ds = DepStats {
        count: 0,
        single_cell: 0,
        le2cell: 0,
        depths_m: vec![],
        areas: vec![],
        vols_norm: vec![],
        floors_norm: vec![],
    };
    let mut stack = Vec::new();
    for s in 0..n {
        if seen[s] || raise[s] <= thr {
            continue;
        }
        let (mut area, mut maxr, mut vol, mut floor) = (0u32, 0f32, 0f64, f32::MAX);
        seen[s] = true;
        stack.push(s);
        while let Some(k) = stack.pop() {
            area += 1;
            maxr = maxr.max(raise[k]);
            vol += raise[k] as f64;
            floor = floor.min(field.data[k]);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if !seen[nk] && raise[nk] > thr {
                            seen[nk] = true;
                            stack.push(nk);
                        }
                    }
                }
            }
        }
        ds.count += 1;
        if area == 1 {
            ds.single_cell += 1;
        }
        if area <= 2 {
            ds.le2cell += 1;
        }
        ds.depths_m.push(maxr * 11302.0);
        ds.areas.push(area);
        ds.vols_norm.push(vol);
        ds.floors_norm.push(floor);
    }
    ds
}

fn dep_summary(label: &str, ds: &DepStats) {
    let land = ds.floors_norm.iter().filter(|&&f| f > 0.5).count();
    let below = ds.count - land;
    let mut d = ds.depths_m.clone();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p = |q: f32| d.get(((d.len() as f32 - 1.0).max(0.0) * q) as usize).copied().unwrap_or(0.0);
    let deep1 = ds.depths_m.iter().filter(|&&x| x >= 1.0).count();
    let deep10 = ds.depths_m.iter().filter(|&&x| x >= 10.0).count();
    let deep50 = ds.depths_m.iter().filter(|&&x| x >= 50.0).count();
    eprintln!(
        "  [{label:<26}] pits {:>6} (land {land}, below {below}) | 1cell {:>5} <=2cell {:>5} | depth med {:>5.1} p90 {:>6.1} max {:>6.0} m | >=1m {deep1} >=10m {deep10} >=50m {deep50}",
        ds.count,
        ds.single_cell,
        ds.le2cell,
        p(0.5),
        p(0.9),
        d.last().copied().unwrap_or(0.0)
    );
}

/// C-1 EXPORT AUDIT — read the SHIPPED rasters (the product, not a reconstruction)
/// and report the closed-depression population, slope shares, local relief and lake
/// stats on the conditioned 8192² export. The export is POST-breach, so land pits
/// should be near-zero and the surviving depressions should coincide with real
/// lakes (legitimate), not fabricated FBM dimples.
///
/// Run: cargo test -p ymir-core --test depression_investigation --release -- --ignored c1_export_audit --nocapture
/// Parse `min_m`/`max_m` of the height layer from a `manifest.json` (cheap scan).
fn manifest_min_max(dir: &str) -> (f32, f32) {
    let m = std::fs::read_to_string(format!("{dir}/manifest.json")).expect("manifest.json");
    let grab = |key: &str| -> f32 {
        // the height layer is the first block carrying both keys.
        let i = m.find(key).unwrap();
        let tail = &m[i + key.len()..];
        let s: String = tail
            .chars()
            .skip_while(|c| *c != '-' && !c.is_ascii_digit())
            .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '.' || *c == 'e')
            .collect();
        s.parse().unwrap()
    };
    (grab("\"min_m\""), grab("\"max_m\""))
}

/// C-1 EXPORT AUDIT — compare the SHIPPED products (the verdict, not a
/// reconstruction), un-conditioned vs conditioned (`_c1`), both climates. Reports
/// slope shares, local relief, and lake population per export. The height field is
/// post-breach (pit count ~0 either way — the interesting signal is the morphology
/// and how many lakes the breach had to open).
///
/// Run: cargo test -p ymir-core --test depression_investigation --release -- --ignored c1_export_audit --nocapture
#[test]
#[ignore]
fn c1_export_audit() {
    let base = concat!(env!("CARGO_MANIFEST_DIR"), "/../../exports");
    let cases = [
        ("humid OFF   ", "seed10481999410520546993_8192.humid.ymir"),
        ("humid  C-1  ", "seed10481999410520546993_8192_c1.humid.ymir"),
        ("arid-hot OFF", "seed10481999410520546993_8192.arid-hot.ymir"),
        ("arid-hot C-1", "seed10481999410520546993_8192_c1.arid-hot.ymir"),
    ];
    let n = 8192usize * 8192;
    eprintln!("\n=== C-1 EXPORT AUDIT (shipped 8192² products, post-breach) ===");
    for (label, name) in cases {
        let dir = format!("{base}/{name}");
        if std::fs::metadata(format!("{dir}/height.u16")).is_err() {
            eprintln!("  [{label}] (absent)");
            continue;
        }
        let (min_m, max_m) = manifest_min_max(&dir);
        let hb = std::fs::read(format!("{dir}/height.u16")).unwrap();
        let mut field = GridF32::new(8192, 8192, 0.0);
        for k in 0..n {
            let u = u16::from_le_bytes([hb[2 * k], hb[2 * k + 1]]) as f32;
            field.data[k] = 0.5 + (min_m + (u / 65535.0) * (max_m - min_m)) / 11302.0;
        }
        let lb = std::fs::read(format!("{dir}/lake_mask.u32")).unwrap();
        let lake_cells = (0..n)
            .filter(|&k| {
                u32::from_le_bytes([lb[4 * k], lb[4 * k + 1], lb[4 * k + 2], lb[4 * k + 3]]) != 0
            })
            .count();
        // distinct lake ids in lakes.json (one `"id"` per lake record).
        let lj = std::fs::read_to_string(format!("{dir}/lakes.json")).unwrap_or_default();
        let n_lakes = lj.matches("\"lake_type\"").count();
        let rj_bytes = std::fs::metadata(format!("{dir}/rivers.json")).map(|m| m.len()).unwrap_or(0);

        let (s15, s30, s45) = slope_shares_hd(&field, 11302.0, 400.0 / 8192.0);
        let lr = local_relief_median_m(&field, 11302.0);
        let land = (0..n).filter(|&k| field.data[k] > 0.5).count();
        eprintln!(
            "  [{label}] land {:>4.1}% | slope>15° {s15:>4.1}% >30° {s30:>4.1}% >45° {s45:>4.2}% | relief11² {lr:>3.0} m | lakes {n_lakes:>3} ({:>4.1}% cells) | rivers {} MB",
            100.0 * land as f64 / n as f64,
            100.0 * lake_cells as f64 / n as f64,
            rj_bytes / 1_000_000
        );
    }
}

#[test]
#[ignore]
fn depression_population_investigation() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise, incise_with_progress};
    use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use ymir_core::tectonics_c1::kinematics::PlateKinematics;
    use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
    use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use ymir_core::terrain::upscale::{FbmUpscaleConfig, upscale_with_fbm};
    const PSEED: u64 = 10481999410520546993;
    let target = 8192usize;
    let ss = SteinSteinParams::default();

    // 1. coarse tectonics (grid 64, n_steps 300 — the Viz-0 production defaults).
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

    // 2. upscale + FBM (NO erosion) — the pre-incision HD field.
    let seed = WorldSeed::new(PSEED);
    let mut cfg = FbmUpscaleConfig::c1_hd_production(target);
    cfg.amplitude_base = 0.04; // the production seed's FBM amplitude
    cfg.sample_origin = [0.09375, 0.578125]; // window offset (from the manifest)
    cfg.sample_size = 1.0;
    cfg.stream_power = None;
    cfg.erosion = None;
    let fbm = upscale_with_fbm(&coarse, 0.5, &seed, &cfg).heightmap;

    let km_per_cell = 400.0f32 / target as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let mut sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    sp.mfd_exponent = Some(2.0);

    eprintln!("\n=== #190 CLOSED-DEPRESSION POPULATION (production terrain, 8192^2) ===");
    eprintln!("\n--- STEP 1: population per stage ---");
    dep_summary("coarse post-isostasy", &count_depressions(&coarse));
    dep_summary("after FBM upscale", &count_depressions(&fbm));

    // 3. STEP 2 — iteration curve: 16 incision passes, count pits after each.
    eprintln!("\n--- STEP 2: iteration curve (relief-v3, per pass) ---");
    let mut sp16 = sp.clone();
    sp16.iterations = 16;
    let final16 = incise_with_progress(&fbm, &sp16, &mut |iter, f| {
        dep_summary(&format!("after incision iter {}", iter + 1), &count_depressions(f));
    });

    // 4. STEP 1 attribution — isolate each sub-process at the production 2 iterations.
    eprintln!("\n--- STEP 1b: stage attribution (2 iters, one process removed) ---");
    let mk = |f: &dyn Fn(&mut StreamPowerConfig)| {
        let mut s = sp.clone();
        s.iterations = 2;
        f(&mut s);
        incise(&fbm, &s)
    };
    dep_summary("full relief-v3 (2 iter)", &count_depressions(&mk(&|_| {})));
    dep_summary("NO MFD (single-flow)", &count_depressions(&mk(&|s| s.mfd_exponent = None)));
    dep_summary("NO talus", &count_depressions(&mk(&|s| s.talus_passes = 0)));
    dep_summary("NO hillslope diffusion", &count_depressions(&mk(&|s| s.diffusion = 0.0)));

    // 5. STEP 3 — size distribution of the FINAL (production 2-iter) field's hollows.
    let prod = mk(&|_| {});
    let ds = count_depressions(&prod);
    eprintln!("\n--- STEP 3: hollow size distribution (production 2-iter eroded, pre-breach) ---");
    dep_summary("PRODUCTION eroded", &ds);
    let mut a = ds.areas.clone();
    a.sort_unstable();
    let ap = |q: f32| a.get(((a.len() as f32 - 1.0).max(0.0) * q) as usize).copied().unwrap_or(0);
    let tot_vol_m3: f64 =
        ds.vols_norm.iter().map(|v| v * 11302.0 * (km_per_cell as f64 * 1000.0).powi(2)).sum();
    eprintln!(
        "  area cells: med {} p90 {} max {} | <=2-cell share {:.0}% | total fill volume {:.3e} m3",
        ap(0.5),
        ap(0.9),
        a.last().copied().unwrap_or(0),
        100.0 * ds.le2cell as f32 / ds.count.max(1) as f32,
        tot_vol_m3
    );

    // 6. VALIDATION against the shipped export (avoid the reconstructed-terrain trap). The export is
    //    POST-breach; compare on land cells the breach mostly leaves alone.
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../exports/seed10481999410520546993_8192.humid.ymir"
    );
    let hb = std::fs::read(format!("{dir}/height.u16")).unwrap();
    let (min_m, max_m) = (-5505.853515625f32, 3917.749267578125f32);
    let (mut match_cells, mut land_cells, mut sumabs) = (0u64, 0u64, 0f64);
    for k in 0..target * target {
        let u = u16::from_le_bytes([hb[2 * k], hb[2 * k + 1]]) as f32;
        let export_m = min_m + (u / 65535.0) * (max_m - min_m);
        let repro2_m = (prod.data[k] - 0.5) * 11302.0;
        if prod.data[k] > 0.5 {
            land_cells += 1;
            let d = (repro2_m - export_m).abs();
            sumabs += d as f64;
            if d < 5.0 {
                match_cells += 1;
            }
        }
    }
    eprintln!(
        "\n--- VALIDATION vs export (2-iter pre-breach vs post-breach export) ---\n  land cells {land_cells}: within 5 m of export {:.1}% | mean abs diff {:.1} m",
        100.0 * match_cells as f64 / land_cells.max(1) as f64,
        sumabs / land_cells.max(1) as f64
    );
    let _ = final16;
}

/// C-1 (closures roadmap §1) — the flow-conditioned FBM measured in the production
/// configuration at BOTH resolutions. Reports the closed-depression population per
/// stage (coarse / post-FBM / post-relief) for the un-conditioned baseline and for
/// the relief-budget conditioning, then sweeps `flow_conditioning` (β) at 8192² to
/// locate the value that brings the post-FBM count to the same ORDER as the 16
/// tectonic depressions. ACCEPTANCE = post-FBM count of order 16 at both resolutions.
///
/// Run: cargo test -p ymir-core --test depression_investigation --release -- --ignored c1_flow_conditioning --nocapture
/// Land-cell slope shares (>15° / >30° / >45°) on an HD field. `full_range_m` is
/// the norm→metre span (11302 for the HD field), `km_per_cell` the horizontal cell.
fn slope_shares_hd(field: &GridF32, full_range_m: f32, km_per_cell: f32) -> (f32, f32, f32) {
    let (w, h) = (field.width, field.height);
    let cell_m = km_per_cell * 1000.0;
    let (mut land, mut c15, mut c30, mut c45) = (0u64, 0u64, 0u64, 0u64);
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if field.data[k] <= 0.5 {
                continue;
            }
            land += 1;
            let gx = (field.data[k + 1] - field.data[k - 1]) * 0.5 * full_range_m / cell_m;
            let gy = (field.data[k + w] - field.data[k - w]) * 0.5 * full_range_m / cell_m;
            let deg = (gx * gx + gy * gy).sqrt().atan().to_degrees();
            if deg > 15.0 {
                c15 += 1;
            }
            if deg > 30.0 {
                c30 += 1;
            }
            if deg > 45.0 {
                c45 += 1;
            }
        }
    }
    let f = land.max(1) as f32;
    (c15 as f32 / f * 100.0, c30 as f32 / f * 100.0, c45 as f32 / f * 100.0)
}

/// Median local relief (m) in an 11×11 window over land cells — the inter-channel
/// roughness a conditioning that over-flattens would collapse.
fn local_relief_median_m(field: &GridF32, full_range_m: f32) -> f32 {
    let (w, h) = (field.width, field.height);
    let r = 5i32;
    let mut rel: Vec<f32> = Vec::new();
    // subsample (every 4th cell) to keep the 8192² pass affordable — a median proxy.
    for y in (r as usize..h - r as usize).step_by(4) {
        for x in (r as usize..w - r as usize).step_by(4) {
            let k = y * w + x;
            if field.data[k] <= 0.5 {
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
            rel.push((mx - mn) * full_range_m);
        }
    }
    if rel.is_empty() {
        return 0.0;
    }
    rel.sort_by(|a, b| a.partial_cmp(b).unwrap());
    rel[rel.len() / 2]
}

/// C-1 shape-metrics regression — the conditioning must not over-flatten the
/// mountains. Reports slope shares + local relief for OFF vs the production β at
/// BOTH resolutions, on the post-relief (product) field.
///
/// Run: cargo test -p ymir-core --test depression_investigation --release -- --ignored c1_shape_metrics --nocapture
#[test]
#[ignore]
fn c1_shape_metrics() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use ymir_core::tectonics_c1::kinematics::PlateKinematics;
    use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
    use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use ymir_core::terrain::upscale::{FbmUpscaleConfig, upscale_with_fbm};
    const PSEED: u64 = 10481999410520546993;
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
    let mk_cfg = |target: usize, beta: f64| {
        let mut cfg = FbmUpscaleConfig::c1_hd_production(target);
        cfg.amplitude_base = 0.04;
        cfg.sample_origin = [0.09375, 0.578125];
        cfg.sample_size = 1.0;
        cfg.stream_power = None;
        cfg.erosion = None;
        cfg.flow_conditioning = beta;
        cfg
    };
    eprintln!("\n=== C-1 SHAPE METRICS (post-relief product field) ===");
    for &target in &[2048usize, 8192usize] {
        let km_per_cell = 400.0f32 / target as f32;
        let cell_km2 = km_per_cell * km_per_cell;
        let mut sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
        sp.mfd_exponent = Some(2.0);
        sp.iterations = 2;
        eprintln!("\n--- {target}² ---");
        for (label, beta) in [("OFF", 0.0), ("β=0.1", 0.1)] {
            let fbm = upscale_with_fbm(&coarse, 0.5, &seed, &mk_cfg(target, beta)).heightmap;
            let eroded = incise(&fbm, &sp);
            let (s15, s30, s45) = slope_shares_hd(&eroded, 11302.0, km_per_cell);
            let lr = local_relief_median_m(&eroded, 11302.0);
            let pits = count_depressions(&eroded).count;
            eprintln!(
                "  [{label:<6}] slope>15° {s15:>5.1}% >30° {s30:>4.1}% >45° {s45:>4.2}% | local relief(11²) med {lr:>5.0} m | pits {pits:>6}"
            );
        }
    }
}

#[test]
#[ignore]
fn c1_flow_conditioning_sweep() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use ymir_core::tectonics_c1::kinematics::PlateKinematics;
    use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
    use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use ymir_core::terrain::upscale::{FbmUpscaleConfig, upscale_with_fbm};
    const PSEED: u64 = 10481999410520546993;
    let ss = SteinSteinParams::default();

    // Coarse tectonics — resolution-independent, built once (the Viz-0 defaults).
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

    // The production FBM config for this seed (window offset from the manifest).
    let mk_cfg = |target: usize, beta: f64| {
        let mut cfg = FbmUpscaleConfig::c1_hd_production(target);
        cfg.amplitude_base = 0.04; // the production seed's FBM amplitude
        cfg.sample_origin = [0.09375, 0.578125];
        cfg.sample_size = 1.0;
        cfg.stream_power = None;
        cfg.erosion = None;
        cfg.flow_conditioning = beta;
        cfg
    };

    eprintln!("\n=== C-1 FLOW-CONDITIONED FBM (production terrain) ===");
    let coarse_pits = count_depressions(&coarse).count;
    eprintln!("coarse post-isostasy (both resolutions share it): {coarse_pits} closed depressions");

    // --- Per-stage table at both resolutions: OFF vs conditioned (β = 1.0). ---
    for &target in &[2048usize, 8192usize] {
        let km_per_cell = 400.0f32 / target as f32;
        let cell_km2 = km_per_cell * km_per_cell;
        let mut sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
        sp.mfd_exponent = Some(2.0);
        sp.iterations = 2;

        eprintln!("\n--- {target}² per-stage population ---");
        for (label, beta) in [("OFF (additive)", 0.0), ("conditioned β=0.1", 0.1)] {
            let fbm = upscale_with_fbm(&coarse, 0.5, &seed, &mk_cfg(target, beta)).heightmap;
            let post_fbm = count_depressions(&fbm);
            let eroded = incise(&fbm, &sp);
            let post_relief = count_depressions(&eroded);
            eprintln!(
                "  [{label:<18}] post-FBM {:>6} | post-relief(2it) {:>6}",
                post_fbm.count, post_relief.count
            );
            dep_summary(&format!("{target} {label} post-FBM   "), &post_fbm);
        }
    }

    // --- β sweep at 8192² (post-FBM only — the acceptance metric). ---
    eprintln!(
        "\n--- 8192² relief-budget sweep (post-FBM closed-depression count vs tectonic 16) ---"
    );
    for beta in [1.0, 0.4, 0.2, 0.1, 0.05, 0.02, 0.01] {
        let fbm = upscale_with_fbm(&coarse, 0.5, &seed, &mk_cfg(8192, beta)).heightmap;
        let ds = count_depressions(&fbm);
        eprintln!("  β = {beta:>4.2} : post-FBM pits {:>6}", ds.count);
    }
}
