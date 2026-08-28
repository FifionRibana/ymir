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
