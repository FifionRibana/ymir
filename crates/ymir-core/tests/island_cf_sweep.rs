//! M1 #190 PART B — continental_fraction sweep. Reduce the landmass by having
//! LESS continental crust (fewer continental plates) instead of drowning the
//! world: sea level stays at the isostatic level (steep hypsometric curve → crisp
//! coasts) while the mass is bounded by construction. `round(num_plates · cf)`
//! continental plates; `= 1` is a single Voronoi cell that cannot span the torus.
//!
//! target_land_fraction = None throughout (the whole point); cc = 1.
//! Run: cargo test -p ymir-core --test island_cf_sweep --release -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::land_topology::evaluate_island;
use ymir_core::tectonics_c1::production_upscale::{c1_coarse_raw_altitude, c1_normalize_coarse};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};

fn slope_at(sorted: &[f32], thr: f32, win: f32) -> f32 {
    let lo = sorted.partition_point(|&v| v < thr - win);
    let hi = sorted.partition_point(|&v| v < thr + win);
    (hi - lo) as f32 / (2.0 * win)
}
fn median(v: &mut [f32]) -> f32 {
    if v.is_empty() {
        return f32::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

struct Row {
    np: usize,
    cf: f64,
    n_cont: usize,
    seed: u64,
    emerged: f32,
    masses: usize,
    largest: f32,
    traverse: f32,
    wx: bool,
    wy: bool,
    clean: bool,
    window_km: f32,
    m_per_cell: f32,
    compact: f32,
    slope: f32,
}

#[test]
#[ignore]
fn sweep_continental_fraction() {
    let grid = 64usize;
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / grid as f64,
        dy: 1.0 / grid as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let ss = SteinSteinParams::default();
    let clo = C1Closures::default();
    let seeds: Vec<u64> = (0..20).collect();
    let margin_km = 25.0f32;

    // (num_plates, [continental_fraction so round(np·cf) ∈ {1, 2}])
    let configs: [(usize, [f64; 2]); 4] = [
        (8, [1.0 / 8.0, 2.0 / 8.0]),
        (12, [1.0 / 12.0, 2.0 / 12.0]),
        (16, [1.0 / 16.0, 2.0 / 16.0]),
        (20, [1.0 / 20.0, 2.0 / 20.0]),
    ];

    let mut rows: Vec<Row> = Vec::new();
    for (np, cfs) in configs {
        for cf in cfs {
            let n_cont = (np as f64 * cf).round() as usize;
            for &seed in &seeds {
                let mut init = Phase2InitParams::default();
                init.num_plates = np;
                init.cluster.seed_cluster_count = 1;
                init.cluster.continental_fraction = cf;
                let mut state = init_c1_state_phase_2_r7(grid, seed, &init);
                let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
                run_with_closures(&mut state, &mut kin, &run, &clo, |_, _| {});
                let raw: GridF32 = c1_coarse_raw_altitude(&state, &run.iso_config, &ss);
                let mut sorted = raw.data.clone();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let slope = slope_at(&sorted, 0.0, 0.02); // sea = raw 0 (tlf None)
                let norm = c1_normalize_coarse(raw, None);
                let e = evaluate_island(&norm, 0.5, margin_km, 1);
                let t = &e.topo;
                rows.push(Row {
                    np,
                    cf,
                    n_cont,
                    seed,
                    emerged: t.emerged_fraction,
                    masses: t.num_landmasses,
                    largest: t.largest_area_km2,
                    traverse: t.bbox_km.0.max(t.bbox_km.1),
                    wx: t.wraps_x,
                    wy: t.wraps_y,
                    clean: e.accepted(),
                    window_km: e.window_km,
                    m_per_cell: e.window_km / 8192.0 * 1000.0,
                    compact: e.compactness,
                    slope,
                });
            }
        }
    }

    // CSV (Part B schema — new file so the option-3 CSV stays intact).
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("island_sweep_partB.csv");
    let mut csv = String::from(
        "num_plates,continental_fraction,n_continental_plates,seed,emerged_frac,num_landmasses,\
         largest_area_km2,traverse_km,wraps_x,wraps_y,border_clean,window_km,m_per_cell,\
         compactness,hypso_slope_at_sea\n",
    );
    for r in &rows {
        csv.push_str(&format!(
            "{},{:.4},{},{},{:.4},{},{:.0},{:.0},{},{},{},{:.0},{:.2},{:.3},{:.0}\n",
            r.np,
            r.cf,
            r.n_cont,
            r.seed,
            r.emerged,
            r.masses,
            r.largest,
            r.traverse,
            r.wx,
            r.wy,
            r.clean,
            r.window_km,
            r.m_per_cell,
            r.compact,
            r.slope
        ));
    }
    std::fs::write(&path, csv).unwrap();
    eprintln!(
        "rows: {} — CSV: {}\n(tlf None throughout; cc 1; margin 25 km; isostatic slope ref ≈ 12350)",
        rows.len(),
        path.display()
    );

    // Per (np, cf) summary.
    for (np, cfs) in configs {
        for cf in cfs {
            let g: Vec<&Row> =
                rows.iter().filter(|r| r.np == np && (r.cf - cf).abs() < 1e-9).collect();
            let n_cont = g[0].n_cont;
            let nonwrap = g.iter().filter(|r| !(r.wx || r.wy)).count();
            let clean: Vec<&&Row> = g.iter().filter(|r| r.clean).collect();
            let mut em: Vec<f32> = g.iter().map(|r| r.emerged * 100.0).collect();
            let mut sl: Vec<f32> = g.iter().map(|r| r.slope).collect();
            let best = clean.iter().max_by(|a, b| a.largest.partial_cmp(&b.largest).unwrap());
            eprint!(
                "np {np:>2} cf {cf:.3} → {n_cont} cont. plate(s): non-wrap {nonwrap:>2}/20, \
                 border-clean {:>2}, median emerged {:>4.1}%, median slope {:>6.0}",
                clean.len(),
                median(&mut em),
                median(&mut sl),
            );
            match best {
                Some(b) => eprintln!(
                    " | best clean seed {} → {:.0} km² ({:.0}%), traverse {:.0} km, {:.1} m/cell, compact {:.2}",
                    b.seed,
                    b.largest,
                    b.emerged * 100.0,
                    b.traverse,
                    b.m_per_cell,
                    b.compact
                ),
                None => eprintln!(" | best clean —"),
            }
        }
    }
}
