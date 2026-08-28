//! M1 #190 — framing-roll acceptance guard (regression).
//! A seed whose largest landmass straddles the torus seam must export as ONE
//! contiguous continent, not two/four fragments split across the border. The roll
//! is applied at the coarse SAMPLING ORIGIN (upscale `sample_origin`); this test
//! reproduces exactly what `run_hd` does and checks the exported grid.
//!
//! Run: cargo test -p ymir-core --test rolled_export_seam --release -- --ignored --nocapture

use std::collections::VecDeque;

use ymir_core::export::vector::coastline_geojson;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::cached_product::cached_c1_eroded;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::land_topology::land_topology;
use ymir_core::tectonics_c1::production_upscale::{c1_coarse_raw_altitude, c1_normalize_coarse};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const SEA: f32 = 0.5;

/// Largest NON-PERIODIC 4-connected land component fraction (the exported grid is
/// NOT a torus — a split continent shows up as several smaller components).
fn largest_component_frac(g: &GridF32) -> f32 {
    let (w, h) = (g.width, g.height);
    let n = w * h;
    let is_land = |k: usize| g.data[k] > SEA;
    let mut seen = vec![false; n];
    let mut best = 0usize;
    for start in 0..n {
        if !is_land(start) || seen[start] {
            continue;
        }
        let mut size = 0usize;
        let mut q = VecDeque::new();
        q.push_back(start);
        seen[start] = true;
        while let Some(k) = q.pop_front() {
            size += 1;
            let (x, y) = (k % w, k / w);
            // NON-periodic neighbours (bounded map, no wrap).
            if x + 1 < w && is_land(k + 1) && !seen[k + 1] {
                seen[k + 1] = true;
                q.push_back(k + 1);
            }
            if x > 0 && is_land(k - 1) && !seen[k - 1] {
                seen[k - 1] = true;
                q.push_back(k - 1);
            }
            if y + 1 < h && is_land(k + w) && !seen[k + w] {
                seen[k + w] = true;
                q.push_back(k + w);
            }
            if y > 0 && is_land(k - w) && !seen[k - w] {
                seen[k - w] = true;
                q.push_back(k - w);
            }
        }
        best = best.max(size);
    }
    best as f32 / n as f32
}

/// Number of land cells on the four map borders (a split continent touches both
/// the left and right — or top and bottom — edges).
fn land_on_borders(g: &GridF32) -> usize {
    let (w, h) = (g.width, g.height);
    let d = &g.data;
    let mut c = 0;
    for x in 0..w {
        if d[x] > SEA {
            c += 1;
        }
        if d[(h - 1) * w + x] > SEA {
            c += 1;
        }
    }
    for y in 0..h {
        if d[y * w] > SEA {
            c += 1;
        }
        if d[y * w + (w - 1)] > SEA {
            c += 1;
        }
    }
    c
}

/// Longest LineString in a coastline FeatureCollection, and whether it is a closed
/// ring (first vertex == last vertex).
fn longest_ring_closed(geojson: &[u8]) -> (usize, bool) {
    let v: serde_json::Value = serde_json::from_slice(geojson).unwrap();
    let mut best_len = 0usize;
    let mut best_closed = false;
    let mut consider = |coords: &[serde_json::Value]| {
        if coords.len() >= 2 && coords.len() > best_len {
            best_len = coords.len();
            best_closed = coords[0] == coords[coords.len() - 1];
        }
    };
    let empty: Vec<serde_json::Value> = Vec::new();
    for feat in v["features"].as_array().unwrap_or(&empty) {
        let geom = &feat["geometry"];
        match geom["type"].as_str() {
            Some("LineString") => {
                if let Some(c) = geom["coordinates"].as_array() {
                    consider(c);
                }
            }
            Some("MultiLineString") => {
                for line in geom["coordinates"].as_array().unwrap_or(&empty) {
                    if let Some(c) = line.as_array() {
                        consider(c);
                    }
                }
            }
            _ => {}
        }
    }
    (best_len, best_closed)
}

#[test]
#[ignore]
fn rolled_export_is_one_contiguous_continent() {
    let grid = 64usize;
    let seed = 42u64; // largest mass straddles the x seam (straddles_x)
    let target = 1024usize; // no erosion → fast; roll/contiguity is resolution-agnostic
    let init = Phase2InitParams::default(); // pre-M1: 8 plates, cc 1
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / grid as f64,
        dy: 1.0 / grid as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let clo = C1Closures::default();
    let ss = SteinSteinParams::default();

    // Preview topology (unrolled coarse, tlf None — matches c1_hd_production).
    let mut state = init_c1_state_phase_2_r7(grid, seed, &init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &clo, |_, _| {});
    let coarse = c1_normalize_coarse(c1_coarse_raw_altitude(&state, &run.iso_config, &ss), None);
    let topo = land_topology(&coarse, SEA);
    eprintln!(
        "preview: largest {:.0} km² ({:.1}%), traverse {:.0} km, straddle x={} y={}, band x={} y={}, centre {:?}",
        topo.largest_area_km2,
        topo.largest_area_frac * 100.0,
        topo.bbox_km.0.max(topo.bbox_km.1),
        topo.straddles_x,
        topo.straddles_y,
        topo.wraps_x,
        topo.wraps_y,
        topo.center_cell,
    );
    assert!(
        topo.straddles_x || topo.straddles_y,
        "seed {seed} must straddle a seam for this guard"
    );

    // The roll run_hd computes: centre the mass, integer coarse cells.
    let g = grid as i64;
    let roll_x = ((topo.center_cell.0 as i64) - g / 2).rem_euclid(g) as usize;
    let roll_y = ((topo.center_cell.1 as i64) - g / 2).rem_euclid(g) as usize;
    let offset = [roll_x as f64 / grid as f64, roll_y as f64 / grid as f64];
    eprintln!(
        "roll ({roll_x},{roll_y}) cells → window_offset_in_torus [{:.3},{:.3}]",
        offset[0], offset[1]
    );
    assert_ne!(offset, [0.0, 0.0], "a straddling seed must have a non-zero framing offset");

    let cache = std::env::temp_dir().join("ymir_rolled_export_seam");
    let _ = std::fs::create_dir_all(&cache);

    let export = |sample_origin: [f64; 2]| -> GridF32 {
        let mut up = FbmUpscaleConfig::c1_hd_production(target);
        up.erosion = None; // speed; contiguity is set by the coarse sampling + roll
        up.sample_origin = sample_origin;
        up.sample_size = 1.0;
        cached_c1_eroded(&cache, seed, grid, &init, &run, &clo, &ss, &up).unwrap()
    };

    // Baseline: NO roll (origin [0,0]) → the straddler is split across the border.
    let unrolled = export([0.0, 0.0]);
    let unrolled_frac = largest_component_frac(&unrolled);

    // Rolled export (what run_hd ships).
    let rolled = export(offset);
    let rolled_frac = largest_component_frac(&rolled);
    let periodic_frac = land_topology(&rolled, SEA).largest_area_frac;

    eprintln!(
        "largest NON-periodic land component: unrolled {:.1}% vs rolled {:.1}% (periodic/torus {:.1}%)",
        unrolled_frac * 100.0,
        rolled_frac * 100.0,
        periodic_frac * 100.0,
    );
    eprintln!(
        "land on map borders: unrolled {} cells, rolled {} cells",
        land_on_borders(&unrolled),
        land_on_borders(&rolled),
    );

    // ACCEPTANCE 1 — the rolled export is ONE contiguous piece: its largest
    // non-periodic map component equals the whole torus continent. The unrolled
    // export is split at the seam (largest piece strictly smaller), and rolling
    // moves all land off the map border.
    assert!(
        rolled_frac >= periodic_frac - 0.005,
        "rolled largest map component {rolled_frac} must equal the torus continent {periodic_frac} (contiguous)"
    );
    assert!(
        unrolled_frac < periodic_frac - 0.005,
        "sanity: the unrolled export must be split at the seam (largest {unrolled_frac} < continent {periodic_frac})"
    );
    assert!(
        land_on_borders(&rolled) < land_on_borders(&unrolled),
        "the roll must move land off the map border"
    );

    // ACCEPTANCE 2 — the coastline's main ring is a CLOSED loop (a split continent
    // yields open curves terminating at opposite borders).
    let (len, closed) = longest_ring_closed(&coastline_geojson(&rolled));
    eprintln!("coastline: longest ring {len} pts, closed = {closed}");
    assert!(closed, "the main coastline ring must be a closed loop on the rolled export");

    // ACCEPTANCE 3 — exported land topology matches the preview (roll-invariant
    // largest area + traverse; FBM shifts them only slightly).
    let rolled_topo = land_topology(&rolled, SEA);
    let area_ratio = rolled_topo.largest_area_km2 / topo.largest_area_km2;
    let trav_ratio =
        rolled_topo.bbox_km.0.max(rolled_topo.bbox_km.1) / topo.bbox_km.0.max(topo.bbox_km.1);
    eprintln!("exported vs preview: area ×{area_ratio:.2}, traverse ×{trav_ratio:.2}");
    assert!((0.8..=1.25).contains(&area_ratio), "largest area must match the preview (±25%)");
    assert!((0.85..=1.15).contains(&trav_ratio), "traverse must match the preview (±15%)");

    // REPORT (do not fix) — water_class histogram on the rolled export.
    let wc = water_class(&rolled, SEA);
    let (mut c0, mut c1, mut c2) = (0usize, 0usize, 0usize);
    for &v in &wc {
        match v {
            0 => c0 += 1,
            1 => c1 += 1,
            2 => c2 += 1,
            _ => {}
        }
    }
    eprintln!(
        "water_class histogram (rolled): 0(land)={c0}  1(ocean)={c1}  2(inland water)={c2}  \
         total={}",
        wc.len()
    );
}
