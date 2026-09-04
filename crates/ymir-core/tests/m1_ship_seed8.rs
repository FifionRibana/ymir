//! M1 #190 end-to-end validation — ship the border-clean island (seed 8) at
//! full HD 8192² and run the WHOLE pipeline on it (tectonics → windowed HD
//! upscale + 64M-droplet erosion → climate → drainage → biomes → the full v1
//! `.ymir` container), then report wall-clock per stage and PEAK working set.
//!
//! Config (from the border-clean sweep, the smallest clean non-wrapping island):
//!   num_plates=20, seed_cluster_count=3, seed=8, target_land_fraction=0.10,
//!   window_km=418 (traverse 368 + 2·25 margin) → 51.0 m/cell at 8192².
//!
//! Heavy + offline: `#[ignore]`. Run:
//!   cargo test -p ymir-core --test m1_ship_seed8 -- --ignored --nocapture

use std::time::Instant;

use ymir_core::climate::biomes::Biome;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::climate::{c1_biomes, c1_climate_windowed};
use ymir_core::export::container::{ContinentMeta, ContinentWriter, Grid};
use ymir_core::export::{height, hydro, vector};
use ymir_core::lakes::connectivity;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::cached_product::{
    cached_c1_drainage_windowed, cached_c1_eroded, coarse_normalized_sweep, eroded_key,
    tectonic_key,
};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::C1DrainageConfig;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::land_topology::land_topology;
use ymir_core::tectonics_c1::production_upscale::C1_DOMAIN_KM;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

/// Peak working-set size (MB) of the current process (Windows: kernel32
/// K32GetProcessMemoryInfo; other platforms return 0).
#[cfg(windows)]
fn peak_ws_mb() -> u64 {
    #[repr(C)]
    #[derive(Default)]
    struct Pmc {
        cb: u32,
        page_fault_count: u32,
        peak_ws: usize,
        ws: usize,
        qppp: usize,
        qpp: usize,
        qpnpp: usize,
        qnpp: usize,
        pagefile: usize,
        peak_pagefile: usize,
    }
    unsafe extern "system" {
        fn GetCurrentProcess() -> isize;
        fn K32GetProcessMemoryInfo(process: isize, counters: *mut Pmc, cb: u32) -> i32;
    }
    let mut pmc = Pmc { cb: std::mem::size_of::<Pmc>() as u32, ..Default::default() };
    unsafe {
        K32GetProcessMemoryInfo(GetCurrentProcess(), &mut pmc, pmc.cb);
    }
    (pmc.peak_ws / (1024 * 1024)) as u64
}
#[cfg(not(windows))]
fn peak_ws_mb() -> u64 {
    0
}

#[test]
#[ignore]
fn ship_seed8_full_pipeline() {
    let grid = 64usize;
    let seed = 8u64;
    let target_size = 8192usize;
    let tlf = 0.10f32;
    let window_km = 418.0f32;
    let lat = 45.0f32;

    let mut init = Phase2InitParams::default();
    init.num_plates = 20;
    init.cluster.seed_cluster_count = 3;
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
    let pp = PrecipParams::default();

    // Window centred on the island's bbox centre (matching evaluate_island).
    let field = coarse_normalized_sweep(seed, grid, &init, &run, &clo, &ss, &[tlf]).remove(0).1;
    let topo = land_topology(&field, 0.5);
    let cx = (topo.bbox_min.0 + topo.bbox_max.0) as f64 / 2.0 / grid as f64;
    let cy = (topo.bbox_min.1 + topo.bbox_max.1) as f64 / 2.0 / grid as f64;
    let wf = (window_km / C1_DOMAIN_KM) as f64;
    let origin = [(cx - wf / 2.0).clamp(0.0, 1.0 - wf), (cy - wf / 2.0).clamp(0.0, 1.0 - wf)];
    eprintln!(
        "island: largest {:.0} km² ({:.0}%), {} masses, wrap x={} y={}; window {window_km} km \
         → {:.1} m/cell; origin [{:.3},{:.3}] size {:.3}",
        topo.largest_area_km2,
        topo.largest_area_frac * 100.0,
        topo.num_landmasses,
        topo.wraps_x,
        topo.wraps_y,
        window_km / target_size as f32 * 1000.0,
        origin[0],
        origin[1],
        wf
    );

    let mut upscale = FbmUpscaleConfig::c1_hd_production(target_size);
    upscale.target_land_fraction = Some(tlf);
    upscale.sample_origin = origin;
    upscale.sample_size = wf;

    let cache = std::env::temp_dir().join("ymir_ship_seed8_cache");
    let _ = std::fs::create_dir_all(&cache);

    // ── Eroded (tectonics → windowed HD upscale → 64M-droplet erosion). ──
    let t = Instant::now();
    let eroded = cached_c1_eroded(
        &cache,
        seed,
        grid,
        &init,
        &run,
        &clo,
        &ss,
        &upscale,
        &ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig::default(),
    )
    .unwrap()
    .heightmap;
    eprintln!(
        "[stage] eroded {}² in {:.0} s (peak {} MB)",
        eroded.width,
        t.elapsed().as_secs_f32(),
        peak_ws_mb()
    );

    // ── Climate (windowed). ──
    let t = Instant::now();
    let climate = c1_climate_windowed(&eroded, &ss, lat, &pp, window_km);
    eprintln!("[stage] climate in {:.0} s (peak {} MB)", t.elapsed().as_secs_f32(), peak_ws_mb());

    // ── Drainage (windowed, cached). ──
    let t = Instant::now();
    let ekey = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss, &upscale);
    let drainage = cached_c1_drainage_windowed(
        &cache,
        &ekey,
        &eroded,
        Some((lat, &pp)),
        &C1DrainageConfig::default(),
        &ss,
        window_km,
    )
    .unwrap();
    eprintln!("[stage] drainage in {:.0} s (peak {} MB)", t.elapsed().as_secs_f32(), peak_ws_mb());

    // ── Biomes. ──
    let biomes = c1_biomes(&eroded, &climate);

    // ── Full .ymir export (every layer). ──
    let t = Instant::now();
    let (w, h) = (eroded.width, eroded.height);
    let n = w * h;
    let height_layer = height::metric_height_u16(&eroded, &ss);
    let meta = ContinentMeta {
        name: format!("seed{seed}_island_{target_size}"),
        seed,
        grid: Grid { width: w, height: h },
        window_km: window_km as f64,
        tectonic_domain_km: C1_DOMAIN_KM as f64,
        window_offset_in_torus: origin,
        latitude_deg: lat as f64,
        latitude_span_deg: window_km as f64 / 111.0, // geographic-span identity (Finding 25)
        geographic_scale_ratio: 1.0,                 // identity (Finding 24)
        stein_stein: ss,
        sea_level_m: 0.0,
        max_elevation_m: height_layer.max_m as f64,
        max_depth_m: height_layer.min_m as f64,
    };
    let dir = std::env::temp_dir().join("ymir_ship_seed8_out").join(format!("{}.ymir", meta.name));
    let mut writer = ContinentWriter::new(&dir, meta).unwrap();
    writer.add_raster_u16("height", &height_layer.codes).unwrap();
    writer
        .set_metric_range("height", height_layer.min_m as f64, height_layer.max_m as f64)
        .unwrap();

    let coastline = vector::coastline_geojson(&eroded);
    writer.add_vector_file("coastline", "coastline.geojson", &coastline).unwrap();
    writer.set_level_m("coastline", 0.0).unwrap();
    let cell_size_m = window_km / w as f32 * 1000.0;
    let threshold = vector::DEFAULT_CLIFF_THRESHOLD_DEG;
    let cliffs = vector::cliffs_geojson(&eroded, &ss, cell_size_m, threshold);
    writer.add_vector_file("cliffs", "cliffs.geojson", &cliffs).unwrap();
    writer.set_slope_threshold_deg("cliffs", threshold as f64).unwrap();

    let temperature: Vec<i16> = climate
        .temperature
        .data
        .iter()
        .map(|&c| (c * 100.0).round().clamp(i16::MIN as f32, i16::MAX as f32) as i16)
        .collect();
    let precipitation: Vec<u16> = climate
        .precipitation
        .data
        .iter()
        .map(|&p| precip_mm_per_year(p).round().clamp(0.0, u16::MAX as f32) as u16)
        .collect();
    let biome: Vec<u8> = biomes.iter().map(|b| b.to_u8()).collect();
    writer.add_raster_i16("temperature", &temperature).unwrap();
    writer.add_raster_u16("precipitation", &precipitation).unwrap();
    writer.add_raster_u8("biome", &biome).unwrap();
    writer.add_raster_u32("lake_mask", &drainage.lake_map).unwrap();
    writer.add_raster_f32("flow_accumulation", &drainage.flow.accumulation.data).unwrap();
    writer
        .add_vector_file(
            "rivers",
            "rivers.json",
            &hydro::rivers_json(&drainage, 400.0 / 8192.0 * (400.0 / 8192.0)),
        )
        .unwrap();
    writer.add_vector_file("lakes", "lakes.json", &hydro::lakes_json(&drainage)).unwrap();
    let water = connectivity::water_class(&eroded, vector::SEA_LEVEL_NORM);
    writer.add_raster_u8("water_class", &water).unwrap();
    let manifest = writer.finish().unwrap();
    eprintln!("[stage] export in {:.0} s → {}", t.elapsed().as_secs_f32(), manifest.display());

    // ── Report: emerged fraction of the exported window + peak RSS. ──
    let sea = vector::SEA_LEVEL_NORM;
    let emerged = eroded.data.iter().filter(|&&v| v > sea).count() as f32 / n as f32;
    let ocean_border = {
        // sanity: the window edge should be ocean (border-clean at HD).
        let mut land_on_edge = 0usize;
        for x in 0..w {
            if eroded.data[x] > sea {
                land_on_edge += 1;
            }
            if eroded.data[(h - 1) * w + x] > sea {
                land_on_edge += 1;
            }
        }
        for y in 0..h {
            if eroded.data[y * w] > sea {
                land_on_edge += 1;
            }
            if eroded.data[y * w + (w - 1)] > sea {
                land_on_edge += 1;
            }
        }
        land_on_edge
    };
    let biome_ocean = biome.iter().filter(|&&b| b == Biome::Ocean.to_u8()).count();
    eprintln!(
        "RESULT: {w}²  emerged {:.1}%  ocean-biome {:.1}%  land cells on window edge {ocean_border}  \
         rivers {}  lakes {}  PEAK RSS {} MB",
        emerged * 100.0,
        biome_ocean as f32 / n as f32 * 100.0,
        drainage.rivers.segments.len(),
        drainage.lakes.len(),
        peak_ws_mb()
    );

    assert_eq!(eroded.width, target_size);
    assert!(emerged > 0.0 && emerged < 0.5, "a window framing one island is mostly ocean");
}
