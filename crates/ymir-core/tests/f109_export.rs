//! ADR 0001 Finding 109 — **the `.ymir` exports the author opens in LL**, seed 1, 8192².
//!
//! The author judges in the viz and then in LL — the Finding 82 path — so the round's deliverable
//! is a file, not a panel. Two of them:
//!
//! * `SHIPPED` — the checkbox OFF, the delivered world;
//! * `AGE_F109_s0.024` — the checkbox ON: `SlopeFloor::Absolute { s_eq: 0.024 }` + A1.
//!
//! ⚠️ **Both, deliberately, though the round asked for one.** Every line of the looking
//! instruction is comparative — "Δ ≥ 1 km +2", "D_L > 5: 32 → 0 %", "R8 0.092 → 0.048" — and an
//! A/B in LL needs two files. The extra cost is one 8192² build.
//!
//! ⚠️ **What these files are NOT.** They are the BREACHED field, which is what `.ymir` has always
//! shipped, while every gate of Findings 95–108 reads the ERODED one. Finding 108 measured the
//! archetype at `Δ 0.02 m` eroded against `623.56 m` breached: the two stages disagree about the
//! deepest object on the continent. `f109_toggle` reports both counts beside each other.
//!
//! ⚠️ The seam is proved BEFORE this bench runs: `f109_toggle` asserts the production path is
//! bit-identical to the Findings 106/108 bench on three seeds. Without that, these files would be
//! a different world from the one the numbers describe.
//!
//! `geo_scale_ratio = 7.5` — the author's practice, as at Finding 68.
//!
//! Run: cargo test -p ymir-core --release --test f109_export -- --ignored --nocapture

use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::climate::{c1_biomes, c1_climate_placed};
use ymir_core::export::container::{ContinentMeta, ContinentWriter, Grid};
use ymir_core::export::{height, hydro, vector};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, DrainageClimate, c1_drainage_windowed};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::coast_metrics::{CoastShape, NECK_KM, coast_shape_thresholds};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::breach_monotone;
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;
const TARGET: usize = 8192;
const CELL_KM: f32 = DOMAIN_KM / TARGET as f32;
const GEO_RATIO: f32 = 7.5;
/// The middle of Finding 106's window, re-read under Finding 108's differential class.
const S_EQ: f32 = 0.024;

#[derive(Clone, Copy, PartialEq)]
enum Variant {
    /// The checkbox OFF: the delivered world. Exported too, because every line of the round's
    /// looking instruction is COMPARATIVE ("Delta >= 1 km +2", "32 -> 0 %") and an A/B needs two
    /// files in LL.
    Shipped,
    /// The checkbox ON: ADR Finding 109, `SlopeFloor::Absolute { s_eq: 0.024 }` + A1.
    Age,
}

fn label(v: Variant) -> &'static str {
    match v {
        Variant::Shipped => "SHIPPED",
        Variant::Age => "AGE_F109_s0.024",
    }
}

fn build(v: Variant) -> GridF32 {
    let ss = SteinSteinParams::default();
    let run_cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run_cfg, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let refr = false;
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: TARGET,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: if refr { 0.0 } else { 0.04 },
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: !refr,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: !refr,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    // ADR Finding 109 -- the ONLY difference between the two exports. It is a CONFIG value, not
    // a knob edit, so this export and the viz reach the identical code path in `upscale_from_c1`,
    // and `f109_toggle` has already proved that path bit-identical to the Findings 106/108 bench.
    if v == Variant::Age {
        cfg.slope_floor = Some(ymir_core::terrain::upscale::SlopeFloor::Absolute { s_eq: S_EQ });
    }
    upscale_from_c1_with_progress(
        &state,
        &run_cfg.iso_config,
        &ss,
        &seed,
        &cfg,
        &edifices,
        &volc,
        Some(&kin),
        &mut |_| {},
        &|| false,
    )
    .0
    .heightmap
}

fn row(name: &str, s: &CoastShape) {
    eprintln!(
        "   {name:<26} {:>8} {:>10.0} {:>9.2} {:>9.2} {:>9.2} {:>8.2} {:>8.3}",
        s.count, s.coast_km, s.med_len_km, s.p90_len_km, s.max_len_km, s.share_pct, s.local_axis_r
    );
}

#[test]
#[ignore]
fn f109_export() {
    let ss = SteinSteinParams::default();
    let pp = PrecipParams::default();
    let out_root = std::path::Path::new("../../exports/f109_age");
    eprintln!("\n==========  Finding 75 · three coastal states, two length scales  ==========");
    eprintln!(
        "   cell = {:.4} km = {:.1} m. KILOMETRE scale: spur ≥ 1.00 km ({:.1} cells), neck \
         {NECK_KM} km. CELL scale: spur ≥ {:.4} km (2.0 cells), neck {:.4} km (1.0 cell).",
        CELL_KM,
        CELL_KM * 1000.0,
        1.0 / CELL_KM,
        2.0 * CELL_KM,
        CELL_KM
    );

    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;

    for v in [Variant::Shipped, Variant::Age] {
        let name = label(v);
        eprintln!("\n╔═══════════════ {name} ═══════════════╗");
        let raw = build(v);
        let (w, h) = (raw.width, raw.height);
        let pre0 = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
        // The export draws its coastline from the CONDITIONED field, so the coast is measured
        // on that one; the eroded-field row is printed beside it because Finding 56's published
        // figures were taken there and the two must stay comparable.
        let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);

        eprintln!(
            "   {:<26} {:>8} {:>10} {:>9} {:>9} {:>9} {:>8} {:>8}",
            "detector / field", "count", "coast km", "med km", "p90 km", "max km", "share%", "R"
        );
        for (fname, f) in [("eroded", &raw), ("breached (EXPORTED)", &field)] {
            let polys = marching_squares(f, SEA);
            let km_scale = coast_shape_thresholds(&polys, CELL_KM, 1.0, NECK_KM);
            let cell_scale = coast_shape_thresholds(&polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
            row(&format!("≥1 km · {fname}"), &km_scale);
            row(&format!("≥2 cells · {fname}"), &cell_scale);
        }

        // ── the export ────────────────────────────────────────────────────────
        let climate = c1_climate_placed(&field, &ss, 45.0, 40.0, &pp, DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
        let drainage = assemble_hd_drainage(
            &field,
            &dclim,
            Some(pre),
            &dcfg,
            &ss,
            DOMAIN_KM,
            GEO_RATIO,
            None,
            false,
        )
        .drainage;
        let biomes = c1_biomes(&field, &climate);
        let hl = height::metric_height_u16(&field, &ss);
        let meta = ContinentMeta {
            name: format!("seed{PSEED}_{TARGET}_{name}"),
            seed: PSEED,
            grid: Grid { width: w, height: h },
            window_km: DOMAIN_KM as f64,
            tectonic_domain_km: DOMAIN_KM as f64,
            window_offset_in_torus: [0.0, 0.578_125],
            latitude_deg: 45.0,
            latitude_span_deg: 40.0,
            geographic_scale_ratio: GEO_RATIO as f64,
            stein_stein: ss,
            sea_level_m: 0.0,
            max_elevation_m: hl.max_m as f64,
            max_depth_m: hl.min_m as f64,
        };
        let dir = out_root.join(format!("{}.ymir", meta.name));
        let mut wr = ContinentWriter::new(&dir, meta).unwrap();
        wr.add_raster_u16("height", &hl.codes).unwrap();
        wr.set_metric_range("height", hl.min_m as f64, hl.max_m as f64).unwrap();
        let coastline = vector::coastline_geojson(&field, None);
        wr.add_vector_file("coastline", "coastline.geojson", &coastline).unwrap();
        wr.set_level_m("coastline", 0.0).unwrap();
        let cell_m = DOMAIN_KM / w as f32 * 1000.0;
        let thr = vector::DEFAULT_CLIFF_THRESHOLD_DEG;
        let cliffs = vector::cliffs_geojson(&field, &ss, cell_m, thr);
        wr.add_vector_file("cliffs", "cliffs.geojson", &cliffs).unwrap();
        wr.set_slope_threshold_deg("cliffs", thr as f64).unwrap();
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
        wr.add_raster_i16("temperature", &temperature).unwrap();
        wr.add_raster_u16("precipitation", &precipitation).unwrap();
        wr.add_raster_u8("biome", &biome).unwrap();
        wr.add_raster_u32("lake_mask", &drainage.lake_map).unwrap();
        wr.add_raster_f32("flow_accumulation", &drainage.flow.accumulation.data).unwrap();
        wr.add_vector_file(
            "rivers",
            "rivers.json",
            &hydro::rivers_json(&drainage, CELL_KM * CELL_KM),
        )
        .unwrap();
        wr.add_vector_file("lakes", "lakes.json", &hydro::lakes_json(&drainage)).unwrap();
        let water = connectivity::water_class(&field, vector::SEA_LEVEL_NORM);
        wr.add_raster_u8("water_class", &water).unwrap();
        let manifest = wr.finish().unwrap();

        // ── the mouth question, from the EXPORT's own numbers ─────────────────
        let mut mouths: Vec<f32> = drainage
            .rivers
            .segments
            .iter()
            .enumerate()
            .filter(|(i, s)| {
                s.downstream.is_none()
                    && drainage.segment_kind[*i]
                        == ymir_core::tectonics_c1::drainage::SegmentKind::Watercourse
                    && !drainage.segment_profile_m[*i].is_empty()
            })
            .map(|(i, _)| *drainage.segment_profile_m[i].last().unwrap())
            .collect();
        mouths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let q = |p: f64| -> f32 {
            if mouths.is_empty() {
                f32::NAN
            } else {
                mouths[((mouths.len() - 1) as f64 * p) as usize]
            }
        };
        // the true SEA mouths: the terminus drains into an ocean cell
        let sea_mouths: Vec<f32> = drainage
            .rivers
            .segments
            .iter()
            .enumerate()
            .filter(|(i, s)| {
                s.downstream.is_none()
                    && drainage.segment_kind[*i]
                        == ymir_core::tectonics_c1::drainage::SegmentKind::Watercourse
            })
            .filter(|(_, s)| {
                let &(x, y) = s.points.last().unwrap();
                water[y as usize * w + x as usize] == 1 || {
                    let k = y as usize * w + x as usize;
                    let d = drainage.flow.direction[k];
                    d != ymir_core::terrain::flow::DIR_NONE && {
                        let nx = (x as i32 + ymir_core::terrain::flow::D8_DX[d as usize])
                            .rem_euclid(w as i32);
                        let ny = (y as i32 + ymir_core::terrain::flow::D8_DY[d as usize])
                            .rem_euclid(h as i32);
                        water[ny as usize * w + nx as usize] == 1
                    }
                }
            })
            .map(|(i, _)| *drainage.segment_profile_m[i].last().unwrap_or(&f32::NAN))
            .filter(|v| v.is_finite())
            .collect();
        let mut sm = sea_mouths.clone();
        sm.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!(
            "   MOUTHS — `sea_level_m` exported = 0.0 (hard-coded, hd.rs:1227)\n     all \
             {} Watercourse termini: p10 {:.1} MEDIAN {:.1} p75 {:.1} p90 {:.1} p95 {:.1} max \
             {:.1} m\n     the {} that actually drain into an OCEAN cell: min {:.1} MEDIAN {:.1} \
             p90 {:.1} max {:.1} m",
            mouths.len(),
            q(0.10),
            q(0.50),
            q(0.75),
            q(0.90),
            q(0.95),
            mouths.last().copied().unwrap_or(f32::NAN),
            sm.len(),
            sm.first().copied().unwrap_or(f32::NAN),
            sm.get(sm.len() / 2).copied().unwrap_or(f32::NAN),
            sm.get((sm.len().max(1) - 1) * 9 / 10).copied().unwrap_or(f32::NAN),
            sm.last().copied().unwrap_or(f32::NAN)
        );
        eprintln!(
            "   LAKES {} | max Watercourse {:.3} m³/s (ratio {GEO_RATIO}) | export → {}",
            drainage.lakes.len(),
            (0..drainage.rivers.segments.len())
                .filter(|&i| drainage.segment_kind[i]
                    == ymir_core::tectonics_c1::drainage::SegmentKind::Watercourse)
                .map(|i| drainage.segment_discharge_m3s[i])
                .fold(0.0f32, f32::max),
            manifest.display()
        );
    }
}
