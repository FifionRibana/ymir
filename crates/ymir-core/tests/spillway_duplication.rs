//! STEP 3a — diagnose the spillway DUPLICATION before touching the typing. If the pairs come
//! from double emission, typing first would create two `Spillway` objects where one belongs
//! and would mask the cause.
//!
//! The pairs are NOT all alike: 88 468 vs 88 449 (near-equal → probably two paths traced from
//! NEIGHBOURING start cells on the same basin, i.e. a col FLAT over several cells) but
//! 62 094 = 62 094 and 8 846 = 8 846 (exactly equal → true duplicates). Both cases may
//! coexist and they call for DIFFERENT remedies.
//!
//! Reported: spillways per basin (is it 1:1?), the duplicated groups with their start cells
//! and elevations, and whether the col is flat across them.
//!
//! Run: cargo test -p ymir-core --test spillway_duplication -- --ignored --nocapture

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, below_sea_basin_lakes, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const SPAN_DEG: f32 = 10.0;
const FULL_RANGE_M: f32 = 11302.0;

fn run(target: usize, lat: f32, label: &str) {
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
    // THE shipped config — this is exactly what `production_hd_config` is for.
    let cfg = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.09375, 0.578125],
        sample_size: 1.0,
        amplitude_base: 0.04,
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: true,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: true,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    let (up, _c) = upscale_from_c1_with_progress(
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
    );
    let raw = up.heightmap;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field = ymir_core::terrain::flow::breach_monotone(
        &raw,
        &pre.flow.filled,
        &pre.lake_map,
        0.5,
        raw.width,
        raw.height,
    );
    let (w, _h) = (field.width, field.height);
    let climate =
        c1_climate_placed(&field, &ss, lat, SPAN_DEG, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let bs = below_sea_basin_lakes(&field, &dclim, &dcfg, &ss, DOMAIN_KM, Some(&pre.lake_map));

    eprintln!(
        "\n================  SPILLWAY DUPLICATION — {label} {lat}° @ {target}²  ================"
    );
    eprintln!("below-sea LAKES: {} | SPILLWAYS: {}", bs.lakes.len(), bs.spillways.len());

    // 1:1 with basins? Group by lake_id.
    let mut per_lake: std::collections::BTreeMap<u32, Vec<usize>> = Default::default();
    for (i, sw) in bs.spillways.iter().enumerate() {
        per_lake.entry(sw.lake_id).or_default().push(i);
    }
    let multi: Vec<(&u32, &Vec<usize>)> = per_lake.iter().filter(|(_, v)| v.len() > 1).collect();
    eprintln!(
        "distinct lake_ids emitting a spillway: {} | ids emitting MORE THAN ONE: {}",
        per_lake.len(),
        multi.len()
    );
    if multi.is_empty() {
        eprintln!("  ⇒ emission is 1 spillway per basin. The duplication is NOT here.");
    }

    // ── THE DECISIVE TEST: RECIPROCAL CHAINING. A → B and B → A is a 2-cycle: both basins
    //    trace a spillway across the SAME col into each other. Physically only one can spill
    //    (or they are one water body sharing a col and should have merged).
    let by_id: std::collections::BTreeMap<u32, &_> =
        bs.spillways.iter().map(|s| (s.lake_id, s)).collect();
    let mut recip: Vec<(u32, u32)> = Vec::new();
    for sw in &bs.spillways {
        if let Some(t) = sw.chained_into {
            if let Some(other) = by_id.get(&t) {
                if other.chained_into == Some(sw.lake_id) && sw.lake_id < t {
                    recip.push((sw.lake_id, t));
                }
            }
        }
    }
    eprintln!(
        "\nRECIPROCAL CHAIN PAIRS (A→B and B→A): {} pair(s) = {} spillways",
        recip.len(),
        recip.len() * 2
    );
    for (a, b) in recip.iter().take(10) {
        let (sa, sb) = (by_id[a], by_id[b]);
        let same_start = sa.points.first() == sb.points.first();
        eprintln!(
            "  {a} ↔ {b} | drainage {:.0} / {:.0} km² | start {:?} vs {:?} | SAME start cell: {same_start}",
            sa.drainage_km2,
            sb.drainage_km2,
            sa.points.first(),
            sb.points.first()
        );
    }
    eprintln!(
        "  basins emitting a spillway but NOT inventoried as a lake: {} (spillways {} − lakes {})",
        bs.spillways.len() as i64 - bs.lakes.len() as i64,
        bs.spillways.len(),
        bs.lakes.len()
    );

    // ── DANGLING-ID CHECK (decision 2): the consumer must never have to resolve a lake_id
    //    absent from lakes.json — neither a spillway's OWN id nor its chain target.
    let alive: std::collections::HashSet<u32> = bs.lakes.iter().map(|l| l.base.id).collect();
    let own_dangling: Vec<u32> =
        bs.spillways.iter().filter(|s| !alive.contains(&s.lake_id)).map(|s| s.lake_id).collect();
    let tgt_dangling: Vec<u32> =
        bs.spillways.iter().filter_map(|s| s.chained_into).filter(|t| !alive.contains(t)).collect();
    eprintln!(
        "\nDANGLING IDS — spillway.lake_id absent from lakes: {} {:?}\n              spillway.chained_into absent from lakes: {} {:?}",
        own_dangling.len(),
        &own_dangling.iter().take(8).collect::<Vec<_>>(),
        tgt_dangling.len(),
        &tgt_dangling.iter().take(8).collect::<Vec<_>>()
    );

    // Group by (rounded) drainage to expose the near-equal pairs the author saw.
    let mut by_q: std::collections::BTreeMap<i64, Vec<usize>> = Default::default();
    for (i, sw) in bs.spillways.iter().enumerate() {
        by_q.entry((sw.drainage_km2 / 50.0).round() as i64).or_default().push(i);
    }
    let dup: Vec<(&i64, &Vec<usize>)> = by_q.iter().filter(|(_, v)| v.len() > 1).collect();
    eprintln!("\nspillway GROUPS with near-equal drainage (±50 km²): {}", dup.len());
    eprintln!(
        "{:>6} | {:>8} | {:>11} | {:>7} | {:>13} | {:>9} | {:>8}",
        "group", "lake_id", "drain km²", "pts", "start cell", "z_start m", "chained"
    );
    for (_k, idxs) in dup.iter().take(10) {
        for &i in idxs.iter() {
            let sw = &bs.spillways[i];
            let p0 = sw.points.first().copied().unwrap_or((0, 0));
            let z0 = sw.profile_m.first().copied().unwrap_or(0.0);
            eprintln!(
                "{:>6} | {:>8} | {:>11.0} | {:>7} | {:>6},{:<6} | {z0:>9.1} | {:>8}",
                i,
                sw.lake_id,
                sw.drainage_km2,
                sw.points.len(),
                p0.0,
                p0.1,
                sw.chained_into.map(|v| v.to_string()).unwrap_or_else(|| "sea".into())
            );
        }
        // Is the col FLAT across the start cells of this group?
        let zs: Vec<f32> = idxs
            .iter()
            .map(|&i| {
                let p = bs.spillways[i].points[0];
                field.data[p.1 as usize * w + p.0 as usize] * FULL_RANGE_M
            })
            .collect();
        let (mn, mx) = (
            zs.iter().cloned().fold(f32::MAX, f32::min),
            zs.iter().cloned().fold(f32::MIN, f32::max),
        );
        // pairwise cell distance between the group's start cells
        let mut maxd = 0.0f32;
        for a in 0..idxs.len() {
            for b in a + 1..idxs.len() {
                let (pa, pb) = (bs.spillways[idxs[a]].points[0], bs.spillways[idxs[b]].points[0]);
                let d = (((pa.0 as f32 - pb.0 as f32).powi(2)
                    + (pa.1 as f32 - pb.1 as f32).powi(2))
                .sqrt())
                .max(0.0);
                maxd = maxd.max(d);
            }
        }
        eprintln!(
            "        └ start-cell elevation spread {:.2} m | max start-cell separation {maxd:.1} cells  ⇒ {}",
            mx - mn,
            if mx - mn < 0.5 && maxd > 1.5 {
                "FLAT COL, neighbouring starts"
            } else if maxd <= 1.5 {
                "same/adjacent cell — true duplicate"
            } else {
                "distinct cols"
            }
        );
    }
    eprintln!(
        "\n(1 spillway per basin + near-equal drainages ⇒ the duplicates are NOT emitted twice by the\n basin loop; look at whether several BASINS share a col, or whether the same physical outflow is\n inventoried as several below-sea components.)"
    );
}

#[test]
#[ignore]
fn spillway_duplication_arid_hot() {
    run(8192, 25.0, "arid-hot");
}
