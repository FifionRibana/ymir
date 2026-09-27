//! ADR 0001 Finding 117-A — **how far each family-1 lake COULD fall, and what holds its base.**
//!
//! Finding 114 split the 29 lakes into family 1 (a step below the sill, 16–20 lakes) and family
//! 2 (below-sea basins, flat to the ocean). Finding 115 showed the retreat wave reaches any step in
//! a fraction of a pass. This bench asks the remaining geometric question: for each family-1 lake,
//! follow its D8 chain to the **terminal base** — the sea, another lake, or a family-2 basin — and
//! report the **drainable**: how far the lake level sits above the lowest profile the base allows.
//!
//! ## ⛔ Two drainables, not one
//!
//! `Σ S_eq·dx` with `S_eq = c·A^(−1/2)` is the **age closure's** floor (c = 0.024, Finding 96/108).
//! In the delivered world that floor does not exist — Findings 61/62: *"with no uplift the only
//! fixed point is base level"*. So every lake gets two columns: **delivered** (`level − base`) and
//! **closed** (`level − base − Σ S_eq·dx` over channel cells `A ≥ A_c`). The closed one is always
//! the smaller.
//!
//! ## ⚠️ The oracle control is NOT runnable
//!
//! The round asks that the lakes Finding 95 drained (24.05 % → 15.92 % of land at 300 passes) be
//! the high-drainable ones. Finding 95 recorded the **area fraction only**; no per-lake inventory
//! survives and the oracle costs 1 h 43 to rebuild (refused at Finding 96). Declared, not faked.
//!
//! Run: cargo test -p ymir-core --release --test f117_drain -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, DrainageClimate, c1_drainage_windowed};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
/// Finding 83's shipped base-level bound: the incision may not cut below `sea + 0.5 m`.
const BASE_LEVEL_M: f32 = 0.5;
const MAX_WALK: usize = 20_000;

#[test]
#[ignore]
fn f117_drain() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 117-A . the drainable, and what holds the base  ==========");

    let f = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let (w, h) = (f.width, f.height);
    let pre = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
    let bf = breach_monotone(&f, &pre.flow.filled, &pre.lake_map, SEA, w, h);
    let climate = c1_climate_placed(&bf, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let pre2 = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
    let dr = assemble_hd_drainage(
        &bf,
        &dclim,
        Some(pre2),
        &dcfg,
        &ss,
        DOMAIN_KM,
        GEO_RATIO,
        None,
        false,
    )
    .drainage;
    // the walk is on the BREACHED field with its own flow: the stage the exported network lives on
    let brd = c1_drainage_windowed(&bf, None, &dcfg, &ss, DOMAIN_KM);
    let m = |k: usize| c1_altitude_norm_to_metres(bf.data[k], &ss);
    let lake_of = |k: usize| dr.lake_map[k];
    let level_of = |id: u32| dr.lakes.iter().find(|l| l.base.id == id).map(|l| l.level_m);
    let is_below_sea = |id: u32| id >= 1_000_001;

    eprintln!(
        "\n   {:>8} {:>8} {:>7} {:>10} {:>9} {:>9} {:>9} {:>9}  {}",
        "area km²",
        "level m",
        "floor m",
        "chain km",
        "base m",
        "ΣS_eq·dx",
        "drain DEL",
        "drain CLO",
        "terminal base"
    );
    let (mut n_sea, mut n_fam2, mut n_lake, mut n_lost) = (0usize, 0usize, 0usize, 0usize);
    let (mut d_del, mut d_clo) = (Vec::new(), Vec::new());
    let mut lakes: Vec<_> = dr.lakes.iter().filter(|l| l.area_km2 >= 1.0).collect();
    lakes.sort_by(|a, b| b.area_km2.total_cmp(&a.area_km2));
    for l in &lakes {
        // family 2 = below-sea basins: their base is the sea by construction (Finding 114)
        if is_below_sea(l.base.id) {
            continue;
        }
        let (ox, oy) = (l.base.outlet.0 as usize, l.base.outlet.1 as usize);
        let mut k = oy * w + ox;
        if k >= w * h {
            continue;
        }
        let (mut km, mut seq_m, mut base_m, mut base_name) = (0f32, 0f32, f32::NAN, String::new());
        for _ in 0..MAX_WALK {
            if bf.data[k] <= SEA {
                base_m = BASE_LEVEL_M;
                base_name = "SEA".into();
                break;
            }
            let id = lake_of(k);
            if id != 0 && id != l.base.id {
                base_m = level_of(id).unwrap_or(m(k));
                base_name = if is_below_sea(id) {
                    format!("FAMILY 2 basin {id}")
                } else {
                    format!("lake {id}")
                };
                break;
            }
            let d = brd.flow.direction[k];
            if d == DIR_NONE {
                break;
            }
            let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
            let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
            let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
            let dx = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
            let a_km2 = brd.flow.accumulation.data[k] * cell_km2;
            if a_km2 >= RELIEF_V1_A_C_KM2 {
                seq_m += S_EQ / a_km2.sqrt() * dx; // the closure's graded profile, channels only
            }
            km += dx / 1000.0;
            k = ny * w + nx;
        }
        if base_m.is_nan() {
            n_lost += 1;
            base_name = "**no base found**".into();
        } else if base_name == "SEA" {
            n_sea += 1;
        } else if base_name.starts_with("FAMILY") {
            n_fam2 += 1;
        } else {
            n_lake += 1;
        }
        let drain_del = l.level_m - base_m;
        let drain_clo = drain_del - seq_m;
        if base_m.is_finite() {
            d_del.push(drain_del);
            d_clo.push(drain_clo);
        }
        eprintln!(
            "   {:>8.2} {:>8.1} {:>7.1} {:>10.1} {:>9.1} {:>9.1} {:>9.1} {:>9.1}  {}",
            l.area_km2,
            l.level_m,
            l.level_m - l.depth_m,
            km,
            base_m,
            seq_m,
            drain_del,
            drain_clo,
            base_name
        );
    }
    let sd = sorted(d_del.clone());
    let sc = sorted(d_clo.clone());
    eprintln!(
        "\n   ⇒ family-1 lakes: **{}** · terminal base SEA **{n_sea}** · FAMILY-2 basin **{n_fam2}** \
         · another lake **{n_lake}** · none found {n_lost}",
        d_del.len() + n_lost
    );
    eprintln!(
        "   drainable, DELIVERED world (level − base): p10/p50/p90 **{:.0} / {:.0} / {:.0} m** · \
         > 100 m: **{}** · < 10 m: **{}**",
        pct(&sd, 0.10),
        pct(&sd, 0.50),
        pct(&sd, 0.90),
        d_del.iter().filter(|&&v| v > 100.0).count(),
        d_del.iter().filter(|&&v| v < 10.0).count()
    );
    eprintln!(
        "   drainable, CLOSED world (− Σ S_eq·dx):     p10/p50/p90 **{:.0} / {:.0} / {:.0} m** · \
         > 100 m: **{}** · < 10 m: **{}**",
        pct(&sc, 0.10),
        pct(&sc, 0.50),
        pct(&sc, 0.90),
        d_clo.iter().filter(|&&v| v > 100.0).count(),
        d_clo.iter().filter(|&&v| v < 10.0).count()
    );
    eprintln!(
        "   ⚠️ oracle control (Finding 95, 24.05 → 15.92 %): NOT RUNNABLE — the dossier holds the \
         area fraction only, no per-lake list; the oracle costs 1 h 43 to rebuild."
    );
    eprintln!(
        "\n==========  end Finding 117-A . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
