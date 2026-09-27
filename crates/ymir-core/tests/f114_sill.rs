//! ADR 0001 Finding 114 — **H-2: why the lakes do not drain. It is not a missing mechanism.**
//!
//! The author: lakes are transient, ageing should EMPTY them through their outlet, and the age
//! closure cannot do it because it raises beds without touching sills (Finding 113). The round
//! asks whether threshold incision is a one-pass closure or a small-N loop.
//!
//! ## ⛔ P0 — "lower the sill to its equilibrium" has no interior solution, and the dossier says
//! so three times
//!
//! `E = K·A^m·S^n` is a rate and vanishes only at `S = 0`. There is **no opposing term**:
//! * **L6522** — *"**no uplift term**, so the only fixed point of `h ← (h + f·h_r)/(1+f)` is
//!   `h = h_r` everywhere"*;
//! * **L236** — *"Stream power ran on a static field with no uplift, so channels graded down to
//!   sea level. Since Ymir's tectonics already did the uplift, **the fix is to LIMIT total
//!   incision, not add U: iters 3→2 and K 3000→1500**"* — the design decision block B would undo;
//! * **L1727 (Finding 39)** — *"In nature an overflow outlet incises and lowers the lake; the
//!   model leaves the outlet at the fill elevation… **a future spillway-incision pass (lower the
//!   sill, partially drain the lake) is the realistic fix if undesired. Left to the author's call
//!   — not folded in here.**"* and **L10550** — *"That is H-2 (threshold incision / lake levels)…
//!   **named here rather than attempted**"*. A rule-11b sweep in the dossier even records
//!   `spillway incision` **0**.
//!
//! ⇒ **Proposed twice, never written**, and "equilibrium" can only mean the receiver's elevation
//! (a monotone trench — what `breach_monotone` already computes for every depression NOT held as
//! a lake) or base level (drain everything).
//!
//! ## ⛔ P1 — and the arithmetic says one pass is already a full relaxation
//!
//! For `n = 1` the shipped update is `h ← (h_o + f·h_r)/(1 + f)` with `f = k·dt·A_km²^m / dist_m`.
//! **One pass therefore removes the fraction `f/(1+f)` of the cell's height above its receiver.**
//! At `k·dt = 4500` and `dist ≈ 48.8 m`, `f ≈ 92·√A`: a 20 km² channel keeps **0.24 %** of its
//! drop. The incision is not weak — Finding 90 measured Courant **3 699**.
//!
//! ⇒ **The lakes that survive do not survive because the incision is too weak. They survive
//! because their sill is never incised at all.** Two candidate reasons, and this bench separates
//! them: the sill sits below the channel head (`A < A_c`, hillslope regime, `continue`), or
//! `depression_floor` (A1) excludes it as a depression cell. That is the measurement H-2 needs
//! before anyone writes a mechanism.
//!
//! Run: cargo test -p ymir-core --release --test f114_sill -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::{RELIEF_V1_A_C_KM2, StreamPowerConfig};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, DrainageClimate, c1_drainage_windowed};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;

#[test]
#[ignore]
fn f114_sill() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 114 . the sills, and why they are never incised  =========");
    eprintln!(
        "   shipped law, read from the config (not retyped): k **{}** · dt **{}** · m **{}** · n \
         **{}** · iterations **{}** · A_c **{:.3} km²** = {:.1} cells · k·dt **{}**",
        sp.k,
        sp.dt,
        sp.m,
        sp.n,
        sp.iterations,
        sp.min_area_cells * cell_km2,
        sp.min_area_cells,
        sp.k * sp.dt
    );

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

    // ── A · the state of the sills, before any mechanism ────────────────────
    eprintln!(
        "\n   ── A · every lake ≥ 1 km², its sill, and whether the incision can SEE it ──\n      \
         {:>9} {:>9} {:>9} {:>9} {:>10} {:>9} {:>8} {:>7}  {}",
        "area km²",
        "sill m",
        "floor m",
        "depth m",
        "A@sill km²",
        "drop ERO",
        "drop BRE",
        "N",
        "why"
    );
    let (mut deep, mut shallow, mut blind, mut one_pass) = (0usize, 0usize, 0usize, 0usize);
    let mut depths = Vec::new();
    let mut rows = Vec::new();
    // ⚠️ RULE 12 — the outlet index comes from the inventory built on the BREACHED field, so
    // reading heights on the ERODED one mixes two stages, and the whole "the sill drop is 0"
    // conclusion rests on that number. Both are therefore read, each on its OWN field with its
    // OWN flow direction, and both are printed.
    let brd = c1_drainage_windowed(&bf, None, &dcfg, &ss, DOMAIN_KM);
    for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
        let (ox, oy) = (l.base.outlet.0 as usize, l.base.outlet.1 as usize);
        let k = oy * w + ox;
        if k >= w * h {
            continue;
        }
        // The accumulation the INCISION would see at the sill, in km², and the drop to its
        // receiver — the two quantities the `n = 1` update is made of.
        let a_km2 = pre.flow.accumulation.data[k] * cell_km2;
        let d = pre.flow.direction[k];
        let (drop_m, dist_m) = if d == DIR_NONE {
            (0.0, CELL_M)
        } else {
            let nx = (ox as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
            let ny = (oy as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
            let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
            (
                ((f.data[k] - f.data[ny * w + nx]).max(0.0)) * n2m,
                if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M },
            )
        };
        // the same quantity on the BREACHED field, its own direction, its own heights
        let db = brd.flow.direction[k];
        let drop_b = if db == DIR_NONE {
            0.0
        } else {
            let nx = (ox as i32 + D8_DX[db as usize]).rem_euclid(w as i32) as usize;
            let ny = (oy as i32 + D8_DY[db as usize]).rem_euclid(h as i32) as usize;
            ((bf.data[k] - bf.data[ny * w + nx]).max(0.0)) * n2m
        };
        let fr = sp.k * sp.dt * a_km2.powf(sp.m) / dist_m; // the `f` of the shipped update
        let per_pass = drop_m * fr / (1.0 + fr);
        let depth = l.depth_m;
        depths.push(depth);
        let seen = a_km2 >= RELIEF_V1_A_C_KM2;
        if !seen {
            blind += 1;
        }
        if depth < 30.0 {
            shallow += 1;
        }
        if depth > 100.0 {
            deep += 1;
        }
        let n_pass = if per_pass <= 0.0 { f32::INFINITY } else { depth / per_pass };
        if seen && n_pass <= 1.0 {
            one_pass += 1;
        }
        rows.push((l.area_km2, l.level_m, depth, a_km2, drop_m, drop_b, n_pass, seen));
    }
    rows.sort_by(|a, b| b.0.total_cmp(&a.0));
    for r in rows.iter() {
        eprintln!(
            "      {:>9.2} {:>9.1} {:>9.1} {:>9.1} {:>10.4} {:>10.3} {:>10.3} {:>7}  {}",
            r.0,
            r.1,
            r.1 - r.2,
            r.2,
            r.3,
            r.4,
            r.5,
            if r.6.is_finite() { format!("{:.1}", r.6) } else { "inf".into() },
            if !r.7 {
                "**A < A_c — the incision NEVER touches this sill**"
            } else if r.4 < 0.5 {
                "sill already level with its receiver"
            } else {
                "incisable"
            }
        );
    }
    let sd = sorted(depths.clone());
    eprintln!(
        "\n   ⇒ **{} lakes ≥ 1 km²** · depth p10/p50/p90 **{:.1} / {:.1} / {:.1} m** · under 30 m \
         **{shallow}** · over 100 m **{deep}**",
        rows.len(),
        pct(&sd, 0.10),
        pct(&sd, 0.50),
        pct(&sd, 0.90)
    );
    eprintln!(
        "   ⛔ **sills the incision CANNOT see (A < A_c = {RELIEF_V1_A_C_KM2} km²): {blind} of \
         {}** = **{:.0} %** · sills that would empty their lake in ONE pass if they were \
         incisable: **{one_pass}**",
        rows.len(),
        100.0 * blind as f32 / rows.len().max(1) as f32
    );
    // ── B · WHERE the drop is, if not at the sill ───────────────────────────
    // "S = 0 at the sill" is only half a finding: a LOCAL law sees nothing there, but a
    // knickpoint-retreat mechanism would see a fall further down. Walking the flow path turns the
    // hypothesis into a number — if the drop is 0 at one cell and large at fifty, the head of the
    // fall is downstream and no local operator can reach it.
    eprintln!(
        "\n   ── B · the profile BELOW each sill: where the drop actually is ──\n      {:>9} \
         {:>8} {:>8} {:>8} {:>8} {:>10}",
        "area km²", "d@1", "d@10", "d@50", "d@200", "max S"
    );
    let mut far = Vec::new();
    for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
        let (ox, oy) = (l.base.outlet.0 as usize, l.base.outlet.1 as usize);
        let mut k = oy * w + ox;
        if k >= w * h {
            continue;
        }
        let z0 = bf.data[k] * n2m;
        let (mut d1, mut d10, mut d50, mut d200, mut smax) = (0f32, 0f32, 0f32, 0f32, 0f32);
        let mut prev = z0;
        for step in 1..=200usize {
            let d = brd.flow.direction[k];
            if d == DIR_NONE {
                break;
            }
            let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
            let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
            let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
            let dist = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
            k = ny * w + nx;
            let z = bf.data[k] * n2m;
            smax = smax.max((prev - z).max(0.0) / dist);
            prev = z;
            let drop = z0 - z;
            match step {
                1 => d1 = drop,
                10 => d10 = drop,
                50 => d50 = drop,
                200 => d200 = drop,
                _ => {}
            }
        }
        far.push(d200);
        eprintln!(
            "      {:>9.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>10.4}",
            l.area_km2, d1, d10, d50, d200, smax
        );
    }
    let sf = sorted(far.clone());
    eprintln!(
        "\n   ⇒ drop 200 cells (9.8 km) below the sill: p10/p50/p90 **{:.1} / {:.1} / {:.1} m** ⇒ \
         **{}**",
        pct(&sf, 0.10),
        pct(&sf, 0.50),
        pct(&sf, 0.90),
        if pct(&sf, 0.50) > 20.0 {
            "the fall IS downstream — a LOCAL law cannot reach it; this is knickpoint retreat"
        } else {
            "**flat downstream too — there is no fall to retreat**"
        }
    );
    eprintln!("\n==========  end Finding 114 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
