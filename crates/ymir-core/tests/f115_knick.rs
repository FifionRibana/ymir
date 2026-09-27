//! ADR 0001 Finding 115 — **knickpoint retreat is NOT missing: it is 3 699× over Courant and has
//! already run. What survives it is a GATE, not a wave.**
//!
//! The round opens "the last major fluvial process absent from the pipeline" and asks whether the
//! retreat is a one-pass closure or an N-step solver. ⛔ **Rule 11 answers before the first
//! measurement, and the number is the dossier's own.**
//!
//! ADR **L14093-14100** (Finding 90, from Finding 44's `timescale_plan()`), at the shipped config:
//!
//! | | |
//! |---|---|
//! | `k` / `dt` / `iterations` | 4 500 / 1 / 2 |
//! | **knickpoint celerity `c = K·A_km²^m`** | **180 618 m/yr** at `A_max` = 1 611 km² |
//! | `dt_max = cell_m / c` | 2.70·10⁻⁴ yr |
//! | **`cfl_iterations` / Courant** | **7 399 / 3 699** |
//!
//! > *"**At the shipped K a knickpoint crosses a 400 km continent in a little over two years**"* —
//! > and `T = iterations · dt` = **2 years**.
//!
//! ⇒ The shipped budget is **one continental crossing**. A step 2–10 km from a sill is reached in
//! **0.011–0.055 yr = 0.6 %–2.8 % of ONE pass**. Blocks B and C are therefore moot: there is no
//! retreat to price at one pass against twenty, because it is already effectively instantaneous.
//!
//! ⚠️ And the same arithmetic falls straight out of the update Finding 114 read off the code:
//! `f = k·dt·A_km²^m / dist_m` is dimensionless, so **`k·dt·√A` is the metres the wave travels per
//! pass** — 4 500·√A. At `A_max` = 1 611 km² that is 180 600 m, the dossier's number to four
//! digits. Two independent routes, same constant.
//!
//! ## What this bench measures instead
//!
//! If the wave crosses a continent in the shipped budget, **a step that survives is a place where
//! the law is switched OFF**. There are exactly four ways for that in `incise`:
//!
//! 1. `area < a_c_cells` → `continue` — the hillslope regime (`A_c` = 0.1 km²);
//! 2. `cfg.depression_floor && flow.filled > field` → `continue` — A1, Finding 96/104;
//! 3. the receiver is already at the cell's own height — nothing to transmit;
//! 4. `base_level_floor` — `h_r` floored at `sea + 0.5 m`, Finding 83.
//!
//! The bench walks each lake's outlet downstream and attributes **every cell of the path** to one
//! of the four. That is the question the round should be asking, and it needs no new mechanism.
//!
//! Run: cargo test -p ymir-core --release --test f115_knick -- --ignored --nocapture

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
const WALK: usize = 200;

#[test]
#[ignore]
fn f115_knick() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    let kdt = sp.k * sp.dt;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 115 . what survives a wave that crosses a continent  ======");
    eprintln!(
        "   celerity from the SHIPPED config: c·dt = k·dt·√A = **{kdt}·√A m per pass** ⇒ at A_c \
         ({RELIEF_V1_A_C_KM2} km²) **{:.0} m**, at 1 km² **{:.0} m**, at 20 km² **{:.0} m** · ADR \
         L14093: 180 618 m/yr at A_max 1 611 km², here **{:.0} m**",
        kdt * RELIEF_V1_A_C_KM2.sqrt(),
        kdt,
        kdt * 20f32.sqrt(),
        kdt * 1611f32.sqrt()
    );

    // ⚠️ METHOD RULE 18 — `depression_floor` is set to TRUE by the age closure's own branch
    // (Finding 109); it is FALSE by default in `relief_v3`. So an A1 count measured only with the
    // closure ON cannot tell "the closure is what stops the lakes draining" from "those cells are
    // in depressions anyway". The two readings are opposite. Both worlds, or neither.
    for (label, knobs) in [
        ("F109 OFF (delivered — depression_floor = false)", Knobs::passes(2)),
        (
            "F109 ON  (s_eq 0.024 — depression_floor = TRUE)",
            Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
        ),
    ] {
        eprintln!("\n########  {label}  ########");
        // A1 is ON only in the closure branch — `Knobs::slope_floor_abs` is what sets it.
        let a1_on = knobs.slope_floor_abs.is_some();
        let f = build_field_seed(knobs, PSEED);
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

        // ⚠️ The walk is on the ERODED field with the ERODED flow — the stage the INCISION acts on
        // (rule 12). Finding 114 showed the breached field gives the same zero at the sill.
        eprintln!(
            "\n   ── every lake's outlet, walked {WALK} cells down the ERODED flow ──\n      {:>9} \
         {:>8} {:>9} {:>9} {:>9} {:>8}  {}",
            "area km²",
            "step m",
            "at cell",
            "A@step",
            "T_col",
            "gate@step",
            "path: why the law is off"
        );
        let (mut fam1, mut fam2) = (0usize, 0usize);
        let mut tcols = Vec::new();
        for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
            let (ox, oy) = (l.base.outlet.0 as usize, l.base.outlet.1 as usize);
            let mut k = oy * w + ox;
            if k >= w * h {
                continue;
            }
            // walk, recording the biggest local drop and WHY each cell is or is not incisable
            let (mut best_s, mut best_at, mut best_drop, mut best_a) = (0f32, 0usize, 0f32, 0f32);
            let (mut n_ac, mut n_dep, mut n_flat, mut n_live, mut t_col) =
                (0usize, 0usize, 0usize, 0usize, 0f32);
            for step in 1..=WALK {
                let d = pre.flow.direction[k];
                if d == DIR_NONE {
                    break;
                }
                let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
                let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
                let dist = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
                let nk = ny * w + nx;
                let a_km2 = pre.flow.accumulation.data[k] * cell_km2;
                let drop = (f.data[k] - f.data[nk]).max(0.0) * n2m;
                let s = drop / dist;
                // the four gates, in the order `incise` applies them
                if a_km2 < RELIEF_V1_A_C_KM2 {
                    n_ac += 1;
                } else if a1_on && pre.flow.filled.data[k] > f.data[k] {
                    n_dep += 1;
                } else if drop <= 0.01 {
                    n_flat += 1;
                } else {
                    n_live += 1;
                }
                // the retreat time to bring THIS cell's step up to the sill, in passes
                t_col += dist / (kdt * a_km2.max(RELIEF_V1_A_C_KM2).sqrt());
                if s > best_s {
                    best_s = s;
                    best_at = step;
                    best_drop = drop;
                    best_a = a_km2;
                }
                k = nk;
            }
            let has_step = best_drop > 1.0;
            if has_step {
                fam1 += 1;
                tcols.push(t_col);
            } else {
                fam2 += 1;
            }
            let gate = if !has_step {
                "—".to_string()
            } else if best_a < RELIEF_V1_A_C_KM2 {
                "**A<A_c**".to_string()
            } else {
                "live".to_string()
            };
            eprintln!(
                "      {:>9.2} {:>8.2} {:>9} {:>9.3} {:>9.4} {:>8}  A<A_c {n_ac} · depression {n_dep} \
             · flat {n_flat} · **incisable {n_live}**",
                l.area_km2,
                best_drop,
                if has_step { best_at.to_string() } else { "—".into() },
                best_a,
                t_col,
                gate
            );
        }
        let st = sorted(tcols.clone());
        eprintln!(
            "\n   ⇒ **{fam1} lakes with a step below the sill · {fam2} with none** · T_col (passes to \
         walk the whole {WALK}-cell path at c = k·dt·√A) p10/p50/p90 **{:.4} / {:.4} / {:.4}**",
            pct(&st, 0.10),
            pct(&st, 0.50),
            pct(&st, 0.90)
        );
        eprintln!(
            "   ⇒ the shipped budget is **{} passes**. The reviewer's estimate was **5–20 k_time**.",
            sp.iterations
        );
    }
    eprintln!("\n==========  end Finding 115 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
