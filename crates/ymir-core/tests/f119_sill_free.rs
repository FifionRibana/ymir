//! ADR 0001 Finding 119 — **separate the floor from the sill: A1 protects the floor, lift it on
//! the outlet only.**
//!
//! Finding 118 measured that A1 (`depression_floor`) freezes the sill of every family-1 lake at
//! its input height (closure OFF cuts 8 of 18 sills by 18–702 m; closure ON cuts ≈ 0). This bench
//! lifts A1 on the sill cell, then on the sill + d@1…d@k, and asks whether the lakes empty while
//! the over-dug class (Finding 108's Δ gate) stays at zero.
//!
//! ## Code read, before the measurement (A)
//!
//! The A1 line is `filled > field` — **strict** — on the flow the incision recomputes on its own
//! surface each pass (`compute_flow(&field, &flow_cfg)`, exact priority-flood fill with no
//! epsilon, `terrain/flow.rs`). A spill cell has `filled == field`, so A1 **cannot fire on the sill
//! itself**. `incise_with_floor` never reads `lake_map`. Any protection of a sill is therefore
//! transmitted through its receiver chain: the implicit sweep lands each cell on its receiver's
//! already-updated height, and a receiver A1 skips keeps its height. The A-diagnostic below reads,
//! on the incision's own routing, which cells of sill → d@10 carry `filled > raw`.
//!
//! ⛔ **The first run refuted that reading — and Findings 114–118's "sill" with it.** It found A1
//! firing on the "sill" **18 of 18** on the eroded-ON surface. The reason is the instrument:
//! `Lake::outlet` is defined in `lakes/detection.rs` as *"lake cell whose D8 direction points to a
//! non-lake cell"* — the LAST WATER CELL, `filled − raw > 1e-6` by definition. Every "S at the
//! sill" since Finding 114 was the slope from a water cell to the spill point above it, clamped at
//! zero. The col used below is that cell's RECEIVER on the inventory's own flow, checked to lie
//! out of the lake; the lake's own cells are refused from every mask (A1 on the floor unchanged).
//!
//! ## Rule 12, declared
//!
//! * **The sills** are the closure-ON breached inventory's (`l.base.outlet`, as Findings 117/118).
//!   That inventory is the OUTPUT of the world being modified: the mask is built from the delivered
//!   field, which a production form could not do without a second pipeline (Finding 105). Bench
//!   seam, not a candidate — declared circularity.
//! * **The paths** d@1…d@k are walked on the incision's own routing (`compute_flow` at
//!   `SEA`, no flat perturbation — NOT the drainage chain's), on TWO surfaces: the pre-incision
//!   field (≈ the pass-1 input; it also carries the bathymetry, which only moves below-sea cells)
//!   and the eroded-ON field (≈ the pass-2 input). The mask is their union. The exact pass-2
//!   surface is not reachable from a bench without a new seam; declared.
//! * **Lake status** is read at the lake's FLOOR cell (argmin of the ON breached field over its
//!   ON footprint) in each world's own assembled inventory: `lake_map == 0` there ⇒ emptied.
//! * **The Δ class** is Finding 109's reading (bodies from `f95_criteria`, `fill_field_m` at the
//!   floor, eroded stage) with the breached stage beside it (what the `.ymir` export ships).
//!
//! ## Worlds (rule 18)
//!
//! OFF (no closure) · ON (shipped closure, A1 everywhere) · ON-free k = 0 / 1 / 3 · and one
//! BRACKET: ON with A1 nowhere (= the all-true mask, by the rule-13 test in `stream_power.rs`) —
//! the ceiling of any A1 lifting, not a candidate.
//!
//! Run: cargo test -p ymir-core --release --test f119_sill_free -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, build_field_with_a1_exempt, fill_field_m,
    land_u16, majority, over_dug_depression, pct, set_dump, sorted, take_bodies, to_mask,
};
use std::sync::Arc;
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, SegmentKind, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, FlowResult, breach_monotone, compute_flow,
};

const DOMAIN_KM: f32 = 400.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const DELIVERED_P50_M: f32 = 424.2;
const DIAG_K: usize = 10;

/// Finding 117-A's terminal base per family-1 lake, keyed by area (km²). A CITATION: it is matched
/// to this run's inventory by area and a mismatch is printed, not hidden.
const F117_BASE: [(f32, &str); 18] = [
    (275.11, "F2"),
    (198.72, "F2"),
    (181.98, "F2"),
    (145.30, "SEA"),
    (138.60, "lake"),
    (87.84, "F2"),
    (75.62, "SEA"),
    (54.44, "lake"),
    (48.83, "SEA"),
    (35.81, "F2"),
    (23.67, "F2"),
    (13.73, "lake"),
    (11.55, "SEA"),
    (8.87, "SEA"),
    (8.05, "lake"),
    (7.58, "F2"),
    (6.22, "F2"),
    (5.72, "lake"),
];

fn hash(f: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for v in &f.data {
        for b in v.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
    }
    h
}

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

/// The incision's OWN routing on a surface: `incise_with_floor` calls exactly this.
fn incision_flow(g: &GridF32) -> FlowResult {
    compute_flow(g, &FlowConfig { sea_level: SEA, ..Default::default() })
}

/// d@1…d@k from `start` down `flow.direction` (the start itself is not included).
fn path(flow: &FlowResult, start: usize, k: usize, w: usize, h: usize) -> Vec<usize> {
    let mut out = Vec::with_capacity(k);
    let mut c = start;
    for _ in 0..k {
        let d = flow.direction[c];
        if d == DIR_NONE {
            break;
        }
        let nx = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        c = ny * w + nx;
        out.push(c);
    }
    out
}

/// Indices j ∈ 0..=k (0 = the sill) whose cell A1 would skip on this surface: `filled > raw`.
fn blockers(g: &GridF32, flow: &FlowResult, sill: usize, p: &[usize]) -> Vec<usize> {
    std::iter::once(sill)
        .chain(p.iter().copied())
        .enumerate()
        .filter(|&(_, c)| flow.filled.data[c] > g.data[c])
        .map(|(j, _)| j)
        .collect()
}

fn fmt_blk(b: &[usize]) -> String {
    if b.is_empty() {
        "—".into()
    } else {
        b.iter().map(|j| j.to_string()).collect::<Vec<_>>().join(",")
    }
}

/// One family-1 lake of the ON inventory, with the cells every world is read at.
///
/// ⛔ **`Lake::outlet` is a WATER cell** (`lakes/detection.rs`: *"lake cell whose D8 direction
/// points to a non-lake cell"*, water = `filled − raw > 1e-6`). The first run of this bench read it
/// as the sill and found A1 firing on it 18 of 18 — by construction. The SILL is its receiver on
/// the inventory's own flow: the first cell out of the lake.
struct Lake {
    id: u32,
    area: f32,
    level: f32,
    /// the last water cell (`Lake::outlet`) — read by Findings 114–118 as "the sill"
    outlet: usize,
    /// the sill: `outlet`'s receiver on the inventory's own flow
    col: usize,
    floor: usize,
    base: &'static str,
}

/// One world's per-lake reading.
#[derive(Clone, Copy)]
struct Row {
    col_e: f32,
    col_b: f32,
    status: &'static str,
    area: f32,
    level: f32,
}

struct Totals {
    lakes_all: usize,
    lakes_fam1: usize,
    lake_pct: f32,
    canyons_e: usize,
    canyons_b: usize,
    scanned: usize,
    coast: i64,
    r8: f32,
}

/// The receiver of `c` on `flow`, or `None` at a base.
fn recv(flow: &FlowResult, c: usize, w: usize, h: usize) -> Option<usize> {
    path(flow, c, 1, w, h).first().copied()
}

#[test]
#[ignore]
fn f119_sill_free() {
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 119 . A1 lifted on the sill only  ==========");

    let on_k = Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };

    // ── the ON world and its inventory (Findings 117/118's) ─────────────────
    let ero_on = build_field_seed(on_k, PSEED);
    let (w, h) = (ero_on.width, ero_on.height);
    let n = w * h;
    // the inventory's own stage, kept for the instrument check below
    let (lakes, own_map, col_checks, bre_on) = {
        let d = c1_drainage_windowed(&ero_on, None, &dcfg, &ss, DOMAIN_KM);
        let bre = breach_monotone(&ero_on, &d.flow.filled, &d.lake_map, SEA, w, h);
        let cl = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dclim =
            DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        let dr = assemble_hd_drainage(
            &bre,
            &dclim,
            Some(d),
            &dcfg,
            &ss,
            DOMAIN_KM,
            GEO_RATIO,
            None,
            false,
        )
        .drainage;
        let mut fam1: Vec<_> =
            dr.lakes.iter().filter(|l| l.area_km2 >= 1.0 && l.base.id < 1_000_001).collect();
        fam1.sort_by(|a, b| b.area_km2.total_cmp(&a.area_km2));
        let mut checks = Vec::new();
        let lakes: Vec<Lake> = fam1
            .iter()
            .filter_map(|l| {
                let outlet = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
                if outlet >= n {
                    return None;
                }
                let col = recv(&dr.flow, outlet, w, h)?;
                let floor = (0..n)
                    .filter(|&c| dr.lake_map[c] == l.base.id)
                    .min_by(|&a, &b| bre.data[a].total_cmp(&bre.data[b]))
                    .unwrap_or(outlet);
                let base = F117_BASE
                    .iter()
                    .find(|(a, _)| (a - l.area_km2).abs() < 0.01)
                    .map_or("?", |(_, b)| *b);
                // instrument check: the outlet is water of THIS lake, the col is not, and the col
                // sits AT the lake level (it is the spill point)
                checks.push((
                    dr.lake_map[outlet] == l.base.id,
                    dr.lake_map[col],
                    l.level_m - m(&bre, outlet),
                    m(&bre, col) - l.level_m,
                ));
                Some(Lake {
                    id: l.base.id,
                    area: l.area_km2,
                    level: l.level_m,
                    outlet,
                    col,
                    floor,
                    base,
                })
            })
            .collect();
        (lakes, dr.lake_map.clone(), checks, bre)
    };
    let unmatched = lakes.iter().filter(|l| l.base == "?").count();
    let nl = lakes.len();
    eprintln!(
        "\n   ON inventory: **{nl}** family-1 lakes ≥ 1 km² · Finding 117 base matched by area on {} \
         ({unmatched} unmatched){}",
        nl - unmatched,
        if nl == 18 && unmatched == 0 { "" } else { " — **inventory differs from F117**" }
    );
    let outlet_water = col_checks.iter().filter(|c| c.0).count();
    let col_other_lake = col_checks.iter().filter(|c| c.1 != 0).count();
    let od = sorted(col_checks.iter().map(|c| c.2).collect());
    let cd = sorted(col_checks.iter().map(|c| c.3).collect());
    eprintln!(
        "   ⛔ INSTRUMENT · `Lake::outlet` is a cell of its own lake: **{outlet_water} of {nl}** · \
         its depth below the lake level p10/p50/p90 **{:.2} / {:.2} / {:.2} m**\
         \n   ⛔ INSTRUMENT · the col (outlet's receiver, inventory flow) is water of another lake: \
         **{col_other_lake} of {nl}** · col − level p10/p50/p90 **{:+.2} / {:+.2} / {:+.2} m**",
        pct(&od, 0.10),
        pct(&od, 0.50),
        pct(&od, 0.90),
        pct(&cd, 0.10),
        pct(&cd, 0.50),
        pct(&cd, 0.90)
    );

    // ── A · the incision's routing on the two surfaces, and who A1 skips ────
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let f_pre = incision_flow(&pre);
    let f_on = incision_flow(&ero_on);
    eprintln!(
        "\n   ── A · cells A1 SKIPS (`filled > raw`, incision routing) on col = 0 … d@{DIAG_K} ──\
         \n   {:>8} {:>5} {:>6} {:>20} {:>20}  {:>9} {:>9} {:>7}",
        "area km²", "base", "out.w", "pre-incision", "eroded ON", "S_col ON", "S_col BR", "knick@"
    );
    let (mut out_on, mut col_pre, mut col_on, mut none_pre, mut none_on) = (0, 0, 0, 0, 0);
    let mut paths: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    let mut skipped: std::collections::HashSet<usize> = std::collections::HashSet::new();
    // the breached field's own drainage, for S at the col on the exported stage
    let f_bre = incision_flow(&bre_on);
    let s_step = |g: &GridF32, f: &FlowResult, c: usize| -> f32 {
        match recv(f, c, w, h) {
            None => 0.0,
            Some(r) => {
                let diag = (c % w != r % w) && (c / w != r / w);
                let len = if diag {
                    CELL_KM * 1000.0 * std::f32::consts::SQRT_2
                } else {
                    CELL_KM * 1000.0
                };
                (m(g, c) - m(g, r)).max(0.0) / len
            }
        }
    };
    for l in &lakes {
        let pp = path(&f_pre, l.col, DIAG_K, w, h);
        let po = path(&f_on, l.col, DIAG_K, w, h);
        let bp = blockers(&pre, &f_pre, l.col, &pp);
        let bo = blockers(&ero_on, &f_on, l.col, &po);
        let cells_p: Vec<usize> = std::iter::once(l.col).chain(pp.iter().copied()).collect();
        let cells_o: Vec<usize> = std::iter::once(l.col).chain(po.iter().copied()).collect();
        skipped.extend(bp.iter().map(|&j| cells_p[j]));
        skipped.extend(bo.iter().map(|&j| cells_o[j]));
        let ow = f_on.filled.data[l.outlet] > ero_on.data[l.outlet];
        out_on += ow as usize;
        col_pre += bp.contains(&0) as usize;
        col_on += bo.contains(&0) as usize;
        none_pre += bp.is_empty() as usize;
        none_on += bo.is_empty() as usize;
        // S at the col itself, and the first cell (from the col, eroded ON) with S > 0
        let s_on = s_step(&ero_on, &f_on, l.col);
        let s_br = s_step(&bre_on, &f_bre, l.col);
        let knick = std::iter::once(l.col)
            .chain(po.iter().copied())
            .position(|c| s_step(&ero_on, &f_on, c) > 1e-4);
        eprintln!(
            "   {:>8.2} {:>5} {:>6} {:>20} {:>20}  {:>9.5} {:>9.5} {:>7}",
            l.area,
            l.base,
            if ow { "water" } else { "DRY" },
            fmt_blk(&bp),
            fmt_blk(&bo),
            s_on,
            s_br,
            knick.map_or(format!(">{DIAG_K}"), |j| j.to_string())
        );
        paths.push((pp, po));
    }
    eprintln!(
        "   (out.w = the last water cell is water on the eroded-ON incision surface · S_col = slope \
         from the col to its receiver, eroded ON / breached ON · knick@ = first index from the col \
         with S > 0, eroded ON)\
         \n   ⇒ the old \"sill\" (last water cell) is water on eroded ON: **{out_on} of {nl}** — A1's \
         job, not a defect\
         \n   ⇒ A1 fires AT THE COL: pre **{col_pre} of {nl}** · eroded ON **{col_on} of {nl}** · no \
         skipped cell on col…d@{DIAG_K}: pre **{none_pre}** · eroded ON **{none_on}**"
    );
    drop(f_pre);
    drop(f_on);
    drop(f_bre);

    // ── the masks: col + d@1…d@k, never a cell of the lake itself ───────────
    let mask_for = |k: usize| -> (Arc<Vec<bool>>, usize) {
        let mut mk = vec![false; n];
        let mut refused = 0usize;
        for (l, (pp, po)) in lakes.iter().zip(&paths) {
            for c in std::iter::once(l.col)
                .chain(pp.iter().take(k).copied())
                .chain(po.iter().take(k).copied())
            {
                if own_map[c] == l.id {
                    refused += 1; // A1 on the floor unchanged: the lake's own cells are never lifted
                } else {
                    mk[c] = true;
                }
            }
        }
        (Arc::new(mk), refused)
    };

    // ── the pre-incision references every world is scored against ──────────
    let auth = spurs(&land_u16(&pre, &ss), w);
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc =
            DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg, &ss, &dc, DOMAIN_KM).0
    };

    // ── one world, read ─────────────────────────────────────────────────────
    let read_world = |g: &GridF32| -> (Vec<Row>, Totals) {
        let d = c1_drainage_windowed(g, None, &dcfg, &ss, DOMAIN_KM);
        let bre = breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, w, h);
        // the Δ class, exactly as Finding 109 read it (eroded climate for both fills)
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate {
            precip_internal: &cl_e.precipitation,
            temperature: &cl_e.temperature,
        };
        let (fill_del, under) = fill_field_m(g, &dcfg, &ss, &dc_e, DOMAIN_KM);
        let (fill_bre, _) = fill_field_m(&bre, &dcfg, &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let cr =
            common::f95_criteria(g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg, cell_km2, n2m, w, h);
        set_dump(false);
        let (mut ce, mut cb) = (0usize, 0usize);
        for b in take_bodies() {
            let floor = *b
                .cells
                .iter()
                .min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c]))
                .expect("non-empty");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                ce += 1;
            }
            if over_dug_depression(fill_bre[floor] - fill_pre[floor], b.rim) {
                cb += 1;
            }
        }
        let _ = under;
        drop(fill_del);
        drop(fill_bre);
        let coast = spurs(&land_u16(&bre, &ss), w) as i64 - auth as i64;
        let land: Vec<bool> = (0..n).map(|k| g.data[k] > SEA).collect();
        let r8 = aniso(g, &land, 16).r8;
        drop(land);
        // the world's OWN assembled inventory (breached climate, as Findings 117/118)
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate {
            precip_internal: &cl_b.precipitation,
            temperature: &cl_b.temperature,
        };
        let dr = assemble_hd_drainage(
            &bre,
            &dc_b,
            Some(d),
            &dcfg,
            &ss,
            DOMAIN_KM,
            GEO_RATIO,
            None,
            false,
        )
        .drainage;
        let rows = lakes
            .iter()
            .map(|l| {
                let id = dr.lake_map[l.floor];
                let found = dr.lakes.iter().find(|x| x.base.id == id);
                let (status, area, level) = match (id, found) {
                    (0, _) => ("EMPTIED", 0.0, f32::NAN),
                    (_, Some(x)) if x.area_km2 < 0.5 * l.area => ("shrunk", x.area_km2, x.level_m),
                    (_, Some(x)) if x.area_km2 > 1.5 * l.area => ("GREW", x.area_km2, x.level_m),
                    (_, Some(x)) => ("same", x.area_km2, x.level_m),
                    (_, None) => ("no-summary", 0.0, f32::NAN),
                };
                Row { col_e: m(g, l.col), col_b: m(&bre, l.col), status, area, level }
            })
            .collect();
        let tot = Totals {
            lakes_all: dr.lakes.iter().filter(|x| x.area_km2 >= 1.0).count(),
            lakes_fam1: dr
                .lakes
                .iter()
                .filter(|x| x.area_km2 >= 1.0 && x.base.id < 1_000_001)
                .count(),
            lake_pct: cr.lake_pct,
            canyons_e: ce,
            canyons_b: cb,
            scanned: cr.scanned,
            coast,
            r8,
        };
        (rows, tot)
    };

    // ── the worlds ──────────────────────────────────────────────────────────
    let h_on = hash(&ero_on);
    let t = Instant::now();
    let (rows_on, tot_on) = read_world(&ero_on);
    eprintln!("   ⏱ ON read in {:.1} s", t.elapsed().as_secs_f64());
    let mut worlds: Vec<(String, Vec<Row>, Totals)> = Vec::new();
    {
        let g = build_field_seed(Knobs::passes(2), PSEED);
        let (r, tt) = read_world(&g);
        worlds.push(("OFF (no closure)".into(), r, tt));
    }
    // ON sits second so every table reads OFF → ON → the gesture
    let mut mask_sizes = Vec::new();
    let mut k0_identical = false;
    let mut gesture: Vec<(String, Vec<Row>, Totals)> = Vec::new();
    for k in [0usize, 1, 3] {
        let (mk, refused) = mask_for(k);
        let cells = mk.iter().filter(|&&b| b).count();
        let effective = (0..n).filter(|&c| mk[c] && skipped.contains(&c)).count();
        mask_sizes.push((k, cells));
        eprintln!(
            "
   mask k = {k}: **{cells}** cells · of which A1 skips on the pre or eroded-ON surface              **{effective}** · lake cells refused **{refused}**"
        );
        let g = build_field_with_a1_exempt(on_k, PSEED, mk);
        let hk = hash(&g);
        let diff = (0..n).filter(|&c| g.data[c] != ero_on.data[c]).count();
        eprintln!(
            "\n   ⛔ k = {k}: mask **{cells}** cells · field vs ON **{}** ({diff} cells differ)",
            if hk == h_on { "BIT-IDENTICAL" } else { "DIFFERENT" }
        );
        if k == 0 {
            k0_identical = hk == h_on;
            if k0_identical {
                eprintln!(
                    "      ⇒ k = 0 is a NO-OP by construction: A1 never fired on any mask cell in either \
                     pass. Its row IS the ON row; not re-read."
                );
                continue;
            }
        }
        let (r, tt) = read_world(&g);
        gesture.push((format!("ON-free k = {k}"), r, tt));
    }
    {
        let g = build_field_seed(
            Knobs { slope_floor_uk: Some(S_EQ), depression_floor: false, ..Knobs::passes(2) },
            PSEED,
        );
        let (r, tt) = read_world(&g);
        gesture.push(("BRACKET A1 nowhere".into(), r, tt));
    }
    let mut all: Vec<(String, Vec<Row>, Totals)> = Vec::new();
    all.append(&mut worlds);
    all.push(("ON (shipped closure)".into(), rows_on.clone(), tot_on));
    all.append(&mut gesture);

    // ── B · per lake, per world ─────────────────────────────────────────────
    let on_idx = 1usize;
    for (name, rows, _) in all.iter() {
        eprintln!(
            "\n   ── B · {name} — per lake (Δ vs ON at the ON COL cell; status at the ON floor cell) ──\
             \n   {:>8} {:>5} {:>9} {:>9} {:>9} {:>9} {:>11} {:>9} {:>9}",
            "area km²",
            "base",
            "col_e m",
            "Δcol_e",
            "Δcol_b",
            "re-raise",
            "status",
            "area",
            "Δlevel"
        );
        for (i, l) in lakes.iter().enumerate() {
            let r = rows[i];
            let o = all[on_idx].1[i];
            eprintln!(
                "   {:>8.2} {:>5} {:>9.2} {:>+9.2} {:>+9.2} {:>+9.2} {:>11} {:>9.2} {:>+9.2}",
                l.area,
                l.base,
                r.col_e,
                r.col_e - o.col_e,
                r.col_b - o.col_b,
                r.col_b - r.col_e,
                r.status,
                r.area,
                r.level - l.level
            );
        }
    }

    // ── the summary, one line per world ─────────────────────────────────────
    eprintln!(
        "\n   ── summary ──\n   {:<22} {:>7} {:>7} {:>7} {:>11} {:>11} {:>10} {:>6} {:>6} {:>7} {:>9} {:>6} {:>7}",
        "world",
        "emptied",
        "shrunk",
        "moved",
        "Δcol_e p50",
        "Δcol_b p50",
        "re-raise50",
        "≥1km²",
        "fam1",
        "lake %",
        "canyons",
        "coast",
        "R8"
    );
    for (name, rows, tt) in &all {
        let emptied = rows.iter().filter(|r| r.status == "EMPTIED").count();
        let shrunk = rows.iter().filter(|r| r.status == "shrunk").count();
        let moved: Vec<f32> = rows
            .iter()
            .zip(&all[on_idx].1)
            .map(|(r, o)| r.col_e - o.col_e)
            .filter(|d| d.abs() > 0.01)
            .collect();
        let de = sorted(rows.iter().zip(&all[on_idx].1).map(|(r, o)| r.col_e - o.col_e).collect());
        let db = sorted(rows.iter().zip(&all[on_idx].1).map(|(r, o)| r.col_b - o.col_b).collect());
        let rr = sorted(rows.iter().map(|r| r.col_b - r.col_e).collect());
        eprintln!(
            "   {name:<22} {emptied:>7} {shrunk:>7} {:>7} {:>+11.2} {:>+11.2} {:>+10.2} {:>6} {:>6} {:>7.2} {:>4}e/{:<3}b {:>+6} {:>7.4}",
            moved.len(),
            pct(&de, 0.50),
            pct(&db, 0.50),
            pct(&rr, 0.50),
            tt.lakes_all,
            tt.lakes_fam1,
            tt.lake_pct,
            tt.canyons_e,
            tt.canyons_b,
            tt.coast,
            tt.r8
        );
    }
    eprintln!(
        "   (moved = cols whose eroded height differs from ON by > 1 cm · canyons = Finding 108 Δ \
         class, eroded / breached stage, of {} bodies scanned on ON)",
        all[on_idx].2.scanned
    );
    eprintln!(
        "   masks: {} · k = 0 identical to ON: **{}**",
        mask_sizes.iter().map(|(k, c)| format!("k={k} {c} cells")).collect::<Vec<_>>().join(" · "),
        if k0_identical { "YES" } else { "no" }
    );
    eprintln!("\n==========  end Finding 119 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 119 — **the "emptied" status, verified by footprint, and where the count went.**
///
/// The main bench reads a lake's status at ONE cell (its ON floor) and reported 3 (k = 1) and 6
/// (k = 3) lakes emptied while the lake count stayed at 29 / 18 and the lake fraction moved by
/// 0.02 points — and one "emptied" lake had its col RISE 0.12 m. A one-cell status cannot carry
/// that claim. This reads, per ON lake, the SHARE of its ON footprint still under water in each
/// world, the lake that covers it, and every lake ≥ 1 km² of the world that touches no ON lake.
///
/// Run: cargo test -p ymir-core --release --test f119_sill_free -- --ignored verify --nocapture
#[test]
#[ignore]
fn f119_verify() {
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 119 . the emptied status, by footprint  ==========");
    let on_k = Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };

    // the assembled inventory of any eroded field: (lake_map, lakes)
    let inventory = |g: &GridF32| {
        let (w, h) = (g.width, g.height);
        let d = c1_drainage_windowed(g, None, &dcfg, &ss, DOMAIN_KM);
        let bre = breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, w, h);
        let cl = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc =
            DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        assemble_hd_drainage(&bre, &dc, Some(d), &dcfg, &ss, DOMAIN_KM, GEO_RATIO, None, false)
            .drainage
    };

    let ero_on = build_field_seed(on_k, PSEED);
    let (w, h) = (ero_on.width, ero_on.height);
    let n = w * h;
    let dr_on = inventory(&ero_on);
    let mut fam1: Vec<_> =
        dr_on.lakes.iter().filter(|l| l.area_km2 >= 1.0 && l.base.id < 1_000_001).collect();
    fam1.sort_by(|a, b| b.area_km2.total_cmp(&a.area_km2));
    let lakes: Vec<(u32, f32, f32, usize)> = fam1
        .iter()
        .filter_map(|l| {
            let outlet = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
            recv(&dr_on.flow, outlet, w, h).map(|col| (l.base.id, l.area_km2, l.level_m, col))
        })
        .collect();
    // Finding 87-A's river: the highest-discharge spillway of the ON inventory (92.97 m³/s, 61 cells)
    let spill =
        |dr: &ymir_core::tectonics_c1::drainage::C1DrainageResult| -> Option<(f32, Vec<usize>)> {
            (0..dr.rivers.segments.len())
                .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway)
                .max_by(|&a, &b| {
                    dr.segment_discharge_m3s[a].total_cmp(&dr.segment_discharge_m3s[b])
                })
                .map(|i| {
                    let pts = &dr.rivers.segments[i].points;
                    (
                        dr.segment_discharge_m3s[i],
                        pts.iter().map(|p| p.1 as usize * w + p.0 as usize).collect(),
                    )
                })
        };
    let river_on = spill(&dr_on).expect("the ON inventory has a spillway");
    eprintln!(
        "   F87-A river on ON: **{:.2} m³/s** over **{}** cells (Finding 87-A: 92.97 m³/s, 61 cells)",
        river_on.0,
        river_on.1.len()
    );
    let on_map = dr_on.lake_map.clone();
    let on_big: std::collections::HashSet<u32> =
        dr_on.lakes.iter().filter(|l| l.area_km2 >= 1.0).map(|l| l.base.id).collect();
    drop(dr_on);

    // the same masks as the main bench: col + d@1…d@k on the pre and eroded-ON incision flows,
    // the lake's own cells refused
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let f_pre = incision_flow(&pre);
    let f_on = incision_flow(&ero_on);
    let paths: Vec<(Vec<usize>, Vec<usize>)> = lakes
        .iter()
        .map(|&(_, _, _, col)| (path(&f_pre, col, 3, w, h), path(&f_on, col, 3, w, h)))
        .collect();
    drop((f_pre, f_on, pre));
    let mask_for = |k: usize| -> Arc<Vec<bool>> {
        let mut mk = vec![false; n];
        for (&(id, _, _, col), (pp, po)) in lakes.iter().zip(&paths) {
            for c in std::iter::once(col)
                .chain(pp.iter().take(k).copied())
                .chain(po.iter().take(k).copied())
            {
                if on_map[c] != id {
                    mk[c] = true;
                }
            }
        }
        Arc::new(mk)
    };

    for (name, g) in [
        ("ON (control: 100 % by construction)", ero_on.clone()),
        ("ON-free k = 1", build_field_with_a1_exempt(on_k, PSEED, mask_for(1))),
        ("ON-free k = 3", build_field_with_a1_exempt(on_k, PSEED, mask_for(3))),
    ] {
        let dr = inventory(&g);
        let level_of: std::collections::HashMap<u32, (f32, f32)> =
            dr.lakes.iter().map(|l| (l.base.id, (l.area_km2, l.level_m))).collect();
        eprintln!(
            "\n   ── {name} ──\n   {:>8} {:>9} {:>8} {:>10} {:>9} {:>9}  {}",
            "ON km²", "wet %", "wet km²", "cover id", "its km²", "Δlevel", "verdict"
        );
        let mut truly_empty = 0usize;
        for &(id, area, level, _) in &lakes {
            let cells: Vec<usize> = (0..n).filter(|&c| on_map[c] == id).collect();
            let wet = cells.iter().filter(|&&c| dr.lake_map[c] != 0).count();
            // the world's lake covering most of the footprint
            let mut votes: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
            for &c in &cells {
                if dr.lake_map[c] != 0 {
                    *votes.entry(dr.lake_map[c]).or_default() += 1;
                }
            }
            let cover = votes.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))).map(|e| *e.0);
            let share = 100.0 * wet as f32 / cells.len().max(1) as f32;
            let (ca, cl) = cover.and_then(|i| level_of.get(&i).copied()).unwrap_or((0.0, f32::NAN));
            let verdict = if wet == 0 {
                truly_empty += 1;
                "EMPTY — no water left on the footprint"
            } else if share < 50.0 {
                "mostly drained"
            } else if cl < level - 0.5 {
                "wet, level LOWER"
            } else {
                "wet"
            };
            eprintln!(
                "   {area:>8.2} {share:>8.1}% {:>8.2} {:>10} {ca:>9.2} {:>+9.2}  {verdict}",
                wet as f32 * CELL_KM * CELL_KM,
                cover.map_or("—".into(), |i| i.to_string()),
                cl - level
            );
        }
        // lakes of this world that touch no ON lake at all
        let mut touch: std::collections::HashMap<u32, (usize, bool, f64, f64)> =
            std::collections::HashMap::new();
        for c in 0..n {
            let id = dr.lake_map[c];
            if id == 0 {
                continue;
            }
            let e = touch.entry(id).or_insert((0, false, 0.0, 0.0));
            e.0 += 1;
            e.1 |= on_map[c] != 0 && on_big.contains(&on_map[c]);
            e.2 += (c % w) as f64;
            e.3 += (c / w) as f64;
        }
        let mut news: Vec<(f32, f32, f64, f64, u32)> = touch
            .iter()
            .filter(|(id, v)| !v.1 && level_of.get(id).is_some_and(|x| x.0 >= 1.0))
            .map(|(&id, v)| {
                let (a, l) = level_of[&id];
                (a, l, v.2 / v.0 as f64 * CELL_KM as f64, v.3 / v.0 as f64 * CELL_KM as f64, id)
            })
            .collect();
        news.sort_by(|a, b| b.0.total_cmp(&a.0));
        eprintln!(
            "   ⇒ footprints with NO water left: **{truly_empty} of {}** · lakes ≥ 1 km² of this world \
             touching NO ON lake ≥ 1 km²: **{}**{}",
            lakes.len(),
            news.len(),
            news.iter()
                .map(|x| format!(
                    "\n      new lake {} · {:.2} km² · level {:.1} m · centroid ({:.1}, {:.1}) km{}",
                    x.4,
                    x.0,
                    x.1,
                    x.2,
                    x.3,
                    if x.4 >= 1_000_001 { " · below-sea family" } else { "" }
                ))
                .collect::<String>()
        );
        // ── negative control: the F87-A river unchanged ──
        let r = spill(&dr);
        let dz = river_on.1.iter().map(|&c| (m(&g, c) - m(&ero_on, c)).abs()).fold(0f32, f32::max);
        eprintln!(
            "   CONTROL · F87-A river: {} · same cells as ON: **{}** · max |Δz| along ON's river              (eroded) **{dz:.3} m**",
            r.as_ref().map_or("no spillway".into(), |x| format!(
                "**{:.2} m³/s** over **{}** cells",
                x.0,
                x.1.len()
            )),
            r.as_ref().is_some_and(|x| x.1 == river_on.1)
        );
    }
    eprintln!(
        "\n==========  end Finding 119 verify . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
