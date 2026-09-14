//! ADR 0001 Finding 83 block B2 — **the first lever: estuaries.**
//!
//! The bound of Finding 80 refuses the incision the right to cut below sea level ANYWHERE. A real
//! coast does not work like that: the trunks reach the sea in drowned valleys. `free_above_km2`
//! releases the bound above a drainage area `A_est` — binary first, a graded version only if the
//! binary one works.
//!
//! **`A_est` is a PROXY and is labelled one in the code.** Ymir has no sea-level history; a real
//! ria is a valley cut at a low stand and flooded at a high one. Drainage area is a correlate of
//! "this valley would have been drowned", not the cause. The sweep is in REAL km² of the
//! simulated grid — the geographic scale ratio (Finding 24) never reaches the incision.
//!
//! Run: cargo test -p ymir-core --release --test f83_estuary -- --ignored --nocapture

mod common;

use common::{
    A_C_CELLS, CELL_KM, CELL_KM2, Knobs, SEA, build_field, land_u16, pct, richardson_line, sorted,
    spectrum, to_mask,
};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{FlowConfig, compute_flow};

/// The sweep, in km² of the simulated grid. 50 km² is 20 976 cells at 8192², 1 000 km² is
/// 419 526 — 0.6 % of the 160 000 km² domain, so only the trunks.
const SWEEP: [f32; 3] = [50.0, 200.0, 1000.0];

/// A spur of the variant is ADDED when the pre-incision mask has no spur within this many cells
/// of its neck midpoint. Declared: at 3 cells two shores of the same headland are not confused.
const MATCH_CELLS: f32 = 3.0;
/// Dilation radius, in cells, of the excursion ARC searched for the channel at its head.
const HEAD_REACH: i32 = 6;

/// 8-connected components of a boolean mask; returns a label per cell (0 = none) and the count.
fn components(mask: &[bool], w: usize) -> (Vec<u32>, u32) {
    let n = mask.len();
    let h = n / w;
    let mut lab = vec![0u32; n];
    let mut next = 0u32;
    let mut stack: Vec<u32> = Vec::new();
    for s in 0..n {
        if !mask[s] || lab[s] != 0 {
            continue;
        }
        next += 1;
        lab[s] = next;
        stack.push(s as u32);
        while let Some(kk) = stack.pop() {
            let k = kk as usize;
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if mask[nk] && lab[nk] == 0 {
                            lab[nk] = next;
                            stack.push(nk as u32);
                        }
                    }
                }
            }
        }
    }
    (lab, next)
}

/// Multi-source BFS distance in cells from every `source` cell.
fn bfs_dist(source: &[bool], w: usize) -> Vec<u16> {
    let n = source.len();
    let h = n / w;
    let mut d = vec![u16::MAX; n];
    let mut q: Vec<u32> = (0..n).filter(|&k| source[k]).map(|k| k as u32).collect();
    for &k in &q {
        d[k as usize] = 0;
    }
    let mut i = 0usize;
    while i < q.len() {
        let k = q[i] as usize;
        i += 1;
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if d[nk] == u16::MAX {
                    d[nk] = d[k] + 1;
                    q.push(nk as u32);
                }
            }
        }
    }
    d
}

fn control(tag: &str, f: &GridF32, pre: &GridF32, n2m: f32) {
    let n = f.data.len();
    let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
    let n_land = land.iter().filter(|&&l| l).count();
    let hh = sorted((0..n).filter(|&k| land[k]).map(|k| (f.data[k] - SEA) * n2m).collect());
    let mean = hh.iter().map(|&x| x as f64).sum::<f64>() / hh.len().max(1) as f64;
    let flow = compute_flow(f, &FlowConfig { sea_level: SEA, ..Default::default() });
    let pits = (0..n).filter(|&k| flow.filled.data[k] > f.data[k] + 1e-7).count();
    let chan = (0..n).filter(|&k| land[k] && flow.accumulation.data[k] >= A_C_CELLS).count();
    // PAIRED work: only cells that are land in BOTH fields (Findings 63-64, 81-B3, 82-C1).
    let (mut wsum, mut wn) = (0.0f64, 0usize);
    for k in 0..n {
        if pre.data[k] > SEA && land[k] {
            wsum += ((pre.data[k] - f.data[k]) * n2m) as f64;
            wn += 1;
        }
    }
    let drowned = (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] <= SEA).count();
    eprintln!(
        "      CONTROL {tag:<22} land {n_land} | hypso mean {mean:.2} m p50 {:.2} | pits {pits} | \
         channel {:.3} % | PAIRED work {:.2} m ({wn} cells) | DROWNED vs pre {drowned}",
        pct(&hh, 0.50),
        100.0 * chan as f64 / n_land.max(1) as f64,
        wsum / wn.max(1) as f64
    );
}

#[test]
#[ignore]
fn f83_estuary() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 83 · B2 — the estuary lever  ==========");
    eprintln!(
        "   A_est PROXY, in km2 of the SIMULATED grid. Cell = {CELL_KM:.4} km, cell area \
         {CELL_KM2:.6} km2, so 50 km2 = {:.0} cells, 200 = {:.0}, 1000 = {:.0}.",
        50.0 / CELL_KM2,
        200.0 / CELL_KM2,
        1000.0 / CELL_KM2
    );

    let pre = build_field(Knobs::no_incision());
    let w = pre.width;
    let n = pre.data.len();
    let pre_land = land_u16(&pre, &ss);
    let pre_polys = marching_squares(&to_mask(&pre_land, w), 0.5);
    let (pre_spurs, _) = coast_spurs(&pre_polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    let pre_base = coast_shape_thresholds(&pre_polys, CELL_KM, 1.0, NECK_KM).count;
    let pre_cells = coast_shape_thresholds(&pre_polys, CELL_KM, 2.0 * CELL_KM, CELL_KM).count;
    eprintln!(
        "\n   PRE-INCISION authority: >= 1 km {pre_base} | >= 2 cells {pre_cells} | {} spurs \
         located for the ADDED test",
        pre_spurs.len()
    );

    // The reference the lever is measured against is the SHIPPED bound (A_est = None).
    let bnd = build_field(Knobs::shipped());
    let bnd_land = land_u16(&bnd, &ss);
    let bnd_sea: Vec<bool> = bnd_land.iter().map(|&l| !l).collect();

    for a_est in [None, Some(SWEEP[0]), Some(SWEEP[1]), Some(SWEEP[2])] {
        let tag = match a_est {
            None => "A_est = none (PRODUCTION)".to_string(),
            Some(a) => format!("A_est = {a:.0} km2"),
        };
        eprintln!("\n╔═══ {tag} ═══╗");
        let f = build_field(Knobs { a_est_km2: a_est, ..Knobs::shipped() });
        let land = land_u16(&f, &ss);
        let polys = marching_squares(&to_mask(&land, w), 0.5);
        let a1 = coast_shape_thresholds(&polys, CELL_KM, 1.0, NECK_KM);
        let a2 = coast_shape_thresholds(&polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
        let sp = spectrum(&polys);
        eprintln!(
            "   spurs: >= 1 km {} (Δ vs pre {:+}) | >= 2 cells {} (**Δ {:+}**) | coast {:.0} km | \
             spectrum peak {} at x{} white",
            a1.count,
            a1.count as i64 - pre_base as i64,
            a2.count,
            a2.count as i64 - pre_cells as i64,
            a1.coast_km,
            sp.map_or("n/a".into(), |s| format!("{}", s.3)),
            sp.map_or("n/a".into(), |s| format!("{:.2}", s.4))
        );
        richardson_line("  Richardson", &land, w, CELL_KM, 10.0);
        control(&tag, &f, &pre, n2m);

        let Some(a_km2) = a_est else {
            continue;
        };

        // ── which excursions are ADDED, and does each carry a channel at its head? ──
        let flow = compute_flow(&f, &FlowConfig { sea_level: SEA, ..Default::default() });
        let (spurs, _) = coast_spurs(&polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
        let mut added = 0usize;
        let mut with_head = 0usize;
        let mut head_areas: Vec<f32> = Vec::new();
        for s in &spurs {
            let near = pre_spurs
                .iter()
                .any(|p| (p.mid.0 - s.mid.0).hypot(p.mid.1 - s.mid.1) <= MATCH_CELLS);
            if near {
                continue;
            }
            added += 1;
            // The head is at the INNER end of the excursion, not at its neck: searching a disc
            // around the neck midpoint alone would miss the channel feeding a long ria. So the
            // search is the whole ARC of the excursion, dilated by `HEAD_REACH`.
            let pl = &polys[s.poly];
            let (lo, hi) = (s.i.min(s.j), s.i.max(s.j).min(pl.len().saturating_sub(1)));
            let mut best = 0.0f32;
            for pt in &pl[lo..=hi] {
                let (cx, cy) = (pt.0.round() as i32, pt.1.round() as i32);
                for dy in -HEAD_REACH..=HEAD_REACH {
                    for dx in -HEAD_REACH..=HEAD_REACH {
                        let (nx, ny) = (cx + dx, cy + dy);
                        if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= w {
                            continue;
                        }
                        let k = ny as usize * w + nx as usize;
                        if f.data[k] > SEA {
                            best = best.max(flow.accumulation.data[k]);
                        }
                    }
                }
            }
            let km2 = best * CELL_KM2;
            head_areas.push(km2);
            if km2 >= a_km2 {
                with_head += 1;
            }
        }
        let ha = sorted(head_areas);
        eprintln!(
            "   ADDED excursions (>= 2 cells, no pre-incision spur within {MATCH_CELLS} cells): \
             **{added}** | with A >= A_est at the head: **{with_head}** ({:.1} %) | head A km2 \
             p10 {:.2} p50 {:.2} p90 {:.2} max {:.1}",
            100.0 * with_head as f64 / added.max(1) as f64,
            pct(&ha, 0.10),
            pct(&ha, 0.50),
            pct(&ha, 0.90),
            ha.last().copied().unwrap_or(0.0)
        );
        eprintln!(
            "      ⚠️ the accumulation here is D8 (`compute_flow`); the GATE reads the MFD \
             accumulation (p = 2) the incision itself builds, which disperses. The two are not \
             the same number and the share above is therefore a LOWER bound on the gate's own view."
        );

        // ── the estuaries themselves: what the released bound actually drowned ──
        let drowned: Vec<bool> = (0..n).map(|k| bnd_land[k] && !land[k]).collect();
        let (lab, nc) = components(&drowned, w);
        let dist_sea = bfs_dist(&bnd_sea, w); // how far inland from the BOUNDED shoreline
        let pre_clamp =
            build_field(Knobs { a_est_km2: a_est, bathymetry_off: true, ..Knobs::shipped() });
        let mut cells = vec![0usize; nc as usize + 1];
        let mut reach = vec![0u16; nc as usize + 1];
        let mut mouth = vec![0usize; nc as usize + 1];
        let mut depth: Vec<Vec<f32>> = vec![Vec::new(); nc as usize + 1];
        for k in 0..n {
            let c = lab[k] as usize;
            if c == 0 {
                continue;
            }
            cells[c] += 1;
            reach[c] = reach[c].max(dist_sea[k]);
            depth[c].push(-(pre_clamp.data[k] - SEA) * n2m);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            let touches = (-1i32..=1).any(|dy| {
                (-1i32..=1).any(|dx| {
                    let (nx, ny) = (x + dx, y + dy);
                    nx >= 0
                        && ny >= 0
                        && (nx as usize) < w
                        && (ny as usize) < w
                        && bnd_sea[ny as usize * w + nx as usize]
                })
            });
            if touches {
                mouth[c] += 1;
            }
        }
        // an ESTUARY, declared: a drowned component of at least 4 cells that touches the sea.
        let est: Vec<usize> =
            (1..=nc as usize).filter(|&c| cells[c] >= 4 && mouth[c] > 0).collect();
        let lens = sorted(est.iter().map(|&c| reach[c] as f32 * CELL_KM).collect::<Vec<f32>>());
        let mouths = sorted(est.iter().map(|&c| mouth[c] as f32 * CELL_KM).collect::<Vec<f32>>());
        let areas = sorted(est.iter().map(|&c| cells[c] as f32 * CELL_KM2).collect::<Vec<f32>>());
        let deps = sorted(
            est.iter()
                .map(|&c| {
                    let s = sorted(depth[c].clone());
                    s[s.len() / 2]
                })
                .collect::<Vec<f32>>(),
        );
        eprintln!(
            "   ESTUARIES (drowned components >= 4 cells touching the bounded sea): **{}** of \
             {nc} drowned components; {} cells drowned in total",
            est.len(),
            (0..n).filter(|&k| drowned[k]).count()
        );
        eprintln!(
            "      length inland km  p50 {:.2}  p90 {:.2}  MAX {:.2}\n      mouth width km    p50 \
             {:.3}  p90 {:.3}  MAX {:.3}\n      area km2          p50 {:.3}  p90 {:.3}  MAX \
             {:.2}\n      depth m (PRE-CLAMP) p50 {:.2}  p90 {:.2}  MAX {:.2}",
            pct(&lens, 0.50),
            pct(&lens, 0.90),
            lens.last().copied().unwrap_or(0.0),
            pct(&mouths, 0.50),
            pct(&mouths, 0.90),
            mouths.last().copied().unwrap_or(0.0),
            pct(&areas, 0.50),
            pct(&areas, 0.90),
            areas.last().copied().unwrap_or(0.0),
            pct(&deps, 0.50),
            pct(&deps, 0.90),
            deps.last().copied().unwrap_or(0.0)
        );
        eprintln!(
            "      REAL rias for the plausibility anchor: Rance ~20 km long, 0.5-2 km wide; Aber \
             Wrac'h ~9 km, 0.3-1 km; Ria de Vigo ~35 km, 3-7 km."
        );
    }
}
