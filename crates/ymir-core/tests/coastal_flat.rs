//! ADR Finding 76 block A — does the coastal fur sit on a SEA-LEVEL FLAT?
//!
//! The candidate mechanism, read out of the code before measuring (`stream_power.rs:410` and
//! `:530`): the ocean is a fixed base node (`receiver[k] = k` when `field <= sea_level`), so a
//! land cell whose D8 receiver is an ocean cell takes `hr` = THAT OCEAN CELL'S HEIGHT, and the
//! implicit update `h ← (h + f·hr + …)/(1+f)` clamped to `[hr, ho]` drives it onto `hr` in ONE
//! step at Courant ≫ 1 (Finding 61). Next iteration it is itself ≤ sea level, so its upstream
//! neighbour is dragged in turn: **the coast eats inland, one cell per iteration, per channel.**
//!
//! If that is right: dent length = flat length, dent count = channels at the coast.
//! **Stop rule 1: if the correlation is below 0.5 the mechanism is false and the seam is not
//! written.**
//!
//! Run: cargo test -p ymir-core --release --test coastal_flat -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;
const TARGET: usize = 8192;
const CELL_KM: f32 = DOMAIN_KM / TARGET as f32;

pub fn build_field(incision: bool) -> GridF32 {
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
    let cfg = production_hd_config(&ProductionHdOpts {
        target_size: TARGET,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
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
    let mut cfg = cfg;
    if !incision {
        // FBM and both closures STAY ON — only the incision is removed. This is NOT the
        // Finding 75 REFERENCE, which also drops the FBM: separating the two is the point.
        cfg.stream_power = None;
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

fn pct(v: &[f32], p: f64) -> f32 {
    if v.is_empty() { f32::NAN } else { v[(((v.len() - 1) as f64) * p) as usize] }
}

#[test]
#[ignore]
fn coastal_flat() {
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 76 block A · is the fur a sea-level flat?  ==========");
    let f = build_field(true);
    let (w, h) = (f.width, f.height);
    let n = w * h;
    // Metres per unit of the normalised field — so every elevation below is REAL metres above
    // the sea, not a norm difference nobody can size.
    let norm_to_m =
        ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres(1.0, &ss)
            - ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres(0.0, &ss);
    let m_above_sea = |k: usize| (f.data[k] - SEA) * norm_to_m;
    eprintln!("   field {w}×{h}, 1 norm unit = {norm_to_m:.1} m, cell = {:.1} m", CELL_KM * 1000.0);

    let flow = compute_flow(&f, &FlowConfig { sea_level: SEA, ..Default::default() });
    let is_land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
    let recv = |k: usize| -> Option<usize> {
        let d = flow.direction[k];
        if d == DIR_NONE {
            return None;
        }
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let nx = (x + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = (y + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        Some(ny * w + nx)
    };

    // ── distance to the sea ALONG THE FLOW, in cells: a BFS upstream from every land cell
    //    whose receiver is ocean. This is `d_i`, the quantity the seam's floor would use.
    let mut dsea = vec![u32::MAX; n];
    let mut queue: Vec<u32> = Vec::new();
    let mut donors: Vec<Vec<u32>> = vec![Vec::new(); n];
    for k in 0..n {
        if is_land[k] {
            if let Some(r) = recv(k) {
                donors[r].push(k as u32);
                if !is_land[r] {
                    dsea[k] = 1;
                    queue.push(k as u32);
                }
            }
        }
    }
    let n_mouth = queue.len();
    let mut qi = 0usize;
    while qi < queue.len() {
        let k = queue[qi] as usize;
        qi += 1;
        for &dk in &donors[k] {
            let d = dk as usize;
            if dsea[d] == u32::MAX {
                dsea[d] = dsea[k] + 1;
                queue.push(dk);
            }
        }
    }
    let reached = dsea.iter().filter(|&&d| d != u32::MAX).count();
    let land_n = is_land.iter().filter(|&&l| l).count();
    eprintln!(
        "   {n_mouth} land cells drain DIRECTLY into an ocean cell; {reached} of {land_n} land \
         cells ({:.1} %) have a flow path to the sea",
        100.0 * reached as f64 / land_n as f64
    );

    // ── A2 · bed elevation by distance from the coast ─────────────────────────
    eprintln!("\n── A2 · bed elevation of a land cell, by flow distance to the sea ──");
    eprintln!(
        "   {:>8} {:>9} {:>10} {:>10} {:>10} {:>10}",
        "d cells", "n", "p10 m", "MEDIAN m", "p90 m", "mean m"
    );
    for d in [1u32, 2, 3, 5, 10, 20] {
        let mut v: Vec<f32> = (0..n).filter(|&k| dsea[k] == d).map(m_above_sea).collect();
        if v.is_empty() {
            eprintln!("   {d:>8} {:>9} — EMPTY (rule 10: not a reading)", 0);
            continue;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mean = v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64;
        eprintln!(
            "   {d:>8} {:>9} {:>10.3} {:>10.3} {:>10.3} {:>10.3}",
            v.len(),
            pct(&v, 0.10),
            pct(&v, 0.50),
            pct(&v, 0.90),
            mean
        );
    }

    // ── A1 · spur length against the length of the sea-level flat it sits on ──
    // The flat: walk UPSTREAM from the mouth cell (the land cell draining into the ocean),
    // always taking the highest-accumulation donor, counting consecutive cells at or below
    // `sea + FLAT_TOL_M`. That is the reach the drag would have levelled.
    const FLAT_TOL_M: f32 = 0.01;
    let polys = marching_squares(&f, SEA);
    let (spurs, _) = coast_spurs(&polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    eprintln!("\n── A1 · spur length vs flat length, over {} spurs ≥ 2 cells ──", spurs.len());
    assert!(spurs.len() >= 100, "rule 10: population {} is too small to read", spurs.len());

    // index the mouth cells so a spur can be attached to the nearest one
    let mouths: Vec<u32> =
        (0..n as u32).filter(|&k| is_land[k as usize] && dsea[k as usize] == 1).collect();
    // grid buckets of 16 cells so "nearest mouth" is O(1) per spur
    const BS: usize = 16;
    let gw = w.div_ceil(BS);
    let mut bucket: Vec<Vec<u32>> = vec![Vec::new(); gw * h.div_ceil(BS)];
    for &mk in &mouths {
        let (x, y) = (mk as usize % w, mk as usize / w);
        bucket[(y / BS) * gw + (x / BS)].push(mk);
    }
    let flat_len = |start: usize| -> usize {
        let mut cur = start;
        let mut len = 0usize;
        for _ in 0..200 {
            if !is_land[cur] || m_above_sea(cur) > FLAT_TOL_M {
                break;
            }
            len += 1;
            match donors[cur].iter().copied().max_by(|&a, &b| {
                flow.accumulation.data[a as usize]
                    .partial_cmp(&flow.accumulation.data[b as usize])
                    .unwrap_or(std::cmp::Ordering::Equal)
            }) {
                Some(d) => cur = d as usize,
                None => break,
            }
        }
        len
    };
    let mut pairs: Vec<(f32, f32)> = Vec::new(); // (spur len in CELLS, flat len in CELLS)
    let mut unmatched = 0usize;
    for sp in &spurs {
        let (mx, my) = (sp.mid.0 as usize, sp.mid.1 as usize);
        let (bx, by) = (mx / BS, my / BS);
        let mut best: Option<(f32, u32)> = None;
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (bx as i32 + dx, by as i32 + dy);
                if nx < 0 || ny < 0 || nx as usize >= gw || (ny as usize) >= h.div_ceil(BS) {
                    continue;
                }
                for &mk in &bucket[ny as usize * gw + nx as usize] {
                    let (px, py) = ((mk as usize % w) as f32, (mk as usize / w) as f32);
                    let d2 = (px - sp.mid.0).powi(2) + (py - sp.mid.1).powi(2);
                    if best.map_or(true, |b| d2 < b.0) {
                        best = Some((d2, mk));
                    }
                }
            }
        }
        match best {
            Some((_, mk)) => {
                pairs.push((sp.len_km / CELL_KM, flat_len(mk as usize) as f32));
            }
            None => unmatched += 1,
        }
    }
    eprintln!(
        "   {} spurs matched to a mouth within 1.5 buckets ({BS} cells); {unmatched} unmatched \
         (no channel mouth nearby — those are NOT channel dents)",
        pairs.len()
    );
    let r = {
        let nn = pairs.len() as f64;
        let (sx, sy): (f64, f64) =
            pairs.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0 as f64, a.1 + p.1 as f64));
        let (mx, my) = (sx / nn, sy / nn);
        let (mut sxy, mut sxx, mut syy) = (0.0f64, 0.0f64, 0.0f64);
        for p in &pairs {
            let (dx, dy) = (p.0 as f64 - mx, p.1 as f64 - my);
            sxy += dx * dy;
            sxx += dx * dx;
            syy += dy * dy;
        }
        (sxy / (sxx.sqrt() * syy.sqrt()), mx, my, sxy / sxx.max(1e-12))
    };
    let mut fl: Vec<f32> = pairs.iter().map(|p| p.1).collect();
    let mut sl: Vec<f32> = pairs.iter().map(|p| p.0).collect();
    fl.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sl.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eprintln!(
        "   FLAT length (cells ≤ sea + {FLAT_TOL_M} m walked upstream): p10 {:.1} MEDIAN {:.1} \
         p90 {:.1} max {:.1} | zero-length {:.1} %",
        pct(&fl, 0.10),
        pct(&fl, 0.50),
        pct(&fl, 0.90),
        fl.last().copied().unwrap_or(0.0),
        100.0 * fl.iter().filter(|&&x| x == 0.0).count() as f64 / fl.len() as f64
    );
    eprintln!(
        "   SPUR length (cells): p10 {:.1} MEDIAN {:.1} p90 {:.1} max {:.1}",
        pct(&sl, 0.10),
        pct(&sl, 0.50),
        pct(&sl, 0.90),
        sl.last().copied().unwrap_or(0.0)
    );
    eprintln!(
        "   ⇒ **r = {:.3}** (mean spur {:.2} cells, mean flat {:.2} cells, slope flat/spur \
         {:.3})",
        r.0, r.1, r.2, r.3
    );

    // ── A3 · the counter-hypothesis: is the zero a CLASSIFICATION, not a flat? ─
    // The drag clamps to `hr`, the RECEIVER's height — not to `sea_level` — so an exact-sea
    // plateau would mean something else entirely.
    eprintln!("\n── A3 · where do the near-shore cells actually sit? ──");
    let near: Vec<usize> = (0..n).filter(|&k| dsea[k] <= 3 && dsea[k] != u32::MAX).collect();
    let sea_bits = SEA.to_bits();
    let (mut exact, mut just_below, mut just_above, mut higher) = (0usize, 0usize, 0usize, 0usize);
    for &k in &near {
        let e = m_above_sea(k);
        if f.data[k].to_bits() == sea_bits {
            exact += 1;
        } else if e <= 0.0 {
            just_below += 1;
        } else if e <= 0.01 {
            just_above += 1;
        } else {
            higher += 1;
        }
    }
    let tot = near.len().max(1);
    eprintln!(
        "   over {} land cells within 3 of the sea: EXACTLY sea_level {exact} ({:.2} %) | \
         (sea−0.01, sea] {just_below} ({:.1} %) | (sea, sea+0.01] {just_above} ({:.1} %) | \
         above sea+0.01 {higher} ({:.1} %)",
        tot,
        100.0 * exact as f64 / tot as f64,
        100.0 * just_below as f64 / tot as f64,
        100.0 * just_above as f64 / tot as f64,
        100.0 * higher as f64 / tot as f64
    );
    // and how many mouth cells sit exactly on their receiver's height — the clamp's signature
    let mut on_recv = 0usize;
    for &mk in &mouths {
        if let Some(r0) = recv(mk as usize) {
            if f.data[mk as usize].to_bits() == f.data[r0].to_bits() {
                on_recv += 1;
            }
        }
    }
    eprintln!(
        "   of {} mouth cells, {on_recv} ({:.1} %) sit EXACTLY on their receiver's height — the \
         `clamp(hr, ho)` signature of the drag",
        mouths.len(),
        100.0 * on_recv as f64 / mouths.len().max(1) as f64
    );

    // ── A4 · STOP RULE 1 DEMANDS AN ANSWER: what ARE the spurs, then? ────────
    // The flat hypothesis is dead, so the excursions are characterised directly. Three
    // questions, each with a yes/no answer: does the excursion enclose LAND or SEA (a tooth or
    // a bay); is its tip a CHANNEL (accumulation above A_c) or hillslope; and how deep is the
    // notch between the neck and the tip.
    eprintln!("\n── A4 · what the spurs ARE (stop rule 1's second half) ──");
    let acc = &flow.accumulation.data;
    let mut coast_acc: Vec<f32> = Vec::new();
    for pl in &polys {
        for &(px, py) in pl.iter() {
            let (x, y) = (px as usize, py as usize);
            if x < w && y < h && is_land[y * w + x] {
                coast_acc.push(acc[y * w + x]);
            }
        }
    }
    coast_acc.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (mut encloses_land, mut encloses_sea) = (0usize, 0usize);
    let (mut tip_acc, mut notch_m) = (Vec::new(), Vec::new());
    let mut tip_is_channel = 0usize;
    // RELIEF_V1_A_C_KM2 / cell_km2 — the channel head the incision itself uses
    let a_c_cells = 0.1f32 / (CELL_KM * CELL_KM);
    for sp in &spurs {
        let pl = &polys[sp.poly];
        let (mut cx, mut cy) = (0.0f64, 0.0f64);
        for p in &pl[sp.i..=sp.j] {
            cx += p.0 as f64;
            cy += p.1 as f64;
        }
        let m = (sp.j - sp.i + 1) as f64;
        let (cix, ciy) = ((cx / m) as usize, (cy / m) as usize);
        if cix < w && ciy < h {
            if is_land[ciy * w + cix] {
                encloses_land += 1;
            } else {
                encloses_sea += 1;
            }
        }
        let mut far = (0.0f32, pl[sp.i]);
        for p in &pl[sp.i..=sp.j] {
            let d = ((p.0 - sp.mid.0).powi(2) + (p.1 - sp.mid.1).powi(2)).sqrt();
            if d > far.0 {
                far = (d, *p);
            }
        }
        let (tx, ty) = (far.1.0 as usize, far.1.1 as usize);
        let (nx, ny) = (sp.mid.0 as usize, sp.mid.1 as usize);
        if tx < w && ty < h && nx < w && ny < h {
            let a = acc[ty * w + tx];
            tip_acc.push(a);
            if a >= a_c_cells {
                tip_is_channel += 1;
            }
            notch_m.push(m_above_sea(ny * w + nx) - m_above_sea(ty * w + tx));
        }
    }
    tip_acc.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    notch_m.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!(
        "   the excursion encloses LAND for {encloses_land} spurs ({:.1} %) and SEA for \
         {encloses_sea} ({:.1} %) ⇒ they are {}",
        100.0 * encloses_land as f64 / spurs.len() as f64,
        100.0 * encloses_sea as f64 / spurs.len() as f64,
        if encloses_land > encloses_sea {
            "TEETH (land slivers pointing seaward)"
        } else {
            "BAYS (sea intrusions)"
        }
    );
    eprintln!(
        "   accumulation at the TIP (cells): p10 {:.0} MEDIAN {:.0} p90 {:.0} | coastline-wide \
         median {:.0} | tips at or above A_c = {a_c_cells:.1} cells: **{tip_is_channel} of {} \
         ({:.1} %)**",
        pct(&tip_acc, 0.10),
        pct(&tip_acc, 0.50),
        pct(&tip_acc, 0.90),
        pct(&coast_acc, 0.50),
        tip_acc.len(),
        100.0 * tip_is_channel as f64 / tip_acc.len().max(1) as f64
    );
    eprintln!(
        "   NOTCH depth (neck − tip, m): p10 {:.2} MEDIAN {:.2} p90 {:.2} — negative means the \
         tip is HIGHER than the neck",
        pct(&notch_m, 0.10),
        pct(&notch_m, 0.50),
        pct(&notch_m, 0.90)
    );
    // ── A5 · FBM against INCISION, at the CELL scale ─────────────────────────
    // Finding 51 concluded the chain is fringe-free up to the incision's input — measured at
    // the KILOMETRE scale, the one that cannot see a 2-to-5-cell tooth. This is the same
    // separation read at the scale the fur actually lives at: FBM and both closures ON, only
    // `stream_power` removed. (Finding 75's REFERENCE drops the FBM too, so it cannot
    // attribute between them.)
    let no_inc = build_field(false);
    let polys_ni = marching_squares(&no_inc, SEA);
    eprintln!("\n── A5 · is the fur the FBM's or the INCISION's? (cell scale) ──");
    eprintln!(
        "   {:<34} {:>10} {:>10} {:>12} {:>9} {:>8}",
        "field", "≥1 km", "≥2 cells", "coast km", "p90 km", "R cell"
    );
    for (nm, pp) in
        [("SHIPPED (FBM + incision)", &polys), ("FBM + closures, NO incision", &polys_ni)]
    {
        let a = coast_shape_thresholds(pp, CELL_KM, 1.0, NECK_KM);
        let b = coast_shape_thresholds(pp, CELL_KM, 2.0 * CELL_KM, CELL_KM);
        eprintln!(
            "   {nm:<34} {:>10} {:>10} {:>12.0} {:>9.2} {:>8.3}",
            a.count, b.count, a.coast_km, a.p90_len_km, b.local_axis_r
        );
    }

    // ── C · the spectrum of spur POSITIONS along the coast ───────────────────
    // Autocorrelation surrogate: the distribution of GAPS between consecutive spur roots, per
    // continuous coast segment > 20 km, against a WHITE baseline — the same number of spurs
    // with positions drawn uniformly over the same length, averaged over DRAWS draws. A real
    // coast has no wavelength; a threshold artefact does.
    const DRAWS: usize = 64;
    let spectrum = |pp: &[Vec<(f32, f32)>], label: &str| {
        let (sp2, _) = coast_spurs(pp, CELL_KM, 2.0 * CELL_KM, CELL_KM);
        let mut by_poly: std::collections::HashMap<usize, Vec<f32>> = Default::default();
        for s in &sp2 {
            by_poly.entry(s.poly).or_default().push(s.root_km);
        }
        let mut gaps: Vec<f32> = Vec::new();
        let mut used = 0usize;
        for (pi, mut roots) in by_poly {
            let pl = &pp[pi];
            let mut len = 0.0f32;
            for i in 1..pl.len() {
                len += ((pl[i].0 - pl[i - 1].0).powi(2) + (pl[i].1 - pl[i - 1].1).powi(2)).sqrt();
            }
            if len * CELL_KM < 20.0 || roots.len() < 8 {
                continue;
            }
            used += 1;
            roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            for wdw in roots.windows(2) {
                gaps.push((wdw[1] - wdw[0]) / CELL_KM); // in CELLS
            }
        }
        if gaps.len() < 32 {
            eprintln!("   {label:<28} population {} gaps — rule 10: not a reading", gaps.len());
            return;
        }
        let mut hist = [0usize; 33];
        for &g in &gaps {
            hist[(g.round() as usize).min(32)] += 1;
        }
        let mut base = [0.0f64; 33];
        let mut st = 0x2545_F491_4F6C_DD1Du64;
        let mut rnd = || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            (st >> 11) as f64 / (1u64 << 53) as f64
        };
        let total_len_cells: f32 = gaps.iter().sum();
        for _ in 0..DRAWS {
            let mut u: Vec<f32> =
                (0..gaps.len() + 1).map(|_| (rnd() as f32) * total_len_cells).collect();
            u.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            for wdw in u.windows(2) {
                base[((wdw[1] - wdw[0]).round() as usize).min(32)] += 1.0;
            }
        }
        for b in base.iter_mut() {
            *b /= DRAWS as f64;
        }
        let (mut peak, mut peak_ratio) = (0usize, 0.0f64);
        for b in 1..=32 {
            let r = hist[b] as f64 / base[b].max(0.5);
            if r > peak_ratio {
                peak_ratio = r;
                peak = b;
            }
        }
        let med = {
            let mut g = gaps.clone();
            g.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            g[g.len() / 2]
        };
        eprintln!(
            "   {label:<28} {used} segments > 20 km, {} gaps, median gap {med:.1} cells | \
             dominant spacing **{peak} cells** at **×{peak_ratio:.2}** the white baseline {}",
            gaps.len(),
            if peak_ratio >= 3.0 {
                "⇒ ABOVE 3× — a wavelength exists"
            } else {
                "⇒ below 3×"
            }
        );
    };
    eprintln!("\n── C · spectrum of spur positions (white baseline, {DRAWS} draws) ──");
    spectrum(&polys, "SHIPPED");
    spectrum(&polys_ni, "no incision");
    eprintln!(
        "\n   ⇒ STOP RULE 1: {}",
        if r.0 >= 0.5 {
            "r ≥ 0.50 — the mechanism stands, the seam may be written."
        } else {
            "r < 0.50 — THE MECHANISM IS FALSE. The seam is NOT written."
        }
    );
}
