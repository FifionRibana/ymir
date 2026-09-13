//! ADR Finding 77 — naming the fur's mechanism: the HOLLOWS, and three sweeps.
//!
//! Finding 76 measured the wrong population and said so is not enough: block A4 took, as a
//! spur's "tip", the polyline sample farthest from the neck midpoint — and since 81.4 % of
//! excursions enclose LAND, that sample is the tip of a land SLIVER, i.e. a ridge, whose
//! accumulation is 1 by construction. **"86 % carry no channel" measured the ridges.** The
//! cells that MAKE the fringe are the hollows between the slivers: ex-land (land in the
//! pre-incision field, sea in the shipped one). No block of Finding 76 read them.
//!
//! Every quantity below carries its population: LAND / SEA / EX-LAND.
//!
//! Run: cargo test -p ymir-core --release --test coastal_attribution -- --ignored --nocapture

mod common;

use common::{A_C_CELLS, CELL_KM, DOMAIN_KM, Knobs, SEA, build_field, pct, sorted, spectrum};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{FlowConfig, compute_flow};

/// One hollow: an 8-connected component of EX-LAND cells.
struct Hollow {
    cells: usize,
    /// deepest point of the hollow in the shipped field, metres below sea (negative)
    depth_m: f32,
    /// how far that same cell stood ABOVE sea before the incision, metres
    pre_m: f32,
    /// pre − shipped at the deepest cell: what the incision actually removed, metres
    cut_m: f32,
    /// bounding box, cells
    long_cells: usize,
    short_cells: usize,
    /// the HEAD: the ex-land cell standing highest in the PRE-incision field — the inland end
    head_acc_cells: f32,
    /// the largest accumulation found anywhere in the hollow (a bank cell sits BESIDE a
    /// channel, so its own A is small while the channel's is not)
    max_acc_cells: f32,
    /// the largest accumulation in the 8-neighbourhood of the hollow — the channel it flanks
    neighbour_acc_cells: f32,
}

fn hollows(pre: &GridF32, shipped: &GridF32, acc: &[f32], norm_to_m: f32) -> Vec<Hollow> {
    let (w, h) = (shipped.width, shipped.height);
    let n = w * h;
    // THE MASK, declared: pre-incision LAND ∧ shipped SEA. Chosen over "the isoline segment
    // between two consecutive spur roots" because it needs no chord geometry and cannot
    // mis-pair; its cost is that it also catches ex-land far from any spur, which is measured
    // and reported rather than assumed away.
    let ex: Vec<bool> = (0..n).map(|k| pre.data[k] > SEA && shipped.data[k] <= SEA).collect();
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    let mut stack: Vec<u32> = Vec::new();
    for s0 in 0..n {
        if !ex[s0] || seen[s0] {
            continue;
        }
        seen[s0] = true;
        stack.push(s0 as u32);
        let (mut cells, mut comp) = (0usize, Vec::new());
        let (mut x0, mut y0, mut x1, mut y1) = (usize::MAX, usize::MAX, 0usize, 0usize);
        while let Some(kk) = stack.pop() {
            let k = kk as usize;
            comp.push(k);
            cells += 1;
            let (x, y) = (k % w, k / w);
            x0 = x0.min(x);
            y0 = y0.min(y);
            x1 = x1.max(x);
            y1 = y1.max(y);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if ex[nk] && !seen[nk] {
                        seen[nk] = true;
                        stack.push(nk as u32);
                    }
                }
            }
        }
        // deepest cell (shipped), head cell (highest PRE), accumulations
        let (mut deep_k, mut deep_v) = (comp[0], f32::MAX);
        let (mut head_k, mut head_v) = (comp[0], f32::MIN);
        let mut max_acc = 0.0f32;
        let mut nb_acc = 0.0f32;
        for &k in &comp {
            if shipped.data[k] < deep_v {
                deep_v = shipped.data[k];
                deep_k = k;
            }
            if pre.data[k] > head_v {
                head_v = pre.data[k];
                head_k = k;
            }
            max_acc = max_acc.max(acc[k]);
            let (x, y) = (k % w, k / w);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if !ex[nk] {
                        nb_acc = nb_acc.max(acc[nk]);
                    }
                }
            }
        }
        let (bw, bh) = (x1 - x0 + 1, y1 - y0 + 1);
        out.push(Hollow {
            cells,
            depth_m: (shipped.data[deep_k] - SEA) * norm_to_m,
            pre_m: (pre.data[deep_k] - SEA) * norm_to_m,
            cut_m: (pre.data[deep_k] - shipped.data[deep_k]) * norm_to_m,
            long_cells: bw.max(bh),
            short_cells: bw.min(bh),
            head_acc_cells: acc[head_k],
            max_acc_cells: max_acc,
            neighbour_acc_cells: nb_acc,
        });
    }
    out
}

fn six(name: &str, polys: &[Vec<(f32, f32)>]) -> (usize, usize, f32, f32) {
    let a = coast_shape_thresholds(polys, CELL_KM, 1.0, NECK_KM);
    let b = coast_shape_thresholds(polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    eprintln!(
        "   {name:<30} {:>8} {:>10} {:>11.0} {:>9.2} {:>8.3}",
        a.count, b.count, a.coast_km, a.p90_len_km, b.local_axis_r
    );
    (a.count, b.count, a.coast_km, a.p90_len_km)
}

#[test]
#[ignore]
fn coastal_attribution() {
    let ss =
        ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams::default();
    let norm_to_m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 77 · the hollows, and three sweeps  ==========");
    eprintln!(
        "   1 norm unit = {norm_to_m:.0} m · cell = {:.1} m · A_c = {A_C_CELLS:.2} cells",
        CELL_KM * 1000.0
    );

    let shipped = build_field(Knobs::shipped());
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (shipped.width, shipped.height);
    let flow = compute_flow(&shipped, &FlowConfig { sea_level: SEA, ..Default::default() });

    // ── A · the hollows ───────────────────────────────────────────────────────
    let hs = hollows(&pre, &shipped, &flow.accumulation.data, norm_to_m);
    eprintln!("\n── A · the HOLLOWS (population EX-LAND: pre-incision land ∧ shipped sea) ──");
    eprintln!("   mask: pre > SEA ∧ shipped ≤ SEA, 8-connected components");
    assert!(!hs.is_empty(), "rule 10: the mask found no hollow at all");
    let dep = sorted(hs.iter().map(|x| x.depth_m).collect());
    let cut = sorted(hs.iter().map(|x| x.cut_m).collect());
    let lng = sorted(hs.iter().map(|x| x.long_cells as f32).collect());
    let sht = sorted(hs.iter().map(|x| x.short_cells as f32).collect());
    let cel = sorted(hs.iter().map(|x| x.cells as f32).collect());
    let hac = sorted(hs.iter().map(|x| x.head_acc_cells).collect());
    let mac = sorted(hs.iter().map(|x| x.max_acc_cells).collect());
    let nac = sorted(hs.iter().map(|x| x.neighbour_acc_cells).collect());
    eprintln!("   population **{}** hollows", hs.len());
    eprintln!("   {:<38} {:>10} {:>10} {:>10}", "column", "p10", "MEDIAN", "p90");
    for (nm, v) in [
        ("depth below sea, m (SHIPPED)", &dep),
        ("the same cell BEFORE incision, m", &sorted(hs.iter().map(|x| x.pre_m).collect())),
        ("INCISED amount pre−shipped, m", &cut),
        ("long side, cells", &lng),
        ("short side, cells", &sht),
        ("area, cells", &cel),
        ("HEAD accumulation, cells", &hac),
        ("MAX accumulation inside, cells", &mac),
        ("MAX accumulation of a NEIGHBOUR", &nac),
    ] {
        eprintln!(
            "   {nm:<38} {:>10.3} {:>10.3} {:>10.3}",
            pct(v, 0.10),
            pct(v, 0.50),
            pct(v, 0.90)
        );
    }
    let ge =
        |v: &[f32]| 100.0 * v.iter().filter(|&&x| x >= A_C_CELLS).count() as f64 / v.len() as f64;
    eprintln!(
        "   share at or above A_c = {A_C_CELLS:.1} cells — HEAD **{:.1} %** · MAX inside {:.1} % \
         · best NEIGHBOUR **{:.1} %**",
        ge(&hac),
        ge(&mac),
        ge(&nac)
    );
    // the contingency table the round asked for
    let (mut sh_sh, mut sh_dp, mut dl_sh, mut dl_dp) = (0usize, 0usize, 0usize, 0usize);
    for x in &hs {
        match (x.head_acc_cells >= A_C_CELLS, x.depth_m > -2.0) {
            (true, true) => sh_sh += 1,
            (true, false) => sh_dp += 1,
            (false, true) => dl_sh += 1,
            (false, false) => dl_dp += 1,
        }
    }
    eprintln!(
        "   contingency (head ≥ A_c) × (depth shallower than −2 m):\n     head ≥ A_c & shallow \
         {sh_sh} | head ≥ A_c & deep {sh_dp} | head < A_c & shallow **{dl_sh}** | head < A_c & \
         deep {dl_dp}"
    );
    // NEGATIVE CONTROL 1 — structural, and declared as such: the mask on (pre, pre) is empty by
    // construction, so it only proves the component finder does not invent cells.
    let none = hollows(&pre, &pre, &flow.accumulation.data, norm_to_m);
    eprintln!(
        "   negative control (mask applied to pre vs pre): **{}** hollows — must be 0, and it is \
         structural: no cell can be land and sea in the same field",
        none.len()
    );
    assert!(none.is_empty(), "the component finder invented cells");
    // NEGATIVE CONTROL 2 — the informative one: does the mask leak inland?
    let polys_sh = marching_squares(&shipped, SEA);
    let (sp_sh, _) =
        ymir_core::terrain::coast_metrics::coast_spurs(&polys_sh, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    let mut near = vec![false; w * h];
    for s in &sp_sh {
        let (cx, cy) = (s.mid.0 as i32, s.mid.1 as i32);
        for dy in -24i32..=24 {
            for dx in -24i32..=24 {
                let (nx, ny) = (cx + dx, cy + dy);
                if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                    near[ny as usize * w + nx as usize] = true;
                }
            }
        }
    }
    let ex_total: usize = hs.iter().map(|x| x.cells).sum();
    eprintln!(
        "   mask extent: {ex_total} ex-land cells in {} hollows; {} spurs ≥ 2 cells on SHIPPED",
        hs.len(),
        sp_sh.len()
    );

    // ── the Finding 76 A4 correction, side by side ────────────────────────────
    eprintln!(
        "\n── A4 CORRECTED · tips vs HEADS ──\n   Finding 76 read the RIDGE tips: 13.7 % at or \
         above A_c. Read on the HOLLOW heads instead: **{:.1} %**. Read on the best neighbour \
         of the hollow (the channel a bank would flank): **{:.1} %**.",
        ge(&hac),
        ge(&nac)
    );

    // ── B · λ against A_c ─────────────────────────────────────────────────────
    eprintln!("\n── B · λ against A_c ──");
    eprintln!(
        "   {:<30} {:>8} {:>10} {:>11} {:>9} {:>8}",
        "variant", "≥1 km", "≥2 cells", "coast km", "p90 km", "R cell"
    );
    let base_polys = marching_squares(&pre, SEA);
    six("PRE-INCISION (the authority)", &base_polys);
    six("SHIPPED (A_c = 0.1 km²)", &polys_sh);
    let mut rows: Vec<(String, Option<(usize, usize, f32, usize, f64)>)> = Vec::new();
    rows.push(("PRE-INCISION".into(), spectrum(&base_polys)));
    rows.push(("A_c 0.1 (SHIPPED)".into(), spectrum(&polys_sh)));
    for a in [0.025f32, 0.4, 1.0] {
        let f = build_field(Knobs { a_c_km2: Some(a), ..Knobs::shipped() });
        let p = marching_squares(&f, SEA);
        six(&format!("A_c = {a} km²"), &p);
        rows.push((format!("A_c {a}"), spectrum(&p)));
    }

    // ── C · λ against the diffusion ───────────────────────────────────────────
    eprintln!("\n── C · λ against the diffusion (F59 §4: the anchored value is 1.28) ──");
    eprintln!(
        "   {:<30} {:>8} {:>10} {:>11} {:>9} {:>8}",
        "variant", "≥1 km", "≥2 cells", "coast km", "p90 km", "R cell"
    );
    for d in [0.0f32, 0.32, 1.28] {
        let sub = if d > 1.0 { Some(8) } else { None };
        let f =
            build_field(Knobs { diffusion: Some(d), diffusion_substeps: sub, ..Knobs::shipped() });
        let p = marching_squares(&f, SEA);
        six(&format!("diffusion = {d}{}", if sub.is_some() { " (substeps 8)" } else { "" }), &p);
        rows.push((format!("diffusion {d}"), spectrum(&p)));
    }

    // ── C-bis · MY candidate: the bank-widening closure and the talus ─────────
    // `lateral_erosion` planes the banks PERPENDICULAR to flow out to a physical half-width
    // `K_lat·A^m`, only where `A ≥ A_c` — it lowers cells that are NOT channels, in a comb at
    // right angles to them. `talus_passes` relaxes to the repose angle everywhere. Neither was
    // ever swept; Finding 55 swept the MFD exponent.
    eprintln!("\n── C-bis · the two operators that move NON-channel cells (never swept) ──");
    eprintln!(
        "   {:<30} {:>8} {:>10} {:>11} {:>9} {:>8}",
        "variant", "≥1 km", "≥2 cells", "coast km", "p90 km", "R cell"
    );
    for (nm, kn) in [
        (
            "lateral_erosion = 0 (banks OFF)",
            Knobs { lateral_erosion: Some(0.0), ..Knobs::shipped() },
        ),
        ("talus_passes = 0", Knobs { talus_passes: Some(0), ..Knobs::shipped() }),
    ] {
        let f = build_field(kn);
        let p = marching_squares(&f, SEA);
        six(nm, &p);
        rows.push((nm.to_string(), spectrum(&p)));
    }

    eprintln!("\n── the SPECTRA, one line per variant ──");
    eprintln!(
        "   {:<30} {:>9} {:>8} {:>12} {:>10} {:>12}",
        "variant", "segments", "gaps", "median gap", "λ cells", "× white"
    );
    for (nm, s) in &rows {
        match s {
            Some((seg, g, med, peak, ratio)) => eprintln!(
                "   {nm:<30} {seg:>9} {g:>8} {med:>12.1} {peak:>10} {ratio:>12.2}{}",
                if *ratio >= 3.0 { "  ⇒ wavelength" } else { "" }
            ),
            None => eprintln!("   {nm:<30} {:>9} — rule 10: under 32 gaps, NOT a reading", "—"),
        }
    }

    // ── D · the class proof: fill the shallow hollows, post-loop ──────────────
    eprintln!(
        "\n── D-post · fill ex-land hollows shallower than a threshold (CLASS TEST, not a remedy) ──"
    );
    eprintln!(
        "   ⚠️ ADR line 18 already records what a deposition term did: \"it deposits in channels \
         faster than it incises … the common cause of the terraces, the missing valleys and the \
         COASTAL SEDIMENT DUMP\". This tests the CLASS; it does not propose it."
    );
    eprintln!(
        "   {:<30} {:>8} {:>10} {:>11} {:>9} {:>8}",
        "threshold", "≥1 km", "≥2 cells", "coast km", "p90 km", "R cell"
    );
    let (pre_a, pre_b, pre_km, _) = (
        coast_shape_thresholds(&base_polys, CELL_KM, 1.0, NECK_KM).count,
        coast_shape_thresholds(&base_polys, CELL_KM, 2.0 * CELL_KM, CELL_KM).count,
        coast_shape_thresholds(&base_polys, CELL_KM, 1.0, NECK_KM).coast_km,
        0,
    );
    for t in [0.0f32, 0.5, 1.0, 2.0] {
        let mut g = shipped.clone();
        let mut lifted = 0usize;
        for k in 0..w * h {
            if pre.data[k] > SEA && shipped.data[k] <= SEA {
                let d = -(shipped.data[k] - SEA) * norm_to_m; // metres below sea, positive
                if d <= t {
                    g.data[k] = SEA + 0.01 / norm_to_m;
                    lifted += 1;
                }
            }
        }
        let p = marching_squares(&g, SEA);
        let (a, b, km, _) = six(&format!("fill ≤ {t} m ({lifted} cells)"), &p);
        eprintln!(
            "        Δ against PRE-INCISION: ≥1 km {:+}, ≥2 cells {:+}, coast {:+.0} km \
             ({:+.1} %)",
            a as i64 - pre_a as i64,
            b as i64 - pre_b as i64,
            km - pre_km,
            100.0 * (km - pre_km) / pre_km
        );
    }
    eprintln!(
        "   (threshold 0 m is the negative control: it lifts only cells already at sea and must \
         leave the counts at SHIPPED's.)"
    );
}

/// ADR Finding 77, second pass — **the hollow mask of the first pass LEAKED**, and this is the
/// corrected instrument.
///
/// `pre land AND shipped sea` catches every cell the incision ever took below sea level, which
/// at 8192 is 298 598 cells in 360 components of median 44 and p90 2 660 cells — the DROWNED
/// VALLEYS, not the 1 831 teeth of median 4 cells. So the first pass's "81.4 % of heads at or
/// above A_c" answers a question nobody asked. The objection was about the hollows BETWEEN two
/// adjacent slivers; this restricts the population to them.
///
/// Restriction, declared: a hollow counts when at least half its cells lie within `ZONE` cells
/// of a coastline sample that BELONGS to a spur of 2 cells or more.
#[test]
#[ignore]
fn coastal_hollows_fringe() {
    let ss =
        ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams::default();
    let norm_to_m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    const ZONE: i32 = 3;
    eprintln!("\n==========  Finding 77 second pass · the FRINGE hollows only  ==========");
    let shipped = build_field(Knobs::shipped());
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (shipped.width, shipped.height);
    let n = w * h;
    let flow = compute_flow(&shipped, &FlowConfig { sea_level: SEA, ..Default::default() });
    let acc = &flow.accumulation.data;

    let polys = marching_squares(&shipped, SEA);
    let (sp, _) =
        ymir_core::terrain::coast_metrics::coast_spurs(&polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    let mut zone = vec![false; n];
    for s in &sp {
        for p in &polys[s.poly][s.i..=s.j] {
            let (px, py) = (p.0 as i32, p.1 as i32);
            for dy in -ZONE..=ZONE {
                for dx in -ZONE..=ZONE {
                    let (nx, ny) = (px + dx, py + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        zone[ny as usize * w + nx as usize] = true;
                    }
                }
            }
        }
    }
    eprintln!(
        "   {} spurs >= 2 cells | spur zone = {} cells (within {ZONE} of a spur's own samples)",
        sp.len(),
        zone.iter().filter(|&&z| z).count()
    );

    let ex: Vec<bool> = (0..n).map(|k| pre.data[k] > SEA && shipped.data[k] <= SEA).collect();
    let depths =
        sorted((0..n).filter(|&k| ex[k]).map(|k| (shipped.data[k] - SEA) * norm_to_m).collect());
    let in_zone: Vec<f32> = sorted(
        (0..n).filter(|&k| ex[k] && zone[k]).map(|k| (shipped.data[k] - SEA) * norm_to_m).collect(),
    );
    eprintln!(
        "\n-- the EX-LAND population, per cell, metres below sea --\n   ALL {} cells: p1 {:.3} \
         p10 {:.3} MEDIAN {:.3} p90 {:.3} max {:.3}\n   IN THE SPUR ZONE {} cells: p1 {:.3} p10 \
         {:.3} MEDIAN {:.3} p90 {:.3} max {:.3}",
        depths.len(),
        pct(&depths, 0.01),
        pct(&depths, 0.10),
        pct(&depths, 0.50),
        pct(&depths, 0.90),
        depths.last().copied().unwrap_or(f32::NAN),
        in_zone.len(),
        pct(&in_zone, 0.01),
        pct(&in_zone, 0.10),
        pct(&in_zone, 0.50),
        pct(&in_zone, 0.90),
        in_zone.last().copied().unwrap_or(f32::NAN)
    );
    let sh = |v: &[f32], t: f32| {
        100.0 * v.iter().filter(|&&d| d > -t).count() as f64 / v.len().max(1) as f64
    };
    eprintln!(
        "   share shallower than -0.5 / -1 / -2 / -5 m -- ALL: {:.2} / {:.2} / {:.2} / {:.2} % | \
         IN ZONE: {:.2} / {:.2} / {:.2} / {:.2} %",
        sh(&depths, 0.5),
        sh(&depths, 1.0),
        sh(&depths, 2.0),
        sh(&depths, 5.0),
        sh(&in_zone, 0.5),
        sh(&in_zone, 1.0),
        sh(&in_zone, 2.0),
        sh(&in_zone, 5.0)
    );

    // WHERE DOES -20 m COME FROM? Every ex-land cell reads exactly -20.000, which is a
    // constant, not a distribution. Either the near-shore ocean IS a -20 m shelf (and the
    // incision merely relaxes land onto its receiver, Finding 61) or something clamps. The
    // ocean's own distribution decides it.
    let shore_ocean: Vec<f32> = sorted(
        (0..n)
            .filter(|&k| {
                shipped.data[k] <= SEA && {
                    let (x, y) = ((k % w) as i32, (k / w) as i32);
                    (-1i32..=1).any(|dy| {
                        (-1i32..=1).any(|dx| {
                            let (nx, ny) = (x + dx, y + dy);
                            nx >= 0
                                && ny >= 0
                                && (nx as usize) < w
                                && (ny as usize) < h
                                && shipped.data[ny as usize * w + nx as usize] > SEA
                        })
                    })
                }
            })
            .map(|k| (shipped.data[k] - SEA) * norm_to_m)
            .collect(),
    );
    let all_ocean = sorted(
        (0..n)
            .filter(|&k| shipped.data[k] <= SEA)
            .map(|k| (shipped.data[k] - SEA) * norm_to_m)
            .collect(),
    );
    let pre_shore: Vec<f32> = sorted(
        (0..n).filter(|&k| pre.data[k] <= SEA).map(|k| (pre.data[k] - SEA) * norm_to_m).collect(),
    );
    eprintln!(
        "
-- where -20 m comes from --
   SHIPPED ocean cells ADJACENT TO LAND, {} of them:          p1 {:.3} p10 {:.3} MEDIAN {:.3} p90 {:.3} max {:.3} m
   SHIPPED all {} ocean cells:          p1 {:.1} MEDIAN {:.1} max {:.3} m
   PRE-INCISION all {} ocean cells: p1 {:.1} MEDIAN          {:.1} max {:.3} m",
        shore_ocean.len(),
        pct(&shore_ocean, 0.01),
        pct(&shore_ocean, 0.10),
        pct(&shore_ocean, 0.50),
        pct(&shore_ocean, 0.90),
        shore_ocean.last().copied().unwrap_or(f32::NAN),
        all_ocean.len(),
        pct(&all_ocean, 0.01),
        pct(&all_ocean, 0.50),
        all_ocean.last().copied().unwrap_or(f32::NAN),
        pre_shore.len(),
        pct(&pre_shore, 0.01),
        pct(&pre_shore, 0.50),
        pre_shore.last().copied().unwrap_or(f32::NAN)
    );

    let mut seen = vec![false; n];
    let mut kept: Vec<(usize, f32, f32, f32, f32)> = Vec::new();
    let mut total = 0usize;
    for s0 in 0..n {
        if !ex[s0] || seen[s0] {
            continue;
        }
        total += 1;
        let mut comp = Vec::new();
        let mut stack = vec![s0 as u32];
        seen[s0] = true;
        while let Some(kk) = stack.pop() {
            let k = kk as usize;
            comp.push(k);
            let (x, y) = (k % w, k / w);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if ex[nk] && !seen[nk] {
                        seen[nk] = true;
                        stack.push(nk as u32);
                    }
                }
            }
        }
        if comp.iter().filter(|&&k| zone[k]).count() * 2 < comp.len() {
            continue;
        }
        let (mut deep, mut head_k, mut head_v, mut nbr) = (f32::MAX, comp[0], f32::MIN, 0.0f32);
        for &k in &comp {
            deep = deep.min((shipped.data[k] - SEA) * norm_to_m);
            if pre.data[k] > head_v {
                head_v = pre.data[k];
                head_k = k;
            }
            let (x, y) = (k % w, k / w);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if !ex[nk] {
                            nbr = nbr.max(acc[nk]);
                        }
                    }
                }
            }
        }
        kept.push((
            comp.len(),
            deep,
            acc[head_k],
            nbr,
            (pre.data[head_k] - shipped.data[head_k]) * norm_to_m,
        ));
    }
    eprintln!(
        "\n-- FRINGE hollows: **{}** of {} components kept (at least half their cells in the \
         spur zone) --",
        kept.len(),
        total
    );
    assert!(!kept.is_empty(), "rule 10: the restriction emptied the population");
    let cells = sorted(kept.iter().map(|x| x.0 as f32).collect());
    let deep = sorted(kept.iter().map(|x| x.1).collect());
    let head = sorted(kept.iter().map(|x| x.2).collect());
    let nbr = sorted(kept.iter().map(|x| x.3).collect());
    let cut = sorted(kept.iter().map(|x| x.4).collect());
    eprintln!("   {:<36} {:>10} {:>10} {:>10}", "column", "p10", "MEDIAN", "p90");
    for (nm, v) in [
        ("area, cells", &cells),
        ("depth below sea, m", &deep),
        ("incised at the head, m", &cut),
        ("HEAD accumulation, cells", &head),
        ("best NEIGHBOUR accumulation, cells", &nbr),
    ] {
        eprintln!(
            "   {nm:<36} {:>10.2} {:>10.2} {:>10.2}",
            pct(v, 0.10),
            pct(v, 0.50),
            pct(v, 0.90)
        );
    }
    let ge =
        |v: &[f32]| 100.0 * v.iter().filter(|&&x| x >= A_C_CELLS).count() as f64 / v.len() as f64;
    eprintln!(
        "   at or above A_c = {A_C_CELLS:.1} cells -- HEAD **{:.1} %** - best NEIGHBOUR **{:.1} \
         %**\n   (Finding 76 read the ridge TIPS: 13.7 %. The first pass read the drowned \
         valleys: 81.4 %. THIS is the number the objection asked for.)",
        ge(&head),
        ge(&nbr)
    );
    let _ = DOMAIN_KM;
}
