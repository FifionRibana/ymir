//! ADR Finding 78 — the 20 m shelf clamp: its intent, its price, and what it is NOT.
//!
//! `BathymetryProfile::shelf_min_depth_m = 20.0` (`bathymetry.rs:49`) clamps every ocean cell to
//! at least 20 m of depth (`:178`, with `min_depth = shelf_min_depth_m.max(1.0)` at `:171`).
//!
//! ORDERING, corrected against the round's premise: `apply_bathymetry_profile` is called at
//! `production_upscale.rs:506` and `incise_lithology` at `:451`. **The clamp runs AFTER the
//! incision**, so it cannot be what the relaxation landed on — it overwrites the incision's
//! output, on sub-sea cells only. The `bathymetry_off` knob produces the field the incision
//! actually made, which is the control that separates the two.
//!
//! Run: cargo test -p ymir-core --release --test coastal_shelf -- --ignored --nocapture

mod common;

use common::{A_C_CELLS, CELL_KM, Knobs, SEA, build_field, pct, sorted, spectrum};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};

struct Row {
    name: String,
    spurs_km: usize,
    spurs_cell: usize,
    coast_km: f32,
    p90: f32,
    r: f32,
}

fn six(name: &str, polys: &[Vec<(f32, f32)>]) -> Row {
    let a = coast_shape_thresholds(polys, CELL_KM, 1.0, NECK_KM);
    let b = coast_shape_thresholds(polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    Row {
        name: name.to_string(),
        spurs_km: a.count,
        spurs_cell: b.count,
        coast_km: a.coast_km,
        p90: a.p90_len_km,
        r: b.local_axis_r,
    }
}

/// THE CONTROL BLOCK, in full. Every variant, every column, so "the terrain moved" is a number
/// and not an impression.
fn control(tag: &str, f: &GridF32, pre: &GridF32, shipped: &GridF32, n2m: f32) {
    let n = f.data.len();
    let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
    let n_land = land.iter().filter(|&&l| l).count();
    // land/sea mask movement against SHIPPED — the bathymetry docstring's own invariant
    let moved = (0..n).filter(|&k| (f.data[k] > SEA) != (shipped.data[k] > SEA)).count();
    let h = sorted((0..n).filter(|&k| land[k]).map(|k| (f.data[k] - SEA) * n2m).collect());
    let mean = h.iter().map(|&x| x as f64).sum::<f64>() / h.len().max(1) as f64;
    let flow = compute_flow(f, &FlowConfig { sea_level: SEA, ..Default::default() });
    let pits = (0..n).filter(|&k| flow.filled.data[k] > f.data[k] + 1e-7).count();
    let chan = (0..n).filter(|&k| land[k] && flow.accumulation.data[k] >= A_C_CELLS).count();
    // erosion work: pre-incision minus this field, paired over cells that are land in BOTH
    let (mut wsum, mut wn) = (0.0f64, 0usize);
    for k in 0..n {
        if pre.data[k] > SEA && land[k] {
            wsum += ((pre.data[k] - f.data[k]) * n2m) as f64;
            wn += 1;
        }
    }
    eprintln!(
        "      CONTROL {tag:<26} land {n_land} ({:.4} %) | mask moved vs SHIPPED **{moved}** | \
         hypso mean {mean:.2} m p50 {:.2} | pits {pits} | channel share {:.3} % | work {:.2} m",
        100.0 * n_land as f64 / n as f64,
        pct(&h, 0.50),
        100.0 * chan as f64 / n_land.max(1) as f64,
        wsum / wn.max(1) as f64
    );
}

/// Coastal step and bank slope — the consumer numbers. Population: LAND cells 8-adjacent to sea.
fn banks(tag: &str, f: &GridF32, n2m: f32) {
    let (w, h) = (f.width, f.height);
    let n = w * h;
    let mut step: Vec<f32> = Vec::new();
    let cell_m = CELL_KM * 1000.0;
    for k in 0..n {
        if f.data[k] <= SEA {
            continue;
        }
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let (mut lo, mut found) = (f32::MAX, false);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if f.data[nk] <= SEA {
                    lo = lo.min(f.data[nk]);
                    found = true;
                }
            }
        }
        if found {
            step.push((f.data[k] - lo) * n2m);
        }
    }
    let step = sorted(step);
    let slope = sorted(step.iter().map(|&s| s / cell_m).collect());
    let below =
        |t: f32| 100.0 * slope.iter().filter(|&&v| v < t).count() as f64 / slope.len() as f64;
    let above =
        |t: f32| 100.0 * slope.iter().filter(|&&v| v > t).count() as f64 / slope.len() as f64;
    eprintln!(
        "      BANKS   {tag:<26} {} shore-land cells | STEP land-sea p10 {:.2} p50 **{:.2}** p90 \
         {:.2} m | local slope under 15 % {:.1} % · over 25 % {:.1} %",
        step.len(),
        pct(&step, 0.10),
        pct(&step, 0.50),
        pct(&step, 0.90),
        below(0.15),
        above(0.25)
    );
}

/// Ex-land hollow depth, the enlarged population (pre-land AND delivered-sea), kept separate
/// from Finding 77's restricted 46 as the round asked.
fn hollow_depth(tag: &str, pre: &GridF32, f: &GridF32, n2m: f32) {
    let n = f.data.len();
    let d = sorted(
        (0..n)
            .filter(|&k| pre.data[k] > SEA && f.data[k] <= SEA)
            .map(|k| (f.data[k] - SEA) * n2m)
            .collect(),
    );
    if d.is_empty() {
        eprintln!("      HOLLOWS {tag:<26} EMPTY — rule 10, not a reading");
        return;
    }
    eprintln!(
        "      HOLLOWS {tag:<26} {} ex-land cells | depth p10 {:.3} p50 **{:.3}** p90 {:.3} max \
         {:.3} m",
        d.len(),
        pct(&d, 0.10),
        pct(&d, 0.50),
        pct(&d, 0.90),
        d.last().copied().unwrap_or(f32::NAN)
    );
}

#[test]
#[ignore]
fn coastal_shelf() {
    let ss =
        ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 78 · the 20 m shelf clamp  ==========");
    eprintln!(
        "   ORDERING: incise_lithology at production_upscale.rs:451, apply_bathymetry_profile at \
         :506 — the clamp runs AFTER the incision."
    );

    let shipped = build_field(Knobs::shipped());
    let pre = build_field(Knobs::no_incision());
    let nobath = build_field(Knobs { bathymetry_off: true, ..Knobs::shipped() });
    let (w, h) = (shipped.width, shipped.height);
    let n = w * h;

    // ── A2 · what do the coastal channel cells relax ONTO? ────────────────────
    eprintln!(
        "\n── A2 · receiver height of COASTAL CHANNEL cells (A >= A_c, receiver is ocean) ──"
    );
    for (tag, f) in [("SHIPPED (post-clamp)", &shipped), ("pre-bathymetry", &nobath)] {
        let flow = compute_flow(f, &FlowConfig { sea_level: SEA, ..Default::default() });
        let mut recv_h: Vec<f32> = Vec::new();
        for k in 0..n {
            if f.data[k] <= SEA || flow.accumulation.data[k] < A_C_CELLS {
                continue;
            }
            let d = flow.direction[k];
            if d == DIR_NONE {
                continue;
            }
            let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
            let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
            let r = ny * w + nx;
            if f.data[r] <= SEA {
                recv_h.push((f.data[r] - SEA) * n2m);
            }
        }
        let v = sorted(recv_h);
        let distinct = {
            let mut u: Vec<u32> = v.iter().map(|x| x.to_bits()).collect();
            u.dedup();
            u.len()
        };
        eprintln!(
            "   {tag:<22} {} cells | p1 {:.3} p10 {:.3} MEDIAN {:.3} p90 {:.3} max {:.3} m | \
             **{distinct} distinct values**",
            v.len(),
            pct(&v, 0.01),
            pct(&v, 0.10),
            pct(&v, 0.50),
            pct(&v, 0.90),
            v.last().copied().unwrap_or(f32::NAN)
        );
    }

    // ── B · the shelf sweep ───────────────────────────────────────────────────
    let pre_polys = marching_squares(&pre, SEA);
    let mut rows: Vec<(Row, Option<(usize, usize, f32, usize, f64)>)> =
        vec![(six("PRE-INCISION (authority)", &pre_polys), spectrum(&pre_polys))];
    eprintln!("\n── B · shelf_min_depth_m sweep ──");
    eprintln!(
        "   `bathymetry.rs:171` takes `.max(1.0)`, so the 0.5 m point is CLAMPED TO 1 m by \
         production code. It is swept anyway and reported as what it is: a duplicate of 1 m."
    );
    for (nm, kn, reuse) in [
        ("shelf 20 m (SHIPPED)", Knobs::shipped(), true),
        ("shelf 5 m", Knobs { shelf_min_depth_m: Some(5.0), ..Knobs::shipped() }, false),
        ("shelf 1 m", Knobs { shelf_min_depth_m: Some(1.0), ..Knobs::shipped() }, false),
        (
            "shelf 0.5 m (->1 by .max)",
            Knobs { shelf_min_depth_m: Some(0.5), ..Knobs::shipped() },
            false,
        ),
        (
            "shelf 200 m (NEG CONTROL)",
            Knobs { shelf_min_depth_m: Some(200.0), ..Knobs::shipped() },
            false,
        ),
    ] {
        let f = if reuse { shipped.clone() } else { build_field(kn) };
        let p = marching_squares(&f, SEA);
        let r = six(nm, &p);
        eprintln!(
            "   {:<28} >=1km {:>6} >=2c {:>6} coast {:>8.0} p90 {:>6.2} R {:>6.3}",
            r.name, r.spurs_km, r.spurs_cell, r.coast_km, r.p90, r.r
        );
        hollow_depth(nm, &pre, &f, n2m);
        banks(nm, &f, n2m);
        control(nm, &f, &pre, &shipped, n2m);
        rows.push((r, spectrum(&p)));
    }

    // ── C · the anchored damper on a low floor ────────────────────────────────
    eprintln!("\n── C · diffusion 1.28 (the F59 section 4 anchored value) on shelf 1 m ──");
    for (nm, kn) in [(
        "shelf 1 + diff 1.28",
        Knobs {
            shelf_min_depth_m: Some(1.0),
            diffusion: Some(1.28),
            diffusion_substeps: Some(8),
            ..Knobs::shipped()
        },
    )] {
        let f = build_field(kn);
        let p = marching_squares(&f, SEA);
        let r = six(nm, &p);
        eprintln!(
            "   {:<28} >=1km {:>6} >=2c {:>6} coast {:>8.0} p90 {:>6.2} R {:>6.3}",
            r.name, r.spurs_km, r.spurs_cell, r.coast_km, r.p90, r.r
        );
        hollow_depth(nm, &pre, &f, n2m);
        banks(nm, &f, n2m);
        control(nm, &f, &pre, &shipped, n2m);
        rows.push((r, spectrum(&p)));

        // ── D · the deposition CLASS test, now that the hollows are shallow ───
        eprintln!("\n── D-post · fill ex-land hollows shallower than a threshold (CLASS TEST) ──");
        for t in [0.0f32, 0.5, 1.0, 2.0] {
            let mut g = f.clone();
            let mut lifted = 0usize;
            for k in 0..n {
                if pre.data[k] > SEA && f.data[k] <= SEA {
                    let d = -(f.data[k] - SEA) * n2m;
                    if d <= t {
                        g.data[k] = SEA + 0.01 / n2m;
                        lifted += 1;
                    }
                }
            }
            let gp = marching_squares(&g, SEA);
            let gr = six(&format!("fill <= {t} m ({lifted} cells)"), &gp);
            eprintln!(
                "   {:<28} >=1km {:>6} >=2c {:>6} coast {:>8.0}",
                gr.name, gr.spurs_km, gr.spurs_cell, gr.coast_km
            );
            rows.push((gr, spectrum(&gp)));
        }
    }

    eprintln!("\n── Δ AGAINST THE PRE-INCISION AUTHORITY (the criterion) ──");
    let (b_km, b_cell, b_coast) = (rows[0].0.spurs_km, rows[0].0.spurs_cell, rows[0].0.coast_km);
    eprintln!(
        "   {:<30} {:>10} {:>12} {:>14} {:>9} {:>9}",
        "variant", "d >=1 km", "d >=2 cells", "d coast km", "lambda", "x white"
    );
    for (r, s) in &rows {
        let (l, x) = match s {
            Some((_, _, _, p, ra)) => (format!("{p}"), format!("{ra:.2}")),
            None => ("-".into(), "rule 10".into()),
        };
        eprintln!(
            "   {:<30} {:>+10} {:>+12} {:>+14.0} {:>9} {:>9}",
            r.name,
            r.spurs_km as i64 - b_km as i64,
            r.spurs_cell as i64 - b_cell as i64,
            r.coast_km - b_coast,
            l,
            x
        );
    }
}
