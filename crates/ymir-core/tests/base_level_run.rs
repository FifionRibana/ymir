//! ADR Finding 80 blocks B0–B5 — the base-level bound at 8192², on five masks.
//!
//! Every quantity carries its POPULATION and its STAGE (pre-incision / delivered / bounded) and
//! its GRID (f32 8192 / u16 8192 / terrain 2048 / ocean 1024).
//!
//! Run: cargo test -p ymir-core --release --test base_level_run -- --ignored --nocapture

mod common;

use common::{A_C_CELLS, CELL_KM, Knobs, SEA, build_field, pct, sorted};
use ymir_core::export::height::metric_height_u16;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{FlowConfig, compute_flow};

const EPS: f32 = 0.5;

/// Land/sea as a field the tracer can read. Tracing at 0.5 gives the boundary staircase.
fn to_mask(land: &[bool], w: usize, h: usize) -> GridF32 {
    let mut m = GridF32::new(w, h, 0.0);
    for k in 0..w * h {
        if land[k] {
            m.data[k] = 1.0;
        }
    }
    m
}

/// The DELIVERED land mask: encode to `height.u16` with the production encoder, decode, and
/// threshold at 0 m. This is the raster Living Landz actually reads.
fn land_u16(f: &GridF32, ss: &SteinSteinParams) -> (Vec<bool>, f32) {
    let hl = metric_height_u16(f, ss);
    let step = (hl.max_m - hl.min_m) / 65535.0;
    let land = hl
        .codes
        .iter()
        .map(|&c| hl.min_m + (c as f32 / 65535.0) * (hl.max_m - hl.min_m) > 0.0)
        .collect();
    (land, step)
}

/// Binary majority downsample by `f`, the author's declared consumer method: threshold at 8192²
/// then resample the MASK, not the heights. Ties (exactly half) go to SEA — declared, because a
/// tie rule that went to land would grow every sliver.
fn majority(land: &[bool], w: usize, f: usize) -> (Vec<bool>, usize) {
    let nw = w / f;
    let mut out = vec![false; nw * nw];
    let half = f * f / 2;
    for y in 0..nw {
        for x in 0..nw {
            let mut c = 0usize;
            for dy in 0..f {
                for dx in 0..f {
                    if land[(y * f + dy) * w + x * f + dx] {
                        c += 1;
                    }
                }
            }
            out[y * nw + x] = c > half;
        }
    }
    (out, nw)
}

struct Six {
    km: usize,
    cell: usize,
    coast: f32,
    p90: f32,
    r: f32,
}

fn six(land: &[bool], w: usize, cell_km: f32) -> Six {
    let polys = marching_squares(&to_mask(land, w, w), 0.5);
    let a = coast_shape_thresholds(&polys, cell_km, 1.0, NECK_KM);
    let b = coast_shape_thresholds(&polys, cell_km, 2.0 * cell_km, cell_km);
    Six { km: a.count, cell: b.count, coast: a.coast_km, p90: a.p90_len_km, r: b.local_axis_r }
}

/// Chebyshev distance to the boundary, in pixels of THIS grid, by multi-source BFS. Declared as
/// Chebyshev (8-connected unit steps), not Euclidean: at d = 0.5 px the band is exactly the
/// boundary-adjacent ring either way, and at 5 px the difference is a corner rounding.
fn band_counts(land: &[bool], w: usize, radii: &[usize]) -> Vec<(usize, usize)> {
    let n = w * w;
    let mut dist = vec![u16::MAX; n];
    let mut q: Vec<u32> = Vec::new();
    for k in 0..n {
        let (x, y) = (k % w, k / w);
        let mut edge = false;
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                if nx >= 0
                    && ny >= 0
                    && (nx as usize) < w
                    && (ny as usize) < w
                    && land[ny as usize * w + nx as usize] != land[k]
                {
                    edge = true;
                }
            }
        }
        if edge {
            dist[k] = 0;
            q.push(k as u32);
        }
    }
    let mut i = 0usize;
    while i < q.len() {
        let k = q[i] as usize;
        i += 1;
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= w {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if dist[nk] == u16::MAX {
                    dist[nk] = dist[k] + 1;
                    q.push(nk as u32);
                }
            }
        }
    }
    radii
        .iter()
        .map(|&r| {
            let sand = (0..n).filter(|&k| land[k] && (dist[k] as usize) < r.max(1)).count();
            let foam = (0..n).filter(|&k| !land[k] && (dist[k] as usize) < r.max(1)).count();
            (sand, foam)
        })
        .collect()
}

/// Land components of a single pixel — the slivers a shader paints entirely as beach.
fn single_pixel_islands(land: &[bool], w: usize) -> usize {
    let n = w * w;
    let mut c = 0usize;
    for k in 0..n {
        if !land[k] {
            continue;
        }
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let mut nb = 0usize;
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let (nx, ny) = (x + dx, y + dy);
                if nx >= 0
                    && ny >= 0
                    && (nx as usize) < w
                    && (ny as usize) < w
                    && land[ny as usize * w + nx as usize]
                {
                    nb += 1;
                }
            }
        }
        if nb == 0 {
            c += 1;
        }
    }
    c
}

fn control(tag: &str, f: &GridF32, pre: &GridF32, n2m: f32) {
    let n = f.data.len();
    let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
    let n_land = land.iter().filter(|&&l| l).count();
    let h = sorted((0..n).filter(|&k| land[k]).map(|k| (f.data[k] - SEA) * n2m).collect());
    let mean = h.iter().map(|&x| x as f64).sum::<f64>() / h.len().max(1) as f64;
    let flow = compute_flow(f, &FlowConfig { sea_level: SEA, ..Default::default() });
    let pits = (0..n).filter(|&k| flow.filled.data[k] > f.data[k] + 1e-7).count();
    let chan = (0..n).filter(|&k| land[k] && flow.accumulation.data[k] >= A_C_CELLS).count();
    let (mut wsum, mut wn) = (0.0f64, 0usize);
    for k in 0..n {
        if pre.data[k] > SEA && land[k] {
            wsum += ((pre.data[k] - f.data[k]) * n2m) as f64;
            wn += 1;
        }
    }
    let drowned = (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] <= SEA).count();
    eprintln!(
        "      CONTROL {tag:<24} land {n_land} | hypso mean {mean:.2} m p50 {:.2} | pits {pits} \
         | channel {:.3} % | work {:.2} m | DROWNED {drowned}",
        pct(&h, 0.50),
        100.0 * chan as f64 / n_land.max(1) as f64,
        wsum / wn.max(1) as f64
    );
}

#[test]
#[ignore]
fn base_level_run() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 80 · the base-level bound, eps = {EPS} m  ==========");

    let pre = build_field(Knobs::no_incision());
    let shipped = build_field(Knobs::shipped());
    let bounded = build_field(Knobs { base_level_m: Some(EPS), ..Knobs::shipped() });
    let w = shipped.width;

    // ── B0 · the consumer view, and the u16 step ──────────────────────────────
    let (_, step) = land_u16(&shipped, &ss);
    eprintln!(
        "\n── B0 · the consumer view ──\n   u16 step = **{step:.4} m** ⇒ eps = {EPS} m is \
         **{:.1} steps** (the bound needs at least 2)",
        EPS / step
    );

    let stages: Vec<(&str, &GridF32)> =
        vec![("PRE-INCISION", &pre), ("DELIVERED", &shipped), ("BASE LEVEL", &bounded)];
    let mut grids: Vec<(&str, Vec<(String, Vec<bool>, usize, f32)>)> = Vec::new();
    for (nm, f) in &stages {
        let l_f32: Vec<bool> = f.data.iter().map(|&v| v > SEA).collect();
        let (l_u16, _) = land_u16(f, &ss);
        let (l_2048, w2) = majority(&l_u16, w, 4);
        let (l_1024, w1) = majority(&l_u16, w, 8);
        grids.push((
            nm,
            vec![
                ("f32 8192".to_string(), l_f32, w, CELL_KM),
                ("u16 8192".to_string(), l_u16, w, CELL_KM),
                ("terrain 2048".to_string(), l_2048, w2, CELL_KM * 4.0),
                ("ocean 1024".to_string(), l_1024, w1, CELL_KM * 8.0),
            ],
        ));
    }

    // ── B1 · Δ on five masks ──────────────────────────────────────────────────
    eprintln!(
        "\n── B1 · the criterion, on every mask (Δ against PRE-INCISION on the same grid) ──"
    );
    eprintln!(
        "   {:<14} {:<14} {:>8} {:>9} {:>10} {:>7} {:>7} | {:>9} {:>11}",
        "grid", "stage", ">=1 km", ">=2 cells", "coast km", "p90", "R", "D >=1 km", "D >=2 cells"
    );
    for gi in 0..4 {
        let base = six(&grids[0].1[gi].1, grids[0].1[gi].2, grids[0].1[gi].3);
        for (si, (snm, cols)) in grids.iter().enumerate() {
            let (gnm, land, gw, ck) = &cols[gi];
            let s = six(land, *gw, *ck);
            eprintln!(
                "   {:<14} {snm:<14} {:>8} {:>9} {:>10.0} {:>7.2} {:>7.3} | {:>+9} {:>+11}",
                if si == 0 { gnm.as_str() } else { "" },
                s.km,
                s.cell,
                s.coast,
                s.p90,
                s.r,
                s.km as i64 - base.km as i64,
                s.cell as i64 - base.cell as i64
            );
        }
    }

    // the XOR of the two consumer shorelines — the "two coasts" cost
    eprintln!("\n   XOR of the terrain (2048) and ocean (1024) shorelines, on the 2048 grid:");
    for (nm, cols) in &grids {
        let (_, l2, w2, _) = &cols[2];
        let (_, l1, w1, _) = &cols[3];
        let mut xor = 0usize;
        for y in 0..*w2 {
            for x in 0..*w2 {
                let a = l2[y * w2 + x];
                let b = l1[(y / 2) * w1 + x / 2];
                if a != b {
                    xor += 1;
                }
            }
        }
        eprintln!(
            "     {nm:<14} {xor} cells disagree ({:.3} % of the 2048 grid)",
            100.0 * xor as f64 / (w2 * w2) as f64
        );
    }

    // ── B1b · WHICH operator drowns the residual ──────────────────────────────
    let n = shipped.data.len();
    let resid: Vec<usize> =
        (0..n).filter(|&k| pre.data[k] > SEA && bounded.data[k] <= SEA).collect();
    eprintln!(
        "\n── B1b · the residual drowned cells under the bound: **{}** (against {} delivered) ──",
        resid.len(),
        (0..n).filter(|&k| pre.data[k] > SEA && shipped.data[k] <= SEA).count()
    );
    if !resid.is_empty() {
        for (nm, kn) in [
            (
                "bound + talus OFF",
                Knobs { base_level_m: Some(EPS), talus_passes: Some(0), ..Knobs::shipped() },
            ),
            (
                "bound + diffusion OFF",
                Knobs { base_level_m: Some(EPS), diffusion: Some(0.0), ..Knobs::shipped() },
            ),
        ] {
            let f = build_field(kn);
            let d = (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] <= SEA).count();
            eprintln!("     {nm:<24} drowned {d}");
        }
    }

    // ── B2 · the price, per pinned cell ───────────────────────────────────────
    let floor_norm = SEA + EPS / n2m;
    let pinned: Vec<usize> = (0..n)
        .filter(|&k| bounded.data[k] > SEA && (bounded.data[k] - floor_norm).abs() < 1e-7)
        .collect();
    let prevented = sorted(
        pinned.iter().map(|&k| (bounded.data[k] - shipped.data[k]) * n2m).collect::<Vec<f32>>(),
    );
    eprintln!(
        "\n── B2 · the price ──\n   cells pinned at sea + eps: **{}** | incision PREVENTED \
         (bounded - delivered, m): p10 {:.3} MEDIAN **{:.3}** p90 {:.3} max {:.2}",
        pinned.len(),
        pct(&prevented, 0.10),
        pct(&prevented, 0.50),
        pct(&prevented, 0.90),
        prevented.last().copied().unwrap_or(f32::NAN)
    );

    // ── B3 · the control block ────────────────────────────────────────────────
    eprintln!("\n── B3 · the control block ──");
    for (nm, f) in &stages {
        control(nm, f, &pre, n2m);
    }

    // ── B4 · the shader metrics, on the consumer grids ────────────────────────
    eprintln!("\n── B4 · shader metrics on the LL grids (band radii in PIXELS of that grid) ──");
    eprintln!(
        "   {:<14} {:<14} {:>10} {:>12} {:>12} {:>14}",
        "grid", "stage", "contour km", "sand d<1 px", "foam d<1 px", "1-px islands"
    );
    for gi in 2..4 {
        for (snm, cols) in &grids {
            let (gnm, land, gw, ck) = &cols[gi];
            let s = six(land, *gw, *ck);
            let b = band_counts(land, *gw, &[1, 5]);
            eprintln!(
                "   {:<14} {snm:<14} {:>10.0} {:>12} {:>12} {:>14}",
                gnm,
                s.coast,
                b[0].0,
                b[0].1,
                single_pixel_islands(land, *gw)
            );
        }
    }
    // the 4-to-8-cell population at 8192, never isolated before
    for (snm, cols) in &grids {
        let (_, l8, _, _) = &cols[1];
        let mut sizes: Vec<usize> = Vec::new();
        let mut seen = vec![false; n];
        for s0 in 0..n {
            if !l8[s0] || seen[s0] {
                continue;
            }
            let mut stack = vec![s0 as u32];
            seen[s0] = true;
            let mut c = 0usize;
            while let Some(kk) = stack.pop() {
                let k = kk as usize;
                c += 1;
                if c > 64 {
                    break;
                }
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                for dy in -1i32..=1 {
                    for dx in -1i32..=1 {
                        let (nx, ny) = (x + dx, y + dy);
                        if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < w {
                            let nk = ny as usize * w + nx as usize;
                            if l8[nk] && !seen[nk] {
                                seen[nk] = true;
                                stack.push(nk as u32);
                            }
                        }
                    }
                }
            }
            sizes.push(c);
        }
        eprintln!(
            "     {snm:<14} 8192 land components of 4-8 cells: **{}** (1-3 cells: {})",
            sizes.iter().filter(|&&c| (4..=8).contains(&c)).count(),
            sizes.iter().filter(|&&c| c <= 3).count()
        );
    }

    // ── B5 · sensitivity to epsilon ───────────────────────────────────────────
    eprintln!("\n── B5 · epsilon sweep (Δ >=2 cells on the u16 8192 mask) ──");
    let base_u16 = six(&grids[0].1[1].1, w, CELL_KM);
    for e in [0.2f32, 0.5, 1.0] {
        let f = if (e - EPS).abs() < 1e-6 {
            bounded.clone()
        } else {
            build_field(Knobs { base_level_m: Some(e), ..Knobs::shipped() })
        };
        let (l_u16, _) = land_u16(&f, &ss);
        let l_f32: Vec<bool> = f.data.iter().map(|&v| v > SEA).collect();
        let su = six(&l_u16, w, CELL_KM);
        let sf = six(&l_f32, w, CELL_KM);
        let fnorm = SEA + e / n2m;
        let pin = (0..n).filter(|&k| f.data[k] > SEA && (f.data[k] - fnorm).abs() < 1e-7).count();
        eprintln!(
            "   eps {e:>4} m ({:>4.1} u16 steps) | Δ(f32) {:>+7} | Δ(u16) {:>+7} | pinned cells \
             {pin}",
            e / step,
            sf.cell as i64 - base_u16.cell as i64,
            su.cell as i64 - base_u16.cell as i64
        );
    }
}

/// ADR Finding 80 B2, CORRECTED. The first pass compared the bounded field to the DELIVERED one,
/// whose bathymetry clamp puts every drowned cell at exactly -1.0 m, so "incision prevented" read
/// a constant 1.500 m: the clamp, not the incision. **Method rule 12 applied to my own
/// instrument.** Both sides are rebuilt with `bathymetry_off`, which is the stage the bound acts
/// at. Also carries the paired hypsometry the unpaired mean of B3 cannot give.
#[test]
#[ignore]
fn base_level_price_preclamp() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!(
        "
==========  Finding 80 B2 corrected · the price, on the PRE-CLAMP field  =========="
    );
    let free = build_field(Knobs { bathymetry_off: true, ..Knobs::shipped() });
    let bound =
        build_field(Knobs { bathymetry_off: true, base_level_m: Some(EPS), ..Knobs::shipped() });
    let pre = build_field(Knobs::no_incision());
    let (w, n) = (free.width, free.data.len());
    let floor_norm = SEA + EPS / n2m;

    // the price: at every cell the bound pinned, how much incision did it refuse?
    let pinned: Vec<usize> =
        (0..n).filter(|&k| (bound.data[k] - floor_norm).abs() < 1e-7).collect();
    let prevented =
        sorted(pinned.iter().map(|&k| (bound.data[k] - free.data[k]) * n2m).collect::<Vec<f32>>());
    eprintln!(
        "   pinned cells **{}** | incision PREVENTED (bound - free, PRE-CLAMP), m: p10 {:.3}          p50 **{:.3}** p90 **{:.3}** p99 {:.2} max {:.2}",
        pinned.len(),
        pct(&prevented, 0.10),
        pct(&prevented, 0.50),
        pct(&prevented, 0.90),
        pct(&prevented, 0.99),
        prevented.last().copied().unwrap_or(f32::NAN)
    );

    // the pinned REACHES: runs of consecutive pinned cells along the flow
    let flow = compute_flow(&bound, &FlowConfig { sea_level: SEA, ..Default::default() });
    let is_pin: Vec<bool> = (0..n).map(|k| (bound.data[k] - floor_norm).abs() < 1e-7).collect();
    let mut seen = vec![false; n];
    let mut runs: Vec<f32> = Vec::new();
    for s0 in 0..n {
        if !is_pin[s0] || seen[s0] {
            continue;
        }
        // walk downstream from an unvisited pinned cell
        let (mut cur, mut len) = (s0, 0usize);
        while is_pin[cur] && !seen[cur] && len < 10_000 {
            seen[cur] = true;
            len += 1;
            let d = flow.direction[cur];
            if d == ymir_core::terrain::flow::DIR_NONE {
                break;
            }
            let nx = ((cur % w) as i32 + ymir_core::terrain::flow::D8_DX[d as usize])
                .rem_euclid(w as i32) as usize;
            let ny = ((cur / w) as i32 + ymir_core::terrain::flow::D8_DY[d as usize])
                .rem_euclid(w as i32) as usize;
            cur = ny * w + nx;
        }
        runs.push(len as f32);
    }
    let runs = sorted(runs);
    eprintln!(
        "   PINNED REACHES (consecutive pinned cells along the flow): {} runs | p50 **{:.0}**          p90 **{:.0}** max {:.0} cells",
        runs.len(),
        pct(&runs, 0.50),
        pct(&runs, 0.90),
        runs.last().copied().unwrap_or(f32::NAN)
    );

    // the PAIRED hypsometry: only cells that are land in BOTH, so the -17.26 m of B3 can be
    // split into "erosion changed" and "290 895 low cells rejoined the population".
    let shipped = build_field(Knobs::shipped());
    let bounded_clamped = build_field(Knobs { base_level_m: Some(EPS), ..Knobs::shipped() });
    let both: Vec<usize> =
        (0..n).filter(|&k| shipped.data[k] > SEA && bounded_clamped.data[k] > SEA).collect();
    let mean = |f: &GridF32, set: &[usize]| {
        set.iter().map(|&k| ((f.data[k] - SEA) * n2m) as f64).sum::<f64>() / set.len() as f64
    };
    let all_s: Vec<usize> = (0..n).filter(|&k| shipped.data[k] > SEA).collect();
    let all_b: Vec<usize> = (0..n).filter(|&k| bounded_clamped.data[k] > SEA).collect();
    eprintln!(
        "
   HYPSOMETRY, decomposed (rule 3):
     UNPAIRED  delivered {:.2} m over {} cells          -> bounded {:.2} m over {} cells  = **{:+.2} m**
     PAIRED (the {} cells that are          land in BOTH)  {:.2} -> {:.2} m  = **{:+.2} m**
     the difference is the {} low          cells the bound gave back, at a mean of {:.2} m",
        mean(&shipped, &all_s),
        all_s.len(),
        mean(&bounded_clamped, &all_b),
        all_b.len(),
        mean(&bounded_clamped, &all_b) - mean(&shipped, &all_s),
        both.len(),
        mean(&shipped, &both),
        mean(&bounded_clamped, &both),
        mean(&bounded_clamped, &both) - mean(&shipped, &both),
        all_b.len() - both.len(),
        {
            let back: Vec<usize> =
                all_b.iter().copied().filter(|&k| shipped.data[k] <= SEA).collect();
            mean(&bounded_clamped, &back)
        }
    );
    let _ = pre;
}
