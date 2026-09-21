//! ADR 0001 Finding 111 — **is the river comb a D8 artefact, a piedmont, or the instability the
//! dossier already confirmed?**
//!
//! ⛔ **P0, stated before any measurement: block B is NOT answerable as the round writes it.**
//! *"Un artefact D8 : θ_réseau pique sur 0/45/90°"* — **every river point lies on the D8 grid**.
//! `extract_rivers` follows `flow_result.direction`, which has exactly eight values, so **every
//! consecutive step is a multiple of 45° by construction**, in the comb, in a dendritic basin, in
//! any world. A test that returns the same answer for both hypotheses discriminates nothing.
//!
//! ⇒ θ is therefore read at **four scales** — the step, and chords of 8, 16 and 32 cells — so the
//! reader watches the quantisation die instead of taking the repair on trust.
//!
//! ⛔ **And rule 11 found a THIRD hypothesis the round does not offer, already CONFIRMED.**
//! ADR L460 (Finding 10): *"**Smith–Bretherton parallel rilling (CONFIRMED)**. On a SMOOTH plane
//! tilted 30° off the grid, no FBM, the incision spontaneously forms regularly-spaced parallel
//! rills running straight DOWNSLOPE (↘, diagonal — **following the slope, NOT the grid axes**; v1
//! and v2 give the identical concentration R = 0.19, **exonerating the solver**)… the classic
//! linear instability of `E = K·A^m·S^n` (m < 1) on smooth slopes, damped only by hillslope
//! diffusion… the regime split EXCLUDES channel cells from diffusion, so the rills are never
//! damped laterally."*
//!
//! So the dossier has already run this round's block B on a control where the grid was excluded by
//! construction, and the rills followed the slope. **That object looks exactly like a piedmont —
//! parallel channels down the slope — and is not one: it is a numerical instability of the law,
//! with no deposition anywhere in `relief_v3`.** It also predicts Finding 110's measurement: the
//! closure FLATTENS the beds (σ 3.89° → 1.28°), and Smith–Bretherton is an instability *of smooth
//! slopes*.
//!
//! ⇒ A fourth measurement is therefore added, the one that separates the three: **cross-slope
//! SPACING**. Rilling is regularly spaced (a wavelength); a piedmont's channels are irregular and
//! **converge downstream** (fans merge); a grid artefact is tied to the cell pitch.
//!
//! ⚠️ Rule 11: `piedmont`, `piémont`, `alluvial`, `fan`, and every real-world reference —
//! **NOTHING FOUND**. Ten rounds of "comb" and the dossier never once asked whether it was a
//! landform. `parallel` 51 hits, earliest relevant **L460** above. `structure tensor` 2 (L14741
//! records it as "NOTHING FOUND" at Finding 97's own grep). `R8` 112, earliest L3916.
//!
//! ⚠️ Calibration carried forward from Finding 97, unchanged: R8 reads **0.0402** on a provably
//! isotropic synthetic and **0.9999** on a 45° striped one.
//!
//! Run: cargo test -p ymir-core --release --test f111_piedmont -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, pct, sorted};
use std::path::Path;
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
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const TILE: usize = 512;
/// Finding 110-B's comb tile: max Δ(R network) between Finding 109 OFF and ON.
const COMB: (usize, usize) = (5120, 4096);
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f111_eye");

/// Circular statistics of a set of axial directions (mod π): `R2`, `R8` and the dominant axis.
fn axial(thetas: &[f32]) -> (f32, f32, f32) {
    if thetas.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let (mut c2, mut s2, mut c8, mut s8) = (0f64, 0f64, 0f64, 0f64);
    for &t in thetas {
        let t = t as f64;
        c2 += (2.0 * t).cos();
        s2 += (2.0 * t).sin();
        c8 += (8.0 * t).cos();
        s8 += (8.0 * t).sin();
    }
    let n = thetas.len() as f64;
    let th = 0.5 * s2.atan2(c2);
    let mut d = th.to_degrees();
    if d < 0.0 {
        d += 180.0;
    }
    (((c2 * c2 + s2 * s2).sqrt() / n) as f32, ((c8 * c8 + s8 * s8).sqrt() / n) as f32, d as f32)
}

/// Chord directions of a polyline at scale `l` cells, as axial angles in radians (mod π).
fn chords(points: &[(u32, u32)], l: usize) -> Vec<f32> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i + l < points.len() {
        let (x0, y0) = (points[i].0 as f32, points[i].1 as f32);
        let (x1, y1) = (points[i + l].0 as f32, points[i + l].1 as f32);
        let (dx, dy) = (x1 - x0, y1 - y0);
        if dx != 0.0 || dy != 0.0 {
            let mut t = dy.atan2(dx);
            if t < 0.0 {
                t += std::f32::consts::PI;
            }
            out.push(t);
        }
        i += l.max(1);
    }
    out
}

fn spearman(a: &[f32], b: &[f32]) -> f32 {
    let rank = |v: &[f32]| -> Vec<f32> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
        let mut r = vec![0f32; v.len()];
        for (p, &i) in idx.iter().enumerate() {
            r[i] = p as f32;
        }
        r
    };
    let (ra, rb) = (rank(a), rank(b));
    let n = a.len() as f32;
    let (ma, mb) = (ra.iter().sum::<f32>() / n, rb.iter().sum::<f32>() / n);
    let (mut num, mut da, mut db) = (0f32, 0f32, 0f32);
    for i in 0..a.len() {
        num += (ra[i] - ma) * (rb[i] - mb);
        da += (ra[i] - ma) * (ra[i] - ma);
        db += (rb[i] - mb) * (rb[i] - mb);
    }
    if da * db <= 0.0 { 0.0 } else { num / (da * db).sqrt() }
}

fn save(g: &GridF32, name: &str) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    let (lo, hi) = g.data.iter().fold((f32::MAX, f32::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    let span = (hi - lo).max(f32::EPSILON);
    let mut o = GridF32::new(g.width, g.height, 0.0);
    for (d, &v) in o.data.iter_mut().zip(g.data.iter()) {
        *d = (v - lo) / span;
    }
    o.save_png_u8(&Path::new(OUT).join(name)).expect("png");
    eprintln!("   wrote {name}");
}

#[test]
#[ignore]
fn f111_piedmont() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 111 . D8 artefact, piedmont, or Smith-Bretherton?  =========");

    let f = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let (w, h) = (f.width, f.height);
    let nt = w.div_ceil(TILE);
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

    // ── P0: the quantisation, at four scales, over the WHOLE network ────────
    eprintln!("\n   ── P0 · θ of the network at four scales, whole continent ──");
    for l in [1usize, 8, 16, 32] {
        let mut th = Vec::new();
        for s in &dr.rivers.segments {
            th.extend(chords(&s.points, l));
        }
        let (r2, r8, t) = axial(&th);
        eprintln!(
            "      chord **{l:>2} cells** ({:>5.2} km): n = {:>7} · R2 **{r2:.4}** · **R8 \
             {r8:.4}** · θ {t:>5.1}°{}",
            l as f32 * CELL_KM,
            th.len(),
            if l == 1 { "  ⇐ 1.0000 expected BY CONSTRUCTION: D8 has 8 directions" } else { "" }
        );
    }

    // ── A: the cloud, every land tile ───────────────────────────────────────
    let mut rows: Vec<(usize, f32, f32, f32, f32, f32, f32)> = Vec::new(); // tile,R2net,θnet,R8terr,θterr,R8net32,sigma
    for ty in 0..nt {
        for tx in 0..nt {
            let (x0, y0) = (tx * TILE, ty * TILE);
            let (cw, ch) = (TILE.min(w - x0), TILE.min(h - y0));
            let mut crop = GridF32::new(cw, ch, 0.0);
            for y in 0..ch {
                for x in 0..cw {
                    crop.data[y * cw + x] = bf.data[(y0 + y) * w + x0 + x];
                }
            }
            let land: Vec<bool> = crop.data.iter().map(|&v| v > SEA).collect();
            let land_share = land.iter().filter(|&&b| b).count() as f32 / land.len() as f32;
            if land_share < 0.9 {
                continue; // declared: land tiles only, ≥ 90 % land
            }
            let mut th32 = Vec::new();
            let mut th1 = Vec::new();
            for s in &dr.rivers.segments {
                let inside: Vec<(u32, u32)> = s
                    .points
                    .iter()
                    .copied()
                    .filter(|&(x, y)| {
                        (x as usize) >= x0
                            && (x as usize) < x0 + cw
                            && (y as usize) >= y0
                            && (y as usize) < y0 + ch
                    })
                    .collect();
                th32.extend(chords(&inside, 32));
                th1.extend(chords(&inside, 1));
            }
            if th32.len() < 32 {
                continue;
            }
            let (r2n, r8n, tn) = axial(&th32);
            let a = aniso(&crop, &land, 16);
            let mut sl = Vec::new();
            for y in 1..ch - 1 {
                for x in 1..cw - 1 {
                    let gx = (crop.data[y * cw + x + 1] - crop.data[y * cw + x - 1]) * n2m
                        / (2.0 * CELL_M);
                    let gy = (crop.data[(y + 1) * cw + x] - crop.data[(y - 1) * cw + x]) * n2m
                        / (2.0 * CELL_M);
                    sl.push((gx * gx + gy * gy).sqrt().atan().to_degrees());
                }
            }
            let m = sl.iter().sum::<f32>() / sl.len().max(1) as f32;
            let sg =
                (sl.iter().map(|v| (v - m) * (v - m)).sum::<f32>() / sl.len().max(1) as f32).sqrt();
            rows.push((ty * nt + tx, r2n, tn, a.r8, a.theta2_deg, r8n, sg));
        }
    }
    let dth = |a: f32, b: f32| -> f32 {
        let d = (a - b).abs() % 180.0;
        if d > 90.0 { 180.0 - d } else { d }
    };
    let rn: Vec<f32> = rows.iter().map(|r| r.1).collect();
    let rt: Vec<f32> = rows.iter().map(|r| r.3).collect();
    eprintln!(
        "\n   ── A · {} land tiles (≥ 90 % land, ≥ 32 chords) ──\n      Spearman R2(network, \
         chord 32) vs R8(terrain) = **{:.3}**",
        rows.len(),
        spearman(&rn, &rt)
    );
    let mut hi: Vec<&(usize, f32, f32, f32, f32, f32, f32)> =
        rows.iter().filter(|r| r.1 >= 0.25).collect();
    hi.sort_by(|a, b| b.1.total_cmp(&a.1));
    let dths: Vec<f32> = hi.iter().map(|r| dth(r.2, r.4)).collect();
    eprintln!(
        "      high-R tiles (R2 network ≥ 0.25): **{}** of {} · |θ_net − θ_terr| p50 **{:.1}°** \
         p90 **{:.1}°** · share within 20° **{:.0} %** · share within 45±10° **{:.0} %**",
        hi.len(),
        rows.len(),
        if dths.is_empty() { 0.0 } else { pct(&sorted(dths.clone()), 0.50) },
        if dths.is_empty() { 0.0 } else { pct(&sorted(dths.clone()), 0.90) },
        100.0 * dths.iter().filter(|&&d| d <= 20.0).count() as f32 / dths.len().max(1) as f32,
        100.0 * dths.iter().filter(|&&d| (35.0..=55.0).contains(&d)).count() as f32
            / dths.len().max(1) as f32
    );
    eprintln!("      the ten most anisotropic tiles:");
    for r in hi.iter().take(10) {
        let (tx, ty) = ((r.0 % nt) * TILE, (r.0 / nt) * TILE);
        eprintln!(
            "         ({tx:>4},{ty:>4}) R2_net **{:.3}** θ_net {:>5.1}° · R8_terr **{:.4}** \
             θ_terr {:>5.1}° ⇒ **|Δθ| {:>5.1}°** · R8_net(32) {:.3} · σ(slope) {:>5.2}°",
            r.1,
            r.2,
            r.3,
            r.4,
            dth(r.2, r.4),
            r.5,
            r.6
        );
    }

    // ── B: the comb tile, and the negative control ──────────────────────────
    let max_wc = (0..dr.rivers.segments.len())
        .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse)
        .max_by(|&a, &b| dr.segment_catchment_cells[a].total_cmp(&dr.segment_catchment_cells[b]));
    let ctrl = max_wc
        .and_then(|i| dr.rivers.segments[i].points.first().copied())
        .map(|(x, y)| ((x as usize / TILE) * TILE, (y as usize / TILE) * TILE))
        .unwrap_or((0, 0));
    for (name, (x0, y0)) in [("THE COMB TILE", COMB), ("the dendritic control", ctrl)] {
        eprintln!("\n   ── B · {name} at ({x0}, {y0}) ──");
        let (cw, ch) = (TILE.min(w - x0), TILE.min(h - y0));
        let mut crop = GridF32::new(cw, ch, 0.0);
        for y in 0..ch {
            for x in 0..cw {
                crop.data[y * cw + x] = bf.data[(y0 + y) * w + x0 + x];
            }
        }
        let land: Vec<bool> = crop.data.iter().map(|&v| v > SEA).collect();
        let a = aniso(&crop, &land, 16);
        for l in [1usize, 8, 16, 32] {
            let mut th = Vec::new();
            for s in &dr.rivers.segments {
                let inside: Vec<(u32, u32)> = s
                    .points
                    .iter()
                    .copied()
                    .filter(|&(x, y)| {
                        (x as usize) >= x0
                            && (x as usize) < x0 + cw
                            && (y as usize) >= y0
                            && (y as usize) < y0 + ch
                    })
                    .collect();
                th.extend(chords(&inside, l));
            }
            let (r2, r8, t) = axial(&th);
            eprintln!(
                "      network chord {l:>2}: n = {:>6} · R2 **{r2:.4}** · R8 **{r8:.4}** · θ \
                 **{t:>5.1}°**",
                th.len()
            );
        }
        eprintln!(
            "      TERRAIN (structure tensor, w = 16): R2 **{:.4}** · R8 **{:.4}** · θ \
             **{:>5.1}°** · C p50 {:.3} · H(θ) {:.4}",
            a.r2, a.r8, a.theta2_deg, a.c_p50, a.h_theta
        );
        save(&crop, &format!("{}_height.png", if x0 == COMB.0 { "comb" } else { "control" }));
        // the network, burnt in, so the eye can see whether the channels MERGE downstream
        let mut net = crop.clone();
        let (lo, hi2) =
            crop.data.iter().fold((f32::MAX, f32::MIN), |(a, b), &v| (a.min(v), b.max(v)));
        for s in &dr.rivers.segments {
            for &(x, y) in &s.points {
                let (x, y) = (x as usize, y as usize);
                if x >= x0 && x < x0 + cw && y >= y0 && y < y0 + ch {
                    net.data[(y - y0) * cw + (x - x0)] = if x0 == COMB.0 { lo } else { hi2 };
                }
            }
        }
        save(&net, &format!("{}_network.png", if x0 == COMB.0 { "comb" } else { "control" }));
    }

    eprintln!("\n==========  end Finding 111 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
