//! ADR 0001 Finding 110-B — **is the river comb the D8 routing exported off a flattened surface?**
//!
//! The author saw, in the viz with Finding 109 ON: **384 → 667 rivers, 60 → 38 lakes** (both
//! targets) and a comb on the drainage layer. Block A is settled by the author — `cross_rill` OFF,
//! the comb survives — so this bench is block B and C only.
//!
//! ⛔ **Rule 11 found the answer already written, at Finding 65**, and it is the earliest hit that
//! matters: *"THREE incompatible answers to 'is there a channel here', and the one that carved the
//! terrain is used by nothing"* — definition 1 (`mfd_accumulation`, `mfd_exponent = 2`) **carves
//! the terrain and nothing else**; definition 2 (**D8**, `extract_rivers` reading
//! `flow_result.direction`) feeds *"rivers.json, drainage density, catchment_km2, drainage_km2,
//! Strahler order, the microscope"* — and the viz drainage layer. At 8192² the two disagree by
//! **5.25×** on the same field. Finding 11 (ADR L583) bounded the MFD remedy deliberately: *"MFD
//! for the stream-power accumulation/incision only, **D8 kept for rivers/lakes (blast radius
//! contained)**"*.
//!
//! ⇒ The comb is on the layer that was knowingly left on D8. **This bench measures that rather
//! than asserting it**, because "it must be D8" is the kind of sentence Findings 101 and 107 were
//! both written to punish.
//!
//! ⚠️ **Declared chain (rule 12).** 8192², the EXPORT chain: breached field, `assemble_hd_drainage`,
//! `head_km2 = RELIEF_V1_A_C_KM2 = 0.1`, `full_tree = false`, `stream_km2 = 20` — **byte-for-byte
//! the viz's own config** (`hd.rs:701-702` plus the `C1DrainageThresholds` defaults), so the layer
//! measured here is the layer the author looked at. ⚠️ But the viz runs at **2048²** by default
//! and this bench at 8192², so **the author's 384 → 667 is NOT reproduced here and is not claimed**;
//! the Δ below is this bench's own.
//!
//! Run: cargo test -p ymir-core --release --test f110_comb -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, pct, sorted};
use std::collections::HashMap;
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
/// Tile side in cells for the R / R8 maps. 512 = 25 km, the Finding 97 tile.
const TILE: usize = 512;

/// One segment, reduced to what block B asks about.
struct Seg {
    order: u8,
    km2: f32,
    kind: SegmentKind,
    /// Mean local slope of the TERRAIN under the bed, degrees.
    slope: f32,
    /// Standard deviation of that slope along the bed — the flatness statistic.
    sigma: f32,
    tile: usize,
}

fn slope_deg_at(f: &GridF32, k: usize, n2m: f32, w: usize, h: usize) -> f32 {
    let (x, y) = ((k % w) as i32, (k / w) as i32);
    let at = |dx: i32, dy: i32| -> f32 {
        let nx = (x + dx).rem_euclid(w as i32) as usize;
        let ny = (y + dy).rem_euclid(h as i32) as usize;
        f.data[ny * w + nx] * n2m
    };
    let gx = (at(1, 0) - at(-1, 0)) / (2.0 * CELL_M);
    let gy = (at(0, 1) - at(0, -1)) / (2.0 * CELL_M);
    (gx * gx + gy * gy).sqrt().atan().to_degrees()
}

fn std_dev(v: &[f32]) -> f32 {
    if v.len() < 2 {
        return 0.0;
    }
    let m = v.iter().sum::<f32>() / v.len() as f32;
    (v.iter().map(|x| (x - m) * (x - m)).sum::<f32>() / v.len() as f32).sqrt()
}

/// Finding 52's axis R, on the STEP directions of the river network: `|mean(e^{2iθ})|` over every
/// consecutive point pair in a tile. 0 = isotropic, 1 = every step parallel. The 2nd harmonic is
/// the right one here — a comb is 2-fold (parallel lines), not the 4-fold D8 pattern R8 catches.
fn river_r_per_tile(segs: &[(Vec<(u32, u32)>, usize)], w: usize, nt: usize) -> Vec<f32> {
    let mut acc = vec![(0.0f64, 0.0f64, 0usize); nt * nt];
    for (pts, _) in segs {
        for p in pts.windows(2) {
            let (x0, y0) = (p[0].0 as f32, p[0].1 as f32);
            let (x1, y1) = (p[1].0 as f32, p[1].1 as f32);
            let (dx, dy) = (x1 - x0, y1 - y0);
            if dx == 0.0 && dy == 0.0 {
                continue;
            }
            let th = dy.atan2(dx);
            let t = (p[0].1 as usize / TILE) * nt + (p[0].0 as usize / TILE);
            if t < acc.len() {
                acc[t].0 += (2.0 * th as f64).cos();
                acc[t].1 += (2.0 * th as f64).sin();
                acc[t].2 += 1;
            }
        }
    }
    let _ = w;
    acc.iter()
        .map(|&(c, s, n)| if n < 64 { 0.0 } else { ((c * c + s * s).sqrt() / n as f64) as f32 })
        .collect()
}

#[test]
#[ignore]
fn f110_comb() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 110-B . the river comb, attributed  ==========");
    eprintln!(
        "   declared: stream_km2 = **{}** (trunk), head_km2 = **{}** (headward extension, main \
         stem only) — the viz's own config",
        dcfg.thresholds.stream_km2, dcfg.thresholds.head_km2
    );

    let mut per_label: Vec<(String, Vec<Seg>, Vec<f32>, GridF32)> = Vec::new();
    for (label, knobs) in [
        ("F109 OFF (delivered)", Knobs::passes(2)),
        ("F109 ON  (s_eq 0.024)", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
    ] {
        eprintln!("\n--------  {label}  --------");
        let f = build_field_seed(knobs, PSEED);
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

        let mut segs = Vec::new();
        let mut pts_for_r = Vec::new();
        for (i, s) in dr.rivers.segments.iter().enumerate() {
            if s.points.len() < 2 {
                continue;
            }
            let sl: Vec<f32> = s
                .points
                .iter()
                .map(|&(x, y)| slope_deg_at(&bf, y as usize * w + x as usize, n2m, w, h))
                .collect();
            let (x0, y0) = s.points[0];
            segs.push(Seg {
                order: s.strahler_order,
                km2: dr.segment_catchment_cells[i] * cell_km2,
                kind: dr.segment_kind[i],
                slope: sl.iter().sum::<f32>() / sl.len() as f32,
                sigma: std_dev(&sl),
                tile: (y0 as usize / TILE) * nt + (x0 as usize / TILE),
            });
            pts_for_r.push((s.points.clone(), i));
        }
        let r = river_r_per_tile(&pts_for_r, w, nt);

        // ── the counts the author is looking at ─────────────────────────────
        let wc = segs.iter().filter(|s| s.kind == SegmentKind::Watercourse).count();
        let sw = segs.len() - wc;
        let mut by_order: HashMap<u8, usize> = HashMap::new();
        for s in segs.iter().filter(|s| s.kind == SegmentKind::Watercourse) {
            *by_order.entry(s.order).or_default() += 1;
        }
        let mut ord: Vec<_> = by_order.iter().map(|(k, v)| (*k, *v)).collect();
        ord.sort();
        eprintln!(
            "   segments **{}** = {wc} Watercourse + {sw} Spillway · by Strahler: {:?}",
            segs.len(),
            ord
        );
        let o1: Vec<&Seg> =
            segs.iter().filter(|s| s.kind == SegmentKind::Watercourse && s.order == 1).collect();
        let o2: Vec<&Seg> =
            segs.iter().filter(|s| s.kind == SegmentKind::Watercourse && s.order >= 2).collect();
        let q =
            |v: Vec<f32>, p: f64| -> f32 { if v.is_empty() { 0.0 } else { pct(&sorted(v), p) } };
        eprintln!(
            "   order 1  (n = {:>4}): catchment p50 **{:>8.2} km²** (p10 {:.2}, p90 {:.2}) · bed \
             slope p50 **{:>5.2}°** · σ(slope) p50 **{:>5.2}°**",
            o1.len(),
            q(o1.iter().map(|s| s.km2).collect(), 0.50),
            q(o1.iter().map(|s| s.km2).collect(), 0.10),
            q(o1.iter().map(|s| s.km2).collect(), 0.90),
            q(o1.iter().map(|s| s.slope).collect(), 0.50),
            q(o1.iter().map(|s| s.sigma).collect(), 0.50),
        );
        eprintln!(
            "   order ≥2 (n = {:>4}): catchment p50 **{:>8.2} km²** (p10 {:.2}, p90 {:.2}) · bed \
             slope p50 **{:>5.2}°** · σ(slope) p50 **{:>5.2}°**",
            o2.len(),
            q(o2.iter().map(|s| s.km2).collect(), 0.50),
            q(o2.iter().map(|s| s.km2).collect(), 0.10),
            q(o2.iter().map(|s| s.km2).collect(), 0.90),
            q(o2.iter().map(|s| s.slope).collect(), 0.50),
            q(o2.iter().map(|s| s.sigma).collect(), 0.50),
        );
        eprintln!(
            "   below A_c ({RELIEF_V1_A_C_KM2} km²): **{}** of {} Watercourse · below stream_km2 \
             ({}): **{}**",
            segs.iter()
                .filter(|s| s.kind == SegmentKind::Watercourse && s.km2 < RELIEF_V1_A_C_KM2)
                .count(),
            wc,
            dcfg.thresholds.stream_km2,
            segs.iter()
                .filter(|s| s.kind == SegmentKind::Watercourse
                    && s.km2 < dcfg.thresholds.stream_km2)
                .count()
        );
        per_label.push((label.to_string(), segs, r, bf));
    }

    // ── the Δ, and the comb tile ────────────────────────────────────────────
    let (off, on) = (&per_label[0], &per_label[1]);
    let w = off.3.width;
    let nt = w.div_ceil(TILE);
    let mut best = (0usize, f32::MIN);
    for t in 0..on.2.len().min(off.2.len()) {
        let d = on.2[t] - off.2[t];
        if on.2[t] > 0.0 && off.2[t] > 0.0 && d > best.1 {
            best = (t, d);
        }
    }
    let (tx, ty) = ((best.0 % nt) * TILE, (best.0 / nt) * TILE);
    eprintln!(
        "\n   ── the comb tile: max Δ(R_rivers) ON − OFF ──\n   tile **({tx}, {ty})** {TILE}² = \
         {:.0} km across · R_rivers **{:.4} → {:.4}** (Δ **{:+.4}**)",
        TILE as f32 * CELL_KM,
        off.2[best.0],
        on.2[best.0],
        best.1
    );
    // ── C1: the terrain's own R8 on that tile ───────────────────────────────
    for (label, _, _, bf) in [off, on] {
        let h = bf.height;
        let mut crop = GridF32::new(TILE.min(w - tx), TILE.min(h - ty), 0.0);
        for y in 0..crop.height {
            for x in 0..crop.width {
                crop.data[y * crop.width + x] = bf.data[(ty + y) * w + tx + x];
            }
        }
        let land: Vec<bool> = crop.data.iter().map(|&v| v > SEA).collect();
        eprintln!(
            "   C1 · TERRAIN R8 on that tile, {label}: **{:.4}** (Finding 109 whole-field: 0.0484)",
            aniso(&crop, &land, 16).r8
        );
    }
    // ── the Δ by Strahler order ─────────────────────────────────────────────
    let count = |s: &Vec<Seg>, o: u8| {
        s.iter().filter(|x| x.kind == SegmentKind::Watercourse && x.order == o).count() as i64
    };
    eprintln!("\n   ── Δ(segments) by Strahler order, Watercourse only ──");
    let mut tot = 0i64;
    for o in 1..=6u8 {
        let (a, b) = (count(&off.1, o), count(&on.1, o));
        if a + b == 0 {
            continue;
        }
        tot += b - a;
        eprintln!("      order {o}: **{a} → {b}** ⇒ **{:+}**", b - a);
    }
    eprintln!(
        "      ⇒ total **{:+}**, of which order 1 is **{:+}** = **{:.0} %**",
        tot,
        count(&on.1, 1) - count(&off.1, 1),
        100.0 * (count(&on.1, 1) - count(&off.1, 1)) as f32 / tot.max(1) as f32
    );

    // ── C2: the Finding 87-A river ──────────────────────────────────────────
    eprintln!(
        "\n   ── C2 · the Finding 87-A river ──\n   ⚠️ It is a **Spillway**, and Finding 87-F \
         records that `strahler_order` is MEANINGLESS on one (*\"the kind is authoritative\"*), so \
         the round's \"Strahler ≥ 3\" is not answerable as posed. The biggest **Watercourse** is \
         reported instead, which is answerable."
    );
    for (label, segs, _, _) in [off, on] {
        if let Some(s) = segs
            .iter()
            .filter(|s| s.kind == SegmentKind::Watercourse)
            .max_by(|a, b| a.km2.total_cmp(&b.km2))
        {
            eprintln!(
                "      {label}: biggest Watercourse **{:.0} km²** · Strahler **{}** · bed slope \
                 {:.2}° · σ {:.2}°",
                s.km2, s.order, s.slope, s.sigma
            );
        }
    }
    eprintln!(
        "\n==========  end Finding 110-B . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}

/// ADR Finding 110-B2 — **is `strahler_order` answerable at all on this network?**
///
/// Block B's question is "order 1 everywhere = non-convergent parallel lines". That reading
/// assumes the order reflects CONVERGENCE. The main bench found the biggest `Watercourse` —
/// **79 669 km²** — reading order **1** with Finding 109 ON and order **5** with it OFF, on the
/// same continent and the same threshold. Two candidate causes, and they are not the same finding:
///
/// 1. `clip_rivers_to_lakes` (drainage.rs:771) splits a segment into runs and copies
///    `strahler_order: s.strahler_order` onto every fragment — so order survives clipping but is
///    no longer a property of the fragment;
/// 2. **`full_tree = false`** extends only the MAIN STEM headward, so the tributary tree is never
///    extracted and a trunk with no *extracted* tributary is order 1 **by Strahler's own rule**.
///
/// If (2) dominates, the order is measuring EXTRACTION, not hierarchy, and block B's 41 % must not
/// be reported as an answer. The discriminator is one run with `full_tree = true`.
///
/// Run: cargo test -p ymir-core --release --test f110_comb -- --ignored full_tree --nocapture
#[test]
#[ignore]
fn f110_full_tree() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 110-B2 . Strahler under `full_tree = true`  ==========");
    for (label, knobs) in [
        ("F109 OFF", Knobs::passes(2)),
        ("F109 ON ", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
    ] {
        let f = build_field_seed(knobs, PSEED);
        let (w, h) = (f.width, f.height);
        for full in [false, true] {
            let mut dcfg = C1DrainageConfig::default();
            dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
            dcfg.thresholds.full_tree = full;
            let pre = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
            let bf = breach_monotone(&f, &pre.flow.filled, &pre.lake_map, SEA, w, h);
            let climate =
                c1_climate_placed(&bf, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
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
            let mut by: HashMap<u8, usize> = HashMap::new();
            let mut big = (0.0f32, 0u8);
            for (i, s) in dr.rivers.segments.iter().enumerate() {
                if dr.segment_kind[i] != SegmentKind::Watercourse {
                    continue;
                }
                *by.entry(s.strahler_order).or_default() += 1;
                let a = dr.segment_catchment_cells[i] * cell_km2;
                if a > big.0 {
                    big = (a, s.strahler_order);
                }
            }
            let mut o: Vec<_> = by.iter().map(|(k, v)| (*k, *v)).collect();
            o.sort();
            eprintln!(
                "   {label} · full_tree = {full:<5} ⇒ {} Watercourse · orders {:?} · biggest \
                 **{:.0} km²** at order **{}**",
                by.values().sum::<usize>(),
                o,
                big.0,
                big.1
            );
        }
    }
    eprintln!("==========  end 110-B2 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
