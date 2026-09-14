//! ADR 0001 Finding 83 block C1 — **one budget, labelled**, and the PAIRED decomposition that
//! Finding 82 owed.
//!
//! Finding 81 read −5.5 % and Finding 82 read −7.5 %. They are not a contradiction: they are two
//! STAGES (breached / eroded) and therefore two land masks and two populations. Both are measured
//! here, side by side, so the retained figure can be named instead of chosen.
//!
//! And the confound Finding 82 named rather than buried is closed: the distance and altitude
//! bands were RECOMPUTED on each field, so a cell could change band when the coast moved by
//! 290 895 cells. Here the bands are assigned ONCE, on the delivered field, and the population is
//! split three ways — cells that are land in BOTH (the paired core), cells CEDED (land delivered,
//! sea bounded) and cells GAINED — so the three add up to the unpaired delta by construction.
//!
//! Run: cargo test -p ymir-core --release --test f83_budget -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, c1_drainage_windowed, potential_evaporation_mm, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;

/// Chebyshev distance to the land/sea boundary, in cells.
fn dist_to_coast(land: &[bool], w: usize) -> Vec<u16> {
    let n = w * w;
    let mut d = vec![u16::MAX; n];
    let mut q: Vec<u32> = Vec::new();
    for k in 0..n {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let mut edge = false;
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
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
            d[k] = 0;
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
                if d[nk] == u16::MAX {
                    d[nk] = d[k] + 1;
                    q.push(nk as u32);
                }
            }
        }
    }
    d
}

/// Runoff surplus `(P − PE)⁺ · cell_km²` for every cell, in km²·mm/yr.
fn surplus(f: &GridF32, ss: &SteinSteinParams, lat: f32, span: f32) -> Vec<f32> {
    let climate = c1_climate_placed(f, ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
    let cell_km2 = CELL_KM * CELL_KM;
    (0..f.data.len())
        .map(|k| {
            let p = precip_mm_per_year(climate.precipitation.data[k]);
            let pe = potential_evaporation_mm(climate.temperature.data[k]);
            (p - pe).max(0.0) * cell_km2
        })
        .collect()
}

/// Breach the raw eroded field the way production does before the drainage reads it.
fn breached(raw: &GridF32, ss: &SteinSteinParams) -> GridF32 {
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let d0 = c1_drainage_windowed(raw, None, &dcfg, ss, DOMAIN_KM);
    breach_monotone(raw, &d0.flow.filled, &d0.lake_map, SEA, raw.width, raw.height)
}

#[test]
#[ignore]
fn f83_budget() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 83 · C1 — one budget, and the PAIRED decomposition  =====");

    let del_e = build_field(Knobs::pre83());
    let bnd_e = build_field(Knobs::shipped());
    let (w, n) = (del_e.width, del_e.data.len());
    let del_b = breached(&del_e, &ss);
    let bnd_b = breached(&bnd_e, &ss);

    let m3s = |v: f64| runoff_km2_to_m3s(v as f32) as f64;

    for (stage, del, bnd) in [
        ("ERODED (Finding 82's stage)", &del_e, &bnd_e),
        ("BREACHED (what the drainage reads, Finding 81's stage)", &del_b, &bnd_b),
    ] {
        eprintln!("\n╔═══ stage: {stage} ═══╗");
        let land_d: Vec<bool> = del.data.iter().map(|&v| v > SEA).collect();
        let land_b: Vec<bool> = bnd.data.iter().map(|&v| v > SEA).collect();
        let core: Vec<usize> = (0..n).filter(|&k| land_d[k] && land_b[k]).collect();
        let ceded: Vec<usize> = (0..n).filter(|&k| land_d[k] && !land_b[k]).collect();
        let gained: Vec<usize> = (0..n).filter(|&k| !land_d[k] && land_b[k]).collect();
        eprintln!(
            "   population: land delivered {} | land bounded {} | PAIRED CORE {} | CEDED {} | \
             GAINED {}",
            land_d.iter().filter(|&&l| l).count(),
            land_b.iter().filter(|&&l| l).count(),
            core.len(),
            ceded.len(),
            gained.len()
        );

        // bands assigned ONCE, on the DELIVERED field — the "before" of the comparison
        let dist = dist_to_coast(&land_d, w);
        let band_d: Vec<u8> = (0..n)
            .map(|k| {
                let dkm = dist[k] as f32 * CELL_KM;
                if dkm < 1.0 {
                    0
                } else if dkm < 5.0 {
                    1
                } else if dkm < 20.0 {
                    2
                } else {
                    3
                }
            })
            .collect();
        let band_a: Vec<u8> = (0..n)
            .map(|k| {
                let alt = (del.data[k] - SEA) * n2m;
                if alt < 10.0 {
                    0
                } else if alt < 100.0 {
                    1
                } else if alt < 500.0 {
                    2
                } else {
                    3
                }
            })
            .collect();

        for (bed, lat, span) in [("humid", 45.0f32, 40.0f32), ("arid-hot", 25.0, 10.0)] {
            let sd = surplus(del, &ss, lat, span);
            let sb = surplus(bnd, &ss, lat, span);
            let tot_d: f64 = (0..n).filter(|&k| land_d[k]).map(|k| sd[k] as f64).sum();
            let tot_b: f64 = (0..n).filter(|&k| land_b[k]).map(|k| sb[k] as f64).sum();
            let core_d: f64 = core.iter().map(|&k| sd[k] as f64).sum();
            let core_b: f64 = core.iter().map(|&k| sb[k] as f64).sum();
            let ced: f64 = ceded.iter().map(|&k| sd[k] as f64).sum();
            let gai: f64 = gained.iter().map(|&k| sb[k] as f64).sum();
            let d_tot = m3s(tot_b) - m3s(tot_d);
            eprintln!(
                "\n   ── {bed} ──\n      UNPAIRED budget {:.1} → {:.1} m³/s (**{:+.1} %**)",
                m3s(tot_d),
                m3s(tot_b),
                100.0 * (tot_b - tot_d) / tot_d
            );
            eprintln!(
                "      decomposition of Δ = {d_tot:+.2} m³/s:  paired core {:+.2}  ·  CEDED \
                 {:+.2}  ·  GAINED {:+.2}   (sum {:+.2}, residual {:.4})",
                m3s(core_b) - m3s(core_d),
                -m3s(ced),
                m3s(gai),
                (m3s(core_b) - m3s(core_d)) - m3s(ced) + m3s(gai),
                ((m3s(core_b) - m3s(core_d)) - m3s(ced) + m3s(gai) - d_tot).abs()
            );
            eprintln!(
                "      the ceded cells carry {:.1} % of |Δ|; the paired core {:.1} %",
                100.0 * m3s(ced) / d_tot.abs().max(1e-9),
                100.0 * (m3s(core_b) - m3s(core_d)).abs() / d_tot.abs().max(1e-9)
            );
            eprintln!(
                "      {:<24} {:>12} {:>12} {:>12} {:>12}",
                "band (fixed, delivered)", "delivered", "bounded", "Δ m³/s", "share of Δcore"
            );
            let dcore = m3s(core_b) - m3s(core_d);
            for (lbl, sel) in
                [("dist 0-1 km", 0u8), ("dist 1-5 km", 1), ("dist 5-20 km", 2), ("dist > 20 km", 3)]
            {
                let a: f64 =
                    core.iter().filter(|&&k| band_d[k] == sel).map(|&k| sd[k] as f64).sum();
                let b: f64 =
                    core.iter().filter(|&&k| band_d[k] == sel).map(|&k| sb[k] as f64).sum();
                eprintln!(
                    "      {lbl:<24} {:>12.2} {:>12.2} {:>12.2} {:>11.1} %",
                    m3s(a),
                    m3s(b),
                    m3s(b) - m3s(a),
                    100.0 * (m3s(b) - m3s(a)) / dcore
                );
            }
            for (lbl, sel) in
                [("alt < 10 m", 0u8), ("alt 10-100 m", 1), ("alt 100-500 m", 2), ("alt > 500 m", 3)]
            {
                let a: f64 =
                    core.iter().filter(|&&k| band_a[k] == sel).map(|&k| sd[k] as f64).sum();
                let b: f64 =
                    core.iter().filter(|&&k| band_a[k] == sel).map(|&k| sb[k] as f64).sum();
                eprintln!(
                    "      {lbl:<24} {:>12.2} {:>12.2} {:>12.2} {:>11.1} %",
                    m3s(a),
                    m3s(b),
                    m3s(b) - m3s(a),
                    100.0 * (m3s(b) - m3s(a)) / dcore
                );
            }
            // where the ceded water was, by the same fixed bands
            let mut cb = [0.0f64; 4];
            for &k in &ceded {
                cb[band_d[k] as usize] += sd[k] as f64;
            }
            eprintln!(
                "      CEDED water by distance band: 0-1 km {:.2} | 1-5 {:.2} | 5-20 {:.2} | >20 \
                 {:.2} m³/s",
                m3s(cb[0]),
                m3s(cb[1]),
                m3s(cb[2]),
                m3s(cb[3])
            );
        }
    }
}
