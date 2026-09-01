//! DIAGNOSTIC (read-only, the shipped arid-hot H-1c export IS the verdict).
//!
//! DEFECT 1 — PHANTOM RIVERS at the head of the discharge ranking. Rivers #1..#6 run 0 m →
//! −20 m over 1 km, Strahler 1, equal width at both ends, draining 88 468 km². A 1 km
//! order-1 course draining 88 000 km² is impossible. Hypothesis to CONFIRM OR REFUTE: they
//! are the stretch of trunk crossing the floor the endorheic shrinkage EXPOSED — formerly
//! clipped away because it lay inside the lake, now dry and therefore kept.
//! The two causes call for different remedies, so they must be separated:
//!   (a) legitimate but MISCHARACTERISED — water does cross the exposed floor, but it is
//!       sheet runoff over a dry lake bed, not a watercourse fit for `rivers.json`;
//!   (b) an ORDERING problem — tracks traced against the OLD footprint, orphaned when the
//!       lake receded.
//! The decisive measurement: do the phantom cells lie on floor EXPOSED by the shrinkage
//! (old lake ∧ ¬new lake)? The OLD footprint is rebuilt deterministically here; the NEW one
//! is read from the export's `lake_mask.u32`.
//!
//! DEFECT 2 — COASTLINE COMB. Same discriminant that settled the river comb: is the barbing
//! in the HEIGHT FIELD (real serration) or only in the EXTRACTED CONTOUR (marching squares
//! on a rectangular grid — the same class of constraint as D8)?
//!
//! Run: cargo test -p ymir-core --test h1c_phantom_and_coast -- --ignored --nocapture

use std::path::Path;

use serde::Deserialize;

const DIR: &str = "../../exports/seed10481999410520546993_8192.h1_c.arid-hot.ymir";
const N: usize = 8192;
const HMIN: f32 = -5535.70458984375;
const HMAX: f32 = 4145.34228515625;
const DOMAIN_KM: f32 = 400.0;
const KM_PER_CELL: f32 = DOMAIN_KM / N as f32;

#[derive(Deserialize)]
struct Rivers {
    segments: Vec<Seg>,
}
#[derive(Deserialize)]
struct Seg {
    points: Vec<[i64; 2]>,
    strahler_order: u32,
    drainage_km2: f32,
    discharge_m3s: f32,
    width_m: f32,
    profile_m: Vec<f32>,
}

fn load_height() -> Vec<f32> {
    let b = std::fs::read(Path::new(DIR).join("height.u16")).unwrap();
    let span = HMAX - HMIN;
    b.chunks_exact(2)
        .map(|c| HMIN + (u16::from_le_bytes([c[0], c[1]]) as f32 / 65535.0) * span)
        .collect()
}

fn load_lake_mask() -> Vec<u32> {
    let b = std::fs::read(Path::new(DIR).join("lake_mask.u32")).unwrap();
    b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}

/// Rebuild the OLD (pre-shrink) lake footprint: the climate-free pre-breach drainage, which
/// is exactly what production carries before H-1c settles the endorheic basins.
fn old_footprint() -> (Vec<u32>, ymir_core::grid::GridF32) {
    use ymir_core::seed::WorldSeed;
    use ymir_core::tectonics::isostasy::IsostasyConfig;
    use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
    use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
    use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
    use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
    use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
    use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use ymir_core::tectonics_c1::kinematics::PlateKinematics;
    use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
    use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use ymir_core::terrain::upscale::FbmUpscaleConfig;

    let ss = SteinSteinParams::default();
    let run_cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state =
        init_c1_state_phase_2_r7(64, 10481999410520546993, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run_cfg, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(10481999410520546993);
    let cell_km2 = KM_PER_CELL * KM_PER_CELL;
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = FbmUpscaleConfig::c1_hd_production(N);
    cfg.amplitude_base = 0.04; // PRODUCTION amplitude (the viz overrides c1_hd_production's 0.16)
    cfg.sample_origin = [0.09375, 0.578125];
    cfg.sample_size = 1.0;
    cfg.erosion = None;
    cfg.lithology = LithologyConfig {
        enabled: true,
        soft_multiplier: 10.0,
        volcanic_multiplier: 3.0,
        rift_age_threshold: 1.0,
    };
    cfg.fracture = FractureConfig {
        enabled: true,
        amplitude: 6.0,
        decay_km: 25.0,
        domain_km: DOMAIN_KM,
        ..Default::default()
    };
    cfg.stream_power = {
        let mut sp = ymir_core::erosion::stream_power::StreamPowerConfig::relief_v3(
            cell_km2,
            ss.depth_scale_m as f32,
        );
        sp.mfd_exponent = Some(2.0);
        sp.iterations = 2;
        Some(sp)
    };
    let (up, _c) = upscale_from_c1_with_progress(
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
    );
    let raw = up.heightmap;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field = ymir_core::terrain::flow::breach_monotone(
        &raw,
        &pre.flow.filled,
        &pre.lake_map,
        0.5,
        raw.width,
        raw.height,
    );
    (pre.lake_map, field)
}

#[test]
#[ignore]
fn phantom_rivers() {
    if !Path::new(DIR).exists() {
        eprintln!("export missing — skip");
        return;
    }
    let rivers: Rivers =
        serde_json::from_slice(&std::fs::read(Path::new(DIR).join("rivers.json")).unwrap())
            .unwrap();
    let segs = &rivers.segments;
    eprintln!("\n=== DEFECT 1 — PHANTOM RIVER CENSUS (arid-hot H-1c export) ===");
    eprintln!("segments: {}", segs.len());

    // Ground the threshold in the data: for ORDER-1 segments, how does drainage_km2
    // distribute? A real headwater drains little; a phantom inherits a trunk catchment.
    let mut o1: Vec<f32> =
        segs.iter().filter(|s| s.strahler_order == 1).map(|s| s.drainage_km2).collect();
    o1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |v: &Vec<f32>, p: f32| v[((v.len() as f32 * p) as usize).min(v.len() - 1)];
    eprintln!(
        "order-1 drainage_km2  p50 {:.1} | p90 {:.0} | p99 {:.0} | max {:.0}   ({} order-1 segments)",
        q(&o1, 0.5),
        q(&o1, 0.9),
        q(&o1, 0.99),
        q(&o1, 1.0),
        o1.len()
    );
    // PHANTOM = order 1, short, and draining far more than a headwater ever should.
    let thr = q(&o1, 0.99).max(1000.0);
    let phantom: Vec<&Seg> = segs
        .iter()
        .filter(|s| s.strahler_order == 1 && s.drainage_km2 >= thr && s.points.len() <= 40)
        .collect();
    eprintln!(
        "\nPHANTOM criterion: order 1 ∧ drainage ≥ {thr:.0} km² ∧ ≤ 40 pts → {} segment(s)",
        phantom.len()
    );
    // Sort by DISCHARGE — the author sees them at the head of that ranking.
    let mut phantom = phantom;
    phantom.sort_by(|a, b| b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap());
    let below_sea =
        phantom.iter().filter(|s| s.profile_m.last().copied().unwrap_or(0.0) < 0.0).count();
    eprintln!(
        "  of which END BELOW SEA LEVEL (z_end < 0): {below_sea} / {}  ← below-sea spillway signature",
        phantom.len()
    );
    eprintln!(
        "{:>4} | {:>5} | {:>10} | {:>9} | {:>8} | {:>9} {:>9} | {:>8}",
        "#", "pts", "drain km²", "Q m³/s", "width m", "z_start", "z_end", "Δz m"
    );
    for (i, s) in phantom.iter().take(12).enumerate() {
        let (z0, z1) = (
            s.profile_m.first().copied().unwrap_or(0.0),
            s.profile_m.last().copied().unwrap_or(0.0),
        );
        eprintln!(
            "{:>4} | {:>5} | {:>10.0} | {:>9.0} | {:>8.0} | {:>9.1} {:>9.1} | {:>8.1}",
            i + 1,
            s.points.len(),
            s.drainage_km2,
            s.discharge_m3s,
            s.width_m,
            z0,
            z1,
            z1 - z0
        );
    }

    // ── THE DECISIVE MEASUREMENT: are their cells on floor EXPOSED by the shrinkage?
    let new_lake = load_lake_mask();
    let (old_lake, _field) = old_footprint();
    let exposed: Vec<bool> = (0..N * N).map(|k| old_lake[k] != 0 && new_lake[k] == 0).collect();
    let n_exposed = exposed.iter().filter(|&&b| b).count();
    eprintln!(
        "\nfootprints: OLD lake {} cells | NEW lake {} cells | EXPOSED floor {} cells ({:.0} km²)",
        old_lake.iter().filter(|&&v| v != 0).count(),
        new_lake.iter().filter(|&&v| v != 0).count(),
        n_exposed,
        n_exposed as f32 * KM_PER_CELL * KM_PER_CELL
    );
    let frac_on = |s: &Seg| {
        let hit = s
            .points
            .iter()
            .filter(|p| {
                let (x, y) =
                    (p[0].rem_euclid(N as i64) as usize, p[1].rem_euclid(N as i64) as usize);
                exposed[y * N + x]
            })
            .count();
        100.0 * hit as f32 / s.points.len().max(1) as f32
    };
    let mut on_exposed = 0usize;
    for s in &phantom {
        if frac_on(s) > 50.0 {
            on_exposed += 1;
        }
    }
    eprintln!(
        "PHANTOMS whose cells lie on EXPOSED FLOOR (>50 % of points): {on_exposed} / {}",
        phantom.len()
    );
    let mean_frac: f32 =
        phantom.iter().map(|s| frac_on(s)).sum::<f32>() / phantom.len().max(1) as f32;
    eprintln!("  mean fraction of phantom points on exposed floor: {mean_frac:.0} %");
    // Control: the same fraction over NORMAL high-order rivers (should be ~0).
    let normal: Vec<&Seg> = segs.iter().filter(|s| s.strahler_order >= 3).collect();
    let ctrl: f32 = normal.iter().map(|s| frac_on(s)).sum::<f32>() / normal.len().max(1) as f32;
    eprintln!("  control — same fraction over order ≥3 trunks: {ctrl:.1} %");
    eprintln!(
        "\n(If the phantoms sit on exposed floor and the trunks do not, the segments are the trunk's\n crossing of the drained lake bed, kept by the clip because those cells are no longer lake.)"
    );
}

#[test]
#[ignore]
fn coastline_comb() {
    if !Path::new(DIR).exists() {
        eprintln!("export missing — skip");
        return;
    }
    let h = load_height();
    // The coastline contour as exported (marching squares at 0 m).
    let gj: serde_json::Value =
        serde_json::from_slice(&std::fs::read(Path::new(DIR).join("coastline.geojson")).unwrap())
            .unwrap();
    // Collect vertices of the longest linestring.
    let mut best: Vec<[f64; 2]> = Vec::new();
    fn walk(v: &serde_json::Value, best: &mut Vec<[f64; 2]>) {
        if let Some(arr) = v.as_array() {
            // a linestring is an array of [x,y] pairs
            if arr.len() > 2
                && arr[0].as_array().map(|p| p.len() == 2 && p[0].is_number()) == Some(true)
            {
                if arr.len() > best.len() {
                    *best = arr
                        .iter()
                        .filter_map(|p| {
                            let q = p.as_array()?;
                            Some([q[0].as_f64()?, q[1].as_f64()?])
                        })
                        .collect();
                }
                return;
            }
            for x in arr {
                walk(x, best);
            }
        } else if let Some(o) = v.as_object() {
            for (_k, x) in o {
                walk(x, best);
            }
        }
    }
    walk(&gj, &mut best);
    eprintln!("\n=== DEFECT 2 — COASTLINE: terrain serration or contour artefact? ===");
    eprintln!("longest coastline ring: {} vertices", best.len());
    if best.len() < 50 {
        eprintln!("not enough vertices — skip");
        return;
    }

    // (1) The CONTOUR's own zigzag: axial concentration of its step directions, and the
    //     turn angle between consecutive steps. Marching squares on a rectangular grid
    //     produces steps on grid axes/diagonals and alternating ±90° turns (the barbs).
    let mut hist = [0u32; 12];
    let (mut sc, mut ss) = (0.0f64, 0.0f64);
    let mut turns = Vec::new();
    let mut prev: Option<(f64, f64)> = None;
    let mut steplen = Vec::new();
    for w2 in best.windows(2) {
        let (dx, dy) = (w2[1][0] - w2[0][0], w2[1][1] - w2[0][1]);
        let m = (dx * dx + dy * dy).sqrt();
        if m < 1e-9 {
            continue;
        }
        steplen.push(m);
        let th = dy.atan2(dx);
        sc += (2.0 * th).cos();
        ss += (2.0 * th).sin();
        let b = (th.rem_euclid(std::f64::consts::PI)) / (std::f64::consts::PI / 12.0);
        hist[(b as usize).min(11)] += 1;
        if let Some((px, py)) = prev {
            let a0 = py.atan2(px);
            let mut d = (th - a0).to_degrees();
            while d > 180.0 {
                d -= 360.0;
            }
            while d < -180.0 {
                d += 360.0;
            }
            turns.push(d.abs());
        }
        prev = Some((dx, dy));
    }
    let nn = steplen.len().max(1) as f64;
    let r_axial = ((sc / nn).powi(2) + (ss / nn).powi(2)).sqrt();
    turns.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let tq = |p: f32| turns[((turns.len() as f32 * p) as usize).min(turns.len() - 1)];
    let sharp = 100.0 * turns.iter().filter(|&&t| t > 80.0).count() as f32 / turns.len() as f32;
    eprintln!(
        "CONTOUR: axial R {r_axial:.2} | 12-bin hist {hist:?}\n  turn angle p50 {:.0}° p90 {:.0}° | turns > 80° : {sharp:.0} %  | mean step {:.2} cells",
        tq(0.5),
        tq(0.9),
        steplen.iter().sum::<f64>() / nn
    );

    // (2) The TERRAIN under the coast: is it serrated, or is it a near-flat shelf where any
    //     discretisation zigzags? Sample the elevation profile perpendicular to the local
    //     contour direction, and the coastal slope.
    let cell_m = KM_PER_CELL * 1000.0;
    let at = |x: i64, y: i64| -> f32 {
        h[(y.rem_euclid(N as i64) as usize) * N + (x.rem_euclid(N as i64) as usize)]
    };
    let mut prof = [0.0f64; 17];
    let mut slopes = Vec::new();
    let mut cnt = 0u64;
    for i in (2..best.len() - 2).step_by(7) {
        let (dx, dy) = (best[i + 1][0] - best[i - 1][0], best[i + 1][1] - best[i - 1][1]);
        let m = (dx * dx + dy * dy).sqrt();
        if m < 1e-9 {
            continue;
        }
        let (nx, ny) = (-dy / m, dx / m); // ⊥ the contour
        let (cx, cy) = (best[i][0], best[i][1]);
        for (bi, t) in (-8i64..=8).enumerate() {
            prof[bi] +=
                at((cx + nx * t as f64).round() as i64, (cy + ny * t as f64).round() as i64) as f64;
        }
        // local coastal slope in degrees over ±4 cells across the shore
        let za = at((cx + nx * 4.0).round() as i64, (cy + ny * 4.0).round() as i64);
        let zb = at((cx - nx * 4.0).round() as i64, (cy - ny * 4.0).round() as i64);
        slopes.push((((za - zb).abs() / (8.0 * cell_m)).atan()).to_degrees());
        cnt += 1;
    }
    slopes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let sq = |p: f32| slopes[((slopes.len() as f32 * p) as usize).min(slopes.len() - 1)];
    eprintln!(
        "TERRAIN ⊥ the shore, mean elevation (m) at offsets −8..+8 cells (×{:.0} m):",
        cell_m
    );
    let line: Vec<String> =
        prof.iter().map(|v| format!("{:>6.1}", v / cnt.max(1) as f64)).collect();
    eprintln!("  {}", line.join(" "));
    eprintln!(
        "  coastal slope across the shore: p10 {:.2}° | p50 {:.2}° | p90 {:.2}°",
        sq(0.1),
        sq(0.5),
        sq(0.9)
    );
    // ── THE DECISIVE LINK: do the BARBS concentrate where the shore is FLAT? Pair each
    //    vertex's turn angle with the local coastal slope and bucket.
    let mut pairs: Vec<(f32, f32)> = Vec::new(); // (slope_deg, turn_deg)
    for i in (2..best.len() - 2).step_by(3) {
        let (dxa, dya) = (best[i][0] - best[i - 1][0], best[i][1] - best[i - 1][1]);
        let (dxb, dyb) = (best[i + 1][0] - best[i][0], best[i + 1][1] - best[i][1]);
        let (ma, mb) = ((dxa * dxa + dya * dya).sqrt(), (dxb * dxb + dyb * dyb).sqrt());
        if ma < 1e-9 || mb < 1e-9 {
            continue;
        }
        let mut d = (dyb.atan2(dxb) - dya.atan2(dxa)).to_degrees();
        while d > 180.0 {
            d -= 360.0;
        }
        while d < -180.0 {
            d += 360.0;
        }
        let (nx, ny) = (-dya / ma, dxa / ma);
        let (cx, cy) = (best[i][0], best[i][1]);
        let za = at((cx + nx * 4.0).round() as i64, (cy + ny * 4.0).round() as i64);
        let zb = at((cx - nx * 4.0).round() as i64, (cy - ny * 4.0).round() as i64);
        let sl = (((za - zb).abs() / (8.0 * cell_m)).atan()).to_degrees();
        pairs.push((sl, d.abs() as f32));
    }
    eprintln!("\nBARBS vs SHORE SLOPE (the decisive link):");
    eprintln!(
        "{:>16} | {:>8} | {:>12} | {:>12}",
        "shore slope", "vertices", "mean turn °", "turns>80° %"
    );
    for (lo, hi, lbl) in [
        (0.0f32, 0.5f32, "  < 0.5° (flat)"),
        (0.5, 2.0, "  0.5 – 2°"),
        (2.0, 5.0, "  2 – 5°"),
        (5.0, 15.0, "  5 – 15°"),
        (15.0, 1e9, "  > 15° (steep)"),
    ] {
        let sel: Vec<f32> =
            pairs.iter().filter(|(s, _)| *s >= lo && *s < hi).map(|(_, t)| *t).collect();
        if sel.is_empty() {
            continue;
        }
        let mean = sel.iter().sum::<f32>() / sel.len() as f32;
        let sharp = 100.0 * sel.iter().filter(|&&t| t > 80.0).count() as f32 / sel.len() as f32;
        eprintln!("{lbl:>16} | {:>8} | {mean:>12.1} | {sharp:>12.1}", sel.len());
    }
    eprintln!(
        "\n(If sharp turns concentrate on the FLAT shore buckets while steep shores stay smooth, the barbing\n is the 0 m contour wandering on ground with no gradient — a MARCHING-SQUARES artefact, not a terrain\n serration. Remedy is geometric (contour smoothing / sub-cell interpolation), far lighter than C-4.)"
    );
}
