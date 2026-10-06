//! ADR Finding 150 — geology v1, round 1 of 3: the inventory of what Ymir already knows about the rock (the coarse
//! tectonic state, the C-2 / C-3 / C-3b closures), its maps, and whether C-3 already changes the relief. Declared in
//! `docs/reports/geology_v1/f150_inventory/f150_declared.md`. No production code.
//!
//! Run: cargo test -p ymir-core --release --test f150_geology -- --ignored --nocapture

mod common;

use common::{CANONICAL_ORIGIN, CELL_KM, DOMAIN_KM, Knobs, PSEED, SEA, aniso, build_world, pct, sorted};
use std::time::Instant;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries};
use ymir_core::tectonics_c1::closures::fracture::{FractureConfig, derive_coarse_density};
use ymir_core::tectonics_c1::closures::lithology::{LithologyConfig, build_coarse_k, stamp_volcanic_k, upscale_k_to_hd};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::place_edifices;
use ymir_core::tectonics_c1::debug_labels::derive_tectonic_labels;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

/// The equilibrium-slope floor of the témoin (`Knobs::slope_floor_abs`), as in `f126_coast.rs`.
const S_EQ: f32 = 0.024;

fn temoin() -> Knobs {
    Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }
}

/// A 1024² north-up PNG from a per-pixel colour function of HD (x, y) (internal y = 0 is south), the sea greyed.
fn save_map(path: &std::path::Path, land: &dyn Fn(usize, usize) -> bool, col: &dyn Fn(usize, usize) -> [u8; 3]) {
    let n = 1024u32;
    let img = image::RgbImage::from_fn(n, n, |px, py| {
        let x = px as usize * 8 + 4;
        let y = (n - 1 - py) as usize * 8 + 4;
        let c = col(x, y);
        if land(x, y) {
            image::Rgb(c)
        } else {
            image::Rgb([(c[0] as u16 / 4 + 60) as u8, (c[1] as u16 / 4 + 62) as u8, (c[2] as u16 / 4 + 70) as u8])
        }
    });
    img.save(path).expect("write the map");
}

/// A blue → yellow ramp for t ∈ [0, 1].
fn ramp(t: f32) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    let a = [40.0, 30.0, 110.0];
    let m = [30.0, 150.0, 140.0];
    let b = [250.0, 225.0, 60.0];
    let (p, q, u) = if t < 0.5 { (a, m, t * 2.0) } else { (m, b, (t - 0.5) * 2.0) };
    [(p[0] + (q[0] - p[0]) * u) as u8, (p[1] + (q[1] - p[1]) * u) as u8, (p[2] + (q[2] - p[2]) * u) as u8]
}

#[test]
#[ignore]
fn f150_geology() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/reports/geology_v1/f150_inventory");
    std::fs::create_dir_all(&out).expect("the report folder");
    eprintln!("\n==========  Finding 150 . geology v1, the inventory (témoin C2 /10 col; C-2, C-3, C-3b ON)  ==========");
    let wd = build_world(temoin(), None, PSEED, None);
    let f = &wd.heightmap;
    let (w, h) = (f.width, f.height);
    let n = w * h;
    let st = &wd.state;
    let (nx, ny) = (st.plate_id.nx(), st.plate_id.ny());
    let nc = nx * ny;
    let land = |x: usize, y: usize| f.data[y * w + x] > SEA;
    // the coarse cell under an HD pixel: the témoin's framing (canonical origin, size 1), nearest
    let (ox, oy) = (CANONICAL_ORIGIN[0] * nx as f64, CANONICAL_ORIGIN[1] * ny as f64);
    let coarse = |x: usize, y: usize| -> (usize, usize) {
        let sx = ox + (x as f64 + 0.5) * nx as f64 / w as f64;
        let sy = oy + (y as f64 + 0.5) * ny as f64 / h as f64;
        ((sx.floor() as i64).rem_euclid(nx as i64) as usize, (sy.floor() as i64).rem_euclid(ny as i64) as usize)
    };
    // ── I: the coarse state, filled?
    let s_vals: Vec<f32> = (0..nc).map(|k| st.s.get(k % nx, k / nx) as f32).collect();
    let age_vals: Vec<f32> = (0..nc).map(|k| st.age.get(k % nx, k / nx) as f32).collect();
    let cont: Vec<bool> = (0..nc).map(|k| st.plate_type.get(k % nx, k / nx) == PlateType::Continental).collect();
    let craton: Vec<bool> = (0..nc).map(|k| st.cratonic_mask.get(k % nx, k / nx)).collect();
    let mut ids: Vec<u16> = (0..nc).map(|k| st.plate_id.get(k % nx, k / nx)).collect();
    let pid = ids.clone();
    ids.sort_unstable();
    ids.dedup();
    let cont_age: Vec<f32> = (0..nc).filter(|&k| cont[k]).map(|k| age_vals[k]).collect();
    eprintln!("\n   I · coarse grid {nx}×{ny} (the témoin's settled state, 300 steps)");
    eprintln!(
        "   I · s (crust thickness, non-dim.): min / p10 / p50 / p90 / max {:.3} / {:.3} / {:.3} / {:.3} / {:.3} · on continental cells p50 {:.3} · on oceanic p50 {:.3}",
        pct(&sorted(s_vals.clone()), 0.0),
        pct(&sorted(s_vals.clone()), 0.1),
        pct(&sorted(s_vals.clone()), 0.5),
        pct(&sorted(s_vals.clone()), 0.9),
        pct(&sorted(s_vals.clone()), 1.0),
        pct(&sorted((0..nc).filter(|&k| cont[k]).map(|k| s_vals[k]).collect()), 0.5),
        pct(&sorted((0..nc).filter(|&k| !cont[k]).map(|k| s_vals[k]).collect()), 0.5)
    );
    eprintln!(
        "   I · age (non-dim.): min / p10 / p50 / p90 / max {:.3} / {:.3} / {:.3} / {:.3} / {:.3} · distinct values {} · continental cells with age < 1 (C-3's rift-soft) {} of {}",
        pct(&sorted(age_vals.clone()), 0.0),
        pct(&sorted(age_vals.clone()), 0.1),
        pct(&sorted(age_vals.clone()), 0.5),
        pct(&sorted(age_vals.clone()), 0.9),
        pct(&sorted(age_vals.clone()), 1.0),
        {
            let mut v: Vec<u32> = age_vals.iter().map(|a| a.to_bits()).collect();
            v.sort_unstable();
            v.dedup();
            v.len()
        },
        cont_age.iter().filter(|&&a| a < 1.0).count(),
        cont_age.len()
    );
    eprintln!(
        "   I · plate_type: continental {} · oceanic {} of {nc} · plate_id: {} distinct plates (num_plates at init {}) · cratonic_mask: {} cells ({:.1} % of the continental)",
        cont.iter().filter(|&&c| c).count(),
        cont.iter().filter(|&&c| !c).count(),
        ids.len(),
        st.num_plates,
        craton.iter().filter(|&&c| c).count(),
        100.0 * craton.iter().filter(|&&c| c).count() as f64 / cont.iter().filter(|&&c| c).count().max(1) as f64
    );
    let info = classify_boundaries(&st.plate_id, &wd.kin);
    let mut bt = [0usize; 4];
    for k in 0..nc {
        let i = match info.boundary_type.get(k % nx, k / nx) {
            BoundaryType::Internal => 0,
            BoundaryType::Convergent => 1,
            BoundaryType::Divergent => 2,
            _ => 3,
        };
        bt[i] += 1;
    }
    let lab = derive_tectonic_labels(st, &wd.kin);
    let cnt = |v: &[bool]| v.iter().filter(|&&b| b).count();
    eprintln!(
        "   I · boundaries: internal {} · convergent {} · divergent {} · transform/other {} · labels: craton {} · rift {} · subduction upper {} · slab {} · collision {} · divergent {}",
        bt[0],
        bt[1],
        bt[2],
        bt[3],
        cnt(&lab.craton),
        cnt(&lab.rift),
        cnt(&lab.subduction_upper),
        cnt(&lab.subduction_slab),
        cnt(&lab.collision),
        cnt(&lab.divergent)
    );
    // ── I: the HD fields (C-3 classes, C-3b density, the edifices)
    let lcfg = LithologyConfig { enabled: true, soft_multiplier: 10.0, volcanic_multiplier: 3.0, rift_age_threshold: 1.0 };
    let k_rift = upscale_k_to_hd(&build_coarse_k(st, &lcfg), w, h, CANONICAL_ORIGIN, 1.0);
    let edifices = place_edifices(st, &wd.kin, &WorldSeed::new(PSEED), DOMAIN_KM, &wd.volc);
    let mut k_volc = vec![1.0f32; n];
    stamp_volcanic_k(&mut k_volc, &edifices, CANONICAL_ORIGIN, 1.0, CELL_KM, w, h, &lcfg);
    let fcfg = FractureConfig { enabled: true, amplitude: 6.0, decay_km: 25.0, domain_km: DOMAIN_KM, ..Default::default() };
    let dens = upscale_k_to_hd(&derive_coarse_density(st, &wd.kin, &fcfg, None), w, h, CANONICAL_ORIGIN, 1.0);
    let soft = |c: usize| k_rift[c] > 1.5;
    let volc = |c: usize| k_volc[c] > 1.5;
    let land_n = (0..n).filter(|&c| f.data[c] > SEA).count();
    let land_soft = (0..n).filter(|&c| f.data[c] > SEA && soft(c)).count();
    let land_volc = (0..n).filter(|&c| f.data[c] > SEA && volc(c)).count();
    let dl: Vec<f32> = (0..n).filter(|&c| f.data[c] > SEA).map(|c| dens[c]).collect();
    eprintln!(
        "   I · HD land {land_n} cells · C-3 rift-soft (K×10) {land_soft} ({:.2} %) · volcaniclastic (K×3) {land_volc} ({:.2} %) · hard {:.2} % · edifices {} · craters {} · C-3b density on land p10 / p50 / p90 {:.3} / {:.3} / {:.3} (K = 1 + 6·density)",
        100.0 * land_soft as f64 / land_n as f64,
        100.0 * land_volc as f64 / land_n as f64,
        100.0 - 100.0 * (land_soft + land_volc) as f64 / land_n as f64,
        edifices.len(),
        wd.craters.len(),
        pct(&sorted(dl.clone()), 0.1),
        pct(&sorted(dl.clone()), 0.5),
        pct(&sorted(dl), 0.9)
    );
    // ── I: the maps
    let landf = |x: usize, y: usize| land(x, y);
    let (smin, smax) = (pct(&sorted(s_vals.clone()), 0.0), pct(&sorted(s_vals.clone()), 1.0));
    save_map(&out.join("map_crust_thickness_s.png"), &landf, &|x, y| {
        let (i, j) = coarse(x, y);
        ramp((s_vals[j * nx + i] - smin) / (smax - smin).max(1e-6))
    });
    let (amin, amax) = (pct(&sorted(age_vals.clone()), 0.0), pct(&sorted(age_vals.clone()), 1.0));
    save_map(&out.join("map_age.png"), &landf, &|x, y| {
        let (i, j) = coarse(x, y);
        ramp((age_vals[j * nx + i] - amin) / (amax - amin).max(1e-6))
    });
    save_map(&out.join("map_plate_type.png"), &landf, &|x, y| {
        let (i, j) = coarse(x, y);
        if cont[j * nx + i] { [200, 160, 90] } else { [60, 110, 170] }
    });
    save_map(&out.join("map_plate_id.png"), &landf, &|x, y| {
        let (i, j) = coarse(x, y);
        let p = pid[j * nx + i] as u32;
        [((p * 97) % 200 + 40) as u8, ((p * 57) % 200 + 40) as u8, ((p * 151) % 200 + 40) as u8]
    });
    save_map(&out.join("map_craton.png"), &landf, &|x, y| {
        let (i, j) = coarse(x, y);
        if craton[j * nx + i] { [205, 178, 100] } else if cont[j * nx + i] { [120, 120, 120] } else { [60, 80, 110] }
    });
    save_map(&out.join("map_tectonic_settings.png"), &landf, &|x, y| {
        let (i, j) = coarse(x, y);
        let k = j * nx + i;
        if lab.subduction_upper[k] {
            [245, 140, 30]
        } else if lab.subduction_slab[k] {
            [70, 120, 230]
        } else if lab.collision[k] {
            [230, 60, 200]
        } else if lab.divergent[k] {
            [50, 220, 220]
        } else if lab.rift[k] {
            [40, 185, 175]
        } else if lab.craton[k] {
            [205, 178, 100]
        } else if lab.continental[k] {
            [140, 140, 140]
        } else {
            [60, 80, 110]
        }
    });
    save_map(&out.join("map_c3_lithology_class.png"), &landf, &|x, y| {
        let c = y * w + x;
        if volc(c) { [175, 70, 210] } else if soft(c) { [40, 185, 175] } else { [150, 140, 120] }
    });
    save_map(&out.join("map_c3b_fracture_density.png"), &landf, &|x, y| ramp(dens[y * w + x]));
    eprintln!("   I · maps written to {}", out.display());
    // ── D: C-3 ON (A, the témoin) against OFF (B)
    let zm = |g: &ymir_core::grid::GridF32, c: usize| c1_altitude_norm_to_metres(g.data[c], &ss);
    let wb = build_world(Knobs { lithology_off: true, ..temoin() }, None, PSEED, None);
    let g = &wb.heightmap;
    let either = |c: usize| f.data[c] > SEA || g.data[c] > SEA;
    let dz: Vec<f32> = (0..n).filter(|&c| either(c)).map(|c| zm(f, c) - zm(g, c)).collect();
    let adz: Vec<f32> = dz.iter().map(|d| d.abs()).collect();
    eprintln!(
        "\n   D · Δz = z(C-3 ON) − z(C-3 OFF) on {} land cells: |Δz| p50 / p90 / p99 / max {:.2} / {:.2} / {:.1} / {:.1} m · signed p1 / p99 {:.1} / {:.1} m · cells |Δz| > 10 m {} ({:.2} %) · > 100 m {}",
        dz.len(),
        pct(&sorted(adz.clone()), 0.5),
        pct(&sorted(adz.clone()), 0.9),
        pct(&sorted(adz.clone()), 0.99),
        pct(&sorted(adz.clone()), 1.0),
        pct(&sorted(dz.clone()), 0.01),
        pct(&sorted(dz.clone()), 0.99),
        adz.iter().filter(|&&d| d > 10.0).count(),
        100.0 * adz.iter().filter(|&&d| d > 10.0).count() as f64 / adz.len().max(1) as f64,
        adz.iter().filter(|&&d| d > 100.0).count()
    );
    // where: on the soft class, on a volcanic footprint, within 5 km of either (a chamfer distance), elsewhere
    let seed_mask: Vec<bool> = (0..n).map(|c| soft(c) || volc(c)).collect();
    let mut dist = vec![f32::INFINITY; n];
    for c in 0..n {
        if seed_mask[c] {
            dist[c] = 0.0;
        }
    }
    let d1 = 1.0f32;
    let d2 = std::f32::consts::SQRT_2;
    for y in 0..h {
        for x in 0..w {
            let c = y * w + x;
            let mut v = dist[c];
            if x > 0 { v = v.min(dist[c - 1] + d1); }
            if y > 0 {
                v = v.min(dist[c - w] + d1);
                if x > 0 { v = v.min(dist[c - w - 1] + d2); }
                if x + 1 < w { v = v.min(dist[c - w + 1] + d2); }
            }
            dist[c] = v;
        }
    }
    for y in (0..h).rev() {
        for x in (0..w).rev() {
            let c = y * w + x;
            let mut v = dist[c];
            if x + 1 < w { v = v.min(dist[c + 1] + d1); }
            if y + 1 < h {
                v = v.min(dist[c + w] + d1);
                if x + 1 < w { v = v.min(dist[c + w + 1] + d2); }
                if x > 0 { v = v.min(dist[c + w - 1] + d2); }
            }
            dist[c] = v;
        }
    }
    let five_km = 5.0 / CELL_KM;
    let (mut on_soft, mut on_volc, mut near, mut else_) = (0usize, 0usize, 0usize, 0usize);
    for c in 0..n {
        if !either(c) || (zm(f, c) - zm(g, c)).abs() <= 10.0 {
            continue;
        }
        if volc(c) {
            on_volc += 1;
        } else if soft(c) {
            on_soft += 1;
        } else if dist[c] <= five_km {
            near += 1;
        } else {
            else_ += 1;
        }
    }
    let tot = (on_soft + on_volc + near + else_).max(1) as f64;
    eprintln!(
        "   D · where the |Δz| > 10 m cells are: on the rift-soft class {:.1} % · on a volcanic footprint {:.1} % · within 5 km of either {:.1} % · elsewhere {:.1} %",
        100.0 * on_soft as f64 / tot,
        100.0 * on_volc as f64 / tot,
        100.0 * near as f64 / tot,
        100.0 * else_ as f64 / tot
    );
    let relief = |g: &ymir_core::grid::GridF32| -> (f32, f32, f32) {
        let landv: Vec<bool> = g.data.iter().map(|&v| v > SEA).collect();
        let alt: Vec<f32> = (0..n).filter(|&c| landv[c]).map(|c| zm(g, c)).collect();
        let mut sig = Vec::new();
        for y in 1..h - 1 {
            for x in (1..w - 1).step_by(7) {
                let k = y * w + x;
                if !landv[k] {
                    continue;
                }
                let (mut s, mut s2) = (0f64, 0f64);
                for dy in 0..3 {
                    for dx in 0..3 {
                        let v = zm(g, (y + dy - 1) * w + x + dx - 1) as f64;
                        s += v;
                        s2 += v * v;
                    }
                }
                sig.push(((s2 / 9.0 - (s / 9.0) * (s / 9.0)).max(0.0)).sqrt() as f32);
            }
        }
        (pct(&sorted(alt), 0.5), aniso(g, &landv, 16).r8, pct(&sorted(sig), 0.5))
    };
    let (ra, r8a, sa) = relief(f);
    let (rb, r8b, sb) = relief(g);
    eprintln!(
        "   D · land altitude p50: ON {ra:.1} m · OFF {rb:.1} m ({:+.2} %) · R8: ON {r8a:.4} · OFF {r8b:.4} ({:+.4}) · local σ 3×3 p50: ON {sa:.2} m · OFF {sb:.2} m",
        100.0 * (ra - rb) / rb.abs().max(1e-6),
        r8a - r8b
    );
    // the Δz map (north up): red = higher with C-3, blue = lower
    save_map(&out.join("map_d_dz_c3_on_minus_off.png"), &landf, &|x, y| {
        let c = y * w + x;
        let d = zm(f, c) - zm(g, c);
        let t = (d.abs() / 100.0).min(1.0);
        if d >= 0.0 {
            [(200.0 + 55.0 * t) as u8, (200.0 - 170.0 * t) as u8, (200.0 - 170.0 * t) as u8]
        } else {
            [(200.0 - 170.0 * t) as u8, (200.0 - 120.0 * t) as u8, (200.0 + 55.0 * t) as u8]
        }
    });
    eprintln!("\n==========  end Finding 150 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
