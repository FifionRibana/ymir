//! ADR Finding 152 — geology v1, round 3 of 3: the build on the témoin and its gates (G-blocks, G-circles, G-nonempty,
//! G-placers, G-veins with its control, the cost, the export sizes), d_ref measured, the maps and the views' captures.
//! G-field and G-rules are permanent tests (`geology::history`, `geology::zoning`) plus the 6/6 guard. Declared in
//! `docs/reports/geology_v1/f152_build/f152_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f152_geology -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, build_world, pct, sorted, viz_hd_lakes_on};
use std::collections::HashSet;
use std::time::Instant;
use ymir_core::climate::precipitation::precip_mm_per_year;
use ymir_core::geology::rocks::{self, ROCK_CLASSES};
use ymir_core::geology::rules::{DEFAULT_RULES_TOML, parse_rules};
use ymir_core::geology::{DensitySource, GeologyInputs, build_geology, export::export_files, render, zone};
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const S_EQ: f32 = 0.024;

fn temoin() -> Knobs {
    Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }
}

/// The share of `cells` (x, y on a grid of `w`) lying within `tol` (coarse units) of a coarse grid line (a coarse
/// coordinate whose fractional part is near 0 or 0.5), either axis.
fn on_grid_lines(cells: &[(usize, usize)], w: usize, origin: [f64; 2], size: f64, nx: usize) -> f64 {
    let per = nx as f64 * size / w as f64; // coarse units per cell
    let tol = per;
    let near = |c: f64| {
        let f = c - c.floor();
        f.min(1.0 - f).min((f - 0.5).abs()) < tol
    };
    let hit = cells.iter().filter(|&&(x, y)| near(origin[0] * nx as f64 + x as f64 * per) || near(origin[1] * nx as f64 + y as f64 * per)).count();
    hit as f64 / cells.len().max(1) as f64
}

/// Q = 4πA/P² of a 4-connected patch (P = its 4-neighbour edges to other cells).
fn isoperimetric(cells: &[usize], inside: &dyn Fn(usize, usize) -> bool, w: usize) -> f64 {
    let mut p = 0usize;
    for &c in cells {
        let (x, y) = (c % w, c / w);
        for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
            let (xx, yy) = (x as i64 + dx, y as i64 + dy);
            if xx < 0 || yy < 0 || !inside(xx as usize, yy as usize) {
                p += 1;
            }
        }
    }
    4.0 * std::f64::consts::PI * cells.len() as f64 / (p as f64 * p as f64).max(1.0)
}

/// Q of a rasterised disc of `area` cells.
fn disc_q(area: usize) -> f64 {
    let r = (area as f64 / std::f64::consts::PI).sqrt();
    let n = (2.0 * r).ceil() as usize + 4;
    let c = n as f64 / 2.0;
    let inside = |x: usize, y: usize| ((x as f64 + 0.5 - c).powi(2) + (y as f64 + 0.5 - c).powi(2)).sqrt() <= r;
    let cells: Vec<usize> = (0..n * n).filter(|&k| inside(k % n, k / n)).collect();
    isoperimetric(&cells, &|x, y| x < n && y < n && inside(x, y), n)
}

/// A 1024² north-up PNG of a south-first RGBA buffer (one pixel of every 8).
fn save_rgba(path: &std::path::Path, rgba: &[u8], w: usize, h: usize) {
    let n = 1024u32;
    let f = w / n as usize;
    let img = image::RgbImage::from_fn(n, n, |px, py| {
        let x = px as usize * f + f / 2;
        let y = (h / f - 1 - py as usize) * f + f / 2;
        let k = (y * w + x) * 4;
        image::Rgb([rgba[k], rgba[k + 1], rgba[k + 2]])
    });
    img.save(path).expect("write the capture");
}

#[test]
#[ignore]
fn f152_geology() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/reports/geology_v1/f152_build");
    std::fs::create_dir_all(&out).expect("the report folder");
    eprintln!("\n==========  Finding 152 . geology v1, the build and its gates (témoin C2 /10 col)  ==========");
    let wd = build_world(temoin(), None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, temoin(), PSEED, 45.0, 40.0);
    let f = &v.conditioned;
    let (w, h) = (f.width, f.height);
    let z_m: Vec<f32> = f.data.iter().map(|&x| c1_altitude_norm_to_metres(x, &ss)).collect();
    let precip_mm: Vec<f32> = v.precipitation.data.iter().map(|&p| precip_mm_per_year(p)).collect();
    let run = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0 / 64.0, dy: 1.0 / 64.0, iso_config: IsostasyConfig::c1_default(), drainage_max_distance: 30 };
    let tb = Instant::now();
    let product = build_geology(&GeologyInputs {
        seed: PSEED,
        grid: 64,
        init: &Phase2InitParams::default(),
        run: &run,
        closures: &C1Closures::default(),
        field: f,
        z_m: &z_m,
        cell_km: CELL_KM,
        drainage: &v.drainage,
        wetland: &v.wetland,
        temp_c: &v.temperature.data,
        precip_mm: &precip_mm,
        sample_origin: wd.cfg.sample_origin,
        sample_size: wd.cfg.sample_size,
        domain_km: 400.0,
        volcanism: Some(&wd.volc),
    });
    let t_build = tb.elapsed().as_secs_f64();
    let ctx = &product.context;
    let (w4, h4) = (ctx.w4, ctx.h4);
    let t = product.timings;
    eprintln!(
        "   build {t_build:.2} s: history {:.2} · rocks {:.2} · context {:.2} · fossil belt {} coarse cells · active belt {} · merges {} · edifices {}",
        t.history_s,
        t.rocks_s,
        t.context_s,
        product.fossil.iter().filter(|&&b| b).count(),
        product.active.iter().filter(|&&b| b).count(),
        product.history.merges.len(),
        product.edifices.len()
    );
    // the rock shares
    let land_n = product.rocks.iter().filter(|&&r| r != rocks::NONE).count();
    let mut counts = [0usize; 9];
    for &r in &product.rocks {
        counts[r as usize] += 1;
    }
    eprintln!("   rocks · share of the non-water cells:");
    for c in &ROCK_CLASSES[1..] {
        eprintln!("      {} ({}): {:.2} %", c.name_fr, c.key, 100.0 * counts[c.id as usize] as f64 / land_n as f64);
    }
    // d_ref, MEASURED: the p90 of d over the land ¼ cells of the belts and the arc zone
    let zone_cells: Vec<f32> = (0..w4 * h4).filter(|&q| ctx.rock[q] != rocks::NONE && (ctx.belt_fossil[q] || ctx.belt_active[q] || ctx.arc_km[q] <= 30.0)).map(|q| ctx.density[q]).collect();
    let d_ref = (pct(&sorted(zone_cells.clone()), 0.9) * 1000.0).round() / 1000.0;
    let mut rules = parse_rules(DEFAULT_RULES_TOML).expect("the shipped rules");
    let file_d_ref = rules.modulation.structural.map(|m| m.d_ref);
    if let Some(m) = rules.modulation.structural.as_mut() {
        m.d_ref = d_ref;
    }
    eprintln!("   d_ref MEASURED = {d_ref:.3} (the p90 of d over {} belt / arc ¼ cells; p50 {:.3}) · the rules file says {file_d_ref:?} · g_min {:?}", zone_cells.len(), pct(&sorted(zone_cells), 0.5), rules.modulation.structural.map(|m| m.g_min));
    // the zoning
    let tz = Instant::now();
    let z = zone(ctx, &rules, DensitySource::Structural);
    let t_zone = tz.elapsed().as_secs_f64();
    let zc = zone(ctx, &rules, DensitySource::C3bOnly);
    eprintln!("   COST · the geology stage {:.2} s (build {t_build:.2} + zoning {t_zone:.2}) · re-zoning alone {t_zone:.2} s, against 249.8 s", t_build + t_zone);
    // ── G-nonempty
    let cell4_km2 = (ctx.cell4_km * ctx.cell4_km) as f64;
    eprintln!("\n   G-nonempty · per resource, the land ¼ cells with favourability ≥ 10 (area), max, placer HD cells:");
    for r in &z.resources {
        let n10 = (0..w4 * h4).filter(|&q| ctx.rock[q] != rocks::NONE && r.fav[q] >= 10).count();
        eprintln!("      {:<22} {:>7} cells ({:>7.1} km²) · max {:>3} · placer cells {}", r.id, n10, n10 as f64 * cell4_km2, r.fav.iter().max().copied().unwrap_or(0), r.placer_cells);
    }
    // ── G-blocks: the rock contacts and the favourability edges against the coarse grid lines
    let (nx, origin, size) = (64usize, wd.cfg.sample_origin, wd.cfg.sample_size);
    let mut contacts: Vec<(usize, usize)> = Vec::new();
    let mut coarse_contacts: Vec<(usize, usize)> = Vec::new();
    let mut land_cells: Vec<(usize, usize)> = Vec::new();
    let coarse_class = |r: u8| matches!(r, rocks::CRATON | rocks::BELT | rocks::RIFT_FILL | rocks::BASEMENT);
    for y in 0..h - 1 {
        for x in 0..w - 1 {
            let c = y * w + x;
            let a = product.rocks[c];
            if a == rocks::NONE {
                continue;
            }
            if (x + y) % 16 == 0 {
                land_cells.push((x, y));
            }
            for b in [product.rocks[c + 1], product.rocks[c + w]] {
                if b != rocks::NONE && b != a {
                    contacts.push((x, y));
                    if coarse_class(a) && coarse_class(b) {
                        coarse_contacts.push((x, y));
                    }
                    break;
                }
            }
        }
    }
    let ctrl = on_grid_lines(&land_cells, w, origin, size, nx);
    let all = on_grid_lines(&contacts, w, origin, size, nx);
    let cc = on_grid_lines(&coarse_contacts, w, origin, size, nx);
    let mut fav_edges: Vec<(usize, usize)> = Vec::new();
    let mut land4: Vec<(usize, usize)> = Vec::new();
    for y in 0..h4 - 1 {
        for x in 0..w4 - 1 {
            let q = y * w4 + x;
            if ctx.rock[q] == rocks::NONE {
                continue;
            }
            land4.push((x, y));
            if z.resources.iter().any(|r| r.fav[q].abs_diff(r.fav[q + 1]) >= 10 || r.fav[q].abs_diff(r.fav[q + w4]) >= 10) {
                fav_edges.push((x, y));
            }
        }
    }
    let ctrl4 = on_grid_lines(&land4, w4, origin, size, nx);
    let fe = on_grid_lines(&fav_edges, w4, origin, size, nx);
    eprintln!(
        "\n   G-blocks · share within 1 cell of a coarse grid line: rock contacts {:.4} ({} cells) · coarse-sourced contacts (craton / belt / rift / basement) {:.4} ({}) · control (land) {:.4} → ratios **{:.2}** / **{:.2}** · favourability edges ¼ {:.4} ({}) · control {:.4} → **{:.2}** · pass ≤ 1.5",
        all,
        contacts.len(),
        cc,
        coarse_contacts.len(),
        ctrl,
        all / ctrl.max(1e-9),
        cc / ctrl.max(1e-9),
        fe,
        fav_edges.len(),
        ctrl4,
        fe / ctrl4.max(1e-9)
    );
    // ── G-circles: each volcanic patch against a disc of its area
    let volc = |c: usize| matches!(product.rocks[c], rocks::BASALTIC_VOLCANIC | rocks::ARC_VOLCANIC);
    let mut seen = vec![false; w * h];
    let mut patches: Vec<(usize, f64, f64)> = Vec::new();
    for c0 in 0..w * h {
        if seen[c0] || !volc(c0) {
            continue;
        }
        let mut stack = vec![c0];
        seen[c0] = true;
        let mut cells = Vec::new();
        while let Some(c) = stack.pop() {
            cells.push(c);
            let (x, y) = (c % w, c / w);
            for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let (xx, yy) = (x as i64 + dx, y as i64 + dy);
                if xx < 0 || yy < 0 || xx as usize >= w || yy as usize >= h {
                    continue;
                }
                let q = yy as usize * w + xx as usize;
                if !seen[q] && volc(q) {
                    seen[q] = true;
                    stack.push(q);
                }
            }
        }
        if cells.len() < 20 {
            continue;
        }
        let q = isoperimetric(&cells, &|x, y| x < w && y < h && volc(y * w + x), w);
        let qd = disc_q(cells.len());
        patches.push((cells.len(), q, q / qd));
    }
    patches.sort_by(|a, b| b.0.cmp(&a.0));
    let worst = patches.iter().map(|p| p.2).fold(0f64, f64::max);
    eprintln!("\n   G-circles · {} volcanic patches (≥ 20 cells) · Q / Q_disc max **{worst:.3}** (pass ≤ 0.8 for every patch) · failing {}", patches.len(), patches.iter().filter(|p| p.2 > 0.8).count());
    for (n, q, r) in &patches {
        eprintln!("      {n:>8} cells · Q {q:.3} · Q / Q_disc {r:.3}{}", if *r > 0.8 { "  ← FAILS" } else { "" });
    }
    // ── G-placers: the ¼ cells a placer raised all hold a river cell of the rule's area
    let parsed = &rules;
    for (k, res) in parsed.resource.iter().enumerate() {
        let Some(pl) = &res.placer else { continue };
        let mut no_placer = rules.clone();
        no_placer.resource[k].placer = None;
        let base = zone(ctx, &no_placer, DensitySource::Structural);
        let raised: Vec<usize> = (0..w4 * h4).filter(|&q| z.resources[k].fav[q] > base.resources[k].fav[q]).collect();
        let river4: HashSet<u32> = ctx.rivers.iter().filter(|rc| rc.area_km2 >= pl.river_min_area_km2).map(|rc| rc.q4).collect();
        let on = raised.iter().filter(|&&q| river4.contains(&(q as u32))).count();
        eprintln!("   G-placers · {}: {} ¼ cells raised by the placer · on a river of ≥ {} km²: **{} ({:.1} %)** · HD placer cells {}", res.id, raised.len(), pl.river_min_area_km2, on, 100.0 * on as f64 / raised.len().max(1) as f64, z.resources[k].placer_cells);
    }
    // ── G-veins: the modulated gold in the belts, and the control (C-3b alone)
    let gi = z.resources.iter().position(|r| r.id == "gold").expect("gold");
    let ti = z.resources.iter().position(|r| r.id == "tin").expect("tin");
    for (name, belt) in [("all belts", (0..w4 * h4).filter(|&q| ctx.rock[q] == rocks::BELT).collect::<Vec<_>>()), ("fossil belts", (0..w4 * h4).filter(|&q| ctx.rock[q] == rocks::BELT && ctx.belt_fossil[q]).collect::<Vec<_>>())] {
        for (rname, ri) in [("gold", gi), ("tin", ti)] {
            let vals: Vec<f32> = belt.iter().map(|&q| z.resources[ri].fav[q] as f32).collect();
            let ctl: Vec<f32> = belt.iter().map(|&q| zc.resources[ri].fav[q] as f32).collect();
            let mx = vals.iter().copied().fold(0f32, f32::max);
            let half = vals.iter().filter(|&&x| x >= 0.5 * mx && mx > 0.0).count();
            let mean = vals.iter().sum::<f32>() / vals.len().max(1) as f32;
            let mean_c = ctl.iter().sum::<f32>() / ctl.len().max(1) as f32;
            eprintln!(
                "   G-veins · {rname} on {name} ({} ¼ cells): p10 / p50 / p90 {:.0} / {:.0} / {:.0} · max {mx:.0} · share ≥ half the max **{:.1} %** · mean {mean:.1} · control C-3b alone: mean {mean_c:.1} (**{:+.1} %**), p50 {:.0}",
                belt.len(),
                pct(&sorted(vals.clone()), 0.1),
                pct(&sorted(vals.clone()), 0.5),
                pct(&sorted(vals.clone()), 0.9),
                100.0 * half as f64 / belt.len().max(1) as f64,
                100.0 * (mean_c - mean) / mean.max(1e-6),
                pct(&sorted(ctl), 0.5)
            );
        }
    }
    // ── the export
    let sha = ymir_core::geology::sha256_hex(DEFAULT_RULES_TOML.as_bytes());
    let files = export_files(&product.rocks, product.w, product.h, &z, &rules, &sha);
    let fav_bytes: usize = files.iter().filter(|(n, _)| n.starts_with("favorabilite_")).map(|(_, b)| b.len()).sum();
    let rock_bytes = files.iter().find(|(n, _)| n == "geologie_roches.png").map_or(0, |(_, b)| b.len());
    eprintln!("\n   EXPORT · geologie_roches.png {:.2} MB · the 12 favourability PNGs **{:.2} MB** · geologie.json {} bytes", rock_bytes as f64 / 1e6, fav_bytes as f64 / 1e6, files.last().map_or(0, |(_, b)| b.len()));
    // ── the maps and the views' captures (render.rs's colours over a grey shaded relief, north up)
    let mut base = vec![0u8; w * h * 4];
    let sl = rocks::slope_field(&z_m, w, h, CELL_KM * 1000.0);
    for c in 0..w * h {
        let g = if product.rocks[c] == rocks::NONE { 70u8 } else { (210.0 - 300.0 * sl[c].min(0.4)) as u8 };
        base[c * 4] = g;
        base[c * 4 + 1] = g;
        base[c * 4 + 2] = if product.rocks[c] == rocks::NONE { 95 } else { g };
        base[c * 4 + 3] = 255;
    }
    let mut rk = base.clone();
    render::rocks_rgba(&product.rocks, &mut rk);
    save_rgba(&out.join("view_rocks.png"), &rk, w, h);
    let mut all_c = base.clone();
    render::chances_rgba(&z, None, w, h, &mut all_c);
    save_rgba(&out.join("view_chances_all.png"), &all_c, w, h);
    for id in ["gold", "tin", "copper", "salt"] {
        let i = z.resources.iter().position(|r| r.id == id).expect("resource");
        let mut b = base.clone();
        render::chances_rgba(&z, Some(i), w, h, &mut b);
        save_rgba(&out.join(format!("view_chances_{id}.png")), &b, w, h);
    }
    let mut ctl = base.clone();
    render::chances_rgba(&zc, Some(gi), w, h, &mut ctl);
    save_rgba(&out.join("view_chances_gold_CONTROL_c3b_only.png"), &ctl, w, h);
    eprintln!("   views written to {}", out.display());
    eprintln!("\n==========  end Finding 152 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
