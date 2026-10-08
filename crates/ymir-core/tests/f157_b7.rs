//! ADR Finding 157-B7 — the production world's defects on the validation crops (OFF), stage by stage: comb teeth
//! (F124-P4), planar walls (F155), axis alignment (terrain R8; network R8 of the 1–3 km² reaches). Declared in
//! `docs/reports/lakes_gorges/f157_close/f157_declared.md`; split out of `f157_b` after it ran out of memory (it held
//! four worlds at once): one stage at a time here, the crops' origins those `f157_viz` / `f157_b` printed (the same
//! rule).
//!
//! Run: cargo test -p ymir-core --release --test f157_b7 -- --ignored --exact f157_b7 --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, build_world, viz_hd_lakes_on};
use std::time::Instant;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageResult, SegmentKind};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction, carve, skeleton};

const S_EQ: f32 = 0.024;
const DOMAIN_KM: f32 = 400.0;
const HALF: usize = 307;
/// The crops' data origins (south-first), as `f157_viz` printed them.
const CROPS: [(&str, usize, usize); 12] = [
    ("lake1", 3184, 2193),
    ("lake2", 3712, 2532),
    ("lake3", 3817, 2789),
    ("lake4", 1913, 3277),
    ("lake7", 1599, 3390),
    ("lake8", 3753, 3701),
    ("lake9", 4024, 3800),
    ("lake10", 2644, 3906),
    ("lake11", 1122, 3980),
    ("lake12", 3161, 4386),
    ("lake19", 3035, 5491),
    ("control", 2688, 5888),
];

fn chords(points: &[(u32, u32)], l: usize) -> Vec<f32> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i + l < points.len() {
        let (dx, dy) = (points[i + l].0 as f32 - points[i].0 as f32, points[i + l].1 as f32 - points[i].1 as f32);
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

fn r8c8(thetas: &[f32]) -> f32 {
    if thetas.is_empty() {
        return f32::NAN;
    }
    let (mut c8, mut s8) = (0f64, 0f64);
    for &t in thetas {
        c8 += (8.0 * t as f64).cos();
        s8 += (8.0 * t as f64).sin();
    }
    ((c8 * c8 + s8 * s8).sqrt() / thetas.len() as f64) as f32
}

#[test]
#[ignore]
fn f157_b7() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let bare = ValleyConstruction { light_k_time_fraction: None, ..off };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    eprintln!("\n==========  Finding 157 B7 . the production defects on the crops (OFF), stage by stage  ==========");
    // per crop: planar % and terrain R8 on a field
    let field_stats = |g: &GridF32| -> Vec<(f64, f32)> {
        let (w, h) = (g.width, g.height);
        let z: Vec<f32> = g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let cm = CELL_KM * 1000.0;
        CROPS
            .iter()
            .map(|&(_, x0, y0)| {
                let (mut land, mut pl) = (0usize, 0usize);
                for y in y0.max(1)..(y0 + 2 * HALF + 1).min(h - 1) {
                    for x in x0.max(1)..(x0 + 2 * HALF + 1).min(w - 1) {
                        let k = y * w + x;
                        if g.data[k] <= SEA {
                            continue;
                        }
                        land += 1;
                        let gx = (z[k + 1] - z[k - 1]) / (2.0 * cm);
                        let gy = (z[k + w] - z[k - w]) / (2.0 * cm);
                        if ((gx * gx + gy * gy).sqrt().atan().to_degrees() - 28.0).abs() <= 0.5 {
                            pl += 1;
                        }
                    }
                }
                let s = 2 * HALF + 1;
                let mut c = GridF32::new(s, s, 0.0);
                for y in 0..s {
                    for x in 0..s {
                        c.data[y * s + x] = g.data[(y0 + y) * w + x0 + x];
                    }
                }
                let lm: Vec<bool> = c.data.iter().map(|&v| v > SEA).collect();
                (100.0 * pl as f64 / land.max(1) as f64, aniso(&c, &lm, 16).r8)
            })
            .collect()
    };
    let net_stats = |dr: &C1DrainageResult, walls: &[bool], w: usize| -> Vec<(usize, usize, f32)> {
        CROPS
            .iter()
            .map(|&(_, x0, y0)| {
                let inside = |&(x, y): &(u32, u32)| (x as usize) >= x0 && (x as usize) <= x0 + 2 * HALF && (y as usize) >= y0 && (y as usize) <= y0 + 2 * HALF;
                let (mut teeth, mut reaches, mut th) = (0usize, 0usize, Vec::new());
                for (i, sg) in dr.rivers.segments.iter().enumerate() {
                    if dr.segment_kind[i] != SegmentKind::Watercourse || sg.points.len() < 2 || !sg.points.iter().any(inside) {
                        continue;
                    }
                    reaches += 1;
                    let on_wall = sg.points.iter().filter(|&&(x, y)| walls[y as usize * w + x as usize]).count();
                    if on_wall as f32 >= 0.8 * sg.points.len() as f32 {
                        teeth += 1;
                    }
                    let (lx, ly) = *sg.points.last().unwrap();
                    let a = dr.flow.accumulation.data[ly as usize * w + lx as usize] * cell_km2;
                    if (1.0..=3.0).contains(&a) {
                        let pts: Vec<(u32, u32)> = sg.points.iter().copied().filter(inside).collect();
                        th.extend(chords(&pts, 32));
                    }
                }
                (teeth, reaches, r8c8(&th))
            })
            .collect()
    };
    // S1 and the construction
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let w = s1.width;
    let st_s1 = field_stats(&s1);
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let (built, mk) = carve(&s1, &sk, &off, &ss);
    drop((sk, s1));
    let walls: Vec<bool> = mk.carved.iter().zip(&mk.floor).map(|(&c, &f)| c && !f).collect();
    drop(mk);
    let st_con = field_stats(&built);
    drop(built);
    eprintln!("   S1 and the construction measured ({:.0} s)", t0.elapsed().as_secs_f64());
    // C1 bare (no light pass)
    let (st_bare, nt_bare) = {
        let wd = build_world(kn(bare), None, PSEED, None);
        let st = field_stats(&wd.heightmap);
        let v = viz_hd_lakes_on(&wd, kn(bare), PSEED, 45.0, 40.0);
        drop(wd);
        (st, net_stats(&v.drainage, &walls, w))
    };
    eprintln!("   C1 bare measured ({:.0} s)", t0.elapsed().as_secs_f64());
    // OFF (production)
    let (st_off, nt_off) = {
        let wd = build_world(kn(off), None, PSEED, None);
        let st = field_stats(&wd.heightmap);
        let v = viz_hd_lakes_on(&wd, kn(off), PSEED, 45.0, 40.0);
        drop(wd);
        (st, net_stats(&v.drainage, &walls, w))
    };
    eprintln!("\n   per crop · planar walls (% of the land at 28° ± 0.5°) / terrain R8, per stage · comb teeth / watercourse reaches and the 1–3 km² network R8, C1 bare and OFF:");
    for (i, (name, x0, y0)) in CROPS.iter().enumerate() {
        eprintln!(
            "      {name:<8} ({x0}, {y0}) · S1 {:.2} % / {:.4} · construction {:.2} % / {:.4} · C1 bare {:.2} % / {:.4} · OFF {:.2} % / {:.4} · teeth C1 bare {}/{} (R8 {:.3}) · OFF {}/{} (R8 {:.3})",
            st_s1[i].0, st_s1[i].1, st_con[i].0, st_con[i].1, st_bare[i].0, st_bare[i].1, st_off[i].0, st_off[i].1, nt_bare[i].0, nt_bare[i].1, nt_bare[i].2, nt_off[i].0, nt_off[i].1, nt_off[i].2
        );
    }
    let mean = |v: &[(f64, f32)]| -> (f64, f64) { (v.iter().map(|x| x.0).sum::<f64>() / v.len() as f64, v.iter().map(|x| x.1 as f64).sum::<f64>() / v.len() as f64) };
    let (a, b, c, d) = (mean(&st_s1), mean(&st_con), mean(&st_bare), mean(&st_off));
    eprintln!(
        "   mean over the 12 crops · planar S1 {:.2} → construction {:.2} → C1 bare {:.2} → OFF {:.2} % · terrain R8 {:.4} → {:.4} → {:.4} → {:.4} · teeth C1 bare {} → OFF {} (of {} / {} reaches)",
        a.0, b.0, c.0, d.0, a.1, b.1, c.1, d.1,
        nt_bare.iter().map(|x| x.0).sum::<usize>(),
        nt_off.iter().map(|x| x.0).sum::<usize>(),
        nt_bare.iter().map(|x| x.1).sum::<usize>(),
        nt_off.iter().map(|x| x.1).sum::<usize>()
    );
    eprintln!("\n==========  end Finding 157 B7 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
