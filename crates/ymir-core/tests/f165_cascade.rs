//! ADR Finding 165 — the continent's physiography: the physics level follows C1's history. The references first
//! (`f165_europe`: Europe, the Alps and Corsica classified, the thresholds frozen before any world of ours is read),
//! then the diagnostic, the history and the measures (`f165_cascade`).
//! Declared in `docs/reports/relief_method/f165_continent/f165_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f165_cascade -- --ignored --exact f165_europe --nocapture
//!      cargo test -p ymir-core --release --test f165_cascade -- --ignored --exact f165_main::f165_cascade --nocapture
//!      cargo test -p ymir-core --release --test f165_cascade -- --ignored --exact f165_main::f165_a2_physics --nocapture

mod common;
#[path = "f165/main.rs"]
mod f165_main;

use common::cascade_bench::{NOTICE, load_grid, root};
use std::path::PathBuf;
use ymir_core::cascade::physio::{CLASS_NAMES, HILL_R7_M, MOUNTAIN, MacroStats, PLATEAU_Z_M, classify, macro_stats};
use ymir_core::grid::GridF32;

const CELL: f32 = 1.5625;

fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f165_continent/images")
}

/// A class map, north up (row 0 of the grid is the south).
fn class_png(classes: &[u8], w: usize, h: usize, scale: u32) -> image::RgbaImage {
    let col = |c: u8| match c {
        0 => [70u8, 110, 170],
        1 => [170, 210, 120],
        2 => [220, 190, 120],
        3 => [150, 160, 70],
        _ => [130, 80, 50],
    };
    let mut img = image::RgbaImage::new(w as u32, h as u32);
    for y in 0..h {
        for x in 0..w {
            let c = col(classes[y * w + x]);
            img.put_pixel(x as u32, (h - 1 - y) as u32, image::Rgba([c[0], c[1], c[2], 255]));
        }
    }
    image::imageops::resize(&img, w as u32 * scale, h as u32 * scale, image::imageops::FilterType::Nearest)
}

fn print_macro(tag: &str, m: &MacroStats) {
    eprintln!(
        "   {tag:<14} plain {:5.1} % · plateau {:5.1} % · hill {:5.1} % · mountain {:5.1} % · land {:.0} km² · mountain components {} · width {:.1} km · width / √land {:.4} · deciles {}",
        100.0 * m.fractions[0],
        100.0 * m.fractions[1],
        100.0 * m.fractions[2],
        100.0 * m.fractions[3],
        m.land_km2,
        m.mtn_components,
        m.mtn_width_km,
        m.mtn_ratio,
        m.deciles.iter().map(|v| format!("{v:.0}")).collect::<Vec<_>>().join("/")
    );
}

pub fn europe_grids() -> (GridF32, GridF32) {
    let e = load_grid(&root().join("data/europe/europe_1920.bin")).expect("data/europe (prep_europe.py)");
    let r = load_grid(&root().join("data/europe/europe_1920_regions.bin")).unwrap();
    (e, r)
}

#[test]
#[ignore]
fn f165_europe() {
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 165 . B — the references classified (thresholds: Kapos 2000 + hill R7 ≥ {HILL_R7_M} m, plateau z ≥ {PLATEAU_Z_M} m)  ==========");
    let (e, reg) = europe_grids();
    let ce = classify(&e.data, e.width, e.height, CELL);
    let me = macro_stats(&e.data, &ce, e.width, e.height, CELL);
    print_macro("Europe", &me);
    for (id, name, want) in [(1.0f32, "Paris Basin", "plain"), (2.0, "Massif Central", "plateau"), (3.0, "Meseta", "plateau"), (4.0, "Corsica (in Europe)", "mountain")] {
        let mut cnt = [0usize; 5];
        for k in 0..ce.len() {
            if reg.data[k] == id && ce[k] != 0 {
                cnt[ce[k] as usize] += 1;
            }
        }
        let tot: usize = cnt[1..].iter().sum();
        let best = (1..5).max_by_key(|&c| cnt[c]).unwrap();
        eprintln!(
            "   {name:<20} {} land cells: plain {:.0} % · plateau {:.0} % · hill {:.0} % · mountain {:.0} % → {} ({} wanted: {})",
            tot,
            100.0 * cnt[1] as f32 / tot.max(1) as f32,
            100.0 * cnt[2] as f32 / tot.max(1) as f32,
            100.0 * cnt[3] as f32 / tot.max(1) as f32,
            100.0 * cnt[4] as f32 / tot.max(1) as f32,
            CLASS_NAMES[best],
            want,
            if CLASS_NAMES[best] == want { "HOLDS" } else { "DOES NOT HOLD" }
        );
    }
    let a = load_grid(&root().join("data/europe/alps_256.bin")).unwrap();
    let ca = classify(&a.data, a.width, a.height, CELL);
    print_macro("Alps (400 km)", &macro_stats(&a.data, &ca, a.width, a.height, CELL));
    let c = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let cc = classify(&c.data, c.width, c.height, CELL);
    print_macro("Corsica", &macro_stats(&c.data, &cc, c.width, c.height, CELL));
    let mtn = |cl: &[u8]| cl.iter().filter(|&&x| x == MOUNTAIN).count() as f32 / cl.iter().filter(|&&x| x != 0).count().max(1) as f32;
    eprintln!("   Corsica mostly mountain: {:.0} % → {}", 100.0 * mtn(&cc), if mtn(&cc) > 0.5 { "HOLDS" } else { "DOES NOT HOLD" });
    class_png(&ce, e.width, e.height, 1).save(out().join("f165_classes_europe.png")).unwrap();
    class_png(&ca, a.width, a.height, 2).save(out().join("f165_classes_alps.png")).unwrap();
    class_png(&cc, c.width, c.height, 4).save(out().join("f165_classes_corse.png")).unwrap();
    eprintln!(
        "   images f165_classes_{{europe,alps,corse}}.png: sea blue, plain light green, plateau tan, hill olive, mountain brown, north up · ETOPO 2022 (NOAA NCEI, public domain) · {NOTICE}"
    );
}
