//! ADR 0001 Finding 101 — name the three basins, attribute the two canyons, show the tile.
//! **No production change.**
//!
//! ⚠️ Two framing points established from the grep before any measurement:
//!
//!  * **The basin count is NOT monotone in the factor.** Finding 100 measured 11 / 11 / 14 / 14 /
//!    13 / 13 for delivered / ×0.3 / ×0.4 / ×0.5 / ×0.6 / ×0.7. A quantity that moves like that as
//!    its driver moves steadily is not responding to the driver, so block A must survive the
//!    possibility that ×0.4 created **no** basin and that marginal bodies crossed a detector.
//!  * **The floors will be set by the bathymetry clamp, not by the terrain.** The earliest
//!    `wc == 2` hit in the dossier (L1570, Finding 38) reports *"floor −19.9 m"* for every
//!    below-sea region — the **old 20 m shelf clamp**, which Finding 79 took to 1 m. Since
//!    `apply_bathymetry_profile` runs AFTER the incision and takes `min_depth =
//!    shelf_min_depth_m.max(1.0)`, today's floors should sit at ≈ −1.0 m **for old and new alike**,
//!    which makes the final-field depth useless as a discriminator. This bench therefore reads the
//!    floors with **`bathymetry_off`** as well, which is the only stage where the terrain's own
//!    depth survives.
//!  * And Finding 90-D supplies a THIRD possible origin the round does not list: **the breach
//!    itself makes 6 100 of 387 710 `wc == 2` cells**, in 169 trenches (p50 39 cells, deepest
//!    point p50 1.85 m below sea).
//!
//! Run: cargo test -p ymir-core --release --test f101_promote -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, aniso, build_field, f95_criteria, pct, slope_deg, sorted, tile};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const ORACLE_P50_M: f32 = 488.1;
const DELIVERED_P50_M: f32 = 424.2;
const KS: f32 = 0.0451;
const TX: usize = 2048;
const TY: usize = 5120;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f101_promote");

fn hillshade(t: &GridF32, n2m: f32) -> GridF32 {
    let (w, h) = (t.width, t.height);
    let (az, alt) = (315f32.to_radians(), 45f32.to_radians());
    let (lx, ly, lz) = (alt.cos() * az.sin(), -alt.cos() * az.cos(), alt.sin());
    let mut g = GridF32::new(w, h, 0.0);
    for y in 0..h {
        for x in 0..w {
            let (gx, gy) = t.gradient_at(x, y);
            let (dx, dy) = (gx * n2m / CELL_M, gy * n2m / CELL_M);
            let nrm = (dx * dx + dy * dy + 1.0).sqrt();
            g.data[y * w + x] = ((-dx * lx - dy * ly + lz) / nrm).clamp(0.0, 1.0);
        }
    }
    g
}

fn sigma_p50(f: &GridF32, land: &[bool], n2m: f32) -> f32 {
    let (w, h) = (f.width, f.height);
    let mut out = Vec::new();
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if !land[k] {
                continue;
            }
            let (mut s, mut s2) = (0f64, 0f64);
            for dy in 0..3 {
                for dx in 0..3 {
                    let v = (f.data[(y + dy - 1) * w + x + dx - 1] * n2m) as f64;
                    s += v;
                    s2 += v * v;
                }
            }
            out.push(((s2 / 9.0 - (s / 9.0) * (s / 9.0)).max(0.0)).sqrt() as f32);
        }
    }
    pct(&sorted(out), 0.50)
}

/// A `wc == 2` component, identified by its BOTTOM CELL (Findings 81 / 85: never by id — ids are
/// scan-order and do not survive a rebuild).
struct Comp {
    bottom: usize,
    cells: usize,
    floor_m: f32,
    /// distance in cells from the component to the nearest sea cell, 0 = touching
    to_sea: usize,
}

fn components(f: &GridF32, wc: &[u8], w: usize, h: usize, ss: &SteinSteinParams) -> Vec<Comp> {
    let n = w * h;
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for s in 0..n {
        if wc[s] != 2 || seen[s] {
            continue;
        }
        let (mut stack, mut cells) = (vec![s], Vec::new());
        seen[s] = true;
        while let Some(k) = stack.pop() {
            cells.push(k);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for d in 0..8 {
                let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if wc[nk] == 2 && !seen[nk] {
                    seen[nk] = true;
                    stack.push(nk);
                }
            }
        }
        let bottom =
            *cells.iter().min_by(|&&a, &&b| f.data[a].partial_cmp(&f.data[b]).unwrap()).unwrap();
        // distance to open sea: BFS-free proxy — the minimum Chebyshev distance to a wc==3/ocean
        // cell within a 64-cell box, capped. Declared as a proxy.
        let (bx, by) = ((bottom % w) as i32, (bottom / w) as i32);
        let mut to_sea = 65usize;
        'r: for r in 1..64i32 {
            for dy in -r..=r {
                for dx in -r..=r {
                    if dx.abs() != r && dy.abs() != r {
                        continue;
                    }
                    let (nx, ny) = (bx + dx, by + dy);
                    if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if wc[nk] != 2 && f.data[nk] <= SEA {
                        to_sea = r as usize;
                        break 'r;
                    }
                }
            }
        }
        out.push(Comp {
            bottom,
            cells: cells.len(),
            floor_m: c1_altitude_norm_to_metres(f.data[bottom], ss),
            to_sea,
        });
    }
    out
}

/// ⚠️ **ADR Finding 107 RENAMED THE CLASS.** What this function computes is the **CUT** column,
/// `pre-incision − delivered`, which is what Findings 87-104 printed as "the canyon class". The
/// gate is now the **over-dug depression**: `filled − raw > 50 m`, the floor below its OWN sill.
/// A valley the drainage runs through has a large cut and a fill of zero, so the two columns are
/// not the same class. This function is kept unrenamed in spirit -- it reproduces a historical
/// number -- but it is `cut_class`, not the gate.
///
/// The cut class of Finding 88-D4, computed on an ARBITRARY stage with a FIXED lake inventory,
/// so "eroded vs breached" is an A/B on one population. ⚠️ This inventory comes straight from
/// `c1_drainage_windowed` on the breached field and is SIMPLER than `f95_criteria`'s (no water
/// balance, no below-sea merge), so its absolute count is not the table's — only the A/B is.
fn cut_class(
    stage: &GridF32,
    pre: &GridF32,
    lake_map: &[u32],
    ids: &[u32],
    cell_km2: f32,
    n2m: f32,
    w: usize,
    h: usize,
) -> (usize, usize, Vec<(u32, f32, f32, usize)>) {
    let n = w * h;
    let m =
        |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &SteinSteinParams::default());
    let (mut klass, mut scanned, mut named) = (0usize, 0usize, Vec::new());
    // ONE pass to bucket the cells by lake: the first version scanned all 67 M cells per lake.
    let want: HashSet<u32> = ids.iter().copied().collect();
    let mut by_lake: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..n {
        let id = lake_map[k];
        if id != 0 && want.contains(&id) {
            by_lake.entry(id).or_default().push(k);
        }
    }
    let mut keys: Vec<u32> = by_lake.keys().copied().collect();
    keys.sort_unstable();
    for id in keys {
        let cells = by_lake[&id].clone();
        if cells.len() as f32 * cell_km2 < 1.0 {
            continue;
        }
        scanned += 1;
        let cuts: Vec<f32> = cells.iter().map(|&k| m(pre, k) - m(stage, k)).collect();
        let inside: HashSet<usize> = cells.iter().copied().collect();
        let mut ring: HashSet<usize> = HashSet::new();
        for &k in &cells {
            for d in 0..8 {
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                let nk = ny * w + nx;
                if !inside.contains(&nk) {
                    ring.insert(nk);
                }
            }
        }
        let rim: Vec<f32> = ring.iter().map(|&k| slope_deg(stage, k, n2m, w, h)).collect();
        let empty = rim.is_empty();
        let (cp, rp) = (pct(&sorted(cuts), 0.50), pct(&sorted(rim), 0.50));
        if cp > 50.0 && !empty && rp > 30.0 {
            klass += 1;
            named.push((id, cp, rp, cells.len()));
        }
    }
    (klass, scanned, named)
}

fn save(g: &GridF32, name: &str) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    g.save_png_u8(&Path::new(OUT).join(name)).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!("   wrote {name}");
}

#[test]
#[ignore]
fn f101_promote() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 101 . naming the basins, attributing the canyons  ==========");

    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let delivered = build_field(Knobs::passes(2));
    let kn =
        |m: f32| Knobs { slope_floor_uk: Some(KS * m), depression_floor: true, ..Knobs::passes(2) };
    let f04 = build_field(kn(0.4));
    let f05 = build_field(kn(0.5));
    // the only stage where the terrain's own depth survives the clamp
    let f04_nb = build_field(Knobs { bathymetry_off: true, ..kn(0.4) });
    let del_nb = build_field(Knobs { bathymetry_off: true, ..Knobs::passes(2) });
    eprintln!("   five fields in {:.1} s", t0.elapsed().as_secs_f64());

    // instrument line: the ×0.4 control must reproduce the sweep
    let d04 = c1_drainage_windowed(&f04, None, &on, &ss, DOMAIN_KM);
    let bf04 = breach_monotone(&f04, &d04.flow.filled, &d04.lake_map, SEA, w, h);
    let c04 = f95_criteria(&f04, &bf04, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
    let land04: Vec<bool> = (0..n).map(|k| f04.data[k] > SEA).collect();
    eprintln!(
        "\n-- INSTRUMENT . the x0.4 control against Finding 100 --\n   relief **{:.1} m** (sweep \
         493.4) . canyons **{}/{}** (sweep 2/29) . basins **{}** (sweep 14) . catchment **{:.0}** \
         (sweep 9780) . R8 **{:.4}** (sweep 0.0563)",
        c04.p50,
        c04.klass,
        c04.scanned,
        c04.below_sea_spillways,
        c04.max_catchment_km2,
        aniso(&f04, &land04, 16).r8
    );

    // ══ A — the three basins, by BOTTOM CELL ═════════════════════════════════
    let wc_del = water_class(&delivered, SEA);
    let wc_04 = water_class(&f04, SEA);
    let cd = components(&delivered, &wc_del, w, h, &ss);
    let c4 = components(&f04, &wc_04, w, h, &ss);
    eprintln!(
        "\n-- A . wc==2 components by BOTTOM CELL (never by id: Findings 81, 85) --\n   delivered \
         **{}** components . x0.4 **{}**",
        cd.len(),
        c4.len()
    );
    // match by bottom cell, then by proximity (a bottom can shift by a cell)
    let dm: HashMap<usize, &Comp> = cd.iter().map(|c| (c.bottom, c)).collect();
    let near = |k: usize| -> Option<&Comp> {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for r in 0..=8i32 {
            for dy in -r..=r {
                for dx in -r..=r {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                        continue;
                    }
                    if let Some(c) = dm.get(&(ny as usize * w + nx as usize)) {
                        return Some(c);
                    }
                }
            }
        }
        None
    };
    let mut matched = 0usize;
    eprintln!("   the components of x0.4 that have NO counterpart in the delivered field:");
    for c in &c4 {
        match near(c.bottom) {
            Some(d) => {
                matched += 1;
                let da = (c.cells as f32 / d.cells.max(1) as f32 - 1.0) * 100.0;
                if da.abs() > 5.0 {
                    eprintln!(
                        "     (matched, but area moved {da:+.0} %: {} -> {} cells at ({}, {}))",
                        d.cells,
                        c.cells,
                        c.bottom % w,
                        c.bottom / w
                    );
                }
            }
            None => {
                let (bx, by) = (c.bottom % w, c.bottom / w);
                // the same cell's altitude in the delivered field and in the NO-BATHYMETRY fields
                eprintln!(
                    "     ⛔ NEW at ({bx}, {by}): **{} cells = {:.3} km²** . floor **{:.3} m** . \
                     distance to sea **{}** cells . at that cell: delivered **{:.3} m**, x0.4 \
                     no-bathymetry **{:.3} m**, delivered no-bathymetry **{:.3} m**, pre-incision \
                     **{:.1} m**",
                    c.cells,
                    c.cells as f32 * cell_km2,
                    c.floor_m,
                    c.to_sea,
                    c1_altitude_norm_to_metres(delivered.data[c.bottom], &ss),
                    c1_altitude_norm_to_metres(f04_nb.data[c.bottom], &ss),
                    c1_altitude_norm_to_metres(del_nb.data[c.bottom], &ss),
                    c1_altitude_norm_to_metres(pre.data[c.bottom], &ss)
                );
            }
        }
    }
    eprintln!(
        "   control: **{matched} of {}** delivered components found again at x0.4 by bottom cell",
        cd.len()
    );
    let fl: Vec<f32> = sorted(c4.iter().map(|c| c.floor_m).collect());
    let fld: Vec<f32> = sorted(cd.iter().map(|c| c.floor_m).collect());
    eprintln!(
        "   ⚠️ FLOORS, the clamp test: x0.4 p10/p50/p90 **{:.3} / {:.3} / {:.3} m** . delivered \
         **{:.3} / {:.3} / {:.3} m** -- Finding 38 reported -19.9 m under the OLD 20 m clamp; \
         production is 1 m since Finding 79",
        pct(&fl, 0.10),
        pct(&fl, 0.50),
        pct(&fl, 0.90),
        pct(&fld, 0.10),
        pct(&fld, 0.50),
        pct(&fld, 0.90)
    );

    // ══ B — the two canyons: eroded or breached ══════════════════════════════
    let ids: Vec<u32> = {
        let mut s: HashSet<u32> = d04.lake_map.iter().copied().filter(|&i| i != 0).collect();
        let mut v: Vec<u32> = s.drain().collect();
        v.sort_unstable();
        v
    };
    let (ke, se, ne) = cut_class(&f04, &pre, &d04.lake_map, &ids, cell_km2, n2m, w, h);
    let (kb, sb, nb) = cut_class(&bf04, &pre, &d04.lake_map, &ids, cell_km2, n2m, w, h);
    eprintln!(
        "\n-- B . the canyon class at x0.4, one FIXED inventory, two stages --\n   ⚠️ inventory \
         from `c1_drainage_windowed` on the breached field, SIMPLER than `f95_criteria`'s (no \
         water balance, no below-sea merge): {} bodies >= 1 km² scanned\n   **ERODED stage: \
         {ke}** . **BREACHED stage: {kb}** (scanned {se} / {sb})",
        se.max(sb)
    );
    for (tag, v) in [("eroded", &ne), ("breached", &nb)] {
        for (id, cp, rp, cells) in v {
            eprintln!(
                "     {tag}: body {id} . {cells} cells ({:.2} km²) . median cut **{cp:.1} m** . \
                 rim p50 **{rp:.1}°**",
                *cells as f32 * cell_km2
            );
        }
    }

    // ══ C — the images ═══════════════════════════════════════════════════════
    eprintln!("\n-- C . the tile ({TX}, {TY}) and the continent --");
    let ld: Vec<bool> = (0..n).map(|k| delivered.data[k] > SEA).collect();
    for (label, f) in [("delivered", &delivered), ("x0.4", &f04), ("x0.5", &f05)] {
        let (g, gl) = tile(f, &ld, TX, TY, 1024);
        let cut: Vec<f32> = (0..n)
            .filter(|&k| pre.data[k] > SEA && f.data[k] > SEA)
            .map(|k| (pre.data[k] - f.data[k]) * n2m)
            .collect();
        eprintln!(
            "   {label:<10} R8 **{:.4}** . sigma **{:.2} m** . median cut per cell **{:.1} m**",
            aniso(&g, &gl, 16).r8,
            sigma_p50(&g, &gl, n2m),
            pct(&sorted(cut), 0.50)
        );
        save(&hillshade(&g, n2m), &format!("{label}_hillshade.png"));
    }
    // the continent-scale water view: the scale at which "France, not Scotland" was asked
    for (label, f) in [("delivered", &delivered), ("x0.4", &f04)] {
        let wc = water_class(f, SEA);
        let d = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
        let mut g = GridF32::new(1024, 1024, 0.0);
        for y in 0..1024 {
            for x in 0..1024 {
                let (mut water, mut land) = (0u32, 0u32);
                for dy in 0..8 {
                    for dx in 0..8 {
                        let k = (y * 8 + dy) * w + x * 8 + dx;
                        if f.data[k] <= SEA {
                            continue;
                        }
                        land += 1;
                        if d.lake_map[k] != 0 || wc[k] == 2 {
                            water += 1;
                        }
                    }
                }
                g.data[y * 1024 + x] =
                    if land == 0 { 0.0 } else { 0.25 + 0.75 * (1.0 - water as f32 / land as f32) };
            }
        }
        save(&g, &format!("continent_water_{label}.png"));
    }
    eprintln!("\n==========  end Finding 101 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
