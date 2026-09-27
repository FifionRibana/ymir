//! ADR 0001 Finding 122 — **the "brush" lakes of the hybrid: attributed before corrected.**
//!
//! The author judged the four Finding 121 states and named one defect: lakes "painted with a
//! brush" (the lower-right one). Block A attributes it: per lake of the C2/10 world, the share of
//! its cells that are a CONSTRUCTED valley floor, and its arms.
//!
//! ## Instrument (rule 12)
//!
//! * **Constructed-floor mask** — `carve` of the PRE build (`Knobs::no_incision`) at the shipped
//!   `F121_AGE_K`, as Finding 121 did for its interfluve and wall masks: the in-pipeline carve is
//!   not returned by the pipeline, and the PRE build also carries the C-2 rims and the bathymetry
//!   (crater and sub-sea cells only). Declared.
//! * **Lakes** — the assembled inventory of the C2/10 world (the export chain, lake_map of every id,
//!   family 2 included): the objects the author's microscope lists.
//! * **Arms** — a morphological opening of each lake by a disk of radius `W` (the law's floor width
//!   at the lake's own trunk cells, else 200 m): the arms are the lake cells the opening removes,
//!   counted as 8-connected components of ≥ 9 cells (declared).
//! * **Frame** — the viz frames seed 1 at torus offset x = 0.09375 (the author's export manifest),
//!   the benches at 0.0: a viz cell `x_v` is the bench cell `x_v + 768`. Both are printed.
//! * **Control** — the delivered (OFF) world's lakes over the SAME mask.
//!
//! Run: cargo test -p ymir-core --release --test f122_brush -- --ignored f122_a --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::collections::VecDeque;
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::valley_construction::{
    CarveMasks, F121_AGE_K, Skeleton, ValleyConstruction, carve, skeleton,
};
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
/// The viz's torus x-offset for seed 1 (author's manifest), in 8192² cells: 0.09375 × 8192.
const VIZ_SHIFT_X: usize = 768;
const MIN_ARM_CELLS: usize = 9;

fn dcfg() -> C1DrainageConfig {
    let mut d = C1DrainageConfig::default();
    d.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    d.thresholds.full_tree = false;
    d
}

/// The export chain's assembled inventory of an eroded field.
fn inventory(g: &GridF32, ss: &SteinSteinParams) -> C1DrainageResult {
    let (w, h) = (g.width, g.height);
    let d = c1_drainage_windowed(g, None, &dcfg(), ss, DOMAIN_KM);
    let bre = breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, w, h);
    let cl = c1_climate_placed(&bre, ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
    assemble_hd_drainage(&bre, &dc, Some(d), &dcfg(), ss, DOMAIN_KM, GEO_RATIO, None, false)
        .drainage
}

/// Chessboard BFS distance, in cells, from `src`, restricted to `inside` (sources need not be).
fn bfs(src: &[usize], inside: &dyn Fn(usize) -> bool, w: usize, h: usize, cap: u32) -> Vec<u32> {
    let mut d = vec![u32::MAX; w * h];
    let mut q = VecDeque::new();
    for &s in src {
        d[s] = 0;
        q.push_back(s);
    }
    while let Some(k) = q.pop_front() {
        if d[k] >= cap {
            continue;
        }
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let nk = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                if d[nk] == u32::MAX && inside(nk) {
                    d[nk] = d[k] + 1;
                    q.push_back(nk);
                }
            }
        }
    }
    d
}

/// One lake read against the construction's masks.
struct LakeRead {
    id: u32,
    km2: f32,
    level_m: f32,
    cx_bench: f32,
    cy: f32,
    floor_cell: usize,
    floor_share: f32,
    carved_share: f32,
    w_m: f32,
    arms: usize,
    arm_km2: f32,
    arm_floor_share: f32,
    arm_trunk_share: f32,
    arm_drowned_floor_share: f32,
    dl: f32,
}

#[allow(clippy::too_many_arguments)]
fn read_lakes(
    dr: &C1DrainageResult,
    g: &GridF32,
    sk: &Skeleton,
    floor: &[bool],
    carved: &[bool],
    trunk_near: &[bool],
    vc: &ValleyConstruction,
    ss: &SteinSteinParams,
) -> Vec<LakeRead> {
    let (w, h) = (g.width, g.height);
    let n = w * h;
    let cell_km2 = CELL_KM * CELL_KM;
    let mut cells_of: std::collections::BTreeMap<u32, Vec<usize>> = Default::default();
    for (k, &id) in dr.lake_map.iter().enumerate().take(n) {
        if id != 0 {
            cells_of.entry(id).or_default().push(k);
        }
    }
    let mut out = Vec::new();
    for l in &dr.lakes {
        let Some(cells) = cells_of.get(&l.base.id) else { continue };
        let id = l.base.id;
        let inl = |k: usize| dr.lake_map[k] == id;
        let nc = cells.len() as f32;
        let floor_share = cells.iter().filter(|&&k| floor[k]).count() as f32 / nc;
        let carved_share = cells.iter().filter(|&&k| carved[k]).count() as f32 / nc;
        // the law's floor width at the lake's own trunk cells
        let ws = sorted(
            cells.iter().filter(|&&k| sk.trunk[k]).map(|&k| vc.width_m(sk.area_km2[k])).collect(),
        );
        let w_m = if ws.is_empty() { 200.0 } else { pct(&ws, 0.5) };
        let r = (w_m / CELL_M).round().max(1.0) as u32;
        // opening by a disk of radius r: shore distance → core → body
        let shore: Vec<usize> = cells
            .iter()
            .copied()
            .filter(|&k| {
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                (-1i32..=1).any(|dy| {
                    (-1i32..=1).any(|dx| {
                        !inl((y + dy).rem_euclid(h as i32) as usize * w
                            + (x + dx).rem_euclid(w as i32) as usize)
                    })
                })
            })
            .collect();
        let ds = bfs(&shore, &inl, w, h, r + 1);
        let core: Vec<usize> = cells.iter().copied().filter(|&k| ds[k] > r).collect();
        let body = bfs(&core, &inl, w, h, r);
        let arm: Vec<usize> = cells.iter().copied().filter(|&k| body[k] > r).collect();
        // components of the arm cells
        let mut seen = std::collections::HashSet::new();
        let (mut arms, mut arm_cells) = (0usize, Vec::new());
        let armset: std::collections::HashSet<usize> = arm.iter().copied().collect();
        for &s in &arm {
            if seen.contains(&s) {
                continue;
            }
            let mut comp = vec![s];
            let mut q = vec![s];
            seen.insert(s);
            while let Some(k) = q.pop() {
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                for dy in -1i32..=1 {
                    for dx in -1i32..=1 {
                        let nk = (y + dy).rem_euclid(h as i32) as usize * w
                            + (x + dx).rem_euclid(w as i32) as usize;
                        if armset.contains(&nk) && seen.insert(nk) {
                            q.push(nk);
                            comp.push(nk);
                        }
                    }
                }
            }
            if comp.len() >= MIN_ARM_CELLS {
                arms += 1;
                arm_cells.extend(comp);
            }
        }
        let na = arm_cells.len().max(1) as f32;
        let arm_floor = arm_cells.iter().filter(|&&k| floor[k]).count();
        let arm_trunk = arm_cells.iter().filter(|&&k| trunk_near[k]).count() as f32 / na;
        let drowned = arm_cells
            .iter()
            .filter(|&&k| floor[k] && !sk.chi_m[k].is_nan() && sk.floor_m(k, vc.age_k) < l.level_m)
            .count() as f32
            / arm_floor.max(1) as f32;
        // shoreline development D_L = perimeter / (2 sqrt(π A)), 4-connected edges
        let mut per = 0usize;
        for &k in cells {
            let (x, y) = (k % w, k / w);
            for (dx, dy) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                let nx = (x as i32 + dx).rem_euclid(w as i32) as usize;
                let ny = (y as i32 + dy).rem_euclid(h as i32) as usize;
                if !inl(ny * w + nx) {
                    per += 1;
                }
            }
        }
        let a = nc * cell_km2;
        let dl = per as f32 * CELL_KM / (2.0 * (std::f32::consts::PI * a).sqrt());
        let (sx, sy) =
            cells.iter().fold((0f64, 0f64), |(a, b), &k| (a + (k % w) as f64, b + (k / w) as f64));
        let floor_cell =
            *cells.iter().min_by(|&&a, &&b| g.data[a].total_cmp(&g.data[b])).expect("non-empty");
        out.push(LakeRead {
            id,
            km2: a,
            level_m: l.level_m,
            cx_bench: (sx / nc as f64) as f32,
            cy: (sy / nc as f64) as f32,
            floor_cell,
            floor_share,
            carved_share,
            w_m,
            arms,
            arm_km2: arm_cells.len() as f32 * cell_km2,
            arm_floor_share: arm_floor as f32 / na,
            arm_trunk_share: arm_trunk,
            arm_drowned_floor_share: drowned,
            dl,
        });
        let _ = ss;
    }
    out.sort_by(|a, b| b.km2.total_cmp(&a.km2));
    out
}

fn print_lakes(name: &str, v: &[LakeRead], w: usize) {
    eprintln!(
        "\n   ── {name} — every inventoried lake ({}) ──\n   {:>8} {:>8} {:>9} {:>14} {:>14} {:>11} {:>7} \
         {:>7} {:>6} {:>5} {:>8} {:>8} {:>8} {:>9} {:>6}",
        v.len(),
        "id",
        "km²",
        "level m",
        "centre bench",
        "centre viz",
        "floor cell",
        "floor%",
        "carved%",
        "W m",
        "arms",
        "arm km²",
        "arm fl%",
        "arm tr%",
        "drowned%",
        "D_L"
    );
    for l in v {
        let vx = (l.cx_bench as i64 - VIZ_SHIFT_X as i64).rem_euclid(w as i64);
        eprintln!(
            "   {:>8} {:>8.2} {:>9.1} {:>6.0},{:<7.0} {:>6},{:<7.0} {:>11} {:>6.1}% {:>6.1}% {:>6.0} {:>5} \
             {:>8.2} {:>7.1}% {:>7.1}% {:>8.1}% {:>6.2}",
            l.id,
            l.km2,
            l.level_m,
            l.cx_bench,
            l.cy,
            vx,
            l.cy,
            l.floor_cell,
            100.0 * l.floor_share,
            100.0 * l.carved_share,
            l.w_m,
            l.arms,
            l.arm_km2,
            100.0 * l.arm_floor_share,
            100.0 * l.arm_trunk_share,
            100.0 * l.arm_drowned_floor_share,
            l.dl
        );
    }
}

#[test]
#[ignore]
fn f122_a() {
    let ss = SteinSteinParams::default();
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 122-A . the brush, attributed  ==========");
    let vc = ValleyConstruction::f121(F121_AGE_K, Some(0.1));
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
    let (_, masks) = carve(&pre, &sk, &vc, &ss);
    let trunk_src: Vec<usize> = (0..n).filter(|&k| sk.trunk[k]).collect();
    let near = bfs(&trunk_src, &|_| true, w, h, 2);
    let trunk_near: Vec<bool> = near.iter().map(|&d| d <= 2).collect();
    drop(near);
    eprintln!(
        "   masks from the PRE carve at k = {F121_AGE_K}: floor {} cells · carved {} cells · trunk {} \
         cells ({:.1} s)",
        masks.floor.iter().filter(|&&b| b).count(),
        masks.carved.iter().filter(|&&b| b).count(),
        trunk_src.len(),
        t0.elapsed().as_secs_f64()
    );
    drop(pre);

    let c2 = build_field_seed(
        Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
        PSEED,
    );
    let dr = inventory(&c2, &ss);
    let v = read_lakes(&dr, &c2, &sk, &masks.floor, &masks.carved, &trunk_near, &vc, &ss);
    print_lakes("C2/10 hybrid", &v, w);
    // the lower-right candidates in the VIZ frame (x right, y down)
    let mut lr: Vec<&LakeRead> = v.iter().filter(|l| l.km2 >= 1.0).collect();
    lr.sort_by(|a, b| {
        let key = |l: &LakeRead| {
            (l.cx_bench as i64 - VIZ_SHIFT_X as i64).rem_euclid(w as i64) as f32 + l.cy
        };
        key(b).total_cmp(&key(a))
    });
    eprintln!(
        "\n   lower-right candidates (viz frame, max x + y, ≥ 1 km²): {}",
        lr.iter()
            .take(3)
            .map(|l| format!("id {} ({:.1} km², floor cell {})", l.id, l.km2, l.floor_cell))
            .collect::<Vec<_>>()
            .join(" · ")
    );
    drop((c2, dr));

    let off = build_field_seed(Knobs::passes(2), PSEED);
    let dro = inventory(&off, &ss);
    let vo = read_lakes(&dro, &off, &sk, &masks.floor, &masks.carved, &trunk_near, &vc, &ss);
    print_lakes("CONTROL · OFF delivered (same masks)", &vo, w);
    let _ = h;
    eprintln!("\n==========  end 122-A . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f122_brush");

/// Hillshade + sea + lakes + rivers of a crop, NORTH UP (the viz mirrors the south-first rows,
/// Finding 27), optionally tinting the construction's carved (walls) and floor cells.
fn render(
    name: &str,
    g: &GridF32,
    dr: &C1DrainageResult,
    (x0, y0, cw, ch): (usize, usize, usize, usize),
    tint: Option<(&[bool], &[bool])>,
    ss: &SteinSteinParams,
) {
    use ymir_core::tectonics_c1::drainage::SegmentKind;
    use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    let (w, h) = (g.width, g.height);
    let m = |k: usize| c1_altitude_norm_to_metres(g.data[k], ss);
    let mut river = std::collections::HashSet::new();
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] == SegmentKind::Watercourse {
            for &(x, y) in &s.points {
                river.insert(y as usize * w + x as usize);
            }
        }
    }
    let (az, alt) = (315f32.to_radians(), 45f32.to_radians());
    let mut img = image::RgbImage::new(cw as u32, ch as u32);
    for yy in 0..ch {
        for xx in 0..cw {
            let (x, y) = ((x0 + xx) % w, (y0 + yy) % h);
            let k = y * w + x;
            let at = |dx: i32, dy: i32| {
                m(((y as i32 + dy).rem_euclid(h as i32) as usize) * w
                    + (x as i32 + dx).rem_euclid(w as i32) as usize)
            };
            let gx = (at(1, 0) - at(-1, 0)) / (2.0 * CELL_M);
            let gy = (at(0, 1) - at(0, -1)) / (2.0 * CELL_M);
            // normal · light, light from the NW (azimuth clockwise from north), data y = NORTH
            let (lx, ly, lz) = (alt.cos() * az.sin(), alt.cos() * az.cos(), alt.sin());
            let shade =
                ((-gx * lx - gy * ly + lz) / (1.0 + gx * gx + gy * gy).sqrt()).clamp(0.0, 1.0);
            let mut c = if g.data[k] <= SEA && dr.lake_map[k] == 0 {
                [150.0f32, 185.0, 220.0]
            } else {
                let v = 55.0 + 200.0 * shade;
                [v, v, v]
            };
            if let Some((carved, floor)) = tint {
                if floor[k] {
                    c = [c[0] * 0.55 + 230.0 * 0.45, c[1] * 0.55 + 120.0 * 0.45, c[2] * 0.55];
                } else if carved[k] {
                    c = [c[0] * 0.75 + 240.0 * 0.25, c[1] * 0.75 + 200.0 * 0.25, c[2] * 0.75];
                }
            }
            if dr.lake_map[k] != 0 {
                c = [
                    c[0] * 0.35 + 40.0 * 0.65,
                    c[1] * 0.35 + 110.0 * 0.65,
                    c[2] * 0.35 + 210.0 * 0.65,
                ];
            }
            if river.contains(&k) {
                c = [20.0, 45.0, 150.0];
            }
            img.put_pixel(
                xx as u32,
                (ch - 1 - yy) as u32,
                image::Rgb([c[0] as u8, c[1] as u8, c[2] as u8]),
            );
        }
    }
    std::fs::create_dir_all(OUT).expect("report dir");
    let p = std::path::Path::new(OUT).join(name);
    img.save(&p).expect("png");
    eprintln!("   wrote {}", p.display());
}

/// ADR Finding 122-B/C — **the remedy measured on the lakes, and the renders for the author.**
///
/// B: the C2/10 world with `basin_base`, read with block A's instrument, against C2/10.
/// C: the tile (2048, 5120) (bench frame; the viz frame is x − 768) in relief for the livré,
/// A1+B2, C1 bare, C2/10 and C2/10 + basin base, and the brush lake 1000011's crop. Images are
/// NORTH UP, as the viz shows them.
///
/// Run: cargo test -p ymir-core --release --test f122_brush -- --ignored f122_bc --nocapture
#[test]
#[ignore]
fn f122_bc() {
    let ss = SteinSteinParams::default();
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 122-B/C . the remedy's lakes, and the tile  ==========");
    let v121 = ValleyConstruction::f121(F121_AGE_K, Some(0.1));
    let v122 = ValleyConstruction::f122(F121_AGE_K, Some(0.1));
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let sk121 = skeleton(&pre, &v121, &ss, DOMAIN_KM);
    let (_, m121) = carve(&pre, &sk121, &v121, &ss);
    let sk122 = skeleton(&pre, &v122, &ss, DOMAIN_KM);
    let (_, m122) = carve(&pre, &sk122, &v122, &ss);
    let trunk_src: Vec<usize> = (0..n).filter(|&k| sk122.trunk[k]).collect();
    let near = bfs(&trunk_src, &|_| true, w, h, 2);
    let trunk_near: Vec<bool> = near.iter().map(|&d| d <= 2).collect();
    drop(near);
    let moved: Vec<f32> = (0..n)
        .filter(|&k| sk122.trunk[k])
        .map(|k| sk122.floor_m(k, F121_AGE_K) - sk121.floor_m(k, F121_AGE_K))
        .filter(|d| d.abs() > 0.01)
        .collect();
    let mv = sorted(moved.clone());
    let q = |p: f64| if mv.is_empty() { 0.0 } else { pct(&mv, p) };
    eprintln!(
        "   basin_base moves the floor on **{} of {} trunk cells** · Δfloor p10/p50/p90 **{:.1} / \
         {:.1} / {:.1} m** · carved {} → {} cells · floor {} → {} cells",
        moved.len(),
        trunk_src.len(),
        q(0.10),
        q(0.50),
        q(0.90),
        m121.carved.iter().filter(|&&b| b).count(),
        m122.carved.iter().filter(|&&b| b).count(),
        m121.floor.iter().filter(|&&b| b).count(),
        m122.floor.iter().filter(|&&b| b).count()
    );
    drop(pre);
    let tile = (2048usize, 5120usize, 1024usize, 1024usize);
    let mut crop: Option<(usize, usize, usize, usize)> = None;
    let c2 = |vc: ValleyConstruction| Knobs {
        valley: Some(vc),
        slope_floor_abs: Some(S_EQ),
        ..Knobs::passes(2)
    };
    let worlds: [(&str, &str, Knobs, Option<(&Skeleton, &ValleyConstruction, &CarveMasks)>); 5] = [
        ("C2/10", "c2_10", c2(v121), Some((&sk121, &v121, &m121))),
        ("C2/10 + basin base", "c2_10_basin", c2(v122), Some((&sk122, &v122, &m122))),
        (
            "C1 bare",
            "c1_bare",
            Knobs { valley: Some(ValleyConstruction::f121(F121_AGE_K, None)), ..Knobs::passes(2) },
            None,
        ),
        ("A1+B2", "a1b2", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, None),
        ("livré", "livre", Knobs::passes(2), None),
    ];
    for (label, file, knobs, masks) in worlds {
        let t = Instant::now();
        let g = build_field_seed(knobs, PSEED);
        let dr = inventory(&g, &ss);
        if let Some((sk, vc, mk)) = masks {
            let v = read_lakes(&dr, &g, sk, &mk.floor, &mk.carved, &trunk_near, vc, &ss);
            print_lakes(label, &v, w);
            if crop.is_none() {
                let id = dr.lake_map[31_634_531];
                let (mut x0, mut x1, mut y0, mut y1) = (w, 0, h, 0);
                for k in 0..n {
                    if id != 0 && dr.lake_map[k] == id {
                        x0 = x0.min(k % w);
                        x1 = x1.max(k % w);
                        y0 = y0.min(k / w);
                        y1 = y1.max(k / w);
                    }
                }
                let mg = 96;
                let (cx0, cy0) = (x0.saturating_sub(mg), y0.saturating_sub(mg));
                crop = Some((cx0, cy0, (x1 + mg).min(w - 1) - cx0, (y1 + mg).min(h - 1) - cy0));
                eprintln!(
                    "   brush crop, lake id {id} (bench frame x0, y0, w, h): {:?}",
                    crop.unwrap()
                );
            }
        }
        render(&format!("tile_2048_5120_{file}.png"), &g, &dr, tile, None, &ss);
        let cr = crop.expect("the C2/10 world runs first");
        render(&format!("brush_1000011_{file}.png"), &g, &dr, cr, None, &ss);
        if let Some((_, _, mk)) = masks {
            render(
                &format!("brush_1000011_{file}_valleys.png"),
                &g,
                &dr,
                cr,
                Some((&mk.carved, &mk.floor)),
                &ss,
            );
        }
        eprintln!("   ⏱ {label} {:.1} s", t.elapsed().as_secs_f64());
    }
    eprintln!("\n==========  end 122-B/C . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 124, part 2 — **PERMANENT: the brush lake is a bowl under the definition.** At the
/// canonical framing, the C2/10 world built with `ValleyConstruction::new` (basin base): the family-2
/// lake at Finding 122's brush lake (bench floor cell (5219, 3861), i.e. viz x 5219 − 768 = 4451)
/// must read **D_L < 3** and hold **no constructed floor cell below its level**. It fails loudly if
/// the definition ever drifts back to the sea base.
///
/// Run: cargo test -p ymir-core --release --test f122_brush -- --ignored f124_lake --nocapture
#[test]
#[ignore]
fn f124_lake_1000011_is_a_bowl() {
    let ss = SteinSteinParams::default();
    let vc = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
    let (_, masks) = carve(&pre, &sk, &vc, &ss);
    drop(pre);
    let trunk_src: Vec<usize> = (0..n).filter(|&k| sk.trunk[k]).collect();
    let near = bfs(&trunk_src, &|_| true, w, h, 2);
    let trunk_near: Vec<bool> = near.iter().map(|&d| d <= 2).collect();
    drop(near);
    let g = build_field_seed(
        Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
        PSEED,
    );
    let dr = inventory(&g, &ss);
    // the family-2 lake nearest the brush lake's floor, viz frame
    let (tx, ty) = (4451usize, 3861usize);
    let lakes = read_lakes(&dr, &g, &sk, &masks.floor, &masks.carved, &trunk_near, &vc, &ss);
    let l = lakes
        .iter()
        .filter(|l| l.id >= 1_000_001 && l.km2 >= 10.0)
        .min_by(|a, b| {
            let d = |l: &LakeRead| {
                let (fx, fy) = ((l.floor_cell % w) as f32, (l.floor_cell / w) as f32);
                (fx - tx as f32).powi(2) + (fy - ty as f32).powi(2)
            };
            d(a).total_cmp(&d(b))
        })
        .expect("a family-2 lake near the brush lake");
    let below = (0..n)
        .filter(|&k| {
            dr.lake_map[k] == l.id
                && masks.floor[k]
                && !sk.chi_m[k].is_nan()
                && sk.floor_m(k, vc.age_k) < l.level_m - 0.5
        })
        .count();
    eprintln!(
        "\n   the brush lake under the definition: id {} · {:.1} km² · level {:.1} m · floor cell ({}, {}) \
         · D_L **{:.2}** · carved share {:.1} % · constructed floor cells below its level **{below}**",
        l.id,
        l.km2,
        l.level_m,
        l.floor_cell % w,
        l.floor_cell / w,
        l.dl,
        100.0 * l.carved_share
    );
    let _ = h;
    assert!(l.dl < 3.0, "the brush lake reads D_L {:.2} ≥ 3: the brush is back", l.dl);
    assert_eq!(below, 0, "{below} constructed floor cells lie below the lake's level");
}
