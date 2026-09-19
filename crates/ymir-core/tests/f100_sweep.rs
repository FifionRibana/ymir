//! ADR 0001 Finding 100 — the 0.3–0.7 sweep, under the proposal template. **No production change.**
//!
//! Finding 99 found the campaign's best row at the LOWEST factor it measured (×0.7: relief +7.7 %
//! of the oracle, coast +2, canyons 3.4 %, erosion 60.6 %) and extrapolated that the oracle's
//! 488.1 m sits below it. This sweeps the range and asks the three thresholds the round names:
//! where the coast leaves +2/+3, where the relief crosses the oracle, where the canyons reopen.
//!
//! ⚠️ **Two columns of the template are load-bearing and are printed on every row**: the drainage
//! integrity (below-sea basins, biggest catchment — the pair that disqualified Finding 99's block
//! C after four criteria had called it a winner), and the erosion as a per-cell INTENSITY against
//! the reference at the same budget (rule 14 as amended at Finding 98).
//!
//! ⚠️ **And an arithmetic correction to the round's own block A**: it asks where the relief
//! "crosses the oracle (ratio 1.000; F99: 0.952 at ×0.7)". Those are different quantities — the
//! ratio Finding 99 printed is against the ×1 CALIBRATION (552.4 m), and the oracle is 488.1 m,
//! i.e. **ratio 0.8836**. Both columns are printed.
//!
//! Run: cargo test -p ymir-core --release --test f100_sweep -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Crit, Knobs, SEA, aniso, build_field, f95_criteria, land_u16, majority, pct, sorted,
    tile, to_mask,
};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const ORACLE_P50_M: f32 = 488.1;
const DELIVERED_P50_M: f32 = 424.2;
/// Finding 97 B's declared calibration, the ×1 point of this sweep.
const KS: f32 = 0.0451;
/// Finding 99's ×1 relief, the denominator of the "ratio to the ×1 calibration" column.
const CAL_P50_M: f32 = 552.4;
/// The tile of maximum Δ(R8), found at Finding 99 and reused so the image is comparable.
const TX: usize = 2048;
const TY: usize = 5120;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f100_sweep");

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

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

/// Local 3×3 standard deviation of altitude, metres. ⚠️ Needed BESIDE R8 for block B: a field with
/// nothing in it scores R8 ≈ 0 too, so R8 alone cannot tell "comb removed" from "texture lost".
/// σ holding while R8 falls = the comb went; σ collapsing with it = the relief went.
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

/// Drainage density, km of channel per km² of land — Finding 99's column, now with Finding 62's
/// mechanism behind it: the median drainage area itself moves out from under `A_c`.
fn density(f: &GridF32, on: &C1DrainageConfig, ss: &SteinSteinParams, cell_km2: f32) -> (f32, f32) {
    let (w, h) = (f.width, f.height);
    let d = c1_drainage_windowed(f, None, on, ss, DOMAIN_KM);
    let (mut len, mut area) = (0f64, 0f64);
    let mut a_all: Vec<f32> = Vec::new();
    for k in 0..w * h {
        if f.data[k] <= SEA {
            continue;
        }
        area += cell_km2 as f64;
        let a = d.flow.accumulation.data[k] * cell_km2;
        if k % 7 == 0 {
            a_all.push(a);
        }
        if a < RELIEF_V1_A_C_KM2 {
            continue;
        }
        let dir = d.flow.direction[k];
        if dir == DIR_NONE {
            continue;
        }
        let diag = D8_DX[dir as usize] != 0 && D8_DY[dir as usize] != 0;
        len += if diag { CELL_M as f64 * 1.4142136 } else { CELL_M as f64 } / 1000.0;
    }
    ((len / area.max(1e-9)) as f32, pct(&sorted(a_all), 0.50))
}

fn save(g: &GridF32, name: &str) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    let p = Path::new(OUT).join(name);
    g.save_png_u8(&p).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!("   wrote {}", p.display());
}

#[test]
#[ignore]
fn f100_sweep() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 100 . the 0.3-0.7 sweep  ==========");
    eprintln!(
        "   k_time HELD at x1 . k_s(x1) = {KS:.4} . oracle {ORACLE_P50_M} m = ratio \
         **{:.4}** of the x1 calibration ({CAL_P50_M} m) -- NOT 1.000",
        ORACLE_P50_M / CAL_P50_M
    );

    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let delivered = build_field(Knobs::passes(2));
    let auth = spurs(&land_u16(&pre, &ss), w);
    let intensity = |f: &GridF32| -> f64 {
        let c: Vec<usize> = (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] > SEA).collect();
        c.iter().map(|&k| ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64).sum::<f64>()
            / c.len().max(1) as f64
    };
    let i_del = intensity(&delivered);
    let (dens_del, a50_del) = density(&delivered, &on, &ss, cell_km2);
    let (dens_pre, a50_pre) = density(&pre, &on, &ss, cell_km2);
    eprintln!(
        "   prerequisites in {:.1} s . delivered cuts {i_del:.1} m/cell . density pre \
         **{dens_pre:.3}** (median A {a50_pre:.4} km2) vs delivered **{dens_del:.3}** (median A \
         **{a50_del:.4}** km2) -- Finding 62's collapse, re-measured",
        t0.elapsed().as_secs_f64()
    );

    // ── the sweep ────────────────────────────────────────────────────────────
    let mut fields: Vec<(f32, GridF32, f64)> = Vec::new();
    for mult in [0.3f32, 0.4, 0.5, 0.6, 0.7] {
        let t = Instant::now();
        let f = build_field(Knobs {
            slope_floor_uk: Some(KS * mult),
            depression_floor: true,
            ..Knobs::passes(2)
        });
        fields.push((mult, f, t.elapsed().as_secs_f64()));
    }

    let mut rows: Vec<(f32, Crit, f64, f32, i64, f32)> = Vec::new();
    for (label, f, cost) in std::iter::once(("delivered".to_string(), &delivered, 0.0))
        .chain(fields.iter().map(|(m, f, c)| (format!("A1+B2 x{m}"), f, *c)))
    {
        let d = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(f, &d.flow.filled, &d.lake_map, SEA, w, h);
        let c = f95_criteria(f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
        let r8 = aniso(f, &land, 16).r8;
        let inten = intensity(f);
        let sp = spurs(&land_u16(&bf, &ss), w) as i64 - auth as i64;
        let (dens, a50) = density(f, &on, &ss, cell_km2);
        eprintln!(
            "   **{label}**: relief **{:.1} m** = ratio **{:.4} of x1** / **{:.4} of ORACLE** \
             ({:+.1} %) . coast **{sp:+}** . canyons **{}/{} ({:.1} %)** . lakes {:.2} % . \
             **erosion {:.1} m/cell = {:.1} %** . R8 **{r8:.4}** . density {dens:.3} (median A \
             {a50:.4}) . ⛔ **basins {} ({:+.0} %) . catchment {:.0} km2 ({:+.0} %)** . {cost:.1} s",
            c.p50,
            c.p50 / CAL_P50_M,
            c.p50 / ORACLE_P50_M,
            100.0 * (c.p50 - ORACLE_P50_M) / ORACLE_P50_M,
            c.klass,
            c.scanned,
            100.0 * c.klass as f32 / c.scanned.max(1) as f32,
            c.lake_pct,
            inten,
            100.0 * inten / i_del,
            c.below_sea_basins,
            100.0 * (c.below_sea_basins as f32 / 11.0 - 1.0),
            c.max_catchment_km2,
            100.0 * (c.max_catchment_km2 / 9779.0 - 1.0)
        );
        let m = if label == "delivered" { 0.0 } else { label[7..].parse().unwrap_or(0.0) };
        rows.push((m, c, inten, r8, sp, dens));
    }

    // ── the three thresholds, read off the sweep ────────────────────────────
    eprintln!("\n-- the three thresholds the round asked for --");
    let sw: Vec<&(f32, Crit, f64, f32, i64, f32)> = rows.iter().filter(|r| r.0 > 0.0).collect();
    let coast_break = sw.iter().rev().find(|r| r.4 > 3).map(|r| r.0);
    let oracle_cross = sw.iter().find(|r| r.1.p50 <= ORACLE_P50_M).map(|r| r.0);
    let canyon_open = sw
        .iter()
        .rev()
        .find(|r| 100.0 * r.1.klass as f32 / r.1.scanned.max(1) as f32 > 6.0)
        .map(|r| r.0);
    let rule14 = sw.iter().rev().find(|r| 100.0 * r.2 / i_del >= 70.0).map(|r| r.0);
    eprintln!(
        "   coast leaves +2/+3 at factor **{}** . relief crosses the ORACLE at **{}** . canyons \
         exceed 6 % at **{}** . rule 14's 70 % reached at **{}**",
        coast_break.map_or("NOT IN 0.3-0.7".into(), |v| format!("x{v}")),
        oracle_cross.map_or("NOT IN 0.3-0.7".into(), |v| format!("x{v}")),
        canyon_open.map_or("NOT IN 0.3-0.7 (holds)".into(), |v| format!("x{v}")),
        rule14.map_or("NOT IN 0.3-0.7".into(), |v| format!("x{v}"))
    );
    let integrity_ok = sw.iter().all(|r| {
        (r.1.below_sea_basins as f32 / 11.0 - 1.0).abs() <= 0.2 || r.1.below_sea_basins <= 15
    });
    eprintln!(
        "   ⛔ DRAINAGE INTEGRITY across the sweep: basins {:?} (delivered 11) . catchments {:?} \
         (delivered 9779) ⇒ **{}**",
        sw.iter().map(|r| r.1.below_sea_basins).collect::<Vec<_>>(),
        sw.iter().map(|r| r.1.max_catchment_km2 as i32).collect::<Vec<_>>(),
        if integrity_ok {
            "HELD -- no candidate disqualified"
        } else {
            "**BROKEN -- disqualified**"
        }
    );

    // ── B — the image, with sigma beside R8 ─────────────────────────────────
    let l_del: Vec<bool> = (0..n).map(|k| delivered.data[k] > SEA).collect();
    let best = sw
        .iter()
        .filter(|r| r.4 <= 5)
        .min_by(|a, b| {
            (a.1.p50 - ORACLE_P50_M).abs().partial_cmp(&(b.1.p50 - ORACLE_P50_M).abs()).unwrap()
        })
        .map(|r| r.0)
        .unwrap_or(0.7);
    eprintln!("\n-- B . the image at ({TX}, {TY}) . best coast-clean row = x{best} --");
    for (label, f) in std::iter::once(("delivered".to_string(), &delivered)).chain(
        fields
            .iter()
            .filter(|(m, _, _)| *m == 0.7 || *m == best)
            .map(|(m, f, _)| (format!("x{m}"), f)),
    ) {
        let (g, gl) = tile(f, &l_del, TX, TY, 1024);
        // the WESTERN THIRD, the region the round asks about
        let (gw, gwl) = tile(f, &l_del, TX, TY, 320);
        eprintln!(
            "   {label:<12} whole tile: R8 **{:.4}** sigma **{:.2} m** . western third (320^2): \
             R8 **{:.4}** sigma **{:.2} m**",
            aniso(&g, &gl, 16).r8,
            sigma_p50(&g, &gl, n2m),
            aniso(&gw, &gwl, 16).r8,
            sigma_p50(&gw, &gwl, n2m)
        );
        save(&hillshade(&g, n2m), &format!("{label}_hillshade.png"));
    }
    eprintln!(
        "\n==========  end Finding 100 . total {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
