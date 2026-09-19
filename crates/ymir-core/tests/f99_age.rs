//! ADR 0001 Finding 99 — is the age the FACTOR of B2? **No production change.**
//!
//!  * **A** — sweep `(U/K)^{1/n}` at `k_time` fixed: does the relief follow the factor while the
//!    coast, the canyons and R8 hold? With Finding 98's A1+B2 ×1 row reproduced as the control.
//!  * **B** — the image and the drainage density, on the tile of **maximum Δ(R8)** (Finding 97's
//!    own correction: selecting by max R8 finds tiles anisotropic in *every* field and is blind to
//!    the difference).
//!  * **C** — ⛔ the round asks for `RELIEF_V3_DIFFUSION` 0.08 → 1.28. **That is ×16, and the
//!    dossier has already run it and disqualified it, twice.** Finding 45 measured the equivalent
//!    patch: mean land altitude 685 → 702 m (+2.5 %, *"small and in the unhelpful direction"*) and
//!    **below-sea basins 43 → 1994 (×46)**, biggest catchment 110 → 48 km². Its verdict: *"the
//!    patch did not change the altitude — it destroyed the drainage"*. The production comment at
//!    `stream_power.rs:937` says why: with `diffuse_channels = true` a 16× Laplacian *"backfills
//!    the channels the incision just cut, and a backfilled channel IS a closed depression"*.
//!    So block C is run **with Finding 45's own column added** (`Crit::below_sea_spillways`), because
//!    measuring only R8 / relief / erosion would reproduce Finding 45's original error exactly.
//!
//! Run: cargo test -p ymir-core --release --test f99_age -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Crit, Knobs, SEA, aniso, build_field, f95_criteria, land_u16, majority, pct, sorted,
    tile, to_mask,
};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::{
    RELIEF_V1_A_C_KM2, RELIEF_V1_K, RELIEF_V3_DIFFUSION, RELIEF_V3_K_MULT,
};
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
/// Finding 97 B's declared calibration: the delivered channels' own Flint intercept.
const KS: f32 = 0.0451;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f99_age");

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

/// Shaded relief (light at 315°, 45°). Reads the SLOPE, so altitude quantisation cannot band it —
/// the renderer Finding 97 made the dossier's reference for texture questions.
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

/// Drainage density: channel length (cells with `A ≥ A_c`, each counted as its own link length)
/// per km² of land, per quadrant. A closure that removed the network would show here.
fn density(f: &GridF32, on: &C1DrainageConfig, ss: &SteinSteinParams, cell_km2: f32) -> [f32; 4] {
    let (w, h) = (f.width, f.height);
    let d = c1_drainage_windowed(f, None, on, ss, DOMAIN_KM);
    let (mut len, mut area) = ([0f64; 4], [0f64; 4]);
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if f.data[k] <= SEA {
                continue;
            }
            let q = (if y < h / 2 { 0 } else { 2 }) + usize::from(x >= w / 2);
            area[q] += cell_km2 as f64;
            if d.flow.accumulation.data[k] * cell_km2 < RELIEF_V1_A_C_KM2 {
                continue;
            }
            let dir = d.flow.direction[k];
            if dir == DIR_NONE {
                continue;
            }
            let diag = D8_DX[dir as usize] != 0 && D8_DY[dir as usize] != 0;
            len[q] += if diag { CELL_M as f64 * 1.4142136 } else { CELL_M as f64 } / 1000.0;
        }
    }
    [0, 1, 2, 3].map(|q| (len[q] / area[q].max(1e-9)) as f32)
}

/// The 1024² tile of maximum **Δ(R8)** between two fields. ⚠️ Finding 97 selected by max R8 and
/// recorded the error: that finds tiles anisotropic in EVERY field and is blind to the difference.
fn max_delta_tile(a: &GridF32, b: &GridF32, la: &[bool], lb: &[bool]) -> (usize, usize, f32, f32) {
    let (w, h) = (a.width, a.height);
    let mut best = (0usize, 0usize, 0f32, 0f32);
    let mut bd = -1f32;
    for by in (0..h).step_by(1024) {
        for bx in (0..w).step_by(1024) {
            let (ga, gla) = tile(a, la, bx, by, 1024);
            let (gb, glb) = tile(b, lb, bx, by, 1024);
            let (ra, rb) = (aniso(&ga, &gla, 16), aniso(&gb, &glb, 16));
            if ra.windows < 64 || rb.windows < 64 {
                continue;
            }
            let d = (ra.r8 - rb.r8).abs();
            if d > bd {
                bd = d;
                best = (bx, by, ra.r8, rb.r8);
            }
        }
    }
    best
}

fn save(g: &GridF32, name: &str) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    let p = Path::new(OUT).join(name);
    g.save_png_u8(&p).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!("   wrote {}", p.display());
}

#[test]
#[ignore]
fn f99_age() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let k0 = RELIEF_V1_K * RELIEF_V3_K_MULT;
    eprintln!("\n==========  Finding 99 . is the age the FACTOR of B2?  ==========");
    eprintln!(
        "   shipped: k = {k0:.0}, k_time = {:.0} (HELD at x1 all round) . k_s calibration \
         {KS:.4} . RELIEF_V3_DIFFUSION = {RELIEF_V3_DIFFUSION}",
        k0 * 2.0
    );

    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let delivered = build_field(Knobs::passes(2));
    let auth = spurs(&land_u16(&pre, &ss), w);
    let cut = |f: &GridF32, k: usize| ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64;
    let land_cells =
        |f: &GridF32| (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] > SEA).count() as f64;
    // rule 14 as amended at Finding 98: per-cell INTENSITY, against the reference at the SAME
    // budget (which is `delivered` throughout, since k_time is held at x1 all round)
    let intensity = |f: &GridF32| -> f64 {
        (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] > SEA).map(|k| cut(f, k)).sum::<f64>()
            / land_cells(f).max(1.0)
    };
    let i_del = intensity(&delivered);
    eprintln!(
        "   prerequisites + delivered in {:.1} s . delivered cuts **{i_del:.1} m per land cell**",
        t0.elapsed().as_secs_f64()
    );

    // ══ A — the factor sweep, k_time held ════════════════════════════════════
    let mut fields: Vec<(String, GridF32, f64)> = Vec::new();
    for mult in [0.7f32, 0.85, 1.0, 1.2] {
        let t = Instant::now();
        let f = build_field(Knobs {
            slope_floor_uk: Some(KS * mult),
            depression_floor: true,
            ..Knobs::passes(2)
        });
        fields.push((format!("A1+B2 factor x{mult}"), f, t.elapsed().as_secs_f64()));
    }
    // ══ C — the diffusion, with Finding 45's column ══════════════════════════
    let t = Instant::now();
    let diff = build_field(Knobs { diffusion: Some(1.28), ..Knobs::passes(2) });
    let t_diff = t.elapsed().as_secs_f64();

    let mut rows: Vec<(String, Crit, f64, f64, f32, i64)> = Vec::new();
    for (label, f, cost) in std::iter::once(("delivered".to_string(), &delivered, 0.0))
        .chain(fields.iter().map(|(l, f, c)| (l.clone(), f, *c)))
        .chain(std::iter::once(("C diffusion 1.28 (x16)".to_string(), &diff, t_diff)))
    {
        let d = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(f, &d.flow.filled, &d.lake_map, SEA, w, h);
        let c = f95_criteria(f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
        let r8 = aniso(f, &land, 16).r8;
        let inten = intensity(f);
        let sp = spurs(&land_u16(&bf, &ss), w) as i64 - auth as i64;
        eprintln!(
            "   **{label}**: relief **{:.1} m ({:+.1} % vs oracle)** . coast **{sp:+}** . canyons \
             **{}/{} ({:.1} %)** . lakes {:.2} % . **erosion {:.1} m/cell = {:.1} % of delivered** \
             . R8 **{r8:.4}** . ⛔ **below-sea basins {} . max catchment {:.0} km2** . {cost:.1} s",
            c.p50,
            100.0 * (c.p50 - ORACLE_P50_M) / ORACLE_P50_M,
            c.klass,
            c.scanned,
            100.0 * c.klass as f32 / c.scanned.max(1) as f32,
            c.lake_pct,
            inten,
            100.0 * inten / i_del,
            c.below_sea_spillways,
            c.max_catchment_km2
        );
        rows.push((label, c, cost, inten, r8, sp));
    }

    // ══ B — the image and the drainage density ═══════════════════════════════
    eprintln!("\n-- B . drainage density (channel km per km2 of land, per quadrant) --");
    for (label, f) in [
        ("pre-incision", &pre),
        ("delivered", &delivered),
        ("A1+B2 x1", &fields[2].1),
        ("A1+B2 x0.7", &fields[0].1),
        ("C diffusion x16", &diff),
    ] {
        let q = density(f, &on, &ss, cell_km2);
        let mean = (q[0] + q[1] + q[2] + q[3]) / 4.0;
        eprintln!(
            "   {label:<18} NW {:.3} NE {:.3} SW {:.3} SE {:.3} . **mean {mean:.3}**",
            q[0], q[1], q[2], q[3]
        );
    }

    let l_del: Vec<bool> = (0..n).map(|k| delivered.data[k] > SEA).collect();
    let l_ab: Vec<bool> = (0..n).map(|k| fields[2].1.data[k] > SEA).collect();
    let (bx, by, r_del, r_ab) = max_delta_tile(&delivered, &fields[2].1, &l_del, &l_ab);
    eprintln!(
        "\n   the tile of maximum DELTA(R8) -- the selection Finding 97 said to use: **({bx}, \
         {by})** . delivered R8 {r_del:.4} vs A1+B2 {r_ab:.4} (delta {:+.4})",
        r_ab - r_del
    );
    // ⚠️ the fourth panel the round asks for is the ORACLE, which does not exist as a field:
    // Finding 95 recorded numbers and rebuilding the 300-pass field costs 1 h 43, refused by the
    // author. The pre-incision field is the declared substitute.
    for (label, f) in [
        ("delivered", &delivered),
        ("a1b2_x1", &fields[2].1),
        ("a1b2_x0.7", &fields[0].1),
        ("pre_incision_SUBSTITUTE_for_oracle", &pre),
        ("c_diffusion_x16", &diff),
    ] {
        let (g, _) = tile(f, &l_del, bx, by, 1024);
        save(&hillshade(&g, n2m), &format!("{label}_hillshade.png"));
    }

    eprintln!("\n-- the sweep, read as a dial --");
    let base = rows.iter().find(|r| r.0.contains("x1")).map(|r| r.1.p50).unwrap_or(1.0);
    for (label, c, _cost, inten, r8, sp) in &rows {
        eprintln!(
            "   {label}: relief {:.1} m = **{:.3} of the x1 calibration** . coast {sp:+} . R8 \
             {r8:.4} . erosion {:.1} % . below-sea basins {}",
            c.p50,
            c.p50 / base,
            100.0 * inten / i_del,
            c.below_sea_spillways
        );
    }
    eprintln!(
        "\n==========  end Finding 99 . total {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
