//! ADR 0001 Finding 107-C — **the panels**, because the À VALIDER list was ids and numbers.
//!
//! The author asked where the bodies this round argues about can be looked at, and the honest
//! answer was **nowhere**: Finding 107 named `1000001`, `660.79 km²`, `rim 30.1°` and gave no
//! coordinate and no image. This bench produces the locator and the zooms.
//!
//! ⚠️ **The footprints are the GATE'S OWN** (`common::take_bodies`), captured inside the scoring
//! loop rather than re-derived. Finding 101-B had to declare that its simpler re-derived
//! inventory gave 26 bodies against the table's 29 and that only the A/B survived; a picture that
//! draws a different set of bodies from the one the gate scored is worse than no picture.
//!
//! ⚠️ **Rule 14b — the renderer is calibrated and stated, not "an image".** Three products per
//! view, because each can lie differently:
//!  * `_hillshade` — Finding 97's shaded relief, azimuth 315°, altitude 45°, reads the SLOPE, so
//!    it cannot band on altitude. This is the one that shows whether a rim is a wall or a ramp.
//!  * `_fill` — the DISCRIMINANT ITSELF, `filled − raw` in metres on a LOCAL ramp whose
//!    metres-per-code is printed. Black is 0 (through-cut); bright is deep below its own sill.
//!    **This is the quantity the gate reads**, drawn rather than tabulated.
//!  * `_mask` — the body footprint the gate scored, so nobody has to guess which hollow is meant.
//!
//! Run: cargo test -p ymir-core --release --test f107_eye -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, f95_criteria, set_dump, take_bodies};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const DELIVERED_P50_M: f32 = 424.2;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f107_eye");
/// Half-width of a zoom, SCALED TO ITS BODY.
///
/// ⚠️ A fixed 512 (50 km across) was the first version, and it made the panels unusable for the
/// small end: body 1 is 33.44 km², about 6 km across, so it filled **1.3 %** of its own frame and
/// read as a dim smudge. The rule below gives roughly three body-widths of context whatever the
/// size — 17 km for a 33 km² lake, 75 km for a 757 km² basin — so every panel frames its subject.
fn pad_for(cells: usize) -> usize {
    (((cells as f32).sqrt() * 1.5) as usize).clamp(96, 768)
}

/// Finding 97's hillshade, verbatim — the calibrated renderer, not a fresh one.
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

fn crop(g: &GridF32, cx: usize, cy: usize, pad: usize) -> GridF32 {
    let (w, h) = (g.width, g.height);
    let (x0, y0) = (cx.saturating_sub(pad), cy.saturating_sub(pad));
    let (x1, y1) = ((cx + pad).min(w), (cy + pad).min(h));
    let mut o = GridF32::new(x1 - x0, y1 - y0, 0.0);
    for y in y0..y1 {
        for x in x0..x1 {
            o.data[(y - y0) * o.width + (x - x0)] = g.data[y * w + x];
        }
    }
    o
}

/// Save on a LOCAL ramp and PRINT its calibration, so a reader can convert a code back to metres.
///
/// ⛔ **The first three renderings did NOT do this, and the calibration line they printed was a
/// fiction.** `GridF32::save_png_u8` clamps to `0..1` and multiplies by 255 — it does not stretch.
/// The fill grids are in NORMALISED altitude, so a hollow 875 m below its sill wrote code **20**
/// instead of 186: the discriminant was rendered at about a thirteenth of its dynamic range while
/// the console claimed `1.539 m/code`. **Rule 14b exists for exactly this** — a renderer that
/// states a calibration it does not perform is worse than an uncalibrated one, because the
/// statement is what a reader trusts. Caught by opening the output, not by reasoning about it.
fn save(g: &GridF32, name: &str, unit: &str, scale: f32) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    let (lo, hi) = g.data.iter().fold((f32::MAX, f32::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    let span = (hi - lo).max(f32::EPSILON);
    let mut out = GridF32::new(g.width, g.height, 0.0);
    for (o, &v) in out.data.iter_mut().zip(g.data.iter()) {
        *o = (v - lo) / span;
    }
    out.save_png_u8(&Path::new(OUT).join(name)).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!(
        "   wrote {name:<44} range **{:.2} .. {:.2} {unit}** ⇒ **{:.3} {unit}/code**",
        lo * scale,
        hi * scale,
        span * scale / 255.0
    );
}

#[test]
#[ignore]
fn f107_eye() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 107-C . the panels  ==========");

    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);

    for (label, knobs) in [
        ("delivered", Knobs::passes(2)),
        (
            "c0.027",
            Knobs { slope_floor_uk: Some(0.027), depression_floor: true, ..Knobs::passes(2) },
        ),
    ] {
        eprintln!("\n--------  {label}  --------");
        let f = build_field_seed(knobs, PSEED);
        let dr = c1_drainage_windowed(&f, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(&f, &dr.flow.filled, &dr.lake_map, SEA, w, h);

        // `filled − raw` on the ERODED field: the discriminant, as a field.
        let mut fill = GridF32::new(w, h, 0.0);
        for k in 0..w * h {
            fill.data[k] = (dr.flow.filled.data[k] - f.data[k]).max(0.0);
        }

        set_dump(true);
        let _ = f95_criteria(&f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        set_dump(false);
        let bodies = take_bodies();

        // ── the LOCATOR: where on the continent, at a glance ───────────────
        let hs = hillshade(&f, n2m);
        let mut loc = hs.clone();
        for b in bodies.iter().filter(|b| b.is_dug || b.is_cut) {
            // over-dug bodies burn WHITE, cut-only bodies burn BLACK, so a body in both classes
            // reads white and a disagreement is visible without a legend.
            let v = if b.is_dug { 1.0 } else { 0.0 };
            for &k in &b.cells {
                loc.data[k] = v;
            }
        }
        save(&loc, &format!("{label}_locator.png"), "shade", 1.0);
        save(&fill, &format!("{label}_fill_global.png"), "m", n2m);

        // ── the ZOOMS: EVERY body in either class, no selection at all ─────
        //
        // ⛔ **Selection was tried twice and failed the same way twice.** By AREA gave six panels
        // of which not one was an ordinary lake — the real class had no image. By area AND deep
        // share still dropped body 1000015, which the À VALIDER list then pointed the author at.
        // **Pointing someone at a file that does not exist is the same failure as giving them no
        // file, and no third selection rule fixes a rule-shaped problem.**
        //
        // ⇒ render EVERY classed body. There are ~20 of them, the zooms are small, and the
        // rendering loop costs nothing against the field builds. The assertion below is rule 13
        // applied to the renderer: **no body this round can discuss without a picture of it.**
        let of_interest: Vec<_> = bodies.iter().filter(|b| b.is_dug || b.is_cut).collect();
        let expected = of_interest.len();
        let mut rendered = 0usize;
        for b in of_interest.iter() {
            let (mut sx, mut sy) = (0usize, 0usize);
            for &k in &b.cells {
                sx += k % w;
                sy += k / w;
            }
            let (cx, cy) = (sx / b.cells.len(), sy / b.cells.len());
            eprintln!(
                "   body **{}** · {:.2} km² · centre **({cx}, {cy})** · cut **{:.2} m** · fill \
                 max **{:.2} m** · **deep share {:.1} %** · rim **{:.1}°** ⇒ cut {} / over-dug {}",
                b.id,
                b.km2,
                b.cut,
                b.fill,
                b.share,
                b.rim,
                if b.is_cut { "YES" } else { "no" },
                if b.is_dug { "YES" } else { "no" }
            );
            let pad = pad_for(b.cells.len());
            eprintln!("     frame **{} km across**", 2 * pad * 400 / 8192);
            let tag = format!("{label}_body{}", b.id);
            save(&crop(&hs, cx, cy, pad), &format!("{tag}_hillshade.png"), "shade", 1.0);
            save(&crop(&fill, cx, cy, pad), &format!("{tag}_fill.png"), "m", n2m);
            let mut mask = GridF32::new(w, h, 0.0);
            for &k in &b.cells {
                mask.data[k] = 1.0;
            }
            save(&crop(&mask, cx, cy, pad), &format!("{tag}_mask.png"), "in/out", 1.0);
            rendered += 1;
            for suffix in ["hillshade", "fill", "mask"] {
                let p = Path::new(OUT).join(format!("{tag}_{suffix}.png"));
                assert!(p.exists(), "{tag}_{suffix}.png was not written");
            }
        }
        // Rule 13 for the renderer: a body the round can name is a body the author can open.
        assert_eq!(rendered, expected, "{label}: {rendered} panels for {expected} classed bodies");
        eprintln!("   ⇒ **{rendered} of {expected} classed bodies rendered**, 3 files each");
    }
    eprintln!(
        "\n==========  end Finding 107-C . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
