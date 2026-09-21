//! ADR 0001 Finding 108 — the class is DIFFERENTIAL: `Δ(fill) > 50 m`, "the incision dug a hole".
//!
//! **Author, 2026-09-20.** The topological definition alone (`filled − raw > 50 m`, Finding 107)
//! counts PRE-EXISTING tectonic hollows. The complete definition is differential:
//!
//! ```text
//! Δ(fill) = (fill)_delivered − (fill)_pre-incision   at the body's FLOOR cell
//! ```
//!
//! ⚠️ **AND THE FILL ITSELF WAS BLIND, which is the addendum that makes this round possible.**
//! The priority flood takes `h ≤ sea` as base level, so it never raises a below-sea cell and
//! `filled − raw` reads **0** there whatever the hollow. That is the whole of Finding 107-D's
//! cleanest-looking number — *"fill p50 exactly 0.00 m on all six below-sea basins"* — and it
//! meant "the instrument cannot see there", not "they are not hollows". For every `wc == 2` body
//! the fill is therefore `spill_level − raw`, the Finding 85 levels
//! (`BasinSummary::spill_level_m`, metres).
//!
//! ⇒ **Finding 107-D's "16 → 12" and "the cut gate has four false positives" were WITHDRAWN
//! before this bench ran**, because they were computed on the blind instrument. Re-read here:
//! the delivered field's below-sea basins are real hollows with real sills, but SHALLOW ones —
//! their col sits **1.5 m** above the floor against 46–127 m before the incision — so the cut
//! class does carry four of them and Δ does drop all four. The conclusion survives; the number
//! it was built on did not.
//!
//! ## Pairing (rule 12 / 17, declared)
//!
//! Δ is read at the body's **floor cell** — `argmin(raw)` over the footprint — and both fills are
//! read **AT THAT SAME CELL**. Never by id: Finding 87's 623.6 m lake is `1000013` in its
//! inventory and `1000015` in Finding 107-D's, and only the floor cell **(3487, 5902)** is stable
//! (ADR L11573). That cell is tracked by name here so the archetype cannot fall out silently.
//!
//! ## The pre-incision STAGE, declared
//!
//! `Knobs::no_incision()` — tectonics + FBM, no stream power, no droplets. The bench ASSERTS it is
//! the Finding 76–80 coastal-authority stage by re-measuring its spur count (**20** on seed 1,
//! Finding 106's `coast authority` line). If that number moves, the stage is not the one the
//! dossier's coastal numbers were read on and every Δ here is against a different reference.
//!
//! Run: cargo test -p ymir-core --release --test f108_delta -- --ignored --nocapture

mod common;

#[allow(unused_imports)]
use common::{
    CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, f95_criteria, land_u16, majority,
    over_dug_depression, pct, set_dump, sorted, take_bodies, to_mask,
};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, below_sea_basin_lakes_infil, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const DELIVERED_P50_M: f32 = 424.2;
const WINDOW: [f32; 3] = [0.021, 0.024, 0.027];
/// ADR L11573 — the floor of the 623.6 m lake, the archetype of the class. Tracked by CELL
/// because its id is inventory-dependent (`1000013` at Finding 87, `1000015` at Finding 107-D).
const ARCHETYPE: (usize, usize) = (3487, 5902);
/// Finding 106's `coast authority` on seed 1's pre-incision field. The stage assertion.
const AUTHORITY_SEED1: usize = 20;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f108_eye");

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Finding 106's spur count, copied VERBATIM so "coast authority 20" is the same number.
fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

fn median(v: &[f32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    let s = sorted(v.to_vec());
    s[(s.len() - 1) / 2]
}

/// Share of lakes ≥ 1 km² with `D_L > 5`. Finding 102's raster bias: a disk reads 1.284.
/// Copied verbatim from `f106_constant.rs` rather than re-derived — same population, same number.
fn dl_share(lake_map: &[u32], w: usize, h: usize, cell_km2: f32) -> f32 {
    let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
    for (k, &id) in lake_map.iter().enumerate().take(w * h) {
        if id != 0 {
            by.entry(id).or_default().push(k);
        }
    }
    let (mut n, mut hi) = (0usize, 0usize);
    for cells in by.values() {
        let a = cells.len() as f32 * cell_km2;
        if a < 1.0 {
            continue;
        }
        let inside: HashSet<usize> = cells.iter().copied().collect();
        let mut per = 0usize;
        for &k in cells {
            let (x, y) = (k % w, k / w);
            for (dx, dy) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                let nx = (x as i32 + dx).rem_euclid(w as i32) as usize;
                let ny = (y as i32 + dy).rem_euclid(h as i32) as usize;
                if !inside.contains(&(ny * w + nx)) {
                    per += 1;
                }
            }
        }
        n += 1;
        let p = per as f32 * (400.0 / w as f32);
        if p / (2.0 * (std::f32::consts::PI * a).sqrt()) > 5.0 {
            hi += 1;
        }
    }
    100.0 * hi as f32 / n.max(1) as f32
}

/// ADR Finding 109 — the fill instrument moved to `common::fill_field_m` so this bench and
/// `f109_toggle` read the SAME one. The documentation of why it exists, and the correction to my
/// own diagnosis of `C1Lake::level_m`, live at the definition site.
fn fill_field(
    f: &GridF32,
    on: &C1DrainageConfig,
    ss: &SteinSteinParams,
    dclim: &DrainageClimate,
) -> (Vec<f32>, Vec<bool>) {
    common::fill_field_m(f, on, ss, dclim, DOMAIN_KM)
}

/// One scored body, with its differential column.
struct Row {
    id: u32,
    km2: f32,
    cut: f32,
    rim: f32,
    floor: usize,
    fill_del: f32,
    fill_pre: f32,
    delta: f32,
    under: bool,
    is_cut: bool,
}

impl Row {
    /// The gate. Rim carries a **declared ±1°**: a body in [29°, 31°] does not decide.
    fn is_delta(&self) -> bool {
        self.delta > 50.0 && self.rim > 30.0
    }
    fn rim_ambiguous(&self) -> bool {
        (29.0..=31.0).contains(&self.rim)
    }
}

#[allow(clippy::too_many_arguments)]
fn rows(
    f: &GridF32,
    pre: &GridF32,
    fill_del: &[f32],
    fill_pre: &[f32],
    under: &[bool],
    ss: &SteinSteinParams,
    on: &C1DrainageConfig,
    cell_km2: f32,
    n2m: f32,
    w: usize,
    h: usize,
) -> Vec<Row> {
    let dr = c1_drainage_windowed(f, None, on, ss, DOMAIN_KM);
    let bf = breach_monotone(f, &dr.flow.filled, &dr.lake_map, SEA, w, h);
    set_dump(true);
    let _ = f95_criteria(f, &bf, pre, DELIVERED_P50_M, ss, on, cell_km2, n2m, w, h);
    set_dump(false);
    take_bodies()
        .into_iter()
        .map(|b| {
            // ⚠️ The FLOOR cell, and both fills read AT IT. Never by id (Finding 81/85).
            let floor = *b
                .cells
                .iter()
                .min_by(|&&a, &&c| f.data[a].total_cmp(&f.data[c]))
                .expect("non-empty body");
            Row {
                id: b.id,
                km2: b.km2,
                cut: b.cut,
                rim: b.rim,
                floor,
                fill_del: fill_del[floor],
                fill_pre: fill_pre[floor],
                delta: fill_del[floor] - fill_pre[floor],
                under: under[floor],
                is_cut: b.is_cut,
            }
        })
        .collect()
}

fn report(rows: &[Row], label: &str, w: usize) -> (usize, usize, usize) {
    let (mut ord, mut sub, mut amb) = (0usize, 0usize, 0usize);
    for r in rows {
        if r.is_delta() {
            if r.under {
                sub += 1;
            } else {
                ord += 1;
            }
        }
        if r.rim_ambiguous() && r.delta > 50.0 {
            amb += 1;
        }
    }
    eprintln!(
        "   ⇒ **{label}**: Δ class **{}** = **{ord} ordinary + {sub} below-sea** of {} scanned \
         (**±{amb}** on the declared ±1° rim) · cut class {}",
        ord + sub,
        rows.len(),
        rows.iter().filter(|r| r.is_cut).count()
    );
    for r in rows.iter().filter(|r| r.is_delta() || r.is_cut) {
        eprintln!(
            "       body {:>7} {:>8.2} km² floor ({:>4},{:>4}) · cut {:>8.2} · fill del {:>8.2} \
             pre {:>8.2} ⇒ **Δ {:>9.2} m** · rim {:>5.1}°{} ⇒ cut {} / Δ {}",
            r.id,
            r.km2,
            r.floor % w,
            r.floor / w,
            r.cut,
            r.fill_del,
            r.fill_pre,
            r.delta,
            r.rim,
            if r.under { " [BELOW-SEA]" } else { "            " },
            if r.is_cut { "YES" } else { " no" },
            if r.is_delta() { "YES" } else { " no" }
        );
    }
    (ord, sub, amb)
}

#[test]
#[ignore]
fn f108_delta() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let t0 = Instant::now();
    let s2 = splitmix64(PSEED);
    let s3 = splitmix64(s2);
    eprintln!("\n==========  Finding 108 . the differential class  ==========");

    for (si, seed) in [PSEED, s2, s3].into_iter().enumerate() {
        eprintln!("\n========  SEED {} ({seed})  ========", si + 1);
        let pre = build_field_seed(Knobs::no_incision(), seed);
        let (w, h) = (pre.width, pre.height);

        // ── the STAGE assertion (rule 12/17) ────────────────────────────────
        let auth = spurs(&land_u16(&pre, &ss), w);
        eprintln!(
            "   PRE-INCISION STAGE: `Knobs::no_incision()` · coast authority **{auth} spurs**{}",
            if si == 0 {
                assert_eq!(
                    auth, AUTHORITY_SEED1,
                    "the pre-incision stage is NOT the Finding 76-80 coastal-authority field"
                );
                " ⇒ **IS the Finding 76-80 authority field** (20, Finding 106)"
            } else {
                ""
            }
        );

        // ⚠️ Rule 12 — EACH field gets its OWN climate. The first run computed it once on the
        // pre-incision field and handed it to every stage, so the two populations being
        // differenced were not built the same way.
        let cl_pre = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_pre = DrainageClimate {
            precip_internal: &cl_pre.precipitation,
            temperature: &cl_pre.temperature,
        };
        let (fill_pre, under_pre) = fill_field(&pre, &on, &ss, &dc_pre);
        eprintln!(
            "   pre-incision fill: **{}** below-sea cells · fill at the archetype floor (3487, \
             5902) = **{:.2} m**",
            under_pre.iter().filter(|&&u| u).count(),
            fill_pre[ARCHETYPE.1 * w + ARCHETYPE.0]
        );

        let delivered = build_field_seed(Knobs::passes(2), seed);
        let mut fields = vec![("delivered".to_string(), delivered)];
        for c in WINDOW {
            fields.push((
                format!("c = {c}"),
                build_field_seed(
                    Knobs { slope_floor_uk: Some(c), depression_floor: true, ..Knobs::passes(2) },
                    seed,
                ),
            ));
        }

        let mut rates: Vec<(String, f32, usize, usize)> = Vec::new();
        for (label, f) in &fields {
            eprintln!("\n--------  seed {} · {label}  --------", si + 1);
            let cl = c1_climate_placed(f, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
            let dclim = DrainageClimate {
                precip_internal: &cl.precipitation,
                temperature: &cl.temperature,
            };
            let (fill_del, under) = fill_field(f, &on, &ss, &dclim);
            let rs = rows(f, &pre, &fill_del, &fill_pre, &under, &ss, &on, cell_km2, n2m, w, h);
            let (ord, sub, amb) = report(&rs, label, w);
            rates.push((
                label.clone(),
                100.0 * (ord + sub) as f32 / rs.len().max(1) as f32,
                ord + sub,
                amb,
            ));

            // ── CONTROLS, on the delivered field only ───────────────────────
            if label == "delivered" {
                let dd = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
                let cutc: Vec<&Row> = rs.iter().filter(|r| r.is_cut).collect();
                let tight = cutc.iter().min_by(|a, b| a.delta.total_cmp(&b.delta));
                eprintln!(
                    "\n   ⛔ **CONTROL 1** — the sixteen of Finding 88-D4 read Δ > 50 m: cut class \
                     **{}**, of which **{}** pass ⇒ **{}** · tightest: body **{}** at **Δ {:.2} \
                     m**{}",
                    cutc.len(),
                    cutc.iter().filter(|r| r.delta > 50.0).count(),
                    if cutc.iter().all(|r| r.delta > 50.0) {
                        "TOTAL INCLUSION"
                    } else {
                        "**FAILS -- the definition excludes members of the sixteen**"
                    },
                    tight.map(|r| r.id).unwrap_or(0),
                    tight.map(|r| r.delta).unwrap_or(0.0),
                    tight.map(|r| if r.under { " [BELOW-SEA]" } else { "" }).unwrap_or("")
                );

                // the ARCHETYPE, by cell, never by id
                let ak = ARCHETYPE.1 * w + ARCHETYPE.0;
                match rs.iter().find(|r| r.floor == ak) {
                    Some(r) => eprintln!(
                        "   ⛔ **THE ARCHETYPE** (623.6 m lake, ADR L11573) is body **{}**, {:.2} \
                         km² · fill del **{:.2}** pre **{:.2}** ⇒ **Δ {:.2} m** · rim {:.1}° ⇒ in \
                         the Δ class: **{}**",
                        r.id,
                        r.km2,
                        r.fill_del,
                        r.fill_pre,
                        r.delta,
                        r.rim,
                        if r.is_delta() { "YES" } else { "**NO**" }
                    ),
                    None => {
                        let owner = rs.iter().find(|r| r.floor / w == ARCHETYPE.1);
                        eprintln!(
                            "   ⚠️ **THE ARCHETYPE'S floor cell (3487, 5902) is not any body's \
                             floor here** · fill del **{:.2} m** pre **{:.2} m** ⇒ Δ **{:.2} m** \
                             (nearest row on that scanline: {:?})",
                            fill_del[ak],
                            fill_pre[ak],
                            fill_del[ak] - fill_pre[ak],
                            owner.map(|r| r.id)
                        );
                    }
                }

                // ⛔ THE ARCHETYPE, READ AT ALL THREE STAGES. Finding 87 splits the 623.6 m depth
                // across them — pre-incision **+12.63 m**, eroded **+0.50 m**, breached
                // **−8.67 m**, *"the operator that carries the cell to −8.67 m is
                // `breach_monotone`"* (ADR L11590) — so "which stage" is not a detail here: if the
                // hollow only exists after the breach, a Δ read on the ERODED field cannot see it
                // and the archetype leaves the class for a reason that is about staging, not
                // provenance.
                let bfa = breach_monotone(f, &dd.flow.filled, &dd.lake_map, SEA, w, h);
                let (fill_bre, _) = fill_field(&bfa, &on, &ss, &dclim);
                let m_at = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
                eprintln!(
                    "   ⛔ **THE ARCHETYPE cell (3487, 5902) AT THREE STAGES** — altitude: \
                     pre-incision **{:.2} m** · eroded **{:.2} m** · breached **{:.2} m** (ADR \
                     L11573: +12.63 / +0.50 / −8.67) · FILL: pre **{:.2}** · eroded **{:.2}** · \
                     breached **{:.2}** ⇒ Δ(eroded) **{:.2} m** · Δ(breached) **{:.2} m**",
                    m_at(&pre, ak),
                    m_at(f, ak),
                    m_at(&bfa, ak),
                    fill_pre[ak],
                    fill_del[ak],
                    fill_bre[ak],
                    fill_del[ak] - fill_pre[ak],
                    fill_bre[ak] - fill_pre[ak]
                );

                eprintln!(
                    "\n   ⛔ **CONTROL 2** — the pre-existing tectonic hollows read Δ < 10 m:"
                );
                let mut big: Vec<&Row> = rs.iter().filter(|r| r.under).collect();
                big.sort_by(|a, b| b.km2.total_cmp(&a.km2));
                for r in big.iter().take(8) {
                    eprintln!(
                        "       body {:>7} {:>8.2} km² · cut {:>7.2} · fill del {:>8.2} pre \
                         {:>8.2} ⇒ **Δ {:>9.2} m** ⇒ **{}**",
                        r.id,
                        r.km2,
                        r.cut,
                        r.fill_del,
                        r.fill_pre,
                        r.delta,
                        if r.delta < 10.0 { "PRE-EXISTING" } else { "**DUG**" }
                    );
                }

                eprintln!("\n   ⛔ **CONTROL 3** — through-cut beds read Δ ≈ 0:");
                // ⚠️ Declared (rule 12): the author's control river is "le 92,97 m³/s", and that
                // number belongs to `f95_criteria`'s chain — spillways emitted from the BREACHED
                // field with a detected lake map, not from the eroded field with `None`. The
                // first run profiled a DIFFERENT spillway (8.30 m³/s, 3 cells) and its verdict was
                // therefore about a different bed. The chain is reproduced here so the discharge
                // printed below is the one the round names, or is visibly not.
                let dr_bf = c1_drainage_windowed(&bfa, Some(&dclim), &on, &ss, DOMAIN_KM);
                let pm = dr_bf.lake_map.clone();
                let bsd =
                    below_sea_basin_lakes_infil(&bfa, &dclim, &on, &ss, DOMAIN_KM, Some(&pm), None);
                for sw in
                    bsd.spillways.iter().max_by(|a, b| a.discharge_m3s.total_cmp(&b.discharge_m3s))
                {
                    let d: Vec<f32> = sw
                        .points
                        .iter()
                        .map(|&(x, y)| {
                            let k = y as usize * w + x as usize;
                            fill_del[k] - fill_pre[k]
                        })
                        .collect();
                    eprintln!(
                        "       Spillway max **{:.2} m³/s** · {} cells · **Δ p50 {:.4} m · max \
                         {:.4} m** ⇒ **{}**",
                        sw.discharge_m3s,
                        d.len(),
                        median(&d),
                        d.iter().copied().fold(f32::MIN, f32::max),
                        if median(&d).abs() <= 1.0 {
                            "NOT in the class"
                        } else {
                            "**IN -- FAILS**"
                        }
                    );
                }
                if let Some(i) = (0..dd.rivers.segments.len()).max_by(|&a, &b| {
                    dd.segment_discharge_m3s[a].total_cmp(&dd.segment_discharge_m3s[b])
                }) {
                    let d: Vec<f32> = dd.rivers.segments[i]
                        .points
                        .iter()
                        .map(|&(x, y)| {
                            let k = y as usize * w + x as usize;
                            fill_del[k] - fill_pre[k]
                        })
                        .collect();
                    eprintln!(
                        "       Watercourse max **{:.2} m³/s** · {} cells · **Δ p50 {:.4} m · max \
                         {:.4} m**",
                        dd.segment_discharge_m3s[i],
                        d.len(),
                        median(&d),
                        d.iter().copied().fold(f32::MIN, f32::max)
                    );
                }
                // ── C: what the CUT was measuring where the two disagree ────
                eprintln!(
                    "\n   ── C · bodies where cut and Δ disagree by > 50 m (cut high, Δ low) ──"
                );
                let mut div: Vec<&Row> = rs.iter().filter(|r| r.cut - r.delta > 50.0).collect();
                div.sort_by(|a, b| (b.cut - b.delta).total_cmp(&(a.cut - a.delta)));
                for r in div.iter().take(6) {
                    let ak = r.floor;
                    let a_km2 = dd.flow.accumulation.data[ak] * cell_km2;
                    eprintln!(
                        "       body {:>7} · cut **{:>8.2}** vs Δ **{:>8.2}** (gap {:>8.2}) · \
                         **catchment at the floor {:>9.2} km²** · receiver {} ⇒ **{}**",
                        r.id,
                        r.cut,
                        r.delta,
                        r.cut - r.delta,
                        a_km2,
                        if dd.flow.direction[ak] == ymir_core::terrain::flow::DIR_NONE {
                            "NONE"
                        } else {
                            "D8"
                        },
                        if a_km2 >= RELIEF_V1_A_C_KM2 {
                            "CHANNEL -- the cut was real incision, a through-cut bed"
                        } else {
                            "hillslope -- not a channel"
                        }
                    );
                }
            }

            // ── the EYE, Finding 106 block C: delivered and c = 0.024 ───────
            if label == "delivered" || label == "c = 0.024" {
                let dd = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
                let land: Vec<bool> = (0..w * h).map(|k| f.data[k] > SEA).collect();
                eprintln!(
                    "   PANEL: D_L>5 **{:.0} %** · Δ class **{}** · R8 **{:.4}**",
                    dl_share(&dd.lake_map, w, h, cell_km2),
                    rs.iter().filter(|r| r.is_delta()).count(),
                    aniso(f, &land, 16).r8
                );
                let mut g = GridF32::new(w, h, 0.0);
                let (az, alt) = (315f32.to_radians(), 45f32.to_radians());
                let (lx, ly, lz) = (alt.cos() * az.sin(), -alt.cos() * az.cos(), alt.sin());
                for y in 0..h {
                    for x in 0..w {
                        let (gx, gy) = f.gradient_at(x, y);
                        let (dx, dy) = (gx * n2m / CELL_M, gy * n2m / CELL_M);
                        g.data[y * w + x] = ((-dx * lx - dy * ly + lz)
                            / (dx * dx + dy * dy + 1.0).sqrt())
                        .clamp(0.0, 1.0);
                    }
                }
                for r in rs.iter().filter(|r| r.is_delta()) {
                    for &k in &[r.floor] {
                        for dy in -6i32..=6 {
                            for dx in -6i32..=6 {
                                let nx = (k % w) as i32 + dx;
                                let ny = (k / w) as i32 + dy;
                                if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                                    g.data[ny as usize * w + nx as usize] = 1.0;
                                }
                            }
                        }
                    }
                }
                std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
                let name = format!("seed{}_{}.png", si + 1, label.replace([' ', '='], ""));
                g.save_png_u8(&Path::new(OUT).join(&name)).expect("png");
                eprintln!("   wrote {name}");
            }
        }

        // ── the window under Δ, with the ±1-body tolerance ──────────────────
        let d_rate = rates[0].1;
        eprintln!(
            "\n   ── the window under Δ, gate `0.25 × delivered` = **{:.1} %** ──",
            0.25 * d_rate
        );
        for (label, rate, n, amb) in rates.iter().skip(1) {
            let gate = 0.25 * d_rate;
            eprintln!(
                "      **{label}**: Δ rate **{rate:.1} %** ({n} bodies, ±{amb} on the rim) ⇒ \
                 **{}**",
                if *rate <= gate + 0.001 {
                    "PASS"
                } else if 100.0 * (*n as f32 - 1.0 - *amb as f32).max(0.0) / 29.0 <= gate {
                    "**PASS within ±1 body**"
                } else {
                    "**FAIL**"
                }
            );
        }
    }
    eprintln!("\n==========  end Finding 108 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// **ADR Finding 108-D, rule 13 — the gate's two PERMANENT controls, both from measurements.**
///
/// These run in the normal suite: no field build, just the gate predicate against numbers this
/// round measured. If the gate is ever loosened or tightened, one of them fails.
#[test]
fn the_over_dug_gate_spares_the_river_and_keeps_the_tightest_body() {
    // NEGATIVE -- Finding 87-A's spillway, 92.97 m3/s over 61 cells (the length reproduced here
    // exactly). Its WORST cell reads delta = 0.0451 m. Paired with a wall-like rim so the test
    // proves the DEPTH clause carries it, not the slope clause: a bed can be steep-sided and
    // still not be a hollow.
    assert!(
        !over_dug_depression(0.0451, 74.9),
        "the gate must spare a through-cut bed however steep its sides"
    );
    // ... and it is not marginal: three orders of magnitude of headroom.
    assert!(0.0451 * 1000.0 < 50.0, "the river's margin is not three orders of magnitude");

    // POSITIVE -- body 51, the tightest UNAMBIGUOUS member of seed 1's class: delta 125.75 m,
    // rim 32.7 deg. Body 30 is deeper in the gate on delta (116.65 m) but its rim is 30.1 deg,
    // inside the declared +-1 deg band, so it is not allowed to anchor anything.
    assert!(over_dug_depression(125.75, 32.7), "the tightest measured member must stay in");
    assert!(
        !(29.0..=31.0).contains(&32.7),
        "the positive control must sit OUTSIDE the declared +-1 deg rim band"
    );

    // The clause that separates Finding 108 from Finding 107: a PRE-EXISTING hollow. Body 1000001
    // is 660.79 km2 with a real sill, and reads delta -75.73 m -- the erosion made it shallower.
    assert!(!over_dug_depression(-75.73, 31.0), "a pre-existing tectonic basin is not over-dug");
    // And the one the round named as the test of that clause: 1000012, delta +1.34 m.
    assert!(!over_dug_depression(1.34, 30.7), "1000012 is pre-existing, not dug");
}

/// ADR Finding 108 block C — **what was the CUT measuring** where it and Δ disagree?
///
/// The main test's block C came out EMPTY on all three seeds: no body has `cut − Δ > 50 m`. That
/// is already an answer, but the round names two bodies (10: cut 1 034 m against fill 2.88 m;
/// 21: 393 against 0.91) and a number beats an inference. Seed 1, delivered only, every body
/// ranked by `cut − Δ` in BOTH directions.
///
/// Run: cargo test -p ymir-core --release --test f108_delta -- --ignored block_c --nocapture
#[test]
#[ignore]
fn f108_block_c() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let cl_pre = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dc_pre = DrainageClimate {
        precip_internal: &cl_pre.precipitation,
        temperature: &cl_pre.temperature,
    };
    let (fill_pre, _) = fill_field(&pre, &on, &ss, &dc_pre);
    let f = build_field_seed(Knobs::passes(2), PSEED);
    let cl = c1_climate_placed(&f, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dclim =
        DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
    let (fill_del, under) = fill_field(&f, &on, &ss, &dclim);
    let dd = c1_drainage_windowed(&f, None, &on, &ss, DOMAIN_KM);
    let mut rs = rows(&f, &pre, &fill_del, &fill_pre, &under, &ss, &on, cell_km2, n2m, w, h);
    rs.sort_by(|a, b| (b.cut - b.delta).total_cmp(&(a.cut - a.delta)));
    eprintln!("\n==========  Finding 108-C . every body by `cut − Δ`, seed 1 delivered  =========");
    eprintln!(
        "   bodies with **cut − Δ > 50 m** (the round's hypothesis): **{}**",
        rs.iter().filter(|r| r.cut - r.delta > 50.0).count()
    );
    for r in rs.iter().take(6).chain(rs.iter().rev().take(3)) {
        let a_km2 = dd.flow.accumulation.data[r.floor] * cell_km2;
        eprintln!(
            "   body {:>7} {:>8.2} km² floor ({:>4},{:>4}){} · **cut {:>9.2}** · fill del \
             {:>9.2} pre {:>9.2} ⇒ **Δ {:>10.2}** · gap **{:>10.2}** · catchment at floor \
             {:>9.3} km² ⇒ {}",
            r.id,
            r.km2,
            r.floor % w,
            r.floor / w,
            if r.under { " [B-S]" } else { "      " },
            r.cut,
            r.fill_del,
            r.fill_pre,
            r.delta,
            r.cut - r.delta,
            a_km2,
            if a_km2 >= RELIEF_V1_A_C_KM2 { "CHANNEL" } else { "hillslope" }
        );
    }
    eprintln!("==========  end 108-C  ==========\n");
}
