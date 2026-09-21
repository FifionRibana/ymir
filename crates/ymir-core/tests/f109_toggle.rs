//! ADR 0001 Finding 109 — **the guards, run BEFORE the checkbox is allowed to exist.**
//!
//! The author judges in the viz and then in LL — the Finding 82 path. This bench is what licenses
//! showing him anything: if the production toggle does not reproduce the bench field **bit for
//! bit**, the world he looks at is not the world Findings 106–108 measured, and nothing should be
//! shown at all.
//!
//! **Finding 83's pattern, invoked deliberately: the hash IS the guard table.** When the field is
//! identical every gate measured on it holds by construction.
//!
//! | guard | what it forbids |
//! |---|---|
//! | `None` == delivered, 3 seeds | a seam that is not inert when off |
//! | `Absolute(0.024)` == the bench at c = 0.024, 3 seeds | a production path reading a different stage — **the test Finding 105's two-pass form FAILED**, 5.9 / 7.8 / 6.6 M cells differing |
//! | `Absolute(0.024)` != delivered, 3 seeds | the branch silently never firing, which would make guard 1 pass and guard 2 vacuous |
//! | the Finding 108 gates on the PRODUCTION output | the criteria chain reading something other than the heightmap |
//!
//! ⚠️ **Why guard 2 can pass here and could not at Finding 105.** The calibrated form measures
//! `k_s` on an intermediate stage no bench knob reproduces (between `:451` and `:506` the pipeline
//! also stamps craters and recomputes slope). **A constant has no calibration pass, therefore no
//! stage to misalign**: `upscale_from_c1` sets two knobs and makes one call to `once`, exactly
//! like the `None` branch.
//!
//! Run: cargo test -p ymir-core --release --test f109_toggle -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, fill_field_m, land_u16, majority,
    over_dug_depression, set_dump, take_bodies, to_mask,
};
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, DrainageClimate, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const DELIVERED_P50_M: f32 = 424.2;
/// The middle of Finding 106's window, re-read at Finding 108 under the differential class.
const S_EQ: f32 = 0.024;

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// FNV-1a over the raw f32 bits — Finding 105's hash, so the two rounds are comparable.
fn hash(f: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for v in &f.data {
        for b in v.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
    }
    h
}

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

fn components(f: &GridF32, w: usize, h: usize) -> usize {
    let n = w * h;
    let wc = water_class(f, SEA);
    let mut seen = vec![false; n];
    let mut c = 0usize;
    for s in 0..n {
        if wc[s] != 2 || seen[s] {
            continue;
        }
        c += 1;
        let mut q = vec![s];
        seen[s] = true;
        while let Some(k) = q.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nx = (x + dx).rem_euclid(w as i32) as usize;
                    let ny = (y + dy).rem_euclid(h as i32) as usize;
                    let nk = ny * w + nx;
                    if wc[nk] == 2 && !seen[nk] {
                        seen[nk] = true;
                        q.push(nk);
                    }
                }
            }
        }
    }
    c
}

/// Finding 102's shoreline-development share, raster bias declared (a disk reads 1.284).
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

#[test]
#[ignore]
fn f109_toggle() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let s2 = splitmix64(PSEED);
    let s3 = splitmix64(s2);
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 109 . the toggle's guards  ==========");

    let mut all_ok = true;
    for (si, seed) in [PSEED, s2, s3].into_iter().enumerate() {
        eprintln!("\n========  SEED {} ({seed})  ========", si + 1);

        let t = Instant::now();
        let delivered = build_field_seed(Knobs::passes(2), seed);
        let t_del = t.elapsed().as_secs_f64();
        let h_del = hash(&delivered);

        // ── GUARD 1: the seam is inert when off ─────────────────────────────
        let off = build_field_seed(Knobs { slope_floor_abs: None, ..Knobs::passes(2) }, seed);
        let g1 = hash(&off) == h_del;
        all_ok &= g1;
        eprintln!(
            "   GUARD 1 · `slope_floor: None` **{:016x}** vs delivered **{h_del:016x}** ⇒ **{}**",
            hash(&off),
            if g1 { "BIT-IDENTICAL" } else { "**DIFFERENT -- the seam is not inert**" }
        );

        // ── GUARD 2: production == bench, bit for bit ───────────────────────
        let t = Instant::now();
        let prod =
            build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, seed);
        let t_prod = t.elapsed().as_secs_f64();
        let bench = build_field_seed(
            Knobs { slope_floor_uk: Some(S_EQ), depression_floor: true, ..Knobs::passes(2) },
            seed,
        );
        let (hp, hb) = (hash(&prod), hash(&bench));
        let g2 = hp == hb;
        let diff = (0..prod.data.len()).filter(|&k| prod.data[k] != bench.data[k]).count();
        all_ok &= g2;
        eprintln!(
            "   GUARD 2 · production **{hp:016x}** vs bench c = {S_EQ} **{hb:016x}** ⇒ **{}** \
             ({diff} cells differ) · Finding 105's two-pass form failed this with 5.9-7.8 M",
            if g2 { "BIT-IDENTICAL" } else { "**DIFFERENT -- DO NOT SHOW** " }
        );

        // ── GUARD 3: the branch actually fires ──────────────────────────────
        // Without this, a `slope_floor` that never reaches the incision would make GUARD 1 pass
        // and GUARD 2 vacuous -- both hashes would simply be `delivered`.
        let g3 = hp != h_del;
        all_ok &= g3;
        eprintln!(
            "   GUARD 3 · production vs delivered ⇒ **{}** ({} cells differ) — the branch fires",
            if g3 { "DIFFERENT, as required" } else { "**IDENTICAL -- the branch never fired**" },
            (0..prod.data.len()).filter(|&k| prod.data[k] != delivered.data[k]).count()
        );
        eprintln!(
            "   ⏱ delivered **{t_del:.1} s** · toggled **{t_prod:.1} s** ⇒ **{:+.1} s** (Finding \
             105's calibrated form cost +56 to +79 s)",
            t_prod - t_del
        );

        // ── the Finding 108 gates, re-read on the PRODUCTION output ─────────
        if si == 0 {
            let (w, h) = (prod.width, prod.height);
            let pre = build_field_seed(Knobs::no_incision(), seed);
            let auth = spurs(&land_u16(&pre, &ss), w);
            let cl_pre =
                c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
            let dc_pre = DrainageClimate {
                precip_internal: &cl_pre.precipitation,
                temperature: &cl_pre.temperature,
            };
            let (fill_pre, _) = fill_field_m(&pre, &on, &ss, &dc_pre, DOMAIN_KM);
            let cl = c1_climate_placed(&prod, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
            let dclim = DrainageClimate {
                precip_internal: &cl.precipitation,
                temperature: &cl.temperature,
            };
            let (fill_del, under) = fill_field_m(&prod, &on, &ss, &dclim, DOMAIN_KM);
            let dr = c1_drainage_windowed(&prod, None, &on, &ss, DOMAIN_KM);
            let bf = breach_monotone(&prod, &dr.flow.filled, &dr.lake_map, SEA, w, h);
            set_dump(true);
            let cr = common::f95_criteria(
                &prod,
                &bf,
                &pre,
                DELIVERED_P50_M,
                &ss,
                &on,
                cell_km2,
                n2m,
                w,
                h,
            );
            set_dump(false);
            let bodies = take_bodies();
            let (mut ord, mut sub) = (0usize, 0usize);
            for b in &bodies {
                let floor = *b
                    .cells
                    .iter()
                    .min_by(|&&a, &&c| prod.data[a].total_cmp(&prod.data[c]))
                    .expect("non-empty");
                if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                    if under[floor] {
                        sub += 1;
                    } else {
                        ord += 1;
                    }
                }
            }
            let land: Vec<bool> = (0..w * h).map(|k| prod.data[k] > SEA).collect();
            eprintln!(
                "\n   ── the Finding 108 gates on the PRODUCTION output (seed 1, s_eq = {S_EQ}) ──\
                 \n   Δ class **{}** = {ord} ordinary + {sub} below-sea of {} scanned · coast \
                 **{:+}** · D_L>5 **{:.0} %** · R8 **{:.4}** · components **{}** · catchment \
                 **{:.0} km²** · lakes {:.2} % · `to_nothing` {}",
                ord + sub,
                cr.scanned,
                spurs(&land_u16(&bf, &ss), w) as i64 - auth as i64,
                dl_share(&dr.lake_map, w, h, cell_km2),
                aniso(&prod, &land, 16).r8,
                components(&prod, w, h),
                cr.max_catchment_km2,
                cr.lake_pct,
                cr.to_nothing
            );
            // ⚠️ The round asks for BOTH stages, not a choice: `f95_criteria` reads the ERODED
            // field, while the `.ymir` export ships the BREACHED one. Finding 108 measured the
            // archetype at Δ 0.02 m eroded against 623.56 m breached, so the two stages disagree
            // about the deepest object on the continent and the export must say which it is.
            let (fill_bre, _) = fill_field_m(&bf, &on, &ss, &dclim, DOMAIN_KM);
            let deep = (0..w * h).filter(|&k| fill_bre[k] > 50.0).count();
            eprintln!(
                "   ── the BREACHED stage, reported beside it (what the `.ymir` export ships) ──\
                 \n   cells more than 50 m below their own sill: eroded **{}** · breached **{}** \
                 ⇒ the export carries the breach's hollows, not the incision's",
                (0..w * h).filter(|&k| fill_del[k] > 50.0).count(),
                deep
            );
        }
    }
    eprintln!(
        "\n==========  end Finding 109 . guards green on three seeds: **{}** . {:.1} s  ==========\n",
        if all_ok { "YES -- the toggle may be shown" } else { "**NO -- SHOW NOTHING**" },
        t0.elapsed().as_secs_f64()
    );
    assert!(all_ok, "a guard failed: the toggle must not be put in front of the author");
}

/// ADR Finding 109 — **the cost, re-measured with nothing else on the machine.**
///
/// The three-seed guard run reported +1.8 / +17.3 / +15.2 s, and I do not believe its own numbers:
/// `f109_export` was compiling concurrently, and the DELIVERED builds — which should be the same
/// work three times — ranged 93.6 to 141.1 s, a spread with no physical cause. My prediction said
/// "more than +10 s means something runs twice and I should find out what"; the honest answer is
/// that the machine was the second thing. Seed 1 only, two builds, alternating, three rounds so
/// the ordering cannot carry the result.
///
/// Run ALONE: cargo test -p ymir-core --release --test f109_toggle -- --ignored cost --nocapture
#[test]
#[ignore]
fn f109_cost() {
    eprintln!("\n==========  Finding 109 . the cost, uncontaminated  ==========");
    let (mut d, mut t) = (Vec::new(), Vec::new());
    for r in 0..3 {
        let a = Instant::now();
        let _ = build_field_seed(Knobs::passes(2), PSEED);
        d.push(a.elapsed().as_secs_f64());
        let b = Instant::now();
        let _ = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        t.push(b.elapsed().as_secs_f64());
        eprintln!(
            "   round {} · delivered **{:.1} s** · toggled **{:.1} s** ⇒ **{:+.1} s**",
            r + 1,
            d[r],
            t[r],
            t[r] - d[r]
        );
    }
    let md = d.iter().sum::<f64>() / 3.0;
    let mt = t.iter().sum::<f64>() / 3.0;
    eprintln!(
        "   ⇒ mean delivered **{md:.1} s** · toggled **{mt:.1} s** ⇒ **{:+.1} s** ({:+.1} %) · \
         spread within delivered **{:.1} s** (the contamination measure)",
        mt - md,
        100.0 * (mt - md) / md,
        d.iter().cloned().fold(f64::MIN, f64::max) - d.iter().cloned().fold(f64::MAX, f64::min)
    );
    eprintln!("==========  end  ==========\n");
}
