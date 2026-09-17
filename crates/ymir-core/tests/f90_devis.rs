//! ADR 0001 Finding 90 block D — the QUOTE for a base-level floor in `breach_monotone`.
//!
//! **Not a patch, not a gate, nothing promoted.** The round asked what such a floor would cost;
//! this measures the population it would act on and the population it would break, on the shipped
//! field, so the decision is taken against numbers.
//!
//! The candidate line is one clause inside `breach_monotone_protected`:
//! `let mut target = (height.data[nb] - EPS).max(sea_level + eps_floor);`
//! — the same bound Finding 83 put on the relaxation TARGET, applied to the carve target.
//!
//! Run: cargo test -p ymir-core --release --test f90_devis -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone};

const DOMAIN_KM: f32 = 400.0;

#[test]
#[ignore]
fn f90_devis() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let eps_floor_norm = 0.5 / n2m; // Finding 83's RELIEF_V3_BASE_LEVEL_M, in norm units
    eprintln!("\n==========  Finding 90 · D — the quote for a floor in the breach  ==========");
    eprintln!(
        "   the candidate bound: `target.max(sea + {:.1} m)` = `max(.., {:.3e} norm)` — the same \
         0.5 m Finding 83 put on the relaxation target",
        0.5,
        SEA + eps_floor_norm
    );

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;

    let raw = build_field(Knobs::passes(2));
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let d = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let br = breach_monotone(&raw, &d.flow.filled, &d.lake_map, SEA, w, h);
    let m = |g: &ymir_core::grid::GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);

    // ── the population the floor would act on ───────────────────────────────
    let drowned: Vec<usize> = (0..n).filter(|&k| raw.data[k] > SEA && br.data[k] <= SEA).collect();
    let below_floor: Vec<usize> =
        drowned.iter().copied().filter(|&k| br.data[k] < SEA + eps_floor_norm).collect();
    eprintln!(
        "\n── D1 · the population ──\n   cells the breach drowns: **{}** ({:.2} km²) · of those, \
         **{} ({:.1} %)** end BELOW the candidate floor and would be clamped",
        drowned.len(),
        drowned.len() as f32 * cell_km2,
        below_floor.len(),
        100.0 * below_floor.len() as f64 / drowned.len().max(1) as f64
    );

    // ── the trenches: connected components of drowned cells ─────────────────
    let dset: std::collections::HashSet<usize> = drowned.iter().copied().collect();
    let mut seen = vec![false; n];
    let mut comps: Vec<(usize, f32)> = Vec::new(); // (size, deepest metres below sea)
    for &s in &drowned {
        if seen[s] {
            continue;
        }
        let mut stack = vec![s];
        seen[s] = true;
        let mut size = 0usize;
        let mut deepest = 0.0f32;
        while let Some(k) = stack.pop() {
            size += 1;
            deepest = deepest.max(-m(&br, k));
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dd in 0..8 {
                let nx = (x + D8_DX[dd]).rem_euclid(w as i32) as usize;
                let ny = (y + D8_DY[dd]).rem_euclid(h as i32) as usize;
                let nk = ny * w + nx;
                if dset.contains(&nk) && !seen[nk] {
                    seen[nk] = true;
                    stack.push(nk);
                }
            }
        }
        comps.push((size, deepest));
    }
    comps.sort_by(|a, b| b.0.cmp(&a.0));
    let sizes: Vec<f32> = comps.iter().map(|c| c.0 as f32).collect();
    let deeps: Vec<f32> = comps.iter().map(|c| c.1).collect();
    eprintln!(
        "── D2 · the trenches ──\n   **{} distinct drowned trenches** · size in cells p50 {:.0} \
         p90 {:.0} max {:.0} · their deepest point below sea p50 **{:.2} m** p90 {:.2} m max \
         **{:.2} m**",
        comps.len(),
        pct(&sorted(sizes.clone()), 0.50),
        pct(&sorted(sizes.clone()), 0.90),
        sizes.first().copied().unwrap_or(0.0),
        pct(&sorted(deeps.clone()), 0.50),
        pct(&sorted(deeps.clone()), 0.90),
        deeps.iter().copied().fold(0.0f32, f32::max)
    );
    let would_fail = comps.iter().filter(|c| c.1 > 0.5).count();
    eprintln!(
        "   of those, **{would_fail} ({:.1} %)** are deeper than the candidate floor, so their \
         carve could NOT reach its target and the pit would fall to the FILL mop-up instead — a \
         flat pond, which is exactly the option Finding 13 rejected (*\"the pit becomes a flat \
         pond rather than a drained channel (the author asked for carve)\"*)",
        100.0 * would_fail as f64 / comps.len().max(1) as f64
    );

    // ── what it would NOT disturb ───────────────────────────────────────────
    let wc_e = water_class(&raw, SEA);
    let wc_b = water_class(&br, SEA);
    let enc_b = (0..n).filter(|&k| wc_b[k] == 2).count();
    let enc_new = (0..n).filter(|&k| wc_b[k] == 2 && wc_e[k] != 2).count();
    eprintln!(
        "── D3 · what it would NOT disturb ──\n   the breach makes {enc_new} of {enc_b} \
         enclosed-below-sea cells (**{:.1} %**), so bounding it moves the below-sea population — \
         and every Finding 71–88 chain figure with it — by **at most {:.1} %**",
        100.0 * enc_new as f64 / enc_b.max(1) as f64,
        100.0 * enc_new as f64 / enc_b.max(1) as f64
    );

    // ── the guard that would have to be re-proved ───────────────────────────
    let pits_e =
        (0..n).filter(|&k| raw.data[k] > SEA && d.flow.filled.data[k] > raw.data[k] + 1e-6).count();
    eprintln!(
        "── D4 · the guard that would have to be re-proved ──\n   land cells the priority-flood \
         considers pits in the ERODED field: **{pits_e}** ({:.2} % of land). Finding 14's \
         permanent non-ignored guard `flow::tests::breach_leaves_no_interior_pit` and the \
         `river_climbs` acceptance test are what a floor puts at risk: the FILL mop-up would still \
         guarantee monotonicity, so the guards would stay GREEN while the terrain silently gained \
         **{would_fail}** flat ponds. **That is the real cost, and no existing guard sees it.**",
        100.0 * pits_e as f64 / (0..n).filter(|&k| raw.data[k] > SEA).count().max(1) as f64
    );
    eprintln!("\n==========  end Finding 90 · D (quote only, nothing promoted)  ==========\n");
}
