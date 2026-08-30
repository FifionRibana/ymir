//! C-2 placement — structural verification on the REAL C1 state (not a table of
//! values: placement correctness is judged on STRUCTURE). Reports, for two seeds:
//!   - ARCS: offset from the O-C margin (mean/std in km) — a tight std means the
//!     edifices form a line parallel to the trench at a consistent offset;
//!   - HOTSPOT CHAINS: age monotonic along the chain and the members colinear
//!     (age increasing opposite to plate motion — the discriminating causal test);
//!   - RIFTS: fraction sitting on Divergent∧Continental cells;
//!   - per-mechanism counts, compared across two seeds (coherent with tectonics,
//!     not a random draw).
//!
//! Run: cargo test -p ymir-core --test volcanism_placement --release -- --ignored --nocapture

use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::boundary_classification::{
    BoundaryType, classify_boundaries, oc_override_seed_mask,
};
use ymir_core::tectonics_c1::closures::volcanism::placement::place_edifices;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, VolcanoSetting};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const TORUS_KM: f32 = 1024.0; // C1_DOMAIN_KM — the coarse tectonic torus span.

fn build_state(seed_val: u64) -> (C1State, PlateKinematics) {
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, seed_val, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});
    (state, kin)
}

fn torus_km_dist(a: (f32, f32), b: (f32, f32)) -> f32 {
    let mut dx = (a.0 - b.0).abs();
    let mut dy = (a.1 - b.1).abs();
    if dx > 0.5 {
        dx = 1.0 - dx;
    }
    if dy > 0.5 {
        dy = 1.0 - dy;
    }
    (dx * dx + dy * dy).sqrt() * TORUS_KM
}

fn report_seed(seed_val: u64) {
    let (state, kin) = build_state(seed_val);
    let seed = WorldSeed::new(seed_val);
    let cfg = VolcanismConfig { enabled: true, domain_km: TORUS_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, TORUS_KM, &cfg);

    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let boundary = classify_boundaries(&state.plate_id, &kin);
    let arc_mask = oc_override_seed_mask(&boundary, &state.plate_id, &state.plate_type);
    let margin_uv: Vec<(f32, f32)> = (0..ny)
        .flat_map(|j| (0..nx).map(move |i| (i, j)))
        .filter(|&(i, j)| arc_mask.get(i, j))
        .map(|(i, j)| ((i as f32 + 0.5) / nx as f32, (j as f32 + 0.5) / ny as f32))
        .collect();

    let arcs: Vec<_> = edifices.iter().filter(|e| e.setting == VolcanoSetting::Arc).collect();
    let rifts: Vec<_> = edifices.iter().filter(|e| e.setting == VolcanoSetting::Rift).collect();
    let hots: Vec<_> = edifices.iter().filter(|e| e.setting == VolcanoSetting::Hotspot).collect();

    eprintln!("\n===== seed {seed_val} =====");
    eprintln!(
        "counts: arcs {} | rifts {} | hotspot edifices {} | O-C margin cells {} | total {}",
        arcs.len(),
        rifts.len(),
        hots.len(),
        margin_uv.len(),
        edifices.len()
    );

    // ARCS: two offset measures.
    //  (a) nearest-neighbour distance (underestimates on a curved/clustered margin);
    //  (b) PERPENDICULAR to the local boundary tangent — the honest offset, robust
    //      to which margin cell is nearest. The tangent at the nearest margin cell
    //      M0 is the direction to its nearest OTHER margin cell; the normal is its
    //      perpendicular, and the perpendicular offset is |(P − M0)·normal|.
    if !arcs.is_empty() && margin_uv.len() >= 2 {
        // periodic component helper (wrap to [-0.5,0.5])
        let wrap = |mut d: f32| {
            if d > 0.5 {
                d -= 1.0;
            }
            if d < -0.5 {
                d += 1.0;
            }
            d
        };
        let mut nn = Vec::new();
        let mut perp = Vec::new();
        for e in &arcs {
            // nearest margin cell M0
            let (m0, _) = margin_uv
                .iter()
                .map(|&m| (m, torus_km_dist(e.center_uv, m)))
                .fold(((0.0, 0.0), f32::MAX), |acc, x| if x.1 < acc.1 { x } else { acc });
            nn.push(torus_km_dist(e.center_uv, m0));
            // tangent = toward nearest OTHER margin cell
            let m1 = margin_uv
                .iter()
                .filter(|&&m| m != m0)
                .map(|&m| (m, torus_km_dist(m0, m)))
                .fold(((0.0f32, 0.0f32), f32::MAX), |acc, x| if x.1 < acc.1 { x } else { acc })
                .0;
            let (tx, ty) = (wrap(m1.0 - m0.0), wrap(m1.1 - m0.1));
            let tl = (tx * tx + ty * ty).sqrt().max(1e-9);
            let (nxn, nyn) = (-ty / tl, tx / tl); // normal to the tangent
            let (px, py) = (wrap(e.center_uv.0 - m0.0), wrap(e.center_uv.1 - m0.1));
            perp.push((px * nxn + py * nyn).abs() * TORUS_KM);
        }
        let stat = |v: &[f32]| {
            let m = v.iter().sum::<f32>() / v.len() as f32;
            let s = (v.iter().map(|x| (x - m).powi(2)).sum::<f32>() / v.len() as f32).sqrt();
            (m, s)
        };
        let (nnm, nns) = stat(&nn);
        let (pm, ps) = stat(&perp);
        // Applied offset magnitude, measured directly foot→edifice (construction).
        use ymir_core::tectonics_c1::closures::volcanism::placement::arc_sites;
        let applied: Vec<f32> = arc_sites(&state, &boundary, TORUS_KM, &cfg)
            .iter()
            .map(|&(foot, pos)| torus_km_dist(foot, pos))
            .collect();
        let (am, asd) = stat(&applied);
        eprintln!(
            "  ARCS: applied foot→edifice {:.0}±{:.0} km | nearest-margin {:.0}±{:.0} | PERP-to-tangent {:.0}±{:.0} (target {:.0})",
            am, asd, nnm, nns, pm, ps, cfg.trench_arc_offset_km
        );
    }

    // RIFTS: fraction sitting on a Divergent∧Continental cell.
    if !rifts.is_empty() {
        let on = rifts
            .iter()
            .filter(|e| {
                let i = ((e.center_uv.0 * nx as f32) as usize).min(nx - 1);
                let j = ((e.center_uv.1 * ny as f32) as usize).min(ny - 1);
                matches!(boundary.boundary_type.get(i, j), BoundaryType::Divergent)
                    && matches!(state.plate_type.get(i, j), PlateType::Continental)
            })
            .count();
        eprintln!(
            "  RIFTS: {}/{} on Divergent∧Continental ({:.0}%)",
            on,
            rifts.len(),
            100.0 * on as f32 / rifts.len() as f32
        );
    }

    // HOTSPOT CHAINS: grouped in contiguous runs of `hotspot_chain_len`. Report the
    // age sequence, monotonicity, and per-step direction (colinearity ⇒ a straight
    // age-progressive chain aligned with plate motion).
    let clen = cfg.hotspot_chain_len;
    for (c, chunk) in hots.chunks(clen).enumerate() {
        if chunk.len() < 2 {
            continue;
        }
        let ages: Vec<f32> = chunk.iter().map(|e| e.age_frac).collect();
        let monotonic = ages.windows(2).all(|w| w[1] > w[0]);
        // Direction of each consecutive step (should be constant ⇒ colinear).
        let steps: Vec<(f32, f32)> = chunk
            .windows(2)
            .map(|w| {
                let mut dx = w[1].center_uv.0 - w[0].center_uv.0;
                let mut dy = w[1].center_uv.1 - w[0].center_uv.1;
                if dx > 0.5 {
                    dx -= 1.0;
                }
                if dx < -0.5 {
                    dx += 1.0;
                }
                if dy > 0.5 {
                    dy -= 1.0;
                }
                if dy < -0.5 {
                    dy += 1.0;
                }
                let l = (dx * dx + dy * dy).sqrt().max(1e-9);
                (dx / l, dy / l)
            })
            .collect();
        // Max angular deviation between steps (0 ⇒ perfectly colinear).
        let d0 = steps[0];
        let max_dev_deg = steps
            .iter()
            .map(|s| (s.0 * d0.0 + s.1 * d0.1).clamp(-1.0, 1.0).acos().to_degrees())
            .fold(0.0f32, f32::max);
        eprintln!(
            "  HOTSPOT chain {c}: ages [{}] monotonic={monotonic} | colinear (max step dev {:.1}°) | active first = {}",
            ages.iter().map(|a| format!("{a:.2}")).collect::<Vec<_>>().join(", "),
            max_dev_deg,
            chunk[0].active
        );
    }
}

#[test]
#[ignore]
fn c2_placement_structure() {
    eprintln!("\n=== C-2 PLACEMENT STRUCTURE (real C1 state) ===");
    report_seed(10481999410520546993);
    report_seed(42);
}
