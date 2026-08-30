//! C-3 lithology-source coverage probe (read-only, the decision measurement). The
//! four fractions of the continent that a causal, NON-NOISE K field could
//! structure, so option 1 (3-class now) vs option 2 (C1 terrane plumbing) is
//! ARITHMETIC, not preference:
//!   - craton BASE (build_phase_1_1_cratonic_mask geometric rule: plate seed x <
//!     nx/2, approximated here by the plate CENTROID x — Voronoi cells cluster on
//!     their seed), NOT the FBM-refined shield (~7 %);
//!   - rift (age ≈ 0, causal, rift-spawned);
//!   - volcanic footprints (C-2 placement, physical basal area);
//!   - residual "generic continental" (no causal structure available).
//! Fractions are area-preserving under bilinear upscale, so the coarse figures are
//! the HD figures; volcanic is computed from physical basal km².
//!
//! Run: cargo test -p ymir-core --test c3_lithology_probe --release -- --ignored --nocapture

use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const DOMAIN_KM: f32 = 400.0;

fn probe(seed_val: u64) {
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, seed_val, &Phase2InitParams::default());
    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let cell_km2 = (DOMAIN_KM / nx as f32).powi(2);
    let x_threshold = nx as f64 / 2.0;

    // Craton BASE measured AT INIT (centroid ≈ seed before advection; on the final
    // state the plates have moved and the centroid proxy degenerates). Fraction of
    // INIT continental cells whose plate's centroid x < nx/2 (the geometric rule).
    let craton_base_pc = {
        let np = state.plate_id.data().iter().copied().max().unwrap_or(0) as usize + 1;
        let (mut sx, mut c) = (vec![0.0f64; np], vec![0u32; np]);
        for j in 0..ny {
            for i in 0..nx {
                let p = state.plate_id.get(i, j) as usize;
                sx[p] += i as f64;
                c[p] += 1;
            }
        }
        let (mut cont0, mut cra0) = (0u32, 0u32);
        for j in 0..ny {
            for i in 0..nx {
                if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                    cont0 += 1;
                    let p = state.plate_id.get(i, j) as usize;
                    if c[p] > 0 && sx[p] / (c[p] as f64) < x_threshold {
                        cra0 += 1;
                    }
                }
            }
        }
        100.0 * cra0 as f32 / cont0.max(1) as f32
    };

    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});

    let (mut cont, mut shield, mut rift) = (0u32, 0u32, 0u32);
    for j in 0..ny {
        for i in 0..nx {
            if !matches!(state.plate_type.get(i, j), PlateType::Continental) {
                continue;
            }
            cont += 1;
            if state.cratonic_mask.get(i, j) {
                shield += 1;
            }
            if state.age.get(i, j) < 1.0 {
                rift += 1;
            }
        }
    }
    let cont_km2 = cont as f32 * cell_km2;

    // Volcanic footprint: sum of edifice basal discs (physical km²), as % of land.
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &WorldSeed::new(seed_val), DOMAIN_KM, &volc);
    let volc_km2: f32 =
        edifices.iter().map(|e| std::f32::consts::PI * (e.basal_diameter_km * 0.5).powi(2)).sum();

    let pc = |n: u32| 100.0 * n as f32 / cont.max(1) as f32;
    let volc_pc = 100.0 * volc_km2 / cont_km2.max(1.0);
    // Genuinely-causal-and-non-noise coverage = rift + volcanic (the craton base is
    // a geometric placeholder; the shield is noise-refined — neither is causal).
    let residual = 100.0 - pc(rift) - volc_pc;
    eprintln!("\n=== seed {seed_val} ===  continental {cont} cells = {cont_km2:.0} km²");
    eprintln!(
        "  craton BASE (geometric PLACEHOLDER, seed_x<nx/2) {craton_base_pc:.0}%  | shield (noise-refined) {:.0}%  | rift (age~0, causal) {:.0}%  | volcanic (causal) {:.1}% ({:.0} km², {} edifices)",
        pc(shield),
        pc(rift),
        volc_pc,
        volc_km2,
        edifices.len(),
    );
    eprintln!(
        "  → residual GENERIC continental (NO causal, non-noise structure) ≈ {residual:.0}%  [craton base excluded: geometric placeholder, not physics]"
    );
}

#[test]
#[ignore]
fn c3_coverage() {
    eprintln!("\n=== C-3 LITHOLOGY COVERAGE (option 1 vs 2, by the numbers) ===");
    probe(10481999410520546993);
    probe(42);
}
