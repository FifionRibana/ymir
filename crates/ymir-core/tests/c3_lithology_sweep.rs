//! C-3 lithology spread SWEEP (closures roadmap §3) — the measurement that CHOOSES
//! the soft↔hard erodibility multiplier, per the author's directive: the spread is a
//! MEASUREMENT, not a prediction. Runs the WHOLE production chain
//! (`upscale_from_c1_with_progress`, the export recipe: relief-v3 stream-power
//! incision, droplets off), lithology OFF as baseline then ON at ×3/×10/×30/×100 soft,
//! and reports PER LITHOLOGY CLASS (hard / rift-soft / volcaniclastic):
//!   - local relief (1 km window, m), median slope (°), steep-cell share (>30°),
//!   - channel incision depth (ridge−channel over channel cells, m) → the W/D handle,
//! plus the GLOBAL invariants that must survive (closed depressions — C-1 must not
//! regress — and the land fraction).
//!
//! Method rule 3 (ADR 0001): the bench reproduces the WHOLE chain, not just the
//! incision — the K field is built and threaded exactly as production does it.
//!
//! Two effects separated BY DESIGN: hard = ×1.0 (the relief-v3 reference), soft
//! ABOVE. So the ~80 % hard bulk is the reference at every multiplier → the hard-class
//! columns should be ~constant across the sweep (no global slowdown to disentangle);
//! ONLY the contrast (soft/volcanic columns) moves. The bench prints the hard column
//! every row so that invariance is visible.
//!
//! Run: cargo test -p ymir-core --test c3_lithology_sweep --release -- --ignored --nocapture

use ymir_core::erosion::stream_power::StreamPowerConfig;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::lithology::{self, LithologyConfig};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{Edifice, VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const FULL_RANGE_M: f32 = 11302.0;

const HARD: u8 = 0;
const SOFT: u8 = 1;
const VOLC: u8 = 2;

/// Per-HD-cell lithology CLASS (0 hard / 1 rift-soft / 2 volcaniclastic), built from
/// the same causal signals as the K field but as crisp classes (so the sweep can
/// group metrics by class without value collisions when soft_mult == volc_mult).
fn class_map(
    state: &C1State,
    edifices: &[Edifice],
    so: [f64; 2],
    ss: f64,
    km_per_hd: f32,
    w: usize,
    h: usize,
) -> Vec<u8> {
    // Coarse rift mask (1.0 continental age<1) → upscale → threshold.
    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let mut rift = GridF32::new(nx, ny, 0.0);
    for j in 0..ny {
        for i in 0..nx {
            if matches!(state.plate_type.get(i, j), PlateType::Continental)
                && (state.age.get(i, j) as f32) < 1.0
            {
                rift.set(i, j, 1.0);
            }
        }
    }
    let up = lithology::upscale_k_to_hd(&rift, w, h, so, ss);
    let mut cls = vec![HARD; w * h];
    for (k, &v) in up.iter().enumerate() {
        if v > 0.5 {
            cls[k] = SOFT;
        }
    }
    // Volcaniclastic footprints stamped over the top (basal discs).
    let (so32, ss32) = ([so[0] as f32, so[1] as f32], ss as f32);
    for e in edifices {
        let fx = (e.center_uv.0 - so32[0]).rem_euclid(1.0) / ss32;
        let fy = (e.center_uv.1 - so32[1]).rem_euclid(1.0) / ss32;
        if fx >= 1.0 || fy >= 1.0 {
            continue;
        }
        let (cx, cy) = (fx * w as f32, fy * h as f32);
        let rb = (e.basal_diameter_km * 0.5) / km_per_hd;
        if rb < 1.0 {
            continue;
        }
        let (i0, i1) =
            ((cx - rb).floor().max(0.0) as usize, ((cx + rb).ceil() as usize).min(w - 1));
        let (j0, j1) =
            ((cy - rb).floor().max(0.0) as usize, ((cy + rb).ceil() as usize).min(h - 1));
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d <= rb {
                    cls[j * w + i] = VOLC;
                }
            }
        }
    }
    cls
}

/// Per-class morphometrics on the eroded field: (cells, median local relief m,
/// median slope °, steep>30° share ‰, median channel incision depth m).
struct ClassStats {
    cells: u64,
    relief_m: f32,
    slope_deg: f32,
    steep_permille: f32,
    incision_m: f32,
}

fn per_class_stats(field: &GridF32, cls: &[u8], class: u8, km_per_cell: f32) -> ClassStats {
    let (w, h) = (field.width, field.height);
    let cell_m = km_per_cell * 1000.0;
    let r = 5i32; // ~ (5·km_per_cell) window; at 2048² ≈ 1 km
    let mut relief = Vec::new();
    let mut slopes = Vec::new();
    let mut incision = Vec::new();
    let (mut cells, mut steep) = (0u64, 0u64);
    // Channel threshold via flow accumulation (99th pct over land).
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let flow =
        compute_flow(field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let mut sorted = flow.accumulation.data.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let chan_thr = sorted[(sorted.len() as f64 * 0.99) as usize];
    for y in r..(h as i32 - r) {
        for x in r..(w as i32 - r) {
            let k = y as usize * w + x as usize;
            if cls[k] != class || field.data[k] <= 0.5 {
                continue;
            }
            cells += 1;
            let (gx, gy) = field.gradient_at(x as usize, y as usize);
            let slope_m = ((gx * gx + gy * gy).sqrt() * FULL_RANGE_M) / cell_m;
            let deg = slope_m.atan().to_degrees();
            slopes.push(deg);
            if deg > 30.0 {
                steep += 1;
            }
            let (mut mn, mut mx) = (f32::MAX, f32::MIN);
            for dy in -r..=r {
                for dx in -r..=r {
                    let v = field.data[(y + dy) as usize * w + (x + dx) as usize];
                    mn = mn.min(v);
                    mx = mx.max(v);
                }
            }
            relief.push((mx - mn) * FULL_RANGE_M);
            if flow.accumulation.data[k] >= chan_thr {
                incision.push((mx - field.data[k]) * FULL_RANGE_M); // ridge − channel
            }
        }
    }
    let med = |mut v: Vec<f32>| -> f32 {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    ClassStats {
        cells,
        relief_m: med(relief),
        slope_deg: med(slopes),
        steep_permille: 1000.0 * steep as f32 / cells.max(1) as f32,
        incision_m: med(incision),
    }
}

fn count_closed_depressions(field: &GridF32) -> usize {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let flow =
        compute_flow(field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let thr = 0.1f32 / FULL_RANGE_M;
    let mut seen = vec![false; n];
    let mut count = 0usize;
    let mut stack = Vec::new();
    for s in 0..n {
        if seen[s] || flow.filled.data[s] - field.data[s] <= thr {
            continue;
        }
        count += 1;
        seen[s] = true;
        stack.push(s);
        while let Some(k) = stack.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if !seen[nk] && flow.filled.data[nk] - field.data[nk] > thr {
                            seen[nk] = true;
                            stack.push(nk);
                        }
                    }
                }
            }
        }
    }
    count
}

fn land_fraction(field: &GridF32) -> f32 {
    let land = field.data.iter().filter(|&&v| v > 0.5).count();
    100.0 * land as f32 / field.data.len() as f32
}

fn run_sweep(target: usize) {
    let ss = SteinSteinParams::default();
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);

    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);

    // The production export config (relief-v3 stream-power, droplets off).
    let km_per_cell = DOMAIN_KM / target as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let mut base_cfg = FbmUpscaleConfig::c1_hd_production(target);
    base_cfg.sample_origin = [0.09375, 0.578125];
    base_cfg.sample_size = 1.0;
    base_cfg.erosion = None;
    let mut sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    sp.mfd_exponent = Some(2.0);
    sp.iterations = 2;
    base_cfg.stream_power = Some(sp);

    let cls = class_map(
        &state,
        &edifices,
        base_cfg.sample_origin,
        base_cfg.sample_size,
        km_per_cell,
        target,
        target,
    );
    let (nh, nsf, nv) = (
        cls.iter().filter(|&&c| c == HARD).count(),
        cls.iter().filter(|&&c| c == SOFT).count(),
        cls.iter().filter(|&&c| c == VOLC).count(),
    );
    eprintln!("\n================  C-3 SWEEP  {target}²  ================");
    eprintln!(
        "class coverage (HD cells): hard {} ({:.1}%) | rift-soft {} ({:.2}%) | volcaniclastic {} ({:.2}%)",
        nh,
        100.0 * nh as f32 / cls.len() as f32,
        nsf,
        100.0 * nsf as f32 / cls.len() as f32,
        nv,
        100.0 * nv as f32 / cls.len() as f32,
    );
    eprintln!(
        "{:>10} | {:>32} | {:>32} | {:>32} | {:>8} {:>6}",
        "soft ×",
        "HARD (relief/slope/steep‰/inc)",
        "SOFT (relief/slope/steep‰/inc)",
        "VOLC (relief/slope/steep‰/inc)",
        "pits",
        "land%"
    );

    for &mult in &[1.0f32, 3.0, 10.0, 30.0, 100.0] {
        let mut cfg = base_cfg.clone();
        cfg.lithology = LithologyConfig {
            enabled: mult > 1.0,
            soft_multiplier: mult,
            // Volcaniclastic held at √(soft) so it tracks the sweep as an intermediate
            // class (Stock & Montgomery: volcaniclastic sits between granite & mudstone).
            volcanic_multiplier: mult.sqrt().max(1.0),
            rift_age_threshold: 1.0,
        };
        let (res, _craters) = upscale_from_c1_with_progress(
            &state,
            &run.iso_config,
            &ss,
            &seed,
            &cfg,
            &edifices,
            &volc,
            &mut |_| {},
            &|| false,
        );
        let field = &res.heightmap;
        let sh = per_class_stats(field, &cls, HARD, km_per_cell);
        let sf = per_class_stats(field, &cls, SOFT, km_per_cell);
        let sv = per_class_stats(field, &cls, VOLC, km_per_cell);
        let fmt = |s: &ClassStats| {
            format!(
                "{:>6.0}/{:>4.1}/{:>5.0}/{:>5.0}",
                s.relief_m, s.slope_deg, s.steep_permille, s.incision_m
            )
        };
        eprintln!(
            "{:>10} | {:>32} | {:>32} | {:>32} | {:>8} {:>6.1}",
            if mult > 1.0 { format!("×{mult:.0}") } else { "OFF".into() },
            fmt(&sh),
            fmt(&sf),
            fmt(&sv),
            count_closed_depressions(field),
            land_fraction(field),
        );
        let _ = (sf.cells, sv.cells, sh.cells);
    }
    eprintln!(
        "(relief/slope/steep‰/incision, all in the class; hard columns ~constant = reference intact,\n only soft/volc move = the contrast. Softer K erodes DOWN → lower relief, wider open valleys.)"
    );
}

#[test]
#[ignore]
fn c3_sweep() {
    run_sweep(2048);
}

#[test]
#[ignore]
fn c3_sweep_8192() {
    run_sweep(8192);
}
