//! Shared bench tooling for the coastal-fringe campaign (ADR Findings 75–77).
//!
//! Cargo does not compile `tests/common/mod.rs` as its own test target, so this is the one
//! place a helper can live without becoming a fourth copy. It exists because Finding 75 already
//! had to name one duplication as debt (`crop_image` against `render_coast_crop`) and Finding 76
//! nearly added a second (the spectrum).
//!
//! Nothing here is production code: these are instruments, and the production path is reached
//! only through `production_hd_config` + `upscale_from_c1_with_progress`.

#![allow(dead_code)]

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::coast_metrics::coast_spurs;
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

pub const PSEED: u64 = 10_481_999_410_520_546_993;
pub const DOMAIN_KM: f32 = 400.0;
pub const SEA: f32 = 0.5;
pub const TARGET: usize = 8192;
pub const CELL_KM: f32 = DOMAIN_KM / TARGET as f32;
pub const CELL_KM2: f32 = CELL_KM * CELL_KM;
/// `RELIEF_V1_A_C_KM2 / cell_km2` at 8192² — the channel head the incision itself uses.
pub const A_C_CELLS: f32 = 0.1 / CELL_KM2;

/// One knob-set for a production build. `None` everywhere = the SHIPPED relief.
#[derive(Clone, Copy, Default)]
pub struct Knobs {
    /// `false` → `stream_power = None`: FBM and both closures stay ON, only the incision goes.
    /// That is NOT Finding 75's REFERENCE, which also drops the FBM.
    pub no_incision: bool,
    pub a_c_km2: Option<f32>,
    pub diffusion: Option<f32>,
    pub diffusion_substeps: Option<usize>,
    pub lateral_erosion: Option<f32>,
    pub talus_passes: Option<usize>,
    /// `BathymetryProfile::shelf_min_depth_m` — the clamp of ADR Finding 78. ⚠️ the production
    /// code takes `shelf_min_depth_m.max(1.0)` (`bathymetry.rs:171`), so a bench point below
    /// 1 m is CLAMPED TO 1 and must be declared, not silently swept.
    pub shelf_min_depth_m: Option<f32>,
    /// Drop the bathymetry re-map entirely — this is the field the INCISION actually produced,
    /// before the sub-sea overwrite. `apply_bathymetry_profile` runs at
    /// `production_upscale.rs:506`, i.e. AFTER `incise_lithology` at :451.
    pub bathymetry_off: bool,
    /// `a_c_slope_law` at the TARGET-GRID calibration of VALIDATION_NOTE section 1 as amended
    /// (`S_ref = 0.1128`). The production gate stays OFF; this is a bench variant.
    pub a_c_law: bool,
}

impl Knobs {
    pub fn shipped() -> Self {
        Self::default()
    }
    pub fn no_incision() -> Self {
        Self { no_incision: true, ..Self::default() }
    }
}

/// Build the 8192² field at the production seed and config, with `k` applied to the
/// stream-power stage only. Everything upstream (tectonics, upscale, FBM, closures) is the
/// production path, so a sweep point differs from SHIPPED in exactly one term.
pub fn build_field(k: Knobs) -> GridF32 {
    let ss = SteinSteinParams::default();
    let run_cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run_cfg, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: TARGET,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: 0.04,
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: true,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: true,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    if k.bathymetry_off {
        cfg.bathymetry = None;
    } else if let Some(d) = k.shelf_min_depth_m {
        if let Some(b) = cfg.bathymetry.as_mut() {
            b.shelf_min_depth_m = d;
        }
    }
    if k.no_incision {
        cfg.stream_power = None;
    } else if let Some(sp) = cfg.stream_power.as_mut() {
        if let Some(a) = k.a_c_km2 {
            sp.min_area_cells = a / CELL_KM2;
        }
        if let Some(d) = k.diffusion {
            sp.diffusion = d;
        }
        if let Some(s) = k.diffusion_substeps {
            sp.diffusion_substeps = s;
        }
        if let Some(l) = k.lateral_erosion {
            sp.lateral_erosion = l;
        }
        if let Some(t) = k.talus_passes {
            sp.talus_passes = t;
        }
        if k.a_c_law {
            sp.a_c_slope_law = Some(ymir_core::erosion::stream_power::ChannelHeadLaw {
                s_ref: 0.1128,
                s_min: ymir_core::erosion::stream_power::CHANNEL_HEAD_S_MIN,
            });
        }
    }
    upscale_from_c1_with_progress(
        &state,
        &run_cfg.iso_config,
        &ss,
        &seed,
        &cfg,
        &edifices,
        &volc,
        Some(&kin),
        &mut |_| {},
        &|| false,
    )
    .0
    .heightmap
}

pub fn pct(v: &[f32], p: f64) -> f32 {
    if v.is_empty() { f32::NAN } else { v[(((v.len() - 1) as f64) * p) as usize] }
}

pub fn sorted(mut v: Vec<f32>) -> Vec<f32> {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v
}

/// The spectrum of Finding 76 block C: gaps between consecutive spur roots along each
/// continuous coast segment longer than 20 km, against a WHITE baseline — the same number of
/// spurs with positions drawn uniformly over the same length, averaged over `DRAWS` draws.
/// A real coast has no wavelength; a threshold or instability artefact does.
///
/// Returns `(segments, gaps, median_gap_cells, peak_cells, peak_over_white)`. A population under
/// 32 gaps returns `None` — rule 10: too few to read, reported as such rather than as a flat
/// spectrum.
pub const DRAWS: usize = 64;

pub fn spectrum(polys: &[Vec<(f32, f32)>]) -> Option<(usize, usize, f32, usize, f64)> {
    let (sp, _) = coast_spurs(polys, CELL_KM, 2.0 * CELL_KM, CELL_KM);
    let mut by_poly: std::collections::HashMap<usize, Vec<f32>> = Default::default();
    for s in &sp {
        by_poly.entry(s.poly).or_default().push(s.root_km);
    }
    let mut gaps: Vec<f32> = Vec::new();
    let mut used = 0usize;
    let mut keys: Vec<usize> = by_poly.keys().copied().collect();
    keys.sort_unstable(); // determinism: a HashMap walk is not an order
    for pi in keys {
        let mut roots = by_poly.remove(&pi).unwrap();
        let pl = &polys[pi];
        let mut len = 0.0f32;
        for i in 1..pl.len() {
            len += ((pl[i].0 - pl[i - 1].0).powi(2) + (pl[i].1 - pl[i - 1].1).powi(2)).sqrt();
        }
        if len * CELL_KM < 20.0 || roots.len() < 8 {
            continue;
        }
        used += 1;
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        for w in roots.windows(2) {
            gaps.push((w[1] - w[0]) / CELL_KM);
        }
    }
    if gaps.len() < 32 {
        return None;
    }
    let mut hist = [0usize; 33];
    for &g in &gaps {
        hist[(g.round() as usize).min(32)] += 1;
    }
    let mut base = [0.0f64; 33];
    let mut st = 0x2545_F491_4F6C_DD1Du64;
    let mut rnd = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let total: f32 = gaps.iter().sum();
    for _ in 0..DRAWS {
        let mut u: Vec<f32> = (0..gaps.len() + 1).map(|_| (rnd() as f32) * total).collect();
        u.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        for w in u.windows(2) {
            base[((w[1] - w[0]).round() as usize).min(32)] += 1.0;
        }
    }
    for b in base.iter_mut() {
        *b /= DRAWS as f64;
    }
    let (mut peak, mut ratio) = (0usize, 0.0f64);
    for b in 1..=32 {
        let r = hist[b] as f64 / base[b].max(0.5);
        if r > ratio {
            ratio = r;
            peak = b;
        }
    }
    let g = sorted(gaps.clone());
    Some((used, gaps.len(), g[g.len() / 2], peak, ratio))
}
