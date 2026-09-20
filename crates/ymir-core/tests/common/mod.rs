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

// ═══════════════════════════════════════════════════════════════════════════════════════════════
// ADR 0001 Finding 88 — RULE 12, as code: a bench that TRACES declares its field and proves by
// hash that it is production's, BEFORE it reads a single reason off it.
//
// Third instance in eight findings, which is why it stops being a habit and becomes a function:
//   * Finding 80-B2 read the CLAMP where it meant the incision;
//   * Finding 85 paired basins BY ID where ids are scan-order (Finding 81's own trap);
//   * Finding 87-D traced on a plain `compute_flow` where production routes on
//     `C1DrainageResult::flow`, which carries `FlatPerturbation` — and reported the WRONG REASON
//     for fifteen of nineteen lakes until it was re-run on the right field.
//
// The cost of the error is always the same shape: a number that looks like a measurement and is
// an artefact of the instrument. `declared_flow` makes the declaration mechanical.
// ═══════════════════════════════════════════════════════════════════════════════════════════════

/// A stable 64-bit digest of every byte of a flow field that a trace can depend on: the D8
/// direction of each cell and its accumulation. FNV-1a, so the value is reproducible across runs
/// and machines (`DefaultHasher` is explicitly not).
pub fn flow_field_hash(flow: &ymir_core::terrain::flow::FlowResult) -> u64 {
    let mut hv: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |b: u8| {
        hv ^= b as u64;
        hv = hv.wrapping_mul(0x0000_0100_0000_01b3);
    };
    for &d in &flow.direction {
        eat(d);
    }
    for &a in &flow.accumulation.data {
        for b in a.to_bits().to_le_bytes() {
            eat(b);
        }
    }
    hv
}

/// Rule 12. Returns the field to trace on, having PRINTED its hash and that of production's, and
/// PANICKED if they differ. `label` names the candidate in the failure message.
///
/// Use it at the top of every bench block that walks `flow.direction`. If the two fields are the
/// same object the hashes trivially agree and the line is still worth its cost: it puts the field's
/// identity in the transcript beside the numbers read off it, which is what Finding 87-D lacked.
pub fn declared_flow<'a>(
    label: &str,
    candidate: &'a ymir_core::terrain::flow::FlowResult,
    production: &ymir_core::terrain::flow::FlowResult,
) -> &'a ymir_core::terrain::flow::FlowResult {
    let (hc, hp) = (flow_field_hash(candidate), flow_field_hash(production));
    println!("  [rule 12] field `{label}` 0x{hc:016x} | production 0x{hp:016x}");
    assert_eq!(
        hc, hp,
        "RULE 12: the bench would trace on `{label}` (0x{hc:016x}) while production routes on          0x{hp:016x}. Finding 87-D paid three runs for exactly this."
    );
    candidate
}

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
    /// ADR Finding 80/83 -- OVERRIDE the base-level epsilon, in metres above sea level.
    /// `None` = the SHIPPED `RELIEF_V3_BASE_LEVEL_M` (0.5 m). Since Finding 83 the bound
    /// is production, so `None` no longer means "no bound" -- use `base_level_off`.
    pub base_level_m: Option<f32>,
    /// ADR Finding 105 -- drop the droplet pass (`cfg.erosion = None`), to isolate the field the
    /// incision itself produces from the one the pipeline delivers.
    pub erosion_off: bool,
    /// ADR Finding 105 -- `FbmUpscaleConfig::slope_floor_factor`: the TWO-PASS closure, the
    /// shipped form. Mutually exclusive with `slope_floor_uk` in practice (that one is the
    /// bench's manual form, with `k_s` supplied from outside).
    pub slope_floor_factor: Option<f32>,
    /// ADR Finding 97 B2 -- `StreamPowerConfig::slope_floor_uk`: the local equilibrium-slope
    /// floor, `h_r_eff = min(h_r + S_eq(A)·dx, h_o)`. `None` = shipped.
    pub slope_floor_uk: Option<f32>,
    /// ADR Finding 96 A1 -- `StreamPowerConfig::depression_floor`: the incision may not cut a
    /// cell below the SPILL LEVEL of the depression it sits in. `false` = shipped.
    pub depression_floor: bool,
    /// ADR Finding 83 -- the A/B control: drop the shipped bound and rebuild the
    /// pre-Finding-83 world (the one with the coastal fringe). This is what the benches
    /// must use for the "DELIVERED (pre-83)" column; `Knobs::shipped()` is now BOUNDED.
    pub base_level_off: bool,
    /// ADR Finding 83-B2 -- `BaseLevelFloor::free_above_km2`, the estuary gate, in km2 of
    /// the SIMULATED grid. **PROXY**, bench-only; production is `None`.
    pub a_est_km2: Option<f32>,
    /// ADR Finding 90-C -- `StreamPowerConfig::k`. Used ONLY with [`Knobs::integrating`], which
    /// holds `k_time = k * dt * iterations` at the shipped 9000 while raising the pass count, so
    /// the sweep changes the DISTRIBUTION of the erosion work and not its quantity. That is what
    /// makes the Finding 44 item-2 experiment well posed (Finding 44: "`K` and the duration are
    /// not separately observable" -- holding `k_time` is the only way to move one dial alone).
    pub k: Option<f32>,
    /// ADR Finding 88-D1 -- `StreamPowerConfig::iterations`. Production ships **2**; the point of
    /// the knob is to read the field AFTER ONE pass, because at Courant 1353 (Finding 61) a cell
    /// relaxes onto its receiver in a single step and "how much of the cut is the first pass" is
    /// the only way to tell a one-shot landing from a rhythm.
    pub iterations: Option<usize>,
}

impl Knobs {
    /// The SHIPPED relief. Since ADR Finding 83 that includes the base-level bound.
    pub fn shipped() -> Self {
        Self::default()
    }
    /// The pre-Finding-83 world: every shipped knob, bound OFF. This is the field every
    /// "delivered" column of Findings 76-82 was measured on.
    pub fn pre83() -> Self {
        Self { base_level_off: true, ..Self::default() }
    }
    pub fn no_incision() -> Self {
        Self { no_incision: true, ..Self::default() }
    }
    /// ADR Finding 88-D1 -- the shipped relief after exactly `it` incision passes.
    pub fn passes(it: usize) -> Self {
        Self { iterations: Some(it), ..Self::default() }
    }
    /// ADR Finding 90-C -- `n` passes at the SHIPPED `k_time` (9000): `k = 9000 / n` with
    /// `dt = 1.0`. The integrated erosion budget is identical to the delivered field's; only the
    /// Courant number falls, by `n / 2`. At `n = 7399` the configuration sits at Courant 1.0,
    /// which is Finding 44's item 2 and the whole point.
    pub fn integrating(n: usize) -> Self {
        Self {
            iterations: Some(n),
            k: Some(ymir_core::erosion::stream_power::SHIPPED_K_TIME / n as f32),
            ..Self::default()
        }
    }
}

/// ADR 0001 Finding 88-D3 — the C-3 erodibility multiplier field and the C-2 edifices, built by
/// the SAME calls `build_field` makes, so block D3 reads the field production incised and does not
/// rebuild one. `None` for the K field means every contributing closure was off (it never is here).
///
/// This repeats the 300-step coarse run, which is the price of not reconstructing: the ADR records
/// six diagnoses misled by a rebuilt terrain.
pub fn production_k_and_edifices()
-> (Option<Vec<f32>>, Vec<ymir_core::tectonics_c1::closures::volcanism::Edifice>) {
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
    let cfg = production_hd_config(&ProductionHdOpts {
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
    let k = ymir_core::tectonics_c1::production_upscale::production_k_field(
        &state,
        Some(&kin),
        &cfg,
        &edifices,
        &volc,
        TARGET,
        TARGET,
    );
    (k, edifices)
}

/// Build the 8192² field at the production seed and config, with `k` applied to the
/// stream-power stage only. Everything upstream (tectonics, upscale, FBM, closures) is the
/// production path, so a sweep point differs from SHIPPED in exactly one term.
pub fn build_field(k: Knobs) -> GridF32 {
    build_field_with_floor(k, None)
}

/// [`build_field`] on an ARBITRARY seed. ADR Finding 102: every measurement of this chantier
/// before that finding was made on `PSEED` alone, and Finding 92 asked for a second seed without
/// getting one. The seed drives both the tectonic init and [`WorldSeed`], as production does.
pub fn build_field_seed(k: Knobs, seed: u64) -> GridF32 {
    build_field_inner(k, None, seed)
}

/// [`build_field`] with ADR Finding 96 B3's external incision floor (norm units, one entry per
/// target cell). `None` is byte-identical. Separate entry point because [`Knobs`] is `Copy` and
/// a grid cannot live in it -- the same constraint ADR Finding 88 hit on `LakeCellInfo`.
pub fn build_field_with_floor(k: Knobs, floor: Option<std::sync::Arc<Vec<f32>>>) -> GridF32 {
    build_field_inner(k, floor, PSEED)
}

fn build_field_inner(k: Knobs, floor: Option<std::sync::Arc<Vec<f32>>>, pseed: u64) -> GridF32 {
    let ss = SteinSteinParams::default();
    let run_cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, pseed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run_cfg, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(pseed);
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
        if let Some(it) = k.iterations {
            sp.iterations = it;
        }
        if let Some(kk) = k.k {
            sp.k = kk;
        }
        if let Some(l) = k.lateral_erosion {
            sp.lateral_erosion = l;
        }
        if let Some(t) = k.talus_passes {
            sp.talus_passes = t;
        }
        sp.depression_floor = k.depression_floor; // ADR Finding 96 A1
        sp.slope_floor_uk = k.slope_floor_uk; // ADR Finding 97 B2
    }
    cfg.slope_floor_factor = k.slope_floor_factor; // ADR Finding 105
    if k.erosion_off {
        cfg.erosion = None; // ADR Finding 105
    }
    if cfg.stream_power.is_some() {
        let sp = cfg.stream_power.as_mut().unwrap();
        // ADR Finding 83 -- `relief_v3` now ARRIVES bounded, so the knobs edit what is
        // there rather than installing it. `base_level_off` wins over the others.
        if k.base_level_off {
            sp.base_level_floor = None;
        } else {
            let mut b = sp.base_level_floor.clone().unwrap_or(
                ymir_core::erosion::stream_power::BaseLevelFloor {
                    epsilon_m: ymir_core::erosion::stream_power::RELIEF_V3_BASE_LEVEL_M,
                    free_above_km2: None,
                },
            );
            if let Some(e) = k.base_level_m {
                b.epsilon_m = e;
            }
            b.free_above_km2 = k.a_est_km2;
            sp.base_level_floor = Some(b);
        }
        if k.a_c_law {
            sp.a_c_slope_law = Some(ymir_core::erosion::stream_power::ChannelHeadLaw {
                s_ref: 0.1128,
                s_min: ymir_core::erosion::stream_power::CHANNEL_HEAD_S_MIN,
            });
        }
    }
    cfg.incision_floor = floor; // ADR Finding 96 B3
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

// ─────────────────────────────────────────────────────────────────────────────
// ADR Finding 83 — the RICHARDSON instrument, shared.
//
// It lives here and not in a bench because two benches need it (B1's baseline and B2's sweep)
// and the ADR has already had to name one instrument duplication as debt (Finding 75).
// ─────────────────────────────────────────────────────────────────────────────

/// Divider start offsets averaged over. A single start point makes the step count depend on where
/// the walk began, which is a property of the instrument and not of the coast.
pub const RICHARDSON_OFFSETS: usize = 8;
/// Log-spaced rulers over the declared range (2 cells to `max_km`).
pub const RICHARDSON_RULERS: usize = 13;
/// A polygon must carry the biggest ruler several times over or its contribution is a rounding
/// error that moves with the ruler. Declared, and FIXED across variants.
pub const RICHARDSON_MIN_POLY_KM: f32 = 20.0;

/// One divider ("compass") walk around a closed polyline: how many rulers of length `eps` (in the
/// polygon's own units, i.e. CELLS) fit around it, starting at vertex `start`.
///
/// CHORD distance, not arc — the dividers' two points are `eps` apart in a straight line, which is
/// Richardson's own instrument. The final partial step counts as a fraction, so the result is
/// continuous in `eps`; without that the curve is a staircase of integers at the large rulers and
/// the regression reads the staircase instead of the coast.
pub fn divider_steps(poly: &[(f32, f32)], eps: f32, start: usize) -> f64 {
    let n = poly.len();
    if n < 4 || eps <= 0.0 {
        return 0.0;
    }
    let at = |i: usize| poly[i % n];
    let p0 = at(start);
    let (mut px, mut py) = p0;
    let mut k = start;
    let mut cur = (px, py);
    let mut steps = 0.0f64;
    let mut advanced = 0usize;
    let e2 = (eps * eps) as f64;
    loop {
        let mut hit: Option<(f32, f32)> = None;
        while advanced < n {
            let b = at(k + 1);
            let (fx, fy) = ((cur.0 - px) as f64, (cur.1 - py) as f64);
            let (dx, dy) = ((b.0 - cur.0) as f64, (b.1 - cur.1) as f64);
            let a = dx * dx + dy * dy;
            if a > 1e-18 {
                let bq = fx * dx + fy * dy;
                let c = fx * fx + fy * fy - e2;
                let disc = bq * bq - a * c;
                if disc >= 0.0 {
                    let t = (-bq + disc.sqrt()) / a;
                    if t > 0.0 && t <= 1.0 {
                        hit = Some((
                            cur.0 + t as f32 * (b.0 - cur.0),
                            cur.1 + t as f32 * (b.1 - cur.1),
                        ));
                        break;
                    }
                }
            }
            cur = b;
            k += 1;
            advanced += 1;
        }
        match hit {
            Some(q) => {
                steps += 1.0;
                px = q.0;
                py = q.1;
                cur = q;
            }
            None => {
                let d = (((p0.0 - px).powi(2) + (p0.1 - py).powi(2)).sqrt() / eps) as f64;
                return steps + d.min(1.0);
            }
        }
    }
}

/// Perimeter of a polyline in CELLS at the finest ruler (its own vertex spacing).
pub fn perimeter_cells(poly: &[(f32, f32)]) -> f32 {
    let mut s = 0.0f32;
    for i in 1..poly.len() {
        s += ((poly[i].0 - poly[i - 1].0).powi(2) + (poly[i].1 - poly[i - 1].1).powi(2)).sqrt();
    }
    s
}

/// The Richardson curve: `(ruler_km, length_km)` for `RICHARDSON_RULERS` log-spaced rulers from
/// 2 cells to `max_km`, averaged over `RICHARDSON_OFFSETS` start points, summed over `polys`.
pub fn richardson(polys: &[&Vec<(f32, f32)>], cell_km: f32, max_km: f32) -> Vec<(f64, f64)> {
    let lo = 2.0f64;
    let hi = (max_km / cell_km) as f64;
    (0..RICHARDSON_RULERS)
        .map(|i| {
            let eps = lo * (hi / lo).powf(i as f64 / (RICHARDSON_RULERS - 1) as f64);
            let mut total = 0.0f64;
            for pl in polys {
                let n = pl.len();
                let mut acc = 0.0f64;
                for o in 0..RICHARDSON_OFFSETS {
                    acc += divider_steps(pl, eps as f32, o * n / RICHARDSON_OFFSETS);
                }
                total += acc / RICHARDSON_OFFSETS as f64;
            }
            (eps * cell_km as f64, total * eps * cell_km as f64)
        })
        .collect()
}

/// Least squares on `log10 L = c + s * log10(eps)`; returns `(D, R2, worst_residual_pct)` with
/// `D = 1 - s`.
///
/// ⚠️ On a nearly FLAT curve (a smooth coast) the total sum of squares is tiny, so `R2` is small
/// for reasons that have nothing to do with linearity. Judge linearity by the residual per cent,
/// and read `R2` only when D is far from 1.
pub fn richardson_fit(pts: &[(f64, f64)]) -> (f64, f64, f64) {
    if pts.len() < 3 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    let x: Vec<f64> = pts.iter().map(|p| p.0.log10()).collect();
    let y: Vec<f64> = pts.iter().map(|p| p.1.max(1e-9).log10()).collect();
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let sxy: f64 = x.iter().zip(&y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    let s = sxy / sxx;
    let c = my - s * mx;
    let sst: f64 = y.iter().map(|b| (b - my) * (b - my)).sum();
    let mut sse = 0.0;
    let mut worst = 0.0f64;
    for (a, b) in x.iter().zip(&y) {
        let r = b - (c + s * a);
        sse += r * r;
        worst = worst.max((10f64.powf(r.abs()) - 1.0) * 100.0);
    }
    (1.0 - s, if sst > 0.0 { 1.0 - sse / sst } else { f64::NAN }, worst)
}

/// The DELIVERED land mask: encode with the production encoder, decode, threshold at 0 m. This is
/// the raster Living Landz reads, and the criterion is written on it.
pub fn land_u16(f: &GridF32, ss: &SteinSteinParams) -> Vec<bool> {
    let hl = ymir_core::export::height::metric_height_u16(f, ss);
    hl.codes
        .iter()
        .map(|&c| hl.min_m + (c as f32 / 65535.0) * (hl.max_m - hl.min_m) > 0.0)
        .collect()
}

/// Land/sea as a field the contour tracer can read; tracing at 0.5 gives the mask staircase.
pub fn to_mask(land: &[bool], w: usize) -> GridF32 {
    let mut m = GridF32::new(w, w, 0.0);
    for k in 0..w * w {
        if land[k] {
            m.data[k] = 1.0;
        }
    }
    m
}

/// Binary majority downsample — the author's declared consumer method (Finding 80-B0): threshold
/// at 8192² then resample the MASK, not the heights. Ties go to SEA, declared.
pub fn majority(land: &[bool], w: usize, f: usize) -> (Vec<bool>, usize) {
    let nw = w / f;
    let mut out = vec![false; nw * nw];
    let half = f * f / 2;
    for y in 0..nw {
        for x in 0..nw {
            let mut c = 0usize;
            for dy in 0..f {
                for dx in 0..f {
                    if land[(y * f + dy) * w + x * f + dx] {
                        c += 1;
                    }
                }
            }
            out[y * nw + x] = c > half;
        }
    }
    (out, nw)
}

/// One Richardson line: traces the mask, keeps the fixed polygon population, and reports D over
/// the whole range, under 1 km, over 1 km, and on the LARGEST polygon alone (a paired population
/// across variants — the whole-population D also moves when the number of islands does).
pub fn richardson_line(tag: &str, land: &[bool], w: usize, cell_km: f32, max_km: f32) {
    use ymir_core::terrain::contour::marching_squares;
    let polys = marching_squares(&to_mask(land, w), 0.5);
    let keep: Vec<&Vec<(f32, f32)>> =
        polys.iter().filter(|p| perimeter_cells(p) * cell_km >= RICHARDSON_MIN_POLY_KM).collect();
    let biggest = polys
        .iter()
        .max_by(|a, b| {
            perimeter_cells(a).partial_cmp(&perimeter_cells(b)).unwrap_or(std::cmp::Ordering::Equal)
        })
        .cloned()
        .unwrap_or_default();
    let all = richardson(&keep, cell_km, max_km);
    let one = richardson(&[&biggest], cell_km, max_km);
    let sub: Vec<(f64, f64)> = all.iter().copied().filter(|p| p.0 <= 1.001).collect();
    let sup: Vec<(f64, f64)> = all.iter().copied().filter(|p| p.0 >= 0.999).collect();
    let (d_all, r2_all, res_all) = richardson_fit(&all);
    let (d_sub, _, _) = richardson_fit(&sub);
    let (d_sup, _, _) = richardson_fit(&sup);
    let (d_one, _, res_one) = richardson_fit(&one);
    eprintln!(
        "   {tag:<26} polys {:>5} (of {:>6}) | **D {d_all:.3}** R2 {r2_all:.4} worst resid \
         {res_all:.2} % | D(<1 km) {d_sub:.3} | D(>1 km) {d_sup:.3} | D(largest, PAIRED) \
         {d_one:.3} resid {res_one:.2} %",
        keep.len(),
        polys.len()
    );
    eprintln!(
        "      L(km) by ruler: {}",
        all.iter().map(|(e, l)| format!("{e:.3}km:{l:.0}")).collect::<Vec<_>>().join("  ")
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// ADR Finding 84 block B1 — the candidate "organic" quantities, and the SHAPES that calibrate
// them. Shared because B1 needs them and Finding 83's Richardson bench already needed two of
// them; a third copy is what the ADR has twice named as debt.
// ─────────────────────────────────────────────────────────────────────────────

/// Fill a closed polygon into a `w x w` binary mask by even-odd scanline. Edges are bucketed by
/// row, so the cost is the perimeter and not `rows x edges`.
pub fn fill_polygon(pts: &[(f64, f64)], w: usize) -> Vec<bool> {
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); w];
    let n = pts.len();
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        let (y0, y1) = if a.1 < b.1 { (a.1, b.1) } else { (b.1, a.1) };
        let lo = (y0.floor().max(0.0)) as usize;
        let hi = (y1.ceil().min((w - 1) as f64)) as usize;
        for row in buckets.iter_mut().take(hi + 1).skip(lo) {
            row.push(i);
        }
    }
    let mut out = vec![false; w * w];
    let mut xs: Vec<f64> = Vec::new();
    for y in 0..w {
        let yc = y as f64 + 0.5;
        xs.clear();
        for &i in &buckets[y] {
            let (a, b) = (pts[i], pts[(i + 1) % n]);
            if (a.1 <= yc) != (b.1 <= yc) {
                xs.push(a.0 + (yc - a.1) / (b.1 - a.1) * (b.0 - a.0));
            }
        }
        xs.sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
        for pair in xs.chunks_exact(2) {
            let (x0, x1) = (pair[0].max(0.0) as usize, (pair[1].min((w - 1) as f64)) as usize);
            for x in x0..=x1 {
                out[y * w + x] = true;
            }
        }
    }
    out
}

/// A Koch snowflake, D = log 4 / log 3 = 1.2619, `generations` rounds from a triangle of `side` px.
pub fn koch_polygon(w: usize, side: f64, generations: usize) -> Vec<(f64, f64)> {
    let c = w as f64 / 2.0;
    let r = side / 3f64.sqrt();
    let mut pts: Vec<(f64, f64)> = (0..3)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / 3.0 - std::f64::consts::FRAC_PI_2;
            (c + r * a.cos(), c + r * a.sin())
        })
        .collect();
    for _ in 0..generations {
        let n = pts.len();
        let mut next = Vec::with_capacity(n * 4);
        for i in 0..n {
            let (a, b) = (pts[i], pts[(i + 1) % n]);
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let p = (a.0 + dx / 3.0, a.1 + dy / 3.0);
            let q = (a.0 + 2.0 * dx / 3.0, a.1 + 2.0 * dy / 3.0);
            let (ux, uy) = (q.0 - p.0, q.1 - p.1);
            let (co, si) = (0.5f64, -(3f64.sqrt() / 2.0));
            let m = (p.0 + ux * co - uy * si, p.1 + ux * si + uy * co);
            next.push(a);
            next.push(p);
            next.push(m);
            next.push(q);
        }
        pts = next;
    }
    pts
}

/// A circle of radius `r` px, sampled finely. True D = 1.000, zero curvature inversions.
pub fn circle_polygon(w: usize, r: f64) -> Vec<(f64, f64)> {
    let c = w as f64 / 2.0;
    let n = 8192usize;
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            (c + r * a.cos(), c + r * a.sin())
        })
        .collect()
}

/// A square rotated 30°, side `s` px. True D = 1.000 with FOUR corners — the control that says
/// whether a curvature instrument fires on corners that are content rather than noise.
pub fn rotated_square_polygon(w: usize, s: f64, deg: f64) -> Vec<(f64, f64)> {
    let c = w as f64 / 2.0;
    let t = deg.to_radians();
    let corners = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)];
    let mut out = Vec::new();
    for i in 0..4 {
        let (ax, ay) = corners[i];
        let (bx, by) = corners[(i + 1) % 4];
        for j in 0..2048 {
            let u = j as f64 / 2048.0;
            let (x, y) = ((ax + u * (bx - ax)) * s, (ay + u * (by - ay)) * s);
            out.push((c + x * t.cos() - y * t.sin(), c + x * t.sin() + y * t.cos()));
        }
    }
    out
}

/// **The PERIODIC negative control**: a smooth circle carrying regular radial teeth of
/// `amp` px at a spacing of `lambda_px` along the circumference. Built to Finding 76's measured
/// fur geometry (amplitude 2–5 cells, λ ≈ 14 cells). A quantity that claims to detect "a
/// generator" must fire on this.
pub fn periodic_polygon(w: usize, r: f64, amp: f64, lambda_px: f64) -> Vec<(f64, f64)> {
    let c = w as f64 / 2.0;
    let teeth = (std::f64::consts::TAU * r / lambda_px).round().max(4.0);
    let n = 65536usize;
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            let saw = if (teeth * a).cos() >= 0.0 { amp } else { 0.0 };
            let rr = r + saw;
            (c + rr * a.cos(), c + rr * a.sin())
        })
        .collect()
}

/// A radial fractional-Brownian boundary, `r(θ) = R₀ (1 + a · Σ k^-(H+1/2) cos(kθ + φ_k))`,
/// with `H = 2 − D`. Deterministic phases from a fixed xorshift.
///
/// ⚠️ **KEPT AS A MEASURED NEGATIVE, not used as the isotropic control.** Built at
/// `r = 3000` px, `D = 1.25`, amplitude 5 %, it MEASURES **D = 1.008** over 100 m – 10 km —
/// and that is correct arithmetic, not a bug: an fBm curve's excursion over an arc `Δs`
/// scales as `Δs^H`, so over 2 px of a 18 850 px circumference it is `150 · (2/18850)^0.75
/// ≈ 0.16` px. **Sub-pixel.** A D = 1.25 fBm lobe of moderate amplitude is SMOOTH at the
/// cell scale, so it cannot calibrate a cell-scale instrument; reaching 2 px of roughness
/// at 2 px of arc would need an amplitude of ~64 % of the radius, i.e. a self-intersecting
/// blob. Use [`random_koch_polygon`] instead, which is rough down to `side / 3^g`.
pub fn fbm_polygon(w: usize, r: f64, d_target: f64, amp_frac: f64, kmax: usize) -> Vec<(f64, f64)> {
    let c = w as f64 / 2.0;
    let hurst = 2.0 - d_target;
    let mut st = 0x9E37_79B9_7F4A_7C15u64;
    let mut rnd = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let phase: Vec<f64> = (0..=kmax).map(|_| rnd() * std::f64::consts::TAU).collect();
    let coef: Vec<f64> =
        (0..=kmax).map(|k| if k == 0 { 0.0 } else { (k as f64).powf(-(hurst + 0.5)) }).collect();
    let norm = coef.iter().map(|x| x * x).sum::<f64>().sqrt() * std::f64::consts::SQRT_2.recip();
    let n = 65536usize;
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            let mut s = 0.0;
            for k in 1..=kmax {
                s += coef[k] * (k as f64 * a + phase[k]).cos();
            }
            let rr = r * (1.0 + amp_frac * s / norm.max(1e-12));
            (c + rr * a.cos(), c + rr * a.sin())
        })
        .collect()
}

/// Windows of `win_km` of arc length along every polygon at least that long, as index ranges.
fn arc_windows(poly: &[(f32, f32)], cell_km: f32, win_km: f32) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let (mut i0, mut acc) = (0usize, 0.0f32);
    for i in 1..poly.len() {
        acc += ((poly[i].0 - poly[i - 1].0).powi(2) + (poly[i].1 - poly[i - 1].1).powi(2)).sqrt()
            * cell_km;
        if acc >= win_km {
            out.push((i0, i));
            i0 = i;
            acc = 0.0;
        }
    }
    out
}

/// **Normal-orientation entropy.** The contour normal's AXIAL angle (mod π) binned into 18 bins
/// per 20 km window; Shannon entropy divided by `ln 18`, so 1.0 is a perfectly isotropic window
/// and 0.0 a window whose normal never turns. Returns the MEDIAN over windows and the count.
///
/// Axial (mod π), not directional (mod 2π): a fringe tooth's two flanks have opposite normals and
/// the same axis, and it is the axis that "perpendicular to the coast everywhere" is about.
pub fn normal_entropy(polys: &[Vec<(f32, f32)>], cell_km: f32, win_km: f32) -> (f64, usize) {
    const BINS: usize = 18;
    let mut vals: Vec<f64> = Vec::new();
    for pl in polys {
        if pl.len() < 32 {
            continue;
        }
        for (a, b) in arc_windows(pl, cell_km, win_km) {
            let mut hist = [0.0f64; BINS];
            let mut tot = 0.0f64;
            let mut i = a + 4;
            while i + 4 < b {
                let (ax, ay) = pl[i - 4];
                let (bx, by) = pl[i + 4];
                let (tx, ty) = (bx - ax, by - ay);
                let len = (tx * tx + ty * ty).sqrt();
                if len > 1e-6 {
                    let ang = (-ty as f64).atan2(tx as f64).rem_euclid(std::f64::consts::PI);
                    let bin = ((ang / std::f64::consts::PI) * BINS as f64) as usize;
                    hist[bin.min(BINS - 1)] += 1.0;
                    tot += 1.0;
                }
                i += 1;
            }
            if tot < BINS as f64 * 2.0 {
                continue; // rule 10: too few samples in the window to read a histogram
            }
            let mut e = 0.0;
            for c in hist {
                if c > 0.0 {
                    let p = c / tot;
                    e -= p * p.ln();
                }
            }
            vals.push(e / (BINS as f64).ln());
        }
    }
    let n = vals.len();
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    (if n == 0 { f64::NAN } else { vals[n / 2] }, n)
}

/// **Curvature radii and inversions.** Signed curvature from the circle through
/// `p[i-stencil], p[i], p[i+stencil]`. Returns
/// `(r_p10, r_p50, r_p90 in CELLS, share of sign inversions closer than `near_cells` of arc,
/// inversion count)`.
///
/// ⚠️ **This instrument has a RESOLUTION CEILING and it is not optional.** A circle of radius `R`
/// has a sagitta of `(2·stencil)²/(8R)` cells over the stencil; when that is below the contour's
/// own sub-cell wobble the fitted radius and the curvature SIGN are both raster noise. At
/// `stencil = 2` a rasterised 3 000-cell circle — curvature of one constant sign, by
/// construction — measures **55.8 %** of its inversions closer than 3 cells and a median radius
/// of **3.3 cells**. Both figures are noise. Read this only against the circle's reading at the
/// same stencil, and treat the circle's value as the FLOOR.
pub fn curvature_stats(
    polys: &[Vec<(f32, f32)>],
    near_cells: f32,
    stencil: usize,
) -> (f32, f32, f32, f64, usize) {
    let s_st = stencil.max(1);
    let mut radii: Vec<f32> = Vec::new();
    let (mut gaps_near, mut gaps_all) = (0usize, 0usize);
    for pl in polys {
        if pl.len() < 4 * s_st + 8 {
            continue;
        }
        let mut last_sign = 0i8;
        let mut arc_since = 0.0f32;
        for i in s_st..pl.len() - s_st {
            let (a, b, c) = (pl[i - s_st], pl[i], pl[i + s_st]);
            let cross = (b.0 - a.0) * (c.1 - b.1) - (b.1 - a.1) * (c.0 - b.0);
            let (ab, bc, ca) = (
                ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt(),
                ((c.0 - b.0).powi(2) + (c.1 - b.1).powi(2)).sqrt(),
                ((a.0 - c.0).powi(2) + (a.1 - c.1).powi(2)).sqrt(),
            );
            arc_since += ab;
            if cross.abs() < 1e-9 {
                continue; // collinear: infinite radius, no sign
            }
            radii.push((ab * bc * ca / (2.0 * cross.abs())).min(1.0e6));
            let sign = if cross > 0.0 { 1i8 } else { -1i8 };
            if last_sign != 0 && sign != last_sign {
                gaps_all += 1;
                if arc_since < near_cells {
                    gaps_near += 1;
                }
                arc_since = 0.0;
            }
            if sign != last_sign {
                last_sign = sign;
            }
        }
    }
    let r = sorted(radii);
    (
        pct(&r, 0.10),
        pct(&r, 0.50),
        pct(&r, 0.90),
        if gaps_all == 0 { 0.0 } else { 100.0 * gaps_near as f64 / gaps_all as f64 },
        gaps_all,
    )
}

/// **Excursion-size distribution**, from the same `coast_spurs` walk the counts use: lengths in
/// CELLS, as `(p10, p50, p90, CV, count)`. A periodic fringe is a delta; a real coast is broad.
pub fn excursion_sizes(polys: &[Vec<(f32, f32)>], cell_km: f32) -> (f32, f32, f32, f64, usize) {
    let (sp, _) = coast_spurs(polys, cell_km, 2.0 * cell_km, cell_km);
    let l = sorted(sp.iter().map(|s| s.len_km / cell_km).collect::<Vec<f32>>());
    if l.is_empty() {
        return (f32::NAN, f32::NAN, f32::NAN, f64::NAN, 0);
    }
    let m = l.iter().map(|&x| x as f64).sum::<f64>() / l.len() as f64;
    let var = l.iter().map(|&x| (x as f64 - m).powi(2)).sum::<f64>() / l.len() as f64;
    (pct(&l, 0.10), pct(&l, 0.50), pct(&l, 0.90), var.sqrt() / m.max(1e-9), l.len())
}

/// **The PERIODIC negative control, built ON THE GRID** and not as a polygon.
///
/// ⚠️ The first attempt built it as a polygon with `periodic_polygon` and it was DEGENERATE: a
/// 2-px tooth at a 14-px period has a neck 7 px wide, so `coast_spurs` (neck ≤ 1 cell) found
/// **zero** excursions on it and R came out 0.000 — the control could not exercise the very
/// instruments it was built to calibrate. Finding 76 measured the fur as excursions **2–5 cells
/// long with a neck at or under one cell**, i.e. narrow SPIKES, and that is what this builds:
/// a filled circle plus one-cell-wide radial spikes of `spike_cells`, spaced `lambda_cells` apart
/// along the circumference.
pub fn spiky_circle_mask(w: usize, r: f64, spike_cells: f64, lambda_cells: f64) -> Vec<bool> {
    let c = w as f64 / 2.0;
    let mut m = vec![false; w * w];
    for y in 0..w {
        for x in 0..w {
            let (dx, dy) = (x as f64 + 0.5 - c, y as f64 + 0.5 - c);
            if dx * dx + dy * dy <= r * r {
                m[y * w + x] = true;
            }
        }
    }
    let teeth = (std::f64::consts::TAU * r / lambda_cells).round().max(4.0) as usize;
    for j in 0..teeth {
        let a = std::f64::consts::TAU * j as f64 / teeth as f64;
        let (ca, sa) = (a.cos(), a.sin());
        let mut t = 0.0f64;
        while t <= spike_cells {
            let (px, py) = (c + (r + t) * ca, c + (r + t) * sa);
            let (ix, iy) = (px.round() as i64, py.round() as i64);
            if ix >= 0 && iy >= 0 && (ix as usize) < w && (iy as usize) < w {
                m[iy as usize * w + ix as usize] = true;
            }
            t += 0.4; // sub-cell stepping so the spike is continuous on the raster
        }
    }
    m
}

/// **The ISOTROPIC negative control**: a RANDOMISED Koch snowflake — same 1/3 subdivision, so the
/// dimension is still `log 4 / log 3 = 1.2619`, but the apex is thrown to a random SIDE of each
/// edge. That removes the periodicity and the preferred orientation while keeping roughness at
/// every scale down to the smallest generation (`side / 3^g`), which is what a cell-scale
/// instrument needs and what a low-amplitude fBm lobe cannot give it (see
/// [`fbm_polygon`]'s note).
pub fn random_koch_polygon(w: usize, side: f64, generations: usize, seed: u64) -> Vec<(f64, f64)> {
    let c = w as f64 / 2.0;
    let r = side / 3f64.sqrt();
    let mut st = seed | 1;
    let mut rnd = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut pts: Vec<(f64, f64)> = (0..3)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / 3.0 - std::f64::consts::FRAC_PI_2;
            (c + r * a.cos(), c + r * a.sin())
        })
        .collect();
    for _ in 0..generations {
        let n = pts.len();
        let mut next = Vec::with_capacity(n * 4);
        for i in 0..n {
            let (a, b) = (pts[i], pts[(i + 1) % n]);
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let p = (a.0 + dx / 3.0, a.1 + dy / 3.0);
            let q = (a.0 + 2.0 * dx / 3.0, a.1 + 2.0 * dy / 3.0);
            let (ux, uy) = (q.0 - p.0, q.1 - p.1);
            // the apex, thrown to one side or the other with equal probability
            let sgn = if rnd() < 0.5 { 1.0f64 } else { -1.0 };
            let (co, si) = (0.5f64, sgn * (3f64.sqrt() / 2.0));
            let m = (p.0 + ux * co - uy * si, p.1 + ux * si + uy * co);
            next.push(a);
            next.push(p);
            next.push(m);
            next.push(q);
        }
        pts = next;
    }
    pts
}

// ══ ADR Finding 95's criteria chain, and ADR Finding 97's anisotropy reader ══
//
// Finding 96 copied the criteria chain verbatim out of `f95_long_run.rs:159-310` and recorded
// that it WAS a copy. Finding 97 needs the same instrument, so it moves here and both benches
// read the same code. `f95_long_run` itself is NOT refactored to call it — that bench costs 5 h
// and must not be touched to serve a later round — so this is still a copy of it, and any drift
// is a defect of THIS file.
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, clip_rivers_to_lakes,
    resolve_exorheic_without_outlet, runoff_accumulation, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY};

const CELL_M: f32 = 400_000.0 / 8192.0;

fn median(v: &[f32]) -> f32 {
    pct(&sorted(v.to_vec()), 0.50)
}

/// Max-of-8 slope in degrees, in metres. Shared by the Finding 96 and 97 benches.
pub fn slope_deg(f: &GridF32, k: usize, n2m: f32, w: usize, h: usize) -> f32 {
    let (x, y) = ((k % w) as i32, (k / w) as i32);
    let mut best = 0.0f32;
    for d in 0..8 {
        let nx = (x + D8_DX[d]).rem_euclid(w as i32) as usize;
        let ny = (y + D8_DY[d]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[d] != 0 && D8_DY[d] != 0;
        let dist: f32 = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        best = best.max((f.data[k] - f.data[ny * w + nx]).abs() * n2m / dist);
    }
    best.atan().to_degrees()
}

/// One reading of the instrument over one population of windows.
#[derive(Clone, Default)]
pub struct Aniso {
    pub windows: usize,
    pub c_p50: f32,
    pub c_p90: f32,
    pub share_hi: f32,
    /// normalised Shannon entropy of the orientation histogram, 36 bins, 1.0 = uniform
    pub h_theta: f32,
    /// |<e^{2iθ}>| -- the dossier's own `axis R` (Finding 84): ONE preferred axis
    pub r2: f32,
    /// |<e^{8iθ}>| -- the D8-lattice harmonic: four axes 45° apart
    pub r8: f32,
    /// r8 weighted by coherence, so incoherent windows cannot vote on orientation
    pub r8_w: f32,
}

/// Structure tensor over NON-OVERLAPPING `win`-sized windows: the population is then exactly
/// "every fully-land window of side `win`", which is declarable and reproducible, and it needs no
/// 8192²-sized smoothing buffer.
///
/// `J = Σ ∇h ∇hᵀ` over the window; `C = (λ₁−λ₂)/(λ₁+λ₂) = d / mean` with
/// `d = √(((Jxx−Jyy)/2)² + Jxy²)`; `θ = ½·atan2(2Jxy, Jxx−Jyy)` is the dominant gradient
/// orientation, defined mod π (axial data).
pub fn aniso(f: &GridF32, land: &[bool], win: usize) -> Aniso {
    let (w, h) = (f.width, f.height);
    let nb = 36usize;
    let mut cs: Vec<f32> = Vec::new();
    let mut hist = vec![0f64; nb];
    let (mut s2r, mut s2i, mut s8r, mut s8i) = (0f64, 0f64, 0f64, 0f64);
    let (mut s8rw, mut s8iw, mut wsum) = (0f64, 0f64, 0f64);
    for by in (0..h.saturating_sub(win)).step_by(win) {
        'win: for bx in (0..w.saturating_sub(win)).step_by(win) {
            let (mut jxx, mut jyy, mut jxy) = (0f64, 0f64, 0f64);
            for y in by..by + win {
                for x in bx..bx + win {
                    // a window touching water is dropped entirely: the coastline's own step is
                    // not a terrain orientation, and Finding 76 measured the coast separately
                    if !land[y * w + x] {
                        continue 'win;
                    }
                    let (gx, gy) = f.gradient_at(x, y);
                    jxx += (gx * gx) as f64;
                    jyy += (gy * gy) as f64;
                    jxy += (gx * gy) as f64;
                }
            }
            let mean = 0.5 * (jxx + jyy);
            if mean <= 1e-20 {
                continue;
            }
            let d = (0.25 * (jxx - jyy) * (jxx - jyy) + jxy * jxy).sqrt();
            let c = (d / mean) as f32;
            let theta = 0.5 * (2.0 * jxy).atan2(jxx - jyy); // (-π/2, π/2]
            let t = if theta < 0.0 { theta + std::f64::consts::PI } else { theta };
            cs.push(c);
            hist[((t / std::f64::consts::PI * nb as f64) as usize).min(nb - 1)] += 1.0;
            s2r += (2.0 * t).cos();
            s2i += (2.0 * t).sin();
            s8r += (8.0 * t).cos();
            s8i += (8.0 * t).sin();
            s8rw += c as f64 * (8.0 * t).cos();
            s8iw += c as f64 * (8.0 * t).sin();
            wsum += c as f64;
        }
    }
    let n = cs.len().max(1) as f64;
    let tot: f64 = hist.iter().sum::<f64>().max(1.0);
    let hh: f64 = hist
        .iter()
        .filter(|&&v| v > 0.0)
        .map(|&v| {
            let p = v / tot;
            -p * p.ln()
        })
        .sum();
    let s = sorted(cs.clone());
    Aniso {
        windows: cs.len(),
        c_p50: pct(&s, 0.50),
        c_p90: pct(&s, 0.90),
        share_hi: 100.0 * cs.iter().filter(|&&c| c > 0.7).count() as f32 / cs.len().max(1) as f32,
        h_theta: (hh / (nb as f64).ln()) as f32,
        r2: ((s2r / n).hypot(s2i / n)) as f32,
        r8: ((s8r / n).hypot(s8i / n)) as f32,
        r8_w: ((s8rw / wsum.max(1e-9)).hypot(s8iw / wsum.max(1e-9))) as f32,
    }
}

/// The worst 1024² tile by `r8` at window 16 -- Finding 76's rule-9 crop trap says a fixed crop
/// flatters whichever candidate it was not chosen on, so the tile is chosen PER FIELD and its
/// coordinates are printed.
pub fn worst_tile(f: &GridF32, land: &[bool]) -> (usize, usize, Aniso) {
    let (w, h) = (f.width, f.height);
    let mut best = (0usize, 0usize, Aniso::default());
    for by in (0..h).step_by(1024) {
        for bx in (0..w).step_by(1024) {
            let mut sub = GridF32::new(1024, 1024, 0.0);
            let mut sl = vec![false; 1024 * 1024];
            for y in 0..1024 {
                for x in 0..1024 {
                    sub.data[y * 1024 + x] = f.data[(by + y) * w + bx + x];
                    sl[y * 1024 + x] = land[(by + y) * w + bx + x];
                }
            }
            let a = aniso(&sub, &sl, 16);
            if a.windows >= 64 && a.r8 > best.2.r8 {
                best = (bx, by, a);
            }
        }
    }
    best
}

/// A `side`-sized tile at a FIXED corner, with its land mask.
pub fn tile(f: &GridF32, land: &[bool], bx: usize, by: usize, side: usize) -> (GridF32, Vec<bool>) {
    let w = f.width;
    let mut g = GridF32::new(side, side, 0.0);
    let mut l = vec![false; side * side];
    for y in 0..side {
        for x in 0..side {
            g.data[y * side + x] = f.data[(by + y) * w + bx + x];
            l[y * side + x] = land[(by + y) * w + bx + x];
        }
    }
    (g, l)
}

/// **Directional semivariogram**, the A3 instrument: `gamma(lag, dir) = mean((h(x+lag·d) − h(x))²)/2`
/// along the four axial D8 directions. Stripes with a wavelength show up twice over — as an
/// ANISOTROPY RATIO between the four directions, and as a "hole effect" (a dip in gamma at the
/// wavelength). Reported normalised by the omnidirectional gamma at the same lag, so the terrain's
/// own relief scale divides out.
pub fn variogram(
    f: &GridF32,
    land: &[bool],
    n2m: f32,
    lags: &[usize],
) -> Vec<(usize, [f32; 4], f32)> {
    let (w, h) = (f.width, f.height);
    let dirs: [(i32, i32); 4] = [(1, 0), (1, 1), (0, 1), (-1, 1)];
    let mut out = Vec::new();
    for &lag in lags {
        let mut g = [0f64; 4];
        for (di, &(dx, dy)) in dirs.iter().enumerate() {
            let (mut s, mut c) = (0f64, 0usize);
            for y in (0..h).step_by(3) {
                for x in (0..w).step_by(3) {
                    let (nx, ny) = (x as i32 + dx * lag as i32, y as i32 + dy * lag as i32);
                    if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                        continue;
                    }
                    let (a, b) = (y * w + x, ny as usize * w + nx as usize);
                    if !land[a] || !land[b] {
                        continue;
                    }
                    let d = ((f.data[b] - f.data[a]) * n2m) as f64;
                    s += d * d;
                    c += 1;
                }
            }
            // the diagonals span 1.41x the distance: normalise to distance so the four are
            // comparable, otherwise the diagonals read rougher for a trivial reason
            let scale = if dx != 0 && dy != 0 { 2.0f64 } else { 1.0 };
            g[di] = 0.5 * s / (c.max(1) as f64) / scale;
        }
        let mean = g.iter().sum::<f64>() / 4.0;
        let ratio =
            if mean > 0.0 { (g.iter().cloned().fold(0.0, f64::max) / mean) as f32 } else { 0.0 };
        out.push((lag, [g[0] as f32, g[1] as f32, g[2] as f32, g[3] as f32], ratio));
    }
    out
}

/// Finding 95's criteria chain, **copied VERBATIM** from `f95_long_run.rs:159-310` (commit
/// `751e0bf`) so table C reads the SAME instruments as the oracle row. Deliberately NOT
/// refactored into a shared helper: `f95_long_run` costs 5 h and must not be touched to serve
/// this round. If the two ever drift, this is a COPY and the drift is a defect of this bench.
///
/// `f` = the eroded field before the breach - `bf` = its breach - `pre` = the pre-incision
/// reference the cut is measured against.
#[allow(clippy::too_many_arguments)]
pub fn f95_criteria(
    f: &GridF32,
    bf: &GridF32,
    pre: &GridF32,
    ref_p50: f32,
    ss: &SteinSteinParams,
    on: &C1DrainageConfig,
    cell_km2: f32,
    n2m: f32,
    w: usize,
    h: usize,
) -> Crit {
    let n = w * h;
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], ss);
    // ── paired relief against the DELIVERED reference ──
    let com: Vec<usize> = (0..n).filter(|&k| f.data[k] > SEA && pre.data[k] > SEA).collect();
    let dd = sorted(com.iter().map(|&k| m(&f, k)).collect::<Vec<f32>>());
    let cut: Vec<f32> = com.iter().map(|&k| m(&pre, k) - m(&f, k)).collect();
    let p50 = pct(&dd, 0.50);
    eprintln!(
        "   RELIEF paired (n = {}): p10/p50/p90 **{:.1} / {:.1} / {:.1} m** · median cut \
         {:.1} m · **{:+.1} % against the delivered {ref_p50:.1} m ⇒ {}**",
        com.len(),
        pct(&dd, 0.10),
        p50,
        pct(&dd, 0.90),
        median(&cut),
        100.0 * (p50 - ref_p50) / ref_p50,
        if (100.0 * (p50 - ref_p50) / ref_p50).abs() <= 10.0 {
            "INSIDE ±10 %"
        } else {
            "**OUTSIDE ±10 %**"
        }
    );

    // ── the lake inventory: fraction, canyon class, integration ──
    let climate = c1_climate_placed(&bf, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let pre_d = c1_drainage_windowed(&f, None, &on, &ss, DOMAIN_KM);
    let mut dr = c1_drainage_windowed(&bf, Some(&dclim), &on, &ss, DOMAIN_KM);
    dr.lakes = pre_d.lakes;
    dr.lake_map = pre_d.lake_map;
    let li = std::mem::take(&mut dr.lakes);
    dr.lakes = apply_lake_water_balance(
        &bf,
        &dr.flow,
        &dclim,
        cell_km2,
        &ss,
        &li,
        &mut dr.lake_map,
        None,
        w,
        h,
    );
    let racc = runoff_accumulation(&bf, &dr.flow, &dclim, cell_km2, None, None, w, h);
    let pm = dr.lake_map.clone();
    let bs = below_sea_basin_lakes_infil(&bf, &dclim, &on, &ss, DOMAIN_KM, Some(&pm), None);
    let mut d_before: HashMap<u32, usize> = HashMap::new();
    for &id in &dr.lake_map {
        if id != 0 && id < 1_000_001 {
            *d_before.entry(id).or_default() += 1;
        }
    }
    let mut d_after: HashMap<u32, usize> = HashMap::new();
    for k in 0..n {
        if bs.lake_map[k] != 0 {
            dr.lake_map[k] = bs.lake_map[k];
        } else if dr.lake_map[k] != 0 && dr.lake_map[k] < 1_000_001 {
            *d_after.entry(dr.lake_map[k]).or_default() += 1;
        }
    }
    let absorbed: HashSet<u32> =
        d_before.keys().copied().filter(|id| d_after.get(id).copied().unwrap_or(0) == 0).collect();
    dr.lakes.retain(|l| l.base.id >= 1_000_001 || !absorbed.contains(&l.base.id));
    dr.lakes.extend(bs.lakes.iter().cloned());
    clip_rivers_to_lakes(&mut dr);
    let inv: HashSet<u32> = dr.lakes.iter().map(|l| l.base.id).collect();
    for sw in &bs.spillways {
        let (lx, ly) = *sw.points.last().unwrap_or(&(0, 0));
        dr.push_segment(SegmentRow {
            segment: ymir_core::terrain::flow::RiverSegment {
                points: sw.points.clone(),
                strahler_order: 1,
                avg_flow: 0.0,
                max_flow: 0.0,
                basin_id: 0,
                upstream: vec![],
                downstream: None,
            },
            drainage_km2: sw.drainage_km2,
            navigability: sw.navigability,
            discharge_m3s: sw.discharge_m3s,
            width_m: sw.width_m,
            profile_m: sw.profile_m.clone(),
            catchment_cells: {
                let kk = ly as usize * w + lx as usize;
                dr.flow.accumulation.data.get(kk).copied().unwrap_or(0.0)
            },
            discharge_profile_m3s: vec![sw.discharge_m3s; sw.points.len()],
            kind: SegmentKind::Spillway,
            source_lake: inv.contains(&sw.lake_id).then_some(sw.lake_id),
        });
    }
    let relab = resolve_exorheic_without_outlet(&mut dr);

    let land_km2 = (0..n).filter(|&k| bf.data[k] > SEA).count() as f32 * cell_km2;
    let water_cells = (0..n).filter(|&k| dr.lake_map[k] != 0).count();
    let water_km2 = water_cells as f32 * cell_km2;
    let dump = DUMP.load(std::sync::atomic::Ordering::Relaxed);
    let mut klass = 0usize;
    let mut cut_class = 0usize;
    let mut scanned = 0usize;
    let mut ty: HashMap<String, usize> = HashMap::new();
    for l in &dr.lakes {
        *ty.entry(format!("{:?}", l.lake_type)).or_default() += 1;
    }
    for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
        let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
        if cells.is_empty() {
            continue;
        }
        scanned += 1;
        // ADR Finding 107 -- TWO columns, and the GATE is the second.
        //
        // **CUT** (`pre-incision − delivered`) is what Findings 87-104 called the canyon class. It
        // cannot tell an over-dug closed depression from a deep valley the drainage cuts straight
        // through, because a through-cut valley has a large cut too.
        //
        // **FILL** (`filled − raw`, the priority flood on the ERODED field) is topological: it is
        // how far the floor lies BELOW ITS OWN SILL. A through-cut valley has `filled == raw` along
        // its bed whatever its depth, so it scores 0 and is not in the class — which is what must
        // hold the day an uplift or a glaciation digs a real one.
        let cuts: Vec<f32> = cells.iter().map(|&k| m(&pre, k) - m(&bf, k)).collect();
        let fills: Vec<f32> =
            cells.iter().map(|&k| (pre_d.flow.filled.data[k] - f.data[k]) * n2m).collect();
        // ⚠️ ADR Finding 107-D -- `max` is the LITERAL reading of "the floor below its
        // OWN sill": over a body flooded to a flat sill `filled` is constant and `raw` is
        // least at the floor, so `max(filled - raw)` IS sill-minus-floor. **But it is a
        // ONE-CELL statistic**, and this campaign has now been bitten three times by a
        // threshold read at the resolution of its own quantity. If a footprint is NOT one
        // closed depression -- a below-sea basin draining to the ocean, say -- `max` reports
        // an unrelated pit that happens to lie inside it.
        //
        // `share` is the control THE FIRST RENDERING DEMANDED: body 1000001 scored 237.97 m
        // and its fill zoom came out BLACK -- no basin, thin filaments. The deep share is the
        // fraction of the body actually lying below its sill, which separates a real hollow
        // (a blob) from a filament the max crossed.
        let fill = fills.iter().copied().fold(0.0f32, f32::max);
        let share = 100.0 * fills.iter().filter(|&&v| v > 50.0).count() as f32 / fills.len() as f32;
        let inside: HashSet<usize> = cells.iter().copied().collect();
        let mut ring: HashSet<usize> = HashSet::new();
        for &k in &cells {
            for d in 0..8 {
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                let nk = ny * w + nx;
                if !inside.contains(&nk) {
                    ring.insert(nk);
                }
            }
        }
        let rim: Vec<f32> = ring.iter().map(|&k| slope_deg(&bf, k, n2m, w, h)).collect();
        let steep = !rim.is_empty() && median(&rim) > 30.0;
        let (is_cut, is_dug) = (median(&cuts) > 50.0 && steep, fill > 50.0 && steep);
        if is_cut {
            cut_class += 1;
        }
        if is_dug {
            klass += 1;
        }
        if dump {
            BODIES.lock().expect("bodies").push(Body {
                id: l.base.id,
                cells: cells.clone(),
                km2: cells.len() as f32 * cell_km2,
                cut: median(&cuts),
                fill,
                share,
                rim: if rim.is_empty() { 0.0 } else { median(&rim) },
                is_cut,
                is_dug,
            });
            eprintln!(
                "     body {:>7} · {:>8.2} km² · cut p50 **{:>8.2} m** · fill max **{:>8.2} \
                 m** · fill p50 {:>7.2} m · **deep share {:>5.1} %** · rim p50 {:>5.1}° ⇒ cut {} / over-dug {}",
                l.base.id,
                cells.len() as f32 * cell_km2,
                median(&cuts),
                fill,
                median(&fills),
                share,
                median(&rim),
                if is_cut { "YES" } else { " no" },
                if is_dug { "YES" } else { " no" },
            );
        }
    }
    if dump {
        // ── the second control: a THROUGH-CUT bed must read fill ~ 0 ──────────
        // Finding 107's discriminant only earns the rename if it clears the river the gate must
        // never kill. `filled == raw` along a bed however deep the valley; `filled - raw` is the
        // depth below the bed's OWN sill, which a through-flowing course does not have.
        for want in [SegmentKind::Spillway, SegmentKind::Watercourse] {
            let Some(i) =
                (0..dr.rivers.segments.len()).filter(|&i| dr.segment_kind[i] == want).max_by(
                    |&a, &b| dr.segment_discharge_m3s[a].total_cmp(&dr.segment_discharge_m3s[b]),
                )
            else {
                continue;
            };
            let pts = &dr.rivers.segments[i].points;
            let fl: Vec<f32> = pts
                .iter()
                .map(|&(x, y)| {
                    let k = y as usize * w + x as usize;
                    (pre_d.flow.filled.data[k] - f.data[k]) * n2m
                })
                .collect();
            let drop = match (pts.first(), pts.last()) {
                (Some(&(ax, ay)), Some(&(bx, by))) => {
                    m(&f, ay as usize * w + ax as usize) - m(&f, by as usize * w + bx as usize)
                }
                _ => 0.0,
            };
            let p50 = median(&fl);
            eprintln!(
                "     {:?} max **{:.2} m³/s** · {} cells · drop **{drop:.1} m** · FILL along                  the bed: p50 **{p50:.4} m** · max **{:.4} m** ⇒ **{}**",
                want,
                dr.segment_discharge_m3s[i],
                fl.len(),
                fl.iter().copied().fold(0.0f32, f32::max),
                if p50 <= 1.0 {
                    "NOT in the class -- the bed is through-cut"
                } else {
                    "**IN the class -- the discriminant FAILS**"
                }
            );
        }
    }
    let q_unres: f64 = dr
        .lakes
        .iter()
        .filter(|l| l.lake_type == LakeType::Unresolved)
        .map(|l| {
            let c: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
            runoff_km2_to_m3s(c.iter().map(|&k| racc[k]).fold(0.0f32, f32::max)) as f64
        })
        .sum();
    let wc = water_class(&bf, SEA);
    eprintln!(
        "   LAKES: **{} bodies · {water_km2:.0} km² = {:.2} % of {land_km2:.0} km² of \
         land** · {:?} · OVER-DUG **{klass} of {scanned}** (old cut class {cut_class}) · relabels {} · \
         `to_nothing` {} · Unresolved inflow {q_unres:.1} m³/s · wc==2 cells {}",
        dr.lakes.len(),
        100.0 * water_km2 / land_km2,
        ty,
        relab.len(),
        bs.termination.to_nothing,
        (0..n).filter(|&k| wc[k] == 2).count()
    );
    Crit {
        p50,
        cut_class,
        below_sea_spillways: bs.spillways.len(),
        max_catchment_km2: dr.segment_drainage_km2.iter().copied().fold(0.0f32, f32::max),
        klass,
        scanned,
        to_nothing: bs.termination.to_nothing,
        lake_pct: 100.0 * water_km2 / land_km2,
        unresolved: dr.lakes.iter().filter(|l| l.lake_type == LakeType::Unresolved).count(),
        q_unres,
    }
}

/// ADR Finding 107 -- when set, `f95_criteria` prints ONE LINE PER BODY (both columns) and the
/// fill profile of the highest-discharge spillway and watercourse. The controls for the gate
/// redefinition read the SAME `dr` the class is scored on, so they cannot drift from the gate.
static DUMP: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Turn the Finding 107 per-body dump on.
pub fn set_dump(on: bool) {
    DUMP.store(on, std::sync::atomic::Ordering::Relaxed);
}

/// ADR Finding 107 — one scored body, with the FOOTPRINT the gate used.
///
/// The author asked where the bodies the round talks about can be LOOKED AT, and the honest
/// answer was "nowhere": the À VALIDER list was ids and numbers with no coordinate. A renderer
/// that re-derives its own inventory would draw a different set of bodies from the one the gate
/// scored (Finding 101-B: 26 against 29), so the footprints are captured HERE, inside the scoring
/// loop, and handed to the renderer verbatim.
pub struct Body {
    pub id: u32,
    pub cells: Vec<usize>,
    pub km2: f32,
    pub cut: f32,
    pub fill: f32,
    /// ADR Finding 107-D -- percentage of the body's cells lying more than 50 m below its sill.
    /// The control on [`Body::fill`], which is a ONE-CELL maximum.
    pub share: f32,
    pub rim: f32,
    pub is_cut: bool,
    pub is_dug: bool,
}

static BODIES: std::sync::Mutex<Vec<Body>> = std::sync::Mutex::new(Vec::new());

/// Drain the bodies captured by the last dumping `f95_criteria` call.
pub fn take_bodies() -> Vec<Body> {
    std::mem::take(&mut BODIES.lock().expect("bodies"))
}

/// Finding 95's criteria, as one row of table C.
pub struct Crit {
    pub p50: f32,
    /// ADR Finding 107 -- the OLD, cut-based count that Findings 87-104 printed as "canyon class".
    /// Kept so the redefinition is auditable rather than silent: every historical number in the
    /// dossier is this column, and `klass` is now the topological one.
    pub cut_class: usize,
    /// ADR Finding 45's DRAINAGE-INTEGRITY column: the number of SPILLWAYS -- the traced outflows
    /// of below-sea basins over their cols, i.e. `bs.spillways.len()`.
    ///
    /// ⚠️ **It does NOT count basins, and Finding 101 had to prove that the hard way.** Finding 99
    /// named it `below_sea_basins`, Finding 100 built an admissibility gate on it (+-20 % or it is
    /// Finding 45 again), and that gate nearly disqualified the campaign's best row on a move from
    /// 11 to 14. Finding 101 then counted the `wc == 2` COMPONENTS by bottom cell: **20 in the
    /// delivered field, 20 in the candidate, 20 of 20 matched, none new**. The terrain did not
    /// change; the ROUTING did. A column that decides admissibility must say what it counts.
    ///
    /// Added at Finding 99 because Finding 45 measured a 16x diffusion, saw the hypsometry move
    /// +17 m ("small and in the unhelpful direction"), reported it as harmless -- and this number
    /// had gone 43 -> 1994. Its own verdict: *"the patch did not change the altitude, it destroyed
    /// the drainage"*. Any row that changes the hillslope term must show this column or it is
    /// repeating that error.
    pub below_sea_spillways: usize,
    /// The largest river catchment, km2 (signified). Finding 45 watched it fall 110 -> 48 km2
    /// while the altitude held: the second half of the same collapse.
    pub max_catchment_km2: f32,
    /// ADR Finding 107 -- **the OVER-DUG DEPRESSION count**, the gate. A body >= 1 km2 whose floor
    /// sits more than 50 m below its OWN sill (`filled - raw`, the priority flood of Findings
    /// 13/14 on the ERODED field) behind walls at p50 > 30 deg.
    ///
    /// ⚠️ It is NOT what Findings 87-104 printed under this name: those scored the **cut**,
    /// `pre-incision - delivered`, which is [`Crit::cut_class`]. A deep valley the drainage runs
    /// straight through has a large cut and a fill of ZERO, and a cut-based gate would kill it --
    /// including the day an uplift or a glaciation digs a real gorge. The discriminant is
    /// topological, and it was free: the chain already computes `filled`.
    pub klass: usize,
    pub scanned: usize,
    pub to_nothing: usize,
    pub lake_pct: f32,
    pub unresolved: usize,
    pub q_unres: f64,
}
