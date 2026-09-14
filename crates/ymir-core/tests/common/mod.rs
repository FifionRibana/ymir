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
    /// ADR Finding 80/83 -- OVERRIDE the base-level epsilon, in metres above sea level.
    /// `None` = the SHIPPED `RELIEF_V3_BASE_LEVEL_M` (0.5 m). Since Finding 83 the bound
    /// is production, so `None` no longer means "no bound" -- use `base_level_off`.
    pub base_level_m: Option<f32>,
    /// ADR Finding 83 -- the A/B control: drop the shipped bound and rebuild the
    /// pre-Finding-83 world (the one with the coastal fringe). This is what the benches
    /// must use for the "DELIVERED (pre-83)" column; `Knobs::shipped()` is now BOUNDED.
    pub base_level_off: bool,
    /// ADR Finding 83-B2 -- `BaseLevelFloor::free_above_km2`, the estuary gate, in km2 of
    /// the SIMULATED grid. **PROXY**, bench-only; production is `None`.
    pub a_est_km2: Option<f32>,
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
