//! ADR Finding 165-C -- **the physics level driven by C1's history** (the author: « Le relief est fait à 64² on
//! devrait avoir toutes ces zones déjà formées »; « quelque chose qui soit physique »).
//!
//! **The record** ([`run_c1_recorded`]): C1's own loop (`run_with_closures_observed`), with a READ-ONLY callback that
//! keeps the raw coarse altitude every `every` steps (plus the state before the first step), the cells classified
//! Convergent at any step (the sutures, F151's fossil belts) and the cells whose `plate_type` changed. C1's outputs are
//! bit-identical with and without it (`the_record_leaves_c1_bit_identical`).
//!
//! **The material frame**: C1's transport is Eulerian upwind on a fixed grid, and with `rigid_continental_crust` the
//! continental crust has zero velocity (`apply_continental_rigidity`, a no-flux boundary), so a continental column
//! never leaves its cell: the Eulerian record at a continental cell IS its material history. Only oceanic crust moves;
//! the cells that changed type are reported ([`C1History::type_changed`]).
//!
//! **The drive** ([`physics_history`]): N_phys steps over T years. U(x, t) = Δh_iso / Δt between the bracketing
//! snapshots (piecewise constant), upsampled like h_iso; the sea is not held (it follows U); the stream power, the
//! talus and the diffusion of the current level; the isostatic rebound of each step's erosion
//! ([`super::Rebound`]). No equilibrium criterion, no deposition. Declared in
//! `docs/reports/relief_method/f165_continent/f165_declared.md` (C).

use super::{
    CascadeConfig, CascadeProgress, Drive, LevelRecord, Rebound, amplify, run_loop, upsample_to,
};
use crate::grid::GridF32;
use crate::tectonics::isostasy::IsostasyConfig;
use crate::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries};
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::production_upscale::{
    c1_coarse_raw_altitude, c1_normalize_coarse, c1_production_altitude_craton, calibrate_to_land_fraction,
};
use crate::tectonics_c1::state::C1State;
pub use crate::tectonics_c1::time_loop::C1Term;
use crate::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures_observed, run_with_closures_terms};
use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::Field2D;
use std::cell::RefCell;

/// C1's recorded history at its own grid (64²).
#[derive(Clone, Debug)]
pub struct C1History {
    /// The C1 step each snapshot was taken after (0 = before the first step).
    pub steps: Vec<usize>,
    pub n_steps: usize,
    /// The raw coarse altitude (`c1_coarse_raw_altitude`) of each snapshot.
    pub raw: Vec<GridF32>,
    /// Cells classified Convergent at any step.
    pub suture: Vec<bool>,
    /// Cells whose `plate_type` differs at the end from the start.
    pub type_changed: Vec<bool>,
}

/// Run C1 (exactly `run_with_closures_observed`) and record its history every `every` steps.
pub fn run_c1_recorded(
    state: &mut C1State,
    kin: &mut PlateKinematics,
    cfg: &C1TimeLoopConfig,
    closures: &C1Closures,
    every: usize,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
) -> C1History {
    let n = state.nx() * state.ny();
    let type0: Vec<_> = state.plate_type.data().to_vec();
    let mut raw = vec![c1_coarse_raw_altitude(state, iso, ss)];
    let mut steps = vec![0usize];
    let mut suture = vec![false; n];
    run_with_closures_observed(state, kin, cfg, closures, |step, st, k| {
        let b = classify_boundaries(&st.plate_id, k);
        for (s, t) in suture.iter_mut().zip(b.boundary_type.data()) {
            if *t == BoundaryType::Convergent {
                *s = true;
            }
        }
        if (step + 1) % every == 0 {
            raw.push(c1_coarse_raw_altitude(st, iso, ss));
            steps.push(step + 1);
        }
    });
    let type_changed = state.plate_type.data().iter().zip(&type0).map(|(a, b)| a != b).collect();
    C1History { steps, n_steps: cfg.n_steps, raw, suture, type_changed }
}

impl C1History {
    /// The snapshots normalised like the final coarse field (`c1_coarse_normalized_altitude`): with a target land
    /// fraction, the FINAL snapshot's sea-level shift is applied to all, so the sea level does not move in time.
    pub fn normalized(&self, target_land_fraction: Option<f32>) -> Vec<GridF32> {
        let shift = target_land_fraction.map_or(0.0, |f| {
            let mut d = self.raw.last().unwrap().data.clone();
            -calibrate_to_land_fraction(&mut d, f)
        });
        self.raw
            .iter()
            .map(|r| {
                let mut g = r.clone();
                for v in &mut g.data {
                    *v += shift;
                }
                c1_normalize_coarse(g, None)
            })
            .collect()
    }
}

/// ADR Finding 165-C2/C3 -- how the uplift is spread in time.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Timing {
    /// Δh_iso between the bracketing snapshots, over their interval.
    History,
    /// Control 3: the whole Δh_iso (final − initial) in the last 10 % of T.
    Late10,
    /// Negative control 1: U ∝ the FINAL h_iso, constant, to equilibrium, no rebound, U₀ calibrated on
    /// [`STEADY_PEAK_M`] -- the current level ([`super::physics_level`], the same loop with `Drive::steady`).
    Steady,
}

/// ADR Finding 166-B -- C1's history with the excluded terms taken out of `s`: the snapshots are h_src = C1's isostasy
/// of s_src = s − Σ(excluded terms' Δs), with C1's `age`, `plate_type` and craton mask at the time.
#[derive(Clone, Debug)]
pub struct C1Sources {
    /// `raw` holds h_src (raw coarse altitude) at the snapshots; the sutures and the type changes as F165.
    pub hist: C1History,
    /// Σ Δs per term over the run (the run's grid, row-major).
    pub terms: Vec<(C1Term, Vec<f64>)>,
    pub excluded: Vec<C1Term>,
    /// s at the start and at the end (the run's).
    pub s0: Vec<f64>,
    pub s_final: Vec<f64>,
    /// The craton mask (the isostasy's buoyancy).
    pub craton: Vec<bool>,
}

/// ADR Finding 167-B -- the thickness → altitude mapping of the record: C1's own (F165–F166), or the physical Airy
/// mapping with an absolute sea level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mapping {
    C1,
    Airy(Airy),
}

/// ADR Finding 167-B -- **Airy isostasy with an absolute sea level** (Whitehead & Clift 2009, p. 3 and eq. 1–2):
/// e = b · L · S̃ − H₀ above sea level, b = (ρm − ρc)/ρm, and water-loaded below it (× ρm/(ρm − ρw)), one continuous
/// monotone function of S̃; H₀ from the anchor e(S̃ = 1) = `anchor_km`. Declared in
/// `docs/reports/relief_method/f167_isostasy/f167_declared.md` (B).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Airy {
    /// km of crust per unit S̃ (C1's convention, S̃ = 1 ↔ 35 km).
    pub km_per_s: f64,
    pub rho_c: f64,
    pub rho_m: f64,
    pub rho_w: f64,
    /// The altitude of S̃ = 1 (km).
    pub anchor_km: f64,
    /// The S̃ the oceanic-plate cells are read at (C1's « phantom oceanic advective spike » is not a thickness); `None`
    /// = every cell's own S̃.
    pub oceanic_s: Option<f64>,
    /// C1's Gaussian blur of its altitude (`IsostasyConfig::altitude_smoothing_sigma`), kept.
    pub sigma: f32,
}

impl Airy {
    /// F167's declaration: ρc 2 800, ρm 3 300, ρw 1 030 kg/m³, 35 km per S̃, S̃ = 1 → +0.45 km, the oceanic cells at S̃ = 0.2,
    /// C1's σ = 0.5.
    pub fn declared() -> Self {
        Self { km_per_s: 35.0, rho_c: 2800.0, rho_m: 3300.0, rho_w: 1030.0, anchor_km: 0.45, oceanic_s: Some(0.2), sigma: 0.5 }
    }
    /// The altitude (km) of a column of thickness S̃.
    pub fn altitude_km(&self, s: f64) -> f64 {
        let b = (self.rho_m - self.rho_c) / self.rho_m;
        let h0 = b * self.km_per_s - self.anchor_km;
        let e = b * self.km_per_s * s - h0;
        if e >= 0.0 { e } else { e * self.rho_m / (self.rho_m - self.rho_w) }
    }
}

/// ADR Finding 167-B -- the Airy mapping of a C1 state's thickness, in C1's raw coarse units (m / `depth_scale_m`), so
/// [`C1History::normalized`] treats it as C1's.
pub fn airy_raw(s: &Field2D, plate_type: &PlateTypeField, a: &Airy, depth_scale_m: f64) -> GridF32 {
    let (nx, ny) = (s.nx(), s.ny());
    let data: Vec<f32> = s
        .data()
        .iter()
        .zip(plate_type.data())
        .map(|(v, t)| {
            let sv = match (a.oceanic_s, t) {
                (Some(o), PlateType::Oceanic) => o,
                _ => *v,
            };
            (a.altitude_km(sv) * 1000.0 / depth_scale_m) as f32
        })
        .collect();
    let g = GridF32::from_vec(nx, ny, data);
    if a.sigma > 0.0 { g.gaussian_blur(a.sigma) } else { g }
}

/// ADR Finding 166-B1 -- run C1 (exactly; `run_with_closures_terms`) and record h_src every `every` steps, the terms in
/// `exclude` taken out of `s`. With nothing excluded the snapshots are [`run_c1_recorded`]'s, bit for bit
/// (`the_sources_with_every_term_are_f165s_history`).
#[allow(clippy::too_many_arguments)]
pub fn run_c1_sources(
    state: &mut C1State,
    kin: &mut PlateKinematics,
    cfg: &C1TimeLoopConfig,
    closures: &C1Closures,
    every: usize,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    exclude: &[C1Term],
) -> C1Sources {
    run_c1_sources_mapped(state, kin, cfg, closures, every, iso, ss, exclude, Mapping::C1)
}

/// ADR Finding 167-B -- [`run_c1_sources`] with the thickness → altitude mapping chosen; `Mapping::C1` is
/// [`run_c1_sources`] exactly (`the_c1_mapping_is_f166s_record`).
#[allow(clippy::too_many_arguments)]
pub fn run_c1_sources_mapped(
    state: &mut C1State,
    kin: &mut PlateKinematics,
    cfg: &C1TimeLoopConfig,
    closures: &C1Closures,
    every: usize,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    exclude: &[C1Term],
    mapping: Mapping,
) -> C1Sources {
    let (nx, ny) = (state.nx(), state.ny());
    let n = nx * ny;
    let s0 = state.s.data().to_vec();
    let type0: Vec<_> = state.plate_type.data().to_vec();
    let removed = RefCell::new(vec![0f64; n]);
    let sums = RefCell::new(C1Term::ALL.iter().map(|&t| (t, vec![0f64; n])).collect::<Vec<_>>());
    let src_alt = |st: &C1State, removed: &[f64]| -> GridF32 {
        if let Mapping::Airy(a) = &mapping {
            let s_src = Field2D::from_vec(nx, ny, st.s.data().iter().zip(removed).map(|(s, r)| s - r).collect());
            return airy_raw(&s_src, &st.plate_type, a, ss.depth_scale_m);
        }
        if exclude.is_empty() {
            return c1_coarse_raw_altitude(st, iso, ss);
        }
        let s_src = Field2D::from_vec(nx, ny, st.s.data().iter().zip(removed).map(|(s, r)| s - r).collect());
        c1_production_altitude_craton(&s_src, &st.age, &st.plate_type, st.cratonic_mask.data(), iso, ss)
    };
    let mut raw = vec![src_alt(state, &removed.borrow())];
    let mut steps = vec![0usize];
    let mut suture = vec![false; n];
    run_with_closures_terms(
        state,
        kin,
        cfg,
        closures,
        |step, st, k| {
            let b = classify_boundaries(&st.plate_id, k);
            for (s, t) in suture.iter_mut().zip(b.boundary_type.data()) {
                if *t == BoundaryType::Convergent {
                    *s = true;
                }
            }
            if (step + 1) % every == 0 {
                raw.push(src_alt(st, &removed.borrow()));
                steps.push(step + 1);
            }
        },
        &mut |term, before, after| {
            let mut sm = sums.borrow_mut();
            let acc = &mut sm.iter_mut().find(|x| x.0 == term).unwrap().1;
            let ex = exclude.contains(&term);
            let mut rm = removed.borrow_mut();
            for k in 0..n {
                let d = after[k] - before[k];
                acc[k] += d;
                if ex {
                    rm[k] += d;
                }
            }
        },
    );
    let type_changed = state.plate_type.data().iter().zip(&type0).map(|(a, b)| a != b).collect();
    C1Sources {
        hist: C1History { steps, n_steps: cfg.n_steps, raw, suture, type_changed },
        terms: sums.into_inner(),
        excluded: exclude.to_vec(),
        s0,
        s_final: state.s.data().to_vec(),
        craton: state.cratonic_mask.data().to_vec(),
    }
}

/// ADR Finding 166-A2 -- the isostatic datum (C1's raw sea level h_sea = h_min + fraction · (p_cap − h_min) of the
/// whole field, `isostasy.rs::compute_isostasy_inner`), recomputed for the diagnosis.
pub fn isostatic_datum(s: &[f64], craton: &[bool], iso: &IsostasyConfig) -> f64 {
    let b = 1.0 - iso.rho_crust as f64 / iso.rho_mantle as f64;
    let bc = iso.craton_rho_crust.map_or(b, |r| 1.0 - r as f64 / iso.rho_mantle as f64);
    let mut h: Vec<f64> = s.iter().zip(craton).map(|(v, c)| v * if *c { bc } else { b }).collect();
    let h_min = h.iter().copied().fold(f64::INFINITY, f64::min);
    h.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let cap = match iso.sea_level_mode {
        crate::tectonics::isostasy::SeaLevelMode::PercentileCapped { cap_percentile } => {
            h[((cap_percentile as f64) * (h.len() - 1) as f64).round() as usize]
        }
        _ => *h.last().unwrap(),
    };
    h_min + iso.sea_level_fraction as f64 * (cap - h_min)
}

/// ADR Finding 166-A2 -- a land cell's altitude (m) at thickness `s` and datum `h_sea`, by C1's land ramp
/// (`peak · (h − h_sea) / (land_ceiling − h_sea)`, the land ceiling the fixed S̃ reference), before the smoothing.
pub fn land_altitude_m(s: f64, craton: bool, h_sea: f64, iso: &IsostasyConfig) -> f64 {
    let b = 1.0 - iso.rho_crust as f64 / iso.rho_mantle as f64;
    let bc = iso.craton_rho_crust.map_or(b, |r| 1.0 - r as f64 / iso.rho_mantle as f64);
    let lc = iso.land_ref_thickness.map_or(2.0, |r| r as f64) * b;
    let h = s * if craton { bc } else { b };
    iso.max_elevation_m as f64 * (h - h_sea) / (lc - h_sea).max(1e-9)
}

/// ADR Finding 166 -- **T of the sources' drive**, calibrated once on the témoin (`f166_cascade`: 7 trials, 5.623 × 10⁶ yr
/// → 2 728 m).
pub const T_F166_YEARS: f64 = 5.623e6;

/// ADR Finding 165-C2 -- **T, the one constant calibrated on our worlds**: C1's history mapped to years, calibrated once
/// on the témoin so its peak at 256² falls in [2 700, 3 000] m (`f165_cascade`: 9 trials, 1.155 × 10⁶ yr → 2 852 m).
pub const T_F165_YEARS: f64 = 1.155e6;

/// The steady control's peak target (m), F160's « type Corse » value.
pub const STEADY_PEAK_M: f32 = 2850.0;

/// The historical level's settings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HistoryRun {
    /// T, the physical duration of C1's history (years): the ONE calibrated constant.
    pub t_years: f64,
    /// N_phys.
    pub steps: usize,
    pub rebound: Option<Rebound>,
    pub timing: Timing,
}

/// The historical level's result.
pub struct HistoryLevel {
    pub rec: LevelRecord,
    /// The initial field (norm), h_iso(t₀) upsampled plus the roughness.
    pub z0: GridF32,
    /// The snapshots at the level (norm), for the viz.
    pub snapshots: Vec<GridF32>,
}

/// ADR Finding 165-C2 -- the physics level at `n` cells driven by C1's history (`snaps`: the normalised 64²
/// snapshots, the first the initial state, the last the final one). The roughness is the current level's
/// (`cfg.roughness_m`, F163). `None` if cancelled.
pub fn physics_history(
    snaps: &[GridF32],
    n: usize,
    cfg: &CascadeConfig,
    run: &HistoryRun,
    progress: &mut dyn FnMut(CascadeProgress),
    cancel: &dyn Fn() -> bool,
) -> Option<HistoryLevel> {
    let n2m = cfg.norm_to_m();
    let up: Vec<GridF32> = snaps.iter().map(|s| upsample_to(s, n)).collect();
    if run.timing == Timing::Steady {
        let (rec, _) = super::physics_level(snaps.last().unwrap(), n, cfg, STEADY_PEAK_M, run.steps, progress, cancel)?;
        return Some(HistoryLevel { z0: rec.upscaled.clone(), rec, snapshots: up });
    }
    let mut z0 = up[0].clone();
    // the roughness of the current level (F163), on the land of the initial state
    if cfg.roughness_m > 0.0 {
        let cell_km = cfg.cell_km(n);
        let g = amplify::noise01(n, cell_km, 8.0 * cell_km, 2, cfg.seed, "cascade_phys_rough");
        let land: Vec<usize> = (0..n * n).filter(|&k| z0.data[k] > cfg.sea_level).collect();
        let mean = land.iter().map(|&k| g[k] as f64).sum::<f64>() / land.len().max(1) as f64;
        let rms = (land.iter().map(|&k| (g[k] as f64 - mean).powi(2)).sum::<f64>() / land.len().max(1) as f64).sqrt();
        let floor = cfg.sea_level + 1.0 / n2m;
        for &k in &land {
            let r = ((g[k] as f64 - mean) / rms.max(1e-12)) as f32 * cfg.roughness_m / n2m;
            z0.data[k] = (z0.data[k] + r).max(floor);
        }
    }
    let intervals = up.len() - 1;
    let dt = run.t_years / run.steps as f64;
    let mut cfg_t = cfg.clone();
    cfg_t.dt_yr = dt as f32;
    // the uplift per step (norm): Δh over the interval × the step's share of it
    let per_step: Vec<Vec<f32>> = match run.timing {
        Timing::History => (0..intervals)
            .map(|j| {
                let f = (intervals as f64 / run.steps as f64) as f32;
                up[j + 1].data.iter().zip(&up[j].data).map(|(b, a)| (b - a) * f).collect()
            })
            .collect(),
        Timing::Late10 => {
            let late = ((run.steps as f64) * 0.1).round().max(1.0) as usize;
            vec![up[intervals].data.iter().zip(&up[0].data).map(|(b, a)| (b - a) / late as f32).collect()]
        }
        Timing::Steady => unreachable!(),
    };
    let zeros = vec![0f32; n * n];
    let steps = run.steps;
    let timing = run.timing;
    let drive = Drive {
        du: Box::new(move |i| {
            let v: &[f32] = match timing {
                Timing::History => &per_step[((i * intervals) / steps).min(intervals - 1)],
                Timing::Late10 => {
                    let late = ((steps as f64) * 0.1).round().max(1.0) as usize;
                    if i + late >= steps { &per_step[0] } else { &zeros }
                }
                Timing::Steady => unreachable!(),
            };
            std::borrow::Cow::Owned(v.to_vec())
        }),
        hold_sea: false,
        eq_stop: false,
        rebound: run.rebound,
    };
    let uplift = GridF32::new(n, n, 0.0);
    let rec = run_loop(z0.clone(), uplift, &cfg_t, n, run.steps, Some(0), &drive, progress, cancel)?;
    Some(HistoryLevel { rec, z0, snapshots: up })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use crate::tectonics_c1::time_loop::run_with_closures;

    fn c1(steps: usize) -> (C1State, PlateKinematics, C1TimeLoopConfig) {
        let cfg = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: steps,
            dx: 1.0 / 64.0,
            dy: 1.0 / 64.0,
            iso_config: IsostasyConfig::c1_default(),
            drainage_max_distance: 30,
        };
        let state = init_c1_state_phase_2_r7(64, 3, &Phase2InitParams::default());
        let kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        (state, kin, cfg)
    }

    /// ADR Finding 165-C1 -- the record is read-only: C1's s, age, plate_id, plate_type and cratonic_mask are
    /// bit-identical with and without it, and it keeps the declared snapshots.
    #[test]
    fn the_record_leaves_c1_bit_identical() {
        let (mut a, mut ka, cfg) = c1(40);
        run_with_closures(&mut a, &mut ka, &cfg, &C1Closures::default(), |_, _| {});
        let (mut b, mut kb, _) = c1(40);
        let h = run_c1_recorded(&mut b, &mut kb, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default());
        assert_eq!(a.s.data(), b.s.data());
        assert_eq!(a.age.data(), b.age.data());
        assert_eq!(a.plate_id.data(), b.plate_id.data());
        assert_eq!(a.plate_type.data(), b.plate_type.data());
        assert_eq!(a.cratonic_mask.data(), b.cratonic_mask.data());
        assert_eq!(h.steps, vec![0, 10, 20, 30, 40]);
        assert_eq!(h.raw.len(), 5);
        assert!(h.suture.iter().any(|&s| s), "some cell is convergent at some step");
    }

    /// ADR Finding 167-B5, the negative control -- the record under C1's mapping is F166's bit for bit; the Airy mapping
    /// follows the declaration (S̃ = 1 → +0.45 km, S̃ = 2 → +5.75 km, S̃ = 0.2 → −5.51 km, continuous at sea level) and
    /// reads the oceanic cells at S̃ = 0.2.
    #[test]
    fn the_c1_mapping_is_f166s_record() {
        let (mut a, mut ka, cfg) = c1(40);
        let x = run_c1_sources(&mut a, &mut ka, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default(), &[C1Term::Erosion]);
        let (mut b, mut kb, _) = c1(40);
        let y = run_c1_sources_mapped(&mut b, &mut kb, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default(), &[C1Term::Erosion], Mapping::C1);
        for (p, q) in x.hist.raw.iter().zip(&y.hist.raw) {
            assert_eq!(p.data, q.data);
        }
        let air = Airy::declared();
        assert!((air.altitude_km(1.0) - 0.45).abs() < 1e-9);
        assert!((air.altitude_km(2.0) - 5.753).abs() < 0.01, "{}", air.altitude_km(2.0));
        assert!((air.altitude_km(0.2) + 5.51).abs() < 0.01, "{}", air.altitude_km(0.2));
        let bl = (air.rho_m - air.rho_c) / air.rho_m * air.km_per_s;
        let s0 = (bl - air.anchor_km) / bl;
        assert!(air.altitude_km(s0 - 1e-9).abs() < 1e-6 && air.altitude_km(s0 + 1e-9).abs() < 1e-6, "continuous at sea level");
        let (mut c, mut kc, _) = c1(40);
        let z = run_c1_sources_mapped(&mut c, &mut kc, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default(), &[C1Term::Erosion], Mapping::Airy(Airy { sigma: 0.0, ..air }));
        let oceanic = c.plate_type.data().iter().position(|t| *t == PlateType::Oceanic).unwrap();
        assert!((z.hist.raw.last().unwrap().data[oceanic] as f64 * 5000.0 / 1000.0 - air.altitude_km(0.2)).abs() < 1e-3);
    }

    /// ADR Finding 166-B4, the negative control -- with no term excluded the sources' record is F165's history bit for
    /// bit (the snapshots), so with F165's T the drive is F165's history.
    #[test]
    fn the_sources_with_every_term_are_f165s_history() {
        let (mut a, mut ka, cfg) = c1(40);
        let h = run_c1_recorded(&mut a, &mut ka, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default());
        let (mut b, mut kb, _) = c1(40);
        let s = run_c1_sources(&mut b, &mut kb, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default(), &[]);
        assert_eq!(h.raw.len(), s.hist.raw.len());
        for (x, y) in h.raw.iter().zip(&s.hist.raw) {
            assert_eq!(x.data, y.data);
        }
        assert_eq!(a.s.data(), b.s.data());
        // excluding the erosion changes the record, and only where the erosion acted
        let (mut c, mut kc, _) = c1(40);
        let e = run_c1_sources(&mut c, &mut kc, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default(), &[C1Term::Erosion]);
        assert_eq!(c.s.data(), a.s.data(), "C1 itself is untouched");
        assert_ne!(e.hist.raw.last().unwrap().data, h.raw.last().unwrap().data);
        let ero = &e.terms.iter().find(|t| t.0 == C1Term::Erosion).unwrap().1;
        assert!(ero.iter().all(|&v| v <= 1e-12), "the erosion only removes");
    }

    /// ADR Finding 165-C3, negative control 1 -- the « history » U ∝ the final h_iso, constant, to equilibrium, no
    /// rebound, IS the current level bit for bit (`the_steady_physics_level_is_pinned` pins that level as it was before
    /// the loop was shared).
    #[test]
    fn the_steady_history_is_the_current_level() {
        let (mut s, mut k, cfg) = c1(40);
        let h = run_c1_recorded(&mut s, &mut k, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default());
        let snaps = h.normalized(None);
        let c = CascadeConfig { mfd_exponent: Some(6.0), roughness_m: 40.0, seed: 7, ..CascadeConfig::declared(5000.0) };
        let a = physics_history(&snaps, 128, &c, &HistoryRun { t_years: 0.0, steps: 120, rebound: None, timing: Timing::Steady }, &mut |_| {}, &|| false).unwrap();
        let (b, _) = super::super::physics_level(snaps.last().unwrap(), 128, &c, STEADY_PEAK_M, 120, &mut |_| {}, &|| false).unwrap();
        assert_eq!(a.rec.z.data, b.z.data);
        assert_eq!(a.rec.steps, b.steps);
    }

    /// ADR Finding 165-C2 -- the historical level runs its N_phys steps (no equilibrium stop), its uplift sums to
    /// Δh_iso (final − initial) under both timings, and the rebound raises an eroded field.
    #[test]
    fn the_history_sums_its_uplift_and_the_rebound_lifts() {
        let (mut s, mut k, cfg) = c1(40);
        let h = run_c1_recorded(&mut s, &mut k, &cfg, &C1Closures::default(), 10, &IsostasyConfig::c1_default(), &SteinSteinParams::default());
        let snaps = h.normalized(None);
        let c = CascadeConfig::declared(5000.0);
        let n = 128;
        let want: Vec<f32> = {
            let (a, b) = (upsample_to(&snaps[0], n), upsample_to(snaps.last().unwrap(), n));
            b.data.iter().zip(&a.data).map(|(x, y)| x - y).collect()
        };
        for timing in [Timing::History, Timing::Late10] {
            let run = HistoryRun { t_years: 2.0e6, steps: 40, rebound: None, timing };
            let l = physics_history(&snaps, n, &c, &run, &mut |_| {}, &|| false).unwrap();
            assert_eq!(l.rec.steps, 40);
            let err = l.rec.sum_uplift.iter().zip(&want).map(|(a, b)| (a - b).abs()).fold(0f32, f32::max);
            assert!(err < 1e-4, "{timing:?}: Σ uplift off by {err}");
        }
        let a = physics_history(&snaps, n, &c, &HistoryRun { t_years: 2.0e6, steps: 40, rebound: None, timing: Timing::History }, &mut |_| {}, &|| false).unwrap();
        let b = physics_history(&snaps, n, &c, &HistoryRun { t_years: 2.0e6, steps: 40, rebound: Some(Rebound::declared()), timing: Timing::History }, &mut |_| {}, &|| false).unwrap();
        let mean = |r: &LevelRecord| r.z.data.iter().map(|&v| v as f64).sum::<f64>();
        assert!(mean(&b.rec) > mean(&a.rec), "the rebound lifts");
        // the filter passes the mean (Ĝ(0) = 1): Σ rebound = (ρc / ρm) · Σ erosion
        let (sr, se) = (b.rec.sum_rebound.iter().map(|&v| v as f64).sum::<f64>(), b.rec.sum_erosion.iter().map(|&v| v as f64).sum::<f64>());
        assert!((sr + 2700.0 / 3300.0 * se).abs() < 0.01 * sr.abs().max(1e-9), "Σ rebound {sr} against Σ erosion {se}");
        assert!((Rebound::declared().alpha_km - 64.4).abs() < 0.5, "α {}", Rebound::declared().alpha_km);
    }
}
