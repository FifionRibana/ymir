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
use crate::tectonics_c1::closures::subduction::OoSubductionParams;
use crate::tectonics_c1::time_loop::{C1Closures, C1Ext, C1TimeLoopConfig, run_with_closures_observed, run_with_closures_terms_ext};
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
    run_c1_sources_mapped_oo(state, kin, cfg, closures, every, iso, ss, exclude, mapping, None)
}

/// ADR Finding 175 -- [`run_c1_sources_mapped`] with the ocean-ocean subduction path (`oo`, option 5); `None` is
/// [`run_c1_sources_mapped`] bit for bit.
#[allow(clippy::too_many_arguments)]
pub fn run_c1_sources_mapped_oo(
    state: &mut C1State,
    kin: &mut PlateKinematics,
    cfg: &C1TimeLoopConfig,
    closures: &C1Closures,
    every: usize,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    exclude: &[C1Term],
    mapping: Mapping,
    oo: Option<OoSubductionParams>,
) -> C1Sources {
    run_c1_sources_mapped_flags(state, kin, cfg, closures, every, iso, ss, exclude, mapping, oo, false, false)
}

/// ADR Finding 176 -- [`run_c1_sources_mapped_oo`] with the accretion « merge at collision » flag and the measured-only variant
/// B (`arcs_emerge`: under Airy, an oceanic cell that has received arc mass is read at 0.2 + its cumulative arc mass, the
/// positive part of the subduction term, instead of 0.2). `(false, false)` is [`run_c1_sources_mapped_oo`] bit for bit.
#[allow(clippy::too_many_arguments)]
pub fn run_c1_sources_mapped_flags(
    state: &mut C1State,
    kin: &mut PlateKinematics,
    cfg: &C1TimeLoopConfig,
    closures: &C1Closures,
    every: usize,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    exclude: &[C1Term],
    mapping: Mapping,
    oo: Option<OoSubductionParams>,
    merge_at_collision: bool,
    arcs_emerge: bool,
) -> C1Sources {
    run_c1_sources_mapped_ext(state, kin, cfg, closures, every, iso, ss, exclude, mapping, C1Ext { oo, merge_at_collision, andean_gap: None, andean_arc: None }, arcs_emerge)
}

/// ADR Finding 178 -- [`run_c1_sources_mapped_flags`] with every C1 flag of [`C1Ext`] (option 5, the merge at collision, the
/// Andean seed). `andean_gap: None` is [`run_c1_sources_mapped_flags`] bit for bit.
#[allow(clippy::too_many_arguments)]
pub fn run_c1_sources_mapped_ext(
    state: &mut C1State,
    kin: &mut PlateKinematics,
    cfg: &C1TimeLoopConfig,
    closures: &C1Closures,
    every: usize,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    exclude: &[C1Term],
    mapping: Mapping,
    ext: C1Ext,
    arcs_emerge: bool,
) -> C1Sources {
    let (nx, ny) = (state.nx(), state.ny());
    let n = nx * ny;
    let arcs = RefCell::new(vec![0f64; n]);
    let s0 = state.s.data().to_vec();
    let type0: Vec<_> = state.plate_type.data().to_vec();
    let removed = RefCell::new(vec![0f64; n]);
    let sums = RefCell::new(C1Term::ALL.iter().map(|&t| (t, vec![0f64; n])).collect::<Vec<_>>());
    let src_alt = |st: &C1State, removed: &[f64]| -> GridF32 {
        if let Mapping::Airy(a) = &mapping {
            let s_src = Field2D::from_vec(nx, ny, st.s.data().iter().zip(removed).map(|(s, r)| s - r).collect());
            if arcs_emerge {
                // F176-B: the arc-loaded oceanic cells read at 0.2 + their cumulative arc mass
                let ar = arcs.borrow();
                let mut pt = st.plate_type.clone();
                let mut sv = s_src.data().to_vec();
                for k in 0..n {
                    let (i, j) = (k % nx, k / nx);
                    if pt.get(i, j) == PlateType::Oceanic && ar[k] > 0.0 {
                        pt.set(i, j, PlateType::Continental);
                        sv[k] = 0.2 + ar[k];
                    }
                }
                return airy_raw(&Field2D::from_vec(nx, ny, sv), &pt, a, ss.depth_scale_m);
            }
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
    run_with_closures_terms_ext(
        state,
        kin,
        cfg,
        closures,
        ext,
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
            if arcs_emerge && term == C1Term::Subduction {
                let mut ar = arcs.borrow_mut();
                for k in 0..n {
                    ar[k] += (after[k] - before[k]).max(0.0);
                }
            }
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
    /// ADR Finding 169-C1 -- the continent's age A (years): C1's step k lasts d_N · r^(N−k), d_N = T/N, r solved so the
    /// steps sum to A (r = 1 exactly at A = T). `None` = linear time (F165–F168).
    pub age_years: Option<f64>,
    /// ADR Finding 172-B1 -- the lowstand: the base level at −Δ during the history's last T_q. `None` = the sea at 0.
    pub sea_low: Option<SeaLow>,
}

/// ADR Finding 172-B -- **the Quaternary lowstand**: the base level −`delta_m` during the last `t_q_years` of the history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeaLow {
    pub delta_m: f64,
    pub t_q_years: f64,
}

impl SeaLow {
    /// The declaration (`f172_declared.md`): Δ = 134 m (Lambeck et al. 2014, PNAS 111, 15296: « a slow fall to −134 m from
    /// 29 to 21 ka BP »), T_q = 2.58 Myr (the base of the Quaternary, ICS).
    pub fn declared() -> Self {
        Self { delta_m: 134.0, t_q_years: 2.58e6 }
    }
    /// Per step: is it in the lowstand (its remaining time, the sum of the later steps' durations, ≤ T_q)?
    pub fn steps_low(&self, durations: &[f64]) -> Vec<bool> {
        let mut rem = 0.0;
        let mut out = vec![false; durations.len()];
        for i in (0..durations.len()).rev() {
            rem += durations[i];
            out[i] = rem <= self.t_q_years;
        }
        out
    }
}

/// ADR Finding 173-B -- **a sea-level curve** s(age) in m above present: a published series (Spratt & Lisiecki 2016's
/// `SeaLev_longPC1`, 0–798 ka every 1 kyr). Between the series' oldest age and `quaternary_kyr` the level is `before_m`
/// (0 = the declaration); older, 0. ADR Finding 174-0: the parsed series is re-zeroed so it reads 0 m today (the stack reads
/// +8.96 m at 0 ka) -- an alignment of reference, not a calibration.
#[derive(Clone, Debug, PartialEq)]
pub struct SeaCurve {
    pub age_kyr: Vec<f64>,
    pub level_m: Vec<f64>,
    pub before_m: f64,
    pub quaternary_kyr: f64,
    /// The sub-step (kyr) of the C1 steps overlapping the series.
    pub sub_kyr: f64,
}

impl SeaCurve {
    /// Parse NOAA's `spratt2016.txt` (study 19982): the `SeaLev_longPC1` column against `age_calkaBP`, the NaN rows dropped,
    /// every level minus the youngest row's (F174-0: 0 m today).
    pub fn from_noaa_spratt2016(text: &str) -> Result<Self, String> {
        let mut col: Option<usize> = None;
        let (mut age, mut lev) = (Vec::new(), Vec::new());
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
            if f.first() == Some(&"age_calkaBP") {
                col = f.iter().position(|s| *s == "SeaLev_longPC1");
                continue;
            }
            let c = col.ok_or("no header with SeaLev_longPC1")?;
            let (Some(a), Some(v)) = (f.first(), f.get(c)) else { continue };
            let (Ok(a), Ok(v)) = (a.parse::<f64>(), v.parse::<f64>()) else { continue };
            if v.is_finite() {
                age.push(a);
                lev.push(v);
            }
        }
        if age.len() < 2 {
            return Err("too few rows".into());
        }
        let now = lev[0];
        for v in &mut lev {
            *v -= now;
        }
        Ok(Self { age_kyr: age, level_m: lev, before_m: 0.0, quaternary_kyr: 2580.0, sub_kyr: 2.0 })
    }
    /// s at `age_kyr` (linear between the series' points).
    pub fn level_at(&self, age_kyr: f64) -> f64 {
        let last = *self.age_kyr.last().unwrap();
        if age_kyr > last {
            return if age_kyr <= self.quaternary_kyr { self.before_m } else { 0.0 };
        }
        let i = self.age_kyr.partition_point(|&a| a <= age_kyr).clamp(1, self.age_kyr.len() - 1);
        let (a0, a1) = (self.age_kyr[i - 1], self.age_kyr[i]);
        let w = if a1 > a0 { (age_kyr - a0) / (a1 - a0) } else { 0.0 };
        self.level_m[i - 1] * (1.0 - w) + self.level_m[i] * w
    }
    /// The series' mean (m).
    pub fn mean(&self) -> f64 {
        self.level_m.iter().sum::<f64>() / self.level_m.len() as f64
    }
}

/// ADR Finding 169-C1 -- the durations (years) of `steps` C1 steps whose last lasts `t_years / steps` and which sum to
/// `age_years`, growing geometrically into the past; and r. At `age_years == t_years`, r = 1 exactly.
pub fn step_durations(t_years: f64, age_years: f64, steps: usize) -> (Vec<f64>, f64) {
    let d_n = t_years / steps as f64;
    let total = |r: f64| (1..=steps).map(|k| d_n * r.powi((steps - k) as i32)).sum::<f64>();
    let r = if age_years == t_years {
        1.0
    } else {
        assert!(age_years > t_years, "the age A must be at least T");
        let (mut lo, mut hi) = (1.0f64, 2.0f64);
        while total(hi) < age_years {
            hi *= 2.0;
        }
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            if total(m) > age_years { hi = m } else { lo = m }
        }
        0.5 * (lo + hi)
    };
    ((1..=steps).map(|k| d_n * r.powi((steps - k) as i32)).collect(), r)
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
    physics_history_rock(snaps, n, cfg, run, None, progress, cancel)
}

/// ADR Finding 171-B/C -- [`physics_history`] with a rock law at the level (`None` = [`physics_history`] exactly).
pub fn physics_history_rock(
    snaps: &[GridF32],
    n: usize,
    cfg: &CascadeConfig,
    run: &HistoryRun,
    rock: Option<&super::Rock>,
    progress: &mut dyn FnMut(CascadeProgress),
    cancel: &dyn Fn() -> bool,
) -> Option<HistoryLevel> {
    physics_history_sea(snaps, n, cfg, run, rock, None, progress, cancel)
}

/// ADR Finding 173-B -- [`physics_history_rock`] with the base level following a sea-level curve (`None` =
/// [`physics_history_rock`] exactly): the C1 steps overlapping the curve's series are cut into ⌈d / sub⌉ equal sub-steps, each
/// with its step's uplift spread uniformly and its base level s at its mid-age. `run.sea_low` is ignored when a curve is given.
#[allow(clippy::too_many_arguments)]
pub fn physics_history_sea(
    snaps: &[GridF32],
    n: usize,
    cfg: &CascadeConfig,
    run: &HistoryRun,
    rock: Option<&super::Rock>,
    curve: Option<&SeaCurve>,
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
        let tor = amplify::torus_at(cfg.torus, n);
        let g = amplify::noise01_t(n, cell_km, 8.0 * cell_km, 2, cfg.seed, "cascade_phys_rough", tor);
        // F169-A1: summed in terrain order
        let land: Vec<usize> = amplify::terrain_order(n, tor).into_iter().filter(|&k| z0.data[k] > cfg.sea_level).collect();
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
    // ADR Finding 173-B: the sea-level curve -- the sub-step plan (C1 step, the uplift's share, dt, the base level)
    if let (Some(cv), Timing::History) = (curve, timing) {
        let d: Vec<f64> = match run.age_years {
            Some(a) => step_durations(run.t_years, a, steps).0,
            None => vec![run.t_years / steps as f64; steps],
        };
        let interval = move |i: usize| ((i * intervals) / steps).min(intervals - 1);
        let mut dsum = vec![0f64; intervals];
        for (i, di) in d.iter().enumerate() {
            dsum[interval(i)] += di;
        }
        let mut rem = vec![0f64; steps];
        let mut acc = 0f64;
        for i in (0..steps).rev() {
            rem[i] = acc;
            acc += d[i];
        }
        let window_yr = cv.age_kyr.last().copied().unwrap_or(0.0) * 1000.0;
        let mut plan: Vec<(usize, f32, f32, f32)> = Vec::new();
        for i in 0..steps {
            let ns = if rem[i] < window_yr { (d[i] / (cv.sub_kyr * 1000.0)).ceil().max(1.0) as usize } else { 1 };
            let share = match run.age_years {
                Some(_) => d[i] / dsum[interval(i)],
                None => intervals as f64 / steps as f64,
            } / ns as f64;
            let dt = d[i] / ns as f64;
            for k in 0..ns {
                let mid_age_kyr = (rem[i] + d[i] - (k as f64 + 0.5) * dt) / 1000.0;
                let sea = cfg.sea_level + (cv.level_at(mid_age_kyr) as f32) / n2m;
                plan.push((i, share as f32, dt as f32, sea));
            }
        }
        let delta: Vec<Vec<f32>> = (0..intervals).map(|j| up[j + 1].data.iter().zip(&up[j].data).map(|(b, a)| b - a).collect()).collect();
        let total = plan.len();
        let plan = std::rc::Rc::new(plan);
        let (p1, p2, p3) = (plan.clone(), plan.clone(), plan.clone());
        let drive = Drive {
            du: Box::new(move |s| {
                let (i, f, _, _) = p1[s.min(total - 1)];
                std::borrow::Cow::Owned(delta[interval(i)].iter().map(|v| v * f).collect())
            }),
            hold_sea: false,
            eq_stop: false,
            rebound: run.rebound,
            dt: Some(Box::new(move |s| p2[s.min(total - 1)].2)),
            rock,
            sea: Some(Box::new(move |s| p3[s.min(total - 1)].3)),
        };
        let uplift = GridF32::new(n, n, 0.0);
        let rec = run_loop(z0.clone(), uplift, &cfg_t, n, total, Some(0), &drive, progress, cancel)?;
        return Some(HistoryLevel { rec, z0, snapshots: up });
    }
    // ADR Finding 172-B1: the lowstand's base level per step (norm)
    let sea_fn = |durations: &[f64]| -> Option<Box<dyn Fn(usize) -> f32>> {
        run.sea_low.map(|sl| {
            let low = sl.steps_low(durations);
            let (hi_s, lo_s) = (cfg.sea_level, cfg.sea_level - (sl.delta_m as f32) / n2m);
            Box::new(move |i: usize| if low[i.min(low.len() - 1)] { lo_s } else { hi_s }) as Box<dyn Fn(usize) -> f32>
        })
    };
    // ADR Finding 169-C1: the continent's age -- each step its own duration, the interval's Δh spread in proportion
    if let (Some(age), Timing::History) = (run.age_years, timing) {
        let (d, _) = step_durations(run.t_years, age, steps);
        let sea = sea_fn(&d);
        let interval = move |i: usize| ((i * intervals) / steps).min(intervals - 1);
        let mut dsum = vec![0f64; intervals];
        for (i, di) in d.iter().enumerate() {
            dsum[interval(i)] += di;
        }
        let delta: Vec<Vec<f32>> = (0..intervals).map(|j| up[j + 1].data.iter().zip(&up[j].data).map(|(b, a)| b - a).collect()).collect();
        let dd = d.clone();
        let drive = Drive {
            du: Box::new(move |i| {
                let j = interval(i);
                let f = (dd[i] / dsum[j]) as f32;
                std::borrow::Cow::Owned(delta[j].iter().map(|v| v * f).collect())
            }),
            hold_sea: false,
            eq_stop: false,
            rebound: run.rebound,
            dt: Some(Box::new(move |i| d[i.min(steps - 1)] as f32)),
            rock,
            sea,
        };
        let uplift = GridF32::new(n, n, 0.0);
        let rec = run_loop(z0.clone(), uplift, &cfg_t, n, run.steps, Some(0), &drive, progress, cancel)?;
        return Some(HistoryLevel { rec, z0, snapshots: up });
    }
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
        dt: None,
        rock,
        sea: sea_fn(&vec![run.t_years / steps as f64; steps]),
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
        let a = physics_history(&snaps, 128, &c, &HistoryRun { t_years: 0.0, steps: 120, rebound: None, timing: Timing::Steady, age_years: None, sea_low: None }, &mut |_| {}, &|| false).unwrap();
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
            let run = HistoryRun { t_years: 2.0e6, steps: 40, rebound: None, timing, age_years: None, sea_low: None };
            let l = physics_history(&snaps, n, &c, &run, &mut |_| {}, &|| false).unwrap();
            assert_eq!(l.rec.steps, 40);
            let err = l.rec.sum_uplift.iter().zip(&want).map(|(a, b)| (a - b).abs()).fold(0f32, f32::max);
            assert!(err < 1e-4, "{timing:?}: Σ uplift off by {err}");
        }
        let a = physics_history(&snaps, n, &c, &HistoryRun { t_years: 2.0e6, steps: 40, rebound: None, timing: Timing::History, age_years: None, sea_low: None }, &mut |_| {}, &|| false).unwrap();
        let b = physics_history(&snaps, n, &c, &HistoryRun { t_years: 2.0e6, steps: 40, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: None, sea_low: None }, &mut |_| {}, &|| false).unwrap();
        let mean = |r: &LevelRecord| r.z.data.iter().map(|&v| v as f64).sum::<f64>();
        assert!(mean(&b.rec) > mean(&a.rec), "the rebound lifts");
        // the filter passes the mean (Ĝ(0) = 1): Σ rebound = (ρc / ρm) · Σ erosion
        let (sr, se) = (b.rec.sum_rebound.iter().map(|&v| v as f64).sum::<f64>(), b.rec.sum_erosion.iter().map(|&v| v as f64).sum::<f64>());
        assert!((sr + 2700.0 / 3300.0 * se).abs() < 0.01 * sr.abs().max(1e-9), "Σ rebound {sr} against Σ erosion {se}");
        assert!((Rebound::declared().alpha_km - 64.4).abs() < 0.5, "α {}", Rebound::declared().alpha_km);
    }
}

/// ADR Finding 169 -- the torus, the terrain frame and the age, tested (permanent).
#[cfg(test)]
mod f169_tests {
    use super::*;
    use crate::cascade::amplify::{AmpConfig, Chain, Variant, Warp, WarpPlace, noise01_t};
    use crate::cascade::{CascadeConfig, roll};

    /// 31 synthetic 64² snapshots (norm, sea 0.5): an island straddling both seams of the torus, rising in time, with a
    /// ridge that grows late.
    fn snaps() -> Vec<GridF32> {
        (0..31)
            .map(|s| {
                let a = s as f32 / 30.0;
                let mut g = GridF32::new(64, 64, 0.44);
                for y in 0..64 {
                    for x in 0..64 {
                        // periodic distance to (60, 3)
                        let dx = ((x as f32 - 60.0 + 32.0).rem_euclid(64.0)) - 32.0;
                        let dy = ((y as f32 - 3.0 + 32.0).rem_euclid(64.0)) - 32.0;
                        let r2 = (dx * dx + 1.6 * dy * dy) / 140.0;
                        let ridge = (-((dx - 0.6 * dy) / 2.5).powi(2)).exp() * (-r2 / 3.0).exp();
                        g.data[y * 64 + x] = 0.44 + (0.06 + 0.02 * a) * (-r2).exp() + 0.05 * a * a * ridge;
                    }
                }
                g
            })
            .collect()
    }

    fn cfg(origin: [usize; 2]) -> CascadeConfig {
        CascadeConfig {
            mfd_exponent: Some(6.0),
            roughness_m: 40.0,
            receiver_tau: Some(0.5),
            seed: 11,
            torus: Some(origin),
            ..CascadeConfig::declared(5000.0)
        }
    }

    fn run(snaps: &[GridF32], origin: [usize; 2], age: Option<f64>) -> (GridF32, GridF32) {
        let run = HistoryRun { t_years: 2.0e6, steps: 30, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: None };
        let l = physics_history(snaps, 128, &cfg(origin), &run, &mut |_| {}, &|| false).unwrap();
        let n2m = cfg(origin).norm_to_m();
        let zm = GridF32 { width: 128, height: 128, data: l.rec.z.data.iter().map(|v| (v - 0.5) * n2m).collect() };
        let amp = AmpConfig {
            budget_scale: 0.01,
            n1_recal_depth: 2,
            n1_pi_m: vec![(256, 3.0)],
            warp: Some(Warp { amp: 0.5, corr: 4.0, place: WarpPlace::All }),
            torus: Some(origin),
            ..AmpConfig::declared(11, Variant::N1, 1.0)
        };
        let mut ch = Chain::new(zm.clone(), amp);
        let out = ch.next_level(&|| false).unwrap().result();
        (zm, out)
    }

    /// F169-A1, the invariance test (declared before measuring: bit for bit): the input rolled by (dx, dy) 64² cells,
    /// with that origin, gives the reference output rolled by (dx, dy)·n/64 at the physics level (128²) and at the
    /// amplification level (256²) -- shifts on one axis, the other, and both edges.
    #[test]
    fn the_torus_cascade_is_translation_invariant_bit_for_bit() {
        let base = snaps();
        let (p0, a0) = run(&base, [0, 0], None);
        for d in [[5usize, 0usize], [0, 23], [37, 50]] {
            let rolled: Vec<GridF32> = base.iter().map(|g| roll(g, d)).collect();
            let (p, a) = run(&rolled, d, None);
            assert_eq!(p.data, roll(&p0, [2 * d[0], 2 * d[1]]).data, "physics level, shift {d:?}");
            assert_eq!(a.data, roll(&a0, [4 * d[0], 4 * d[1]]).data, "amplification level, shift {d:?}");
        }
    }

    /// F169-C1 (permanent): A = T is linear time bit for bit (r = 1), and an older continent is not.
    #[test]
    fn the_age_at_t_is_linear_time_bit_for_bit() {
        let base = snaps();
        let (p_lin, _) = run(&base, [0, 0], None);
        let (p_t, _) = run(&base, [0, 0], Some(2.0e6));
        assert_eq!(p_lin.data, p_t.data);
        let (p_old, _) = run(&base, [0, 0], Some(2.0e7));
        assert_ne!(p_lin.data, p_old.data);
        let (d, r) = step_durations(5.623e6, 3.0e8, 300);
        assert!((d.iter().sum::<f64>() / 3.0e8 - 1.0).abs() < 1e-9 && (d[299] - 5.623e6 / 300.0).abs() < 1e-6, "r {r}");
        assert!((r - 1.019306).abs() < 1e-5, "r {r}");
        assert_eq!(step_durations(5.623e6, 5.623e6, 300).1, 1.0);
    }

    /// ADR Finding 171-B4 and C4, the negative controls (permanent): a rock law with ρ_K = 1 is F170 bit for bit, and a
    /// heterogeneous K of σ_lnK = 0 is K1 bit for bit, at the physics level and at an amplification level; with ρ_K = 10
    /// the cover's exhaustion changes the field.
    #[test]
    fn the_rock_law_reduces_to_f170_and_k1_bit_for_bit() {
        use crate::cascade::Rock;
        use crate::cascade::amplify::KHetero;
        let base = snaps();
        let run = HistoryRun { t_years: 2.0e6, steps: 30, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: Some(2.0e7), sea_low: None };
        let c = cfg([0, 0]);
        let ref_ = physics_history(&base, 128, &c, &run, &mut |_| {}, &|| false).unwrap();
        let cover = vec![50.0f32; 128 * 128];
        let r1 = Rock { rho_k: 1.0, cover_m: cover.clone(), ln_k: None };
        let a = physics_history_rock(&base, 128, &c, &run, Some(&r1), &mut |_| {}, &|| false).unwrap();
        assert_eq!(a.rec.z.data, ref_.rec.z.data, "rho_k = 1");
        let r10 = Rock { rho_k: 10.0, cover_m: cover.clone(), ln_k: None };
        let b = physics_history_rock(&base, 128, &c, &run, Some(&r10), &mut |_| {}, &|| false).unwrap();
        assert_ne!(b.rec.z.data, ref_.rec.z.data, "rho_k = 10 changes the field");
        assert!(b.rec.e_cum_m.iter().any(|&e| e > 50.0), "the cover is exhausted somewhere");
        let k0 = KHetero { sigma: 0.0, ..KHetero::declared() };
        let r10z = Rock { ln_k: Some(k0.band(128, 400.0 / 128.0, 0, k0.resolved(400.0 / 128.0), 11, Some([0, 0]))), ..r10.clone() };
        let bz = physics_history_rock(&base, 128, &c, &run, Some(&r10z), &mut |_| {}, &|| false).unwrap();
        assert_eq!(bz.rec.z.data, b.rec.z.data, "sigma = 0 at the physics level");
        // the amplification level: σ = 0 is the one-k chain
        let n2m = c.norm_to_m();
        let zm = GridF32 { width: 128, height: 128, data: b.rec.z.data.iter().map(|v| (v - 0.5) * n2m).collect() };
        let amp = |kh: Option<KHetero>| AmpConfig { budget_scale: 0.01, n1_recal_depth: 2, warp: Some(Warp { amp: 0.5, corr: 4.0, place: WarpPlace::All }), torus: Some([0, 0]), k_hetero: kh, ..AmpConfig::declared(11, Variant::N1, 1.0) };
        let one = Chain::new(zm.clone(), amp(None)).next_level(&|| false).unwrap().result();
        let mut ch = Chain::new(zm.clone(), amp(Some(k0))).with_ln_k(vec![0.0; 128 * 128]);
        let zero = ch.next_level(&|| false).unwrap().result();
        assert_eq!(one.data, zero.data, "sigma = 0 at the amplification level");
        let mut ch = Chain::new(zm, amp(Some(KHetero::declared()))).with_ln_k(KHetero::declared().band(128, 400.0 / 128.0, 0, 3, 11, Some([0, 0])));
        let het = ch.next_level(&|| false).unwrap();
        assert_ne!(one.data, het.result().data, "sigma = 0.5 changes the level");
        let l = het.ln_k.as_ref().unwrap();
        let sd = (l.iter().map(|v| (*v as f64).powi(2)).sum::<f64>() / l.len() as f64).sqrt();
        assert!(sd > 0.3 && sd < 0.6, "ln M's spread {sd}");
    }

    /// ADR Finding 172-B3, the negative control (permanent): a lowstand of Δ = 0 is F171's K2 bit for bit; Δ = 134 m
    /// changes the field; the lowstand covers the steps whose remaining time is ≤ T_q.
    #[test]
    fn the_lowstand_at_zero_is_f171_bit_for_bit() {
        use crate::cascade::Rock;
        use crate::cascade::amplify::KHetero;
        let base = snaps();
        let c = cfg([0, 0]);
        let kh = KHetero::declared();
        let rock = Rock { rho_k: 10.0, cover_m: vec![80.0; 128 * 128], ln_k: Some(kh.band(128, 400.0 / 128.0, 0, kh.resolved(400.0 / 128.0), 11, Some([0, 0]))) };
        for age in [None, Some(2.0e7)] {
            let run = |sl: Option<SeaLow>| HistoryRun { t_years: 2.0e6, steps: 30, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: sl };
            let a = physics_history_rock(&base, 128, &c, &run(None), Some(&rock), &mut |_| {}, &|| false).unwrap();
            let z = physics_history_rock(&base, 128, &c, &run(Some(SeaLow { delta_m: 0.0, t_q_years: 1.0e6 })), Some(&rock), &mut |_| {}, &|| false).unwrap();
            assert_eq!(a.rec.z.data, z.rec.z.data, "delta = 0, age {age:?}");
            let l = physics_history_rock(&base, 128, &c, &run(Some(SeaLow { delta_m: 134.0, t_q_years: 1.0e6 })), Some(&rock), &mut |_| {}, &|| false).unwrap();
            assert_ne!(a.rec.z.data, l.rec.z.data, "delta = 134 m changes the field");
        }
        let low = SeaLow { delta_m: 134.0, t_q_years: 2.58e6 }.steps_low(&vec![5.623e6 / 300.0; 300]);
        assert_eq!(low.iter().filter(|&&b| b).count(), 137, "the last 2.58 Myr of 5.623 Myr in 300 steps");
        assert!(low[299] && !low[0]);
    }

    /// ADR Finding 173-B4 (permanent): the constant lowstand through the new entry point (no curve) is F172 bit for bit; a
    /// curve parses, interpolates, falls back to `before_m` then 0; a non-zero curve changes the field and adds sub-steps.
    #[test]
    fn the_sea_curve_entry_point_keeps_f172_and_sub_steps() {
        use crate::cascade::Rock;
        let base = snaps();
        let c = cfg([0, 0]);
        let rock = Rock { rho_k: 10.0, cover_m: vec![80.0; 128 * 128], ln_k: None };
        let run = HistoryRun { t_years: 2.0e6, steps: 30, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: None, sea_low: Some(SeaLow { delta_m: 134.0, t_q_years: 1.0e6 }) };
        let a = physics_history_rock(&base, 128, &c, &run, Some(&rock), &mut |_| {}, &|| false).unwrap();
        let b = physics_history_sea(&base, 128, &c, &run, Some(&rock), None, &mut |_| {}, &|| false).unwrap();
        assert_eq!(a.rec.z.data, b.rec.z.data, "the constant lowstand, through the new entry point");
        let text = "# x
age_calkaBP	SeaLev_shortPC1	SeaLev_longPC1
0	1	10
1	NaN	-20
2	NaN	NaN
3	0	-40
";
        let cv = SeaCurve::from_noaa_spratt2016(text).unwrap();
        assert_eq!(cv.age_kyr, vec![0.0, 1.0, 3.0]);
        // F174-0: re-zeroed on the youngest row (10 m): 0 today, every level shifted by the same −10 m
        assert_eq!(cv.level_m, vec![0.0, -30.0, -50.0]);
        assert_eq!(cv.level_at(0.0), 0.0);
        assert!((cv.level_at(0.5) + 15.0).abs() < 1e-9 && (cv.level_at(2.0) + 40.0).abs() < 1e-9);
        assert_eq!(cv.level_at(100.0), 0.0);
        let cv2 = SeaCurve { before_m: -50.0, ..cv.clone() };
        assert_eq!(cv2.level_at(100.0), -50.0);
        assert_eq!(cv2.level_at(3000.0), 0.0);
        // a curve over the last 200 kyr of a 2 Myr history (30 steps of 66.7 kyr): the 3 last steps cut into 34 sub-steps each
        let big = SeaCurve { age_kyr: vec![0.0, 200.0], level_m: vec![-100.0, -100.0], before_m: 0.0, quaternary_kyr: 2580.0, sub_kyr: 2.0 };
        let run0 = HistoryRun { sea_low: None, ..run };
        let plain = physics_history_rock(&base, 128, &c, &run0, Some(&rock), &mut |_| {}, &|| false).unwrap();
        let curved = physics_history_sea(&base, 128, &c, &run0, Some(&rock), Some(&big), &mut |_| {}, &|| false).unwrap();
        assert_eq!(curved.rec.steps, 27 + 3 * 34, "the sub-steps (the steps ending 0, 66.7 and 133.3 kyr before the end)");
        assert_ne!(plain.rec.z.data, curved.rec.z.data);
    }

    /// F169-A1: the torus noise follows the terrain (an origin rolls it exactly) and wraps (no seam).
    #[test]
    fn the_torus_noise_follows_the_terrain_and_wraps() {
        let n = 64;
        let a = noise01_t(n, 6.25, 25.0, 3, 5, "t", Some([0, 0]));
        let b = noise01_t(n, 6.25, 25.0, 3, 5, "t", Some([9, 40]));
        let ga = GridF32 { width: n, height: n, data: a.clone() };
        assert_eq!(b, roll(&ga, [9, 40]).data);
        // across the seam the step is like any neighbour step
        let step = |x0: usize, x1: usize| (0..n).map(|y| (a[y * n + x0] - a[y * n + x1]).abs()).sum::<f32>() / n as f32;
        let inner: f32 = (0..n - 1).map(|x| step(x, x + 1)).sum::<f32>() / (n - 1) as f32;
        assert!(step(n - 1, 0) < 3.0 * inner, "seam {} vs {}", step(n - 1, 0), inner);
    }

    /// F169-A1, a diagnostic kept for whoever breaks the invariance: compares each amplification stage of the reference
    /// and of a rolled input (the first differing stage names the unanchored source; it found the warp's f32 sampling).
    #[test]
    #[ignore]
    fn f169_debug_amp_stages() {
        let base = snaps();
        let stage = |sn: &[GridF32], o: [usize; 2]| {
            let (zm, _) = run(sn, o, None);
            let amp = AmpConfig { budget_scale: 0.01, n1_recal_depth: 2, n1_pi_m: vec![(256, 3.0)], warp: Some(Warp { amp: 0.5, corr: 4.0, place: WarpPlace::All }), torus: Some(o), ..AmpConfig::declared(11, Variant::N1, 1.0) };
            let mut ch = Chain::new(zm, amp.clone());
            let b = amp.budget(256);
            let mut v: Vec<(&str, Vec<f32>)> = Vec::new();
            let l = ch.start_level();
            v.push(("upscaled", l.upscaled.data.clone()));
            v.push(("ocean", l.ocean.iter().map(|&b| b as u8 as f32).collect()));
            v.push(("d5", l.d5.iter().map(|&b| b as u8 as f32).collect()));
            v.push(("area", l.area.clone()));
            v.push(("s0", l.s0.clone()));
            v.push(("pi", l.pi.clone()));
            v.push(("z0", l.z.data.clone()));
            v.push(("prev_work", crate::cascade::upsample2(&l.prev_work).data));
            v.push(("prev2", l.prev2_work.as_ref().map(|g| crate::cascade::upsample_to(g, 256).data).unwrap_or_default()));
            l.erode(b.erosion, &amp, &|| false);
            v.push(("erode", l.z.data.clone()));
            l.talus(b.talus, &amp, &|| false);
            v.push(("talus", l.z.data.clone()));
            l.deposit(b.deposit, &amp, &|| false);
            v.push(("deposit", l.z.data.clone()));
            v
        };
        let a = stage(&base, [0, 0]);
        let d = [5usize, 0usize];
        let rolled: Vec<GridF32> = base.iter().map(|g| roll(g, d)).collect();
        let b = stage(&rolled, d);
        for ((name, x), (_, y)) in a.iter().zip(&b) {
            if x.is_empty() { continue; }
            let xr = roll(&GridF32 { width: 256, height: 256, data: x.clone() }, [4 * d[0], 4 * d[1]]).data;
            let diff: Vec<usize> = (0..xr.len()).filter(|&k| xr[k].to_bits() != y[k].to_bits()).collect();
            let maxd = diff.iter().map(|&k| (xr[k] - y[k]).abs()).fold(0f32, f32::max);
            eprintln!("{name}: {} cells differ, max {maxd}, first {:?}", diff.len(), diff.first().map(|&k| (k % 256, k / 256)));
        }
    }
}
