//! ADR Finding 168 -- **the cascade's C1 parameter profile** (gated: the cascade and its benches only; C1's
//! production default is unchanged).
//!
//! « Production » is `Phase2InitParams::default()`. A « continent » profile overrides ONE parameter, the initial
//! continental fraction (`ContinentalClusterParams::continental_fraction`, a fraction of the plates: round(8 f)
//! continental plates out of 8). It is a design target, not a physical calibration: L1 and L2 were set once on the
//! témoin for ~40 % and 55–60 % land at the end of the history-driven physics level, then frozen for every seed.
//! Declared in `docs/reports/relief_method/f168_continent_size/f168_declared.md`.

use crate::cascade::history::{Airy, C1Term, Mapping, run_c1_sources_mapped};
use crate::cascade::physio::{Circumnavigation, circumnavigation};
use crate::tectonics::isostasy::IsostasyConfig;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::init_r7::clustering::start_plate_order;
use crate::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7_start};
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};

/// The « continent L1 » initial continental fraction (F168-B1, the témoin's trial nearest 40 % land).
pub const CONTINENT_L1: f64 = 0.375;
/// The « continent L2 » initial continental fraction (F168-B1, the témoin's trial inside 55–60 % land).
pub const CONTINENT_L2: f64 = 0.625;

/// The C1 parameter profile the cascade's history is run with.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum C1Profile {
    /// C1's production default.
    #[default]
    Production,
    /// Production with the initial continental fraction overridden.
    Continent(f64),
    /// ADR Finding 169-A4 -- the continent profile with its BFS grown from a SELECTED start plate ([`select_continent`]).
    Selected { fraction: f64, start: u16 },
}

impl C1Profile {
    pub const L1: Self = Self::Continent(CONTINENT_L1);
    pub const L2: Self = Self::Continent(CONTINENT_L2);

    /// The C1 init parameters of this profile.
    pub fn init_params(self) -> Phase2InitParams {
        let mut p = Phase2InitParams::default();
        match self {
            Self::Production => {}
            Self::Continent(f) | Self::Selected { fraction: f, .. } => p.cluster.continental_fraction = f,
        }
        p
    }

    /// The BFS start plate (`None` = C1's RNG draw).
    pub fn start(self) -> Option<u16> {
        match self {
            Self::Selected { start, .. } => Some(start),
            _ => None,
        }
    }

    /// C1's initial state under this profile.
    pub fn init_state(self, grid: usize, seed: u64) -> C1State {
        init_c1_state_phase_2_r7_start(grid, seed, &self.init_params(), self.start())
    }
}

/// ADR Finding 169-A4 -- the land range (share of the torus) the selection accepts.
pub const SELECT_LAND: (f32, f32) = (0.50, 0.62);
/// ADR Finding 169-A4 -- the fractions tried, in preference order (the larger land first).
pub const SELECT_FRACTIONS: [f64; 2] = [CONTINENT_L2, 0.5];

/// One candidate of the selection.
#[derive(Clone, Copy, Debug)]
pub struct Candidate {
    pub fraction: f64,
    pub start: u16,
    /// Its rank in C1's RNG draw order (0 = today's start plate).
    pub rng_rank: usize,
    /// The end land (64², h_src > 0 under Airy), share of the torus.
    pub land: f32,
    pub circ: Circumnavigation,
}

impl Candidate {
    pub fn qualifies(&self) -> bool {
        (SELECT_LAND.0..=SELECT_LAND.1).contains(&self.land) && !self.circ.wraps && self.circ.sea_r.is_some()
    }
    pub fn profile(&self) -> C1Profile {
        C1Profile::Selected { fraction: self.fraction, start: self.start }
    }
}

/// The selection's record: every candidate and the chosen one.
#[derive(Clone, Debug)]
pub struct Selection {
    pub candidates: Vec<Candidate>,
    pub chosen: Option<usize>,
}

/// The C1 run of the cascade's history (the benches' and the viz's).
pub fn c1_run_config() -> C1TimeLoopConfig {
    C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    }
}

/// ADR Finding 169-A4 -- **the continental plates' selection**, an initial condition: for each fraction of
/// [`SELECT_FRACTIONS`] and each start plate, C1 runs to its end under the Airy record (the sources' drive), and the
/// candidate qualifies if its end land is within [`SELECT_LAND`] and it is circumnavigable (F168's test, 64²). The chosen
/// one has the largest r*; ties go to the earlier fraction, then to C1's RNG draw order. Deterministic from the seed.
pub fn select_continent(seed: u64) -> Selection {
    let cfg = c1_run_config();
    let ss = SteinSteinParams::default();
    let mut candidates = Vec::new();
    for &fraction in &SELECT_FRACTIONS {
        let params = C1Profile::Continent(fraction).init_params();
        for (rng_rank, &start) in start_plate_order(&params.cluster, params.num_plates).iter().enumerate() {
            let mut st = init_c1_state_phase_2_r7_start(64, seed, &params, Some(start));
            let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
            let src = run_c1_sources_mapped(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &ss, &[C1Term::Erosion], Mapping::Airy(Airy::declared()));
            let last = src.hist.normalized(None).pop().unwrap();
            let land: Vec<bool> = last.data.iter().map(|&v| v > 0.5).collect();
            let share = land.iter().filter(|&&b| b).count() as f32 / land.len() as f32;
            candidates.push(Candidate { fraction, start, rng_rank, land: share, circ: circumnavigation(&land, 64) });
        }
    }
    // the candidates are in preference order, so a strict `>` keeps the earlier one on a tie
    let mut chosen: Option<usize> = None;
    for (i, c) in candidates.iter().enumerate() {
        if c.qualifies() && chosen.is_none_or(|j| c.circ.sea_r > candidates[j].circ.sea_r) {
            chosen = Some(i);
        }
    }
    Selection { candidates, chosen }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::isostasy::IsostasyConfig;
    use crate::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
    use crate::tectonics_c1::kinematics::PlateKinematics;
    use crate::tectonics_c1::state::C1State;
    use crate::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use crate::tectonics_v2::boundaries::plate_type::PlateType;

    fn run(p: C1Profile, seed: u64) -> C1State {
        let cfg = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: 30,
            dx: 1.0 / 64.0,
            dy: 1.0 / 64.0,
            iso_config: IsostasyConfig::c1_default(),
            drainage_max_distance: 30,
        };
        let mut st = init_c1_state_phase_2_r7(64, seed, &p.init_params());
        let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
        run_with_closures(&mut st, &mut kin, &cfg, &C1Closures::default(), |_, _| {});
        st
    }

    /// F168-B3, the negative control (permanent): the continent profile at production's fraction (0.29) is C1's
    /// production run bit for bit; L1 and L2 add continental plates.
    #[test]
    fn the_continent_profile_at_production_fraction_is_c1_bit_for_bit() {
        let prod = run(C1Profile::Production, 42);
        let l0 = run(C1Profile::Continent(0.29), 42);
        let bits = |s: &C1State| s.s.data().iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&prod), bits(&l0));
        assert_eq!(prod.plate_id.data(), l0.plate_id.data());
        assert_eq!(prod.plate_type.data(), l0.plate_type.data());
        let cont = |p: C1Profile| {
            let st = init_c1_state_phase_2_r7(64, 42, &p.init_params());
            st.plate_type.data().iter().filter(|t| **t == PlateType::Continental).count()
        };
        let (c0, c1, c2) = (cont(C1Profile::Production), cont(C1Profile::L1), cont(C1Profile::L2));
        assert!(c0 < c1 && c1 < c2, "continental cells {c0} {c1} {c2}");
    }

    /// F169-A4, the negative control (permanent): a start plate equal to C1's first RNG draw is today's init bit for bit,
    /// for the production fraction and for L2.
    #[test]
    fn the_start_plate_at_the_rng_draw_is_c1_bit_for_bit() {
        for seed in [42u64, 9] {
            for p in [C1Profile::Production, C1Profile::L2] {
                let params = p.init_params();
                let first = start_plate_order(&params.cluster, params.num_plates)[0];
                let a = crate::tectonics_c1::init_r7::init_c1_state_phase_2_r7(64, seed, &params);
                let b = init_c1_state_phase_2_r7_start(64, seed, &params, Some(first));
                assert_eq!(a.plate_type.data(), b.plate_type.data());
                assert_eq!(a.s.data().iter().map(|v| v.to_bits()).collect::<Vec<_>>(), b.s.data().iter().map(|v| v.to_bits()).collect::<Vec<_>>());
                assert_eq!(a.cratonic_mask.data(), b.cratonic_mask.data());
            }
        }
    }
}
