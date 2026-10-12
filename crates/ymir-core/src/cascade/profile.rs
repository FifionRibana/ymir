//! ADR Finding 168 -- **the cascade's C1 parameter profile** (gated: the cascade and its benches only; C1's
//! production default is unchanged).
//!
//! « Production » is `Phase2InitParams::default()`. A « continent » profile overrides ONE parameter, the initial
//! continental fraction (`ContinentalClusterParams::continental_fraction`, a fraction of the plates: round(8 f)
//! continental plates out of 8). It is a design target, not a physical calibration: L1 and L2 were set once on the
//! témoin for ~40 % and 55–60 % land at the end of the history-driven physics level, then frozen for every seed.
//! Declared in `docs/reports/relief_method/f168_continent_size/f168_declared.md`.

use crate::cascade::history::{Airy, C1Term, Mapping, run_c1_sources_mapped};
use crate::cascade::margin::{self, Margin};
use crate::cascade::physio::{Circumnavigation, circumnavigation};
use crate::tectonics::isostasy::IsostasyConfig;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::init_r7::clustering::start_plate_order;
use crate::tectonics_c1::init_r7::{ClusterGrowth, Phase2InitParams, init_c1_state_phase_2_r7_growth, init_c1_state_phase_2_r7_start};
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
    /// ADR Finding 169-A4 / 170-A1 -- a SELECTED continent ([`select_continent`], [`select_continent_plates`]): `plates`
    /// Voronoi plates (production's 8, or more), `k` of them continental, grown from `start` by C1's BFS or by the
    /// compact family.
    Selected { plates: usize, k: usize, start: u16, compact: bool },
}

impl C1Profile {
    pub const L1: Self = Self::Continent(CONTINENT_L1);
    pub const L2: Self = Self::Continent(CONTINENT_L2);

    /// The C1 init parameters of this profile.
    pub fn init_params(self) -> Phase2InitParams {
        let mut p = Phase2InitParams::default();
        match self {
            Self::Production => {}
            Self::Continent(f) => p.cluster.continental_fraction = f,
            Self::Selected { plates, k, .. } => {
                p.num_plates = plates;
                p.cluster.continental_fraction = k as f64 / plates as f64;
            }
        }
        p
    }

    /// How the continental cluster grows.
    pub fn growth(self) -> ClusterGrowth {
        match self {
            Self::Selected { start, compact: false, .. } => ClusterGrowth::Bfs(start),
            Self::Selected { start, compact: true, .. } => ClusterGrowth::Compact(start),
            _ => ClusterGrowth::Rng,
        }
    }

    /// C1's initial state under this profile.
    pub fn init_state(self, grid: usize, seed: u64) -> C1State {
        init_c1_state_phase_2_r7_growth(grid, seed, &self.init_params(), self.growth())
    }

    /// ADR Finding 174-B -- C1's initial state under this profile with the margin applied ([`margin::apply`]);
    /// [`Margin::off`] is [`Self::init_state`] bit for bit.
    pub fn init_state_margin(self, grid: usize, seed: u64, m: &Margin) -> C1State {
        let mut st = self.init_state(grid, seed);
        margin::apply(&mut st, m, &self.init_params(), seed);
        st
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
    /// F170: the plate count, the continental count, the family.
    pub plates: usize,
    pub k: usize,
    pub compact: bool,
    pub start: u16,
    /// Its rank in C1's RNG draw order (0 = today's start plate).
    pub rng_rank: usize,
    /// The end land (64², h_src > 0 under Airy), share of the torus.
    pub land: f32,
    pub circ: Circumnavigation,
}

impl Candidate {
    pub fn qualifies(&self) -> bool {
        self.qualifies_min(0)
    }
    /// F170-A2: within the land range and circumnavigable with r* ≥ `min_r`.
    pub fn qualifies_min(&self, min_r: usize) -> bool {
        (SELECT_LAND.0..=SELECT_LAND.1).contains(&self.land) && !self.circ.wraps && self.circ.sea_r.is_some_and(|r| r >= min_r)
    }
    pub fn profile(&self) -> C1Profile {
        C1Profile::Selected { plates: self.plates, k: self.k, start: self.start, compact: self.compact }
    }
}

/// The selection's record: every candidate and the chosen one.
#[derive(Clone, Debug)]
pub struct Selection {
    pub candidates: Vec<Candidate>,
    pub chosen: Option<usize>,
    /// The distinct C1 runs (F170: candidates with the same plate set share one).
    pub runs: usize,
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
            let k = ((params.num_plates as f64 * fraction).round() as usize).clamp(1, params.num_plates);
            candidates.push(Candidate { fraction, plates: params.num_plates, k, compact: false, start, rng_rank, land: share, circ: circumnavigation(&land, 64) });
        }
    }
    // the candidates are in preference order, so a strict `>` keeps the earlier one on a tie
    let mut chosen: Option<usize> = None;
    for (i, c) in candidates.iter().enumerate() {
        if c.qualifies() && chosen.is_none_or(|j| c.circ.sea_r > candidates[j].circ.sea_r) {
            chosen = Some(i);
        }
    }
    let runs = candidates.len();
    Selection { candidates, chosen, runs }
}

/// ADR Finding 171-A1 -- the continent profile's plate counts: 16 by default, then the fallbacks.
pub const PLATES_FALLBACK: [usize; 3] = [16, 20, 24];

/// F170-A1: the continental counts tried at `plates`: ⌈0.35 N⌉ … ⌊0.75 N⌋.
pub fn k_range(plates: usize) -> Vec<usize> {
    let (lo, hi) = ((0.35 * plates as f64).ceil() as usize, (0.75 * plates as f64).floor() as usize);
    (lo..=hi).collect()
}

/// ADR Finding 171-A1 -- **the default selection**: 16 plates, then 20, then 24 if nothing qualifies (r* ≥ 1). Returns the
/// retained plate count (`None` when even 24 fails) and every selection tried.
pub fn select_continent_default(seed: u64) -> (Option<usize>, Vec<(usize, Selection)>) {
    select_continent_default_margin(seed, &Margin::off())
}

/// ADR Finding 174-B1 -- [`select_continent_default`] with every candidate's C1 run starting from the margined state;
/// [`Margin::off`] is [`select_continent_default`] bit for bit.
pub fn select_continent_default_margin(seed: u64, m: &Margin) -> (Option<usize>, Vec<(usize, Selection)>) {
    let mut tried = Vec::new();
    for &n in &PLATES_FALLBACK {
        let s = select_continent_plates_margin(seed, n, &k_range(n), 1, m);
        let ok = s.chosen.is_some();
        tried.push((n, s));
        if ok {
            return (Some(n), tried);
        }
    }
    (None, tried)
}

/// The end land (64², h_src > 0 under the Airy sources' record) and its circumnavigation for one C1 initial state.
fn end_land(mut st: C1State) -> (f32, Circumnavigation) {
    let cfg = c1_run_config();
    let ss = SteinSteinParams::default();
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let src = run_c1_sources_mapped(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &ss, &[C1Term::Erosion], Mapping::Airy(Airy::declared()));
    let last = src.hist.normalized(None).pop().unwrap();
    let land: Vec<bool> = last.data.iter().map(|&v| v > 0.5).collect();
    (land.iter().filter(|&&b| b).count() as f32 / land.len() as f32, circumnavigation(&land, 64))
}

/// ADR Finding 170-A1/A2 -- **the selection with `plates` plates**: for every k in `ks`, every start plate and both families
/// (C1's BFS, the compact growth), C1 runs to its end (candidates with the same plate set share one run). A candidate
/// qualifies within [`SELECT_LAND`] and circumnavigable with r* ≥ `min_r`. The chosen one has the largest r*; ties go to the
/// start plate's rank in C1's RNG draw order, then BFS before compact, then the larger k. Deterministic from the seed.
pub fn select_continent_plates(seed: u64, plates: usize, ks: &[usize], min_r: usize) -> Selection {
    select_continent_plates_margin(seed, plates, ks, min_r, &Margin::off())
}

/// ADR Finding 174-B1 -- [`select_continent_plates`] from the margined initial states (the memo keyed by the plate set,
/// as before: the margin is a function of the plate set and the seed).
pub fn select_continent_plates_margin(seed: u64, plates: usize, ks: &[usize], min_r: usize, m: &Margin) -> Selection {
    let base = Phase2InitParams { num_plates: plates, ..Phase2InitParams::default() };
    let order = start_plate_order(&base.cluster, plates);
    let mut memo: std::collections::HashMap<Vec<bool>, (f32, Circumnavigation)> = std::collections::HashMap::new();
    let mut candidates = Vec::new();
    // the preference order: RNG rank, then BFS before compact, then the larger k
    for (rng_rank, &start) in order.iter().enumerate() {
        for compact in [false, true] {
            let mut ks_desc = ks.to_vec();
            ks_desc.sort_unstable_by(|a, b| b.cmp(a));
            for k in ks_desc {
                let p = C1Profile::Selected { plates, k, start, compact };
                let mut st = p.init_state(64, seed);
                // the plate set: which plates are continental
                let mut set = vec![false; plates];
                for (pid, ty) in st.plate_id.data().iter().zip(st.plate_type.data()) {
                    if *ty == crate::tectonics_v2::boundaries::plate_type::PlateType::Continental {
                        set[*pid as usize] = true;
                    }
                }
                let (land, circ) = *memo.entry(set).or_insert_with(|| {
                    margin::apply(&mut st, m, &p.init_params(), seed);
                    end_land(st)
                });
                candidates.push(Candidate { fraction: k as f64 / plates as f64, plates, k, compact, start, rng_rank, land, circ });
            }
        }
    }
    let mut chosen: Option<usize> = None;
    for (i, c) in candidates.iter().enumerate() {
        if c.qualifies_min(min_r) && chosen.is_none_or(|j| c.circ.sea_r > candidates[j].circ.sea_r) {
            chosen = Some(i);
        }
    }
    Selection { candidates, chosen, runs: memo.len() }
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
    /// F170-A5 (permanent): at 12 and 16 plates, the BFS family from C1's first RNG draw is C1's own init at that plate
    /// count bit for bit; the compact family grows a connected group of k plates.
    #[test]
    fn more_plates_keep_the_rng_init_and_the_compact_family_is_connected() {
        for plates in [12usize, 16] {
            let base = Phase2InitParams { num_plates: plates, ..Phase2InitParams::default() };
            let first = start_plate_order(&base.cluster, plates)[0];
            let k = plates / 2;
            let p = C1Profile::Selected { plates, k, start: first, compact: false };
            let params = p.init_params();
            let a = crate::tectonics_c1::init_r7::init_c1_state_phase_2_r7(64, 42, &params);
            let b = p.init_state(64, 42);
            assert_eq!(a.plate_type.data(), b.plate_type.data());
            assert_eq!(a.s.data().iter().map(|v| v.to_bits()).collect::<Vec<_>>(), b.s.data().iter().map(|v| v.to_bits()).collect::<Vec<_>>());
            let c = C1Profile::Selected { plates, k, start: first, compact: true }.init_state(64, 42);
            let mut set = vec![false; plates];
            for (pid, ty) in c.plate_id.data().iter().zip(c.plate_type.data()) {
                if *ty == PlateType::Continental {
                    set[*pid as usize] = true;
                }
            }
            assert_eq!(set.iter().filter(|&&s| s).count(), k, "{plates} plates");
            assert!(set[first as usize]);
        }
    }

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
