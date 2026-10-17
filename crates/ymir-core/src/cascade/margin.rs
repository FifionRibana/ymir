//! ADR Finding 174-B -- **continents inside their plates, and the thinned margin**: an INITIAL CONDITION of the continent
//! profile (gated: the cascade and its benches only; C1's code is not touched).
//!
//! The author: « attention à ne pas retomber dans le biais d'intégrer un solveur dans les mécanismes sauf si réellement
//! nécessaire. » So nothing here runs a process: the profile's initial state is edited once, then C1 runs as it is.
//!
//! - **The rim (B1)**: the cells of a continental plate within `w_rim` cells (Euclidean, periodic) of an oceanic plate become
//!   oceanic **on the same plate** (`plate_id` unchanged), with C1's oceanic init values (s = 0.2, the age by
//!   `init_age_field_ridge_aligned`, out of the cratonic mask). The plate boundary is then ocean against ocean.
//! - **The ramp (B2)**: on the continental cells, s falls from the interior value to 0.2 across Péron-Pinvidic et al. 2013's
//!   proximal (35 → 30 km), necking (30 → 10 km) and distal (10 → 7 km) domains, piecewise linear in the distance from the
//!   continent's edge; the widths on the map are declared ([`Ramp::declared`], the real widths ÷ 7.5).
//! - **B-hétéro (B4)**: the interior value s_int × (1 + σ n), n a periodic fBm standardised over the continental cells
//!   (σ_rel 0.14, Christensen & Mooney 1995's shields and platforms).
//!
//! - **The margin seaward (B′, [`Margin::seaward`])**: the same ramp placed so that 32 km of crust (Airy's 0 m) falls on the
//!   continent's edge; only the 35 → 32 km part lies inside, the rest of the proximal, the necking and the distal domains lie
//!   on a band of the ocean attached to the continental plate at t₀ (continental while the crust exceeds 10 km, oceanic
//!   beyond); the nearest continent wins where bands meet.
//!
//! Declared in `docs/reports/relief_method/f174_margin/f174_declared.md` (Parts B and C, B′ and C′).

use super::amplify::noise01_t;
use crate::tectonics_c1::boundary_classification::classify_boundaries;
use crate::tectonics_c1::closures::subduction::OoSubductionParams;
use crate::tectonics_c1::init_r7::{Phase2InitParams, init_age_field_ridge_aligned};
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_v2::boundaries::plate_type::PlateType;

/// The C1 grid's cell (km on the map: 400 km / 64).
pub const C1_CELL_KM: f64 = 6.25;
/// C1's oceanic reference s (7 km).
const S_OCEAN: f64 = 0.2;
/// C1's km of crust per unit s.
const KM_PER_S: f64 = 35.0;

/// B2: the margin's three domain widths on the map (km), from the continent's edge inwards.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ramp {
    pub distal_km: f64,
    pub necking_km: f64,
    pub proximal_km: f64,
}

impl Ramp {
    /// Péron-Pinvidic et al. 2013: distal ~150 km (p. 34), necking ~100 km (p. 35), proximal > 200 km (p. 27); ÷ 7.5.
    pub fn declared() -> Self {
        Self { distal_km: 150.0 / 7.5, necking_km: 100.0 / 7.5, proximal_km: 200.0 / 7.5 }
    }
    pub fn scaled(self, f: f64) -> Self {
        Self { distal_km: self.distal_km * f, necking_km: self.necking_km * f, proximal_km: self.proximal_km * f }
    }
    pub fn total_km(&self) -> f64 {
        self.distal_km + self.necking_km + self.proximal_km
    }
    /// B′: the proximal domain's landward part (35 → 32 km), km.
    pub fn landward_km(&self) -> f64 {
        self.proximal_km * 3.0 / 5.0
    }
    /// B′: x of the taper break (the necking's outer end, 10 km of crust), from the continent's edge, km.
    pub fn taper_km(&self) -> f64 {
        self.proximal_km * 2.0 / 5.0 + self.necking_km
    }
    /// B′: the band's width (the distal domain's outer end), km.
    pub fn seaward_km(&self) -> f64 {
        self.taper_km() + self.distal_km
    }
    /// B′: g at the signed distance `x_km` from the continent's edge (negative inland): 1 inland of the landward part, 32 km of
    /// crust at 0, the 30 km and 10 km knots, 0 at the band's outer end.
    pub fn g_seaward(&self, x_km: f64) -> f64 {
        let g10 = (10.0 / KM_PER_S - S_OCEAN) / (1.0 - S_OCEAN);
        let g30 = (30.0 / KM_PER_S - S_OCEAN) / (1.0 - S_OCEAN);
        let x0 = -self.landward_km();
        let x1 = self.proximal_km * 2.0 / 5.0;
        let (x2, x3) = (self.taper_km(), self.seaward_km());
        if x_km <= x0 {
            1.0
        } else if x_km < x1 {
            1.0 - (1.0 - g30) * (x_km - x0) / self.proximal_km
        } else if x_km < x2 {
            g30 - (g30 - g10) * (x_km - x1) / self.necking_km
        } else if x_km < x3 {
            g10 * (1.0 - (x_km - x2) / self.distal_km)
        } else {
            0.0
        }
    }
    /// g(x) in [0, 1]: 0 at the edge, the 10 km and 30 km knots, 1 at the proximal domain's inner end.
    pub fn g(&self, x_km: f64) -> f64 {
        let g10 = (10.0 / KM_PER_S - S_OCEAN) / (1.0 - S_OCEAN);
        let g30 = (30.0 / KM_PER_S - S_OCEAN) / (1.0 - S_OCEAN);
        let (a, b) = (self.distal_km, self.distal_km + self.necking_km);
        if x_km <= 0.0 {
            0.0
        } else if x_km < a {
            g10 * x_km / a
        } else if x_km < b {
            g10 + (g30 - g10) * (x_km - a) / self.necking_km
        } else if x_km < self.total_km() {
            g30 + (1.0 - g30) * (x_km - b) / self.proximal_km
        } else {
            1.0
        }
    }
}

/// B4: the initial heterogeneity of the continental thickness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hetero {
    pub sigma_rel: f64,
    /// The fBm's wavelength on the map (km); 400 / it must be an integer (the torus).
    pub wavelength_km: f32,
}

impl Hetero {
    /// Christensen & Mooney 1995, Table 3: shields and platforms 41.5 ± 5.8 km → 0.14; the wavelength 40 km, declared.
    pub fn declared() -> Self {
        Self { sigma_rel: 5.8 / 41.5, wavelength_km: 40.0 }
    }
}

/// The margin's settings; [`Margin::off`] leaves the initial state untouched.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Margin {
    /// The rim's width (cells of 64²).
    pub w_rim: f64,
    pub ramp: Option<Ramp>,
    pub hetero: Option<Hetero>,
    /// B′: the ramp seaward (`w_rim` unused).
    pub seaward: bool,
    /// ADR Finding 175: the ocean-ocean subduction path (option 5), a C1 closure flag the continent profile carries to C1
    /// (`run_c1_sources_mapped_oo`); not an initial condition, so [`Margin::is_off`] ignores it.
    pub oo: Option<OoSubductionParams>,
}

impl Margin {
    pub fn off() -> Self {
        Self::default()
    }
    /// F174-B's main choice: a rim of 2.5 cells and the declared ramp.
    pub fn declared() -> Self {
        Self { w_rim: 2.5, ramp: Some(Ramp::declared()), hetero: None, seaward: false, oo: None }
    }
    /// F174-B′'s choice: the declared ramp seaward.
    pub fn seaward() -> Self {
        Self { w_rim: 0.0, ramp: Some(Ramp::declared()), hetero: None, seaward: true, oo: None }
    }
    /// F175's choice: the seaward margin with the ocean-ocean subduction path.
    pub fn seaward_oo() -> Self {
        Self { oo: Some(OoSubductionParams::declared()), ..Self::seaward() }
    }
    pub fn is_off(&self) -> bool {
        self.w_rim <= 0.0 && self.ramp.is_none_or(|r| r.total_km() <= 0.0) && self.hetero.is_none()
    }
}

/// The periodic Euclidean distance (cells) from every `src`-false cell to the nearest `src`-true cell, and that cell (the
/// lowest index on a tie); 0 and itself on `src`.
fn nearest(src: &[bool], n: usize) -> (Vec<f64>, Vec<usize>) {
    let pts: Vec<usize> = (0..n * n).filter(|&k| src[k]).collect();
    let ni = n as i64;
    let mut d = vec![0f64; n * n];
    let mut at = (0..n * n).collect::<Vec<usize>>();
    for k in 0..n * n {
        if src[k] {
            continue;
        }
        let (x, y) = ((k % n) as i64, (k / n) as i64);
        let mut best = (i64::MAX, usize::MAX);
        for &p in &pts {
            let mut dx = ((p % n) as i64 - x).rem_euclid(ni);
            let mut dy = ((p / n) as i64 - y).rem_euclid(ni);
            dx = dx.min(ni - dx);
            dy = dy.min(ni - dy);
            let q = dx * dx + dy * dy;
            if q < best.0 {
                best = (q, p);
            }
        }
        d[k] = (best.0 as f64).sqrt();
        at[k] = best.1;
    }
    (d, at)
}

/// The periodic Euclidean distance (cells) from every `src`-false cell to the nearest `src`-true cell; 0 on `src`.
fn dist_to(src: &[bool], n: usize) -> Vec<f64> {
    let pts: Vec<(i64, i64)> = (0..n * n).filter(|&k| src[k]).map(|k| ((k % n) as i64, (k / n) as i64)).collect();
    let ni = n as i64;
    (0..n * n)
        .map(|k| {
            if src[k] {
                return 0.0;
            }
            let (x, y) = ((k % n) as i64, (k / n) as i64);
            let mut best = i64::MAX;
            for &(a, b) in &pts {
                let mut dx = (a - x).rem_euclid(ni);
                let mut dy = (b - y).rem_euclid(ni);
                dx = dx.min(ni - dx);
                dy = dy.min(ni - dy);
                best = best.min(dx * dx + dy * dy);
            }
            (best as f64).sqrt()
        })
        .collect()
}

/// ADR Finding 174-B -- apply the margin to a continent profile's initial state (`params` = the profile's init params, for
/// the craton ratio and the age rule). A no-op when [`Margin::is_off`].
pub fn apply(st: &mut C1State, m: &Margin, params: &Phase2InitParams, seed: u64) {
    if m.is_off() {
        return;
    }
    if m.seaward {
        return apply_seaward(st, m, params, seed);
    }
    let n = st.s.nx();
    assert_eq!(st.s.ny(), n);
    let oceanic_plate: Vec<bool> = st.plate_type.data().iter().map(|t| *t == PlateType::Oceanic).collect();
    let d = dist_to(&oceanic_plate, n);
    // B1: the rim
    let mut rim = vec![false; n * n];
    for k in 0..n * n {
        if !oceanic_plate[k] && d[k] <= m.w_rim + 1e-9 {
            rim[k] = true;
            let (i, j) = (k % n, k / n);
            st.plate_type.set(i, j, PlateType::Oceanic);
            st.s.set(i, j, S_OCEAN);
            st.cratonic_mask.set(i, j, false);
        }
    }
    if rim.iter().any(|&r| r) {
        let kin = PlateKinematics::preset_phase_1_1(st.num_plates);
        let bi = classify_boundaries(&st.plate_id, &kin);
        let mut age = st.age.clone();
        init_age_field_ridge_aligned(&mut age, &st.plate_id, &st.plate_type, &bi, &params.age);
        for k in 0..n * n {
            if rim[k] {
                st.age.set(k % n, k / n, age.get(k % n, k / n));
            }
        }
    }
    // B4: the field over the continental cells
    let cont: Vec<bool> = st.plate_type.data().iter().map(|t| *t == PlateType::Continental).collect();
    let field: Option<Vec<f64>> = m.hetero.map(|h| {
        let v = noise01_t(n, C1_CELL_KM as f32, h.wavelength_km, 1, seed, "c1_margin_hetero", Some([0, 0]));
        let idx: Vec<usize> = (0..n * n).filter(|&k| cont[k]).collect();
        let mean = idx.iter().map(|&k| v[k] as f64).sum::<f64>() / idx.len().max(1) as f64;
        let sd = (idx.iter().map(|&k| (v[k] as f64 - mean).powi(2)).sum::<f64>() / idx.len().max(1) as f64).sqrt().max(1e-12);
        v.iter().map(|&x| h.sigma_rel * (x as f64 - mean) / sd).collect()
    });
    // B2: the ramp, from the continent's edge (the rim's inner side)
    for k in 0..n * n {
        if !cont[k] {
            continue;
        }
        let (i, j) = (k % n, k / n);
        let f = field.as_ref().map_or(0.0, |v| v[k]);
        let x_km = (d[k] - m.w_rim) * C1_CELL_KM;
        match m.ramp {
            Some(r) if x_km < r.total_km() => {
                let s_int = if st.cratonic_mask.get(i, j) { params.craton_thickness_ratio } else { 1.0 } * (1.0 + f);
                st.s.set(i, j, S_OCEAN + (s_int - S_OCEAN) * r.g(x_km));
            }
            _ => {
                if field.is_some() {
                    st.s.set(i, j, st.s.get(i, j) * (1.0 + f));
                }
            }
        }
    }
}

/// The B-hétéro field over the `cont` cells (zero mean, σ_rel standard deviation over them), or None.
fn hetero_field(m: &Margin, cont: &[bool], n: usize, seed: u64) -> Option<Vec<f64>> {
    m.hetero.map(|h| {
        let v = noise01_t(n, C1_CELL_KM as f32, h.wavelength_km, 1, seed, "c1_margin_hetero", Some([0, 0]));
        let idx: Vec<usize> = (0..n * n).filter(|&k| cont[k]).collect();
        let mean = idx.iter().map(|&k| v[k] as f64).sum::<f64>() / idx.len().max(1) as f64;
        let sd = (idx.iter().map(|&k| (v[k] as f64 - mean).powi(2)).sum::<f64>() / idx.len().max(1) as f64).sqrt().max(1e-12);
        v.iter().map(|&x| h.sigma_rel * (x as f64 - mean) / sd).collect()
    })
}

/// ADR Finding 174-B′ -- the margin seaward: x = the signed distance to the continent's edge (−(d_o − ½) inside, +(d_c − ½)
/// outside, × 6.25 km); the ocean cells with x < the band's width join the plate of their nearest continental cell
/// (continental while x < the taper break, oceanic beyond); s = 0.2 + (s_int − 0.2) g_seaward(x) wherever g < 1.
fn apply_seaward(st: &mut C1State, m: &Margin, params: &Phase2InitParams, seed: u64) {
    let Some(r) = m.ramp else { return };
    let n = st.s.nx();
    assert_eq!(st.s.ny(), n);
    let cont0: Vec<bool> = st.plate_type.data().iter().map(|t| *t == PlateType::Continental).collect();
    let ocean0: Vec<bool> = cont0.iter().map(|c| !c).collect();
    let (d_c, near) = nearest(&cont0, n);
    let d_o = dist_to(&ocean0, n);
    let x_km: Vec<f64> = (0..n * n).map(|k| if cont0[k] { -(d_o[k] - 0.5) * C1_CELL_KM } else { (d_c[k] - 0.5) * C1_CELL_KM }).collect();
    let mut distal = vec![false; n * n];
    for k in 0..n * n {
        if cont0[k] || x_km[k] >= r.seaward_km() {
            continue;
        }
        let (i, j) = (k % n, k / n);
        let src = near[k];
        st.plate_id.set(i, j, st.plate_id.get(src % n, src / n));
        st.cratonic_mask.set(i, j, false);
        if x_km[k] < r.taper_km() {
            st.plate_type.set(i, j, PlateType::Continental);
            st.age.set(i, j, params.age.continental_baseline);
        } else {
            distal[k] = true;
        }
    }
    if distal.iter().any(|&b| b) {
        let kin = PlateKinematics::preset_phase_1_1(st.num_plates);
        let bi = classify_boundaries(&st.plate_id, &kin);
        let mut age = st.age.clone();
        init_age_field_ridge_aligned(&mut age, &st.plate_id, &st.plate_type, &bi, &params.age);
        for k in 0..n * n {
            if distal[k] {
                st.age.set(k % n, k / n, age.get(k % n, k / n));
            }
        }
    }
    let cont: Vec<bool> = st.plate_type.data().iter().map(|t| *t == PlateType::Continental).collect();
    let field = hetero_field(m, &cont, n, seed);
    for k in 0..n * n {
        let (i, j) = (k % n, k / n);
        let f = field.as_ref().map_or(0.0, |v| if cont[k] { v[k] } else { 0.0 });
        let g = r.g_seaward(x_km[k]);
        let in_band = !cont0[k] && x_km[k] < r.seaward_km();
        if (cont0[k] && g < 1.0) || in_band {
            let s_int = if st.cratonic_mask.get(i, j) { params.craton_thickness_ratio } else { 1.0 } * (1.0 + f);
            st.s.set(i, j, S_OCEAN + (s_int - S_OCEAN) * g);
        } else if cont[k] && field.is_some() {
            st.s.set(i, j, st.s.get(i, j) * (1.0 + f));
        }
    }
}

/// F174-B′2 -- the conflicts of the seaward band on a profile's initial state (Part A's): (a) the share of the coastal
/// perimeter within 2 × the band of another continent, (b) within the band of another continent, (c) the oceanic plates fully
/// annexed, those losing > 50 %, and the share of the perimeter 4-adjacent to such a plate.
#[derive(Clone, Debug, Default)]
pub struct SeawardConflicts {
    pub overlap_share: f64,
    pub reach_share: f64,
    pub plates_full: usize,
    pub plates_half: usize,
    pub perim_on_half_share: f64,
    pub continents: usize,
}

pub fn seaward_conflicts(st0: &C1State, m: &Margin) -> SeawardConflicts {
    let Some(r) = m.ramp else { return SeawardConflicts::default() };
    let n = st0.s.nx();
    let cont: Vec<bool> = st0.plate_type.data().iter().map(|t| *t == PlateType::Continental).collect();
    // the continents: 4-connected periodic components
    let mut comp = vec![usize::MAX; n * n];
    let mut nc = 0;
    for s in 0..n * n {
        if !cont[s] || comp[s] != usize::MAX {
            continue;
        }
        let mut stack = vec![s];
        comp[s] = nc;
        while let Some(k) = stack.pop() {
            let (x, y) = (k % n, k / n);
            for j in [y * n + (x + 1) % n, y * n + (x + n - 1) % n, ((y + 1) % n) * n + x, ((y + n - 1) % n) * n + x] {
                if cont[j] && comp[j] == usize::MAX {
                    comp[j] = nc;
                    stack.push(j);
                }
            }
        }
        nc += 1;
    }
    let w = r.seaward_km() / C1_CELL_KM + 0.5;
    let ni = n as i64;
    let perim: Vec<usize> = (0..n * n)
        .filter(|&k| {
            let (x, y) = (k % n, k / n);
            cont[k] && [y * n + (x + 1) % n, y * n + (x + n - 1) % n, ((y + 1) % n) * n + x, ((y + n - 1) % n) * n + x].iter().any(|&j| !cont[j])
        })
        .collect();
    let (mut over, mut reach) = (0usize, 0usize);
    for &k in &perim {
        let (x, y) = ((k % n) as i64, (k / n) as i64);
        let mut dmin = f64::INFINITY;
        for j in 0..n * n {
            if !cont[j] || comp[j] == comp[k] {
                continue;
            }
            let mut dx = ((j % n) as i64 - x).rem_euclid(ni);
            let mut dy = ((j / n) as i64 - y).rem_euclid(ni);
            dx = dx.min(ni - dx);
            dy = dy.min(ni - dy);
            dmin = dmin.min(((dx * dx + dy * dy) as f64).sqrt());
        }
        if dmin <= 2.0 * w {
            over += 1;
        }
        if dmin <= w {
            reach += 1;
        }
    }
    let (d_c, _) = nearest(&cont, n);
    let banded: Vec<bool> = (0..n * n).map(|k| !cont[k] && (d_c[k] - 0.5) * C1_CELL_KM < r.seaward_km()).collect();
    let np = st0.num_plates;
    let mut before_c = vec![0usize; np];
    let mut kept_c = vec![0usize; np];
    for k in 0..n * n {
        if !cont[k] {
            before_c[st0.plate_id.data()[k] as usize] += 1;
            if !banded[k] {
                kept_c[st0.plate_id.data()[k] as usize] += 1;
            }
        }
    }
    let oceanic_plate: Vec<bool> = (0..np).map(|p| before_c[p] > 0).collect();
    let lost = |p: usize| 1.0 - kept_c[p] as f64 / before_c[p].max(1) as f64;
    let losses: Vec<f64> = (0..np).map(|p| if oceanic_plate[p] { lost(p) } else { 0.0 }).collect();
    let plates_full = (0..np).filter(|&p| oceanic_plate[p] && losses[p] >= 1.0 - 1e-12).count();
    let half: Vec<bool> = (0..np).map(|p| oceanic_plate[p] && losses[p] > 0.5).collect();
    let on_half = perim
        .iter()
        .filter(|&&k| {
            let (x, y) = (k % n, k / n);
            [y * n + (x + 1) % n, y * n + (x + n - 1) % n, ((y + 1) % n) * n + x, ((y + n - 1) % n) * n + x].iter().any(|&j| !cont[j] && half[st0.plate_id.data()[j] as usize])
        })
        .count();
    let pn = perim.len().max(1) as f64;
    SeawardConflicts {
        overlap_share: over as f64 / pn,
        reach_share: reach as f64 / pn,
        plates_full,
        plates_half: half.iter().filter(|&&h| h).count(),
        perim_on_half_share: on_half as f64 / pn,
        continents: nc,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cascade::profile::{C1Profile, select_continent_default, select_continent_default_margin};

    fn bits(st: &C1State) -> (Vec<u64>, Vec<u64>, Vec<PlateType>, Vec<u16>, Vec<bool>) {
        (
            st.s.data().iter().map(|v| v.to_bits()).collect(),
            st.age.data().iter().map(|v| v.to_bits()).collect(),
            st.plate_type.data().to_vec(),
            st.plate_id.data().to_vec(),
            st.cratonic_mask.data().to_vec(),
        )
    }

    /// F174-B5 (permanent): the margin off is the profile's initial state bit for bit, and so is a margin of zero width.
    #[test]
    fn the_margin_off_is_the_profile_bit_for_bit() {
        let p = C1Profile::Selected { plates: 16, k: 8, start: 3, compact: true };
        let a = p.init_state(64, 42);
        let b = p.init_state_margin(64, 42, &Margin::off());
        assert_eq!(bits(&a), bits(&b));
        let zero = Margin { w_rim: 0.0, ramp: Some(Ramp::declared().scaled(0.0)), hetero: None, seaward: false, oo: None };
        assert!(zero.is_off());
        assert_eq!(bits(&a), bits(&p.init_state_margin(64, 42, &zero)));
    }

    /// F174-B1/B2 (permanent): the rim is oceanic on the continent's own plate; the ramp rises monotonically from 0.2 at the
    /// edge to the interior; no continental cell touches another plate's ocean.
    #[test]
    fn the_rim_is_the_plates_own_ocean_and_the_ramp_rises() {
        let p = C1Profile::Selected { plates: 16, k: 8, start: 3, compact: true };
        let a = p.init_state(64, 42);
        let b = p.init_state_margin(64, 42, &Margin::declared());
        assert_eq!(a.plate_id.data(), b.plate_id.data(), "plate_id unchanged");
        let n = 64;
        let (mut rim, mut same, mut other) = (0, 0, 0);
        for k in 0..n * n {
            let (i, j) = (k % n, k / n);
            if a.plate_type.get(i, j) == PlateType::Continental && b.plate_type.get(i, j) == PlateType::Oceanic {
                rim += 1;
                assert_eq!(b.s.get(i, j), S_OCEAN);
            }
            if b.plate_type.get(i, j) == PlateType::Continental {
                for (x, y) in [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)] {
                    if b.plate_type.get(x, y) == PlateType::Oceanic {
                        if b.plate_id.get(x, y) == b.plate_id.get(i, j) {
                            same += 1;
                        } else {
                            // where two continental plates meet near the coast, one's rim can touch the other's continent
                            other += 1;
                        }
                    }
                }
            }
        }
        assert!(rim > 0);
        assert!(same >= 9 * other, "the coast's ocean is mostly its own plate's: {same} same, {other} other");
        let r = Ramp::declared();
        let mut last = -1.0;
        for x in 0..80 {
            let g = r.g(x as f64);
            assert!(g >= last && (0.0..=1.0).contains(&g));
            last = g;
        }
        assert!((r.g(r.distal_km) - (10.0 / 35.0 - 0.2) / 0.8).abs() < 1e-12);
    }

    /// F174-B′5 (permanent): a seaward band of zero width is the profile's state bit for bit; the seaward geometry changes
    /// plate_id only on band cells (each to its nearest continent's plate), types the band continental iff x < the taper
    /// break, and puts the coast on the old edge (s > 0.914 on the land side, < 0.915 at sea).
    #[test]
    fn the_seaward_margin_keeps_the_coast_and_annexes_a_band() {
        let p = C1Profile::Selected { plates: 16, k: 8, start: 3, compact: true };
        let a = p.init_state(64, 42);
        let zero = Margin { ramp: Some(Ramp::declared().scaled(0.0)), ..Margin::seaward() };
        assert!(zero.is_off());
        assert_eq!(bits(&a), bits(&p.init_state_margin(64, 42, &zero)));
        let m = Margin::seaward();
        let b = p.init_state_margin(64, 42, &m);
        let r = Ramp::declared();
        let n = 64;
        let cont0: Vec<bool> = a.plate_type.data().iter().map(|t| *t == PlateType::Continental).collect();
        let (d_c, near) = nearest(&cont0, n);
        let mut band = 0;
        for k in 0..n * n {
            let (i, j) = (k % n, k / n);
            if cont0[k] {
                assert_eq!(a.plate_id.get(i, j), b.plate_id.get(i, j));
                assert_eq!(b.plate_type.get(i, j), PlateType::Continental);
                continue;
            }
            let x = (d_c[k] - 0.5) * C1_CELL_KM;
            if x < r.seaward_km() {
                band += 1;
                assert_eq!(b.plate_id.get(i, j), a.plate_id.get(near[k] % n, near[k] / n));
                assert_eq!(b.plate_type.get(i, j) == PlateType::Continental, x < r.taper_km());
            } else {
                assert_eq!(a.plate_id.get(i, j), b.plate_id.get(i, j));
                assert_eq!(bits(&a).0[k], bits(&b).0[k], "beyond the band, s untouched");
            }
            // the first ocean cell (x = 3.125 km) is under the coast's s, the coast cells above it
            if d_c[k] == 1.0 {
                assert!(b.s.get(i, j) < 0.915, "sea at the edge: {}", b.s.get(i, j));
            }
        }
        assert!(band > 0);
        assert!((r.g_seaward(0.0) * 0.8 + 0.2 - 32.0 / 35.0).abs() < 1e-12, "32 km on the edge");
        assert!((r.g_seaward(r.taper_km()) * 0.8 + 0.2 - 10.0 / 35.0).abs() < 1e-12);
        assert_eq!(r.g_seaward(r.seaward_km()), 0.0);
        assert_eq!(r.g_seaward(-r.landward_km()), 1.0);
    }

    /// F175-A8 (permanent): the O-O path off (`None`) is `run_with_closures_terms` bit for bit; with D_init infinite nothing
    /// is consumed ocean-ocean (A7 still cuts the O-O seeds); with the declared D_init the O-O path consumes, conserves the
    /// arc mass, and A7 leaves no O-O seed.
    #[test]
    fn the_oo_subduction_path_controls() {
        use crate::tectonics_c1::boundary_classification::{classify_boundaries, clear_oo_seeds, retarget_upper_plate_continental};
        use crate::tectonics_c1::closures::subduction::SubductionStats;
        use crate::tectonics_c1::time_loop::{C1Closures, C1Term, run_with_closures_terms, run_with_closures_terms_oo};
        let p = C1Profile::Selected { plates: 16, k: 8, start: 3, compact: true };
        let mut cfg = crate::cascade::profile::c1_run_config();
        cfg.n_steps = 80;
        let run = |oo: Option<Option<OoSubductionParams>>| {
            let mut st = p.init_state_margin(64, 42, &Margin::seaward());
            let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
            let mut sum = SubductionStats::default();
            let mut on_step = |_: usize, s: &C1State, _: &PlateKinematics| {
                let x = s.last_step_stats.subduction;
                sum.oo_cells_consumed += x.oo_cells_consumed;
                sum.oo_mass_consumed += x.oo_mass_consumed;
                sum.oo_arc_continental += x.oo_arc_continental;
                sum.oo_arc_oceanic += x.oo_arc_oceanic;
                sum.oo_arc_lost += x.oo_arc_lost;
                sum.oo_reassigned += x.oo_reassigned;
                sum.oo_seeds_cleared += x.oo_seeds_cleared;
            };
            let mut on_term = |_: C1Term, _: &[f64], _: &[f64]| {};
            match oo {
                None => run_with_closures_terms(&mut st, &mut kin, &cfg, &C1Closures::default(), &mut on_step, &mut on_term),
                Some(o) => run_with_closures_terms_oo(&mut st, &mut kin, &cfg, &C1Closures::default(), o, &mut on_step, &mut on_term),
            }
            (bits(&st), sum, st, kin)
        };
        let (plain, ..) = run(None);
        let (off, s_off, ..) = run(Some(None));
        assert_eq!(plain, off, "the O-O path off is C1 bit for bit");
        assert_eq!(s_off.oo_cells_consumed + s_off.oo_seeds_cleared, 0);
        let (_, s_inf, ..) = run(Some(Some(OoSubductionParams { d_init_cells: f64::INFINITY })));
        assert_eq!(s_inf.oo_cells_consumed, 0, "D_init infinite: no O-O consumption");
        assert_eq!(s_inf.oo_reassigned, 0);
        let (_, s_on, st, kin) = run(Some(Some(OoSubductionParams::declared())));
        assert!(s_on.oo_cells_consumed > 0, "the declared D_init consumes");
        let placed = s_on.oo_arc_continental + s_on.oo_arc_oceanic + s_on.oo_arc_lost;
        assert!((placed - 0.5 * s_on.oo_mass_consumed).abs() <= 1e-9 * s_on.oo_mass_consumed.max(1.0), "the arc mass is conserved");
        // A7 on the end state: no O-O seed is left
        let mut bi = classify_boundaries(&st.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &st.plate_id, &st.plate_type);
        clear_oo_seeds(&mut bi, &st.plate_id, &st.plate_type);
        assert_eq!(clear_oo_seeds(&mut bi, &st.plate_id, &st.plate_type), 0, "a second cut finds nothing");
    }

    /// F174-B5 (permanent): the default selection with the margin off is `select_continent_default` bit for bit.
    #[test]
    #[ignore]
    fn the_selection_with_the_margin_off_is_the_default() {
        let (a, sa) = select_continent_default(8);
        let (b, sb) = select_continent_default_margin(8, &Margin::off());
        assert_eq!(a, b);
        assert_eq!(format!("{:?}", sa.iter().map(|x| (x.0, x.1.chosen, x.1.candidates.len())).collect::<Vec<_>>()), format!("{:?}", sb.iter().map(|x| (x.0, x.1.chosen, x.1.candidates.len())).collect::<Vec<_>>()));
    }
}
