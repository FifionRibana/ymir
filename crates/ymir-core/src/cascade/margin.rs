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
//! Declared in `docs/reports/relief_method/f174_margin/f174_declared.md` (Parts B and C).

use super::amplify::noise01_t;
use crate::tectonics_c1::boundary_classification::classify_boundaries;
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
}

impl Margin {
    pub fn off() -> Self {
        Self::default()
    }
    /// F174-B's main choice: a rim of 2.5 cells and the declared ramp.
    pub fn declared() -> Self {
        Self { w_rim: 2.5, ramp: Some(Ramp::declared()), hetero: None }
    }
    pub fn is_off(&self) -> bool {
        self.w_rim <= 0.0 && self.ramp.is_none_or(|r| r.total_km() <= 0.0) && self.hetero.is_none()
    }
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
        let zero = Margin { w_rim: 0.0, ramp: Some(Ramp::declared().scaled(0.0)), hetero: None };
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
