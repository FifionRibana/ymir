//! ADR Finding 168 -- **the cascade's C1 parameter profile** (gated: the cascade and its benches only; C1's
//! production default is unchanged).
//!
//! « Production » is `Phase2InitParams::default()`. A « continent » profile overrides ONE parameter, the initial
//! continental fraction (`ContinentalClusterParams::continental_fraction`, a fraction of the plates: round(8 f)
//! continental plates out of 8). It is a design target, not a physical calibration: L1 and L2 were set once on the
//! témoin for ~40 % and 55–60 % land at the end of the history-driven physics level, then frozen for every seed.
//! Declared in `docs/reports/relief_method/f168_continent_size/f168_declared.md`.

use crate::tectonics_c1::init_r7::Phase2InitParams;

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
}

impl C1Profile {
    pub const L1: Self = Self::Continent(CONTINENT_L1);
    pub const L2: Self = Self::Continent(CONTINENT_L2);

    /// The C1 init parameters of this profile.
    pub fn init_params(self) -> Phase2InitParams {
        let mut p = Phase2InitParams::default();
        if let Self::Continent(f) = self {
            p.cluster.continental_fraction = f;
        }
        p
    }
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
}
