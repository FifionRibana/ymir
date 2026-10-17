//! F175 diagnostic: a hash of C1's production run (64², 300 steps, the default closures), to compare two code states.
use std::hash::{Hash, Hasher};
use ymir_core::cascade::profile::c1_run_config;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{C1Closures, run_with_closures};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

#[test]
#[ignore]
fn f175_c1_hash() {
    for seed in [42u64, 8, 10481999410520546993] {
        let cfg = c1_run_config();
        let mut st = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
        run_with_closures(&mut st, &mut kin, &cfg, &C1Closures::default(), |_, _| {});
        let mut h = std::collections::hash_map::DefaultHasher::new();
        for v in st.s.data() {
            v.to_bits().hash(&mut h);
        }
        for v in st.age.data() {
            v.to_bits().hash(&mut h);
        }
        st.plate_id.data().hash(&mut h);
        for t in st.plate_type.data() {
            (*t == PlateType::Continental).hash(&mut h);
        }
        eprintln!("seed {seed}: C1 hash {:016x}", h.finish());
    }
}
