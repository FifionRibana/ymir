//! ADR Finding 151-F / 152-B1 — the tectonic HISTORY of the coarse grid, recorded during the C1 loop through the
//! read-only observer (`run_with_closures_observed`): how long each cell was a convergent / collision / active-margin
//! boundary, and the SUTURES of the plate merges. It feeds the fossil belts and the structural density of the zoning.
//!
//! The geology runs its own coarse pass (≈ 0.2 s at 64²): the terrain's cached tectonic run is never touched, so the
//! history cannot change the field.

use std::collections::{HashMap, HashSet};

use crate::grid::GridF32;
use crate::tectonics_c1::closures::fracture::{FractureConfig, derive_coarse_density};
use crate::tectonics_c1::debug_labels::{CoarseTectonicLabels, derive_tectonic_labels};
use crate::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures_observed};
use crate::tectonics_v2::boundaries::plate_type::PlateType;

/// A collision lasting at least this many steps leaves a fossil belt (PROXY, F151: ≈ 7 Ma at ≈ 0.67 Ma per step).
pub const FOSSIL_COLLISION_STEPS: u32 = 10;
/// The active belt: a cell within this many coarse cells of a convergent boundary at the end (PROXY, F151).
pub const ACTIVE_BELT_DILATION: i64 = 2;
/// The structural density's decay from the sutures and the past collisions (km): C-3b's own decay (25 km).
pub const STRUCTURAL_DECAY_KM: f32 = 25.0;

/// A plate merge (an accretion): the plate that vanished, the one it went into, and its suture.
#[derive(Clone, Debug)]
pub struct Merge {
    pub step: usize,
    pub vanished: u16,
    pub into: u16,
    pub cells: usize,
    pub continental_share: f32,
}

/// The coarse grid's tectonic history, all row-major `nx·ny`.
#[derive(Clone, Debug)]
pub struct TectonicHistory {
    pub nx: usize,
    pub ny: usize,
    pub convergent_steps: Vec<u32>,
    pub collision_steps: Vec<u32>,
    pub upper_plate_steps: Vec<u32>,
    pub last_collision_step: Vec<u32>,
    /// A cell of the two plates touching each other on the step before they merged.
    pub suture: Vec<bool>,
    pub merges: Vec<Merge>,
}

/// Run the coarse tectonics with the history recorded. Returns the SAME settled state and kinematics as
/// `run_coarse_tectonics` (the observer is read-only; F151 measured the state bit-identical).
#[must_use]
pub fn record_history(
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
) -> (C1State, PlateKinematics, TectonicHistory) {
    let mut state = init_c1_state_phase_2_r7(grid, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    let (nx, ny) = (state.nx(), state.ny());
    let nc = nx * ny;
    let mut h = TectonicHistory {
        nx,
        ny,
        convergent_steps: vec![0; nc],
        collision_steps: vec![0; nc],
        upper_plate_steps: vec![0; nc],
        last_collision_step: vec![0; nc],
        suture: vec![false; nc],
        merges: Vec::new(),
    };
    let mut prev_id: Vec<u16> = (0..nc).map(|c| state.plate_id.get(c % nx, c / nx)).collect();
    let mut prev_cont: Vec<bool> = (0..nc).map(|c| state.plate_type.get(c % nx, c / nx) == PlateType::Continental).collect();
    run_with_closures_observed(&mut state, &mut kin, run, closures, |step, s, k| {
        let lab = derive_tectonic_labels(s, k);
        for c in 0..nc {
            if lab.subduction_upper[c] || lab.subduction_slab[c] || lab.collision[c] {
                h.convergent_steps[c] += 1;
            }
            if lab.collision[c] {
                h.collision_steps[c] += 1;
                h.last_collision_step[c] = step as u32;
            }
            if lab.subduction_upper[c] {
                h.upper_plate_steps[c] += 1;
            }
        }
        let now: Vec<u16> = (0..nc).map(|c| s.plate_id.get(c % nx, c / nx)).collect();
        let ids_now: HashSet<u16> = now.iter().copied().collect();
        let mut vanished: Vec<u16> = prev_id.iter().copied().collect::<HashSet<u16>>().difference(&ids_now).copied().collect();
        vanished.sort_unstable();
        for a in vanished {
            let mut votes: HashMap<u16, usize> = HashMap::new();
            for c in 0..nc {
                if prev_id[c] == a {
                    *votes.entry(now[c]).or_insert(0) += 1;
                }
            }
            let Some((&b, _)) = votes.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))) else { continue };
            let mut cells = 0usize;
            let mut cont = 0usize;
            for c in 0..nc {
                if prev_id[c] != a && prev_id[c] != b {
                    continue;
                }
                let other = if prev_id[c] == a { b } else { a };
                let (i, j) = ((c % nx) as i64, (c / nx) as i64);
                let touches = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)]
                    .iter()
                    .any(|(di, dj)| prev_id[((j + dj).rem_euclid(ny as i64) as usize) * nx + (i + di).rem_euclid(nx as i64) as usize] == other);
                if touches {
                    h.suture[c] = true;
                    cells += 1;
                    cont += prev_cont[c] as usize;
                }
            }
            h.merges.push(Merge { step, vanished: a, into: b, cells, continental_share: cont as f32 / cells.max(1) as f32 });
        }
        prev_id = now;
        prev_cont = (0..nc).map(|c| s.plate_type.get(c % nx, c / nx) == PlateType::Continental).collect();
    });
    (state, kin, h)
}

impl TectonicHistory {
    /// The FOSSIL belt (F151's PROXY rule): a collision of ≥ [`FOSSIL_COLLISION_STEPS`], or a suture cell on
    /// continental crust (at the end), and not convergent at the end.
    #[must_use]
    pub fn fossil_belt(&self, state: &C1State, end: &CoarseTectonicLabels) -> Vec<bool> {
        let nx = self.nx;
        (0..self.nx * self.ny)
            .map(|c| {
                let cont = state.plate_type.get(c % nx, c / nx) == PlateType::Continental;
                let conv_now = end.subduction_upper[c] || end.subduction_slab[c] || end.collision[c];
                !conv_now && (self.collision_steps[c] >= FOSSIL_COLLISION_STEPS || (self.suture[c] && cont))
            })
            .collect()
    }

    /// The ACTIVE belt: within [`ACTIVE_BELT_DILATION`] coarse cells of a convergent boundary at the end.
    #[must_use]
    pub fn active_belt(&self, end: &CoarseTectonicLabels) -> Vec<bool> {
        let (nx, ny) = (self.nx, self.ny);
        let mut m = vec![false; nx * ny];
        for c in 0..nx * ny {
            if end.subduction_upper[c] || end.subduction_slab[c] || end.collision[c] {
                let (i, j) = ((c % nx) as i64, (c / nx) as i64);
                for dj in -ACTIVE_BELT_DILATION..=ACTIVE_BELT_DILATION {
                    for di in -ACTIVE_BELT_DILATION..=ACTIVE_BELT_DILATION {
                        m[((j + dj).rem_euclid(ny as i64) as usize) * nx + (i + di).rem_euclid(nx as i64) as usize] = true;
                    }
                }
            }
        }
        m
    }

    /// The structural density of the ZONING (never the incision): `max(C-3b density, exp(−d / 25 km))`, d the
    /// distance to the sutures and the past collision cells. With `c3b_only`, C-3b's density alone (G-veins' control).
    #[must_use]
    pub fn structural_density(&self, state: &C1State, kin: &PlateKinematics, fcfg: &FractureConfig, c3b_only: bool) -> GridF32 {
        let (nx, ny) = (self.nx, self.ny);
        let mut d = derive_coarse_density(state, kin, fcfg, None);
        if c3b_only {
            return d;
        }
        let seeds: Vec<usize> = (0..nx * ny).filter(|&c| self.suture[c] || self.collision_steps[c] > 0).collect();
        if seeds.is_empty() {
            return d;
        }
        let km_per_cell = fcfg.domain_km / nx as f32;
        for c in 0..nx * ny {
            let (i, j) = ((c % nx) as i64, (c / nx) as i64);
            let mut best = f32::INFINITY;
            for &s in &seeds {
                let (si, sj) = ((s % nx) as i64, (s / nx) as i64);
                let dx = (i - si).rem_euclid(nx as i64).min((si - i).rem_euclid(nx as i64)) as f32;
                let dy = (j - sj).rem_euclid(ny as i64).min((sj - j).rem_euclid(ny as i64)) as f32;
                best = best.min((dx * dx + dy * dy).sqrt());
            }
            let s = (-(best * km_per_cell) / STRUCTURAL_DECAY_KM).exp();
            let (ii, jj) = (c % nx, c / nx);
            let k = jj * nx + ii;
            d.data[k] = d.data[k].max(s);
        }
        d
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::isostasy::IsostasyConfig;

    /// ADR Finding 152 G-field — building the zoning's structural density leaves C-3b's own density (the one the
    /// incision reads) unchanged, and the recorded run settles on the same state as the plain run.
    #[test]
    fn the_structural_density_never_touches_c3b() {
        let run = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: 40,
            dx: 1.0 / 32.0,
            dy: 1.0 / 32.0,
            iso_config: IsostasyConfig::c1_default(),
            drainage_max_distance: 30,
        };
        let (state, kin, hist) = record_history(7, 32, &Phase2InitParams::default(), &run, &C1Closures::default());
        let mut plain = init_c1_state_phase_2_r7(32, 7, &Phase2InitParams::default());
        let mut pk = PlateKinematics::preset_phase_1_1(plain.num_plates);
        crate::tectonics_c1::time_loop::run_with_closures(&mut plain, &mut pk, &run, &C1Closures::default(), |_, _| {});
        let bits = |s: &C1State| -> Vec<u64> { (0..s.nx() * s.ny()).map(|c| s.s.get(c % s.nx(), c / s.nx()).to_bits()).collect() };
        assert_eq!(bits(&state), bits(&plain), "the observed run is the plain run");
        let fcfg = FractureConfig { enabled: true, amplitude: 6.0, decay_km: 25.0, domain_km: 400.0, ..Default::default() };
        let before = derive_coarse_density(&state, &kin, &fcfg, None);
        let s = hist.structural_density(&state, &kin, &fcfg, false);
        let after = derive_coarse_density(&state, &kin, &fcfg, None);
        assert_eq!(before.data, after.data, "C-3b's density is untouched");
        assert!(s.data.iter().zip(&before.data).all(|(a, b)| a >= b), "the structural density is ≥ C-3b's everywhere");
        assert_eq!(hist.structural_density(&state, &kin, &fcfg, true).data, before.data, "the control is C-3b alone");
    }
}
