//! Post-run statistics on `S̃` by plate type and boundary flag.
//!
//! Step 5 introduces type-heterogeneous cells, so domain-level
//! `mean(S̃)` is no longer a useful descriptor. The harness instead
//! reports the per-class means and standard deviations below, used
//! both to feed the calibration loop (target
//! `s_oceanic_mean ∈ [0.18, 0.22]`) and to populate the physics
//! report's acceptance criteria.
//!
//! "Interior" continental cells are those with no boundary flag —
//! the `None` variant — distinguishing them from collision-band and
//! rift cells whose `S̃` evolves away from the reference 1.0.

use super::boundary_flag::{BoundaryFlag, BoundaryFlagField};
use super::plate_type::{PlateType, PlateTypeField};
use super::super::field::Field2D;

/// Mean + standard deviation aggregate.
#[derive(Clone, Copy, Debug, Default)]
pub struct MeanStd {
    pub mean: f64,
    pub std: f64,
    pub count: usize,
}

fn mean_std_of<'a, I: Iterator<Item = &'a f64>>(it: I) -> MeanStd {
    let values: Vec<f64> = it.copied().collect();
    let n = values.len();
    if n == 0 {
        return MeanStd::default();
    }
    let mean: f64 = values.iter().sum::<f64>() / n as f64;
    let var: f64 = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n as f64;
    MeanStd {
        mean,
        std: var.sqrt(),
        count: n,
    }
}

/// Mean + std of `S̃` over oceanic cells.
pub fn s_oceanic(s: &Field2D, plate_types: &PlateTypeField) -> MeanStd {
    let it = s
        .data()
        .iter()
        .zip(plate_types.data().iter())
        .filter_map(|(v, &t)| match t {
            PlateType::Oceanic => Some(v),
            _ => None,
        });
    mean_std_of(it)
}

/// Mean + std of `S̃` over continental cells with boundary_flag ==
/// None (strict interior). The acceptance target at Step 5 is this
/// band sitting in `[0.9, 1.1]`.
pub fn s_continental_interior(
    s: &Field2D,
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
) -> MeanStd {
    let mut it_iter = s
        .data()
        .iter()
        .zip(plate_types.data().iter())
        .zip(flags.data().iter())
        .filter_map(|((v, &t), &f)| match (t, f) {
            (PlateType::Continental, BoundaryFlag::None) => Some(v),
            _ => None,
        });
    let owned: Vec<&f64> = it_iter.by_ref().collect();
    mean_std_of(owned.into_iter())
}

/// Mean of `S̃` over continental collision cells. Reported as
/// telemetry only; the value grows with time (orogen thickening) and
/// has no Step-5 acceptance threshold.
pub fn s_continental_collision_mean(
    s: &Field2D,
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
) -> Option<f64> {
    let mut sum = 0.0_f64;
    let mut count = 0usize;
    for ((v, &t), &f) in s
        .data()
        .iter()
        .zip(plate_types.data().iter())
        .zip(flags.data().iter())
    {
        if matches!(
            (t, f),
            (PlateType::Continental, BoundaryFlag::ContinentalCollision)
        ) {
            sum += v;
            count += 1;
        }
    }
    if count == 0 {
        None
    } else {
        Some(sum / count as f64)
    }
}

/// Count of distinct boundary-mechanism variants with a non-zero
/// `Q_something` cell in the run. Subduction and OceanicSubduction
/// both count as the "subduction" mechanism — one unit, not two.
pub fn boundary_type_diversity(
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
    rates_nonzero: BoundaryMechanismActive,
) -> u32 {
    let mut subduction_present = false;
    let mut spread_present = false;
    let mut coll_present = false;
    let mut rift_v_present = false;
    let mut has_oceanic = false;
    let mut has_continental = false;

    for (&t, &f) in plate_types.data().iter().zip(flags.data().iter()) {
        match t {
            PlateType::Oceanic => has_oceanic = true,
            PlateType::Continental => has_continental = true,
        }
        match (t, f) {
            (_, BoundaryFlag::Subduction | BoundaryFlag::OceanicSubduction) => {
                subduction_present = true;
            }
            (PlateType::Oceanic, BoundaryFlag::Rift) => spread_present = true,
            (PlateType::Continental, BoundaryFlag::Rift) => rift_v_present = true,
            (PlateType::Continental, BoundaryFlag::ContinentalCollision) => {
                coll_present = true;
            }
            _ => {}
        }
    }

    let _ = has_oceanic;
    let _ = has_continental;

    let mut n = 0u32;
    if subduction_present && rates_nonzero.sub {
        n += 1;
    }
    if spread_present && rates_nonzero.spread {
        n += 1;
    }
    if coll_present && rates_nonzero.coll_v {
        n += 1;
    }
    if rift_v_present && rates_nonzero.rift_v {
        n += 1;
    }
    n
}

/// Helper to signal which rate coefficients are non-zero in the run.
/// `Q_arc` is not a distinct diversity unit — it re-publishes subducted
/// mass, so if `sub` is counted, arc comes with it.
#[derive(Clone, Copy, Debug)]
pub struct BoundaryMechanismActive {
    pub sub: bool,
    pub spread: bool,
    pub coll_v: bool,
    pub rift_v: bool,
}

impl From<&super::boundary_flag::BoundaryRates> for BoundaryMechanismActive {
    fn from(r: &super::boundary_flag::BoundaryRates) -> Self {
        Self {
            sub: r.k_sub != 0.0,
            spread: r.k_spread != 0.0,
            coll_v: r.k_coll_v != 0.0,
            rift_v: r.k_rift_v != 0.0,
        }
    }
}

impl BoundaryMechanismActive {
    /// Closed-mode constructor: in Step 6 Closed mode, rate
    /// coefficients (k_arc, k_spread, k_coll_v, k_rift_v) are
    /// typically zeroed because creation is driven by the
    /// recycling budget fractions, not per-cell rates. Use the
    /// fractions to decide which mechanisms are actually active.
    ///
    /// Note: `k_sub` drives the drain budget even in Closed mode,
    /// so it is taken from rates. The four creation channels
    /// (arc/spread/coll_v/rift_v) come from the config fractions.
    pub fn from_closed_mode(
        rates: &super::boundary_flag::BoundaryRates,
        recycling: &super::super::recycling::RecyclingConfig,
    ) -> Self {
        Self {
            sub: rates.k_sub != 0.0,
            spread: recycling.spread_fraction != 0.0,
            coll_v: recycling.coll_v_fraction != 0.0,
            rift_v: recycling.rift_v_fraction != 0.0,
        }
    }
}

/// Interface-cell mask: a cell is "on the oceanic/continental
/// interface" if its plate type differs from at least one 4-neighbour.
/// Returned as a boolean field the caller can iterate against `S̃`
/// gradient magnitudes.
pub fn interface_mask(
    plate_types: &PlateTypeField,
    idx_x: &super::super::field::PeriodicIndex,
    idx_y: &super::super::field::PeriodicIndex,
) -> Vec<bool> {
    let nx = plate_types.nx();
    let ny = plate_types.ny();
    let mut m = vec![false; nx * ny];
    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);
            let t = plate_types.get(i, j);
            let neighbours = [
                plate_types.get(ip, j),
                plate_types.get(im, j),
                plate_types.get(i, jp),
                plate_types.get(i, jm),
            ];
            if neighbours.iter().any(|&n| n != t) {
                m[j * nx + i] = true;
            }
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::field::PeriodicIndex;

    #[test]
    fn s_oceanic_stats_ignore_continental_cells() {
        let nx = 3;
        let ny = 1;
        let mut s = Field2D::new(nx, ny);
        s.set(0, 0, 0.2);
        s.set(1, 0, 0.2);
        s.set(2, 0, 1.0);
        let mut pt = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        pt.set(2, 0, PlateType::Continental);
        let st = s_oceanic(&s, &pt);
        assert_eq!(st.count, 2);
        assert!((st.mean - 0.2).abs() < 1e-14);
        assert!(st.std < 1e-14);
    }

    #[test]
    fn interior_continental_excludes_boundary_flags() {
        let nx = 3;
        let ny = 1;
        let mut s = Field2D::new(nx, ny);
        s.set(0, 0, 1.0);
        s.set(1, 0, 1.2);
        s.set(2, 0, 0.9);
        let pt = PlateTypeField::filled(nx, ny, PlateType::Continental);
        let mut fl = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        fl.set(1, 0, BoundaryFlag::ContinentalCollision);
        let st = s_continental_interior(&s, &pt, &fl);
        assert_eq!(st.count, 2);
        assert!((st.mean - 0.95).abs() < 1e-14);
    }

    #[test]
    fn diversity_counts_mechanisms_once() {
        // 3×1 with one oceanic subduction + one continental
        // collision cell. Two mechanisms active → diversity = 2.
        let nx = 3;
        let ny = 1;
        let mut pt = PlateTypeField::filled(nx, ny, PlateType::Continental);
        pt.set(0, 0, PlateType::Oceanic);
        let mut fl = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        fl.set(0, 0, BoundaryFlag::OceanicSubduction);
        fl.set(1, 0, BoundaryFlag::ContinentalCollision);
        let active = BoundaryMechanismActive {
            sub: true,
            spread: true,
            coll_v: true,
            rift_v: true,
        };
        let d = boundary_type_diversity(&pt, &fl, active);
        assert_eq!(d, 2);
    }

    #[test]
    fn interface_mask_fires_on_oceanic_continental_boundary() {
        let nx = 4;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut pt = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        for i in 0..nx {
            pt.set(i, 3, PlateType::Continental);
        }
        let m = interface_mask(&pt, &idx_x, &idx_y);
        // Row j=3 is continental, row j=2 is oceanic; interface on
        // both rows (across j=2↔3 faces, plus wrap across j=3↔0).
        for i in 0..nx {
            assert!(m[3 * nx + i], "continental row should be interface");
            assert!(m[2 * nx + i], "row adjacent to continental should be interface");
            assert!(m[0 * nx + i], "wrap neighbour of continental should be interface");
        }
        // Row j=1 is fully oceanic with oceanic neighbours → not interface.
        for i in 0..nx {
            assert!(!m[1 * nx + i], "oceanic interior should not be interface");
        }
    }
}
