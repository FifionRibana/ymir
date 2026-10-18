//! ADR Finding 180-A2 -- **the per-seed mass bound on Davis-Suppe's increments** (Dahlen 1990, p. 59, eq. 1: the wedge's
//! cross-section grows by exactly hV per unit length, whatever the erosion; Willett et al. 2001 p. 456, 459).
//!
//! Each Davis-Suppe cell is attached to its nearest seed (Davis-Suppe's metric: 8-connected, √2 diagonals, within the plate, up
//! to its reach). Per seed σ and per step, the demand R_σ is the sum of its cells' increments and the supply
//! S_σ = (v_rel · n̂) · dt · 64 · s_in is the section hV of one cell length of trench; an incoming or trench cell feeding k seeds
//! gives each S / k. If R_σ > S_σ every increment of σ is multiplied by S_σ / R_σ (the wedge's shape kept, no iteration); a seed
//! with R_σ ≤ S_σ keeps its increments bit for bit. s_in is the incoming cell's whole s (an upper bound: C1 has no sediment
//! layer). Declared in `docs/reports/relief_method/f180_andean_mass/f180_declared.md`.

use std::collections::BinaryHeap;

use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::voronoi::PlateIdField;

use super::kinematics::PlateKinematics;
use super::state::BoolField;

/// The bound's record of one step.
#[derive(Clone, Debug, Default)]
pub struct DsBoundStats {
    /// The seeds, and those whose increments were scaled.
    pub seeds: usize,
    pub bounded: usize,
    /// The mass removed per cell (s units), row-major (empty when the bound is off).
    pub removed: Vec<f64>,
    /// (seed cell, scaled) for every seed.
    pub seed_bounded: Vec<(usize, bool)>,
    /// The largest Σ increments / S_σ after the bound (≤ 1 when the bound holds; 0 with no seed).
    pub max_ratio: f64,
}

/// ADR Finding 180-A2 -- bound the step's Davis-Suppe increments (`s_after` − `s_before`) seed by seed. `upper` holds every seed
/// of the step (the boundary seeds and the Andean ones); `andean_feeds` maps an Andean seed to its trench cell and `partners`
/// gives a trench cell's sinking partner and v_rel · n̂. `supply_factor` multiplies every S_σ (1 declared; ∞ = no bound).
#[allow(clippy::too_many_arguments)]
pub fn apply_ds_mass_bound(
    s_before: &[f64],
    s_after: &mut Field2D,
    plate_id: &PlateIdField,
    upper: &BoolField,
    kinematics: &PlateKinematics,
    dt: f64,
    max_distance: f64,
    andean_feeds: &[(usize, usize)],
    partners: &[Option<(usize, f64)>],
    supply_factor: f64,
) -> DsBoundStats {
    let (nx, ny) = (plate_id.nx(), plate_id.ny());
    let n = nx * ny;
    let pid = |k: usize| plate_id.get(k % nx, k / nx);
    let cells_per_step = dt * nx as f64;
    // the seeds' supplies, shared among the seeds a feeding cell feeds
    let seeds: Vec<usize> = (0..n).filter(|&k| upper.get(k % nx, k / nx)).collect();
    let mut is_andean = vec![false; n];
    let mut feed_of = vec![usize::MAX; n];
    for &(sd, tr) in andean_feeds {
        is_andean[sd] = true;
        feed_of[sd] = tr;
    }
    // (feeding cell, raw supply) per seed
    let mut src: Vec<(usize, f64)> = Vec::with_capacity(seeds.len());
    for &k in &seeds {
        if is_andean[k] {
            let tr = feed_of[k];
            match partners.get(tr).copied().flatten() {
                Some((q, dot)) => src.push((tr, dot * cells_per_step * s_before[q])),
                None => src.push((usize::MAX, 0.0)),
            }
        } else {
            let (x, y) = (k % nx, k / nx);
            let (vxc, vyc) = kinematics.velocities[pid(k) as usize];
            let mut best: Option<(usize, f64)> = None;
            for (q, nx_n, ny_n) in [
                (y * nx + (x + 1) % nx, 1.0, 0.0),
                (y * nx + (x + nx - 1) % nx, -1.0, 0.0),
                (((y + 1) % ny) * nx + x, 0.0, 1.0),
                (((y + ny - 1) % ny) * nx + x, 0.0, -1.0),
            ] {
                if pid(q) == pid(k) {
                    continue;
                }
                let (vxn, vyn) = kinematics.velocities[pid(q) as usize];
                let dot = (vxc - vxn) * nx_n + (vyc - vyn) * ny_n;
                if dot > 0.0 && best.is_none_or(|b| dot > b.1) {
                    best = Some((q, dot));
                }
            }
            match best {
                Some((q, dot)) => src.push((n + q, dot * cells_per_step * s_before[q])),
                None => src.push((usize::MAX, 0.0)),
            }
        }
    }
    let mut share: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    for &(f, _) in &src {
        if f != usize::MAX {
            *share.entry(f).or_insert(0) += 1;
        }
    }
    let supply: Vec<f64> = src.iter().map(|&(f, s)| if f == usize::MAX { 0.0 } else { supply_factor * s / share[&f] as f64 }).collect();
    // the nearest seed of every cell (Davis-Suppe's metric, within the plate, up to its reach)
    let mut dist = vec![f64::INFINITY; n];
    let mut owner = vec![usize::MAX; n];
    let mut heap = BinaryHeap::new();
    for (si, &k) in seeds.iter().enumerate() {
        dist[k] = 0.0;
        owner[k] = si;
        heap.push(std::cmp::Reverse((0u64, si, k)));
    }
    let diag = std::f64::consts::SQRT_2;
    while let Some(std::cmp::Reverse((dk, si, k))) = heap.pop() {
        let d = dk as f64 / 1e9;
        if owner[k] != si || d > dist[k] + 1e-12 || d >= max_distance {
            continue;
        }
        let (x, y) = (k % nx, k / nx);
        for (dx, dy, c) in [(1i64, 0i64, 1.0), (-1, 0, 1.0), (0, 1, 1.0), (0, -1, 1.0), (1, 1, diag), (1, -1, diag), (-1, 1, diag), (-1, -1, diag)] {
            let j = ((y as i64 + dy).rem_euclid(ny as i64) as usize) * nx + (x as i64 + dx).rem_euclid(nx as i64) as usize;
            if pid(j) != pid(k) {
                continue;
            }
            let nd = d + c;
            if nd < dist[j] - 1e-12 {
                dist[j] = nd;
                owner[j] = si;
                heap.push(std::cmp::Reverse(((nd * 1e9) as u64, si, j)));
            }
        }
    }
    // the demand per seed, the scaling
    let mut demand = vec![0f64; seeds.len()];
    for k in 0..n {
        let inc = s_after.data()[k] - s_before[k];
        if inc > 0.0 && owner[k] != usize::MAX {
            demand[owner[k]] += inc;
        }
    }
    let mut stats = DsBoundStats { seeds: seeds.len(), removed: vec![0.0; n], ..Default::default() };
    let scale: Vec<Option<f64>> = (0..seeds.len()).map(|si| (demand[si] > supply[si]).then(|| supply[si] / demand[si])).collect();
    for k in 0..n {
        let inc = s_after.data()[k] - s_before[k];
        if inc <= 0.0 || owner[k] == usize::MAX {
            continue;
        }
        if let Some(f) = scale[owner[k]] {
            let new = s_before[k] + inc * f;
            stats.removed[k] = inc - inc * f;
            s_after.set(k % nx, k / nx, new);
        }
    }
    let mut after = vec![0f64; seeds.len()];
    for k in 0..n {
        let inc = s_after.data()[k] - s_before[k];
        if inc > 0.0 && owner[k] != usize::MAX {
            after[owner[k]] += inc;
        }
    }
    for si in 0..seeds.len() {
        let b = scale[si].is_some();
        stats.bounded += b as usize;
        stats.seed_bounded.push((seeds[si], b));
        if after[si] > 0.0 {
            let r = if supply[si] > 0.0 { after[si] / supply[si] } else { f64::INFINITY };
            stats.max_ratio = stats.max_ratio.max(r);
        }
    }
    stats
}
