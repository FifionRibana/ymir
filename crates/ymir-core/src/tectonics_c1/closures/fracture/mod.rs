//! C-3b inherited structure (closures roadmap §3b) — an ISOTROPIC erodibility
//! modulation by tectonic FRACTURE DENSITY, derived causally from the tectonic state.
//! Intact cratonic interior retains its relief (K = reference ×1); tectonically active
//! belts near plate contacts and accretion sutures are more fractured and erode more.
//!
//! ## Scope — density only (the orientation was measured out)
//!
//! C-3b began as a DIRECTIONAL fracture field (density + orientation) meant to align
//! valleys on the tectonic fabric via anisotropic incision. That mechanism was BUILT
//! and MEASURED and it does not work: with the global erosion rate held constant
//! (mean-preserving anisotropic K), fabric alignment |flow·strike| still FELL
//! (0.639 → 0.536, both resolutions) and closed depressions rose (1001 → 2146) — the
//! incision RATE cannot re-route a receiver fixed by topography, so it cannot align
//! valleys. The true mechanism is flow ROUTING, but (a) its blast radius covers
//! `compute_flow` → C-1, rivers, lakes, the whole hydro chain, and (b) C1's directional
//! field is too poor to feed it (constant per-plate velocities, no strain, no
//! deformation history → a continent-wide uniform grain = an artefact), and (c) the
//! premise is weak (the Appalachian trellis is FOLDED STRATA, out of scope, not
//! fractures). So C-3b ships the DENSITY only; the orientation limitation and its
//! specification are recorded in ADR 0001 (C-3b).
//!
//! ## Literature (density → erodibility)
//!
//! Molnar 2007 (*Tectonics, fracturing of rock, and erosion*): tectonics erodes mostly
//! by FRACTURING rock — fragmentation → plucking, avenues for water. Clarke & Burbank
//! 2011; Zondervan et al.: fracturing spans ~1–2 orders of K, and homogenises the
//! inter-lithology contrast by ~1 order. Domain of validity: bedrock rivers, brittle
//! upper crust (< ~10 km), detachment-limited — the relief-v3 production regime.
//!
//! ## What is derived (causally, never noise / geometry)
//!
//! Density = proximity to real tectonic CONTACTS: the dynamic boundary classification
//! (`classify_boundaries`, from `plate_id` + kinematics — NOT `cratonic_mask`, the
//! FBM-noise-refined field C-3 rejected, and NOT the geometric craton placeholder) and
//! the accretion SUTURES (Phase B). `density → 0` far from every contact, so the intact
//! cratonic interior EMERGES at K = 1 (the reference — global-slowdown nil by
//! construction, the C-3 design that must survive), rising toward 1 at a contact.
//! Limit: a mobile belt far from any CURRENT boundary reads as craton (C1 records no
//! past belts) — a weaker but CAUSAL contrast, preferred over a broad noise-derived one.

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;
use crate::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries};
use crate::tectonics_c1::closures::lithology::upscale_k_to_hd;
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;

/// Version of the fracture derivation. ⚠️ BUMP on any change that moves the field, so
/// a fracture-enabled eroded cache invalidates. Added to the eroded key only when on.
pub const FRACTURE_ALGO: u32 = 2;

/// C-3b fracture config. `enabled = false` (default) → no field, pre-C-3b byte-identical.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct FractureConfig {
    /// Master switch. Default `false` (OFF → byte-identical).
    pub enabled: bool,
    /// Erodibility amplitude at a tectonic contact: `K = 1 + amplitude · density`,
    /// `density ∈ [0,1]` (0 = intact craton → K = 1 reference). Swept, not fixed
    /// (Molnar's ~1–2 orders is the outer bound).
    pub amplitude: f32,
    /// Extra density multiplier ON suture bands (Phase B; accretion welds are linear
    /// zones of weakness). `density_suture = min(1, suture_multiplier · contact_density)`.
    pub suture_multiplier: f32,
    /// Density decay length (km) from a tectonic contact: `density = exp(-dist/decay)`.
    /// Sets how wide the fractured belt around a boundary is.
    pub decay_km: f32,
    /// Geometric domain span (km) — sets km per coarse cell. NEVER geo_scale_ratio.
    pub domain_km: f32,
}

impl Default for FractureConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            amplitude: 4.0,
            suture_multiplier: 2.0,
            // Narrow orogenic fracture belt (~40–80 km): keeps the intact craton the
            // MAJORITY (~53 % of land) so it reads as the ×1 reference, with the
            // fractured belt a minority (measured; see c3b_fracture_sweep).
            decay_km: 25.0,
            domain_km: 400.0,
        }
    }
}

impl FractureConfig {
    /// `true` when OFF — for serde `skip_serializing_if` so a disabled config drops out
    /// of the eroded cache key (byte-identical to pre-C-3b).
    #[must_use]
    pub fn is_disabled(&self) -> bool {
        !self.enabled
    }
}

/// Octile (8-neighbour) distance in CELLS to the nearest `true` cell of `seed`,
/// periodic, capped at `cap`. Small BFS-relaxation (64² → sub-ms).
fn distance_to_mask(seed: &[bool], nx: usize, ny: usize, cap: f32) -> Vec<f32> {
    let diag = std::f32::consts::SQRT_2;
    let mut dist = vec![cap; nx * ny];
    let mut front = Vec::new();
    for k in 0..nx * ny {
        if seed[k] {
            dist[k] = 0.0;
            front.push(k);
        }
    }
    let nbr: [(i32, i32, f32); 8] = [
        (1, 0, 1.0),
        (-1, 0, 1.0),
        (0, 1, 1.0),
        (0, -1, 1.0),
        (1, 1, diag),
        (1, -1, diag),
        (-1, 1, diag),
        (-1, -1, diag),
    ];
    while let Some(k) = front.pop() {
        let (x, y) = ((k % nx) as i32, (k / nx) as i32);
        let dk = dist[k];
        for &(dx, dy, c) in &nbr {
            let nxp = (x + dx).rem_euclid(nx as i32) as usize;
            let nyp = (y + dy).rem_euclid(ny as i32) as usize;
            let nk = nyp * nx + nxp;
            let nd = dk + c;
            if nd < dist[nk] && nd < cap {
                dist[nk] = nd;
                front.push(nk);
            }
        }
    }
    dist
}

/// Coarse fracture DENSITY field ∈ [0,1]: `exp(-dist_to_contact_km / decay)`. `1` at a
/// plate boundary (or suture, Phase B), `→ 0` in the intact interior (K = reference).
/// `suture_mask` (Phase B; `None` for now) adds seed cells with the suture weight.
#[must_use]
pub fn derive_coarse_density(
    state: &C1State,
    kin: &PlateKinematics,
    cfg: &FractureConfig,
    suture_mask: Option<&[bool]>,
) -> GridF32 {
    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let info = classify_boundaries(&state.plate_id, kin);
    // Contacts that FRACTURE the basement: convergent orogens + transform faults
    // (compression + shear). Divergent boundaries are rifts — that softness is C-3's
    // domain (rift-soft), not a fractured-hard belt — so they are excluded here, which
    // also keeps the belt a minority so the intact craton emerges as the majority.
    let mut contact = vec![false; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            if matches!(
                info.boundary_type.get(i, j),
                BoundaryType::Convergent | BoundaryType::Transform
            ) {
                contact[j * nx + i] = true;
            }
        }
    }
    let km_per_cell = cfg.domain_km / nx as f32;
    let cap = (nx.max(ny) as f32) * 0.75;
    let dist = distance_to_mask(&contact, nx, ny, cap);
    // Suture distance (Phase B): a separate, closer-weighted contribution.
    let sdist = suture_mask.map(|m| distance_to_mask(m, nx, ny, cap));

    let mut density = GridF32::new(nx, ny, 0.0);
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let d_km = dist[k] * km_per_cell;
            let mut dens = (-d_km / cfg.decay_km.max(1e-3)).exp();
            if let Some(sd) = &sdist {
                let sd_km = sd[k] * km_per_cell;
                let sdens = cfg.suture_multiplier * (-sd_km / cfg.decay_km.max(1e-3)).exp();
                dens = dens.max(sdens);
            }
            density.set(i, j, dens.clamp(0.0, 1.0));
        }
    }
    density
}

/// Build the HD per-cell erodibility MULTIPLIER `K = 1 + amplitude · density`,
/// registered to the terrain via `(sample_origin, sample_size)`. `1.0` in the intact
/// interior (the hard-basement reference), rising toward contacts / sutures.
#[must_use]
pub fn build_hd_density_k(
    state: &C1State,
    kin: &PlateKinematics,
    cfg: &FractureConfig,
    suture_mask: Option<&[bool]>,
    dst_w: usize,
    dst_h: usize,
    sample_origin: [f64; 2],
    sample_size: f64,
) -> Vec<f32> {
    let density = derive_coarse_density(state, kin, cfg, suture_mask);
    let dens_hd = upscale_k_to_hd(&density, dst_w, dst_h, sample_origin, sample_size);
    dens_hd.iter().map(|&d| 1.0 + cfg.amplitude * d.clamp(0.0, 1.0)).collect()
}
