//! Geology v1 (ADR Findings 150–152): the ROCK grid of a world and the FAVOURABILITY of its resources.
//!
//! - The rock grid is fixed for a world. The favourabilities come from a declarative rules file (TOML); reloading the
//!   file reruns the zoning only.
//! - The values are RELATIVE favourabilities (0–100), never deposits nor probabilities: « Favorabilité dans Ymir,
//!   rareté dans Living Landz ».
//! - The stage is READ-ONLY: it reads the eroded world and writes only its own buffers (G-field). Its structural density
//!   never reaches the incision.
//! - Spec: `docs/reports/geology_v1/f151_spec/spec_geologie_v1.md`; format: `docs/geology_format.md`.

pub mod export;
pub mod history;
pub mod render;
pub mod rocks;
pub mod rules;
pub mod zoning;

use std::collections::HashSet;
use std::time::Instant;

use crate::grid::GridF32;
use crate::tectonics_c1::closures::fracture::FractureConfig;
use crate::seed::WorldSeed;
use crate::tectonics_c1::closures::volcanism::{Edifice, VolcanismConfig, place_edifices};
use crate::tectonics_c1::debug_labels::{CoarseTectonicLabels, derive_tectonic_labels};
use crate::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, C1DrainageResult, LakeType};
use crate::tectonics_c1::init_r7::Phase2InitParams;
use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};

pub use history::TectonicHistory;
pub use rules::{GeologyRules, parse_rules, sha256_hex};
pub use zoning::{DensitySource, Zoning, ZoningContext, zone};

/// The world's inputs to the geology stage.
pub struct GeologyInputs<'a> {
    pub seed: u64,
    pub grid: usize,
    pub init: &'a Phase2InitParams,
    pub run: &'a C1TimeLoopConfig,
    pub closures: &'a C1Closures,
    /// The final (conditioned) field, normalised, and the metres per HD cell.
    pub field: &'a GridF32,
    pub z_m: &'a [f32],
    pub cell_km: f32,
    pub drainage: &'a C1DrainageResult,
    pub wetland: &'a [u8],
    pub temp_c: &'a [f32],
    pub precip_mm: &'a [f32],
    pub sample_origin: [f64; 2],
    pub sample_size: f64,
    pub domain_km: f32,
    /// C-2's config: the edifices are placed from the history's settled state (`None`: no volcanism).
    pub volcanism: Option<&'a VolcanismConfig>,
}

/// The seconds of each part of the stage.
#[derive(Clone, Copy, Debug, Default)]
pub struct GeologyTimings {
    pub history_s: f64,
    pub rocks_s: f64,
    pub context_s: f64,
    pub zoning_s: f64,
}

/// The geology of a world: the rock grid, the zoning context (kept for re-zoning), the history and the zoning.
#[derive(Clone, Debug)]
pub struct GeologyProduct {
    pub w: usize,
    pub h: usize,
    pub rocks: Vec<u8>,
    pub context: ZoningContext,
    pub history: TectonicHistory,
    /// Coarse masks (nx·ny): the fossil and active belts.
    pub fossil: Vec<bool>,
    pub active: Vec<bool>,
    pub edifices: Vec<Edifice>,
    pub timings: GeologyTimings,
}

/// ADR Finding 154-S — what the SUBSTRATUM needs from the tectonics: the coarse history's end state, its masks
/// (craton, belt = fossil or active, rift) and the edifices. It needs C1 alone (no relief, no drainage), so it is
/// available from the construction stage on.
pub struct TectonicSources {
    pub state: C1State,
    pub kin: PlateKinematics,
    pub history: TectonicHistory,
    pub end: CoarseTectonicLabels,
    pub nx: usize,
    pub ny: usize,
    pub fossil: Vec<bool>,
    pub active: Vec<bool>,
    pub craton: Vec<bool>,
    pub belt: Vec<bool>,
    pub edifices: Vec<Edifice>,
}

/// ADR Finding 154-S — run the coarse history and place the edifices (C-2's config; `None`: no volcanism).
#[must_use]
pub fn tectonic_sources(seed: u64, grid: usize, init: &Phase2InitParams, run: &C1TimeLoopConfig, closures: &C1Closures, volcanism: Option<&VolcanismConfig>) -> TectonicSources {
    let (state, kin, history) = history::record_history(seed, grid, init, run, closures);
    let end = derive_tectonic_labels(&state, &kin);
    let fossil = history.fossil_belt(&state, &end);
    let active = history.active_belt(&end);
    let (nx, ny) = (state.nx(), state.ny());
    let craton: Vec<bool> = (0..nx * ny).map(|c| state.cratonic_mask.get(c % nx, c / nx)).collect();
    let belt: Vec<bool> = fossil.iter().zip(&active).map(|(a, b)| *a || *b).collect();
    let edifices = match volcanism {
        Some(v) if v.enabled => place_edifices(&state, &kin, &WorldSeed::new(seed), v.domain_km, v),
        _ => Vec::new(),
    };
    TectonicSources { state, kin, history, end, nx, ny, fossil, active, craton, belt, edifices }
}

impl TectonicSources {
    /// ADR Finding 154-S — the substratum on the HD grid (`w·h`), read on the relief `z_m` (m).
    #[must_use]
    pub fn substratum(&self, z_m: &[f32], w: usize, h: usize, cell_km: f32, sample_origin: [f64; 2], sample_size: f64) -> Vec<u8> {
        rocks::build_substratum(&rocks::SubstratumInputs {
            w,
            h,
            z_m,
            cell_km,
            nx: self.nx,
            ny: self.ny,
            craton: &self.craton,
            belt: &self.belt,
            rift: &self.end.rift,
            sample_origin,
            sample_size,
            edifices: &self.edifices,
        })
    }
}

/// Build the rock grid and the zoning context (no rules involved).
#[must_use]
pub fn build_geology(inp: &GeologyInputs) -> GeologyProduct {
    let t = Instant::now();
    let src = tectonic_sources(inp.seed, inp.grid, inp.init, inp.run, inp.closures, inp.volcanism);
    let (nx, ny) = (src.nx, src.ny);
    let fcfg = FractureConfig { enabled: true, amplitude: 6.0, decay_km: 25.0, domain_km: inp.domain_km, ..Default::default() };
    let density = src.history.structural_density(&src.state, &src.kin, &fcfg, false);
    let density_c3b = src.history.structural_density(&src.state, &src.kin, &fcfg, true);
    let history_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let (w, h) = (inp.field.width, inp.field.height);
    let dr = inp.drainage;
    let sea: Vec<bool> = (0..w * h).map(|c| dr.lake_map[c] == 0 && inp.field.data[c] <= C1_SEA_LEVEL_NORM).collect();
    let endorheic: HashSet<u32> = dr.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).map(|l| l.base.id).collect();
    let slope = rocks::slope_field(inp.z_m, w, h, inp.cell_km * 1000.0);
    // ADR Finding 154-S -- the surface over the substratum (the substratum read on the final relief: F152's grid)
    let substratum = src.substratum(inp.z_m, w, h, inp.cell_km, inp.sample_origin, inp.sample_size);
    let rock_grid = rocks::apply_surface(
        &substratum,
        &rocks::SurfaceInputs {
            w,
            h,
            sea: &sea,
            lake_map: &dr.lake_map,
            endorheic: &endorheic,
            acc: &dr.flow.accumulation.data,
            cell_km: inp.cell_km,
            precip_mm: inp.precip_mm,
        },
        &slope,
    );
    let rocks_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let upper: Vec<bool> = src.end.subduction_upper.clone();
    let context = zoning::build_context(&zoning::ContextInputs {
        w,
        h,
        cell_km: inp.cell_km,
        rocks: &rock_grid,
        slope: &slope,
        acc: &dr.flow.accumulation.data,
        dir: &dr.flow.direction,
        lake_map: &dr.lake_map,
        endorheic: &endorheic,
        wetland: inp.wetland,
        temp_c: inp.temp_c,
        precip_mm: inp.precip_mm,
        sea: &sea,
        nx,
        ny,
        fossil: &src.fossil,
        active: &src.active,
        upper_plate: &upper,
        density: &density,
        density_c3b: &density_c3b,
        sample_origin: inp.sample_origin,
        sample_size: inp.sample_size,
        edifices: &src.edifices,
    });
    let context_s = t.elapsed().as_secs_f64();
    let TectonicSources { history, fossil, active, edifices, .. } = src;
    GeologyProduct { w, h, rocks: rock_grid, context, history, fossil, active, edifices, timings: GeologyTimings { history_s, rocks_s, context_s, zoning_s: 0.0 } }
}

/// Zone a product with a rules file, timed.
#[must_use]
pub fn zone_timed(p: &GeologyProduct, rules: &GeologyRules) -> (Zoning, f64) {
    let t = Instant::now();
    let z = zone(&p.context, rules, DensitySource::Structural);
    (z, t.elapsed().as_secs_f64())
}
