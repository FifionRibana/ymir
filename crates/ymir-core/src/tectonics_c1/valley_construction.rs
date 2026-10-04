//! ADR 0001 Finding 121 — **valley construction: the SHAPE stage of the hybrid** (construct the
//! mature form along the trunks, then let a light incision supply the texture).
//!
//! **A bench / viz seam. `None` in production, byte-identical, nothing promoted.**
//!
//! ## Why this is not Option A again
//!
//! The C1 design document (`docs/c1_lightweight_dynamic_tectonics.md` §3.1) rejected "static with
//! invented rules — geometric templates modulated by procedural noise" for three reasons, copied
//! verbatim so the hybrid answers each:
//!
//! 1. *"cannot produce coherent chronological diversity (young vs old mountain ranges on the same
//!    continent)"*;
//! 2. *"parameter calibration drifts indefinitely with no physical anchor"*;
//! 3. *"'looks invented' is a real failure mode for users who can identify imposed geometry by
//!    eye"*.
//!
//! And ADR Finding 96 built the mature form as a FIELD (χ everywhere) and rejected it: 7.1× rougher
//! than the delivered field, cliffs across every divide, D8's axes stamped into the relief. This
//! module differs on each count: the valleys are laid ONLY along trunks of `A ≥ a_min_km2` (Finding
//! 120: 10 km²), the interfluves keep the tectonic + FBM field (`min(field, valley)`, so no divide
//! is ever built), and the cross-section is laid by the distance to a SMOOTHED skeleton, not to the
//! D8 cells.
//!
//! ## The laws, labelled
//!
//! | law | form | value | status |
//! |---|---|---|---|
//! | floor profile | `z = z_base + k·χ`, `χ = ∫ (A0/A)^0.5 dx`, A0 = 1 km², A clamped at `a_c_km2` | θ = 0.5 | ANCHORED — Harel et al. 2016 eqs. (4)–(7), PDF p. 15 |
//! | age `k` | one number, m/m at A0 | calibrated so the floor surface's paired p50 = 488 m (Finding 95's oracle; Finding 96's own calibration) | **PROXY** on a target |
//! | valley floor width | `W = a·A^0.3`, A in km², W in m | exponent 0.3 | exponent ANCHORED — Clubb et al. 2022 p. 437 |
//! | width coefficient `a` | | `W(10 km²) = 200 m` | **PROXY** (Clubb Table 2 p. 452 is unit-inconsistent, Finding 120) |
//! | walls | planar, 28° | | ANCHORED — Whipple & Tucker 1999 Table 1 p. 17,663 (colluvial slopes) |
//! | base | sea + 0.5 m | | ANCHORED — ADR Finding 83 |
//! | skeleton smoothing | moving average over ±`smooth_m` | 250 m | **PROXY** — against Finding 96's D8 stamping |
//! | trunk threshold | `A ≥ 10 km²` (4 194 cells at 8192²) | | DECISION — Finding 120 |
//!
//! Everything is in physical units: km², m, degrees. Nothing is expressed in cells or in Strahler
//! order (the author's amendment to Finding 120).

use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, C1DrainageConfig, c1_drainage_windowed};
use crate::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, c1_metres_to_altitude_norm,
};
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, breach_monotone, compute_flow};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// **PROXY — the age `k` calibrated by `f121_hybrid`** (seed 1, 8192², 400 km, at the CANONICAL
/// framing since Finding 124 — 0.07186 at Findings 121–123's framing): the floor surface
/// `base + k·χ` over the paired land of the pre-incision field has p50 = **488.1 m**, Finding 95's
/// oracle. Finding 96 calibrated its χ field on the same oracle and found 0.0719 with `A` clamped
/// at `A_c`; this skeleton gives 0.07186. A calibration on a TARGET, not a measurement of an age.
pub const F121_AGE_K: f32 = 0.07183;

/// The construction's parameters, all in physical units. See the module table for each label.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ValleyConstruction {
    /// Trunk threshold, km² of D8 drained area (in-domain, never the geo-ratio-signified area).
    /// DECISION (Finding 120).
    pub a_min_km2: f32,
    /// The age `k` of `z = z_base + k·χ` (m per m of χ, A0 = 1 km²). PROXY.
    pub age_k: f32,
    /// The drained area χ is clamped at (km²): Finding 96 B1′, the fluvial law is not applied
    /// below its own channel threshold. ANCHORED (Whipple & Tucker Table 1, via the anchoring
    /// report).
    pub a_c_km2: f32,
    /// `a` in `W = a·A^b` (m, A in km²). PROXY.
    pub width_coef_m: f32,
    /// `b` in `W = a·A^b`. ANCHORED (Clubb 2022 p. 437).
    pub width_exp: f32,
    /// Wall slope, degrees. ANCHORED (Whipple & Tucker 1999 p. 17,663).
    pub wall_deg: f32,
    /// Base level above the sea, m. ANCHORED (Finding 83).
    pub base_m: f32,
    /// Half-length of the skeleton smoothing, m. PROXY.
    pub smooth_m: f32,
    /// `None` = the BARE construction (no incision at all). `Some(f)` = ONE incision pass at
    /// `f × k_time` of the shipped stream power, AFTER the construction (the hybrid).
    pub light_k_time_fraction: Option<f32>,
    /// **ADR Finding 122-B — χ from the col of a closed depression.** `false` = Finding 121 (χ from
    /// the sea everywhere, byte-identical). `true`: a trunk whose D8 path enters a closed depression
    /// of the input field — an inland below-sea basin (`water_class` 2, Finding 85's spill level) or
    /// a lake the pre-incision drainage holds flat — takes that depression's SPILL LEVEL as its base
    /// and counts χ from the depression's col, so a constructed floor never lies below the water
    /// the depression will hold. Outside closed depressions nothing changes. Finding 122-A attributed
    /// the author's "brush" lake to exactly that: valleys built down to sea + 0.5 m inside a basin
    /// whose col stays at 459 m.
    #[serde(default = "basin_base_default")]
    pub basin_base: bool,
    /// **ADR Finding 124-B1 — break the plane.** `None` = the planar 28° wall (Findings 121–123,
    /// byte-identical). `Some` gives the wall a hillslope PROFILE (concave foot, straight colluvial
    /// segment, convex crest — the form Whipple & Tucker 1999 p. 17,663 describe for colluvial
    /// slopes) and puts back the tectonic + FBM detail the wall replaced. Finding 123's rule 18: the
    /// planar walls create the comb's teeth (4.6 % have an A1+B2 counterpart).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_profile: Option<WallProfile>,
    /// **ADR Finding 124-5 (E) — W(k), the valley widens with age.** `None` = Findings 121–124's
    /// width, independent of k (byte-identical). `Some(γ)`: `W = a₀·(k/k₀)^γ·A^b` with
    /// `k₀ = F121_AGE_K`, so the age knob that deepens the floors also widens them. **γ is a PROXY**:
    /// lateral widening is a matter of TIME, and Clubb 2022's `a` is a snapshot of present valleys,
    /// not a rate — nothing anchors γ. `Some(0.0)` is `None`'s width.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub width_age_gamma: Option<f32>,
    /// **ADR Finding 126-B — mur ↔ mer, the coast clause.** `None` = Findings 121–125: nothing
    /// floors a WALL cell, so the profile's foot, its crest and its signed detail can lay a wall
    /// below the sea (Finding 126-A: the bathymetry then reposes it at −1.00 m, and an enclosed one
    /// fires Finding 38's invariant). `Some(ε)`: no wall cell is laid below `sea + ε` m — the floors'
    /// own base (Finding 83's `base_m`) transposed to the walls, so a wall that reaches the coast
    /// stops at the coast. The clause's other half ("nothing raises a coastal cell above the
    /// original terrain") needs no code: `carve` only lowers (`out = min(field, V)`, asserted).
    /// DECISION, on the ANCHORED `base_m`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_sea_floor_m: Option<f32>,
    /// **ADR Finding 126-D — the skeleton smoothed in units of the valley's width.** `None` = the
    /// ±`smooth_m` moving average (Findings 121–125: a half-window pinned at each end, and a line
    /// shorter than `2·half + 3` samples left exactly as D8 drew it). `Some(c)`: at every sample the
    /// half-length is `c·W(A)` of THAT sample — a law in W, never in cells — and the window shrinks
    /// symmetrically toward the ends, so only the two endpoints stay pinned and every line is
    /// smoothed. **PROXY**.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub smooth_w: Option<f32>,
    /// **ADR Finding 127-B — the off-grid tracé.** `None` = Findings 120–126: every trunk polyline is
    /// its D8 cell chain. `Some(t)`: the MEDIAN LINE of each segment with `a_min_km2 ≤ A <
    /// t.below_km2` is retraced by continuous steepest descent on the INPUT field (interpolated by
    /// `t.interp`), half a cell per step, from the segment's source down to the entry of its D8
    /// receiver's corridor (± W/2). The accumulation, the junctions and the choice of trunks stay D8
    /// (the concentration, Finding 112). Where the gradient is below `t.flat_slope` the step falls back
    /// on the D8 path, declared and counted; a trace that reaches another corridor, the sea, or runs out
    /// of steps keeps its D8 geometry, counted. Not Finding 113's snap (that moved an EXPORTED river
    /// toward a thalweg): this replaces the skeleton's geometry before the construction digs, at a
    /// continuous angle. **PROXY** (`flat_slope`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skeleton_trace: Option<SkeletonTrace>,
    /// **ADR Finding 128-A — the confluence clause.** `false` = Findings 121–127: every cell is laid
    /// from its NEAREST sample, with a minimum across different lines only among the 8 neighbours'
    /// samples. At a confluence a tributary's sample can then be nearer to a trunk-floor cell than
    /// the trunk's own, and its higher floor dams the trunk (Finding 127-C: col +164 / +374 m). `true`:
    /// on the FLOOR BAND of a primitive (within `W/2` of one of its samples), the LINE with the largest
    /// drained area covering the cell lays it, from ITS nearest covering sample. The trunk wins its
    /// band. **Not a minimum** (Finding 121 refuted the minimum of the cones, which flattens the long
    /// profile), and never within one line (a line's own downstream sample never takes over its
    /// upstream band). DECISION.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub trunk_band: bool,
    /// **ADR Finding 128-C — the skeleton's directions by D8-LTD** (Orlandini et al. 2003 with λ = 1,
    /// the parameter-free form of Orlandini, Moretti & Gavioli 2014). `false` = the D8 of
    /// `compute_flow`. `true`: every land cell of the breached input takes one of the two D8 pointers
    /// of its steepest Tarboton facet, whichever keeps the CUMULATIVE transverse deviation along its
    /// path the smallest. The deviation is inherited from the donor with the largest area. Where no
    /// facet descends strictly (a flat), the pointer is `compute_flow`'s and the deviation is reset,
    /// declared and counted. **The skeleton only**: the hydrology is never routed on it.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub ltd_directions: bool,
    /// **ADR Finding 130-A2 — the confluence clause RESTRICTED to downstream of the junction.** Only
    /// with [`Self::trunk_band`]. A larger line takes a smaller line's cell only when the larger line's
    /// nearest covering sample lies AT OR DOWNSTREAM OF the point where the smaller line's downstream
    /// chain of lines joins the larger line ([`Skeleton::line_parent`]); where that chain never meets
    /// the larger line, there is no junction and no takeover. A cell no sample owns keeps the full
    /// clause's rule. `false` = Finding 128's full clause. DECISION, gated.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub trunk_band_downstream: bool,
    /// **ADR Finding 131-A — the CONCORDANT confluence** (Playfair 1802: a tributary joins the trunk at the
    /// level of the trunk's bed). Only with [`Self::trunk_band`]. The cells Finding 130-A0 located — on a
    /// larger line's floor band, UPSTREAM of the junction of the smaller line (their owner) with it — are
    /// laid at the floor of the larger line AT THE JUNCTION, instead of from its nearest covering sample.
    /// No other cell changes: every other takeover of the clause is off. `false` = off. DECISION, gated.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub trunk_band_concordant: bool,
    /// **ADR Finding 132-P4 — "a present lake is a base level".** The author's decision (2026-09-29): *« Un lac
    /// présent est un niveau de base pour les rivières qui s'y jettent ; un lac vidé ne l'est plus. »* The base
    /// follows the lake's existence, not the breach's state (`breach_monotone` conditions the drainage, it does
    /// not empty a lake). With `Some`, a χ path that enters a lake's footprint stops there: the lake's cells take
    /// χ = 0 and base = the lake's surface (its col), and every cell upstream integrates χ from the lake with that
    /// base. `None` = off (a breached lake passes χ through, Findings 122-131), byte-identical. Which lakes count
    /// as present is [`LakeBase`]'s choice — see there for the circularity it names. DECISION, gated.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lake_base: Option<LakeBase>,
    /// **ADR Finding 140 — the gorge's retreat by age** (`spec_gorge_age_v3.md`). With `Some` (and
    /// `lake_base = InputLakesAndBasins`, `basin_base`), each present lake ≥ 1 km² stands at `L(r_lake)`
    /// between its input spill and its floor, its outlet leaves through a gorge hung from that level (the
    /// invariant on the samples' floor), its rim keeps `≥ L(r_lake)`, and its falls are computed
    /// ([`Skeleton::gorge_falls`]). `None` = off, byte-identical (skipped in the serialised key). DECISION,
    /// gated.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gorge_retreat: Option<GorgeRetreat>,
}

/// ADR Finding 140 — the gorge's retreat parameters (`spec_gorge_age_v3.md` § 1).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct GorgeRetreat {
    /// The world's retreat, from the age selector by map (a): ×0.7 → 0, ×1 → 1, ×1.4 → 2. DECISION.
    pub r_world: f32,
    /// The per-lake exponent: `r_lake = min(2, r_world·(A/A_ref)^p)`, p ∈ {0, 0.25, 0.5}. DECISION (the author:
    /// default 0.5). The link to the knickpoint celerity (F115, `c = K·A^m`) is ANCHORED in form, not in value.
    pub p: f32,
    /// The gorge's steepness factor: `S = min(tan 28°, max(m·S_loi(A), S_req))`. DECISION, m = 10.
    pub m: f32,
    /// The reference outlet area, km². A DECLARED CONSTANT: 418.7 (MEASURED, F139-R: the median outlet area of
    /// the témoin's 14 D8 lakes), never the current world's, so a lake's retreat does not depend on the others.
    pub a_ref_km2: f32,
    /// ADR Finding 142 — spec v4's scope: only the bodies that are input lakes standing alone (their ring touches
    /// no other body's closed-depression cell and no cell under the sea), the construction's proxy for the 14 D8
    /// lakes. The other bodies keep the extended lake base. `false` = v3 (every body ≥ 1 km²).
    #[serde(default)]
    pub d8_scope: bool,
    /// ADR Finding 142 — spec v4's minimal plain ([`gorge_plain`]): a drained bowl's residual hollows filled to
    /// their spill. `false` = v3 (none).
    #[serde(default)]
    pub plain: bool,
    /// ADR Finding 143 — the EXACT scope: the lowest cells of the bodies to keep (`u32::MAX` pads), a per-world list
    /// measured on ON's final lakes (the témoin's 14 D8 lakes, by footprint). It replaces the proxy test. The
    /// construction knows no final lake, so no construction-side criterion gives the exact set (F143-S). `None` = the
    /// proxy (v4). Skipped in the key when `None`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scope_lows: Option<[u32; 16]>,
    /// ADR Finding 143-C1 — BENCH OPTION: φ = 0 (no head fall). Skipped in the key when `false`.
    #[serde(default, skip_serializing_if = "is_false")]
    pub phi_zero: bool,
    /// ADR Finding 143-C2 — BENCH OPTION: the light pass does not modify the designed geometry (the kept bodies'
    /// footprints, their rings, and the gorge's corridor: its path cells and their radius-2 neighbourhood); those
    /// cells are restored after it (`production_upscale.rs`). Skipped in the key when `false`.
    #[serde(default, skip_serializing_if = "is_false")]
    pub freeze_design: bool,
}

fn is_false(b: &bool) -> bool {
    !*b
}

impl GorgeRetreat {
    /// Spec v3's values: m = 10, A_ref = 418.7 km².
    pub fn v3(r_world: f32, p: f32) -> Self {
        Self { r_world, p, m: 10.0, a_ref_km2: 418.7, d8_scope: false, plain: false, scope_lows: None, phi_zero: false, freeze_design: false }
    }

    /// Spec v4: v3 with the D8 scope and the minimal plain (`spec_gorge_age_v4.md`).
    pub fn v4(r_world: f32, p: f32) -> Self {
        Self { d8_scope: true, plain: true, ..Self::v3(r_world, p) }
    }
}

/// ADR Finding 140 — one fall the gorge construction lays, tagged (its drop exceeds `H_f`).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct GorgeFall {
    /// `"waterfall_head"` (at the lip) or `"waterfall_shortage"` (at the next base, beyond the 28° cap).
    pub kind: GorgeFallKind,
    /// The cell (x, y) of the fall's top.
    pub x: u32,
    pub y: u32,
    /// The lake's outflow cell (x, y): the lake cell just upstream of the col.
    pub exit_x: u32,
    pub exit_y: u32,
    /// The drop, m.
    pub height_m: f32,
}

/// ADR Finding 140 — the two causes of a fall (the author's « les deux »).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GorgeFallKind {
    WaterfallHead,
    WaterfallShortage,
}

/// ADR Finding 140 — a present lake as the gorge construction saw it (diagnostics, benches).
#[derive(Clone, Debug, Default)]
pub struct GorgeBody {
    /// The body's lowest cell (its key) and size.
    pub low: u32,
    pub cells: u32,
    /// Its outflow cell and the col (the outflow's receiver), `u32::MAX` if none.
    pub exit: u32,
    pub col: u32,
    /// Levels (m): the input spill, the bed sill (minimax from the lowest cell to the outflow), the floor.
    pub l_in: f32,
    pub l_bed: f32,
    pub l_floor: f32,
    /// The col's drained area (km²), the lake's retreat and its level `L(r_lake)`.
    pub a_out_km2: f32,
    pub r_lake: f32,
    pub level: f32,
    /// The head fall's share, the gorge's slope, the gorge's drop D_g (m), the head fall and the shortage (m).
    pub phi: f32,
    pub slope: f32,
    pub d_g: f32,
    pub head_fall_m: f32,
    pub shortage_m: f32,
    /// The outlet path's length (km) and the number of cells the gorge raised.
    pub path_km: f32,
    pub gorge_cells: u32,
    /// ADR Finding 142 — an input lake (else a closed-depression component); its outer ring touches a
    /// closed-depression cell of another body (a merged water body), or a cell under the sea.
    pub input_lake: bool,
    pub touches_depression: bool,
    pub touches_below_sea: bool,
}

/// ADR Finding 140 — φ, the head fall's share: 0 with probability 1/3, otherwise U[0.1, 0.5], from splitmix64 of
/// `key` (the body's lowest cell: the construction knows no seed and no final lake id; deviation 3 of the spec).
fn gorge_phi(key: u64) -> f32 {
    let mut z = key.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    let u = (z >> 11) as f64 / (1u64 << 53) as f64;
    if u < 1.0 / 3.0 { 0.0 } else { (0.1 + 0.4 * (u - 1.0 / 3.0) / (2.0 / 3.0)) as f32 }
}

/// ADR Finding 132-P4 — which lakes are "present" when the construction integrates χ.
///
/// ⚠️ **A circularity, named.** The lakes of the world's `lakes.json` exist only AFTER the construction and the
/// light pass (the HD assembly reads the delivered field). When χ integrates, the only lakes known are the
/// construction INPUT's pre-drainage lakes. The fixed point (a second pass on a first pass's lakes) needs a lake
/// set from outside this config and is not built.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum LakeBase {
    /// Every lake of the construction input's pre-drainage (`c1_drainage_windowed`'s `lake_map`) is present.
    InputLakes,
    /// **ADR Finding 133-F** — [`Self::InputLakes`], and every `basin_base` closed depression (the open-ocean
    /// flood's spill above the breached input) is ALSO a present lake, filled to its col. The author's decision
    /// (2026-09-30): *« Les lacs des bassins sous la mer sont des lacs présents. »* χ stops at the depression's
    /// SHORE (the first of its cells a path enters) with χ = 0 and base = the col, instead of walking to the col —
    /// which is also where Finding 132's residual was born (that walk pushed cells with no lake test). The level is
    /// the col, known when χ integrates; a basin lake held BELOW its col by an endorheic balance would be known only
    /// after the water balance on the delivered field — circular, named, not built. Needs `basin_base`.
    InputLakesAndBasins,
}

/// ADR Finding 127-B — the interpolant whose gradient the tracé descends.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum TraceInterp {
    /// Catmull-Rom bicubic (Keys, a = −0.5), with its analytic gradient.
    Bicubic,
    /// Bilinear: its gradient is piecewise constant per cell — the negative control (axis artefacts).
    Bilinear,
}

/// ADR Finding 127-B — the tracé's parameters.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SkeletonTrace {
    pub interp: TraceInterp,
    /// The retraced segments: `A < below_km2` (km², domain). DECISION (Finding 126: the 1–10 km²
    /// sub-valleys carry B2's R8).
    pub below_km2: f32,
    /// Below this slope (m/m) a step falls back on the D8 path. **PROXY**.
    pub flat_slope: f32,
}

impl SkeletonTrace {
    /// Finding 127's setting: the 1–10 km² segments, a 1 m/km flat threshold.
    pub fn f127(interp: TraceInterp) -> Self {
        Self { interp, below_km2: 10.0, flat_slope: 1e-3 }
    }
}

/// ADR Finding 127-B — what the tracé did.
#[derive(Clone, Debug, Default)]
pub struct TraceStats {
    /// Segments retraced (≥ 2 own cells below `below_km2`).
    pub segments: usize,
    /// ... that entered their D8 receiver's corridor first (kept).
    pub same_receiver: usize,
    /// ... that entered ANOTHER line's corridor, or the sea, first (their D8 geometry kept).
    pub other_receiver: usize,
    /// ... that ran out of steps (their D8 geometry kept).
    pub lost: usize,
    /// Continuous half-cell steps, and D8 fallback steps on flats (all retraced segments).
    pub steps: usize,
    pub flat_steps: usize,
    /// Every traced sample's distance to its own D8 path, in W of the nearest D8 cell (kept segments).
    pub dev_w: Vec<f32>,
    /// DIAGNOSTIC (Finding 127-B0): every retraced path as the descent drew it, whatever its outcome
    /// (`true` = it reached its own receiver and was kept). Nothing reads it to build.
    pub paths: Vec<(bool, Vec<(f32, f32)>)>,
}

/// ADR Finding 124-B1 — the wall's profile and its texture. Every value is a **PROXY**.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct WallProfile {
    /// Concave FOOT: the slope ramps linearly from 0 at the floor edge to the colluvial slope over
    /// this horizontal length (m). PROXY.
    pub foot_m: f32,
    /// Convex CREST: a smooth minimum of the wall and the terrain, of this height scale (m). PROXY.
    pub crest_m: f32,
    /// Radius (m) of the box blur that defines the field's local TREND; the detail put back on the
    /// wall is `field − trend`. PROXY.
    pub detail_radius_m: f32,
    /// Gain of the detail put back on the wall (1 = the replaced field's own detail). PROXY.
    pub detail_gain: f32,
    /// **ADR Finding 127-A — bruit ↔ pied, the foot clause.** `None` = Findings 124–126: the detail
    /// is put back on the whole wall, foot included (Finding 126-C: under B1 the teeth live in the
    /// concave foot, 89.6 % of their mouths). `Some(τ)`: the detail's amplitude is ZERO on the foot
    /// (`u < foot_m`), FULL beyond `u = (1 + τ)·foot_m`, with a smoothstep between — the transition
    /// is declared as a fraction τ of the foot's own length. With `foot_m` 0 (no foot) it is a no-op.
    /// **PROXY** (τ).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub foot_quiet: Option<f32>,
}

impl WallProfile {
    /// Finding 124's first setting: a 300 m foot, a 20 m crest, the detail below 500 m at gain 1.
    pub fn f124() -> Self {
        Self { foot_m: 300.0, crest_m: 20.0, detail_radius_m: 500.0, detail_gain: 1.0, foot_quiet: None }
    }
}

impl ValleyConstruction {
    /// **THE DEFINITION** (ADR Finding 124, the author's decision D: "uniquement leur bassin"):
    /// Finding 121's laws with Finding 122-B's `basin_base` — χ from the col of every closed
    /// depression. What the viz's "Vallées construites" states build.
    pub fn new(age_k: f32, light_k_time_fraction: Option<f32>) -> Self {
        Self { basin_base: true, ..Self::f121(age_k, light_k_time_fraction) }
    }

    /// Finding 121's HISTORICAL values (χ from the sea everywhere), kept so Findings 121–123 can be
    /// reproduced and so the viz can show the old definition beside the new one.
    pub fn f121(age_k: f32, light_k_time_fraction: Option<f32>) -> Self {
        Self {
            a_min_km2: 10.0,
            age_k,
            a_c_km2: 0.1,
            // W(10 km²) = 200 m ⇒ a = 200 / 10^0.3
            width_coef_m: 200.0 / 10f32.powf(0.3),
            width_exp: 0.3,
            wall_deg: 28.0,
            base_m: 0.5,
            smooth_m: 250.0,
            light_k_time_fraction,
            basin_base: false,
            wall_profile: None,
            width_age_gamma: None,
            wall_sea_floor_m: None,
            smooth_w: None,
            skeleton_trace: None,
            trunk_band: false,
            ltd_directions: false,
            trunk_band_downstream: false,
            trunk_band_concordant: false,
            lake_base: None,
            gorge_retreat: None,
        }
    }

    /// Finding 122-B's name for [`Self::new`], kept for the benches that used it.
    pub fn f122(age_k: f32, light_k_time_fraction: Option<f32>) -> Self {
        Self::new(age_k, light_k_time_fraction)
    }

    /// Valley floor width (m) at drained area `a_km2`: `a·A^b`, times `(k/k₀)^γ` under Finding
    /// 124-5's W(k).
    pub fn width_m(&self, a_km2: f32) -> f32 {
        let age = self.width_age_gamma.map_or(1.0, |g| (self.age_k / F121_AGE_K).max(0.0).powf(g));
        self.width_coef_m * age * a_km2.max(0.0).powf(self.width_exp)
    }
}

/// A config serialised before Finding 124 carries no `basin_base`: it reads the definition.
fn basin_base_default() -> bool {
    true
}

/// The skeleton a construction stands on: the D8 network of the BREACHED input field (Finding
/// 120's A0 instrument, reproduced exactly), χ along it, and the smoothed trunk polylines.
pub struct Skeleton {
    pub width: usize,
    pub height: usize,
    pub cell_m: f32,
    /// D8 drained area, km², per cell.
    pub area_km2: Vec<f32>,
    /// χ (m) of each land cell along its D8 path to its base; NaN on sea cells.
    pub chi_m: Vec<f32>,
    /// Altitude (m) of each land cell's base: `base_m` above the sea, or the altitude of a land
    /// terminal (a D8 cell with no receiver).
    pub base_alt_m: Vec<f32>,
    /// Trunk cells: land, `A ≥ a_min_km2`.
    pub trunk: Vec<bool>,
    /// Smoothed, densified trunk polylines: (x, y) in continuous cell units (cell centres at
    /// +0.5, UNWRAPPED along each line), floor altitude (m), half floor width (m).
    pub polylines: Vec<Vec<(f32, f32, f32, f32)>>,
    /// ADR Finding 127-B — what the off-grid tracé did (`None` when it is off).
    pub trace_stats: Option<TraceStats>,
    /// The D8 pointer of every cell the skeleton was built on (`compute_flow`'s, or D8-LTD's under
    /// Finding 128-C); `DIR_NONE` at a base.
    pub direction: Vec<u8>,
    /// ADR Finding 128-C — land cells where D8-LTD found no strictly descending facet and fell back on
    /// `compute_flow`'s pointer (0 when LTD is off).
    pub ltd_flat_cells: usize,
    /// ADR Finding 130-A — per polyline, the line it ends on and the index, in THAT line's polyline, of
    /// the junction (the parent's sample laid from the cell the line ends on). A line traced from a
    /// head stops on the first cell an earlier line already owns, so its last cell IS the junction.
    /// `None`: the line ends at the sea or a base, or the skeleton is traced off-grid (Finding 127-B).
    pub line_parent: Vec<Option<(u32, u32)>>,
    /// ADR Finding 140 — the rim clamp: every outer-ring cell of a present lake keeps `≥` this (m); NaN elsewhere.
    /// `None` when the gorge is off.
    pub rim_floor_m: Option<Vec<f32>>,
    /// ADR Finding 140 — the tagged falls the gorge construction laid (empty when off).
    pub gorge_falls: Vec<GorgeFall>,
    /// ADR Finding 140 — the present lakes it treated, and each cell's body (`u32::MAX` none), for the benches.
    pub gorge_bodies: Vec<GorgeBody>,
    pub gorge_body_of: Option<Vec<u32>>,
    /// ADR Finding 143 — the gorge's designed path: each kept body's outlet cells from the col to where the gorge
    /// meets the law below (empty when off).
    pub gorge_path: Vec<u32>,
}

impl Skeleton {
    /// Floor altitude (m) of a land cell at age `k`.
    pub fn floor_m(&self, k: usize, age_k: f32) -> f32 {
        self.base_alt_m[k] + age_k * self.chi_m[k]
    }
}

fn recv(dir: &[u8], k: usize, w: usize, h: usize) -> Option<usize> {
    let d = dir[k];
    if d == DIR_NONE {
        return None;
    }
    let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
    let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
    Some(ny * w + nx)
}

/// Build the skeleton of `field` (normalised altitude). `domain_km` is the span of the grid.
///
/// The flow is Finding 120's: `c1_drainage_windowed` → `breach_monotone` → the drainage chain's
/// own `compute_flow` on the breached field (same sea level, same flat perturbation).
pub fn skeleton(
    field: &GridF32,
    vc: &ValleyConstruction,
    ss: &SteinSteinParams,
    domain_km: f32,
) -> Skeleton {
    skeleton_patched(field, vc, ss, domain_km, None)
}

/// ADR Finding 131-P — the pointers the skeleton will stand on, handed to a bench before the areas are
/// accumulated: `(breached field, the pre-drainage's lake map, &mut pointers)`. The areas are then recounted
/// on the patched pointers. Diagnostic only; `None` is [`skeleton`], byte-identical.
pub type SkeletonPatch<'a> = &'a dyn Fn(&GridF32, &[u32], &mut Vec<u8>);

/// [`skeleton`] with an optional [`SkeletonPatch`].
pub fn skeleton_patched(
    field: &GridF32,
    vc: &ValleyConstruction,
    ss: &SteinSteinParams,
    domain_km: f32,
    patch: Option<SkeletonPatch>,
) -> Skeleton {
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let cell_m = domain_km * 1000.0 / w as f32;
    let cell_km2 = (cell_m / 1000.0) * (cell_m / 1000.0);
    let dcfg = C1DrainageConfig::default();
    let (bf, lake_map, lake_level) = {
        let pre = c1_drainage_windowed(field, None, &dcfg, ss, domain_km);
        let bf = breach_monotone(field, &pre.flow.filled, &pre.lake_map, C1_SEA_LEVEL_NORM, w, h);
        // ADR Finding 132-P4 -- the input's lakes and their surface (m), kept only when asked for
        let level: Option<Vec<f32>> = vc.lake_base.map(|_| {
            (0..n)
                .map(|k| {
                    if pre.lake_map[k] != 0 {
                        c1_altitude_norm_to_metres(pre.flow.filled.data[k], ss)
                    } else {
                        f32::NAN
                    }
                })
                .collect()
        });
        let keep = patch.is_some() || vc.lake_base.is_some();
        (bf, keep.then_some(pre.lake_map), level)
    };
    let flow = compute_flow(
        &bf,
        &FlowConfig {
            sea_level: C1_SEA_LEVEL_NORM,
            flat_perturbation: dcfg.flat_perturbation.clone(),
            dinf: dcfg.dinf,
        },
    );
    let land: Vec<bool> = bf.data.iter().map(|&v| v > C1_SEA_LEVEL_NORM).collect();
    // ADR Finding 128-C -- the skeleton's own directions (D8, or D8-LTD on the same breached surface)
    let (direction, accumulation, ltd_flat_cells) = if vc.ltd_directions {
        let zm: Vec<f32> = bf.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
        let sink: Vec<bool> = land.iter().map(|&l| !l).collect();
        let (d, flats) = ltd_directions(&zm, &sink, &flow.direction, w, h, cell_m);
        let acc = accumulate_cells(&d, &sink, w, h);
        (d, acc, flats)
    } else {
        (flow.direction, flow.accumulation.data, 0)
    };
    // ADR Finding 131-P -- a bench's patch of the pointers, the areas recounted on them
    let (direction, accumulation) = match patch {
        Some(f) => {
            let mut d = direction;
            f(&bf, lake_map.as_deref().expect("kept for the patch"), &mut d);
            let sink: Vec<bool> = land.iter().map(|&l| !l).collect();
            let acc = accumulate_cells(&d, &sink, w, h);
            (d, acc)
        }
        None => (direction, accumulation),
    };
    let area_km2: Vec<f32> = accumulation.iter().map(|&a| a * cell_km2).collect();
    drop(accumulation);
    let dir = &direction;
    // ADR Finding 122-B -- the closed depressions and their spill levels: a priority flood seeded
    // on the OPEN OCEAN only, so an inland below-sea basin and a held lake both rise to their col.
    let spill: Option<Vec<f32>> = vc.basin_base.then(|| ocean_flood(&bf));
    let eps = 0.01 / c1_altitude_norm_to_metres(1.0, ss).max(1.0); // 1 cm, in norm units
    let dep = |k: usize| -> bool { spill.as_ref().is_some_and(|f| f[k] > bf.data[k] + eps) };
    let spill_m = |k: usize| -> f32 {
        c1_altitude_norm_to_metres(spill.as_ref().expect("basin_base")[k], ss)
    };

    // ADR Finding 140 -- the gorge's retreat: the present lakes as bodies and their level L(r_lake), known before the
    // χ walk so that the lake stops take it as their base
    let gorge = match (vc.gorge_retreat, vc.lake_base, lake_map.as_deref(), lake_level.as_deref()) {
        (Some(g), Some(LakeBase::InputLakesAndBasins), Some(lm), Some(lv)) if vc.basin_base => Some(gorge_bodies(
            field, lm, lv, &land, dir, &area_km2, &dep, &spill_m, w, h, cell_km2, ss, vc, &g,
        )),
        _ => None,
    };
    let gorge_level = |c: usize, fallback: f32| -> f32 {
        match &gorge {
            Some((bo, bs)) if bo[c] != u32::MAX => bs[bo[c] as usize].level,
            _ => fallback,
        }
    };

    // χ, integrated from each base upstream (memoised walk down the receivers)
    let mut chi = vec![f32::NAN; n];
    let mut base = vec![f32::NAN; n];
    let link_m = |a: usize, b: usize| -> f32 {
        let diag = (a % w != b % w) && (a / w != b / w);
        if diag { cell_m * std::f32::consts::SQRT_2 } else { cell_m }
    };
    let mut path: Vec<usize> = Vec::new();
    for s in 0..n {
        if !land[s] || !chi[s].is_nan() {
            continue;
        }
        path.clear();
        let mut c = s;
        let (b, mut x);
        loop {
            if !chi[c].is_nan() {
                b = base[c];
                x = chi[c];
                break;
            }
            // ADR Finding 132-P4 -- a present lake is a base level: χ stops at its footprint
            if let (Some(lm), Some(lv)) = (lake_map.as_deref(), lake_level.as_deref())
                && lm[c] != 0
            {
                b = gorge_level(c, lv[c]); // ADR Finding 140 -- the lake at L(r_lake)
                x = 0.0;
                chi[c] = 0.0;
                base[c] = b;
                break;
            }
            // ADR Finding 122-B -- inside a closed depression: walk to its col; the base is the
            // spill level and χ is counted from the col.
            if land[c] && dep(c) {
                let lvl = spill_m(c);
                // ADR Finding 133-F -- a basin lake is a present lake: χ stops at its shore
                if vc.lake_base == Some(LakeBase::InputLakesAndBasins) {
                    b = gorge_level(c, lvl); // ADR Finding 140 -- the lake at L(r_lake)
                    x = 0.0;
                    chi[c] = 0.0;
                    base[c] = b;
                    break;
                }
                let mut terminal = false;
                loop {
                    path.push(c);
                    match recv(dir, c, w, h) {
                        None => {
                            terminal = true;
                            break;
                        }
                        Some(r) if !land[r] || !dep(r) => break, // the col (or the sea)
                        Some(r) => c = r,
                    }
                    if path.len() > n {
                        panic!("D8 cycle in a closed depression of the valley skeleton");
                    }
                }
                if terminal {
                    let t = path.pop().expect("pushed");
                    chi[t] = 0.0;
                    base[t] = lvl;
                }
                b = lvl;
                x = 0.0;
                break;
            }
            match recv(dir, c, w, h) {
                None => {
                    // a land terminal: its own altitude is the base
                    b = c1_altitude_norm_to_metres(bf.data[c], ss);
                    x = 0.0;
                    chi[c] = 0.0;
                    base[c] = b;
                    break;
                }
                Some(r) => {
                    path.push(c);
                    if !land[r] {
                        // ADR Finding 122-B -- an INLAND below-sea basin is not the sea: its base
                        // is its spill level. The open ocean keeps Finding 83's `base_m`.
                        b = if dep(r) { spill_m(r) } else { vc.base_m };
                        x = 0.0;
                        // the link c → r is integrated in the unwind below, with r = sea
                        break;
                    }
                    if path.len() > n {
                        panic!("D8 cycle in the valley skeleton");
                    }
                    c = r;
                }
            }
        }
        for i in (0..path.len()).rev() {
            let c = path[i];
            let r = recv(dir, c, w, h).expect("a path cell has a receiver");
            let a = area_km2[c].max(vc.a_c_km2);
            x += (1.0 / a).sqrt() * link_m(c, r);
            chi[c] = x;
            base[c] = b;
        }
    }

    // ADR Finding 140 -- the gorge: the invariant on the outlet path's floor, the rim clamp, the falls
    let mut gorge_path: Vec<u32> = Vec::new();
    let (rim_floor_m, gorge_falls, gorge_bodies_out, gorge_body_of) = match gorge {
        Some((body_of, mut bodies)) => {
            let g = vc.gorge_retreat.expect("the gate");
            let tan28 = 28f32.to_radians().tan();
            let nb8 = |c: usize, k: usize| -> usize {
                (((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize) * w
                    + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize
            };
            // the rim: every outer-ring cell (not in another body) keeps >= the lake's level
            let mut rf = vec![f32::NAN; n];
            for c in 0..n {
                let bi = body_of[c];
                if bi == u32::MAX {
                    continue;
                }
                let lvl = bodies[bi as usize].level;
                for k in 0..8 {
                    let m = nb8(c, k);
                    if body_of[m] == u32::MAX && land[m] {
                        rf[m] = if rf[m].is_nan() { lvl } else { rf[m].max(lvl) };
                    }
                }
            }
            let mut falls = Vec::new();
            for b in bodies.iter_mut() {
                if b.col == u32::MAX {
                    continue;
                }
                let lv = b.level;
                let mut path = vec![b.col as usize];
                let mut s = vec![0f32];
                let mut c = b.col as usize;
                while path.len() < 20_000 {
                    let Some(r) = recv(dir, c, w, h) else { break };
                    if !land[r] || body_of[r] != u32::MAX {
                        break;
                    }
                    s.push(s.last().expect("one") + link_m(c, r));
                    path.push(r);
                    c = r;
                }
                // the law below along the path, frozen before the gorge raises any base
                let lawv: Vec<f32> = path.iter().map(|&c| base[c] + vc.age_k * chi[c]).collect();
                let s_req = |top: f32| -> f32 {
                    if top <= lawv[0] {
                        return 0.0;
                    }
                    (1..path.len()).map(|i| (top - lawv[i]) / s[i]).fold(f32::INFINITY, f32::min).max(0.0)
                };
                let sl = vc.age_k * b.a_out_km2.max(vc.a_c_km2).powf(-0.5);
                let slope_for = |top: f32| -> f32 {
                    let r = s_req(top);
                    (g.m * sl).max(if r.is_finite() { r } else { 0.0 }).min(tan28)
                };
                let fit = |top: f32, sg: f32| (0..path.len()).find(|&i| top - sg * s[i] <= lawv[i]);
                let s_end = *s.last().expect("one");
                let s0 = slope_for(lv);
                let dg = match fit(lv, s0) {
                    Some(i) => lv - lawv[i],
                    None => s0 * s_end,
                };
                let top = lv - b.phi * dg;
                let sg = slope_for(top);
                let mut raised = 0u32;
                let mut end = path.len();
                for (i, &c) in path.iter().enumerate() {
                    let z = if i == 0 { lv } else { top - sg * s[i] };
                    if i > 0 && z <= lawv[i] {
                        end = i;
                        break;
                    }
                    if chi[c].is_finite() && z - vc.age_k * chi[c] > base[c] {
                        base[c] = z - vc.age_k * chi[c];
                        raised += 1;
                    }
                }
                let hf = (1.5 * sg * 100.0).max(10.0);
                let head = if path.len() > 1 { lv - (top - sg * s[1]) } else { 0.0 };
                let short = if end == path.len() && path.len() > 1 {
                    ((top - sg * s_end) - lawv[path.len() - 1]).max(0.0)
                } else {
                    0.0
                };
                let (ex, cx) = (b.exit as usize, b.col as usize);
                if head > hf {
                    falls.push(GorgeFall {
                        kind: GorgeFallKind::WaterfallHead,
                        x: (cx % w) as u32,
                        y: (cx / w) as u32,
                        exit_x: (ex % w) as u32,
                        exit_y: (ex / w) as u32,
                        height_m: head,
                    });
                }
                if short > hf {
                    let e = path[path.len() - 1];
                    falls.push(GorgeFall {
                        kind: GorgeFallKind::WaterfallShortage,
                        x: (e % w) as u32,
                        y: (e / w) as u32,
                        exit_x: (ex % w) as u32,
                        exit_y: (ex / w) as u32,
                        height_m: short,
                    });
                }
                gorge_path.extend(path[..end].iter().map(|&c| c as u32));
                (b.slope, b.d_g, b.head_fall_m, b.shortage_m, b.path_km, b.gorge_cells) =
                    (sg, dg, head, short, s_end / 1000.0, raised);
            }
            (Some(rf), falls, bodies, Some(body_of))
        }
        None => (None, Vec::new(), Vec::new(), None),
    };

    // trunks and their polylines
    let trunk: Vec<bool> = (0..n).map(|k| land[k] && area_km2[k] >= vc.a_min_km2).collect();
    let mut donors = vec![0u8; n];
    for k in 0..n {
        if trunk[k]
            && let Some(r) = recv(dir, k, w, h)
            && trunk[r]
        {
            donors[r] = donors[r].saturating_add(1);
        }
    }
    let mut visited = vec![false; n];
    let mut raw_lines: Vec<Vec<usize>> = Vec::new();
    // ADR Finding 130-A -- the line (index in `raw_lines`) and position that own each visited cell
    let mut owner: Vec<(u32, u32)> = vec![(u32::MAX, 0); n];
    for head in 0..n {
        if !trunk[head] || donors[head] != 0 || visited[head] {
            continue;
        }
        let mut line = Vec::new();
        let mut c = head;
        loop {
            line.push(c);
            if visited[c] {
                break; // a junction: the line ends ON the parent's cell
            }
            visited[c] = true;
            owner[c] = (raw_lines.len() as u32, line.len() as u32 - 1);
            match recv(dir, c, w, h) {
                Some(r) if trunk[r] => c = r,
                _ => break,
            }
        }
        if line.len() >= 2 {
            raw_lines.push(line);
        }
    }
    let (polylines, trace_stats, line_parent) = match vc.skeleton_trace {
        None => {
            let (polylines, dense_at): (Vec<_>, Vec<_>) = raw_lines
                .iter()
                .map(|line| {
                    let pts = line_samples(line, w, h, &base, &chi, &area_km2, vc);
                    let sm = smooth_positions(pts, line, &area_km2, vc, cell_m);
                    let at = densify_index(&sm);
                    (densify(sm), at)
                })
                .unzip();
            // ADR Finding 130-A -- each line's parent and the junction's index in the parent's polyline
            let line_parent = raw_lines
                .iter()
                .enumerate()
                .map(|(li, line)| {
                    let (pl, pos) = owner[*line.last().expect("len >= 2")];
                    (pl != u32::MAX && pl as usize != li)
                        .then(|| (pl, dense_at[pl as usize][pos as usize] as u32))
                })
                .collect();
            (polylines, None, line_parent)
        }
        Some(tr) => {
            let (p, st) =
                traced_polylines(&raw_lines, field, ss, cell_m, &base, &chi, &area_km2, vc, tr);
            let np = p.len();
            (p, Some(st), vec![None; np])
        }
    };
    Skeleton {
        width: w,
        height: h,
        cell_m,
        area_km2,
        chi_m: chi,
        base_alt_m: base,
        trunk,
        polylines,
        trace_stats,
        direction,
        ltd_flat_cells,
        line_parent,
        rim_floor_m,
        gorge_falls,
        gorge_bodies: gorge_bodies_out,
        gorge_body_of,
        gorge_path,
    }
}

/// ADR Finding 140 — the present lakes ≥ 1 km² as bodies (the input lakes and the closed-depression components),
/// each with its levels, outflow, col, outlet area, retreat and level `L(r_lake)`. Returns each cell's body
/// (`u32::MAX` none) and the bodies.
#[allow(clippy::too_many_arguments)]
fn gorge_bodies(
    field: &GridF32,
    lm: &[u32],
    lv: &[f32],
    land: &[bool],
    dir: &[u8],
    area_km2: &[f32],
    dep: &dyn Fn(usize) -> bool,
    spill_m: &dyn Fn(usize) -> f32,
    w: usize,
    h: usize,
    cell_km2: f32,
    ss: &SteinSteinParams,
    vc: &ValleyConstruction,
    g: &GorgeRetreat,
) -> (Vec<u32>, Vec<GorgeBody>) {
    use std::collections::HashMap;
    let n = w * h;
    let nb8 = |c: usize, k: usize| -> usize {
        (((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize) * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize
    };
    let z = |c: usize| c1_altitude_norm_to_metres(field.data[c], ss);
    // label: the input lakes by id, then the depression components (8-connected)
    let mut label = vec![u32::MAX; n];
    let mut members: Vec<(Vec<usize>, bool)> = Vec::new(); // (cells, is an input lake)
    let mut by_id: HashMap<u32, u32> = HashMap::new();
    for c in 0..n {
        if lm[c] != 0 {
            let id = *by_id.entry(lm[c]).or_insert_with(|| {
                members.push((Vec::new(), true));
                (members.len() - 1) as u32
            });
            label[c] = id;
            members[id as usize].0.push(c);
        }
    }
    for s0 in 0..n {
        if label[s0] != u32::MAX || !land[s0] || !dep(s0) {
            continue;
        }
        let id = members.len() as u32;
        members.push((vec![s0], false));
        label[s0] = id;
        let mut i = 0;
        while i < members[id as usize].0.len() {
            let c = members[id as usize].0[i];
            i += 1;
            for k in 0..8 {
                let m = nb8(c, k);
                if label[m] == u32::MAX && land[m] && lm[m] == 0 && dep(m) {
                    label[m] = id;
                    members[id as usize].0.push(m);
                }
            }
        }
    }
    let mut body_of = vec![u32::MAX; n];
    let mut bodies = Vec::new();
    for (id, (cells, is_lake)) in members.iter().enumerate() {
        if (cells.len() as f32) * cell_km2 < 1.0 {
            continue;
        }
        // ADR Finding 142 -- what the ring touches: another body's depression cell, or the sea
        let (mut t_dep, mut t_sea) = (false, false);
        for &c in cells {
            for k in 0..8 {
                let m = nb8(c, k);
                if label[m] == id as u32 {
                    continue;
                }
                t_sea |= !land[m];
                t_dep |= land[m] && lm[m] == 0 && dep(m);
            }
        }
        let low = *cells.iter().min_by(|&&a, &&b| z(a).total_cmp(&z(b))).expect("non-empty");
        match g.scope_lows {
            // ADR Finding 143 -- the exact scope: the listed bodies only
            Some(list) => {
                if !list.contains(&(low as u32)) {
                    continue;
                }
            }
            // ADR Finding 142 -- spec v4's scope: an input lake standing alone (the D8 lakes' proxy)
            None => {
                if g.d8_scope && !(*is_lake && !t_dep && !t_sea) {
                    continue;
                }
            }
        }
        let bi = bodies.len() as u32;
        for &c in cells {
            body_of[c] = bi;
        }
        let l_in = if *is_lake { lv[cells[0]] } else { spill_m(low) };
        let l_floor = z(low);
        // the outflow: the largest-area cell whose receiver leaves the body
        let mut exit: Option<usize> = None;
        for &c in cells {
            if let Some(r) = recv(dir, c, w, h)
                && label[r] != id as u32
                && exit.is_none_or(|e| area_km2[c] > area_km2[e])
            {
                exit = Some(c);
            }
        }
        let col = exit.and_then(|e| recv(dir, e, w, h));
        // the bed sill: minimax from the lowest cell to the outflow, within the body
        let l_bed = match exit {
            Some(e) => {
                let mut best: HashMap<usize, f32> = HashMap::new();
                let mut heap = BinaryHeap::new();
                best.insert(low, z(low));
                heap.push(Item(z(low), low));
                let mut found = l_in;
                while let Some(Item(d, c)) = heap.pop() {
                    if best.get(&c).is_some_and(|&b| d > b) {
                        continue;
                    }
                    if c == e {
                        found = d;
                        break;
                    }
                    for k in 0..8 {
                        let m = nb8(c, k);
                        if label[m] != id as u32 {
                            continue;
                        }
                        let dm = d.max(z(m));
                        if best.get(&m).is_none_or(|&b| dm < b) {
                            best.insert(m, dm);
                            heap.push(Item(dm, m));
                        }
                    }
                }
                found.min(l_in)
            }
            None => l_in,
        };
        let l_floor = l_floor.min(l_bed);
        let a_out = col.map_or(0.0, |c| area_km2[c]);
        let r_lake = if a_out > 0.0 { (g.r_world * (a_out / g.a_ref_km2).powf(g.p)).min(2.0) } else { 0.0 };
        let level = if r_lake <= 1.0 { l_in - r_lake * (l_in - l_bed) } else { l_bed - (r_lake - 1.0) * (l_bed - l_floor) };
        bodies.push(GorgeBody {
            low: low as u32,
            cells: cells.len() as u32,
            exit: exit.map_or(u32::MAX, |e| e as u32),
            col: col.map_or(u32::MAX, |c| c as u32),
            l_in,
            l_bed,
            l_floor,
            a_out_km2: a_out,
            r_lake,
            level,
            phi: if g.phi_zero { 0.0 } else { gorge_phi(low as u64) },
            input_lake: *is_lake,
            touches_depression: t_dep,
            touches_below_sea: t_sea,
            ..Default::default()
        });
    }
    let _ = vc;
    (body_of, bodies)
}

/// ADR Finding 142 — one body's minimal plain ([`gorge_plain`]).
#[derive(Clone, Debug, Default)]
pub struct GorgePlain {
    /// The body (index in [`Skeleton::gorge_bodies`]) and its retreat.
    pub body: u32,
    pub r_lake: f32,
    /// The footprint's spill (m): its lowest ring cell, the surface of the lake that is left.
    pub spill_m: f32,
    /// Cells raised by more than 1 cm, their area (km²), the deepest fill (m) and the deposited volume (km³).
    pub cells: u32,
    pub area_km2: f32,
    pub max_depth_m: f32,
    pub volume_km3: f64,
}

/// ADR Finding 142 — **the minimal plain** (`spec_gorge_age_v4.md` § 2): in every body the gorge has partly or wholly
/// drained (`r_lake > min_r`; spec v4: 1), the input footprint's cells that are not under the lake left are filled to
/// their spill. A priority-flood fill restricted to the footprint, seeded on its ring (every outside neighbour at its
/// height), so the outlet is the lake left or the gorge's head. A cell whose fill level is the footprint's spill (the
/// lowest ring cell) is the lake's bed, under water: it is left as it is. A cell whose fill level stands higher is a
/// residual hollow: it is raised to that level.
///
/// ⚠️ **The construction's FIRST term that RAISES the terrain**: a deposit (the sediments of the lake's life), not a
/// cut. [`carve`] still only lowers; this runs after it, on its output, before the light pass
/// (`production_upscale.rs`). Its volume is published per body (rule 14). Works in normalised units (the metres map is
/// affine), so the fill copies its spill's exact value and leaves no round-trip micro-pit.
pub fn gorge_plain(field: &mut GridF32, sk: &Skeleton, ss: &SteinSteinParams, min_r: f32) -> Vec<GorgePlain> {
    use std::collections::HashSet;
    let Some(bo) = sk.gorge_body_of.as_deref() else {
        return Vec::new();
    };
    let (w, h) = (sk.width, sk.height);
    assert_eq!((field.width, field.height), (w, h), "skeleton and field sizes differ");
    let n = w * h;
    let cell_km2 = (sk.cell_m / 1000.0) * (sk.cell_m / 1000.0);
    let m_per_norm = c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);
    let eps = 0.01 / m_per_norm; // 1 cm
    let nb8 = |c: usize, k: usize| -> usize {
        (((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize) * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize
    };
    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); sk.gorge_bodies.len()];
    for c in 0..n {
        if bo[c] != u32::MAX {
            cells[bo[c] as usize].push(c);
        }
    }
    let mut out = Vec::new();
    for (bi, b) in sk.gorge_bodies.iter().enumerate() {
        if b.r_lake <= min_r || cells[bi].is_empty() {
            continue;
        }
        let me = bi as u32;
        let mut heap = BinaryHeap::new();
        let mut ring: HashSet<usize> = HashSet::new();
        let mut spill = f32::INFINITY;
        for &c in &cells[bi] {
            for k in 0..8 {
                let m = nb8(c, k);
                if bo[m] != me && ring.insert(m) {
                    spill = spill.min(field.data[m]);
                    heap.push(Item(field.data[m], m));
                }
            }
        }
        let mut done: HashSet<usize> = HashSet::with_capacity(cells[bi].len());
        let mut rec = GorgePlain { body: me, r_lake: b.r_lake, spill_m: c1_altitude_norm_to_metres(spill, ss), ..Default::default() };
        while let Some(Item(lev, c)) = heap.pop() {
            for k in 0..8 {
                let m = nb8(c, k);
                if bo[m] != me || !done.insert(m) {
                    continue;
                }
                let z = field.data[m];
                let l = lev.max(z);
                if l > z && l > spill + eps {
                    let d = (l - z) * m_per_norm;
                    field.data[m] = l;
                    if d > 0.01 {
                        rec.cells += 1;
                    }
                    rec.max_depth_m = rec.max_depth_m.max(d);
                    rec.volume_km3 += d as f64 * cell_km2 as f64 * 1e-3;
                }
                heap.push(Item(l, m));
            }
        }
        rec.area_km2 = rec.cells as f32 * cell_km2;
        out.push(rec);
    }
    out
}

/// ADR Finding 143-C2 — the gorge's designed geometry, for the bench option `freeze_design`: the kept bodies'
/// footprints, their rings (8-neighbours), and the gorge's corridor (the path cells and their Chebyshev radius-2
/// neighbourhood). Empty grid when the gate is off.
pub fn gorge_design_mask(sk: &Skeleton) -> Vec<bool> {
    let (w, h) = (sk.width, sk.height);
    let n = w * h;
    let mut m = vec![false; n];
    let Some(bo) = sk.gorge_body_of.as_deref() else {
        return m;
    };
    let at = |c: usize, dx: i32, dy: i32| -> usize {
        (((c / w) as i32 + dy).rem_euclid(h as i32) as usize) * w + ((c % w) as i32 + dx).rem_euclid(w as i32) as usize
    };
    for c in 0..n {
        if bo[c] != u32::MAX {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    m[at(c, dx, dy)] = true;
                }
            }
        }
    }
    for &c in &sk.gorge_path {
        for dy in -2..=2 {
            for dx in -2..=2 {
                m[at(c as usize, dx, dy)] = true;
            }
        }
    }
    m
}

/// ADR Finding 128-C — D8-LTD pointers (Orlandini et al. 2003 §2, with λ = 1; Orlandini, Moretti &
/// Gavioli 2014 eqs. (1)–(6)) on the surface `zm` (metres, torus), `sink` cells taking no pointer.
///
/// Cells are processed in descending elevation ([O03] §2.2). For each one:
/// - its theoretical direction is Tarboton's steepest facet (the continuous angle `r` from the facet's
///   cardinal toward its diagonal pointer, clamped to the facet);
/// - each candidate step `q` (the cardinal, or the diagonal) has the signed transverse deviation
///   `q × t` from that direction;
/// - the pointer kept is the one whose CUMULATIVE deviation (the inherited one + its own) is the
///   smaller in absolute value.
///
/// The inherited deviation is the one conveyed by the donor with the largest area. A candidate that
/// is not strictly lower is refused. With no strictly descending facet, the cell takes `fallback`'s
/// pointer and conveys a zero deviation; those cells are counted (the second return). Ties go to the
/// cardinal, as [O14] eq. (5) does.
pub fn ltd_directions(
    zm: &[f32],
    sink: &[bool],
    fallback: &[u8],
    w: usize,
    h: usize,
    cell_m: f32,
) -> (Vec<u8>, usize) {
    let (dir, st) = ltd_directions_stats(zm, sink, fallback, w, h, cell_m);
    (dir, st.flats)
}

/// ADR Finding 129-B2 — what D8-LTD's choices rested on.
#[derive(Clone, Debug, Default)]
pub struct LtdStats {
    /// Land cells with no strictly descending facet (the fallback pointer).
    pub flats: usize,
    /// Choices where the two cumulative deviations tie EXACTLY in f32 (`|c1| == |c2|`).
    pub exact_ties: usize,
    /// Choices the tie rule decided: `||c1| − |c2|| ≤ TIE` (exact ties included).
    pub ties: usize,
    /// The cells of `ties`, for locating them.
    pub tie_cells: Vec<u32>,
}

/// [`ltd_directions`], with what its choices rested on (Finding 129-B2).
pub fn ltd_directions_stats(
    zm: &[f32],
    sink: &[bool],
    fallback: &[u8],
    w: usize,
    h: usize,
    cell_m: f32,
) -> (Vec<u8>, LtdStats) {
    ltd_directions_masked(zm, sink, None, fallback, w, h, cell_m)
}

/// ADR Finding 130 — [`ltd_directions_stats`] on any elevation type (`f32` gives the same pointers
/// bit for bit: differences are taken exactly in `f64`, then rounded to `f32` as before), with an
/// optional `active` mask: a cell outside it is never processed and never a receiver (a CLOSED
/// boundary, Orlandini et al. 2003's synthetic valley; B), and `f64` elevations carry a flat
/// resolution too fine for `f32` (P).
pub fn ltd_directions_masked<T: Copy + Into<f64>>(
    z: &[T],
    sink: &[bool],
    active: Option<&[bool]>,
    fallback: &[u8],
    w: usize,
    h: usize,
    cell_m: f32,
) -> (Vec<u8>, LtdStats) {
    use std::f32::consts::{FRAC_PI_4, SQRT_2};
    let zm = |k: usize| -> f64 { z[k].into() };
    let on = |k: usize| active.is_none_or(|a| a[k]);
    let n = w * h;
    let nb = |c: usize, k: usize| -> usize {
        let x = ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
        let y = ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize;
        y * w + x
    };
    // the eight facets: (cardinal pointer, diagonal pointer) in D8 order (even = cardinal)
    const FACETS: [(usize, usize); 8] = [(0, 1), (0, 7), (2, 1), (2, 3), (4, 3), (4, 5), (6, 5), (6, 7)];
    let mut order: Vec<u32> = (0..n as u32).filter(|&c| !sink[c as usize] && on(c as usize)).collect();
    order.sort_by(|&a, &b| zm(b as usize).total_cmp(&zm(a as usize)).then(a.cmp(&b)));
    let mut dir = vec![DIR_NONE; n];
    let mut area = vec![1f32; n];
    let mut best_area = vec![0f32; n];
    let mut best_dev = vec![0f32; n];
    let mut flats = 0usize;
    let (mut exact_ties, mut ties, mut tie_cells) = (0usize, 0usize, Vec::new());
    for &c in &order {
        let c = c as usize;
        let e0 = zm(c);
        // Tarboton's steepest facet
        let mut best: Option<(f32, usize, usize, f32)> = None; // (slope, cardinal, diagonal, r)
        for &(kc, kd) in &FACETS {
            let (e1, e2) = (zm(nb(c, kc)), zm(nb(c, kd)));
            let (s1, s2) = (((e0 - e1) as f32) / cell_m, ((e1 - e2) as f32) / cell_m);
            let r = s2.atan2(s1);
            let (r, s) = if r < 0.0 {
                (0.0, s1)
            } else if r > FRAC_PI_4 {
                (FRAC_PI_4, ((e0 - e2) as f32) / (SQRT_2 * cell_m))
            } else {
                (r, (s1 * s1 + s2 * s2).sqrt())
            };
            if s > 0.0 && best.is_none_or(|b| s > b.0) {
                best = Some((s, kc, kd, r));
            }
        }
        let inherited = best_dev[c];
        let chosen = best.and_then(|(_, kc, kd, r)| {
            let (q1, q2) = ((D8_DX[kc] as f32, D8_DY[kc] as f32), (D8_DX[kd] as f32, D8_DY[kd] as f32));
            // the theoretical direction: from the cardinal toward the diagonal side by r
            let p = (q2.0 - q1.0, q2.1 - q1.1);
            let t = (r.cos() * q1.0 + r.sin() * p.0, r.cos() * q1.1 + r.sin() * p.1);
            let cross = |q: (f32, f32)| q.0 * t.1 - q.1 * t.0;
            let (c1, c2) = (inherited + cross(q1), inherited + cross(q2));
            let (n1, n2) = (nb(c, kc), nb(c, kd));
            let ok = |m: usize| on(m) && (sink[m] || zm(m) < e0);
            // ADR Finding 128-C0 -- [O14] eq. (5)'s "≤" gives an EXACT tie to the cardinal. Rounding
            // must not break it cell by cell (a slope of 1:4 ties at every fourth step, and f32 noise
            // then dephased neighbouring paths into 137 confluences on a 30² plane): ties within
            // 1e-4 cell go to the cardinal.
            const TIE: f32 = 1e-4;
            let cardinal = c1.abs() <= c2.abs() + TIE;
            // ADR Finding 129-B2 -- whether the tie rule decided (and whether the tie was exact)
            let tie = (c1.abs() - c2.abs()).abs() <= TIE;
            let exact = c1.abs() == c2.abs();
            let first = if cardinal { (kc, n1, c1) } else { (kd, n2, c2) };
            let second = if cardinal { (kd, n2, c2) } else { (kc, n1, c1) };
            let pick = if ok(first.1) {
                Some(first)
            } else if ok(second.1) {
                Some(second)
            } else {
                None
            };
            pick.map(|p| (p, tie, exact))
        });
        let chosen = chosen.map(|(p, tie, exact)| {
            if tie {
                ties += 1;
                tie_cells.push(c as u32);
            }
            if exact {
                exact_ties += 1;
            }
            p
        });
        let (k, r, conveyed) = match chosen {
            Some((k, m, dev)) => (k as u8, m, dev),
            None => {
                flats += 1;
                let k = fallback[c];
                if k == DIR_NONE {
                    continue;
                }
                (k, nb(c, k as usize), 0.0)
            }
        };
        dir[c] = k;
        if area[c] > best_area[r] {
            best_area[r] = area[c];
            best_dev[r] = conveyed;
        }
        area[r] += area[c];
    }
    (dir, LtdStats { flats, exact_ties, ties, tie_cells })
}

/// Drained area in CELLS along `dir` (each cell counts 1), by Kahn's order on the receiver graph;
/// `sink` cells receive but do not pass on.
pub fn accumulate_cells(dir: &[u8], sink: &[bool], w: usize, h: usize) -> Vec<f32> {
    let n = w * h;
    let recv = |c: usize| -> Option<usize> {
        let k = dir[c];
        if k == DIR_NONE || sink[c] {
            return None;
        }
        let x = ((c % w) as i32 + D8_DX[k as usize]).rem_euclid(w as i32) as usize;
        let y = ((c / w) as i32 + D8_DY[k as usize]).rem_euclid(h as i32) as usize;
        Some(y * w + x)
    };
    let mut indeg = vec![0u32; n];
    for c in 0..n {
        if let Some(r) = recv(c) {
            indeg[r] += 1;
        }
    }
    let mut acc = vec![1f32; n];
    let mut stack: Vec<usize> = (0..n).filter(|&c| indeg[c] == 0).collect();
    while let Some(c) = stack.pop() {
        if let Some(r) = recv(c) {
            acc[r] += acc[c];
            indeg[r] -= 1;
            if indeg[r] == 0 {
                stack.push(r);
            }
        }
    }
    acc
}

/// A trunk line's samples: its D8 cells UNWRAPPED along the torus (cell centres at +0.5), each with
/// its floor altitude (m) and half floor width (m).
fn line_samples(
    line: &[usize],
    w: usize,
    h: usize,
    base: &[f32],
    chi: &[f32],
    area_km2: &[f32],
    vc: &ValleyConstruction,
) -> Vec<(f32, f32, f32, f32)> {
    let mut pts: Vec<(f32, f32, f32, f32)> = Vec::with_capacity(line.len());
    let (mut px, mut py) = (0f32, 0f32);
    for (i, &c) in line.iter().enumerate() {
        let (mut x, mut y) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
        if i > 0 {
            while x - px > w as f32 / 2.0 {
                x -= w as f32;
            }
            while px - x > w as f32 / 2.0 {
                x += w as f32;
            }
            while y - py > h as f32 / 2.0 {
                y -= h as f32;
            }
            while py - y > h as f32 / 2.0 {
                y += h as f32;
            }
        }
        px = x;
        py = y;
        let zf = base[c] + vc.age_k * chi[c];
        pts.push((x, y, zf, 0.5 * vc.width_m(area_km2[c])));
    }
    pts
}

/// ADR Finding 130-A — the index, in [`densify`]'s output, of each input sample (same step rule).
fn densify_index(sm: &[(f32, f32, f32, f32)]) -> Vec<usize> {
    let mut at = Vec::with_capacity(sm.len());
    let mut k = 0usize;
    for i in 0..sm.len() {
        at.push(k);
        k += 1;
        if i + 1 < sm.len() {
            let (a, b) = (sm[i], sm[i + 1]);
            let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
            let steps = (len / 0.5).ceil() as usize;
            k += steps.saturating_sub(1);
        }
    }
    at
}

/// Densify a polyline to ≤ 0.5 cell between samples (linear in position, floor and width).
fn densify(sm: Vec<(f32, f32, f32, f32)>) -> Vec<(f32, f32, f32, f32)> {
    let mut dense = Vec::with_capacity(sm.len() * 2);
    for i in 0..sm.len() {
        dense.push(sm[i]);
        if i + 1 < sm.len() {
            let (a, b) = (sm[i], sm[i + 1]);
            let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
            let steps = (len / 0.5).ceil() as usize;
            for s in 1..steps {
                let t = s as f32 / steps as f32;
                dense.push((
                    a.0 + t * (b.0 - a.0),
                    a.1 + t * (b.1 - a.1),
                    a.2 + t * (b.2 - a.2),
                    a.3 + t * (b.3 - a.3),
                ));
            }
        }
    }
    dense
}

/// ADR Finding 127-B — the value (m) and gradient (m/m) of `zm` (metres) at continuous cell
/// coordinates `(x, y)` (cell centres at +0.5), on the torus.
pub fn interp_grad(
    zm: &[f32],
    w: usize,
    h: usize,
    cell_m: f32,
    how: TraceInterp,
    x: f32,
    y: f32,
) -> (f32, f32, f32) {
    let at = |ix: i64, iy: i64| {
        zm[iy.rem_euclid(h as i64) as usize * w + ix.rem_euclid(w as i64) as usize]
    };
    let (fx, fy) = (x - 0.5, y - 0.5);
    let (x0, y0) = (fx.floor(), fy.floor());
    let (tx, ty) = (fx - x0, fy - y0);
    let (ix, iy) = (x0 as i64, y0 as i64);
    match how {
        TraceInterp::Bilinear => {
            let (z00, z10, z01, z11) = (at(ix, iy), at(ix + 1, iy), at(ix, iy + 1), at(ix + 1, iy + 1));
            let v = (z00 * (1.0 - tx) + z10 * tx) * (1.0 - ty) + (z01 * (1.0 - tx) + z11 * tx) * ty;
            let gx = ((z10 - z00) * (1.0 - ty) + (z11 - z01) * ty) / cell_m;
            let gy = ((z01 - z00) * (1.0 - tx) + (z11 - z10) * tx) / cell_m;
            (v, gx, gy)
        }
        TraceInterp::Bicubic => {
            let cubic = |t: f32| -> ([f32; 4], [f32; 4]) {
                let (t2, t3) = (t * t, t * t * t);
                (
                    [
                        0.5 * (-t3 + 2.0 * t2 - t),
                        0.5 * (3.0 * t3 - 5.0 * t2 + 2.0),
                        0.5 * (-3.0 * t3 + 4.0 * t2 + t),
                        0.5 * (t3 - t2),
                    ],
                    [
                        0.5 * (-3.0 * t2 + 4.0 * t - 1.0),
                        0.5 * (9.0 * t2 - 10.0 * t),
                        0.5 * (-9.0 * t2 + 8.0 * t + 1.0),
                        0.5 * (3.0 * t2 - 2.0 * t),
                    ],
                )
            };
            let ((wx, dwx), (wy, dwy)) = (cubic(tx), cubic(ty));
            let (mut v, mut gx, mut gy) = (0f32, 0f32, 0f32);
            for (j, (&wyj, &dwyj)) in wy.iter().zip(&dwy).enumerate() {
                for (i, (&wxi, &dwxi)) in wx.iter().zip(&dwx).enumerate() {
                    let z = at(ix - 1 + i as i64, iy - 1 + j as i64);
                    v += wxi * wyj * z;
                    gx += dwxi * wyj * z;
                    gy += wxi * dwyj * z;
                }
            }
            (v, gx / cell_m, gy / cell_m)
        }
    }
}

enum Outcome {
    Same,
    Other,
    Lost,
}

/// One retraced line: its polyline, its outcome (`None` = not retraced), continuous steps, flat
/// steps, and its samples' deviations from its D8 path (in W).
type Traced = (Vec<(f32, f32, f32, f32)>, Option<Outcome>, usize, usize, Vec<f32>, Vec<(f32, f32)>);

/// ADR Finding 127-B — the polylines with every `a_min`–`below_km2` segment's median line retraced
/// off the D8 lattice, and what the tracé did.
#[allow(clippy::too_many_arguments)]
fn traced_polylines(
    raw_lines: &[Vec<usize>],
    field: &GridF32,
    ss: &SteinSteinParams,
    cell_m: f32,
    base: &[f32],
    chi: &[f32],
    area_km2: &[f32],
    vc: &ValleyConstruction,
    tr: SkeletonTrace,
) -> (Vec<Vec<(f32, f32, f32, f32)>>, TraceStats) {
    use rayon::prelude::*;
    use std::collections::HashSet;
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let zm: Vec<f32> = field.data.par_iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
    let sea_m = c1_altitude_norm_to_metres(C1_SEA_LEVEL_NORM, ss);
    // the owner of every trunk cell: the first line to visit it (a later line ENDS on it: a junction)
    let mut owner = vec![u32::MAX; n];
    for (li, l) in raw_lines.iter().enumerate() {
        for &c in l {
            if owner[c] == u32::MAX {
                owner[c] = li as u32;
            }
        }
    }
    // the junction a line ends on: (the owning line, the index of that cell in it)
    let parent: Vec<Option<(u32, u32)>> = raw_lines
        .iter()
        .enumerate()
        .map(|(li, l)| {
            let last = *l.last().expect("a line has ≥ 2 cells");
            let o = owner[last];
            (o != li as u32).then(|| {
                let j = raw_lines[o as usize]
                    .iter()
                    .position(|&c| c == last)
                    .expect("the owner holds its cell");
                (o, j as u32)
            })
        })
        .collect();
    let d8_poly = |line: &[usize]| {
        let pts = line_samples(line, w, h, base, chi, area_km2, vc);
        densify(smooth_positions(pts, line, area_km2, vc, cell_m))
    };
    let torus = |a: f32, b: f32, len: usize| -> f32 {
        let mut d = (a - b).rem_euclid(len as f32);
        if d > len as f32 / 2.0 {
            d -= len as f32;
        }
        d
    };
    let near = |c: usize, px: f32, py: f32| -> (f32, f32) {
        let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
        (px + torus(cx, px, w), py + torus(cy, py, h))
    };
    let results: Vec<Traced> = raw_lines
        .par_iter()
        .enumerate()
        .map(|(li, line)| {
            let me = li as u32;
            let own = |c: usize| owner[c] == me;
            let k_end = line.iter().take_while(|&&c| own(c) && area_km2[c] < tr.below_km2).count();
            if k_end < 2 {
                return (d8_poly(line), None, 0, 0, Vec::new(), Vec::new());
            }
            // the receiver: this line's own continuation (A reaches below_km2), the line it joins,
            // or its end (a terminal, the sea)
            let recv_own = k_end < line.len() && own(line[k_end]);
            let recv_parent = if recv_own { None } else { parent[li].map(|(p, _)| p) };
            let prefix: HashSet<usize> = line[..k_end].iter().copied().collect();
            let upstream = |mut l: u32| -> bool {
                for _ in 0..100_000 {
                    match parent[l as usize] {
                        Some((p, j)) if p == me => return (j as usize) < k_end,
                        Some((p, _)) => l = p,
                        None => return false,
                    }
                }
                false
            };
            let c0 = line[0];
            let mut p = ((c0 % w) as f32 + 0.5, (c0 / w) as f32 + 0.5);
            let mut pts: Vec<(f32, f32)> = vec![p];
            let (mut prog, mut steps, mut flat) = (0usize, 0usize, 0usize);
            let mut dev: Vec<f32> = Vec::new();
            let max_steps = 16 * k_end + 64;
            let r_max = 6i64;
            let last = line[line.len() - 1];
            let outcome = loop {
                if steps + flat > max_steps {
                    break Outcome::Lost;
                }
                let (cx, cy) = (p.0.rem_euclid(w as f32), p.1.rem_euclid(h as f32));
                let (ix, iy) = (cx.floor() as i64, cy.floor() as i64);
                let (mut hit_recv, mut hit_other) = (false, false);
                for dy in -r_max..=r_max {
                    for dx in -r_max..=r_max {
                        let c = (iy + dy).rem_euclid(h as i64) as usize * w
                            + (ix + dx).rem_euclid(w as i64) as usize;
                        let o = owner[c];
                        if o == u32::MAX || (o == me && prefix.contains(&c)) {
                            continue;
                        }
                        let (ccx, ccy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
                        let (ddx, ddy) = (torus(cx, ccx, w), torus(cy, ccy, h));
                        if (ddx * ddx + ddy * ddy).sqrt() * cell_m > 0.5 * vc.width_m(area_km2[c]) {
                            continue;
                        }
                        if (recv_own && o == me) || recv_parent == Some(o) {
                            hit_recv = true;
                        } else if o != me && !upstream(o) {
                            hit_other = true;
                        }
                    }
                }
                if hit_recv {
                    break Outcome::Same;
                }
                if hit_other && steps + flat >= 2 {
                    break Outcome::Other;
                }
                let (zp, gx, gy) = interp_grad(&zm, w, h, cell_m, tr.interp, cx, cy);
                if !recv_own && recv_parent.is_none() {
                    // the end of a line with no receiver: its terminal cell, or the sea
                    let (lx, ly) = near(last, p.0, p.1);
                    let dl = ((lx - p.0).powi(2) + (ly - p.1).powi(2)).sqrt() * cell_m;
                    if zp <= sea_m || dl <= (0.5 * vc.width_m(area_km2[last])).max(cell_m) {
                        break Outcome::Same;
                    }
                } else if zp <= sea_m {
                    break Outcome::Other;
                }
                // progress on the own D8 prefix: the nearest of its next cells
                let hi = (prog + 40).min(k_end);
                let (mut bj, mut bd) = (prog, f32::INFINITY);
                for (j, &c) in line[prog..hi].iter().enumerate() {
                    let (qx, qy) = near(c, p.0, p.1);
                    let d = (qx - p.0).powi(2) + (qy - p.1).powi(2);
                    if d < bd {
                        bd = d;
                        bj = prog + j;
                    }
                }
                prog = bj;
                dev.push(bd.sqrt() * cell_m / vc.width_m(area_km2[line[prog]]).max(1e-3));
                let gn = (gx * gx + gy * gy).sqrt();
                if gn < tr.flat_slope {
                    // a flat: one step along the D8 path, onto the prefix cell after the nearest one
                    let j = (prog + 1).min(k_end - 1);
                    p = near(line[j], p.0, p.1);
                    pts.push(p);
                    prog = j;
                    flat += 1;
                    if j == k_end - 1 {
                        break Outcome::Same; // the D8 path's next cell IS the receiver
                    }
                    continue;
                }
                p = (p.0 - 0.5 * gx / gn, p.1 - 0.5 * gy / gn);
                pts.push(p);
                steps += 1;
            };
            match outcome {
                Outcome::Same => {
                    // floor and width by arclength FRACTION along the D8 prefix (monotone; both ends
                    // keep their D8 values)
                    let d8s = line_samples(&line[..k_end], w, h, base, chi, area_km2, vc);
                    let cum = |q: &[(f32, f32)]| -> Vec<f32> {
                        let mut a = vec![0f32; q.len()];
                        for i in 1..q.len() {
                            a[i] = a[i - 1]
                                + ((q[i].0 - q[i - 1].0).powi(2) + (q[i].1 - q[i - 1].1).powi(2))
                                    .sqrt();
                        }
                        a
                    };
                    let a8 = cum(&d8s.iter().map(|s| (s.0, s.1)).collect::<Vec<_>>());
                    let at = cum(&pts);
                    let l8 = *a8.last().expect("≥ 2 samples");
                    let lt = at.last().copied().unwrap_or(0.0).max(1e-6);
                    let mut out: Vec<(f32, f32, f32, f32)> = Vec::with_capacity(pts.len() + line.len());
                    for (i, q) in pts.iter().enumerate() {
                        let s = at[i] / lt * l8;
                        let k = a8.partition_point(|&v| v <= s).clamp(1, a8.len() - 1);
                        let t = ((s - a8[k - 1]) / (a8[k] - a8[k - 1]).max(1e-6)).clamp(0.0, 1.0);
                        let (za, zb) = (d8s[k - 1], d8s[k]);
                        out.push((q.0, q.1, za.2 + t * (zb.2 - za.2), za.3 + t * (zb.3 - za.3)));
                    }
                    let end = *out.last().expect("a trace has a sample");
                    let mut tail = if recv_own {
                        let rest = &line[k_end..];
                        let rs = line_samples(rest, w, h, base, chi, area_km2, vc);
                        smooth_positions(rs, rest, area_km2, vc, cell_m)
                    } else {
                        line_samples(&line[line.len() - 1..], w, h, base, chi, area_km2, vc)
                    };
                    let sx = end.0 + torus(tail[0].0, end.0, w) - tail[0].0;
                    let sy = end.1 + torus(tail[0].1, end.1, h) - tail[0].1;
                    for q in &mut tail {
                        q.0 += sx;
                        q.1 += sy;
                    }
                    out.extend(tail);
                    (densify(out), Some(Outcome::Same), steps, flat, dev, pts)
                }
                o => (d8_poly(line), Some(o), steps, flat, Vec::new(), pts),
            }
        })
        .collect();
    let mut st = TraceStats::default();
    let mut polys = Vec::with_capacity(results.len());
    for (poly, o, steps, flat, dev, path) in results {
        polys.push(poly);
        let Some(o) = o else { continue };
        st.segments += 1;
        st.steps += steps;
        st.flat_steps += flat;
        st.paths.push((matches!(o, Outcome::Same), path));
        match o {
            Outcome::Same => {
                st.same_receiver += 1;
                st.dev_w.extend(dev);
            }
            Outcome::Other => st.other_receiver += 1,
            Outcome::Lost => st.lost += 1,
        }
    }
    (polys, st)
}

/// The moving average of a trunk line's POSITIONS (floor altitude and width untouched). `line` is
/// the line's D8 cells, `pts` their unwrapped samples, one per cell.
///
/// `smooth_w` `None` (Findings 121–125, byte-identical): a fixed half-window of `smooth_m`, pinned
/// over `half` samples at each end, and a line shorter than `2·half + 3` samples left as D8 drew
/// it. ADR Finding 126-D, `Some(c)`: the half-window at sample `i` is `c·W(A_i)` metres, shrunk to
/// `min(i, len − 1 − i)` so the window stays symmetric; only the two endpoints are pinned.
fn smooth_positions(
    pts: Vec<(f32, f32, f32, f32)>,
    line: &[usize],
    area_km2: &[f32],
    vc: &ValleyConstruction,
    cell_m: f32,
) -> Vec<(f32, f32, f32, f32)> {
    let mean = |pts: &[(f32, f32, f32, f32)], i: usize, hk: usize| -> (f32, f32) {
        let (mut sx, mut sy) = (0f32, 0f32);
        for p in &pts[i - hk..=i + hk] {
            sx += p.0;
            sy += p.1;
        }
        let d = (2 * hk + 1) as f32;
        (sx / d, sy / d)
    };
    if let Some(c) = vc.smooth_w {
        let len = pts.len();
        let mut o = pts.clone();
        for i in 1..len.saturating_sub(1) {
            let want = (c * vc.width_m(area_km2[line[i]]) / cell_m).round().max(0.0) as usize;
            let hk = want.min(i).min(len - 1 - i);
            if hk > 0 {
                (o[i].0, o[i].1) = mean(&pts, i, hk);
            }
        }
        return o;
    }
    // moving average of the POSITION only, endpoints pinned
    let half_k = (vc.smooth_m / cell_m).round().max(0.0) as usize;
    if half_k > 0 && pts.len() >= 2 * half_k + 3 {
        let mut o = pts.clone();
        for i in half_k..pts.len() - half_k {
            (o[i].0, o[i].1) = mean(&pts, i, half_k);
        }
        o
    } else {
        pts
    }
}

/// ADR Finding 122-B -- a priority flood seeded on the OPEN-OCEAN cells only (`water_class` 1):
/// every other cell rises to the lowest col on its way to the ocean, so an inland below-sea basin
/// and a closed land depression both read their spill level. Not the drainage chain's fill, which
/// takes every `h <= sea` cell as base level and never sees an inland basin (Finding 85).
pub fn ocean_flood(bf: &GridF32) -> Vec<f32> {
    use crate::lakes::connectivity::{WATER_CLASS_OCEAN, water_class};
    let (w, h) = (bf.width, bf.height);
    let n = w * h;
    let wc = water_class(bf, C1_SEA_LEVEL_NORM);
    let mut f = vec![f32::INFINITY; n];
    let mut heap = BinaryHeap::new();
    for k in 0..n {
        if wc[k] == WATER_CLASS_OCEAN {
            f[k] = bf.data[k];
            heap.push(Item(bf.data[k], k));
        }
    }
    while let Some(Item(v, c)) = heap.pop() {
        if v > f[c] {
            continue;
        }
        let (x, y) = ((c % w) as i32, (c / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let nb = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                let nv = bf.data[nb].max(v);
                if nv < f[nb] {
                    f[nb] = nv;
                    heap.push(Item(nv, nb));
                }
            }
        }
    }
    // a grid with no open ocean (a test field): nothing is closed
    for (v, &b) in f.iter_mut().zip(bf.data.iter()) {
        if !v.is_finite() {
            *v = b;
        }
    }
    f
}

/// ADR Finding 124-B1 -- a separable box blur of radius `r` cells on the torus (the field's trend).
fn box_blur_torus(z: &[f32], w: usize, h: usize, r: i64) -> Vec<f32> {
    let pass = |src: &[f32], horizontal: bool| -> Vec<f32> {
        let mut out = vec![0f32; w * h];
        let (len, lines) = if horizontal { (w, h) } else { (h, w) };
        let at = |line: usize, i: i64| -> usize {
            let i = i.rem_euclid(len as i64) as usize;
            if horizontal { line * w + i } else { i * w + line }
        };
        let win = (2 * r + 1) as f64;
        for line in 0..lines {
            let mut acc = 0f64;
            for i in -r..=r {
                acc += src[at(line, i)] as f64;
            }
            for i in 0..len as i64 {
                out[at(line, i)] = (acc / win) as f32;
                acc += src[at(line, i + r + 1)] as f64 - src[at(line, i - r)] as f64;
            }
        }
        out
    };
    pass(&pass(z, true), false)
}

/// ADR Finding 127-A -- the detail's weight at `u` metres beyond the floor edge: 0 on the foot
/// (`u < foot_m`), 1 beyond `(1 + tau)·foot_m`, a smoothstep between. No foot (`foot_m ≤ 0`): 1.
fn foot_weight(u: f32, foot_m: f32, tau: f32) -> f32 {
    if foot_m <= 0.0 {
        return 1.0;
    }
    let t = ((u - foot_m) / (tau.max(1e-6) * foot_m)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// What the carve did, per cell.
pub struct CarveMasks {
    /// The cell was lowered by the construction.
    pub carved: Vec<bool>,
    /// The cell lies on a valley FLOOR (within `W/2` of the skeleton) and was lowered.
    pub floor: Vec<bool>,
}

#[derive(PartialEq)]
struct Item(f32, usize);
impl Eq for Item {}
impl PartialOrd for Item {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Item {
    // min-heap on the key (distance), deterministic tie-break on the cell index
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.total_cmp(&self.0).then_with(|| o.1.cmp(&self.1))
    }
}

/// Lay the valleys: `out = min(field, V)`, with `V(x) = z_floor(p) + max(0, |x − p| − W(p)/2)·tan θ`
/// and `p` the skeleton sample NEAREST to `x`. **It can only lower a cell** (asserted in the tests).
///
/// ⚠️ **Nearest, not lowest.** Taking the minimum of the cones over every sample lets a DOWNSTREAM
/// sample, whose floor is lower, flatten the floor over half a valley width upstream of it — the
/// long profile would no longer be the χ law. So the cross-section is laid from the nearest sample
/// (vector propagation of the nearest source by Euclidean distance), and the minimum of two cones
/// is taken only between samples of DIFFERENT polylines, where it keeps the surface continuous
/// across the line where two valleys meet. Propagation stops where a wall reaches the terrain, so
/// a valley never cuts through a ridge into the next basin.
pub fn carve(
    field: &GridF32,
    sk: &Skeleton,
    vc: &ValleyConstruction,
    ss: &SteinSteinParams,
) -> (GridF32, CarveMasks) {
    let (out, masks, _) = carve_diag(field, sk, vc, ss);
    (out, masks)
}

/// ADR Finding 130-A0 — what [`carve`] decided, per cell (diagnostic, benches only).
pub struct CarveDiag {
    /// The nearest sample (index over every polyline, in order) that reached the cell; `u32::MAX` if none.
    pub who: Vec<u32>,
    /// Under the confluence clause: the sample that laid the cell when a larger line took it.
    pub banded: Option<Vec<u32>>,
    /// The polyline of each sample, and the sample's index within its polyline.
    pub line_of: Vec<u32>,
    pub pos_of: Vec<u32>,
}

/// [`carve`], with what it decided per cell ([`CarveDiag`]). Byte-identical output.
pub fn carve_diag(
    field: &GridF32,
    sk: &Skeleton,
    vc: &ValleyConstruction,
    ss: &SteinSteinParams,
) -> (GridF32, CarveMasks, CarveDiag) {
    let (w, h) = (sk.width, sk.height);
    assert_eq!((field.width, field.height), (w, h), "skeleton and field sizes differ");
    let n = w * h;
    let tan = vc.wall_deg.to_radians().tan();
    let zfield: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
    let mut src: Vec<(f32, f32, f32, f32)> = Vec::new();
    let mut line_of: Vec<u32> = Vec::new();
    let mut pos_of: Vec<u32> = Vec::new();
    for (li, l) in sk.polylines.iter().enumerate() {
        for (pi, &p) in l.iter().enumerate() {
            src.push(p);
            line_of.push(li as u32);
            pos_of.push(pi as u32);
        }
    }
    let mut line_start: Vec<usize> = Vec::with_capacity(sk.polylines.len());
    let mut acc_len = 0usize;
    for l in &sk.polylines {
        line_start.push(acc_len);
        acc_len += l.len();
    }
    // ADR Finding 130-A2 -- the junction index on `big` of `small`'s downstream chain of lines
    let junction_on = |small: u32, big: u32| -> Option<u32> {
        let mut l = small;
        for _ in 0..sk.polylines.len() {
            match sk.line_parent.get(l as usize).copied().flatten() {
                Some((p, j)) if p == big => return Some(j),
                Some((p, _)) => l = p,
                None => return None,
            }
        }
        None
    };
    // (distance m, cone value m, on the floor?, u = distance beyond the floor edge m)
    let foot = vc.wall_profile.map_or(0.0, |p| p.foot_m.max(0.0));
    let geo = |s: usize, c: usize| -> (f32, f32, bool, f32) {
        let (sx, sy, zf, hw) = src[s];
        let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
        let mut dx = (cx - sx).rem_euclid(w as f32);
        if dx > w as f32 / 2.0 {
            dx -= w as f32;
        }
        let mut dy = (cy - sy).rem_euclid(h as f32);
        if dy > h as f32 / 2.0 {
            dy -= h as f32;
        }
        let d = (dx * dx + dy * dy).sqrt() * sk.cell_m;
        let u = (d - hw).max(0.0);
        // ADR Finding 124-B1 -- the concave foot: slope ramping 0 → tan over `foot` metres
        let rise = if foot > 0.0 && u < foot {
            tan * u * u / (2.0 * foot)
        } else {
            tan * (u - 0.5 * foot)
        };
        (d, zf + rise, d <= hw, u)
    };
    let mut bestd = vec![f32::INFINITY; n];
    let mut who = vec![u32::MAX; n];
    let mut heap = BinaryHeap::new();
    for (i, sp) in src.iter().enumerate() {
        let cx = (sp.0.floor() as i64).rem_euclid(w as i64) as usize;
        let cy = (sp.1.floor() as i64).rem_euclid(h as i64) as usize;
        let c = cy * w + cx;
        let (d, v, _, _) = geo(i, c);
        if d < bestd[c] && v < zfield[c] {
            bestd[c] = d;
            who[c] = i as u32;
            heap.push(Item(d, c));
        }
    }
    while let Some(Item(d, c)) = heap.pop() {
        if d > bestd[c] {
            continue;
        }
        let s = who[c] as usize;
        let (x, y) = ((c % w) as i32, (c / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let nb = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                let (dd, vv, _, _) = geo(s, nb);
                if dd < bestd[nb] && vv < zfield[nb] {
                    bestd[nb] = dd;
                    who[nb] = s as u32;
                    heap.push(Item(dd, nb));
                }
            }
        }
    }
    // ADR Finding 128-A -- the confluence clause. Per cell of any primitive's FLOOR BAND: the line
    // with the largest covering half-width (the largest drained area, W being monotone in A) and, of
    // THAT line, the nearest covering sample.
    let band: Option<(Vec<u32>, Vec<f32>)> = vc.trunk_band.then(|| {
        let mut b_line = vec![u32::MAX; n];
        let mut b_hw = vec![0f32; n];
        let mut b_smp = vec![u32::MAX; n];
        let mut b_d = vec![f32::INFINITY; n];
        for (i, &(sx, sy, _zf, hw)) in src.iter().enumerate() {
            let r = hw / sk.cell_m;
            let li = line_of[i];
            for yy in (sy - r).floor() as i64..=(sy + r).ceil() as i64 {
                for xx in (sx - r).floor() as i64..=(sx + r).ceil() as i64 {
                    let d = ((xx as f32 + 0.5 - sx).powi(2) + (yy as f32 + 0.5 - sy).powi(2)).sqrt() * sk.cell_m;
                    if d > hw {
                        continue;
                    }
                    let c = yy.rem_euclid(h as i64) as usize * w + xx.rem_euclid(w as i64) as usize;
                    if b_line[c] == li {
                        b_hw[c] = b_hw[c].max(hw);
                        if d < b_d[c] {
                            b_d[c] = d;
                            b_smp[c] = i as u32;
                        }
                    } else if b_line[c] == u32::MAX || hw > b_hw[c] {
                        (b_line[c], b_hw[c], b_smp[c], b_d[c]) = (li, hw, i as u32, d);
                    }
                }
            }
        }
        (b_smp, b_hw)
    });
    let trend: Option<Vec<f32>> = vc.wall_profile.filter(|p| p.detail_gain != 0.0).map(|p| {
        let r = (p.detail_radius_m / sk.cell_m).round().max(1.0) as i64;
        box_blur_torus(&zfield, w, h, r)
    });
    let sea_m = c1_altitude_norm_to_metres(C1_SEA_LEVEL_NORM, ss);
    let mut out = field.clone();
    let mut carved = vec![false; n];
    let mut floor = vec![false; n];
    let mut banded_by: Option<Vec<u32>> = band.as_ref().map(|_| vec![u32::MAX; n]);
    for c in 0..n {
        // ADR Finding 128-A -- a LARGER line's floor band takes the cell from a smaller line's sample
        let banded = band.as_ref().and_then(|(smp, bhw)| {
            let s = smp[c];
            let wins = s != u32::MAX
                && (who[c] == u32::MAX
                    || (line_of[s as usize] != line_of[who[c] as usize]
                        && bhw[c] > src[who[c] as usize].3
                        // ADR Finding 130-A2 -- only at or downstream of the junction
                        && (!vc.trunk_band_downstream
                            || junction_on(line_of[who[c] as usize], line_of[s as usize])
                                .is_some_and(|j| pos_of[s as usize] >= j))));
            wins.then_some(s as usize)
        });
        // ADR Finding 131-A -- the concordant confluence: only A0's cells (upstream of the owner's junction
        // with the larger line), laid at the larger line's floor AT the junction; no other takeover
        let (banded, concordant_floor) = match (banded, vc.trunk_band_concordant) {
            (Some(s), true) => {
                let o = who[c];
                let j = if o == u32::MAX { None } else { junction_on(line_of[o as usize], line_of[s]) };
                match j {
                    Some(j) if pos_of[s] < j => {
                        (Some(s), Some(src[line_start[line_of[s] as usize] + j as usize].2))
                    }
                    _ => (None, None),
                }
            }
            (b, _) => (b, None),
        };
        if who[c] == u32::MAX && banded.is_none() {
            continue;
        }
        if let (Some(b), Some(s)) = (banded_by.as_mut(), banded) {
            b[c] = s as u32;
        }
        let me = banded.unwrap_or(who[c] as usize);
        let (x, y) = ((c % w) as i32, (c / w) as i32);
        let (_, mut v, mut on_floor, mut u_sel) = geo(me, c);
        if let Some(zj) = concordant_floor {
            (v, on_floor, u_sel) = (zj, true, 0.0);
        }
        // continuity across the line where two DIFFERENT valleys meet (not where the clause decided:
        // the larger line lays its band, with no minimum)
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let nb = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                let o = who[nb];
                if banded.is_none() && o != u32::MAX && line_of[o as usize] != line_of[me] {
                    let (_, vv, f, uu) = geo(o as usize, c);
                    if vv < v {
                        v = vv;
                        on_floor = f;
                        u_sel = uu;
                    }
                }
            }
        }
        // ADR Finding 124-B1 -- the convex crest (a smooth minimum with the terrain) and the detail
        if let Some(p) = vc.wall_profile
            && !on_floor
            && v < zfield[c]
        {
            let k = p.crest_m.max(1e-3);
            v -= k * (1.0 + (-(zfield[c] - v) / k).exp()).ln();
            if let Some(t) = &trend {
                // ADR Finding 127-A -- bruit ↔ pied: no detail on the concave foot
                let quiet = p.foot_quiet.map_or(1.0, |tau| foot_weight(u_sel, p.foot_m, tau));
                v += p.detail_gain * quiet * (zfield[c] - t[c]);
            }
        }
        // ADR Finding 126-B -- mur ↔ mer: no WALL cell is laid below sea + ε
        if let Some(eps) = vc.wall_sea_floor_m
            && !on_floor
        {
            v = v.max(sea_m + eps);
        }
        // ADR Finding 140 -- the rim: a present lake's outer ring keeps >= its level L(r_lake)
        if let Some(rf) = sk.rim_floor_m.as_deref()
            && rf[c].is_finite()
        {
            v = v.max(rf[c]);
        }
        if v < zfield[c] {
            let nv = c1_metres_to_altitude_norm(v, ss);
            if nv < out.data[c] {
                out.data[c] = nv;
                carved[c] = true;
                floor[c] = on_floor;
            }
        }
    }
    (out, CarveMasks { carved, floor }, CarveDiag { who, banded: banded_by, line_of, pos_of })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Finding 124-5 (E) — W(k) is gated and scales as `(k/k₀)^γ`: `None` and `Some(0)` give
    /// Findings 121–124's width at every age; `Some(1)` at k × 1.4 gives 1.4× that width, `Some(0.5)`
    /// √1.4×; at k₀ every γ gives the same width.
    #[test]
    fn the_width_widens_with_age_only_when_asked() {
        let at = |mult: f32, g: Option<f32>| {
            ValleyConstruction { width_age_gamma: g, ..ValleyConstruction::new(F121_AGE_K * mult, None) }
                .width_m(10.0)
        };
        let w0 = ValleyConstruction::new(F121_AGE_K, None).width_m(10.0);
        assert!((w0 - 200.0).abs() < 1e-3, "W(10 km²) = 200 m, got {w0}");
        for mult in [0.7f32, 1.0, 1.4] {
            assert_eq!(at(mult, None), w0, "no W(k): the width ignores the age");
            assert!((at(mult, Some(0.0)) - w0).abs() < 1e-4, "γ = 0 is the ungated width");
        }
        assert!((at(1.4, Some(1.0)) / w0 - 1.4).abs() < 1e-4, "γ = 1: W ∝ k");
        assert!((at(1.4, Some(0.5)) / w0 - 1.4f32.sqrt()).abs() < 1e-4, "γ = 0.5: W ∝ √k");
        assert!((at(1.0, Some(1.0)) - w0).abs() < 1e-3, "at k₀ every γ gives W₀");
        // gated in the serialised config (the cache digest): absent unless asked
        let js = serde_json::to_string(&ValleyConstruction::new(F121_AGE_K, None)).unwrap();
        assert!(!js.contains("width_age_gamma"), "{js}");
    }

    /// A 96×96 V-shaped basin draining north into a sea band: the flow converges on the centre
    /// column, which becomes the trunk.
    fn basin(ss: &SteinSteinParams) -> GridF32 {
        let (w, h) = (96usize, 96usize);
        let mut d = vec![0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let m = if y < 6 {
                    -50.0
                } else {
                    40.0 + 6.0 * (y as f32 - 6.0) + 9.0 * (x as f32 - 48.0).abs()
                };
                d[y * w + x] = c1_metres_to_altitude_norm(m, ss);
            }
        }
        GridF32 { width: w, height: h, data: d }
    }

    /// ADR Finding 124-B1, rule 13 — the wall profile only LOWERS, never raises, and changes the
    /// walls (negative control: it must not be a no-op) while leaving the floor cells alone.
    #[test]
    fn the_wall_profile_lowers_only_and_bends_the_walls() {
        let ss = SteinSteinParams::default();
        let f = basin(&ss);
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.a_min_km2 = 0.5;
        vc.smooth_m = 0.0;
        let sk = skeleton(&f, &vc, &ss, 4.8);
        let (plain, m0) = carve(&f, &sk, &vc, &ss);
        vc.wall_profile = Some(WallProfile::f124());
        let (bent, m1) = carve(&f, &sk, &vc, &ss);
        for k in 0..f.data.len() {
            assert!(bent.data[k] <= f.data[k], "cell {k} RAISED by the wall profile");
        }
        let walls_changed = (0..f.data.len())
            .filter(|&k| m0.carved[k] && !m0.floor[k] && bent.data[k] != plain.data[k])
            .count();
        assert!(
            walls_changed > 20,
            "negative control: the profile must bend the walls ({walls_changed})"
        );
        for k in 0..f.data.len() {
            if m0.floor[k] && m1.floor[k] {
                assert_eq!(bent.data[k], plain.data[k], "floor cell {k} moved");
            }
        }
    }

    /// ADR Finding 126-B, rule 13 — the coast clause keeps every WALL cell at or above `sea + ε`,
    /// and changes nothing else: every cell the unclamped carve left above `sea + ε` is bit-identical.
    /// Negative control FIRST: without the clause, the profile and its signed detail lay at least one
    /// wall cell below the sea (a 60 m dimple on a wall near the sea band, which the detail puts back).
    #[test]
    fn the_coast_clause_keeps_every_wall_above_the_sea() {
        let ss = SteinSteinParams::default();
        let mut f = basin(&ss);
        let w = f.width;
        for y in 0..f.height {
            for x in 0..w {
                let d2 = ((x as f32 - 52.0).powi(2) + (y as f32 - 9.0).powi(2)) / 4.0;
                if d2 < 9.0 {
                    let k = y * w + x;
                    let m = c1_altitude_norm_to_metres(f.data[k], &ss) - 60.0 * (-d2).exp();
                    f.data[k] = c1_metres_to_altitude_norm(m, &ss);
                }
            }
        }
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.a_min_km2 = 0.5;
        vc.smooth_m = 0.0;
        vc.wall_profile = Some(WallProfile::f124());
        let sk = skeleton(&f, &vc, &ss, 4.8);
        let (free, m0) = carve(&f, &sk, &vc, &ss);
        let sea = c1_altitude_norm_to_metres(C1_SEA_LEVEL_NORM, &ss);
        let land_in = |k: usize| f.data[k] > C1_SEA_LEVEL_NORM;
        let drowned = (0..f.data.len())
            .filter(|&k| land_in(k) && m0.carved[k] && !m0.floor[k] && free.data[k] <= C1_SEA_LEVEL_NORM)
            .count();
        assert!(drowned > 0, "negative control: without the clause a wall must drown ({drowned})");
        let eps = 0.5;
        let (held, m1) = carve(&f, &sk, &ValleyConstruction { wall_sea_floor_m: Some(eps), ..vc }, &ss);
        for k in 0..f.data.len() {
            if land_in(k) && m1.carved[k] && !m1.floor[k] {
                let z = c1_altitude_norm_to_metres(held.data[k], &ss);
                assert!(z >= sea + eps - 1e-3, "wall cell {k} at {z} m, below sea + ε");
            }
            assert!(held.data[k] <= f.data[k], "cell {k} RAISED by the clause");
            if c1_altitude_norm_to_metres(free.data[k], &ss) >= sea + eps {
                assert_eq!(held.data[k], free.data[k], "cell {k} moved though it was above sea + ε");
            }
        }
    }

    /// ADR Finding 126-D, rule 13 — the width law smooths a line the fixed half-window leaves as
    /// D8 drew it. A 9-sample D8 staircase: under `smooth_m` 250 m at 50 m cells (half 5, so a line
    /// needs 13 samples) it is untouched — the negative control; under `smooth_w` 2 (W = 100 m at
    /// 1 km², a 200 m = 4-cell half-window) its interior moves and its two endpoints do not. A
    /// straight line stays exactly straight under both.
    #[test]
    fn the_width_law_smooths_a_line_the_fixed_window_leaves_raw() {
        let vc = ValleyConstruction::f121(0.02, None); // W(1 km²) = 200/10^0.3 ≈ 100 m
        let cell_m = 50.0;
        let stair: Vec<(f32, f32, f32, f32)> = (0..9)
            .map(|i| (0.5 + (i / 2) as f32, 0.5 + i as f32, 10.0, 50.0))
            .collect();
        let line: Vec<usize> = (0..9).collect();
        let area = vec![1.0f32; 9];
        let fixed = smooth_positions(stair.clone(), &line, &area, &vc, cell_m);
        assert_eq!(fixed, stair, "negative control: the fixed window must leave a 9-sample line raw");
        let wlaw = smooth_positions(stair.clone(), &line, &area, &ValleyConstruction { smooth_w: Some(2.0), ..vc }, cell_m);
        assert_eq!((wlaw[0], wlaw[8]), (stair[0], stair[8]), "the endpoints stay pinned");
        let moved = (1..8).filter(|&i| wlaw[i].0 != stair[i].0).count();
        assert!(moved >= 3, "the width law must smooth the interior ({moved} moved)");
        for i in 0..9 {
            assert_eq!((wlaw[i].2, wlaw[i].3), (stair[i].2, stair[i].3), "only the position moves");
        }
        let straight: Vec<(f32, f32, f32, f32)> = (0..30).map(|i| (0.5, 0.5 + i as f32, 1.0, 50.0)).collect();
        let l30: Vec<usize> = (0..30).collect();
        let a30 = vec![5.0f32; 30];
        for v in [vc, ValleyConstruction { smooth_w: Some(2.0), ..vc }] {
            let s = smooth_positions(straight.clone(), &l30, &a30, &v, cell_m);
            for (a, b) in s.iter().zip(&straight) {
                assert!((a.0 - b.0).abs() < 1e-5 && (a.1 - b.1).abs() < 1e-4, "a straight line bent");
            }
        }
    }

    /// ADR Finding 127-A, rule 13 — the foot clause silences the detail on the concave foot only.
    /// The weight is 0 on the foot, 1 beyond `(1 + τ)·foot`, monotone between. Negative control
    /// first: with the clause the carve CHANGES some wall cells. It never raises a cell, and without
    /// a foot (B1b's planar wall) it is a bit-exact no-op.
    #[test]
    fn the_foot_clause_silences_the_detail_on_the_foot_only() {
        assert_eq!(foot_weight(0.0, 300.0, 0.5), 0.0);
        assert_eq!(foot_weight(300.0, 300.0, 0.5), 0.0);
        assert_eq!(foot_weight(450.0, 300.0, 0.5), 1.0);
        assert_eq!(foot_weight(900.0, 300.0, 0.5), 1.0);
        assert_eq!(foot_weight(10.0, 0.0, 0.5), 1.0, "no foot: full detail");
        let mut prev = 0.0;
        for i in 0..100 {
            let v = foot_weight(250.0 + 3.0 * i as f32, 300.0, 0.5);
            assert!(v >= prev, "the weight must be monotone");
            prev = v;
        }
        let ss = SteinSteinParams::default();
        let f = basin(&ss);
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.a_min_km2 = 0.5;
        vc.smooth_m = 0.0;
        let sk = skeleton(&f, &vc, &ss, 4.8);
        let prof = WallProfile::f124();
        let with = |wp: WallProfile| carve(&f, &sk, &ValleyConstruction { wall_profile: Some(wp), ..vc }, &ss);
        let (loud, m0) = with(prof);
        let (quiet, _) = with(WallProfile { foot_quiet: Some(0.5), ..prof });
        let changed = (0..f.data.len())
            .filter(|&k| m0.carved[k] && !m0.floor[k] && loud.data[k] != quiet.data[k])
            .count();
        assert!(changed > 0, "negative control: the clause must change the foot's walls ({changed})");
        for k in 0..f.data.len() {
            assert!(quiet.data[k] <= f.data[k], "cell {k} RAISED by the foot clause");
        }
        let planar = WallProfile { foot_m: 0.0, crest_m: 0.0, ..prof };
        let (a, _) = with(planar);
        let (b, _) = with(WallProfile { foot_quiet: Some(0.5), ..planar });
        assert_eq!(a.data, b.data, "no foot: the clause must be a bit-exact no-op");
    }

    /// A rounded V valley whose axis runs at 22.5° (between two D8 directions) down to a sea, 192²
    /// at 50 m: D8 draws its trunk as a staircase of E and SE steps.
    fn oblique_valley(ss: &SteinSteinParams) -> GridF32 {
        let n = 192usize;
        let th = 22.5f32.to_radians();
        let mut d = vec![0f32; n * n];
        for y in 0..n {
            for x in 0..n {
                let (u, v) = (x as f32 - 16.0, y as f32 - 16.0);
                let along = u * th.cos() + v * th.sin();
                let across = -u * th.sin() + v * th.cos();
                let m = if along > 150.0 {
                    -30.0
                } else {
                    300.0 - 2.0 * along + 8.0 * (across * across + 9.0).sqrt()
                };
                d[y * n + x] = c1_metres_to_altitude_norm(m, ss);
            }
        }
        GridF32 { width: n, height: n, data: d }
    }

    /// Per-cell-of-arclength chord angles of the longest polyline, degrees.
    fn chord_angles(sk: &Skeleton) -> Vec<f32> {
        let l = sk.polylines.iter().max_by_key(|l| l.len()).expect("a trunk");
        let mut out = Vec::new();
        let (mut i0, mut arc) = (0usize, 0f32);
        for i in 1..l.len() {
            arc += ((l[i].0 - l[i - 1].0).powi(2) + (l[i].1 - l[i - 1].1).powi(2)).sqrt();
            if arc >= 1.0 {
                out.push((l[i].1 - l[i0].1).atan2(l[i].0 - l[i0].0).to_degrees());
                i0 = i;
                arc = 0.0;
            }
        }
        out
    }

    /// ADR Finding 127-B, rule 13 — the off-grid tracé leaves the D8 lattice and keeps its receiver.
    /// On the 22.5° valley, every one-cell chord of the D8 trunk is at 0° or 45°, i.e. 22.5° off the
    /// axis (the negative control). The traced trunk (bicubic, the whole line retraced) stays within
    /// a few degrees of the axis, and every traced segment reaches its own receiver.
    #[test]
    fn the_trace_leaves_the_lattice_and_keeps_its_receiver() {
        let ss = SteinSteinParams::default();
        let f = oblique_valley(&ss);
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.a_min_km2 = 0.5;
        vc.smooth_m = 0.0;
        let dev = |a: &[f32]| a.iter().map(|&t| (t - 22.5).abs()).sum::<f32>() / a.len().max(1) as f32;
        let d8 = skeleton(&f, &vc, &ss, 9.6);
        assert!(d8.trace_stats.is_none(), "the tracé is off by default");
        let a8 = chord_angles(&d8);
        assert!(a8.len() > 40, "the trunk must be long ({})", a8.len());
        assert!(dev(&a8) > 15.0, "negative control: D8 must sit on its lattice ({:.1}°)", dev(&a8));
        let tr = SkeletonTrace { below_km2: 1.0e6, ..SkeletonTrace::f127(TraceInterp::Bicubic) };
        let tk = skeleton(&f, &ValleyConstruction { skeleton_trace: Some(tr), ..vc }, &ss, 9.6);
        let st = tk.trace_stats.as_ref().expect("the tracé's stats");
        assert!(st.segments >= 1, "a segment must be retraced");
        assert_eq!(st.same_receiver, st.segments, "every trace must reach its receiver: {st:?}");
        let at = chord_angles(&tk);
        assert!(dev(&at) < 8.0, "the tracé must leave the lattice ({:.1}°)", dev(&at));
    }

    /// A hand-made skeleton on a flat 500 m field, 64² at 50 m. A TRUNK along y = 32.5 (half-width 100 m,
    /// floor 100 → 131 m upstream) and a TRIBUTARY down x = 30.5 into it (half-width 25 m, floor 300 m at
    /// its source → 140 m at the trunk's axis), whose last samples lie INSIDE the trunk's floor band.
    fn confluence() -> (GridF32, Skeleton, SteinSteinParams) {
        let ss = SteinSteinParams::default();
        let n = 64usize;
        let f = GridF32 { width: n, height: n, data: vec![c1_metres_to_altitude_norm(500.0, &ss); n * n] };
        let trunk: Vec<(f32, f32, f32, f32)> = (0..119).map(|i| {
            let x = 2.5 + 0.5 * i as f32;
            (x, 32.5, 100.0 + 0.5 * (61.5 - x), 100.0)
        }).collect();
        let trib: Vec<(f32, f32, f32, f32)> = (0..61).map(|i| {
            let y = 2.5 + 0.5 * i as f32;
            (30.5, y, 300.0 - 160.0 * (y - 2.5) / 30.0, 25.0)
        }).collect();
        let sk = Skeleton {
            width: n,
            height: n,
            cell_m: 50.0,
            area_km2: vec![0.0; n * n],
            chi_m: vec![0.0; n * n],
            base_alt_m: vec![0.0; n * n],
            trunk: vec![false; n * n],
            polylines: vec![trunk, trib],
            trace_stats: None,
            direction: vec![DIR_NONE; n * n],
            ltd_flat_cells: 0,
            // the tributary (line 1) ends on the trunk's sample at x = 30.5 (index 56)
            line_parent: vec![None, Some((0, 56))],
            rim_floor_m: None,
            gorge_falls: Vec::new(),
            gorge_bodies: Vec::new(),
            gorge_body_of: None,
            gorge_path: Vec::new(),
        };
        (f, sk, ss)
    }

    /// ADR Finding 128-A, rule 13 — the confluence clause lets the trunk lay its whole floor band. Negative
    /// control FIRST: without it, the tributary's nearer samples lay some trunk-band cells ABOVE the trunk's
    /// floor (a dam). With it, every trunk-band cell lies at the trunk's own floor, cells away from the
    /// junction are bit-identical, and the tributary's own floor upstream is untouched.
    #[test]
    fn the_confluence_clause_lets_the_trunk_lay_its_band() {
        let (f, sk, ss) = confluence();
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.wall_deg = 28.0;
        let (free, _) = carve(&f, &sk, &vc, &ss);
        let (held, mh) = carve(&f, &sk, &ValleyConstruction { trunk_band: true, ..vc }, &ss);
        let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
        let trunk_floor = |x: usize| 100.0 + 0.5 * (61.5 - (x as f32 + 0.5));
        let in_band = |x: usize, y: usize| ((y as f32 + 0.5) - 32.5).abs() * 50.0 <= 100.0 && (3..61).contains(&x);
        let mut dammed = 0;
        for y in 0..64 {
            for x in 0..64 {
                if in_band(x, y) && m(&free, y * 64 + x) > trunk_floor(x) + 1.0 {
                    dammed += 1;
                }
            }
        }
        assert!(dammed > 0, "negative control: without the clause the tributary must dam the trunk ({dammed})");
        for y in 0..64 {
            for x in 0..64 {
                let k = y * 64 + x;
                if in_band(x, y) {
                    assert!(m(&held, k) <= trunk_floor(x) + 0.3, "({x},{y}) {} m above the trunk's floor", m(&held, k));
                    assert!(mh.floor[k], "({x},{y}) must be the trunk's floor");
                }
                if (x as i32 - 30).abs() > 4 || y < 20 {
                    assert_eq!(held.data[k], free.data[k], "({x},{y}) moved away from the junction");
                }
                assert!(held.data[k] <= f.data[k], "({x},{y}) raised");
            }
        }
    }

    /// ADR Finding 130-A2, rule 13 — the RESTRICTED clause takes a smaller line's band cells only at or
    /// downstream of the junction. Negative control FIRST: the full clause re-lays cells whose nearest
    /// trunk sample lies UPSTREAM of the tributary's junction (index 56). Restricted: those cells are
    /// bit-identical to the construction without any clause; the ones at or downstream of the junction
    /// are laid exactly as by the full clause.
    #[test]
    fn the_restricted_clause_takes_the_band_only_downstream_of_the_junction() {
        let (f, sk, ss) = confluence();
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.wall_deg = 28.0;
        let (free, _) = carve(&f, &sk, &vc, &ss);
        let (full, _, df) = carve_diag(&f, &sk, &ValleyConstruction { trunk_band: true, ..vc }, &ss);
        let restricted = ValleyConstruction { trunk_band: true, trunk_band_downstream: true, ..vc };
        let (rest, _) = carve(&f, &sk, &restricted, &ss);
        let taken = df.banded.expect("the clause is on");
        let (mut upstream_moved, mut upstream, mut downstream) = (0, 0, 0);
        for k in 0..64 * 64 {
            let s = taken[k];
            // a trunk (line 0) sample took the cell from the tributary (line 1)
            if s == u32::MAX || df.line_of[s as usize] != 0 || df.who[k] == u32::MAX || df.line_of[df.who[k] as usize] != 1 {
                continue;
            }
            if df.pos_of[s as usize] < 56 {
                upstream += 1;
                upstream_moved += (full.data[k] != free.data[k]) as usize;
                assert_eq!(rest.data[k], free.data[k], "cell {k}: upstream of the junction the restricted clause acted");
            } else {
                downstream += 1;
                assert_eq!(rest.data[k], full.data[k], "cell {k}: downstream the restricted clause is not the full one");
            }
        }
        assert!(upstream_moved > 0, "negative control: the full clause must re-lay cells upstream of the junction ({upstream})");
        assert!(downstream > 0, "the fixture must have band cells downstream of the junction");
        for k in 0..64 * 64 {
            assert!(rest.data[k] <= f.data[k], "cell {k} raised above the field");
        }
    }

    /// ADR Finding 131-A, rule 13 — the CONCORDANT confluence lays the tributary's cells in the trunk's band,
    /// upstream of the junction, at the trunk's floor AT the junction (115.5 m, x = 30.5). Negative control
    /// FIRST: without any clause some of them dam the trunk; the full clause lays them higher than the
    /// junction's floor (from samples upstream). Concordant: every one of them at the junction's floor, and
    /// every other cell bit-identical to the construction without a clause.
    #[test]
    fn the_concordant_clause_lays_the_tributary_at_the_junction_floor() {
        let (f, sk, ss) = confluence();
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.wall_deg = 28.0;
        let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
        let (free, _) = carve(&f, &sk, &vc, &ss);
        let (full, _, df) = carve_diag(&f, &sk, &ValleyConstruction { trunk_band: true, ..vc }, &ss);
        let (conc, mc) = carve(&f, &sk, &ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..vc }, &ss);
        let zj = 100.0 + 0.5 * (61.5 - 30.5);
        let taken = df.banded.expect("the clause is on");
        let (mut a0, mut dammed, mut full_higher) = (0, 0, 0);
        for k in 0..64 * 64 {
            let s = taken[k];
            let a0_cell = s != u32::MAX
                && df.line_of[s as usize] == 0
                && df.who[k] != u32::MAX
                && df.line_of[df.who[k] as usize] == 1
                && df.pos_of[s as usize] < 56;
            if a0_cell {
                a0 += 1;
                dammed += (m(&free, k) > zj + 1.0) as usize;
                full_higher += (m(&full, k) > zj + 0.2) as usize;
                assert!((m(&conc, k) - zj).abs() < 0.2, "cell {k}: {} m, not the junction's floor {zj}", m(&conc, k));
                assert!(mc.floor[k], "cell {k} must be a floor cell");
            } else {
                assert_eq!(conc.data[k], free.data[k], "cell {k}: outside A0's cells the concordant clause acted");
            }
        }
        assert!(a0 > 0 && dammed > 0, "negative control: A0's cells exist ({a0}) and some dam the trunk without a clause ({dammed})");
        assert!(full_higher > 0, "negative control: the full clause lays some of them above the junction's floor ({full_higher})");
    }

    /// ADR Finding 132-P4, rule 13 — a present lake is the base level of its catchment. A plane falls 4 m per cell
    /// to the sea on 128² cells of 400 m; a bowl 150 m deep (a closed depression above the 5 km² lake floor) is set in it. Negative control FIRST: without the field, the breached bowl
    /// passes χ through, so a cell upstream takes the SEA's base. With it, that cell's base is the lake's surface
    /// and its χ counts from the lake only; a cell downstream of the lake keeps its base and χ bit for bit.
    #[test]
    fn a_present_lake_is_the_base_level_of_its_catchment() {
        let ss = SteinSteinParams::default();
        let n = 128usize;
        let z_m = |x: f32, y: f32| -> f32 {
            if y >= 120.0 {
                return -50.0;
            }
            let d = ((x - 64.0).powi(2) + (y - 50.0).powi(2)).sqrt();
            let bowl = if d < 12.0 { 150.0 * (1.0 - (d / 12.0).powi(2)) } else { 0.0 };
            5.0 + 4.0 * (120.0 - y) - bowl
        };
        let f = GridF32 {
            width: n,
            height: n,
            data: (0..n * n).map(|k| c1_metres_to_altitude_norm(z_m((k % n) as f32 + 0.5, (k / n) as f32 + 0.5), &ss)).collect(),
        };
        let vc = ValleyConstruction::new(F121_AGE_K, None);
        let off = skeleton(&f, &vc, &ss, 51.2);
        let on = skeleton(&f, &ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..vc }, &ss, 51.2);
        let (up, down) = (20 * n + 64, 100 * n + 64);
        assert!(
            (off.base_alt_m[up] - vc.base_m).abs() < 1e-3,
            "negative control: without the field the upstream cell takes the sea's base ({} m)",
            off.base_alt_m[up]
        );
        assert!(on.base_alt_m[up] > 150.0, "with the field its base is the lake's surface, not the sea ({} m)", on.base_alt_m[up]);
        assert!(on.chi_m[up] < off.chi_m[up], "its χ counts from the lake only ({} vs {})", on.chi_m[up], off.chi_m[up]);
        assert_eq!(on.base_alt_m[down].to_bits(), off.base_alt_m[down].to_bits(), "downstream of the lake nothing moves");
        assert_eq!(on.chi_m[down].to_bits(), off.chi_m[down].to_bits(), "downstream of the lake nothing moves");
    }

    /// ADR Finding 140, rule 13 — the gorge's retreat. A plane falls 4 m per cell to the sea on 128² cells of 400 m;
    /// a bowl 150 m deep sits in it (a present lake). Negative controls FIRST:
    /// (a) without the gorge, the extended lake base still lets the outlet valley, laid from the sea, cut the rim
    /// more than 10 m below the lake's spill (F136-K7's notch);
    /// (b) the gate without the extended lake base is inert, bit for bit.
    /// With the gorge at r = 0 the rim keeps ≥ the spill (the lake full); at r = 2 the lake's level is its floor
    /// and the col is cut down to it (the lake drained). The gate is absent from the serialised key when off.
    #[test]
    fn the_gorge_keeps_the_rim_at_the_lakes_level_and_retreats_with_age() {
        let ss = SteinSteinParams::default();
        let n = 128usize;
        let z_m = |x: f32, y: f32| -> f32 {
            if y >= 120.0 {
                return -50.0;
            }
            let d = ((x - 64.0).powi(2) + (y - 50.0).powi(2)).sqrt();
            let bowl = if d < 12.0 { 150.0 * (1.0 - (d / 12.0).powi(2)) } else { 0.0 };
            5.0 + 4.0 * (120.0 - y) - bowl
        };
        let f = GridF32 {
            width: n,
            height: n,
            data: (0..n * n).map(|k| c1_metres_to_altitude_norm(z_m((k % n) as f32 + 0.5, (k / n) as f32 + 0.5), &ss)).collect(),
        };
        let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
        let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, None) };
        let gorge = |r: f32| ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(r, 0.0)), ..ext };
        // the key: absent when off
        let js = serde_json::to_string(&ext).unwrap();
        assert!(!js.contains("gorge_retreat"), "{js}");
        assert!(serde_json::to_string(&gorge(1.0)).unwrap().contains("gorge_retreat"), "negative control: present when on");
        // (b) the gate needs the extended lake base: without it, bit-identical
        let plain = ValleyConstruction::new(F121_AGE_K, None);
        let (a, _) = carve(&f, &skeleton(&f, &plain, &ss, 51.2), &plain, &ss);
        let with = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.0)), ..plain };
        let (b, _) = carve(&f, &skeleton(&f, &with, &ss, 51.2), &with, &ss);
        assert!(a.data.iter().zip(b.data.iter()).all(|(x, y)| x.to_bits() == y.to_bits()), "the gate alone must be inert");
        // the body and its ring, from the gorge's own diagnostics
        let sk0 = skeleton(&f, &gorge(0.0), &ss, 51.2);
        assert_eq!(sk0.gorge_bodies.len(), 1, "one present lake ≥ 1 km²");
        let body = sk0.gorge_bodies[0].clone();
        assert!(body.col != u32::MAX, "the lake has an outflow");
        let bo = sk0.gorge_body_of.as_ref().expect("diagnostics");
        let ring: Vec<usize> = (0..n * n)
            .filter(|&k| {
                bo[k] == u32::MAX
                    && f.data[k] > C1_SEA_LEVEL_NORM
                    && [(-1i32, -1i32), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)].iter().any(|&(dx, dy)| {
                        let (x, y) = ((k % n) as i32 + dx, (k / n) as i32 + dy);
                        x >= 0 && y >= 0 && x < n as i32 && y < n as i32 && bo[y as usize * n + x as usize] == 0
                    })
            })
            .collect();
        let ring_min = |g: &GridF32| ring.iter().map(|&k| m(g, k)).fold(f32::INFINITY, f32::min);
        // (a) the negative control: without the gorge the rim is notched
        let (no_gorge, _) = carve(&f, &skeleton(&f, &ext, &ss, 51.2), &ext, &ss);
        assert!(
            ring_min(&no_gorge) < body.l_in - 10.0,
            "negative control: without the gorge the outlet cuts the rim ({} m against the spill {} m)",
            ring_min(&no_gorge),
            body.l_in
        );
        // r = 0: the lake full, the rim kept
        assert!((body.level - body.l_in).abs() < 1e-3, "r = 0: L = L_in");
        let (young, _) = carve(&f, &sk0, &gorge(0.0), &ss);
        assert!(ring_min(&young) >= body.l_in - 0.5, "r = 0: the rim keeps ≥ the spill ({} m against {} m)", ring_min(&young), body.l_in);
        // r = 2: drained, the col cut to the floor
        let sk2 = skeleton(&f, &gorge(2.0), &ss, 51.2);
        let b2 = &sk2.gorge_bodies[0];
        assert!((b2.level - b2.l_floor).abs() < 1e-3, "r = 2: L = L_floor");
        let (old, _) = carve(&f, &sk2, &gorge(2.0), &ss);
        assert!(m(&old, b2.col as usize) <= b2.l_floor + 1.0, "r = 2: the col is cut to the floor ({} m against {} m)", m(&old, b2.col as usize), b2.l_floor);
        assert!(ring_min(&old) >= b2.l_floor - 0.5, "r = 2: no ring cell under the drained level");
    }

    /// ADR Finding 142, rule 13 — spec v4's minimal plain. The bowl of the gorge test (150 m deep, a present lake)
    /// carries a second, separate hollow inside its footprint (a 40 m dimple 8 cells from the centre). Drained (r = 2):
    /// - **negative control first**: without the plain the dimple stays a closed hollow outside the lake (a dry run of
    ///   the fill would raise cells);
    /// - with the plain: the fill raises cells, only up (never down), only inside the footprint, by the volume it
    ///   reports; run again, it raises nothing (no closed hollow left outside the lake).
    #[test]
    fn the_minimal_plain_fills_a_drained_bowls_residual_hollows_and_only_raises() {
        let ss = SteinSteinParams::default();
        let n = 128usize;
        let z_m = |x: f32, y: f32| -> f32 {
            if y >= 120.0 {
                return -50.0;
            }
            let d = ((x - 64.0).powi(2) + (y - 50.0).powi(2)).sqrt();
            let bowl = if d < 12.0 { 150.0 * (1.0 - (d / 12.0).powi(2)) } else { 0.0 };
            let e = ((x - 56.0).powi(2) + (y - 50.0).powi(2)).sqrt();
            let dimple = if e < 2.5 { 40.0 * (1.0 - (e / 2.5).powi(2)) } else { 0.0 };
            5.0 + 4.0 * (120.0 - y) - bowl - dimple
        };
        let f = GridF32 {
            width: n,
            height: n,
            data: (0..n * n).map(|k| c1_metres_to_altitude_norm(z_m((k % n) as f32 + 0.5, (k / n) as f32 + 0.5), &ss)).collect(),
        };
        let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, None) };
        let vc = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(2.0, 0.0)), ..ext };
        let sk = skeleton(&f, &vc, &ss, 51.2);
        assert_eq!(sk.gorge_bodies.len(), 1, "the bowl is in v4's scope (an input lake standing alone)");
        assert!(sk.gorge_bodies[0].r_lake > 1.0, "drained");
        let (carved, _) = carve(&f, &sk, &vc, &ss);
        // negative control: the dimple is a residual hollow after the construction alone
        let mut dry = carved.clone();
        let rec0 = gorge_plain(&mut dry, &sk, &ss, 1.0);
        assert_eq!(rec0.len(), 1);
        assert!(rec0[0].cells > 0 && rec0[0].max_depth_m > 1.0, "negative control: a hollow is left without the plain ({:?})", rec0[0]);
        // the plain
        let mut out = carved.clone();
        let rec = gorge_plain(&mut out, &sk, &ss, 1.0);
        let bo = sk.gorge_body_of.as_ref().expect("diagnostics");
        let cell_km2 = (sk.cell_m / 1000.0).powi(2) as f64;
        let mut vol = 0f64;
        for k in 0..n * n {
            assert!(out.data[k] >= carved.data[k], "the plain only raises (cell {k})");
            if out.data[k] != carved.data[k] {
                assert_eq!(bo[k], 0, "the plain acts only inside the footprint (cell {k})");
                vol += ((c1_altitude_norm_to_metres(out.data[k], &ss) - c1_altitude_norm_to_metres(carved.data[k], &ss)) as f64) * cell_km2 * 1e-3;
            }
        }
        assert!((vol - rec[0].volume_km3).abs() <= 1e-6 + 1e-3 * vol, "the reported volume is the raise ({vol} vs {})", rec[0].volume_km3);
        let again = gorge_plain(&mut out.clone(), &sk, &ss, 1.0);
        assert_eq!(again[0].cells, 0, "after the plain no closed hollow is left outside the lake");
        // below the threshold (r ≤ 1) the plain does nothing
        let young = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(1.0, 0.0)), ..ext };
        let sky = skeleton(&f, &young, &ss, 51.2);
        assert!(gorge_plain(&mut carved.clone(), &sky, &ss, 1.0).is_empty(), "r_lake ≤ 1: no plain");
    }

    /// ADR Finding 142, rule 13 — spec v4's scope. The plane of the gorge test carries an above-sea bowl standing
    /// alone, and a shallow bowl opening onto the fringe of a below-sea basin (a closed depression of the open-ocean
    /// flood). Negative control first: v3 treats a body that is not an input lake standing alone. v4 keeps only the
    /// lone bowl; a dropped body's cells keep the extended lake base, bit for bit as without the gorge.
    #[test]
    fn v4_keeps_only_the_input_lakes_standing_alone() {
        let ss = SteinSteinParams::default();
        let n = 128usize;
        let z_m = |x: f32, y: f32| -> f32 {
            if y >= 120.0 {
                return -50.0;
            }
            let bowl = |cx: f32, cy: f32, r: f32, depth: f32| {
                let d = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
                if d < r { depth * (1.0 - (d / r).powi(2)) } else { 0.0 }
            };
            5.0 + 4.0 * (120.0 - y) - bowl(32.0, 40.0, 10.0, 120.0) - bowl(96.0, 60.0, 14.0, 500.0) - bowl(96.0, 46.0, 7.0, 60.0)
        };
        let f = GridF32 {
            width: n,
            height: n,
            data: (0..n * n).map(|k| c1_metres_to_altitude_norm(z_m((k % n) as f32 + 0.5, (k / n) as f32 + 0.5), &ss)).collect(),
        };
        let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, None) };
        let v3 = skeleton(&f, &ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.0)), ..ext }, &ss, 51.2);
        let v4 = skeleton(&f, &ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(1.0, 0.0)), ..ext }, &ss, 51.2);
        let lone = 40 * n + 32;
        let lone_ok = |b: &GorgeBody| b.input_lake && !b.touches_depression && !b.touches_below_sea;
        assert!(v3.gorge_bodies.iter().any(|b| !lone_ok(b)), "negative control: v3 treats a body that is not a lone input lake ({:?})", v3.gorge_bodies);
        assert!(v3.gorge_body_of.as_ref().unwrap()[lone] != u32::MAX, "v3 treats the lone bowl");
        let bo4 = v4.gorge_body_of.as_ref().unwrap();
        assert!(bo4[lone] != u32::MAX, "v4 keeps the lone bowl");
        assert!(v4.gorge_bodies.iter().all(lone_ok), "v4 keeps only lone input lakes");
        assert!(v4.gorge_bodies.len() < v3.gorge_bodies.len(), "v4 drops the merged bodies");
        // a dropped body's cells keep the extended lake base (as without the gorge)
        let none = skeleton(&f, &ext, &ss, 51.2);
        let bo3 = v3.gorge_body_of.as_ref().unwrap();
        for k in 0..n * n {
            if bo3[k] != u32::MAX && bo4[k] == u32::MAX {
                assert_eq!(v4.base_alt_m[k].to_bits(), none.base_alt_m[k].to_bits(), "a dropped body's cell keeps the extended base (cell {k})");
            }
        }
        // the key: v4 is not v3
        let js3 = serde_json::to_string(&ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.5)), ..ext }).unwrap();
        let js4 = serde_json::to_string(&ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(1.0, 0.5)), ..ext }).unwrap();
        assert_ne!(js3, js4, "v4 must move the key");
    }

    /// ADR Finding 143, rule 13 — the exact scope. Two above-sea bowls stand alone on the plane of the gorge test.
    /// Negative control first: v4's proxy keeps both. With `scope_lows` listing only the first bowl's lowest cell,
    /// only it is kept; the other bowl's cells keep the extended lake base bit for bit (as without the gorge). The
    /// list moves the key, and v4 without it serialises as before (no `scope_lows`). The bench options φ = 0 and the
    /// frozen design move the key too and are absent when off.
    #[test]
    fn the_exact_scope_keeps_only_the_listed_bodies() {
        let ss = SteinSteinParams::default();
        let n = 128usize;
        let z_m = |x: f32, y: f32| -> f32 {
            if y >= 120.0 {
                return -50.0;
            }
            let bowl = |cx: f32, cy: f32, r: f32, depth: f32| {
                let d = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
                if d < r { depth * (1.0 - (d / r).powi(2)) } else { 0.0 }
            };
            5.0 + 4.0 * (120.0 - y) - bowl(32.0, 40.0, 10.0, 120.0) - bowl(96.0, 40.0, 10.0, 120.0)
        };
        let f = GridF32 {
            width: n,
            height: n,
            data: (0..n * n).map(|k| c1_metres_to_altitude_norm(z_m((k % n) as f32 + 0.5, (k / n) as f32 + 0.5), &ss)).collect(),
        };
        let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, None) };
        let proxy = skeleton(&f, &ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(1.0, 0.0)), ..ext }, &ss, 51.2);
        assert_eq!(proxy.gorge_bodies.len(), 2, "negative control: the proxy keeps both lone bowls ({:?})", proxy.gorge_bodies);
        let (a, b) = (40 * n + 32, 40 * n + 96);
        let bo_p = proxy.gorge_body_of.as_ref().unwrap();
        let low_a = proxy.gorge_bodies[bo_p[a] as usize].low;
        let mut list = [u32::MAX; 16];
        list[0] = low_a;
        let exact = GorgeRetreat { scope_lows: Some(list), ..GorgeRetreat::v4(1.0, 0.0) };
        let sk = skeleton(&f, &ValleyConstruction { gorge_retreat: Some(exact), ..ext }, &ss, 51.2);
        assert_eq!(sk.gorge_bodies.len(), 1, "the exact scope keeps the listed body only");
        assert_eq!(sk.gorge_bodies[0].low, low_a);
        let bo = sk.gorge_body_of.as_ref().unwrap();
        assert!(bo[a] != u32::MAX && bo[b] == u32::MAX, "bowl A kept, bowl B dropped");
        let none = skeleton(&f, &ext, &ss, 51.2);
        for k in 0..n * n {
            if bo_p[k] != u32::MAX && bo[k] == u32::MAX {
                assert_eq!(sk.base_alt_m[k].to_bits(), none.base_alt_m[k].to_bits(), "a dropped body's cell keeps the extended base (cell {k})");
            }
        }
        // the key
        let js = |g: GorgeRetreat| serde_json::to_string(&ValleyConstruction { gorge_retreat: Some(g), ..ext }).unwrap();
        let v4 = js(GorgeRetreat::v4(1.0, 0.0));
        assert!(!v4.contains("scope_lows") && !v4.contains("phi_zero") && !v4.contains("freeze_design"), "v4 serialises as before: {v4}");
        assert_ne!(js(exact), v4, "the list moves the key");
        assert_ne!(js(GorgeRetreat { phi_zero: true, ..GorgeRetreat::v4(1.0, 0.0) }), v4, "φ = 0 moves the key");
        assert_ne!(js(GorgeRetreat { freeze_design: true, ..GorgeRetreat::v4(1.0, 0.0) }), v4, "the frozen design moves the key");
        // φ = 0: no head fall's share
        let z = skeleton(&f, &ValleyConstruction { gorge_retreat: Some(GorgeRetreat { phi_zero: true, ..GorgeRetreat::v4(1.0, 0.0) }), ..ext }, &ss, 51.2);
        assert!(z.gorge_bodies.iter().all(|b| b.phi == 0.0), "φ = 0 on every body");
    }

    /// ADR Finding 133-F, rule 13 — a below-sea basin's lake is a present lake: χ stops at its shore. A plane
    /// falls 4 m per cell to the sea on 128² cells of 400 m; an inland bowl whose floor lies 100 m BELOW the sea
    /// is set in it (a closed depression of the open-ocean flood). Negative control FIRST: under Finding 132's
    /// `InputLakes`, a cell upstream integrates χ across the bowl's slope down to the water (the col walk), so a
    /// cell on that slope carries χ > 0. Extended: every land cell of the depression has χ = 0 and the col as its
    /// base, the upstream cell's χ is smaller, and a cell outside the bowl's catchment is bit-identical.
    #[test]
    fn a_below_sea_basin_lake_stops_chi_at_its_shore() {
        let ss = SteinSteinParams::default();
        let n = 128usize;
        let z_m = |x: f32, y: f32| -> f32 {
            if y >= 120.0 {
                return -50.0;
            }
            let d = ((x - 64.0).powi(2) + (y - 50.0).powi(2)).sqrt();
            let bowl = if d < 20.0 { 450.0 * (1.0 - (d / 20.0).powi(2)) } else { 0.0 };
            5.0 + 4.0 * (120.0 - y) - bowl
        };
        let f = GridF32 {
            width: n,
            height: n,
            data: (0..n * n).map(|k| c1_metres_to_altitude_norm(z_m((k % n) as f32 + 0.5, (k / n) as f32 + 0.5), &ss)).collect(),
        };
        let vc = ValleyConstruction::new(F121_AGE_K, None);
        let f132 = skeleton(&f, &ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..vc }, &ss, 51.2);
        let ext = skeleton(&f, &ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..vc }, &ss, 51.2);
        // a land cell on the bowl's slope, below the col: the depression
        let slope = (50 * n) + 64 - 14 * n; // (64, 36)
        assert!(z_m(64.5, 36.5) > 0.0, "the slope cell is land");
        assert!(f132.chi_m[slope] > 0.0, "negative control: under InputLakes the col walk gives the slope χ > 0 ({})", f132.chi_m[slope]);
        assert_eq!(ext.chi_m[slope], 0.0, "extended: the basin lake's land cell has χ = 0");
        let up = 5 * n + 64;
        assert!(ext.chi_m[up] < f132.chi_m[up], "the upstream cell counts χ from the shore ({} vs {})", ext.chi_m[up], f132.chi_m[up]);
        assert!((ext.base_alt_m[up] - f132.base_alt_m[up]).abs() < 1e-3, "both take the col as base");
        let away = 110 * n + 10;
        assert_eq!(ext.chi_m[away].to_bits(), f132.chi_m[away].to_bits(), "outside the bowl's catchment nothing moves");
    }

    /// ADR Finding 128-C, rule 13 — D8-LTD follows a planar slope off the lattice. On a plane whose
    /// slope points 14.04° off north (a 1:4 lateral:longitudinal gradient, Paik 2008's plane) every D8
    /// path runs due north, 14° off (the negative control). LTD paths' start-to-end direction stays
    /// within 1.5° of the slope, and each cell has exactly one pointer.
    #[test]
    fn the_ltd_tree_follows_a_planar_slope_off_the_lattice() {
        let n = 80usize;
        let th = (0.25f32).atan(); // 14.04° off the north axis
        let zm: Vec<f32> = (0..n * n)
            .map(|k| {
                let (x, y) = ((k % n) as f32, (k / n) as f32);
                // descending toward north (−y) and slightly toward +x
                1000.0 + (y * th.cos() - x * th.sin()) * 10.0
            })
            .collect();
        let sink: Vec<bool> = (0..n * n).map(|k| {
            let (x, y) = (k % n, k / n);
            x == 0 || y == 0 || x == n - 1 || y == n - 1
        }).collect();
        // the D8 fallback / reference: the steepest of the eight neighbours
        let d8: Vec<u8> = (0..n * n)
            .map(|c| {
                if sink[c] {
                    return DIR_NONE;
                }
                let (x, y) = ((c % n) as i32, (c / n) as i32);
                let mut best = (DIR_NONE, 0f32);
                for k in 0..8 {
                    let m = (y + D8_DY[k]) as usize * n + (x + D8_DX[k]) as usize;
                    let s = (zm[c] - zm[m]) / if k % 2 == 1 { 1.414 } else { 1.0 };
                    if s > best.1 {
                        best = (k as u8, s);
                    }
                }
                best.0
            })
            .collect();
        let (ltd, flats) = ltd_directions(&zm, &sink, &d8, n, n, 10.0);
        assert_eq!(flats, 0, "a plane has no flat");
        let heading = |dir: &[u8], c0: usize| -> Option<f32> {
            let mut c = c0;
            let mut steps = 0;
            while !sink[c] && steps < 400 {
                let k = dir[c] as usize;
                c = ((c / n) as i32 + D8_DY[k]) as usize * n + ((c % n) as i32 + D8_DX[k]) as usize;
                steps += 1;
            }
            let (dx, dy) = ((c % n) as f32 - (c0 % n) as f32, (c0 / n) as f32 - (c / n) as f32);
            (steps >= 40).then(|| dx.atan2(dy).to_degrees())
        };
        let mean_dev = |dir: &[u8]| {
            let v: Vec<f32> = (0..n * n)
                .filter(|&c| c / n > 60 && !sink[c])
                .filter_map(|c| heading(dir, c))
                .map(|a| (a - th.to_degrees()).abs())
                .collect();
            assert!(v.len() > 100, "the paths must be long ({})", v.len());
            v.iter().sum::<f32>() / v.len() as f32
        };
        let (dev8, devl) = (mean_dev(&d8), mean_dev(&ltd));
        assert!(dev8 > 10.0, "negative control: D8 must run on its lattice ({dev8:.2}°)");
        assert!(devl < 1.5, "D8-LTD must follow the slope ({devl:.2}°)");
        let acc = accumulate_cells(&ltd, &sink, n, n);
        let land_total = (0..n * n).filter(|&c| !sink[c]).count() as f32;
        let into_sinks: f32 = (0..n * n)
            .filter(|&c| !sink[c])
            .filter(|&c| {
                let k = ltd[c] as usize;
                sink[((c / n) as i32 + D8_DY[k]) as usize * n + ((c % n) as i32 + D8_DX[k]) as usize]
            })
            .map(|c| acc[c])
            .sum();
        assert_eq!(into_sinks, land_total, "every land cell drains exactly once (one pointer each)");
    }

    /// ADR Finding 122-B, rule 13 — with `basin_base`, no constructed floor upstream of an inland
    /// below-sea basin lies below that basin's spill level; without it (Finding 121), the valleys
    /// are built down toward sea + 0.5 m inside the basin's rim (the negative control).
    #[test]
    fn basin_base_keeps_the_floors_above_the_water_of_a_closed_basin() {
        let ss = SteinSteinParams::default();
        let mut f = basin(&ss);
        let w = f.width;
        // an inland basin on the trunk: a paraboloid dipping ~100 m below the sea, radius 8 cells
        for y in 0..f.height {
            for x in 0..w {
                let d2 = ((x as f32 - 48.0).powi(2) + (y as f32 - 50.0).powi(2)) / 64.0;
                if d2 < 1.0 {
                    let k = y * w + x;
                    let m = c1_altitude_norm_to_metres(f.data[k], &ss) - 400.0 * (1.0 - d2);
                    f.data[k] = c1_metres_to_altitude_norm(m, &ss);
                }
            }
        }
        let run = |basin: bool| {
            let mut vc = ValleyConstruction::f121(0.02, None);
            vc.a_min_km2 = 0.5;
            vc.smooth_m = 0.0;
            vc.basin_base = basin;
            let sk = skeleton(&f, &vc, &ss, 4.8);
            carve(&f, &sk, &vc, &ss)
        };
        // the lowest carved altitude upstream of the basin (beyond its rim)
        let low_up = |out: &GridF32, m: &CarveMasks| -> f32 {
            (60 * w..f.data.len())
                .filter(|&k| m.carved[k])
                .map(|k| c1_altitude_norm_to_metres(out.data[k], &ss))
                .fold(f32::INFINITY, f32::min)
        };
        let (o_off, m_off) = run(false);
        let (o_on, m_on) = run(true);
        let (lo_off, lo_on) = (low_up(&o_off, &m_off), low_up(&o_on, &m_on));
        // the basin's rim on the trunk: 40 + 6·(42 − 6) m ≈ 256 m, its lowest col
        assert!(
            lo_off < 150.0,
            "negative control: F121 must build below the basin's col, got {lo_off}"
        );
        assert!(lo_on > 240.0, "basin_base must keep the floors above the col, got {lo_on}");
        for k in 0..f.data.len() {
            assert!(o_on.data[k] <= f.data[k], "cell {k} RAISED with basin_base");
        }
    }

    /// ADR Finding 121, rule 13 — the construction only LOWERS, lays the floor the law states on
    /// the trunk, and leaves the field untouched beyond the walls.
    #[test]
    fn the_valley_construction_lowers_only_and_lays_its_floor() {
        let ss = SteinSteinParams::default();
        let f = basin(&ss);
        // 96 cells over 4.8 km: 50 m cells, 0.0025 km² per cell
        let domain_km = 4.8;
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.a_min_km2 = 0.5; // 200 cells: the centre column carries far more
        vc.smooth_m = 0.0; // a straight trunk: no smoothing needed, exact geometry
        let sk = skeleton(&f, &vc, &ss, domain_km);
        let n_trunk = sk.trunk.iter().filter(|&&b| b).count();
        assert!(n_trunk > 20, "the basin must carry a trunk, got {n_trunk} cells");
        let (out, masks) = carve(&f, &sk, &vc, &ss);
        // 1. never raises
        for k in 0..f.data.len() {
            assert!(out.data[k] <= f.data[k], "cell {k} was RAISED");
        }
        assert!(masks.carved.iter().any(|&b| b), "negative control: the construction must bite");
        // 2. on the trunk, where the law's floor is below the field, the cell IS the floor
        let mut checked = 0;
        for k in 0..f.data.len() {
            if sk.trunk[k] {
                let fl = sk.floor_m(k, vc.age_k);
                let z = c1_altitude_norm_to_metres(f.data[k], &ss);
                if fl < z - 1.0 {
                    let got = c1_altitude_norm_to_metres(out.data[k], &ss);
                    assert!(
                        (got - fl).abs() < 0.05,
                        "trunk cell {k}: {got} m against floor {fl} m"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 10, "the floor clause must be exercised, got {checked}");
        // 3. far from the trunk (the basin's edges), the field is untouched
        let w = f.width;
        for y in 10..90 {
            for x in [0usize, 1, 94, 95] {
                assert_eq!(out.data[y * w + x], f.data[y * w + x], "({x},{y}) was touched");
            }
        }
        // 4. the walls: no carved cell is steeper than the wall law allows toward the floor
        let tan = vc.wall_deg.to_radians().tan();
        for y in 20..80 {
            for x in 1..w - 1 {
                let k = y * w + x;
                if masks.carved[k] && !masks.floor[k] {
                    let (a, b) = (
                        c1_altitude_norm_to_metres(out.data[k - 1], &ss),
                        c1_altitude_norm_to_metres(out.data[k + 1], &ss),
                    );
                    let g = (a - b).abs() / (2.0 * sk.cell_m);
                    assert!(g <= tan * 1.05 + 1e-3, "({x},{y}) wall gradient {g} > tan 28°");
                }
            }
        }
    }
}
