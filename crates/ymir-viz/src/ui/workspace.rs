//! HD workspace — the egui interface refined to the Claude Design mock (UI
//! rewrite step d2: polish). Builds on d1's functional skeleton (generate →
//! view HD layers → inspect, non-blocking) and matches the mock's layout:
//! a top bar (brand + Standard/Expert), a left control panel (progressive
//! disclosure), a central map with a pipeline frieze + contextual legend, and
//! a right inspector with grouped sections. Copper is the CHROME; the data
//! palettes stay canonical.
//!
//! Scope notes (honest divergences from the mock):
//! - The frieze shows the FIVE layers the HD product (step b) actually carries
//!   (Relief / Drainage / Précip. / Température / Biomes). The mock's separate
//!   Tectonique + Érosion nodes need the coarse snapshot / a pre-erosion field
//!   that the HD flow does not surface — deferred.
//! - Expert parameter sliders are EXPOSED (progressive disclosure) but not yet
//!   wired to the engine (only seed / resolution / latitude drive generation);
//!   wiring the rest is a follow-up. They are tagged as such.
//! - Live per-node frieze animation during compute is step e (here the frieze
//!   is the selector + post-generation state).
//! - Fonts approximate the mock (egui defaults; Space Grotesk / IBM Plex are
//!   not bundled).

use std::path::PathBuf;
use std::sync::Arc;

use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};
use egui::Color32 as C;

use crate::bridge::c1::{
    C1RunSpec, C1SolverBridge, CacheRegime, CellInspection, HdParams, HdPhase, HdResult, HdState,
    PreviewShape, RiverCellMap, inspect_cell,
};
use crate::visualization::colormap::hypsometric_bipolar;
use ymir_core::climate::biomes::Biome;
use ymir_core::climate::precipitation::{precip_mm_per_year, wind_zonal_dir};
use ymir_core::climate::temperature::sea_level_temperature;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{LakeType, Navigability, SegmentKind, UnresolvedReason};
use ymir_core::tectonics_c1::land_topology::domain_metrics;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

// ── Palette (the mock's exact hex) ───────────────────────────────────────
const COPPER: C = C::from_rgb(0xB8, 0x73, 0x33);
const COPPER_BRIGHT: C = C::from_rgb(0xC9, 0x85, 0x3F);
const BRONZE: C = C::from_rgb(0x9A, 0x73, 0x50);
const INK: C = C::from_rgb(0x1A, 0x12, 0x0A); // text on copper
const PANEL: C = C::from_rgb(0x1B, 0x1B, 0x1B);
const PANEL2: C = C::from_rgb(0x16, 0x16, 0x16);
const VIEWPORT: C = C::from_rgb(0x0A, 0x0A, 0x0A);
const FIELD: C = C::from_rgb(0x11, 0x11, 0x11);
const TEXT: C = C::from_rgb(0xE0, 0xE0, 0xE0);
const TEXT_BRIGHT: C = C::from_rgb(0xEA, 0xEA, 0xEA);
const DIM: C = C::from_rgb(0x8A, 0x8A, 0x8A);
const DIM2: C = C::from_rgb(0x6E, 0x6E, 0x6E);
const GREEN: C = C::from_rgb(0x5F, 0xA5, 0x5F);

pub struct HdWorkspacePlugin;

impl Plugin for HdWorkspacePlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<WorkspaceState>();
        app.add_systems(EguiPrimaryContextPass, draw_workspace);
    }
}

/// Which sink a watercourse terminates in (Finding 28 — inspection microscope).
#[derive(Clone, Copy, PartialEq)]
enum Sink {
    Sea,
    ExoLake,
    EndoLake,
    /// A below-sea ENCLOSED cell (water_class == 2) that is NOT an inventoried lake. Distinct
    /// from Sea: tagging it "→ mer" would contradict water_class (the authority) and reinstate the
    /// very altitude proxy Finding 31 removed.
    ///
    /// ⚠️ **ADR Finding 91-D: this was labelled "évaporatif" and the condition cannot support the
    /// word.** It reads `water_class` alone — never `net_evap`, never `a_eq`, never `lake_type` —
    /// so it says "enclosed below sea and in no inventory", which is Finding 89's `to_nothing`,
    /// not a regime. Finding 92-C removed the population that made it common; what can still
    /// reach it is a mouth on a cell no body covers.
    SubSeaSink,
    Unknown,
}

/// An aggregated WATERCOURSE (a trunk + its tributaries), assembled UI-side from the
/// exported per-segment drainage — NOT recomputed. `segments` are every segment of one river
/// SYSTEM (ADR Finding 124-3: each reach in exactly one); `trunk` is the source→sink main stem
/// (the max-area climb, Finding 123); `catchment_km2` is its share of the geometric partition
/// of the land (an AREA, signified). Everything read from `C1DrainageResult`.
#[derive(Clone)]
struct Watercourse {
    segments: Vec<usize>,
    trunk: Vec<usize>, // source → sink
    order: u8,
    discharge_m3s: f32,
    width_mouth_m: f32,
    width_source_m: f32,
    catchment_km2: f32,
    mouth_xy: (u32, u32),
    sink: Sink,
    /// If the sink is a lake, its `lake_map` id (for the clickable chain).
    sink_lake_id: Option<u32>,
    length_km: f32,
    tributaries: usize,
    /// What the trunk mouth IS (`segment_kind`). `Spillway` = the outflow of a closed
    /// below-sea basin over its col — a real flow with NO hierarchy, so `order` is
    /// meaningless for it and the list must not present it as a river.
    kind: SegmentKind,
    /// For a `Spillway`, the below-sea basin it drains — `None` when that basin sits
    /// below the lake-inventory floor (absent from `drainage.lakes`).
    source_lake_id: Option<u32>,
}

/// A navigation jump requested by a detail panel (the clickable hydrological chain).
#[derive(Clone, Copy)]
enum NavAction {
    Lake(usize),
    River(usize),
}

/// Which entity list the inspection panel is showing.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum InspectTab {
    #[default]
    Rivers,
    Lakes,
    Volcanoes,
}

/// ADR Finding 141 -- the gorge's retreat toggle (Finding 140) stays hidden until its crash is fixed.
const GORGE_TOGGLE_VISIBLE: bool = false;

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
enum HdLayer {
    #[default]
    Relief,
    Drainage,
    Precipitation,
    Temperature,
    Biomes,
    /// ADR Finding 152/153 -- the rock classes (a frieze step, like the biomes). The favourabilities are a MASK
    /// (`WorkspaceState::chance_overlay`), drawn over any step.
    Rocks,
}

impl HdLayer {
    fn label(self) -> &'static str {
        match self {
            HdLayer::Relief => "Relief",
            HdLayer::Drainage => "Drainage",
            HdLayer::Precipitation => "Précipitation",
            HdLayer::Temperature => "Température",
            HdLayer::Biomes => "Biomes",
            HdLayer::Rocks => "Roches",
        }
    }
    fn desc(self) -> &'static str {
        match self {
            HdLayer::Relief => "Élévation hypsométrique (relief + bathymétrie), après érosion.",
            HdLayer::Drainage => "Réseau hydrographique — navigabilité des rivières et lacs.",
            HdLayer::Precipitation => {
                "Précipitation annuelle — bandes ITCZ / subtropicale / mid-latitude."
            }
            HdLayer::Temperature => "Température de surface — gradient latitudinal + lapse rate.",
            HdLayer::Biomes => "Classification de Whittaker (température × précipitation).",
            HdLayer::Rocks => "Classes de roche (géologie v1) : fixes pour un monde.",
        }
    }
}

/// The frieze nodes (suite e): `(worker phase, short label, viewable layer)`.
/// The layer is `Some` only when the HD product RETAINS a displayable layer for
/// that node → clickable/selectable. `None` = a pure build-progress stage (no
/// retained layer to view): Tectonique (coarse S̃/plate/craton, not in the HD
/// product) and Relief (pre-erosion upscale, not retained). The Érosion node
/// carries the Relief layer (the eroded heightmap = the final relief). Précip.
/// and Temp. share the Climate phase.
const FRIEZE: [(HdPhase, &str, Option<HdLayer>); 8] = [
    (HdPhase::Tectonic, "Tectonique", None),
    (HdPhase::Relief, "Relief", None),
    (HdPhase::Erosion, "Érosion", Some(HdLayer::Relief)),
    (HdPhase::Drainage, "Drainage", Some(HdLayer::Drainage)),
    (HdPhase::Climate, "Précip.", Some(HdLayer::Precipitation)),
    (HdPhase::Climate, "Temp.", Some(HdLayer::Temperature)),
    (HdPhase::Biomes, "Biomes", Some(HdLayer::Biomes)),
    // ADR Finding 153 -- the rocks are a step of the pipeline, a view like the biomes
    (HdPhase::Geology, "Roches", Some(HdLayer::Rocks)),
];

/// Visual state of a frieze node, derived from the live HD event stream.
#[derive(Clone, Copy, PartialEq)]
enum NodeVis {
    Pending,
    Running,
    Done { cached: bool },
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Standard,
    Expert,
}

#[derive(Resource)]
struct WorkspaceState {
    seed: u64,
    resolution: usize,
    latitude: f32,
    mode: Mode,
    layer: HdLayer,
    /// River-network overlay on the map (any layer): polylines, thickness by
    /// Strahler order, colour by navigability — ORPHAN reaches (no downstream,
    /// not ending at a sink) drawn RED so they stand out for validation.
    river_overlay: bool,
    /// ADR Finding 146 -- the "Rivières LL" layer: the smoothed main stems of `rivers_ll.json`, filtered by length
    /// (km), catchment (km²) and Strahler order (the selection is the author's; nothing is fixed).
    rll_overlay: bool,
    rll_min_len_km: f32,
    rll_min_area_km2: f32,
    rll_min_strahler: u8,
    /// ADR Findings 147-F / 149-T -- « Masquer faisceaux et micro-rivières » hides the fusion candidates (the
    /// lower-discharge river of a parallel pair) and the micro-rivers into a lake (< 1 km); a separate box hides the
    /// spillways that retrace watercourses. Display only.
    rll_hide_fusion: bool,
    rll_hide_retrace: bool,
    /// ADR Finding 152 -- the geology rules file, the resource the « Chances » view shows (`None`: all, the highest
    /// favourability wins), and a re-zoning done by « Recharger les règles » (zoning, rules, seconds).
    geology_rules_path: String,
    chance_only: Option<usize>,
    /// ADR Finding 153 -- the « Chances » mask (over any step).
    chance_overlay: bool,
    zoning_override: Option<(ymir_core::geology::Zoning, ymir_core::geology::rules::LoadedRules, f64)>,
    inspector_open: bool,
    // Expert params (exposed, wiring deferred — tagged in the UI).
    climat_open: bool,
    relief_open: bool,
    drainage_open: bool,
    craton_density: f32,
    shield_frac: f32,
    shelf_width: f32,
    erosion: f32,
    channel_jitter: f32,
    dinf: bool,
    /// Physical size (km) the tectonic domain represents — the domain IS the map
    /// (no crop). Pure relabelling of the 64² pattern: `m_per_px = domain_km ·
    /// 1000 / resolution`. Moving it recomputes NO tectonics (see `m_per_px`).
    domain_km: f32,
    /// Finding 24 — geographic scale ratio: the map DRAWS `domain_km` but the exported
    /// hydrology SIGNIFIES `domain_km · ratio`. Pure presentation multiplier on the
    /// river width / discharge / navigability; NOTHING physical (terrain, climate,
    /// biomes) sees it. Instant (post-process, no re-solve). 1.0 = identity.
    geo_scale_ratio: f32,
    /// Finding 25 — CLIMATIC latitude span (°) centred on `latitude`, INDEPENDENT of
    /// `domain_km`. Widen it to cross several belts (tundra↔desert) on a small island.
    /// REAL physics (temperature, wind, precip → biomes) — re-runs the climate. At the
    /// geographic value `domain_km/111` the run is the byte-identical windowed path.
    latitude_span_deg: f32,
    /// Framing offset (integer COARSE cells) applied to the coarse field for both
    /// the preview and the export (auto ± manual pan), stored mod grid. A cyclic
    /// rotation on the torus — no bounds, wraps at the edges. This is what
    /// `run_hd` writes to `window_offset_in_torus`.
    offset_cells: [i64; 2],
    /// Last auto-computed framing offset (cells) from the preview — the target of
    /// "snap to auto" and the reference for [`panned`].
    auto_offset_cells: [i64; 2],
    /// The offset has been manually panned away from the auto value.
    panned: bool,
    /// Seed text edited but not yet submitted — debounced: the coarse recompute
    /// fires on focus loss / Enter, not per keystroke (seeds are 20-digit).
    seed_dirty: bool,
    /// The initial auto-preview has been requested (so it fires once at startup).
    bootstrapped: bool,
    /// Opt-in: write a v1 `.ymir` delivery container after the biome phase.
    export_ymir: bool,
    /// EXPERIMENTAL (ADR 0001): use routed stream-power incision instead of the
    /// droplet pass for the next HD run. Off by default (production unchanged).
    stream_power: bool,
    /// EXPERIMENTAL (ADR 0001, Finding 7): with stream-power on, use `relief_v2` — the
    /// two bounding closures (nonlinear hillslope diffusion + lateral widening).
    closures: bool,
    /// EXPERIMENTAL (ADR 0001, Finding 9): cross-rill diffusion (diffuse everywhere) to
    /// damp the Smith–Bretherton comb; `cross_rill_d` is its coefficient.
    cross_rill: bool,
    cross_rill_d: f32,
    /// ADR 0001 Finding 109 -- "Age du continent": the absolute slope floor toggle.
    /// `false` = the delivered world, byte-identical.
    age_closure: bool,
    /// Its constant, restricted to Finding 106's MEASURED window {0.021, 0.024, 0.027}.
    age_s_eq: f32,
    /// ADR 0001 Finding 121 -- "Vallées construites": 0 = off (the delivered world), 1 = C1 bare
    /// construction, 2 = C2 hybrid at k_time/10, 3 = C2 hybrid at k_time/3, 4 = C2 /10 with
    /// Finding 122-B's `basin_base` (χ from the col of each closed depression).
    valley_mode: usize,
    /// Index into the three measured ages {0.7, 1, 1.4} × `F121_AGE_K`.
    valley_age: usize,
    /// ADR Finding 124-5 (E) -- W(k)'s γ, index into {0, 0.5, 1}: `W = a₀·(k/k₀)^γ·A^0.3`. γ = 0 is
    /// the ungated width (`None`, the bench's digest). PROXY.
    valley_width_gamma: usize,
    /// ADR Finding 132-P4 -- "Lac = niveau de base": a present lake (the construction input's pre-drainage
    /// lakes, `LakeBase::InputLakes`) stops χ at its shore. Off = `None`, the guarded digests.
    valley_lake_base: bool,
    /// ADR Finding 135-W4 -- "Mur ↔ mer" (Finding 126-B's clause): no wall cell laid below sea + 0.5 m
    /// (`wall_sea_floor_m = Some(0.5)`). Off = `None`, the guarded digests. Gated candidate, nothing promoted.
    valley_wall_sea: bool,
    /// ADR Finding 140 -- "Recul de la gorge (âge)": the gorge's retreat, r read from the age selector by map (a).
    /// Off = `None`, the guarded digests. Gated, nothing promoted.
    valley_gorge: bool,
    /// ADR Finding 140 -- the per-lake exponent p, an index into {0, 0.25, 0.5} (the author: default 0.5).
    valley_gorge_p: usize,
    /// EXPERIMENTAL (ADR 0001, Finding 11): MFD incision — dendritic valleys, no solver.
    mfd: bool,
    mfd_p: f32,
    /// ADR 0001 Finding 83 -- the base-level bound, now the PRODUCTION default, kept here
    /// as the A/B control the author asked not to remove. `true` (default) = production;
    /// UNCHECKING sets `HdParams::base_level_off` and reproduces the pre-Finding-83 world.
    base_level: bool,
    /// EXPERIMENTAL: FBM amplitude_base for the striation ladder (0.16 = production).
    /// EXPERIMENTAL (C-2, closures roadmap §2): inject volcanic edifices (arcs /
    /// hotspot chains / rifts) derived from the tectonic state. ON by default since ADR Finding 151 (the author's
    /// production is the guard's world). Arc volcanism is barely visible on seeds with
    /// few O-C margins — the hotspot chains and rifts carry the render.
    volcanism: bool,
    /// EXPERIMENTAL (C-3, closures roadmap §3): lithological heterogeneity — a
    /// causal per-cell erodibility multiplier (hard basement ×1.0, rift-soft ×10,
    /// volcaniclastic ×3). ON by default since ADR Finding 151. Only bites
    /// with stream-power on.
    lithology: bool,
    /// EXPERIMENTAL (C-3b, closures roadmap §3b): inherited structure — fracture-density
    /// erodibility (intact craton ×1 reference, fractured belts near plate contacts
    /// erode more). ON by default since ADR Finding 151. Only bites with
    /// stream-power on.
    fracture: bool,
    /// EXPERIMENTAL (H-1): infiltration — a causal permeability field (lithology matrix +
    /// fracture density, double porosity) removes part of the precipitation surplus from
    /// the SURFACE balance, so over-supplied basins can flip endorheic. Off by default.
    infiltration: bool,
    /// Destination directory for the `.ymir` container (`<dir>/<name>.ymir/`).
    export_dir: String,
    // Derived / cache.
    current: Option<Arc<HdResult>>,
    /// Fast tectonic-shape preview (coarse continent) shown before the HD run.
    preview: Option<Arc<PreviewShape>>,
    preview_tex: Option<egui::TextureHandle>,
    /// Framing offset the `preview_tex` was built at — rebuild the texture when the
    /// user pans so the drawn frame follows the offset.
    preview_tex_offset: [i64; 2],
    /// Sub-cell drag accumulator (pixels) for smooth drag-to-pan on the preview.
    pan_accum: [f32; 2],
    /// Zoom factor for the tectonic PREVIEW (1× = fit; >1 shows scroll bars to pan
    /// the magnified view — distinct from `zoom`, the HD map zoom).
    preview_zoom: f32,
    /// Pan CENTRE (uv, [0,1]²) of the zoomed HD result view — driven by the X/Y
    /// scroll sliders; clamped so the zoom window stays inside the map.
    map_pan: [f32; 2],
    river_map: Option<RiverCellMap>,
    /// Aggregated watercourses (rebuilt with each HD result). Assembly of exported data.
    watercourses: Option<Arc<Vec<Watercourse>>>,
    /// Selected entity (index into `watercourses` / `drainage.lakes`) for the microscope.
    selected_river: Option<usize>,
    selected_lake: Option<usize>,
    selected_volcano: Option<usize>,
    /// segment index → watercourse index (for hover-select on the map).
    seg_to_wc: Vec<usize>,
    inspect_tab: InspectTab,
    texture: Option<egui::TextureHandle>,
    tex_layer: Option<HdLayer>,
    /// Overlay state the cached `texture` was built with (rebuild when it flips).
    tex_overlay: bool,
    /// Tectonic-debug overlays (the microscope) currently enabled.
    overlays: TectonicOverlays,
    /// Overlay settings the cached `texture` was built with (rebuild on change).
    tex_overlays: TectonicOverlays,
    /// ADR Finding 135-V -- how the Relief layer is drawn: 0 hypsometry (the shipped layer), 1 hillshade,
    /// 2 the signed difference against `diff_ref`. A VIEW: it reads the world and never writes it.
    relief_view: usize,
    /// ADR Finding 157-V -- the lakes' outline on the hillshade (a view; on by default, can be hidden).
    lake_outline: bool,
    /// ADR Finding 135-V -- the reference world of the difference view, stored by the author.
    diff_ref: Option<Arc<HdResult>>,
    /// ADR Finding 135-V -- the difference's saturation, an index into `DIFF_SAT_M`.
    diff_sat: usize,
    /// ADR Finding 135-V -- (cells with |Δz| ≥ the saturation, cells with Δz ≠ 0) of the last difference texture.
    diff_stats: Option<(usize, usize)>,
    hover: Option<CellInspection>,
    hover_xy: Option<(usize, usize)>,
    /// Active map tool: SELECT (microscope — hover/click inspects an entity) or PAN
    /// (hand — click-drag moves the view). Wheel/CTRL-wheel work in either tool.
    tool: MapTool,
    /// Pipeline frieze collapsed (hidden) to give the canvas more room.
    pipeline_collapsed: bool,
    /// Bird-eye minimap shown over the canvas.
    minimap: bool,
    /// Show map symbols (C-2 volcano markers) over the canvas.
    show_symbols: bool,
    /// Microscope SELECTION drawer (left) open.
    micro_open: bool,
    /// The selected river's profile PATH cells (source→sink), parallel to the plotted profile —
    /// so a hover on the long profile can place a matching point on the map highlight.
    profile_path: Vec<(u32, u32)>,
    /// Fraction (0=source, 1=sink) currently hovered on the long profile, or `None`.
    profile_hover: Option<f32>,
    /// Pending "frame this element" request from a microscope double-click: (uv centre, uv extent).
    /// The map consumes it (pan to the centre, zoom so the extent fits) then clears it.
    focus_target: Option<([f32; 2], f32)>,
    /// Absolute on-screen size (px) of the FULL map square — the camera scale. `0` = fit on the
    /// next frame. Kept absolute (not derived from the panel-dependent available size) so resizing
    /// the side panels does NOT rescale the canvas.
    map_px: f32,
}

/// Map interaction tool (mock: two buttons — select / hand).
#[derive(Clone, Copy, PartialEq, Eq)]
enum MapTool {
    Select,
    Pan,
}

impl Default for WorkspaceState {
    fn default() -> Self {
        Self {
            seed: 42, // restored pre-M1 default (island preset is opt-in, not default)
            resolution: 2048,
            latitude: 45.0,
            mode: Mode::Standard,
            layer: HdLayer::Relief,
            river_overlay: false,
            rll_overlay: false,
            rll_min_len_km: 0.0,
            rll_min_area_km2: 0.0,
            rll_min_strahler: 1,
            rll_hide_fusion: false,
            rll_hide_retrace: false,
            geology_rules_path: "geology_rules.toml".to_string(),
            chance_only: None,
            chance_overlay: false,
            zoning_override: None,
            inspector_open: true,
            climat_open: true,
            relief_open: false,
            drainage_open: false,
            craton_density: 2.85,
            shield_frac: 0.45,
            shelf_width: 140.0,
            erosion: 0.5,
            channel_jitter: 0.3,
            dinf: false,
            domain_km: 1024.0,
            geo_scale_ratio: 1.0,
            latitude_span_deg: 1024.0 / 111.0, // geographic span at the default domain (identity)
            offset_cells: [0, 0],
            auto_offset_cells: [0, 0],
            panned: false,
            seed_dirty: false,
            bootstrapped: false,
            export_ymir: false,
            // relief-v3 is the PRODUCTION default (#190): the author has generated + validated every
            // map through stream-power + closures + MFD for many rounds. Leaving it opt-in silently
            // reverted a fresh run to the old droplet terrain. FBM amplitude defaults to 0.04 (the
            // production seed's value — fine relief detail); the toggle stays in Expert mode only.
            stream_power: true,
            closures: true,
            cross_rill: false,
            age_closure: false, // ADR Finding 109 -- off ships
            age_s_eq: 0.024,    // the middle of the Finding 106 window
            valley_mode: 0,     // ADR Finding 121 -- off ships
            valley_age: 1,      // k × 1
            valley_width_gamma: 0, // ADR Finding 124-5 -- γ 0: the width ignores the age
            valley_lake_base: false, // ADR Finding 132-P4 -- off ships
            valley_wall_sea: false,  // ADR Finding 135-W4 -- off ships
            valley_gorge: false,     // ADR Finding 140 -- off ships
            valley_gorge_p: 2,       // ADR Finding 140 -- p = 0.5, the author's default
            cross_rill_d: 0.40,
            mfd: true,
            mfd_p: 2.0,
            base_level: true,    // ADR Finding 83 -- production ships the bound
            // ADR Finding 151-0 — the author's production is the guard's world: C-2, C-3, C-3b ON at startup
            volcanism: true,
            lithology: true,
            fracture: true,
            infiltration: false, // H-1 opt-in (Expert); production byte-identical
            export_dir: "exports".to_string(),
            current: None,
            preview: None,
            preview_tex: None,
            preview_tex_offset: [0, 0],
            pan_accum: [0.0, 0.0],
            preview_zoom: 1.0,
            river_map: None,
            watercourses: None,
            selected_river: None,
            selected_lake: None,
            selected_volcano: None,
            seg_to_wc: Vec::new(),
            inspect_tab: InspectTab::Rivers,
            texture: None,
            tex_layer: None,
            tex_overlay: false,
            overlays: TectonicOverlays::default(),
            tex_overlays: TectonicOverlays::default(),
            relief_view: 0, // ADR Finding 135-V -- the shipped layer
            lake_outline: true,
            diff_ref: None,
            diff_sat: 2,
            diff_stats: None,
            hover: None,
            hover_xy: None,
            map_pan: [0.5, 0.5],
            tool: MapTool::Select,
            pipeline_collapsed: false,
            minimap: true,
            show_symbols: true,
            micro_open: true,
            map_px: 0.0,
            profile_path: Vec::new(),
            profile_hover: None,
            focus_target: None,
        }
    }
}

/// Convert a DATA-cell bounding box to a texture-uv (centre, extent) for the map camera, honouring
/// the north-up mirror (data cell y → texture row h-1-y). `extent` is the larger uv side.
fn bbox_to_uv(minx: u32, maxx: u32, miny: u32, maxy: u32, w: usize, h: usize) -> ([f32; 2], f32) {
    let (wf, hf) = (w as f32, h as f32);
    let ux0 = minx as f32 / wf;
    let ux1 = (maxx + 1) as f32 / wf;
    let uy0 = (h as u32 - 1 - maxy) as f32 / hf; // min data-y is the BOTTOM texture row
    let uy1 = (h as u32 - miny) as f32 / hf;
    let cx = (ux0 + ux1) * 0.5;
    let cy = (uy0 + uy1) * 0.5;
    let ext = (ux1 - ux0).max(uy1 - uy0).max(1.0 / wf);
    ([cx, cy], ext)
}

/// Bounding box (uv centre + extent) of a watercourse's cells, for double-click framing.
fn river_bbox_uv(hd: &HdResult, wc: &Watercourse) -> Option<([f32; 2], f32)> {
    let (mut minx, mut miny, mut maxx, mut maxy) = (u32::MAX, u32::MAX, 0u32, 0u32);
    let mut any = false;
    for &s in &wc.segments {
        for &(x, y) in &hd.drainage.rivers.segments[s].points {
            minx = minx.min(x);
            maxx = maxx.max(x);
            miny = miny.min(y);
            maxy = maxy.max(y);
            any = true;
        }
    }
    any.then(|| bbox_to_uv(minx, maxx, miny, maxy, hd.width, hd.height))
}

/// Bounding box (uv centre + extent) of a lake's footprint, for double-click framing.
fn lake_bbox_uv(hd: &HdResult, lake_id: u32) -> Option<([f32; 2], f32)> {
    let (w, h) = (hd.width, hd.height);
    let (mut minx, mut miny, mut maxx, mut maxy) = (u32::MAX, u32::MAX, 0u32, 0u32);
    let mut any = false;
    for y in 0..h {
        for x in 0..w {
            if hd.drainage.lake_map[y * w + x] == lake_id {
                minx = minx.min(x as u32);
                maxx = maxx.max(x as u32);
                miny = miny.min(y as u32);
                maxy = maxy.max(y as u32);
                any = true;
            }
        }
    }
    any.then(|| bbox_to_uv(minx, maxx, miny, maxy, w, h))
}

fn seg_label(ui: &mut egui::Ui, text: &str, active: bool) -> egui::Response {
    let (bg, fg) = if active { (COPPER, INK) } else { (C::TRANSPARENT, DIM) };
    ui.add(
        egui::Button::new(egui::RichText::new(text).color(fg).size(11.0))
            .fill(bg)
            .corner_radius(4.0)
            .min_size(egui::vec2(0.0, 22.0)),
    )
}

/// A framed segmented control (equal-width buttons in a recessed track), as the
/// mock's resolution / flow selectors. Returns the clicked index.
fn seg_row(ui: &mut egui::Ui, labels: &[&str], active: usize) -> Option<usize> {
    let mut clicked = None;
    egui::Frame::default()
        .fill(FIELD)
        .stroke(egui::Stroke::new(1.0, C::from_rgb(0x2e, 0x2e, 0x2e)))
        .inner_margin(3)
        .corner_radius(6)
        .show(ui, |ui| {
            let n = labels.len() as f32;
            let w = ((ui.available_width() - 3.0 * (n - 1.0)) / n).max(10.0);
            ui.spacing_mut().item_spacing.x = 3.0;
            ui.horizontal(|ui| {
                for (i, l) in labels.iter().enumerate() {
                    let (bg, fg) = if i == active { (COPPER, INK) } else { (C::TRANSPARENT, DIM) };
                    if ui
                        .add(
                            egui::Button::new(egui::RichText::new(*l).color(fg).size(12.0))
                                .fill(bg)
                                .corner_radius(4.0)
                                .min_size(egui::vec2(w, 24.0)),
                        )
                        .clicked()
                    {
                        clicked = Some(i);
                    }
                }
            });
        });
    clicked
}

/// A left/right-panel content block padded to the mock's 16 px gutter, with
/// per-block top/bottom padding. Separators between blocks stay full-bleed.
fn block<R>(
    ui: &mut egui::Ui,
    gutter: i8,
    top: i8,
    bottom: i8,
    add: impl FnOnce(&mut egui::Ui) -> R,
) -> R {
    egui::Frame::default()
        .inner_margin(egui::Margin { left: gutter, right: gutter, top, bottom })
        .show(ui, add)
        .inner
}

/// A small uppercase field label ("GRAINE", "RÉSOLUTION").
fn field_label(ui: &mut egui::Ui, t: &str) {
    ui.label(egui::RichText::new(t).color(DIM2).size(10.0));
}

/// An expert parameter row (the mock's pattern): a label (left) + its value in
/// copper mono (right), then a full-width slider below (no inline value).
fn param_slider(
    ui: &mut egui::Ui,
    label: &str,
    v: &mut f32,
    range: std::ops::RangeInclusive<f32>,
    decimals: usize,
) {
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(label).color(C::from_rgb(0x9a, 0x9a, 0x9a)).size(11.0));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(
                egui::RichText::new(format!("{:.*}", decimals, *v))
                    .color(COPPER_BRIGHT)
                    .monospace()
                    .size(11.0),
            );
        });
    });
    ui.add_space(2.0);
    ui.spacing_mut().slider_width = ui.available_width();
    ui.add(egui::Slider::new(v, range).show_value(false));
    ui.add_space(10.0);
}

/// A sub-heading inside an expert section (e.g. "ISOSTASIE").
fn sub_heading(ui: &mut egui::Ui, t: &str) {
    ui.add_space(4.0);
    ui.label(egui::RichText::new(t).color(C::from_rgb(0x5e, 0x5e, 0x5e)).size(9.5));
    ui.add_space(6.0);
}

/// Custom latitude slider matching the mock: a warm→polar gradient track, a
/// copper thumb, and tick marks + labels aligned under 0/30/45/60/75/90.
fn latitude_slider(ui: &mut egui::Ui, lat: &mut f32) {
    let w = ui.available_width();
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(w, 40.0), egui::Sense::click_and_drag());
    let p = ui.painter_at(rect);
    let x0 = rect.left();
    let x1 = rect.right();
    let span = (x1 - x0).max(1.0);
    let track_y = rect.top() + 6.0;

    // Gradient track (equator warm → polar white).
    let stops = [
        (0.0, C::from_rgb(0xE1, 0x78, 0x46)),
        (0.34, C::from_rgb(0x6E, 0xBE, 0x6E)),
        (0.67, C::from_rgb(0x5A, 0x8C, 0xCD)),
        (1.0, C::from_rgb(0xE1, 0xEB, 0xF8)),
    ];
    let steps = 64;
    for s in 0..steps {
        let t = s as f32 / (steps - 1) as f32;
        let x = x0 + span * t;
        p.rect_filled(
            egui::Rect::from_min_size(
                egui::pos2(x, track_y),
                egui::vec2(span / steps as f32 + 1.0, 4.0),
            ),
            0.0,
            grad_at(&stops, t),
        );
    }

    // Drag / click to set.
    if (resp.dragged() || resp.clicked()) {
        if let Some(pos) = resp.interact_pointer_pos() {
            let f = ((pos.x - x0) / span).clamp(0.0, 1.0);
            *lat = (f * 90.0).round();
        }
    }

    // Thumb.
    let tx = x0 + (*lat / 90.0).clamp(0.0, 1.0) * span;
    let tc = egui::pos2(tx, track_y + 2.0);
    p.circle_filled(tc, 7.5, C::from_rgb(0x2a, 0x21, 0x18));
    p.circle_filled(tc, 6.0, COPPER);

    // Ticks + labels.
    let ticks = [
        (0.0, "équat."),
        (30.0, "subtrop."),
        (45.0, "tempéré"),
        (60.0, "subpol."),
        (90.0, "polaire"),
    ];
    for (deg, label) in ticks {
        let tickx = x0 + (deg / 90.0) * span;
        p.line_segment(
            [egui::pos2(tickx, track_y + 10.0), egui::pos2(tickx, track_y + 15.0)],
            egui::Stroke::new(1.0, C::from_rgb(0x3a, 0x3a, 0x3a)),
        );
        let (align, lx) = if deg == 0.0 {
            (egui::Align2::LEFT_TOP, tickx)
        } else if deg == 90.0 {
            (egui::Align2::RIGHT_TOP, tickx)
        } else {
            (egui::Align2::CENTER_TOP, tickx)
        };
        p.text(
            egui::pos2(lx, track_y + 17.0),
            align,
            format!("{deg:.0}°"),
            egui::FontId::monospace(8.5),
            DIM2,
        );
        p.text(
            egui::pos2(lx, track_y + 27.0),
            align,
            label,
            egui::FontId::proportional(8.0),
            C::from_rgb(0x55, 0x55, 0x55),
        );
    }
}

fn section_header(ui: &mut egui::Ui, open: &mut bool, title: &str, badge: Option<(&str, C)>) {
    ui.horizontal(|ui| {
        let arrow = if *open { "▼" } else { "▶" };
        if ui
            .add(
                egui::Label::new(egui::RichText::new(arrow).color(COPPER).size(10.0))
                    .sense(egui::Sense::click()),
            )
            .clicked()
        {
            *open = !*open;
        }
        if ui
            .add(
                egui::Label::new(egui::RichText::new(title).color(TEXT_BRIGHT).strong().size(12.5))
                    .sense(egui::Sense::click()),
            )
            .clicked()
        {
            *open = !*open;
        }
        if let Some((b, col)) = badge {
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(egui::RichText::new(b).color(col).size(9.0));
            });
        }
    });
}

fn draw_workspace(
    mut contexts: EguiContexts,
    bridge: Res<C1SolverBridge>,
    mut ws: ResMut<WorkspaceState>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    // Latch a freshly-completed HD result (identity via Arc pointer).
    if let HdState::Completed { result, .. } = &bridge.hd {
        let is_new = ws.current.as_ref().map(|c| !Arc::ptr_eq(c, result)).unwrap_or(true);
        if is_new {
            ws.current = Some(result.clone());
            ws.zoning_override = None; // ADR Finding 152 -- a new world: its own zoning
            GUARD_REFUSES.store(
                result.bench_guard.refuses_numbers() || result.lake_guard.refuses_numbers(),
                std::sync::atomic::Ordering::Relaxed,
            );
            ws.river_map = Some(RiverCellMap::from_drainage(&result.drainage));
            let wcs = aggregate_watercourses(result, ws.domain_km, ws.geo_scale_ratio);
            let mut seg2wc = vec![usize::MAX; result.drainage.rivers.segments.len()];
            for (wi, wc) in wcs.iter().enumerate() {
                for &s in &wc.segments {
                    seg2wc[s] = wi;
                }
            }
            ws.watercourses = Some(Arc::new(wcs));
            ws.seg_to_wc = seg2wc;
            ws.selected_river = None;
            ws.selected_lake = None;
            ws.texture = None;
            ws.tex_layer = None;
            ws.hover = None;
            // Default zoom 100% (fit) on a fresh map — refit next frame, recentred.
            ws.map_px = 0.0;
            ws.map_pan = [0.5, 0.5];
        }
    }
    // Latch a fresh tectonic-shape preview → show it (and drop any stale HD map
    // so the preview is what's on screen until the user hits Générer).
    if let Some(preview) = &bridge.preview {
        let is_new = ws.preview.as_ref().map(|p| !Arc::ptr_eq(p, preview)).unwrap_or(true);
        if is_new {
            ws.auto_offset_cells = preview.auto_offset_cells;
            // A fresh preview follows a seed change (which cleared `panned`), so the
            // framing snaps to the new auto offset. A manual pan set `panned`, which
            // survives a mere refresh of the same seed.
            if !ws.panned {
                ws.offset_cells = preview.auto_offset_cells;
            }
            ws.preview = Some(preview.clone());
            ws.preview_tex = None;
            ws.current = None;
            ws.texture = None;
            ws.hover = None;
        }
    }
    let hd_running = matches!(bridge.hd, HdState::Running { .. });

    // One-time auto-preview at startup: the default seed shows a framed continent
    // without a click (coarse pass only — never the HD pipeline).
    if !ws.bootstrapped && ws.preview.is_none() && !hd_running {
        request_preview(&mut ws, &bridge, true);
        ctx.request_repaint();
    }

    // Keep animating (waiter pulse) + polling the worker while a run is live.
    if hd_running {
        ctx.request_repaint();
    }

    top_bar(ctx, &bridge, &mut ws);
    left_panel(ctx, &bridge, &mut ws, hd_running);
    microscope_drawer(ctx, &mut ws);
    right_panel(ctx, &mut ws);
    central_panel(ctx, &bridge, &mut ws);
}

/// Microscope SELECTION as a LEFT drawer (second left panel), resizable, collapsible to a thin
/// re-open strip. Holds the entity list; the selected element's DETAIL shows in the right zone.
fn microscope_drawer(ctx: &egui::Context, ws: &mut WorkspaceState) {
    if ws.current.is_none() {
        return;
    }
    if !ws.micro_open {
        egui::SidePanel::left("micro_closed")
            .exact_width(28.0)
            .frame(egui::Frame::default().fill(PANEL2).inner_margin(4))
            .show(ctx, |ui| {
                ui.vertical_centered(|ui| {
                    if ui
                        .button(egui::RichText::new("⌖").size(15.0))
                        .on_hover_text("Ouvrir le microscope")
                        .clicked()
                    {
                        ws.micro_open = true;
                    }
                });
            });
        return;
    }
    egui::SidePanel::left("microscope")
        .default_width(248.0)
        .width_range(180.0..=420.0)
        .resizable(true)
        .frame(egui::Frame::default().fill(PANEL).inner_margin(10))
        .show(ctx, |ui| microscope_list(ui, ws));
}

// ── Top bar ──────────────────────────────────────────────────────────────
fn top_bar(ctx: &egui::Context, bridge: &C1SolverBridge, ws: &mut WorkspaceState) {
    egui::TopBottomPanel::top("topbar")
        .exact_height(34.0)
        .frame(egui::Frame::default().fill(PANEL2).inner_margin(egui::Margin::symmetric(12, 0)))
        .show(ctx, |ui| {
            ui.horizontal_centered(|ui| {
                ui.spacing_mut().item_spacing.x = 6.0;
                ui.label(egui::RichText::new("◇").color(COPPER).size(13.0));
                ui.label(
                    egui::RichText::new("YMIR")
                        .color(C::from_rgb(0xD8, 0xD8, 0xD8))
                        .strong()
                        .size(12.0),
                );
                ui.add_space(14.0);
                for m in ["Fichier", "Génération", "Couches", "Vue", "Aide"] {
                    ui.label(
                        egui::RichText::new(m).color(C::from_rgb(0x7d, 0x7d, 0x7d)).size(12.0),
                    );
                    ui.add_space(8.0);
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    // Mode toggle (Standard / Expert).
                    if seg_label(ui, "Expert", ws.mode == Mode::Expert).clicked() {
                        ws.mode = Mode::Expert;
                    }
                    if seg_label(ui, "Standard", ws.mode == Mode::Standard).clicked() {
                        ws.mode = Mode::Standard;
                    }
                    ui.label(egui::RichText::new("MODE").color(DIM2).size(10.0));
                    ui.add_space(12.0);
                    // Stats.
                    let cells = ws.current.as_ref().map(|c| c.width * c.height).unwrap_or(0);
                    let stat = format!(
                        "{}² px  ·  {:.0} km  ·  {:.0} m/px  ·  {} cellules  ·  seed {}",
                        ws.resolution,
                        ws.domain_km,
                        m_per_px(ws.domain_km, ws.resolution),
                        cells,
                        ws.seed
                    );
                    ui.label(egui::RichText::new(stat).color(DIM2).monospace().size(11.0));
                    // M1 island-continent verdict from the land-topology diagnostics.
                    if let Some(cur) = ws.current.as_ref() {
                        let t = &cur.land_topology;
                        let island = t.num_landmasses > 0 && !t.wraps_x && !t.wraps_y;
                        let (verdict, col) = if island {
                            ("île ✓", C::from_rgb(0x6a, 0xb0, 0x6a))
                        } else {
                            ("continent tore ✗", C::from_rgb(0xc0, 0x7a, 0x4a))
                        };
                        ui.add_space(12.0);
                        let land = format!(
                            "terre {:.0}%  ·  masse max {:.0} km²  ·  {}",
                            t.largest_area_frac * 100.0,
                            t.largest_area_km2,
                            verdict,
                        );
                        ui.label(egui::RichText::new(land).color(col).monospace().size(11.0));
                    }
                    let _ = bridge;
                });
            });
        });
}

// ── Left control panel ───────────────────────────────────────────────────
fn left_panel(
    ctx: &egui::Context,
    bridge: &C1SolverBridge,
    ws: &mut WorkspaceState,
    hd_running: bool,
) {
    egui::SidePanel::left("controls")
        .exact_width(296.0)
        .frame(egui::Frame::default().fill(PANEL).inner_margin(0))
        .show(ctx, |ui| {
            ui.spacing_mut().item_spacing = egui::vec2(8.0, 8.0);

            // GENERATE + progress pinned at the bottom (declared first so it
            // reserves the bottom strip; the central area fills the rest).
            egui::TopBottomPanel::bottom("generate_bar")
                .frame(egui::Frame::default().fill(C::from_rgb(0x19, 0x19, 0x19)).inner_margin(egui::Margin { left: 16, right: 16, top: 13, bottom: 15 }))
                .show_inside(ui, |ui| {
                    let label = if hd_running {
                        "GÉNÉRATION EN COURS…"
                    } else if ws.current.is_some() {
                        "RÉGÉNÉRER"
                    } else {
                        "GÉNÉRER"
                    };
                    let btn = egui::Button::new(egui::RichText::new(label).color(INK).strong().size(13.0))
                        .fill(COPPER)
                        .corner_radius(7.0)
                        .min_size(egui::vec2(ui.available_width(), 38.0));
                    if ui.add_enabled(!hd_running, btn).clicked() {
                        let spec = C1RunSpec { seed: ws.seed, ..C1RunSpec::default() };
                        // `Some(dir)` = export enabled; empty path falls back to
                        // "exports" so the checkbox alone is enough to opt in.
                        let export_dir = if ws.export_ymir {
                            let d = ws.export_dir.trim();
                            Some(PathBuf::from(if d.is_empty() { "exports" } else { d }))
                        } else {
                            None
                        };
                        // Ship exactly the framing the user sees (auto ± manual pan),
                        // as a normalised offset in [0,1) — cells / grid.
                        let grid = ws.preview.as_ref().map(|p| p.coarse.width as i64).unwrap_or(64);
                        let off = [
                            ws.offset_cells[0].rem_euclid(grid) as f64 / grid as f64,
                            ws.offset_cells[1].rem_euclid(grid) as f64 / grid as f64,
                        ];
                        let params = HdParams {
                            target_size: ws.resolution,
                            latitude_deg: ws.latitude,
                            domain_km: ws.domain_km,
                            manual_offset: Some(off),
                            stream_power: ws.stream_power,
                            closures: ws.closures,
                            cross_rill: ws.cross_rill,
                            // ADR Finding 121 -- the hybrid (modes 2/3) carries A1+B2 as its
                            // finish, exactly as the bench does, whatever the age box says.
                            slope_floor_s_eq: if ws.valley_mode >= 2 {
                                Some(0.024)
                            } else {
                                ws.age_closure.then_some(ws.age_s_eq)
                            },
                            valley_construction: {
                                use ymir_core::tectonics_c1::valley_construction::{
                                    F121_AGE_K, ValleyConstruction,
                                };
                                let k = F121_AGE_K * [0.7f32, 1.0, 1.4][ws.valley_age.min(2)];
                                match ws.valley_mode {
                                    // ADR Finding 124 -- `new` IS the definition (basin base, the
                                    // author's decision D); mode 4 keeps Finding 121's sea base
                                    1 => Some(ValleyConstruction::new(k, None)),
                                    2 => Some(ValleyConstruction::new(k, Some(0.1))),
                                    3 => Some(ValleyConstruction::new(k, Some(1.0 / 3.0))),
                                    4 => Some(ValleyConstruction::f121(k, Some(0.1))),
                                    _ => None,
                                }
                                // ADR Finding 124-5 (E) -- W(k); γ 0 stays `None` (same width, and
                                // the digest the guard knows)
                                .map(|vc| ValleyConstruction {
                                    width_age_gamma: [None, Some(0.5), Some(1.0)]
                                        [ws.valley_width_gamma.min(2)],
                                    // ADR Finding 132-P4 -- off stays `None` (the guarded digest)
                                    lake_base: (ws.valley_lake_base || ws.valley_gorge).then_some(
                                        // ADR Finding 133 -- the extended variant (basin lakes are present lakes)
                                        ymir_core::tectonics_c1::valley_construction::LakeBase::InputLakesAndBasins,
                                    ),
                                    // ADR Finding 135-W4 -- mur ↔ mer, off stays `None` (the guarded digest)
                                    wall_sea_floor_m: ws.valley_wall_sea.then_some(0.5),
                                    // ADR Finding 140 -- the gorge's retreat: r from the age selector by map (a)
                                    // (×0.7 → 0, ×1 → 1, ×1.4 → 2), p from its selector
                                    gorge_retreat: ws.valley_gorge.then(|| {
                                        ymir_core::tectonics_c1::valley_construction::GorgeRetreat::v4(
                                            [0.0f32, 1.0, 2.0][ws.valley_age.min(2)],
                                            [0.0f32, 0.25, 0.5][ws.valley_gorge_p.min(2)],
                                        )
                                    }),
                                    ..vc
                                })
                            },
                            cross_rill_d: ws.cross_rill_d,
                            mfd: ws.mfd,
                            mfd_p: ws.mfd_p,
                            base_level_m: None, // shipped epsilon (Finding 83)
                            base_level_off: !ws.base_level,
                            geo_scale_ratio: ws.geo_scale_ratio,
                            // At the geographic span → None (byte-identical windowed climate);
                            // otherwise the explicit span drives the placed (per-belt) climate.
                            latitude_span_deg: {
                                let geo = ws.domain_km / 111.0;
                                if (ws.latitude_span_deg - geo).abs() < 0.05 {
                                    None
                                } else {
                                    Some(ws.latitude_span_deg)
                                }
                            },
                            export_dir,
                            // C-2 opt-in. `domain_km` is overwritten with the run's
                            // geometric span in the bridge, so only `enabled` matters here.
                            volcanism: ws.volcanism.then(|| {
                                ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig {
                                    enabled: true,
                                    ..Default::default()
                                }
                            }),
                            // C-3 opt-in. Chosen multipliers (see the sweep report):
                            // soft rift ×10, volcaniclastic ×3. `enabled` gates it.
                            lithology: ws.lithology.then(|| {
                                ymir_core::tectonics_c1::closures::lithology::LithologyConfig {
                                    enabled: true,
                                    soft_multiplier: 10.0,
                                    volcanic_multiplier: 3.0,
                                    rift_age_threshold: 1.0,
                                }
                            }),
                            // C-3b opt-in. Chosen from the sweep: amplitude ×6 (visible
                            // over ×4 without doubling the lake-population pressure that
                            // ×8 would add), narrow orogenic belt (decay 25 km).
                            fracture: ws.fracture.then(|| {
                                ymir_core::tectonics_c1::closures::fracture::FractureConfig {
                                    enabled: true,
                                    amplitude: 6.0,
                                    decay_km: 25.0,
                                    ..Default::default()
                                }
                            }),
                            // H-1 opt-in. Defaults anchored on the BFI upper range
                            // (f_cap 0.7) and the rainfall supply rate (k_ref).
                            infiltration: ws.infiltration.then(|| {
                                ymir_core::tectonics_c1::closures::infiltration::InfiltrationConfig {
                                    enabled: true,
                                    ..Default::default()
                                }
                            }),
                            // Debug microscope overlay — derive the tectonic labels so
                            // the overlay panel has data (one cheap coarse pass).
                            emit_tectonic_labels: true,
                            geology_rules: Some(std::path::PathBuf::from(&ws.geology_rules_path)),
                        };
                        let _ = bridge.submit_hd(spec, params);
                    }
                    // Manual refresh of the coarse preview (auto-preview fires on
                    // seed change; this re-runs the same seed, keeping any pan).
                    ui.add_space(6.0);
                    let prev_btn = egui::Button::new(
                        egui::RichText::new("RAFRAÎCHIR L'APERÇU").color(BRONZE).strong().size(11.5),
                    )
                    .fill(C::from_rgb(0x24, 0x22, 0x1e))
                    .stroke(egui::Stroke::new(1.0, BRONZE))
                    .corner_radius(6.0)
                    .min_size(egui::vec2(ui.available_width(), 28.0));
                    if ui.add_enabled(!hd_running, prev_btn).clicked() {
                        request_preview(ws, bridge, false);
                    }
                    ui.add_enabled_ui(!hd_running, |ui| {
                        // relief-v3 is the PRODUCTION recipe and the default (#190). In STANDARD mode it
                        // is FORCED ON — only the tuning (FBM amplitude, MFD p) is exposed, so a run can
                        // never silently fall back to the old droplet/relief-v1 terrain. EXPERT mode keeps
                        // the full enable/disable toggles for experimentation.
                        if ws.mode == Mode::Expert {
                            ui.checkbox(
                                &mut ws.stream_power,
                                egui::RichText::new("Incision stream-power (exp.)").color(DIM2).size(11.0),
                            )
                            .on_hover_text(
                                "Remplace l'érosion par gouttelettes par l'incision stream-power \
                                 routée + versants (ADR 0001). Off = ancien pipeline gouttelettes.",
                            );
                            if ws.stream_power {
                                ui.checkbox(
                                    &mut ws.closures,
                                    egui::RichText::new("Closures relief-v2 (exp.)").color(DIM2).size(11.0),
                                )
                                .on_hover_text(
                                    "Les deux fermetures qui BORNENT le relief (ADR 0001, Finding 7): \
                                     diffusion de versant NON LINÉAIRE à pente critique + élargissement \
                                     latéral. Off = relief-v1 (slits 1 px).",
                                );
                                ui.checkbox(
                                    &mut ws.base_level,
                                    egui::RichText::new("Niveau de base (F83, production)")
                                        .color(DIM2)
                                        .size(11.0),
                                )
                                .on_hover_text(
                                    "ADR 0001 Finding 83 — le niveau de base EST la production \
                                     depuis le 2026-09-15: l'incision ne creuse plus sous \
                                     sea + 0,5 m, et la frange côtière est à l'autorité \
                                     pré-incision (Δ éperons ≥ 2 cellules = 0 aux quatre \
                                     résolutions).\n\n\
                                     DÉCOCHER = LE MONDE D'AVANT (contrôle A/B): la frange \
                                     revient, cuvettes sous-marines 20 → 62, budget +5,8 %, \
                                     largeurs de rivière +13 %.\n\n\
                                     ⚠ Lacs et rivières sont recalculés dans les deux sens, \
                                     pas seulement le raster (Findings 81-82).",
                                );
                                // ADR 0001 Finding 109 — the AGE closure, on the Finding 83
                                // toggle pattern: the box reaches `eroded_key`, so lakes, rivers,
                                // biomes and spillways are RE-DERIVED rather than redrawn from the
                                // other terrain's cache entry.
                                ui.checkbox(
                                    &mut ws.age_closure,
                                    egui::RichText::new("Âge du continent — constante (F109)")
                                        .color(DIM2)
                                        .size(11.0),
                                )
                                .on_hover_text(
                                    "ADR 0001 Findings 106/108/109 — le plancher de pente ABSOLU: \
                                     S_eq = s_eq·A^(-1/2), plus A1 (exclusion des cellules en \
                                     dépression au moment de l'incision, F104). UNE passe, AUCUNE \
                                     calibration, aucun étage à nommer.\n\n\
                                     COCHER = le continent a VIEILLI: la classe sur-creusée tombe \
                                     12 → 0 (seed 1), la part de lacs D_L > 5 va de 32 % à 0 %, R8 \
                                     de 0,092 à 0,048.\n\n\
                                     DÉCOCHER = le livré, AU BIT (garde sur trois seeds).\n\n\
                                     ⚠ Lacs, rivières, biomes et spillways sont RECALCULÉS dans \
                                     les deux sens — la clé de cache contient s_eq — pas seulement \
                                     le raster.\n\n\
                                     ⚠ RIEN N'EST PROMU: relief_v3 est inchangé, le défaut est \
                                     décoché.",
                                );
                                if ws.age_closure {
                                    // ⚠️ The measured WINDOW and nothing else. Outside
                                    // {0.021, 0.024, 0.027} no row was measured, on any seed,
                                    // under either class definition — a free slider here would
                                    // invite exactly the extrapolation Findings 106-108 refused.
                                    let sv = [0.021f32, 0.024, 0.027];
                                    let sc = sv
                                        .iter()
                                        .position(|&v| (v - ws.age_s_eq).abs() < 1e-6)
                                        .unwrap_or(1);
                                    ui.horizontal(|ui| {
                                        ui.label(
                                            egui::RichText::new("s_eq").color(DIM2).size(11.0),
                                        )
                                        .on_hover_text(
                                            "La FENÊTRE mesurée (F106, relue au F108 sous la \
                                             classe différentielle) — non vide sur les trois \
                                             seeds. 0,024 en est le milieu. Rien n'est offert en \
                                             dehors: aucune ligne n'y a jamais été mesurée.",
                                        );
                                        if let Some(i) =
                                            seg_row(ui, &["0.021", "0.024", "0.027"], sc)
                                        {
                                            ws.age_s_eq = sv[i];
                                        }
                                    });
                                }
                                // ADR 0001 Finding 121 — the constructed valleys. Four states
                                // with the age box: livré / A1+B2 / C1 nue / C2 hybride.
                                ui.horizontal(|ui| {
                                    ui.label(
                                        egui::RichText::new("Vallées construites (F121)")
                                            .color(DIM2)
                                            .size(11.0),
                                    )
                                    .on_hover_text(
                                        "ADR 0001 Finding 121 — la FORME construite le long des \
                                         troncs du réseau pré-incision (A ≥ 10 km²), la TEXTURE \
                                         laissée à une incision légère.\n\n\
                                         off = le livré, AU BIT.\n\
                                         C1 nue = vallées posées (fond χ, largeur W = a·A^0,3, \
                                         versants 28°), AUCUNE incision.\n\
                                         C2 /10, /3 = C1 puis UNE passe de stream power à \
                                         k_time/10 ou k_time/3, avec A1+B2 (s_eq 0,024) comme \
                                         finition — forcé, comme au banc.\n\
                                         Depuis le F124, la base au col est la DÉFINITION \
                                         (χ recompté depuis le col de chaque dépression fermée: \
                                         les vallées ne descendent plus sous l'eau d'un bassin). \
                                         « C2 /10 niveau mer » montre l'ancienne définition (F121).\n\n\
                                         Lois: fond ANCRÉ (Harel 2016), exposant de largeur \
                                         ANCRÉ (Clubb 2022), 28° ANCRÉ (Whipple & Tucker 1999); \
                                         l'âge k et le coefficient de largeur sont des PROXY.\n\n\
                                         ⚠ Le monde vu ici est celui du banc f121_hybrid à seed 1, \
                                         8192², 400 km. À une autre résolution, AUCUN verdict \
                                         (F89-C5).\n\n\
                                         ⚠ RIEN N'EST PROMU: le défaut est off.",
                                    );
                                    // ADR Finding 123 — a combo box: a five-button row was
                                    // clipped to three states in a narrow panel.
                                    let names =
                                        // ADR Finding 124 -- STABLE names: the state the author
                                        // validated keeps its name when the default changes
                                        ["off", "C1 nue", "C2 /10 col (défaut)", "C2 /3", "C2 /10 niveau mer"];
                                    egui::ComboBox::from_id_salt("valley_mode")
                                        .selected_text(names[ws.valley_mode.min(4)])
                                        .show_ui(ui, |ui| {
                                            for (i, nm) in names.iter().enumerate() {
                                                ui.selectable_value(&mut ws.valley_mode, i, *nm);
                                            }
                                        });
                                });
                                if ws.valley_mode > 0 {
                                    ui.horizontal(|ui| {
                                        ui.label(
                                            egui::RichText::new("âge k").color(DIM2).size(11.0),
                                        )
                                        .on_hover_text(
                                            "Les trois âges MESURÉS au bloc C3 du F121: \
                                             k × {0,7 ; 1 ; 1,4}, k = 0,07183 (PROXY calé sur \
                                             l'oracle F95, 488 m). Rien d'autre n'est offert.",
                                        );
                                        if let Some(i) =
                                            seg_row(ui, &["×0.7", "×1", "×1.4"], ws.valley_age)
                                        {
                                            ws.valley_age = i;
                                        }
                                    });
                                    ui.horizontal(|ui| {
                                        ui.label(
                                            egui::RichText::new("W(k) γ · PROXY")
                                                .color(DIM2)
                                                .size(11.0),
                                        )
                                        .on_hover_text(
                                            "ADR Finding 124-5 (E) — la vallée s'élargit avec \
                                             l'âge: W = a₀·(k/k₀)^γ·A^0,3, k₀ = 0,07183.\n\n\
                                             γ est un PROXY: l'élargissement latéral est une \
                                             affaire de TEMPS, et le a de Clubb 2022 est un \
                                             instantané des vallées actuelles, pas un taux — rien \
                                             n'ancre γ.\n\n\
                                             γ 0 = la largeur des F121–F124 (indépendante de k). \
                                             γ 0,5 et 1 ne sont gardés par aucun banc: le badge \
                                             dira « non gardé ».",
                                        );
                                        if let Some(i) =
                                            seg_row(ui, &["0", "0,5", "1"], ws.valley_width_gamma)
                                        {
                                            ws.valley_width_gamma = i;
                                        }
                                    });
                                    ui.checkbox(
                                        &mut ws.valley_lake_base,
                                        egui::RichText::new("Lac = niveau de base (F133)").color(DIM2).size(11.0),
                                    )
                                    .on_hover_text(
                                        "ADR Findings 132-P4 / 133-F — décisions de l'auteur.\n\n\
                                         2026-09-29 : « Un lac présent est un niveau de base pour les rivières \
                                         qui s'y jettent ; un lac vidé ne l'est plus. »\n\
                                         2026-09-30 : « Les lacs des bassins sous la mer sont des lacs présents. »\n\n\
                                         χ s'arrête au rivage de chaque lac présent : l'amont compte χ depuis le lac, \
                                         avec pour base la surface du lac (son col).\n\n\
                                         ⚠ Circularité : les lacs de lakes.json n'existent qu'après la \
                                         construction. « Présent » = les lacs du pré-drainage du champ d'ENTRÉE \
                                         de la construction et les dépressions fermées de la base au col, \
                                         remplies à leur col (LakeBase::InputLakesAndBasins). Le point fixe \
                                         n'est pas construit.\n\n\
                                         ⚠ Aucun banc ne garde cet état : le badge dira « non gardé ». \
                                         RIEN N'EST PROMU.",
                                    );
                                    ui.checkbox(
                                        &mut ws.valley_wall_sea,
                                        egui::RichText::new("Mur ↔ mer (F126-B)").color(DIM2).size(11.0),
                                    )
                                    .on_hover_text(
                                        "ADR Findings 126-B / 135-W4 — la clause « mur ↔ mer » : aucune cellule de \
                                         MUR de vallée n'est posée sous la mer + 0,5 m (le fond, lui, peut \
                                         l'être). Elle portait tous les bancs des F127–F134 ; le viz ne la posait \
                                         jamais.\n\n\
                                         Candidate gatée. Aucun banc ne garde cet état : le badge dira « non \
                                         gardé ». RIEN N'EST PROMU.",
                                    );
                                    // ADR Finding 141 -- HIDDEN until Finding 140's crash is fixed (the gorge at ×1, p 0.5 fires
                                    // Finding 38's invariant in the HD assembly): `valley_gorge` stays false, the world is
                                    // the guarded one
                                    if GORGE_TOGGLE_VISIBLE {
                                    ui.checkbox(
                                        &mut ws.valley_gorge,
                                        egui::RichText::new("Recul de la gorge (âge)").color(DIM2).size(11.0),
                                    )
                                    .on_hover_text(
                                        "ADR Finding 140 — spec v3. Chaque lac présent (≥ 1 km²) se tient à un \
                                         niveau L(r) entre son déversoir et son fond ; son exutoire part en \
                                         GORGE accrochée à ce niveau (pente 10 × la loi gradée, raidie si la \
                                         place manque, plafond 28°) ; son rebord garde ≥ L(r) ; une chute de \
                                         tête (un lac sur trois sans) est taggée.\n\n\
                                         r suit le sélecteur d'âge : ×0,7 → 0 (lac plein), ×1 → 1 (lac sur \
                                         son seuil de fond), ×1,4 → 2 (vidé). Recul par lac : \
                                         r_lac = r · (A/418,7 km²)^p.\n\n\
                                         Active la base des lacs étendue. Aucun banc ne garde cet état : le \
                                         badge dira « non gardé ». RIEN N'EST PROMU.",
                                    );
                                    if ws.valley_gorge {
                                        ui.horizontal(|ui| {
                                            ui.label(egui::RichText::new("recul par lac p").color(DIM2).size(11.0))
                                                .on_hover_text(
                                                    "p = 0 : tous les lacs reculent ensemble. p = 0,5 : les \
                                                     grands exutoires reculent plus vite (la célérité d'une \
                                                     rupture croît avec l'aire, F115). L'auteur : défaut 0,5.",
                                                );
                                            if let Some(i) = seg_row(ui, &["0", "0,25", "0,5"], ws.valley_gorge_p) {
                                                ws.valley_gorge_p = i;
                                            }
                                        });
                                    }
                                    }
                                }
                                ui.checkbox(
                                    &mut ws.cross_rill,
                                    egui::RichText::new("cross_rill (GS non convergé, F9)").color(DIM2).size(11.0),
                                )
                                .on_hover_text(
                                    "STEP 2b (ADR 0001, Finding 9) — `cross_rill`: diffusion de \
                                     versant PARTOUT (on abandonne l'exclusion des chenaux du split \
                                     de régime) pour amortir le rilling de Smith–Bretherton.\n\n\
                                     ⚠ RENOMMÉE au F109. Elle s'appelait « Anti-peigne », ce qui \
                                     décrivait une INTENTION et non le code: l'ADR L15258 enregistre \
                                     que cette case lie la branche de versant GAUSS-SEIDEL NON \
                                     CONVERGÉE (F9), une troisième chose. Le peigne, lui, est traité \
                                     à sa cause par relief-v3 (MFD).",
                                );
                                if ws.cross_rill {
                                    let dv = [0.25f32, 0.40, 0.55];
                                    let dc = dv
                                        .iter()
                                        .position(|&v| (v - ws.cross_rill_d).abs() < 1e-4)
                                        .unwrap_or(1);
                                    ui.horizontal(|ui| {
                                        ui.label(egui::RichText::new("D").color(DIM2).size(11.0));
                                        if let Some(i) = seg_row(ui, &["0.25", "0.40", "0.55"], dc) {
                                            ws.cross_rill_d = dv[i];
                                        }
                                    });
                                }
                                ui.checkbox(
                                    &mut ws.mfd,
                                    egui::RichText::new("relief-v3 (MFD+talus+breach)").color(DIM2).size(11.0),
                                )
                                .on_hover_text(
                                    "ADR 0001 Findings 11/14/15 — la chaîne relief-v3 complète: incision \
                                     MULTI-FLUX + gradage TALUS + K×3, conditionnée par BREACH (profils \
                                     monotones, lacs à plat). Le correctif final recommandé.",
                                );
                            }
                            // C-2 volcanism (independent of the relief chain).
                            ui.checkbox(
                                &mut ws.volcanism,
                                egui::RichText::new("Volcanisme (C-2)").color(DIM2).size(11.0),
                            )
                            .on_hover_text(
                                "Closures roadmap §2 — édifices dérivés de la tectonique: arcs \
                                 insulaires (subduction), chaînes de points chauds (âge croissant), \
                                 rifts. Lacs de cratère acides si l'édifice dégaze. Le volcanisme \
                                 d'arc est peu visible sur les seeds à peu de marges O-C — regarder \
                                 les chaînes de hotspot et les rifts. Le log [HD volcanism] indique \
                                 ce qui a été produit.",
                            );
                            // C-3 lithology (needs stream-power on to bite).
                            ui.checkbox(
                                &mut ws.lithology,
                                egui::RichText::new("Lithologie (C-3)").color(DIM2).size(11.0),
                            )
                            .on_hover_text(
                                "Closures roadmap §3 — hétérogénéité lithologique CAUSALE (jamais du \
                                 bruit): érodabilité par cellule, socle dur ×1 (référence), rift tendre \
                                 ×10, footprints volcaniclastiques ×3. La roche tendre s'érode plus → \
                                 vallées ouvertes; le socle dur retient le relief. N'agit qu'avec la \
                                 stream-power. Le log [HD C-3 lithology] confirme.",
                            );
                            // C-3b inherited structure (needs stream-power on to bite).
                            ui.checkbox(
                                &mut ws.fracture,
                                egui::RichText::new("Structure héritée (C-3b)").color(DIM2).size(11.0),
                            )
                            .on_hover_text(
                                "Closures roadmap §3b — densité de fracturation CAUSALE (proximité des \
                                 contacts tectoniques, jamais du bruit ni cratonic_mask): érodabilité \
                                 isotrope, craton intact ×1 (référence, garde son relief), ceintures \
                                 fracturées près des orogènes/transformants s'érodent/disséquent plus. \
                                 L'orientation a été écartée (non atteignable avec C1 — voir ADR). N'agit \
                                 qu'avec la stream-power. Le log [HD C-3b fracture] confirme.",
                            );
                            // H-1 infiltration (agit sur le BILAN hydrique, pas la géométrie).
                            ui.checkbox(
                                &mut ws.infiltration,
                                egui::RichText::new("Infiltration (H-1)").color(DIM2).size(11.0),
                            )
                            .on_hover_text(
                                "H-1 — premier terme SOUTERRAIN. Un champ de perméabilité CAUSAL \
                                 (matrice lithologique + densité de fracturation, double porosité — \
                                 la fracturation gagne 5-6 ordres, cf. Heath) retire une fraction du \
                                 surplus de précipitation du bilan de SURFACE. Des bassins \
                                 sur-alimentés peuvent basculer en endoréique et des niveaux baisser. \
                                 Ne touche ni la géométrie ni le routage. Pas de terme de pente (la \
                                 littérature ne le soutient pas). Le log [HD H-1 infiltration] confirme.",
                            );
                        } else {
                            // STANDARD — the production recipe, forced on. No fallback toggle.
                            ws.stream_power = true;
                            ws.closures = true;
                            ws.mfd = true;
                            ws.cross_rill = false;
                            ui.label(
                                egui::RichText::new("Relief : stream-power + relief-v3 (défaut production)")
                                    .color(DIM2)
                                    .size(11.0),
                            )
                            .on_hover_text(
                                "Recette de production (ADR 0001): incision stream-power routée + closures \
                                 relief-v2 + MFD/talus/breach relief-v3. Réglages avancés en mode Expert.",
                            );
                        }
                        // Tuning exposed in BOTH modes (only when the recipe is on).
                        if ws.stream_power {
                            // ADR Finding 117 -- the "FBM amp" selector that stood here was INERT:
                            // since C-1 the relief-budget cap binds at every cell and
                            // `amplitude_base` changes nothing (ADR "The DEAD KNOB"; Finding 116
                            // re-proved it byte-identical at 8x). A dead control in the author's
                            // panel is a defect dressed as configuration. The real lever is
                            // `flow_conditioning` (beta), which carries two roles and must be
                            // split before it is exposed -- not this round.
                            if ws.mfd {
                                let pv = [4.0f32, 2.0, 1.1];
                                let pc = pv
                                    .iter()
                                    .position(|&v| (v - ws.mfd_p).abs() < 1e-4)
                                    .unwrap_or(1);
                                ui.horizontal(|ui| {
                                    ui.label(egui::RichText::new("p").color(DIM2).size(11.0));
                                    if let Some(i) = seg_row(ui, &["4", "2", "1.1"], pc) {
                                        ws.mfd_p = pv[i];
                                    }
                                });
                            }
                        }
                        ui.checkbox(
                            &mut ws.export_ymir,
                            egui::RichText::new("Exporter .ymir").color(DIM2).size(11.0),
                        );
                        ui.add_enabled_ui(ws.export_ymir, |ui| {
                            ui.horizontal(|ui| {
                                ui.label(egui::RichText::new("Dossier").color(DIM2).size(11.0));
                                ui.add(
                                    egui::TextEdit::singleline(&mut ws.export_dir)
                                        .hint_text("exports")
                                        .desired_width(f32::INFINITY),
                                );
                            });
                        });
                    });
                    progress_block(ui, &bridge.hd, bridge);
                });

            // Logo block (full-bleed tinted header).
            egui::Frame::default()
                .fill(C::from_rgb(0x1f, 0x1c, 0x18))
                .inner_margin(egui::Margin { left: 16, right: 16, top: 15, bottom: 14 })
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new("⬡").color(COPPER_BRIGHT).size(30.0));
                        ui.add_space(4.0);
                        ui.vertical(|ui| {
                            ui.spacing_mut().item_spacing.y = 3.0;
                            ui.label(egui::RichText::new("YMIR").color(TEXT_BRIGHT).strong().size(19.0));
                            ui.label(egui::RichText::new("CONTINENT GENERATOR").color(BRONZE).size(9.5));
                        });
                    });
                });
            ui.separator();

            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.add_enabled_ui(!hd_running, |ui| {
                    // Seed + resolution.
                    block(ui, 16, 14, 16, |ui| {
                        field_label(ui, "GRAINE");
                        ui.add_space(1.0);
                        ui.horizontal(|ui| {
                            let mut s = ws.seed.to_string();
                            let resp = ui.add(
                                egui::TextEdit::singleline(&mut s)
                                    .desired_width(ui.available_width() - 42.0)
                                    .font(egui::TextStyle::Monospace),
                            );
                            if resp.changed() {
                                if let Ok(v) = s.parse::<u64>() {
                                    if v != ws.seed {
                                        ws.seed = v;
                                        ws.seed_dirty = true; // submit on focus loss
                                    }
                                }
                            }
                            // Debounce: one coarse recompute when editing completes
                            // (Enter / focus loss), never per keystroke.
                            if resp.lost_focus() && ws.seed_dirty {
                                request_preview(ws, bridge, true);
                            }
                            if ui
                                .add(egui::Button::new(egui::RichText::new("⟳").color(COPPER_BRIGHT)).min_size(egui::vec2(34.0, 0.0)))
                                .on_hover_text("Graine aléatoire")
                                .clicked()
                            {
                                ws.seed = ws.seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                                request_preview(ws, bridge, true); // immediate
                            }
                        });
                        ui.add_space(10.0);
                        field_label(ui, "RÉSOLUTION");
                        ui.add_space(1.0);
                        let cur = [512usize, 1024, 2048, 4096, 8192].iter().position(|&r| r == ws.resolution).unwrap_or(1);
                        if let Some(i) = seg_row(ui, &["512", "1024", "2048", "4096", "8192"], cur) {
                            ws.resolution = [512, 1024, 2048, 4096, 8192][i];
                        }
                        ui.add_space(10.0);
                        // Domain size (km): the domain IS the map (no crop). The whole
                        // torus renders at `resolution`, so m_per_px = domain_km·1000/
                        // resolution (the single source of truth, `m_per_px`). Moving
                        // this relabels the SAME coarse field — NO tectonic recompute.
                        field_label(ui, "DOMAINE (km)");
                        ui.add_space(1.0);
                        let mpp = m_per_px(ws.domain_km, ws.resolution);
                        ui.horizontal(|ui| {
                            ui.add(egui::Slider::new(&mut ws.domain_km, 200.0..=1200.0).step_by(1.0).show_value(false));
                            ui.label(egui::RichText::new(format!("{:.0} km · {:.0} m/px", ws.domain_km, mpp)).color(COPPER_BRIGHT).monospace().size(12.0));
                        });
                        // Compact domain-as-map readout for the current seed + domain,
                        // recomputed instantly from the coarse preview field (cheap 64²
                        // — no tectonic pass). Shown once a preview exists.
                        domain_readout(ui, ws, mpp);
                        ui.add_space(10.0);
                        // Finding 24 — geographic scale ratio (hydrology only). Instant:
                        // a pure post-process on river width/discharge/navigability. The
                        // terrain, climate and biomes stay on the real domain_km.
                        field_label(ui, "ÉCHELLE GÉOGRAPHIQUE (hydrologie)");
                        ui.add_space(1.0);
                        ui.horizontal(|ui| {
                            ui.add(egui::Slider::new(&mut ws.geo_scale_ratio, 1.0..=15.0).step_by(0.1).show_value(false));
                            ui.label(egui::RichText::new(format!("×{:.1} · signifie {:.0} km", ws.geo_scale_ratio, ws.domain_km * ws.geo_scale_ratio)).color(COPPER_BRIGHT).monospace().size(12.0));
                        });
                        ui.label(egui::RichText::new("largeur/débit des rivières uniquement — n'affecte ni le relief ni les biomes").color(DIM).size(10.0));
                        ui.add_space(10.0);
                        // Manual framing pan (cyclic; composes on top of the auto offset).
                        framing_controls(ui, ws);
                        ui.add_space(8.0);
                        if ui
                            .add(egui::Button::new(
                                egui::RichText::new("Réinitialiser").color(DIM).size(11.0),
                            ))
                            .on_hover_text("Cadrage auto + domaine par défaut, graine inchangée")
                            .clicked()
                        {
                            reset_current_seed(ws);
                        }
                    });
                    ui.separator();

                    // Climat section (latitude — the wired climate knob).
                    let mut climat_open = ws.climat_open;
                    block(ui, 16, 12, if climat_open { 4 } else { 12 }, |ui| {
                        section_header(ui, &mut climat_open, "Climat", Some(("EXPRESSIF", C::from_rgb(0x5a, 0x7a, 0x5a))));
                    });
                    if climat_open {
                        block(ui, 16, 0, 16, |ui| {
                            ui.horizontal(|ui| {
                                ui.label(egui::RichText::new("Latitude du continent").color(C::from_rgb(0xbd, 0xbd, 0xbd)).size(11.5));
                                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                    ui.label(egui::RichText::new(format!("{:.0}° ({})", ws.latitude, zone_name(ws.latitude))).color(COPPER_BRIGHT).monospace().size(13.0));
                                });
                            });
                            ui.add_space(8.0);
                            latitude_slider(ui, &mut ws.latitude);
                            // Finding 25 — climatic latitude SPAN, independent of the
                            // physical size. Widen to cross belts (real physics → re-run).
                            ui.add_space(10.0);
                            let geo = ws.domain_km / 111.0;
                            let tag = if (ws.latitude_span_deg - geo).abs() < 0.05 { "géographique" } else { "élargie" };
                            ui.horizontal(|ui| {
                                ui.label(egui::RichText::new("Étendue latitudinale").color(C::from_rgb(0xbd, 0xbd, 0xbd)).size(11.5));
                                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                    ui.label(egui::RichText::new(format!("{:.1}° ({tag})", ws.latitude_span_deg)).color(COPPER_BRIGHT).monospace().size(13.0));
                                });
                            });
                            ui.add(egui::Slider::new(&mut ws.latitude_span_deg, 1.0..=40.0).step_by(0.1).show_value(false));
                            ui.label(egui::RichText::new("élargir pour croiser plusieurs ceintures (toundra ↔ désert)").color(DIM).size(10.0));
                            ui.add_space(8.0);
                            // TASK 1 — placement widget: shows the map rectangle over the
                            // latitude gradient + wind belts, live as the sliders move.
                            latitude_placement_widget(ui, ws.latitude, ws.latitude_span_deg);
                        });
                    }
                    ui.separator();

                    match ws.mode {
                        Mode::Standard => {
                            block(ui, 16, 14, 16, |ui| {
                                egui::Frame::default().fill(C::from_rgb(0x1a, 0x16, 0x11)).stroke(egui::Stroke::new(1.0, C::from_rgb(0x2c, 0x24, 0x17))).inner_margin(11).corner_radius(6).show(ui, |ui| {
                                    ui.label(egui::RichText::new("Plaques, isostasie, érosion, bathymétrie & drainage — réglages fins dans le mode Expert.").color(C::from_rgb(0x8a, 0x7a, 0x60)).size(10.5));
                                });
                                ui.add_space(9.0);
                                let b = egui::Button::new(egui::RichText::new("Passer en Expert →").color(COPPER_BRIGHT).size(11.5)).fill(C::from_rgb(0x22, 0x1c, 0x14)).min_size(egui::vec2(ui.available_width(), 30.0));
                                if ui.add(b).clicked() {
                                    ws.mode = Mode::Expert;
                                }
                            });
                        }
                        Mode::Expert => {
                            block(ui, 16, 12, 16, |ui| {
                                expert_sections(ui, ws);
                            });
                        }
                    }
                    ui.separator();
                });
            });
        });
}

fn expert_sections(ui: &mut egui::Ui, ws: &mut WorkspaceState) {
    ui.label(
        egui::RichText::new("Ces réglages sont EXPOSÉS ; câblage moteur à suivre.")
            .color(C::from_rgb(0x6a, 0x5a, 0x40))
            .size(9.5)
            .italics(),
    );
    ui.add_space(4.0);

    let mut relief_open = ws.relief_open;
    section_header(ui, &mut relief_open, "Relief & closures", Some(("AVANCÉ", DIM2)));
    ws.relief_open = relief_open;
    if ws.relief_open {
        ui.add_space(6.0);
        sub_heading(ui, "ISOSTASIE");
        param_slider(ui, "Densité cratonique", &mut ws.craton_density, 2.6..=3.1, 2);
        param_slider(ui, "Fraction bouclier", &mut ws.shield_frac, 0.0..=1.0, 2);
        sub_heading(ui, "BATHYMÉTRIE & ÉROSION");
        param_slider(ui, "Largeur plateau cont.", &mut ws.shelf_width, 40.0..=320.0, 0);
        param_slider(ui, "Intensité d'érosion", &mut ws.erosion, 0.0..=1.0, 2);
    }
    ui.add_space(6.0);
    let mut drainage_open = ws.drainage_open;
    section_header(ui, &mut drainage_open, "Drainage", Some(("AVANCÉ", DIM2)));
    ws.drainage_open = drainage_open;
    if ws.drainage_open {
        ui.add_space(6.0);
        param_slider(ui, "Perturbation du tracé", &mut ws.channel_jitter, 0.0..=1.0, 2);
        ui.label(
            egui::RichText::new("Mode d'écoulement")
                .color(C::from_rgb(0x9a, 0x9a, 0x9a))
                .size(11.0),
        );
        ui.add_space(3.0);
        let cur = if ws.dinf { 1 } else { 0 };
        if let Some(i) = seg_row(ui, &["D8", "D∞"], cur) {
            ws.dinf = i == 1;
        }
    }
}

fn progress_block(ui: &mut egui::Ui, hd: &HdState, _bridge: &C1SolverBridge) {
    match hd {
        HdState::Idle => {}
        HdState::Running { current, done, .. } => {
            ui.add_space(8.0);
            for r in done {
                ui.label(
                    egui::RichText::new(format!(
                        "✓ {} — {} ({:.1}s)",
                        r.phase.label(),
                        r.regime.label(),
                        r.elapsed.as_secs_f32()
                    ))
                    .color(DIM)
                    .size(10.5),
                );
            }
            if let Some(p) = current {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(
                        egui::RichText::new(format!("{} …", p.label())).color(TEXT).size(11.5),
                    );
                });
            }
            ui.horizontal(|ui| {
                ui.label(
                    egui::RichText::new("● Interface active pendant le calcul")
                        .color(GREEN)
                        .size(10.0),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui
                        .add(
                            egui::Label::new(
                                egui::RichText::new("Annuler")
                                    .color(C::from_rgb(0xcc, 0x77, 0x77))
                                    .size(11.0),
                            )
                            .sense(egui::Sense::click()),
                        )
                        .clicked()
                    {
                        _bridge.request_cancel();
                    }
                });
            });
        }
        HdState::Completed { done, total, .. } => {
            ui.add_space(8.0);
            ui.label(
                egui::RichText::new(format!(
                    "⚡ Rendu complet · {} étapes · {:.1}s",
                    done.len(),
                    total.as_secs_f32()
                ))
                .color(C::from_rgb(0x7a, 0x7a, 0x7a))
                .size(10.5),
            );
        }
        HdState::Failed { error } => {
            ui.add_space(8.0);
            ui.label(
                egui::RichText::new(format!("Échec : {error}"))
                    .color(C::from_rgb(0xE1, 0x78, 0x46))
                    .size(11.0),
            );
        }
    }
}

// ── Right inspector ──────────────────────────────────────────────────────
fn right_panel(ctx: &egui::Context, ws: &mut WorkspaceState) {
    if !ws.inspector_open {
        egui::SidePanel::right("inspector_closed")
            .exact_width(34.0)
            .frame(egui::Frame::default().fill(PANEL))
            .show(ctx, |ui| {
                if ui
                    .add(
                        egui::Label::new(
                            egui::RichText::new("‹\nINSPECTION").color(DIM).size(11.0),
                        )
                        .sense(egui::Sense::click()),
                    )
                    .clicked()
                {
                    ws.inspector_open = true;
                }
            });
        return;
    }
    egui::SidePanel::right("inspector")
        .default_width(320.0)
        .width_range(240.0..=620.0)
        .resizable(true)
        .frame(egui::Frame::default().fill(PANEL).inner_margin(0))
        .show(ctx, |ui| {
            ui.spacing_mut().item_spacing = egui::vec2(8.0, 6.0);
            block(ui, 15, 13, 13, |ui| {
                ui.horizontal(|ui| {
                    ui.label(
                        egui::RichText::new("Inspection").color(TEXT_BRIGHT).strong().size(12.5),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("›").on_hover_text("Replier").clicked() {
                            ws.inspector_open = false;
                        }
                    });
                });
            });
            ui.separator();
            let hover = ws.hover; // copy before the &mut borrow below
            let geo = ws.hover_xy.and_then(|(x, y)| geo_cell(ws, x, y));
            egui::ScrollArea::vertical().show(ui, |ui| {
                // Primary: the selected element's detail (microscope).
                block(ui, 15, 14, 14, |ui| microscope_detail(ui, ws));
                // Secondary: the cell currently under the cursor.
                ui.separator();
                match hover {
                    Some(c) => {
                        group_title(ui, "CELLULE SURVOLÉE");
                        block(ui, 8, 14, 14, |ui| inspection(ui, &c, geo.as_ref()));
                    }
                    None => {
                        ui.add_space(10.0);
                        ui.label(
                            egui::RichText::new(if ws.current.is_some() {
                                "Survolez la carte pour inspecter une cellule."
                            } else {
                                "Générez un continent."
                            })
                            .color(C::from_rgb(0x7a, 0x7a, 0x7a))
                            .size(11.5),
                        );
                    }
                }
            });
        });
}

/// ADR 0001 Finding 88-C — the `Unresolved` reason in full, for the cell panel.
fn unresolved_reason_fr(r: Option<UnresolvedReason>) -> &'static str {
    match r {
        Some(UnresolvedReason::NoExteriorNeighbour) => "l'empreinte n'a aucun voisin extérieur",
        Some(UnresolvedReason::SaddleHasNoLowerNeighbour) => {
            "le col n'a aucun voisin extérieur plus bas (pas d'échappée, F37c)"
        }
        Some(UnresolvedReason::ReturnsIntoOwnFootprint) => {
            "la descente retombe dans sa propre empreinte (boucle, F37c)"
        }
        Some(UnresolvedReason::NoD8Direction) => "la descente s'arrête sur un plat (DIR_NONE)",
        Some(UnresolvedReason::NoSinkInOneGrid) => "aucun puits en une grille de pas",
        Some(UnresolvedReason::OutletOffGrid) => "exutoire hors grille",
        Some(UnresolvedReason::NoOutletReach) => {
            "aucun tronçon exporté ne part de son empreinte (F86)"
        }
        None => "raison non renseignée",
    }
}

/// The same, abbreviated for the lake list and the river source chip.
fn unresolved_reason_short_fr(r: Option<UnresolvedReason>) -> &'static str {
    match r {
        Some(UnresolvedReason::NoExteriorNeighbour) => "⚠ aucun voisin extérieur",
        Some(UnresolvedReason::SaddleHasNoLowerNeighbour) => "⚠ col sans échappée",
        Some(UnresolvedReason::ReturnsIntoOwnFootprint) => "⚠ boucle sur son empreinte",
        Some(UnresolvedReason::NoD8Direction) => "⚠ plat (DIR_NONE)",
        Some(UnresolvedReason::NoSinkInOneGrid) => "⚠ aucun puits",
        Some(UnresolvedReason::OutletOffGrid) => "⚠ exutoire hors grille",
        Some(UnresolvedReason::NoOutletReach) => "⚠ aucun tronçon exporté",
        None => "⚠ sans exutoire",
    }
}

/// ADR Finding 123 — set when the displayed world FAILED the identity guard: every key/value number
/// of the microscope is then replaced by "≠ banc" (the author's rule: the viz refuses to display a
/// number when the viz and bench hashes differ).
static GUARD_REFUSES: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn kv(ui: &mut egui::Ui, k: &str, v: String) {
    let v = if GUARD_REFUSES.load(std::sync::atomic::Ordering::Relaxed) {
        "≠ banc".to_string()
    } else {
        v
    };
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(k).color(DIM).size(12.0));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(
                egui::RichText::new(v).color(C::from_rgb(0xe6, 0xe6, 0xe6)).monospace().size(12.0),
            );
        });
    });
    ui.separator();
}

fn group_title(ui: &mut egui::Ui, t: &str) {
    ui.add_space(6.0);
    ui.label(egui::RichText::new(t).color(BRONZE).size(9.5));
}

/// ADR Finding 153 -- the geology of the hovered cell: its rock class and its resources (strongest first), read from the
/// active zoning (the reloaded rules if any).
struct GeoCell {
    rock: u8,
    resources: Vec<ymir_core::geology::zoning::CellResource>,
}

fn geo_cell(ws: &WorkspaceState, x: usize, y: usize) -> Option<GeoCell> {
    use ymir_core::geology::zoning::{CHANCE_FACTOR, DensitySource, explain};
    let g = ws.current.as_ref()?.geology.as_ref()?;
    let (z, rules) = match &ws.zoning_override {
        Some((z, l, _)) => (z, &l.rules),
        None => (&g.zoning, &g.rules.rules),
    };
    let ctx = &g.product.context;
    let q = (y / CHANCE_FACTOR).min(ctx.h4 - 1) * ctx.w4 + (x / CHANCE_FACTOR).min(ctx.w4 - 1);
    Some(GeoCell { rock: *g.product.rocks.get(y * g.product.w + x)?, resources: explain(ctx, rules, z, q, DensitySource::Structural) })
}

fn inspection(ui: &mut egui::Ui, c: &CellInspection, geo: Option<&GeoCell>) {
    let [br, bg, bb] = c.biome.color();
    egui::Frame::default()
        .fill(FIELD)
        .stroke(egui::Stroke::new(1.0, C::from_rgb(0x2a, 0x2a, 0x2a)))
        .inner_margin(egui::Margin::symmetric(9, 5))
        .corner_radius(5)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(egui::RichText::new("⬛").color(C::from_rgb(br, bg, bb)).size(11.0));
                ui.label(
                    egui::RichText::new(format!("x {} · y {}", c.x, c.y))
                        .color(C::from_rgb(0xbb, 0xbb, 0xbb))
                        .monospace()
                        .size(11.0),
                );
            });
        });

    group_title(ui, "GÉOLOGIE");
    kv(ui, "Altitude", format!("{:+.0} m", c.altitude_m));
    if let Some(d) = c.depth_m {
        kv(ui, "Profondeur", format!("−{d:.0} m"));
    }
    // ADR Finding 151-0 — the coarse tectonic state under the cell (6.25 km cells), when the labels were built
    match c.tectonic {
        Some((continental, craton, s)) => {
            kv(ui, "Épaisseur crustale / S̃", format!("{s:.2} (sans dim.)"));
            kv(ui, "Type de plaque", if continental { "continentale".into() } else { "océanique".into() });
            kv(ui, "Craton", if craton { "oui".into() } else { "non".into() });
        }
        None => {
            kv(ui, "Épaisseur crustale / S̃", "—".into());
            kv(ui, "Type de plaque", "—".into());
            kv(ui, "Craton", "—".into());
        }
    }

    group_title(ui, "CLIMAT");
    kv(ui, "Température", format!("{:.1} °C", c.temperature_c));
    kv(ui, "Précipitation", format!("{:.0} mm/an", c.precip_mm));

    group_title(ui, "HYDROLOGIE");
    let drainage = match &c.river {
        Some(r) => match r.navigability {
            Navigability::Ship => "Rivière — navire".to_string(),
            Navigability::Barge => "Rivière — chaland".to_string(),
            Navigability::SmallBoat => "Rivière — barque".to_string(),
            Navigability::NonNavigable => "Rivière".to_string(),
        },
        None => match &c.lake {
            Some(l) => match l.lake_type {
                LakeType::Exorheic => "Lac (exoréique)".to_string(),
                LakeType::Endorheic => "Lac (endoréique)".to_string(),
                LakeType::CraterAcidic => "Lac de cratère (acide)".to_string(),
                LakeType::CraterNeutral => "Lac de cratère (eau douce)".to_string(),
                // ADR Finding 86 — the balance says it overflows and the trace does not resolve.
                // Finding 88-C: say WHICH failure, because two different defects wore these words.
                LakeType::Unresolved => format!(
                    "⚠ Lac SANS EXUTOIRE tracé — {}",
                    unresolved_reason_fr(l.unresolved_reason)
                ),
            },
            None => "—".to_string(),
        },
    };
    kv(ui, "Drainage", drainage);
    if let Some(r) = &c.river {
        kv(ui, "Débit", format!("{:.1} m³/s", r.discharge_m3s));
        kv(ui, "Largeur chenal", format!("{:.0} m", r.width_m));
        kv(ui, "Bassin versant", format!("{:.0} km²", r.drainage_km2));
    }
    kv(ui, "Ruissellement", format!("{:.0} mm/an", c.runoff_mm));

    group_title(ui, "BIOME");
    egui::Frame::default()
        .fill(FIELD)
        .stroke(egui::Stroke::new(1.0, C::from_rgb(0x2a, 0x2a, 0x2a)))
        .inner_margin(egui::Margin::symmetric(12, 10))
        .corner_radius(6)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(egui::RichText::new("⬛").color(C::from_rgb(br, bg, bb)).size(16.0));
                ui.add_space(2.0);
                ui.label(egui::RichText::new(french_biome(c.biome)).color(TEXT_BRIGHT).size(13.0));
            });
        });

    // ADR Finding 153 -- the rock and the resources of the cell
    if let Some(g) = geo {
        let rc = &ymir_core::geology::rocks::ROCK_CLASSES[(g.rock as usize).min(8)];
        group_title(ui, "ROCHES");
        egui::Frame::default()
            .fill(FIELD)
            .stroke(egui::Stroke::new(1.0, C::from_rgb(0x2a, 0x2a, 0x2a)))
            .inner_margin(egui::Margin::symmetric(12, 10))
            .corner_radius(6)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new("⬛").color(C::from_rgb(rc.color[0], rc.color[1], rc.color[2])).size(16.0));
                    ui.add_space(2.0);
                    ui.label(egui::RichText::new(rc.name_fr).color(TEXT_BRIGHT).size(13.0));
                });
                if g.rock != 0 {
                    kv(ui, "Dureté", ["tendre", "moyenne", "dure"][rc.hardness.min(2) as usize].into());
                }
            });
        group_title(ui, "RESSOURCES");
        ui.label(egui::RichText::new("favorabilité relative 0–100, lue sur la grille ¼ (195 m), comme l'export").color(DIM).size(9.5));
        if g.resources.is_empty() {
            ui.label(egui::RichText::new("aucune").color(DIM).size(11.0));
        }
        for r in &g.resources {
            ui.horizontal(|ui| {
                ui.label(egui::RichText::new("■").color(C::from_rgb(r.color[0], r.color[1], r.color[2])).size(12.0));
                ui.label(egui::RichText::new(&r.name_fr).color(TEXT_BRIGHT).size(11.5));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(egui::RichText::new(format!("{}", r.fav)).color(COPPER_BRIGHT).monospace().size(11.5));
                });
            });
            let why = match &r.origin {
                ymir_core::geology::zoning::FavOrigin::Rule { index, source, modulation } => {
                    let m = modulation.map_or(String::new(), |g| format!(" · modulation structurale × {g:.2}"));
                    format!("règle {} : {}{m}", index + 1, source.clone().unwrap_or_else(|| "(sans source)".into()))
                }
                ymir_core::geology::zoning::FavOrigin::Placer { source } => format!("placer : {}", source.clone().unwrap_or_else(|| "(sans source)".into())),
            };
            ui.label(egui::RichText::new(why).color(DIM).size(9.5));
        }
    }
}

// ── Central: frieze + map ────────────────────────────────────────────────
fn central_panel(ctx: &egui::Context, bridge: &C1SolverBridge, ws: &mut WorkspaceState) {
    egui::CentralPanel::default()
        .frame(egui::Frame::default().fill(C::from_rgb(0x0f, 0x0f, 0x0f)).inner_margin(0))
        .show(ctx, |ui| {
            // Pipeline frieze — a top strip (#161616) matching the mock. Collapsible (a thin bar
            // remains to re-expand it), so the canvas can take the full height.
            egui::TopBottomPanel::top("frieze_strip")
                .frame(egui::Frame::default().fill(PANEL2).inner_margin(egui::Margin {
                    left: 16,
                    right: 16,
                    top: 8,
                    bottom: if ws.pipeline_collapsed { 8 } else { 14 },
                }))
                .show_inside(ui, |ui| {
                    ui.horizontal(|ui| {
                        let chevron = if ws.pipeline_collapsed { "▸" } else { "▾" };
                        if ui
                            .add(
                                egui::Button::new(
                                    egui::RichText::new(format!("{chevron} PIPELINE"))
                                        .color(DIM2)
                                        .size(9.5),
                                )
                                .frame(false),
                            )
                            .on_hover_text("Replier / déplier le pipeline")
                            .clicked()
                        {
                            ws.pipeline_collapsed = !ws.pipeline_collapsed;
                        }
                        if let HdState::Running { current: Some(p), progress, .. } = &bridge.hd {
                            ui.with_layout(
                                egui::Layout::right_to_left(egui::Align::Center),
                                |ui| {
                                    let pct = progress
                                        .map(|(d, t)| {
                                            format!(" {}%", if t > 0 { d * 100 / t } else { 0 })
                                        })
                                        .unwrap_or_default();
                                    ui.label(
                                        egui::RichText::new(format!("{}{} …", p.label(), pct))
                                            .color(COPPER_BRIGHT)
                                            .size(10.5),
                                    );
                                    ui.spinner();
                                },
                            );
                        }
                    });
                    if !ws.pipeline_collapsed {
                        frieze(ui, &bridge.hd, ws);
                    }
                });
            // Map viewport (#0A0A0A) — the microscope entity list is now the LEFT drawer and its
            // detail the RIGHT zone, so the canvas takes the full central area.
            egui::CentralPanel::default()
                .frame(egui::Frame::default().fill(VIEWPORT).inner_margin(0))
                .show_inside(ui, |ui| map(ui, ws));
        });
}

/// Node state for `phase` derived from the live HD event stream (step e).
fn node_vis(hd: &HdState, has_result: bool, phase: HdPhase) -> NodeVis {
    match hd {
        HdState::Running { current, done, .. } => {
            if let Some(r) = done.iter().find(|r| r.phase == phase) {
                NodeVis::Done { cached: r.regime == CacheRegime::Hit }
            } else if *current == Some(phase) {
                NodeVis::Running
            } else {
                NodeVis::Pending
            }
        }
        HdState::Completed { done, .. } => {
            let cached = done
                .iter()
                .find(|r| r.phase == phase)
                .map(|r| r.regime == CacheRegime::Hit)
                .unwrap_or(false);
            NodeVis::Done { cached }
        }
        _ => {
            if has_result {
                NodeVis::Done { cached: false }
            } else {
                NodeVis::Pending
            }
        }
    }
}

/// The canonical frieze node highlighted for a viewed layer (selector mode):
/// Relief → the "Relief" node (not Tectonique/Érosion, which share that layer).
/// The frieze index that OWNS a layer (the selection ring lands here). Derived
/// from the table so it always points at the selectable node carrying `layer`
/// — notably Relief → Érosion (index 2), not the non-selectable Relief node.
fn canonical_node(layer: HdLayer) -> usize {
    FRIEZE.iter().position(|(_, _, v)| *v == Some(layer)).unwrap_or(0)
}

fn frieze(ui: &mut egui::Ui, hd: &HdState, ws: &mut WorkspaceState) {
    // Header (PIPELINE label + running indicator) is drawn by the collapsible wrapper.
    ui.add_space(2.0);
    let full = ui.available_width();
    let (rect, _) = ui.allocate_exact_size(egui::vec2(full, 44.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let n = FRIEZE.len();
    let margin = 30.0;
    let y = rect.top() + 12.0;
    let x0 = rect.left() + margin;
    let x1 = rect.right() - margin;
    let node_x = |i: usize| x0 + (x1 - x0) * (i as f32 / (n - 1) as f32);
    let time = ui.ctx().input(|i| i.time) as f32;

    // Three ORTHOGONAL states, each on its own visual dimension (they overlay):
    //   1. completeness  → node FILL + the copper line (persists; never retracts)
    //   2. running       → an animated PROGRESS ARC (distinct from a filled dot)
    //   3. selected      → an OUTER RING (navigation; independent of the above)
    let states: Vec<NodeVis> =
        FRIEZE.iter().map(|(p, _, _)| node_vis(hd, ws.current.is_some(), *p)).collect();
    let (running_idx, frac) = match hd {
        HdState::Running { current, progress, .. } => (
            current.and_then(|c| FRIEZE.iter().position(|(p, _, _)| *p == c)),
            progress.map(|(d, t)| if t > 0 { d as f32 / t as f32 } else { 0.0 }),
        ),
        _ => (None, None),
    };
    let selected_idx = ws.current.is_some().then(|| canonical_node(ws.layer));

    const R: f32 = 6.0;
    let dark = C::from_rgb(0x1a, 0x1a, 0x1a);

    // Dimension 1 — completeness line: grey base + copper up to the last DONE
    // node (stays full after the build; independent of selection).
    painter.line_segment(
        [egui::pos2(x0, y), egui::pos2(x1, y)],
        egui::Stroke::new(2.0, C::from_rgb(0x2e, 0x2e, 0x2e)),
    );
    if let Some(last_done) = states.iter().rposition(|s| matches!(s, NodeVis::Done { .. })) {
        painter.line_segment(
            [egui::pos2(x0, y), egui::pos2(node_x(last_done), y)],
            egui::Stroke::new(2.0, COPPER),
        );
    }

    // Muted copper for NON-selectable nodes (Tectonique / Relief): pure build
    // stages with no retained layer to view. They still fill / run (progression)
    // but read as "technical", carry no selection ring, and take no click.
    let muted_done = C::from_rgb(0x63, 0x4d, 0x33);
    for (i, (_, label, view_layer)) in FRIEZE.iter().enumerate() {
        let center = egui::pos2(node_x(i), y);
        let vis = states[i];
        let selectable = view_layer.is_some();
        let selected = selectable && selected_idx == Some(i);

        // Dimension 1/2 — the dot by completeness / running.
        match vis {
            NodeVis::Done { cached } => {
                painter.circle_filled(center, R, if selectable { COPPER } else { muted_done });
                if cached {
                    let ring = if selectable {
                        C::from_rgb(0xE9, 0xDC, 0xC4)
                    } else {
                        C::from_rgb(0x8a, 0x7a, 0x60)
                    };
                    painter.circle_stroke(center, R - 2.0, egui::Stroke::new(1.5, ring));
                }
            }
            NodeVis::Pending => {
                painter.circle_filled(center, R, C::from_rgb(0x20, 0x20, 0x20));
                painter.circle_stroke(
                    center,
                    R,
                    egui::Stroke::new(1.5, C::from_rgb(0x3f, 0x3f, 0x3f)),
                );
            }
            NodeVis::Running => {
                // Distinct from "done": a dark disc + an animated progress ARC.
                painter.circle_filled(center, R, dark);
                painter.circle_stroke(
                    center,
                    R + 2.0,
                    egui::Stroke::new(2.0, C::from_rgb(0x3a, 0x2e, 0x20)),
                );
                let start = -std::f32::consts::FRAC_PI_2;
                let (a0, a1) = if running_idx == Some(i) {
                    match frac {
                        Some(f) => (start, start + f.clamp(0.0, 1.0) * std::f32::consts::TAU),
                        None => (time * 2.4, time * 2.4 + 1.8), // waiter: spinning arc
                    }
                } else {
                    (start, start) // not the live node (shouldn't happen)
                };
                // Progression stays visible even on non-selectable nodes, just muted.
                let arc = if selectable { COPPER_BRIGHT } else { muted_done };
                stroke_arc(&painter, center, R + 2.0, a0, a1, egui::Stroke::new(2.5, arc));
            }
        }

        // Dimension 3 — selection: an outer ring (selectable nodes only).
        if selected {
            painter.circle_stroke(center, R + 4.5, egui::Stroke::new(2.0, COPPER_BRIGHT));
        }

        // Label: non-selectable nodes read muted across every state; selectable
        // nodes brighten on selection / running / completion.
        let lab_col = if !selectable {
            match vis {
                NodeVis::Running => C::from_rgb(0x8a, 0x74, 0x58),
                NodeVis::Done { .. } => C::from_rgb(0x74, 0x6a, 0x5a),
                NodeVis::Pending => C::from_rgb(0x4a, 0x44, 0x3a),
            }
        } else if selected {
            COPPER_BRIGHT
        } else {
            match vis {
                NodeVis::Running => COPPER_BRIGHT,
                NodeVis::Done { .. } => C::from_rgb(0xc4, 0xc4, 0xc4),
                NodeVis::Pending => C::from_rgb(0x5a, 0x5a, 0x5a),
            }
        };
        painter.text(
            egui::pos2(node_x(i), y + 18.0),
            egui::Align2::CENTER_TOP,
            *label,
            egui::FontId::proportional(9.5),
            lab_col,
        );

        // Click a node to view its layer — selectable nodes only (a result must
        // exist to show). Non-selectable nodes get no click / hover / pointer.
        if let (true, Some(layer)) = (ws.current.is_some(), *view_layer) {
            let hit = egui::Rect::from_center_size(center, egui::vec2((x1 - x0) / n as f32, 44.0));
            let resp = ui.interact(hit, ui.id().with(("frieze", i)), egui::Sense::click());
            if resp.hovered() {
                ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
            }
            if resp.clicked() {
                ws.layer = layer;
            }
        }
    }
}

/// Stroke a circular arc (egui has no arc primitive) as short segments.
fn stroke_arc(
    painter: &egui::Painter,
    center: egui::Pos2,
    radius: f32,
    a0: f32,
    a1: f32,
    stroke: egui::Stroke,
) {
    let steps = (((a1 - a0).abs() / 0.22).ceil() as usize).max(1);
    let mut prev = None;
    for k in 0..=steps {
        let t = a0 + (a1 - a0) * (k as f32 / steps as f32);
        let p = center + egui::vec2(t.cos(), t.sin()) * radius;
        if let Some(pp) = prev {
            painter.line_segment([pp, p], stroke);
        }
        prev = Some(p);
    }
}

/// Coarse continent → relief heightmap image (no FBM detail; the shape only).
fn preview_color_image(coarse: &ymir_core::grid::GridF32) -> egui::ColorImage {
    let (w, h) = (coarse.width, coarse.height);
    let mut rgba = vec![0u8; w * h * 4];
    for k in 0..w * h {
        let [r, g, b] = relief_color(coarse.data[k]);
        rgba[k * 4..k * 4 + 4].copy_from_slice(&[r, g, b, 255]);
    }
    // North-up view (Finding 27): mirror the south-first coarse field so the preview
    // matches the generated HD map (both show north at the top).
    flip_rows_rgba(&mut rgba, w, h);
    egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba)
}

/// Render the fast tectonic-shape preview: coarse continent + the export-window
/// box + a one-line island verdict, so a seed can be judged before the HD run.
/// THE single source of truth for the export scale (M1 #190). The domain IS the
/// map (no crop), so a pixel is `domain_km · 1000 / resolution` metres — every
/// scale readout (slider label, top bar, preview chip) derives from here, so they
/// cannot disagree by construction.
fn m_per_px(domain_km: f32, resolution: usize) -> f32 {
    domain_km * 1000.0 / resolution.max(1) as f32
}

/// Sea-level norm on the coarse field (sea = 0.5 in the C1 vertical contract).
const SEA_NORM: f32 = 0.5;

const OK_GREEN: C = C::from_rgb(0x6a, 0xb0, 0x6a);
const WARN_ORANGE: C = C::from_rgb(0xc0, 0x7a, 0x4a);

/// Cyclically shift a periodic coarse field so output cell `(x,y)` reads torus
/// cell `(x+rx, y+ry) mod (w,h)`. Matches the HD export's `sample_origin =
/// offset/grid`, so the preview frames exactly what ships. Pure relabelling — no
/// recompute; `rem_euclid` makes any integer offset valid (cyclic, no clamp).
fn roll_grid(src: &GridF32, rx: i64, ry: i64) -> GridF32 {
    let (w, h) = (src.width, src.height);
    let (wi, hi) = (w as i64, h as i64);
    let mut out = vec![0.0f32; w * h];
    for y in 0..h {
        let sy = ((y as i64 + ry).rem_euclid(hi)) as usize;
        for x in 0..w {
            let sx = ((x as i64 + rx).rem_euclid(wi)) as usize;
            out[y * w + x] = src.data[sy * w + sx];
        }
    }
    GridF32::from_vec(w, h, out)
}

/// The coarse preview field with the current framing offset applied (the frame the
/// user sees and the export ships). `None` until the first preview lands.
fn framed_coarse(ws: &WorkspaceState) -> Option<GridF32> {
    ws.preview.as_ref().map(|p| roll_grid(&p.coarse, ws.offset_cells[0], ws.offset_cells[1]))
}

/// Params for a coarse preview: NO manual offset (the preview computes the auto
/// offset; the UI applies the roll at render time). Never runs the HD pipeline.
fn preview_params(ws: &WorkspaceState) -> HdParams {
    HdParams {
        base_level_m: None,    // the preview never runs the HD pipeline
        base_level_off: false, // idem -- and `false` is production, not "no bound"
        target_size: ws.resolution,
        latitude_deg: ws.latitude,
        domain_km: ws.domain_km,
        manual_offset: None,
        stream_power: false, // preview is coarse-only; SP applies to the HD run
        closures: false,
        cross_rill: false,
        slope_floor_s_eq: None, // ADR Finding 109
        valley_construction: None, // ADR Finding 121
        cross_rill_d: 0.40,
        mfd: false,
        mfd_p: 2.0,
        geo_scale_ratio: 1.0,
        latitude_span_deg: None,
        export_dir: None,
        volcanism: None,
        lithology: None,
        fracture: None,
        infiltration: None,
        emit_tectonic_labels: false, // preview is coarse-only; overlays need the HD run
        geology_rules: None,
    }
}

/// Request a fresh CHEAP coarse preview for `ws.seed` (never the HD pipeline).
/// `reset_pan` = true for a seed change (snap the framing back to auto + drop the
/// stale HD result when the preview lands); false for a manual refresh of the same
/// seed (keep any manual pan).
fn request_preview(ws: &mut WorkspaceState, bridge: &C1SolverBridge, reset_pan: bool) {
    if reset_pan {
        ws.panned = false;
        ws.current = None; // a new shape invalidates the previous HD result
    }
    ws.seed_dirty = false;
    ws.bootstrapped = true;
    let spec = C1RunSpec { seed: ws.seed, ..C1RunSpec::default() };
    let _ = bridge.submit_preview(spec, preview_params(ws));
}

/// Reset the framing/scale for the CURRENT seed (Reset button): snap the offset to
/// auto, restore the default domain, clear the HD result. Seed + preview untouched.
fn reset_current_seed(ws: &mut WorkspaceState) {
    ws.offset_cells = ws.auto_offset_cells;
    ws.panned = false;
    ws.domain_km = 1024.0;
    ws.current = None;
    ws.texture = None;
    ws.tex_layer = None;
    ws.hover = None;
}

/// Manual framing controls: two cyclic cell offsets (x, y) with an auto/manual
/// indicator and a snap-to-auto button. Panning recomputes metrics instantly (pure
/// relabelling of the same coarse field — no tectonic recompute, no HD run).
fn framing_controls(ui: &mut egui::Ui, ws: &mut WorkspaceState) {
    let grid = ws.preview.as_ref().map(|p| p.coarse.width as i64).unwrap_or(64);
    field_label(ui, "CADRAGE (cellules, cyclique)");
    ui.add_space(1.0);
    ui.horizontal(|ui| {
        let (mut x, mut y) = (ws.offset_cells[0], ws.offset_cells[1]);
        let rx = ui.add(egui::DragValue::new(&mut x).prefix("x ").speed(0.15));
        let ry = ui.add(egui::DragValue::new(&mut y).prefix("y ").speed(0.15));
        if rx.changed() || ry.changed() {
            ws.offset_cells = [x.rem_euclid(grid), y.rem_euclid(grid)];
            ws.panned = ws.offset_cells != ws.auto_offset_cells;
        }
        let (tag, col) = if ws.panned { ("manuel", WARN_ORANGE) } else { ("auto", OK_GREEN) };
        ui.label(egui::RichText::new(tag).color(col).monospace().size(11.0));
        if ws.panned
            && ui
                .add(egui::Button::new(egui::RichText::new("↺ auto").color(BRONZE).size(10.5)))
                .on_hover_text("Revenir au cadrage automatique")
                .clicked()
        {
            ws.offset_cells = ws.auto_offset_cells;
            ws.panned = false;
        }
    });
}

/// Compact domain-as-map readout under the DOMAINE slider: m/px band flag, then
/// (once a coarse preview exists) the largest-mass traverse + bbox %, the per-side
/// ocean margins, and the geometric seed verdict. Recomputed from the coarse field
/// each frame via `domain_metrics` — a 64² pass, no tectonic recompute — so it
/// tracks the slider instantly. All figures scale with `ws.domain_km`.
fn domain_readout(ui: &mut egui::Ui, ws: &WorkspaceState, mpp: f32) {
    ui.add_space(6.0);
    let in_band = (30.0..=50.0).contains(&mpp);
    let (flag, fcol) = if in_band {
        ("dans la cible 30–50 m", OK_GREEN)
    } else {
        ("hors cible 30–50 m", WARN_ORANGE)
    };
    ui.label(
        egui::RichText::new(format!("{mpp:.1} m/px — {flag}")).color(fcol).monospace().size(10.5),
    );

    let Some(coarse) = framed_coarse(ws) else {
        ui.label(egui::RichText::new("· aperçu en cours…").color(DIM2).monospace().size(10.0));
        return;
    };
    let ss = SteinSteinParams::default();
    let m = domain_metrics(&coarse, SEA_NORM, &ss, ws.domain_km, ws.resolution);
    let traverse = m.extent_km.0.max(m.extent_km.1);
    ui.label(
        egui::RichText::new(format!(
            "masse max {traverse:.0} km · bbox {:.0}×{:.0}% du domaine",
            m.bbox_frac_x * 100.0,
            m.bbox_frac_y * 100.0,
        ))
        .color(DIM)
        .monospace()
        .size(10.5),
    );
    ui.label(
        egui::RichText::new(format!(
            "marges océan  N{:.0} S{:.0} E{:.0} W{:.0} km",
            m.margin_n_km, m.margin_s_km, m.margin_e_km, m.margin_w_km,
        ))
        .color(DIM)
        .monospace()
        .size(10.5),
    );
    let (vtxt, vcol) = if m.verdict_pass {
        ("verdict ✓ île centrée entourée d'océan".to_string(), OK_GREEN)
    } else {
        (format!("verdict ✗ {}", m.verdict_reason), WARN_ORANGE)
    };
    ui.label(egui::RichText::new(vtxt).color(vcol).monospace().size(10.5));
}

fn render_preview(ui: &mut egui::Ui, ws: &mut WorkspaceState, preview: &PreviewShape) {
    let grid = preview.coarse.width as i64;
    // Rolled (framed) coarse field: rebuild the texture whenever the offset changes
    // (a pan) so the drawn frame follows the offset. Metrics use the same field.
    let framed = roll_grid(&preview.coarse, ws.offset_cells[0], ws.offset_cells[1]);
    if ws.preview_tex.is_none() || ws.preview_tex_offset != ws.offset_cells {
        let img = preview_color_image(&framed);
        ws.preview_tex =
            Some(ui.ctx().load_texture("tecto_preview", img, egui::TextureOptions::NEAREST));
        ws.preview_tex_offset = ws.offset_cells;
    }
    let handle = ws.preview_tex.as_ref().unwrap().clone();
    let ss = SteinSteinParams::default();
    let m = domain_metrics(&framed, SEA_NORM, &ss, ws.domain_km, ws.resolution);
    let mpp = m_per_px(ws.domain_km, ws.resolution);

    // Fixed toolbar: zoom controls + the verdict/metrics line (stays put while the
    // magnified image scrolls beneath it).
    let (vtag, vcol) = if m.verdict_pass {
        ("île centrée ✓", OK_GREEN)
    } else {
        ("non conforme ✗", WARN_ORANGE)
    };
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("Zoom").color(DIM2).size(11.0));
        if ui.small_button("−").clicked() {
            ws.preview_zoom = (ws.preview_zoom / 1.25).max(1.0);
        }
        ui.label(egui::RichText::new(format!("{:.0}%", ws.preview_zoom * 100.0)).color(COPPER_BRIGHT).monospace().size(11.0));
        if ui.small_button("+").clicked() {
            ws.preview_zoom = (ws.preview_zoom * 1.25).min(12.0);
        }
        if ui.small_button("⟲").clicked() {
            ws.preview_zoom = 1.0;
        }
        ui.label(
            egui::RichText::new(format!(
                "masse {:.0} km · bbox {:.0}×{:.0}% · N{:.0} S{:.0} E{:.0} W{:.0} km · {:.0} m/px · {vtag}",
                m.extent_km.0.max(m.extent_km.1),
                m.bbox_frac_x * 100.0,
                m.bbox_frac_y * 100.0,
                m.margin_n_km,
                m.margin_s_km,
                m.margin_e_km,
                m.margin_w_km,
                mpp,
            ))
            .color(vcol)
            .monospace()
            .size(11.0),
        );
    });

    // Scroll bars ("sliders") appear once the zoomed image exceeds the viewport.
    // `drag_to_scroll(false)` keeps content-drag for RECADRAGE; the scroll bars pan
    // the view. Wheel over the image zooms.
    let avail = ui.available_size();
    let side = avail.x.min(avail.y).max(1.0);
    let img_side = side * ws.preview_zoom;
    egui::ScrollArea::both().auto_shrink([false, false]).drag_to_scroll(false).show(ui, |ui| {
        let resp = ui
            .add(
                egui::Image::new(egui::load::SizedTexture::new(handle.id(), egui::vec2(img_side, img_side)))
                    .sense(egui::Sense::drag()),
            )
            .on_hover_text("Glisser pour recadrer (rotation cyclique du tore) · molette pour zoomer · barres pour défiler");
        let r = resp.rect;

        // Wheel-to-zoom when hovering the image.
        if resp.hovered() {
            let scroll = ui.input(|i| i.raw_scroll_delta.y);
            if scroll.abs() > 0.0 {
                ws.preview_zoom = (ws.preview_zoom * (1.0 + scroll * 0.001)).clamp(1.0, 12.0);
            }
        }

        // Drag-to-recadrage: shift the framing by whole coarse cells (cyclic). Uses
        // the ZOOMED image size for pixels-per-cell so the feel is scale-correct.
        if resp.dragged() {
            let ppc = (img_side / grid as f32).max(1.0);
            let d = resp.drag_delta();
            ws.pan_accum[0] += d.x;
            ws.pan_accum[1] += d.y;
            let step_x = (ws.pan_accum[0] / ppc).trunc() as i64;
            let step_y = (ws.pan_accum[1] / ppc).trunc() as i64;
            if step_x != 0 || step_y != 0 {
                ws.pan_accum[0] -= step_x as f32 * ppc;
                ws.pan_accum[1] -= step_y as f32 * ppc;
                ws.offset_cells[0] = (ws.offset_cells[0] - step_x).rem_euclid(grid);
                // North-up view (Finding 27): the preview is mirrored vertically, so a
                // downward drag must roll the frame the opposite way to follow the cursor.
                ws.offset_cells[1] = (ws.offset_cells[1] + step_y).rem_euclid(grid);
                ws.panned = ws.offset_cells != ws.auto_offset_cells;
            }
        }
        if resp.drag_stopped() {
            ws.pan_accum = [0.0, 0.0];
        }

        // Largest-mass bounding box (map frame, from the ocean margins), drawn on the
        // image so it scrolls/zooms with it. Row 0 (y=0 = south) is at the top.
        let dk = ws.domain_km.max(1.0);
        let box_rect = egui::Rect::from_min_size(
            egui::pos2(
                r.min.x + (m.margin_w_km / dk) * r.width(),
                r.min.y + (m.margin_s_km / dk) * r.height(),
            ),
            egui::vec2(m.bbox_frac_x * r.width(), m.bbox_frac_y * r.height()),
        );
        ui.painter_at(r).rect_stroke(
            box_rect,
            2.0,
            egui::Stroke::new(2.0, COPPER_BRIGHT),
            egui::StrokeKind::Middle,
        );
    });
}

fn map(ui: &mut egui::Ui, ws: &mut WorkspaceState) {
    if ws.current.is_none() {
        if let Some(preview) = ws.preview.clone() {
            render_preview(ui, ws, &preview);
        } else {
            ui.centered_and_justified(|ui| {
                ui.label(
                    egui::RichText::new(
                        "Aucun continent — « Aperçu tectonique » pour voir la forme, puis « Générer ».",
                    )
                    .color(DIM),
                );
            });
        }
        return;
    }
    // (Re)build the texture on layer/result/overlay change.
    if ws.texture.is_none()
        || ws.tex_layer != Some(ws.layer)
        || ws.tex_overlay != ws.river_overlay
        || ws.tex_overlays != ws.overlays
    {
        let layer = ws.layer;
        let overlay = ws.river_overlay;
        let overlays = ws.overlays;
        let t_tex = std::time::Instant::now();
        let (img, stats) = {
            let hd = ws.current.as_ref().unwrap();
            let rm = ws.river_map.as_ref().unwrap();
            // ADR Finding 135-V -- the relief view
            let relief = match (ws.relief_view, &ws.diff_ref) {
                (1, _) => ReliefView::Shade { outline: ws.lake_outline },
                (2, Some(r)) => ReliefView::Diff(r, DIFF_SAT_M[ws.diff_sat.min(2)]),
                _ => ReliefView::Hypso,
            };
            let geo = hd.geology.as_ref().map(|g| (ws.zoning_override.as_ref().map_or(&g.zoning, |o| &o.0), ws.chance_only, g.product.rocks.as_slice(), ws.chance_overlay));
            layer_color_image(hd, layer, rm, overlay, overlays, relief, geo)
        };
        ws.diff_stats = stats;
        ws.texture = Some(ui.ctx().load_texture("hd_map", img, egui::TextureOptions::NEAREST));
        ws.tex_layer = Some(ws.layer);
        ws.tex_overlay = ws.river_overlay;
        ws.tex_overlays = ws.overlays;
        eprintln!(
            "[HD timing] map texture rebuild {:.1}s (UI thread)",
            t_tex.elapsed().as_secs_f32()
        );
    }
    let handle = ws.texture.as_ref().unwrap().clone();

    // ── Camera (item 11): an ABSOLUTE on-screen scale `map_px` (the full-map size in px) + a uv
    //    pan centre. The scale is FIXED (not derived from the panel-dependent available size), so
    //    resizing the side panels does NOT rescale the canvas. The map is painted freely and
    //    clipped to the viewport (fills it); drag pans, CTRL+wheel zooms on the cursor, wheel
    //    scrolls, SHIFT+wheel scrolls horizontally.
    let hd = ws.current.clone().unwrap();
    let (w, h) = (hd.width, hd.height);
    let vp = ui.available_rect_before_wrap();
    let resp = ui.interact(vp, ui.id().with("map_canvas"), egui::Sense::click_and_drag());
    let fit = vp.width().min(vp.height()).max(64.0);
    if ws.map_px <= 0.0 {
        ws.map_px = fit; // fit on the first frame
    }
    ws.map_pan[0] = ws.map_pan[0].clamp(0.0, 1.0);
    ws.map_pan[1] = ws.map_pan[1].clamp(0.0, 1.0);

    // ── Navigation input. ──
    let pan_drag =
        (ws.tool == MapTool::Pan && resp.dragged()) || resp.dragged_by(egui::PointerButton::Middle);
    if pan_drag {
        let d = resp.drag_delta();
        ws.map_pan[0] -= d.x / ws.map_px; // drag right -> view slides left
        ws.map_pan[1] -= d.y / ws.map_px;
    }
    if resp.hovered() {
        let (scroll, zoom_delta, mods) =
            ui.ctx().input(|i| (i.smooth_scroll_delta, i.zoom_delta(), i.modifiers));
        let zf = if (zoom_delta - 1.0).abs() > 1e-3 {
            zoom_delta
        } else if (mods.ctrl || mods.command) && scroll.y != 0.0 {
            1.0 + scroll.y / 400.0
        } else {
            1.0
        };
        if (zf - 1.0).abs() > 1e-4 {
            let new_px = (ws.map_px * zf).clamp(fit * 0.5, fit * 16.0);
            if let Some(pos) = resp.hover_pos() {
                // keep the uv under the cursor fixed while zooming.
                let map_min = vp.center() - egui::vec2(ws.map_pan[0], ws.map_pan[1]) * ws.map_px;
                let cur_u = (pos.x - map_min.x) / ws.map_px;
                let cur_v = (pos.y - map_min.y) / ws.map_px;
                ws.map_pan[0] = cur_u + (vp.center().x - pos.x) / new_px;
                ws.map_pan[1] = cur_v + (vp.center().y - pos.y) / new_px;
            }
            ws.map_px = new_px;
        } else if !mods.ctrl && !mods.command && (scroll.x != 0.0 || scroll.y != 0.0) {
            ws.map_pan[0] -= scroll.x / ws.map_px;
            ws.map_pan[1] -= scroll.y / ws.map_px;
        }
        let icon = if pan_drag {
            egui::CursorIcon::Grabbing
        } else if ws.tool == MapTool::Pan {
            egui::CursorIcon::Grab
        } else {
            egui::CursorIcon::Crosshair
        };
        ui.ctx().set_cursor_icon(icon);
    }
    // Consume a pending microscope double-click "frame element": centre on it + fit its extent
    // (with a margin), clamped to the zoom range.
    if let Some((center, ext)) = ws.focus_target.take() {
        ws.map_pan = center;
        let target = (vp.width().min(vp.height()) * 0.8) / ext.max(1e-4);
        ws.map_px = target.clamp(fit * 0.5, fit * 16.0);
    }
    ws.map_pan[0] = ws.map_pan[0].clamp(0.0, 1.0);
    ws.map_pan[1] = ws.map_pan[1].clamp(0.0, 1.0);

    // ── Transform for this frame + paint the (clipped) map. ──
    let map_px = ws.map_px;
    let map_min = vp.center() - egui::vec2(ws.map_pan[0], ws.map_pan[1]) * map_px;
    let map_rect = egui::Rect::from_min_size(map_min, egui::vec2(map_px, map_px));
    let to_screen = |x: u32, y: u32| -> egui::Pos2 {
        let ncx = (x as f32 + 0.5) / w as f32;
        let ncy = ((h as u32 - 1 - y) as f32 + 0.5) / h as f32;
        map_min + egui::vec2(ncx, ncy) * map_px
    };
    let pnt = ui.painter_at(vp);
    pnt.image(
        handle.id(),
        map_rect,
        egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
        C::WHITE,
    );

    overlay_chip(ui, vp, ws.layer);
    // ADR Finding 152 -- the geology views' legends
    let geo_items: Vec<(C, String, String)> = match (ws.layer, ws.current.as_ref().and_then(|h| h.geology.as_ref())) {
        (HdLayer::Rocks, Some(_)) => ymir_core::geology::rocks::ROCK_CLASSES[1..]
            .iter()
            .map(|r| (C::from_rgb(r.color[0], r.color[1], r.color[2]), r.name_fr.to_string(), ["tendre", "moyenne", "dure"][r.hardness.min(2) as usize].to_string()))
            .collect(),
        _ => Vec::new(),
    };
    legend_box(ui, vp, ws.layer, (ws.relief_view, ws.diff_sat, ws.diff_stats, ws.diff_ref.is_some()), geo_items);
    // ADR Finding 153 -- the « Chances » mask's own legend (top left of the map)
    if ws.chance_overlay
        && let Some(g) = ws.current.as_ref().and_then(|h| h.geology.as_ref())
    {
        let z = ws.zoning_override.as_ref().map_or(&g.zoning, |o| &o.0);
        let items: Vec<(C, String, String)> = z
            .resources
            .iter()
            .enumerate()
            .filter(|(i, _)| ws.chance_only.is_none_or(|o| o == *i))
            .map(|(_, r)| (C::from_rgb(r.color[0], r.color[1], r.color[2]), r.name_fr.clone(), format!("max {}", r.fav.iter().max().copied().unwrap_or(0))))
            .collect();
        chances_legend(ui, vp, &items, ws.chance_only.is_none());
    }

    // Hover -> cell -> inspect (only when the cursor is over the map, not the letterbox).
    ws.hover = None;
    ws.hover_xy = None;
    if let Some(pos) = resp.hover_pos() {
        let fx = (pos.x - map_min.x) / map_px;
        let fy = (pos.y - map_min.y) / map_px;
        if (0.0..1.0).contains(&fx) && (0.0..1.0).contains(&fy) {
            let cx = ((fx * w as f32) as usize).min(w - 1);
            // The texture is north-up mirrored (Finding 27); the DATA cell is the mirror row.
            let tex_row = ((fy * h as f32) as usize).min(h - 1);
            let cy = (h - 1) - tex_row;
            if let Some(rm) = ws.river_map.as_ref() {
                ws.hover = Some(inspect_cell(&hd, rm, cx, cy));
                ws.hover_xy = Some((cx, cy));
            }
            // Hover-select only with the SELECT tool (the hand pans without touching selection).
            if ws.tool == MapTool::Select {
                match ws.inspect_tab {
                    InspectTab::Rivers => {
                        if let Some(info) = ws.river_map.as_ref().and_then(|rm| rm.at(cx, cy)) {
                            if let Some(&wi) = ws.seg_to_wc.get(info.segment) {
                                if wi != usize::MAX {
                                    ws.selected_river = Some(wi);
                                }
                            }
                        }
                    }
                    InspectTab::Lakes => {
                        let id = hd.drainage.lake_map[cy * hd.width + cx];
                        if id != 0 {
                            if let Some(li) = hd.drainage.lakes.iter().position(|l| l.base.id == id)
                            {
                                ws.selected_lake = Some(li);
                            }
                        }
                    }
                    // Volcanoes are picked from the microscope list, not the map hover.
                    InspectTab::Volcanoes => {}
                }
            }
            // Reticle: faint crosshair through the cell + a copper cell box.
            let s = to_screen(cx as u32, cy as u32);
            let cell = (map_px / w as f32).max(8.0);
            let faint = C::from_rgba_unmultiplied(0xC9, 0x85, 0x3F, 90);
            pnt.line_segment(
                [egui::pos2(s.x, vp.top()), egui::pos2(s.x, vp.bottom())],
                egui::Stroke::new(1.0, faint),
            );
            pnt.line_segment(
                [egui::pos2(vp.left(), s.y), egui::pos2(vp.right(), s.y)],
                egui::Stroke::new(1.0, faint),
            );
            pnt.rect_stroke(
                egui::Rect::from_center_size(s, egui::vec2(cell, cell)),
                0.0,
                egui::Stroke::new(1.5, COPPER_BRIGHT),
                egui::StrokeKind::Middle,
            );
            // Hover coord readout (top-right of the viewport).
            let galley = pnt.layout_no_wrap(
                format!("x {cx}  y {cy}"),
                egui::FontId::monospace(10.5),
                C::from_rgb(0x9a, 0x9a, 0x9a),
            );
            let bg = egui::Rect::from_min_size(
                egui::pos2(vp.right() - 8.0 - galley.size().x - 12.0, vp.top() + 8.0),
                galley.size() + egui::vec2(12.0, 8.0),
            );
            pnt.rect_filled(bg, 6.0, C::from_rgba_unmultiplied(18, 18, 18, 210));
            pnt.galley(egui::pos2(bg.left() + 6.0, bg.top() + 4.0), galley, TEXT);
        }
    }

    // Selection highlight (Finding 28) — painter-drawn.
    if ws.inspect_tab == InspectTab::Rivers {
        if let Some(wc) =
            ws.selected_river.and_then(|i| ws.watercourses.as_ref().and_then(|v| v.get(i)))
        {
            // Main TRUNK in orange, the secondary sources (tributaries) in yellow. Draw tributaries
            // first (underneath), then the trunk on top so the main stem reads clearly.
            let orange = C::from_rgb(255, 150, 40);
            let yellow = C::from_rgb(255, 220, 90);
            let trunk: std::collections::HashSet<usize> = wc.trunk.iter().copied().collect();
            let draw = |pnt: &egui::Painter, s: usize, col: C, wdt: f32| {
                let pts: Vec<egui::Pos2> = hd.drainage.rivers.segments[s]
                    .points
                    .iter()
                    .map(|&(px, py)| to_screen(px, py))
                    .collect();
                if pts.len() >= 2 {
                    pnt.add(egui::Shape::line(pts, egui::Stroke::new(wdt, col)));
                }
            };
            for &s in &wc.segments {
                if !trunk.contains(&s) {
                    draw(&pnt, s, yellow, 1.6);
                }
            }
            for &s in &wc.trunk {
                draw(&pnt, s, orange, 2.4);
            }
            pnt.circle_filled(to_screen(wc.mouth_xy.0, wc.mouth_xy.1), 4.0, sink_label(wc.sink).1);
            // Profile-hover tracking point: mirror the long-profile cursor onto the river highlight
            // at the same length-from-source (Finding 29 UI).
            if let Some(fx) = ws.profile_hover {
                let np = ws.profile_path.len();
                if np >= 2 {
                    let i = ((fx * (np - 1) as f32).round() as usize).min(np - 1);
                    let (px, py) = ws.profile_path[i];
                    let sp = to_screen(px, py);
                    pnt.circle_filled(sp, 5.0, C::WHITE);
                    pnt.circle_stroke(sp, 5.0, egui::Stroke::new(1.5, COPPER_BRIGHT));
                }
            }
        }
    } else if let Some(lk) = ws.selected_lake.and_then(|i| hd.drainage.lakes.get(i)) {
        let id = lk.base.id;
        let cyan = C::from_rgba_unmultiplied(120, 230, 240, 220);
        for y in 0..h {
            for x in 0..w {
                if hd.drainage.lake_map[y * w + x] == id {
                    let edge = [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)].iter().any(|&(dx, dy)| {
                        let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                        nx < 0
                            || ny < 0
                            || nx as usize >= w
                            || ny as usize >= h
                            || hd.drainage.lake_map[ny as usize * w + nx as usize] != id
                    });
                    if edge {
                        pnt.circle_filled(to_screen(x as u32, y as u32), 1.5, cyan);
                    }
                }
            }
        }
    }

    // ADR Finding 146 -- the "Rivières LL" layer: the smoothed main stems passing the filters, stroke ∝ the bankfull width
    if ws.rll_overlay {
        let to_screen_f = |x: f32, y: f32| -> egui::Pos2 {
            map_min + egui::vec2(x / w as f32, (h as f32 - y) / h as f32) * map_px
        };
        for r in hd.rivers_ll.iter().filter(|r| rll_keep(ws, r)) {
            let pts: Vec<egui::Pos2> = r.points.iter().map(|p| to_screen_f(p[0], p[1])).collect();
            if pts.len() >= 2 {
                let wdt = (0.6 + 0.35 * r.width_m.max(1.0).log2()).clamp(0.8, 3.5);
                pnt.add(egui::Shape::line(pts, egui::Stroke::new(wdt, C::from_rgb(0x4a, 0xa8, 0xf0))));
            }
        }
    }

    // C-2 volcano markers — shown when the "Symboles" toggle is on (default). A
    // triangle at each crater: red = active/degassing, grey = extinct. The selected
    // one is enlarged with a white ring. A selected volcano stays ringed even with
    // symbols off, so the microscope selection is never lost.
    for (i, vo) in hd
        .volcanoes
        .iter()
        .enumerate()
        .filter(|(i, _)| ws.show_symbols || ws.selected_volcano == Some(*i))
    {
        let c = to_screen(vo.center_px.0 as u32, vo.center_px.1 as u32);
        let sel = ws.selected_volcano == Some(i);
        let s = if sel { 9.0 } else { 6.0 };
        let fill =
            if vo.active { C::from_rgb(0xe0, 0x50, 0x30) } else { C::from_rgb(0x8a, 0x8a, 0x8a) };
        let tri = vec![
            egui::pos2(c.x, c.y - s),
            egui::pos2(c.x - s * 0.9, c.y + s * 0.7),
            egui::pos2(c.x + s * 0.9, c.y + s * 0.7),
        ];
        pnt.add(egui::Shape::convex_polygon(
            tri,
            fill,
            egui::Stroke::new(1.5, C::from_rgb(0x20, 0x14, 0x10)),
        ));
        if sel {
            pnt.circle_stroke(c, s + 3.0, egui::Stroke::new(2.0, C::WHITE));
        }
    }

    // ADR Finding 140 -- the gorge's tagged falls: a cyan diamond on the Relief layer (with "Symboles"); hovering
    // shows its kind, height and lake
    if ws.show_symbols && ws.layer == HdLayer::Relief {
        let hover = ui.input(|i| i.pointer.hover_pos());
        for (f, lk) in &hd.gorge_falls {
            let c = to_screen(f.x, f.y);
            let s = 6.0;
            let dia = vec![
                egui::pos2(c.x, c.y - s),
                egui::pos2(c.x + s, c.y),
                egui::pos2(c.x, c.y + s),
                egui::pos2(c.x - s, c.y),
            ];
            pnt.add(egui::Shape::convex_polygon(
                dia,
                C::from_rgb(0x40, 0xd0, 0xf0),
                egui::Stroke::new(1.5, C::from_rgb(0x10, 0x20, 0x30)),
            ));
            if hover.is_some_and(|h| h.distance(c) < 9.0) {
                let kind = match f.kind {
                    ymir_core::tectonics_c1::valley_construction::GorgeFallKind::WaterfallHead => "chute de tête",
                    ymir_core::tectonics_c1::valley_construction::GorgeFallKind::WaterfallShortage => {
                        "chute de manque de place"
                    }
                };
                let lake = lk.map_or("—".to_string(), |l| l.to_string());
                pnt.text(
                    egui::pos2(c.x + 10.0, c.y - 8.0),
                    egui::Align2::LEFT_TOP,
                    format!("{kind} · {:.0} m · lac {lake}", f.height_m),
                    egui::FontId::proportional(11.0),
                    C::WHITE,
                );
            }
        }
    }

    // Floating tool bar (top-left), minimap (bottom-left), zoom controls (bottom-right).
    canvas_toolbar(ui, vp, ws);
    if ws.minimap {
        let uv = egui::Rect::from_min_max(
            egui::pos2((vp.left() - map_min.x) / map_px, (vp.top() - map_min.y) / map_px),
            egui::pos2((vp.right() - map_min.x) / map_px, (vp.bottom() - map_min.y) / map_px),
        );
        minimap_overlay(ui, vp, ws, &hd, uv);
    }
    zoom_controls(ui, vp, ws);
    if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
        ws.selected_river = None;
        ws.selected_lake = None;
        ws.selected_volcano = None;
    }
}

/// ADR Finding 146 -- does a Living Landz river pass the layer's selection filters? A spillway has no Strahler
/// order and passes that filter.
fn rll_keep(ws: &WorkspaceState, r: &ymir_core::export::rivers_ll::RiverLl) -> bool {
    r.length_km >= ws.rll_min_len_km
        && r.catchment_km2 >= ws.rll_min_area_km2
        && r.strahler.is_none_or(|s| s >= ws.rll_min_strahler)
        && !(ws.rll_hide_fusion && (r.fusion_candidate || r.micro_lake_inflow))
        && !(ws.rll_hide_retrace && r.retraces_spillway)
}

/// Floating canvas tool bar (mock): SELECT (microscope) / PAN (hand) tools, plus the hydro
/// overlay and minimap toggles. Painted as a small pill at the canvas top-left.
fn canvas_toolbar(ui: &mut egui::Ui, rect: egui::Rect, ws: &mut WorkspaceState) {
    let bar = egui::Rect::from_min_size(
        egui::pos2(rect.left() + 10.0, rect.top() + 10.0),
        egui::vec2(0.0, 0.0),
    );
    let mut child = ui.new_child(
        egui::UiBuilder::new()
            .max_rect(egui::Rect::from_min_size(bar.min, egui::vec2(520.0, 34.0))),
    );
    child.horizontal(|ui| {
        egui::Frame::default()
            .fill(C::from_rgba_unmultiplied(18, 18, 18, 235))
            .stroke(egui::Stroke::new(1.0, C::from_rgb(0x2e, 0x2e, 0x2e)))
            .corner_radius(6.0)
            .inner_margin(egui::Margin::symmetric(6, 4))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    let sel = seg_label(ui, "⌖ Sélection", ws.tool == MapTool::Select);
                    if sel.clicked() {
                        ws.tool = MapTool::Select;
                    }
                    sel.on_hover_text("Survol/clic inspecte une entité (microscope).");
                    let pan = seg_label(ui, "✋ Main", ws.tool == MapTool::Pan);
                    if pan.clicked() {
                        ws.tool = MapTool::Pan;
                    }
                    pan.on_hover_text("Clic-glisser déplace la vue (pan).");
                    ui.separator();
                    if ui
                        .checkbox(&mut ws.river_overlay, egui::RichText::new("Réseau").size(11.0))
                        .on_hover_text("Réseau hydro (épaisseur ∝ Strahler ; orphelins en rouge)")
                        .changed()
                    {
                        ws.texture = None;
                    }
                    ui.checkbox(&mut ws.minimap, egui::RichText::new("Minimap").size(11.0));
                    // ADR Finding 146 -- the rivers for Living Landz, with their selection filters
                    ui.checkbox(&mut ws.rll_overlay, egui::RichText::new("Rivières LL").size(11.0))
                        .on_hover_text("Rivières pour Living Landz : cours principaux lissés dans leur vallée (rivers_ll.json)");
                    if ws.rll_overlay {
                        tool_menu(ui, egui::RichText::new("filtres ▾").size(11.0), |ui| {
                            ui.set_min_width(300.0);
                            ui.add(egui::Slider::new(&mut ws.rll_min_len_km, 0.0..=100.0).logarithmic(true).text("longueur min (km)"));
                            ui.add(egui::Slider::new(&mut ws.rll_min_area_km2, 0.0..=5000.0).logarithmic(true).text("aire min (km²)"));
                            ui.add(egui::Slider::new(&mut ws.rll_min_strahler, 1..=8).text("Strahler min"));
                            // ADR Findings 147-F / 149-T -- the annotations, with the counts and km each box hides
                            let count_km = |f: &dyn Fn(&ymir_core::export::rivers_ll::RiverLl) -> bool| -> (usize, f32) {
                                ws.current.as_ref().map_or((0, 0.0), |hd| {
                                    let v: Vec<_> = hd.rivers_ll.iter().filter(|r| f(r)).collect();
                                    (v.len(), v.iter().map(|r| r.length_km).sum())
                                })
                            };
                            let (n_mask, km_mask) = count_km(&|r| r.fusion_candidate || r.micro_lake_inflow);
                            let (n_retrace, km_retrace) = count_km(&|r| r.retraces_spillway);
                            ui.checkbox(&mut ws.rll_hide_fusion, format!("Masquer faisceaux et micro-rivières ({n_mask} · {km_mask:.0} km)"))
                                .on_hover_text(
                                    "Faisceaux : la rivière de plus faible débit d'une paire parallèle (à ≤ 2 cellules sur > 2 km, sans confluence commune). \
                                     Micro-rivières : moins de 1 km et finissant dans un lac. Affichage seul : le monde et les exports sont inchangés.",
                                );
                            ui.checkbox(&mut ws.rll_hide_retrace, format!("masquer les déversoirs retracés ({n_retrace} · {km_retrace:.0} km)"))
                                .on_hover_text("Un déversoir dont ≥ 50 % du tracé repasse sur des cours d'eau (dessiné sur leurs sommets)");
                            if let Some(hd) = ws.current.as_ref() {
                                let kept = hd.rivers_ll.iter().filter(|r| rll_keep(ws, r)).count();
                                let km: f32 = hd.rivers_ll.iter().filter(|r| rll_keep(ws, r)).map(|r| r.length_km).sum();
                                ui.label(format!("{kept} rivières retenues sur {} · {km:.0} km", hd.rivers_ll.len()));
                            }
                        });
                    }
                    ui.checkbox(&mut ws.show_symbols, egui::RichText::new("Symboles").size(11.0))
                        .on_hover_text(
                            "Marqueurs sur la carte (volcans : ▲ rouge actif / gris éteint)",
                        );
                    // ADR Finding 153 -- the « Chances » mask (over any step): the resource, the rules file and its reload. The rocks
                    // are the frieze's « Roches » step.
                    if ws.current.as_ref().is_some_and(|h| h.geology.is_some()) {
                        ui.separator();
                        if ui.checkbox(&mut ws.chance_overlay, egui::RichText::new("Chances").size(11.0)).on_hover_text("Favorabilité des ressources par-dessus l'étape affichée (opacité = favorabilité)").changed() {
                            ws.texture = None;
                        }
                        tool_menu(ui, egui::RichText::new("⛏ ressources ▾").size(11.0), |ui| {
                            ui.set_min_width(320.0);
                            let g = ws.current.as_ref().and_then(|h| h.geology.clone());
                            if let Some(g) = g {
                                if ui.checkbox(&mut ws.chance_overlay, "Afficher les chances").changed() {
                                    ws.texture = None;
                                }
                                ui.separator();
                                let names: Vec<(String, [u8; 3])> = ws.zoning_override.as_ref().map_or(&g.zoning, |o| &o.0).resources.iter().map(|r| (r.name_fr.clone(), r.color)).collect();
                                if ui.radio(ws.chance_only.is_none(), "Toutes (la plus forte l'emporte)").clicked() {
                                    ws.chance_only = None;
                                    ws.texture = None;
                                }
                                for (i, (n, col)) in names.iter().enumerate() {
                                    let txt = egui::RichText::new(format!("■ {n}")).color(C::from_rgb(col[0], col[1], col[2]));
                                    if ui.radio(ws.chance_only == Some(i), txt).clicked() {
                                        ws.chance_only = Some(i);
                                        ws.texture = None;
                                    }
                                }
                                ui.separator();
                                ui.horizontal(|ui| {
                                    ui.label("fichier");
                                    ui.text_edit_singleline(&mut ws.geology_rules_path);
                                });
                                if ui.button("Recharger les règles").clicked() {
                                    let loaded = ymir_core::geology::rules::load_rules(Some(std::path::Path::new(&ws.geology_rules_path)));
                                    let (z, secs) = ymir_core::geology::zone_timed(&g.product, &loaded.rules);
                                    ws.zoning_override = Some((z, loaded, secs));
                                    ws.texture = None;
                                }
                                let (label, sha, secs, err) = match &ws.zoning_override {
                                    Some((_, l, t)) => (l.label.clone(), l.sha256.clone(), *t, l.error.clone()),
                                    None => (g.rules.label.clone(), g.rules.sha256.clone(), g.zoning_s, g.rules.error.clone()),
                                };
                                ui.label(egui::RichText::new(format!("règles : {label} · {} · re-zonage {secs:.2} s", &sha[..12])).weak());
                                if let Some(e) = err {
                                    ui.label(egui::RichText::new(format!("⚠ fichier refusé, règles intégrées utilisées : {e}")).color(WARN_ORANGE));
                                }
                                ui.label(egui::RichText::new(format!("étage géologie {:.2} s par monde", g.stage_s)).weak());
                            }
                        });
                    }
                    // ADR Finding 135-V -- the relief view (a view of the world, never the world)
                    if ws.layer == HdLayer::Relief && ws.current.is_some() {
                        ui.separator();
                        let name = ["Hypsométrie", "Ombrage", "Différence"][ws.relief_view.min(2)];
                        tool_menu(ui, egui::RichText::new(format!("🗻 {name} ▾")).size(11.0), |ui| {
                            ui.set_min_width(280.0);
                            let views = [
                                ("Hypsométrie", "La couche livrée : l'altitude en couleur."),
                                (
                                    "Ombrage",
                                    "Ombrage du relief (lumière azimut 315°, hauteur 45°) calculé sur le \
                                     champ conditionné en mètres, multiplié à la couleur hypsométrique.",
                                ),
                                (
                                    "Différence",
                                    "Monde courant − référence, en mètres, sur le champ conditionné : \
                                     rouge = le courant est plus haut, bleu = plus bas, blanc = égal. \
                                     Échelle linéaire, SATURÉE au-delà de ± la valeur choisie.",
                                ),
                            ];
                            for (i, (lbl, hint)) in views.iter().enumerate() {
                                if ui.radio(ws.relief_view == i, *lbl).on_hover_text(*hint).clicked() {
                                    ws.relief_view = i;
                                    ws.texture = None;
                                }
                            }
                            // ADR Finding 157-V -- the lakes' outline on the hillshade
                            if ws.relief_view == 1
                                && ui
                                    .checkbox(&mut ws.lake_outline, "Contour des lacs")
                                    .on_hover_text(
                                        "Le contour des lacs du monde (le masque final de run_hd) sur l'ombrage : \
                                         un trait cyan côté lac doublé d'un trait sombre côté berge, lisible sur \
                                         les versants clairs comme sombres. Une vue : le monde ne change pas.",
                                    )
                                    .changed()
                            {
                                ws.texture = None;
                            }
                            ui.separator();
                            if ui
                                .button("Mémoriser ce monde comme référence")
                                .on_hover_text(
                                    "La vue Différence affiche (monde courant − référence). Mémorisez un \
                                     état (p. ex. le toggle OFF), puis générez l'autre.",
                                )
                                .clicked()
                            {
                                ws.diff_ref = ws.current.clone();
                                ws.texture = None;
                            }
                            let txt = match (&ws.diff_ref, &ws.current) {
                                (Some(r), Some(c)) if Arc::ptr_eq(r, c) => {
                                    "Référence = le monde courant (Δz ≡ 0)".to_string()
                                }
                                (Some(r), _) => format!("Référence mémorisée · {}", guard_short(&r.bench_guard)),
                                (None, _) => "Aucune référence mémorisée".to_string(),
                            };
                            ui.label(egui::RichText::new(txt).size(10.5).color(DIM));
                            ui.label(egui::RichText::new("Saturation").size(10.5).color(DIM2));
                            if let Some(i) = seg_row(ui, &["±10 m", "±100 m", "±1000 m"], ws.diff_sat) {
                                ws.diff_sat = i;
                                ws.texture = None;
                            }
                        });
                    }
                    // ── Microscope tectonique — un MENU compact (toggles empilés) pour
                    //    superposer les couches causales sur la vue courante. Disponible
                    //    quand les labels ont été dérivés (run HD ; pas l'aperçu coarse).
                    let has_labels =
                        ws.current.as_ref().and_then(|h| h.tectonic.as_ref()).is_some();
                    if has_labels {
                        ui.separator();
                        let n_on = [
                            ws.overlays.rift,
                            ws.overlays.subduction,
                            ws.overlays.collision,
                            ws.overlays.craton,
                            ws.overlays.divergent,
                            ws.overlays.lithology,
                        ]
                        .iter()
                        .filter(|b| **b)
                        .count();
                        let btn = if n_on > 0 {
                            format!("🔬 Tectonique ({n_on}) ▾")
                        } else {
                            "🔬 Tectonique ▾".to_string()
                        };
                        tool_menu(ui, egui::RichText::new(btn).size(11.0), |ui| {
                            ui.set_min_width(180.0);
                            let mut chk =
                                |ui: &mut egui::Ui, on: &mut bool, label, hint, col: [u8; 3]| {
                                    let txt = egui::RichText::new(label)
                                        .size(12.0)
                                        .color(C::from_rgb(col[0], col[1], col[2]));
                                    if ui.checkbox(on, txt).on_hover_text(hint).changed() {
                                        ws.texture = None;
                                    }
                                };
                            chk(ui, &mut ws.overlays.rift, "Rift (tendre C-3)", "Croûte jeune (age≈0), la classe tendre C-3 — teal", [40, 185, 175]);
                            chk(ui, &mut ws.overlays.subduction, "Subduction", "Marge continentale chevauchante (orange) + slab océanique plongeant (bleu)", [245, 140, 30]);
                            chk(ui, &mut ws.overlays.collision, "Collision / accrétion", "Convergence continent-continent (empreinte orogène) — magenta", [230, 60, 200]);
                            chk(ui, &mut ws.overlays.craton, "Craton (socle dur)", "Bouclier cratonique dur (intérieur continental ancien) — or", [205, 178, 100]);
                            chk(ui, &mut ws.overlays.divergent, "Divergent", "Frontières divergentes (axes d'accrétion) — cyan", [50, 220, 220]);
                            chk(ui, &mut ws.overlays.lithology, "Litho C-3 (classes)", "Classes lithologiques : rift tendre (teal) + footprints volcaniclastiques (violet)", [175, 70, 210]);
                            ui.separator();
                            ui.horizontal(|ui| {
                                ui.label(egui::RichText::new("opacité").size(11.0).color(DIM2));
                                if ui
                                    .add(
                                        egui::Slider::new(&mut ws.overlays.opacity, 0.15..=0.9)
                                            .show_value(false),
                                    )
                                    .changed()
                                {
                                    ws.texture = None;
                                }
                            });
                        })
                        .response
                        .on_hover_text("Superpose des couches tectoniques causales sur la vue courante");
                    }
                });
            });
    });
}

/// Bird-eye minimap: the whole map thumbnail with the current view window outlined.
fn minimap_overlay(
    ui: &mut egui::Ui,
    rect: egui::Rect,
    ws: &WorkspaceState,
    hd: &HdResult,
    uv: egui::Rect,
) {
    let _ = hd;
    let side = 132.0f32;
    let mm = egui::Rect::from_min_size(
        egui::pos2(rect.left() + 10.0, rect.bottom() - side - 10.0),
        egui::vec2(side, side),
    );
    let p = ui.painter_at(rect);
    p.rect_filled(mm, 4.0, C::from_rgba_unmultiplied(10, 10, 10, 220));
    if let Some(tex) = ws.texture.as_ref() {
        // Full map thumbnail (uv 0..1) inside the minimap frame.
        let inner = mm.shrink(3.0);
        egui::Image::new(egui::load::SizedTexture::new(tex.id(), inner.size()))
            .uv(egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)))
            .paint_at(ui, inner);
        // View window rectangle (uv → minimap coords).
        let inner_sz = inner.size();
        let win = egui::Rect::from_min_max(
            egui::pos2(inner.left() + uv.min.x * inner_sz.x, inner.top() + uv.min.y * inner_sz.y),
            egui::pos2(inner.left() + uv.max.x * inner_sz.x, inner.top() + uv.max.y * inner_sz.y),
        );
        p.rect_stroke(win, 0.0, egui::Stroke::new(1.5, COPPER_BRIGHT), egui::StrokeKind::Middle);
    }
    p.rect_stroke(
        mm,
        4.0,
        egui::Stroke::new(1.0, C::from_rgb(0x2e, 0x2e, 0x2e)),
        egui::StrokeKind::Middle,
    );
}

fn overlay_chip(ui: &mut egui::Ui, rect: egui::Rect, layer: HdLayer) {
    let p = ui.painter_at(rect);
    let pos = egui::pos2(rect.left() + 12.0, rect.top() + 12.0);
    let title = layer.label();
    let desc = layer.desc();
    let tg = p.layout_no_wrap(title.to_string(), egui::FontId::proportional(13.0), TEXT_BRIGHT);
    let dg = p.layout(desc.to_string(), egui::FontId::proportional(10.5), DIM, 250.0);
    let w = tg.size().x.max(dg.size().x) + 26.0;
    let h = tg.size().y + dg.size().y + 18.0;
    let box_rect = egui::Rect::from_min_size(pos, egui::vec2(w, h));
    p.rect_filled(box_rect, 7.0, C::from_rgba_unmultiplied(18, 18, 18, 210));
    p.circle_filled(egui::pos2(pos.x + 11.0, pos.y + 11.0), 3.5, COPPER_BRIGHT);
    p.galley(egui::pos2(pos.x + 20.0, pos.y + 6.0), tg, TEXT_BRIGHT);
    p.galley(egui::pos2(pos.x + 13.0, pos.y + 22.0), dg, DIM);
}

/// `relief` = (the relief view, the saturation index, the difference's stats, a reference is stored), Finding 135-V.
fn legend_box(
    ui: &mut egui::Ui,
    rect: egui::Rect,
    layer: HdLayer,
    relief: (usize, usize, Option<(usize, usize)>, bool),
    geo_items: Vec<(C, String, String)>,
) {
    let diff_view = layer == HdLayer::Relief && relief.0 == 2 && relief.3;
    let p = ui.painter_at(rect);
    let items: Vec<(C, String, String)> = match layer {
        HdLayer::Relief => vec![], // scale below
        HdLayer::Drainage => vec![
            (C::from_rgb(0x5A, 0xA0, 0xF0), "Barque", "small-boat"),
            (C::from_rgb(0x28, 0x6E, 0xE6), "Chaland", "barge"),
            (C::from_rgb(0x14, 0x46, 0xC8), "Navire", "ship"),
            (C::from_rgb(0x1E, 0x5A, 0xB4), "Lac", ""),
        ],
        HdLayer::Precipitation => vec![
            (C::from_rgb(0xE1, 0xC8, 0x8C), "Désert", "<250"),
            (C::from_rgb(0xC8, 0xC3, 0x6E), "Steppe", "250–500"),
            (C::from_rgb(0x96, 0xB4, 0x5A), "Tempéré-sec", "500–800"),
            (C::from_rgb(0x50, 0x96, 0xC8), "Océanique", "800–1500"),
            (C::from_rgb(0x1E, 0x5A, 0xC8), "Très humide", ">1500"),
        ],
        HdLayer::Temperature => vec![
            (C::from_rgb(0xE1, 0xEB, 0xF8), "Polaire", "<−5°"),
            (C::from_rgb(0x5A, 0x8C, 0xCD), "Boréal", "−5–5°"),
            (C::from_rgb(0x6E, 0xBE, 0x6E), "Tempéré", "5–20°"),
            (C::from_rgb(0xE1, 0x78, 0x46), "Chaud", ">20°"),
        ],
        HdLayer::Biomes => (0..10).map(|i| (biome_hex(i), biome_fr(i), "")).collect(),
        HdLayer::Rocks => Vec::new(),
    }
    .into_iter()
    .map(|(c, a, b): (C, &str, &str)| (c, a.to_string(), b.to_string()))
    .chain(geo_items)
    .collect();
    let title = match layer {
        HdLayer::Relief if diff_view => "DIFFÉRENCE COURANT − RÉFÉRENCE (M)",
        HdLayer::Relief if relief.0 == 1 => "RELIEF — OMBRAGE (AZ 315°, H 45°)",
        HdLayer::Relief => "RELIEF — HYPSOMÉTRIE",
        HdLayer::Drainage => "DRAINAGE — NAVIGABILITÉ",
        HdLayer::Precipitation => "PRÉCIPITATION (MM/AN)",
        HdLayer::Temperature => "TEMPÉRATURE",
        HdLayer::Biomes => "BIOMES",
        HdLayer::Rocks => "ROCHES",
    };
    let row_h = 16.0;
    let rows = if diff_view { 3.0 } else if items.is_empty() { 2.0 } else { items.len() as f32 };
    let bh = 26.0 + rows * row_h;
    let bw = 210.0;
    let bpos = egui::pos2(rect.left() + 12.0, rect.bottom() - bh - 12.0);
    let box_rect = egui::Rect::from_min_size(bpos, egui::vec2(bw, bh));
    p.rect_filled(box_rect, 8.0, C::from_rgba_unmultiplied(18, 18, 18, 220));
    p.text(
        egui::pos2(bpos.x + 12.0, bpos.y + 9.0),
        egui::Align2::LEFT_TOP,
        title,
        egui::FontId::proportional(9.5),
        C::from_rgb(0x7a, 0x7a, 0x7a),
    );
    if diff_view {
        // ADR Finding 135-V -- the signed scale, its saturation, the saturated cells
        let sat = DIFF_SAT_M[relief.1.min(2)];
        let gy = bpos.y + 26.0;
        let gx0 = bpos.x + 12.0;
        let gw = bw - 24.0;
        let steps = 40;
        for s in 0..steps {
            let t = 2.0 * s as f32 / (steps - 1) as f32 - 1.0;
            let [r, g, b] = diff_rgb(t);
            p.rect_filled(
                egui::Rect::from_min_size(
                    egui::pos2(gx0 + gw * (t + 1.0) / 2.0, gy),
                    egui::vec2(gw / steps as f32 + 1.0, 9.0),
                ),
                0.0,
                C::from_rgb(r, g, b),
            );
        }
        for (f, lbl) in [(0.0, format!("−{sat:.0}")), (0.5, "0".to_string()), (1.0, format!("+{sat:.0}"))] {
            p.text(
                egui::pos2(gx0 + gw * f, gy + 12.0),
                if f == 0.0 {
                    egui::Align2::LEFT_TOP
                } else if f == 1.0 {
                    egui::Align2::RIGHT_TOP
                } else {
                    egui::Align2::CENTER_TOP
                },
                lbl,
                egui::FontId::proportional(8.0),
                DIM,
            );
        }
        let (ns, nd) = relief.2.unwrap_or((0, 0));
        p.text(
            egui::pos2(gx0, gy + 26.0),
            egui::Align2::LEFT_TOP,
            format!("saturé au-delà de ±{sat:.0} m · {ns} cellules saturées / {nd} ≠ 0"),
            egui::FontId::proportional(8.0),
            DIM,
        );
    } else if layer == HdLayer::Relief {
        // Gradient scale.
        let gy = bpos.y + 26.0;
        let gx0 = bpos.x + 12.0;
        let gw = bw - 24.0;
        let stops = [
            (0.0, C::from_rgb(0x0A, 0x19, 0x46)),
            (0.3, C::from_rgb(0x8C, 0xC3, 0xD7)),
            (0.45, C::from_rgb(0x46, 0x82, 0x46)),
            (0.75, C::from_rgb(0xC8, 0xB4, 0x78)),
            (1.0, C::from_rgb(0xE6, 0xE6, 0xEB)),
        ];
        let steps = 40;
        for s in 0..steps {
            let t = s as f32 / (steps - 1) as f32;
            let col = grad_at(&stops, t);
            let x = gx0 + gw * t;
            p.rect_filled(
                egui::Rect::from_min_size(
                    egui::pos2(x, gy),
                    egui::vec2(gw / steps as f32 + 1.0, 9.0),
                ),
                0.0,
                col,
            );
        }
        for (i, lbl) in ["Abysse", "Plateau", "Plaine", "Collines", "Sommets"].iter().enumerate() {
            let x = gx0 + gw * (i as f32 / 4.0);
            p.text(
                egui::pos2(x, gy + 12.0),
                egui::Align2::LEFT_TOP,
                *lbl,
                egui::FontId::proportional(8.0),
                DIM,
            );
        }
    } else {
        for (i, (col, lbl, sub)) in items.iter().enumerate() {
            let ry = bpos.y + 26.0 + i as f32 * row_h;
            p.rect_filled(
                egui::Rect::from_min_size(
                    egui::pos2(bpos.x + 12.0, ry + 1.0),
                    egui::vec2(11.0, 11.0),
                ),
                3.0,
                *col,
            );
            p.text(
                egui::pos2(bpos.x + 28.0, ry),
                egui::Align2::LEFT_TOP,
                lbl,
                egui::FontId::proportional(11.0),
                C::from_rgb(0xcf, 0xcf, 0xcf),
            );
            if !sub.is_empty() {
                p.text(
                    egui::pos2(bpos.x + bw - 12.0, ry),
                    egui::Align2::RIGHT_TOP,
                    sub,
                    egui::FontId::monospace(9.5),
                    C::from_rgb(0x77, 0x77, 0x77),
                );
            }
        }
    }
}

/// ADR Finding 153 -- the « Chances » mask's legend, top left of the map.
fn chances_legend(ui: &mut egui::Ui, rect: egui::Rect, items: &[(C, String, String)], all: bool) {
    let p = ui.painter_at(rect);
    let row_h = 15.0;
    let bh = 40.0 + items.len() as f32 * row_h;
    let bw = 230.0;
    let bpos = egui::pos2(rect.left() + 12.0, rect.top() + 12.0);
    p.rect_filled(egui::Rect::from_min_size(bpos, egui::vec2(bw, bh)), 8.0, C::from_rgba_unmultiplied(18, 18, 18, 220));
    p.text(egui::pos2(bpos.x + 12.0, bpos.y + 8.0), egui::Align2::LEFT_TOP, "CHANCES — OPACITÉ = FAVORABILITÉ", egui::FontId::proportional(9.5), C::from_rgb(0x7a, 0x7a, 0x7a));
    p.text(
        egui::pos2(bpos.x + 12.0, bpos.y + 21.0),
        egui::Align2::LEFT_TOP,
        if all { "toutes : la plus forte l'emporte · grille ¼ (195 m)" } else { "une ressource · grille ¼ (195 m)" },
        egui::FontId::proportional(8.5),
        DIM,
    );
    for (i, (col, lbl, sub)) in items.iter().enumerate() {
        let ry = bpos.y + 36.0 + i as f32 * row_h;
        p.rect_filled(egui::Rect::from_min_size(egui::pos2(bpos.x + 12.0, ry + 1.0), egui::vec2(10.0, 10.0)), 3.0, *col);
        p.text(egui::pos2(bpos.x + 28.0, ry), egui::Align2::LEFT_TOP, lbl, egui::FontId::proportional(10.5), C::from_rgb(0xcf, 0xcf, 0xcf));
        p.text(egui::pos2(bpos.x + bw - 12.0, ry), egui::Align2::RIGHT_TOP, sub, egui::FontId::monospace(9.0), C::from_rgb(0x77, 0x77, 0x77));
    }
}

/// ADR Finding 153 -- a toolbar menu that stays open while its items are used: it closes on a click OUTSIDE it or on
/// Escape (egui's default closes on any click, so a checkbox closed its menu).
fn tool_menu<R>(ui: &mut egui::Ui, text: egui::RichText, add: impl FnOnce(&mut egui::Ui) -> R) -> egui::InnerResponse<Option<R>> {
    use egui::containers::menu::{MenuButton, MenuConfig};
    let (resp, inner) = MenuButton::new(text).config(MenuConfig::new().close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)).ui(ui, add);
    egui::InnerResponse::new(inner.map(|i| i.inner), resp)
}

fn zoom_controls(ui: &mut egui::Ui, rect: egui::Rect, ws: &mut WorkspaceState) {
    let w = 150.0;
    let area = egui::Rect::from_min_size(
        egui::pos2(rect.right() - w - 12.0, rect.bottom() - 38.0),
        egui::vec2(w, 26.0),
    );
    let mut child = ui.new_child(
        egui::UiBuilder::new()
            .max_rect(area)
            .layout(egui::Layout::left_to_right(egui::Align::Center)),
    );
    // The zoom is the absolute camera scale `map_px` relative to the fit size for THIS viewport.
    let fit = rect.width().min(rect.height()).max(64.0);
    let cur = if ws.map_px <= 0.0 { fit } else { ws.map_px };
    child.horizontal(|ui| {
        if ui.small_button("−").clicked() {
            ws.map_px = (cur / 1.25).max(fit * 0.5);
        }
        ui.label(
            egui::RichText::new(format!("{:.0}%", cur / fit * 100.0))
                .monospace()
                .color(C::from_rgb(0xbb, 0xbb, 0xbb))
                .size(11.0),
        );
        if ui.small_button("+").clicked() {
            ws.map_px = (cur * 1.25).min(fit * 16.0);
        }
        if ui.small_button("⤢").on_hover_text("Recadrer").clicked() {
            ws.map_px = fit;
            ws.map_pan = [0.5, 0.5];
        }
    });
}

// ── Helpers: colours, zones, biome names ─────────────────────────────────
fn zone_name(l: f32) -> &'static str {
    if l < 10.0 {
        "équatoriale"
    } else if l < 30.0 {
        "subtropicale"
    } else if l < 45.0 {
        "tempérée"
    } else if l < 60.0 {
        "subpolaire"
    } else {
        "polaire"
    }
}

fn grad_at(stops: &[(f32, C)], t: f32) -> C {
    if t <= stops[0].0 {
        return stops[0].1;
    }
    for w in stops.windows(2) {
        let (t0, c0) = w[0];
        let (t1, c1) = w[1];
        if t <= t1 {
            let f = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0);
            let lerp = |a: u8, b: u8| (a as f32 + (b as f32 - a as f32) * f) as u8;
            return C::from_rgb(lerp(c0.r(), c1.r()), lerp(c0.g(), c1.g()), lerp(c0.b(), c1.b()));
        }
    }
    stops[stops.len() - 1].1
}

fn biome_hex(i: usize) -> C {
    const H: [[u8; 3]; 10] = [
        [0x1E, 0x32, 0x5A],
        [0xC8, 0xCD, 0xD7],
        [0x46, 0x6E, 0x5A],
        [0xC8, 0xC3, 0x6E],
        [0x50, 0xA0, 0x50],
        [0x28, 0x6E, 0x46],
        [0xE1, 0xC8, 0x8C],
        [0xBE, 0xAF, 0x5A],
        [0x78, 0xAF, 0x46],
        [0x14, 0x6E, 0x32],
    ];
    let c = H[i.min(9)];
    C::from_rgb(c[0], c[1], c[2])
}

fn biome_fr(i: usize) -> &'static str {
    [
        "Océan",
        "Toundra",
        "Taïga",
        "Steppe",
        "Forêt tempérée",
        "Forêt pluviale tempérée",
        "Désert",
        "Savane",
        "Forêt tropicale saisonnière",
        "Forêt tropicale",
    ][i.min(9)]
}

fn french_biome(b: Biome) -> &'static str {
    match b {
        Biome::Ocean => "Océan",
        Biome::Tundra => "Toundra",
        Biome::BorealForest => "Taïga",
        Biome::TemperateGrassland => "Steppe",
        Biome::TemperateForest => "Forêt tempérée",
        Biome::TemperateRainforest => "Forêt pluviale tempérée",
        Biome::Desert => "Désert",
        Biome::Savanna => "Savane",
        Biome::TropicalSeasonalForest => "Forêt tropicale saisonnière",
        Biome::TropicalRainforest => "Forêt tropicale",
        Biome::Lake => "Lac (eau intérieure)",
        Biome::Wetland => "Zone humide (fluente, douce)",
    }
}

/// Tectonic-debug OVERLAY toggles (the "microscope"): semi-transparent causal
/// layers the author superimposes on any base view to SEE where each closure acts.
/// Registered to the terrain via the HD window (`hd.sample_origin`/`sample_size`).
#[derive(Clone, Copy, PartialEq)]
struct TectonicOverlays {
    /// Rift-soft (continental `age≈0`) — the C-3 soft class, teal.
    rift: bool,
    /// Subduction: overriding continental margin (orange) + subducting slab (blue).
    subduction: bool,
    /// Continental collision / accretion footprint (C-C convergent), magenta.
    collision: bool,
    /// Hard cratonic shield, gold.
    craton: bool,
    /// Divergent (spreading) boundaries, cyan.
    divergent: bool,
    /// C-3 lithology CLASS: rift-soft (teal) + volcaniclastic footprints (violet).
    lithology: bool,
    /// Blend strength (0 = base only, 1 = full tint).
    opacity: f32,
}

impl Default for TectonicOverlays {
    fn default() -> Self {
        Self {
            rift: false,
            subduction: false,
            collision: false,
            craton: false,
            divergent: false,
            lithology: false,
            opacity: 0.55,
        }
    }
}

impl TectonicOverlays {
    fn any(&self) -> bool {
        self.rift
            || self.subduction
            || self.collision
            || self.craton
            || self.divergent
            || self.lithology
    }
}

/// Blend the enabled tectonic overlays into the data-space RGBA (BEFORE the north-up
/// flip, so it registers with the coarse labels). Mask layers are nearest-sampled
/// from `hd.tectonic`; the volcaniclastic footprints are rasterised from `hd.edifices`
/// (the same basal discs the C-3 K field stamps).
fn blend_tectonic_overlays(rgba: &mut [u8], hd: &HdResult, ov: TectonicOverlays) {
    let Some(lab) = hd.tectonic.as_ref() else { return };
    let (w, h) = (hd.width, hd.height);
    let (nx, ny) = (lab.nx, lab.ny);
    let a = ov.opacity.clamp(0.0, 1.0);
    let sx0 = hd.sample_origin[0] * nx as f64;
    let sy0 = hd.sample_origin[1] * ny as f64;
    let dx = hd.sample_size * nx as f64 / w as f64;
    let dy = hd.sample_size * ny as f64 / h as f64;
    let blend = |rgba: &mut [u8], k: usize, tint: [u8; 3]| {
        for c in 0..3 {
            let base = rgba[k * 4 + c] as f32;
            rgba[k * 4 + c] = (base * (1.0 - a) + tint[c] as f32 * a) as u8;
        }
    };
    // Mask layers, priority low→high (later overrides at a shared cell).
    for j in 0..h {
        let cj = ((sy0 + j as f64 * dy).round() as isize).rem_euclid(ny as isize) as usize;
        for i in 0..w {
            let ci = ((sx0 + i as f64 * dx).round() as isize).rem_euclid(nx as isize) as usize;
            let ck = cj * nx + ci;
            let k = j * w + i;
            let mut tint: Option<[u8; 3]> = None;
            if ov.craton && lab.craton[ck] {
                tint = Some([205, 178, 100]); // gold shield
            }
            if ov.divergent && lab.divergent[ck] {
                tint = Some([50, 220, 220]); // cyan
            }
            if ov.collision && lab.collision[ck] {
                tint = Some([230, 60, 200]); // magenta
            }
            if ov.subduction && lab.subduction_slab[ck] {
                tint = Some([70, 130, 225]); // subducting slab, blue
            }
            if ov.subduction && lab.subduction_upper[ck] {
                tint = Some([245, 140, 30]); // overriding margin, orange
            }
            if (ov.rift || ov.lithology) && lab.rift[ck] {
                tint = Some([40, 185, 175]); // rift-soft, teal
            }
            if let Some(t) = tint {
                blend(rgba, k, t);
            }
        }
    }
    // Volcaniclastic footprints (only for the lithology class overlay) — rasterise the
    // edifice basal discs, exactly as `lithology::stamp_volcanic_k` does at HD.
    if ov.lithology {
        let (so, ss) =
            ([hd.sample_origin[0] as f32, hd.sample_origin[1] as f32], hd.sample_size as f32);
        for e in &hd.edifices {
            let fx = (e.center_uv.0 - so[0]).rem_euclid(1.0) / ss;
            let fy = (e.center_uv.1 - so[1]).rem_euclid(1.0) / ss;
            if fx >= 1.0 || fy >= 1.0 {
                continue;
            }
            let (cx, cy) = (fx * w as f32, fy * h as f32);
            let rb = (e.basal_diameter_km * 0.5) / hd.km_per_cell;
            if rb < 1.0 {
                continue;
            }
            let (i0, i1) =
                ((cx - rb).floor().max(0.0) as usize, ((cx + rb).ceil() as usize).min(w - 1));
            let (j0, j1) =
                ((cy - rb).floor().max(0.0) as usize, ((cy + rb).ceil() as usize).min(h - 1));
            for j in j0..=j1 {
                for i in i0..=i1 {
                    let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                    if d <= rb {
                        blend(rgba, j * w + i, [175, 70, 210]); // volcaniclastic, violet
                    }
                }
            }
        }
    }
}

/// ADR Finding 135-V -- how the Relief layer is drawn. A VIEW: it reads the world and never writes it.
#[derive(Clone, Copy)]
enum ReliefView<'a> {
    /// The shipped hypsometric layer.
    Hypso,
    /// The hypsometric colour times a hillshade; `outline` draws the lakes' outline over it (ADR Finding 157-V).
    Shade { outline: bool },
    /// The current world minus a reference world (m), saturated at the given magnitude.
    Diff(&'a HdResult, f32),
}

/// ADR Finding 135-V -- the difference view's saturations (m).
const DIFF_SAT_M: [f32; 3] = [10.0, 100.0, 1000.0];

/// ADR Finding 135-V -- the signed palette, `t` in [−1, 1]: blue (lower) → white → red (higher).
fn diff_rgb(t: f32) -> [u8; 3] {
    let t = t.clamp(-1.0, 1.0);
    if t >= 0.0 {
        let v = (255.0 * (1.0 - t)) as u8;
        [255, v, v]
    } else {
        let v = (255.0 * (1.0 + t)) as u8;
        [v, v, 255]
    }
}

/// A short label for a guard verdict (the difference view's reference line).
fn guard_short(g: &ymir_core::tectonics_c1::bench_guard::GuardStatus) -> String {
    use ymir_core::tectonics_c1::bench_guard::GuardStatus;
    match g {
        GuardStatus::Match { label } => format!("= banc ({label})"),
        GuardStatus::Mismatch { label, .. } => format!("≠ banc ({label})"),
        GuardStatus::NoReference { viz } => format!("non gardé (champ {viz})"),
    }
}

/// ADR Finding 135-V -- a Lambertian hillshade of the CONDITIONED field in metres: light from azimuth 315°
/// (north-west), altitude 45°, normalised so that flat ground reads 1. The field is stored south-first, so data
/// `y + 1` is north.
fn hillshade(hd: &HdResult) -> Vec<f32> {
    use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
    use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    let ss = SteinSteinParams::default();
    let (w, h) = (hd.width, hd.height);
    let z: Vec<f32> = hd.eroded.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
    let cell = (hd.km_per_cell * 1000.0).max(1e-3);
    let (az, alt) = (315f32.to_radians(), 45f32.to_radians());
    let l = [alt.cos() * az.sin(), alt.cos() * az.cos(), alt.sin()];
    let mut s = vec![1.0f32; w * h];
    for y in 0..h {
        let (ym, yp) = (y.saturating_sub(1), (y + 1).min(h - 1));
        for x in 0..w {
            let (xm, xp) = (x.saturating_sub(1), (x + 1).min(w - 1));
            let dzdx = (z[y * w + xp] - z[y * w + xm]) / (cell * (xp - xm).max(1) as f32);
            let dzdy = (z[yp * w + x] - z[ym * w + x]) / (cell * (yp - ym).max(1) as f32);
            let nn = (dzdx * dzdx + dzdy * dzdy + 1.0).sqrt();
            s[y * w + x] = ((-dzdx * l[0] - dzdy * l[1] + l[2]) / nn).max(0.0) / l[2];
        }
    }
    s
}

// ── Layer → RGBA image (canonical palettes) ──────────────────────────────
/// The second value is the difference view's (saturated cells, cells ≠ 0), `None` for every other view.
fn layer_color_image(
    hd: &HdResult,
    layer: HdLayer,
    river_map: &RiverCellMap,
    overlay: bool,
    tectonic: TectonicOverlays,
    relief: ReliefView<'_>,
    geo: Option<(&ymir_core::geology::Zoning, Option<usize>, &[u8], bool)>,
) -> (egui::ColorImage, Option<(usize, usize)>) {
    use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
    use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    let (w, h) = (hd.width, hd.height);
    let mut rgba = vec![0u8; w * h * 4];
    // ADR Finding 135-V -- the relief view's inputs (a reference of another size falls back to hypsometry)
    let relief = match relief {
        ReliefView::Diff(r, _) if (r.width, r.height) != (w, h) => ReliefView::Hypso,
        v => v,
    };
    let shade = match (layer, relief) {
        (HdLayer::Relief, ReliefView::Shade { .. } | ReliefView::Diff(..)) => Some(hillshade(hd)),
        _ => None,
    };
    let ss = SteinSteinParams::default();
    let (mut n_sat, mut n_diff) = (0usize, 0usize);
    for k in 0..w * h {
        let c = match layer {
            HdLayer::Relief => match relief {
                ReliefView::Hypso => relief_color(hd.eroded.data[k]),
                ReliefView::Shade { .. } => {
                    let [r, g, b] = relief_color(hd.eroded.data[k]);
                    let f = (0.25 + 0.75 * shade.as_ref().map_or(1.0, |s| s[k])).min(1.35);
                    [(r as f32 * f).min(255.0) as u8, (g as f32 * f).min(255.0) as u8, (b as f32 * f).min(255.0) as u8]
                }
                ReliefView::Diff(reference, sat) => {
                    let d = c1_altitude_norm_to_metres(hd.eroded.data[k], &ss)
                        - c1_altitude_norm_to_metres(reference.eroded.data[k], &ss);
                    if d != 0.0 {
                        n_diff += 1;
                    }
                    if d.abs() >= sat {
                        n_sat += 1;
                    }
                    let [r, g, b] = diff_rgb(d / sat);
                    // a faint shade of the current world, for orientation only
                    let f = (0.8 + 0.2 * shade.as_ref().map_or(1.0, |s| s[k])).min(1.0);
                    [(r as f32 * f) as u8, (g as f32 * f) as u8, (b as f32 * f) as u8]
                }
            },
            HdLayer::Precipitation => precip_color(precip_mm_per_year(hd.precipitation.data[k])),
            HdLayer::Temperature => temp_color(hd.temperature.data[k]),
            HdLayer::Biomes => {
                let [r, g, b] = hd.biomes[k].color();
                [r, g, b]
            }
            HdLayer::Drainage => drainage_color(hd, river_map, k),
            // ADR Finding 152 -- the geology views draw over the hypsometric relief (water stays as it is)
            HdLayer::Rocks => relief_color(hd.eroded.data[k]),
        };
        rgba[k * 4] = c[0];
        rgba[k * 4 + 1] = c[1];
        rgba[k * 4 + 2] = c[2];
        rgba[k * 4 + 3] = 255;
    }
    if let Some((zoning, only, rocks, chances)) = geo {
        if layer == HdLayer::Rocks {
            ymir_core::geology::render::rocks_rgba(rocks, &mut rgba);
        }
        // ADR Finding 153 -- the « Chances » mask, over any step
        if chances {
            ymir_core::geology::render::chances_rgba(zoning, only, w, h, &mut rgba);
        }
    }
    // ADR Finding 157-V -- the lakes' outline on the hillshade (before the rivers, which stay on top)
    if layer == HdLayer::Relief && matches!(relief, ReliefView::Shade { outline: true }) {
        draw_lake_outline(&mut rgba, &lake_outline_mask(&hd.drainage.lake_map, w, h));
    }
    if overlay {
        draw_river_overlay(&mut rgba, hd);
    }
    if tectonic.any() {
        blend_tectonic_overlays(&mut rgba, hd, tectonic);
    }
    // Finding 27 — the field is stored y=0=SOUTH (the export/LL contract); the VIEWER
    // shows NORTH-UP, so mirror the texture rows here. The DATA, the export and the
    // inspector's cell lookup all stay south-first; only the pixels are flipped.
    flip_rows_rgba(&mut rgba, w, h);
    let stats = matches!((layer, relief), (HdLayer::Relief, ReliefView::Diff(..))).then_some((n_sat, n_diff));
    (egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba), stats)
}

/// ADR Finding 157-V -- the lakes' outline: 1 on a lake cell with a 4-neighbour outside its lake (the inner line), 2 on a
/// non-lake cell with a 4-neighbour in a lake (the outer line), 0 elsewhere. Two lakes side by side each get their inner
/// line; the map's border closes nothing.
fn lake_outline_mask(lake_map: &[u32], w: usize, h: usize) -> Vec<u8> {
    let mut m = vec![0u8; w * h];
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            let id = lake_map[k];
            let nb = [(x > 0).then(|| k - 1), (x + 1 < w).then(|| k + 1), (y > 0).then(|| k - w), (y + 1 < h).then(|| k + w)];
            if id != 0 {
                if nb.iter().flatten().any(|&j| lake_map[j] != id) {
                    m[k] = 1;
                }
            } else if nb.iter().flatten().any(|&j| lake_map[j] != 0) {
                m[k] = 2;
            }
        }
    }
    m
}

/// ADR Finding 157-V -- paint [`lake_outline_mask`]: a bright cyan inner line on the water, a dark navy outer line on the
/// shore, so the outline reads on a light slope and on a dark one alike.
fn draw_lake_outline(rgba: &mut [u8], mask: &[u8]) {
    for (k, &v) in mask.iter().enumerate() {
        let c = match v {
            1 => [0u8, 229, 255],
            2 => [10, 22, 48],
            _ => continue,
        };
        rgba[k * 4..k * 4 + 3].copy_from_slice(&c);
    }
}

/// Vertically mirror a row-major RGBA buffer in place: display row `j` ↔ data row
/// `h-1-j`. Turns the south-first stored field into a north-up image (Finding 27).
fn flip_rows_rgba(rgba: &mut [u8], w: usize, h: usize) {
    for j in 0..h / 2 {
        let (a, b) = (j * w * 4, (h - 1 - j) * w * 4);
        for k in 0..w * 4 {
            rgba.swap(a + k, b + k);
        }
    }
}

/// TASK 1 — the latitude-placement widget (Finding 28). A vertical −90…+90° strip (north up,
/// matching the map) painted with the thermal gradient AND the wind BELTS as bands (trade
/// easterlies / westerlies / polar easterlies, with a direction arrow each), and the current
/// map rectangle `[centre−span/2, centre+span/2]` drawn over it. Makes the CONSEQUENCE of the
/// span visible — how many belts it crosses — before generating. Pure function of centre/span;
/// no generation, no data read.
fn latitude_placement_widget(ui: &mut egui::Ui, centre: f32, span: f32) {
    let width = ui.available_width();
    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, 168.0), egui::Sense::hover());
    let p = ui.painter_at(rect);
    let h = rect.height();
    let y_of = |lat: f32| rect.top() + (90.0 - lat) / 180.0 * h; // +90 top, −90 bottom
    // Thermal gradient background (40 bands): cold blue → warm red via sea-level T.
    let bands = 40usize;
    for b in 0..bands {
        let lat_hi = 90.0 - (b as f32 / bands as f32) * 180.0;
        let lat_lo = 90.0 - ((b + 1) as f32 / bands as f32) * 180.0;
        let t = sea_level_temperature((lat_hi + lat_lo) * 0.5);
        let f = ((t + 25.0) / 52.0).clamp(0.0, 1.0); // -25..27 → 0..1
        let col = C::from_rgb(
            (40.0 + 200.0 * f) as u8,
            (90.0 + 60.0 * (1.0 - (f - 0.5).abs() * 2.0)) as u8,
            (200.0 - 170.0 * f) as u8,
        );
        p.rect_filled(
            egui::Rect::from_min_max(
                egui::pos2(rect.left(), y_of(lat_hi)),
                egui::pos2(rect.right(), y_of(lat_lo)),
            ),
            0.0,
            col,
        );
    }
    // Belt boundaries (±30, ±60) + a labelled arrow per belt.
    let faint = C::from_rgba_unmultiplied(255, 255, 255, 60);
    for lat in [-60.0f32, -30.0, 0.0, 30.0, 60.0] {
        let y = y_of(lat);
        p.line_segment(
            [egui::pos2(rect.left(), y), egui::pos2(rect.right(), y)],
            egui::Stroke::new(1.0, faint),
        );
    }
    let belt = |mid: f32, name: &str| {
        let dir = wind_zonal_dir(mid);
        let arrow = if dir > 0 { "→" } else { "←" };
        p.text(
            egui::pos2(rect.left() + 5.0, y_of(mid)),
            egui::Align2::LEFT_CENTER,
            format!("{arrow} {name}"),
            egui::FontId::proportional(9.5),
            C::from_rgba_unmultiplied(255, 255, 255, 200),
        );
    };
    belt(75.0, "polaires");
    belt(45.0, "westerlies");
    belt(15.0, "alizés");
    belt(-15.0, "alizés");
    belt(-45.0, "westerlies");
    belt(-75.0, "polaires");
    // Equator marker.
    p.text(
        egui::pos2(rect.right() - 4.0, y_of(0.0)),
        egui::Align2::RIGHT_CENTER,
        "0°",
        egui::FontId::monospace(9.0),
        C::from_rgba_unmultiplied(255, 255, 255, 150),
    );
    // The map rectangle at [centre−span/2, centre+span/2].
    let (y0, y1) = (y_of(centre + span / 2.0), y_of(centre - span / 2.0));
    let map_rect = egui::Rect::from_min_max(
        egui::pos2(rect.left() + 1.0, y0),
        egui::pos2(rect.right() - 1.0, y1),
    );
    p.rect_filled(map_rect, 2.0, C::from_rgba_unmultiplied(0xC9, 0x85, 0x3F, 40));
    p.rect_stroke(map_rect, 2.0, egui::Stroke::new(2.0, COPPER_BRIGHT), egui::StrokeKind::Middle);
    p.text(
        egui::pos2(map_rect.center().x, y0 - 1.0),
        egui::Align2::CENTER_BOTTOM,
        format!("{:.0}°", centre + span / 2.0),
        egui::FontId::monospace(9.0),
        COPPER_BRIGHT,
    );
    p.text(
        egui::pos2(map_rect.center().x, y1 + 1.0),
        egui::Align2::CENTER_TOP,
        format!("{:.0}°", centre - span / 2.0),
        egui::FontId::monospace(9.0),
        COPPER_BRIGHT,
    );
}

/// Classify where a mouth cell `(x,y)` drains → (sink kind, the sink lake id if any). Reads
/// `lake_map`, lake types, and `water_class` — no recompute. NOTE (Finding 31): the SEA label
/// is derived ONLY from `water_class == 1` (the authority — flood-fill from the domain borders),
/// NEVER from altitude. A low-altitude cell that is an enclosed below-sea basin (class 2) is NOT
/// the sea; testing altitude here mislabelled rivers draining into sub-threshold basins as
/// "outlet = sea". Lake membership (any marked basin, regardless of inventory size) wins first.
fn classify_sink(hd: &HdResult, wc: &[u8], x: u32, y: u32) -> (Sink, Option<u32>) {
    let (w, h) = (hd.width, hd.height);
    let endo: std::collections::HashSet<u32> = hd
        .drainage
        .lakes
        .iter()
        .filter(|l| l.lake_type == LakeType::Endorheic)
        .map(|l| l.base.id)
        .collect();
    let mut sea = false;
    for dy in -1i32..=1 {
        for dx in -1i32..=1 {
            let (nx, ny) = (x as i32 + dx, y as i32 + dy);
            if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                let k = ny as usize * w + nx as usize;
                let id = hd.drainage.lake_map[k];
                if id != 0 {
                    return (
                        if endo.contains(&id) { Sink::EndoLake } else { Sink::ExoLake },
                        Some(id),
                    );
                }
                if wc.get(k).copied() == Some(1) {
                    sea = true; // ONLY water_class==1 is the ocean
                }
            }
        }
    }
    // Finding 36 — the mouth ends in an ENCLOSED below-sea cell (class 2) that is NOT an inventoried
    // lake: a sub-sea evaporative sink, NOT the sea. Reading class 2 as "→ mer" (or from the −20 m
    // altitude) is exactly the proxy this thread removed; water_class is the authority. The mouth's
    // OWN cell decides — a class-1 neighbour still wins (a true coastal mouth touches the ocean).
    let mk = y as usize * w + x as usize;
    if !sea && wc.get(mk).copied() == Some(2) {
        return (Sink::SubSeaSink, None);
    }
    (if sea { Sink::Sea } else { Sink::Unknown }, None)
}

/// Aggregate the exported river SEGMENTS into browsable WATERCOURSES (Finding 28). A
/// watercourse = every segment that drains to the same terminal (a `downstream == None`
/// mouth after the lake-clip); its TRUNK is the source→sink main stem, chosen by taking
/// the highest-discharge tributary at each confluence. Pure assembly of `C1DrainageResult`
/// (discharge / width / catchment / points) — nothing is recomputed. `ratio` (the current
/// geographic-scale slider) signifies the trunk length, consistent with the exported
/// width/discharge. Sorted by discharge (biggest rivers first).
/// `chain_lakes`: link an inflow reach to its exorheic lake's outlet reach so an entry is a
/// river SYSTEM (ADR Finding 45). `false` reproduces the pre-Finding-45 behaviour, where an
/// entry was a reach between two water bodies — kept so the bench can measure both on one run.
///
/// **HISTORICAL since ADR Finding 124-3** (Findings 45–123's rule, kept for the benches that
/// measure it): the shipped aggregation is [`aggregate_watercourses`].
fn aggregate_watercourses_opt(
    hd: &HdResult,
    domain_km: f32,
    ratio: f32,
    chain_lakes: bool,
) -> Vec<Watercourse> {
    let d = &hd.drainage;
    let wc = ymir_core::lakes::connectivity::water_class(&hd.eroded, 0.5);
    let segs = &d.rivers.segments;
    let n = segs.len();
    let q = |i: usize| d.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
    // ── CHAIN THROUGH AN EXORHEIC LAKE (ADR Finding 45). A reach that ends at a lake shore
    //    gets `downstream = None`, so 87–92 % of terminals are lake INFLOWS and an entry was
    //    a REACH BETWEEN TWO WATER BODIES rather than a river SYSTEM — which is why the head
    //    of the discharge sort showed a fragment. Water physically continues through an
    //    EXORHEIC lake (it has an outflow), so for aggregation the inflow reach is linked to
    //    that lake's outlet reach. NOT for an endorheic lake (the water dies there — a true
    //    terminus) nor for a below-sea basin (its outflow is a typed `Spillway`, no hierarchy).
    let wc_class = ymir_core::lakes::connectivity::water_class(&hd.eroded, 0.5);
    // lake id → the outlet reach: a segment whose FIRST cell borders that lake. Several can
    // qualify on a wide shore; take the highest-discharge one (the main outflow).
    let mut outlet_of: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for (i, sg) in segs.iter().enumerate() {
        let (sx, sy) = sg.points[0];
        if let (Sink::ExoLake, Some(id)) = classify_sink(hd, &wc_class, sx, sy) {
            let better = outlet_of.get(&id).map_or(true, |&j| q(i) > q(j));
            if better {
                outlet_of.insert(id, i);
            }
        }
    }
    // The continuation of a terminal reach, if it dies on an exorheic lake that has an outlet.
    let across_lake = |i: usize| -> Option<usize> {
        let (mx, my) = *segs[i].points.last()?;
        match classify_sink(hd, &wc_class, mx, my) {
            (Sink::ExoLake, Some(id)) => outlet_of.get(&id).copied().filter(|&o| o != i),
            _ => None,
        }
        .filter(|_| chain_lakes)
    };

    // root(i): follow downstream to the terminal (mouth), crossing exorheic lakes; guard
    // against any cycle (a lake link could close one — two lakes spilling into each other).
    let mut root = vec![0usize; n];
    for i in 0..n {
        let mut j = i;
        let mut seen = 0usize;
        loop {
            seen += 1;
            if seen > n {
                break;
            }
            match segs[j].downstream {
                Some(k) if k < n => j = k,
                _ => match across_lake(j) {
                    Some(o) => j = o,
                    None => break,
                },
            }
        }
        root[i] = j;
    }
    let mut groups: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for i in 0..n {
        groups.entry(root[i]).or_default().push(i);
    }
    let km_per_cell = if hd.width > 0 { domain_km / hd.width as f32 } else { 0.0 };
    // Longest upstream path (in points) per segment — memoised, so the trunk is the
    // geographic main stem (river length), not the short high-discharge branch.
    let mut pathlen = vec![0u32; n];
    {
        let mut order: Vec<usize> = (0..n).collect();
        // process shallowest-first: a segment's pathlen needs its upstreams' — but the
        // DAG has no cycles, so a simple relaxation to a fixed point over `order` (few
        // passes) suffices; bound the passes by the deepest Strahler order.
        order.sort_by_key(|&i| segs[i].strahler_order);
        for _ in 0..16 {
            let mut changed = false;
            for &i in &order {
                let up = segs[i]
                    .upstream
                    .iter()
                    .copied()
                    .filter(|&u| u < n)
                    .map(|u| pathlen[u])
                    .max()
                    .unwrap_or(0);
                let v = up + segs[i].points.len() as u32;
                if v != pathlen[i] {
                    pathlen[i] = v;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    // ADR Finding 123-C2 -- the drained area of a reach (signified km², the exported array).
    let area_of = |i: usize| d.segment_drainage_km2.get(i).copied().unwrap_or(0.0);
    // Two areas within 1 % (PROXY) are a TIE -- an area a clipped fragment inherited from its parent
    // (Findings 42/93) or two branches of one catchment differing by a few cells -- and the tie goes
    // to the longest-path rule this replaces. Without it, 3 then 2 delivered trunks lost their lake
    // crossing on ties of 0.00-0.05 %; with it, 0 (f123_identity_control).
    let by_area = |a: usize, b: usize| -> std::cmp::Ordering {
        let (x, y) = (area_of(a), area_of(b));
        if (x - y).abs() <= 0.01 * x.max(y) {
            pathlen[a].cmp(&pathlen[b])
        } else {
            x.total_cmp(&y)
        }
    };
    // The longest reach flowing INTO the lake an outlet reach starts on — the trunk's
    // continuation upstream. Needs `pathlen`, hence defined after it.
    let inflow_to_lake_of = |i: usize| -> Option<usize> {
        let (sx, sy) = segs[i].points[0];
        if !chain_lakes {
            return None;
        }
        let id = match classify_sink(hd, &wc_class, sx, sy) {
            (Sink::ExoLake, Some(id)) => id,
            _ => return None,
        };
        segs
            .iter()
            .enumerate()
            .filter(|&(j, sg)| {
                j != i
                    && sg.downstream.is_none()
                    && sg.points.last().is_some_and(|&(mx, my)| {
                        matches!(classify_sink(hd, &wc_class, mx, my), (Sink::ExoLake, Some(k)) if k == id)
                    })
            })
            .max_by(|&(a, _), &(b, _)| {
                if TRUNK_BY_MAX_AREA {
                    by_area(a, b)
                } else {
                    pathlen[a].cmp(&pathlen[b])
                }
            })
            .map(|(j, _)| j)
    };

    let mut out = Vec::with_capacity(groups.len());
    // DETERMINISM (a core project invariant, and it was broken here): `groups` is a HashMap,
    // whose iteration order is randomised per process. Combined with a stable sort that TIES on
    // discharge, the head of the list between two equal-discharge entries was decided by hash
    // order — so "the head is the trunk" would have been luck, not a property. Iterate by root
    // index instead.
    // Which reaches END on a cell that another reach also ends on (the shared-mouth case).
    let mut terminus_count: std::collections::HashMap<(u32, u32), usize> = Default::default();
    for sg in segs.iter().filter(|sg| sg.downstream.is_none()) {
        *terminus_count.entry(*sg.points.last().unwrap()).or_insert(0) += 1;
    }
    let shared_terminus: Vec<bool> = segs
        .iter()
        .map(|sg| {
            sg.downstream.is_none()
                && terminus_count.get(sg.points.last().unwrap()).copied().unwrap_or(0) > 1
        })
        .collect();

    let mut group_list: Vec<(usize, Vec<usize>)> = groups.into_iter().collect();
    group_list.sort_by_key(|(r, _)| *r);
    for (r, members) in group_list {
        // Trunk: from the mouth, climb the LONGEST tributary at each confluence.
        let mut trunk = Vec::new();
        let mut cur = r;
        for _ in 0..=n {
            trunk.push(cur);
            let ups = segs[cur].upstream.iter().copied().filter(|&u| u < n);
            let next = if TRUNK_BY_MAX_AREA {
                ups.max_by(|&a, &b| by_area(a, b))
            } else {
                ups.max_by_key(|&u| pathlen[u])
            };
            match next {
                Some(next) => cur = next,
                // Symmetric to `across_lake`: at the head of an outlet reach, continue up the
                // LONGEST reach that flows INTO the same lake, so the trunk is the whole main
                // stem rather than the segment below the lake.
                None => match inflow_to_lake_of(cur) {
                    Some(up) => cur = up,
                    None => break,
                },
            }
        }
        let source = *trunk.last().unwrap();
        trunk.reverse(); // source → sink
        let pts: usize = trunk.iter().map(|&s| segs[s].points.len()).sum();
        let (mx, my) = *segs[r].points.last().unwrap();
        let (sink, sink_lake_id) = classify_sink(hd, &wc, mx, my);
        out.push(Watercourse {
            trunk,
            tributaries: members.len().saturating_sub(1),
            segments: members,
            order: segs[r].strahler_order,
            discharge_m3s: if shared_terminus[r] {
                d.segment_discharge_profile_m3s
                    .get(r)
                    .filter(|p| p.len() >= 2)
                    .map(|p| p[p.len() - 2])
                    .unwrap_or_else(|| q(r))
            } else {
                q(r)
            },
            width_mouth_m: d.segment_width_m.get(r).copied().unwrap_or(0.0),
            width_source_m: d.segment_width_m.get(source).copied().unwrap_or(0.0),
            // SHARED-MOUTH ATTRIBUTION (ADR Finding 47, the author's call): a system reports
            // ITS OWN catchment, never the union. When the root's terminus is also another
            // reach's terminus, that cell's accumulation is the sum of both branches — 46 of 46
            // duplicated-terminal groups at 2048² were each claiming the whole. Read the last
            // point of the root reach that is NOT the shared confluence instead; with a
            // per-point discharge profile that is one index back. A reach of a single point at
            // a shared mouth has nothing unshared to read and keeps the cell value (flagged in
            // the ADR as the residual case, ~0 at both resolutions).
            catchment_km2: {
                let own = if shared_terminus[r] {
                    d.segment_discharge_profile_m3s
                        .get(r)
                        .filter(|p| p.len() >= 2)
                        .map(|p| p[p.len() - 2])
                } else {
                    None
                };
                match own {
                    Some(q) => q / ymir_core::tectonics_c1::drainage::runoff_km2_to_m3s(300.0),
                    None => d.segment_drainage_km2.get(r).copied().unwrap_or(0.0),
                }
            },
            mouth_xy: (mx, my),
            sink,
            sink_lake_id,
            length_km: pts as f32 * km_per_cell * ratio,
            kind: d.segment_kind.get(r).copied().unwrap_or(SegmentKind::Watercourse),
            source_lake_id: d.segment_source_lake.get(r).copied().flatten(),
        });
    }
    // WATERCOURSES FIRST, then discharge. A spillway carries the whole catchment of its
    // closed basin, so on discharge alone it outranks every real river (Finding 41: the
    // author's list opened on six spillways, "the first real river is #7"). The kind is
    // the primary key so the head of the list is a river again; spillways stay browsable
    // at the tail, where their discharge is honest but their rank no longer misleads.
    out.sort_by(|a, b| {
        let key = |w: &Watercourse| u8::from(w.kind == SegmentKind::Spillway);
        key(a)
            .cmp(&key(b))
            .then(
                b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap_or(std::cmp::Ordering::Equal),
            )
            // TIE-BREAK, deterministic and meaningful: a river SYSTEM outranks a stump that
            // shares its mouth cell and therefore reads the same discharge (ADR Finding 46 —
            // two disjoint reaches ending on one cell both read that cell's accumulation, which
            // is the union of their catchments). More segments = the actual system. Then the
            // mouth cell, so the order is total and reproducible.
            .then(b.segments.len().cmp(&a.segments.len()))
            .then(a.mouth_xy.cmp(&b.mouth_xy))
    });
    out
}

/// **ADR Finding 123-C2 — the trunk is the path of MAX drained area** from the mouth to the source,
/// at every confluence and at every exorheic-lake crossing; two areas within 1 % (PROXY) tie, and
/// the tie goes to the longest-path climb (in points) this replaces. That climb followed a long
/// MINOR branch wherever one existed. Finding 93's control, measured before/after in
/// `f123_identity_control` (trunks of one reach, lengths, lake crossings): 0 trunks lose a lake
/// crossing on the delivered world or on C2/10, the one-reach share is unchanged.
const TRUNK_BY_MAX_AREA: bool = true;

/// ADR Finding 124-3 — the water BODIES: every lake id and every connected enclosed below-sea
/// component (`water_class == 2`, 8-connectivity, torus), UNIFIED where they overlap — a lake
/// flooding a below-sea basin up to its spill level is one body. `body[k]` is the canonical body of
/// cell k, `u32::MAX` on land and on the open sea. PORTED in `ymir-core/tests/f124_export.rs`.
fn water_bodies(
    d: &ymir_core::tectonics_c1::drainage::C1DrainageResult,
    wc: &[u8],
    w: usize,
    h: usize,
) -> Vec<u32> {
    let n = w * h;
    let max_lake = d.lake_map.iter().copied().max().unwrap_or(0) as usize;
    let mut comp = vec![u32::MAX; n];
    let mut nc = 0u32;
    for s in 0..n {
        if wc[s] != 2 || comp[s] != u32::MAX {
            continue;
        }
        comp[s] = nc;
        let mut st = vec![s];
        while let Some(k) = st.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nk = (y + dy).rem_euclid(h as i32) as usize * w
                        + (x + dx).rem_euclid(w as i32) as usize;
                    if wc[nk] == 2 && comp[nk] == u32::MAX {
                        comp[nk] = nc;
                        st.push(nk);
                    }
                }
            }
        }
        nc += 1;
    }
    let mut p: Vec<usize> = (0..max_lake + 1 + nc as usize).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut c = x;
        while p[c] != r {
            let nx = p[c];
            p[c] = r;
            c = nx;
        }
        r
    }
    for k in 0..n {
        if d.lake_map[k] != 0 && comp[k] != u32::MAX {
            let a = find(&mut p, d.lake_map[k] as usize);
            let b = find(&mut p, max_lake + 1 + comp[k] as usize);
            if a != b {
                p[a] = b;
            }
        }
    }
    (0..n)
        .map(|k| {
            if d.lake_map[k] != 0 {
                find(&mut p, d.lake_map[k] as usize) as u32
            } else if comp[k] != u32::MAX {
                find(&mut p, max_lake + 1 + comp[k] as usize) as u32
            } else {
                u32::MAX
            }
        })
        .collect()
}

/// The water body a point touches (its 3×3), if any.
fn body_near(body: &[u32], w: usize, h: usize, x: u32, y: u32) -> Option<u32> {
    for dy in -1i32..=1 {
        for dx in -1i32..=1 {
            let k = (y as i32 + dy).rem_euclid(h as i32) as usize * w
                + (x as i32 + dx).rem_euclid(w as i32) as usize;
            if body[k] != u32::MAX {
                return Some(body[k]);
            }
        }
    }
    None
}

/// A reach's points without the confluence cell it shares with its receiver (Finding 124-3).
fn own_points(s: &ymir_core::terrain::flow::RiverSegment) -> &[(u32, u32)] {
    if s.downstream.is_some() && s.points.len() >= 2 {
        &s.points[..s.points.len() - 1]
    } else {
        &s.points[..]
    }
}

/// **ADR Finding 124-3 — THE shipped aggregation: one object per river SYSTEM, and every cell of
/// the network in exactly one.** Three gestures on top of Findings 45 / 47 / 123:
/// - **one reach, one object** — systems whose reaches share a cell are ONE object (union of their
///   terminals): the common downstream is one trunk. The rule it replaces grouped by terminal, and
///   a spillway running down a watercourse's valley to the same mouth made two objects claiming the
///   same cells (550 cells, 19 pairs on C2/10);
/// - **the chain crosses every body that overflows** — a reach dying in a water body that has a
///   SPILLWAY is chained to it, as Finding 45 chains an exorheic lake to its outlet; the climb back
///   up is the exact inverse of that chaining (a climb rebuilt separately walked into reaches the
///   chaining gave to another object);
/// - **the basin is an AREA** — `catchment_km2` is the object's share of the GEOMETRIC partition of
///   the land (every land cell follows D8 to the first river cell or water body it meets), so the
///   objects' basins are disjoint and sum to at most the land. It was the runoff-equivalent area
///   (`runoff / 300 mm`), 1.5× the area on average, and its sum read 1.30–1.40× the land.
/// The trunk is Finding 123's max-area climb; since Finding 124-3 the core reads each reach's area
/// at its OWN end (the confluence cell is its receiver's), so the climb no longer ties the teeth of
/// a comb to their collector. PORTED in `ymir-core/tests/f124_export.rs::aggregate_v2`.
fn aggregate_watercourses(hd: &HdResult, domain_km: f32, ratio: f32) -> Vec<Watercourse> {
    let d = &hd.drainage;
    let (w, h) = (hd.width, hd.height);
    let wc = ymir_core::lakes::connectivity::water_class(&hd.eroded, 0.5);
    let body = water_bodies(d, &wc, w, h);
    let segs = &d.rivers.segments;
    let n = segs.len();
    let q = |i: usize| d.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
    let is_spill = |i: usize| d.segment_kind.get(i) == Some(&SegmentKind::Spillway);
    // a spillway's body: where it starts, else the lake it names
    let mut body_of_lake: std::collections::HashMap<u32, u32> = Default::default();
    for k in 0..w * h {
        if d.lake_map[k] != 0 {
            body_of_lake.entry(d.lake_map[k]).or_insert(body[k]);
        }
    }
    let mut outlet_of: std::collections::HashMap<u32, usize> = Default::default();
    let mut spill_of: std::collections::HashMap<u32, usize> = Default::default();
    for (i, sg) in segs.iter().enumerate() {
        let (sx, sy) = sg.points[0];
        if is_spill(i) {
            let b = body_near(&body, w, h, sx, sy).or_else(|| {
                d.segment_source_lake.get(i).copied().flatten().and_then(|id| body_of_lake.get(&id).copied())
            });
            if let Some(b) = b {
                if spill_of.get(&b).is_none_or(|&j| q(i) > q(j)) {
                    spill_of.insert(b, i);
                }
            }
        } else if let (Sink::ExoLake, Some(id)) = classify_sink(hd, &wc, sx, sy) {
            if outlet_of.get(&id).is_none_or(|&j| q(i) > q(j)) {
                outlet_of.insert(id, i);
            }
        }
    }
    let end_body = |i: usize| -> Option<u32> {
        let &(x, y) = segs[i].points.last()?;
        body_near(&body, w, h, x, y)
    };
    // the continuation of a terminal reach: the spillway of the body it dies in, else (Finding 45)
    // the outlet of the exorheic lake it dies in
    let across = |i: usize| -> Option<usize> {
        if let Some(&s) = end_body(i).and_then(|b| spill_of.get(&b)) {
            return Some(s).filter(|&o| o != i);
        }
        let &(mx, my) = segs[i].points.last()?;
        match classify_sink(hd, &wc, mx, my) {
            (Sink::ExoLake, Some(id)) => outlet_of.get(&id).copied().filter(|&o| o != i),
            _ => None,
        }
    };
    let mut root = vec![0usize; n];
    for i in 0..n {
        let (mut j, mut seen) = (i, 0usize);
        loop {
            seen += 1;
            if seen > n {
                break; // a cycle through bodies: stop, never loop
            }
            match segs[j].downstream {
                Some(k) if k < n => j = k,
                _ => match across(j) {
                    Some(o) => j = o,
                    None => break,
                },
            }
        }
        root[i] = j;
    }
    // one reach, one object: the terminals of reaches sharing a cell are unified
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut c = x;
        while p[c] != r {
            let nx = p[c];
            p[c] = r;
            c = nx;
        }
        r
    }
    let mut owner: std::collections::HashMap<usize, usize> = Default::default();
    for i in 0..n {
        for &(x, y) in own_points(&segs[i]) {
            let k = y as usize * w + x as usize;
            match owner.get(&k) {
                Some(&o) => {
                    let (a, b) = (find(&mut parent, root[o]), find(&mut parent, root[i]));
                    if a != b {
                        parent[a] = b;
                    }
                }
                None => {
                    owner.insert(k, i);
                }
            }
        }
    }
    let mut groups: std::collections::HashMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        let g = find(&mut parent, root[i]);
        groups.entry(g).or_default().push(i);
    }
    let km_per_cell = if w > 0 { domain_km / w as f32 } else { 0.0 };
    let mut pathlen = vec![0u32; n];
    {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|&i| segs[i].strahler_order);
        for _ in 0..16 {
            let mut changed = false;
            for &i in &order {
                let up = segs[i]
                    .upstream
                    .iter()
                    .copied()
                    .filter(|&u| u < n)
                    .map(|u| pathlen[u])
                    .max()
                    .unwrap_or(0);
                let v = up + segs[i].points.len() as u32;
                if v != pathlen[i] {
                    pathlen[i] = v;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    // Finding 123-C2's max-area climb and its 1 % tie (PROXY), on the reach's OWN area (core, F124-3)
    let area_of = |i: usize| d.segment_drainage_km2.get(i).copied().unwrap_or(0.0);
    let by_area = |a: usize, b: usize| -> std::cmp::Ordering {
        let (x, y) = (area_of(a), area_of(b));
        if (x - y).abs() <= 0.01 * x.max(y) { pathlen[a].cmp(&pathlen[b]) } else { x.total_cmp(&y) }
    };
    // the climb across a body: the exact inverse of `across`
    let mut chained_to: std::collections::HashMap<usize, Vec<usize>> = Default::default();
    for j in 0..n {
        if segs[j].downstream.is_none() {
            if let Some(o) = across(j) {
                chained_to.entry(o).or_default().push(j);
            }
        }
    }
    // Finding 47's shared-mouth attribution: a terminal on a cell another terminal ends on reads
    // one point back
    let mut terminus_count: std::collections::HashMap<(u32, u32), usize> = Default::default();
    for sg in segs.iter().filter(|sg| sg.downstream.is_none()) {
        *terminus_count.entry(*sg.points.last().unwrap()).or_insert(0) += 1;
    }
    let own_q = |r: usize| -> f32 {
        let shared = terminus_count.get(segs[r].points.last().unwrap()).copied().unwrap_or(0) > 1;
        if shared {
            d.segment_discharge_profile_m3s
                .get(r)
                .filter(|p| p.len() >= 2)
                .map(|p| p[p.len() - 2])
                .unwrap_or_else(|| q(r))
        } else {
            q(r)
        }
    };
    let mut group_list: Vec<(usize, Vec<usize>)> = groups.into_iter().collect();
    group_list.sort_by_key(|(r, _)| *r); // determinism: never the HashMap's order
    let mut out = Vec::with_capacity(group_list.len());
    let mut seg_obj = vec![u32::MAX; n];
    for (o, (_, members)) in group_list.iter().enumerate() {
        for &s in members {
            seg_obj[s] = o as u32;
        }
    }
    for (_, members) in group_list {
        let mut roots: Vec<usize> = members.iter().map(|&i| root[i]).collect();
        roots.sort_unstable();
        roots.dedup();
        let rep = *roots.iter().max_by(|&&a, &&b| by_area(a, b)).unwrap();
        let mut trunk = Vec::new();
        let mut cur = rep;
        let mut visited = std::collections::HashSet::new();
        while visited.insert(cur) {
            trunk.push(cur);
            let next = segs[cur].upstream.iter().copied().filter(|&u| u < n).max_by(|&a, &b| by_area(a, b));
            match next {
                Some(nx) => cur = nx,
                None => match chained_to.get(&cur).and_then(|v| v.iter().copied().max_by(|&a, &b| by_area(a, b))) {
                    Some(up) => cur = up,
                    None => break,
                },
            }
        }
        let source = *trunk.last().unwrap();
        trunk.reverse(); // source → sink
        let path_cells = trunk_length_cells(segs, &trunk, w, h);
        let (mx, my) = *segs[rep].points.last().unwrap();
        let (sink, sink_lake_id) = classify_sink(hd, &wc, mx, my);
        // a system is a RIVER as soon as one of its reaches is a watercourse (a spillway running
        // down a valley, or a basin fed by rivers); its order is then a watercourse's
        let river = members.iter().any(|&s| !is_spill(s));
        let order = if !is_spill(rep) {
            segs[rep].strahler_order
        } else {
            members.iter().filter(|&&s| !is_spill(s)).map(|&s| segs[s].strahler_order).max().unwrap_or(1)
        };
        out.push(Watercourse {
            trunk,
            tributaries: members.len().saturating_sub(1),
            segments: members,
            order,
            discharge_m3s: roots.iter().map(|&r| own_q(r)).sum(),
            width_mouth_m: d.segment_width_m.get(rep).copied().unwrap_or(0.0),
            width_source_m: d.segment_width_m.get(source).copied().unwrap_or(0.0),
            catchment_km2: 0.0, // the geometric partition, below
            mouth_xy: (mx, my),
            sink,
            sink_lake_id,
            length_km: path_cells * km_per_cell * ratio,
            kind: if river { SegmentKind::Watercourse } else { SegmentKind::Spillway },
            source_lake_id: d.segment_source_lake.get(rep).copied().flatten(),
        });
    }
    // ── the geometric partition: every land cell takes the object of the first river cell or water
    //    body its D8 path meets; a body belongs to the object of the reach LEAVING it (outlet or
    //    spillway), else of the largest reach dying in it. The open sea and the coastal fringe no
    //    river captures belong to no object.
    let cells = catchment_partition(d, &wc, &body, &seg_obj, out.len());
    let cell_km2 = km_per_cell * (domain_km / h.max(1) as f32);
    for (o, x) in out.iter_mut().enumerate() {
        x.catchment_km2 = cells[o] as f32 * cell_km2 * ratio * ratio;
    }
    out.sort_by(|a, b| {
        let key = |w: &Watercourse| u8::from(w.kind == SegmentKind::Spillway);
        key(a)
            .cmp(&key(b))
            .then(b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap_or(std::cmp::Ordering::Equal))
            .then(b.segments.len().cmp(&a.segments.len()))
            .then(a.mouth_xy.cmp(&b.mouth_xy))
    });
    out
}

/// ADR Finding 124-3 — the LENGTH of a trunk, in cells: the path through its reaches' points in
/// order, each step Euclidean on the torus (1 or √2 on D8, a crossed water body by its chord), the
/// confluence cell a reach shares with its receiver counted ONCE. It was `points × cell`, which
/// counted every confluence twice (+7.6 % on a 61 km trunk) and a diagonal step as one cell.
fn trunk_length_cells(segs: &[ymir_core::terrain::flow::RiverSegment], trunk: &[usize], w: usize, h: usize) -> f32 {
    let (mut len, mut prev): (f32, Option<(u32, u32)>) = (0.0, None);
    for &s in trunk {
        for &p in &segs[s].points {
            if let Some(q) = prev {
                if q != p {
                    let dx = (p.0 as i64 - q.0 as i64).rem_euclid(w as i64);
                    let dy = (p.1 as i64 - q.1 as i64).rem_euclid(h as i64);
                    let (dx, dy) = (dx.min(w as i64 - dx) as f32, dy.min(h as i64 - dy) as f32);
                    len += (dx * dx + dy * dy).sqrt();
                }
            }
            prev = Some(p);
        }
    }
    len
}

/// ADR Finding 124-3 — cells per object of the GEOMETRIC partition of the land (see
/// `aggregate_watercourses`). Counted in CELLS (u64): 67 M f32 increments of a 0.0024 km² cell lose
/// the sum. PORTED in `ymir-core/tests/f124_export.rs::geometric_partition`.
fn catchment_partition(
    d: &ymir_core::tectonics_c1::drainage::C1DrainageResult,
    wc: &[u8],
    body: &[u32],
    seg_obj: &[u32],
    n_obj: usize,
) -> Vec<u64> {
    let (w, h) = (d.width, d.height);
    let n = w * h;
    let segs = &d.rivers.segments;
    const UNK: u32 = u32::MAX;
    const NONE: u32 = u32::MAX - 1;
    let mut obj = vec![UNK; n];
    // a body: the object leaving it (outlet / spillway, largest discharge), else the largest dying in it
    let mut body_obj: std::collections::HashMap<u32, (bool, f32, u32)> = Default::default();
    for (i, s) in segs.iter().enumerate() {
        let a = d.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
        let (sx, sy) = s.points[0];
        if let Some(b) = body_near(body, w, h, sx, sy) {
            let e = body_obj.entry(b).or_insert((true, a, seg_obj[i]));
            if !e.0 || a > e.1 {
                *e = (true, a, seg_obj[i]);
            }
        }
        if s.downstream.is_none() {
            let &(mx, my) = s.points.last().unwrap();
            if let Some(b) = body_near(body, w, h, mx, my) {
                let e = body_obj.entry(b).or_insert((false, a, seg_obj[i]));
                if !e.0 && a > e.1 {
                    *e = (false, a, seg_obj[i]);
                }
            }
        }
    }
    for (i, s) in segs.iter().enumerate() {
        for &(x, y) in own_points(s) {
            let k = y as usize * w + x as usize;
            if obj[k] == UNK {
                obj[k] = seg_obj[i];
            }
        }
    }
    for k in 0..n {
        if obj[k] == UNK {
            if body[k] != u32::MAX {
                obj[k] = body_obj.get(&body[k]).map_or(NONE, |e| e.2);
            } else if wc[k] == 1 {
                obj[k] = NONE;
            }
        }
    }
    let dir = &d.flow.direction;
    let mut path = Vec::new();
    for s in 0..n {
        if obj[s] != UNK {
            continue;
        }
        path.clear();
        let mut cur = s;
        let res = loop {
            if obj[cur] != UNK {
                break obj[cur];
            }
            path.push(cur);
            let dd = dir[cur];
            if dd == DIR_NONE || path.len() > n {
                break NONE;
            }
            let (x, y) = ((cur % w) as i32, (cur / w) as i32);
            cur = (y + D8_DY[dd as usize]).rem_euclid(h as i32) as usize * w
                + (x + D8_DX[dd as usize]).rem_euclid(w as i32) as usize;
        };
        for &k in &path {
            obj[k] = res;
        }
    }
    let mut cells = vec![0u64; n_obj];
    for k in 0..n {
        if wc[k] != 1 && (obj[k] as usize) < n_obj {
            cells[obj[k] as usize] += 1;
        }
    }
    cells
}

fn sink_label(s: Sink) -> (&'static str, C) {
    match s {
        Sink::Sea => ("→ mer", C::from_rgb(0x5a, 0x9a, 0xd0)),
        Sink::ExoLake => ("→ lac exoréique", C::from_rgb(0x5a, 0x9a, 0xd0)),
        Sink::EndoLake => ("→ bassin endoréique", C::from_rgb(0x3a, 0xb0, 0xa0)),
        // ADR Finding 91-D / 92-D — the condition is `!sea && water_class == 2`: an enclosed
        // below-sea cell that NO water body covers. It never reads `net_evap`, `a_eq` or
        // `lake_type`, so "évaporatif" asserted a process it cannot see — and in a humid bed
        // `net_evap = 0 ⇒ a_eq = ∞` (Finding 39), so nothing there evaporates. Since Finding 92-C
        // every enclosed component carries a body, so what remains in this branch is a mouth on a
        // cell outside every inventory: Finding 89's `to_nothing`, named as what it is.
        Sink::SubSeaSink => {
            ("→ cuvette sous-marine SANS exutoire tracé", C::from_rgb(0xd0, 0x6a, 0x5a))
        }
        Sink::Unknown => ("→ ?", WARN_ORANGE),
    }
}

/// Microscope SELECTION drawer (left, Finding 28): entity tabs (rivers as watercourses / lakes)
/// + a selectable list. Click a selected row again to DESELECT it. Read-only.
fn microscope_list(ui: &mut egui::Ui, ws: &mut WorkspaceState) {
    let hd = match ws.current.clone() {
        Some(h) => h,
        None => return,
    };
    let wcs = ws.watercourses.clone();
    {
        // ADR Finding 123 — the identity guard, stated before any number.
        use ymir_core::tectonics_c1::bench_guard::GuardStatus;
        let (txt, col, tip) = match &hd.bench_guard {
            GuardStatus::Match { label } => (
                format!("= banc ({label})"),
                C::from_rgb(0x5a, 0xc0, 0x7a),
                "Ce monde est BIT À BIT celui qu'un banc a mesuré à ce réglage: ses nombres se \
                 lisent contre le Finding correspondant."
                    .to_string(),
            ),
            GuardStatus::Mismatch { label, bench, viz } => (
                format!("≠ banc ({label}) — nombres masqués"),
                C::from_rgb(0xe0, 0x60, 0x50),
                format!(
                    "Même réglage que le banc « {label} », mais le champ DIFFÈRE (banc {bench}, \
                     viz {viz}). Aucun nombre n'est affiché: ce monde n'est pas le monde mesuré."
                ),
            ),
            GuardStatus::NoReference { viz } => (
                "non gardé — aucun banc à ce réglage".to_string(),
                C::from_rgb(0xd0, 0xb0, 0x50),
                format!(
                    "Aucun banc n'a mesuré ce réglage exact (seed, résolution, cadrage, états). Les \
                     nombres sont ceux de CE monde et ne se comparent à aucun Finding. Champ {viz}."
                ),
            ),
        };
        ui.label(egui::RichText::new(txt).color(col).size(10.5)).on_hover_text(tip);
        // ADR Finding 134 — the lakes, guarded apart: the final lakes (after the crater pass) against
        // the bench's, under the final drainage's key (field + drainage + climate). One climate is
        // guarded, 45° / 40°.
        let (txt, col) = match &hd.lake_guard {
            GuardStatus::Match { .. } => ("lacs = banc".to_string(), C::from_rgb(0x5a, 0xc0, 0x7a)),
            GuardStatus::Mismatch { bench, viz, .. } => (
                format!("lacs ≠ banc (banc {bench}, viz {viz}) — nombres masqués"),
                C::from_rgb(0xe0, 0x60, 0x50),
            ),
            GuardStatus::NoReference { .. } => (
                "lacs non gardés — un seul climat gardé (45° / 40°)".to_string(),
                C::from_rgb(0xd0, 0xb0, 0x50),
            ),
        };
        ui.label(egui::RichText::new(txt).color(col).size(10.5)).on_hover_text(
            "Les lacs finals (après la passe des lacs de cratère) comparés à ceux du banc, sous la \
             clé du drainage final (champ + drainage + climat). Un seul climat est gardé: \
             latitude 45°, étendue 40°.",
        );
    }
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("MICROSCOPE").color(TEXT_BRIGHT).strong().size(11.5));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui.button("‹").on_hover_text("Replier le microscope").clicked() {
                ws.micro_open = false;
            }
        });
    });
    ui.horizontal(|ui| {
        if ui
            .selectable_label(
                ws.inspect_tab == InspectTab::Rivers,
                format!("Rivières ({})", wcs.as_ref().map(|w| w.len()).unwrap_or(0)),
            )
            .clicked()
        {
            ws.inspect_tab = InspectTab::Rivers;
        }
        if ui
            .selectable_label(
                ws.inspect_tab == InspectTab::Lakes,
                format!("Lacs ({})", hd.drainage.lakes.len()),
            )
            .clicked()
        {
            ws.inspect_tab = InspectTab::Lakes;
        }
        if !hd.volcanoes.is_empty()
            && ui
                .selectable_label(
                    ws.inspect_tab == InspectTab::Volcanoes,
                    format!("Volcans ({})", hd.volcanoes.len()),
                )
                .clicked()
        {
            ws.inspect_tab = InspectTab::Volcanoes;
        }
    });
    ui.separator();
    egui::ScrollArea::vertical().id_salt("entity_list").show(ui, |ui| match ws.inspect_tab {
        InspectTab::Rivers => {
            if let Some(wcs) = &wcs {
                for (i, wc) in wcs.iter().enumerate() {
                    let (sl, sc) = sink_label(wc.sink);
                    // A spillway has no Strahler order (emitted as 1 regardless of
                    // catchment) and no tributaries — show what it IS instead.
                    let txt = if wc.kind == SegmentKind::Spillway {
                        format!(
                            "#{:<3} déversoir · {:.0} m³/s · bassin {} · {}",
                            i + 1,
                            wc.discharge_m3s,
                            wc.source_lake_id
                                .map(|id| id.to_string())
                                .unwrap_or_else(|| "hors inventaire".into()),
                            sl
                        )
                    } else {
                        format!(
                            "#{:<3} S{} · {:.0} m³/s · {} trib · {}",
                            i + 1,
                            wc.order,
                            wc.discharge_m3s,
                            wc.tributaries,
                            sl
                        )
                    };
                    let r = ui.selectable_label(
                        ws.selected_river == Some(i),
                        egui::RichText::new(txt).color(sc).size(11.0),
                    );
                    if r.double_clicked() {
                        // Double-click: select AND frame the river on the map.
                        ws.selected_river = Some(i);
                        if let Some(bb) = river_bbox_uv(&hd, wc) {
                            ws.focus_target = Some(bb);
                        }
                    } else if r.clicked() {
                        // Single click toggles (re-click deselects).
                        ws.selected_river = (ws.selected_river != Some(i)).then_some(i);
                    }
                }
            }
        }
        InspectTab::Lakes => {
            let mut idx: Vec<usize> = (0..hd.drainage.lakes.len()).collect();
            idx.sort_by(|&a, &b| {
                hd.drainage.lakes[b]
                    .area_km2
                    .partial_cmp(&hd.drainage.lakes[a].area_km2)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            for i in idx {
                let lk = &hd.drainage.lakes[i];
                let (ty, tc) = match lk.lake_type {
                    LakeType::Exorheic => ("exoréique", C::from_rgb(0x5a, 0x9a, 0xd0)),
                    LakeType::Endorheic => ("endoréique (salé)", C::from_rgb(0x3a, 0xb0, 0xa0)),
                    LakeType::CraterAcidic => ("cratère (acide)", C::from_rgb(0xc9, 0xc0, 0x3a)),
                    LakeType::CraterNeutral => {
                        ("cratère (eau douce)", C::from_rgb(0x6a, 0x8a, 0xc0))
                    }
                    // ADR Finding 86 — drawn in warning red, because it is a disagreement between
                    // the water balance and the traced network, not a kind of lake.
                    LakeType::Unresolved => (
                        unresolved_reason_short_fr(lk.unresolved_reason),
                        C::from_rgb(0xd0, 0x6a, 0x5a),
                    ),
                };
                let txt = format!(
                    "#{} · {:.0} km² · {:.0} m · {}",
                    lk.base.id, lk.area_km2, lk.level_m, ty
                );
                let id = lk.base.id;
                let r = ui.selectable_label(
                    ws.selected_lake == Some(i),
                    egui::RichText::new(txt).color(tc).size(11.0),
                );
                if r.double_clicked() {
                    ws.selected_lake = Some(i);
                    if let Some(bb) = lake_bbox_uv(&hd, id) {
                        ws.focus_target = Some(bb);
                    }
                } else if r.clicked() {
                    ws.selected_lake = (ws.selected_lake != Some(i)).then_some(i);
                }
            }
        }
        InspectTab::Volcanoes => {
            let km_per_cell = ws.domain_km / hd.width as f32;
            // Active (degassing) first, then by crater size.
            let mut idx: Vec<usize> = (0..hd.volcanoes.len()).collect();
            idx.sort_by(|&a, &b| {
                let (va, vb) = (&hd.volcanoes[a], &hd.volcanoes[b]);
                vb.active.cmp(&va.active).then(
                    vb.radius_px.partial_cmp(&va.radius_px).unwrap_or(std::cmp::Ordering::Equal),
                )
            });
            for i in idx {
                let v = &hd.volcanoes[i];
                let (state, col) = if v.active {
                    ("ACTIF (dégazage)", C::from_rgb(0xe0, 0x62, 0x40))
                } else {
                    ("éteint", C::from_rgb(0x9a, 0x9a, 0x9a))
                };
                let dia_km = 2.0 * v.radius_px * km_per_cell;
                let txt = format!("#{:<2} {state} · cratère Ø {dia_km:.1} km", i + 1);
                let r = ui.selectable_label(
                    ws.selected_volcano == Some(i),
                    egui::RichText::new(txt).color(col).size(11.0),
                );
                if r.double_clicked() {
                    ws.selected_volcano = Some(i);
                    // Frame the volcano. center_px is in DATA space; the camera uv is
                    // TEXTURE space (north-up flipped), so reuse bbox_to_uv (which flips
                    // y via h-1-y) with a small box around the crater centre — passing
                    // data y directly would land at the vertically mirrored spot.
                    let (cx, cy) = (v.center_px.0 as u32, v.center_px.1 as u32);
                    let r_px = (v.radius_px.max(4.0) * 2.0) as u32;
                    let ([bu, bv], ext) = bbox_to_uv(
                        cx.saturating_sub(r_px),
                        (cx + r_px).min(hd.width as u32 - 1),
                        cy.saturating_sub(r_px),
                        (cy + r_px).min(hd.height as u32 - 1),
                        hd.width,
                        hd.height,
                    );
                    ws.focus_target = Some(([bu, bv], ext));
                } else if r.clicked() {
                    ws.selected_volcano = (ws.selected_volcano != Some(i)).then_some(i);
                }
            }
        }
    });
}

/// Microscope DETAIL (right inspection zone, Finding 28): the selected river long profile OR
/// lake sheet. The clickable hydrological chain (Finding 29) jumps between entities.
fn microscope_detail(ui: &mut egui::Ui, ws: &mut WorkspaceState) {
    let hd = match ws.current.clone() {
        Some(h) => h,
        None => return,
    };
    let wcs = ws.watercourses.clone();
    // Geo-ratio note: displayed discharges/areas are ×ratio² compressed vs the real balance.
    if (ws.geo_scale_ratio - 1.0).abs() > 0.01 {
        ui.label(
            egui::RichText::new(format!(
                "débits ×{:.0} compressés (échelle {:.1})",
                ws.geo_scale_ratio * ws.geo_scale_ratio,
                ws.geo_scale_ratio
            ))
            .color(WARN_ORANGE)
            .size(10.0),
        );
    }
    let nav = match ws.inspect_tab {
        InspectTab::Rivers => match (ws.selected_river, &wcs) {
            (Some(i), Some(w)) if i < w.len() => river_profile_panel(ui, ws, &hd, &w[i]),
            _ => {
                ui.add_space(20.0);
                ui.label(
                    egui::RichText::new("Sélectionnez un cours d'eau (microscope à gauche).")
                        .color(DIM)
                        .size(11.0),
                );
                None
            }
        },
        InspectTab::Lakes => match ws.selected_lake {
            Some(i) if i < hd.drainage.lakes.len() => {
                lake_sheet_panel(ui, &hd, i, ws.domain_km, wcs.as_deref().map(|v| v.as_slice()))
            }
            _ => {
                ui.add_space(20.0);
                ui.label(
                    egui::RichText::new("Sélectionnez un lac (microscope à gauche).")
                        .color(DIM)
                        .size(11.0),
                );
                None
            }
        },
        InspectTab::Volcanoes => match ws.selected_volcano.and_then(|i| hd.volcanoes.get(i)) {
            Some(v) => {
                ui.add_space(8.0);
                let (state, col) = if v.active {
                    ("volcan ACTIF (dégazage)", C::from_rgb(0xe0, 0x62, 0x40))
                } else {
                    ("volcan éteint", C::from_rgb(0xb0, 0xb0, 0xb0))
                };
                ui.label(egui::RichText::new(state).color(col).strong().size(13.0));
                let dia_km = 2.0 * v.radius_px * (ws.domain_km / hd.width as f32);
                for line in [
                    format!("cratère : Ø {dia_km:.1} km"),
                    if v.active {
                        "eau de cratère : acide (pH < 2) si un lac s'y forme — sinon relief sec"
                    } else {
                        "eau de cratère : eau douce si un lac s'y forme — sinon relief sec"
                    }
                    .to_string(),
                ] {
                    ui.label(egui::RichText::new(line).color(DIM).size(11.0));
                }
                None
            }
            None => {
                ui.add_space(20.0);
                ui.label(
                    egui::RichText::new("Sélectionnez un volcan (microscope à gauche).")
                        .color(DIM)
                        .size(11.0),
                );
                None
            }
        },
    };
    // Clickable-chain jump (Finding 29): follow the hydrological graph in one click.
    match nav {
        Some(NavAction::Lake(i)) => {
            ws.selected_lake = Some(i);
            ws.inspect_tab = InspectTab::Lakes;
        }
        Some(NavAction::River(i)) => {
            ws.selected_river = Some(i);
            ws.inspect_tab = InspectTab::Rivers;
        }
        None => {}
    }
}

/// TASK 3 — river long profile: bed elevation source→sink (from `profile_m`) + the figures
/// the inspector already carries. The SINK is marked explicitly; any climb is flagged (this
/// doubles as the monotonicity inspector). Reads `segment_profile_m` / discharge / width.
fn river_profile_panel(
    ui: &mut egui::Ui,
    ws: &mut WorkspaceState,
    hd: &HdResult,
    wc: &Watercourse,
) -> Option<NavAction> {
    let d = &hd.drainage;
    let (sl, sc) = sink_label(wc.sink);
    let mut nav = None;
    // SOURCE basin (Finding 29, upstream side): if the watercourse's source cell borders a lake, that
    // lake OVERFLOWS into this river — expose it, clickable, to walk the hydrological chain UPSTREAM.
    let source_lake = wc
        .trunk
        .first()
        .and_then(|&s| d.rivers.segments.get(s))
        .and_then(|seg| seg.points.first())
        .and_then(|&(sx, sy)| {
            let (w, h) = (hd.width, hd.height);
            let mut found = 0u32;
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let (x, y) = (sx as i32 + dx, sy as i32 + dy);
                    if x >= 0 && y >= 0 && (x as usize) < w && (y as usize) < h {
                        let id = d.lake_map[y as usize * w + x as usize];
                        if id != 0 {
                            found = id;
                        }
                    }
                }
            }
            (found != 0).then_some(found)
        })
        .and_then(|id| hd.drainage.lakes.iter().position(|l| l.base.id == id));
    ui.horizontal_wrapped(|ui| {
        let kv = |ui: &mut egui::Ui, k: &str, v: String| {
            ui.label(egui::RichText::new(format!("{k} ")).color(DIM).size(10.5));
            ui.label(egui::RichText::new(v).color(COPPER_BRIGHT).monospace().size(11.0));
            ui.add_space(10.0);
        };
        kv(ui, "Débit", format!("{:.0} m³/s", wc.discharge_m3s));
        kv(ui, "Bassin", format!("{:.0} km²", wc.catchment_km2));
        kv(ui, "Longueur", format!("{:.0} km", wc.length_km));
        kv(ui, "Ordre", format!("S{}", wc.order));
        kv(ui, "Largeur embouchure", format!("{:.0} m", wc.width_mouth_m));
        // Finding 37 — honest label: the first exported point sits at `stream_km2` (~20 km²), NOT the
        // true source, until the upstream extension (head_km2 = A_c) is enabled. Then it IS the source.
        kv(ui, "Largeur au 1er point", format!("{:.0} m", wc.width_source_m));
        // Clickable SOURCE basin (Finding 29, upstream): the lake that overflows into this river.
        if let Some(li) = source_lake {
            let lk = &hd.drainage.lakes[li];
            let (ty, tc) = match lk.lake_type {
                LakeType::Exorheic => ("exoréique", C::from_rgb(0x5a, 0x9a, 0xd0)),
                LakeType::Endorheic => ("endoréique", C::from_rgb(0x3a, 0xb0, 0xa0)),
                LakeType::CraterAcidic => ("cratère (acide)", C::from_rgb(0xc9, 0xc0, 0x3a)),
                LakeType::CraterNeutral => ("cratère (eau douce)", C::from_rgb(0x6a, 0x8a, 0xc0)),
                // ADR Finding 86
                LakeType::Unresolved => (
                    unresolved_reason_short_fr(lk.unresolved_reason),
                    C::from_rgb(0xd0, 0x6a, 0x5a),
                ),
            };
            ui.label(egui::RichText::new("Source ").color(DIM).size(10.5));
            if ui
                .add(
                    egui::Button::new(
                        egui::RichText::new(format!("⬅ lac #{} ({ty})", lk.base.id))
                            .color(tc)
                            .strong()
                            .size(11.5),
                    )
                    .frame(true),
                )
                .on_hover_text("Aller au bassin source (remonter la chaîne)")
                .clicked()
            {
                nav = Some(NavAction::Lake(li));
            }
            ui.add_space(10.0);
        }
        // Clickable SINK (Finding 29): jump to the lake it drains into. A SUB-threshold basin
        // (Finding 32) is a real sink but not in the inventory list → labelled honestly, no jump.
        ui.label(egui::RichText::new("Exutoire ").color(DIM).size(10.5));
        match wc.sink_lake_id.and_then(|id| hd.drainage.lakes.iter().position(|l| l.base.id == id))
        {
            Some(li) => {
                if ui
                    .add(
                        egui::Button::new(
                            egui::RichText::new(format!("{sl} ⮕")).color(sc).strong().size(11.5),
                        )
                        .frame(true),
                    )
                    .on_hover_text("Aller au lac")
                    .clicked()
                {
                    nav = Some(NavAction::Lake(li));
                }
            }
            None if wc.sink_lake_id.is_some() => {
                ui.label(
                    egui::RichText::new(format!("{sl} (sous-seuil, non listé)"))
                        .color(sc)
                        .size(11.0),
                );
            }
            None => {
                ui.label(egui::RichText::new(sl).color(sc).strong().size(11.5));
            }
        }
    });
    // Bed profile source → sink. For a normal river we WALK the flow field from the main-stem
    // headwater (the breached field is monotone, so this reflects the true bed — no dependence on
    // the clip's segment links, which would inject phantom steps at junctions; reads
    // `flow.direction` + `eroded`, both exported). But a below-sea SPILLWAY is the one watercourse
    // whose bed is NOT a flow-field streamline: it crosses a divide the flow field itself routes
    // BACK into the lake, so re-walking gives a flat lake-level profile (the "49 m → 49 m" the author
    // saw). A spillway is an appended segment (`max_flow == 0`); plot its exported `profile_m` (the
    // real spillway bed, descending to the sink) instead. Finding 39.
    let (w, h) = (hd.width, hd.height);
    let ss = SteinSteinParams::default();
    let is_spillway =
        wc.trunk.first().and_then(|&s| d.rivers.segments.get(s)).is_some_and(|s| s.max_flow == 0.0);
    let mut elev: Vec<f32> = Vec::new();
    // Cells parallel to `elev` (source→sink) so a profile-hover can place a point on the map.
    let mut path: Vec<(u32, u32)> = Vec::new();
    if is_spillway {
        for &s in &wc.trunk {
            if let (Some(pr), Some(seg)) = (d.segment_profile_m.get(s), d.rivers.segments.get(s)) {
                elev.extend(pr.iter().copied());
                path.extend(seg.points.iter().copied());
            }
        }
    } else if let Some(&(sx, sy)) =
        wc.trunk.first().and_then(|&s| d.rivers.segments[s].points.first())
    {
        let mut k = sy as usize * w + sx as usize;
        for _ in 0..(w + h) {
            elev.push(c1_altitude_norm_to_metres(hd.eroded.data[k], &ss));
            path.push(((k % w) as u32, (k / w) as u32));
            if hd.eroded.data[k] <= 0.5 || d.lake_map[k] != 0 {
                break; // reached the sink (sea or lake)
            }
            let dir = d.flow.direction[k];
            if dir == DIR_NONE {
                break;
            }
            let (nx, ny) =
                ((k % w) as i32 + D8_DX[dir as usize], (k / w) as i32 + D8_DY[dir as usize]);
            if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                break;
            }
            k = ny as usize * w + nx as usize;
        }
    }
    // Share the path with the map (it draws the tracking point on the highlight). Truncated to the
    // longer of the two so `path`/`elev` stay parallel; cleared hover until the plot is hovered.
    ws.profile_path = path;
    ws.profile_hover = None;
    if elev.len() < 2 {
        ui.label(egui::RichText::new("Profil indisponible.").color(DIM).size(11.0));
        return nav;
    }
    let (mut lo, mut hi) = (f32::MAX, f32::MIN);
    for &e in &elev {
        lo = lo.min(e);
        hi = hi.max(e);
    }
    let mut max_climb = 0.0f32;
    for w in elev.windows(2) {
        max_climb = max_climb.max(w[1] - w[0]);
    }
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("Profil en long (source → exutoire)").color(DIM).size(10.5));
        let (mtxt, mcol) = if max_climb > 1.0 {
            (format!("⚠ remontée +{max_climb:.0} m"), WARN_ORANGE)
        } else {
            ("monotone ✓".to_string(), OK_GREEN)
        };
        ui.label(egui::RichText::new(mtxt).color(mcol).size(10.5));
    });
    let (rect, plot_resp) = ui.allocate_exact_size(
        egui::vec2(ui.available_width().max(120.0), 120.0),
        egui::Sense::hover(),
    );
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 3.0, C::from_rgb(0x14, 0x14, 0x18));
    let span_e = (hi - lo).max(1.0);
    let n = elev.len();
    let x_of = |i: usize| rect.left() + 6.0 + (i as f32 / (n - 1) as f32) * (rect.width() - 12.0);
    let y_of = |e: f32| rect.bottom() - 6.0 - ((e - lo) / span_e) * (rect.height() - 16.0);
    // sea level (0 m) reference line if in range.
    if lo <= 0.0 && hi >= 0.0 {
        let y0 = y_of(0.0);
        p.line_segment(
            [egui::pos2(rect.left(), y0), egui::pos2(rect.right(), y0)],
            egui::Stroke::new(1.0, C::from_rgba_unmultiplied(0x5a, 0x9a, 0xd0, 90)),
        );
    }
    let pts: Vec<egui::Pos2> = (0..n).map(|i| egui::pos2(x_of(i), y_of(elev[i]))).collect();
    p.add(egui::Shape::line(pts.clone(), egui::Stroke::new(1.5, COPPER_BRIGHT)));
    // Source (start) + sink (end) markers.
    p.circle_filled(pts[0], 3.0, C::from_rgb(0x9a, 0x9a, 0x9a));
    p.circle_filled(pts[n - 1], 4.0, sc);
    p.text(
        pts[0] + egui::vec2(2.0, -4.0),
        egui::Align2::LEFT_BOTTOM,
        "source",
        egui::FontId::proportional(9.0),
        DIM,
    );
    p.text(
        pts[n - 1] + egui::vec2(-2.0, -4.0),
        egui::Align2::RIGHT_BOTTOM,
        sl,
        egui::FontId::proportional(9.0),
        sc,
    );
    p.text(
        egui::pos2(rect.left() + 4.0, rect.top() + 2.0),
        egui::Align2::LEFT_TOP,
        format!("{hi:.0} m"),
        egui::FontId::monospace(9.0),
        DIM,
    );
    p.text(
        egui::pos2(rect.left() + 4.0, rect.bottom() - 2.0),
        egui::Align2::LEFT_BOTTOM,
        format!("{lo:.0} m"),
        egui::FontId::monospace(9.0),
        DIM,
    );
    // Hover readout: abscissa (length from the source), a vertical cursor bar, a point tracking the
    // profile, and the interpolated bed height at that abscissa.
    let (x0, xw) = (rect.left() + 6.0, (rect.width() - 12.0).max(1.0));
    if let Some(pos) = plot_resp.hover_pos() {
        let fx = ((pos.x - x0) / xw).clamp(0.0, 1.0);
        let fi = fx * (n - 1) as f32;
        let i0 = fi.floor() as usize;
        let i1 = (i0 + 1).min(n - 1);
        let t = fi - i0 as f32;
        let e = elev[i0] * (1.0 - t) + elev[i1] * t;
        ws.profile_hover = Some(fx); // the map draws a matching point on the river highlight
        let sx = x0 + fx * xw;
        let sy = y_of(e);
        let length_km = fx * wc.length_km;
        // vertical cursor bar + point on the profile.
        p.line_segment(
            [egui::pos2(sx, rect.top()), egui::pos2(sx, rect.bottom())],
            egui::Stroke::new(1.0, C::from_rgba_unmultiplied(0xC9, 0x85, 0x3F, 120)),
        );
        p.circle_filled(egui::pos2(sx, sy), 3.5, COPPER_BRIGHT);
        // readout chip (top, following the cursor side).
        let txt = format!("{length_km:.0} km · {e:.0} m");
        let galley = p.layout_no_wrap(txt, egui::FontId::monospace(10.0), TEXT);
        let pad = egui::vec2(10.0, 6.0);
        let left = (sx + 8.0).min(rect.right() - galley.size().x - pad.x);
        let bgr = egui::Rect::from_min_size(
            egui::pos2(left.max(rect.left() + 2.0), rect.top() + 3.0),
            galley.size() + pad,
        );
        p.rect_filled(bgr, 4.0, C::from_rgba_unmultiplied(18, 18, 18, 225));
        p.galley(bgr.min + pad * 0.5, galley, TEXT);
    }
    nav
}

/// TASK 4 — lake sheet (Finding 28): the figures `C1Lake` carries + a couple computed
/// UI-side from `lake_map` (shore length, inlet count), and the endorheic CONSEQUENCE spelt
/// out (the content Living Landz should consume). Avg depth / inflow / evaporation are NOT
/// exported (the water balance discards them) — stated as such, not faked.
fn lake_sheet_panel(
    ui: &mut egui::Ui,
    hd: &HdResult,
    lake_idx: usize,
    domain_km: f32,
    wcs: Option<&[Watercourse]>,
) -> Option<NavAction> {
    let mut nav = None;
    let lk = &hd.drainage.lakes[lake_idx];
    let (w, h) = (hd.width, hd.height);
    let km_per_cell = if w > 0 { domain_km / w as f32 } else { 0.0 };
    // Shore length: count lake-boundary edges (a lake cell with a non-lake 4-neighbour).
    let id = lk.base.id;
    let (mut shore_cells, mut inlets) = (0usize, 0usize);
    for y in 0..h {
        for x in 0..w {
            if hd.drainage.lake_map[y * w + x] == id {
                for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx < 0
                        || ny < 0
                        || nx as usize >= w
                        || ny as usize >= h
                        || hd.drainage.lake_map[ny as usize * w + nx as usize] != id
                    {
                        shore_cells += 1;
                    }
                }
            }
        }
    }
    // Inlets: river reaches whose mouth is WITHIN 2 CELLS of this lake (Finding 34 — tolerate the
    // small gap left when the footprint grew after the tracks were traced; a track ending a cell
    // short still counts as a tributary).
    for s in &hd.drainage.rivers.segments {
        let &(mx, my) = s.points.last().unwrap();
        let mut adj = false;
        for dy in -2i32..=2 {
            for dx in -2i32..=2 {
                let (nx, ny) = (mx as i32 + dx, my as i32 + dy);
                if nx >= 0
                    && ny >= 0
                    && (nx as usize) < w
                    && (ny as usize) < h
                    && hd.drainage.lake_map[ny as usize * w + nx as usize] == id
                {
                    adj = true;
                }
            }
        }
        if adj {
            inlets += 1;
        }
    }
    let shore_km = shore_cells as f32 * km_per_cell;
    let endo = lk.lake_type == LakeType::Endorheic;
    let kv = |ui: &mut egui::Ui, k: &str, v: String| {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(k).color(DIM).size(11.0));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(egui::RichText::new(v).color(COPPER_BRIGHT).monospace().size(11.5));
            });
        });
    };
    ui.label(egui::RichText::new(format!("Lac #{id}")).color(TEXT_BRIGHT).strong().size(12.5));
    ui.add_space(2.0);
    kv(ui, "Type", if endo { "endoréique".into() } else { "exoréique".into() });
    kv(ui, "Superficie", format!("{:.1} km²", lk.area_km2));
    kv(ui, "Rive (périmètre)", format!("{shore_km:.0} km"));
    kv(ui, "Niveau", format!("{:.0} m", lk.level_m));
    kv(ui, "Profondeur max", format!("{:.0} m", lk.depth_m));
    kv(ui, "Affluents", format!("{inlets}"));
    kv(ui, "Exutoire", if endo { "aucun (fermé)".into() } else { "oui → aval".into() });
    ui.add_space(6.0);
    let (msg, col) = if endo {
        (
            "Bassin FERMÉ : l'eau s'évapore et CONCENTRE le sel → sel exploitable, pas de poisson, eau non potable, pas d'agriculture riveraine.",
            C::from_rgb(0x3a, 0xb0, 0xa0),
        )
    } else {
        (
            "Lac de PASSAGE : eau douce, exutoire vers l'aval — poisson, eau potable, berges cultivables.",
            C::from_rgb(0x5a, 0x9a, 0xd0),
        )
    };
    egui::Frame::default()
        .fill(C::from_rgb(0x17, 0x1c, 0x1b))
        .inner_margin(8)
        .corner_radius(5)
        .show(ui, |ui| {
            ui.label(egui::RichText::new(msg).color(col).size(10.5));
        });
    ui.add_space(6.0);
    // Clickable hydrological chain (Finding 29): inlets (upstream) + the outlet (downstream).
    // If the lake is EXORHEIC but no outlet reach is found, that inconsistency (H2) shows here.
    if let Some(wcs) = wcs {
        let adj = |mx: u32, my: u32| -> bool {
            (-2i32..=2).any(|dy| {
                (-2i32..=2).any(|dx| {
                    let (nx, ny) = (mx as i32 + dx, my as i32 + dy);
                    nx >= 0
                        && ny >= 0
                        && (nx as usize) < w
                        && (ny as usize) < h
                        && hd.drainage.lake_map[ny as usize * w + nx as usize] == id
                })
            })
        };
        ui.label(egui::RichText::new("CHAÎNE HYDROLOGIQUE").color(TEXT_BRIGHT).strong().size(10.5));
        // Inlets: watercourses whose MOUTH sits on this lake.
        let inflow_wcs: Vec<usize> = wcs
            .iter()
            .enumerate()
            .filter(|(_, wc)| adj(wc.mouth_xy.0, wc.mouth_xy.1))
            .map(|(i, _)| i)
            .collect();
        ui.horizontal_wrapped(|ui| {
            ui.label(
                egui::RichText::new(format!("Affluents ({}) :", inflow_wcs.len()))
                    .color(DIM)
                    .size(10.0),
            );
            for &wi in inflow_wcs.iter().take(12) {
                if ui
                    .add(
                        egui::Button::new(egui::RichText::new(format!("#{}", wi + 1)).size(10.0))
                            .frame(true),
                    )
                    .clicked()
                {
                    nav = Some(NavAction::River(wi));
                }
            }
        });
        // Outlet: the watercourse that DRAINS this lake — its SOURCE borders the lake, its MOUTH does
        // NOT (a real outlet leaves; it does not also enter), and among those it carries the most
        // discharge (the main stem, not a brook that merely rises next to the shore). An ENDORHEIC
        // lake is CLOSED → no outlet (Finding 37b: otherwise a phantom carve reach touching the basin
        // on both sides showed as an "outlet" #154/#926 that was ALSO an inlet).
        let outlet = if endo {
            None
        } else {
            wcs.iter()
                .enumerate()
                .filter(|(_, wc)| {
                    // The outlet REACH leaves the lake — a SEGMENT (anywhere in the watercourse, not
                    // just its far-upstream trunk source) whose FIRST point borders the shore — while
                    // the watercourse's MOUTH is elsewhere (it drains AWAY, it does not also enter).
                    // trunk.first alone missed it: a normal lake's outlet reach is a MIDDLE segment of
                    // a watercourse whose source is a tributary far upstream (#397 on #4).
                    let leaves = wc.segments.iter().any(|&s| {
                        hd.drainage.rivers.segments[s]
                            .points
                            .first()
                            .map(|&(sx, sy)| adj(sx, sy))
                            .unwrap_or(false)
                    });
                    leaves && !adj(wc.mouth_xy.0, wc.mouth_xy.1)
                })
                .max_by(|a, b| {
                    a.1.discharge_m3s
                        .partial_cmp(&b.1.discharge_m3s)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
        };
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new("Exutoire :").color(DIM).size(10.0));
            match outlet {
                Some(wi) => {
                    if ui
                        .add(
                            egui::Button::new(
                                egui::RichText::new(format!("→ rivière #{} ⮕", wi + 1)).size(10.0),
                            )
                            .frame(true),
                        )
                        .clicked()
                    {
                        nav = Some(NavAction::River(wi));
                    }
                }
                None if !endo => {
                    ui.label(
                        egui::RichText::new(
                            "⚠ AUCUN exutoire tracé — lac exoréique sans sortie (incohérence H2)",
                        )
                        .color(WARN_ORANGE)
                        .size(10.0),
                    );
                }
                None => {
                    ui.label(
                        egui::RichText::new("aucun (bassin fermé, cohérent)").color(DIM).size(10.0),
                    );
                }
            }
        });
    }
    ui.add_space(4.0);
    ui.label(egui::RichText::new("Profondeur moyenne · apport · évaporation : non exportés (bilan hydrique non surfacé — à ajouter côté données).").color(DIM).size(9.5));
    nav
}

/// Draw the river network onto `rgba` (row-major, w×h): each segment as a polyline,
/// stroke thickness by Strahler order, colour by navigability — and any ORPHAN reach
/// (no downstream link and NOT ending at a sink: sea or lake) in RED, so a broken /
/// dangling segment is immediately visible for validation. Lakes are tinted first so
/// rivers read against them. Pure raster pass — no effect on the exported data.
fn draw_river_overlay(rgba: &mut [u8], hd: &HdResult) {
    let (w, h) = (hd.width, hd.height);
    let d = &hd.drainage;
    // TASK 5 — draw EVERY water body the export carries (detect_lakes AND below-sea basins,
    // both in `lake_map`) so a channel never appears to stop in the void; distinguish by
    // lake_type: exorheic (through-flow) deep blue, endorheic (closed, saline) teal.
    let endorheic: std::collections::HashSet<u32> = d
        .lakes
        .iter()
        .filter(|lk| lk.lake_type == LakeType::Endorheic)
        .map(|lk| lk.base.id)
        .collect();
    for k in 0..w * h {
        let id = d.lake_map[k];
        if id != 0 {
            let c = if endorheic.contains(&id) { [30, 150, 140] } else { [30, 90, 180] };
            rgba[k * 4] = c[0];
            rgba[k * 4 + 1] = c[1];
            rgba[k * 4 + 2] = c[2];
        }
    }
    // Sea membership at the unified C1 sea level (eroded is normalised, sea = 0.5).
    let at_sink = |x: u32, y: u32| -> bool {
        let (x, y) = (x as i32, y as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx >= 0 && ny >= 0 && nx < w as i32 && ny < h as i32 {
                    let k = ny as usize * w + nx as usize;
                    if hd.eroded.data[k] <= 0.5 || d.lake_map[k] != 0 {
                        return true;
                    }
                }
            }
        }
        false
    };
    let put = |rgba: &mut [u8], x: i32, y: i32, c: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as usize) < w && (y as usize) < h {
            let k = (y as usize * w + x as usize) * 4;
            rgba[k] = c[0];
            rgba[k + 1] = c[1];
            rgba[k + 2] = c[2];
        }
    };
    for (i, seg) in d.rivers.segments.iter().enumerate() {
        let &(lx, ly) = seg.points.last().unwrap();
        let orphan = seg.downstream.is_none() && !at_sink(lx, ly);
        let nav = d.segment_navigability.get(i).copied().unwrap_or(Navigability::NonNavigable);
        let col = if orphan {
            [230, 30, 30]
        } else {
            match nav {
                Navigability::Ship => [20, 70, 200],
                Navigability::Barge => [40, 110, 230],
                Navigability::SmallBoat => [90, 160, 240],
                Navigability::NonNavigable => [120, 150, 190],
            }
        };
        // Thickness (half-width in cells) grows with Strahler order; orphans get a minimum of 1 so a
        // single dangling cell is still visible. Finding 39 — a below-sea SPILLWAY (`max_flow == 0`)
        // is a lake's OUTLET: if the lake is drawn, its outlet must be too, even when it is a 2–3 cell
        // stub to the adjacent sea. Give it the same minimum so it never vanishes.
        let is_spillway = seg.max_flow == 0.0;
        let r = ((seg.strahler_order as i32 - 2).max(0))
            .min(2)
            .max(orphan as i32)
            .max(is_spillway as i32);
        for &(px, py) in &seg.points {
            for dy in -r..=r {
                for dx in -r..=r {
                    put(rgba, px as i32 + dx, py as i32 + dy, col);
                }
            }
        }
    }
}

fn relief_color(norm: f32) -> [u8; 3] {
    let [r, g, b, _] = hypsometric_bipolar(norm, 0.5);
    [r, g, b]
}
fn precip_color(mm: f32) -> [u8; 3] {
    if mm < 250.0 {
        [225, 200, 140]
    } else if mm < 500.0 {
        [200, 195, 110]
    } else if mm < 800.0 {
        [150, 180, 90]
    } else if mm < 1500.0 {
        [80, 150, 200]
    } else {
        [30, 90, 200]
    }
}
fn temp_color(t: f32) -> [u8; 3] {
    if t < -5.0 {
        [225, 235, 248]
    } else if t < 5.0 {
        [90, 140, 205]
    } else if t < 20.0 {
        [110, 190, 110]
    } else {
        [225, 120, 70]
    }
}
fn drainage_color(hd: &HdResult, river_map: &RiverCellMap, k: usize) -> [u8; 3] {
    let w = hd.width;
    let (x, y) = (k % w, k / w);
    let [br, bg, bb] = relief_color(hd.eroded.data[k]);
    let dim = |v: u8| (v as f32 * 0.55) as u8;
    let mut col = [dim(br), dim(bg), dim(bb)];
    if hd.drainage.lake_map[k] != 0 {
        col = [30, 90, 180];
    }
    if let Some(info) = river_map.at(x, y) {
        col = match info.navigability {
            Navigability::Ship => [20, 70, 200],
            Navigability::Barge => [40, 110, 230],
            Navigability::SmallBoat => [90, 160, 240],
            Navigability::NonNavigable => [90, 120, 160],
        };
    }
    col
}

#[cfg(test)]
mod spillway_typing_bench {
    //! STEP 3b — the PRACTICAL SYMPTOM of the spillway mischaracterisation, measured on the
    //! production path. The author's report was "the first real river is #7": the microscope
    //! Rivières list, sorted by discharge, opened on six SPILLWAYS — an order-1 segment
    //! draining 88 468 km² outranks every genuine river because it carries the whole
    //! catchment of a closed below-sea basin.
    //!
    //! This bench lives HERE, next to `aggregate_watercourses`, so it exercises the SAME
    //! function the microscope calls (method rule 3: reproduce the production chain, not a
    //! reconstruction of it). It goes through `run_hd` end to end, so the segments, their
    //! kinds and the sort are all the shipped ones.
    //!
    //! Run: cargo test -p ymir-viz --release spillway_sort -- --ignored --nocapture

    use super::*;
    use crate::bridge::c1::C1RunSpec;
    use crate::bridge::c1::events::C1Event;
    use crate::bridge::c1::hd::{HdParams, run_hd};
    use crossbeam_channel::bounded;
    use std::sync::atomic::AtomicBool;
    use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
    use ymir_core::tectonics_c1::closures::infiltration::InfiltrationConfig;
    use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
    use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;

    const PSEED: u64 = 10_481_999_410_520_546_993;
    const DOMAIN_KM: f32 = 400.0;

    fn head_of_list(target: usize, lat: f32, label: &str) {
        let spec = C1RunSpec { seed: PSEED, ..C1RunSpec::default() };
        let params = HdParams {
            target_size: target,
            latitude_deg: lat,
            latitude_span_deg: Some(10.0),
            domain_km: DOMAIN_KM,
            manual_offset: None,
            // ADR Finding 83 — the bound SHIPS; these bench literals predate it.
            base_level_m: None,
            base_level_off: false,
            // relief-v3 = the shipped triple (stream_power + closures + mfd, cross_rill off).
            stream_power: true,
            closures: true,
            cross_rill: false,
            slope_floor_s_eq: None, // ADR Finding 109
            valley_construction: None, // ADR Finding 121
            cross_rill_d: 0.40,
            mfd: true,
            mfd_p: 2.0,
            geo_scale_ratio: 1.0,
            export_dir: None,
            volcanism: Some(VolcanismConfig {
                enabled: true,
                domain_km: DOMAIN_KM,
                ..Default::default()
            }),
            lithology: Some(LithologyConfig { enabled: true, ..Default::default() }),
            fracture: Some(FractureConfig {
                enabled: true,
                amplitude: 6.0,
                decay_km: 25.0,
                domain_km: DOMAIN_KM,
                ..Default::default()
            }),
            infiltration: Some(InfiltrationConfig { enabled: true, ..Default::default() }),
            emit_tectonic_labels: false,
            geology_rules: None,
        };
        let (tx, rx) = bounded(256);
        let cancel = Arc::new(AtomicBool::new(false));
        std::thread::spawn(move || run_hd(&spec, &params, &tx, &cancel));
        let mut result = None;
        while let Ok(e) = rx.recv() {
            match e {
                C1Event::HdCompleted { result: r, .. } => {
                    result = Some(r);
                    break;
                }
                C1Event::HdFailed { .. } => panic!("HD run failed"),
                _ => {}
            }
        }
        let hd = result.expect("HdCompleted");
        let raw_sw =
            hd.drainage.segment_kind.iter().filter(|k| **k == SegmentKind::Spillway).count();
        eprintln!(
            "RAW segments: {} | segment_kind len {} | tagged Spillway: {}",
            hd.drainage.rivers.segments.len(),
            hd.drainage.segment_kind.len(),
            raw_sw
        );
        let wcs = aggregate_watercourses(&hd, DOMAIN_KM, 1.0);
        let n_sw = wcs.iter().filter(|w| w.kind == SegmentKind::Spillway).count();
        let first_sw = wcs.iter().position(|w| w.kind == SegmentKind::Spillway);
        let first_river = wcs.iter().position(|w| w.kind == SegmentKind::Watercourse);

        eprintln!("\n=========  DISCHARGE SORT HEAD — {label} {lat}° @ {target}²  =========");
        eprintln!(
            "entries: {} | watercourses: {} | spillways: {}",
            wcs.len(),
            wcs.len() - n_sw,
            n_sw
        );
        eprintln!(
            "first WATERCOURSE at rank {:?} | first SPILLWAY at rank {:?}",
            first_river.map(|i| i + 1),
            first_sw.map(|i| i + 1)
        );
        eprintln!("\ntop 10:");
        for (i, w) in wcs.iter().take(10).enumerate() {
            eprintln!(
                "  #{:<3} {:>10} · Q {:>9.0} m³/s · A {:>9.0} km² · S{} · {} trib · {}",
                i + 1,
                if w.kind == SegmentKind::Spillway { "SPILLWAY" } else { "river" },
                w.discharge_m3s,
                w.catchment_km2,
                w.order,
                w.tributaries,
                sink_label(w.sink).0
            );
        }
        if n_sw > 0 {
            eprintln!("\nspillway block (first 8 of {n_sw}):");
            for (i, w) in
                wcs.iter().enumerate().filter(|(_, w)| w.kind == SegmentKind::Spillway).take(8)
            {
                eprintln!(
                    "  #{:<4} Q {:>9.0} m³/s · A {:>9.0} km² · source_lake {:?}",
                    i + 1,
                    w.discharge_m3s,
                    w.catchment_km2,
                    w.source_lake_id
                );
            }
            let dangling = wcs
                .iter()
                .filter(|w| w.kind == SegmentKind::Spillway && w.source_lake_id.is_none())
                .count();
            eprintln!(
                "\nspillways with NO inventoried source lake: {dangling} / {n_sw} (null by \
                 convention — a real basin below the lake-inventory floor)"
            );
        }
        // THE assertion: a spillway can no longer head the list while any river exists.
        if first_river.is_some() {
            assert_eq!(first_river, Some(0), "the head of the discharge sort must be a river");
        }
    }

    #[test]
    #[ignore]
    fn spillway_sort_arid_hot_2048() {
        head_of_list(2048, 25.0, "arid-hot");
    }

    #[test]
    #[ignore]
    fn spillway_sort_arid_hot_8192() {
        head_of_list(8192, 25.0, "arid-hot");
    }
}

#[cfg(test)]
mod network_fragmentation_bench {
    //! DIAGNOSTIC ONLY — the network FRAGMENTS as resolution increases (largest assembled
    //! catchment 1087 km² at 2048² against 110 km² at 8192², with FEWER list entries on the
    //! finer grid: 792 against 1359). The opposite of the expected behaviour, and it hits the
    //! MEASURING INSTRUMENT: the microscope assembles watercourses by following `downstream`
    //! to a terminal, so if the assembly breaks at high resolution, every measurement that
    //! reads that list describes truncated fragments.
    //!
    //! Four questions, four measurements, no fix:
    //!   1. WHERE the assembly breaks — terminals that are neither sea, lake nor sub-sea sink,
    //!      and whether the raster still has a receiver there.
    //!   2. ASSEMBLY or NETWORK — the mouth's area read from the flow-accumulation RASTER
    //!      (independent of segmentation) beside the assembled figure.
    //!   3. WHY FEWER ENTRIES — the km2-to-cell threshold conversion and the network's extent.
    //!   4. Is `clip_rivers_to_lakes` involved — how many terminations it creates at each grid.
    //!
    //! Run: cargo test -p ymir-viz --release network_fragmentation -- --ignored --nocapture

    use super::*;
    use crate::bridge::c1::C1RunSpec;
    use crate::bridge::c1::events::C1Event;
    use crate::bridge::c1::hd::{HdParams, run_hd};
    use crossbeam_channel::bounded;
    use std::sync::atomic::AtomicBool;
    use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
    use ymir_core::tectonics_c1::closures::infiltration::InfiltrationConfig;
    use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
    use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;
    use ymir_core::tectonics_c1::drainage::C1DrainageConfig;

    const PSEED: u64 = 10_481_999_410_520_546_993;
    const DOMAIN_KM: f32 = 400.0;

    fn production_run(target: usize, lat: f32) -> Arc<HdResult> {
        let spec = C1RunSpec { seed: PSEED, ..C1RunSpec::default() };
        let params = HdParams {
            target_size: target,
            latitude_deg: lat,
            latitude_span_deg: Some(10.0),
            domain_km: DOMAIN_KM,
            manual_offset: None,
            // ADR Finding 83 — the bound SHIPS; these bench literals predate it.
            base_level_m: None,
            base_level_off: false,
            stream_power: true,
            closures: true,
            cross_rill: false,
            slope_floor_s_eq: None, // ADR Finding 109
            valley_construction: None, // ADR Finding 121
            cross_rill_d: 0.40,
            mfd: true,
            mfd_p: 2.0,
            geo_scale_ratio: 1.0,
            export_dir: None,
            volcanism: Some(VolcanismConfig {
                enabled: true,
                domain_km: DOMAIN_KM,
                ..Default::default()
            }),
            lithology: Some(LithologyConfig { enabled: true, ..Default::default() }),
            fracture: Some(FractureConfig {
                enabled: true,
                amplitude: 6.0,
                decay_km: 25.0,
                domain_km: DOMAIN_KM,
                ..Default::default()
            }),
            infiltration: Some(InfiltrationConfig { enabled: true, ..Default::default() }),
            emit_tectonic_labels: false,
            geology_rules: None,
        };
        let (tx, rx) = bounded(256);
        let cancel = Arc::new(AtomicBool::new(false));
        std::thread::spawn(move || run_hd(&spec, &params, &tx, &cancel));
        loop {
            match rx.recv().expect("worker channel") {
                C1Event::HdCompleted { result, .. } => return result,
                C1Event::HdFailed { error } => panic!("HD run failed: {error}"),
                _ => {}
            }
        }
    }

    /// The numbers question 2 needs across resolutions: normalised mouth position + both
    /// readings of its catchment (raster-geometric and assembled-effective).
    struct MouthProbe {
        uv: (f64, f64),
        geom_km2: f32,
        effective_km2: f32,
    }

    fn report(target: usize, lat: f32) -> MouthProbe {
        let hd = production_run(target, lat);
        let d = &hd.drainage;
        let (w, h) = (hd.width, hd.height);
        let n = w * h;
        let cell_km2 = (DOMAIN_KM / w as f32) * (DOMAIN_KM / h as f32);
        let wc = ymir_core::lakes::connectivity::water_class(&hd.eroded, 0.5);
        let segs = &d.rivers.segments;
        let wcs = aggregate_watercourses(&hd, DOMAIN_KM, 1.0);

        eprintln!("\n============  NETWORK FRAGMENTATION — {lat} deg @ {target}^2  ============");
        eprintln!(
            "cell = {:.4} km2 ({:.0} m); segments {} | watercourse entries {} | river cells {}",
            cell_km2,
            DOMAIN_KM * 1000.0 / w as f32,
            segs.len(),
            wcs.len(),
            segs.iter().map(|s| s.points.len()).sum::<usize>()
        );

        // 1. TERMINAL AUDIT — where does the chain stop, and is that legitimate?
        let mut term_sea = 0usize;
        let mut term_lake = 0usize;
        let mut term_subsea = 0usize;
        let mut term_none = 0usize;
        let mut term_none_max_km2 = 0.0f32;
        let mut term_none_with_receiver = 0usize;
        let mut terminals = 0usize;
        for (i, s) in segs.iter().enumerate() {
            if s.downstream.is_some() {
                continue;
            }
            terminals += 1;
            let (mx, my) = *s.points.last().unwrap();
            let (sink, _) = classify_sink(&hd, &wc, mx, my);
            match sink {
                Sink::Sea => term_sea += 1,
                Sink::ExoLake | Sink::EndoLake => term_lake += 1,
                Sink::SubSeaSink => term_subsea += 1,
                Sink::Unknown => {
                    term_none += 1;
                    let a = d.segment_drainage_km2.get(i).copied().unwrap_or(0.0);
                    term_none_max_km2 = term_none_max_km2.max(a);
                    // Does the RASTER still route out of this cell? A `None` here while the
                    // D8 field has a receiver is an assembly break, not a real terminus.
                    let k = my as usize * w + mx as usize;
                    if d.flow.direction[k] != DIR_NONE {
                        term_none_with_receiver += 1;
                    }
                }
            }
        }
        eprintln!(
            "\n1. TERMINALS ({terminals} of {} segments)\n   to sea {term_sea} | to lake \
             {term_lake} | to sub-sea sink {term_subsea} | NONE OF THESE {term_none}\n   of \
             those: {} still have a D8 receiver in the raster (assembly break, not a \
             terminus); largest {term_none_max_km2:.0} km2",
            segs.len(),
            term_none_with_receiver
        );
        // Isolated fragments: no upstream AND no downstream — a reach chained to nothing.
        let isolated =
            segs.iter().filter(|s| s.upstream.is_empty() && s.downstream.is_none()).count();
        eprintln!(
            "   isolated fragments (no upstream AND no downstream): {isolated} ({:.1} % of segments)",
            100.0 * isolated as f32 / segs.len().max(1) as f32
        );

        // Follow the chain down from the largest reach and say where it stops.
        let big = (0..segs.len())
            .max_by(|&a, &b| {
                d.segment_drainage_km2[a]
                    .partial_cmp(&d.segment_drainage_km2[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(0);
        let mut hops = 0usize;
        let mut cur = big;
        while let Some(nx) = segs[cur].downstream {
            if nx >= segs.len() || hops > segs.len() {
                break;
            }
            cur = nx;
            hops += 1;
        }
        let (ex, ey) = *segs[cur].points.last().unwrap();
        let (esink, _) = classify_sink(&hd, &wc, ex, ey);
        eprintln!(
            "   largest reach (#{big}, {:.0} km2 eff) -> {hops} hop(s) -> terminal #{cur} at \
             ({ex},{ey}) {} | D8 receiver {} | raster area there {:.0} km2",
            d.segment_drainage_km2[big],
            sink_label(esink).0,
            if d.flow.direction[ey as usize * w + ex as usize] == DIR_NONE {
                "NONE"
            } else {
                "yes"
            },
            d.flow.accumulation.data[ey as usize * w + ex as usize] * cell_km2
        );

        // 2. RASTER vs ASSEMBLY at the mouth. The raster accumulation is a GEOMETRIC cell
        //    count, untouched by segmentation; `segment_drainage_km2` is an EFFECTIVE area
        //    derived from the runoff accumulation. Reading both separates the two causes.
        let mut best_k = 0usize;
        let mut best_acc = -1.0f32;
        for k in 0..n {
            let a = d.flow.accumulation.data[k];
            if a > best_acc && wc[k] == ymir_core::lakes::connectivity::WATER_CLASS_LAND {
                best_acc = a;
                best_k = k;
            }
        }
        let (bx, by) = ((best_k % w) as u32, (best_k / w) as u32);
        // The assembled figure at that same cell: the segment whose points contain it.
        let at_cell = segs
            .iter()
            .enumerate()
            .find(|(_, s)| s.points.iter().any(|&(x, y)| x == bx && y == by))
            .map(|(i, _)| d.segment_drainage_km2[i]);
        eprintln!(
            "\n2. MOUTH — max raster accumulation on land at ({bx},{by}) = {:.0} km2 \
             GEOMETRIC\n   assembled effective area of the segment covering that cell: {}\n   \
             list head: {:.0} km2 effective, {:.1} m3/s",
            best_acc * cell_km2,
            at_cell.map(|a| format!("{a:.0} km2")).unwrap_or_else(|| "NOT ON THE NETWORK".into()),
            wcs.first().map(|x| x.catchment_km2).unwrap_or(0.0),
            wcs.first().map(|x| x.discharge_m3s).unwrap_or(0.0)
        );
        let max_seg_eff = d.segment_drainage_km2.iter().copied().fold(0.0f32, f32::max);
        let max_geom_km2 =
            d.flow.accumulation.data.iter().copied().fold(0.0f32, f32::max) * cell_km2;
        eprintln!(
            "   max over ALL segments: {max_seg_eff:.0} km2 effective | max raster accumulation \
             anywhere: {max_geom_km2:.0} km2 geometric\n   that max-accumulation cell: lake_map \
             = {} | water_class = {}",
            d.lake_map[best_k], wc[best_k]
        );

        // 2b. THE decisive pair: at every SEA-terminating mouth, the GEOMETRIC area (raster) and
        //     the EFFECTIVE area (runoff-derived) side by side. If geometric holds across
        //     resolutions while effective collapses, the topology is intact and the collapse is
        //     the runoff/climate integration — not the network, not the segmentation.
        let mut sea_mouths: Vec<(f32, f32)> = Vec::new(); // (geometric, effective)
        for (i, s) in segs.iter().enumerate() {
            if s.downstream.is_some() {
                continue;
            }
            let (mx, my) = *s.points.last().unwrap();
            if classify_sink(&hd, &wc, mx, my).0 != Sink::Sea {
                continue;
            }
            let k = my as usize * w + mx as usize;
            sea_mouths.push((
                d.flow.accumulation.data[k] * cell_km2,
                d.segment_drainage_km2.get(i).copied().unwrap_or(0.0),
            ));
        }
        sea_mouths.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!("\n2b. SEA MOUTHS ({} of them), by GEOMETRIC area — top 5:", sea_mouths.len());
        for (g, e) in sea_mouths.iter().take(5) {
            eprintln!(
                "   geometric {g:>8.0} km2 | effective {e:>8.0} km2 | effective/geometric {:.3}",
                e / g.max(1e-6)
            );
        }
        let tot_g: f32 = sea_mouths.iter().map(|x| x.0).sum();
        let tot_e: f32 = sea_mouths.iter().map(|x| x.1).sum();
        eprintln!(
            "   summed over all sea mouths: geometric {tot_g:.0} km2 | effective {tot_e:.0} km2 \
             | ratio {:.3}",
            tot_e / tot_g.max(1e-6)
        );
        let dry_big = sea_mouths.iter().filter(|(g, e)| *g >= 100.0 && *e <= 0.0).count();
        eprintln!(
            "   sea mouths with geometric >= 100 km2 but effective == 0: {dry_big} (a mouth with \
             a real catchment and NO discharge)"
        );

        // 2c. Is that zero the ENDORHEIC MASK or genuine aridity? Re-accumulate runoff with and
        //     without the endorheic sink mask (infiltration omitted in BOTH — it is a
        //     discriminant, not a production figure). If unmasked is large where masked is 0,
        //     a lake declared endorheic is killing a route that geometrically reaches the sea.
        let dclim = ymir_core::tectonics_c1::drainage::DrainageClimate {
            precip_internal: &hd.precipitation,
            temperature: &hd.temperature,
        };
        let mut endo_mask = vec![false; n];
        for lk in d.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic) {
            for k in 0..n {
                if d.lake_map[k] == lk.base.id {
                    endo_mask[k] = true;
                }
            }
        }
        let unmasked = ymir_core::tectonics_c1::drainage::runoff_accumulation(
            &hd.eroded, &d.flow, &dclim, cell_km2, None, None, w, h,
        );
        let masked = ymir_core::tectonics_c1::drainage::runoff_accumulation(
            &hd.eroded,
            &d.flow,
            &dclim,
            cell_km2,
            Some(&endo_mask),
            None,
            w,
            h,
        );
        let mut killed = 0usize;
        let mut killed_max_g = 0.0f32;
        for s in segs.iter().filter(|s| s.downstream.is_none()) {
            let (mx, my) = *s.points.last().unwrap();
            if classify_sink(&hd, &wc, mx, my).0 != Sink::Sea {
                continue;
            }
            let k = my as usize * w + mx as usize;
            let g = d.flow.accumulation.data[k] * cell_km2;
            if g >= 100.0 && masked[k] <= 0.0 && unmasked[k] > 0.0 {
                killed += 1;
                killed_max_g = killed_max_g.max(g);
            }
        }
        eprintln!(
            "   of those, killed by the ENDORHEIC MASK (unmasked runoff > 0, masked == 0): \
             {killed}, largest {killed_max_g:.0} km2 geometric\n   endorheic lake cells: {} | \
             endorheic lakes: {}",
            endo_mask.iter().filter(|&&b| b).count(),
            d.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).count()
        );

        // 2d. DUPLICATE MOUTHS — the sea-mouth table shows the same geometric area twice over
        //     (1459/1459, 738/738). Count distinct terminal cells against terminals.
        let distinct_terminal_cells: std::collections::HashSet<(u32, u32)> = segs
            .iter()
            .filter(|s| s.downstream.is_none())
            .map(|s| *s.points.last().unwrap())
            .collect();
        eprintln!(
            "\n2d. DUPLICATE TERMINALS — {terminals} terminals for {} distinct last cells \
             ({} duplicated)",
            distinct_terminal_cells.len(),
            terminals - distinct_terminal_cells.len()
        );

        // 2e. `segment_drainage_km2` is runoff_accum / 300 mm — a runoff-EQUIVALENT area, so
        //     effective/geometric IS the catchment's mean net runoff over 300 mm/yr. If that
        //     mean moves with resolution, the "catchment" moves with it while the topology
        //     stands still. Read the climate terms directly, over land.
        let mut n_land = 0usize;
        let mut sum_p = 0.0f64;
        let mut sum_pe = 0.0f64;
        let mut sum_net = 0.0f64;
        let mut n_wet = 0usize;
        for k in 0..n {
            if hd.eroded.data[k] <= 0.5 {
                continue;
            }
            n_land += 1;
            let p = ymir_core::climate::precipitation::precip_mm_per_year(hd.precipitation.data[k]);
            let pe =
                ymir_core::tectonics_c1::drainage::potential_evaporation_mm(hd.temperature.data[k]);
            sum_p += p as f64;
            sum_pe += pe as f64;
            sum_net += (p - pe).max(0.0) as f64;
            if p > pe {
                n_wet += 1;
            }
        }
        let nl = n_land.max(1) as f64;
        eprintln!(
            "\n2e. CLIMATE TERMS over land ({n_land} cells)\n   mean precip {:.0} mm/yr | mean PE \
             {:.0} mm/yr | mean net runoff max(0,p-pe) {:.1} mm/yr\n   cells with p > pe: {n_wet} \
             ({:.2} %) | net runoff / 300 mm reference = {:.3}",
            sum_p / nl,
            sum_pe / nl,
            sum_net / nl,
            100.0 * n_wet as f64 / nl,
            sum_net / nl / 300.0
        );

        // 2f. The transport is resolution-invariant on FIXED geometry (see
        //     `tests/precip_resolution_invariance.rs`: 391 mm/yr at every grid from 512² to
        //     8192²), so a production difference in mean precipitation must come from the
        //     TERRAIN. Measure the terrain the climate reads: altitude, roughness, temperature.
        let mut sum_alt = 0.0f64;
        let mut sum_t = 0.0f64;
        let mut sum_slope = 0.0f64;
        let m_per_cell = DOMAIN_KM * 1000.0 / w as f32;
        let ss = ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams::default();
        let to_m = |nrm: f32| {
            ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres(nrm, &ss)
        };
        for k in 0..n {
            if hd.eroded.data[k] <= 0.5 {
                continue;
            }
            sum_alt += to_m(hd.eroded.data[k]) as f64;
            sum_t += hd.temperature.data[k] as f64;
            let (x, y) = (k % w, k / w);
            let (gx, gy) = (
                hd.eroded.get((x + 1) as i32, y as i32) - hd.eroded.get(x as i32 - 1, y as i32),
                hd.eroded.get(x as i32, (y + 1) as i32) - hd.eroded.get(x as i32, y as i32 - 1),
            );
            // Central differences → gradient magnitude in normalised units per 2 cells.
            let dz_m = to_m(0.5 + (gx * gx + gy * gy).sqrt() / 2.0) - to_m(0.5);
            sum_slope += (dz_m / m_per_cell).atan().to_degrees() as f64;
        }
        // HYPSOMETRY — does the whole distribution shift, or only the peaks? That decides
        // whether the cause is a global vertical scale or added fine relief.
        let mut alts: Vec<f32> =
            (0..n).filter(|&k| hd.eroded.data[k] > 0.5).map(|k| to_m(hd.eroded.data[k])).collect();
        alts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let q = |f: f32| alts[((alts.len() as f32 - 1.0) * f) as usize];
        eprintln!(
            "\n2f. THE TERRAIN THE CLIMATE READS (land)\n   mean altitude {:.0} m | mean \
             temperature {:.2} C | mean slope {:.2} deg | emerged {:.2} %\n   hypsometry p10 \
             {:.0} | p50 {:.0} | p90 {:.0} | p99 {:.0} | max {:.0} m",
            sum_alt / nl,
            sum_t / nl,
            sum_slope / nl,
            100.0 * nl / n as f64,
            q(0.10),
            q(0.50),
            q(0.90),
            q(0.99),
            alts[alts.len() - 1]
        );
        // Audit line — the same statistics in RAW NORMALISED units, so the metric figures can
        // be checked against the linear `c1_altitude_norm_to_metres` map without trusting it.
        let mut nrm: Vec<f32> =
            (0..n).filter(|&k| hd.eroded.data[k] > 0.5).map(|k| hd.eroded.data[k]).collect();
        nrm.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!(
            "   raw normalised (sea level 0.5): mean {:.5} | p50 {:.5} | p90 {:.5} | max {:.5}",
            nrm.iter().map(|&v| v as f64).sum::<f64>() / nrm.len() as f64,
            nrm[nrm.len() / 2],
            nrm[nrm.len() * 9 / 10],
            nrm[nrm.len() - 1]
        );

        // 3. THRESHOLDS — km2 converted to cells at THIS resolution, and how much of the
        //    raster clears them. A finer grid must have MORE cells above a physical threshold.
        let th = C1DrainageConfig::default().thresholds;
        let stream_cells = (th.stream_km2 / cell_km2).max(1.0);
        let head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
        let head_cells = (head_km2 / cell_km2).max(1.0);
        let above_stream = (0..n).filter(|&k| d.flow.accumulation.data[k] >= stream_cells).count();
        let above_head = (0..n).filter(|&k| d.flow.accumulation.data[k] >= head_cells).count();
        eprintln!(
            "\n3. THRESHOLDS (physical km2, converted per resolution)\n   stream {:.0} km2 = \
             {stream_cells:.1} cells -> {above_stream} cells clear it ({:.0} km2 of channel)\n   \
             head {head_km2:.2} km2 = {head_cells:.1} cells -> {above_head} cells clear it\n   \
             land cells: {} | river cells actually exported: {}",
            th.stream_km2,
            above_stream as f32 * cell_km2,
            (0..n).filter(|&k| wc[k] == ymir_core::lakes::connectivity::WATER_CLASS_LAND).count(),
            segs.iter().map(|s| s.points.len()).sum::<usize>()
        );

        // 4. CLIPPING — a run that ends at a lake shore gets `downstream = None`. Count the
        //    terminations the clip creates, and the lake footprint the network must cross.
        let lake_cells = d.lake_map.iter().filter(|&&id| id != 0).count();
        let starts_at_lake = segs
            .iter()
            .filter(|s| {
                let (sx, sy) = s.points[0];
                (-1i32..=1).any(|dy| {
                    (-1i32..=1).any(|dx| {
                        let (nx, ny) = (sx as i32 + dx, sy as i32 + dy);
                        nx >= 0
                            && ny >= 0
                            && (nx as usize) < w
                            && (ny as usize) < h
                            && d.lake_map[ny as usize * w + nx as usize] != 0
                    })
                })
            })
            .count();
        eprintln!(
            "\n4. CLIPPING\n   lake cells {lake_cells} ({:.0} km2, {:.2} % of the grid) | lakes \
             {}\n   terminations at a lake shore: {term_lake} | segments STARTING at a lake \
             shore: {starts_at_lake}",
            lake_cells as f32 * cell_km2,
            100.0 * lake_cells as f32 / n as f32,
            d.lakes.len()
        );
        MouthProbe {
            uv: (bx as f64 / w as f64, by as f64 / h as f64),
            geom_km2: best_acc * cell_km2,
            effective_km2: at_cell.unwrap_or(0.0),
        }
    }

    /// ADR Finding 45 — MICROSCOPE SEMANTICS, before and after, on ONE production run so the
    /// comparison is exact (both aggregations read the same `HdResult`).
    ///
    /// The defect: a reach ending at a lake shore gets `downstream = None`, so 87–92 % of
    /// terminals were lake INFLOWS and a list entry was a REACH BETWEEN TWO WATER BODIES, not
    /// a river SYSTEM. The head of the discharge sort therefore showed a fragment. Chaining an
    /// inflow reach to its EXORHEIC lake's outlet reach makes an entry a system again — not for
    /// endorheic lakes (the water dies there) nor for below-sea basins (typed `Spillway`).
    ///
    /// Run: cargo test -p ymir-viz --release microscope_semantics -- --ignored --nocapture
    fn semantics_before_after(target: usize, lat: f32) {
        let hd = production_run(target, lat);
        eprintln!("\n=====  MICROSCOPE SEMANTICS — {lat} deg @ {target}^2  =====");
        for (label, chain) in
            [("BEFORE (reach between water bodies)", false), ("AFTER (river system)", true)]
        {
            let wcs = aggregate_watercourses_opt(&hd, DOMAIN_KM, 1.0, chain);
            let rivers: Vec<&Watercourse> =
                wcs.iter().filter(|w| w.kind == SegmentKind::Watercourse).collect();
            let n_sw = wcs.len() - rivers.len();
            eprintln!(
                "\n{label}\n  entries {} (rivers {}, spillways {})",
                wcs.len(),
                rivers.len(),
                n_sw
            );
            eprintln!("  head of the discharge sort:");
            for (i, w) in rivers.iter().take(5).enumerate() {
                eprintln!(
                    "    #{:<3} Q {:>8.1} m3/s | A {:>8.0} km2 | S{} | {:>4} trib | {:>4} seg | \
                     {:>7.0} km | {}",
                    i + 1,
                    w.discharge_m3s,
                    w.catchment_km2,
                    w.order,
                    w.tributaries,
                    w.segments.len(),
                    w.length_km,
                    sink_label(w.sink).0
                );
            }
            let biggest = rivers.iter().map(|w| w.segments.len()).max().unwrap_or(0);
            let longest = rivers.iter().map(|w| w.length_km).fold(0.0f32, f32::max);
            eprintln!(
                "  largest system: {biggest} segments | longest trunk {longest:.0} km | \
                 mean segments/entry {:.1}",
                rivers.iter().map(|w| w.segments.len()).sum::<usize>() as f32
                    / rivers.len().max(1) as f32
            );
        }
    }

    /// ADR Finding 46 — the fragment area-inheritance fix, plus the RULE-7 CONTROL BLOCK
    /// (hypsometry in metres, closed-depression / below-sea count, network extent) because
    /// this touches a field every consumer reads, whatever the subject of the round.
    ///
    /// Run: cargo test -p ymir-viz --release fragment_areas -- --ignored --nocapture
    fn fragment_report(target: usize, lat: f32) {
        let hd = production_run(target, lat);
        let d = &hd.drainage;
        let (w, h) = (hd.width, hd.height);
        let n = w * h;
        let cell_km2 = (DOMAIN_KM / w as f32) * (DOMAIN_KM / h as f32);
        let segs = &d.rivers.segments;
        let wc = ymir_core::lakes::connectivity::water_class(&hd.eroded, 0.5);

        eprintln!("\n=====  FRAGMENT AREAS + RULE-7 CONTROL — {lat} deg @ {target}^2  =====");

        // ── the defect's own metric: fragments (no upstream, no downstream) and what they claim
        let frags: Vec<usize> = (0..segs.len())
            .filter(|&i| segs[i].upstream.is_empty() && segs[i].downstream.is_none())
            .collect();
        let frag_max = frags.iter().map(|&i| d.segment_drainage_km2[i]).fold(0.0f32, f32::max);
        let global_max = d.segment_drainage_km2.iter().copied().fold(0.0f32, f32::max);
        let frag_over_half: usize =
            frags.iter().filter(|&&i| d.segment_drainage_km2[i] > 0.5 * global_max).count();
        eprintln!(
            "\nFRAGMENTS (no upstream AND no downstream): {} of {} segments ({:.1} %)\n  \
             largest fragment area {frag_max:.0} km2 | largest area anywhere {global_max:.0} km2\n  \
             fragments claiming > half the largest area: {frag_over_half}  <-- the defect's \
             signature (an inherited trunk figure on a 1-cell reach)",
            frags.len(),
            segs.len(),
            100.0 * frags.len() as f32 / segs.len().max(1) as f32
        );

        // ── duplicated terminals
        let terminals: Vec<usize> =
            (0..segs.len()).filter(|&i| segs[i].downstream.is_none()).collect();
        let distinct: std::collections::HashSet<(u32, u32)> =
            terminals.iter().map(|&i| *segs[i].points.last().unwrap()).collect();
        // A duplicate is only a READING problem when the two rows claim the same area.
        let mut by_cell: std::collections::HashMap<(u32, u32), Vec<usize>> = Default::default();
        for &i in &terminals {
            by_cell.entry(*segs[i].points.last().unwrap()).or_default().push(i);
        }
        let same_area_pairs = by_cell
            .values()
            .filter(|v| v.len() > 1)
            .filter(|v| {
                let a = d.segment_drainage_km2[v[0]];
                v.iter().all(|&j| (d.segment_drainage_km2[j] - a).abs() < 0.5)
            })
            .count();
        eprintln!(
            "DUPLICATED TERMINALS: {} terminals for {} distinct last cells ({} duplicated), of \
             which {same_area_pairs} groups still claim the SAME area",
            terminals.len(),
            distinct.len(),
            terminals.len() - distinct.len()
        );

        // ── DUPLICATE TERMINAL ANATOMY. My first attribution (inherited area on a middle
        //    fragment) was refuted by the measurement: the head stump reaches the MOUTH, so its
        //    full area is CORRECT. The defect is that two segments end on the same mouth cell.
        //    Print the biggest such group so the cause is read, not guessed.
        if let Some((cell, group)) = by_cell.iter().filter(|(_, v)| v.len() > 1).max_by(|a, b| {
            d.segment_drainage_km2[a.1[0]]
                .partial_cmp(&d.segment_drainage_km2[b.1[0]])
                .unwrap_or(std::cmp::Ordering::Equal)
        }) {
            eprintln!(
                "
ANATOMY of the largest duplicated terminal, at cell {cell:?}:"
            );
            for &i in group {
                let sg = &segs[i];
                eprintln!(
                    "  seg #{i:<6} pts {:<5} S{} | A {:>8.0} km2 | Q {:>7.1} | up {:?} | down {:?}
                         first {:?} last {:?}",
                    sg.points.len(),
                    sg.strahler_order,
                    d.segment_drainage_km2[i],
                    d.segment_discharge_m3s[i],
                    sg.upstream.len(),
                    sg.downstream,
                    sg.points.first(),
                    sg.points.last()
                );
            }
            // Is one a suffix of the other (the clip emitting a run twice), or are they two
            // genuinely different paths converging on the same cell?
            if group.len() == 2 {
                let (x, y) = (&segs[group[0]].points, &segs[group[1]].points);
                let shared = x.iter().filter(|p| y.contains(p)).count();
                eprintln!(
                    "  shared points: {shared} of {}/{} -> {}",
                    x.len(),
                    y.len(),
                    if shared == x.len().min(y.len()) {
                        "ONE IS A SUBSET of the other (the same run emitted twice)"
                    } else if shared == 0 {
                        "DISJOINT paths converging on one cell (two tributaries, legitimate)"
                    } else {
                        "PARTIAL overlap"
                    }
                );
            }
        }

        // ── the head of the discharge sort (the author's test)
        let wcs = aggregate_watercourses(&hd, DOMAIN_KM, 1.0);
        let rivers: Vec<&Watercourse> =
            wcs.iter().filter(|x| x.kind == SegmentKind::Watercourse).collect();
        eprintln!(
            "\nMICROSCOPE — {} entries ({} rivers, {} spillways); head of the discharge sort:",
            wcs.len(),
            rivers.len(),
            wcs.len() - rivers.len()
        );
        for (i, x) in rivers.iter().take(5).enumerate() {
            eprintln!(
                "  #{:<3} Q {:>8.1} m3/s | A {:>8.0} km2 | S{} | {:>4} trib | {:>4} seg | \
                 {:>6.0} km | {}",
                i + 1,
                x.discharge_m3s,
                x.catchment_km2,
                x.order,
                x.tributaries,
                x.segments.len(),
                x.length_km,
                sink_label(x.sink).0
            );
        }

        // ── RULE-7 CONTROL BLOCK: the observables this change is NOT supposed to move.
        let ss = ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams::default();
        let to_m = |v: f32| {
            ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres(v, &ss)
        };
        let mut alts: Vec<f32> =
            (0..n).filter(|&k| hd.eroded.data[k] > 0.5).map(|k| to_m(hd.eroded.data[k])).collect();
        alts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let q = |f: f32| alts[((alts.len() as f32 - 1.0) * f) as usize];
        let below_sea_basins = wcs.iter().filter(|x| x.kind == SegmentKind::Spillway).count();
        let river_cells: usize = segs.iter().map(|s| s.points.len()).sum();
        eprintln!(
            "\nRULE-7 CONTROL BLOCK (must be unchanged — this touches per-segment arrays only)\n  \
             hypsometry (m): mean {:.0} | p10 {:.0} | p50 {:.0} | p90 {:.0} | max {:.0}\n  \
             below-sea basins (spillway entries): {below_sea_basins}\n  \
             network extent: {river_cells} river cells = {:.0} km of channel\n  \
             lake cells: {} ({:.0} km2) | inventoried lakes: {}",
            alts.iter().map(|&v| v as f64).sum::<f64>() / alts.len().max(1) as f64,
            q(0.10),
            q(0.50),
            q(0.90),
            alts[alts.len() - 1],
            river_cells as f32 * (DOMAIN_KM / w as f32),
            d.lake_map.iter().filter(|&&x| x != 0).count(),
            d.lake_map.iter().filter(|&&x| x != 0).count() as f32 * cell_km2,
            d.lakes.len()
        );
        let _ = wc;
    }

    /// ADR Finding 47 — the navigability question, in three columns so the choice is made on
    /// numbers: the thresholds the CODE ships, the ones the author PROPOSED (and which were
    /// never applied), and ones re-anchored on the measured GEOMETRIC distribution.
    ///
    /// Run: cargo test -p ymir-viz --release navigability_columns -- --ignored --nocapture
    fn navigability_report(target: usize, lat: f32) {
        let hd = production_run(target, lat);
        let d = &hd.drainage;
        let (w, h) = (hd.width, hd.height);
        let cell_km2 = (DOMAIN_KM / w as f32) * (DOMAIN_KM / h as f32);
        let segs = &d.rivers.segments;
        let wc = ymir_core::lakes::connectivity::water_class(&hd.eroded, 0.5);

        eprintln!("\n=====  NAVIGABILITY — {lat} deg @ {target}^2  =====");

        // The three quantities, per reach: geometric catchment, runoff-equivalent area, discharge.
        let geo: Vec<f32> =
            (0..segs.len()).map(|i| d.segment_catchment_cells[i] * cell_km2).collect();
        let runoff: Vec<f32> = d.segment_drainage_km2.clone();
        let disch: Vec<f32> = d.segment_discharge_m3s.clone();

        // At the MOUTHS (sea-terminating reaches) — the population the calibration was anchored on.
        let mouths: Vec<usize> = (0..segs.len())
            .filter(|&i| {
                segs[i].downstream.is_none() && {
                    let (mx, my) = *segs[i].points.last().unwrap();
                    classify_sink(&hd, &wc, mx, my).0 == Sink::Sea
                }
            })
            .collect();
        let pct = |v: &mut Vec<f32>, f: f32| {
            v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            v[((v.len() as f32 - 1.0) * f).max(0.0) as usize]
        };
        let mut mg: Vec<f32> = mouths.iter().map(|&i| geo[i]).collect();
        let mut mr: Vec<f32> = mouths.iter().map(|&i| runoff[i]).collect();
        let mut md: Vec<f32> = mouths.iter().map(|&i| disch[i]).collect();
        eprintln!(
            "\nAT THE {} SEA MOUTHS — the two 'areas' side by side\n  \
             GEOMETRIC km2  p50 {:>8.0} | p90 {:>8.0} | max {:>8.0}\n  \
             RUNOFF-EQ km2  p50 {:>8.0} | p90 {:>8.0} | max {:>8.0}   <- what navigability reads\n  \
             DISCHARGE m3/s p50 {:>8.2} | p90 {:>8.2} | max {:>8.2}",
            mouths.len(),
            pct(&mut mg, 0.50),
            pct(&mut mg, 0.90),
            pct(&mut mg, 1.0),
            pct(&mut mr, 0.50),
            pct(&mut mr, 0.90),
            pct(&mut mr, 1.0),
            pct(&mut md, 0.50),
            pct(&mut md, 0.90),
            pct(&mut md, 1.0),
        );

        // THREE COLUMNS. Each threshold set applied to the quantity it was written for.
        let sets: [(&str, [f32; 4]); 3] = [
            ("code today (never changed)", [20.0, 500.0, 5_000.0, 50_000.0]),
            ("author's proposal (parked)", [10.0, 100.0, 1_000.0, 8_000.0]),
            ("re-anchored on GEOMETRIC", {
                // Anchored the way the author anchored his: ship just under the measured max,
                // then decade steps down, so every class is populated.
                let mx = pct(&mut mg.clone(), 1.0).max(1.0);
                [mx / 1000.0, mx / 100.0, mx / 10.0, mx / 1.5]
            }),
        ];
        eprintln!(
            "\nCLASS SPLIT over all {} reaches (non-nav / small boat / barge / ship)",
            segs.len()
        );
        for (label, th) in &sets {
            for (qlabel, q) in [("runoff-eq", &runoff), ("geometric", &geo)] {
                let mut c = [0usize; 4];
                for &v in q.iter() {
                    let k = if v >= th[3] {
                        3
                    } else if v >= th[2] {
                        2
                    } else if v >= th[1] {
                        1
                    } else {
                        0
                    };
                    c[k] += 1;
                }
                eprintln!(
                    "  {:<28} on {:<10} [{:>6.0}/{:>6.0}/{:>6.0}/{:>6.0}] -> non-nav {:>6} | \
                     boat {:>5} | barge {:>4} | ship {:>4}",
                    label, qlabel, th[0], th[1], th[2], th[3], c[0], c[1], c[2], c[3]
                );
            }
        }

        // What a DISCHARGE-based set would need, anchored on real rivers.
        eprintln!(
            "\nIF THE THRESHOLDS WERE IN m3/s (the quantity a boat actually cares about)\n  \
             measured mouth discharge: p50 {:.2} | p90 {:.2} | max {:.2} m3/s\n  \
             Earth anchors: Thames ~65 m3/s (barge), Seine ~560 (ship-capable), Rhine ~2300.\n  \
             So this continent's LARGEST river ({:.1} m3/s) is {:.1}x smaller than the Thames.",
            pct(&mut md.clone(), 0.50),
            pct(&mut md.clone(), 0.90),
            pct(&mut md.clone(), 1.0),
            pct(&mut md.clone(), 1.0),
            65.0 / pct(&mut md.clone(), 1.0).max(1e-6)
        );
        let _ = h;
    }

    /// Is the absence of barge/ship rivers a THRESHOLD problem or a CONTINENT problem? The
    /// arid-hot band (25 deg, span 10) is the driest production config. Measure the same
    /// continent under the humid band and see whether discharge moves by the order needed.
    #[test]
    #[ignore]
    fn navigability_humid_2048() {
        navigability_report(2048, 45.0);
    }

    #[test]
    #[ignore]
    fn navigability_columns_2048() {
        navigability_report(2048, 25.0);
    }

    #[test]
    #[ignore]
    fn navigability_columns_8192() {
        navigability_report(8192, 25.0);
    }

    #[test]
    #[ignore]
    fn fragment_areas_2048() {
        fragment_report(2048, 25.0);
    }

    #[test]
    #[ignore]
    fn fragment_areas_8192() {
        fragment_report(8192, 25.0);
    }

    #[test]
    #[ignore]
    fn microscope_semantics_arid_hot_2048() {
        semantics_before_after(2048, 25.0);
    }

    #[test]
    #[ignore]
    fn microscope_semantics_arid_hot_8192() {
        semantics_before_after(8192, 25.0);
    }

    #[test]
    #[ignore]
    fn network_fragmentation_arid_hot() {
        let a = report(2048, 25.0);
        let b = report(8192, 25.0);
        eprintln!(
            "\n============  CROSS-RESOLUTION  ============\n\
             mouth uv        2048^2 ({:.4},{:.4})   8192^2 ({:.4},{:.4})\n\
             GEOMETRIC km2   2048^2 {:.0}   8192^2 {:.0}   (ratio {:.2})\n\
             EFFECTIVE km2   2048^2 {:.0}   8192^2 {:.0}   (ratio {:.2})",
            a.uv.0,
            a.uv.1,
            b.uv.0,
            b.uv.1,
            a.geom_km2,
            b.geom_km2,
            b.geom_km2 / a.geom_km2.max(1e-6),
            a.effective_km2,
            b.effective_km2,
            b.effective_km2 / a.effective_km2.max(1e-6),
        );
    }
}

/// ADR Finding 123 — **the identity guard, end to end through `run_hd`**: the viz's own pipeline, with
/// the HdParams the workspace sends for seed 1 at 8192² (auto framing, C-2/C-3/C-3b/H-1 on,
/// latitude 45°, span 40°, ratio 7.5), must read `Match` against the bench's reference for the
/// delivered world and for "C2 /10 col (défaut)" (the definition, basin base). A `NoReference` here means the viz's configuration is not
/// the bench's; a `Mismatch` means the same configuration yields another world.
///
/// Run: cargo test -p ymir-viz --release f123_viz_guard -- --ignored --nocapture
/// ADR Finding 151-0 — the viz starts as the author's production, which is the 6/6 guard's world: C-2, C-3 and C-3b
/// checked, H-1 not. The guard's literal (`f123_viz_guard::run`) and these defaults must not drift apart.
#[cfg(test)]
mod lake_outline {
    use super::*;

    /// ADR Finding 157-V, rule 13 -- the outline: a 3 x 3 lake in a 7 x 7 map gets its 8 border cells as the inner
    /// line and the 12 shore cells around it (4-adjacent) as the outer line; its centre is not outlined. Negative
    /// control first: a map without a lake has no outline. Two touching lakes each get their inner line.
    #[test]
    fn the_outline_rings_each_lake_inside_and_out() {
        let (w, h) = (7usize, 7usize);
        assert!(lake_outline_mask(&vec![0u32; w * h], w, h).iter().all(|&v| v == 0), "negative control: no lake, no outline");
        let mut lm = vec![0u32; w * h];
        for y in 2..5 {
            for x in 2..5 {
                lm[y * w + x] = 7;
            }
        }
        let m = lake_outline_mask(&lm, w, h);
        assert_eq!(m.iter().filter(|&&v| v == 1).count(), 8, "the inner line");
        assert_eq!(m.iter().filter(|&&v| v == 2).count(), 12, "the outer line, 4-adjacent shore cells");
        assert_eq!(m[3 * w + 3], 0, "the lake's centre is not outlined");
        // two lakes side by side
        let mut two = vec![0u32; w * h];
        two[3 * w + 2] = 1;
        two[3 * w + 3] = 2;
        let m2 = lake_outline_mask(&two, w, h);
        assert_eq!((m2[3 * w + 2], m2[3 * w + 3]), (1, 1), "each lake gets its own inner line");
        // the paint leaves non-outline pixels alone
        let mut rgba = vec![100u8; w * h * 4];
        draw_lake_outline(&mut rgba, &m);
        assert_eq!(&rgba[0..4], &[100, 100, 100, 100]);
        assert_eq!(&rgba[(2 * w + 2) * 4..(2 * w + 2) * 4 + 3], &[0, 229, 255]);
    }
}

#[cfg(test)]
mod startup_defaults {
    #[test]
    fn the_workspace_starts_with_the_guards_closures() {
        let ws = super::WorkspaceState::default();
        assert!(ws.volcanism, "C-2 volcanism is on in the guard's world");
        assert!(ws.lithology, "C-3 lithology is on in the guard's world");
        assert!(ws.fracture, "C-3b fracture is on in the guard's world");
        assert!(!ws.infiltration, "H-1 infiltration is off in the guard's world");
    }
}

#[cfg(test)]
mod f123_viz_guard {
    use super::*;
    use crate::bridge::c1::C1RunSpec;
    use crate::bridge::c1::events::C1Event;
    use crate::bridge::c1::hd::{HdParams, run_hd};
    use crossbeam_channel::bounded;
    use std::sync::atomic::AtomicBool;
    use ymir_core::tectonics_c1::bench_guard::GuardStatus;
    use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
    use ymir_core::tectonics_c1::closures::infiltration::InfiltrationConfig;
    use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
    use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;
    use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

    const PSEED: u64 = 10_481_999_410_520_546_993;

    fn run(
        label: &str,
        slope: Option<f32>,
        valley: Option<ValleyConstruction>,
    ) -> (GuardStatus, GuardStatus) {
        let spec = C1RunSpec { seed: PSEED, ..C1RunSpec::default() };
        // the workspace's own literal (workspace.rs, the "Générer HD" handler), seed 1 auto-framed
        let params = HdParams {
            target_size: 8192,
            latitude_deg: 45.0,
            domain_km: 400.0,
            manual_offset: Some([6.0 / 64.0, 37.0 / 64.0]),
            stream_power: true,
            closures: true,
            cross_rill: false,
            slope_floor_s_eq: slope,
            valley_construction: valley,
            cross_rill_d: 0.40,
            mfd: true,
            mfd_p: 2.0,
            base_level_m: None,
            base_level_off: false,
            geo_scale_ratio: 7.5,
            latitude_span_deg: Some(40.0),
            export_dir: None,
            volcanism: Some(VolcanismConfig { enabled: true, ..Default::default() }),
            lithology: Some(LithologyConfig {
                enabled: true,
                soft_multiplier: 10.0,
                volcanic_multiplier: 3.0,
                rift_age_threshold: 1.0,
            }),
            fracture: Some(FractureConfig {
                enabled: true,
                amplitude: 6.0,
                decay_km: 25.0,
                ..Default::default()
            }),
            infiltration: Some(InfiltrationConfig { enabled: true, ..Default::default() }),
            emit_tectonic_labels: false,
            geology_rules: None,
        };
        let (tx, rx) = bounded(256);
        let cancel = Arc::new(AtomicBool::new(false));
        std::thread::spawn(move || run_hd(&spec, &params, &tx, &cancel));
        while let Ok(e) = rx.recv() {
            match e {
                C1Event::HdCompleted { result, .. } => {
                    eprintln!("   {label}: field {:?} · lakes {:?}", result.bench_guard, result.lake_guard);
                    return (result.bench_guard.clone(), result.lake_guard.clone());
                }
                C1Event::HdFailed { error } => panic!("HD run failed: {error}"),
                _ => {}
            }
        }
        panic!("the HD worker hung up")
    }

    /// ADR Finding 125 — all SIX states the menu offers, each with the params the workspace sends
    /// for it (age k × 1, W(k) γ 0; the hybrid modes force the closure at 0.024, "C1 nue" follows
    /// the age box, off by default), must read `Match`.
    #[test]
    #[ignore]
    fn f123_viz_guard() {
        eprintln!("\n==========  Finding 123/125 . the guard through run_hd, six states  ==========");
        let k = F121_AGE_K;
        let states: [(&str, Option<f32>, Option<ValleyConstruction>); 6] = [
            ("livré", None, None),
            ("A1+B2", Some(0.024), None),
            ("C1 nue", None, Some(ValleyConstruction::new(k, None))),
            ("C2 /10 col (défaut)", Some(0.024), Some(ValleyConstruction::new(k, Some(0.1)))),
            ("C2 /3", Some(0.024), Some(ValleyConstruction::new(k, Some(1.0 / 3.0)))),
            ("C2 /10 niveau mer", Some(0.024), Some(ValleyConstruction::f121(k, Some(0.1)))),
        ];
        let mut verdicts = Vec::new();
        for (label, slope, valley) in states {
            verdicts.push((label, run(label, slope, valley)));
        }
        let say = |v: &GuardStatus| if matches!(v, GuardStatus::Match { .. }) { "= banc" } else { "≠ / non gardé" };
        for (label, (f, l)) in &verdicts {
            eprintln!("   {label:<22} field {} · lakes {}", say(f), say(l));
        }
        // ADR Finding 134 — six states, field AND lakes
        for (label, (f, l)) in verdicts {
            assert!(matches!(f, GuardStatus::Match { .. }), "{label} field: {f:?}");
            assert!(matches!(l, GuardStatus::Match { .. }), "{label} lakes: {l:?}");
        }
    }
}

/// ADR Finding 133 (visual validation) — a READ-ONLY bench on the viz's own path: `run_hd` itself (its protected
/// breach, its climate, its drainage), OFF and ON (the extended `lake_base`), with `f123_viz_guard`'s literal. The
/// lake counts, the new lakes' crops, the 3 based lakes with the largest upstream |Δz|, a control crop certified by
/// zero changed cells, every lake per crop, and the badges. Declared before the run
/// (`f133v_declared.md`). Test-only: no production code changes.
///
/// Run: cargo test -p ymir-viz --release f133v_viz_path -- --ignored --nocapture
#[cfg(test)]
mod f133v_bench {
    use super::*;
    use crate::bridge::c1::C1RunSpec;
    use crate::bridge::c1::events::C1Event;
    use crate::bridge::c1::hd::{HdParams, HdResult, run_hd};
    use crossbeam_channel::bounded;
    use std::collections::{HashMap, HashSet};
    use std::sync::atomic::AtomicBool;
    use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
    use ymir_core::tectonics_c1::closures::infiltration::InfiltrationConfig;
    use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
    use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
    use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;
    use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, LakeBase, ValleyConstruction};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

    const PSEED: u64 = 10_481_999_410_520_546_993;

    fn run(valley: ValleyConstruction) -> Arc<HdResult> {
        run_with(valley, None)
    }

    /// [`run`] with the `.ymir` export written to `export` (ADR Finding 157-B5).
    fn run_with(valley: ValleyConstruction, export: Option<std::path::PathBuf>) -> Arc<HdResult> {
        let spec = C1RunSpec { seed: PSEED, ..C1RunSpec::default() };
        // f123_viz_guard's literal (the workspace's own), for "C2 /10 col (défaut)"
        let params = HdParams {
            target_size: 8192,
            latitude_deg: 45.0,
            domain_km: 400.0,
            manual_offset: Some([6.0 / 64.0, 37.0 / 64.0]),
            stream_power: true,
            closures: true,
            cross_rill: false,
            slope_floor_s_eq: Some(0.024),
            valley_construction: Some(valley),
            cross_rill_d: 0.40,
            mfd: true,
            mfd_p: 2.0,
            base_level_m: None,
            base_level_off: false,
            geo_scale_ratio: 7.5,
            latitude_span_deg: Some(40.0),
            export_dir: export,
            volcanism: Some(VolcanismConfig { enabled: true, ..Default::default() }),
            lithology: Some(LithologyConfig {
                enabled: true,
                soft_multiplier: 10.0,
                volcanic_multiplier: 3.0,
                rift_age_threshold: 1.0,
            }),
            fracture: Some(FractureConfig { enabled: true, amplitude: 6.0, decay_km: 25.0, ..Default::default() }),
            infiltration: Some(InfiltrationConfig { enabled: true, ..Default::default() }),
            emit_tectonic_labels: false,
            geology_rules: None,
        };
        let (tx, rx) = bounded(256);
        let cancel = Arc::new(AtomicBool::new(false));
        std::thread::spawn(move || run_hd(&spec, &params, &tx, &cancel));
        while let Ok(e) = rx.recv() {
            match e {
                C1Event::HdCompleted { result, .. } => return result,
                C1Event::HdFailed { error } => panic!("HD run failed: {error}"),
                _ => {}
            }
        }
        panic!("the HD worker hung up")
    }

    /// ADR Finding 157-B2 / B3 / B5 — the lake base alone, on `run_hd` itself (the guard's C2 /10 col literal), OFF and ON
    /// extended: the lake counts and badges, the new lakes' crops (F133v's rule, 615 cells), the control crop (zero
    /// changed cells, else |Δz| ≤ 1 m declared), the lakes per crop, the hillshade crops with the lakes' outline
    /// (`layer_color_image`'s own buffer, north up), and the two `.ymir` exports for Living Landz. Declared in
    /// `docs/reports/lakes_gorges/f157_close/f157_declared.md`.
    ///
    /// Run: cargo test -p ymir-viz --release f157_viz -- --ignored --nocapture
    #[test]
    #[ignore]
    fn f157_viz() {
        use ymir_core::tectonics_c1::bench_guard::{field_hash, lake_fingerprint};
        let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        let out = root.join("docs/reports/lakes_gorges/f157_close");
        let ex = root.join("exports/f157");
        for d in ["off", "on"] {
            std::fs::create_dir_all(ex.join(d)).expect("the export folder");
        }
        eprintln!("\n==========  Finding 157 B . the lake base alone on the viz path (run_hd), OFF / ON extended  ==========");
        let ss = SteinSteinParams::default();
        let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
        let t = std::time::Instant::now();
        let off = run_with(base, Some(ex.join("off")));
        let t_off = t.elapsed().as_secs_f64();
        let t = std::time::Instant::now();
        let on = run_with(ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..base }, Some(ex.join("on")));
        let t_on = t.elapsed().as_secs_f64();
        let (w, h) = (off.width, off.height);
        let n = w * h;
        let cell_km2 = off.km_per_cell * off.km_per_cell;
        let (lm_off, lm_on) = (&off.drainage.lake_map, &on.drainage.lake_map);
        eprintln!("   lakes on the viz path: OFF **{}** · ON **{}** (F135-W: 26 / 34)", off.drainage.lakes.len(), on.drainage.lakes.len());
        eprintln!(
            "   field hashes: OFF {:016x} · ON {:016x} · badges: field OFF {:?} · ON {:?} · lakes OFF {:?} · ON {:?}",
            field_hash(&off.eroded),
            field_hash(&on.eroded),
            off.bench_guard,
            on.bench_guard,
            off.lake_guard,
            on.lake_guard
        );
        eprintln!("   lake fingerprints: OFF {} · ON {}", lake_fingerprint(lm_off, &off.drainage.lakes), lake_fingerprint(lm_on, &on.drainage.lakes));
        eprintln!("   run_hd wall time (the cache may HIT): OFF {t_off:.0} s · ON {t_on:.0} s");
        for d in ["off", "on"] {
            let mut files: Vec<String> = std::fs::read_dir(ex.join(d))
                .map(|r| r.flatten().map(|e| format!("{} ({:.1} MB)", e.file_name().to_string_lossy(), e.metadata().map(|m| m.len() as f64 / 1e6).unwrap_or(0.0))).collect())
                .unwrap_or_default();
            files.sort();
            eprintln!("   export {d}: {} · {:?}", ex.join(d).canonicalize().map(|p| p.display().to_string()).unwrap_or_default(), files);
        }
        let zo: Vec<f32> = off.eroded.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let zn: Vec<f32> = on.eroded.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let diff: Vec<bool> = (0..n).map(|i| off.eroded.data[i] != on.eroded.data[i]).collect();
        eprintln!("   conditioned cells that differ OFF vs ON: {}", diff.iter().filter(|&&b| b).count());
        let half = 307usize;
        let crop = |cx: usize, cy: usize| -> (usize, usize) { (cx.clamp(half, w - 1 - half) - half, cy.clamp(half, h - 1 - half) - half) };
        let lakes_in = |lm: &[u32], (x0, y0): (usize, usize)| -> usize {
            let mut st = HashSet::new();
            for y in y0..y0 + 2 * half + 1 {
                for x in x0..x0 + 2 * half + 1 {
                    let l = lm[y * w + x];
                    if l != 0 {
                        st.insert(l);
                    }
                }
            }
            st.len()
        };
        // the crops: the new lakes of ON
        let mut cells_of: HashMap<u32, Vec<usize>> = HashMap::new();
        for i in 0..n {
            if lm_on[i] != 0 {
                cells_of.entry(lm_on[i]).or_default().push(i);
            }
        }
        let mut crops: Vec<(String, usize, usize)> = Vec::new();
        let mut ids: Vec<&u32> = cells_of.keys().collect();
        ids.sort();
        for &id in ids {
            let cells = &cells_of[&id];
            let under = cells.iter().filter(|&&i| lm_off[i] != 0).count();
            if (under as f32) < 0.5 * cells.len() as f32 {
                let cx = (cells.iter().map(|&i| (i % w) as f64).sum::<f64>() / cells.len() as f64).round() as usize;
                let cy = (cells.iter().map(|&i| (i / w) as f64).sum::<f64>() / cells.len() as f64).round() as usize;
                crops.push((format!("lake{id}"), cx, cy));
                eprintln!("      new lake {id} ({:.1} km²) centre ({cx}, {cy})", cells.len() as f32 * cell_km2);
            }
        }
        // the control crop: zero changed cells, else every cell |Δz| ≤ 1 m (declared), the most land
        let sat = |f: &dyn Fn(usize) -> bool| -> Vec<u64> {
            let mut s = vec![0u64; (w + 1) * (h + 1)];
            for y in 0..h {
                for x in 0..w {
                    s[(y + 1) * (w + 1) + x + 1] = f(y * w + x) as u64 + s[y * (w + 1) + x + 1] + s[(y + 1) * (w + 1) + x] - s[y * (w + 1) + x];
                }
            }
            s
        };
        let rect = |s: &[u64], x0: usize, y0: usize| {
            let (x1, y1) = (x0 + 2 * half, y0 + 2 * half);
            s[(y1 + 1) * (w + 1) + x1 + 1] + s[y0 * (w + 1) + x0] - s[y0 * (w + 1) + x1 + 1] - s[(y1 + 1) * (w + 1) + x0]
        };
        let sl = sat(&|i| off.eroded.data[i] > 0.5);
        let mut control: Option<(String, usize, usize)> = None;
        for (rule, test) in [("zero changed cells", Box::new(|i: usize| diff[i]) as Box<dyn Fn(usize) -> bool>), ("|Δz| ≤ 1 m (declared threshold)", Box::new(|i: usize| (zn[i] - zo[i]).abs() > 1.0))] {
            let sd = sat(&*test);
            let mut best: Option<(u64, usize, usize)> = None;
            for y0 in (0..h - 2 * half).step_by(64) {
                for x0 in (0..w - 2 * half).step_by(64) {
                    if rect(&sd, x0, y0) == 0 {
                        let land = rect(&sl, x0, y0);
                        if best.is_none_or(|b| land > b.0) {
                            best = Some((land, x0, y0));
                        }
                    }
                }
            }
            match best {
                Some((land, x0, y0)) => {
                    eprintln!("   control crop ({rule}): origin ({x0}, {y0}) · land {:.0} % of the window", 100.0 * land as f64 / ((2 * half + 1) * (2 * half + 1)) as f64);
                    control = Some(("control".to_string(), x0 + half, y0 + half));
                    break;
                }
                None => eprintln!("   control crop ({rule}): NONE at 615 × 615"),
            }
        }
        if let Some(c) = control {
            crops.push(c);
        }
        // the lakes per crop, and the hillshade with the outline (the viz's own buffer), north up
        let ov = TectonicOverlays::default();
        let (img_off, _) = layer_color_image(&off, HdLayer::Relief, &RiverCellMap::from_drainage(&off.drainage), false, ov, ReliefView::Shade { outline: true }, None);
        let (img_on, _) = layer_color_image(&on, HdLayer::Relief, &RiverCellMap::from_drainage(&on.drainage), false, ov, ReliefView::Shade { outline: true }, None);
        let save = |img: &egui::ColorImage, (x0, y0): (usize, usize), name: &str| {
            let side = (2 * half + 1) as u32;
            let mut im = image::RgbaImage::new(side, side);
            // the image is north up: data row y is image row h − 1 − y
            for yy in 0..side as usize {
                for xx in 0..side as usize {
                    let (x, y) = (x0 + xx, h - 1 - (y0 + 2 * half - yy));
                    let c = img.pixels[y * w + x];
                    im.put_pixel(xx as u32, yy as u32, image::Rgba([c.r(), c.g(), c.b(), 255]));
                }
            }
            im.save(out.join(name)).expect("png");
        };
        eprintln!("   the crops (615², data origin x0, y0 south-first; images north up):");
        for (name, cx, cy) in &crops {
            let o = crop(*cx, *cy);
            let ch = (o.1..o.1 + 2 * half + 1).flat_map(|y| (o.0..o.0 + 2 * half + 1).map(move |x| y * w + x)).filter(|&k| diff[k]).count();
            eprintln!("      {name:<12} centre ({cx}, {cy}) · origin {o:?} · lakes OFF **{}** · ON **{}** · changed cells {ch}", lakes_in(lm_off, o), lakes_in(lm_on, o));
            save(&img_off, o, &format!("f157_{name}_off.png"));
            save(&img_on, o, &format!("f157_{name}_on.png"));
        }
        eprintln!("\n==========  end Finding 157 B (viz)  ==========\n");
    }

    #[test]
    #[ignore]
    fn f133v_viz_path() {
        eprintln!("\n==========  F133 visual validation . the viz path (run_hd), OFF / ON extended  ==========");
        let ss = SteinSteinParams::default();
        let k = F121_AGE_K;
        let base = ValleyConstruction::new(k, Some(0.1));
        let off = run(base);
        let on = run(ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..base });
        let (w, h) = (off.width, off.height);
        let n = w * h;
        let cell_km = off.km_per_cell;
        let cell_km2 = cell_km * cell_km;
        let (lm_off, lm_on) = (&off.drainage.lake_map, &on.drainage.lake_map);
        // 1 · counts, 4 · badges
        eprintln!(
            "   1 · lakes on the viz path: OFF **{}** · ON **{}** (the bench assembly: 25 / 33)",
            off.drainage.lakes.len(),
            on.drainage.lakes.len()
        );
        eprintln!("   4 · badges: OFF {:?} · ON {:?}", off.bench_guard, on.bench_guard);
        let zo: Vec<f32> = off.eroded.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let zn: Vec<f32> = on.eroded.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let diff: Vec<bool> = (0..n).map(|i| off.eroded.data[i] != on.eroded.data[i]).collect();
        eprintln!("   cells whose conditioned z differs OFF vs ON: {}", diff.iter().filter(|&&b| b).count());
        let half = 307i64;
        let crop = |cx: i64, cy: i64| -> (usize, usize, usize, usize) {
            let (x0, x1) = ((cx - half).max(0) as usize, ((cx + half) as usize).min(w - 1));
            let (y0, y1) = ((cy - half).max(0) as usize, ((cy + half) as usize).min(h - 1));
            (x0, x1, y0, y1)
        };
        let lakes_in = |lm: &[u32], b: (usize, usize, usize, usize)| -> usize {
            let mut s = HashSet::new();
            for y in b.2..=b.3 {
                for x in b.0..=b.1 {
                    let l = lm[y * w + x];
                    if l != 0 {
                        s.insert(l);
                    }
                }
            }
            s.len()
        };
        let row = |name: String, cx: i64, cy: i64| {
            let b = crop(cx, cy);
            let ch = (b.2..=b.3).flat_map(|y| (b.0..=b.1).map(move |x| (x, y))).filter(|&(x, y)| diff[y * w + x]).count();
            eprintln!(
                "      {name:<40} centre ({cx}, {cy}) · crop x {}–{} · y {}–{} · lakes OFF **{}** · ON **{}** · cells with Δz ≠ 0: {ch}",
                b.0,
                b.1,
                b.2,
                b.3,
                lakes_in(lm_off, b),
                lakes_in(lm_on, b)
            );
        };
        // 2(a) · the new lakes of ON
        eprintln!("   2a · the lakes new in ON (< 50 % of their cells under any OFF lake):");
        let mut cells_of: HashMap<u32, Vec<usize>> = HashMap::new();
        for i in 0..n {
            if lm_on[i] != 0 {
                cells_of.entry(lm_on[i]).or_default().push(i);
            }
        }
        let mut news: Vec<(u32, f64, f64, usize)> = Vec::new();
        for (&id, cells) in &cells_of {
            let under = cells.iter().filter(|&&i| lm_off[i] != 0).count();
            if (under as f32) < 0.5 * cells.len() as f32 {
                let cx = cells.iter().map(|&i| (i % w) as f64).sum::<f64>() / cells.len() as f64;
                let cy = cells.iter().map(|&i| (i / w) as f64).sum::<f64>() / cells.len() as f64;
                news.push((id, cx, cy, cells.len()));
            }
        }
        news.sort_by(|a, b| a.0.cmp(&b.0));
        for &(id, cx, cy, c) in &news {
            row(format!("new lake {id} ({:.1} km²)", c as f32 * cell_km2), cx.round() as i64, cy.round() as i64);
        }
        // 2(b) · the based lakes with the largest upstream |Δz|: the OFF flow path to the first ON lake
        let dir = &off.drainage.flow.direction;
        let mut first = vec![u32::MAX; n];
        let mut path = Vec::new();
        for s in 0..n {
            if first[s] != u32::MAX {
                continue;
            }
            path.clear();
            let mut c = s;
            let v;
            loop {
                if first[c] != u32::MAX {
                    v = first[c];
                    break;
                }
                if lm_on[c] != 0 {
                    v = lm_on[c];
                    break;
                }
                path.push(c);
                let d = dir[c];
                if d == DIR_NONE {
                    v = 0;
                    break;
                }
                let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
                if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
                    v = 0;
                    break;
                }
                c = y as usize * w + x as usize;
            }
            for &p in &path {
                first[p] = v;
            }
        }
        let mut per: HashMap<u32, (f64, Vec<f32>)> = HashMap::new();
        for i in 0..n {
            if diff[i] && lm_on[i] == 0 && first[i] != 0 {
                let dz = (zn[i] - zo[i]).abs();
                let e = per.entry(first[i]).or_insert((0.0, Vec::new()));
                e.0 += dz as f64 * cell_km2 as f64 * 1e-3;
                e.1.push(dz);
            }
        }
        let mut ranked: Vec<(u32, f64, f32, usize)> = per
            .into_iter()
            .map(|(id, (s, mut v))| {
                v.sort_by(f32::total_cmp);
                let p90 = v[((v.len() - 1) as f32 * 0.9) as usize];
                (id, s, p90, v.len())
            })
            .collect();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
        eprintln!("   2b · the based lakes (ON final lakes) with the largest upstream |Δz|, statistic Σ|Δz| × cell area:");
        for &(id, s, p90, c) in ranked.iter().take(3) {
            let cells = &cells_of[&id];
            let cx = cells.iter().map(|&i| (i % w) as f64).sum::<f64>() / cells.len() as f64;
            let cy = cells.iter().map(|&i| (i / w) as f64).sum::<f64>() / cells.len() as f64;
            eprintln!("      lake {id}: Σ|Δz|·A **{s:.2} km³** over {c} upstream cells · p90 |Δz| {p90:.1} m");
            row(format!("based lake {id}"), cx.round() as i64, cy.round() as i64);
        }
        // 2(c) · the control crop: zero changed cells, the most land
        let sat = |f: &dyn Fn(usize) -> bool| -> Vec<u64> {
            let mut s = vec![0u64; (w + 1) * (h + 1)];
            for y in 0..h {
                for x in 0..w {
                    s[(y + 1) * (w + 1) + x + 1] =
                        f(y * w + x) as u64 + s[y * (w + 1) + x + 1] + s[(y + 1) * (w + 1) + x] - s[y * (w + 1) + x];
                }
            }
            s
        };
        let sd = sat(&|i| diff[i]);
        let sl = sat(&|i| off.eroded.data[i] > 0.5);
        let rect = |s: &[u64], x0: usize, y0: usize, x1: usize, y1: usize| {
            s[(y1 + 1) * (w + 1) + x1 + 1] + s[y0 * (w + 1) + x0] - s[y0 * (w + 1) + x1 + 1] - s[(y1 + 1) * (w + 1) + x0]
        };
        let side = 615usize;
        let mut best: Option<(u64, usize, usize)> = None;
        for y0 in (0..h - side).step_by(64) {
            for x0 in (0..w - side).step_by(64) {
                let (x1, y1) = (x0 + side - 1, y0 + side - 1);
                if rect(&sd, x0, y0, x1, y1) == 0 {
                    let land = rect(&sl, x0, y0, x1, y1);
                    if best.is_none_or(|b| land > b.0) {
                        best = Some((land, x0, y0));
                    }
                }
            }
        }
        match best {
            Some((land, x0, y0)) => {
                eprintln!(
                    "   2c · the control crop: 0 changed cells (certified), land {:.0} % of the window",
                    100.0 * land as f64 / (side * side) as f64
                );
                row("control (no Δz)".to_string(), (x0 + side / 2) as i64, (y0 + side / 2) as i64);
            }
            None => eprintln!("   2c · NO 615 × 615 window with zero changed cells"),
        }
        eprintln!("\n==========  end F133 visual validation  ==========\n");
    }

    /// The table of passage: the lakes the viz path lists beyond the bench assembly's count, by type.
    ///
    /// Run: cargo test -p ymir-viz --release f133v_passage -- --ignored --nocapture
    #[test]
    #[ignore]
    fn f133v_passage() {
        use ymir_core::tectonics_c1::drainage::LakeType;
        let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
        for (label, vc) in [("OFF", base), ("ON extended", ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..base })] {
            let r = run(vc);
            let w = r.width;
            let mut by_type: HashMap<String, usize> = HashMap::new();
            for l in &r.drainage.lakes {
                *by_type.entry(format!("{:?}", l.lake_type)).or_insert(0) += 1;
            }
            eprintln!("   {label}: {} lakes · by type {:?}", r.drainage.lakes.len(), by_type);
            for l in r.drainage.lakes.iter().filter(|l| matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
                let cells: Vec<usize> = (0..r.drainage.lake_map.len()).filter(|&i| r.drainage.lake_map[i] == l.base.id).collect();
                let (cx, cy) = if cells.is_empty() {
                    (f64::NAN, f64::NAN)
                } else {
                    (
                        cells.iter().map(|&i| (i % w) as f64).sum::<f64>() / cells.len() as f64,
                        cells.iter().map(|&i| (i / w) as f64).sum::<f64>() / cells.len() as f64,
                    )
                };
                eprintln!(
                    "      crater lake {} ({:?}) · {:.2} km² · level {:.1} m · centre ({cx:.0}, {cy:.0})",
                    l.base.id, l.lake_type, l.area_km2, l.level_m
                );
            }
            let active = r.volcanoes.iter().filter(|c| c.active).count();
            eprintln!("      active craters: {active}");
        }
    }

    /// ADR Finding 135-V — the layers' own buffers (the SAME `layer_color_image` the viz draws), cropped on lake 2:
    /// the difference ON extended − OFF at ±1000 m and ±100 m, the hillshade of OFF and of ON. NOT a screen grab.
    /// The worlds' field hashes are printed: a view must not change them.
    ///
    /// Run: F135_DIR=<dir> cargo test -p ymir-viz --release f135v_relief_capture -- --ignored --nocapture
    #[test]
    #[ignore]
    fn f135v_relief_capture() {
        use ymir_core::tectonics_c1::bench_guard::field_hash;
        let dir = std::env::var("F135_DIR").expect("F135_DIR");
        let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
        let off = run(base);
        let on = run(ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..base });
        let (h0, h1) = (field_hash(&off.eroded), field_hash(&on.eroded));
        let (w, h) = (on.width, on.height);
        let (cx, cy) = (4019usize, h - 1 - 2839usize); // lake 2, data (4019, 2839), shown north-up
        let half = 307usize;
        let save = |img: &egui::ColorImage, name: &str| {
            let side = (2 * half + 1) as u32;
            let mut out = image::RgbaImage::new(side, side);
            for yy in 0..side as usize {
                for xx in 0..side as usize {
                    let (x, y) = (cx + xx - half, cy + yy - half);
                    let c = img.pixels[y * w + x];
                    out.put_pixel(xx as u32, yy as u32, image::Rgba([c.r(), c.g(), c.b(), 255]));
                }
            }
            let p = std::path::Path::new(&dir).join(name);
            out.save(&p).expect("png");
            eprintln!("   wrote {}", p.display());
        };
        let rm_off = RiverCellMap::from_drainage(&off.drainage);
        let rm_on = RiverCellMap::from_drainage(&on.drainage);
        let ov = TectonicOverlays::default();
        for (sat, name) in [(1000.0f32, "v_diff_lake2_1000m.png"), (100.0, "v_diff_lake2_100m.png")] {
            let (img, st) = layer_color_image(&on, HdLayer::Relief, &rm_on, false, ov, ReliefView::Diff(&off, sat), None);
            eprintln!("   difference ON − OFF at ±{sat} m: (saturated, ≠ 0) = {st:?}");
            save(&img, name);
        }
        let (img, _) = layer_color_image(&off, HdLayer::Relief, &rm_off, false, ov, ReliefView::Shade { outline: false }, None);
        save(&img, "v_shade_lake2_OFF.png");
        let (img, _) = layer_color_image(&on, HdLayer::Relief, &rm_on, false, ov, ReliefView::Shade { outline: false }, None);
        save(&img, "v_shade_lake2_ON.png");
        let (img, _) = layer_color_image(&off, HdLayer::Relief, &rm_off, false, ov, ReliefView::Hypso, None);
        save(&img, "v_hypso_lake2_OFF.png");
        assert_eq!(field_hash(&off.eroded), h0, "a view wrote the world");
        assert_eq!(field_hash(&on.eroded), h1, "a view wrote the world");
        eprintln!("   field hashes unchanged by the views: OFF {h0:016x} · ON {h1:016x} · guard OFF {:?}", off.bench_guard);
    }

    /// ADR Finding 134, item 0.3 — the viz path's lakes, lake by lake, OFF and ON (extended lake base),
    /// in `bench_guard::lake_listing`'s format, written to `F134_LAKES_DIR`: the bench's
    /// `f134_lake_guard` writes the same two files from its own assembly, and the identity is a diff.
    ///
    /// Run: F134_LAKES_DIR=<dir> cargo test -p ymir-viz --release f134v_lake_listing -- --ignored --nocapture
    #[test]
    #[ignore]
    fn f134v_lake_listing() {
        use ymir_core::tectonics_c1::bench_guard::{lake_fingerprint, lake_listing};
        let dir = std::env::var("F134_LAKES_DIR").expect("F134_LAKES_DIR");
        let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
        for (name, vc) in [
            ("viz_lakes_OFF.txt", base),
            ("viz_lakes_ON.txt", ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..base }),
        ] {
            let r = run(vc);
            let fp = lake_fingerprint(&r.drainage.lake_map, &r.drainage.lakes);
            eprintln!("   {name}: {} lakes · {fp} · field {:?} · lakes {:?}", r.drainage.lakes.len(), r.bench_guard, r.lake_guard);
            let mut txt = format!("# viz run_hd · lake guard {:?} · {fp}\n", r.lake_guard);
            for l in lake_listing(&r.drainage.lake_map, &r.drainage.lakes) {
                txt.push_str(&l);
                txt.push('\n');
            }
            std::fs::write(std::path::Path::new(&dir).join(name), txt).expect("write listing");
        }
    }
}

/// ADR Finding 124-3 — the shipped aggregation's three gestures on a hand-built world (permanent,
/// fast). 12×6: the open sea is column 11, an enclosed below-sea basin sits at (2,4). W0 runs along
/// row 1 to the sea; S1 is the basin's SPILLWAY, running down W0's valley to the same mouth (the
/// Watercourse × Spillway pair that shared 297 cells on C2/10); I2 dies in the basin; W3 runs along
/// row 5 on its own. Expected: two objects — {W0, S1, I2}, a RIVER — and {W3}; every reach in
/// exactly one; the basins disjoint and summing to at most the land.
#[cfg(test)]
mod f124_objects {
    use super::*;
    use ymir_core::terrain::flow::{FlowResult, RiverNetwork, RiverSegment};

    #[test]
    fn one_reach_one_object_and_the_basins_are_a_partition() {
        let (w, h) = (12usize, 6usize);
        let n = w * h;
        let mut z = vec![0.8f32; n];
        for y in 0..h {
            z[y * w + 11] = 0.2; // the open sea (touches the grid edge)
        }
        z[4 * w + 2] = 0.3; // an enclosed below-sea basin
        let eroded = GridF32::from_vec(w, h, z);
        let seg = |points: Vec<(u32, u32)>| RiverSegment {
            points,
            strahler_order: 1,
            avg_flow: 1.0,
            max_flow: 1.0,
            basin_id: 0,
            upstream: vec![],
            downstream: None,
        };
        let w0 = seg((5..=10).map(|x| (x, 1)).collect());
        let s1 = seg([(3, 4), (4, 3), (5, 2)].into_iter().chain((6..=10).map(|x| (x, 1))).collect());
        let i2 = seg(vec![(0, 4), (1, 4)]);
        let w3 = seg((0..=10).map(|x| (x, 5)).collect());
        let segs = vec![w0, s1, i2, w3];
        let m = segs.len();
        let mut direction = vec![DIR_NONE; n];
        for y in 0..h {
            for x in 0..11 {
                if (x, y) != (2, 4) {
                    direction[y * w + x] = 2; // east, to the sea
                }
            }
        }
        let drainage = ymir_core::tectonics_c1::drainage::C1DrainageResult {
            flow: FlowResult {
                filled: GridF32::new(w, h, 0.0),
                direction,
                accumulation: GridF32::new(w, h, 1.0),
                basins: vec![0; n],
                num_basins: 1,
            },
            segment_drainage_km2: vec![10.0, 50.0, 5.0, 8.0],
            segment_catchment_cells: vec![1.0; m],
            segment_navigability: vec![Navigability::NonNavigable; m],
            segment_discharge_m3s: vec![1.0, 5.0, 0.5, 0.8],
            segment_width_m: vec![1.0; m],
            segment_profile_m: segs.iter().map(|s| vec![0.0; s.points.len()]).collect(),
            segment_discharge_profile_m3s: segs.iter().map(|s| vec![1.0; s.points.len()]).collect(),
            segment_kind: vec![
                SegmentKind::Watercourse,
                SegmentKind::Spillway,
                SegmentKind::Watercourse,
                SegmentKind::Watercourse,
            ],
            segment_source_lake: vec![None; m],
            rivers: RiverNetwork { segments: segs },
            lakes: vec![],
            lake_map: vec![0; n],
            width: w,
            height: h,
        };
        let hd = HdResult {
            width: w,
            height: h,
            land_topology: ymir_core::tectonics_c1::land_topology::land_topology(&eroded, 0.5),
            eroded,
            temperature: GridF32::new(w, h, 10.0),
            precipitation: GridF32::new(w, h, 0.0),
            drainage,
            biomes: vec![Biome::Ocean; n],
            volcanoes: vec![],
            tectonic: None,
            edifices: vec![],
            sample_origin: [0.0, 0.0],
            sample_size: 1.0,
            bench_guard: ymir_core::tectonics_c1::bench_guard::GuardStatus::NoReference { viz: String::new() },
            lake_guard: ymir_core::tectonics_c1::bench_guard::GuardStatus::NoReference { viz: String::new() },
            gorge_falls: Vec::new(),
            rivers_ll: Vec::new(),
            geology: None,
            km_per_cell: 1.0,
        };
        let wcs = aggregate_watercourses(&hd, w as f32, 1.0);
        assert_eq!(wcs.len(), 2, "two systems: {:?}", wcs.iter().map(|x| &x.segments).collect::<Vec<_>>());
        let mut owner = vec![usize::MAX; m];
        for (o, x) in wcs.iter().enumerate() {
            for &s in &x.segments {
                assert_eq!(owner[s], usize::MAX, "reach {s} in two objects");
                owner[s] = o;
            }
        }
        assert!(owner.iter().all(|&o| o != usize::MAX), "every reach in an object");
        assert_eq!(owner[0], owner[1], "the spillway down W0's valley is W0's system");
        assert_eq!(owner[2], owner[1], "the reach dying in the basin is chained to its spillway");
        assert_ne!(owner[3], owner[0], "W3 is its own system");
        let merged = &wcs[owner[0]];
        assert_eq!(merged.kind, SegmentKind::Watercourse, "a system with a watercourse is a river");
        // the trunk climbs from the mouth to the spillway (largest own area) and across the basin
        assert_eq!(merged.trunk, vec![2, 1], "trunk source → mouth: I2, then the spillway");
        // the viz's cell area: (domain / w) · (domain / h), 1 × 2 km here
        let cell_km2 = (w as f32 / w as f32) * (w as f32 / h as f32);
        let land = (0..n).filter(|&k| k % w != 11).count() as f32 * cell_km2;
        let sum: f32 = wcs.iter().map(|x| x.catchment_km2).sum();
        assert!(sum <= land + 1e-3, "the basins are a partition: Σ {sum} > land {land}");
        assert!(wcs.iter().all(|x| x.catchment_km2 > 0.0), "every system drains an area");
    }
}
