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
use ymir_core::tectonics_c1::drainage::{LakeType, Navigability};
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
    /// A below-sea ENCLOSED cell (water_class == 2) that is NOT an inventoried lake — a sub-sea
    /// evaporative sink (Finding 36). Distinct from Sea: tagging it "→ mer" would contradict
    /// water_class (the authority) and reinstate the very altitude proxy Finding 31 removed.
    SubSeaSink,
    Unknown,
}

/// An aggregated WATERCOURSE (a trunk + its tributaries), assembled UI-side from the
/// exported per-segment drainage — NOT recomputed. `segments` are every segment that
/// drains to the same terminal; `trunk` is the source→sink main stem (the max-discharge
/// choice at each confluence). Everything read from `C1DrainageResult`.
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
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
enum HdLayer {
    #[default]
    Relief,
    Drainage,
    Precipitation,
    Temperature,
    Biomes,
}

impl HdLayer {
    fn label(self) -> &'static str {
        match self {
            HdLayer::Relief => "Relief",
            HdLayer::Drainage => "Drainage",
            HdLayer::Precipitation => "Précipitation",
            HdLayer::Temperature => "Température",
            HdLayer::Biomes => "Biomes",
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
const FRIEZE: [(HdPhase, &str, Option<HdLayer>); 7] = [
    (HdPhase::Tectonic, "Tectonique", None),
    (HdPhase::Relief, "Relief", None),
    (HdPhase::Erosion, "Érosion", Some(HdLayer::Relief)),
    (HdPhase::Drainage, "Drainage", Some(HdLayer::Drainage)),
    (HdPhase::Climate, "Précip.", Some(HdLayer::Precipitation)),
    (HdPhase::Climate, "Temp.", Some(HdLayer::Temperature)),
    (HdPhase::Biomes, "Biomes", Some(HdLayer::Biomes)),
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
    /// EXPERIMENTAL (ADR 0001, Finding 11): MFD incision — dendritic valleys, no solver.
    mfd: bool,
    mfd_p: f32,
    /// EXPERIMENTAL: FBM amplitude_base for the striation ladder (0.16 = production).
    fbm_amplitude: f32,
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
    /// segment index → watercourse index (for hover-select on the map).
    seg_to_wc: Vec<usize>,
    inspect_tab: InspectTab,
    texture: Option<egui::TextureHandle>,
    tex_layer: Option<HdLayer>,
    /// Overlay state the cached `texture` was built with (rebuild when it flips).
    tex_overlay: bool,
    hover: Option<CellInspection>,
    hover_xy: Option<(usize, usize)>,
    zoom: f32,
    /// Active map tool: SELECT (microscope — hover/click inspects an entity) or PAN
    /// (hand — click-drag moves the view). Wheel/CTRL-wheel work in either tool.
    tool: MapTool,
    /// Pipeline frieze collapsed (hidden) to give the canvas more room.
    pipeline_collapsed: bool,
    /// Bird-eye minimap shown over the canvas.
    minimap: bool,
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
            cross_rill_d: 0.40,
            mfd: true,
            mfd_p: 2.0,
            fbm_amplitude: 0.04, // #190 — the production seed's amplitude (fine relief detail)
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
            seg_to_wc: Vec::new(),
            inspect_tab: InspectTab::Rivers,
            texture: None,
            tex_layer: None,
            tex_overlay: false,
            hover: None,
            hover_xy: None,
            zoom: 1.0,
            map_pan: [0.5, 0.5],
            tool: MapTool::Select,
            pipeline_collapsed: false,
            minimap: true,
        }
    }
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
    right_panel(ctx, &mut ws);
    central_panel(ctx, &bridge, &mut ws);
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
                            cross_rill_d: ws.cross_rill_d,
                            mfd: ws.mfd,
                            mfd_p: ws.mfd_p,
                            fbm_amplitude: Some(ws.fbm_amplitude as f64),
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
                                    &mut ws.cross_rill,
                                    egui::RichText::new("Anti-peigne (diffusion partout)").color(DIM2).size(11.0),
                                )
                                .on_hover_text(
                                    "STEP 2b (ADR 0001, Finding 9): diffusion de versant PARTOUT pour \
                                     amortir le rilling de Smith–Bretherton. Solveur GS non convergé.",
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
                            let vals = [0.16f32, 0.08, 0.04, 0.02, 0.01];
                            let cur =
                                vals.iter().position(|&v| (v - ws.fbm_amplitude).abs() < 1e-4).unwrap_or(2);
                            ui.horizontal(|ui| {
                                ui.label(egui::RichText::new("FBM amp").color(DIM2).size(11.0));
                                if let Some(i) = seg_row(ui, &["0.16", "0.08", "0.04", "0.02", "0.01"], cur) {
                                    ws.fbm_amplitude = vals[i];
                                }
                            });
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
        .exact_width(282.0)
        .frame(egui::Frame::default().fill(PANEL).inner_margin(0))
        .show(ctx, |ui| {
            ui.spacing_mut().item_spacing = egui::vec2(8.0, 6.0);
            block(ui, 15, 13, 13, |ui| {
                ui.horizontal(|ui| {
                    ui.label(
                        egui::RichText::new("Inspection de cellule")
                            .color(TEXT_BRIGHT)
                            .strong()
                            .size(12.5),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("›").on_hover_text("Replier").clicked() {
                            ws.inspector_open = false;
                        }
                    });
                });
            });
            ui.separator();
            match &ws.hover {
                Some(c) => {
                    let c = *c;
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        block(ui, 15, 14, 14, |ui| inspection(ui, &c));
                    });
                }
                None => {
                    ui.add_space(48.0);
                    ui.vertical_centered(|ui| {
                        ui.label(
                            egui::RichText::new("⊹")
                                .color(C::from_rgb(0x4a, 0x4a, 0x4a))
                                .size(22.0),
                        );
                        ui.add_space(6.0);
                        ui.label(
                            egui::RichText::new(if ws.current.is_some() {
                                "Survolez la carte pour inspecter\nles grandeurs d'une cellule."
                            } else {
                                "Générez un continent."
                            })
                            .color(C::from_rgb(0x7a, 0x7a, 0x7a))
                            .size(12.0),
                        );
                    });
                }
            }
        });
}

fn kv(ui: &mut egui::Ui, k: &str, v: String) {
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

fn inspection(ui: &mut egui::Ui, c: &CellInspection) {
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
    kv(ui, "Épaisseur crustale / S̃", "—".into()); // coarse — deferred
    kv(ui, "Type de plaque", "—".into());
    kv(ui, "Craton", "—".into());

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
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                let pct = progress
                                    .map(|(d, t)| format!(" {}%", if t > 0 { d * 100 / t } else { 0 }))
                                    .unwrap_or_default();
                                ui.label(
                                    egui::RichText::new(format!("{}{} …", p.label(), pct))
                                        .color(COPPER_BRIGHT)
                                        .size(10.5),
                                );
                                ui.spinner();
                            });
                        }
                    });
                    if !ws.pipeline_collapsed {
                        frieze(ui, &bridge.hd, ws);
                    }
                });
            // Inspection microscope (Finding 28) — entity lists + long profile / lake sheet.
            if ws.current.is_some() {
                egui::TopBottomPanel::bottom("inspection_dock")
                    .resizable(true)
                    .default_height(210.0)
                    .frame(egui::Frame::default().fill(PANEL).inner_margin(8))
                    .show_inside(ui, |ui| inspection_dock(ui, ws));
            }
            // Map viewport (#0A0A0A).
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
        target_size: ws.resolution,
        latitude_deg: ws.latitude,
        domain_km: ws.domain_km,
        manual_offset: None,
        stream_power: false, // preview is coarse-only; SP applies to the HD run
        closures: false,
        cross_rill: false,
        cross_rill_d: 0.40,
        mfd: false,
        mfd_p: 2.0,
        fbm_amplitude: None,
        geo_scale_ratio: 1.0,
        latitude_span_deg: None,
        export_dir: None,
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
    if ws.texture.is_none() || ws.tex_layer != Some(ws.layer) || ws.tex_overlay != ws.river_overlay
    {
        let layer = ws.layer;
        let overlay = ws.river_overlay;
        let t_tex = std::time::Instant::now();
        let img = {
            let hd = ws.current.as_ref().unwrap();
            let rm = ws.river_map.as_ref().unwrap();
            layer_color_image(hd, layer, rm, overlay)
        };
        ws.texture = Some(ui.ctx().load_texture("hd_map", img, egui::TextureOptions::NEAREST));
        ws.tex_layer = Some(ws.layer);
        ws.tex_overlay = ws.river_overlay;
        eprintln!("[HD timing] map texture rebuild {:.1}s (UI thread)", t_tex.elapsed().as_secs_f32());
    }
    let handle = ws.texture.as_ref().unwrap().clone();

    let z = ws.zoom.max(1.0);
    let uv_half = 0.5 / z;
    // Clamp the pan CENTRE so the zoom window stays inside the map; recentre at fit.
    if z > 1.0 {
        let (lo, hi) = (uv_half, 1.0 - uv_half);
        ws.map_pan[0] = ws.map_pan[0].clamp(lo, hi);
        ws.map_pan[1] = ws.map_pan[1].clamp(lo, hi);
    } else {
        ws.map_pan = [0.5, 0.5];
    }
    let avail = ui.available_size();
    let side = avail.x.min(avail.y).max(1.0);
    // Centre the square map in the canvas (fills the shorter dimension; pan/zoom to explore).
    ui.add_space(((avail.y - side) * 0.5).max(0.0));
    ui.horizontal(|ui| {
        ui.add_space(((avail.x - side) * 0.5).max(0.0));
        let (cx, cy) = (ws.map_pan[0], ws.map_pan[1]);
        let uv = egui::Rect::from_min_max(
            egui::pos2(cx - uv_half, cy - uv_half),
            egui::pos2(cx + uv_half, cy + uv_half),
        );
        // The canvas senses click+drag (pan with the hand tool / middle button) and wheel.
        let resp = ui.add(
            egui::Image::new(egui::load::SizedTexture::new(handle.id(), egui::vec2(side, side)))
                .uv(uv)
                .sense(egui::Sense::click_and_drag()),
        );

        // ── Canvas navigation (mock): drag-to-pan, CTRL+wheel zoom around the cursor,
        //    wheel = vertical scroll, SHIFT+wheel = horizontal (egui routes it to scroll.x).
        let uv_span = 2.0 * uv_half;
        let pan_drag = (ws.tool == MapTool::Pan && resp.dragged())
            || resp.dragged_by(egui::PointerButton::Middle);
        if pan_drag && z > 1.0 {
            let d = resp.drag_delta();
            ws.map_pan[0] -= d.x / side * uv_span; // drag right → view slides left
            ws.map_pan[1] -= d.y / side * uv_span;
        }
        if resp.hovered() {
            let (scroll, zoom_delta, mods) =
                ui.ctx().input(|i| (i.smooth_scroll_delta, i.zoom_delta(), i.modifiers));
            // CTRL+wheel (or pinch) → zoom. egui usually routes CTRL+wheel into zoom_delta; the
            // scroll.y fallback covers the case it does not.
            let zf = if (zoom_delta - 1.0).abs() > 1e-3 {
                zoom_delta
            } else if (mods.ctrl || mods.command) && scroll.y != 0.0 {
                1.0 + scroll.y / 400.0
            } else {
                1.0
            };
            if (zf - 1.0).abs() > 1e-4 {
                let new_z = (z * zf).clamp(1.0, 12.0);
                if let Some(p) = resp.hover_pos() {
                    // keep the uv under the cursor fixed while zooming.
                    let rel = p - resp.rect.min;
                    let cur_u =
                        uv.min.x + (rel.x / resp.rect.width()).clamp(0.0, 1.0) * uv_span;
                    let cur_v =
                        uv.min.y + (rel.y / resp.rect.height()).clamp(0.0, 1.0) * uv_span;
                    let k = 1.0 - z / new_z;
                    ws.map_pan[0] += (cur_u - ws.map_pan[0]) * k;
                    ws.map_pan[1] += (cur_v - ws.map_pan[1]) * k;
                }
                ws.zoom = new_z;
            } else if z > 1.0 && !mods.ctrl && !mods.command && (scroll.x != 0.0 || scroll.y != 0.0)
            {
                ws.map_pan[0] -= scroll.x / side * uv_span;
                ws.map_pan[1] -= scroll.y / side * uv_span;
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

        // Overlays (painter-drawn boxes, matching the mock's absolute divs).
        let hd = ws.current.clone().unwrap();
        overlay_chip(ui, resp.rect, ws.layer);
        legend_box(ui, resp.rect, ws.layer);

        // Hover → cell → inspect.
        ws.hover = None;
        ws.hover_xy = None;
        if let Some(pos) = resp.hover_pos() {
            let rel = pos - resp.rect.min;
            let fx =
                uv.min.x + (rel.x / resp.rect.width()).clamp(0.0, 0.999) * (uv.max.x - uv.min.x);
            let fy =
                uv.min.y + (rel.y / resp.rect.height()).clamp(0.0, 0.999) * (uv.max.y - uv.min.y);
            let cx = ((fx * hd.width as f32) as usize).min(hd.width - 1);
            // The texture is north-up mirrored (Finding 27); the DATA cell is south-first,
            // so the hovered data row is the mirror of the texture row.
            let tex_row = ((fy * hd.height as f32) as usize).min(hd.height - 1);
            let cy = (hd.height - 1) - tex_row;
            if let Some(rm) = ws.river_map.as_ref() {
                ws.hover = Some(inspect_cell(&hd, rm, cx, cy));
                ws.hover_xy = Some((cx, cy));
            }
            // Hover-select (Finding 28): only with the SELECT (microscope) tool — the hand tool
            // pans without disturbing the current selection. Selects in the ACTIVE tab.
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
                            if let Some(li) =
                                hd.drainage.lakes.iter().position(|l| l.base.id == id)
                            {
                                ws.selected_lake = Some(li);
                            }
                        }
                    }
                }
            }

            // Reticle: faint crosshair through the cell + a copper cell box.
            let pnt = ui.painter_at(resp.rect);
            let uvw = uv.max.x - uv.min.x;
            let uvh = uv.max.y - uv.min.y;
            let ncx = (cx as f32 + 0.5) / hd.width as f32;
            // Screen v of the data cell = the texture (north-up) row, i.e. its mirror.
            let ncy = ((hd.height - 1 - cy) as f32 + 0.5) / hd.height as f32;
            let sx = resp.rect.left() + ((ncx - uv.min.x) / uvw) * resp.rect.width();
            let sy = resp.rect.top() + ((ncy - uv.min.y) / uvh) * resp.rect.height();
            let bw = ((1.0 / hd.width as f32 / uvw) * resp.rect.width()).max(8.0);
            let bh = ((1.0 / hd.height as f32 / uvh) * resp.rect.height()).max(8.0);
            let faint = C::from_rgba_unmultiplied(0xC9, 0x85, 0x3F, 90);
            pnt.line_segment(
                [egui::pos2(sx, resp.rect.top()), egui::pos2(sx, resp.rect.bottom())],
                egui::Stroke::new(1.0, faint),
            );
            pnt.line_segment(
                [egui::pos2(resp.rect.left(), sy), egui::pos2(resp.rect.right(), sy)],
                egui::Stroke::new(1.0, faint),
            );
            pnt.rect_stroke(
                egui::Rect::from_center_size(egui::pos2(sx, sy), egui::vec2(bw, bh)),
                0.0,
                egui::Stroke::new(1.5, COPPER_BRIGHT),
                egui::StrokeKind::Middle,
            );

            // Hover coord readout (top-right) — the DATA cell (south-first indices).
            let txt = format!("x {cx}  y {cy}");
            let anchor = egui::pos2(resp.rect.right() - 8.0, resp.rect.top() + 8.0);
            let galley = ui.painter().layout_no_wrap(
                txt,
                egui::FontId::monospace(10.5),
                C::from_rgb(0x9a, 0x9a, 0x9a),
            );
            let bg = egui::Rect::from_min_size(
                egui::pos2(anchor.x - galley.size().x - 10.0, anchor.y),
                galley.size() + egui::vec2(12.0, 8.0),
            );
            ui.painter().rect_filled(bg, 6.0, C::from_rgba_unmultiplied(18, 18, 18, 210));
            ui.painter().galley(egui::pos2(bg.left() + 6.0, bg.top() + 4.0), galley, TEXT);
        }

        // Selection highlight (Finding 28) — painter-drawn (cheap, no texture rebake). The
        // north-up flip is honoured: a data cell's screen v uses its mirror row.
        {
            let (w, h) = (hd.width, hd.height);
            let uvw = uv.max.x - uv.min.x;
            let uvh = uv.max.y - uv.min.y;
            let to_screen = |x: u32, y: u32| -> egui::Pos2 {
                let ncx = (x as f32 + 0.5) / w as f32;
                let ncy = ((h as u32 - 1 - y) as f32 + 0.5) / h as f32;
                egui::pos2(
                    resp.rect.left() + ((ncx - uv.min.x) / uvw) * resp.rect.width(),
                    resp.rect.top() + ((ncy - uv.min.y) / uvh) * resp.rect.height(),
                )
            };
            let pnt = ui.painter_at(resp.rect);
            if ws.inspect_tab == InspectTab::Rivers {
                if let Some(wc) =
                    ws.selected_river.and_then(|i| ws.watercourses.as_ref().and_then(|w| w.get(i)))
                {
                    let gold = C::from_rgb(255, 210, 70);
                    for &s in &wc.segments {
                        let pts: Vec<egui::Pos2> = hd.drainage.rivers.segments[s]
                            .points
                            .iter()
                            .map(|&(px, py)| to_screen(px, py))
                            .collect();
                        if pts.len() >= 2 {
                            pnt.add(egui::Shape::line(pts, egui::Stroke::new(2.0, gold)));
                        }
                    }
                    // Mark the mouth.
                    pnt.circle_filled(
                        to_screen(wc.mouth_xy.0, wc.mouth_xy.1),
                        4.0,
                        sink_label(wc.sink).1,
                    );
                }
            } else if let Some(lk) = ws.selected_lake.and_then(|i| hd.drainage.lakes.get(i)) {
                let id = lk.base.id;
                let cyan = C::from_rgba_unmultiplied(120, 230, 240, 220);
                for y in 0..h {
                    for x in 0..w {
                        if hd.drainage.lake_map[y * w + x] == id {
                            // boundary cells only (cheap outline).
                            let edge =
                                [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)].iter().any(|&(dx, dy)| {
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
        }

        // Floating tool bar (top-left): select / hand + the overlay & minimap toggles.
        canvas_toolbar(ui, resp.rect, ws);
        // Bird-eye minimap (bottom-left) showing the current view window.
        if ws.minimap {
            minimap_overlay(ui, resp.rect, ws, &hd, uv);
        }
        // Zoom controls (bottom-right).
        zoom_controls(ui, resp.rect, ws);
        // ESC clears the current selection (deselect).
        if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
            ws.selected_river = None;
            ws.selected_lake = None;
        }
    });
}

/// Floating canvas tool bar (mock): SELECT (microscope) / PAN (hand) tools, plus the hydro
/// overlay and minimap toggles. Painted as a small pill at the canvas top-left.
fn canvas_toolbar(ui: &mut egui::Ui, rect: egui::Rect, ws: &mut WorkspaceState) {
    let bar = egui::Rect::from_min_size(
        egui::pos2(rect.left() + 10.0, rect.top() + 10.0),
        egui::vec2(0.0, 0.0),
    );
    let mut child = ui.new_child(
        egui::UiBuilder::new().max_rect(egui::Rect::from_min_size(bar.min, egui::vec2(360.0, 34.0))),
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
        p.rect_stroke(
            win,
            0.0,
            egui::Stroke::new(1.5, COPPER_BRIGHT),
            egui::StrokeKind::Middle,
        );
    }
    p.rect_stroke(mm, 4.0, egui::Stroke::new(1.0, C::from_rgb(0x2e, 0x2e, 0x2e)), egui::StrokeKind::Middle);
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

fn legend_box(ui: &mut egui::Ui, rect: egui::Rect, layer: HdLayer) {
    let p = ui.painter_at(rect);
    let items: Vec<(C, &str, &str)> = match layer {
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
    };
    let title = match layer {
        HdLayer::Relief => "RELIEF — HYPSOMÉTRIE",
        HdLayer::Drainage => "DRAINAGE — NAVIGABILITÉ",
        HdLayer::Precipitation => "PRÉCIPITATION (MM/AN)",
        HdLayer::Temperature => "TEMPÉRATURE",
        HdLayer::Biomes => "BIOMES",
    };
    let row_h = 16.0;
    let rows = if items.is_empty() { 2.0 } else { items.len() as f32 };
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
    if layer == HdLayer::Relief {
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
                *lbl,
                egui::FontId::proportional(11.0),
                C::from_rgb(0xcf, 0xcf, 0xcf),
            );
            if !sub.is_empty() {
                p.text(
                    egui::pos2(bpos.x + bw - 12.0, ry),
                    egui::Align2::RIGHT_TOP,
                    *sub,
                    egui::FontId::monospace(9.5),
                    C::from_rgb(0x77, 0x77, 0x77),
                );
            }
        }
    }
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
    child.horizontal(|ui| {
        if ui.small_button("−").clicked() {
            ws.zoom = (ws.zoom / 1.25).max(1.0);
        }
        ui.label(
            egui::RichText::new(format!("{:.0}%", ws.zoom * 100.0))
                .monospace()
                .color(C::from_rgb(0xbb, 0xbb, 0xbb))
                .size(11.0),
        );
        if ui.small_button("+").clicked() {
            ws.zoom = (ws.zoom * 1.25).min(9.0);
        }
        if ui.small_button("⤢").on_hover_text("Recadrer").clicked() {
            ws.zoom = 1.0;
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

// ── Layer → RGBA image (canonical palettes) ──────────────────────────────
fn layer_color_image(
    hd: &HdResult,
    layer: HdLayer,
    river_map: &RiverCellMap,
    overlay: bool,
) -> egui::ColorImage {
    let (w, h) = (hd.width, hd.height);
    let mut rgba = vec![0u8; w * h * 4];
    for k in 0..w * h {
        let c = match layer {
            HdLayer::Relief => relief_color(hd.eroded.data[k]),
            HdLayer::Precipitation => precip_color(precip_mm_per_year(hd.precipitation.data[k])),
            HdLayer::Temperature => temp_color(hd.temperature.data[k]),
            HdLayer::Biomes => {
                let [r, g, b] = hd.biomes[k].color();
                [r, g, b]
            }
            HdLayer::Drainage => drainage_color(hd, river_map, k),
        };
        rgba[k * 4] = c[0];
        rgba[k * 4 + 1] = c[1];
        rgba[k * 4 + 2] = c[2];
        rgba[k * 4 + 3] = 255;
    }
    if overlay {
        draw_river_overlay(&mut rgba, hd);
    }
    // Finding 27 — the field is stored y=0=SOUTH (the export/LL contract); the VIEWER
    // shows NORTH-UP, so mirror the texture rows here. The DATA, the export and the
    // inspector's cell lookup all stay south-first; only the pixels are flipped.
    flip_rows_rgba(&mut rgba, w, h);
    egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba)
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
fn aggregate_watercourses(hd: &HdResult, domain_km: f32, ratio: f32) -> Vec<Watercourse> {
    let d = &hd.drainage;
    let wc = ymir_core::lakes::connectivity::water_class(&hd.eroded, 0.5);
    let segs = &d.rivers.segments;
    let n = segs.len();
    let q = |i: usize| d.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
    // root(i): follow downstream to the terminal (mouth); guard against any cycle.
    let mut root = vec![0usize; n];
    for i in 0..n {
        let mut j = i;
        for _ in 0..=n {
            match segs[j].downstream {
                Some(k) if k < n => j = k,
                _ => break,
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
    let mut out = Vec::with_capacity(groups.len());
    for (r, members) in groups {
        // Trunk: from the mouth, climb the LONGEST tributary at each confluence.
        let mut trunk = Vec::new();
        let mut cur = r;
        for _ in 0..=n {
            trunk.push(cur);
            match segs[cur].upstream.iter().copied().filter(|&u| u < n).max_by_key(|&u| pathlen[u])
            {
                Some(next) => cur = next,
                None => break,
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
            discharge_m3s: q(r),
            width_mouth_m: d.segment_width_m.get(r).copied().unwrap_or(0.0),
            width_source_m: d.segment_width_m.get(source).copied().unwrap_or(0.0),
            catchment_km2: d.segment_drainage_km2.get(r).copied().unwrap_or(0.0),
            mouth_xy: (mx, my),
            sink,
            sink_lake_id,
            length_km: pts as f32 * km_per_cell * ratio,
        });
    }
    out.sort_by(|a, b| {
        b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap_or(std::cmp::Ordering::Equal)
    });
    out
}

fn sink_label(s: Sink) -> (&'static str, C) {
    match s {
        Sink::Sea => ("→ mer", C::from_rgb(0x5a, 0x9a, 0xd0)),
        Sink::ExoLake => ("→ lac exoréique", C::from_rgb(0x5a, 0x9a, 0xd0)),
        Sink::EndoLake => ("→ bassin endoréique", C::from_rgb(0x3a, 0xb0, 0xa0)),
        Sink::SubSeaSink => ("→ puits sous-marin (évaporatif)", C::from_rgb(0x3a, 0xb0, 0xa0)),
        Sink::Unknown => ("→ ?", WARN_ORANGE),
    }
}

/// TASK 2/3/4 — the inspection DOCK (Finding 28): entity lists (rivers as watercourses,
/// lakes) on the left with select-and-highlight, and the selected entity's detail (river
/// long profile OR lake sheet) on the right. A MICROSCOPE — read-only, nothing editable.
fn inspection_dock(ui: &mut egui::Ui, ws: &mut WorkspaceState) {
    let hd = match ws.current.clone() {
        Some(h) => h,
        None => return,
    };
    let wcs = ws.watercourses.clone();
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("MICROSCOPE").color(TEXT_BRIGHT).strong().size(11.5));
        ui.add_space(10.0);
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
        // Finding 32 note: displayed discharges/areas are geo-ratio COMPRESSED (×ratio² on
        // catchment), while the water balance reasons on REAL quantities. Flag it so the two
        // are not read as homogeneous.
        if (ws.geo_scale_ratio - 1.0).abs() > 0.01 {
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    egui::RichText::new(format!(
                        "débits ×{:.0} compressés (échelle {:.1})",
                        ws.geo_scale_ratio * ws.geo_scale_ratio,
                        ws.geo_scale_ratio
                    ))
                    .color(WARN_ORANGE)
                    .size(10.0),
                );
            });
        }
    });
    ui.separator();
    let mut nav: Option<NavAction> = None;
    ui.columns(2, |cols| {
        // ── Left: entity list ──
        egui::ScrollArea::vertical().id_salt("entity_list").show(&mut cols[0], |ui| {
            match ws.inspect_tab {
                InspectTab::Rivers => {
                    if let Some(wcs) = &wcs {
                        for (i, wc) in wcs.iter().enumerate() {
                            let (sl, sc) = sink_label(wc.sink);
                            let txt = format!(
                                "#{:<3} S{} · {:.0} m³/s · {} trib · {}",
                                i + 1,
                                wc.order,
                                wc.discharge_m3s,
                                wc.tributaries,
                                sl
                            );
                            if ui
                                .selectable_label(
                                    ws.selected_river == Some(i),
                                    egui::RichText::new(txt).color(sc).size(11.0),
                                )
                                .clicked()
                            {
                                ws.selected_river = Some(i);
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
                            LakeType::Endorheic => {
                                ("endoréique (salé)", C::from_rgb(0x3a, 0xb0, 0xa0))
                            }
                        };
                        let txt = format!(
                            "#{} · {:.0} km² · {:.0} m · {}",
                            lk.base.id, lk.area_km2, lk.level_m, ty
                        );
                        if ui
                            .selectable_label(
                                ws.selected_lake == Some(i),
                                egui::RichText::new(txt).color(tc).size(11.0),
                            )
                            .clicked()
                        {
                            ws.selected_lake = Some(i);
                        }
                    }
                }
            }
        });
        // ── Right: detail ──
        let out = egui::ScrollArea::vertical().id_salt("entity_detail").show(
            &mut cols[1],
            |ui| -> Option<NavAction> {
                match ws.inspect_tab {
                    InspectTab::Rivers => match (ws.selected_river, &wcs) {
                        (Some(i), Some(w)) if i < w.len() => river_profile_panel(ui, &hd, &w[i]),
                        _ => {
                            ui.add_space(20.0);
                            ui.label(
                                egui::RichText::new("Sélectionnez un cours d'eau.")
                                    .color(DIM)
                                    .size(11.0),
                            );
                            None
                        }
                    },
                    InspectTab::Lakes => match ws.selected_lake {
                        Some(i) if i < hd.drainage.lakes.len() => lake_sheet_panel(
                            ui,
                            &hd,
                            i,
                            ws.domain_km,
                            wcs.as_deref().map(|v| v.as_slice()),
                        ),
                        _ => {
                            ui.add_space(20.0);
                            ui.label(
                                egui::RichText::new("Sélectionnez un lac.").color(DIM).size(11.0),
                            );
                            None
                        }
                    },
                }
            },
        );
        nav = out.inner;
    });
    // Apply the clickable-chain jump (Finding 29): follow the hydrological graph in one click.
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
fn river_profile_panel(ui: &mut egui::Ui, hd: &HdResult, wc: &Watercourse) -> Option<NavAction> {
    let d = &hd.drainage;
    let (sl, sc) = sink_label(wc.sink);
    let mut nav = None;
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
    if is_spillway {
        for &s in &wc.trunk {
            if let Some(pr) = d.segment_profile_m.get(s) {
                elev.extend(pr.iter().copied());
            }
        }
    } else if let Some(&(sx, sy)) =
        wc.trunk.first().and_then(|&s| d.rivers.segments[s].points.first())
    {
        let mut k = sy as usize * w + sx as usize;
        for _ in 0..(w + h) {
            elev.push(c1_altitude_norm_to_metres(hd.eroded.data[k], &ss));
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
    let (rect, _) = ui.allocate_exact_size(
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
        let r =
            ((seg.strahler_order as i32 - 2).max(0)).min(2).max(orphan as i32).max(is_spillway as i32);
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
