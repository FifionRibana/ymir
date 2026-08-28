//! UI module — single-engine (C1) after the v2 sunset.
//!
//! The legacy v2 panels (top bar, parameter editor, metrics dashboard, workflow
//! panel) and the per-phase pipeline toolbar were removed with the v2 engine.
//! The C1 engine renders its own egui windows (control panel + cell inspector in
//! `visualization::c1_plugin`). This module now only installs the egui plugin and
//! the shared copper/bronze theme; the dockable 3-zone layout + pipeline frieze
//! land in the upcoming UI rewrite.

use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPlugin, EguiPrimaryContextPass, egui};

pub mod workspace;

pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(EguiPlugin::default());
        app.add_systems(EguiPrimaryContextPass, configure_egui_style);
    }
}

/// The Ymir dark theme — the Claude Design mock's palette (copper `#B87333` /
/// bright `#C9853F`, panels `#1B1B1B`, ink `#0d0d0d`, text `#E0E0E0`). Applied
/// once; the workspace (step d2) builds its chrome on top.
fn configure_egui_style(mut contexts: EguiContexts, mut done: Local<bool>) {
    if *done {
        return;
    }
    let Ok(ctx) = contexts.ctx_mut() else { return };
    *done = true;

    use egui::Color32 as C;
    let copper = C::from_rgb(0xB8, 0x73, 0x33);
    let copper_bright = C::from_rgb(0xC9, 0x85, 0x3F);
    let bronze = C::from_rgb(0x9A, 0x72, 0x4A);
    let panel = C::from_rgb(0x1B, 0x1B, 0x1B);
    let panel2 = C::from_rgb(0x16, 0x16, 0x16);
    let field = C::from_rgb(0x11, 0x11, 0x11);
    let border = C::from_rgb(0x2A, 0x2A, 0x2A);
    let text = C::from_rgb(0xE0, 0xE0, 0xE0);

    let mut style = (*ctx.style()).clone();
    let v = &mut style.visuals;
    v.dark_mode = true;
    v.override_text_color = Some(text);
    v.panel_fill = panel;
    v.window_fill = panel;
    v.extreme_bg_color = field; // text edits / sliders trough
    v.faint_bg_color = panel2;
    v.window_stroke = egui::Stroke::new(1.0, border);
    v.widgets.noninteractive.bg_fill = panel;
    v.widgets.noninteractive.bg_stroke = egui::Stroke::new(1.0, border);
    v.widgets.inactive.bg_fill = C::from_rgb(0x20, 0x20, 0x20);
    v.widgets.inactive.weak_bg_fill = C::from_rgb(0x20, 0x20, 0x20);
    v.widgets.hovered.bg_fill = bronze;
    v.widgets.hovered.weak_bg_fill = C::from_rgb(0x2a, 0x2a, 0x2a);
    v.widgets.active.bg_fill = copper;
    v.widgets.active.weak_bg_fill = copper;
    v.selection.bg_fill = copper.linear_multiply(0.4);
    v.selection.stroke = egui::Stroke::new(1.0, copper_bright);
    v.widgets.active.fg_stroke = egui::Stroke::new(1.0, copper_bright);
    let r = egui::CornerRadius::same(5);
    v.widgets.noninteractive.corner_radius = r;
    v.widgets.inactive.corner_radius = r;
    v.widgets.hovered.corner_radius = r;
    v.widgets.active.corner_radius = r;

    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(8.0, 4.0);
    ctx.set_style(style);
}
