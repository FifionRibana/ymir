//! UI module — single-engine (C1) after the v2 sunset.
//!
//! The legacy v2 panels (top bar, parameter editor, metrics dashboard, workflow
//! panel) and the per-phase pipeline toolbar were removed with the v2 engine.
//! The C1 engine renders its own egui windows (control panel + cell inspector in
//! `visualization::c1_plugin`). This module now only installs the egui plugin and
//! the shared copper/bronze theme; the dockable 3-zone layout + pipeline frieze
//! land in the upcoming UI rewrite.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts, EguiPlugin, EguiPrimaryContextPass};

pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(EguiPlugin::default());
        app.add_systems(EguiPrimaryContextPass, configure_egui_style);
    }
}

/// Copper/bronze accent theme over egui's default dark visuals. Reusable by the
/// upcoming UI rewrite (the brand accents: copper `#B87333`, bronze `#A0724A`).
fn configure_egui_style(mut contexts: EguiContexts, mut done: Local<bool>) {
    if *done {
        return;
    }
    let Ok(ctx) = contexts.ctx_mut() else { return };
    *done = true;

    let mut style = (*ctx.style()).clone();
    let copper = egui::Color32::from_rgb(0xB8, 0x73, 0x33);
    let bronze = egui::Color32::from_rgb(0xA0, 0x72, 0x4A);

    style.visuals.widgets.active.bg_fill = copper;
    style.visuals.widgets.hovered.bg_fill = bronze;
    style.visuals.selection.bg_fill = copper.linear_multiply(0.4);
    style.visuals.selection.stroke = egui::Stroke::new(1.0, copper);
    style.visuals.widgets.active.fg_stroke = egui::Stroke::new(1.0, copper);

    ctx.set_style(style);
}
