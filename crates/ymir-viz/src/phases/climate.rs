//! Climate phase — temperature + precipitation. **Stub.**
//!
//! `ymir_core::climate::{temperature, precipitation, biomes}` are
//! placeholder modules with no compute logic. The pre-Phase-8h
//! sunset binary already exposed Climate as a coming-soon entry,
//! and we re-introduce the same shape here so the navigation
//! toolbar has a consistent surface — the panel surfaces the M3
//! roadmap pointer rather than throwing a "not implemented"
//! error.

use bevy_egui::egui;

pub fn draw_section(ui: &mut egui::Ui) {
    ui.add_space(4.0);
    ui.heading("Climate (stub)");
    ui.label(
        egui::RichText::new(
            "Temperature + precipitation grids are scheduled for \
             milestone M3. The compute modules in \
             `ymir_core::climate::*` are placeholders and will be \
             populated alongside the panel sliders.",
        )
        .small()
        .weak(),
    );
    ui.add_space(6.0);
    ui.label(
        egui::RichText::new(
            "Planned inputs:\n\
             • Hydrology heightmap (already populated)\n\
             • Latitude band\n\
             • Lapse rate (°C / km)\n\
             • Atmospheric water-content baseline\n\
             \n\
             Planned outputs:\n\
             • `Temperature: GridF32` (°C)\n\
             • `Precipitation: GridF32` (mm/yr)",
        )
        .small()
        .weak(),
    );
}
