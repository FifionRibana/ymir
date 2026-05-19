//! Biome phase — Whittaker classification. **Stub.**
//!
//! Depends on the Climate phase (also stubbed). Surfaces the M3
//! roadmap pointer in the right panel; no compute / render systems
//! ship in this commit.

use bevy_egui::egui;

pub fn draw_section(ui: &mut egui::Ui) {
    ui.add_space(4.0);
    ui.heading("Biome (stub)");
    ui.label(
        egui::RichText::new(
            "Biome classification (Whittaker envelope on temperature × \
             precipitation) is scheduled for milestone M3. Depends on \
             the Climate phase landing first.",
        )
        .small()
        .weak(),
    );
    ui.add_space(6.0);
    ui.label(
        egui::RichText::new(
            "Planned inputs:\n\
             • Temperature grid (°C)\n\
             • Precipitation grid (mm/yr)\n\
             \n\
             Planned outputs:\n\
             • `BiomeId: u8` per cell (tundra, taiga, temperate \
               forest, savanna, desert, …)\n\
             • Color-coded biome map for the central sprite",
        )
        .small()
        .weak(),
    );
}
