use bevy_egui::egui;

use crate::state::TerrainStats;

pub fn draw(ui: &mut egui::Ui, stats: &TerrainStats) {
    egui::CollapsingHeader::new("Statistics").default_open(true).show(ui, |ui| {
        stat_row(ui, "Grid size", &format!("{}x{}", stats.grid_width, stats.grid_height));
        stat_row(ui, "Resolution", &format!("{:.0} m/px", stats.meters_per_pixel));
        stat_row(ui, "Peak altitude", &format!("{:.0} m", stats.peak_altitude));
        stat_row(ui, "Min altitude", &format!("{:.0} m", stats.min_altitude));
        stat_row(ui, "Land ratio", &format!("{:.1}%", stats.land_ratio * 100.0));
        stat_row(ui, "River segments", &format!("{}", stats.river_segments));
        stat_row(ui, "Lakes", &format!("{}", stats.lake_count));
    });
}

fn stat_row(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(label).color(egui::Color32::GRAY));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.monospace(value);
        });
    });
}
