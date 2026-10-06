//! ADR Finding 152-B6 — the geology's export: `geologie_roches.png` (u8 class codes, full resolution),
//! `favorabilite_<id>.png` (u8 0–100, the ¼ grid) and `geologie.json` (the legend, the files, the rules' SHA-256).
//! Every PNG is 8-bit grey, row 0 = the SOUTHERNMOST row (the container's orientation invariant). Format:
//! `docs/geology_format.md`.

use serde::Serialize;

use super::rocks::ROCK_CLASSES;
use super::rules::GeologyRules;
use super::zoning::{CHANCE_FACTOR, Zoning};

/// The geology export's own format version.
pub const GEOLOGY_FORMAT_VERSION: &str = "0.1.0";

/// An 8-bit grey PNG (lossless).
#[must_use]
pub fn png_u8(g: &[u8], w: usize, h: usize) -> Vec<u8> {
    use image::ImageEncoder;
    let mut out = Vec::new();
    image::codecs::png::PngEncoder::new(&mut out).write_image(g, w as u32, h as u32, image::ExtendedColorType::L8).expect("PNG encoding of an in-memory u8 grid");
    out
}

#[derive(Serialize)]
struct ClassEntry {
    id: u8,
    key: &'static str,
    name_fr: &'static str,
    color: String,
    /// 2 hard, 1 medium, 0 soft.
    hardness: u8,
}

#[derive(Serialize)]
struct ResourceEntry {
    id: String,
    name_fr: String,
    color: String,
    file: String,
}

#[derive(Serialize)]
struct RulesInfo {
    format: String,
    version: String,
    sha256: String,
}

#[derive(Serialize)]
struct Manifest {
    format_version: &'static str,
    value: &'static str,
    orientation: &'static str,
    rules: RulesInfo,
    rocks: RocksInfo,
    favourability: FavInfo,
}

#[derive(Serialize)]
struct RocksInfo {
    file: &'static str,
    width: usize,
    height: usize,
    classes: Vec<ClassEntry>,
}

#[derive(Serialize)]
struct FavInfo {
    width: usize,
    height: usize,
    factor: usize,
    resources: Vec<ResourceEntry>,
}

fn hex(c: [u8; 3]) -> String {
    format!("#{:02X}{:02X}{:02X}", c[0], c[1], c[2])
}

/// The files of the export: `(file name, bytes)`, the manifest last.
#[must_use]
pub fn export_files(rocks: &[u8], w: usize, h: usize, z: &Zoning, rules: &GeologyRules, rules_sha256: &str) -> Vec<(String, Vec<u8>)> {
    let mut files = vec![("geologie_roches.png".to_string(), png_u8(rocks, w, h))];
    let mut resources = Vec::new();
    for r in &z.resources {
        let file = format!("favorabilite_{}.png", r.id);
        files.push((file.clone(), png_u8(&r.fav, z.w4, z.h4)));
        resources.push(ResourceEntry { id: r.id.clone(), name_fr: r.name_fr.clone(), color: hex(r.color), file });
    }
    let m = Manifest {
        format_version: GEOLOGY_FORMAT_VERSION,
        value: "relative favourability 0-100 per resource (not a probability, not a deposit)",
        orientation: "row-major, y = 0 = the SOUTH edge: row 0 of each PNG is the southernmost row",
        rules: RulesInfo { format: rules.format.clone(), version: rules.version.clone(), sha256: rules_sha256.to_string() },
        rocks: RocksInfo {
            file: "geologie_roches.png",
            width: w,
            height: h,
            classes: ROCK_CLASSES.iter().map(|c| ClassEntry { id: c.id, key: c.key, name_fr: c.name_fr, color: hex(c.color), hardness: c.hardness }).collect(),
        },
        favourability: FavInfo { width: z.w4, height: z.h4, factor: CHANCE_FACTOR, resources },
    };
    files.push(("geologie.json".to_string(), serde_json::to_vec_pretty(&m).expect("the geology manifest serialises")));
    files
}
