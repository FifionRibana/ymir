//! LAKE ABUNDANCE diagnostic (read-only, the export IS the verdict). Subject 2: why do
//! humid maps fill every basin to its sill (northern-Scotland, not France/Spain)? Reads
//! the shipped `lakes.json` and reports the population by size, depth and type, the
//! exorheic (would-drain) fraction, and estimates what a sill-incision pass would drain.
//! Design question, not a defect — measure, do not propose.
//!
//! Run: cargo test -p ymir-core --test lake_abundance_diagnosis -- --ignored --nocapture

use std::path::Path;

use serde::Deserialize;

const DIR: &str = "../../exports/seed10481999410520546993_8192.ymir";

#[derive(Deserialize)]
struct Lake {
    base: Base,
    depth_m: f32,
    area_km2: f32,
    lake_type: String,
}
#[derive(Deserialize)]
struct Base {
    id: u32,
    max_depth: f32,
    shallow: bool,
}

fn pct(v: &mut [f32], p: f32) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() as f32 * p) as usize).min(v.len().saturating_sub(1))]
}

#[test]
#[ignore]
fn lake_abundance() {
    if !Path::new(DIR).exists() {
        eprintln!("export missing — skip");
        return;
    }
    let lakes: Vec<Lake> =
        serde_json::from_slice(&std::fs::read(Path::new(DIR).join("lakes.json")).unwrap()).unwrap();
    eprintln!("\n=== LAKE ABUNDANCE — export {DIR} ===");
    eprintln!("lakes: {}", lakes.len());

    // Type split.
    let mut types = std::collections::BTreeMap::new();
    for l in &lakes {
        *types.entry(l.lake_type.clone()).or_insert(0u32) += 1;
    }
    eprintln!("by type: {types:?}");
    let below_sea = lakes.iter().filter(|l| l.base.id >= 1_000_001).count();
    let shallow = lakes.iter().filter(|l| l.base.shallow).count();
    eprintln!(
        "below-sea basin ids (≥1_000_001): {below_sea}  | shallow (max_depth<min_depth, FBM-residue-like): {shallow}"
    );

    // Size + depth distributions, and where the water AREA concentrates.
    let total_area: f32 = lakes.iter().map(|l| l.area_km2).sum();
    let mut areas: Vec<f32> = lakes.iter().map(|l| l.area_km2).collect();
    let mut depths: Vec<f32> = lakes.iter().map(|l| l.depth_m).collect();
    eprintln!("\ntotal lake area: {total_area:.0} km²");
    eprintln!(
        "area_km2 p50/p90/p99/max: {:.2}/{:.1}/{:.0}/{:.0}",
        pct(&mut areas.clone(), 0.5),
        pct(&mut areas.clone(), 0.9),
        pct(&mut areas.clone(), 0.99),
        pct(&mut areas, 1.0)
    );
    eprintln!(
        "depth_m  p50/p90/p99/max: {:.1}/{:.0}/{:.0}/{:.0}",
        pct(&mut depths.clone(), 0.5),
        pct(&mut depths.clone(), 0.9),
        pct(&mut depths.clone(), 0.99),
        pct(&mut depths, 1.0)
    );

    // Water-area concentration: what fraction of total lake area is in the biggest 10 %.
    let mut by_area: Vec<f32> = lakes.iter().map(|l| l.area_km2).collect();
    by_area.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let top10 = (lakes.len() as f32 * 0.1).ceil() as usize;
    let top_area: f32 = by_area.iter().take(top10).sum();
    eprintln!(
        "the largest {top10} lakes hold {:.0}% of the lake area ({:.0} of {:.0} km²)",
        100.0 * top_area / total_area,
        top_area,
        total_area
    );

    // Sill-incision estimate. An EXORHEIC lake has an outlet → its sill can incise (that
    // is how post-glacial lakes vanish); the frozen-sill model keeps it as a permanent
    // lake. An ENDORHEIC lake has no outlet (arid) and legitimately persists. So the
    // fraction a sill-incision pass would DRAIN ≈ the exorheic share; the residue is
    // endorheic + crater.
    let exo = lakes.iter().filter(|l| l.lake_type == "Exorheic").count();
    let endo = lakes.iter().filter(|l| l.lake_type == "Endorheic").count();
    let crater = lakes.iter().filter(|l| l.lake_type.starts_with("Crater")).count();
    let exo_area: f32 =
        lakes.iter().filter(|l| l.lake_type == "Exorheic").map(|l| l.area_km2).sum();
    eprintln!("\nsill-incision estimate (exorheic lakes = have an outlet → would incise/drain):");
    eprintln!(
        "  exorheic {exo} ({:.0}% by count, {:.0}% by area) | endorheic {endo} | crater {crater}",
        100.0 * exo as f32 / lakes.len() as f32,
        100.0 * exo_area / total_area
    );
    eprintln!(
        "  → a sill-incision pass would drain ~{:.0}% of the lake AREA; only the endorheic + crater\n    residue ({} lakes) legitimately persists. THIS is the France↔Scotland dial (climate sets the\n    exorheic share via net_evap; sill incision decides whether exorheic basins stay full).",
        100.0 * exo_area / total_area,
        endo + crater
    );
    eprintln!(
        "\nBLAST RADIUS (do NOT implement): sill incision partially drains exorheic lakes → the outlet\n becomes a river, lake footprints shrink/vanish, levels & regimes recompute, river mouths and the\n coastline move, biomes shift. It reopens the whole hydro chain stabilised over ~15 passes\n (drainage + lakes + rivers + coast + biomes). The dial belongs here, not in the FBM noise."
    );
}
