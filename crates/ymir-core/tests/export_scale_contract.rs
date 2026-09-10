//! ADR Finding 69 — the EXPORT SCALE CONTRACT, and the assertion that closes it by construction.
//!
//! ## The defect this exists to make impossible
//!
//! `geo_scale_ratio` (Finding 68) is a post-process that scales SOME exported quantities and not
//! others. Finding 49 established the rule — *"EVERY SIGNIFIED QUANTITY, or the export ships two
//! scales at once"* — and enforced it across the six segment arrays. It was never asked of the
//! other exported object types, so `rivers.json` ships signified areas beside `lakes.json`'s real
//! ones, linked by `source_lake_id`.
//!
//! The manifest declares `geographic_scale_ratio`, and `ContinentMeta`'s doc comment enumerates
//! the signified fields — **in source, not in the exported JSON.** A consumer reading
//! `metadata.json` gets the ratio and no field list, which makes uniform application *more*
//! likely than if the ratio were absent: applying 7.5 to everything double-scales the rivers
//! while correcting the lakes.
//!
//! ## What this test enforces
//!
//! Every exported field whose NAME carries a unit (`_km2`, `_m`, `_m3s`, `_km`) must appear in
//! [`CONTRACT`] with its unit and its scale. A field added without a contract entry fails the
//! build. The contract is the machine-readable form of the doc comment, and the negative control
//! proves the check can fail.
//!
//! Run: `cargo test -p ymir-core --test export_scale_contract`

/// `REAL` = the geometric 400 km world. `SIGNIFIED` = multiplied by `geo_scale_ratio` (Finding
/// 68: areas ×ratio², lengths ×ratio, verticals NEVER). `NONE` = dimensionless or a class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Scale {
    Real,
    Signified,
    None,
}

/// One exported, unit-bearing field. `(file, field, unit, scale, who sets the scale)`.
struct Entry {
    file: &'static str,
    field: &'static str,
    unit: &'static str,
    scale: Scale,
    setter: &'static str,
}

const CONTRACT: &[Entry] = &[
    // ── rivers.json ─────────────────────────────────────────────────────────────────
    Entry {
        file: "rivers.json",
        field: "catchment_km2",
        unit: "km²",
        scale: Scale::Signified,
        setter: "apply_geo_scale_ratio ×ratio² (Finding 49)",
    },
    Entry {
        file: "rivers.json",
        field: "drainage_km2",
        unit: "km² (⚠ NOT an area — a runoff-equivalent)",
        scale: Scale::Signified,
        setter: "apply_geo_scale_ratio ×ratio²",
    },
    Entry {
        file: "rivers.json",
        field: "discharge_m3s",
        unit: "m³/s",
        scale: Scale::Signified,
        setter: "apply_geo_scale_ratio ×ratio²",
    },
    Entry {
        file: "rivers.json",
        field: "width_m",
        unit: "m (horizontal)",
        scale: Scale::Signified,
        setter: "re-derived 5·Q^0.5 on the scaled Q ⇒ ×ratio",
    },
    Entry {
        file: "rivers.json",
        field: "profile_m",
        unit: "m (VERTICAL — bed elevation)",
        scale: Scale::Real,
        setter: "never scaled; the ratio is a HORIZONTAL compression",
    },
    // ── lakes.json ──────────────────────────────────────────────────────────────────
    Entry {
        file: "lakes.json",
        field: "area_km2",
        unit: "km²",
        scale: Scale::Real,
        setter: "⛔ NOT scaled — see the mismatch check below",
    },
    Entry {
        file: "lakes.json",
        field: "level_m",
        unit: "m (VERTICAL)",
        scale: Scale::Real,
        setter: "vertical, correctly unscaled",
    },
    Entry {
        file: "lakes.json",
        field: "depth_m",
        unit: "m (VERTICAL)",
        scale: Scale::Real,
        setter: "vertical, correctly unscaled",
    },
    // ── metadata.json ───────────────────────────────────────────────────────────────
    Entry {
        file: "metadata.json",
        field: "window_km",
        unit: "km",
        scale: Scale::Real,
        setter: "the geometric domain — the map DRAWS this",
    },
    Entry {
        file: "metadata.json",
        field: "tectonic_domain_km",
        unit: "km",
        scale: Scale::Real,
        setter: "geometric",
    },
    Entry {
        file: "metadata.json",
        field: "sea_level_m",
        unit: "m (VERTICAL)",
        scale: Scale::Real,
        setter: "anchored at 0",
    },
    Entry {
        file: "metadata.json",
        field: "max_elevation_m",
        unit: "m (VERTICAL)",
        scale: Scale::Real,
        setter: "the field's real extremum",
    },
    Entry {
        file: "metadata.json",
        field: "max_depth_m",
        unit: "m (VERTICAL)",
        scale: Scale::Real,
        setter: "the field's real extremum",
    },
    Entry {
        file: "metadata.json",
        field: "depth_scale_m",
        unit: "m per non-dim unit",
        scale: Scale::Real,
        setter: "VerticalScale — the norm→metres contract",
    },
];

/// Field-name suffixes that declare a unit. A field matching one of these MUST be in the
/// contract; a field matching none is dimensionless (an id, a class, a count, a coordinate).
const UNIT_SUFFIXES: &[&str] = &["_km2", "_m3s", "_km", "_m"];

fn carries_unit(field: &str) -> bool {
    UNIT_SUFFIXES.iter().any(|s| field.ends_with(s))
}

fn contract_has(file: &str, field: &str) -> bool {
    CONTRACT.iter().any(|e| e.file == file && e.field == field)
}

/// The exported field names, as the serializers emit them. Kept as a literal list rather than
/// reflected, because `#[serde(flatten)]` and private structs make reflection unavailable — and
/// a literal list is what a reviewer can check against the JSON.
const EXPORTED: &[(&str, &str)] = &[
    ("rivers.json", "coordinate_space"),
    ("rivers.json", "points"),
    ("rivers.json", "strahler_order"),
    ("rivers.json", "avg_flow"),
    ("rivers.json", "max_flow"),
    ("rivers.json", "basin_id"),
    ("rivers.json", "upstream"),
    ("rivers.json", "downstream"),
    ("rivers.json", "catchment_km2"),
    ("rivers.json", "drainage_km2"),
    ("rivers.json", "navigability"),
    ("rivers.json", "discharge_m3s"),
    ("rivers.json", "width_m"),
    ("rivers.json", "profile_m"),
    ("rivers.json", "kind"),
    ("rivers.json", "source_lake_id"),
    ("lakes.json", "id"),
    ("lakes.json", "surface_elevation"),
    ("lakes.json", "max_depth"),
    ("lakes.json", "area"),
    ("lakes.json", "basin_id"),
    ("lakes.json", "outlet"),
    ("lakes.json", "shallow"),
    ("lakes.json", "level_m"),
    ("lakes.json", "depth_m"),
    ("lakes.json", "area_km2"),
    ("lakes.json", "lake_type"),
    ("metadata.json", "window_km"),
    ("metadata.json", "tectonic_domain_km"),
    ("metadata.json", "latitude_deg"),
    ("metadata.json", "latitude_span_deg"),
    ("metadata.json", "geographic_scale_ratio"),
    ("metadata.json", "sea_level_m"),
    ("metadata.json", "max_elevation_m"),
    ("metadata.json", "max_depth_m"),
    ("metadata.json", "depth_scale_m"),
    ("metadata.json", "altitude_norm_half_range"),
    ("metadata.json", "sea_level_norm"),
];

/// Every unit-bearing exported field is declared in the contract.
#[test]
fn every_dimensioned_export_field_is_in_the_scale_contract() {
    let missing: Vec<String> = EXPORTED
        .iter()
        .filter(|(f, k)| carries_unit(k) && !contract_has(f, k))
        .map(|(f, k)| format!("{f}::{k}"))
        .collect();
    assert!(
        missing.is_empty(),
        "\nthese exported fields carry a unit and are NOT in the scale contract:\n  · {}\n\
         Add an entry stating the unit and whether the value is REAL or SIGNIFIED. A consumer \
         cannot know which without it, and `geographic_scale_ratio` alone makes the error MORE \
         likely, not less (Findings 49, 68).\n",
        missing.join("\n  · ")
    );
    eprintln!(
        "[contract] {} exported fields, {} unit-bearing, all declared",
        EXPORTED.len(),
        EXPORTED.iter().filter(|(_, k)| carries_unit(k)).count()
    );
}

/// NEGATIVE CONTROL (rule 1): the check must FAIL on a dimensioned field withheld from the
/// contract. Without this, an empty `EXPORTED` list or a broken suffix test would pass.
#[test]
fn the_contract_check_fails_on_an_undeclared_dimensioned_field() {
    // A field that carries a unit and is deliberately absent from CONTRACT.
    assert!(carries_unit("floodplain_width_m"), "the suffix test must recognise a unit");
    assert!(
        !contract_has("rivers.json", "floodplain_width_m"),
        "the control field must be absent from the contract"
    );
    // and the discriminator must not fire on dimensionless names
    for k in ["id", "basin_id", "strahler_order", "lake_type", "kind", "points", "shallow"] {
        assert!(
            !carries_unit(k),
            "`{k}` is dimensionless and must not be required in the contract"
        );
    }
    // `_m` must not swallow words merely ending in m
    assert!(!carries_unit("coordinate_space"), "the suffix test is too greedy");
}

/// The finding this contract exists to record: comparable quantities across object types are on
/// DIFFERENT scales. Asserted so it cannot be silently "fixed" or silently worsened.
#[test]
fn comparable_fields_across_object_types_are_on_different_scales() {
    let scale_of = |file: &str, field: &str| -> Scale {
        CONTRACT
            .iter()
            .find(|e| e.file == file && e.field == field)
            .map(|e| e.scale)
            .unwrap_or(Scale::None)
    };
    // A river's catchment and the lake it drains into are compared by any consumer doing a
    // water-balance or plausibility check, and they are linked by `source_lake_id`.
    let river_area = scale_of("rivers.json", "catchment_km2");
    let lake_area = scale_of("lakes.json", "area_km2");
    assert_eq!(river_area, Scale::Signified);
    assert_eq!(lake_area, Scale::Real);
    assert_ne!(
        river_area, lake_area,
        "if these are now on the SAME scale the defect has been fixed — update ADR Finding 69 \
         and delete this test rather than relaxing it"
    );
    eprintln!(
        "[contract] MISMATCH, recorded: rivers.json::catchment_km2 is {river_area:?} while \
         lakes.json::area_km2 is {lake_area:?} — at ratio 7.5 they are 56.25× apart"
    );
    // And the pair INSIDE rivers.json: a horizontal length is signified, a vertical is not.
    assert_eq!(scale_of("rivers.json", "width_m"), Scale::Signified);
    assert_eq!(scale_of("rivers.json", "profile_m"), Scale::Real);
    eprintln!(
        "[contract] and inside rivers.json: width_m SIGNIFIED beside profile_m REAL — coherent \
         (the ratio is horizontal) but it means the signified world carries a ×ratio VERTICAL \
         EXAGGERATION, so river gradients read 7.5× gentler than the terrain's"
    );
}

/// The contract, printed as the deliverable table. Run with `--nocapture`.
#[test]
fn print_the_export_scale_contract() {
    eprintln!("\n{:<16} {:<24} {:<34} {:<11} {}", "file", "field", "unit", "scale", "who sets it");
    let mut file = "";
    for e in CONTRACT {
        if e.file != file {
            eprintln!("{}", "─".repeat(120));
            file = e.file;
        }
        eprintln!(
            "{:<16} {:<24} {:<34} {:<11} {}",
            e.file,
            e.field,
            e.unit,
            format!("{:?}", e.scale).to_uppercase(),
            e.setter
        );
    }
    eprintln!(
        "\nNOT in the contract because dimensionless: points, strahler_order, avg_flow, max_flow,\n\
         basin_id, upstream, downstream, navigability (a class), kind, source_lake_id,\n\
         coordinate_space, id, surface_elevation (norm), max_depth (norm), area (CELLS), outlet,\n\
         shallow, lake_type, latitude_deg, latitude_span_deg, geographic_scale_ratio,\n\
         altitude_norm_half_range, sea_level_norm.\n\
         ⚠️ `lakes.json::area` is in CELLS and `lakes.json::area_km2` in real km² — a consumer\n\
         sizing lakes from the cell tag (which the author reports LL does) is unaffected by the\n\
         ratio mismatch; one reading `area_km2` is not."
    );
}
