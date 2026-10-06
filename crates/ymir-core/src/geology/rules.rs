//! ADR Finding 152-B4 — the declarative RULES FILE (TOML), read by Ymir. One block per resource; a rule's conditions
//! all hold (AND) and give `chance`; a resource's favourability is the MAX of its rules; placers then travel downstream
//! along the rivers. The format is `docs/geology_format.md`'s.
//!
//! Unknown keys are refused (`deny_unknown_fields`): a typo in the author's file must fail loudly, not silently drop a
//! condition.

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use super::rocks::rock_id;

/// The rules file shipped with Ymir (the v0 proposal of Finding 151, with F152's `modulate`).
pub const DEFAULT_RULES_TOML: &str = include_str!("../../data/geology_rules_v0.toml");
/// The format tag a rules file must carry.
pub const RULES_FORMAT: &str = "ymir-geology-rules";

/// `[modulation.structural]`: `g = g_min + (1 − g_min)·min(d / d_ref, 1)`.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct StructuralModulation {
    pub d_ref: f32,
    pub g_min: f32,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Modulation {
    pub structural: Option<StructuralModulation>,
}

/// One rule: every `Some` condition must hold.
#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    pub rock: Option<Vec<String>>,
    pub belt: Option<String>,
    pub arc_within_km: Option<f32>,
    pub volcano_within_km: Option<f32>,
    pub volcano_active: Option<bool>,
    pub volcano_setting: Option<String>,
    pub wetland: Option<bool>,
    pub slope_min: Option<f32>,
    pub slope_max: Option<f32>,
    pub valley_floor: Option<bool>,
    pub lake: Option<bool>,
    pub endorheic_lake_within_km: Option<f32>,
    pub river_min_area_km2: Option<f32>,
    pub precip_min_mm: Option<f32>,
    pub precip_max_mm: Option<f32>,
    pub temp_min_c: Option<f32>,
    pub temp_max_c: Option<f32>,
    pub coast_within_km: Option<f32>,
    pub chance: u8,
    /// `"structural"`: the favourability is multiplied by the structural modulation g(d).
    pub modulate: Option<String>,
    pub source: Option<String>,
}

/// A placer: from the resource's own zones, downstream along the rivers.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Placer {
    pub source_min_chance: u8,
    pub river_min_area_km2: f32,
    pub start_chance: u8,
    pub decay_per_km: f32,
    pub max_slope: f32,
    pub source: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Resource {
    pub id: String,
    pub name_fr: String,
    pub color: String,
    #[serde(default)]
    pub rule: Vec<Rule>,
    pub placer: Option<Placer>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GeologyRules {
    pub format: String,
    pub version: String,
    pub notes: Option<String>,
    #[serde(default)]
    pub modulation: Modulation,
    pub resource: Vec<Resource>,
}

/// The SHA-256 of the rules file's bytes, in hex (the export's header).
#[must_use]
pub fn sha256_hex(bytes: &[u8]) -> String {
    let d = Sha256::digest(bytes);
    d.iter().map(|b| format!("{b:02x}")).collect()
}

/// `#RRGGBB` → RGB.
pub fn parse_color(s: &str) -> Result<[u8; 3], String> {
    let t = s.trim_start_matches('#');
    if t.len() != 6 {
        return Err(format!("colour '{s}': expected #RRGGBB"));
    }
    let p = |i: usize| u8::from_str_radix(&t[i..i + 2], 16).map_err(|e| format!("colour '{s}': {e}"));
    Ok([p(0)?, p(2)?, p(4)?])
}

/// Parse and validate a rules file.
pub fn parse_rules(text: &str) -> Result<GeologyRules, String> {
    let r: GeologyRules = toml::from_str(text).map_err(|e| format!("rules file: {e}"))?;
    if r.format != RULES_FORMAT {
        return Err(format!("rules file: format '{}' is not '{RULES_FORMAT}'", r.format));
    }
    let mut ids = std::collections::HashSet::new();
    for res in &r.resource {
        if !ids.insert(res.id.clone()) {
            return Err(format!("resource '{}' defined twice", res.id));
        }
        parse_color(&res.color).map_err(|e| format!("resource '{}': {e}", res.id))?;
        for (k, rule) in res.rule.iter().enumerate() {
            let at = format!("resource '{}', rule {}", res.id, k + 1);
            if rule.chance > 100 {
                return Err(format!("{at}: chance {} > 100", rule.chance));
            }
            for key in rule.rock.iter().flatten() {
                if rock_id(key).is_none() {
                    return Err(format!("{at}: unknown rock '{key}'"));
                }
            }
            if let Some(b) = &rule.belt
                && b != "fossil"
                && b != "active"
            {
                return Err(format!("{at}: belt '{b}' is not 'fossil' or 'active'"));
            }
            if let Some(s) = &rule.volcano_setting
                && !matches!(s.as_str(), "arc" | "hotspot" | "rift")
            {
                return Err(format!("{at}: volcano_setting '{s}' is not arc / hotspot / rift"));
            }
            if let Some(m) = &rule.modulate {
                if m != "structural" {
                    return Err(format!("{at}: modulate '{m}' is not 'structural'"));
                }
                if r.modulation.structural.is_none() {
                    return Err(format!("{at}: modulate = \"structural\" needs a [modulation.structural] block"));
                }
            }
        }
        if let Some(p) = &res.placer
            && (p.start_chance > 100 || p.source_min_chance > 100)
        {
            return Err(format!("resource '{}': placer chances must be ≤ 100", res.id));
        }
    }
    Ok(r)
}

/// A loaded rules file: the rules, the SHA-256 of the bytes they came from, where they came from, and the error that
/// made the loader fall back to the shipped rules (if any).
#[derive(Clone, Debug)]
pub struct LoadedRules {
    pub rules: GeologyRules,
    pub sha256: String,
    pub label: String,
    pub error: Option<String>,
}

/// Load the rules from `path` when it exists, else the shipped [`DEFAULT_RULES_TOML`]. A file that fails to parse falls
/// back to the shipped rules WITH its error, so the viz can say so instead of silently zoning with other rules.
#[must_use]
pub fn load_rules(path: Option<&std::path::Path>) -> LoadedRules {
    let shipped = || LoadedRules {
        rules: parse_rules(DEFAULT_RULES_TOML).expect("the shipped rules parse (tested)"),
        sha256: sha256_hex(DEFAULT_RULES_TOML.as_bytes()),
        label: "v0 intégré".to_string(),
        error: None,
    };
    let Some(p) = path.filter(|p| p.exists()) else { return shipped() };
    match std::fs::read(p).map_err(|e| e.to_string()).and_then(|b| {
        let text = String::from_utf8(b.clone()).map_err(|e| e.to_string())?;
        parse_rules(&text).map(|r| (r, sha256_hex(&b)))
    }) {
        Ok((rules, sha256)) => LoadedRules { rules, sha256, label: p.display().to_string(), error: None },
        Err(e) => LoadedRules { error: Some(format!("{}: {e}", p.display())), ..shipped() },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The shipped file parses, has the author's 12 resources, and the vein rules carry `modulate`.
    #[test]
    fn the_default_rules_parse() {
        let r = parse_rules(DEFAULT_RULES_TOML).expect("the shipped rules parse");
        let ids: Vec<&str> = r.resource.iter().map(|x| x.id.as_str()).collect();
        assert_eq!(ids, ["iron", "gold", "silver", "copper", "tin", "stone", "clay", "sand_gravel", "salt", "peat", "sulphur_obsidian", "gems"]);
        let modulated: usize = r.resource.iter().flat_map(|x| x.rule.iter()).filter(|x| x.modulate.is_some()).count();
        assert!(modulated >= 6, "the vein rules are modulated ({modulated})");
        assert!(r.modulation.structural.is_some());
    }

    /// The negative controls: an unknown key, an unknown rock and a missing modulation block are refused.
    #[test]
    fn a_bad_rules_file_is_refused() {
        let base = "format = \"ymir-geology-rules\"\nversion = \"0.1.0\"\n[[resource]]\nid = \"x\"\nname_fr = \"X\"\ncolor = \"#112233\"\n[[resource.rule]]\nchance = 10\n";
        assert!(parse_rules(base).is_ok());
        assert!(parse_rules(&format!("{base}rokc = [\"craton\"]\n")).is_err(), "a typo is refused");
        assert!(parse_rules(&format!("{base}rock = [\"granite\"]\n")).is_err(), "an unknown rock is refused");
        assert!(parse_rules(&format!("{base}modulate = \"structural\"\n")).is_err(), "modulate needs its block");
    }
}
