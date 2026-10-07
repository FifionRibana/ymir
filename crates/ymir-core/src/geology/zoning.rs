//! ADR Finding 152-B4/B5/B7 — the ZONING: the per-world CONTEXT on the ¼ grid (built once, kept with the world) and
//! the evaluation of a rules file on it into one u8 FAVOURABILITY grid per resource (0–100, relative, not a probability:
//! « Favorabilité dans Ymir, rareté dans Living Landz »). Reloading the rules reruns [`zone`] only.

use std::collections::HashSet;

use rayon::prelude::*;

use super::rocks::{self, NONE, RockClass};
use super::rules::{GeologyRules, Rule, parse_color};
use crate::grid::GridF32;
use crate::tectonics_c1::closures::lithology::upscale_k_to_hd;
use crate::tectonics_c1::closures::volcanism::{Edifice, VolcanoSetting};
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

/// The chance grid is the HD grid divided by this (DECISION, F151: one chance cell ≈ 9 hexes).
pub const CHANCE_FACTOR: usize = 4;
/// The placers' river network: HD cells with at least this REAL drained area (km²); a rule's smaller area is clamped.
pub const PLACER_NETWORK_KM2: f32 = 1.0;

/// An edifice on the ¼ grid.
#[derive(Clone, Copy, Debug)]
pub struct EdificeAt {
    pub x4: f32,
    pub y4: f32,
    pub rb_km: f32,
    pub active: bool,
    pub setting: VolcanoSetting,
}

/// One HD river cell of the placers' network, upstream before downstream (increasing accumulation).
#[derive(Clone, Copy, Debug)]
pub struct RiverCell {
    pub q4: u32,
    /// Its receiver's position in the list (`u32::MAX`: none).
    pub recv: u32,
    pub step_km: f32,
    pub slope: f32,
    pub area_km2: f32,
}

/// The per-world context on the ¼ grid.
#[derive(Clone, Debug)]
pub struct ZoningContext {
    pub w4: usize,
    pub h4: usize,
    pub cell4_km: f32,
    pub rock: Vec<u8>,
    pub slope: Vec<f32>,
    pub area_km2: Vec<f32>,
    pub lake: Vec<bool>,
    pub wetland: Vec<bool>,
    pub valley_floor: Vec<bool>,
    pub temp_c: Vec<f32>,
    pub precip_mm: Vec<f32>,
    pub coast_km: Vec<f32>,
    pub endo_km: Vec<f32>,
    pub arc_km: Vec<f32>,
    pub belt_fossil: Vec<bool>,
    pub belt_active: Vec<bool>,
    /// The structural density (sutures + past collisions + C-3b) and C-3b's alone (G-veins' control).
    pub density: Vec<f32>,
    pub density_c3b: Vec<f32>,
    pub edifices: Vec<EdificeAt>,
    pub rivers: Vec<RiverCell>,
}

/// The HD and coarse inputs of the context.
pub struct ContextInputs<'a> {
    pub w: usize,
    pub h: usize,
    pub cell_km: f32,
    pub rocks: &'a [u8],
    pub slope: &'a [f32],
    pub acc: &'a [f32],
    pub dir: &'a [u8],
    pub lake_map: &'a [u32],
    pub endorheic: &'a HashSet<u32>,
    pub wetland: &'a [u8],
    pub temp_c: &'a [f32],
    pub precip_mm: &'a [f32],
    pub sea: &'a [bool],
    pub nx: usize,
    pub ny: usize,
    pub fossil: &'a [bool],
    pub active: &'a [bool],
    pub upper_plate: &'a [bool],
    pub density: &'a GridF32,
    pub density_c3b: &'a GridF32,
    pub sample_origin: [f64; 2],
    pub sample_size: f64,
    pub edifices: &'a [Edifice],
}

/// Two-pass chamfer distance (cells) to the `true` cells.
#[must_use]
pub fn chamfer(m: &[bool], w: usize, h: usize) -> Vec<f32> {
    let mut d: Vec<f32> = m.iter().map(|&b| if b { 0.0 } else { f32::INFINITY }).collect();
    let r2 = std::f32::consts::SQRT_2;
    for y in 0..h {
        for x in 0..w {
            let c = y * w + x;
            let mut v = d[c];
            if x > 0 {
                v = v.min(d[c - 1] + 1.0);
            }
            if y > 0 {
                v = v.min(d[c - w] + 1.0);
                if x > 0 {
                    v = v.min(d[c - w - 1] + r2);
                }
                if x + 1 < w {
                    v = v.min(d[c - w + 1] + r2);
                }
            }
            d[c] = v;
        }
    }
    for y in (0..h).rev() {
        for x in (0..w).rev() {
            let c = y * w + x;
            let mut v = d[c];
            if x + 1 < w {
                v = v.min(d[c + 1] + 1.0);
            }
            if y + 1 < h {
                v = v.min(d[c + w] + 1.0);
                if x + 1 < w {
                    v = v.min(d[c + w + 1] + r2);
                }
                if x > 0 {
                    v = v.min(d[c + w - 1] + r2);
                }
            }
            d[c] = v;
        }
    }
    d
}

/// Build the context (once per world).
#[must_use]
pub fn build_context(inp: &ContextInputs) -> ZoningContext {
    let (w, h, f) = (inp.w, inp.h, CHANCE_FACTOR);
    let (w4, h4) = (w / f, h / f);
    let n4 = w4 * h4;
    let cell4_km = inp.cell_km * f as f32;
    let cell_km2 = inp.cell_km * inp.cell_km;
    let flood = rocks::floodplain(inp.slope, inp.acc, inp.cell_km);
    // the block statistics
    let blocks: Vec<(u8, f32, f32, bool, bool, bool, f32, f32, bool, bool)> = (0..n4)
        .into_par_iter()
        .map(|q| {
            let (bx, by) = (q % w4, q / w4);
            let mut votes = [0u32; 9];
            let (mut s, mut amax, mut t, mut p) = (0f32, 0f32, 0f32, 0f32);
            let (mut lake, mut wet, mut vf, mut endo, mut sea) = (false, false, false, false, 0u32);
            for dy in 0..f {
                for dx in 0..f {
                    let c = (by * f + dy) * w + bx * f + dx;
                    votes[inp.rocks[c] as usize] += 1;
                    s += inp.slope[c];
                    amax = amax.max(inp.acc[c] * cell_km2);
                    t += inp.temp_c[c];
                    p += inp.precip_mm[c];
                    lake |= inp.lake_map[c] != 0;
                    endo |= inp.lake_map[c] != 0 && inp.endorheic.contains(&inp.lake_map[c]);
                    wet |= inp.wetland[c] != 0;
                    vf |= flood[c];
                    sea += inp.sea[c] as u32;
                }
            }
            let k = (f * f) as f32;
            // the majority of the land classes (the water class only when nothing else)
            let rock = (1..9).max_by_key(|&i| (votes[i], std::cmp::Reverse(i))).filter(|&i| votes[i] > 0).map_or(NONE, |i| i as u8);
            (rock, s / k, amax, lake, wet, vf, t / k, p / k, endo, sea * 2 > (f * f) as u32)
        })
        .collect();
    let rock: Vec<u8> = blocks.iter().map(|b| b.0).collect();
    let sea4: Vec<bool> = blocks.iter().map(|b| b.9).collect();
    let endo4: Vec<bool> = blocks.iter().map(|b| b.8).collect();
    // the coarse fields on the ¼ grid
    let fossil_s = rocks::smooth_coarse(inp.fossil, inp.nx, inp.ny, inp.sample_origin, inp.sample_size);
    let active_s = rocks::smooth_coarse(inp.active, inp.nx, inp.ny, inp.sample_origin, inp.sample_size);
    let belt_fossil: Vec<bool> = (0..n4).map(|q| rocks::at_smooth(&fossil_s, q % w4, q / w4, w4, h4) > 0.5).collect();
    let belt_active: Vec<bool> = (0..n4).map(|q| rocks::at_smooth(&active_s, q % w4, q / w4, w4, h4) > 0.5).collect();
    // the upper-plate seeds (nearest coarse cell) → the arc distance
    let (ox, oy) = (inp.sample_origin[0] * inp.nx as f64, inp.sample_origin[1] * inp.ny as f64);
    let upper4: Vec<bool> = (0..n4)
        .map(|q| {
            let sx = ox + (q % w4) as f64 * inp.sample_size * inp.nx as f64 / w4 as f64;
            let sy = oy + (q / w4) as f64 * inp.sample_size * inp.ny as f64 / h4 as f64;
            let ci = (sx.round() as i64).rem_euclid(inp.nx as i64) as usize;
            let cj = (sy.round() as i64).rem_euclid(inp.ny as i64) as usize;
            inp.upper_plate[cj * inp.nx + ci]
        })
        .collect();
    let km = |d: Vec<f32>| -> Vec<f32> { d.into_iter().map(|v| v * cell4_km).collect() };
    let coast_km = km(chamfer(&sea4, w4, h4));
    let endo_km = km(chamfer(&endo4, w4, h4));
    let arc_km = km(chamfer(&upper4, w4, h4));
    let density = upscale_k_to_hd(inp.density, w4, h4, inp.sample_origin, inp.sample_size);
    let density_c3b = upscale_k_to_hd(inp.density_c3b, w4, h4, inp.sample_origin, inp.sample_size);
    // the edifices on the ¼ grid
    let (so, ss) = ([inp.sample_origin[0] as f32, inp.sample_origin[1] as f32], inp.sample_size as f32);
    let edifices: Vec<EdificeAt> = inp
        .edifices
        .iter()
        .filter_map(|e| {
            let fx = (e.center_uv.0 - so[0]).rem_euclid(1.0) / ss;
            let fy = (e.center_uv.1 - so[1]).rem_euclid(1.0) / ss;
            (fx < 1.0 && fy < 1.0).then(|| EdificeAt { x4: fx * w4 as f32, y4: fy * h4 as f32, rb_km: e.basal_diameter_km * 0.5, active: e.active, setting: e.setting })
        })
        .collect();
    // the placers' river network, upstream before downstream
    let a_min = PLACER_NETWORK_KM2 / cell_km2;
    let mut cells: Vec<usize> = (0..w * h).filter(|&c| inp.acc[c] >= a_min && !inp.sea[c]).collect();
    cells.sort_by(|&a, &b| inp.acc[a].total_cmp(&inp.acc[b]).then(a.cmp(&b)));
    let mut pos = vec![u32::MAX; w * h];
    for (i, &c) in cells.iter().enumerate() {
        pos[c] = i as u32;
    }
    let rivers: Vec<RiverCell> = cells
        .iter()
        .map(|&c| {
            let d = inp.dir[c];
            let (recv, step) = if d == DIR_NONE {
                (u32::MAX, 0.0)
            } else {
                let x = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
                let y = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
                let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
                (pos[y * w + x], inp.cell_km * if diag { std::f32::consts::SQRT_2 } else { 1.0 })
            };
            let q4 = ((c / w) / f).min(h4 - 1) * w4 + ((c % w) / f).min(w4 - 1);
            RiverCell { q4: q4 as u32, recv, step_km: step, slope: inp.slope[c], area_km2: inp.acc[c] * cell_km2 }
        })
        .collect();
    drop(pos);
    ZoningContext {
        w4,
        h4,
        cell4_km,
        rock,
        slope: blocks.iter().map(|b| b.1).collect(),
        area_km2: blocks.iter().map(|b| b.2).collect(),
        lake: blocks.iter().map(|b| b.3).collect(),
        wetland: blocks.iter().map(|b| b.4).collect(),
        valley_floor: blocks.iter().map(|b| b.5).collect(),
        temp_c: blocks.iter().map(|b| b.6).collect(),
        precip_mm: blocks.iter().map(|b| b.7).collect(),
        coast_km,
        endo_km,
        arc_km,
        belt_fossil,
        belt_active,
        density,
        density_c3b,
        edifices,
        rivers,
    }
}

/// Which density the structural modulation reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DensitySource {
    /// Sutures + past collisions + C-3b (production).
    Structural,
    /// C-3b alone (G-veins' control).
    C3bOnly,
}

/// One resource's favourability grid.
#[derive(Clone, Debug)]
pub struct ResourceGrid {
    pub id: String,
    pub name_fr: String,
    pub color: [u8; 3],
    /// u8 0–100 on the ¼ grid.
    pub fav: Vec<u8>,
    /// The HD river cells a placer set.
    pub placer_cells: usize,
}

/// A zoning: the favourabilities of every resource of a rules file.
#[derive(Clone, Debug)]
pub struct Zoning {
    pub w4: usize,
    pub h4: usize,
    pub resources: Vec<ResourceGrid>,
}

fn rule_holds(r: &Rule, ctx: &ZoningContext, q: usize, rock_ids: &[u8]) -> bool {
    if !rock_ids.is_empty() && !rock_ids.contains(&ctx.rock[q]) {
        return false;
    }
    if let Some(b) = &r.belt
        && !(if b == "fossil" { ctx.belt_fossil[q] } else { ctx.belt_active[q] })
    {
        return false;
    }
    if let Some(d) = r.arc_within_km
        && ctx.arc_km[q] > d
    {
        return false;
    }
    if r.volcano_within_km.is_some() || r.volcano_active.is_some() || r.volcano_setting.is_some() {
        let d = r.volcano_within_km.unwrap_or(0.0);
        let (x, y) = ((q % ctx.w4) as f32 + 0.5, (q / ctx.w4) as f32 + 0.5);
        let ok = ctx.edifices.iter().any(|e| {
            if r.volcano_active.is_some_and(|a| a != e.active) {
                return false;
            }
            if let Some(s) = &r.volcano_setting {
                let want = match s.as_str() {
                    "arc" => VolcanoSetting::Arc,
                    "hotspot" => VolcanoSetting::Hotspot,
                    _ => VolcanoSetting::Rift,
                };
                if e.setting != want {
                    return false;
                }
            }
            let dist_km = ((x - e.x4).powi(2) + (y - e.y4).powi(2)).sqrt() * ctx.cell4_km;
            (dist_km - e.rb_km).max(0.0) <= d
        });
        if !ok {
            return false;
        }
    }
    let flag = |want: Option<bool>, have: bool| want.is_none_or(|w| w == have);
    if !flag(r.wetland, ctx.wetland[q]) || !flag(r.lake, ctx.lake[q]) || !flag(r.valley_floor, ctx.valley_floor[q]) {
        return false;
    }
    let s = ctx.slope[q];
    if r.slope_min.is_some_and(|v| s < v) || r.slope_max.is_some_and(|v| s > v) {
        return false;
    }
    if r.endorheic_lake_within_km.is_some_and(|d| ctx.endo_km[q] > d) || r.coast_within_km.is_some_and(|d| ctx.coast_km[q] > d) {
        return false;
    }
    if r.river_min_area_km2.is_some_and(|a| ctx.area_km2[q] < a) {
        return false;
    }
    let (p, t) = (ctx.precip_mm[q], ctx.temp_c[q]);
    if r.precip_min_mm.is_some_and(|v| p < v) || r.precip_max_mm.is_some_and(|v| p > v) || r.temp_min_c.is_some_and(|v| t < v) || r.temp_max_c.is_some_and(|v| t > v) {
        return false;
    }
    true
}

/// The structural modulation g(d) (`[modulation.structural]`).
#[inline]
#[must_use]
pub fn modulation(d: f32, d_ref: f32, g_min: f32) -> f32 {
    g_min + (1.0 - g_min) * (d / d_ref.max(1e-6)).min(1.0)
}

/// The best rule of a resource at a ¼ cell: its value (before rounding), its index, and the structural modulation g it
/// carried (if modulated). The FIRST rule of maximal value wins (the file's order), as in [`zone`].
fn best_rule(res: &super::rules::Resource, rock_sets: &[Vec<u8>], ctx: &ZoningContext, q: usize, dens: &[f32], sm: Option<super::rules::StructuralModulation>) -> (f32, Option<(usize, Option<f32>)>) {
    let mut best = 0f32;
    let mut which = None;
    for (k, (r, rs)) in res.rule.iter().zip(rock_sets).enumerate() {
        if r.chance as f32 <= best || !rule_holds(r, ctx, q, rs) {
            continue;
        }
        let mut v = r.chance as f32;
        let mut g = None;
        if r.modulate.is_some()
            && let Some(m) = sm
        {
            let gg = modulation(dens[q], m.d_ref, m.g_min);
            v *= gg;
            g = Some(gg);
        }
        if v > best {
            best = v;
            which = Some((k, g));
        }
    }
    (best, which)
}

fn rock_sets_of(res: &super::rules::Resource) -> Vec<Vec<u8>> {
    res.rule.iter().map(|r| r.rock.iter().flatten().filter_map(|k| rocks::rock_id(k)).collect()).collect()
}

/// Where a cell's favourability came from.
#[derive(Clone, Debug, PartialEq)]
pub enum FavOrigin {
    /// A rule of the file (its index, its `source`, and the structural modulation g it carried).
    Rule { index: usize, source: Option<String>, modulation: Option<f32> },
    /// A placer (downstream along the rivers).
    Placer { source: Option<String> },
}

/// One resource at one cell, as the inspector shows it.
#[derive(Clone, Debug)]
pub struct CellResource {
    pub id: String,
    pub name_fr: String,
    pub color: [u8; 3],
    /// Exactly the zoning grid's (and so the export's) value.
    pub fav: u8,
    pub origin: FavOrigin,
}

/// ADR Finding 153-I — the resources of non-zero favourability at the ¼ cell `q`, strongest first. The favourability
/// is READ from the zoning grid (so it is the export's); only the attribution is recomputed, with the same
/// evaluation as [`zone`]: the rule whose rounded value equals the grid's, else the placer.
#[must_use]
pub fn explain(ctx: &ZoningContext, rules: &GeologyRules, z: &Zoning, q: usize, density: DensitySource) -> Vec<CellResource> {
    let dens = match density {
        DensitySource::Structural => &ctx.density,
        DensitySource::C3bOnly => &ctx.density_c3b,
    };
    let sm = rules.modulation.structural;
    let mut out: Vec<CellResource> = rules
        .resource
        .iter()
        .zip(&z.resources)
        .filter_map(|(res, grid)| {
            let fav = *grid.fav.get(q)?;
            if fav == 0 {
                return None;
            }
            let (best, which) = best_rule(res, &rock_sets_of(res), ctx, q, dens, sm);
            let rule_val = best.round().clamp(0.0, 100.0) as u8;
            let origin = match which {
                Some((k, g)) if rule_val == fav => FavOrigin::Rule { index: k, source: res.rule[k].source.clone(), modulation: g },
                _ => FavOrigin::Placer { source: res.placer.as_ref().and_then(|p| p.source.clone()) },
            };
            Some(CellResource { id: res.id.clone(), name_fr: res.name_fr.clone(), color: grid.color, fav, origin })
        })
        .collect();
    out.sort_by(|a, b| b.fav.cmp(&a.fav));
    out
}

/// Evaluate a rules file on a context.
#[must_use]
pub fn zone(ctx: &ZoningContext, rules: &GeologyRules, density: DensitySource) -> Zoning {
    let n4 = ctx.w4 * ctx.h4;
    let dens = match density {
        DensitySource::Structural => &ctx.density,
        DensitySource::C3bOnly => &ctx.density_c3b,
    };
    let sm = rules.modulation.structural;
    let resources = rules
        .resource
        .iter()
        .map(|res| {
            let rock_sets = rock_sets_of(res);
            let mut fav: Vec<u8> = (0..n4)
                .into_par_iter()
                .map(|q| {
                    if ctx.rock[q] == NONE && !ctx.lake[q] {
                        return 0u8;
                    }
                    best_rule(res, &rock_sets, ctx, q, dens, sm).0.round().clamp(0.0, 100.0) as u8
                })
                .collect();
            // the placers: downstream along the rivers from the resource's own zones
            let mut placer_cells = 0usize;
            if let Some(p) = &res.placer {
                let src = fav.clone();
                let mut val = vec![f32::NEG_INFINITY; ctx.rivers.len()];
                let min_area = p.river_min_area_km2.max(PLACER_NETWORK_KM2);
                for (i, rc) in ctx.rivers.iter().enumerate() {
                    if rc.area_km2 >= min_area && src[rc.q4 as usize] >= p.source_min_chance {
                        val[i] = p.start_chance as f32;
                    }
                }
                for i in 0..ctx.rivers.len() {
                    let rc = ctx.rivers[i];
                    if val[i] > 0.0 && rc.recv != u32::MAX {
                        let carried = val[i] - p.decay_per_km * rc.step_km;
                        let r = rc.recv as usize;
                        if carried > val[r] {
                            val[r] = carried;
                        }
                    }
                }
                for (i, rc) in ctx.rivers.iter().enumerate() {
                    if val[i] > 0.0 && rc.area_km2 >= min_area && rc.slope < p.max_slope {
                        placer_cells += 1;
                        let q = rc.q4 as usize;
                        fav[q] = fav[q].max(val[i].round().clamp(0.0, 100.0) as u8);
                    }
                }
            }
            ResourceGrid { id: res.id.clone(), name_fr: res.name_fr.clone(), color: parse_color(&res.color).unwrap_or([255, 0, 255]), fav, placer_cells }
        })
        .collect();
    Zoning { w4: ctx.w4, h4: ctx.h4, resources }
}

/// The rock classes' legend (re-exported for the viz and the export).
#[must_use]
pub fn rock_legend() -> &'static [RockClass] {
    &rocks::ROCK_CLASSES
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geology::rules::parse_rules;

    /// A synthetic context: a 16×16 ¼ grid, craton on the left half, belt on the right.
    fn ctx() -> ZoningContext {
        let (w4, h4) = (16usize, 16usize);
        let n = w4 * h4;
        let rock: Vec<u8> = (0..n).map(|q| if q % w4 < 8 { rocks::CRATON } else { rocks::BELT }).collect();
        ZoningContext {
            w4,
            h4,
            cell4_km: 0.2,
            rock,
            slope: vec![0.1; n],
            area_km2: vec![0.5; n],
            lake: vec![false; n],
            wetland: vec![false; n],
            valley_floor: vec![false; n],
            temp_c: vec![10.0; n],
            precip_mm: vec![800.0; n],
            coast_km: vec![50.0; n],
            endo_km: vec![50.0; n],
            arc_km: vec![50.0; n],
            belt_fossil: (0..n).map(|q| q % w4 >= 8).collect(),
            belt_active: vec![false; n],
            density: (0..n).map(|q| (q % w4) as f32 / 32.0).collect(),
            density_c3b: vec![0.0; n],
            edifices: Vec::new(),
            rivers: Vec::new(),
        }
    }

    const RULES: &str = "format = \"ymir-geology-rules\"\nversion = \"0.1.0\"\n[modulation.structural]\nd_ref = 0.5\ng_min = 0.2\n\
        [[resource]]\nid = \"iron\"\nname_fr = \"Fer\"\ncolor = \"#B5562F\"\n[[resource.rule]]\nrock = [\"craton\"]\nchance = 40\n\
        [[resource]]\nid = \"gold\"\nname_fr = \"Or\"\ncolor = \"#E6C229\"\n[[resource.rule]]\nrock = [\"belt\"]\nchance = 35\nmodulate = \"structural\"\n";

    /// ADR Finding 152 G-rules — changing ONE rule's chance changes that resource's grid only. The negative control:
    /// the unchanged file gives bit-identical grids.
    #[test]
    fn changing_one_rule_changes_only_its_resource() {
        let c = ctx();
        let a = zone(&c, &parse_rules(RULES).unwrap(), DensitySource::Structural);
        let same = zone(&c, &parse_rules(RULES).unwrap(), DensitySource::Structural);
        assert_eq!(a.resources[0].fav, same.resources[0].fav);
        assert_eq!(a.resources[1].fav, same.resources[1].fav);
        let b = zone(&c, &parse_rules(&RULES.replace("chance = 40", "chance = 60")).unwrap(), DensitySource::Structural);
        assert_ne!(a.resources[0].fav, b.resources[0].fav, "iron's grid changes");
        assert_eq!(a.resources[1].fav, b.resources[1].fav, "gold's grid does not");
    }

    /// The structural modulation: g(0) = g_min, g(d ≥ d_ref) = 1; the modulated gold rises with the density, and with
    /// C-3b alone (0 here) it stays at g_min × chance.
    #[test]
    fn the_structural_modulation_makes_districts() {
        assert!((modulation(0.0, 0.5, 0.2) - 0.2).abs() < 1e-6);
        assert!((modulation(1.0, 0.5, 0.2) - 1.0).abs() < 1e-6);
        let c = ctx();
        let z = zone(&c, &parse_rules(RULES).unwrap(), DensitySource::Structural);
        let gold = &z.resources[1].fav;
        assert!(gold[15] > gold[8], "more favourable where the density is higher ({} > {})", gold[15], gold[8]);
        let zc = zone(&c, &parse_rules(RULES).unwrap(), DensitySource::C3bOnly);
        assert_eq!(zc.resources[1].fav[15], 7, "C-3b alone at 0: 35 × 0.2 = 7");
    }

    /// ADR Finding 153-I — the inspector shows exactly the export's values: at several cells (a modulated rule, a plain
    /// rule, a placer, an empty cell), `explain`'s favourability equals the zoning grid's and the decoded PNG's, and its
    /// attribution is the rule or the placer that gave it.
    #[test]
    fn the_inspector_shows_the_exports_values() {
        let mut c = ctx();
        c.rivers = (0..6).map(|i| RiverCell { q4: (9 + i) as u32, recv: if i < 5 { (i + 1) as u32 } else { u32::MAX }, step_km: 1.0, slope: 0.001, area_km2: 10.0 }).collect();
        c.wetland[11] = true;
        let text = RULES.replace("modulate = \"structural\"\n", "modulate = \"structural\"\nsource = \"orogenic\"\n[[resource.rule]]\nwetland = true\nchance = 100\nsource = \"wet\"\n")
            + "[resource.placer]\nsource_min_chance = 100\nriver_min_area_km2 = 5.0\nstart_chance = 50\ndecay_per_km = 10.0\nmax_slope = 0.01\nsource = \"gravels\"\n";
        let rules = parse_rules(&text).unwrap();
        let z = zone(&c, &rules, DensitySource::Structural);
        let files = crate::geology::export::export_files(&vec![1u8; 64 * 64], 64, 64, &z, &rules, "test");
        for q in [0usize, 9, 11, 12, 13, 15, 200] {
            let cells = explain(&c, &rules, &z, q, DensitySource::Structural);
            for (k, grid) in z.resources.iter().enumerate() {
                let png = &files.iter().find(|(n, _)| n == &format!("favorabilite_{}.png", grid.id)).unwrap().1;
                let img = image::load_from_memory(png).unwrap().to_luma8();
                let exported = img.as_raw()[q];
                assert_eq!(exported, grid.fav[q], "the PNG holds the grid's value at {q}");
                let shown = cells.iter().find(|r| r.id == grid.id).map_or(0, |r| r.fav);
                assert_eq!(shown, exported, "resource {k} at cell {q}: the inspector shows the export's value");
            }
        }
        let at = |q: usize, id: &str| explain(&c, &rules, &z, q, DensitySource::Structural).into_iter().find(|r| r.id == id).unwrap().origin;
        assert!(matches!(at(15, "gold"), FavOrigin::Rule { modulation: Some(_), .. }), "a modulated rule carries its g");
        assert!(matches!(at(0, "iron"), FavOrigin::Rule { modulation: None, .. }), "a plain rule");
        assert_eq!(at(11, "gold"), FavOrigin::Rule { index: 1, source: Some("wet".into()), modulation: None }, "the source cell's own rule");
        assert_eq!(at(12, "gold"), FavOrigin::Placer { source: Some("gravels".into()) }, "downstream of the source: the placer (40 > the rule's 28)");
        assert!(matches!(at(13, "gold"), FavOrigin::Rule { index: 0, .. }), "a tie (placer 30 = rule 30) goes to the rule, as zone's MAX keeps it");
    }

    /// The placers travel DOWNSTREAM only: a source upstream reaches the cells below it, never the ones above.
    #[test]
    fn placers_go_downstream_only() {
        let mut c = ctx();
        // a river of 6 HD cells along row 0, ¼ cells 0..6, flowing right; the source is the belt cell 9 → none here:
        // put the river on the belt side: q4 = 8..14, upstream first
        c.rivers = (0..6).map(|i| RiverCell { q4: (9 + i) as u32, recv: if i < 5 { (i + 1) as u32 } else { u32::MAX }, step_km: 1.0, slope: 0.001, area_km2: 10.0 }).collect();
        let text = RULES.to_string() + "[resource.placer]\nsource_min_chance = 100\nriver_min_area_km2 = 5.0\nstart_chance = 50\ndecay_per_km = 10.0\nmax_slope = 0.01\n";
        // the source: make cell 11's gold chance 100 by a rule on a wetland there only
        c.wetland[11] = true;
        let text = text.replace("modulate = \"structural\"\n", "modulate = \"structural\"\n[[resource.rule]]\nwetland = true\nchance = 100\n");
        let z = zone(&c, &parse_rules(&text).unwrap(), DensitySource::Structural);
        let gold = &z.resources[1].fav;
        assert!(z.resources[1].placer_cells > 0);
        assert!(gold[12] >= 40 && gold[13] >= 30, "downstream of the source: {} {}", gold[12], gold[13]);
        let upstream = zone(&c, &parse_rules(&text.replace("start_chance = 50", "start_chance = 0")).unwrap(), DensitySource::Structural);
        assert_eq!(gold[9], upstream.resources[1].fav[9], "the cells upstream of the source are untouched");
    }
}
