//! ADR block A — **is the "hypsometry diverges" blocker actually "erosion does a different
//! amount of work per grid"?** Diagnostic only: this bench MEASURES, it asserts nothing.
//!
//! ## Why it exists: the published control figure is wrong twice
//!
//! Findings 56/56b publish "un-eroded reference 865 m" beside hypsometry means at BOTH
//! resolutions. Two defects:
//!
//! 1. **It is a 2048² number.** It comes from a hard-coded header string in
//!    `coastal_comb_levers::channel_head_law_verdict`, sourced from `comb_levers_control_block`,
//!    which pins `let target = 2048usize`. A control quoted without its population.
//! 2. **It is not "un-eroded".** `Lever { reference: true }` turns off FOUR things at once —
//!    `amplitude_base = 0.0` (no FBM), `lithology.enabled = false`, `fracture.enabled = false`
//!    AND `stream_power = None`. So `865 − 282 = 583 m` is not the work of erosion; it is the
//!    sum of four levers, and attributing it to incision is method rule 3's defect (decompose a
//!    measured factor before quoting it).
//!
//! The true control did not exist in the repository. This builds it: the production
//! configuration with the INCISION alone removed, everything else identical.
//!
//! ## The three fields
//!
//! | variant | FBM | lithology | fracture | stream power |
//! |---|---|---|---|---|
//! | `shipped` | on | on | on | **on** |
//! | `no-incision` — the TRUE control | on | on | on | **off** |
//! | `coarse` — what "865" actually is | **off** | **off** | **off** | **off** |
//!
//! ## The three measurements
//!
//! - **A1** mean hypsometry (f64 accumulation — an f32 sum saturates over 1.8 M altitudes and
//!   read 13 m low, Finding 56b) and emerged fraction, per variant per resolution.
//! - **A2** slope quantiles of the un-eroded field, DIMENSIONLESS, in two flavours: the
//!   central-difference gradient (comparable with the earlier blocks) and the **D8
//!   steepest-descent gradient, which is the quantity `incise` actually reads** (`s_now` in
//!   `stream_power.rs`). Attributing an incision's slope sensitivity to a Sobel gradient the
//!   incision never sees would be measuring a different object.
//! - **A3** erosion work = `no-incision − shipped`, in metres, paired per cell over the cells
//!   that are land in BOTH fields, plus the difference of the two means for comparison.
//!
//! Run: cargo test -p ymir-core --release --test hypsometry_work_attribution -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, upscale_from_c1_with_progress,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;

#[derive(Clone, Copy, PartialEq)]
enum Variant {
    /// Production, exactly as shipped.
    Shipped,
    /// Production with the INCISION alone removed — the control that was missing.
    NoIncision,
    /// What the published "un-eroded reference 865 m" actually is: no FBM, no lithology, no
    /// fracture AND no incision.
    Coarse,
}

impl Variant {
    fn label(self) -> &'static str {
        match self {
            Variant::Shipped => "shipped (production)",
            Variant::NoIncision => "NO INCISION (true control)",
            Variant::Coarse => "coarse (what '865' is)",
        }
    }
}

/// A full build spec, so A2bis and C1 can perturb one knob at a time from the same builder.
#[derive(Clone, Copy)]
struct Spec {
    v: Variant,
    /// `None` = the bench default (0.04, matching `coastal_comb_levers`); `Some(a)` overrides.
    amplitude_base: Option<f64>,
    /// `None` = `FbmUpscaleConfig::default()`'s fixed 7; `Some(n)` overrides. The NEGATIVE
    /// CONTROL for A2bis — the instrument must be able to see a roughness change at all.
    octaves: Option<usize>,
    /// `None` = the shipped `RELIEF_V1_A_C_KM2 = 0.1`; `Some(a)` overrides (C1's sweep).
    a_c_km2: Option<f32>,
}

impl Spec {
    fn of(v: Variant) -> Self {
        Spec { v, amplitude_base: None, octaves: None, a_c_km2: None }
    }
}

/// Effective geometry of a build, read back from the config that was actually constructed —
/// block C's arithmetic verified by measurement rather than by recomputation.
struct Geometry {
    cell_km: f32,
    cell_km2: f32,
    min_area_cells: f32,
    octaves: usize,
}

fn terrain_spec(target: usize, sp: Spec) -> (GridF32, Geometry) {
    build(target, sp)
}

fn terrain(target: usize, v: Variant) -> GridF32 {
    build(target, Spec::of(v)).0
}

fn build(target: usize, sp: Spec) -> (GridF32, Geometry) {
    let v = sp.v;
    let ss = SteinSteinParams::default();
    let run_cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run_cfg, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let coarse = v == Variant::Coarse;
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: sp.amplitude_base.unwrap_or(if coarse { 0.0 } else { 0.04 }),
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: !coarse,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: !coarse,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    if let Some(n) = sp.octaves {
        cfg.octaves = n;
    }
    if v != Variant::Shipped {
        cfg.stream_power = None; // the ONE lever that separates `NoIncision` from `Shipped`
    } else if let Some(a) = sp.a_c_km2 {
        // C1's sweep. `min_area_cells` is a CELL COUNT derived from a km² constant, which is
        // the whole point of block C: the same physical threshold is 2.6 cells at 2048² and
        // 41.9 at 8192².
        let cell_km2 = (DOMAIN_KM / target as f32).powi(2);
        if let Some(spw) = cfg.stream_power.as_mut() {
            spw.min_area_cells = a / cell_km2;
        }
    }
    let geo = {
        let cell_km = DOMAIN_KM / target as f32;
        Geometry {
            cell_km,
            cell_km2: cell_km * cell_km,
            // READ BACK from the config, not recomputed — if the two ever disagree, this is
            // the number that governs the incision.
            min_area_cells: cfg.stream_power.as_ref().map_or(f32::NAN, |x| x.min_area_cells),
            octaves: cfg.octaves,
        }
    };
    let f = upscale_from_c1_with_progress(
        &state,
        &run_cfg.iso_config,
        &ss,
        &seed,
        &cfg,
        &edifices,
        &volc,
        Some(&kin),
        &mut |_| {},
        &|| false,
    )
    .0
    .heightmap;
    (f, geo)
}

fn q(v: &mut Vec<f32>, f: f64) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if v.is_empty() {
        return f32::NAN;
    }
    v[(((v.len() - 1) as f64 * f).round() as usize).min(v.len() - 1)]
}

/// A1 — mean over land in **f64** (an f32 running sum saturates past ~10⁹ and read 13 m low at
/// 8192², Finding 56b) plus the emerged fraction.
fn hypsometry(field: &GridF32, ss: &SteinSteinParams) -> (f64, f32, f32, f32, usize, f32) {
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, ss);
    let mut alt: Vec<f32> = field.data.iter().filter(|&&x| x > SEA).map(|&x| to_m(x)).collect();
    let n = alt.len();
    let mean = alt.iter().map(|&x| x as f64).sum::<f64>() / n.max(1) as f64;
    let land_pct = 100.0 * n as f32 / field.data.len() as f32;
    (mean, q(&mut alt.clone(), 0.10), q(&mut alt.clone(), 0.50), q(&mut alt, 0.90), n, land_pct)
}

/// A2 — the two slope definitions, both DIMENSIONLESS (rise / run, i.e. tan of the angle).
///
/// `central`: the magnitude of the central-difference gradient — the quantity the earlier
/// blocks reported. `d8`: the steepest descent to a D8 neighbour, `(h − h_min) / dist`, which
/// is EXACTLY what `incise` computes as `s_now` and therefore the only one that can attribute
/// an incision's slope sensitivity.
fn slopes(field: &GridF32, ss: &SteinSteinParams, target: usize) -> (Vec<f32>, Vec<f32>) {
    let (w, h) = (field.width, field.height);
    let m_per_cell = DOMAIN_KM / target as f32 * 1000.0;
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, ss);
    let step = if w > 4096 { 4 } else { 2 };
    let (mut central, mut d8) = (Vec::new(), Vec::new());
    const DX: [i32; 8] = [1, -1, 0, 0, 1, 1, -1, -1];
    const DY: [i32; 8] = [0, 0, 1, -1, 1, -1, 1, -1];
    for y in (1..h - 1).step_by(step) {
        for x in (1..w - 1).step_by(step) {
            let k = y * w + x;
            if field.data[k] <= SEA {
                continue;
            }
            let gx = 0.5 * (to_m(field.data[k + 1]) - to_m(field.data[k - 1]));
            let gy = 0.5 * (to_m(field.data[k + w]) - to_m(field.data[k - w]));
            central.push((gx * gx + gy * gy).sqrt() / m_per_cell);
            // Steepest D8 descent — `s_now` in `stream_power::incise`.
            let ho = to_m(field.data[k]);
            let mut best = 0.0f32;
            for i in 0..8 {
                let (nx, ny) = (x as i32 + DX[i], y as i32 + DY[i]);
                let hn = to_m(field.data[ny as usize * w + nx as usize]);
                let dist = if i < 4 { m_per_cell } else { m_per_cell * std::f32::consts::SQRT_2 };
                let s = (ho - hn) / dist;
                if s > best {
                    best = s;
                }
            }
            d8.push(best);
        }
    }
    (central, d8)
}

fn report_slopes(tag: &str, s: &mut Vec<f32>) {
    let deg = |x: f32| x.atan().to_degrees();
    let n = s.len().max(1) as f32;
    let gt30 = 100.0 * s.iter().filter(|&&x| deg(x) > 30.0).count() as f32 / n;
    let gt45 = 100.0 * s.iter().filter(|&&x| deg(x) > 45.0).count() as f32 / n;
    let (p10, p50, p90, p99) = (
        q(&mut s.clone(), 0.10),
        q(&mut s.clone(), 0.50),
        q(&mut s.clone(), 0.90),
        q(&mut s.clone(), 0.99),
    );
    eprintln!(
        "    {tag:<9} p10 {p10:.4} ({:.2}°) | p50 {p50:.4} ({:.2}°) | p90 {p90:.4} ({:.2}°) | \
         p99 {p99:.4} ({:.2}°) | >30° {gt30:5.2} % | >45° {gt45:5.2} % | n {}",
        deg(p10),
        deg(p50),
        deg(p90),
        deg(p99),
        s.len()
    );
}

#[test]
#[ignore]
fn hypsometry_work_attribution() {
    let ss = SteinSteinParams::default();
    eprintln!(
        "\n==========  BLOCK A — is the blocker 'hypsometry diverges' or 'erosion does \
         different WORK per grid'?  =========="
    );
    eprintln!(
        "\nThe published control '865 m un-eroded' is a 2048² figure AND is not un-eroded: the \
         `reference`\nlever also removes the FBM, the lithology and the fracture field. The true \
         control is built here."
    );

    for target in [2048usize, 8192] {
        let km = DOMAIN_KM / target as f32;
        eprintln!("\n╔══════ {target}²  ({:.1} m/cell) ══════╗", km * 1000.0);
        let mut fields: Vec<(Variant, GridF32)> = Vec::new();
        for v in [Variant::Shipped, Variant::NoIncision, Variant::Coarse] {
            fields.push((v, terrain(target, v)));
        }

        eprintln!("\n── A1 · hypsometry over land (mean in f64) ──");
        eprintln!(
            "  {:<28} {:>9} {:>9} {:>9} {:>9} {:>12}",
            "variant", "mean m", "p10 m", "p50 m", "p90 m", "land %"
        );
        for (v, f) in &fields {
            let (mean, p10, p50, p90, _n, land) = hypsometry(f, &ss);
            eprintln!(
                "  {:<28} {mean:>9.1} {p10:>9.0} {p50:>9.0} {p90:>9.0} {land:>11.2} %",
                v.label()
            );
        }

        // Are `coarse` and `no-incision` the SAME FIELD? Their statistics came out identical to
        // four figures, which is either a coincidence or a mechanism. The mechanism, if it is
        // one: `amplitude_base` is the recorded DEAD KNOB (the C-1 relief budget cap binds
        // everywhere in production, so the configured amplitude never applies), and lithology
        // and fracture are ERODIBILITY multipliers consumed only by the stream power — which is
        // off in both variants. Checked by field equality rather than by trusting the doc
        // comment that says so (method rule 5: re-derive the mechanism).
        {
            let a = &fields.iter().find(|(v, _)| *v == Variant::NoIncision).unwrap().1;
            let b = &fields.iter().find(|(v, _)| *v == Variant::Coarse).unwrap().1;
            let diff = a.data.iter().zip(b.data.iter()).filter(|(x, y)| x != y).count();
            let maxd =
                a.data.iter().zip(b.data.iter()).map(|(x, y)| (x - y).abs()).fold(0.0f32, f32::max);
            eprintln!(
                "\n── control identity · no-incision vs coarse (FBM + lithology + fracture ON \
                 vs OFF, both without incision) ──\n  differing cells {diff} of {} | max |Δ| \
                 {maxd:.6} norm",
                a.data.len()
            );
        }

        // A2bis — the POPULATION IS NAMED on every row. The first version of this table said
        // only "non érodé", and the row it reported was the FBM-bearing `NoIncision` build;
        // the `Coarse` row was skipped by a `continue`. Reporting a slope quantile without its
        // variant is the same defect as the "865 m" control quoted without its resolution, so
        // both rows are printed here and labelled.
        eprintln!("\n── A2bis · slope of each PRE-INCISION field, dimensionless (tan) ──");
        for (v, f) in &fields {
            if *v == Variant::Shipped {
                continue;
            }
            eprintln!("  [{}]", v.label());
            let (mut c, mut d) = slopes(f, &ss, target);
            report_slopes("central", &mut c);
            report_slopes("D8 (s_now)", &mut d);
        }
        eprintln!("  [shipped (production) — for scale]");
        {
            let f = &fields.iter().find(|(v, _)| *v == Variant::Shipped).unwrap().1;
            let (mut c, mut d) = slopes(f, &ss, target);
            report_slopes("central", &mut c);
            report_slopes("D8 (s_now)", &mut d);
        }

        eprintln!("\n── A3 · erosion WORK = no-incision − shipped, metres ──");
        let shipped = &fields.iter().find(|(v, _)| *v == Variant::Shipped).unwrap().1;
        let control = &fields.iter().find(|(v, _)| *v == Variant::NoIncision).unwrap().1;
        let to_m = |x: f32| c1_altitude_norm_to_metres(x, &ss);
        let mut diff: Vec<f32> = Vec::new();
        let (mut sum, mut n_both) = (0.0f64, 0usize);
        let (mut only_ctrl, mut only_ship) = (0usize, 0usize);
        for k in 0..shipped.data.len() {
            let (a, b) = (control.data[k], shipped.data[k]);
            match (a > SEA, b > SEA) {
                (true, true) => {
                    let d = to_m(a) - to_m(b);
                    sum += d as f64;
                    n_both += 1;
                    diff.push(d);
                }
                (true, false) => only_ctrl += 1,
                (false, true) => only_ship += 1,
                _ => {}
            }
        }
        let paired_mean = sum / n_both.max(1) as f64;
        let (mc, _, _, _, _, _) = hypsometry(control, &ss);
        let (msh, _, _, _, _, _) = hypsometry(shipped, &ss);
        eprintln!(
            "  PAIRED (land in both, n = {n_both})   mean {paired_mean:8.1} m | p10 {:.1} | \
             p50 {:.1} | p90 {:.1} | max {:.1}",
            q(&mut diff.clone(), 0.10),
            q(&mut diff.clone(), 0.50),
            q(&mut diff.clone(), 0.90),
            q(&mut diff.clone(), 1.00)
        );
        eprintln!(
            "  DIFFERENCE OF MEANS                mean {:8.1} m  ({mc:.1} − {msh:.1})",
            mc - msh
        );
        eprintln!(
            "  mask disagreement: land only in the control {only_ctrl} | only in shipped \
             {only_ship}  (erosion moving cells across sea level)"
        );
        let eroded_frac =
            100.0 * diff.iter().filter(|&&d| d > 1.0).count() as f32 / diff.len().max(1) as f32;
        eprintln!(
            "  share of common land LOWERED by more than 1 m by the incision: {eroded_frac:.2} %"
        );
    }
}

/// A2bis' NEGATIVE CONTROL, plus block C's arithmetic read back from the built config.
///
/// Two controls, because the obvious one may legitimately move nothing. `flow_conditioning =
/// 0.1` caps the FBM's total downslope rise via `flow_budget_divisor = nscale · Σ(p·l)^o`, and
/// `Σ` runs over the octave count — so adding an octave RESCALES the amplitude to hold the same
/// total rise, and octave-invariance would be a correct property of the code rather than a blind
/// instrument. `amplitude_base × 4` is therefore the second control: if NEITHER moves the
/// quantiles, the measurement is void and A2bis proves nothing.
#[test]
#[ignore]
fn a2bis_roughness_negative_control() {
    let ss = SteinSteinParams::default();
    eprintln!(
        "\n==========  A2bis — NEGATIVE CONTROL: can this instrument see roughness?  =========="
    );
    for target in [2048usize, 8192] {
        eprintln!("\n╔══════ {target}² ══════╗");
        for (tag, sp) in [
            ("fbm, octaves 7 (shipped)", Spec::of(Variant::NoIncision)),
            ("fbm, octaves 8 (+1)", Spec { octaves: Some(8), ..Spec::of(Variant::NoIncision) }),
            (
                "fbm, amplitude_base x4",
                Spec { amplitude_base: Some(0.16), ..Spec::of(Variant::NoIncision) },
            ),
        ] {
            let (f, geo) = terrain_spec(target, sp);
            let (mean, _, p50, p90, _, land) = hypsometry(&f, &ss);
            let (_, mut d) = slopes(&f, &ss, target);
            eprintln!(
                "  [{tag}]  octaves {} | mean {mean:.1} m | alt p50 {p50:.0} p90 {p90:.0} | \
                 land {land:.2} %",
                geo.octaves
            );
            report_slopes("D8", &mut d);
        }
    }
}

/// Block C — the cell arithmetic READ BACK from the built config, then the `A_c` sweep.
///
/// The control block runs on every sweep point and is deliberately CLIMATE-FREE, so no climate
/// bed has to be chosen: the hypsometry, the closed-depression population, the enclosed
/// below-sea regions (a `water_class` geometry, not a water balance) and the network extent are
/// all functions of the height field alone. The rule "never one bed" applies when climate enters
/// the chain; here it does not, and that is why one pass suffices.
#[test]
#[ignore]
fn c1_channel_head_area_sweep() {
    use ymir_core::lakes::connectivity::water_class;
    use ymir_core::terrain::flow::{FlowConfig, compute_flow, mfd_accumulation};
    let ss = SteinSteinParams::default();
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, &ss);
    eprintln!("\n==========  BLOCK C — A_c is physically dimensioned but SUB-CELL  ==========");

    // ── C · the arithmetic, from the config that was actually constructed ────────────
    eprintln!("\n── C · effective geometry, read back from the built StreamPowerConfig ──");
    eprintln!(
        "  {:<8} {:>12} {:>14} {:>18} {:>18}",
        "grid", "cell km", "cell km²", "A_c=0.1 km² cells", "A_c=1.0 km² cells"
    );
    for target in [2048usize, 8192] {
        let (_, g1) = terrain_spec(target, Spec::of(Variant::Shipped));
        let (_, g2) =
            terrain_spec(target, Spec { a_c_km2: Some(1.0), ..Spec::of(Variant::Shipped) });
        eprintln!(
            "  {target:<8} {:>12.6} {:>14.7} {:>18.4} {:>18.4}",
            g1.cell_km, g1.cell_km2, g1.min_area_cells, g2.min_area_cells
        );
    }

    // ── C1 · the sweep ──────────────────────────────────────────────────────────────
    let mut work: Vec<(usize, f32, f64)> = Vec::new();
    for target in [2048usize, 8192] {
        eprintln!("\n╔══════ {target}² ══════╗");
        // The pre-incision control, once per grid — the reference every work figure is against.
        let (control, _) = terrain_spec(target, Spec::of(Variant::NoIncision));
        let (mc, _, _, _, _, land_c) = hypsometry(&control, &ss);
        eprintln!("  control (no incision): mean {mc:.1} m | land {land_c:.2} %");
        for a_c in [0.1f32, 0.4, 1.0] {
            let (f, g) =
                terrain_spec(target, Spec { a_c_km2: Some(a_c), ..Spec::of(Variant::Shipped) });
            let (mean, _, _, _, _, land) = hypsometry(&f, &ss);
            // paired work over the cells that are land in BOTH
            let (mut sum, mut n_both, mut n_moved) = (0.0f64, 0usize, 0usize);
            for k in 0..f.data.len() {
                if control.data[k] > SEA && f.data[k] > SEA {
                    let d = to_m(control.data[k]) - to_m(f.data[k]);
                    sum += d as f64;
                    n_both += 1;
                    if d > 1.0 {
                        n_moved += 1;
                    }
                }
            }
            let w = sum / n_both.max(1) as f64;
            let moved = 100.0 * n_moved as f32 / n_both.max(1) as f32;
            work.push((target, a_c, w));
            // ── control block, climate-free ──
            let flow = compute_flow(&f, &FlowConfig { sea_level: SEA, ..Default::default() });
            let pits = flow
                .filled
                .data
                .iter()
                .zip(f.data.iter())
                .filter(|(a, b)| **a > **b + 1e-6)
                .count();
            let wc = water_class(&f, SEA);
            let below = wc.iter().filter(|&&c| c == 2).count();
            let acc = mfd_accumulation(&flow.filled, &flow.direction, SEA, 2.0, f.width, f.height);
            let head_cells = 0.1 / g.cell_km2; // the SHIPPED channel head, held fixed
            let chan = acc.data.iter().filter(|&&x| x >= head_cells).count();
            eprintln!(
                "  A_c {a_c:>4.1} km² ({:>7.2} cells)  mean {mean:>7.1} m | work {w:>7.1} m | \
                 >1 m {moved:>6.2} % | land {land:>6.2} %\n      control: pits {pits:>9} \
                 ({:>7.0} km²) | below-sea cells {below:>8} ({:>7.0} km²) | network {chan:>9} \
                 cells ({:>8.0} km)",
                g.min_area_cells,
                pits as f32 * g.cell_km2,
                below as f32 * g.cell_km2,
                chan as f32 * g.cell_km
            );
        }
    }
    eprintln!("\n── C1 · work ratio 2048² / 8192², per A_c ──");
    for a_c in [0.1f32, 0.4, 1.0] {
        let lo = work.iter().find(|(t, a, _)| *t == 2048 && *a == a_c).map(|x| x.2).unwrap_or(0.0);
        let hi = work.iter().find(|(t, a, _)| *t == 8192 && *a == a_c).map(|x| x.2).unwrap_or(0.0);
        eprintln!(
            "  A_c {a_c:>4.1} km²   2048² {lo:>7.1} m   8192² {hi:>7.1} m   ratio {:>6.2}",
            lo / hi.max(1e-9)
        );
    }
}
