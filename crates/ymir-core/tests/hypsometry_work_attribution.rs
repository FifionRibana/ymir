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
    /// D1 — the HILLSLOPE ANCHOR, borrowed by hand: multiply `diffusion` by this factor, on the
    /// shipped LINEAR branch, and raise `diffusion_substeps` to hold the stability bound the
    /// config documents (`diffusion / substeps <= 0.2`). `Some((HILLSLOPE_REF_CELL_M/cell_m)^2)`
    /// is numerically what "borrowing the rescaling that only exists in the NONLINEAR branch"
    /// means. Diagnostic sign test; nothing in production changes.
    diffusion_scale: Option<f32>,
    /// Block 4 — set `diffusion` ABSOLUTELY (0.0 switches the hillslope stage off). Distinct
    /// from `diffusion_scale`, which multiplies it.
    diffusion_abs: Option<f32>,
    /// Block 1b — `Some(None)` = pure D8, `Some(Some(p))` = MFD with exponent `p`, `None` =
    /// the shipped value. DIAG only; production is never touched.
    mfd: Option<Option<f32>>,
    /// Iteration-convergence block — override `iterations`. `iterations` is a DURATION dial
    /// (Finding 44: the observable is `k_time = K·dt·iterations`), so this sweeps the erosion
    /// budget, not a numerical tolerance.
    iterations: Option<usize>,
}

impl Spec {
    fn of(v: Variant) -> Self {
        Spec {
            v,
            amplitude_base: None,
            octaves: None,
            a_c_km2: None,
            diffusion_scale: None,
            diffusion_abs: None,
            mfd: None,
            iterations: None,
        }
    }
}

/// Effective geometry of a build, read back from the config that was actually constructed —
/// block C's arithmetic verified by measurement rather than by recomputation.
struct Geometry {
    cell_km: f32,
    cell_km2: f32,
    min_area_cells: f32,
    octaves: usize,
    /// Read back, because the E inventory quoted 0.05 (relief-v1) where production ships
    /// relief-v3's 0.08 — a figure that must come from the config, not from a doc comment.
    diffusion: f32,
    diffusion_substeps: usize,
    k: f32,
    dt: f32,
    iterations: usize,
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
    if let Some(it) = sp.iterations {
        if let Some(spw) = cfg.stream_power.as_mut() {
            spw.iterations = it;
        }
    }
    if let Some(d) = sp.diffusion_abs {
        if let Some(spw) = cfg.stream_power.as_mut() {
            spw.diffusion = d;
        }
    }
    if let Some(p) = sp.mfd {
        if let Some(spw) = cfg.stream_power.as_mut() {
            spw.mfd_exponent = p;
        }
    }
    if let Some(scale) = sp.diffusion_scale {
        if let Some(spw) = cfg.stream_power.as_mut() {
            spw.diffusion *= scale;
            // Hold `diffusion / substeps <= 0.2`, the bound the field's own doc states.
            while spw.diffusion / spw.diffusion_substeps as f32 > 0.2 {
                spw.diffusion_substeps *= 2;
            }
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
            diffusion: cfg.stream_power.as_ref().map_or(f32::NAN, |x| x.diffusion),
            diffusion_substeps: cfg.stream_power.as_ref().map_or(0, |x| x.diffusion_substeps),
            k: cfg.stream_power.as_ref().map_or(f32::NAN, |x| x.k),
            dt: cfg.stream_power.as_ref().map_or(f32::NAN, |x| x.dt),
            iterations: cfg.stream_power.as_ref().map_or(0, |x| x.iterations),
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

/// Blocks A1/A2/B1/C1/D1/E — the intensive / extensive decomposition, on a DEFINED denominator.
///
/// ## The denominator, and why this one
///
/// The decomposition was first built on "land lowered by more than 1 m", a proxy with an
/// arbitrary threshold. The regime AUTHORITY is the criterion the incision actually applies:
/// `accumulation >= min_area_cells`. It is evaluated on the **PRE-INCISION** field, on purpose:
///
/// - it is the INPUT partition, so it cannot be moved by the outcome being measured;
/// - evaluating it on the eroded field would make the denominator a function of the numerator,
///   which is exactly what made the "network extent" control column unusable on an `A_c` sweep
///   (read at a FIXED threshold on the RESULTING field).
///
/// Run: cargo test -p ymir-core --release --test hypsometry_work_attribution -- --ignored --nocapture decomposition
#[test]
#[ignore]
fn intensive_extensive_decomposition() {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow, mfd_accumulation};
    let ss = SteinSteinParams::default();
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, &ss);
    eprintln!("\n==========  A1/A2/B1/C1/D1/E — decomposition on the REGIME criterion  ==========");

    struct Pre {
        field: GridF32,
        acc_cells: GridF32,
        pits: usize,
        cell_km2: f32,
    }
    let mut pre: Vec<(usize, Pre)> = Vec::new();
    for target in [2048usize, 8192] {
        let (field, g) = terrain_spec(target, Spec::of(Variant::NoIncision));
        let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
        let pits = flow
            .filled
            .data
            .iter()
            .zip(field.data.iter())
            .filter(|(a, b)| **a > **b + 1e-6)
            .count();
        let acc_cells =
            mfd_accumulation(&flow.filled, &flow.direction, SEA, 2.0, field.width, field.height);
        pre.push((target, Pre { field, acc_cells, pits, cell_km2: g.cell_km2 }));
    }

    // ── C1 · the accumulation field itself, never measured before ───────────────────
    eprintln!("\n── C1 · PRE-INCISION accumulation over land, km² (then A^m, m = 0.5) ──");
    eprintln!(
        "  {:<8} {:>10} {:>10} {:>10} {:>10} {:>10} {:>11} {:>12}",
        "grid", "p10", "p25", "p50", "p75", "p90", "p99", "= 1 cell"
    );
    for (target, p) in &pre {
        let mut a: Vec<f32> = (0..p.field.data.len())
            .filter(|&k| p.field.data[k] > SEA)
            .map(|k| p.acc_cells.data[k] * p.cell_km2)
            .collect();
        let one = 100.0 * a.iter().filter(|&&x| x <= p.cell_km2 * 1.001).count() as f32
            / a.len().max(1) as f32;
        let qs: Vec<f32> =
            [0.10, 0.25, 0.50, 0.75, 0.90, 0.99].iter().map(|&f| q(&mut a.clone(), f)).collect();
        eprintln!(
            "  {target:<8} {:>10.6} {:>10.6} {:>10.6} {:>10.6} {:>10.5} {:>11.4} {:>11.2} %",
            qs[0], qs[1], qs[2], qs[3], qs[4], qs[5], one
        );
        eprintln!(
            "  {:<8} {:>10.6} {:>10.6} {:>10.6} {:>10.6} {:>10.5} {:>11.4}   (A^0.5)",
            "",
            qs[0].sqrt(),
            qs[1].sqrt(),
            qs[2].sqrt(),
            qs[3].sqrt(),
            qs[4].sqrt(),
            qs[5].sqrt()
        );
    }
    eprintln!("\n── E · closed depressions on the PRE-INCISION field ──");
    for (target, p) in &pre {
        eprintln!(
            "  {target}²: pre-incision pits {:>9} ({:>7.0} km²)",
            p.pits,
            p.pits as f32 * p.cell_km2
        );
    }

    // ── A1 / A2 / B1 / D1 ───────────────────────────────────────────────────────────
    eprintln!("\n── A1/A2/B1/D1 · work conditioned on the REGIME (acc_pre >= A_c) ──");
    let cases: [(&str, f32, bool, Option<f32>); 6] = [
        ("A_c 0.025 km2   [A2]", 0.025, false, None),
        ("A_c 0.100 km2   [shipped]", 0.1, false, None),
        ("A_c 0.400 km2", 0.4, false, None),
        ("A_c 1.000 km2", 1.0, false, None),
        // B1 — `A_c` held at 2.6214 CELLS on both grids: 0.1 km² at 2048², 0.00625 at 8192².
        ("A_c 2.62 CELLS  [B1]", 2.6214, true, None),
        // D1 — the anchor borrowed by hand. The scale is 1 at 2048² (it IS the reference cell)
        // and (195.31/48.83)^2 = 16 at 8192².
        ("A_c 0.100 + ANCHOR [D1]", 0.1, false, Some(0.0)),
    ];
    let mut rows: Vec<(String, usize, f64, f64, f64, f32)> = Vec::new();
    for (label, val, in_cells, anchor) in cases {
        for (target, p) in &pre {
            let a_c_km2 = if in_cells { val * p.cell_km2 } else { val };
            let a_c_cells = a_c_km2 / p.cell_km2;
            let scale = anchor.map(|_| {
                let cell_m = DOMAIN_KM / *target as f32 * 1000.0;
                (195.3125f32 / cell_m).powi(2)
            });
            let (f, g) = build(
                *target,
                Spec {
                    a_c_km2: Some(a_c_km2),
                    diffusion_scale: scale,
                    ..Spec::of(Variant::Shipped)
                },
            );
            let (mut sc, mut nc, mut sh, mut nh, mut n_moved) =
                (0.0f64, 0usize, 0.0f64, 0usize, 0usize);
            for k in 0..f.data.len() {
                if p.field.data[k] <= SEA || f.data[k] <= SEA {
                    continue;
                }
                let d = (to_m(p.field.data[k]) - to_m(f.data[k])) as f64;
                if d > 1.0 {
                    n_moved += 1;
                }
                if p.acc_cells.data[k] >= a_c_cells {
                    sc += d;
                    nc += 1;
                } else {
                    sh += d;
                    nh += 1;
                }
            }
            let n = nc + nh;
            let ext = 100.0 * nc as f64 / n.max(1) as f64;
            let (ic, ih) = (sc / nc.max(1) as f64, sh / nh.max(1) as f64);
            let total = (sc + sh) / n.max(1) as f64;
            eprintln!(
                "  [{label}] {target}²  A_c {a_c_km2:.6} km2 = {a_c_cells:>8.2} cells | \
                 diffusion {:.3}/{} substeps\n      extensive(channel) {ext:>6.2} % | intensive \
                 CHANNEL {ic:>7.1} m | intensive HILLSLOPE {ih:>7.1} m | total {total:>7.1} m | \
                 >1 m {:>6.2} %",
                g.diffusion,
                g.diffusion_substeps,
                100.0 * n_moved as f32 / n.max(1) as f32
            );
            rows.push((label.to_string(), *target, ic, ih, total, ext as f32));
        }
    }
    eprintln!("\n── the INTENSIVE ratios, 2048sq / 8192sq ──");
    eprintln!(
        "  {:<28} {:>10} {:>10} {:>12} {:>18}",
        "case", "CHANNEL", "hillslope", "total(raw)", "extensive lo/hi %"
    );
    for (label, _, _, _) in cases.iter().map(|(l, a, b, c)| (l, a, b, c)) {
        let g = |t: usize| rows.iter().find(|r| r.0 == *label && r.1 == t).cloned();
        if let (Some(lo), Some(hi)) = (g(2048), g(8192)) {
            eprintln!(
                "  {label:<28} {:>10.2} {:>10.2} {:>12.2} {:>8.2} / {:<8.2}",
                lo.2 / hi.2.abs().max(1e-9),
                lo.3 / hi.3.abs().max(1e-9),
                lo.4 / hi.4.abs().max(1e-9),
                lo.5,
                hi.5
            );
        }
    }
}

/// A2bis' negative control at the SCALE OF THE CLAIM, not at the scale of the biggest lever.
///
/// The first control only fired on FBM on/off — a change of 674 433 cells and 37 m — while the
/// reading it licensed lives at 3–4 % on the slope tails. A control that only proves sensitivity
/// to the largest available change does not license a per-cent reading. This one perturbs the
/// PRE-INCISION field with a deterministic checkerboard sized to move the D8 p90 slope by
/// roughly 5 %, and asserts the instrument reports it.
///
/// It is a synthetic perturbation on purpose: it isolates the instrument's sensitivity from any
/// question about what the pipeline would do with a rougher input.
#[test]
#[ignore]
fn a2bis_control_sees_a_five_percent_roughness_shift() {
    let ss = SteinSteinParams::default();
    let target = 2048usize;
    let (base, _) = terrain_spec(target, Spec::of(Variant::NoIncision));
    let m_per_cell = DOMAIN_KM / target as f32 * 1000.0;
    let norm_to_m = 2.0 * 1.13 * ss.depth_scale_m as f32;

    let (_, mut d0) = slopes(&base, &ss, target);
    let p90_0 = q(&mut d0, 0.90);
    // A one-cell-wavelength perturbation adds a gradient of about `2·amp/cell_m` where it is
    // resolved, so to lift p90 by ~5 % the amplitude is set from the measured p90 itself.
    let target_delta = 0.05 * p90_0;
    let amp_m = target_delta * m_per_cell / 2.0;
    let amp_norm = amp_m / norm_to_m;
    let mut bumped = base.clone();
    let w = bumped.width;
    for k in 0..bumped.data.len() {
        if bumped.data[k] > SEA {
            let (x, y) = (k % w, k / w);
            let sign = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
            bumped.data[k] += sign * amp_norm;
        }
    }
    let (_, mut d1) = slopes(&bumped, &ss, target);
    let p90_1 = q(&mut d1, 0.90);
    let rel = (p90_1 - p90_0) / p90_0;
    eprintln!(
        "[A2bis control] D8 p90 {p90_0:.4} → {p90_1:.4} ({:+.2} %) for a {amp_m:.2} m \
         one-cell perturbation",
        100.0 * rel
    );
    assert!(
        rel.abs() > 0.02,
        "the slope instrument does not register a perturbation designed to move p90 by 5 % \
         (measured {:+.2} %). The 3–4 % reading in A2bis is then unlicensed and must not be \
         used to argue that the pre-incision field is invariant.",
        100.0 * rel
    );
}

/// Blocks 1a / 1b / 2 / 4 — is the MFD accumulation field's divergence a failure to CONVERGE,
/// or a per-cell PARTITION of a conserved quantity?
///
/// The quantile comparison of Finding 58 compares two different cell populations (11 M against
/// 627 k). The outlet of a matched basin compares the same object, and the basin's cell count is
/// the physical area both grids must agree on. Run first, because it can invalidate the
/// quantile reading.
///
/// Run: cargo test -p ymir-core --release --test hypsometry_work_attribution -- --ignored --nocapture accumulation_convergence
#[test]
#[ignore]
fn accumulation_convergence() {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow, mfd_accumulation};
    let ss = SteinSteinParams::default();
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, &ss);
    eprintln!("\n==========  1a / 1b / 2 / 4 — the accumulation operator  ==========");

    struct G {
        target: usize,
        cell_km2: f32,
        field: GridF32,
        filled: GridF32,
        dir: Vec<u8>,
        basins: Vec<u32>,
        acc: GridF32,
        land: usize,
    }
    let mut gs: Vec<G> = Vec::new();
    for target in [2048usize, 8192] {
        let (field, g) = terrain_spec(target, Spec::of(Variant::NoIncision));
        let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
        let acc =
            mfd_accumulation(&flow.filled, &flow.direction, SEA, 2.0, field.width, field.height);
        let land = field.data.iter().filter(|&&x| x > SEA).count();
        gs.push(G {
            target,
            cell_km2: g.cell_km2,
            field,
            filled: flow.filled,
            dir: flow.direction,
            basins: flow.basins,
            acc,
            land,
        });
    }
    for g in &gs {
        eprintln!(
            "  {}²: land {} cells = {:.0} km² (both grids must drain the same physical area)",
            g.target,
            g.land,
            g.land as f32 * g.cell_km2
        );
    }

    // ── 1a · outlet-matched basins ──────────────────────────────────────────────────
    //
    // MATCHING METHOD, declared: for each of the largest basins at 2048² (by CELL COUNT, the
    // physical area), take its outlet — the land cell of that basin with maximum accumulation —
    // and convert to normalised domain coordinates. At 8192², search a window of ±8 HD cells
    // (±390 m, i.e. ±2 coarse cells) around that position for the maximum accumulation, and
    // report the DISPLACEMENT of the match. A match landing on the window edge is flagged: it
    // means the tolerance, not the topography, chose the point.
    const WIN_HD: i32 = 8;
    eprintln!(
        "\n── 1a · outlet-matched basins (matching: ±{WIN_HD} HD cells = ±{:.0} m around the \
         2048² outlet's normalised position) ──",
        WIN_HD as f32 * 400_000.0 / 8192.0
    );
    let lo = &gs[0];
    let hi = &gs[1];
    let mut by_basin: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for k in 0..lo.field.data.len() {
        if lo.field.data[k] > SEA && lo.basins[k] != 0 {
            *by_basin.entry(lo.basins[k]).or_default() += 1;
        }
    }
    let mut ranked: Vec<(u32, usize)> = by_basin.into_iter().collect();
    ranked.sort_by(|a, b| b.1.cmp(&a.1));
    // and the same census at HD, for the cell-count areas
    let mut hi_cells: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for k in 0..hi.field.data.len() {
        if hi.field.data[k] > SEA && hi.basins[k] != 0 {
            *hi_cells.entry(hi.basins[k]).or_default() += 1;
        }
    }
    eprintln!(
        "  {:>3} {:>12} {:>12} {:>8} {:>12} {:>12} {:>8} {:>9}",
        "#",
        "acc_lo km²",
        "acc_hi km²",
        "ratio",
        "cells_lo km²",
        "cells_hi km²",
        "ratio",
        "shift m"
    );
    let (mut acc_ratios, mut area_ratios) = (Vec::new(), Vec::new());
    for (rank, (bid, ncells)) in ranked.iter().take(10).enumerate() {
        // outlet at 2048²: max accumulation within the basin
        let (mut best_k, mut best_a) = (usize::MAX, -1.0f32);
        for k in 0..lo.field.data.len() {
            if lo.basins[k] == *bid && lo.field.data[k] > SEA && lo.acc.data[k] > best_a {
                best_a = lo.acc.data[k];
                best_k = k;
            }
        }
        if best_k == usize::MAX {
            continue;
        }
        let (ox, oy) = (best_k % lo.field.width, best_k / lo.field.width);
        // normalised position → HD cell
        let fx = (ox as f32 + 0.5) / lo.field.width as f32;
        let fy = (oy as f32 + 0.5) / lo.field.height as f32;
        let (cx, cy) = ((fx * hi.field.width as f32) as i32, (fy * hi.field.height as f32) as i32);
        let (mut hk, mut ha) = (usize::MAX, -1.0f32);
        for dy in -WIN_HD..=WIN_HD {
            for dx in -WIN_HD..=WIN_HD {
                let (x, y) = (cx + dx, cy + dy);
                if x < 0 || y < 0 || x as usize >= hi.field.width || y as usize >= hi.field.height {
                    continue;
                }
                let k = y as usize * hi.field.width + x as usize;
                if hi.field.data[k] > SEA && hi.acc.data[k] > ha {
                    ha = hi.acc.data[k];
                    hk = k;
                }
            }
        }
        if hk == usize::MAX {
            eprintln!("  {:>3}  NO LAND in the match window — unmatched", rank + 1);
            continue;
        }
        let (hx, hy) = (hk % hi.field.width, hk / hi.field.width);
        let shift_m = (((hx as i32 - cx).pow(2) + (hy as i32 - cy).pow(2)) as f32).sqrt()
            * 400_000.0
            / 8192.0;
        let edge = (hx as i32 - cx).abs() == WIN_HD || (hy as i32 - cy).abs() == WIN_HD;
        let (a_lo, a_hi) = (best_a * lo.cell_km2, ha * hi.cell_km2);
        let c_lo = *ncells as f32 * lo.cell_km2;
        let c_hi = hi_cells.get(&hi.basins[hk]).copied().unwrap_or(0) as f32 * hi.cell_km2;
        acc_ratios.push(a_lo / a_hi.max(1e-9));
        area_ratios.push(c_lo / c_hi.max(1e-9));
        eprintln!(
            "  {:>3} {:>12.2} {:>12.2} {:>8.3} {:>12.0} {:>12.0} {:>8.3} {:>9.0}{}",
            rank + 1,
            a_lo,
            a_hi,
            a_lo / a_hi.max(1e-9),
            c_lo,
            c_hi,
            c_lo / c_hi.max(1e-9),
            shift_m,
            if edge { "  ⚠ AT WINDOW EDGE" } else { "" }
        );
    }
    let med = |v: &mut Vec<f32>| q(v, 0.5);
    eprintln!(
        "  MEDIAN over the matched basins:  outlet accumulation ratio {:.3} | basin cell-count \
         area ratio {:.3}",
        med(&mut acc_ratios.clone()),
        med(&mut area_ratios.clone())
    );

    // ── 1b · the exponent sweep, on the SAME pre-incision fields (no rebuild) ───────
    eprintln!("\n── 1b · MFD exponent sweep — accumulation quantiles in km², and the ratio ──");
    eprintln!(
        "  {:<12} {:>10} {:>10} {:>10} {:>10} {:>10} {:>11}",
        "scheme", "p10", "p25", "p50", "p75", "p90", "p99"
    );
    for label in ["D8", "p = 2", "p = 4", "p = 8"] {
        let mut qq: Vec<Vec<f32>> = Vec::new();
        for g in &gs {
            let a = match label {
                "D8" => {
                    // the plain D8 accumulation `compute_flow` already produced
                    let f = compute_flow(
                        &g.field,
                        &FlowConfig { sea_level: SEA, ..Default::default() },
                    );
                    f.accumulation
                }
                other => {
                    let p: f32 = other.trim_start_matches("p = ").parse().unwrap();
                    mfd_accumulation(&g.filled, &g.dir, SEA, p, g.field.width, g.field.height)
                }
            };
            let mut v: Vec<f32> = (0..g.field.data.len())
                .filter(|&k| g.field.data[k] > SEA)
                .map(|k| a.data[k] * g.cell_km2)
                .collect();
            qq.push(
                [0.10, 0.25, 0.50, 0.75, 0.90, 0.99]
                    .iter()
                    .map(|&f| q(&mut v.clone(), f))
                    .collect(),
            );
        }
        for (i, g) in gs.iter().enumerate() {
            eprintln!(
                "  {:<12} {:>10.6} {:>10.6} {:>10.6} {:>10.6} {:>10.5} {:>11.4}   [{}²]",
                label, qq[i][0], qq[i][1], qq[i][2], qq[i][3], qq[i][4], qq[i][5], g.target
            );
        }
        eprintln!(
            "  {:<12} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>11.2}   RATIO lo/hi",
            "",
            qq[0][0] / qq[1][0],
            qq[0][1] / qq[1][1],
            qq[0][2] / qq[1][2],
            qq[0][3] / qq[1][3],
            qq[0][4] / qq[1][4],
            qq[0][5] / qq[1][5]
        );
    }

    // ── 2 · the threshold that EQUALISES the channel-regime share, both directions ──
    //
    // Pure arithmetic on the pre-incision accumulation — no terrain rebuild. Note this uses
    // each grid's OWN land mask rather than the land-in-both mask of Finding 58's table, which
    // is why the shares differ slightly from it.
    eprintln!("\n── 2 · the A_c that equalises the channel-regime share ──");
    let share = |g: &G, a_c_km2: f32| -> f32 {
        let cells = a_c_km2 / g.cell_km2;
        let n = (0..g.field.data.len())
            .filter(|&k| g.field.data[k] > SEA && g.acc.data[k] >= cells)
            .count();
        100.0 * n as f32 / g.land.max(1) as f32
    };
    let bisect = |g: &G, target_share: f32| -> f32 {
        let (mut a, mut b) = (1e-4f32, 100.0f32);
        for _ in 0..40 {
            let m = (a * b).sqrt(); // geometric bisection: the quantity spans decades
            if share(g, m) > target_share {
                a = m;
            } else {
                b = m;
            }
        }
        (a * b).sqrt()
    };
    let (s_lo, s_hi) = (share(lo, 0.1), share(hi, 0.1));
    eprintln!("  shipped A_c = 0.1 km²: channel share {s_lo:.2} % at 2048², {s_hi:.2} % at 8192²");
    let up = bisect(lo, s_hi);
    let down = bisect(hi, s_lo);
    eprintln!(
        "  RAISING 2048² to match HD's {s_hi:.2} %:  A_c = {up:.5} km² = {:.2} cells  → factor \
         {:.2}× on the shipped 0.1",
        up / lo.cell_km2,
        up / 0.1
    );
    eprintln!(
        "  LOWERING 8192² to match 2048²'s {s_lo:.2} %:  A_c = {down:.5} km² = {:.2} cells  → \
         factor {:.2}× on the shipped 0.1",
        down / hi.cell_km2,
        0.1 / down
    );
    eprintln!(
        "  the two directions disagree by {:.2}× — that asymmetry IS the result, not either value",
        // BUG FIXED: this multiplied the two factors (giving 0.57) instead of taking their
        // ratio. The asymmetry is 8.17 / 4.67 = 1.75x.
        (0.1 / down) / (up / 0.1)
    );

    // ── 4 · orientation: is the hillslope-conditioned work diffusion, or incision? ──
    eprintln!("\n── 4 · hillslope-conditioned work with the DIFFUSION SWITCHED OFF ──");
    eprintln!(
        "  `incise` does `continue` on a sub-threshold cell, so it never writes one; but the \
         accumulation is\n  recomputed at EACH of the 2 iterations, so a cell labelled hillslope \
         on the PRE-INCISION field can\n  still pass the gate later. That residual is what this \
         separates."
    );
    for g in &gs {
        let a_c_cells = 0.1 / g.cell_km2;
        for (tag, diff) in [("shipped (diffusion 0.08)", None), ("diffusion 0.0", Some(0.0f32))] {
            let (f, gg) = build(
                g.target,
                Spec { a_c_km2: Some(0.1), diffusion_abs: diff, ..Spec::of(Variant::Shipped) },
            );
            let (mut sh, mut nh, mut sc, mut nc) = (0.0f64, 0usize, 0.0f64, 0usize);
            for k in 0..f.data.len() {
                if g.field.data[k] <= SEA || f.data[k] <= SEA {
                    continue;
                }
                let d = (to_m(g.field.data[k]) - to_m(f.data[k])) as f64;
                if g.acc.data[k] >= a_c_cells {
                    sc += d;
                    nc += 1;
                } else {
                    sh += d;
                    nh += 1;
                }
            }
            eprintln!(
                "  {}²  [{tag}, read back {:.3}]  intensive HILLSLOPE {:>7.1} m | intensive \
                 CHANNEL {:>7.1} m",
                g.target,
                gg.diffusion,
                sh / nh.max(1) as f64,
                sc / nc.max(1) as f64
            );
        }
    }
}

/// BLOCKING BLOCK — of what state is the delivered field the state?
///
/// `iterations = 2` is fixed, `dt = 1.0` is a unit placeholder (Finding 44), and the accumulation
/// is recomputed at each iteration, so the channel/hillslope partition is reorganised between
/// them (Finding 59 §4). The question this answers: does the erosion STABILISE over the iteration
/// count, or is the shipped field a snapshot of an unconverged trajectory?
///
/// The answer changes what the resolution-convergence programme is comparing, so it is measured
/// before anything else this round.
///
/// Note what CANNOT converge here by construction: `E = K·A^m·S^n` carries **no uplift term**, so
/// the only fixed point of the implicit update `h ← (h + f·h_r)/(1+f)` is `h = h_r` everywhere,
/// i.e. base level. The config says so itself ("iters=3 planed them toward base level"). So this
/// sweep measures a DURATION, and "convergence" would mean planation.
///
/// Run: cargo test -p ymir-core --release --test hypsometry_work_attribution -- --ignored --nocapture iteration_convergence
#[test]
#[ignore]
fn iteration_convergence() {
    use ymir_core::lakes::connectivity::water_class;
    use ymir_core::terrain::flow::{FlowConfig, compute_flow, mfd_accumulation};
    let ss = SteinSteinParams::default();
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, &ss);
    eprintln!(
        "\n==========  BLOCKING — is the delivered field converged in ITERATIONS?  =========="
    );

    for target in [2048usize, 8192] {
        eprintln!("\n╔══════ {target}² ══════╗");
        // the pre-incision control, and its own accumulation + pit baseline
        let (pre, g0) = terrain_spec(target, Spec::of(Variant::NoIncision));
        let flow0 = compute_flow(&pre, &FlowConfig { sea_level: SEA, ..Default::default() });
        let acc0 =
            mfd_accumulation(&flow0.filled, &flow0.direction, SEA, 2.0, pre.width, pre.height);
        let pits0 =
            flow0.filled.data.iter().zip(pre.data.iter()).filter(|(a, b)| **a > **b + 1e-6).count();
        let a_c_cells = 0.1 / g0.cell_km2;
        let n_in =
            (0..pre.data.len()).filter(|&k| pre.data[k] > SEA && acc0.data[k] >= a_c_cells).count();
        let land0 = pre.data.iter().filter(|&&x| x > SEA).count();
        let (m0, _, _, _, _, l0) = hypsometry(&pre, &ss);
        eprintln!(
            "  control (no incision): mean {m0:.1} m | land {l0:.2} % | INPUT channel share \
             {:.2} % | pits {pits0}",
            100.0 * n_in as f32 / land0.max(1) as f32
        );

        let mut prev_work: Option<f64> = None;
        for iters in [1usize, 2, 4, 8] {
            let (f, g) = terrain_spec(
                target,
                Spec { iterations: Some(iters), ..Spec::of(Variant::Shipped) },
            );
            // paired work over land-in-both, as Finding 58
            let (mut sum, mut n_both) = (0.0f64, 0usize);
            for k in 0..f.data.len() {
                if pre.data[k] > SEA && f.data[k] > SEA {
                    sum += (to_m(pre.data[k]) - to_m(f.data[k])) as f64;
                    n_both += 1;
                }
            }
            let work = sum / n_both.max(1) as f64;
            let (mean, _, _, _, _, land) = hypsometry(&f, &ss);
            // OUTPUT channel share: the accumulation of the DELIVERED field, which is the
            // partition the last iteration actually acted on.
            let fl = compute_flow(&f, &FlowConfig { sea_level: SEA, ..Default::default() });
            let acc1 = mfd_accumulation(&fl.filled, &fl.direction, SEA, 2.0, f.width, f.height);
            let land1 = f.data.iter().filter(|&&x| x > SEA).count();
            let n_out =
                (0..f.data.len()).filter(|&k| f.data[k] > SEA && acc1.data[k] >= a_c_cells).count();
            let pits =
                fl.filled.data.iter().zip(f.data.iter()).filter(|(a, b)| **a > **b + 1e-6).count();
            let below = water_class(&f, SEA).iter().filter(|&&c| c == 2).count();
            let d = prev_work.map(|p| work - p);
            eprintln!(
                "  iters {iters:>2} (k_time = K·dt·it = {:.0})  work {work:>7.1} m{}  | mean \
                 {mean:>7.1} m | land {land:>6.2} % | OUT channel {:>6.2} % (IN {:>6.2} %)\n      \
                 control: pits {pits:>9} (pre {pits0}, ×{:.2}) | below-sea {below:>8} cells",
                g.k * g.dt * g.iterations as f32,
                match d {
                    Some(x) => format!("  Δ {x:+7.1}"),
                    None => "          ".to_string(),
                },
                100.0 * n_out as f32 / land1.max(1) as f32,
                100.0 * n_in as f32 / land0.max(1) as f32,
                pits as f32 / pits0.max(1) as f32
            );
            prev_work = Some(work);
        }
    }
    eprintln!(
        "\n  READING: the ratio 2048²/8192² at each iteration count is the corollary — if it \
         moves,\n  the blocker's ×2.42 is partly a compute-budget artefact rather than a model \
         property."
    );
}
