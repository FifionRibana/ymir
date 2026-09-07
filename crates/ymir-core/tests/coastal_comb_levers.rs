//! ADR Finding 55 — the coastal comb: characterise the APRON, then test whether the levers that
//! cured the hillslope comb (Findings 8/9/11) reach it.
//!
//! The nuance that explains why nobody saw this coming: the Smith–Bretherton comb WAS solved on
//! the hillslopes with MFD. This is the same instability in a DIFFERENT REGIME — very low
//! gradient, close to base level, on the apron the incision itself manufactures (Finding 51,
//! stage 4). MFD disperses flow; on a near-flat apron dispersion may not suffice to prevent
//! parallel channelisation. So this is not an MFD regression, it is MFD in a regime it was never
//! tested in — and the levers are therefore a hypothesis, not a fix.
//!
//! PRIMARY METRICS, and the only two being optimised:
//!   • indentation DENSITY per 100 km — target ≈ 1 (production 25.8–29.8)
//!   • length p90 in km            — target ≈ 24 (production 3.5–4.6)
//! ⚠️ Spurs must get LONGER and ~25× FEWER. Everything else (coastline length, spacing CV,
//! share) is REPORTED ONLY; a remedy that SHORTENS spurs moves AWAY from the target.
//!
//! Run: cargo test -p ymir-core --release --test coastal_comb_levers -- --ignored --nocapture

use ymir_core::erosion::stream_power::StreamPowerConfig;
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
use ymir_core::terrain::coast_metrics::{
    CoastShape, coast_shape, densest_window, render_coast_crop,
};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;
const CROP: usize = 640;
/// The apron: land within this altitude of sea level — the terrain the coastline runs through.
const APRON_M: f32 = 50.0;

/// One lever setting. `None` fields keep the shipped value.
#[derive(Clone, Copy)]
struct Lever {
    label: &'static str,
    mfd_p: Option<Option<f32>>, // Some(None) = D8, Some(Some(p)) = MFD with exponent p
    a_c_km2: Option<f32>,
    diffusion: Option<f32>,
    /// ADR Finding 56 — the slope-dependent channel head.
    a_c_law: bool,
    /// Coarse reference: no FBM, no incision.
    reference: bool,
}

impl Lever {
    const fn shipped() -> Self {
        Lever {
            label: "SHIPPED",
            mfd_p: None,
            a_c_km2: None,
            diffusion: None,
            a_c_law: false,
            reference: false,
        }
    }
}

fn build(target: usize, lv: Lever) -> GridF32 {
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
    let cell_km2 = (DOMAIN_KM / target as f32).powi(2);
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: if lv.reference { 0.0 } else { 0.04 },
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: !lv.reference,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: !lv.reference,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    if lv.reference {
        cfg.stream_power = None;
    } else if let Some(sp) = cfg.stream_power.as_mut() {
        if let Some(p) = lv.mfd_p {
            sp.mfd_exponent = p;
        }
        if let Some(a) = lv.a_c_km2 {
            sp.min_area_cells = a / cell_km2;
        }
        if let Some(d) = lv.diffusion {
            sp.diffusion = d;
        }
        if lv.a_c_law {
            sp.a_c_slope_law = Some(ymir_core::erosion::stream_power::ChannelHeadLaw {
                s_ref: ymir_core::erosion::stream_power::CHANNEL_HEAD_S_REF,
                s_min: ymir_core::erosion::stream_power::CHANNEL_HEAD_S_MIN,
            });
        }
    }
    let (up, _) = upscale_from_c1_with_progress(
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
    );
    up.heightmap
}

fn slope_deg(g: &GridF32, x: usize, y: usize, m_per_cell: f32, ss: &SteinSteinParams) -> f32 {
    let at = |a: usize, b: usize| {
        c1_altitude_norm_to_metres(g.data[b.min(g.height - 1) * g.width + a.min(g.width - 1)], ss)
    };
    let gx = 0.5 * (at(x + 1, y) - at(x.saturating_sub(1), y));
    let gy = 0.5 * (at(x, y + 1) - at(x, y.saturating_sub(1)));
    ((gx * gx + gy * gy).sqrt() / m_per_cell).atan().to_degrees()
}

fn row(label: &str, s: &CoastShape) {
    let d = |ok: bool| if ok { ' ' } else { '!' };
    eprintln!(
        "{label:<26} {:>9.1}{} {:>8.2}{} | {:>8.2} {:>8.2} {:>8.0} {:>7.2} {:>7.3}",
        s.density_per_100km,
        d(s.density_per_100km < 5.0),
        s.p90_len_km,
        d(s.p90_len_km > 10.0),
        s.med_len_km,
        s.max_len_km,
        s.coast_km,
        s.spacing_cv,
        s.local_axis_r
    );
}

fn header() {
    eprintln!(
        "{:<26} {:>21} | {:>44}",
        "", "PRIMARY (optimised)", "REPORTED ONLY — never optimise these"
    );
    eprintln!(
        "{:<26} {:>10} {:>10} | {:>8} {:>8} {:>8} {:>7} {:>7}",
        "lever", "dens/100km", "p90 km", "med km", "max km", "coast km", "CV", "localR"
    );
    eprintln!(
        "{:<26} {:>10} {:>10} | {:>8} {:>8} {:>8} {:>7} {:>7}",
        "  TARGET (coarse ref)", "~1", "~24", "~3", "~46", "~1600", "~0.9", "~0.33"
    );
}

fn levers() -> Vec<Lever> {
    vec![
        Lever::shipped(),
        // MFD — the lever that cured the hillslope comb. Does it reach the apron?
        Lever { label: "MFD off (D8)", mfd_p: Some(None), ..Lever::shipped() },
        Lever { label: "MFD p = 1.1 (dispersed)", mfd_p: Some(Some(1.1)), ..Lever::shipped() },
        Lever { label: "MFD p = 4.0 (concentrated)", mfd_p: Some(Some(4.0)), ..Lever::shipped() },
        // A_c — does the comb spacing follow the channel-head threshold?
        Lever { label: "A_c = 0.01 km2 (x0.1)", a_c_km2: Some(0.01), ..Lever::shipped() },
        Lever { label: "A_c = 1.0 km2 (x10)", a_c_km2: Some(1.0), ..Lever::shipped() },
        Lever { label: "A_c = 10.0 km2 (x100)", a_c_km2: Some(10.0), ..Lever::shipped() },
        // Hillslope diffusion — currently too weak to matter (Finding 45b). Does strength help?
        Lever { label: "diffusion x10 (0.8)", diffusion: Some(0.8), ..Lever::shipped() },
        Lever { label: "diffusion x100 (8.0)", diffusion: Some(8.0), ..Lever::shipped() },
        Lever { label: "REFERENCE coarse", reference: true, ..Lever::shipped() },
    ]
}

fn run(target: usize, render: bool) {
    let ss = SteinSteinParams::default();
    let km = DOMAIN_KM / target as f32;
    let m_per_cell = km * 1000.0;
    eprintln!("\n==========  COASTAL COMB LEVERS — {target}² ({m_per_cell:.0} m/cell)  ==========");

    let shipped = build(target, Lever::shipped());
    let shipped_polys = marching_squares(&shipped, SEA);

    // ── APRON CHARACTERISATION, on the shipped config.
    {
        let (w, h) = (shipped.width, shipped.height);
        let mut slopes: Vec<f32> = Vec::new();
        let step = if w > 4096 { 4 } else { 2 };
        for y in (1..h - 1).step_by(step) {
            for x in (1..w - 1).step_by(step) {
                let v = shipped.data[y * w + x];
                if v <= SEA || c1_altitude_norm_to_metres(v, &ss) > APRON_M {
                    continue;
                }
                slopes.push(slope_deg(&shipped, x, y, m_per_cell, &ss));
            }
        }
        slopes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = slopes.len().max(1);
        let q = |f: f32| slopes[(((n - 1) as f64 * f as f64).round() as usize).min(n - 1)];
        // Apron WIDTH: its area divided by the coastline length is the mean inland distance.
        let apron_km2 = slopes.len() as f32 * (step * step) as f32 * km * km;
        let sh = coast_shape(&shipped_polys, km);
        // Comb spacing at the coast, from the spur density.
        let spacing_km = if sh.count > 0 { sh.coast_km / sh.count as f32 } else { 0.0 };
        // What A_c implies for channel-head spacing: sqrt of the critical area.
        let a_c_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
        eprintln!(
            "\nAPRON (land within {APRON_M:.0} m of sea level)\n  \
             area {apron_km2:.0} km² | mean width inland {:.2} km (area / coast length)\n  \
             slope p10 {:.2}° | p50 {:.2}° | p90 {:.2}° | share under 0.5°: {:.1} %\n\
             \nCOMB SPACING\n  \
             observed {spacing_km:.3} km between spurs\n  \
             cell size {km:.3} km  → ratio {:.1} cells\n  \
             sqrt(A_c) {:.3} km (A_c = {a_c_km2} km²) → ratio {:.2}",
            apron_km2 / sh.coast_km.max(1e-6),
            q(0.10),
            q(0.50),
            q(0.90),
            100.0 * slopes.iter().filter(|&&s| s < 0.5).count() as f32 / n as f32,
            spacing_km / km,
            a_c_km2.sqrt(),
            spacing_km / a_c_km2.sqrt(),
        );
    }

    // ── THE LEVER SWEEP.
    eprintln!();
    header();
    let (ox, oy) = densest_window(&shipped_polys, target, target, CROP);
    let out = std::path::Path::new("../../exports/coastal_fringes");
    if render {
        std::fs::create_dir_all(out).expect("out dir");
        eprintln!("{:<26} crop {CROP}×{CROP} at ({ox},{oy}), FIXED for every panel", "");
    }
    for lv in levers() {
        let field = if lv.label == "SHIPPED" { shipped.clone() } else { build(target, lv) };
        let polys = marching_squares(&field, SEA);
        let s = coast_shape(&polys, km);
        row(lv.label, &s);
        if render {
            let name = lv.label.replace([' ', '=', '.', '(', ')', '/'], "_");
            let _ = render_coast_crop(
                &field,
                &polys,
                SEA,
                ox,
                oy,
                CROP,
                &out.join(format!("lever_{name}_{target}.png")),
            );
        }
    }
    eprintln!(
        "\n'!' marks a PRIMARY metric outside the target band (density < 5 per 100 km, p90 > 10 km)."
    );
}

#[test]
#[ignore]
fn comb_levers_2048() {
    run(2048, false);
}

#[test]
#[ignore]
fn comb_levers_8192() {
    run(8192, true);
}

/// ADR rule 7 — the observables a coastline metric CANNOT see. `A_c × 100` renders identically
/// to the un-incised reference, which is either a cure or a switched-off erosion; only the
/// inland relief can tell them apart. Hypsometry in metres, fluvial network extent, and the
/// closed-depression population, for the levers that moved the primary metrics.
#[test]
#[ignore]
fn comb_levers_control_block() {
    let ss = SteinSteinParams::default();
    let target = 2048usize;
    let km = DOMAIN_KM / target as f32;
    let cell_km2 = km * km;
    eprintln!(
        "
=====  RULE-7 CONTROL BLOCK — what the coastline metric cannot see  =====

         A remedy that fixes the coast by disabling the erosion is not a remedy. The question is
         whether the INLAND relief survives.
"
    );
    eprintln!(
        "{:<26} {:>9} {:>9} {:>9} {:>12} {:>11}",
        "lever", "mean m", "p50 m", "p90 m", "channel km", "land %"
    );
    for lv in [
        Lever::shipped(),
        Lever { label: "A_c = 1.0 km2 (x10)", a_c_km2: Some(1.0), ..Lever::shipped() },
        Lever { label: "A_c = 10.0 km2 (x100)", a_c_km2: Some(10.0), ..Lever::shipped() },
        Lever { label: "diffusion x10 (0.8)", diffusion: Some(0.8), ..Lever::shipped() },
        Lever { label: "REFERENCE coarse", reference: true, ..Lever::shipped() },
    ] {
        let f = build(target, lv);
        let n = f.data.len();
        let mut alt: Vec<f32> = f
            .data
            .iter()
            .filter(|&&v| v > SEA)
            .map(|&v| c1_altitude_norm_to_metres(v, &ss))
            .collect();
        alt.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let m = alt.len().max(1);
        let q = |x: f32| alt[(((m - 1) as f64 * x as f64).round() as usize).min(m - 1)];
        // Fluvial network extent: cells whose accumulation clears the SHIPPED channel head,
        // measured the same way for every lever so the comparison is like-for-like.
        let flow = ymir_core::terrain::flow::compute_flow(
            &f,
            &ymir_core::terrain::flow::FlowConfig::default(),
        );
        let head_cells = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2 / cell_km2;
        let chan = flow.accumulation.data.iter().filter(|&&a| a >= head_cells).count();
        eprintln!(
            "{:<26} {:>9.0} {:>9.0} {:>9.0} {:>12.0} {:>10.1}%",
            lv.label,
            alt.iter().map(|&v| v as f64).sum::<f64>() / m as f64,
            q(0.5),
            q(0.9),
            chan as f32 * km,
            100.0 * m as f32 / n as f32
        );
    }
    eprintln!(
        "
⇒ If `A_c × 100` collapses the channel network and flattens the hypsometry towards the
           un-incised reference, it is not a cure — it is erosion switched off, and the uniform
           value must NOT ship. That would leave a SLOPE-DEPENDENT A_c as the only form that can
           raise the threshold on the apron while leaving the hillslopes calibrated."
    );
}

/// ADR Finding 56 — the LAW, judged on the three things together: primary metrics, the rule-7
/// control block proving the incision is alive, and a panel. None is sufficient alone: the
/// uniform `A_c × 100` passed the metrics AND the panel while having switched erosion off.
#[test]
#[ignore]
fn channel_head_law_verdict() {
    let ss = SteinSteinParams::default();
    for target in [2048usize, 8192] {
        let km = DOMAIN_KM / target as f32;
        eprintln!(
            "
==========  CHANNEL-HEAD LAW — {target}²  =========="
        );
        header();
        let mut fields = Vec::new();
        for lv in [
            Lever::shipped(),
            Lever { label: "A_c(S) LAW", a_c_law: true, ..Lever::shipped() },
            Lever { label: "REFERENCE coarse", reference: true, ..Lever::shipped() },
        ] {
            let f = build(target, lv);
            let polys = marching_squares(&f, SEA);
            row(lv.label, &coast_shape(&polys, km));
            fields.push((lv.label, f, polys));
        }
        // RULE-7 CONTROL BLOCK — the incision must still be alive.
        eprintln!(
            "
RULE-7 CONTROL BLOCK  (shipped 282/161/787 m · un-eroded reference 865/679/1836 m)"
        );
        eprintln!("{:<26} {:>9} {:>9} {:>9} {:>10}", "", "mean m", "p50 m", "p90 m", "land %");
        for (label, f, _) in &fields {
            let mut a: Vec<f32> = f
                .data
                .iter()
                .filter(|&&v| v > SEA)
                .map(|&v| c1_altitude_norm_to_metres(v, &ss))
                .collect();
            a.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            let m = a.len().max(1);
            let q = |p: f32| a[(((m - 1) as f64 * p as f64).round() as usize).min(m - 1)];
            eprintln!(
                "{label:<26} {:>9.0} {:>9.0} {:>9.0} {:>9.1}%",
                a.iter().map(|&v| v as f64).sum::<f64>() / m as f64,
                q(0.5),
                q(0.9),
                100.0 * m as f32 / f.data.len() as f32
            );
        }
        if target == 8192 {
            let out = std::path::Path::new("../../exports/coastal_fringes");
            let _ = std::fs::create_dir_all(out);
            let (ox, oy) = densest_window(&fields[0].2, target, target, CROP);
            for (label, f, polys) in &fields {
                let n = label.replace([' ', '(', ')', '/'], "_");
                let _ = render_coast_crop(
                    f,
                    polys,
                    SEA,
                    ox,
                    oy,
                    CROP,
                    &out.join(format!("law_{n}_8192.png")),
                );
            }
            eprintln!(
                "
panels: exports/coastal_fringes/law_*_8192.png (crop fixed at ({ox},{oy}))"
            );
        }
    }
}
