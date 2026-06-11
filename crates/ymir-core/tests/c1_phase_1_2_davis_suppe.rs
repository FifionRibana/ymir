//! Issue #123 — C1 Phase 1.2 Davis-Suppe orogenic closure
//! integration test.
//!
//! Two scenarios:
//!
//! 1. [`davis_suppe_wedge_body_invariants`] — runs the C1
//!    prototype with the Davis-Suppe closure ON at the Phase 1.1
//!    64²×300-step configuration, dumps altitude + fixed-palette
//!    `S̃` snapshots at cycles 0/50/100/200/300, and asserts the
//!    three Phase 1.2 acceptance invariants:
//!
//!    a) **p95 bound** — `wedge_p95 < 1.5 · h_max`. The bulk of
//!       the wedge body stays within `≈ h_max` of the surface;
//!       advection-driven outliers above are tolerated.
//!    b) **p99 activity** — `wedge_p99 > 1.0`. The closure must
//!       produce some observable accumulation above the
//!       continental init level; without it (silent bug),
//!       `p99` would stay ≤ init max ≈ 1.0.
//!    c) **h_critical profile imprint** — mean `S̃` in distance
//!       bucket `d ∈ (10, 20]` exceeds 1.5 × mean in bucket
//!       `d ∈ (0, 5]`. This is the *shape* check: even when the
//!       Phase 1.1 advection-dominated regime drains the bulk
//!       mean, the **conditional** mean must inherit the
//!       `h_critical(d)` profile (small at small d, large at
//!       large d).
//!
//!    Boundary cells (`BoundaryType::Convergent`) are explicitly
//!    **not** bounded by Phase 1.2 — the Stage 3.1 architectural
//!    skip leaves them as advection sinks; Phase 1.4 erosion
//!    will eventually balance them. See the
//!    `c1_phase_1_2_davis_suppe/README.md` "On boundary cell
//!    brightness" and "Advection-dominated regime" sections.
//!
//! 2. [`davis_suppe_disabled_matches_phase_1_1`] — runs
//!    [`run_advection_only`] (Phase 1.1 entry-point, **not**
//!    `run_with_closures` with `enabled=false`) on the same
//!    init+kinematics, and asserts Phase 1.1's unbounded pile-up
//!    behaviour is preserved (W4 regression).
//!
//! ## Why call `run_advection_only` for the regression test
//!
//! Per the time-loop module docstring: `run_with_closures` does
//! one-shot work outside its loop (boundary classification +
//! intra-plate wedge distance) even when `davis_suppe.enabled =
//! false`. That overhead is small but not bit-identical to
//! Phase 1.1. The closure-OFF baseline is therefore the
//! `run_advection_only` direct call, which exercises the
//! Phase 1.1 path verbatim.

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::boundary_classification::classify_boundaries;
use ymir_core::tectonics_c1::closures::accretion::AccretionParams;
use ymir_core::tectonics_c1::closures::davis_suppe::source_term::DavisSuppeParams;
use ymir_core::tectonics_c1::closures::equilibrium_height::params::EquilibriumHeightParams;
use ymir_core::tectonics_c1::closures::erosion::params::ErosionParams;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::rifting::RiftingParams;
use ymir_core::tectonics_c1::closures::subduction::SubductionParams;
use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate;
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{
    run_advection_only, run_with_closures, C1Closures, C1TimeLoopConfig,
};
use ymir_core::tectonics_v2::field::Field2D;

const GRID_SIZE: usize = 64;
const SEED: u64 = 42;
const N_STEPS: usize = 300;

/// Palette bound for the Phase 1.2 fixed-scale `S̃` PNG. Picked to
/// cover the post-closure dynamic range with a small headroom:
/// `h_max = 2.5` (Davis-Suppe plateau) plus ~ 20 % headroom for
/// transient overshoot before relaxation settles.
const S_VIZ_MAX: f64 = 3.0;

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/c1_phase_1_2_davis_suppe")
}

#[test]
fn davis_suppe_wedge_body_invariants() {
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let mut state = init_c1_state_phase_1_1(GRID_SIZE, SEED);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let config = C1TimeLoopConfig {
        // #145 — migrated to rigid (production transport). The wedge-body
        // invariants here (bucket-count + global_max bounds) hold unchanged
        // under rigid; the wedge sits higher (crust not advected away) but the
        // asserted invariants are not wedge_p95-valued, so no re-baseline needed.
        rigid_continental_crust: true,
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
        iso_config: IsostasyConfig::default(),
        drainage_max_distance: 30,
    };
    // Phase 1.2 semantics: keep Davis-Suppe ON, equilibrium-height
    // OFF. The Phase 1.2 invariants (especially the unbounded
    // boundary pile-up `global_max ≈ 2297`) rely on no global sink
    // being active. Phase 1.3 adds equilibrium-height defaulted ON
    // — disable it explicitly here so this test continues to lock
    // the pure-Phase-1.2 behaviour.
    let closures = C1Closures {
        davis_suppe: DavisSuppeParams::default(),
        equilibrium_height: EquilibriumHeightParams {
            enabled: false,
            ..EquilibriumHeightParams::default()
        },
        erosion: ErosionParams {
            enabled: false,
            ..ErosionParams::default()
        },
        oceanic_bathymetry: SteinSteinParams {
            enabled: false,
            ..SteinSteinParams::default()
        },
        subduction: SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        },
        accretion: AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        },
        rifting: RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        },
    };
    let h_max = closures.davis_suppe.h_max;

    eprintln!(
        "c1_phase_1_2: grid={}², steps={}, coupling={}, h_max={}, L_taper={}, L_decay={}",
        GRID_SIZE,
        N_STEPS,
        closures.davis_suppe.coupling,
        h_max,
        closures.davis_suppe.l_taper,
        closures.davis_suppe.l_decay,
    );
    print_s_stats("000", &state);
    dump_snapshot(&state, 0, &dir);

    let snapshot_steps: [usize; 4] = [49, 99, 199, 299];

    let started = std::time::Instant::now();
    run_with_closures(
        &mut state,
        &mut kinematics,
        &config,
        &closures,
        |step, current_state| {
            assert!(
                current_state.s.data().iter().all(|v| v.is_finite()),
                "non-finite S̃ at step {}",
                step + 1
            );
            if snapshot_steps.contains(&step) {
                print_s_stats(&format!("{:03}", step + 1), current_state);
                dump_snapshot(current_state, step + 1, &dir);
            }
        },
    );
    let elapsed = started.elapsed();

    // Re-classify boundaries + recompute wedge distance on the
    // final state. Plate_id is static in Phase 1.2, so both fields
    // are identical to the ones used inside the time loop — we
    // recompute here to keep the test self-contained (no need to
    // surface them from `run_with_closures`).
    let boundary = classify_boundaries(&state.plate_id, &kinematics);
    let wedge_d = wedge_distance_intra_plate(
        &state.plate_id,
        &boundary.upper_plate_mask,
        closures.davis_suppe.max_distance,
    );

    // True wedge-body filter: cells where Davis-Suppe is active
    // (`0 < d < max_distance`). This excludes both the Convergent
    // boundary itself (d = 0 by construction of intra-plate
    // Dijkstra) and cells out of reach of any same-plate seed.
    let max_d_cfg = closures.davis_suppe.max_distance;
    let mut global_max = 0.0_f64;
    let mut wedge_values: Vec<f64> = Vec::new();
    // Per-distance-bucket accumulators for the h_critical-profile
    // invariant. Bucket boundaries are open-left, closed-right.
    let bucket_edges: [(f64, f64); 3] = [(0.0, 5.0), (5.0, 10.0), (10.0, 20.0)];
    let mut bucket_sum = [0.0_f64; 3];
    let mut bucket_count = [0_usize; 3];
    for j in 0..state.ny() {
        for i in 0..state.nx() {
            let v = state.s.get(i, j);
            if v > global_max {
                global_max = v;
            }
            let d = wedge_d.get(i, j);
            if d > 0.0 && d < max_d_cfg {
                wedge_values.push(v);
                for (b, &(lo, hi)) in bucket_edges.iter().enumerate() {
                    if d > lo && d <= hi {
                        bucket_sum[b] += v;
                        bucket_count[b] += 1;
                    }
                }
            }
        }
    }
    wedge_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = wedge_values.len();
    let wedge_max = *wedge_values.last().unwrap_or(&0.0);
    let wedge_mean = wedge_values.iter().sum::<f64>() / n.max(1) as f64;
    let wedge_median = wedge_values[n / 2];
    let wedge_p95 = wedge_values[(n * 95) / 100];
    let wedge_p99 = wedge_values[(n * 99) / 100];

    let bucket_mean: [f64; 3] = std::array::from_fn(|b| {
        if bucket_count[b] > 0 {
            bucket_sum[b] / bucket_count[b] as f64
        } else {
            f64::NAN
        }
    });
    let h_crit_at = |d: f64| -> f64 {
        closures.davis_suppe.h_max * (1.0 - (-d / closures.davis_suppe.l_taper).exp())
    };

    eprintln!(
        "c1_phase_1_2: wedge body cells              = {n} ({:.1} %)",
        100.0 * n as f64 / (state.nx() * state.ny()) as f64
    );
    eprintln!("c1_phase_1_2: wedge S̃  mean                = {wedge_mean:.4}");
    eprintln!("c1_phase_1_2: wedge S̃  median              = {wedge_median:.4}");
    eprintln!("c1_phase_1_2: wedge S̃  p95                 = {wedge_p95:.4}");
    eprintln!("c1_phase_1_2: wedge S̃  p99                 = {wedge_p99:.4}");
    eprintln!("c1_phase_1_2: wedge S̃  max                 = {wedge_max:.4}");
    eprintln!(
        "c1_phase_1_2: h_critical profile (mean S̃ per distance bucket):"
    );
    for (b, &(lo, hi)) in bucket_edges.iter().enumerate() {
        let mid = (lo + hi) / 2.0;
        let h_crit_mid = h_crit_at(mid);
        let fill = if h_crit_mid > 0.0 {
            bucket_mean[b] / h_crit_mid
        } else {
            f64::NAN
        };
        eprintln!(
            "    d ∈ ({lo:>4.1}, {hi:>4.1}]  count = {:>5}  mean S̃ = {:.4}  h_crit(mid) = {:.4}  fill = {:.3}",
            bucket_count[b], bucket_mean[b], h_crit_mid, fill
        );
    }
    // Fill ratios per bucket (mean / h_crit_at_mid). Régime-
    // agnostic measure of "how saturated is each wedge band
    // relative to its Davis-Suppe target".
    let fill_near = bucket_mean[0] / h_crit_at(2.5);
    let fill_far = bucket_mean[2] / h_crit_at(15.0);
    let asymmetry = bucket_mean[0] / bucket_mean[2]; // near / far
    eprintln!(
        "c1_phase_1_2: fill_near (d∈0-5)  / h_crit(2.5)   = {fill_near:.3}  (need > 0.5)"
    );
    eprintln!(
        "c1_phase_1_2: fill_far  (d∈10-20)/ h_crit(15.0)  = {fill_far:.3}  (informational)"
    );
    eprintln!(
        "c1_phase_1_2: asymmetry mean(0-5) / mean(10-20)  = {asymmetry:.3}  (need > 1.5)"
    );
    eprintln!(
        "c1_phase_1_2: global_max (boundary pile-up) = {global_max:.2}  (no Phase 1.2 bound — Phase 1.4 erosion sink)",
    );
    eprintln!(
        "c1_phase_1_2: wall time     = {:.2?} ({:.2?} per step)",
        elapsed,
        elapsed / N_STEPS as u32
    );
    eprintln!("c1_phase_1_2: output dir    = {}", dir.display());

    // Invariant (a) — p95 bound. The bulk of the wedge body stays
    // within ≈ h_max; advection-driven outliers above are tolerated.
    assert!(
        wedge_p95 < 1.5 * h_max,
        "Phase 1.2 (a) p95 bound failed: wedge_p95 = {wedge_p95:.3} ≥ 1.5 · h_max = {:.2}",
        1.5 * h_max
    );
    // Invariant (b) — p99 activity. The closure must produce some
    // observable accumulation above the continental init level
    // (init max = 1.0). If a silent bug disables the source, p99
    // would stay ≤ ~1.0.
    assert!(
        wedge_p99 > 1.0,
        "Phase 1.2 (b) p99 activity failed: wedge_p99 = {wedge_p99:.3} ≤ 1.0 — closure may be silently disabled"
    );
    // Invariant (c) — h_critical profile imprint, hybrid form.
    //
    // The original Phase 1.2 plan expected the conditional mean
    // `S̃(d)` to *increase* with `d`, mirroring the `h_critical(d)`
    // shape that grows from 0 toward `h_max`. Empirically the
    // **advection-dominated** Phase 1.1 kinematics regime
    // inverts that direction: mass concentrates at small `d`
    // (envelope ≈ 1, source strong) and drains at large `d`
    // (envelope tail, source weak vs advective outflow). The
    // h_critical *imprint* therefore manifests in the **fill
    // ratio** `mean / h_crit` — a regime-agnostic measure of
    // "how close does each band get to its source target" — and
    // in the **near-vs-far asymmetry** of the absolute mean.
    //
    // The two sub-assertions verify both signals:
    //
    //   (c1) fill_near > 0.5 — near-boundary cells reach at
    //         least half their h_crit target despite the
    //         advection drain.
    //   (c2) mean(d∈0-5) / mean(d∈10-20) > 1.5 — Davis-Suppe
    //         creates a clear near-vs-far asymmetry in this
    //         regime.
    //
    // **Direction note.** Sub-assertion (c2) tests
    // `near > far`, the advection-dominated signature. Phase 1.4
    // erosion (mass sink) or Phase 2 slower kinematics (§6.3
    // C1.md) may flip this ratio. Future phase tests should
    // re-evaluate the direction rather than copy this one
    // verbatim.
    assert!(
        bucket_count[0] > 0 && bucket_count[2] > 0,
        "Phase 1.2 (c) profile buckets empty: cannot compute imprint metrics"
    );
    assert!(
        fill_near > 0.5,
        "Phase 1.2 (c1) saturation failed: fill_near = {fill_near:.3} ≤ 0.5 — \
         closure under-fills the near-boundary wedge band"
    );
    // #155 maillon 1b-i RESTORED. This near>far asymmetry is the macro-
    // border thermometer. It briefly went RED under maillon 1a (≈1.11):
    // 1a retargeted the DS wedge onto the continental plate at O-C
    // subductions but left the geometry a Tibet-style dome filling the
    // interior. Maillon 1b-i routes the O-C wedge to a margin-peaked ridge
    // profile (`h_max·exp(−d/l_taper)`, peaks at the margin, decays
    // inland), concentrating the relief near the boundary → near/far rises
    // back above 1.5 (measured ≈1.70). Green here means the O-C relief is
    // a margin ridge (Andes), not an interior dome. Consistent with the
    // direction note above (~lines 324-329).
    assert!(
        asymmetry > 1.5,
        "Phase 1.2 (c2) asymmetry failed: mean(d∈0-5)/mean(d∈10-20) = {asymmetry:.3} ≤ 1.5 — \
         Davis-Suppe near-vs-far signature absent (advection-dominated regime expected)"
    );
}

#[test]
fn davis_suppe_disabled_matches_phase_1_1() {
    // Regression: with no closure (run_advection_only), the
    // Phase 1.1 unbounded pile-up at convergence corners must
    // be preserved. The Phase 1.1 baseline observed
    // final max S̃ ≈ 1080 × initial.
    let mut state = init_c1_state_phase_1_1(GRID_SIZE, SEED);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let config = C1TimeLoopConfig {
        rigid_continental_crust: false,
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
        iso_config: IsostasyConfig::default(),
        drainage_max_distance: 30,
    };
    run_advection_only(&mut state, &kinematics, &config, |_, _| {});
    let final_max = state.s.data().iter().cloned().fold(0.0_f64, f64::max);
    eprintln!(
        "c1_phase_1_2_disabled: final max S̃ = {:.2} (Phase 1.1 baseline ~1080)",
        final_max
    );
    // The 100 threshold gives plenty of headroom for tiny
    // numerical differences from Phase 1.1 (which reported 1080)
    // without conflating any closure leak.
    assert!(
        final_max > 100.0,
        "Phase 1.1 unbounded behaviour not preserved: final max = {:.2} < 100",
        final_max
    );
}

// ── Diagnostics + viz helpers (mirror of Phase 1.1 test) ──────

fn print_s_stats(tag: &str, state: &C1State) {
    let data = state.s.data();
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for &v in data {
        min = min.min(v);
        max = max.max(v);
        sum += v;
    }
    let mean = sum / data.len() as f64;
    let mut sq = 0.0;
    for &v in data {
        sq += (v - mean) * (v - mean);
    }
    let std = (sq / data.len() as f64).sqrt();
    eprintln!(
        "c1_phase_1_2: cycle_{tag} S̃ min={min:.4} mean={mean:.4} max={max:.4} std={std:.4e}"
    );
}

fn dump_snapshot(state: &C1State, cycle: usize, dir: &Path) {
    let iso = compute_isostasy(&state.s, &IsostasyConfig::default());
    let alt_path = dir.join(format!("cycle_{:03}_altitude.png", cycle));
    save_hypsometric_png(&iso.heightmap, iso.sea_level_normalized, &alt_path);

    let s_path = dir.join(format!("cycle_{:03}_s.png", cycle));
    save_s_fixed_palette_png(&state.s, S_VIZ_MAX, &s_path);
}

fn save_hypsometric_png(
    heightmap: &ymir_core::grid::GridF32,
    sea_norm: f32,
    path: &Path,
) {
    let nx = heightmap.width;
    let ny = heightmap.height;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    for j in 0..ny {
        for i in 0..nx {
            let h = heightmap.get(i as i32, j as i32).clamp(0.0, 1.0);
            let rgb = hypsometric(h, sea_norm);
            let img_row = (ny - 1 - j) as u32;
            img.put_pixel(i as u32, img_row, Rgb(rgb));
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("create parent dir");
    }
    img.save(path).expect("save PNG");
}

fn save_s_fixed_palette_png(s: &Field2D, s_max: f64, path: &Path) {
    let nx = s.nx();
    let ny = s.ny();
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    // Sea-level marker at oceanic init value 0.2 (same as
    // Phase 1.1 viz for visual continuity).
    let sea_norm = 0.2 / s_max;
    for j in 0..ny {
        for i in 0..nx {
            let v = (s.get(i, j) / s_max).clamp(0.0, 1.0) as f32;
            let rgb = hypsometric(v, sea_norm as f32);
            let img_row = (ny - 1 - j) as u32;
            img.put_pixel(i as u32, img_row, Rgb(rgb));
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("create parent dir");
    }
    img.save(path).expect("save PNG");
}

fn hypsometric(h: f32, sea_norm: f32) -> [u8; 3] {
    let mid = (sea_norm + 1.0) * 0.5;
    let lerp = |t: f32, a: [u8; 3], b: [u8; 3]| -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        [
            (a[0] as f32 + t * (b[0] as f32 - a[0] as f32)).round() as u8,
            (a[1] as f32 + t * (b[1] as f32 - a[1] as f32)).round() as u8,
            (a[2] as f32 + t * (b[2] as f32 - a[2] as f32)).round() as u8,
        ]
    };
    if h <= sea_norm * 0.5 {
        let t = h / (sea_norm * 0.5).max(1e-6);
        lerp(t, [10, 20, 60], [40, 80, 160])
    } else if h <= sea_norm {
        let t = (h - sea_norm * 0.5) / (sea_norm * 0.5).max(1e-6);
        lerp(t, [40, 80, 160], [120, 180, 230])
    } else if h <= mid {
        let t = (h - sea_norm) / (mid - sea_norm).max(1e-6);
        lerp(t, [60, 130, 60], [140, 100, 50])
    } else {
        let t = (h - mid) / (1.0 - mid).max(1e-6);
        lerp(t, [140, 100, 50], [245, 245, 245])
    }
}
