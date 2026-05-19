//! Calibration record — comparison of two `L_plate` normalisation
//! schemes for the cratonic factor, across an ensemble of Voronoï
//! seeds. Run via:
//!
//! ```text
//! cargo test --release -p ymir-core \
//!     --test v2_cratonic_normalization_probe -- --ignored --nocapture
//! ```
//!
//! The output table is the empirical evidence cited in the §4.10
//! patch and the `factor.rs` module docstring for the choice of
//! `L_plate = max BFS depth` over `L_plate = sqrt(plate_area)`.
//!
//! Scheme A (literal reading of issue text): `L_plate =
//! sqrt(plate_area)`, `d_mid = 0.5 (1 − sqrt(Cr))`. This was the
//! initial implementation and is preserved here as the comparison
//! baseline.
//!
//! Scheme B (geometric reading, current `factor.rs`):
//! `L_plate = max BFS depth in plate`, `d_mid = 1 − sqrt(Cr)`.
//!
//! Reports the cratonic_cell_fraction / (Cr · continental_fraction)
//! ratio (1.0 ideal) for each scheme, averaged over seeds. Step 9
//! locked in Scheme B with mean ratio ≈ 1.13 across 31 non-
//! degenerate seeds.

use std::collections::VecDeque;

use ymir_core::tectonics_v2::boundaries::PlateType;
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::voronoi::{VoronoiConfig, VoronoiPlates, generate_voronoi};

fn build_distances(plates: &VoronoiPlates, nx: usize, ny: usize, retained: &[bool]) -> Vec<u32> {
    let n = nx * ny;
    let mut dist: Vec<u32> = vec![u32::MAX; n];
    let mut queue: VecDeque<u32> = VecDeque::with_capacity(n);
    for j in 0..ny {
        for i in 0..nx {
            let pid = plates.plate_id.get(i, j) as usize;
            if !retained[pid] {
                let k = (j * nx + i) as u32;
                dist[k as usize] = 0;
                queue.push_back(k);
            }
        }
    }
    if queue.is_empty() {
        return vec![n as u32; n];
    }
    while let Some(k) = queue.pop_front() {
        let k = k as usize;
        let i = k % nx;
        let j = k / nx;
        let d = dist[k];
        let next = d.saturating_add(1);
        let neighbours = [
            (if i + 1 < nx { i + 1 } else { 0 }, j),
            (if i == 0 { nx - 1 } else { i - 1 }, j),
            (i, if j + 1 < ny { j + 1 } else { 0 }),
            (i, if j == 0 { ny - 1 } else { j - 1 }),
        ];
        for (ni, nj) in neighbours {
            let nk = nj * nx + ni;
            if dist[nk] > next {
                dist[nk] = next;
                queue.push_back(nk as u32);
            }
        }
    }
    dist
}

fn smoothstep(low: f64, high: f64, x: f64) -> f64 {
    if x <= low {
        return 0.0;
    }
    if x >= high {
        return 1.0;
    }
    let span = high - low;
    let t = (x - low) / span;
    t * t * (3.0 - 2.0 * t)
}

fn measure(seed: u64, cr: f64, scheme: char) -> (f64, f64, f64) {
    let nx = 64;
    let ny = 64;
    let vcfg = VoronoiConfig::default();
    let plates = generate_voronoi(nx, ny, &vcfg, seed);

    let mut areas = vec![0u32; plates.num_plates];
    for j in 0..ny {
        for i in 0..nx {
            areas[plates.plate_id.get(i, j) as usize] += 1;
        }
    }
    let domain_area = (nx * ny) as f64;
    let area_threshold = 0.10 * domain_area;
    let retained: Vec<bool> = (0..plates.num_plates)
        .map(|p| {
            matches!(plates.per_plate_type[p], PlateType::Continental)
                && (areas[p] as f64) >= area_threshold
        })
        .collect();
    let dist = build_distances(&plates, nx, ny, &retained);

    let plate_max_dist: Vec<u32> = {
        let mut v = vec![0u32; plates.num_plates];
        for j in 0..ny {
            for i in 0..nx {
                let pid = plates.plate_id.get(i, j) as usize;
                let d = dist[j * nx + i];
                if retained[pid] && d > v[pid] {
                    v[pid] = d;
                }
            }
        }
        v
    };

    let smoothing_width = 0.05;
    let factor = {
        let mut f = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let pid = plates.plate_id.get(i, j) as usize;
                if !retained[pid] {
                    continue;
                }
                let l_plate = match scheme {
                    'A' => (areas[pid] as f64).sqrt().max(1.0),
                    'B' => (plate_max_dist[pid] as f64).max(1.0),
                    _ => unreachable!(),
                };
                let d_norm = dist[j * nx + i] as f64 / l_plate;
                let d_mid = match scheme {
                    'A' => 0.5 * (1.0 - cr.sqrt()),
                    'B' => 1.0 - cr.sqrt(),
                    _ => unreachable!(),
                };
                let half = 0.5 * smoothing_width;
                let v = smoothstep(d_mid - half, d_mid + half, d_norm);
                f.set(i, j, v);
            }
        }
        f
    };

    let mut continental = 0usize;
    for j in 0..ny {
        for i in 0..nx {
            if matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                continental += 1;
            }
        }
    }
    let cratonic = factor.data().iter().filter(|&&v| v > 0.5).count();
    let total = (nx * ny) as f64;
    let cont_frac = continental as f64 / total;
    let crat_frac = cratonic as f64 / total;
    let expected = cr * cont_frac;
    let ratio = if expected > 0.0 { crat_frac / expected } else { 0.0 };
    (cont_frac, crat_frac, ratio)
}

#[test]
#[ignore]
fn compare_normalization_schemes() {
    let cr: f64 = 0.3;
    let mut seeds: Vec<u64> = (1..=30).collect();
    seeds.extend([42u64, 100, 123, 200, 256, 314, 500, 777, 1024, 9999]);
    println!();
    println!(
        "Scheme A: L_plate = sqrt(area), d_mid = 0.5(1 - sqrt(Cr)) = {:.4}",
        0.5 * (1.0 - cr.sqrt())
    );
    println!("Scheme B: L_plate = max BFS dist, d_mid = 1 - sqrt(Cr)    = {:.4}", 1.0 - cr.sqrt());
    println!();
    println!(
        "{:>5} | {:>10} | {:>10} {:>10} {:>10} | {:>10} {:>10} {:>10}",
        "seed", "cont_frac", "A:crat", "A:expect", "A:ratio", "B:crat", "B:expect", "B:ratio"
    );
    let mut sum_a = 0.0;
    let mut sum_b = 0.0;
    let mut na = 0usize;
    let mut nb = 0usize;
    for seed in &seeds {
        let (cf, ca, ra) = measure(*seed, cr, 'A');
        let (_cf, cb, rb) = measure(*seed, cr, 'B');
        let expected = cr * cf;
        println!(
            "{:>5} | {:>10.4} | {:>10.4} {:>10.4} {:>10.3} | {:>10.4} {:>10.4} {:>10.3}",
            seed, cf, ca, expected, ra, cb, expected, rb
        );
        if ra > 0.0 {
            sum_a += ra;
            na += 1;
        }
        if rb > 0.0 {
            sum_b += rb;
            nb += 1;
        }
    }
    println!();
    if na > 0 {
        println!("Scheme A mean ratio over {} non-degenerate seeds = {:.3}", na, sum_a / na as f64);
    }
    if nb > 0 {
        println!("Scheme B mean ratio over {} non-degenerate seeds = {:.3}", nb, sum_b / nb as f64);
    }
}
