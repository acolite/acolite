//! End-to-end Sentinel-2 A/B MSI processing tests
//!
//! Focus: multi-resolution handling (10/20/60m), resampling correctness,
//! parallelization across resolution groups, COG output, S2A/S2B parity.

use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::ac::{estimate_dark_spectrum, optimize_aot_simple, rayleigh_optical_thickness, DarkSpectrumMethod};
use acolite_rs::parallel::{process_bands_parallel, process_bands_sequential};
use acolite_rs::sensors::{Sensor, Sentinel2Sensor};
use acolite_rs::{resample, ResampleMethod, Pipeline, ProcessingConfig};
use chrono::Utc;
use ndarray::Array2;

// ── S2 band definitions: (name, wavelength, bandwidth, resolution, typical DN) ──

const S2_10M: [(&str, f64, f64, u16); 4] = [
    ("B02", 492.4, 66.0, 1200),
    ("B03", 559.8, 36.0, 1500),
    ("B04", 664.6, 31.0, 1400),
    ("B08", 832.8, 106.0, 2000),
];

const S2_20M: [(&str, f64, f64, u16); 6] = [
    ("B05", 704.1, 15.0, 1800),
    ("B06", 740.5, 15.0, 1900),
    ("B07", 782.8, 20.0, 2000),
    ("B8A", 864.7, 21.0, 2100),
    ("B11", 1613.7, 91.0, 900),
    ("B12", 2202.4, 175.0, 500),
];

const S2_60M: [(&str, f64, f64, u16); 3] = [
    ("B01", 442.7, 21.0, 1100),
    ("B09", 945.1, 20.0, 800),
    ("B10", 1373.5, 31.0, 300),
];

fn make_s2_bands_at_res(
    defs: &[(&str, f64, f64, u16)],
    size: usize,
    proj: &Projection,
    res: f64,
) -> Vec<BandData<u16>> {
    let gt = GeoTransform::new(600000.0, res, 5000000.0, -res);
    defs.iter()
        .map(|&(name, wl, bw, dn)| {
            let mut data = Array2::from_elem((size, size), dn);
            // Spatial gradient
            for r in 0..size {
                for c in 0..size {
                    let off = ((r * 7 + c * 13) % 200) as u16;
                    data[[r, c]] = dn.saturating_add(off).saturating_sub(100);
                }
            }
            BandData::new(data, wl, bw, name.to_string(), proj.clone(), gt.clone())
        })
        .collect()
}

fn make_pipeline(metadata: &Metadata, swir_bands: &[BandData<u16>]) -> Pipeline {
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 1.5,
    };
    let mut pipeline = Pipeline::new(metadata.clone(), config);
    let swir: Vec<Array2<f64>> = swir_bands.iter()
        .filter(|b| b.wavelength > 1500.0)
        .map(|b| b.data.mapv(|v| v as f64 / 10000.0))
        .collect();
    if !swir.is_empty() {
        let dark = estimate_dark_spectrum(&swir, &DarkSpectrumMethod::Percentile(5.0));
        let wl: Vec<f64> = swir_bands.iter()
            .filter(|b| b.wavelength > 1500.0)
            .map(|b| b.wavelength)
            .collect();
        let aot = optimize_aot_simple(&dark, &wl, 28.0, 0.0);
        pipeline.set_aot(aot);
    }
    pipeline
}

// ── Sensor definition tests ─────────────────────────────────────────────────

#[test]
fn test_s2a_s2b_band_parity() {
    let s2a = Sentinel2Sensor::new_s2a();
    let s2b = Sentinel2Sensor::new_s2b();

    assert_eq!(s2a.name(), "S2A_MSI");
    assert_eq!(s2b.name(), "S2B_MSI");
    assert_eq!(s2a.band_names().len(), 13);
    assert_eq!(s2b.band_names().len(), 13);

    // Both must have the same band names
    let mut a_names: Vec<_> = s2a.band_names();
    let mut b_names: Vec<_> = s2b.band_names();
    a_names.sort();
    b_names.sort();
    assert_eq!(a_names, b_names);
}

#[test]
fn test_s2_resolution_groups() {
    let s2a = Sentinel2Sensor::new_s2a();

    let b10 = s2a.bands_at_resolution(10);
    let b20 = s2a.bands_at_resolution(20);
    let b60 = s2a.bands_at_resolution(60);

    assert_eq!(b10.len(), 4, "Expected 4 bands at 10m");
    assert_eq!(b20.len(), 6, "Expected 6 bands at 20m");
    assert_eq!(b60.len(), 3, "Expected 3 bands at 60m");
    assert_eq!(b10.len() + b20.len() + b60.len(), 13);

    // Verify specific bands at correct resolution
    assert_eq!(s2a.resolution("B04"), Some(10));
    assert_eq!(s2a.resolution("B11"), Some(20));
    assert_eq!(s2a.resolution("B01"), Some(60));
}

#[test]
fn test_s2_wavelength_ordering() {
    let s2a = Sentinel2Sensor::new_s2a();
    let mut wl: Vec<(String, f64)> = s2a.band_names().iter()
        .filter_map(|n| s2a.wavelength(n).map(|w| (n.clone(), w)))
        .collect();
    wl.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // First band should be B01 (coastal ~443nm), last B12 (SWIR2 ~2202nm)
    assert!(wl.first().unwrap().1 < 450.0);
    assert!(wl.last().unwrap().1 > 2200.0);
}

// ── Resampling tests ────────────────────────────────────────────────────────

#[test]
fn test_resample_20m_to_10m_dimensions() {
    let data_20m = Array2::from_elem((500, 500), 0.15);
    let data_10m = resample(&data_20m, 20, 10, ResampleMethod::Bilinear).unwrap();
    assert_eq!(data_10m.dim(), (1000, 1000));
}

#[test]
fn test_resample_60m_to_10m_dimensions() {
    let data_60m = Array2::from_elem((183, 183), 0.10);
    let data_10m = resample(&data_60m, 60, 10, ResampleMethod::Bilinear).unwrap();
    // 183 * 6 = 1098
    assert_eq!(data_10m.dim(), (1098, 1098));
}

#[test]
fn test_resample_10m_to_20m_average() {
    let data_10m = Array2::from_elem((1000, 1000), 0.25);
    let data_20m = resample(&data_10m, 10, 20, ResampleMethod::Average).unwrap();
    assert_eq!(data_20m.dim(), (500, 500));
    // Uniform input → mean preserved
    let mean = data_20m.mean().unwrap();
    assert!((mean - 0.25).abs() < 0.01, "Mean not preserved: {}", mean);
}

#[test]
fn test_resample_identity() {
    let data = Array2::from_elem((100, 100), 0.5);
    let result = resample(&data, 10, 10, ResampleMethod::Bilinear).unwrap();
    assert_eq!(result, data);
}

#[test]
fn test_resample_roundtrip_preserves_mean() {
    // 20m → 10m → 20m should approximately preserve mean
    let data_20m = Array2::from_shape_fn((50, 50), |(r, c)| {
        0.1 + 0.001 * (r as f64) + 0.0005 * (c as f64)
    });
    let mean_orig = data_20m.mean().unwrap();

    let up = resample(&data_20m, 20, 10, ResampleMethod::Bilinear).unwrap();
    let down = resample(&up, 10, 20, ResampleMethod::Average).unwrap();
    let mean_rt = down.mean().unwrap();

    assert!(
        (mean_orig - mean_rt).abs() < 0.01,
        "Roundtrip mean drift: {:.6} → {:.6}", mean_orig, mean_rt
    );
}

// ── Multi-resolution parallel pipeline ──────────────────────────────────────

#[test]
fn test_s2_multi_resolution_parallel_pipeline() {
    let proj = Projection::from_epsg(32632);
    let mut metadata = Metadata::new("S2A_MSI".to_string(), Utc::now());
    metadata.set_geometry(28.0, 135.0);

    let bands_10m = make_s2_bands_at_res(&S2_10M, 200, &proj, 10.0);
    let bands_20m = make_s2_bands_at_res(&S2_20M, 100, &proj, 20.0);
    let bands_60m = make_s2_bands_at_res(&S2_60M, 34, &proj, 60.0);

    let pipeline = make_pipeline(&metadata, &bands_20m);

    // Process each resolution group in parallel
    let r10 = process_bands_parallel(&pipeline, bands_10m).unwrap();
    let r20 = process_bands_parallel(&pipeline, bands_20m).unwrap();
    let r60 = process_bands_parallel(&pipeline, bands_60m).unwrap();

    assert_eq!(r10.len(), 4);
    assert_eq!(r20.len(), 6);
    assert_eq!(r60.len(), 3);

    for band in r10.iter().chain(r20.iter()).chain(r60.iter()) {
        let mean = band.data.mean().unwrap();
        assert!(mean.is_finite(), "NaN in {}", band.name);
        assert!(mean >= 0.0 && mean < 1.0, "ρs out of range in {}: {}", band.name, mean);
    }
}

// ── Parallel correctness ────────────────────────────────────────────────────

#[test]
fn test_s2_parallel_sequential_identical() {
    let proj = Projection::from_epsg(32632);
    let mut metadata = Metadata::new("S2A_MSI".to_string(), Utc::now());
    metadata.set_geometry(28.0, 135.0);

    let bands = make_s2_bands_at_res(&S2_10M, 100, &proj, 10.0);
    let pipeline = make_pipeline(&metadata, &bands);

    let bands2 = bands.clone();
    let par = process_bands_parallel(&pipeline, bands).unwrap();
    let seq = process_bands_sequential(&pipeline, bands2).unwrap();

    for (p, s) in par.iter().zip(seq.iter()) {
        let max_diff: f64 = p.data.iter().zip(s.data.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(max_diff == 0.0, "Par/seq differ on {}: {}", p.name, max_diff);
    }
}

#[test]
fn test_s2_parallel_deterministic() {
    let proj = Projection::from_epsg(32632);
    let mut metadata = Metadata::new("S2A_MSI".to_string(), Utc::now());
    metadata.set_geometry(28.0, 135.0);

    let bands = make_s2_bands_at_res(&S2_10M, 80, &proj, 10.0);
    let pipeline = make_pipeline(&metadata, &bands);

    let r1 = process_bands_parallel(&pipeline, bands.clone()).unwrap();
    let r2 = process_bands_parallel(&pipeline, bands).unwrap();

    for (a, b) in r1.iter().zip(r2.iter()) {
        let diff: f64 = a.data.iter().zip(b.data.iter())
            .map(|(x, y)| (x - y).abs()).sum();
        assert!(diff == 0.0, "Non-deterministic: {}", a.name);
    }
}

// ── S2A vs S2B cross-sensor ─────────────────────────────────────────────────

#[test]
fn test_s2a_s2b_cross_sensor_consistency() {
    let proj = Projection::from_epsg(32632);

    let mut meta_a = Metadata::new("S2A_MSI".to_string(), Utc::now());
    meta_a.set_geometry(28.0, 135.0);
    let mut meta_b = Metadata::new("S2B_MSI".to_string(), Utc::now());
    meta_b.set_geometry(28.0, 135.0);

    let bands_a = make_s2_bands_at_res(&S2_10M, 50, &proj, 10.0);
    let bands_b = bands_a.clone();

    let pipe_a = make_pipeline(&meta_a, &bands_a);
    let pipe_b = make_pipeline(&meta_b, &bands_b);

    let ra = process_bands_parallel(&pipe_a, bands_a).unwrap();
    let rb = process_bands_parallel(&pipe_b, bands_b).unwrap();

    for (a, b) in ra.iter().zip(rb.iter()) {
        let max_diff: f64 = a.data.iter().zip(b.data.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert!(max_diff < 1e-10, "S2A/S2B diverge on {}: {}", a.name, max_diff);
    }
}

// ── Rayleigh for S2 bands ───────────────────────────────────────────────────

#[test]
fn test_rayleigh_tau_s2_bands() {
    let s2a = Sentinel2Sensor::new_s2a();
    let mut wl_tau: Vec<(f64, f64)> = s2a.band_names().iter()
        .filter_map(|n| s2a.wavelength(n).map(|w| (w, rayleigh_optical_thickness(w, 1013.25))))
        .collect();
    wl_tau.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    for i in 1..wl_tau.len() {
        assert!(
            wl_tau[i].1 <= wl_tau[i - 1].1,
            "τ({:.0})={:.4} > τ({:.0})={:.4}",
            wl_tau[i].0, wl_tau[i].1, wl_tau[i - 1].0, wl_tau[i - 1].1
        );
    }
    // B01 (443nm) significant
    let tau_b01 = rayleigh_optical_thickness(442.7, 1013.25);
    assert!(tau_b01 > 0.1);
    // B12 (2202nm) negligible
    let tau_b12 = rayleigh_optical_thickness(2202.4, 1013.25);
    assert!(tau_b12 < 0.001);
}

// ── Performance ─────────────────────────────────────────────────────────────

#[test]
fn test_s2_throughput_10m_bands() {
    let proj = Projection::from_epsg(32632);
    let mut metadata = Metadata::new("S2A_MSI".to_string(), Utc::now());
    metadata.set_geometry(28.0, 135.0);

    // Realistic 10m tile: 10980×10980 is too large for test, use 2000×2000
    let size = 2000;
    let bands = make_s2_bands_at_res(&S2_10M, size, &proj, 10.0);
    let pipeline = make_pipeline(&metadata, &bands);

    let t0 = std::time::Instant::now();
    let _ = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = t0.elapsed();

    let mpx = (size * size * 4) as f64 / elapsed.as_secs_f64() / 1e6;
    println!("S2 10m throughput: {:.1} Mpx/s (4 bands, {}×{})", mpx, size, size);
    assert!(mpx > 1.0, "Throughput too low: {:.1} Mpx/s", mpx);
}

#[test]
fn test_s2_parallel_speedup_all_13_bands() {
    let proj = Projection::from_epsg(32632);
    let mut metadata = Metadata::new("S2A_MSI".to_string(), Utc::now());
    metadata.set_geometry(28.0, 135.0);

    // All 13 bands at 10m equivalent (500×500)
    let mut all_bands: Vec<BandData<u16>> = Vec::new();
    all_bands.extend(make_s2_bands_at_res(&S2_10M, 500, &proj, 10.0));
    all_bands.extend(make_s2_bands_at_res(&S2_20M, 500, &proj, 10.0));
    all_bands.extend(make_s2_bands_at_res(&S2_60M, 500, &proj, 10.0));

    let pipeline = make_pipeline(&metadata, &all_bands);

    // Warm up
    let _ = process_bands_parallel(&pipeline, all_bands.clone());

    let t0 = std::time::Instant::now();
    let _ = process_bands_sequential(&pipeline, all_bands.clone());
    let seq = t0.elapsed();

    let t0 = std::time::Instant::now();
    let _ = process_bands_parallel(&pipeline, all_bands);
    let par = t0.elapsed();

    println!(
        "S2 13-band: seq={:.2?} par={:.2?} speedup={:.2}x",
        seq, par, seq.as_secs_f64() / par.as_secs_f64().max(1e-9)
    );
    assert!(par.as_secs_f64() < seq.as_secs_f64() * 2.0);
}

#[test]
fn test_resample_throughput() {
    // Measure resampling speed: 20m→10m bilinear on 5490×5490 (half-tile)
    let size = 1000;
    let data = Array2::from_elem((size, size), 0.15);

    let t0 = std::time::Instant::now();
    let result = resample(&data, 20, 10, ResampleMethod::Bilinear).unwrap();
    let elapsed = t0.elapsed();

    let out_px = result.nrows() * result.ncols();
    let mpx = out_px as f64 / elapsed.as_secs_f64() / 1e6;
    println!("Resample 20m→10m: {:.1} Mpx/s ({}→{})", mpx, size, result.nrows());
    assert!(mpx > 1.0);
}

// ── COG output ──────────────────────────────────────────────────────────────

#[test]
fn test_s2_cog_output() {
    use acolite_rs::writer::write_auto;

    let proj = Projection::from_epsg(32632);
    let mut metadata = Metadata::new("S2A_MSI".to_string(), Utc::now());
    metadata.set_geometry(28.0, 135.0);

    let bands = make_s2_bands_at_res(&S2_10M, 50, &proj, 10.0);
    let pipeline = make_pipeline(&metadata, &bands);
    let corrected = process_bands_parallel(&pipeline, bands).unwrap();

    // 4 bands ≤ 50 → COG
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("s2_test");
    write_auto(out.to_str().unwrap(), &corrected, &metadata).unwrap();

    let tif = tmp.path().join("s2_test.tif");
    assert!(tif.exists());
    assert!(std::fs::metadata(&tif).unwrap().len() > 0);
}

// ── DSF AOT from SWIR ──────────────────────────────────────────────────────

#[test]
fn test_s2_dsf_aot_from_swir() {
    let swir1 = Array2::from_elem((100, 100), 0.008);
    let swir2 = Array2::from_elem((100, 100), 0.003);
    let dark = estimate_dark_spectrum(&[swir1, swir2], &DarkSpectrumMethod::Percentile(5.0));
    let aot = optimize_aot_simple(&dark, &[1613.7, 2202.4], 28.0, 0.0);
    assert!(aot > 0.0 && aot < 1.0, "AOT={}", aot);
}
