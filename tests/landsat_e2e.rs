//! End-to-end Landsat 8/9 OLI processing tests
//!
//! Focus: parallelization correctness, performance scaling, COG output, and
//! cross-sensor consistency between L8 and L9.

use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::ac::{
    estimate_dark_spectrum, optimize_aot_simple, rayleigh_optical_thickness,
    earth_sun_distance, dn_to_radiance, dn_to_reflectance, DarkSpectrumMethod,
};
use acolite_rs::parallel::{process_bands_parallel, process_bands_sequential};
use acolite_rs::sensors::{LandsatSensor, Sensor};
use acolite_rs::{Pipeline, ProcessingConfig};
use chrono::Utc;
use ndarray::Array2;

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Landsat OLI band definitions: (name, wavelength, bandwidth, typical DN)
const L8_BANDS: [(&str, f64, f64, u16); 7] = [
    ("B1", 443.0, 16.0, 1200),   // Coastal
    ("B2", 482.0, 60.0, 1400),   // Blue
    ("B3", 561.0, 57.0, 1600),   // Green
    ("B4", 655.0, 37.0, 1500),   // Red
    ("B5", 865.0, 28.0, 2000),   // NIR
    ("B6", 1609.0, 85.0, 1000),  // SWIR1
    ("B7", 2201.0, 187.0, 600),  // SWIR2
];

fn make_landsat_bands(
    sensor_name: &str,
    size: usize,
) -> (Metadata, Vec<BandData<u16>>) {
    let mut metadata = Metadata::new(sensor_name.to_string(), Utc::now());
    metadata.set_geometry(35.0, 140.0); // typical mid-latitude
    metadata.view_zenith = Some(0.0);   // nadir

    let proj = Projection::from_epsg(32610); // UTM 10N
    let gt = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);

    let bands = L8_BANDS
        .iter()
        .map(|&(name, wl, bw, dn)| {
            // Add spatial gradient to simulate real scene
            let mut data = Array2::from_elem((size, size), dn);
            for r in 0..size {
                for c in 0..size {
                    let offset = ((r * 7 + c * 13) % 500) as u16;
                    data[[r, c]] = dn.saturating_add(offset).saturating_sub(250);
                }
            }
            BandData::new(data, wl, bw, name.to_string(), proj.clone(), gt.clone())
        })
        .collect();

    (metadata, bands)
}

fn make_pipeline(metadata: Metadata, bands: &[BandData<u16>]) -> Pipeline {
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 1.5,
        ..ProcessingConfig::default()
    };
    let mut pipeline = Pipeline::new(metadata, config);

    // DSF AOT from SWIR dark spectrum
    let swir: Vec<Array2<f64>> = bands[5..7]
        .iter()
        .map(|b| b.data.mapv(|v| v as f64 / 10000.0))
        .collect();
    let dark = estimate_dark_spectrum(&swir, &DarkSpectrumMethod::Percentile(5.0));
    let aot = optimize_aot_simple(&dark, &[1609.0, 2201.0], 35.0, 0.0);
    pipeline.set_aot(aot);
    pipeline
}

// ── Sensor definition tests ─────────────────────────────────────────────────

#[test]
fn test_l8_l9_sensor_band_parity() {
    let l8 = LandsatSensor::new_l8();
    let l9 = LandsatSensor::new_l9();

    assert_eq!(l8.name(), "L8_OLI");
    assert_eq!(l9.name(), "L9_OLI");

    // L8 and L9 OLI share identical band definitions
    for name in &["B1", "B2", "B3", "B4", "B5", "B6", "B7"] {
        let wl8 = l8.wavelength(name).expect(name);
        let wl9 = l9.wavelength(name).expect(name);
        assert_eq!(wl8, wl9, "Wavelength mismatch for {}", name);

        let bw8 = l8.bandwidth(name).expect(name);
        let bw9 = l9.bandwidth(name).expect(name);
        assert_eq!(bw8, bw9, "Bandwidth mismatch for {}", name);
    }
    assert_eq!(l8.band_names().len(), 7);
    assert_eq!(l9.band_names().len(), 7);
}

// ── Calibration tests ───────────────────────────────────────────────────────

#[test]
fn test_dn_to_radiance_linearity() {
    let dn = Array2::from_shape_fn((10, 10), |(r, c)| (r * 10 + c) as u16 * 100);
    let mult = 0.012;
    let add = -60.0;
    let rad = dn_to_radiance(&dn, mult, add);

    // Verify linearity: rad = dn * mult + add
    for r in 0..10 {
        for c in 0..10 {
            let expected = dn[[r, c]] as f64 * mult + add;
            assert!((rad[[r, c]] - expected).abs() < 1e-10);
        }
    }
}

#[test]
fn test_dn_to_reflectance_sun_angle() {
    let dn = Array2::from_elem((5, 5), 10000u16);
    let mult = 2e-5;
    let add = -0.1;

    // Higher sun elevation → lower reflectance (more energy)
    let rho_high = dn_to_reflectance(&dn, mult, add, 60.0); // high sun
    let rho_low = dn_to_reflectance(&dn, mult, add, 30.0);  // low sun
    assert!(
        rho_high[[0, 0]] < rho_low[[0, 0]],
        "Higher sun should give lower reflectance"
    );
}

#[test]
fn test_earth_sun_distance_annual_cycle() {
    // Perihelion ~Jan 3 (doy ~3), aphelion ~Jul 4 (doy ~185)
    let d_peri = earth_sun_distance(3);
    let d_aph = earth_sun_distance(185);
    assert!(d_peri < d_aph, "Perihelion distance should be less than aphelion");
    assert!((d_peri - 0.983).abs() < 0.02);
    assert!((d_aph - 1.017).abs() < 0.02);
}

// ── Parallel vs sequential correctness ──────────────────────────────────────

#[test]
fn test_landsat_parallel_sequential_identical() {
    let (metadata, bands) = make_landsat_bands("L8_OLI", 100);
    let pipeline = make_pipeline(metadata, &bands);

    let bands2 = bands.clone();
    let par = process_bands_parallel(&pipeline, bands).unwrap();
    let seq = process_bands_sequential(&pipeline, bands2).unwrap();

    assert_eq!(par.len(), seq.len());
    for (p, s) in par.iter().zip(seq.iter()) {
        assert_eq!(p.wavelength, s.wavelength);
        let max_diff: f64 = p.data.iter().zip(s.data.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(max_diff == 0.0, "Band {} diff={}", p.name, max_diff);
    }
}

#[test]
fn test_landsat_parallel_deterministic() {
    let (metadata, bands) = make_landsat_bands("L9_OLI", 80);
    let pipeline = make_pipeline(metadata, &bands);

    let r1 = process_bands_parallel(&pipeline, bands.clone()).unwrap();
    let r2 = process_bands_parallel(&pipeline, bands).unwrap();

    for (a, b) in r1.iter().zip(r2.iter()) {
        let diff: f64 = a.data.iter().zip(b.data.iter())
            .map(|(x, y)| (x - y).abs()).sum();
        assert!(diff == 0.0, "Non-deterministic: band {}", a.name);
    }
}

// ── Full pipeline correctness ───────────────────────────────────────────────

#[test]
fn test_landsat8_e2e_pipeline() {
    let (metadata, bands) = make_landsat_bands("L8_OLI", 200);
    let pipeline = make_pipeline(metadata, &bands);

    let corrected = process_bands_parallel(&pipeline, bands).unwrap();
    assert_eq!(corrected.len(), 7);

    for band in &corrected {
        assert_eq!(band.shape(), (200, 200));
        let mean = band.data.mean().unwrap();
        assert!(mean.is_finite(), "NaN in {}", band.name);
        assert!(mean >= 0.0, "Negative ρs in {}: {}", band.name, mean);
        assert!(mean < 1.0, "ρs > 1 in {}: {}", band.name, mean);
    }

    // Spectral ordering preserved
    for i in 1..corrected.len() {
        assert!(corrected[i].wavelength > corrected[i - 1].wavelength);
    }
}

#[test]
fn test_landsat9_e2e_pipeline() {
    let (metadata, bands) = make_landsat_bands("L9_OLI", 200);
    let pipeline = make_pipeline(metadata, &bands);

    let corrected = process_bands_parallel(&pipeline, bands).unwrap();
    assert_eq!(corrected.len(), 7);

    for band in &corrected {
        let mean = band.data.mean().unwrap();
        assert!(mean.is_finite() && mean >= 0.0 && mean < 1.0);
    }
}

#[test]
fn test_l8_l9_cross_sensor_consistency() {
    // Same input DN → same output (sensors share calibration in this simplified model)
    let (meta8, bands8) = make_landsat_bands("L8_OLI", 50);
    let (meta9, bands9) = make_landsat_bands("L9_OLI", 50);

    let pipe8 = make_pipeline(meta8, &bands8);
    let pipe9 = make_pipeline(meta9, &bands9);

    let r8 = process_bands_parallel(&pipe8, bands8).unwrap();
    let r9 = process_bands_parallel(&pipe9, bands9).unwrap();

    for (a, b) in r8.iter().zip(r9.iter()) {
        assert_eq!(a.wavelength, b.wavelength);
        let max_diff: f64 = a.data.iter().zip(b.data.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        // Same input + same physics → identical output
        assert!(max_diff < 1e-10, "L8/L9 diverge on {}: {}", a.name, max_diff);
    }
}

// ── Rayleigh physics for Landsat bands ──────────────────────────────────────

#[test]
fn test_rayleigh_tau_landsat_bands() {
    let taus: Vec<(f64, f64)> = L8_BANDS
        .iter()
        .map(|&(_, wl, _, _)| (wl, rayleigh_optical_thickness(wl, 1013.25)))
        .collect();

    // Must decrease with wavelength
    for i in 1..taus.len() {
        assert!(
            taus[i].1 < taus[i - 1].1,
            "τ({:.0})={:.4} >= τ({:.0})={:.4}",
            taus[i].0, taus[i].1, taus[i - 1].0, taus[i - 1].1
        );
    }
    // B1 (443nm) should have significant Rayleigh
    assert!(taus[0].1 > 0.1, "τ(443) too low: {}", taus[0].1);
    // B7 (2201nm) should have negligible Rayleigh
    assert!(taus[6].1 < 0.001, "τ(2201) too high: {}", taus[6].1);
}

// ── Performance scaling ─────────────────────────────────────────────────────

#[test]
fn test_parallel_speedup_scaling() {
    // Verify parallel is at least not slower than sequential for non-trivial sizes
    let (metadata, bands) = make_landsat_bands("L8_OLI", 500);
    let pipeline = make_pipeline(metadata, &bands);

    // Warm up
    let _ = process_bands_parallel(&pipeline, bands.clone());

    // Sequential timing
    let t0 = std::time::Instant::now();
    let _ = process_bands_sequential(&pipeline, bands.clone());
    let seq_time = t0.elapsed();

    // Parallel timing
    let t0 = std::time::Instant::now();
    let _ = process_bands_parallel(&pipeline, bands.clone());
    let par_time = t0.elapsed();

    println!(
        "Performance: seq={:.2?} par={:.2?} speedup={:.2}x",
        seq_time, par_time,
        seq_time.as_secs_f64() / par_time.as_secs_f64().max(1e-9)
    );

    // Parallel should not be more than 2x slower (accounts for small workload overhead)
    assert!(
        par_time.as_secs_f64() < seq_time.as_secs_f64() * 2.0,
        "Parallel significantly slower: {:.2?} vs {:.2?}", par_time, seq_time
    );
}

#[test]
fn test_throughput_megapixels_per_second() {
    let size = 1000;
    let (metadata, bands) = make_landsat_bands("L8_OLI", size);
    let pipeline = make_pipeline(metadata, &bands);
    let nbands = bands.len();

    let t0 = std::time::Instant::now();
    let _ = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = t0.elapsed();

    let mpx = (size * size * nbands) as f64 / elapsed.as_secs_f64() / 1e6;
    println!("Throughput: {:.1} Mpx/s ({} bands, {}×{})", mpx, nbands, size, size);

    // Minimum acceptable throughput: 1 Mpx/s (very conservative for CI)
    assert!(mpx > 1.0, "Throughput too low: {:.1} Mpx/s", mpx);
}

// ── COG output ──────────────────────────────────────────────────────────────

#[test]
fn test_landsat_cog_output() {
    use acolite_rs::writer::{cog_available, write_auto};

    if !cog_available() {
        eprintln!("Skipping COG test: gdal-support not available");
        return;
    }

    let (metadata, bands) = make_landsat_bands("L8_OLI", 50);
    let pipeline = make_pipeline(metadata.clone(), &bands);
    let corrected = process_bands_parallel(&pipeline, bands).unwrap();

    // 7 bands ≤ 50 → write_auto selects COG
    assert!(corrected.len() <= 50);

    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("landsat_test");
    write_auto(out.to_str().unwrap(), &corrected, &metadata).unwrap();

    // Should produce per-band .tif files (COG)
    let tif = tmp.path().join("landsat_test_B1.tif");
    assert!(tif.exists(), "COG/GeoTIFF not created");
    assert!(std::fs::metadata(&tif).unwrap().len() > 0);
}

// ── DSF AOT estimation ─────────────────────────────────────────────────────

#[test]
fn test_dsf_aot_from_swir_dark_spectrum() {
    // Simulate dark water pixels in SWIR
    let swir1 = Array2::from_elem((100, 100), 0.005); // very dark
    let swir2 = Array2::from_elem((100, 100), 0.002);

    let dark = estimate_dark_spectrum(&[swir1, swir2], &DarkSpectrumMethod::Percentile(5.0));
    let aot = optimize_aot_simple(&dark, &[1609.0, 2201.0], 35.0, 0.0);

    assert!(aot > 0.0 && aot < 1.0, "AOT={} out of range", aot);
}

#[test]
fn test_dsf_aot_increases_with_turbidity() {
    // Brighter dark spectrum → higher AOT
    let clear = estimate_dark_spectrum(&[Array2::from_elem((50, 50), 0.005), Array2::from_elem((50, 50), 0.002)], &DarkSpectrumMethod::Percentile(5.0));
    let turbid = estimate_dark_spectrum(&[Array2::from_elem((50, 50), 0.05), Array2::from_elem((50, 50), 0.03)], &DarkSpectrumMethod::Percentile(5.0));

    let aot_clear = optimize_aot_simple(&clear, &[1609.0, 2201.0], 35.0, 0.0);
    let aot_turbid = optimize_aot_simple(&turbid, &[1609.0, 2201.0], 35.0, 0.0);

    assert!(
        aot_turbid > aot_clear,
        "Turbid AOT ({}) should exceed clear ({})", aot_turbid, aot_clear
    );
}
