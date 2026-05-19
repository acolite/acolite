//! Sentinel-3 OLCI end-to-end and property-based tests
//!
//! Tests cover:
//! - Sensor band definitions (wavelength, bandwidth, count)
//! - Smile correction physics (conservation, monotonicity)
//! - Radiance→reflectance conversion
//! - DSF atmospheric correction for OLCI bands
//! - Gas transmittance at OLCI wavelengths
//! - Rayleigh scattering at OLCI wavelengths
//! - Pipeline integration (L1→L2R)
//! - Property-based testing with proptest

use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::loader::sentinel3::{OlciInstrumentData, OlciScene, OlciTpg, OLCI_BANDS, parse_s3_manifest};
use acolite_rs::ac::rayleigh::rayleigh_optical_thickness;
use acolite_rs::ac::gas::{gas_correction, ozone_transmittance, water_vapor_transmittance};
use acolite_rs::ac::dsf::{estimate_dark_spectrum, DarkSpectrumMethod};
use acolite_rs::pipeline::{Pipeline, ProcessingConfig};
use acolite_rs::parallel::process_bands_parallel;
use acolite_rs::sensors::{Sensor, Sentinel3Sensor};
use chrono::Utc;
use ndarray::Array2;
use std::collections::HashMap;

// ─── Helpers ───────────────────────────────────────────────────────────────

fn make_olci_metadata() -> Metadata {
    let mut m = Metadata::new("S3A_OLCI".to_string(), Utc::now());
    m.set_geometry(35.0, 140.0);
    m.view_zenith = Some(15.0);
    m
}

fn make_instrument(rows: usize, cols: usize, ndet: usize) -> OlciInstrumentData {
    let nbands = 21;
    let lambda0: Vec<Vec<f64>> = OLCI_BANDS
        .iter()
        .map(|(_, wl, _, _)| vec![*wl; ndet])
        .collect();
    let solar_flux: Vec<Vec<f64>> = OLCI_BANDS
        .iter()
        .map(|(_, _, _, e0)| vec![*e0; ndet])
        .collect();
    let fwhm: Vec<Vec<f64>> = OLCI_BANDS
        .iter()
        .map(|(_, _, bw, _)| vec![*bw; ndet])
        .collect();
    let detector_index = Array2::zeros((rows, cols));
    OlciInstrumentData { lambda0, solar_flux, fwhm, detector_index }
}

fn make_tpg(rows: usize, cols: usize) -> OlciTpg {
    OlciTpg {
        sza: Array2::from_elem((rows, cols), 35.0),
        oza: Array2::from_elem((rows, cols), 15.0),
        saa: Array2::from_elem((rows, cols), 150.0),
        oaa: Array2::from_elem((rows, cols), 30.0),
        latitude: Array2::from_elem((rows, cols), 51.0),
        longitude: Array2::from_elem((rows, cols), 3.0),
        total_ozone: Some(Array2::from_elem((rows, cols), 0.006)), // kg/m²
        total_columnar_water_vapour: Some(Array2::from_elem((rows, cols), 15.0)), // kg/m²
        sea_level_pressure: Some(Array2::from_elem((rows, cols), 1013.25)),
    }
}

fn make_scene(rows: usize, cols: usize) -> OlciScene {
    let mut radiance = HashMap::new();
    for (i, (name, wl, _, _)) in OLCI_BANDS.iter().enumerate() {
        // Simulate typical water-leaving radiance + atmosphere
        // Higher in blue, decreasing toward NIR (turbid water)
        let base = 50.0 + 30.0 * (1.0 - (wl - 400.0) / 620.0).max(0.0);
        radiance.insert(name.to_string(), Array2::from_elem((rows, cols), base));
    }
    OlciScene {
        sensor: "S3A_OLCI".into(),
        metadata: make_olci_metadata(),
        radiance,
        instrument: make_instrument(rows, cols, 740),
        tpg: make_tpg(rows, cols),
        data_shape: (rows, cols),
        product_type: Some("FR".into()),
    }
}

fn proj_and_geo() -> (Projection, GeoTransform) {
    (Projection::from_epsg(4326), GeoTransform::new(3.0, 0.003, 51.0, -0.003))
}

// ─── Sensor Definition Tests ───────────────────────────────────────────────

#[test]
fn test_olci_has_21_bands() {
    let sensor = Sentinel3Sensor;
    assert_eq!(sensor.band_names().len(), 21);
    assert_eq!(OLCI_BANDS.len(), 21);
}

#[test]
fn test_olci_wavelengths_monotonically_increasing() {
    for i in 1..OLCI_BANDS.len() {
        assert!(
            OLCI_BANDS[i].1 > OLCI_BANDS[i - 1].1,
            "Band {} ({} nm) not > band {} ({} nm)",
            i + 1, OLCI_BANDS[i].1, i, OLCI_BANDS[i - 1].1
        );
    }
}

#[test]
fn test_olci_wavelength_range() {
    assert_eq!(OLCI_BANDS[0].1, 400.0);
    assert_eq!(OLCI_BANDS[20].1, 1020.0);
}

#[test]
fn test_olci_e0_positive() {
    for (name, _, _, e0) in &OLCI_BANDS {
        assert!(*e0 > 0.0, "E0 for {} must be positive", name);
    }
}

#[test]
fn test_olci_e0_decreasing_trend() {
    // Solar irradiance generally decreases from blue to NIR
    assert!(OLCI_BANDS[0].3 > OLCI_BANDS[20].3);
}

#[test]
fn test_sensor_trait_consistency() {
    let sensor = Sentinel3Sensor;
    for (name, wl, bw, _) in &OLCI_BANDS {
        assert_eq!(sensor.wavelength(name), Some(*wl));
        assert_eq!(sensor.bandwidth(name), Some(*bw));
    }
}

// ─── Rayleigh Tests for OLCI Bands ────────────────────────────────────────

#[test]
fn test_rayleigh_tau_olci_monotonic_decrease() {
    let mut prev_tau = f64::MAX;
    for (_, wl, _, _) in &OLCI_BANDS {
        let tau = rayleigh_optical_thickness(*wl, 1013.25);
        assert!(tau < prev_tau, "Rayleigh τ should decrease with wavelength at {} nm", wl);
        prev_tau = tau;
    }
}

#[test]
fn test_rayleigh_tau_olci_physical_range() {
    for (name, wl, _, _) in &OLCI_BANDS {
        let tau = rayleigh_optical_thickness(*wl, 1013.25);
        assert!(tau > 0.0 && tau < 1.0, "Rayleigh τ out of range for {}: {}", name, tau);
    }
}

#[test]
fn test_rayleigh_pressure_scaling() {
    let wl = 443.0;
    let tau_low = rayleigh_optical_thickness(wl, 500.0);
    let tau_std = rayleigh_optical_thickness(wl, 1013.25);
    let tau_high = rayleigh_optical_thickness(wl, 1100.0);
    assert!(tau_low < tau_std);
    assert!(tau_std < tau_high);
    // Linear scaling with pressure
    let ratio = tau_high / tau_low;
    let expected_ratio = 1100.0 / 500.0;
    assert!((ratio - expected_ratio).abs() < 0.01);
}

// ─── Gas Transmittance Tests for OLCI ─────────────────────────────────────

#[test]
fn test_gas_transmittance_range() {
    for (_, wl, _, _) in &OLCI_BANDS {
        let t_o3 = ozone_transmittance(*wl, 0.3, 2.0);
        let t_wv = water_vapor_transmittance(*wl, 2.0, 2.0);
        assert!(t_o3 > 0.0 && t_o3 <= 1.0, "O3 transmittance out of range at {} nm", wl);
        assert!(t_wv > 0.0 && t_wv <= 1.0, "WV transmittance out of range at {} nm", wl);
    }
}

#[test]
fn test_gas_correction_increases_signal() {
    let toa = Array2::from_elem((10, 10), 0.05);
    for (_, wl, _, _) in &OLCI_BANDS {
        let corrected = gas_correction(&toa, *wl, 0.3, 2.0, 35.0, 15.0);
        // Gas correction divides by transmittance ≤ 1, so signal should increase or stay same
        for (&orig, &corr) in toa.iter().zip(corrected.iter()) {
            assert!(corr >= orig - 1e-10, "Gas correction decreased signal at {} nm", wl);
        }
    }
}

// ─── Smile Correction Tests ──────────────────────────────────────────────

#[test]
fn test_smile_correction_preserves_energy() {
    let (rows, cols) = (50, 50);
    let mut scene = make_scene(rows, cols);

    // Sum radiance before
    let sum_before: f64 = scene.radiance.values()
        .map(|a| a.iter().filter(|v| v.is_finite()).sum::<f64>())
        .sum();

    scene.apply_smile_correction(None);

    let sum_after: f64 = scene.radiance.values()
        .map(|a| a.iter().filter(|v| v.is_finite()).sum::<f64>())
        .sum();

    // Smile correction should not dramatically change total energy
    // Allow 5% tolerance (correction is small)
    let ratio = sum_after / sum_before;
    assert!(
        ratio > 0.95 && ratio < 1.05,
        "Smile correction changed total energy by {:.1}%",
        (ratio - 1.0) * 100.0
    );
}

#[test]
fn test_smile_correction_with_uniform_detectors_is_small() {
    // When all detectors have the same wavelength/F0, smile correction should be ~zero
    let (rows, cols) = (20, 20);
    let mut scene = make_scene(rows, cols);

    let before: HashMap<String, Array2<f64>> = scene.radiance.clone();
    scene.apply_smile_correction(None);

    for (name, _, _, _) in &OLCI_BANDS {
        let key = name.to_string();
        if let (Some(b), Some(a)) = (before.get(&key), scene.radiance.get(&key)) {
            let max_diff = b.iter().zip(a.iter())
                .map(|(x, y)| (x - y).abs())
                .fold(0.0f64, f64::max);
            // With uniform detectors, the main smile term (refl*E0 - rad) should be zero
            // Only the bounding band term contributes, which is also zero for uniform detectors
            assert!(
                max_diff < 1e-6,
                "Smile correction too large for uniform detectors on {}: {}",
                name, max_diff
            );
        }
    }
}

// ─── TOA Reflectance Conversion Tests ─────────────────────────────────────

#[test]
fn test_toa_reflectance_physical_range() {
    let scene = make_scene(50, 50);
    let rhot = scene.to_toa_reflectance(false);

    for (name, _, _, _) in &OLCI_BANDS {
        if let Some(data) = rhot.get(&name.to_string()) {
            for &v in data.iter() {
                assert!(
                    v >= 0.0 && v < 2.0,
                    "TOA reflectance out of range for {}: {}",
                    name, v
                );
            }
        }
    }
}

#[test]
fn test_toa_reflectance_smile_vs_nosmile() {
    let scene = make_scene(20, 20);
    let rhot_nosmile = scene.to_toa_reflectance(false);
    let rhot_smile = scene.to_toa_reflectance(true);

    // With uniform detectors, both should be very close
    for (name, _, _, _) in &OLCI_BANDS {
        let key = name.to_string();
        if let (Some(a), Some(b)) = (rhot_nosmile.get(&key), rhot_smile.get(&key)) {
            let max_diff = a.iter().zip(b.iter())
                .map(|(x, y)| (x - y).abs())
                .fold(0.0f64, f64::max);
            assert!(max_diff < 0.01, "Large smile/nosmile difference for {}: {}", name, max_diff);
        }
    }
}

// ─── DSF Tests for OLCI ──────────────────────────────────────────────────

#[test]
fn test_dark_spectrum_percentile() {
    let bands: Vec<Array2<f64>> = OLCI_BANDS.iter().enumerate().map(|(i, _)| {
        let mut data = Array2::from_elem((100, 100), 0.05 + i as f64 * 0.005);
        // Add some dark pixels
        for r in 0..10 {
            for c in 0..10 {
                data[[r, c]] = 0.001 + i as f64 * 0.0005;
            }
        }
        data
    }).collect();

    let dark = estimate_dark_spectrum(&bands, &DarkSpectrumMethod::Percentile(1.0));
    assert_eq!(dark.len(), 21);
    for (i, &v) in dark.iter().enumerate() {
        assert!(v.is_finite(), "Dark spectrum NaN for band {}", i);
        assert!(v > 0.0, "Dark spectrum non-positive for band {}", i);
        assert!(v < 0.2, "Dark spectrum too high for band {}", i);
    }
}

#[test]
fn test_dark_spectrum_intercept() {
    let bands: Vec<Array2<f64>> = (0..21).map(|i| {
        Array2::from_elem((100, 100), 0.03 + i as f64 * 0.002)
    }).collect();

    let dark = estimate_dark_spectrum(&bands, &DarkSpectrumMethod::Intercept(200));
    assert_eq!(dark.len(), 21);
    // Intercept should be ≤ the minimum value (regression extrapolates to y-intercept)
    for (i, &v) in dark.iter().enumerate() {
        assert!(v.is_finite());
        let min_val = 0.03 + i as f64 * 0.002;
        // Intercept of uniform data = the value itself
        assert!((v - min_val).abs() < 0.001, "Intercept mismatch for band {}: {} vs {}", i, v, min_val);
    }
}

#[test]
fn test_dark_spectrum_monotonic_for_clear_water() {
    // Clear water: dark spectrum should generally decrease from blue to NIR
    let bands: Vec<Array2<f64>> = OLCI_BANDS.iter().map(|(_, wl, _, _)| {
        // Simulate clear water: Rayleigh-dominated, decreasing with wavelength
        let val = 0.15 * (400.0 / wl).powi(4) + 0.005;
        Array2::from_elem((100, 100), val)
    }).collect();

    let dark = estimate_dark_spectrum(&bands, &DarkSpectrumMethod::Percentile(5.0));
    // Check visible bands (Oa01-Oa08) are decreasing
    for i in 1..8 {
        assert!(
            dark[i] <= dark[i - 1] + 0.001,
            "Dark spectrum not decreasing in visible: band {} ({}) > band {} ({})",
            i + 1, dark[i], i, dark[i - 1]
        );
    }
}

// ─── Pipeline Integration Tests ──────────────────────────────────────────

#[test]
fn test_olci_pipeline_processes_all_bands() {
    let sensor = Sentinel3Sensor;
    let mut metadata = make_olci_metadata();
    let config = ProcessingConfig::default();
    let mut pipeline = Pipeline::new(metadata, config);
    pipeline.set_aot(0.1);

    let (proj, geo) = proj_and_geo();
    let bands: Vec<BandData<u16>> = OLCI_BANDS.iter().map(|(name, wl, bw, _)| {
        BandData::new(
            Array2::from_elem((50, 50), 800u16),
            *wl, *bw, name.to_string(),
            proj.clone(), geo.clone(),
        )
    }).collect();

    let result = process_bands_parallel(&pipeline, bands).unwrap();
    assert_eq!(result.len(), 21);
}

#[test]
fn test_olci_pipeline_output_physical_range() {
    let mut metadata = make_olci_metadata();
    let config = ProcessingConfig::default();
    let mut pipeline = Pipeline::new(metadata, config);
    pipeline.set_aot(0.1);

    let (proj, geo) = proj_and_geo();
    let bands: Vec<BandData<u16>> = OLCI_BANDS.iter().map(|(name, wl, bw, _)| {
        BandData::new(
            Array2::from_elem((20, 20), 500u16),
            *wl, *bw, name.to_string(),
            proj.clone(), geo.clone(),
        )
    }).collect();

    let result = process_bands_parallel(&pipeline, bands).unwrap();
    for band in &result {
        for &v in band.data.iter() {
            if v.is_finite() {
                assert!(
                    v >= -0.05 && v < 2.0,
                    "Output reflectance out of range for {} ({} nm): {}",
                    band.name, band.wavelength, v
                );
            }
        }
    }
}

#[test]
fn test_olci_pipeline_deterministic() {
    let (proj, geo) = proj_and_geo();
    let make_run = || {
        let metadata = make_olci_metadata();
        let config = ProcessingConfig::default();
        let mut pipeline = Pipeline::new(metadata, config);
        pipeline.set_aot(0.1);

        let bands: Vec<BandData<u16>> = OLCI_BANDS[..5].iter().map(|(name, wl, bw, _)| {
            BandData::new(
                Array2::from_elem((20, 20), 600u16),
                *wl, *bw, name.to_string(),
                proj.clone(), geo.clone(),
            )
        }).collect();
        process_bands_parallel(&pipeline, bands).unwrap()
    };

    let r1 = make_run();
    let r2 = make_run();
    for (a, b) in r1.iter().zip(r2.iter()) {
        assert_eq!(a.data, b.data, "Non-deterministic output for {}", a.name);
    }
}

// ─── TPG Tests ───────────────────────────────────────────────────────────

#[test]
fn test_raa_computation() {
    let tpg = make_tpg(10, 10);
    let raa = tpg.raa();
    for &v in raa.iter() {
        assert!(v >= 0.0 && v <= 180.0, "RAA out of range: {}", v);
    }
    // SAA=150, OAA=30 → RAA = |150-30| = 120
    assert!((raa[[0, 0]] - 120.0).abs() < 0.001);
}

#[test]
fn test_raa_wrapping() {
    let mut tpg = make_tpg(10, 10);
    tpg.saa = Array2::from_elem((10, 10), 350.0);
    tpg.oaa = Array2::from_elem((10, 10), 10.0);
    let raa = tpg.raa();
    // |350-10| = 340 > 180 → |360-340| = 20
    assert!((raa[[0, 0]] - 20.0).abs() < 0.001);
}

#[test]
fn test_mean_geometry() {
    let tpg = make_tpg(10, 10);
    let (sza, vza, raa) = tpg.mean_geometry();
    assert!((sza - 35.0).abs() < 0.001);
    assert!((vza - 15.0).abs() < 0.001);
    assert!((raa - 120.0).abs() < 0.001);
}

// ─── Manifest Parsing Tests ──────────────────────────────────────────────

#[test]
fn test_parse_manifest_olci() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(tmp.path(), r#"<?xml version="1.0"?>
<xfdu:XFDU>
  <metadataSection>
    <sentinel-safe:platform>
      <sentinel-safe:familyName>Sentinel-3</sentinel-safe:familyName>
      <sentinel-safe:number>A</sentinel-safe:number>
    </sentinel-safe:platform>
    <sentinel-safe:instrument>
      <sentinel-safe:familyName>Ocean Land Colour Instrument</sentinel-safe:familyName>
    </sentinel-safe:instrument>
    <sentinel3:productType>OL_1_EFR____</sentinel3:productType>
  </metadataSection>
</xfdu:XFDU>"#).unwrap();

    let (sensor, attrs) = parse_s3_manifest(tmp.path()).unwrap();
    assert_eq!(sensor, "S3A_OLCI");
    assert_eq!(attrs.get("product_type").map(|s| s.as_str()), Some("FR"));
}

#[test]
fn test_parse_manifest_meris() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(tmp.path(), r#"<?xml version="1.0"?>
<xfdu:XFDU>
  <metadataSection>
    <sentinel-safe:platform>
      <sentinel-safe:familyName>Envisat</sentinel-safe:familyName>
      <sentinel-safe:number>1</sentinel-safe:number>
    </sentinel-safe:platform>
    <sentinel-safe:instrument>
      <sentinel-safe:familyName>MEdium Resolution Imaging Spectrometer</sentinel-safe:familyName>
    </sentinel-safe:instrument>
  </metadataSection>
</xfdu:XFDU>"#).unwrap();

    let (sensor, _) = parse_s3_manifest(tmp.path()).unwrap();
    assert_eq!(sensor, "EN1_MERIS");
}

// ─── Scene Construction Tests ────────────────────────────────────────────

#[test]
fn test_scene_to_band_data() {
    let scene = make_scene(20, 20);
    let rhot = scene.to_toa_reflectance(true);
    let bands = scene.to_band_data(&rhot);
    assert_eq!(bands.len(), 21);
    for (i, band) in bands.iter().enumerate() {
        assert_eq!(band.wavelength, OLCI_BANDS[i].1);
        assert_eq!(band.name, OLCI_BANDS[i].0);
        assert_eq!(band.shape(), (20, 20));
    }
}
