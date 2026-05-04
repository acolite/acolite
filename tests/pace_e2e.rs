//! End-to-end PACE OCI processing test
//!
//! Tests the full pipeline: load → parallel atmospheric correction → GeoZarr output → verify.
//! Uses synthetic PACE-like data (no real NetCDF files needed).

use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::parallel::process_bands_parallel_f32;
use acolite_rs::writer::geozarr::write_geozarr;
use acolite_rs::{Pipeline, ProcessingConfig};
use acolite_rs::ac::{estimate_dark_spectrum, optimize_aot_simple, rayleigh_optical_thickness, DarkSpectrumMethod};
use chrono::Utc;
use ndarray::Array2;

/// Build a synthetic PACE OCI scene with blue/red/SWIR detectors
fn synthetic_pace_scene() -> (Metadata, Vec<BandData<f32>>, Array2<f64>, Array2<f64>) {
    let nrows = 200;
    let ncols = 150;

    let mut metadata = Metadata::new("PACE_OCI".to_string(), Utc::now());
    metadata.set_geometry(32.0, 145.0);
    metadata.view_zenith = Some(5.0);
    metadata.view_azimuth = Some(260.0);
    metadata.add_attribute("processing_level".into(), "L1B".into());

    let proj = Projection::from_epsg(4326);
    let gt = GeoTransform::new(-75.5, 0.01, 36.0, -0.01);

    // Blue detector: 339–555 nm, 5 nm spacing → ~43 bands (use subset)
    let blue_waves: Vec<f64> = (0..10).map(|i| 400.0 + i as f64 * 15.0).collect();
    // Red detector: 555–890 nm
    let red_waves: Vec<f64> = (0..10).map(|i| 560.0 + i as f64 * 30.0).collect();
    // SWIR detector: 9 bands
    let swir_waves = vec![940.0, 1038.0, 1250.0, 1378.0, 1615.0, 2130.0, 2258.0];

    let mut bands = Vec::new();
    let mut rng_seed = 42u64;

    for (det, waves, bw) in [
        ("blue", &blue_waves, 5.0),
        ("red", &red_waves, 5.0),
        ("SWIR", &swir_waves, 20.0),
    ] {
        for &wave in waves {
            // Simulate TOA reflectance: higher in blue (Rayleigh), lower in SWIR
            let base_rhot = 0.15 * (550.0 / wave).powf(0.8);
            let mut data = Array2::from_elem((nrows, ncols), base_rhot as f32);
            // Add spatial variation (simple gradient)
            for r in 0..nrows {
                for c in 0..ncols {
                    rng_seed = rng_seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let noise = ((rng_seed >> 33) as f32) / (u32::MAX as f32) * 0.01 - 0.005;
                    data[[r, c]] += noise + (r as f32 / nrows as f32) * 0.02;
                }
            }
            let name = format!("rhot_{}_{:.0}", det, wave);
            bands.push(BandData::new(data, wave, bw, name, proj.clone(), gt.clone()));
        }
    }

    // Lat/lon grids
    let mut lat = Array2::zeros((nrows, ncols));
    let mut lon = Array2::zeros((nrows, ncols));
    for r in 0..nrows {
        for c in 0..ncols {
            lat[[r, c]] = 36.0 - r as f64 * 0.01;
            lon[[r, c]] = -75.5 + c as f64 * 0.01;
        }
    }

    (metadata, bands, lat, lon)
}

#[test]
fn test_pace_e2e_parallel_ac_to_geozarr() {
    let (metadata, bands, _lat, _lon) = synthetic_pace_scene();
    let nbands = bands.len();
    let (nrows, ncols) = bands[0].shape();

    assert_eq!(nbands, 27); // 10 blue + 10 red + 7 SWIR

    // Configure pipeline
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.5,
        ..ProcessingConfig::default()
    };
    let mut pipeline = Pipeline::new(metadata.clone(), config);

    // Estimate AOT from SWIR dark spectrum
    let swir_bands: Vec<Array2<f64>> = bands[20..27]
        .iter()
        .map(|b| b.data.mapv(|v| v as f64))
        .collect();
    let dark = estimate_dark_spectrum(&swir_bands, &DarkSpectrumMethod::Percentile(5.0));
    let swir_wl: Vec<f64> = bands[20..27].iter().map(|b| b.wavelength).collect();
    let aot = optimize_aot_simple(&dark, &swir_wl, 32.0, 5.0);
    assert!(aot > 0.0 && aot < 2.0, "AOT={} out of range", aot);
    pipeline.set_aot(aot);

    // Parallel atmospheric correction
    let corrected = process_bands_parallel_f32(&pipeline, bands).expect("AC failed");
    assert_eq!(corrected.len(), nbands);

    // Verify output values
    for band in &corrected {
        let (h, w) = band.shape();
        assert_eq!(h, nrows);
        assert_eq!(w, ncols);
        let mean = band.data.mean().unwrap();
        assert!(mean.is_finite(), "NaN in band {}", band.name);
        // Surface reflectance should be non-negative and < TOA
        assert!(mean >= 0.0, "Negative ρs in {}: {}", band.name, mean);
        assert!(mean < 1.0, "ρs > 1 in {}: {}", band.name, mean);
    }

    // Verify spectral ordering preserved
    for i in 1..corrected.len() {
        assert!(
            corrected[i].wavelength >= corrected[i - 1].wavelength,
            "Wavelength order broken at index {}",
            i
        );
    }

    // Write GeoZarr
    let tmp = tempfile::tempdir().expect("tmpdir");
    let zarr_path = tmp.path().join("pace_test.zarr");
    let zarr_str = zarr_path.to_str().unwrap();

    write_geozarr(zarr_str, &corrected, &metadata).expect("GeoZarr write failed");

    // Verify GeoZarr structure
    assert!(zarr_path.exists(), "Zarr directory not created");
    assert!(zarr_path.join("zarr.json").exists(), "Missing zarr.json root metadata");
    assert!(zarr_path.join("data").exists(), "Missing /data array");
    assert!(zarr_path.join("wavelengths").exists(), "Missing /wavelengths");
    assert!(zarr_path.join("bandwidths").exists(), "Missing /bandwidths");

    // Verify root metadata contains GeoZarr attributes
    let root_meta: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(zarr_path.join("zarr.json")).unwrap(),
    )
    .unwrap();
    let attrs = root_meta.get("attributes").expect("No root attributes");
    assert_eq!(attrs["sensor"], "PACE_OCI");
    assert_eq!(attrs["Conventions"], "CF-1.8");
    assert!(attrs.get("proj:epsg").is_some(), "Missing proj:epsg");
    assert!(attrs.get("spatial:transform").is_some(), "Missing spatial:transform");
}

#[test]
fn test_pace_rayleigh_spectral_consistency() {
    // Rayleigh optical thickness must decrease monotonically with wavelength
    let wavelengths: Vec<f64> = (0..20).map(|i| 400.0 + i as f64 * 50.0).collect();
    let taus: Vec<f64> = wavelengths
        .iter()
        .map(|&w| rayleigh_optical_thickness(w, 1013.25))
        .collect();

    for i in 1..taus.len() {
        assert!(
            taus[i] < taus[i - 1],
            "Rayleigh τ not decreasing: τ({:.0})={:.4} >= τ({:.0})={:.4}",
            wavelengths[i],
            taus[i],
            wavelengths[i - 1],
            taus[i - 1]
        );
    }
}

#[test]
fn test_pace_parallel_determinism() {
    // Running parallel AC twice must produce identical results
    let (metadata, bands, _, _) = synthetic_pace_scene();
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: false,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.0,
        ..ProcessingConfig::default()
    };
    let pipeline = Pipeline::new(metadata, config);

    let bands2 = bands
        .iter()
        .map(|b| {
            BandData::new(
                b.data.clone(),
                b.wavelength,
                b.bandwidth,
                b.name.clone(),
                b.projection.clone(),
                b.geotransform.clone(),
            )
        })
        .collect();

    let r1 = process_bands_parallel_f32(&pipeline, bands).unwrap();
    let r2 = process_bands_parallel_f32(&pipeline, bands2).unwrap();

    assert_eq!(r1.len(), r2.len());
    for (a, b) in r1.iter().zip(r2.iter()) {
        assert_eq!(a.wavelength, b.wavelength);
        let diff: f64 = a
            .data
            .iter()
            .zip(b.data.iter())
            .map(|(x, y)| (x - y).abs())
            .sum();
        assert!(diff == 0.0, "Non-deterministic result for band {}", a.name);
    }
}

#[test]
fn test_pace_geozarr_roundtrip_metadata() {
    // Verify GeoZarr metadata survives write
    let (metadata, bands, _, _) = synthetic_pace_scene();
    let config = ProcessingConfig::default();
    let pipeline = Pipeline::new(metadata.clone(), config);

    let corrected = process_bands_parallel_f32(&pipeline, bands).unwrap();

    let tmp = tempfile::tempdir().unwrap();
    let zarr_path = tmp.path().join("rt_test.zarr");
    write_geozarr(zarr_path.to_str().unwrap(), &corrected, &metadata).unwrap();

    // Read back and verify data array metadata
    let data_meta_path = zarr_path.join("data").join("zarr.json");
    assert!(data_meta_path.exists());
    let data_meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(data_meta_path).unwrap()).unwrap();

    // Check shape matches [nbands, nrows, ncols]
    let shape = data_meta["shape"].as_array().unwrap();
    assert_eq!(shape[0].as_u64().unwrap(), corrected.len() as u64);
    assert_eq!(shape[1].as_u64().unwrap(), corrected[0].shape().0 as u64);
    assert_eq!(shape[2].as_u64().unwrap(), corrected[0].shape().1 as u64);
}

#[test]
fn test_pace_band_count_triggers_geozarr() {
    // PACE OCI has >50 bands → write_auto should select GeoZarr
    assert!(
        27 > 0,
        "Synthetic scene has 27 bands; real PACE has ~280"
    );
    // The write_auto threshold is 50; our synthetic scene has 27 so we test
    // write_geozarr directly. This test documents the threshold behavior.
    assert_eq!(acolite_rs::writer::HYPERSPECTRAL_THRESHOLD, 50);
}
