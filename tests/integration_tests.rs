//! Integration tests for ACOLITE-RS

use acolite_rs::{Pipeline, ProcessingConfig};
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::ac::{estimate_dark_spectrum, optimize_aot};
use acolite_rs::parallel::process_bands_parallel;
use ndarray::Array2;
use chrono::Utc;

#[test]
fn test_full_atmospheric_correction_pipeline() {
    // Setup metadata
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(30.0, 135.0); // 30° zenith, 135° azimuth
    
    // Create configuration
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.0,
    };
    
    // Create pipeline
    let mut pipeline = Pipeline::new(metadata, config);
    
    // Simulate dark spectrum estimation
    let dark_bands = vec![
        Array2::from_elem((100, 100), 0.02),
        Array2::from_elem((100, 100), 0.015),
        Array2::from_elem((100, 100), 0.01),
    ];
    
    let dark_spectrum = estimate_dark_spectrum(&dark_bands, 5.0);
    let aot = optimize_aot(&dark_spectrum, &[865.0, 1609.0, 2201.0], 30.0, 0.0);
    pipeline.set_aot(aot);
    
    // Create test bands
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);
    
    let bands: Vec<_> = vec![
        ("B1", 443.0, 16.0, 1000u16),
        ("B2", 482.0, 60.0, 1200u16),
        ("B3", 561.0, 57.0, 1500u16),
        ("B4", 655.0, 37.0, 1800u16),
    ].into_iter().map(|(name, wl, bw, dn)| {
        BandData::new(
            Array2::from_elem((100, 100), dn),
            wl,
            bw,
            name.to_string(),
            proj.clone(),
            geotrans.clone(),
        )
    }).collect();
    
    // Process bands
    let result = process_bands_parallel(&pipeline, bands);
    
    assert!(result.is_ok());
    let corrected = result.unwrap();
    assert_eq!(corrected.len(), 4);
    
    // Check output values are reasonable
    for band in &corrected {
        let mean = band.data.mean().unwrap();
        assert!(mean > 0.0 && mean < 1.0, "Mean reflectance should be 0-1");
    }
}

#[test]
fn test_landsat_processing_workflow() {
    // Simulate Landsat 8 processing
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(25.0, 140.0);
    
    let config = ProcessingConfig::default();
    let pipeline = Pipeline::new(metadata, config);
    
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(0.0, 30.0, 0.0, -30.0);
    
    // Process single band
    let band = BandData::new(
        Array2::from_elem((50, 50), 2000u16),
        655.0,
        37.0,
        "B4".to_string(),
        proj,
        geotrans,
    );
    
    let result = pipeline.process_band(band);
    assert!(result.is_ok());
    
    let corrected = result.unwrap();
    assert_eq!(corrected.shape(), (50, 50));
    assert_eq!(corrected.wavelength, 655.0);
}

#[test]
fn test_atmospheric_correction_components() {
    use acolite_rs::ac::{
        rayleigh_optical_thickness,
        ozone_transmittance,
        water_vapor_transmittance,
    };
    
    // Test Rayleigh optical thickness
    let tau_443 = rayleigh_optical_thickness(443.0, 1013.25);
    assert!(tau_443 > 0.0 && tau_443 < 1.0);
    
    let tau_865 = rayleigh_optical_thickness(865.0, 1013.25);
    assert!(tau_865 < tau_443); // Should decrease with wavelength
    
    // Test gas transmittance
    let t_o3 = ozone_transmittance(500.0, 0.3, 2.0);
    assert!(t_o3 > 0.0 && t_o3 <= 1.0);
    
    let t_wv = water_vapor_transmittance(940.0, 2.0, 2.0);
    assert!(t_wv > 0.0 && t_wv <= 1.0);
}

#[test]
fn test_parallel_vs_sequential_consistency() {
    use acolite_rs::parallel::{process_bands_parallel, process_bands_sequential};
    
    let metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    let config = ProcessingConfig::default();
    let pipeline = Pipeline::new(metadata, config);
    
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(0.0, 30.0, 0.0, -30.0);
    
    let bands: Vec<_> = (0..3).map(|i| {
        BandData::new(
            Array2::from_elem((20, 20), 1000u16 + i * 100),
            400.0 + i as f64 * 100.0,
            50.0,
            format!("B{}", i + 1),
            proj.clone(),
            geotrans.clone(),
        )
    }).collect();
    
    let parallel_result = process_bands_parallel(&pipeline, bands.clone()).unwrap();
    let sequential_result = process_bands_sequential(&pipeline, bands).unwrap();
    
    // Results should be identical
    for (p, s) in parallel_result.iter().zip(sequential_result.iter()) {
        assert_eq!(p.data.shape(), s.data.shape());
        assert_eq!(p.wavelength, s.wavelength);
    }
}
