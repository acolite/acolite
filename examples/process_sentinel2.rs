//! Example: Process Sentinel-2 scene with multi-resolution

use acolite_rs::{Pipeline, ProcessingConfig, ResampleMethod, resample};
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::sensors::{Sensor, Sentinel2Sensor};
use acolite_rs::ac::{estimate_dark_spectrum, optimize_aot};
use acolite_rs::parallel::process_bands_parallel;
use ndarray::Array2;
use chrono::Utc;

fn main() {
    println!("ACOLITE-RS: Sentinel-2 Processing Example\n");
    
    // 1. Setup Sentinel-2A sensor
    let s2a = Sentinel2Sensor::new_s2a();
    println!("✓ Sensor: {}", s2a.name());
    println!("  Bands: {}", s2a.band_names().len());
    
    // Show resolution groups
    for res in [10, 20, 60] {
        let bands = s2a.bands_at_resolution(res);
        println!("  {}m: {} bands", res, bands.len());
    }
    
    // 2. Setup metadata
    let mut metadata = Metadata::new(s2a.name().to_string(), Utc::now());
    metadata.set_geometry(25.0, 140.0);
    println!("\n✓ Metadata configured");
    
    // 3. Create processing configuration
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.0,
    };
    println!("✓ Processing config: Full AC chain");
    
    // 4. Create pipeline
    let mut pipeline = Pipeline::new(metadata, config);
    
    // 5. Simulate dark spectrum from 20m SWIR bands
    println!("\n→ Estimating AOT from SWIR bands (20m)...");
    let dark_bands = vec![
        Array2::from_elem((500, 500), 0.015), // B11
        Array2::from_elem((500, 500), 0.01),  // B12
    ];
    
    let dark_spectrum = estimate_dark_spectrum(&dark_bands, 5.0);
    let aot = optimize_aot(&dark_spectrum, &[1613.7, 2202.4], 25.0, 0.0);
    pipeline.set_aot(aot);
    println!("  AOT = {:.4}", aot);
    
    // 6. Create synthetic 10m bands
    println!("\n→ Processing 10m bands...");
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(600000.0, 10.0, 5000000.0, -10.0);
    
    let bands_10m: Vec<_> = vec![
        ("B02", 492.4, 66.0, 1200u16),   // Blue
        ("B03", 559.8, 36.0, 1500u16),   // Green
        ("B04", 664.6, 31.0, 1800u16),   // Red
        ("B08", 832.8, 106.0, 2200u16),  // NIR
    ].into_iter().map(|(name, wl, bw, dn)| {
        BandData::new(
            Array2::from_elem((1000, 1000), dn),
            wl,
            bw,
            name.to_string(),
            proj.clone(),
            geotrans.clone(),
        )
    }).collect();
    
    let start = std::time::Instant::now();
    let result_10m = process_bands_parallel(&pipeline, bands_10m).unwrap();
    let elapsed_10m = start.elapsed();
    
    println!("  ✓ Processed {} bands (1000×1000) in {:.2?}", result_10m.len(), elapsed_10m);
    
    // 7. Process 20m bands
    println!("\n→ Processing 20m bands...");
    let geotrans_20m = GeoTransform::new(600000.0, 20.0, 5000000.0, -20.0);
    
    let bands_20m: Vec<_> = vec![
        ("B05", 704.1, 15.0, 1900u16),   // Red Edge 1
        ("B06", 740.5, 15.0, 2000u16),   // Red Edge 2
        ("B8A", 864.7, 21.0, 2100u16),   // NIR narrow
    ].into_iter().map(|(name, wl, bw, dn)| {
        BandData::new(
            Array2::from_elem((500, 500), dn),
            wl,
            bw,
            name.to_string(),
            proj.clone(),
            geotrans_20m.clone(),
        )
    }).collect();
    
    let start = std::time::Instant::now();
    let result_20m = process_bands_parallel(&pipeline, bands_20m).unwrap();
    let elapsed_20m = start.elapsed();
    
    println!("  ✓ Processed {} bands (500×500) in {:.2?}", result_20m.len(), elapsed_20m);
    
    // 8. Demonstrate resampling 20m to 10m
    println!("\n→ Resampling 20m bands to 10m...");
    let start = std::time::Instant::now();
    
    let resampled: Vec<_> = result_20m.iter().map(|band| {
        let resampled_data = resample(&band.data, 20, 10, ResampleMethod::Bilinear).unwrap();
        (band.name.clone(), resampled_data.dim())
    }).collect();
    
    let elapsed_resample = start.elapsed();
    
    for (name, shape) in &resampled {
        println!("  {} resampled to {}×{}", name, shape.0, shape.1);
    }
    println!("  ✓ Resampled in {:.2?}", elapsed_resample);
    
    // 9. Display results
    println!("\n→ Results (10m bands):");
    for band in &result_10m {
        let mean = band.data.mean().unwrap();
        println!("  {}: ρ_mean={:.4} ({:.1}nm)", band.name, mean, band.wavelength);
    }
    
    println!("\n→ Results (20m bands):");
    for band in &result_20m {
        let mean = band.data.mean().unwrap();
        println!("  {}: ρ_mean={:.4} ({:.1}nm)", band.name, mean, band.wavelength);
    }
    
    println!("\n✓ Sentinel-2 processing complete!");
    println!("  Total time: {:.2?}", elapsed_10m + elapsed_20m + elapsed_resample);
}
