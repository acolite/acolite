//! Example: Process Landsat scene

use acolite_rs::{Pipeline, ProcessingConfig};
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::ac::{estimate_dark_spectrum, optimize_aot};
use acolite_rs::parallel::process_bands_parallel;
use ndarray::Array2;
use chrono::Utc;

fn main() {
    println!("ACOLITE-RS: Landsat Processing Example\n");
    
    // 1. Setup metadata
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(30.0, 135.0);
    println!("✓ Metadata configured");
    
    // 2. Create processing configuration
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.0,
    };
    println!("✓ Processing config: Rayleigh + Gas + Aerosol");
    
    // 3. Create pipeline
    let mut pipeline = Pipeline::new(metadata, config);
    
    // 4. Estimate AOT using dark spectrum (NIR/SWIR bands)
    println!("\n→ Estimating AOT from dark spectrum...");
    let dark_bands = vec![
        Array2::from_elem((100, 100), 0.02),
        Array2::from_elem((100, 100), 0.015),
        Array2::from_elem((100, 100), 0.01),
    ];
    
    let dark_spectrum = estimate_dark_spectrum(&dark_bands, 5.0);
    let aot = optimize_aot(&dark_spectrum, &[865.0, 1609.0, 2201.0], 30.0, 0.0);
    pipeline.set_aot(aot);
    println!("  AOT = {:.4}", aot);
    
    // 5. Create synthetic bands (in real use, read from file)
    println!("\n→ Creating test bands...");
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);
    
    let bands: Vec<_> = vec![
        ("B1", 443.0, 16.0, 1000u16),
        ("B2", 482.0, 60.0, 1200u16),
        ("B3", 561.0, 57.0, 1500u16),
        ("B4", 655.0, 37.0, 1800u16),
        ("B5", 865.0, 28.0, 2000u16),
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
    
    println!("  Created {} bands (100×100 pixels)", bands.len());
    
    // 6. Process bands in parallel
    println!("\n→ Processing bands (parallel)...");
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = start.elapsed();
    
    println!("  ✓ Processed in {:.2?}", elapsed);
    
    // 7. Display results
    println!("\n→ Results:");
    for band in &result {
        let mean = band.data.mean().unwrap();
        let min = band.data.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max = band.data.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        
        println!("  {}: ρ_mean={:.4}, ρ_min={:.4}, ρ_max={:.4}", 
                 band.name, mean, min, max);
    }
    
    println!("\n✓ Processing complete!");
}
