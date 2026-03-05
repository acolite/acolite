//! Example: Full workflow with LUTs, NetCDF, and Zarr

use acolite_rs::{Pipeline, ProcessingConfig, ZarrWriter, ZarrReader, NetCdfWriter};
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::ac::{LutManager, estimate_dark_spectrum, optimize_aot};
use acolite_rs::parallel::process_bands_parallel;
use acolite_rs::simd::{parallel_sum, parallel_minmax};
use ndarray::Array2;
use chrono::Utc;
use tempfile::TempDir;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: Full Workflow Example\n");
    
    // 1. Setup LUT manager
    println!("→ Setting up LUT manager...");
    let lut_manager = LutManager::new(None);
    lut_manager.ensure_dir().unwrap();
    
    // Load Rayleigh LUT (will use default if not present)
    let rayleigh_lut = lut_manager.load_rayleigh_lut().unwrap();
    println!("  ✓ Loaded Rayleigh LUT with {} wavelengths", rayleigh_lut.wavelengths.len());
    
    // 2. Setup metadata
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(30.0, 135.0);
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
    
    // 4. Create pipeline
    let mut pipeline = Pipeline::new(metadata.clone(), config);
    
    // 5. Estimate AOT
    println!("\n→ Estimating AOT...");
    let dark_bands = vec![
        Array2::from_elem((100, 100), 0.02),
        Array2::from_elem((100, 100), 0.015),
        Array2::from_elem((100, 100), 0.01),
    ];
    
    let dark_spectrum = estimate_dark_spectrum(&dark_bands, 5.0);
    let aot = optimize_aot(&dark_spectrum, &[865.0, 1609.0, 2201.0], 30.0, 0.0);
    pipeline.set_aot(aot);
    println!("  AOT = {:.4}", aot);
    
    // 6. Create test bands
    println!("\n→ Creating test bands...");
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
    
    println!("  Created {} bands (100×100 pixels)", bands.len());
    
    // 7. Process bands
    println!("\n→ Processing bands (parallel)...");
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Processed in {:.2?}", elapsed);
    
    // 8. SIMD operations for statistics
    println!("\n→ Computing statistics (SIMD)...");
    for band in &result {
        let sum = parallel_sum(&band.data);
        let (min, max) = parallel_minmax(&band.data);
        let mean = sum / (band.data.nrows() * band.data.ncols()) as f64;
        
        println!("  {}: mean={:.4}, min={:.4}, max={:.4}", 
                 band.name, mean, min, max);
    }
    
    // 9. Write to Zarr (intermediate storage)
    println!("\n→ Writing to Zarr...");
    let temp_dir = TempDir::new().unwrap();
    let zarr_path = temp_dir.path().join("intermediate.zarr");
    
    let zarr_writer = ZarrWriter::new(&zarr_path);
    let start = std::time::Instant::now();
    zarr_writer.write_bands(&result).unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Wrote {} bands to Zarr in {:.2?}", result.len(), elapsed);
    
    // 10. Read back from Zarr
    println!("\n→ Reading from Zarr...");
    let zarr_reader = ZarrReader::new(&zarr_path);
    let start = std::time::Instant::now();
    let read_band = zarr_reader.read_band("B1").unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Read band B1 ({}×{}) in {:.2?}", 
             read_band.nrows(), read_band.ncols(), elapsed);
    
    // Verify data integrity
    let original_mean = result[0].data.mean().unwrap();
    let read_mean = read_band.mean().unwrap();
    println!("  Data integrity: original={:.6}, read={:.6}, diff={:.2e}",
             original_mean, read_mean, (original_mean - read_mean).abs());
    
    // 11. Write to NetCDF (final output)
    println!("\n→ Writing to NetCDF...");
    let nc_path = temp_dir.path().join("output.nc");
    let nc_writer = NetCdfWriter::new(nc_path.to_str().unwrap().to_string());
    nc_writer.write_bands(&result, &metadata).unwrap();
    println!("  ✓ Wrote to NetCDF (placeholder)");
    
    println!("\n✓ Full workflow complete!");
    println!("  - LUT management: ✓");
    println!("  - Atmospheric correction: ✓");
    println!("  - SIMD statistics: ✓");
    println!("  - Zarr I/O: ✓");
    println!("  - NetCDF output: ✓");
}
