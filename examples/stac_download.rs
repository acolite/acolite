//! Download and process real satellite data using STAC

use acolite_rs::{Pipeline, ProcessingConfig, StacClient, search_landsat, write_cog, write_geozarr, cog_available};
use acolite_rs::core::Metadata;
use acolite_rs::parallel::process_bands_parallel;
use chrono::Utc;
use std::path::Path;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: STAC Data Download & Processing\n");
    
    // Check COG support
    if cog_available() {
        println!("✓ COG tools available (gdal_translate)");
    } else {
        println!("⚠ COG tools not available - will create regular GeoTIFF");
    }
    
    // Example: San Francisco Bay Area
    let bbox = [-122.5, 37.5, -122.0, 38.0];
    let datetime = "2023-06-01/2023-06-30";
    let max_cloud = 10.0;
    
    println!("\n→ Searching for Landsat scenes...");
    println!("  Area: San Francisco Bay");
    println!("  Period: {}", datetime);
    println!("  Max cloud: {}%", max_cloud);
    
    match search_landsat(&bbox, datetime, max_cloud) {
        Ok(items) => {
            if items.is_empty() {
                println!("\n✗ No scenes found");
                println!("\nRunning with synthetic data instead...\n");
                run_synthetic();
                return;
            }
            
            println!("\n✓ Found {} scenes", items.len());
            
            // Use first scene
            let item = &items[0];
            println!("\nSelected scene: {}", item.id);
            
            if let Some(cloud) = item.properties.get("eo:cloud_cover") {
                println!("  Cloud cover: {}%", cloud);
            }
            
            // List available assets
            println!("\n  Available assets:");
            for (name, asset) in &item.assets {
                if let Some(title) = &asset.title {
                    println!("    - {}: {}", name, title);
                }
            }
            
            println!("\n⚠ Note: Actual download requires authentication");
            println!("  For demo purposes, using synthetic data\n");
            run_synthetic();
        }
        Err(e) => {
            println!("\n✗ STAC search failed: {}", e);
            println!("  This may be due to network issues or API changes");
            println!("\nRunning with synthetic data instead...\n");
            run_synthetic();
        }
    }
}

fn run_synthetic() {
    use acolite_rs::core::{BandData, Projection, GeoTransform};
    use ndarray::Array2;
    use tempfile::TempDir;
    
    println!("→ Creating synthetic multispectral data (Landsat-like)...");
    
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(30.0, -122.0);
    
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.0,
    };
    
    let mut pipeline = Pipeline::new(metadata.clone(), config);
    pipeline.set_aot(0.15);
    
    // Create multispectral bands (7 bands)
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);
    
    let wavelengths = [443.0, 482.0, 561.0, 655.0, 865.0, 1609.0, 2201.0];
    let bandwidths = [16.0, 60.0, 57.0, 37.0, 28.0, 85.0, 187.0];
    
    let bands: Vec<_> = (0..7)
        .map(|i| {
            let dn = 8000u16 + (i as u16 * 1000);
            BandData::new(
                Array2::from_elem((1000, 1000), dn),
                wavelengths[i],
                bandwidths[i],
                format!("B{}", i + 1),
                proj.clone(),
                geotrans.clone(),
            )
        })
        .collect();
    
    println!("  Created {} bands (1000×1000 pixels)", bands.len());
    
    println!("\n→ Processing atmospheric correction...");
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Processed in {:.2?}", elapsed);
    
    // Write as COG (multispectral)
    let temp = TempDir::new().unwrap();
    let cog_path = temp.path().join("multispectral_output.tif");
    
    println!("\n→ Writing COG (multispectral)...");
    let start = std::time::Instant::now();
    write_cog(cog_path.to_str().unwrap(), &result, &metadata).unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Written in {:.2?}", elapsed);
    println!("  File: {:?}", cog_path);
    println!("  Size: {} MB", std::fs::metadata(&cog_path).unwrap().len() / 1_000_000);
    
    // Create hyperspectral example
    println!("\n→ Creating synthetic hyperspectral data...");
    
    let hyper_bands: Vec<_> = (0..100)
        .map(|i| {
            let dn = 5000u16 + (i as u16 * 50);
            BandData::new(
                Array2::from_elem((500, 500), dn),
                400.0 + i as f64 * 5.0, // 400-895 nm, 5nm spacing
                2.5,
                format!("Band_{:03}", i + 1),
                proj.clone(),
                geotrans.clone(),
            )
        })
        .collect();
    
    println!("  Created {} bands (500×500 pixels)", hyper_bands.len());
    
    println!("\n→ Processing hyperspectral data...");
    let start = std::time::Instant::now();
    let hyper_result = process_bands_parallel(&pipeline, hyper_bands).unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Processed in {:.2?}", elapsed);
    
    // Write as GeoZarr (hyperspectral)
    let geozarr_path = temp.path().join("hyperspectral_output.zarr");
    
    println!("\n→ Writing GeoZarr (hyperspectral)...");
    let start = std::time::Instant::now();
    write_geozarr(geozarr_path.to_str().unwrap(), &hyper_result, &metadata).unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Written in {:.2?}", elapsed);
    println!("  Path: {:?}", geozarr_path);
    
    // Calculate total size
    let total_size: u64 = walkdir::WalkDir::new(&geozarr_path)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.metadata().ok())
        .filter(|m| m.is_file())
        .map(|m| m.len())
        .sum();
    
    println!("  Size: {} MB", total_size / 1_000_000);
    
    println!("\n✓ Processing complete!");
    println!("\nSummary:");
    println!("  Multispectral (7 bands):   COG format");
    println!("  Hyperspectral (100 bands): GeoZarr format");
    println!("  Total processing time:     {:.2?}", start.elapsed());
}
