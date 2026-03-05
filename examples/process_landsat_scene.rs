//! Process Landsat scene from GeoTIFF files

use acolite_rs::{Pipeline, ProcessingConfig};
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::io::{read_geotiff_band, write_geotiff_multiband};
use acolite_rs::parallel::process_bands_parallel;
use chrono::Utc;
use std::path::Path;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: Landsat Scene Processing from GeoTIFF\n");
    
    // Check if scene directory provided
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: {} <scene_dir>", args[0]);
        println!("\nExample scene structure:");
        println!("  scene_dir/");
        println!("    LC08_L1TP_..._B1.TIF");
        println!("    LC08_L1TP_..._B2.TIF");
        println!("    ...");
        println!("    LC08_L1TP_..._MTL.txt");
        println!("\nRunning with synthetic data instead...\n");
        run_synthetic();
        return;
    }
    
    let scene_dir = &args[1];
    if let Err(e) = process_scene(scene_dir) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn process_scene(scene_dir: &str) -> acolite_rs::Result<()> {
    println!("→ Processing Landsat scene: {}", scene_dir);
    
    // Find band files
    let band_files: Vec<_> = (1..=7)
        .filter_map(|i| {
            let pattern = format!("*_B{}.TIF", i);
            // In real implementation, would glob for files
            None::<String>
        })
        .collect();
    
    if band_files.is_empty() {
        return Err(acolite_rs::AcoliteError::Processing(
            "No band files found".to_string()
        ));
    }
    
    println!("  Found {} bands", band_files.len());
    
    // Setup metadata
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(30.0, 140.0);
    
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.0,
    };
    
    let pipeline = Pipeline::new(metadata.clone(), config);
    
    // Read bands - convert f64 to u16
    println!("\n→ Reading bands...");
    let wavelengths = [443.0, 482.0, 561.0, 655.0, 865.0, 1609.0, 2201.0];
    let bandwidths = [16.0, 60.0, 57.0, 37.0, 28.0, 85.0, 187.0];
    
    let mut bands = Vec::new();
    for (i, file) in band_files.iter().enumerate() {
        let band_f64 = read_geotiff_band(file, 1)?;
        println!("  Band {}: {}×{} pixels", i + 1, band_f64.data.ncols(), band_f64.data.nrows());
        
        // Convert f64 to u16
        let data_u16 = band_f64.data.mapv(|v| v as u16);
        let band_u16 = BandData::new(
            data_u16,
            wavelengths[i],
            bandwidths[i],
            format!("B{}", i + 1),
            band_f64.projection.clone(),
            band_f64.geotransform.clone(),
        );
        bands.push(band_u16);
    }
    
    // Process
    println!("\n→ Processing atmospheric correction...");
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands)?;
    let elapsed = start.elapsed();
    println!("  ✓ Processed in {:.2?}", elapsed);
    
    // Write output
    let output_path = Path::new(scene_dir).join("acolite_output.tif");
    println!("\n→ Writing output: {:?}", output_path);
    write_geotiff_multiband(output_path.to_str().unwrap(), &result, &metadata)?;
    
    println!("\n✓ Processing complete!");
    Ok(())
}

fn run_synthetic() {
    use ndarray::Array2;
    
    println!("→ Creating synthetic Landsat scene...");
    
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(30.0, -120.0);
    
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
    
    // Create synthetic bands
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);
    
    let wavelengths = [443.0, 482.0, 561.0, 655.0, 865.0, 1609.0, 2201.0];
    let bandwidths = [16.0, 60.0, 57.0, 37.0, 28.0, 85.0, 187.0];
    
    let bands: Vec<_> = (0..7)
        .map(|i| {
            let dn = 8000u16 + (i as u16 * 1000);
            
            BandData::new(
                Array2::from_elem((500, 500), dn),
                wavelengths[i],
                bandwidths[i],
                format!("B{}", i + 1),
                proj.clone(),
                geotrans.clone(),
            )
        })
        .collect();
    
    println!("  Created {} bands (500×500 pixels)", bands.len());
    
    println!("\n→ Processing atmospheric correction...");
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = start.elapsed();
    println!("  ✓ Processed in {:.2?}", elapsed);
    
    // Write output
    use tempfile::TempDir;
    let temp = TempDir::new().unwrap();
    let output_path = temp.path().join("synthetic_output.tif");
    
    println!("\n→ Writing output: {:?}", output_path);
    write_geotiff_multiband(output_path.to_str().unwrap(), &result, &metadata).unwrap();
    
    println!("  File size: {} bytes", std::fs::metadata(&output_path).unwrap().len());
    
    println!("\n✓ Synthetic processing complete!");
    println!("\nStatistics:");
    for band in &result {
        let mean = band.data.mean().unwrap();
        let std = band.data.std(0.0);
        println!("  {} ({:.1} nm): mean={:.4}, std={:.4}", 
                 band.name, band.wavelength, mean, std);
    }
}
