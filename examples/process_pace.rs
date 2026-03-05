//! Process PACE OCI data and output to GeoZarr

use acolite_rs::{Pipeline, ProcessingConfig, PaceOciSensor, write_geozarr};
use acolite_rs::sensors::Sensor;
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::parallel::process_bands_parallel;
use ndarray::Array2;
use chrono::Utc;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: PACE OCI Processing → GeoZarr\n");
    
    let sensor = PaceOciSensor;
    println!("Sensor: {}", sensor.name());
    println!("Bands: {}", sensor.band_names().len());
    
    // Check for real PACE data
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 {
        println!("\nNote: Real PACE NetCDF-4 reader not yet implemented");
        println!("Download PACE data from: https://oceandata.sci.gsfc.nasa.gov/");
        println!("Look for: PACE_OCI.*.L1B.*.nc files\n");
    }
    
    println!("Running with synthetic hyperspectral data...\n");
    process_synthetic();
}

fn process_synthetic() {
    use tempfile::TempDir;
    
    let sensor = PaceOciSensor;
    let band_names = sensor.band_names();
    let nbands = band_names.len();
    
    println!("→ Creating synthetic PACE OCI data...");
    println!("  Bands: {} hyperspectral", nbands);
    println!("  Wavelength range: 315-2260 nm");
    
    let mut metadata = Metadata::new(sensor.name().to_string(), Utc::now());
    metadata.set_geometry(35.0, -75.0); // Atlantic Ocean
    
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.5,
    };
    
    let mut pipeline = Pipeline::new(metadata.clone(), config);
    pipeline.set_aot(0.12); // Typical ocean AOT
    
    // Create hyperspectral bands
    let proj = Projection::from_epsg(4326); // Geographic for ocean
    let geotrans = GeoTransform::new(-75.0, 0.01, 35.0, -0.01); // 1km ~ 0.01 deg
    
    println!("\n→ Generating {} hyperspectral bands...", nbands);
    let bands: Vec<_> = band_names.iter()
        .enumerate()
        .map(|(i, name)| {
            let wl = sensor.wavelength(name).unwrap();
            let bw = sensor.bandwidth(name).unwrap();
            
            // Simulate ocean reflectance spectrum
            let base_dn = if wl < 500.0 {
                3000u16 + (i as u16 * 10) // Lower in blue
            } else if wl < 700.0 {
                4000u16 + (i as u16 * 5) // Peak in green
            } else {
                2000u16 + (i as u16 * 3) // Low in red/NIR
            };
            
            BandData::new(
                Array2::from_elem((200, 200), base_dn),
                wl,
                bw,
                name.clone(),
                proj.clone(),
                geotrans.clone(),
            )
        })
        .collect();
    
    println!("  Created {} bands (200×200 pixels)", bands.len());
    
    println!("\n→ Processing atmospheric correction...");
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = start.elapsed();
    
    println!("  ✓ Processed {} bands in {:.2?}", result.len(), elapsed);
    println!("  Throughput: {:.0} bands/sec", result.len() as f64 / elapsed.as_secs_f64());
    
    // Write to GeoZarr
    let temp = TempDir::new().unwrap();
    let output_path = temp.path().join("pace_oci_output.zarr");
    
    println!("\n→ Writing to GeoZarr...");
    let start = std::time::Instant::now();
    write_geozarr(output_path.to_str().unwrap(), &result, &metadata).unwrap();
    let elapsed = start.elapsed();
    
    println!("  ✓ Written in {:.2?}", elapsed);
    println!("  Path: {:?}", output_path);
    
    // Calculate size
    let total_size: u64 = walkdir::WalkDir::new(&output_path)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.metadata().ok())
        .filter(|m| m.is_file())
        .map(|m| m.len())
        .sum();
    
    println!("  Size: {:.1} MB", total_size as f64 / 1_000_000.0);
    println!("  Compression: {:.1}%", 
             100.0 * (1.0 - total_size as f64 / (result.len() * 200 * 200 * 8) as f64));
    
    // Show sample wavelengths and reflectances
    println!("\n→ Sample ocean reflectances:");
    for (i, band) in result.iter().enumerate().step_by(20) {
        let mean = band.data.mean().unwrap();
        println!("  {:.0} nm: ρs = {:.4}", band.wavelength, mean);
    }
    
    println!("\n✓ PACE OCI processing complete!");
    println!("\nGeoZarr Structure:");
    println!("  pace_oci_output.zarr/");
    println!("    ├── .zattrs (GeoZarr metadata)");
    println!("    ├── data/ ({} bands)", result.len());
    println!("    │   ├── .zarray");
    println!("    │   └── *.0.0 (chunks)");
    println!("    └── wavelengths/");
    println!("        ├── .zarray");
    println!("        └── 0 (wavelength array)");
    
    println!("\nTo process real PACE data:");
    println!("  1. Download from https://oceandata.sci.gsfc.nasa.gov/");
    println!("  2. Search for PACE_OCI L1B files");
    println!("  3. Implement NetCDF-4 reader");
    println!("  4. Run: cargo run --example process_pace <file.nc>");
}
