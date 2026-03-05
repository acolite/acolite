//! Sentinel-3 OLCI processing example

use acolite_rs::{Pipeline, ProcessingConfig, Sentinel3Sensor};
use acolite_rs::sensors::Sensor;
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::ac::LutManager;
use acolite_rs::parallel::process_bands_parallel;
use ndarray::Array2;
use chrono::Utc;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: Sentinel-3 OLCI Processing\n");
    
    // Setup LUT manager and download if needed
    let lut_manager = LutManager::new(None);
    let lut_url = "https://raw.githubusercontent.com/acolite/acolite_luts/main/ACOLITE/Shared/rayleigh_lut.json";
    
    println!("→ Checking Rayleigh LUT...");
    if let Err(e) = lut_manager.download_lut("rayleigh_lut.json", lut_url) {
        println!("  Warning: {}", e);
    }
    
    let sensor = Sentinel3Sensor;
    println!("✓ Sensor: {}", sensor.name());
    println!("  Bands: {}", sensor.band_names().len());
    
    // Setup metadata
    let mut metadata = Metadata::new(sensor.name().to_string(), Utc::now());
    metadata.set_geometry(35.0, 140.0);
    
    // Processing config
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.5,
    };
    
    let mut pipeline = Pipeline::new(metadata, config);
    pipeline.set_aot(0.15);
    
    // Create test bands (first 7 visible bands)
    let proj = Projection::from_epsg(32633);
    let geotrans = GeoTransform::new(300000.0, 300.0, 4500000.0, -300.0);
    
    let bands: Vec<_> = ["Oa01", "Oa02", "Oa03", "Oa04", "Oa05", "Oa06", "Oa07"]
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let wl = sensor.wavelength(name).unwrap();
            let bw = sensor.bandwidth(name).unwrap();
            let dn = 800u16 + (i as u16 * 200);
            
            BandData::new(
                Array2::from_elem((200, 200), dn),
                wl,
                bw,
                name.to_string(),
                proj.clone(),
                geotrans.clone(),
            )
        })
        .collect();
    
    println!("\n→ Processing {} bands (200×200 pixels)...", bands.len());
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = start.elapsed();
    
    println!("✓ Processed in {:.2?}\n", elapsed);
    
    for band in &result {
        let mean = band.data.mean().unwrap();
        println!("  {} ({:.1} nm): ρs = {:.4}", band.name, band.wavelength, mean);
    }
}
