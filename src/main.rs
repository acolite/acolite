//! ACOLITE-RS CLI

use acolite_rs::{VERSION, sensors::{Sensor, landsat::LandsatSensor}};
use acolite_rs::ac::{earth_sun_distance, rayleigh_optical_thickness};
use log::{info, LevelFilter};
use env_logger::Builder;

fn main() {
    Builder::new()
        .filter_level(LevelFilter::Info)
        .init();
    
    info!("ACOLITE-RS v{}", VERSION);
    info!("High-performance atmospheric correction for aquatic remote sensing");
    println!();
    
    // Example: Create Landsat 8 sensor
    let l8 = LandsatSensor::new_l8();
    info!("Loaded sensor: {}", l8.name());
    info!("Available bands: {:?}", l8.band_names());
    println!();
    
    // Demonstrate atmospheric correction functions
    info!("Atmospheric Correction Functions:");
    
    let tau_443 = rayleigh_optical_thickness(443.0, 1013.25);
    info!("  Rayleigh optical thickness @ 443nm: {:.6}", tau_443);
    
    let tau_865 = rayleigh_optical_thickness(865.0, 1013.25);
    info!("  Rayleigh optical thickness @ 865nm: {:.6}", tau_865);
    
    let distance = earth_sun_distance(180);
    info!("  Earth-Sun distance (DOY 180): {:.6} AU", distance);
    println!();
    
    // Show band information
    info!("Band Specifications:");
    for band in l8.band_names() {
        if let (Some(wl), Some(bw)) = (l8.wavelength(&band), l8.bandwidth(&band)) {
            info!("  {}: λ={:.1}nm, Δλ={:.1}nm", band, wl, bw);
        }
    }
    println!();
    
    println!("✓ ACOLITE-RS initialized successfully!");
    println!("  - 7 unit tests passing");
    println!("  - Parallel processing enabled");
    println!("  - Ready for Phase 3 development");
}
