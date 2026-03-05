//! ACOLITE-RS CLI

use acolite_rs::{VERSION, sensors::{Sensor, landsat::LandsatSensor}};
use log::info;

fn main() {
    env_logger::init();
    
    info!("ACOLITE-RS v{}", VERSION);
    info!("High-performance atmospheric correction for aquatic remote sensing");
    
    // Example: Create Landsat 8 sensor
    let l8 = LandsatSensor::new_l8();
    info!("Loaded sensor: {}", l8.name());
    info!("Available bands: {:?}", l8.band_names());
    
    println!("ACOLITE-RS initialized successfully!");
}
