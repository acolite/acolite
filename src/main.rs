//! ACOLITE-RS CLI

use acolite_rs::{VERSION, sensors::{Sensor, landsat::LandsatSensor}};
use acolite_rs::ac::{earth_sun_distance, rayleigh_optical_thickness};

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    log::info!("ACOLITE-RS v{}", VERSION);

    let l8 = LandsatSensor::new_l8();
    log::info!("Sensor: {} ({} bands)", l8.name(), l8.band_names().len());

    for band in l8.band_names() {
        if let (Some(wl), Some(bw)) = (l8.wavelength(&band), l8.bandwidth(&band)) {
            log::info!("  {}: λ={:.1}nm Δλ={:.1}nm τ_r={:.4}", 
                band, wl, bw, rayleigh_optical_thickness(wl, 1013.25));
        }
    }

    log::info!("Earth-Sun distance (DOY 180): {:.6} AU", earth_sun_distance(180));
    println!("\n✓ ACOLITE-RS ready");
}
