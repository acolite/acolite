//! Process PACE data at full resolution and output to GeoTIFF

use acolite_rs::{
    PaceOciSensor, Pipeline, ProcessingConfig,
};
use acolite_rs::sensors::Sensor;
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::io::pace_nc_real::{get_pace_dimensions, read_pace_band, read_pace_geolocation};
use chrono::Utc;
use std::path::Path;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: Full Resolution PACE Processing\n");
    
    let nc_path = "/tmp/acolite_testset/PACE_OCI.20240305T030038.L1B.V3.nc";
    
    if !Path::new(nc_path).exists() {
        eprintln!("✗ File not found: {}", nc_path);
        eprintln!("  Run download_testset first");
        std::process::exit(1);
    }
    
    println!("Input: {}\n", nc_path);
    
    // Get dimensions
    println!("→ Reading metadata...");
    let (num_bands, scans, pixels) = match get_pace_dimensions(nc_path) {
        Ok(dims) => {
            println!("  Dimensions: {} bands × {} scans × {} pixels", dims.0, dims.1, dims.2);
            println!("  Total pixels: {} million", (dims.1 * dims.2) / 1_000_000);
            dims
        }
        Err(e) => {
            eprintln!("✗ Failed: {}", e);
            std::process::exit(1);
        }
    };
    
    // Read geolocation
    println!("\n→ Reading geolocation...");
    let (lat, lon) = match read_pace_geolocation(nc_path) {
        Ok(geo) => {
            let lat_min = geo.0.iter().filter(|x| x.is_finite()).fold(f32::INFINITY, |a, &b| a.min(b));
            let lat_max = geo.0.iter().filter(|x| x.is_finite()).fold(f32::NEG_INFINITY, |a, &b| a.max(b));
            let lon_min = geo.1.iter().filter(|x| x.is_finite()).fold(f32::INFINITY, |a, &b| a.min(b));
            let lon_max = geo.1.iter().filter(|x| x.is_finite()).fold(f32::NEG_INFINITY, |a, &b| a.max(b));
            
            println!("  Latitude: {:.2}° to {:.2}°", lat_min, lat_max);
            println!("  Longitude: {:.2}° to {:.2}°", lon_min, lon_max);
            geo
        }
        Err(e) => {
            eprintln!("✗ Failed: {}", e);
            std::process::exit(1);
        }
    };
    
    // Process subset of bands at full resolution
    let bands_to_process = vec![10, 20, 30, 40, 50]; // Sample across spectrum
    println!("\n→ Processing {} bands at full resolution...", bands_to_process.len());
    
    let sensor = PaceOciSensor;
    let band_names = sensor.band_names();
    
    let center_lat = lat[[scans/2, pixels/2]];
    let center_lon = lon[[scans/2, pixels/2]];
    
    let mut metadata = Metadata::new(sensor.name().to_string(), Utc::now());
    metadata.set_geometry(center_lat as f64, center_lon as f64);
    
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
    pipeline.set_aot(0.12);
    
    let start = std::time::Instant::now();
    
    for (i, &band_idx) in bands_to_process.iter().enumerate() {
        print!("  Band {}/{}: ", i+1, bands_to_process.len());
        
        let band_start = std::time::Instant::now();
        
        // Read band
        let band_data = match read_pace_band(nc_path, band_idx) {
            Ok(data) => data,
            Err(e) => {
                println!("✗ Read failed: {}", e);
                continue;
            }
        };
        
        let read_time = band_start.elapsed();
        
        // Convert to u16 for processing (scale 0-1 to 0-10000)
        let band_u16 = band_data.mapv(|x: f32| (x.clamp(0.0, 1.0) * 10000.0) as u16);
        
        let name = &band_names[band_idx];
        let wl = sensor.wavelength(name).unwrap();
        let bw = sensor.bandwidth(name).unwrap();
        
        let proj = Projection::from_epsg(4326);
        let geotrans = GeoTransform::new(center_lon as f64, 0.01, center_lat as f64, -0.01);
        
        let band_obj = BandData::new(band_u16, wl, bw, name.clone(), proj, geotrans);
        
        // Process
        let proc_start = std::time::Instant::now();
        let processed = match pipeline.process_band(band_obj) {
            Ok(p) => p,
            Err(e) => {
                println!("✗ Process failed: {}", e);
                continue;
            }
        };
        let proc_time = proc_start.elapsed();
        
        // Write GeoTIFF
        let output_path = format!("/tmp/acolite_testset/pace_band_{:03}.tif", band_idx);
        let write_start = std::time::Instant::now();
        
        match acolite_rs::io::write_geotiff_band(&output_path, &processed, &metadata) {
            Ok(_) => {
                let write_time = write_start.elapsed();
                
                let file_size = std::fs::metadata(&output_path).unwrap().len();
                
                println!("✓ {} nm (read: {:.1}s, proc: {:.1}s, write: {:.1}s, size: {:.1} MB)",
                    wl as u32,
                    read_time.as_secs_f32(),
                    proc_time.as_secs_f32(),
                    write_time.as_secs_f32(),
                    file_size as f32 / 1_000_000.0
                );
            }
            Err(e) => {
                println!("✗ GeoTIFF write failed: {}", e);
            }
        }
    }
    
    let total_time = start.elapsed();
    
    println!("\n✓ Processing complete!");
    println!("  Total time: {:.1}s", total_time.as_secs_f32());
    println!("  Time per band: {:.1}s", total_time.as_secs_f32() / bands_to_process.len() as f32);
    println!("  Output: /tmp/acolite_testset/pace_band_*.tif");
}
