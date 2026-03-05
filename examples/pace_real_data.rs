//! Download and process real PACE OCI data from NASA DAAC

use acolite_rs::{
    Pipeline, ProcessingConfig, PaceOciSensor, 
    EarthdataAuth, search_pace_data, download_pace_file,
    write_geozarr
};
use acolite_rs::sensors::Sensor;
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::io::read_pace_band_simple;
use acolite_rs::parallel::process_bands_parallel;
use chrono::Utc;
use std::path::Path;
use ndarray::Array2;

fn downsample_array(arr: &Array2<f32>, target_h: usize, target_w: usize) -> Array2<f64> {
    let (h, w) = arr.dim();
    let step_h = (h / target_h).max(1);
    let step_w = (w / target_w).max(1);
    
    Array2::from_shape_fn((target_h, target_w), |(i, j)| {
        let src_i = (i * step_h).min(h - 1);
        let src_j = (j * step_w).min(w - 1);
        arr[[src_i, src_j]] as f64
    })
}

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: Real PACE Data Integration\n");
    
    // Load Earthdata credentials
    println!("→ Loading Earthdata credentials...");
    let auth = match EarthdataAuth::load() {
        Ok(auth) => {
            println!("  ✓ Loaded credentials for: {}", auth.username);
            auth
        }
        Err(e) => {
            eprintln!("✗ Failed to load credentials: {}", e);
            eprintln!("\nCreate ~/.easi-workflows-auth.conf with:");
            eprintln!("[earthdata]");
            eprintln!("    user: YOUR_EARTHDATA_USERNAME");
            eprintln!("    password: YOUR_EARTHDATA_PASSWORD");
            eprintln!("\nGet credentials at: https://urs.earthdata.nasa.gov/");
            std::process::exit(1);
        }
    };
    
    // Search for PACE data
    println!("\n→ Searching for PACE OCI data...");
    let bbox = [110.0, -45.0, 155.0, -10.0]; // Australia
    let start_date = "2024-02-01";
    let end_date = "2025-01-31";
    
    println!("  Area: Australia (110°E-155°E, 10°S-45°S)");
    println!("  Period: {} to {}", start_date, end_date);
    
    let urls = match search_pace_data(&bbox, start_date, end_date) {
        Ok(urls) => {
            if urls.is_empty() {
                println!("\n✗ No PACE data found for this region/time");
                println!("  Try different dates or region");
                std::process::exit(0);
            }
            println!("  ✓ Found {} granules", urls.len());
            urls
        }
        Err(e) => {
            eprintln!("\n✗ Search failed: {}", e);
            eprintln!("  Check network connection and CMR availability");
            std::process::exit(1);
        }
    };
    
    // Show available files
    println!("\nAvailable PACE files:");
    for (i, url) in urls.iter().enumerate().take(5) {
        let filename = url.split('/').last().unwrap_or("unknown");
        println!("  {}. {}", i + 1, filename);
    }
    
    // Download first file
    let download_url = &urls[0];
    let filename = download_url.split('/').last().unwrap_or("pace_data.nc");
    let output_dir = std::env::temp_dir().join("acolite_pace");
    std::fs::create_dir_all(&output_dir).unwrap();
    let nc_path = output_dir.join(filename);
    
    println!("\n→ Downloading PACE data via S3...");
    println!("  File: {}", filename);
    
    if nc_path.exists() {
        println!("  ✓ File already exists (using cached)");
    } else {
        // Get S3 credentials (not needed for TEA but validates auth)
        use acolite_rs::{get_s3_credentials, download_from_s3_with_token};
        
        match get_s3_credentials(&auth) {
            Ok(_) => {
                println!("  ✓ Authentication validated");
                
                // Convert HTTPS URL to S3 URL
                let s3_url = if download_url.contains("obdaac-tea.earthdatacloud.nasa.gov") {
                    download_url.replace(
                        "https://obdaac-tea.earthdatacloud.nasa.gov/",
                        "s3://"
                    )
                } else {
                    format!("s3://ob-cumulus-prod-public/{}", filename)
                };
                
                println!("  S3 path: {}", s3_url);
                
                match download_from_s3_with_token(&s3_url, nc_path.to_str().unwrap(), &auth) {
                    Ok(_) => {
                        println!("  ✓ Download complete");
                    }
                    Err(e) => {
                        eprintln!("\n✗ S3 download failed: {}", e);
                        std::process::exit(1);
                    }
                }
            }
            Err(e) => {
                eprintln!("\n✗ Failed to get S3 credentials: {}", e);
                eprintln!("  Trying direct HTTPS download...");
                
                match download_pace_file(download_url, nc_path.to_str().unwrap(), &auth) {
                    Ok(_) => println!("  ✓ Download complete"),
                    Err(e) => {
                        eprintln!("\n✗ Download failed: {}", e);
                        std::process::exit(1);
                    }
                }
            }
        }
    }
    
    // Process the data
    println!("\n→ Processing PACE OCI data...");
    process_pace_file(nc_path.to_str().unwrap());
}

fn process_pace_file(nc_path: &str) {
    use tempfile::TempDir;
    use acolite_rs::io::pace_nc_real::{get_pace_dimensions, read_pace_band, read_pace_geolocation};
    
    let sensor = PaceOciSensor;
    
    println!("  Getting dimensions...");
    let (num_bands, scans, pixels) = match get_pace_dimensions(nc_path) {
        Ok(dims) => {
            println!("  ✓ {} bands, {} scans, {} pixels", dims.0, dims.1, dims.2);
            dims
        }
        Err(e) => {
            println!("\n✗ Failed to read dimensions: {}", e);
            return;
        }
    };
    
    println!("  Reading geolocation...");
    let (lat, lon) = match read_pace_geolocation(nc_path) {
        Ok(geo) => {
            println!("  ✓ Geolocation: {} × {}", geo.0.shape()[0], geo.0.shape()[1]);
            println!("  Lat range: {:.2} to {:.2}", 
                geo.0.iter().filter(|x| x.is_finite()).fold(f32::INFINITY, |a, &b| a.min(b)),
                geo.0.iter().filter(|x| x.is_finite()).fold(f32::NEG_INFINITY, |a, &b| a.max(b)));
            println!("  Lon range: {:.2} to {:.2}",
                geo.1.iter().filter(|x| x.is_finite()).fold(f32::INFINITY, |a, &b| a.min(b)),
                geo.1.iter().filter(|x| x.is_finite()).fold(f32::NEG_INFINITY, |a, &b| a.max(b)));
            geo
        }
        Err(e) => {
            println!("\n✗ Failed to read geolocation: {}", e);
            return;
        }
    };
    
    // Read first 10 bands
    let num_to_read = 10.min(num_bands);
    println!("  Reading first {} bands...", num_to_read);
    
    let mut bands = Vec::new();
    for i in 0..num_to_read {
        match read_pace_band(nc_path, i) {
            Ok(band) => {
                let valid = band.iter().filter(|x| x.is_finite() && **x > 0.0).count();
                println!("    Band {}: {} valid pixels", i, valid);
                bands.push(band);
            }
            Err(e) => {
                println!("    Band {}: failed - {}", i, e);
                break;
            }
        }
    }
    
    if bands.is_empty() {
        println!("\n✗ No bands could be read");
        return;
    }
    
    println!("  ✓ Read {} bands ({} × {})", bands.len(), bands[0].shape()[0], bands[0].shape()[1]);
    
    let mut metadata = Metadata::new(sensor.name().to_string(), Utc::now());
    
    // Use actual geolocation center
    let center_lat = lat[[scans/2, pixels/2]];
    let center_lon = lon[[scans/2, pixels/2]];
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
    
    // Convert to BandData for processing (downsample for demo)
    let band_names = sensor.band_names();
    let bands_f64: Vec<_> = bands.iter()
        .enumerate()
        .take(10)
        .filter_map(|(i, band)| {
            let name = &band_names[i];
            let wl = sensor.wavelength(name)?;
            let bw = sensor.bandwidth(name)?;
            
            // Downsample to 200x200 for demo
            let downsampled = downsample_array(band, 200, 200);
            
            let proj = Projection::from_epsg(4326);
            let geotrans = GeoTransform::new(center_lon as f64, 0.01, center_lat as f64, -0.01);
            
            Some(BandData::new(downsampled, wl, bw, name.clone(), proj, geotrans))
        })
        .collect();
    
    if bands_f64.is_empty() {
        println!("\n✗ No bands for processing");
        return;
    }
    
    // Convert f64 to u16 for processing
    let bands: Vec<_> = bands_f64.iter()
        .map(|band| {
            use acolite_rs::core::BandData;
            let data_u16 = band.data.mapv(|v| (v * 10000.0) as u16);
            BandData::new(
                data_u16,
                band.wavelength,
                band.bandwidth,
                band.name.clone(),
                band.projection.clone(),
                band.geotransform.clone(),
            )
        })
        .collect();
    
    println!("  ✓ Read {} bands", bands.len());
    
    println!("\n→ Applying atmospheric correction...");
    let start = std::time::Instant::now();
    let result = match process_bands_parallel(&pipeline, bands) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("✗ Processing failed: {}", e);
            return;
        }
    };
    let elapsed = start.elapsed();
    
    println!("  ✓ Processed {} bands in {:.2?}", result.len(), elapsed);
    
    // Write to GeoZarr
    let temp = TempDir::new().unwrap();
    let output_path = temp.path().join("pace_real_output.zarr");
    
    println!("\n→ Writing to GeoZarr...");
    let start = std::time::Instant::now();
    match write_geozarr(output_path.to_str().unwrap(), &result, &metadata) {
        Ok(_) => {
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
        }
        Err(e) => {
            eprintln!("✗ Write failed: {}", e);
            return;
        }
    }
    
    // Show sample reflectances
    println!("\n→ Sample ocean reflectances:");
    for band in result.iter().take(10) {
        let mean = band.data.mean().unwrap();
        println!("  {:.0} nm: ρs = {:.4}", band.wavelength, mean);
    }
    
    println!("\n✓ Real PACE data processing complete!");
    println!("\nWorkflow:");
    println!("  1. ✓ Authenticated with NASA Earthdata");
    println!("  2. ✓ Searched CMR for PACE data");
    println!("  3. ✓ Downloaded NetCDF file");
    println!("  4. ✓ Read hyperspectral bands");
    println!("  5. ✓ Applied atmospheric correction");
    println!("  6. ✓ Wrote GeoZarr output");
}
