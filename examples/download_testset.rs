//! Download and process multi-sensor test set over Australia

use acolite_rs::{
    EarthdataAuth, search_landsat_data, search_pace_data,
    download_from_s3_with_token, get_s3_credentials,
};
use std::path::PathBuf;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: Multi-Sensor Test Set - Australia\n");
    
    // Load credentials
    let auth = match EarthdataAuth::load() {
        Ok(auth) => {
            println!("✓ Loaded credentials: {}", auth.username);
            auth
        }
        Err(e) => {
            eprintln!("✗ Failed to load credentials: {}", e);
            std::process::exit(1);
        }
    };
    
    // Validate auth
    match get_s3_credentials(&auth) {
        Ok(_) => println!("✓ Authentication validated\n"),
        Err(e) => {
            eprintln!("✗ Authentication failed: {}", e);
            std::process::exit(1);
        }
    }
    
    // Australia bounding box
    let bbox = [110.0, -45.0, 155.0, -10.0]; // West, South, East, North
    let start_date = "2024-01-01";
    let end_date = "2024-12-31";
    
    let output_dir = PathBuf::from("/tmp/acolite_testset");
    std::fs::create_dir_all(&output_dir).unwrap();
    
    println!("Region: Australia ({:.1}°E-{:.1}°E, {:.1}°S-{:.1}°S)", 
        bbox[0], bbox[2], -bbox[3], -bbox[1]);
    println!("Period: {} to {}", start_date, end_date);
    println!("Output: {}\n", output_dir.display());
    
    // Search Landsat
    println!("→ Searching Landsat 8/9...");
    let landsat_urls = match search_landsat_data(&bbox, start_date, end_date) {
        Ok(urls) => {
            println!("  ✓ Found {} scenes", urls.len());
            urls
        }
        Err(e) => {
            println!("  ✗ Search failed: {}", e);
            Vec::new()
        }
    };
    
    // Search PACE
    println!("\n→ Searching PACE OCI...");
    let pace_urls = match search_pace_data(&bbox, start_date, end_date) {
        Ok(urls) => {
            println!("  ✓ Found {} granules", urls.len());
            urls
        }
        Err(e) => {
            println!("  ✗ Search failed: {}", e);
            Vec::new()
        }
    };
    
    // Download samples
    println!("\n→ Downloading test samples...");
    
    // Download 1 Landsat scene
    if !landsat_urls.is_empty() {
        let url = &landsat_urls[0];
        let filename = url.split('/').last().unwrap_or("landsat.tar");
        println!("\n  Landsat: {}", filename);
        
        let s3_url = convert_to_s3(url);
        let output_path = output_dir.join(filename);
        
        if output_path.exists() {
            println!("    ✓ Already downloaded");
        } else {
            match download_from_s3_with_token(&s3_url, output_path.to_str().unwrap(), &auth) {
                Ok(_) => println!("    ✓ Downloaded"),
                Err(e) => println!("    ✗ Failed: {}", e),
            }
        }
    }
    
    // Download 1 PACE granule
    if !pace_urls.is_empty() {
        let url = &pace_urls[0];
        let filename = url.split('/').last().unwrap_or("pace.nc");
        println!("\n  PACE: {}", filename);
        
        let s3_url = convert_to_s3(url);
        let output_path = output_dir.join(filename);
        
        if output_path.exists() {
            println!("    ✓ Already downloaded");
        } else {
            match download_from_s3_with_token(&s3_url, output_path.to_str().unwrap(), &auth) {
                Ok(_) => println!("    ✓ Downloaded"),
                Err(e) => println!("    ✗ Failed: {}", e),
            }
        }
    }
    
    println!("\n✓ Test set ready: {}", output_dir.display());
    println!("\nNext steps:");
    println!("  1. Process Landsat with full resolution");
    println!("  2. Process PACE with full resolution");
    println!("  3. Output to COG format");
    println!("  4. Compare performance vs Python ACOLITE");
}

fn convert_to_s3(url: &str) -> String {
    if url.starts_with("s3://") {
        return url.to_string();
    }
    
    // Convert HTTPS URLs to S3
    if url.contains("lpdaac.earthdatacloud.nasa.gov") {
        url.replace("https://data.lpdaac.earthdatacloud.nasa.gov/", "s3://")
    } else if url.contains("obdaac-tea.earthdatacloud.nasa.gov") {
        url.replace("https://obdaac-tea.earthdatacloud.nasa.gov/", "s3://")
    } else {
        url.to_string()
    }
}
