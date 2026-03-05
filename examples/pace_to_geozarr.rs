//! Convert PACE NetCDF directly to GeoZarr (fast path)

use acolite_rs::io::pace_zarr::{pace_to_geozarr, pace_to_geozarr_corrected};
use std::path::Path;

fn main() {
    env_logger::init();
    
    println!("ACOLITE-RS: PACE to GeoZarr Direct Conversion\n");
    
    let nc_path = "/tmp/acolite_testset/PACE_OCI.20240305T030038.L1B.V3.nc";
    
    if !Path::new(nc_path).exists() {
        eprintln!("✗ File not found: {}", nc_path);
        eprintln!("  Run download_testset first");
        std::process::exit(1);
    }
    
    println!("Input: {}\n", nc_path);
    
    // Test 1: Direct conversion (no AC)
    println!("→ Test 1: Direct conversion (TOA reflectance)");
    let zarr_toa = "/tmp/acolite_testset/pace_toa.zarr";
    let bands_subset = vec![0, 10, 20, 30, 40, 50, 60, 70, 80, 90]; // 10 bands
    
    let start = std::time::Instant::now();
    match pace_to_geozarr(nc_path, zarr_toa, &bands_subset) {
        Ok(_) => {
            let elapsed = start.elapsed();
            println!("  ✓ Complete in {:.1}s", elapsed.as_secs_f32());
            println!("  Output: {}", zarr_toa);
        }
        Err(e) => {
            eprintln!("  ✗ Failed: {}", e);
        }
    }
    
    // Test 2: With atmospheric correction
    println!("\n→ Test 2: With atmospheric correction");
    let zarr_rhos = "/tmp/acolite_testset/pace_rhos.zarr";
    let bands_full = (0..20).collect::<Vec<_>>(); // 20 bands
    
    let start = std::time::Instant::now();
    match pace_to_geozarr_corrected(nc_path, zarr_rhos, &bands_full) {
        Ok(_) => {
            let elapsed = start.elapsed();
            println!("  ✓ Complete in {:.1}s", elapsed.as_secs_f32());
            println!("  Output: {}", zarr_rhos);
        }
        Err(e) => {
            eprintln!("  ✗ Failed: {}", e);
        }
    }
    
    // Test 3: All blue bands
    println!("\n→ Test 3: All blue bands (119 bands)");
    let zarr_all = "/tmp/acolite_testset/pace_all_blue.zarr";
    let bands_all = (0..119).collect::<Vec<_>>();
    
    let start = std::time::Instant::now();
    match pace_to_geozarr(nc_path, zarr_all, &bands_all) {
        Ok(_) => {
            let elapsed = start.elapsed();
            println!("  ✓ Complete in {:.1}s", elapsed.as_secs_f32());
            println!("  Output: {}", zarr_all);
            println!("  Throughput: {:.1} bands/sec", 119.0 / elapsed.as_secs_f32());
        }
        Err(e) => {
            eprintln!("  ✗ Failed: {}", e);
        }
    }
    
    println!("\n✓ GeoZarr conversion complete!");
    println!("\nPerformance comparison:");
    println!("  Old method: ~6.5 sec/band (Python bridge per band)");
    println!("  New method: Direct NetCDF→Zarr (batch processing)");
    println!("\nVerify with:");
    println!("  python3 -c 'import zarr; z=zarr.open(\"{}\"); print(z.tree())'", zarr_all);
}
