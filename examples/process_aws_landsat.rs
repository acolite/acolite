//! Download and process Landsat from AWS S3 (Requester Pays)
//! Usage: AWS_PROFILE=your-profile cargo run --example process_aws_landsat

use acolite_rs::{
    Metadata, BandData, Projection, GeoTransform,
    Pipeline, ProcessingConfig, write_cog, cog_available,
    aws_profile,
    loader::source::download::download_s3_requester_pays,
    loader::geotiff::read_geotiff_band,
    parallel::process_bands_parallel,
};
use std::path::PathBuf;
use std::time::Instant;
use chrono::Utc;
use ndarray::Array2;

const SCENE_ID: &str = "LC08_L1TP_014027_20181208_20200830_02_T1";
const S3_BASE: &str = "s3://usgs-landsat/collection02/level-1/standard/oli-tirs/2018/014/027";

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    println!("ACOLITE-RS: AWS S3 Landsat Processing\n");

    let profile = aws_profile().unwrap_or_else(|| {
        eprintln!("Error: AWS_PROFILE not set");
        eprintln!("Usage: AWS_PROFILE=your-profile cargo run --example process_aws_landsat");
        std::process::exit(1);
    });
    println!("Scene: {}\n", SCENE_ID);

    let cache = PathBuf::from("/tmp/acolite_aws");
    std::fs::create_dir_all(&cache).unwrap();

    // Download bands
    println!("→ Downloading bands...");
    let wavelengths = [(1,443.0,16.0),(2,482.0,60.0),(3,561.0,57.0),(4,655.0,37.0),
                       (5,865.0,28.0),(6,1609.0,85.0),(7,2201.0,187.0)];
    let mut band_files = Vec::new();
    for &(num, _, _) in &wavelengths {
        let s3 = format!("{}/{}/{}_B{}.TIF", S3_BASE, SCENE_ID, SCENE_ID, num);
        let local = cache.join(format!("B{}.TIF", num));
        download_s3_requester_pays(&s3, &local, &profile).unwrap();
        band_files.push(local);
    }

    // Read bands
    println!("\n→ Reading bands...");
    let start = Instant::now();
    let bands: Vec<_> = band_files.iter().enumerate().map(|(i, path)| {
        let mut b = read_geotiff_band(path).unwrap();
        b.wavelength = wavelengths[i].1;
        b.bandwidth = wavelengths[i].2;
        b.name = format!("B{}", wavelengths[i].0);
        b
    }).collect();
    println!("  ✓ {} bands in {:.2}s ({}×{})", bands.len(), start.elapsed().as_secs_f64(),
        bands[0].width(), bands[0].height());

    // Process
    let mut metadata = Metadata::new("LANDSAT_8_OLI".into(), Utc::now());
    metadata.set_geometry(30.0, 135.0);
    let pipeline = Pipeline::new(metadata.clone(), ProcessingConfig::default());

    println!("\n→ Processing AC...");
    let start = Instant::now();
    let corrected = process_bands_parallel(&pipeline, bands).unwrap();
    let px = corrected[0].width() * corrected[0].height() * corrected.len();
    println!("  ✓ {:.2}s ({:.1} Mpx/s)", start.elapsed().as_secs_f64(),
        px as f64 / start.elapsed().as_secs_f64() / 1e6);

    // Write
    let out = PathBuf::from("/tmp/acolite_output/landsat_corrected.tif");
    std::fs::create_dir_all(out.parent().unwrap()).unwrap();
    println!("\n→ Writing COG...");
    let start = Instant::now();
    write_cog(out.to_str().unwrap(), &corrected, &metadata).unwrap();
    let sz = std::fs::metadata(&out).unwrap().len();
    println!("  ✓ {:.1} MB in {:.2}s ({})", sz as f64/1e6, start.elapsed().as_secs_f64(),
        if cog_available() { "COG" } else { "GeoTIFF" });
    println!("\n✓ Output: {}", out.display());
}
