//! Sentinel-3 OLCI processing example
//!
//! Usage:
//!   cargo run --release --features netcdf --example process_sentinel3 -- --scene /path/to/S3A_OL_1_EFR.SEN3
//!   cargo run --release --features full-io --example process_sentinel3 -- --scene /path/to/S3A_OL_1_EFR.SEN3 --data-dir ./data
//!   cargo run --release --features full-io --example process_sentinel3 -- --scene /path/to/S3A_OL_1_EFR.SEN3 --limit -35.0,138.0,-34.5,138.7

use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::loader::sentinel3::{OlciScene, OLCI_BANDS};
use acolite_rs::parallel::process_bands_parallel;
use acolite_rs::pipeline::{Pipeline, ProcessingConfig};
use acolite_rs::sensors::{Sensor, Sentinel3Sensor};
use chrono::Utc;
use ndarray::Array2;
use std::collections::HashMap;

fn main() {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    let scene_dir = get_arg(&args, "--scene");
    let output_dir = get_arg(&args, "--output").unwrap_or_else(|| "/tmp/s3_output".into());
    let limit: Option<[f64; 4]> = get_arg(&args, "--limit").map(|s| {
        let v: Vec<f64> = s.split(',').map(|x| x.parse().expect("invalid --limit")).collect();
        assert!(v.len() == 4, "--limit needs 4 values: south,west,north,east");
        [v[0], v[1], v[2], v[3]]
    });
    let data_dir = get_arg(&args, "--data-dir");

    println!("ACOLITE-RS: Sentinel-3 OLCI Processing\n");

    if let Some(ref _dir) = scene_dir {
        #[cfg(feature = "netcdf")]
        {
            process_real(std::path::Path::new(_dir), &output_dir, limit.as_ref(), data_dir.as_deref());
        }
        #[cfg(not(feature = "netcdf"))]
        {
            println!("Real scene processing requires --features netcdf");
            println!("Falling back to synthetic demo...\n");
            process_synthetic(&output_dir);
        }
    } else {
        println!("No --scene provided, running synthetic demo\n");
        process_synthetic(&output_dir);
    }
}

fn process_synthetic(output_dir: &str) {
    let sensor = Sentinel3Sensor;
    println!("Sensor: {} ({} bands)", sensor.name(), sensor.band_names().len());

    let mut metadata = Metadata::new(sensor.name().to_string(), Utc::now());
    metadata.set_geometry(35.0, 140.0);
    metadata.view_zenith = Some(15.0);

    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        apply_glint: false,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.5,
        pressure: 1013.25,
        wind: 2.0,
    };

    let mut pipeline = Pipeline::new(metadata, config);
    pipeline.set_aot(0.15);

    let proj = Projection::from_epsg(4326);
    let geotrans = GeoTransform::new(3.0, 0.003, 51.0, -0.003);

    let size = 300;
    let bands: Vec<BandData<u16>> = OLCI_BANDS
        .iter()
        .map(|(name, wl, bw, _)| {
            let base = (1200.0 * (400.0 / wl).powf(2.0) + 300.0) as u16;
            BandData::new(
                Array2::from_elem((size, size), base),
                *wl, *bw, name.to_string(),
                proj.clone(), geotrans.clone(),
            )
        })
        .collect();

    println!("\n→ Processing {} bands ({}×{} pixels)...", bands.len(), size, size);
    let start = std::time::Instant::now();
    let result = process_bands_parallel(&pipeline, bands).unwrap();
    let elapsed = start.elapsed();

    let total_pixels = size * size * result.len();
    let mpx_per_sec = total_pixels as f64 / elapsed.as_secs_f64() / 1e6;

    println!("✓ Processed in {:.2?} ({:.1} Mpx/s)\n", elapsed, mpx_per_sec);

    println!("{:<6} {:>8} {:>10}", "Band", "λ (nm)", "ρs mean");
    println!("{}", "-".repeat(28));
    for band in &result {
        let mean = band.data.iter().filter(|v| v.is_finite()).sum::<f64>()
            / band.data.len() as f64;
        println!("{:<6} {:>8.1} {:>10.4}", band.name, band.wavelength, mean);
    }
}

#[cfg(feature = "netcdf")]
fn process_real(
    scene_dir: &std::path::Path,
    output_dir: &str,
    limit: Option<&[f64; 4]>,
    data_dir: Option<&str>,
) {
    use acolite_rs::loader::sentinel3::load_olci_scene;
    use acolite_rs::loader::sentinel3_l2r::{process_olci_l2r, OlciProcessingConfig};

    println!("Loading scene: {}", scene_dir.display());
    if let Some(lim) = limit {
        println!("  Limit: [{:.2}, {:.2}, {:.2}, {:.2}]", lim[0], lim[1], lim[2], lim[3]);
    }

    let mut scene = match load_olci_scene(scene_dir, limit) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Failed to load scene: {}", e);
            return;
        }
    };

    println!("  Sensor: {}", scene.sensor);
    println!("  Shape: {}×{}", scene.data_shape.0, scene.data_shape.1);
    println!("  Bands loaded: {}", scene.radiance.len());

    let mut olci_config = OlciProcessingConfig::default();
    if let Some(dd) = data_dir {
        olci_config.data_dir = Some(std::path::PathBuf::from(dd));
    } else {
        // Try to find data dir automatically
        let mut p = std::env::current_dir().unwrap_or_default();
        for _ in 0..5 {
            let d = p.join("data");
            if d.join("RSR").exists() && d.join("LUT").exists() {
                olci_config.data_dir = Some(d);
                break;
            }
            if !p.pop() { break; }
        }
    }

    if olci_config.data_dir.is_some() {
        println!("  Data dir: {:?} (full DSF)", olci_config.data_dir.as_ref().unwrap());
    } else {
        println!("  No data dir found — using simplified Rayleigh correction");
    }

    let pipeline_config = ProcessingConfig::default();

    println!("\n→ Running atmospheric correction...");
    let start = std::time::Instant::now();
    let result = match process_olci_l2r(&mut scene, &olci_config, &pipeline_config) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Processing failed: {}", e);
            return;
        }
    };
    let elapsed = start.elapsed();

    println!("✓ Processed in {:.2?}", elapsed);
    println!("  AOT: {:.4} (model: {})", result.aot, result.model_name);

    println!("\n{:<6} {:>8} {:>10}", "Band", "λ (nm)", "ρs mean");
    println!("{}", "-".repeat(28));
    for (name, wl, _, _) in &OLCI_BANDS {
        if let Some(data) = result.rhos.get(&name.to_string()) {
            let valid: Vec<f64> = data.iter().filter(|v| v.is_finite() && **v > 0.0).copied().collect();
            let mean = if valid.is_empty() { f64::NAN } else { valid.iter().sum::<f64>() / valid.len() as f64 };
            println!("{:<6} {:>8.1} {:>10.4}", name, wl, mean);
        }
    }
}

fn get_arg(args: &[String], flag: &str) -> Option<String> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
