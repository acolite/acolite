//! Download and process PACE OCI L1B data from OBDAAC
//!
//! Usage:
//!   # Search by location + date range
//!   cargo run --example process_pace --features netcdf -- --lat 36.0 --lon -75.5 --start 2024-07-01 --end 2024-07-02
//!
//!   # Process a local file
//!   cargo run --example process_pace --features netcdf -- --file /path/to/PACE_OCI.*.L1B.*.nc
//!
//!   # Search + download + process
//!   cargo run --example process_pace --features netcdf -- --lat 36.0 --lon -75.5 --start 2024-07-01 --end 2024-07-02 --download

use acolite_rs::{
    Credentials, Pipeline, ProcessingConfig, write_auto,
    search_pace_l1b, CmrGranule,
};
use acolite_rs::loader::source::download::download_file;
use std::path::{Path, PathBuf};

fn main() {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    println!("ACOLITE-RS: PACE OCI Processing\n");

    let file = get_arg(&args, "--file");
    let lat = get_arg(&args, "--lat").map(|s| s.parse::<f64>().expect("invalid --lat"));
    let lon = get_arg(&args, "--lon").map(|s| s.parse::<f64>().expect("invalid --lon"));
    let start = get_arg(&args, "--start");
    let end = get_arg(&args, "--end");
    let do_download = args.contains(&"--download".to_string());
    let output_dir = get_arg(&args, "--output").unwrap_or_else(|| "/tmp/acolite_pace".into());

    // Determine input file(s)
    let input_files: Vec<PathBuf> = if let Some(f) = file {
        vec![PathBuf::from(f)]
    } else if let (Some(lat), Some(lon), Some(start), Some(end)) = (lat, lon, start.as_deref(), end.as_deref()) {
        println!("→ Searching OBDAAC for PACE OCI L1B...");
        println!("  Location: ({}, {})", lat, lon);
        println!("  Date range: {} to {}", start, end);

        let granules = search_pace_l1b(lat, lon, start, end)
            .expect("CMR search failed");

        if granules.is_empty() {
            println!("  No granules found.");
            return;
        }

        println!("  Found {} granule(s):", granules.len());
        for g in &granules {
            println!("    {} → {}", g.id, g.urls.first().unwrap_or(&"(no url)".into()));
        }

        if do_download {
            download_granules(&granules, &output_dir)
        } else {
            println!("\n  Add --download to download and process.");
            println!("  Or use --file <path> to process a local file.");
            return;
        }
    } else {
        eprintln!("Usage:");
        eprintln!("  --file <path>                    Process local PACE L1B file");
        eprintln!("  --lat <lat> --lon <lon>          Search location");
        eprintln!("  --start <date> --end <date>      Date range (YYYY-MM-DD)");
        eprintln!("  --download                       Download found granules");
        eprintln!("  --output <dir>                   Output directory (default: /tmp/acolite_pace)");
        return;
    };

    // Process each file
    for input in &input_files {
        println!("\n→ Processing {:?}", input);
        process_pace_file(input, &output_dir);
    }
}

fn download_granules(granules: &[CmrGranule], output_dir: &str) -> Vec<PathBuf> {
    let creds = Credentials::load().expect(
        "EarthData credentials required. Set EARTHDATA_u/EARTHDATA_p env vars or configure ~/.netrc"
    );
    println!("\n→ Downloading with EarthData credentials ({:?})", creds.source);

    std::fs::create_dir_all(output_dir).expect("create output dir");
    let mut files = Vec::new();

    for g in granules {
        if let Some(url) = g.urls.first() {
            let filename = url.rsplit('/').next().unwrap_or(&g.id);
            let dest = PathBuf::from(output_dir).join(filename);
            println!("  Downloading {} ...", filename);
            download_file(url, &dest, &creds).expect("download failed");
            files.push(dest);
        }
    }
    files
}

#[cfg(feature = "netcdf")]
fn process_pace_file(path: &Path, output_dir: &str) {
    use acolite_rs::load_pace_l1b;
    use acolite_rs::parallel::process_bands_parallel_f32;

    let start = std::time::Instant::now();
    let scene = load_pace_l1b(path, None).expect("Failed to load PACE L1B");
    let load_time = start.elapsed();

    let (nrows, ncols) = (scene.lat.nrows(), scene.lat.ncols());
    println!("  Sensor: {}", scene.metadata.sensor);
    println!("  Bands: {}", scene.bands.len());
    println!("  Size: {}×{}", nrows, ncols);
    println!("  Loaded in {:.2?}", load_time);

    // Show wavelength range
    if let (Some(first), Some(last)) = (scene.bands.first(), scene.bands.last()) {
        println!("  Wavelength range: {:.0}–{:.0} nm", first.wavelength, last.wavelength);
    }

    // Process atmospheric correction
    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 2.5,
    };
    let mut pipeline = Pipeline::new(scene.metadata.clone(), config);
    pipeline.set_aot(0.12);

    println!("\n→ Atmospheric correction...");
    let start = std::time::Instant::now();
    let result = process_bands_parallel_f32(&pipeline, scene.bands).expect("AC failed");
    let ac_time = start.elapsed();
    let mpx = (nrows * ncols * result.len()) as f64 / ac_time.as_secs_f64() / 1e6;
    println!("  ✓ {} bands in {:.2?} ({:.1} Mpx/s)", result.len(), ac_time, mpx);

    // Write output — auto-selects GeoZarr for hyperspectral (>50 bands)
    std::fs::create_dir_all(output_dir).expect("create output dir");
    let stem = path.file_stem().unwrap_or_default().to_string_lossy();
    let out_path = format!("{}/{}_corrected", output_dir, stem);

    println!("\n→ Writing output (auto-format)...");
    let start = std::time::Instant::now();
    write_auto(&out_path, &result, &scene.metadata).expect("Write failed");
    println!("  ✓ Written in {:.2?}", start.elapsed());

    // Summary
    println!("\n→ Sample reflectances:");
    for band in result.iter().step_by(result.len().max(1) / 10) {
        let mean = band.data.mean().unwrap_or(0.0);
        println!("  {:.0} nm: ρs = {:.4}", band.wavelength, mean);
    }
    println!("\n✓ PACE OCI processing complete!");
}

#[cfg(not(feature = "netcdf"))]
fn process_pace_file(_path: &Path, _output_dir: &str) {
    eprintln!("Error: netcdf feature required. Build with: cargo run --features netcdf --example process_pace");
}

fn get_arg(args: &[String], flag: &str) -> Option<String> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
