//! Process Sentinel-2 scene with full LUT-based DSF atmospheric correction
//!
//! Usage:
//!   cargo run --release --features full-io --example process_sentinel2_ac -- --file /path/to/S2*.SAFE
//!   cargo run --release --features full-io --example process_sentinel2_ac -- --file /path/to/S2*.SAFE --output /tmp/out --res 20

use acolite_rs::core::BandData;
use acolite_rs::ac::{gas_lut, dsf};
use acolite_rs::loader::sentinel2::{load_sentinel2_scene, S2_AC_BANDS};
use acolite_rs::writer::write_auto;
use ndarray::Array2;
use std::path::Path;

fn get_arg(args: &[String], flag: &str) -> Option<String> {
    args.iter().position(|a| a == flag).and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    env_logger::init();
    let args: Vec<String> = std::env::args().collect();
    println!("ACOLITE-RS: Sentinel-2 Processing (LUT-based DSF)\n");

    let file = get_arg(&args, "--file");
    let output_dir = get_arg(&args, "--output").unwrap_or_else(|| "/tmp/acolite_s2".into());
    let target_res: u32 = get_arg(&args, "--res").and_then(|s| s.parse().ok()).unwrap_or(20);

    if let Some(safe_dir) = file {
        process_real(Path::new(&safe_dir), &output_dir, target_res);
    } else {
        println!("Usage: process_sentinel2_ac --file /path/to/S2*.SAFE [--output /tmp/out] [--res 20]");
    }
}

fn process_real(safe_dir: &Path, output_dir: &str, target_res: u32) {
    println!("→ Loading scene: {:?} (target {}m)", safe_dir, target_res);
    let start = std::time::Instant::now();

    let scene = load_sentinel2_scene(safe_dir, target_res).expect("Failed to load S2 scene");
    let load_time = start.elapsed();
    println!("  ✓ Loaded {} bands in {:.2?}", scene.bands.len(), load_time);

    let metadata = &scene.metadata;
    let sensor_lut = &scene.sensor_lut;

    if scene.bands.is_empty() || scene.bands[0].data.is_empty() {
        eprintln!("No data");
        return;
    }
    let shape = scene.bands[0].data.dim();
    println!("  Working with {} bands, {}×{} pixels", scene.bands.len(), shape.0, shape.1);

    // Find data directory
    let data_dir = find_data_dir(safe_dir);
    println!("  Data dir: {:?}", data_dir);

    // ── Step 1: Compute gas transmittance ──
    let sza = metadata.sun_zenith;
    let vza = metadata.view_zenith.unwrap_or(0.0);
    let uoz = 0.3;
    let uwv = 1.5;

    let gas_t = gas_lut::compute_gas_transmittance(&data_dir, sensor_lut, sza, vza, uoz, uwv)
        .expect("Gas transmittance computation failed");
    println!("  Gas transmittance computed for {} bands", gas_t.tt_gas.len());

    // ── Step 2: Load aerosol LUTs ──
    #[cfg(feature = "full-io")]
    let luts = {
        let pressures = vec![500.0, 750.0, 1013.0, 1100.0];
        acolite_rs::ac::aerlut::load_acolite_luts(&data_dir, sensor_lut, &pressures)
            .expect("Failed to load aerosol LUTs")
    };
    #[cfg(feature = "full-io")]
    println!("  Loaded {} aerosol models", luts.len());

    // ── Step 3: Convert DN to TOA reflectance ──
    // S2 L1C: DN / 10000 = TOA reflectance (already divided by cos(sza) in L1C)
    let cos_sza = sza.to_radians().cos();
    println!("  Calibration: DN/10000, cos(sza)={:.4}", cos_sza);

    // Filter to AC bands only
    let ac_band_set: std::collections::HashSet<&str> = S2_AC_BANDS.iter().copied().collect();

    let toa_bands: Vec<(String, String, f64, Array2<f64>)> = scene.bands.iter()
        .filter_map(|b| {
            if !ac_band_set.contains(b.name.as_str()) { return None; }
            scene.lut_band_map.iter().find(|(rn, _)| *rn == b.name)
                .map(|(_, lut_bn)| {
                    let toa = b.data.mapv(|v| {
                        let dn = v as f64;
                        if dn == 0.0 { return f64::NAN; }
                        dn / 10000.0
                    });
                    (b.name.clone(), lut_bn.to_string(), b.wavelength, toa)
                })
        })
        .collect();

    println!("  {} bands for AC", toa_bands.len());

    // ── Step 4: Tiled DSF atmospheric correction ──
    let ac_start = std::time::Instant::now();

    let band_names_lut: Vec<String> = toa_bands.iter().map(|(_, ln, _, _)| ln.clone()).collect();
    let wavelengths: Vec<f64> = toa_bands.iter().map(|(_, _, wl, _)| *wl).collect();
    let tt_gas_vec: Vec<f64> = toa_bands.iter().map(|(_, ln, _, _)| {
        *gas_t.tt_gas.get(ln).unwrap_or(&1.0)
    }).collect();

    // Gas-correct TOA before dark spectrum extraction
    let toa_gc_arrays: Vec<Array2<f64>> = toa_bands.iter().enumerate().map(|(i, (_, _, _, toa))| {
        let tt = tt_gas_vec[i];
        toa.mapv(|v| if v.is_finite() && v > 0.0 { v / tt } else { f64::NAN })
    }).collect();

    let mut dsf_config = dsf::DsfConfig::default();
    dsf_config.dark_method = dsf::DarkSpectrumMethod::Intercept(200);

    let pressure = 1013.0;
    let vaa = metadata.view_azimuth.unwrap_or(100.0);
    let saa = metadata.sun_azimuth;
    let mut raa = (saa - vaa).abs();
    if raa > 180.0 { raa = (360.0 - raa).abs(); }
    let thv = vza.max(0.001);
    let ths = sza;
    println!("  Geometry: SZA={:.2}, VZA={:.2}, RAA={:.2}, pressure={:.0}", ths, thv, raa, pressure);

    #[cfg(feature = "full-io")]
    let tiled_result = dsf::optimize_aot_tiled(
        &luts, &toa_gc_arrays, &band_names_lut, &wavelengths, &tt_gas_vec,
        pressure, raa, thv, ths, &dsf_config, (200, 200),
    );
    #[cfg(feature = "full-io")]
    {
        let vals: Vec<f64> = tiled_result.aot_grid.iter().flatten().copied().filter(|v| v.is_finite()).collect();
        let mean_aot = if vals.is_empty() { 0.0 } else { vals.iter().sum::<f64>() / vals.len() as f64 };
        println!("  Tiled DSF: model={}, tiles={}×{}, mean AOT={:.4}",
            tiled_result.model_name, tiled_result.ni, tiled_result.nj, mean_aot);
    }

    // ── Step 5: Apply atmospheric correction ──
    #[cfg(feature = "full-io")]
    let selected_lut = &luts[tiled_result.model_idx];

    let mut result_bands: Vec<BandData<f64>> = Vec::new();
    for (i, (rust_name, lut_bn, wl, toa)) in toa_bands.iter().enumerate() {
        let tt = tt_gas_vec[i];

        #[cfg(feature = "full-io")]
        let corrected = dsf::dsf_correct_band_tiled(
            toa, selected_lut, lut_bn, tt, &tiled_result, pressure, raa, thv, ths,
        );
        #[cfg(not(feature = "full-io"))]
        let corrected = toa.clone();

        let orig = scene.bands.iter().find(|b| b.name == *rust_name).unwrap();
        result_bands.push(BandData::new(
            corrected, *wl, orig.bandwidth, rust_name.clone(),
            orig.projection.clone(), orig.geotransform.clone(),
        ));
    }

    let ac_time = ac_start.elapsed();
    println!("  ✓ AC in {:.2?}", ac_time);

    // ── Step 6: Write output ──
    std::fs::create_dir_all(output_dir).ok();
    let stem = safe_dir.file_name().unwrap().to_string_lossy().replace(".SAFE", "");
    let out_path = format!("{}/{}_corrected", output_dir, stem);

    let write_start = std::time::Instant::now();
    write_auto(&out_path, &result_bands, metadata).expect("Write failed");
    let write_time = write_start.elapsed();

    let total = start.elapsed();
    println!("\n→ Results:");
    for band in &result_bands {
        let valid: Vec<f64> = band.data.iter().copied().filter(|v| v.is_finite()).collect();
        let mean = if valid.is_empty() { 0.0 } else { valid.iter().sum::<f64>() / valid.len() as f64 };
        println!("  {}: ρs_mean={:.6}, λ={:.0}nm", band.name, mean, band.wavelength);
    }
    println!("\n  Load: {:.2?}, AC: {:.2?}, Write: {:.2?}, Total: {:.2?}", load_time, ac_time, write_time, total);
    println!("  Output: {}", out_path);
}

fn find_data_dir(scene_dir: &Path) -> std::path::PathBuf {
    let candidates = [
        std::path::PathBuf::from("data"),
        scene_dir.join("../../data").canonicalize().unwrap_or_default(),
    ];
    for c in &candidates {
        if c.join("RSR").exists() && c.join("LUT").exists() {
            return c.clone();
        }
    }
    let mut p = std::env::current_dir().unwrap_or_default();
    for _ in 0..5 {
        let d = p.join("data");
        if d.join("RSR").exists() && d.join("LUT").exists() {
            return d;
        }
        if !p.pop() { break; }
    }
    std::path::PathBuf::from("data")
}
