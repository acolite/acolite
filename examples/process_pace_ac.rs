//! Process PACE OCI with full LUT-based DSF atmospheric correction
//!
//! Usage:
//!   cargo run --release --features full-io --example process_pace_ac -- --file /path/to/PACE_OCI.*.L1B.*.nc
//!   cargo run --release --features full-io --example process_pace_ac -- --file /path/to/PACE_OCI.*.L1B.*.nc --limit 35.0,-75.0,36.0,-74.0
//!   cargo run --release --features full-io --example process_pace_ac -- --file /path/to/PACE_OCI.*.L1B.*.nc --model auto --aot-mode tiled

use acolite_rs::ac::{dsf, gas_lut};
use acolite_rs::core::BandData;
use acolite_rs::loader::pace::load_pace_l1b;
use acolite_rs::writer::write_auto;
use ndarray::Array2;
use rayon::prelude::*;
use std::path::Path;

fn get_arg(args: &[String], flag: &str) -> Option<String> {
    args.iter().position(|a| a == flag).and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    env_logger::init();
    let args: Vec<String> = std::env::args().collect();
    println!("ACOLITE-RS: PACE OCI Processing (LUT-based DSF)\n");

    let file = get_arg(&args, "--file");
    let output_dir = get_arg(&args, "--output").unwrap_or_else(|| "/tmp/acolite_pace_ac".into());
    let model = get_arg(&args, "--model").unwrap_or_else(|| "auto".into());
    let aot_mode = get_arg(&args, "--aot-mode").unwrap_or_else(|| "fixed".into());
    let limit: Option<[f64; 4]> = get_arg(&args, "--limit").map(|s| {
        let v: Vec<f64> = s.split(',').map(|x| x.parse().expect("invalid --limit value")).collect();
        assert!(v.len() == 4, "--limit needs 4 values: south,west,north,east");
        [v[0], v[1], v[2], v[3]]
    });

    if let Some(f) = file {
        process(Path::new(&f), &output_dir, &model, &aot_mode, limit.as_ref());
    } else {
        println!("Usage: process_pace_ac --file /path/to/PACE_OCI.*.L1B.*.nc [--output /tmp/out] [--model auto|MOD1|MOD2] [--aot-mode fixed|tiled] [--limit s,w,n,e]");
    }
}

fn process(path: &Path, output_dir: &str, model: &str, aot_mode: &str, limit: Option<&[f64; 4]>) {
    let total_start = std::time::Instant::now();

    // ── Load scene ──
    println!("→ Loading PACE OCI scene: {:?}", path);
    let start = std::time::Instant::now();
    let scene = load_pace_l1b(path, limit).expect("Failed to load PACE L1B");
    let load_time = start.elapsed();

    let (nrows, ncols) = (scene.lat.nrows(), scene.lat.ncols());
    println!("  Sensor: {}", scene.metadata.sensor);
    println!("  Bands: {}", scene.bands.len());
    println!("  Size: {}×{}", nrows, ncols);
    println!("  Loaded in {:.2?}", load_time);

    if scene.bands.is_empty() {
        eprintln!("No bands loaded");
        return;
    }

    // ── Find data directory ──
    let data_dir = find_data_dir(path);
    println!("  Data dir: {:?}", data_dir);

    // ── Geometry (scene-average for DSF) ──
    let sza = scene.metadata.sun_zenith;
    let vza = scene.metadata.view_zenith.unwrap_or(0.0);
    let saa = scene.metadata.sun_azimuth;
    let vaa = scene.metadata.view_azimuth.unwrap_or(100.0);
    let mut raa = (saa - vaa).abs();
    if raa > 180.0 { raa = (360.0 - raa).abs(); }
    let thv = vza.max(0.001);
    let ths = sza;
    let pressure = 1013.0;
    let uoz = 0.3;
    let uwv = 1.5;
    println!("  Geometry: SZA={:.2}, VZA={:.2}, RAA={:.2}", ths, thv, raa);

    // ── Collect band info ──
    let wavelengths: Vec<f64> = scene.bands.iter().map(|b| b.wavelength).collect();
    let bandwidths: Vec<f64> = scene.bands.iter().map(|b| b.bandwidth).collect();

    // ── Gas transmittance (per-band, LUT-based matching Python ACOLITE) ──
    println!("\n→ Computing gas transmittance (LUT-based)...");

    let hyper_tg = gas_lut::compute_gas_transmittance_hyper(
        &data_dir, &wavelengths, &bandwidths, ths, thv, pressure, uoz, uwv,
    ).expect("Failed to compute gas transmittance");

    let tt_gas_vec = hyper_tg.tt_gas;
    println!("  Gas transmittance computed for {} bands (LUT-based)", tt_gas_vec.len());

    // ── Load generic aerosol LUTs ──
    println!("\n→ Loading generic aerosol LUTs...");
    let start = std::time::Instant::now();
    let pressures = vec![500.0, 750.0, 1013.0, 1100.0];
    let luts = if model == "auto" {
        acolite_rs::ac::aerlut::load_generic_luts(&data_dir, &pressures)
            .expect("Failed to load generic LUTs")
    } else {
        let full_name = format!("ACOLITE-LUT-202110-{}", model);
        let lut = acolite_rs::ac::aerlut::load_generic_lut(
            &data_dir.join("LUT"), &full_name, &pressures,
        ).expect("Failed to load generic LUT");
        vec![lut]
    };
    println!("  Loaded {} model(s) in {:.2?}", luts.len(), start.elapsed());

    // ── Convert f32 TOA to f64 and gas-correct (parallel) ──
    println!("\n→ Gas-correcting TOA reflectance...");
    let toa_gc_bands: Vec<Array2<f64>> = scene.bands.par_iter().enumerate().map(|(i, band)| {
        let tt = tt_gas_vec[i];
        band.data.mapv(|v| {
            let v = v as f64;
            if v.is_finite() && v > 0.0 { v / tt } else { f64::NAN }
        })
    }).collect();

    // ── DSF atmospheric correction ──
    println!("\n→ DSF atmospheric correction ({} mode)...", aot_mode);
    let ac_start = std::time::Instant::now();

    let mut dsf_config = dsf::DsfConfig::default();
    dsf_config.dark_method = dsf::DarkSpectrumMethod::Intercept(1000);
    // Match Python ACOLITE defaults: use full spectral range, skip bands with tgas < 0.85
    dsf_config.wave_range = (400.0, 2500.0);
    dsf_config.min_tgas_aot = 0.85;
    if model != "auto" {
        dsf_config.fixed_model = Some(model.to_string());
    }
    if aot_mode == "fixed" {
        dsf_config.mode = dsf::DsfMode::Fixed;
    }

    let (dsf_result_fixed, tiled_result) = match dsf_config.mode {
        dsf::DsfMode::Fixed => {
            let result = dsf::optimize_aot_fixed_generic(
                &luts, &toa_gc_bands, &wavelengths, &bandwidths, &tt_gas_vec,
                pressure, raa, thv, ths, &dsf_config,
            );
            println!("  Fixed DSF: model={}, AOT={:.4}, RMSD={:.6}",
                result.model_name, result.aot, result.rmsd);
            (Some(result), None)
        }
        dsf::DsfMode::Tiled(tr, tc) => {
            let tr = dsf::optimize_aot_tiled_generic(
                &luts, &toa_gc_bands, &wavelengths, &bandwidths, &tt_gas_vec,
                pressure, raa, thv, ths, &dsf_config, (tr, tc),
            );
            let vals: Vec<f64> = tr.aot_grid.iter().flatten().copied().filter(|v| v.is_finite()).collect();
            let mean_aot = if vals.is_empty() { 0.0 } else { vals.iter().sum::<f64>() / vals.len() as f64 };
            println!("  Tiled DSF: model={}, tiles={}×{}, mean AOT={:.4}",
                tr.model_name, tr.ni, tr.nj, mean_aot);
            (None, Some(tr))
        }
    };

    // ── Apply correction to all bands (parallel with rayon) ──
    let selected_lut = if let Some(ref r) = dsf_result_fixed {
        &luts[r.model_idx]
    } else {
        &luts[tiled_result.as_ref().unwrap().model_idx]
    };

    let result_bands: Vec<BandData<f64>> = scene.bands.par_iter().enumerate().map(|(i, band)| {
        let tt = tt_gas_vec[i];
        let toa_f64 = band.data.mapv(|v| v as f64);

        let corrected = if let Some(ref r) = dsf_result_fixed {
            dsf::dsf_correct_band_generic(
                &toa_f64, selected_lut, band.wavelength, band.bandwidth,
                tt, r.aot, pressure, raa, thv, ths,
            )
        } else {
            dsf::dsf_correct_band_tiled_generic(
                &toa_f64, selected_lut, band.wavelength, band.bandwidth,
                tt, tiled_result.as_ref().unwrap(), pressure, raa, thv, ths,
            )
        };

        BandData::new(
            corrected, band.wavelength, band.bandwidth, band.name.clone(),
            band.projection.clone(), band.geotransform.clone(),
        )
    }).collect();

    let ac_time = ac_start.elapsed();
    let mpx = (nrows * ncols * result_bands.len()) as f64 / ac_time.as_secs_f64() / 1e6;
    println!("  ✓ {} bands in {:.2?} ({:.1} Mpx/s)", result_bands.len(), ac_time, mpx);

    // ── Write output ──
    std::fs::create_dir_all(output_dir).ok();
    let stem = path.file_stem().unwrap_or_default().to_string_lossy();
    let out_path = format!("{}/{}_L2R", output_dir, stem);

    println!("\n→ Writing output...");
    let write_start = std::time::Instant::now();
    write_auto(&out_path, &result_bands, &scene.metadata).expect("Write failed");
    let write_time = write_start.elapsed();

    // ── Summary ──
    let total = total_start.elapsed();
    println!("\n→ Sample reflectances:");
    let step = result_bands.len().max(1) / 10;
    for band in result_bands.iter().step_by(step.max(1)) {
        let valid: Vec<f64> = band.data.iter().copied().filter(|v| v.is_finite()).collect();
        let mean = if valid.is_empty() { 0.0 } else { valid.iter().sum::<f64>() / valid.len() as f64 };
        println!("  {:.0} nm: ρs = {:.6}", band.wavelength, mean);
    }

    println!("\n  Load: {:.2?}, AC: {:.2?}, Write: {:.2?}, Total: {:.2?}", load_time, ac_time, write_time, total);
    println!("  Output: {}", out_path);
    println!("\n✓ PACE OCI DSF processing complete!");
}

fn find_data_dir(scene_path: &Path) -> std::path::PathBuf {
    let mut p = std::env::current_dir().unwrap_or_default();
    for _ in 0..5 {
        let d = p.join("data");
        if d.join("RSR").exists() && d.join("LUT").exists() {
            return d;
        }
        if !p.pop() { break; }
    }
    // Try relative to scene
    if let Some(parent) = scene_path.parent() {
        let d = parent.join("../../data").canonicalize().unwrap_or_default();
        if d.join("RSR").exists() { return d; }
    }
    std::path::PathBuf::from("data")
}
