//! Process Landsat scene with full LUT-based DSF atmospheric correction
//!
//! Usage:
//!   cargo run --release --features full-io --example process_landsat -- --file /path/to/LC08_*
//!   cargo run --release --features full-io --example process_landsat -- --file /path/to/LC08_* --output /tmp/out --limit -35.0,138.0,-34.5,138.7
//!   cargo run --release --features full-io --example process_landsat -- --file /path/to/LC08_* --model MOD1 --aot-mode fixed

use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::ac::gas_lut;
use acolite_rs::ac::dsf;
use acolite_rs::writer::write_auto;
use ndarray::Array2;
use std::path::Path;

fn get_arg(args: &[String], flag: &str) -> Option<String> {
    args.iter().position(|a| a == flag).and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    env_logger::init();
    let args: Vec<String> = std::env::args().collect();
    println!("ACOLITE-RS: Landsat Processing (LUT-based DSF)\n");

    let file = get_arg(&args, "--file");
    let output_dir = get_arg(&args, "--output").unwrap_or_else(|| "/tmp/acolite_landsat".into());
    let model = get_arg(&args, "--model").unwrap_or_else(|| "auto".into());
    let aot_mode = get_arg(&args, "--aot-mode").unwrap_or_else(|| "tiled".into());
    let limit: Option<[f64; 4]> = get_arg(&args, "--limit").map(|s| {
        let v: Vec<f64> = s.split(',').map(|x| x.parse().expect("invalid --limit value")).collect();
        assert!(v.len() == 4, "--limit needs 4 values: south,west,north,east");
        [v[0], v[1], v[2], v[3]]
    });

    if let Some(scene_dir) = file {
        process_real(Path::new(&scene_dir), &output_dir, &model, &aot_mode, limit.as_ref());
    } else {
        println!("Usage: process_landsat --file /path/to/LC08_* [--output /tmp/out] [--limit s,w,n,e] [--model auto|MOD1|MOD2] [--aot-mode tiled|fixed]");
    }
}

fn process_real(scene_dir: &Path, output_dir: &str, model: &str, aot_mode: &str, limit: Option<&[f64; 4]>) {
    use acolite_rs::{load_landsat_scene_limit, sensors::parse_mtl};

    println!("→ Loading scene: {:?}", scene_dir);
    if let Some(lim) = limit {
        println!("  Limit: [{:.3}, {:.3}, {:.3}, {:.3}]", lim[0], lim[1], lim[2], lim[3]);
    }
    let start = std::time::Instant::now();

    let (bands, _metadata) = load_landsat_scene_limit(scene_dir, limit).expect("Failed to load scene");
    let load_time = start.elapsed();
    println!("  ✓ Loaded {} bands in {:.2?}", bands.len(), load_time);

    // Parse MTL for proper metadata
    let mtl_path = find_mtl(scene_dir);
    let metadata = if let Some(ref p) = mtl_path {
        match parse_mtl(p) {
            Ok(m) => { println!("  Sun zenith={:.1}°", m.sun_zenith); m },
            Err(e) => { eprintln!("  MTL parse warning: {}", e); _metadata },
        }
    } else {
        _metadata
    };

    if bands.is_empty() || bands[0].data.is_empty() {
        eprintln!("No data");
        return;
    }
    let shape = bands[0].data.dim();
    println!("  Working with {} bands, {}×{} pixels", bands.len(), shape.0, shape.1);

    // Determine sensor for LUT loading
    let scene_name = scene_dir.file_name().unwrap().to_string_lossy();
    let (sensor_lut, lut_band_map) = if scene_name.starts_with("LC09") {
        ("L9_OLI", vec![("B1","1"), ("B2","2"), ("B3","3"), ("B4","4"), ("B5","5"), ("B6","6"), ("B7","7")])
    } else {
        ("L8_OLI", vec![("B1","1"), ("B2","2"), ("B3","3"), ("B4","4"), ("B5","5"), ("B6","6"), ("B7","7")])
    };

    // Find data directory (relative to repo root)
    let data_dir = find_data_dir(scene_dir);
    println!("  Data dir: {:?}", data_dir);

    // ── Step 1: Compute gas transmittance ──
    let sza = metadata.sun_zenith;
    let vza = metadata.view_zenith.unwrap_or(0.0);
    let uoz = 0.3;  // default ozone (cm-atm)
    let uwv = 1.5;  // default water vapour (g/cm²)

    let gas_t = gas_lut::compute_gas_transmittance(&data_dir, sensor_lut, sza, vza, uoz, uwv)
        .expect("Gas transmittance computation failed");
    println!("  Gas transmittance computed for {} bands", gas_t.tt_gas.len());

    // ── Step 2: Load aerosol LUTs ──
    #[cfg(feature = "full-io")]
    let luts = {
        let pressures = vec![500.0, 750.0, 1013.0, 1100.0];
        if model == "auto" {
            acolite_rs::ac::aerlut::load_acolite_luts(&data_dir, sensor_lut, &pressures)
                .expect("Failed to load aerosol LUTs")
        } else {
            acolite_rs::ac::aerlut::load_single_lut(&data_dir, sensor_lut, model, &pressures)
                .expect("Failed to load aerosol LUT")
        }
    };
    #[cfg(feature = "full-io")]
    println!("  Loaded {} aerosol model(s)", luts.len());

    // ── Step 3: Convert DN to TOA reflectance ──
    // Landsat Collection 2: rhot = DN * MULT + ADD, then / cos(sza)
    let (refl_mult, refl_add) = parse_reflectance_coeffs(scene_dir);
    let cos_sza = sza.to_radians().cos();
    println!("  Calibration: mult={}, add={}, cos(sza)={:.4}", refl_mult, refl_add, cos_sza);

    let toa_bands: Vec<(String, String, f64, Array2<f64>)> = bands.iter()
        .filter_map(|b| {
            let rust_name = &b.name;
            lut_band_map.iter().find(|(rn, _)| *rn == rust_name)
                .map(|(_, lut_bn)| {
                    let toa = b.data.mapv(|v| {
                        let dn = v as f64;
                        if dn == 0.0 { return f64::NAN; } // fill value
                        (dn * refl_mult + refl_add) / cos_sza
                    });
                    (rust_name.clone(), lut_bn.to_string(), b.wavelength, toa)
                })
        })
        .collect();

    // ── Step 4: Extract dark spectrum and estimate AOT ──
    let ac_start = std::time::Instant::now();

    // Prepare arrays for DSF
    let band_names_lut: Vec<String> = toa_bands.iter().map(|(_, ln, _, _)| ln.clone()).collect();
    let wavelengths: Vec<f64> = toa_bands.iter().map(|(_, _, wl, _)| *wl).collect();
    let tt_gas_vec: Vec<f64> = toa_bands.iter().map(|(_, ln, _, _)| {
        *gas_t.tt_gas.get(ln).unwrap_or(&1.0)
    }).collect();

    // Extract dark spectrum using intercept method (matching Python ACOLITE default)
    // Python gas-corrects TOA BEFORE extracting dark spectrum
    let toa_gc_arrays: Vec<Array2<f64>> = toa_bands.iter().enumerate().map(|(i, (_, _, _, toa))| {
        let tt = tt_gas_vec[i];
        toa.mapv(|v| if v.is_finite() && v > 0.0 { v / tt } else { f64::NAN })
    }).collect();
    let mut dsf_config = dsf::DsfConfig::default();
    // Use intercept method matching Python's dsf_spectrum_option=intercept
    dsf_config.dark_method = dsf::DarkSpectrumMethod::Intercept(200);
    if model != "auto" {
        dsf_config.fixed_model = Some(model.to_string());
    }
    if aot_mode == "fixed" {
        dsf_config.mode = dsf::DsfMode::Fixed;
    }

    // Compute relative azimuth from MTL corner coordinates
    let pressure = 1013.0;
    let azi = compute_raa_from_mtl(&metadata);
    let thv = vza.max(0.001);
    let ths = sza;
    println!("  Geometry: SZA={:.2}, VZA={:.2}, RAA={:.2}, pressure={:.0}", ths, thv, azi, pressure);

    // AOT estimation based on configured mode
    #[cfg(feature = "full-io")]
    let (dsf_result_fixed, tiled_result) = match dsf_config.mode {
        dsf::DsfMode::Fixed => {
            let result = dsf::optimize_aot_fixed(
                &luts, &toa_gc_arrays, &band_names_lut, &wavelengths, &tt_gas_vec,
                pressure, azi, thv, ths, &dsf_config,
            );
            println!("  Fixed DSF: model={}, AOT={:.4}, RMSD={:.6}",
                result.model_name, result.aot, result.rmsd);
            (Some(result), None)
        }
        dsf::DsfMode::Tiled(tr, tc) => {
            let tr = dsf::optimize_aot_tiled(
                &luts, &toa_gc_arrays, &band_names_lut, &wavelengths, &tt_gas_vec,
                pressure, azi, thv, ths, &dsf_config, (tr, tc),
            );
            let vals: Vec<f64> = tr.aot_grid.iter().flatten().copied().filter(|v| v.is_finite()).collect();
            let mean_aot = if vals.is_empty() { 0.0 } else { vals.iter().sum::<f64>() / vals.len() as f64 };
            println!("  Tiled DSF: model={}, tiles={}×{}, mean AOT={:.4}",
                tr.model_name, tr.ni, tr.nj, mean_aot);
            (None, Some(tr))
        }
    };

    // ── Step 5: Apply atmospheric correction ──
    #[cfg(feature = "full-io")]
    let selected_lut = if let Some(ref r) = dsf_result_fixed {
        &luts[r.model_idx]
    } else {
        &luts[tiled_result.as_ref().unwrap().model_idx]
    };

    let mut result_bands: Vec<BandData<f64>> = Vec::new();
    for (i, (rust_name, lut_bn, wl, toa)) in toa_bands.iter().enumerate() {
        let tt = tt_gas_vec[i];

        #[cfg(feature = "full-io")]
        let corrected = if let Some(ref r) = dsf_result_fixed {
            dsf::dsf_correct_band(toa, selected_lut, lut_bn, tt, r.aot, pressure, azi, thv, ths)
        } else {
            dsf::dsf_correct_band_tiled(
                toa, selected_lut, lut_bn, tt, tiled_result.as_ref().unwrap(), pressure, azi, thv, ths,
            )
        };
        #[cfg(not(feature = "full-io"))]
        let corrected = toa.clone();

        let orig = &bands[i];
        result_bands.push(BandData::new(
            corrected,
            *wl,
            orig.bandwidth,
            rust_name.clone(),
            orig.projection.clone(),
            orig.geotransform.clone(),
        ));
    }

    let ac_time = ac_start.elapsed();
    println!("  ✓ AC in {:.2?}", ac_time);

    // ── Step 6: Write output ──
    std::fs::create_dir_all(output_dir).ok();
    let stem = scene_dir.file_name().unwrap().to_string_lossy();
    let out_path = format!("{}/{}_corrected", output_dir, stem);

    let write_start = std::time::Instant::now();
    write_auto(&out_path, &result_bands, &metadata).expect("Write failed");
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

fn find_mtl(scene_dir: &Path) -> Option<String> {
    std::fs::read_dir(scene_dir).ok()?.find_map(|e| {
        let p = e.ok()?.path();
        let name = p.file_name()?.to_string_lossy().to_string();
        if name.ends_with("_MTL.txt") { Some(p.to_string_lossy().to_string()) } else { None }
    })
}

/// Compute relative azimuth angle from MTL corner coordinates.
/// Uses geodesic forward azimuth (matching Python ACOLITE's azimuth_two_points).
/// RAA = |SAA - VAA|, clamped to [0, 180].
fn compute_raa_from_mtl(metadata: &acolite_rs::core::Metadata) -> f64 {
    let get = |key: &str| -> Option<f64> {
        metadata.attributes.get(key).and_then(|v| v.parse().ok())
    };
    let saa = metadata.sun_azimuth;

    // Compute VAA from corner coordinates (nadir line direction)
    if let (Some(ul_lat), Some(ul_lon), Some(ur_lat), Some(ur_lon),
            Some(ll_lat), Some(ll_lon), Some(lr_lat), Some(lr_lon)) = (
        get("CORNER_UL_LAT_PRODUCT"), get("CORNER_UL_LON_PRODUCT"),
        get("CORNER_UR_LAT_PRODUCT"), get("CORNER_UR_LON_PRODUCT"),
        get("CORNER_LL_LAT_PRODUCT"), get("CORNER_LL_LON_PRODUCT"),
        get("CORNER_LR_LAT_PRODUCT"), get("CORNER_LR_LON_PRODUCT"),
    ) {
        // Nadir line: midpoint of top-right edge to midpoint of left-bottom edge
        // (matching Python's image_corners nadir computation)
        let top_lat = (ul_lat + ur_lat) / 2.0;
        let top_lon = (ul_lon + ur_lon) / 2.0;
        let bot_lat = (ll_lat + lr_lat) / 2.0;
        let bot_lon = (ll_lon + lr_lon) / 2.0;

        // Geodesic forward azimuth (matching Python's azimuth_two_points)
        let lat1r = top_lat.to_radians();
        let lat2r = bot_lat.to_radians();
        let dlonr = (bot_lon - top_lon).to_radians();
        let x = dlonr.sin() * lat2r.cos();
        let y = lat1r.cos() * lat2r.sin() - lat1r.sin() * lat2r.cos() * dlonr.cos();
        let vaa = (x.atan2(y).to_degrees() + 360.0) % 360.0;

        let mut raa = (saa - vaa).abs();
        if raa > 180.0 { raa = (360.0 - raa).abs(); }
        println!("  VAA={:.2} (geodesic), SAA={:.2}, RAA={:.2}", vaa, saa, raa);
        return raa;
    }
    println!("  Warning: could not compute RAA from corners, using 150°");
    150.0
}

/// Parse reflectance calibration coefficients from MTL file
fn parse_reflectance_coeffs(scene_dir: &Path) -> (f64, f64) {
    if let Some(mtl_path) = find_mtl(scene_dir) {
        if let Ok(content) = std::fs::read_to_string(&mtl_path) {
            let mut mult = 2.0e-5;
            let mut add = -0.1;
            for line in content.lines() {
                let line = line.trim();
                if line.starts_with("REFLECTANCE_MULT_BAND_1") {
                    if let Some(v) = line.split('=').nth(1) {
                        mult = v.trim().parse().unwrap_or(mult);
                    }
                }
                if line.starts_with("REFLECTANCE_ADD_BAND_1") {
                    if let Some(v) = line.split('=').nth(1) {
                        add = v.trim().parse().unwrap_or(add);
                    }
                }
            }
            return (mult, add);
        }
    }
    (2.0e-5, -0.1) // Collection 2 defaults
}

/// Find the ACOLITE data directory by searching upward from the scene dir
fn find_data_dir(scene_dir: &Path) -> std::path::PathBuf {
    // Try common locations
    let candidates = [
        // Relative to CWD
        std::path::PathBuf::from("data"),
        // Relative to scene dir (up several levels)
        scene_dir.join("../../data").canonicalize().unwrap_or_default(),
    ];
    for c in &candidates {
        if c.join("RSR").exists() && c.join("LUT").exists() {
            return c.clone();
        }
    }
    // Fallback: try the repo root
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
