//! Sentinel-3 OLCI L1→L2R atmospheric correction pipeline
//!
//! Implements the full processing chain matching Python ACOLITE:
//! 1. Load .SEN3 bundle (radiance, TPGs, instrument data)
//! 2. Apply smile correction (optional)
//! 3. Convert radiance → TOA reflectance
//! 4. Gas transmittance correction (O3, H2O, O2, NO2)
//! 5. DSF atmospheric correction (LUT-based or simplified fallback)
//! 6. Output surface reflectance ρs = Rrs·π

use crate::ac::dsf::{DsfConfig, DsfMode};
use crate::ac::gas::gas_correction;
use crate::ac::rayleigh::rayleigh_optical_thickness;
use crate::core::{BandData, Metadata};
use crate::loader::sentinel3::{OlciScene, OLCI_BANDS};
use crate::pipeline::{Pipeline, ProcessingConfig};
use crate::Result;
use ndarray::Array2;
use std::collections::HashMap;

/// Configuration for S3 OLCI processing
#[derive(Debug, Clone)]
pub struct OlciProcessingConfig {
    pub smile_correction: bool,
    pub smile_correction_tgas: bool,
    pub use_supplied_ancillary: bool,
    pub ancillary_data: bool,
    pub dsf: DsfConfig,
    pub output_lt: bool,
    /// Path to data directory containing RSR/ and LUT/ for full DSF
    pub data_dir: Option<std::path::PathBuf>,
}

impl Default for OlciProcessingConfig {
    fn default() -> Self {
        Self {
            smile_correction: true,
            smile_correction_tgas: true,
            use_supplied_ancillary: true,
            ancillary_data: true, // match Python default — use defaults unless ancillary_data=false
            dsf: DsfConfig::default(),
            output_lt: false,
            data_dir: None,
        }
    }
}

/// Result of OLCI L2R processing
pub struct OlciL2rResult {
    pub rhos: HashMap<String, Array2<f64>>,
    pub metadata: Metadata,
    pub aot: f64,
    pub model_name: String,
}

/// Process an OLCI scene through the atmospheric correction pipeline.
pub fn process_olci_l2r(
    scene: &mut OlciScene,
    config: &OlciProcessingConfig,
    pipeline_config: &ProcessingConfig,
) -> Result<OlciL2rResult> {
    let (sza, vza, raa) = scene.tpg.mean_geometry();
    let mu = (sza * std::f64::consts::PI / 180.0).cos();

    // Ancillary: use TPG values only when use_supplied_ancillary=true AND ancillary_data=false
    // (matching Python ACOLITE behavior). Otherwise use pipeline defaults.
    let use_tpg_anc = config.use_supplied_ancillary && !config.ancillary_data;
    let uoz = if use_tpg_anc {
        scene.tpg.total_ozone.as_ref()
            .map(|o| nanmean(o) / 0.02141419) // kg/m² → cm·atm
            .unwrap_or(pipeline_config.ozone)
    } else {
        pipeline_config.ozone
    };
    let uwv = if use_tpg_anc {
        scene.tpg.total_columnar_water_vapour.as_ref()
            .map(|w| nanmean(w) / 10.0) // kg/m² → g/cm²
            .unwrap_or(pipeline_config.water_vapor)
    } else {
        pipeline_config.water_vapor
    };
    let pressure = if use_tpg_anc {
        scene.tpg.sea_level_pressure.as_ref()
            .map(|p| nanmean(p))
            .unwrap_or(pipeline_config.pressure)
    } else {
        pipeline_config.pressure
    };

    // Step 1: Smile correction
    if config.smile_correction {
        let tt_gas = if config.smile_correction_tgas {
            let airmass = 1.0 / (sza.to_radians().cos()) + 1.0 / (vza.to_radians().cos());
            let mut tg = HashMap::new();
            for (name, wl, _, _) in &OLCI_BANDS {
                let t_o3 = crate::ac::gas::ozone_transmittance(*wl, uoz, airmass);
                let t_wv = crate::ac::gas::water_vapor_transmittance(*wl, uwv, airmass);
                tg.insert(name.to_string(), t_o3 * t_wv);
            }
            Some(tg)
        } else {
            None
        };
        scene.apply_smile_correction(tt_gas.as_ref());
    }

    // Step 2: Radiance → TOA reflectance
    let rhot = scene.to_toa_reflectance(config.smile_correction);

    // Step 3: Try full LUT-based DSF if data_dir is available
    #[cfg(feature = "full-io")]
    if let Some(ref data_dir) = config.data_dir {
        return process_olci_dsf(scene, &rhot, data_dir, config, sza, vza, raa, pressure, uoz, uwv);
    }

    // Fallback: simplified Rayleigh-only correction
    process_olci_simplified(scene, &rhot, sza, vza, raa, pressure, uoz, uwv, mu)
}

/// Full LUT-based DSF atmospheric correction for OLCI
#[cfg(feature = "full-io")]
fn process_olci_dsf(
    scene: &OlciScene,
    rhot: &HashMap<String, Array2<f64>>,
    data_dir: &std::path::Path,
    config: &OlciProcessingConfig,
    sza: f64, vza: f64, raa: f64, pressure: f64, uoz: f64, uwv: f64,
) -> Result<OlciL2rResult> {
    use crate::ac::{dsf, gas_lut, aerlut};

    let thv = vza.max(0.001);
    let ths = sza;

    // Collect band info
    let wavelengths: Vec<f64> = OLCI_BANDS.iter().map(|(_, wl, _, _)| *wl).collect();
    let bandwidths: Vec<f64> = OLCI_BANDS.iter().map(|(_, _, bw, _)| *bw).collect();

    // Load actual sensor RSR if available
    let sensor_tag = if scene.sensor.contains("S3A") { "S3A_OLCI" } else { "S3B_OLCI" };
    let rsr_path = data_dir.join("RSR").join(format!("{}.txt", sensor_tag));
    let sensor_rsr = aerlut::SensorRsr::load(&rsr_path).ok();
    if sensor_rsr.is_some() {
        log::info!("Using actual {} RSR for LUT convolution", sensor_tag);
    }

    // Gas transmittance — use sensor RSR if available for accurate convolution
    let band_names_str: Vec<&str> = OLCI_BANDS.iter().map(|(n, _, _, _)| *n).collect();
    let tt_gas_vec = if let Some(ref rsr) = sensor_rsr {
        let hyper_tg = gas_lut::compute_gas_transmittance_sensor(
            data_dir, &band_names_str, &rsr.bands, ths, thv, pressure, uoz, uwv,
        ).map_err(|e| crate::AcoliteError::Processing(format!("Gas LUT: {}", e)))?;
        hyper_tg.tt_gas
    } else {
        let hyper_tg = gas_lut::compute_gas_transmittance_hyper(
            data_dir, &wavelengths, &bandwidths, ths, thv, pressure, uoz, uwv,
        ).map_err(|e| crate::AcoliteError::Processing(format!("Gas LUT: {}", e)))?;
        hyper_tg.tt_gas
    };

    log::info!("DSF: sza={:.4} vza={:.4} raa={:.4} p={:.1} uoz={:.3} uwv={:.3}",
        ths, thv, raa, pressure, uoz, uwv);

    // Load generic aerosol LUTs
    let pressures = vec![500.0, 750.0, 1013.0, 1100.0];
    let luts = aerlut::load_generic_luts(data_dir, &pressures)
        .map_err(|e| crate::AcoliteError::Processing(format!("Aerosol LUT: {}", e)))?;

    // Gas-correct TOA
    let toa_gc_bands: Vec<Array2<f64>> = OLCI_BANDS.iter().enumerate().map(|(i, (name, _, _, _))| {
        let tt = tt_gas_vec[i];
        rhot.get(&name.to_string()).map(|toa| {
            toa.mapv(|v| if v.is_finite() && v > 0.0 { v / tt } else { f64::NAN })
        }).unwrap_or_else(|| Array2::zeros((1, 1)))
    }).collect();

    // DSF optimization
    let mut dsf_config = config.dsf.clone();
    if dsf_config.wave_range == (400.0, 2500.0) {
        dsf_config.wave_range = (400.0, 900.0); // OLCI only goes to ~1020nm
    }
    // Force fixed mode when sensor RSR is available (matching Python default for OLCI)
    if sensor_rsr.is_some() {
        dsf_config.mode = DsfMode::Fixed;
    }

    // Build per-band RSR references for sensor-RSR-aware DSF
    let band_names: Vec<String> = OLCI_BANDS.iter().map(|(n, _, _, _)| n.to_string()).collect();
    let band_rsr_data: Vec<Option<(&[f64], &[f64])>> = if let Some(ref rsr) = sensor_rsr {
        band_names.iter().map(|name| {
            rsr.bands.get(name).map(|(w, r)| (w.as_slice(), r.as_slice()))
        }).collect()
    } else {
        vec![None; OLCI_BANDS.len()]
    };

    let (dsf_result, tiled_result) = match dsf_config.mode {
        DsfMode::Fixed => {
            let result = if sensor_rsr.is_some() {
                dsf::optimize_aot_fixed_sensor_rsr(
                    &luts, &toa_gc_bands, &wavelengths, &bandwidths, &tt_gas_vec,
                    &band_rsr_data, pressure, raa, thv, ths, &dsf_config,
                )
            } else {
                dsf::optimize_aot_fixed_generic(
                    &luts, &toa_gc_bands, &wavelengths, &bandwidths, &tt_gas_vec,
                    pressure, raa, thv, ths, &dsf_config,
                )
            };
            (Some(result), None)
        }
        DsfMode::Tiled(tr, tc) => {
            let tr = dsf::optimize_aot_tiled_generic(
                &luts, &toa_gc_bands, &wavelengths, &bandwidths, &tt_gas_vec,
                pressure, raa, thv, ths, &dsf_config, (tr, tc),
            );
            (None, Some(tr))
        }
    };

    let selected_lut = if let Some(ref r) = dsf_result {
        &luts[r.model_idx]
    } else {
        &luts[tiled_result.as_ref().unwrap().model_idx]
    };

    let aot = dsf_result.as_ref().map(|r| r.aot)
        .unwrap_or_else(|| {
            let vals: Vec<f64> = tiled_result.as_ref().unwrap().aot_grid.iter()
                .flatten().copied().filter(|v| v.is_finite()).collect();
            if vals.is_empty() { 0.0 } else { vals.iter().sum::<f64>() / vals.len() as f64 }
        });
    let model_name = dsf_result.as_ref().map(|r| r.model_name.clone())
        .unwrap_or_else(|| tiled_result.as_ref().unwrap().model_name.clone());

    // Apply correction to all bands
    let mut rhos = HashMap::new();
    for (i, (name, wl, bw, _)) in OLCI_BANDS.iter().enumerate() {
        if let Some(toa) = rhot.get(&name.to_string()) {
            let tt = tt_gas_vec[i];
            let rsr_band = band_rsr_data.get(i).and_then(|r| *r);
            let corrected = if let Some(ref r) = dsf_result {
                dsf::dsf_correct_band_generic_rsr(toa, selected_lut, *wl, *bw, rsr_band, tt, r.aot, pressure, raa, thv, ths)
            } else {
                dsf::dsf_correct_band_tiled_generic(
                    toa, selected_lut, *wl, *bw, tt, tiled_result.as_ref().unwrap(), pressure, raa, thv, ths,
                )
            };
            rhos.insert(name.to_string(), corrected);
        }
    }

    Ok(OlciL2rResult {
        rhos,
        metadata: scene.metadata.clone(),
        aot,
        model_name,
    })
}

/// Simplified Rayleigh-only correction (no LUT required)
fn process_olci_simplified(
    scene: &OlciScene,
    rhot: &HashMap<String, Array2<f64>>,
    sza: f64, vza: f64, _raa: f64, pressure: f64, uoz: f64, uwv: f64, mu: f64,
) -> Result<OlciL2rResult> {
    let mut rhos = HashMap::new();

    for (name, wl, _bw, _) in &OLCI_BANDS {
        if let Some(toa) = rhot.get(&name.to_string()) {
            let toa_gc = gas_correction(toa, *wl, uoz, uwv, sza, vza);
            let tau_ray = rayleigh_optical_thickness(*wl, pressure);
            let rho_ray = tau_ray * 0.75 * (1.0 + mu * mu);
            let rhos_band = toa_gc.mapv(|v| {
                let v_minus_ray = v - rho_ray;
                if v_minus_ray > 0.0 { v_minus_ray } else { 0.0 }
            });
            rhos.insert(name.to_string(), rhos_band);
        }
    }

    Ok(OlciL2rResult {
        rhos,
        metadata: scene.metadata.clone(),
        aot: 0.1,
        model_name: "simplified".into(),
    })
}

fn nanmean(arr: &Array2<f64>) -> f64 {
    let (sum, count) = arr.iter().fold((0.0, 0u64), |(s, c), &v| {
        if v.is_finite() { (s + v, c + 1) } else { (s, c) }
    });
    if count > 0 { sum / count as f64 } else { f64::NAN }
}
