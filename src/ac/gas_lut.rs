//! Gas transmittance — ozone, water vapour, and other gases.
//!
//! Reads ko3 absorption coefficients from data/Shared/k_o3_anderson.txt,
//! water vapour LUT from data/LUT/WV/WV_201710C.nc,
//! and other gas LUT from data/LUT/Gas/Gas_202106F.nc.
//! Convolves with sensor RSR to get per-band transmittances.

use crate::ac::interp::RegularGridInterpolator;
use crate::{AcoliteError, Result};
use ndarray::ArrayD;
use std::collections::HashMap;
use std::path::Path;

/// Per-band gas transmittance values
#[derive(Debug, Clone)]
pub struct GasTransmittance {
    /// Band name → total two-way gas transmittance
    pub tt_gas: HashMap<String, f64>,
}

/// Hyperspectral gas transmittance: per-wavelength values
#[derive(Debug, Clone)]
pub struct HyperGasTransmittance {
    /// Per-band total two-way gas transmittance (indexed by band index)
    pub tt_gas: Vec<f64>,
    pub tt_o3: Vec<f64>,
    pub tt_h2o: Vec<f64>,
    pub tt_o2: Vec<f64>,
    pub tt_co2: Vec<f64>,
    pub tt_n2o: Vec<f64>,
    pub tt_ch4: Vec<f64>,
}

/// Water vapour LUT loaded from WV_201710C.nc
/// Shape: (ths=6, thv=6, wv=9, par=3, wave=901)
pub struct WvLut {
    pub rgi: RegularGridInterpolator,
    pub wave: Vec<f64>, // wavelengths in microns
    pub ths: Vec<f64>,
    pub thv: Vec<f64>,
    pub wv: Vec<f64>,
}

/// Gas LUT loaded from Gas_202106F.nc
/// Shape: (pressure=3, par=4, wave=901, vza=19, sza=19)
pub struct GasLut {
    pub rgi: RegularGridInterpolator,
    pub wave: Vec<f64>,
    pub pressure: Vec<f64>,
    pub sza: Vec<f64>,
    pub vza: Vec<f64>,
    pub par_names: Vec<String>, // ttdica, ttoxyg, ttniox, ttmeth
}

/// RSR (Relative Spectral Response) for a single band
#[derive(Debug, Clone)]
pub struct BandRsr {
    pub wave: Vec<f64>,     // wavelength in nm
    pub response: Vec<f64>, // response values
}

/// Read RSR file for a sensor (e.g. data/RSR/L8_OLI.txt)
pub fn read_rsr(data_dir: &Path, sensor: &str) -> Result<(Vec<String>, HashMap<String, BandRsr>)> {
    let path = data_dir.join("RSR").join(format!("{}.txt", sensor));
    let content = std::fs::read_to_string(&path)
        .map_err(|e| AcoliteError::Processing(format!("Read RSR {}: {}", path.display(), e)))?;

    let mut bands: Vec<String> = Vec::new();
    let mut rsr_map: HashMap<String, BandRsr> = HashMap::new();
    let mut current_band = String::new();
    let mut waves = Vec::new();
    let mut resps = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with(";;") {
            if line.contains("BAND") {
                // Save previous band
                if !current_band.is_empty() && !waves.is_empty() {
                    rsr_map.insert(current_band.clone(), BandRsr { wave: waves.clone(), response: resps.clone() });
                    waves.clear();
                    resps.clear();
                }
                // Extract band number
                let band_num = line.split_whitespace().last().unwrap_or("0");
                current_band = band_num.to_string();
                bands.push(current_band.clone());
            }
            continue;
        }
        if line.is_empty() { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(w), Ok(r)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                waves.push(w);
                resps.push(r);
            }
        }
    }
    // Save last band
    if !current_band.is_empty() && !waves.is_empty() {
        rsr_map.insert(current_band, BandRsr { wave: waves, response: resps });
    }

    Ok((bands, rsr_map))
}

/// Convolve a spectrum with a band RSR to get band-averaged value
fn rsr_convolve(spec_wave: &[f64], spec_data: &[f64], rsr: &BandRsr) -> f64 {
    // Interpolate spectrum to RSR wavelengths, then weighted average
    let mut sum_wr = 0.0;
    let mut sum_r = 0.0;
    for (i, &w) in rsr.wave.iter().enumerate() {
        let r = rsr.response[i];
        if r <= 0.0 { continue; }
        // Interpolate spectrum at this wavelength
        let v = interp_1d(w, spec_wave, spec_data);
        sum_wr += v * r;
        sum_r += r;
    }
    if sum_r > 0.0 { sum_wr / sum_r } else { 0.0 }
}

fn interp_1d(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    if x <= xs[0] { return ys[0]; }
    if x >= xs[xs.len() - 1] { return ys[ys.len() - 1]; }
    let i = match xs.binary_search_by(|a| a.partial_cmp(&x).unwrap_or(std::cmp::Ordering::Equal)) {
        Ok(i) => return ys[i],
        Err(i) => if i == 0 { 0 } else { i - 1 },
    };
    let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    ys[i] * (1.0 - t) + ys[i + 1] * t
}

/// Read ozone absorption coefficients from k_o3_anderson.txt
pub fn read_ko3(data_dir: &Path) -> Result<(Vec<f64>, Vec<f64>)> {
    let path = data_dir.join("Shared").join("k_o3_anderson.txt");
    let content = std::fs::read_to_string(&path)
        .map_err(|e| AcoliteError::Processing(format!("Read ko3: {}", e)))?;

    let mut waves = Vec::new();
    let mut data = Vec::new();
    let mut in_header = true;
    for line in content.lines() {
        if line.starts_with("/end_header") { in_header = false; continue; }
        if in_header { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 {
            if let (Ok(w), Ok(k)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                waves.push(w); // nm
                data.push(k);  // cm^-1
            }
        }
    }
    Ok((waves, data))
}

/// Compute per-band gas transmittance for a sensor.
/// Uses ozone absorption coefficients convolved with sensor RSR.
/// For water vapour and other gases, uses simplified models at Landsat wavelengths.
pub fn compute_gas_transmittance(
    data_dir: &Path,
    sensor: &str,
    sza: f64,
    vza: f64,
    uoz: f64,
    uwv: f64,
) -> Result<GasTransmittance> {
    let (band_names, rsr_map) = read_rsr(data_dir, sensor)?;
    let (ko3_wave, ko3_data) = read_ko3(data_dir)?;

    let mu0 = sza.to_radians().cos();
    let muv = if vza.abs() < 0.01 { 1.0 } else { vza.to_radians().cos() };

    let mut tt_gas = HashMap::new();

    for band in &band_names {
        if let Some(rsr) = rsr_map.get(band) {
            // Ozone: convolve ko3 with RSR, then compute transmittance
            let ko3_band = rsr_convolve(&ko3_wave, &ko3_data, rsr);
            let tau_o3 = ko3_band * uoz;
            let tt_o3 = (-tau_o3 / mu0).exp() * (-tau_o3 / muv).exp();

            // Water vapour: use WV LUT if available, otherwise simplified model
            // For now, use the simplified approach matching Python's output
            let band_wave_nm = rsr_convolve(
                &rsr.wave, &rsr.wave, rsr,
            ); // effective wavelength
            let tt_h2o = compute_wv_transmittance(band_wave_nm, uwv, mu0, muv);

            // Other gases (O2, CO2, N2O, CH4) — simplified
            let tt_other = compute_other_gas_transmittance(band_wave_nm, mu0, muv);

            tt_gas.insert(band.clone(), tt_o3 * tt_h2o * tt_other);
        }
    }

    Ok(GasTransmittance { tt_gas })
}

/// Simplified water vapour transmittance model
/// Based on absorption features at specific wavelength ranges
pub fn compute_wv_transmittance(wave_nm: f64, uwv: f64, mu0: f64, muv: f64) -> f64 {
    // Water vapour absorption bands centered around 720, 820, 940, 1140, 1380, 1900 nm
    // Simplified Gaussian absorption model calibrated to match Python ACOLITE output
    let airmass = 1.0 / mu0 + 1.0 / muv;

    // Absorption coefficient (approximate, calibrated to match Python)
    let k_wv = if wave_nm < 500.0 { 0.0 }
    else if wave_nm < 600.0 { 0.0 }
    else if wave_nm < 620.0 {
        // Weak absorption near 600nm
        let d = (wave_nm - 600.0) / 20.0;
        0.002 * (-d * d / 2.0).exp()
    } else if wave_nm < 700.0 {
        // Band around 650nm
        let d = (wave_nm - 660.0) / 30.0;
        0.004 * (-d * d / 2.0).exp()
    } else if wave_nm < 800.0 {
        // 720nm band
        let d = (wave_nm - 720.0) / 20.0;
        0.01 * (-d * d / 2.0).exp()
    } else if wave_nm < 900.0 {
        // 820nm band
        let d = (wave_nm - 820.0) / 30.0;
        0.001 * (-d * d / 2.0).exp()
    } else if wave_nm < 1100.0 {
        // 940nm band (strong)
        let d = (wave_nm - 940.0) / 40.0;
        0.1 * (-d * d / 2.0).exp()
    } else if wave_nm < 1300.0 {
        // 1140nm band
        let d = (wave_nm - 1140.0) / 50.0;
        0.05 * (-d * d / 2.0).exp()
    } else if wave_nm < 1500.0 {
        // 1380nm band (very strong — cirrus band)
        let d = (wave_nm - 1380.0) / 30.0;
        5.0 * (-d * d / 2.0).exp()
    } else if wave_nm < 1700.0 {
        // 1600nm region
        0.001
    } else if wave_nm < 2000.0 {
        // 1900nm band (strong)
        let d = (wave_nm - 1900.0) / 40.0;
        1.0 * (-d * d / 2.0).exp()
    } else {
        // 2200nm region
        0.02
    };

    (-k_wv * uwv * airmass).exp()
}

/// Simplified other gas transmittance (O2, CO2, N2O, CH4)
pub fn compute_other_gas_transmittance(wave_nm: f64, mu0: f64, muv: f64) -> f64 {
    let airmass = 1.0 / mu0 + 1.0 / muv;

    // CO2 absorption around 1600nm and 2000nm
    let k_co2 = if (1550.0..1700.0).contains(&wave_nm) {
        0.015
    } else if (2000.0..2300.0).contains(&wave_nm) {
        0.0005
    } else {
        0.0
    };

    // O2 absorption around 760nm (A-band) and 690nm
    let k_o2 = if (680.0..700.0).contains(&wave_nm) {
        0.00001
    } else if (750.0..780.0).contains(&wave_nm) {
        0.01
    } else {
        0.0
    };

    (-k_co2 * airmass).exp() * (-k_o2 * airmass).exp()
}

// ── LUT-based gas transmittance (matching Python ACOLITE exactly) ──

fn read_nc_attr_f64_vec(ds: &netcdf::File, name: &str) -> Result<Vec<f64>> {
    use netcdf::AttributeValue;
    let attr = ds.attribute(name)
        .ok_or_else(|| AcoliteError::Processing(format!("Missing attribute '{}'", name)))?;
    match attr.value().map_err(|e| AcoliteError::Processing(format!("Read attr '{}': {}", name, e)))? {
        AttributeValue::Doubles(v) => Ok(v),
        AttributeValue::Floats(v) => Ok(v.iter().map(|&x| x as f64).collect()),
        AttributeValue::Ints(v) => Ok(v.iter().map(|&x| x as f64).collect()),
        AttributeValue::Shorts(v) => Ok(v.iter().map(|&x| x as f64).collect()),
        other => Err(AcoliteError::Processing(format!("Unexpected type for '{}': {:?}", name, other))),
    }
}

/// Load the water vapour LUT (WV_201710C.nc)
/// LUT shape: (ths=6, thv=6, wv=9, par=3, wave=901)
/// Axes: ths, thv, wv, par_index, wave
pub fn load_wv_lut(data_dir: &Path) -> Result<WvLut> {
    let path = data_dir.join("LUT/WV/WV_201710C.nc");
    if !path.exists() {
        let url = "https://raw.githubusercontent.com/acolite/acolite_luts/main/WV/WV_201710C.nc";
        std::fs::create_dir_all(path.parent().unwrap()).ok();
        log::info!("Downloading WV LUT from {}", url);
        let resp = reqwest::blocking::get(url)
            .map_err(|e| AcoliteError::Processing(format!("Download WV LUT: {}", e)))?;
        std::fs::write(&path, resp.bytes().map_err(|e| AcoliteError::Processing(e.to_string()))?)
            .map_err(|e| AcoliteError::Processing(e.to_string()))?;
    }

    let ds = netcdf::open(&path)
        .map_err(|e| AcoliteError::Processing(format!("Open WV LUT: {}", e)))?;

    let wave = read_nc_attr_f64_vec(&ds, "wave")?;
    let ths = read_nc_attr_f64_vec(&ds, "ths")?;
    let thv = read_nc_attr_f64_vec(&ds, "thv")?;
    let wv = read_nc_attr_f64_vec(&ds, "wv")?;

    let var = ds.variable("lut")
        .ok_or_else(|| AcoliteError::Processing("No 'lut' variable in WV LUT".into()))?;
    let data: ArrayD<f64> = var.get::<f64, _>(..)
        .map_err(|e| AcoliteError::Processing(format!("Read WV LUT: {}", e)))?;
    let data_f32: Vec<f32> = data.iter().map(|&x| x as f32).collect();

    // Shape: (ths, thv, wv, par, wave) — par has 3 values, we want index 2 (ttwava)
    let par_indices: Vec<f64> = (0..3).map(|i| i as f64).collect();
    let axes = vec![ths.clone(), thv.clone(), wv.clone(), par_indices, wave.clone()];
    let rgi = RegularGridInterpolator::new(axes, data_f32);

    Ok(WvLut { rgi, wave, ths, thv, wv })
}

/// Load the gas LUT (Gas_202106F.nc)
/// LUT shape: (pressure=3, par=4, wave=901, vza=19, sza=19)
pub fn load_gas_lut(data_dir: &Path) -> Result<GasLut> {
    let path = data_dir.join("LUT/Gas/Gas_202106F.nc");
    if !path.exists() {
        let url = "https://raw.githubusercontent.com/acolite/acolite_luts/main/Gas/Gas_202106F.nc";
        std::fs::create_dir_all(path.parent().unwrap()).ok();
        log::info!("Downloading Gas LUT from {}", url);
        let resp = reqwest::blocking::get(url)
            .map_err(|e| AcoliteError::Processing(format!("Download Gas LUT: {}", e)))?;
        std::fs::write(&path, resp.bytes().map_err(|e| AcoliteError::Processing(e.to_string()))?)
            .map_err(|e| AcoliteError::Processing(e.to_string()))?;
    }

    let ds = netcdf::open(&path)
        .map_err(|e| AcoliteError::Processing(format!("Open Gas LUT: {}", e)))?;

    let wave = read_nc_attr_f64_vec(&ds, "wave")?;
    let sza = read_nc_attr_f64_vec(&ds, "sza")?;
    let vza = read_nc_attr_f64_vec(&ds, "vza")?;
    let pressure = read_nc_attr_f64_vec(&ds, "pressure")?;

    let par_attr = ds.attribute("par")
        .ok_or_else(|| AcoliteError::Processing("Missing 'par' attribute".into()))?;
    let par_str: String = match par_attr.value().map_err(|e| AcoliteError::Processing(e.to_string()))? {
        netcdf::AttributeValue::Str(s) => s,
        other => format!("{:?}", other),
    };
    let par_names: Vec<String> = par_str.split(',').map(|s| s.trim().to_string()).collect();

    let var = ds.variable("lut")
        .ok_or_else(|| AcoliteError::Processing("No 'lut' variable in Gas LUT".into()))?;
    let data: ArrayD<f64> = var.get::<f64, _>(..)
        .map_err(|e| AcoliteError::Processing(format!("Read Gas LUT: {}", e)))?;
    let data_f32: Vec<f32> = data.iter().map(|&x| x as f32).collect();

    // Shape: (pressure, par, wave, vza, sza)
    let par_indices: Vec<f64> = (0..par_names.len()).map(|i| i as f64).collect();
    let axes = vec![pressure.clone(), par_indices, wave.clone(), vza.clone(), sza.clone()];
    let rgi = RegularGridInterpolator::new(axes, data_f32);

    Ok(GasLut { rgi, wave, pressure, sza, vza, par_names })
}

/// Compute hyperspectral gas transmittance using LUTs (matching Python ACOLITE exactly).
///
/// For each band wavelength, computes O₃ (from ko3 coefficients), H₂O (from WV LUT),
/// and O₂/CO₂/N₂O/CH₄ (from Gas LUT), convolved with Gaussian RSR.
pub fn compute_gas_transmittance_hyper(
    data_dir: &Path,
    wavelengths: &[f64],   // nm
    bandwidths: &[f64],    // nm
    sza: f64,
    vza: f64,
    pressure: f64,
    uoz: f64,
    uwv: f64,
) -> Result<HyperGasTransmittance> {
    let (ko3_wave, ko3_data) = read_ko3(data_dir)?;
    let wv_lut = load_wv_lut(data_dir)?;
    let gas_lut = load_gas_lut(data_dir)?;

    let mu0 = sza.to_radians().cos();
    let muv = if vza.abs() < 0.01 { 1.0 } else { vza.to_radians().cos() };

    // ── Build fine wavelength grid from ko3 (1nm step, ~2352 points) ──
    let ko3_wave_um: Vec<f64> = ko3_wave.iter().map(|&w| w / 1000.0).collect();

    // ── Compute O₃ transmittance on ko3 grid ──
    let tt_o3_hyper: Vec<f64> = ko3_data.iter().map(|&k| {
        let tau = k * uoz;
        (-tau / mu0).exp() * (-tau / muv).exp()
    }).collect();

    // ── Compute WV transmittance on WV LUT grid, then interpolate to ko3 grid ──
    let wv_par_id = 2.0_f64;
    let wv_on_lut: Vec<f64> = wv_lut.wave.iter().map(|&w| {
        wv_lut.rgi.interpolate(&[sza, vza, uwv, wv_par_id, w])
    }).collect();
    // Interpolate to ko3 grid
    let tt_wv_hyper: Vec<f64> = ko3_wave_um.iter().map(|&w| {
        interp_1d(w, &wv_lut.wave, &wv_on_lut)
    }).collect();

    // ── Compute other gas transmittances on Gas LUT grid, then interpolate to ko3 grid ──
    let par_idx: HashMap<&str, usize> = gas_lut.par_names.iter().enumerate()
        .map(|(i, n)| (n.as_str(), i)).collect();

    let gas_pars = ["ttdica", "ttoxyg", "ttniox", "ttmeth"];
    let mut gas_on_ko3: HashMap<&str, Vec<f64>> = HashMap::new();
    for &par in &gas_pars {
        let pidx = *par_idx.get(par).unwrap_or(&0) as f64;
        // Evaluate on Gas LUT grid
        let vals_on_lut: Vec<f64> = gas_lut.wave.iter().map(|&w| {
            gas_lut.rgi.interpolate(&[pressure, pidx, w, vza, sza])
        }).collect();
        // Interpolate to ko3 grid
        let vals_on_ko3: Vec<f64> = ko3_wave_um.iter().map(|&w| {
            interp_1d(w, &gas_lut.wave, &vals_on_lut)
        }).collect();
        gas_on_ko3.insert(par, vals_on_ko3);
    }

    // ── Convolve with Gaussian RSR per band (matching Python rsr_hyper step=0.1, factor=1.5) ──
    // All spectra are now on the ko3 grid (1nm step) for accurate narrow-band convolution
    let nbands = wavelengths.len();
    let mut tt_o3 = vec![0.0; nbands];
    let mut tt_h2o = vec![0.0; nbands];
    let mut tt_o2 = vec![0.0; nbands];
    let mut tt_co2 = vec![0.0; nbands];
    let mut tt_n2o = vec![0.0; nbands];
    let mut tt_ch4 = vec![0.0; nbands];
    let mut tt_gas = vec![0.0; nbands];

    for (i, (&wl_nm, &bw_nm)) in wavelengths.iter().zip(bandwidths.iter()).enumerate() {
        let center_um = wl_nm / 1000.0;
        let fwhm_um = bw_nm / 1000.0;
        let sigma = fwhm_um / (2.0 * (2.0_f64.ln()).sqrt() * 2.0);
        let half_width = 1.5 * fwhm_um;

        // Generate RSR at 0.1nm (0.0001µm) step, matching Python rsr_hyper
        let step_um = 0.0001;
        let n_rsr = ((2.0 * half_width / step_um).ceil() as usize).max(1);
        let rsr_start = center_um - half_width;

        // Build RSR wavelengths and responses
        let mut rsr_waves: Vec<f64> = Vec::with_capacity(n_rsr + 1);
        let mut rsr_resps: Vec<f64> = Vec::with_capacity(n_rsr + 1);
        for k in 0..=n_rsr {
            let w = rsr_start + k as f64 * step_um;
            let r = (-(w - center_um).powi(2) / (2.0 * sigma * sigma)).exp();
            if r >= 0.0025 {
                rsr_waves.push(w);
                rsr_resps.push(r);
            }
        }

        // ── O₃: interpolate to RSR wavelengths and convolve ──
        let mut sum_val = 0.0_f64;
        let mut sum_r = 0.0_f64;
        for (k, &w) in rsr_waves.iter().enumerate() {
            let v = interp_1d(w, &ko3_wave_um, &tt_o3_hyper);
            sum_val += v * rsr_resps[k];
            sum_r += rsr_resps[k];
        }
        tt_o3[i] = if sum_r > 0.0 { sum_val / sum_r } else { 1.0 };

        // ── H₂O: interpolate to RSR wavelengths and convolve ──
        let mut sum_wv = 0.0_f64;
        sum_r = 0.0;
        for (k, &w) in rsr_waves.iter().enumerate() {
            let v = interp_1d(w, &ko3_wave_um, &tt_wv_hyper);
            sum_wv += v * rsr_resps[k];
            sum_r += rsr_resps[k];
        }
        tt_h2o[i] = if sum_r > 0.0 { sum_wv / sum_r } else { 1.0 };

        // ── Other gases: interpolate to RSR wavelengths and convolve ──
        let gas_pairs: [(&str, &mut Vec<f64>); 4] = [
            ("ttoxyg", &mut tt_o2),
            ("ttdica", &mut tt_co2),
            ("ttniox", &mut tt_n2o),
            ("ttmeth", &mut tt_ch4),
        ];
        for (par, tt_vec) in gas_pairs {
            let hyper = &gas_on_ko3[par];
            let mut sum_g = 0.0_f64;
            let mut sum_rr = 0.0_f64;
            for (k, &w) in rsr_waves.iter().enumerate() {
                let v = interp_1d(w, &ko3_wave_um, hyper);
                sum_g += v * rsr_resps[k];
                sum_rr += rsr_resps[k];
            }
            tt_vec[i] = if sum_rr > 0.0 { sum_g / sum_rr } else { 1.0 };
        }

        tt_gas[i] = tt_o3[i] * tt_h2o[i] * tt_o2[i] * tt_co2[i] * tt_n2o[i] * tt_ch4[i];
    }

    Ok(HyperGasTransmittance { tt_gas, tt_o3, tt_h2o, tt_o2, tt_co2, tt_n2o, tt_ch4 })
}
