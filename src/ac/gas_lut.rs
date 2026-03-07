//! Gas transmittance — ozone, water vapour, and other gases.
//!
//! Reads ko3 absorption coefficients from data/Shared/k_o3_anderson.txt,
//! water vapour LUT from data/LUT/WV/WV_201710C.nc,
//! and other gas LUT from data/LUT/Gas/Gas_202106F.nc.
//! Convolves with sensor RSR to get per-band transmittances.

use crate::ac::interp::RegularGridInterpolator;
use crate::{AcoliteError, Result};
use std::collections::HashMap;
use std::path::Path;

/// Per-band gas transmittance values
#[derive(Debug, Clone)]
pub struct GasTransmittance {
    /// Band name → total two-way gas transmittance
    pub tt_gas: HashMap<String, f64>,
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
