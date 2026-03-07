//! Aerosol LUT reader — loads 6SV sensor-specific NetCDF LUTs
//!
//! LUT structure per band: (par, azi, thv, ths, wnd, tau) → f32
//! Parameters: utott, dtott, astot, romix, etc.
//! After loading, we build RegularGridInterpolator per (band, parameter).

use crate::ac::interp::RegularGridInterpolator;
use crate::{AcoliteError, Result};
use std::collections::HashMap;
use std::path::Path;

/// Indices into the parameter dimension of the LUT
#[derive(Debug, Clone)]
pub struct LutParIndex {
    pub romix: usize,
    pub astot: usize,
    pub dutott: usize, // virtual: utott * dtott
}

/// Metadata axes for a LUT
#[derive(Debug, Clone)]
pub struct LutMeta {
    pub par_names: Vec<String>,
    pub azi: Vec<f64>,
    pub thv: Vec<f64>,
    pub ths: Vec<f64>,
    pub wnd: Vec<f64>,
    pub tau: Vec<f64>,
}

/// A loaded aerosol model with multiple pressures stacked.
pub struct AerosolLut {
    pub name: String,
    pub meta: LutMeta,
    pub pressures: Vec<f64>,
    pub ipd: LutParIndex,
    pub band_names: Vec<String>,
    /// Key: (band_name, par_key), Value: RGI over [pressure, azi, thv, ths, tau]
    band_rgi: HashMap<(String, usize), RegularGridInterpolator>,
}

// Parameter keys used as second element of the HashMap key
const PAR_ROMIX: usize = 0;
const PAR_ASTOT: usize = 1;
const PAR_DUTOTT: usize = 2;

impl AerosolLut {
    pub fn romix(&self, band: &str, pressure: f64, azi: f64, thv: f64, ths: f64, tau: f64) -> f64 {
        self.band_rgi.get(&(band.to_string(), PAR_ROMIX))
            .map(|rgi| rgi.interpolate(&[pressure, azi, thv, ths, tau]))
            .unwrap_or(f64::NAN)
    }
    pub fn astot(&self, band: &str, pressure: f64, azi: f64, thv: f64, ths: f64, tau: f64) -> f64 {
        self.band_rgi.get(&(band.to_string(), PAR_ASTOT))
            .map(|rgi| rgi.interpolate(&[pressure, azi, thv, ths, tau]))
            .unwrap_or(f64::NAN)
    }
    pub fn dutott(&self, band: &str, pressure: f64, azi: f64, thv: f64, ths: f64, tau: f64) -> f64 {
        self.band_rgi.get(&(band.to_string(), PAR_DUTOTT))
            .map(|rgi| rgi.interpolate(&[pressure, azi, thv, ths, tau]))
            .unwrap_or(f64::NAN)
    }
}

/// Load sensor-specific aerosol LUTs for a given model across multiple pressures.
#[cfg(feature = "full-io")]
pub fn load_sensor_lut(
    lut_dir: &Path,
    base_name: &str,
    sensor: &str,
    pressures: &[f64],
) -> Result<AerosolLut> {
    use netcdf::AttributeValue;

    // Derive the LUT subdirectory: "ACOLITE-LUT-202110-MOD2" → "ACOLITE-LUT-202110"
    let parts: Vec<&str> = base_name.rsplitn(2, '-').collect();
    let lut_subdir = if parts.len() == 2 { parts[1] } else { base_name };
    let sensor_dir = lut_dir.join(lut_subdir).join(sensor);

    // Read first pressure file for metadata
    let first_id = format!("{}-{:04}mb_{}", base_name, pressures[0] as i32, sensor);
    let first_path = sensor_dir.join(format!("{}.nc", first_id));
    if !first_path.exists() {
        return Err(AcoliteError::Processing(format!("LUT not found: {}", first_path.display())));
    }

    let nc = netcdf::open(&first_path)
        .map_err(|e| AcoliteError::Processing(format!("NetCDF open: {}", e)))?;

    // Read par attribute
    let par_str: String = match nc.attribute("par").and_then(|a| a.value().ok()) {
        Some(AttributeValue::Str(s)) => s,
        _ => return Err(AcoliteError::Processing("Missing 'par' attribute".into())),
    };
    let par_names: Vec<String> = par_str.split(',').map(|s| s.to_string()).collect();

    let azi = read_f64_attr_nc(&nc, "azi")?;
    let thv = read_f64_attr_nc(&nc, "thv")?;
    let ths = read_f64_attr_nc(&nc, "ths")?;
    let tau = read_f64_attr_nc(&nc, "tau")?;
    let wnd = read_f64_attr_nc(&nc, "wnd").unwrap_or_else(|_| vec![2.0]);

    // Find parameter indices
    let find_par = |name: &str| -> Result<usize> {
        par_names.iter().position(|p| p == name)
            .ok_or_else(|| AcoliteError::Processing(format!("'{}' not in LUT par", name)))
    };
    let romix_idx = find_par("romix")?;
    let astot_idx = find_par("astot")?;
    let utott_idx = find_par("utott")?;
    let dtott_idx = find_par("dtott")?;

    // Get band names
    let band_names: Vec<String> = nc.variables().map(|v| v.name().to_string()).collect();

    let npar = par_names.len();
    let nazi = azi.len();
    let nthv = thv.len();
    let nths = ths.len();
    let ntau = tau.len();
    let npres = pressures.len();

    drop(nc);

    // For each band, accumulate [pressure, azi, thv, ths, tau] data for romix, astot, dutott
    let vol = npres * nazi * nthv * nths * ntau;
    let mut band_romix: HashMap<String, Vec<f32>> = HashMap::new();
    let mut band_astot: HashMap<String, Vec<f32>> = HashMap::new();
    let mut band_dutott: HashMap<String, Vec<f32>> = HashMap::new();
    for b in &band_names {
        band_romix.insert(b.clone(), Vec::with_capacity(vol));
        band_astot.insert(b.clone(), Vec::with_capacity(vol));
        band_dutott.insert(b.clone(), Vec::with_capacity(vol));
    }

    for &pres in pressures {
        let lutid = format!("{}-{:04}mb_{}", base_name, pres as i32, sensor);
        let path = sensor_dir.join(format!("{}.nc", lutid));
        let nc = netcdf::open(&path)
            .map_err(|e| AcoliteError::Processing(format!("NetCDF open {}: {}", path.display(), e)))?;

        for band in &band_names {
            let var = nc.variable(band)
                .ok_or_else(|| AcoliteError::Processing(format!("Band '{}' not in LUT", band)))?;
            // Shape: (npar, nazi, nthv, nths, nwnd=1, ntau)
            let raw: ndarray::ArrayD<f32> = var.get::<f32, _>(..)
                .map_err(|e| AcoliteError::Processing(format!("Read {}: {}", band, e)))?;
            let raw = raw.into_raw_vec();

            let stride_par = nazi * nthv * nths * 1 * ntau;
            let stride_azi = nthv * nths * 1 * ntau;
            let stride_thv = nths * 1 * ntau;
            let stride_ths = 1 * ntau;

            let romix_v = band_romix.get_mut(band)
                .ok_or_else(|| AcoliteError::Processing(format!("Missing romix for band {}", band)))?;
            let astot_v = band_astot.get_mut(band)
                .ok_or_else(|| AcoliteError::Processing(format!("Missing astot for band {}", band)))?;
            let dutott_v = band_dutott.get_mut(band)
                .ok_or_else(|| AcoliteError::Processing(format!("Missing dutott for band {}", band)))?;

            for ia in 0..nazi {
                for iv in 0..nthv {
                    for is_ in 0..nths {
                        for it in 0..ntau {
                            let base = ia * stride_azi + iv * stride_thv + is_ * stride_ths + it;
                            romix_v.push(raw[romix_idx * stride_par + base]);
                            astot_v.push(raw[astot_idx * stride_par + base]);
                            let u = raw[utott_idx * stride_par + base];
                            let d = raw[dtott_idx * stride_par + base];
                            dutott_v.push(u * d);
                        }
                    }
                }
            }
        }
    }

    // Build interpolators
    let pres_axis = pressures.to_vec();
    let mut band_rgi: HashMap<(String, usize), RegularGridInterpolator> = HashMap::new();

    for band in &band_names {
        let axes = vec![pres_axis.clone(), azi.clone(), thv.clone(), ths.clone(), tau.clone()];
        band_rgi.insert(
            (band.clone(), PAR_ROMIX),
            RegularGridInterpolator::new(axes.clone(), band_romix.remove(band)
                .ok_or_else(|| AcoliteError::Processing(format!("Missing romix data for {}", band)))?),
        );
        band_rgi.insert(
            (band.clone(), PAR_ASTOT),
            RegularGridInterpolator::new(axes.clone(), band_astot.remove(band)
                .ok_or_else(|| AcoliteError::Processing(format!("Missing astot data for {}", band)))?),
        );
        band_rgi.insert(
            (band.clone(), PAR_DUTOTT),
            RegularGridInterpolator::new(axes, band_dutott.remove(band)
                .ok_or_else(|| AcoliteError::Processing(format!("Missing dutott data for {}", band)))?),
        );
    }

    Ok(AerosolLut {
        name: base_name.to_string(),
        meta: LutMeta { par_names, azi, thv, ths, wnd, tau },
        pressures: pressures.to_vec(),
        ipd: LutParIndex { romix: PAR_ROMIX, astot: PAR_ASTOT, dutott: PAR_DUTOTT },
        band_names,
        band_rgi,
    })
}

#[cfg(feature = "full-io")]
fn read_f64_attr_nc(nc: &netcdf::File, name: &str) -> Result<Vec<f64>> {
    use netcdf::AttributeValue;
    let attr = nc.attribute(name)
        .ok_or_else(|| AcoliteError::Processing(format!("Missing attr '{}'", name)))?;
    let val = attr.value()
        .map_err(|e| AcoliteError::Processing(format!("Read attr '{}': {}", name, e)))?;
    match val {
        AttributeValue::Double(v) => Ok(vec![v]),
        AttributeValue::Float(v) => Ok(vec![v as f64]),
        AttributeValue::Doubles(v) => Ok(v),
        AttributeValue::Floats(v) => Ok(v.iter().map(|&x| x as f64).collect()),
        AttributeValue::Int(v) => Ok(vec![v as f64]),
        AttributeValue::Ints(v) => Ok(v.iter().map(|&x| x as f64).collect()),
        _ => Err(AcoliteError::Processing(format!("Unexpected type for '{}'", name))),
    }
}

/// Load both MOD1 and MOD2 LUTs for a sensor
#[cfg(feature = "full-io")]
pub fn load_acolite_luts(
    data_dir: &Path,
    sensor: &str,
    pressures: &[f64],
) -> Result<Vec<AerosolLut>> {
    let lut_dir = data_dir.join("LUT");
    let mut luts = Vec::new();
    for model in &["ACOLITE-LUT-202110-MOD1", "ACOLITE-LUT-202110-MOD2"] {
        match load_sensor_lut(&lut_dir, model, sensor, pressures) {
            Ok(lut) => {
                log::info!("Loaded LUT {} for {}", model, sensor);
                luts.push(lut);
            }
            Err(e) => log::warn!("Could not load {}: {}", model, e),
        }
    }
    if luts.is_empty() {
        return Err(AcoliteError::Processing("No aerosol LUTs loaded".into()));
    }
    Ok(luts)
}

/// Load a single aerosol model LUT (e.g. "MOD1" or "MOD2")
#[cfg(feature = "full-io")]
pub fn load_single_lut(
    data_dir: &Path,
    sensor: &str,
    model: &str,
    pressures: &[f64],
) -> Result<Vec<AerosolLut>> {
    let lut_dir = data_dir.join("LUT");
    let full_name = format!("ACOLITE-LUT-202110-{}", model);
    let lut = load_sensor_lut(&lut_dir, &full_name, sensor, pressures)?;
    log::info!("Loaded single LUT {} for {}", full_name, sensor);
    Ok(vec![lut])
}
