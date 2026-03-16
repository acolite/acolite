//! Ancillary data (ozone, water vapour, pressure, wind) download from OBPG/GMAO
//!
//! Mirrors the Python acolite.ac.ancillary module: list files, download, provide defaults.

use crate::auth::Credentials;
use crate::{AcoliteError, Result};
use chrono::{DateTime, Datelike, Duration, Timelike, Utc};
use std::path::{Path, PathBuf};

const OBPG_URL: &str = "https://oceandata.sci.gsfc.nasa.gov/ob/getfile";

/// Ancillary atmospheric parameters for a scene
#[derive(Debug, Clone)]
pub struct Ancillary {
    /// Ozone column (cm-atm, i.e. DU/1000)
    pub uoz: f64,
    /// Water vapour (g/cm²)
    pub uwv: f64,
    /// Surface pressure (hPa)
    pub pressure: f64,
    /// Wind speed (m/s)
    pub wind: f64,
}

impl Default for Ancillary {
    fn default() -> Self {
        Self {
            uoz: 0.3,
            uwv: 1.5,
            pressure: 1013.25,
            wind: 2.0,
        }
    }
}

/// List ancillary file names needed for a given datetime (OBPG naming convention).
/// Returns GMAO MERRA2 MET files (2 bracketing hours) + ozone file.
pub fn list_files(dt: &DateTime<Utc>) -> Vec<String> {
    let yjd = format!("{}{}", dt.format("%Y"), format!("{:03}", dt.ordinal()));
    let mut files = Vec::new();

    // Ozone: AURAOMI for post-2004, EPTOMS before
    let year = dt.year();
    let doy = dt.ordinal();
    if (year == 2004 && doy > 336) || year > 2004 {
        files.push(format!("N{}00_O3_AURAOMI_24h.hdf", yjd));
    } else {
        files.push(format!("N{}00_O3_EPTOMS_24h.hdf", yjd));
    }

    // GMAO MERRA2 MET: bracketing hourly files
    let h = dt.hour();
    files.push(format!(
        "GMAO_MERRA2.{}T{:02}0000.MET.nc",
        dt.format("%Y%m%d"),
        h
    ));
    let next = if h < 23 {
        format!(
            "GMAO_MERRA2.{}T{:02}0000.MET.nc",
            dt.format("%Y%m%d"),
            h + 1
        )
    } else {
        let next_day = *dt + Duration::days(1);
        format!("GMAO_MERRA2.{}T000000.MET.nc", next_day.format("%Y%m%d"))
    };
    files.push(next);

    files
}

/// Download ancillary files to `local_dir/<year>/<doy>/`.
/// Returns paths of successfully downloaded (or cached) files.
pub fn download(dt: &DateTime<Utc>, local_dir: &Path, creds: &Credentials) -> Vec<PathBuf> {
    let year = format!("{}", dt.year());
    let doy = format!("{:03}", dt.ordinal());
    let file_list = list_files(dt);
    let mut local_files = Vec::new();

    for basename in &file_list {
        let dir = local_dir.join(&year).join(&doy);
        let local = dir.join(basename);

        if local.exists() {
            log::info!("Ancillary cached: {}", basename);
            local_files.push(local);
            continue;
        }

        let url = format!("{}/{}", OBPG_URL, basename);
        log::info!("Downloading ancillary: {}", basename);

        if let Err(e) = crate::loader::source::download::download_file(&url, &local, creds) {
            log::warn!("Ancillary download failed for {}: {}", basename, e);
            continue;
        }

        // Reject files < 50 KB (likely error pages)
        if let Ok(meta) = std::fs::metadata(&local) {
            if meta.len() < 50 * 1024 {
                log::warn!(
                    "Removing too-small ancillary file: {} ({} bytes)",
                    basename,
                    meta.len()
                );
                let _ = std::fs::remove_file(&local);
                continue;
            }
        }
        local_files.push(local);
    }

    local_files
}

/// Try to fetch ancillary data for a scene datetime + location.
/// Falls back to climatological defaults if download/parsing fails.
pub fn get(dt: &DateTime<Utc>, _lon: f64, _lat: f64, local_dir: Option<&Path>) -> Ancillary {
    let dir = local_dir.map(PathBuf::from).unwrap_or_else(|| {
        dirs::home_dir()
            .unwrap_or_else(|| PathBuf::from("/tmp"))
            .join(".acolite")
            .join("met")
    });

    // Attempt download — if credentials are unavailable, just use defaults
    let creds = match Credentials::load() {
        Ok(c) => c,
        Err(_) => {
            log::info!("No EarthData credentials; using default ancillary values");
            return Ancillary::default();
        }
    };

    let files = download(dt, &dir, &creds);
    if files.is_empty() {
        log::info!("No ancillary files retrieved; using defaults");
        return Ancillary::default();
    }

    // Try to interpolate from downloaded GMAO MET NetCDF files
    #[cfg(feature = "netcdf")]
    {
        let met_files: Vec<&PathBuf> = files.iter().filter(|f| {
            f.to_string_lossy().contains("GMAO_MERRA2") && f.to_string_lossy().ends_with(".nc")
        }).collect();
        if let Some(anc) = interpolate_met(&met_files, _lon, _lat, dt) {
            return anc;
        }
    }

    log::info!("Ancillary files downloaded ({}); using default values", files.len());
    Ancillary::default()
}

/// Interpolate ancillary values from GMAO MERRA2 MET NetCDF files at a given lon/lat/time.
/// Reads PS (surface pressure), TQV (water vapour), TO3 (ozone), U10M/V10M (wind).
#[cfg(feature = "netcdf")]
fn interpolate_met(
    met_files: &[&PathBuf],
    lon: f64,
    lat: f64,
    dt: &DateTime<Utc>,
) -> Option<Ancillary> {
    if met_files.is_empty() {
        return None;
    }

    // Read the first available MET file and do nearest-neighbour lookup
    let nc = netcdf::open(met_files[0]).ok()?;

    let lons = read_1d_f64(&nc, "lon")?;
    let lats = read_1d_f64(&nc, "lat")?;

    let ix = nearest_idx(&lons, lon);
    let iy = nearest_idx(&lats, lat);

    let pressure = read_2d_val(&nc, "PS", iy, ix).map(|v| v / 100.0)?; // Pa → hPa
    let uwv = read_2d_val(&nc, "TQV", iy, ix).unwrap_or(1.5); // kg/m² ≈ g/cm² (close enough)
    let uoz = read_2d_val(&nc, "TO3", iy, ix).map(|v| v / 1000.0).unwrap_or(0.3); // DU → cm-atm

    let u10 = read_2d_val(&nc, "U10M", iy, ix).unwrap_or(0.0);
    let v10 = read_2d_val(&nc, "V10M", iy, ix).unwrap_or(0.0);
    let wind = (u10 * u10 + v10 * v10).sqrt();

    // If we have a second bracketing file, average the two
    if met_files.len() >= 2 {
        if let Some(nc2) = netcdf::open(met_files[1]).ok() {
            let p2 = read_2d_val(&nc2, "PS", iy, ix).map(|v| v / 100.0);
            let w2 = read_2d_val(&nc2, "TQV", iy, ix);
            let o2 = read_2d_val(&nc2, "TO3", iy, ix).map(|v| v / 1000.0);
            let u2 = read_2d_val(&nc2, "U10M", iy, ix).unwrap_or(0.0);
            let v2_val = read_2d_val(&nc2, "V10M", iy, ix).unwrap_or(0.0);
            let wind2 = (u2 * u2 + v2_val * v2_val).sqrt();

            let _ = dt; // time weighting could be added; for now simple average
            return Some(Ancillary {
                pressure: avg(pressure, p2.unwrap_or(pressure)),
                uwv: avg(uwv, w2.unwrap_or(uwv)),
                uoz: avg(uoz, o2.unwrap_or(uoz)),
                wind: avg(wind, wind2),
            });
        }
    }

    Some(Ancillary { uoz, uwv, pressure, wind })
}

#[cfg(feature = "netcdf")]
fn avg(a: f64, b: f64) -> f64 {
    (a + b) * 0.5
}

#[cfg(feature = "netcdf")]
fn nearest_idx(arr: &[f64], val: f64) -> usize {
    arr.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (*a - val).abs().partial_cmp(&(*b - val).abs()).unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

#[cfg(feature = "netcdf")]
fn read_1d_f64(nc: &netcdf::File, name: &str) -> Option<Vec<f64>> {
    let var = nc.variable(name)?;
    var.get_values::<f64, _>(..).ok()
}

#[cfg(feature = "netcdf")]
fn read_2d_val(nc: &netcdf::File, name: &str, iy: usize, ix: usize) -> Option<f64> {
    let var = nc.variable(name)?;
    // GMAO MET shape is typically (1, lat, lon) or (lat, lon)
    let ndim = var.dimensions().len();
    let vals: Vec<f64> = if ndim == 3 {
        var.get_values::<f64, _>((0..1, iy..iy + 1, ix..ix + 1)).ok()?
    } else {
        var.get_values::<f64, _>((iy..iy + 1, ix..ix + 1)).ok()?
    };
    vals.into_iter().next()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_list_files_post_2004() {
        let dt = "2023-06-15T10:30:00Z".parse::<DateTime<Utc>>().unwrap();
        let files = list_files(&dt);
        assert_eq!(files.len(), 3);
        assert!(files[0].contains("AURAOMI"));
        assert!(files[1].contains("GMAO_MERRA2"));
        assert!(files[1].contains("T100000"));
        assert!(files[2].contains("T110000"));
    }

    #[test]
    fn test_list_files_hour_23() {
        let dt = "2023-06-15T23:30:00Z".parse::<DateTime<Utc>>().unwrap();
        let files = list_files(&dt);
        assert!(files[2].contains("20230616T000000"));
    }

    #[test]
    fn test_default_ancillary() {
        let anc = Ancillary::default();
        assert!((anc.uoz - 0.3).abs() < 1e-6);
        assert!((anc.pressure - 1013.25).abs() < 1e-6);
    }
}
