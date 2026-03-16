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

    // For now return defaults — full HDF/NetCDF interpolation requires pyhdf-equivalent
    // which is out of scope. The download infrastructure is in place for when
    // netcdf-based interpolation is added.
    log::info!(
        "Ancillary files downloaded ({}); using default interpolation values",
        files.len()
    );
    Ancillary::default()
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
