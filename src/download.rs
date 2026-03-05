//! Generic sensor data download from NASA DAAC

use crate::{Result, AcoliteError};
use std::path::PathBuf;

/// Sensor-specific DAAC configuration
#[derive(Debug, Clone)]
pub struct DaacConfig {
    pub collection: String,
    pub daac_url: String,
    pub file_pattern: String,
}

impl DaacConfig {
    /// Landsat 8/9 Collection 2 Level 1
    pub fn landsat8() -> Self {
        Self {
            collection: "LANDSAT_OT_C2_L1".to_string(),
            daac_url: "https://data.lpdaac.earthdatacloud.nasa.gov".to_string(),
            file_pattern: "LC08_*.tar".to_string(),
        }
    }
    
    /// Landsat 9 Collection 2 Level 1
    pub fn landsat9() -> Self {
        Self {
            collection: "LANDSAT_OT_C2_L1".to_string(),
            daac_url: "https://data.lpdaac.earthdatacloud.nasa.gov".to_string(),
            file_pattern: "LC09_*.tar".to_string(),
        }
    }
    
    /// Sentinel-2 L1C
    pub fn sentinel2() -> Self {
        Self {
            collection: "SENTINEL-2".to_string(),
            daac_url: "https://datapool.asf.alaska.edu".to_string(),
            file_pattern: "S2*_MSIL1C_*.zip".to_string(),
        }
    }
    
    /// PACE OCI L1B
    pub fn pace() -> Self {
        Self {
            collection: "PACE_OCI_L1B_SCI".to_string(),
            daac_url: "https://obdaac-tea.earthdatacloud.nasa.gov".to_string(),
            file_pattern: "PACE_OCI.*.L1B.*.nc".to_string(),
        }
    }
}

/// Search for sensor data
pub fn search_sensor_data(
    sensor: &str,
    bbox: &[f64; 4],
    start_date: &str,
    end_date: &str,
) -> Result<Vec<String>> {
    let config = match sensor.to_lowercase().as_str() {
        "landsat8" | "landsat-8" => DaacConfig::landsat8(),
        "landsat9" | "landsat-9" => DaacConfig::landsat9(),
        "sentinel2" | "sentinel-2" => DaacConfig::sentinel2(),
        "pace" => DaacConfig::pace(),
        _ => return Err(AcoliteError::Processing(format!("Unknown sensor: {}", sensor))),
    };
    
    crate::earthdata::search_cmr_collection(
        &config.collection,
        bbox,
        start_date,
        end_date,
    )
}

/// Download sensor data file
pub async fn download_sensor_data(
    url: &str,
    output_dir: &str,
    auth: &crate::EarthdataAuth,
) -> Result<PathBuf> {
    let filename = url.split('/').last().unwrap_or("data");
    let output_path = PathBuf::from(output_dir).join(filename);
    
    // Use S3 streaming if it's an S3 URL
    if url.starts_with("s3://") {
        crate::earthdata::download_from_s3_blocking(
            url,
            output_path.to_str().unwrap(),
            auth,
        )?;
    } else {
        // Convert HTTPS to S3 if possible
        let s3_url = if url.contains("earthdatacloud.nasa.gov") {
            url.replace("https://", "s3://")
                .replace(".earthdatacloud.nasa.gov/", "/")
        } else {
            url.to_string()
        };
        
        if s3_url.starts_with("s3://") {
            crate::earthdata::download_from_s3_blocking(
                &s3_url,
                output_path.to_str().unwrap(),
                auth,
            )?;
        } else {
            // Fallback to HTTPS
            crate::earthdata::download_pace_file(
                url,
                output_path.to_str().unwrap(),
                auth,
            )?;
        }
    }
    
    Ok(output_path)
}
