//! Landsat scene loader — reads L1 GeoTIFF bands + MTL metadata

use crate::{Result, AcoliteError};
use crate::core::{BandData, Metadata};
use crate::loader::geotiff::read_geotiff_band;
use std::path::{Path, PathBuf};

/// Landsat OLI band definitions
const LANDSAT_BANDS: [(u8, f64, f64); 7] = [
    (1, 443.0, 16.0),   // Coastal
    (2, 482.0, 60.0),   // Blue
    (3, 561.0, 57.0),   // Green
    (4, 655.0, 37.0),   // Red
    (5, 865.0, 28.0),   // NIR
    (6, 1609.0, 85.0),  // SWIR1
    (7, 2201.0, 187.0), // SWIR2
];

/// Load a complete Landsat scene from a directory
pub fn load_landsat_scene(scene_dir: &Path) -> Result<(Vec<BandData<u16>>, Metadata)> {
    let band_numbers: Vec<u8> = LANDSAT_BANDS.iter().map(|b| b.0).collect();
    let bands = load_landsat_bands(scene_dir, &band_numbers)?;

    // TODO: Parse MTL for real metadata (sun angles, datetime, etc.)
    let metadata = Metadata::new(
        detect_spacecraft(scene_dir).to_string(),
        chrono::Utc::now(),
    );

    Ok((bands, metadata))
}

/// Load specific Landsat bands by number
pub fn load_landsat_bands(scene_dir: &Path, band_numbers: &[u8]) -> Result<Vec<BandData<u16>>> {
    band_numbers
        .iter()
        .map(|&num| {
            let path = find_band_file(scene_dir, num)?;
            log::info!("Reading band {}: {:?}", num, path);

            let mut band = read_geotiff_band(&path)?;

            // Assign wavelength from band number
            if let Some(&(_, wl, bw)) = LANDSAT_BANDS.iter().find(|(n, _, _)| *n == num) {
                band.wavelength = wl;
                band.bandwidth = bw;
            }
            band.name = format!("B{}", num);

            Ok(band)
        })
        .collect()
}

/// Find band file in scene directory
fn find_band_file(dir: &Path, band_num: u8) -> Result<PathBuf> {
    let pattern = format!("_B{}", band_num);

    for entry in std::fs::read_dir(dir).map_err(AcoliteError::Io)? {
        let path = entry.map_err(AcoliteError::Io)?.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if name.contains(&pattern) && (name.ends_with(".TIF") || name.ends_with(".tif")) {
                return Ok(path);
            }
        }
    }

    Err(AcoliteError::Processing(format!(
        "Band {} not found in {:?}",
        band_num, dir
    )))
}

/// Detect spacecraft from directory name
fn detect_spacecraft(dir: &Path) -> &str {
    let name = dir
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("");

    if name.starts_with("LC09") {
        "LANDSAT_9_OLI"
    } else if name.starts_with("LC08") {
        "LANDSAT_8_OLI"
    } else {
        "LANDSAT_OLI"
    }
}
