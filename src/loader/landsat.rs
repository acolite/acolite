//! Landsat scene loader — reads L1 GeoTIFF bands + MTL metadata

use crate::core::{BandData, GeoTransform, Metadata};
use crate::loader::geotiff::read_geotiff_band;
use crate::sensors::mtl::parse_mtl;
use crate::{AcoliteError, Result};
use std::path::{Path, PathBuf};

/// Landsat OLI band definitions: (band_num, wavelength_nm, bandwidth_nm)
const LANDSAT_BANDS: [(u8, f64, f64); 7] = [
    (1, 443.0, 16.0),   // Coastal
    (2, 482.0, 60.0),   // Blue
    (3, 561.0, 57.0),   // Green
    (4, 655.0, 37.0),   // Red
    (5, 865.0, 28.0),   // NIR
    (6, 1609.0, 85.0),  // SWIR1
    (7, 2201.0, 187.0), // SWIR2
];

/// Reflectance calibration coefficients parsed from MTL
#[derive(Debug, Clone)]
pub struct ReflectanceCoeffs {
    /// Per-band (band_num, mult, add)
    pub bands: Vec<(u8, f64, f64)>,
}

impl Default for ReflectanceCoeffs {
    fn default() -> Self {
        // Collection 2 defaults
        Self {
            bands: (1..=7).map(|b| (b, 2.0e-5, -0.1)).collect(),
        }
    }
}

/// Load a complete Landsat scene: bands + metadata from MTL.
/// If `limit` is provided as `[south, west, north, east]` in geographic coords,
/// bands are spatially subset after loading.
pub fn load_landsat_scene(scene_dir: &Path) -> Result<(Vec<BandData<u16>>, Metadata)> {
    load_landsat_scene_limit(scene_dir, None)
}

/// Load Landsat scene with optional geographic limit `[south, west, north, east]`.
///
/// When a limit is provided and the `gdal-support` feature is enabled, uses
/// GDAL windowed reads to only load the pixels within the limit from disk.
/// This dramatically reduces memory usage for subsetted processing.
pub fn load_landsat_scene_limit(
    scene_dir: &Path,
    limit: Option<&[f64; 4]>,
) -> Result<(Vec<BandData<u16>>, Metadata)> {
    let metadata = if let Some(mtl_path) = find_mtl(scene_dir) {
        parse_mtl(&mtl_path).unwrap_or_else(|e| {
            log::warn!("MTL parse failed ({}), using fallback metadata", e);
            fallback_metadata(scene_dir)
        })
    } else {
        log::warn!("No MTL file found, using fallback metadata");
        fallback_metadata(scene_dir)
    };

    let band_numbers: Vec<u8> = LANDSAT_BANDS.iter().map(|b| b.0).collect();

    // If limit is provided, compute the pixel window FIRST, then read only that window
    let bands = if let Some(lim) = limit {
        load_landsat_bands_windowed(scene_dir, &band_numbers, lim)?
    } else {
        load_landsat_bands(scene_dir, &band_numbers)?
    };

    Ok((bands, metadata))
}

/// Load Landsat bands with a geographic limit applied as a windowed read.
///
/// Reads only the geotransform from the first band to compute the pixel window,
/// then reads only the required pixels from each band file.
fn load_landsat_bands_windowed(
    scene_dir: &Path,
    band_numbers: &[u8],
    limit: &[f64; 4],
) -> Result<Vec<BandData<u16>>> {
    // First, determine the pixel window from the first band's metadata
    let first_band_path = find_band_file(scene_dir, band_numbers[0])?;
    let (rows, cols, gt, proj, wkt_str) = read_geotiff_metadata(&first_band_path)?;

    let wkt_ref = if wkt_str.is_empty() { None } else { Some(wkt_str.as_str()) };

    let window = crate::loader::latlon_limit_to_pixel_subset(
        gt.x_origin, gt.pixel_width, gt.y_origin, gt.pixel_height,
        rows, cols, limit, wkt_ref,
    );

    match window {
        Some((r0, c0, nr, nc)) => {
            log::info!(
                "Windowed read: {}×{} pixels (from {}×{} full scene, offset [{}, {}])",
                nr, nc, rows, cols, r0, c0
            );
            // Read only the window from each band
            band_numbers
                .iter()
                .map(|&num| {
                    let path = find_band_file(scene_dir, num)?;
                    log::info!("Reading band {} (window {}×{}): {:?}", num, nr, nc, path);
                    let mut band = read_geotiff_band_with_window(&path, (c0, r0, nc, nr))?;
                    if let Some(&(_, wl, bw)) = LANDSAT_BANDS.iter().find(|(n, _, _)| *n == num) {
                        band.wavelength = wl;
                        band.bandwidth = bw;
                    }
                    band.name = format!("B{}", num);
                    Ok(band)
                })
                .collect()
        }
        None => {
            log::warn!("Limit {:?} does not intersect scene, loading full extent", limit);
            load_landsat_bands(scene_dir, band_numbers)
        }
    }
}

/// Read only the geotransform, dimensions, and projection from a GeoTIFF without loading raster data.
#[cfg(feature = "gdal-support")]
fn read_geotiff_metadata(path: &Path) -> Result<(usize, usize, GeoTransform, crate::core::Projection, String)> {
    use gdal::Dataset;
    let ds = Dataset::open(path)
        .map_err(|e| AcoliteError::Gdal(format!("Cannot open {:?}: {}", path, e)))?;
    let (w, h) = ds.raster_size();
    let gt_raw = ds.geo_transform()
        .map_err(|e| AcoliteError::Gdal(format!("GeoTransform: {}", e)))?;
    let wkt = ds.projection();
    let proj = crate::core::Projection::from_wkt(wkt.clone());
    let gt = GeoTransform::new(gt_raw[0], gt_raw[1], gt_raw[3], gt_raw[5]);
    Ok((h as usize, w as usize, gt, proj, wkt))
}

#[cfg(not(feature = "gdal-support"))]
fn read_geotiff_metadata(path: &Path) -> Result<(usize, usize, GeoTransform, crate::core::Projection, String)> {
    // Without GDAL, load the full band to get metadata (fallback)
    let band = crate::loader::geotiff::read_geotiff_band(path)?;
    let (rows, cols) = band.data.dim();
    Ok((rows, cols, band.geotransform, band.projection, String::new()))
}

/// Read a windowed subset of a GeoTIFF band.
/// Window is (x_offset, y_offset, x_size, y_size) in pixel coordinates.
#[cfg(feature = "gdal-support")]
fn read_geotiff_band_with_window(path: &Path, window: (usize, usize, usize, usize)) -> Result<BandData<u16>> {
    crate::loader::geotiff::read_geotiff_band_window(&path.to_string_lossy(), window)
}

#[cfg(not(feature = "gdal-support"))]
fn read_geotiff_band_with_window(path: &Path, window: (usize, usize, usize, usize)) -> Result<BandData<u16>> {
    // Without GDAL, load full band then crop (fallback — still uses more memory)
    let full_band = crate::loader::geotiff::read_geotiff_band(path)?;
    let (x_off, y_off, x_size, y_size) = window;
    use ndarray::s;
    let data = full_band.data.slice(s![y_off..y_off+y_size, x_off..x_off+x_size]).to_owned();
    let gt = GeoTransform::new(
        full_band.geotransform.x_origin + x_off as f64 * full_band.geotransform.pixel_width,
        full_band.geotransform.pixel_width,
        full_band.geotransform.y_origin + y_off as f64 * full_band.geotransform.pixel_height,
        full_band.geotransform.pixel_height,
    );
    Ok(BandData::new(data, 0.0, 0.0, String::new(), full_band.projection, gt))
}

/// Load specific Landsat bands by number
pub fn load_landsat_bands(scene_dir: &Path, band_numbers: &[u8]) -> Result<Vec<BandData<u16>>> {
    band_numbers
        .iter()
        .map(|&num| {
            let path = find_band_file(scene_dir, num)?;
            log::info!("Reading band {}: {:?}", num, path);

            let mut band = read_geotiff_band(&path)?;

            if let Some(&(_, wl, bw)) = LANDSAT_BANDS.iter().find(|(n, _, _)| *n == num) {
                band.wavelength = wl;
                band.bandwidth = bw;
            }
            band.name = format!("B{}", num);

            Ok(band)
        })
        .collect()
}

/// Parse reflectance calibration coefficients from MTL
pub fn parse_reflectance_coeffs(scene_dir: &Path) -> ReflectanceCoeffs {
    let mtl_path = match find_mtl(scene_dir) {
        Some(p) => p,
        None => return ReflectanceCoeffs::default(),
    };
    let content = match std::fs::read_to_string(&mtl_path) {
        Ok(c) => c,
        Err(_) => return ReflectanceCoeffs::default(),
    };

    let mut bands = Vec::new();
    for b in 1u8..=7 {
        let mult_key = format!("REFLECTANCE_MULT_BAND_{}", b);
        let add_key = format!("REFLECTANCE_ADD_BAND_{}", b);
        let mult = extract_mtl_value(&content, &mult_key).unwrap_or(2.0e-5);
        let add = extract_mtl_value(&content, &add_key).unwrap_or(-0.1);
        bands.push((b, mult, add));
    }
    ReflectanceCoeffs { bands }
}

/// Find the *_MTL.txt file in a scene directory
pub fn find_mtl(scene_dir: &Path) -> Option<String> {
    std::fs::read_dir(scene_dir).ok()?.find_map(|e| {
        let p = e.ok()?.path();
        let name = p.file_name()?.to_string_lossy().to_string();
        if name.ends_with("_MTL.txt") {
            Some(p.to_string_lossy().to_string())
        } else {
            None
        }
    })
}

fn extract_mtl_value(content: &str, key: &str) -> Option<f64> {
    content.lines().find_map(|line| {
        let line = line.trim();
        if line.starts_with(key) {
            line.split('=').nth(1)?.trim().parse().ok()
        } else {
            None
        }
    })
}

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

fn detect_spacecraft(dir: &Path) -> &str {
    let name = dir.file_name().and_then(|n| n.to_str()).unwrap_or("");
    if name.starts_with("LC09") {
        "LANDSAT_9_OLI"
    } else if name.starts_with("LC08") {
        "LANDSAT_8_OLI"
    } else {
        "LANDSAT_OLI"
    }
}

fn fallback_metadata(scene_dir: &Path) -> Metadata {
    Metadata::new(detect_spacecraft(scene_dir).to_string(), chrono::Utc::now())
}
