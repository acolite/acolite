//! Generic GeoTIFF reader using GDAL (supports BigTIFF, COG, /vsicurl/)
//!
//! When the `gdal-support` feature is enabled, uses GDAL for reading which
//! supports BigTIFF, Cloud Optimized GeoTIFF (COG), and remote reading via
//! /vsicurl/. Falls back to the pure-Rust `tiff` crate when GDAL is unavailable.

use crate::core::{BandData, GeoTransform, Projection};
use crate::{AcoliteError, Result};
use ndarray::Array2;
use std::path::Path;

/// Read a single u16 band from a GeoTIFF file (local or remote via /vsicurl/).
///
/// Supports:
/// - Classic TIFF and BigTIFF formats
/// - Cloud Optimized GeoTIFF (COG) with internal tiling
/// - Remote reading via GDAL /vsicurl/ (pass URL as path)
/// - Spatial subsetting via `window` parameter
#[cfg(feature = "gdal-support")]
pub fn read_geotiff_band(path: &Path) -> Result<BandData<u16>> {
    read_geotiff_band_gdal(&path.to_string_lossy())
}

/// Read a GeoTIFF band from any GDAL-compatible path string.
/// Accepts local paths, /vsicurl/ URLs, or /vsis3/ paths.
#[cfg(feature = "gdal-support")]
pub fn read_geotiff_band_gdal(path: &str) -> Result<BandData<u16>> {
    use gdal::Dataset;

    let ds = Dataset::open(path)
        .map_err(|e| AcoliteError::Gdal(format!("Cannot open '{}': {}", path, e)))?;
    let rb = ds
        .rasterband(1)
        .map_err(|e| AcoliteError::Gdal(format!("Cannot read band 1: {}", e)))?;
    let (w, h) = ds.raster_size();
    let gt = ds
        .geo_transform()
        .map_err(|e| AcoliteError::Gdal(format!("Cannot read geotransform: {}", e)))?;

    // Read the full band as u16
    let buf = rb
        .read_as::<u16>((0, 0), (w, h), (w, h), None)
        .map_err(|e| AcoliteError::Gdal(format!("Cannot read raster data: {}", e)))?;

    let data = Array2::from_shape_vec((h as usize, w as usize), buf.data().to_vec())
        .map_err(|e| AcoliteError::Processing(format!("Array reshape: {}", e)))?;

    let projection = Projection::from_wkt(ds.projection());

    Ok(BandData::new(
        data,
        0.0,
        0.0,
        extract_name(path),
        projection,
        GeoTransform::new(gt[0], gt[1], gt[3], gt[5]),
    ))
}

/// Read a spatial subset of a GeoTIFF band (windowed read).
///
/// `window` is (x_offset, y_offset, x_size, y_size) in pixel coordinates.
/// This is efficient for COGs as only the needed tiles are fetched.
#[cfg(feature = "gdal-support")]
pub fn read_geotiff_band_window(
    path: &str,
    window: (usize, usize, usize, usize),
) -> Result<BandData<u16>> {
    use gdal::Dataset;

    let ds = Dataset::open(path)
        .map_err(|e| AcoliteError::Gdal(format!("Cannot open '{}': {}", path, e)))?;
    let rb = ds
        .rasterband(1)
        .map_err(|e| AcoliteError::Gdal(format!("Cannot read band 1: {}", e)))?;
    let gt = ds
        .geo_transform()
        .map_err(|e| AcoliteError::Gdal(format!("Cannot read geotransform: {}", e)))?;

    let (x_off, y_off, x_size, y_size) = window;

    let buf = rb
        .read_as::<u16>(
            (x_off as isize, y_off as isize),
            (x_size, y_size),
            (x_size, y_size),
            None,
        )
        .map_err(|e| AcoliteError::Gdal(format!("Cannot read window: {}", e)))?;

    let data = Array2::from_shape_vec((y_size, x_size), buf.data().to_vec())
        .map_err(|e| AcoliteError::Processing(format!("Array reshape: {}", e)))?;

    // Adjust geotransform for the window offset
    let new_gt = GeoTransform::new(
        gt[0] + x_off as f64 * gt[1],
        gt[1],
        gt[3] + y_off as f64 * gt[5],
        gt[5],
    );

    Ok(BandData::new(
        data,
        0.0,
        0.0,
        extract_name(path),
        Projection::from_wkt(ds.projection()),
        new_gt,
    ))
}

/// Read a GeoTIFF band directly from a URL using GDAL /vsicurl/.
///
/// This enables direct reading of COG files from USGS LandsatLook, AWS S3, etc.
/// without downloading the entire file first.
///
/// Example:
/// ```no_run
/// # use acolite_rs::loader::geotiff::read_geotiff_band_url;
/// let band = read_geotiff_band_url(
///     "https://landsatlook.usgs.gov/data/.../LC09_..._B1.TIF",
///     None, // full extent
/// );
/// ```
#[cfg(feature = "gdal-support")]
pub fn read_geotiff_band_url(url: &str, window: Option<(usize, usize, usize, usize)>) -> Result<BandData<u16>> {
    let vsicurl_path = if url.starts_with("/vsicurl/") || url.starts_with("/vsis3/") {
        url.to_string()
    } else if url.starts_with("http") {
        format!("/vsicurl/{}", url)
    } else {
        url.to_string()
    };

    if let Some(w) = window {
        read_geotiff_band_window(&vsicurl_path, w)
    } else {
        read_geotiff_band_gdal(&vsicurl_path)
    }
}

/// Stub when GDAL is not available — uses pure-Rust tiff crate.
/// NOTE: This does NOT support BigTIFF (files > 4 GB) — Landsat Collection 2
/// COGs are BigTIFF and require the `gdal-support` feature.
#[cfg(not(feature = "gdal-support"))]
pub fn read_geotiff_band(path: &Path) -> Result<BandData<u16>> {
    use std::fs::File;
    use std::io::BufReader;
    use tiff::decoder::Decoder;

    let file = File::open(path).map_err(AcoliteError::Io)?;
    let mut decoder = Decoder::new(BufReader::new(file))
        .map_err(|e| AcoliteError::Processing(format!(
            "TIFF decode error for {:?}: {}. Landsat Collection 2 uses BigTIFF which requires the 'gdal-support' feature.",
            path.file_name().unwrap_or_default(), e
        )))?;

    let (w, h) = decoder
        .dimensions()
        .map_err(|e| AcoliteError::Processing(format!("TIFF dims: {}", e)))?;

    let image = decoder
        .read_image()
        .map_err(|e| AcoliteError::Processing(format!("TIFF read: {}", e)))?;

    let pixels: Vec<u16> = match image {
        tiff::decoder::DecodingResult::U16(v) => v,
        tiff::decoder::DecodingResult::U8(v) => v.into_iter().map(|x| x as u16).collect(),
        _ => {
            return Err(AcoliteError::Processing(
                "Unsupported TIFF pixel type".into(),
            ))
        }
    };

    let data = Array2::from_shape_vec((h as usize, w as usize), pixels)
        .map_err(|e| AcoliteError::Processing(format!("Array shape: {}", e)))?;

    Ok(BandData::new(
        data,
        0.0,
        0.0,
        extract_name(&path.to_string_lossy()),
        Projection::from_epsg(32610),
        GeoTransform::new(0.0, 30.0, 0.0, -30.0),
    ))
}

/// Extract a short name from a path or URL
fn extract_name(path: &str) -> String {
    path.rsplit('/')
        .next()
        .unwrap_or(path)
        .trim_end_matches(".TIF")
        .trim_end_matches(".tif")
        .to_string()
}
