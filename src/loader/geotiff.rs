//! Generic GeoTIFF reader using GDAL

use crate::{Result, AcoliteError};
use crate::core::{BandData, Projection, GeoTransform};
use gdal::Dataset;
use ndarray::Array2;
use std::path::Path;

/// Read a single u16 band from GeoTIFF
pub fn read_geotiff_band(path: &Path) -> Result<BandData<u16>> {
    let ds = Dataset::open(path)
        .map_err(|e| AcoliteError::Gdal(format!("Open: {}", e)))?;
    let rb = ds.rasterband(1)
        .map_err(|e| AcoliteError::Gdal(format!("Band: {}", e)))?;
    let (w, h) = ds.raster_size();
    let gt = ds.geo_transform()
        .map_err(|e| AcoliteError::Gdal(format!("GT: {}", e)))?;

    let mut data = Array2::zeros((h as usize, w as usize));
    rb.read_into_slice((0, 0), (w, h), (w, h),
        data.as_slice_mut()
            .ok_or_else(|| AcoliteError::Processing("Non-contiguous array".into()))?,
        None)
        .map_err(|e| AcoliteError::Gdal(format!("Read: {}", e)))?;

    Ok(BandData::new(
        data, 0.0, 0.0,
        path.file_stem().unwrap_or_default().to_string_lossy().to_string(),
        Projection::from_wkt(ds.projection()),
        GeoTransform::new(gt[0], gt[1], gt[3], gt[5]),
    ))
}
