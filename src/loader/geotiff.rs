//! Generic GeoTIFF reader using GDAL

use crate::core::{BandData, GeoTransform, Projection};
use crate::{AcoliteError, Result};
use ndarray::Array2;
use std::path::Path;

/// Read a single u16 band from GeoTIFF
#[cfg(feature = "gdal-support")]
pub fn read_geotiff_band(path: &Path) -> Result<BandData<u16>> {
    use gdal::Dataset;
    let ds = Dataset::open(path).map_err(|e| AcoliteError::Gdal(format!("Open: {}", e)))?;
    let rb = ds
        .rasterband(1)
        .map_err(|e| AcoliteError::Gdal(format!("Band: {}", e)))?;
    let (w, h) = ds.raster_size();
    let gt = ds
        .geo_transform()
        .map_err(|e| AcoliteError::Gdal(format!("GT: {}", e)))?;

    let mut data = Array2::zeros((h as usize, w as usize));
    rb.read_into_slice(
        (0, 0),
        (w, h),
        (w, h),
        data.as_slice_mut()
            .ok_or_else(|| AcoliteError::Processing("Non-contiguous array".into()))?,
        None,
    )
    .map_err(|e| AcoliteError::Gdal(format!("Read: {}", e)))?;

    Ok(BandData::new(
        data,
        0.0,
        0.0,
        path.file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string(),
        Projection::from_wkt(ds.projection()),
        GeoTransform::new(gt[0], gt[1], gt[3], gt[5]),
    ))
}

/// Stub when GDAL is not available — reads a zero-filled band
#[cfg(not(feature = "gdal-support"))]
pub fn read_geotiff_band(path: &Path) -> Result<BandData<u16>> {
    // Use tiff crate as fallback for basic reading
    use std::fs::File;
    use std::io::BufReader;
    use tiff::decoder::Decoder;

    let file = File::open(path).map_err(AcoliteError::Io)?;
    let mut decoder = Decoder::new(BufReader::new(file))
        .map_err(|e| AcoliteError::Processing(format!("TIFF decode: {}", e)))?;

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
        path.file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string(),
        Projection::from_epsg(32610),
        GeoTransform::new(0.0, 30.0, 0.0, -30.0),
    ))
}
