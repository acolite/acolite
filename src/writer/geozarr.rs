//! GeoZarr writer for hyperspectral data
//!
//! Writes bands as a 3D Zarr V3 array (band, y, x) with GeoZarr/CF-convention
//! metadata for CRS and spatial transforms. Optimized for hyperspectral sensors
//! (PACE OCI, PRISMA, DESIS, EnMAP, EMIT, etc.) where per-band COG is impractical.

use crate::core::{BandData, Metadata};
use crate::{AcoliteError, Result};
use std::sync::Arc;
use zarrs::array::{ArrayBuilder, data_type};
use zarrs::array::codec::GzipCodec;
use zarrs::group::GroupBuilder;
use zarrs::filesystem::FilesystemStore;

/// Write bands as GeoZarr (Zarr V3 + geospatial metadata)
///
/// Creates a directory-based Zarr V3 store with:
/// - `/data` — 3D array (band, y, x) of f32 reflectance
/// - `/wavelengths` — 1D array of band center wavelengths
/// - `/bandwidths` — 1D array of band widths
/// - Group-level attributes with sensor, datetime, CRS, and geotransform
pub fn write_geozarr(
    output_path: &str,
    bands: &[BandData<f64>],
    metadata: &Metadata,
) -> Result<()> {
    if bands.is_empty() {
        return Err(AcoliteError::Processing("No bands to write".into()));
    }

    let (nrows, ncols) = bands[0].shape();
    let nbands = bands.len();
    let gt = &bands[0].geotransform;

    log::info!("Writing GeoZarr: {} bands, {}×{} → {}", nbands, nrows, ncols, output_path);

    let store_path = std::path::Path::new(output_path);
    if store_path.exists() {
        std::fs::remove_dir_all(store_path).map_err(AcoliteError::Io)?;
    }
    let store = Arc::new(
        FilesystemStore::new(store_path)
            .map_err(|e| AcoliteError::Processing(format!("Zarr store: {}", e)))?
    );

    // Root group with GeoZarr attributes
    let mut attrs = serde_json::Map::new();
    attrs.insert("sensor".into(), metadata.sensor.clone().into());
    attrs.insert("datetime".into(), metadata.datetime.to_rfc3339().into());
    attrs.insert("sun_zenith".into(), metadata.sun_zenith.into());
    attrs.insert("sun_azimuth".into(), metadata.sun_azimuth.into());
    attrs.insert("Conventions".into(), "CF-1.8".into());
    if let Some(ref wkt) = bands[0].projection.wkt {
        attrs.insert("proj:wkt2".into(), wkt.clone().into());
    }
    if let Some(epsg) = bands[0].projection.epsg {
        attrs.insert("proj:epsg".into(), epsg.into());
    }
    attrs.insert("spatial:transform".into(), serde_json::json!([
        gt.pixel_width, gt.x_rotation, gt.x_origin,
        gt.y_rotation, gt.pixel_height, gt.y_origin
    ]));

    let group = GroupBuilder::new()
        .attributes(attrs)
        .build(store.clone(), "/")
        .map_err(|e| AcoliteError::Processing(format!("Zarr group: {}", e)))?;
    group.store_metadata()
        .map_err(|e| AcoliteError::Processing(format!("Zarr group metadata: {}", e)))?;

    let chunk_y = 512.min(nrows as u64);
    let chunk_x = 512.min(ncols as u64);

    // 3D data array (band, y, x)
    let mut data_attrs = serde_json::Map::new();
    data_attrs.insert("long_name".into(), "surface_reflectance".into());
    data_attrs.insert("units".into(), "dimensionless".into());
    data_attrs.insert("_ARRAY_DIMENSIONS".into(), serde_json::json!(["band", "y", "x"]));

    let data_array = ArrayBuilder::new(
        vec![nbands as u64, nrows as u64, ncols as u64],
        vec![1, chunk_y, chunk_x],
        data_type::float32(),
        f32::NAN,
    )
    .bytes_to_bytes_codecs(vec![Arc::new(
        GzipCodec::new(5).map_err(|e| AcoliteError::Processing(format!("Gzip: {}", e)))?,
    )])
    .dimension_names(["band", "y", "x"].into())
    .attributes(data_attrs)
    .build(store.clone(), "/data")
    .map_err(|e| AcoliteError::Processing(format!("Zarr data array: {}", e)))?;
    data_array.store_metadata()
        .map_err(|e| AcoliteError::Processing(format!("Zarr data metadata: {}", e)))?;

    // Write each band
    for (bi, band) in bands.iter().enumerate() {
        let chunk_data: Vec<f32> = band.data.iter().map(|&v| v as f32).collect();
        data_array.store_array_subset(
            &[bi as u64..bi as u64 + 1, 0..nrows as u64, 0..ncols as u64],
            &chunk_data,
        ).map_err(|e| AcoliteError::Processing(format!("Write band {}: {}", bi, e)))?;
    }

    // Wavelengths 1D array
    write_1d_f32(&store, "/wavelengths", "center_wavelength", "nm",
        &bands.iter().map(|b| b.wavelength as f32).collect::<Vec<_>>())?;

    // Bandwidths 1D array
    write_1d_f32(&store, "/bandwidths", "bandwidth", "nm",
        &bands.iter().map(|b| b.bandwidth as f32).collect::<Vec<_>>())?;

    log::info!("GeoZarr written: {}", output_path);
    Ok(())
}

fn write_1d_f32(
    store: &Arc<FilesystemStore>,
    path: &str,
    long_name: &str,
    units: &str,
    data: &[f32],
) -> Result<()> {
    let n = data.len() as u64;
    let mut attrs = serde_json::Map::new();
    attrs.insert("long_name".into(), long_name.into());
    attrs.insert("units".into(), units.into());
    attrs.insert("_ARRAY_DIMENSIONS".into(), serde_json::json!(["band"]));

    let array = ArrayBuilder::new(vec![n], vec![n], data_type::float32(), f32::NAN)
        .dimension_names(["band"].into())
        .attributes(attrs)
        .build(store.clone(), path)
        .map_err(|e| AcoliteError::Processing(format!("Zarr {}: {}", path, e)))?;
    array.store_metadata()
        .map_err(|e| AcoliteError::Processing(format!("Zarr {} metadata: {}", path, e)))?;
    array.store_chunk(&[0], data)
        .map_err(|e| AcoliteError::Processing(format!("Write {}: {}", path, e)))?;
    Ok(())
}
