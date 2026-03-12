//! Cloud Optimized GeoTIFF (COG) writer using GDAL Rust bindings + rayon parallelism

use crate::core::{BandData, Metadata};
use crate::{AcoliteError, Result};
use gdal::raster::Buffer;
use gdal::spatial_ref::SpatialRef;
use gdal::DriverManager;
use rayon::prelude::*;

/// Write bands as parallel per-band Cloud Optimized GeoTIFFs via GDAL.
pub fn write_cog(output_path: &str, bands: &[BandData<f64>], _metadata: &Metadata) -> Result<()> {
    if bands.is_empty() {
        return Err(AcoliteError::Processing("No bands to write".into()));
    }

    let base = output_path.trim_end_matches(".tif");

    let results: Vec<Result<()>> = bands
        .par_iter()
        .map(|band| {
            let band_path = format!("{}_{}.tif", base, band.name);
            write_band_cog(&band_path, band)
        })
        .collect();

    let mut errors = Vec::new();
    for r in results {
        if let Err(e) = r {
            errors.push(e.to_string());
        }
    }
    if !errors.is_empty() {
        return Err(AcoliteError::Processing(format!(
            "COG write errors: {}",
            errors.join("; ")
        )));
    }

    log::info!("Wrote {} per-band COG files (parallel)", bands.len());
    Ok(())
}

/// Write a single band as COG: MEM dataset → create_copy to COG driver.
fn write_band_cog(path: &str, band: &BandData<f64>) -> Result<()> {
    let (height, width) = band.data.dim();

    // Create in-memory dataset
    let mem_driver = DriverManager::get_driver_by_name("MEM")
        .map_err(|e| AcoliteError::Processing(format!("MEM driver: {}", e)))?;
    let mut mem_ds = mem_driver
        .create_with_band_type::<f32, _>("", width, height, 1)
        .map_err(|e| AcoliteError::Processing(format!("MEM create: {}", e)))?;

    // Set geotransform
    let gt = &band.geotransform;
    mem_ds
        .set_geo_transform(&[
            gt.x_origin,
            gt.pixel_width,
            gt.x_rotation,
            gt.y_origin,
            gt.y_rotation,
            gt.pixel_height,
        ])
        .map_err(|e| AcoliteError::Processing(format!("set_geo_transform: {}", e)))?;

    // Set projection
    if let Some(epsg) = band.projection.epsg {
        if let Ok(srs) = SpatialRef::from_epsg(epsg as u32) {
            let _ = mem_ds.set_spatial_ref(&srs);
        }
    } else if let Some(ref wkt) = band.projection.wkt {
        if !wkt.is_empty() {
            if let Ok(srs) = SpatialRef::from_wkt(wkt) {
                let _ = mem_ds.set_spatial_ref(&srs);
            }
        }
    }

    // Write pixel data (f64 → f32)
    let data: Vec<f32> = band.data.iter().map(|&v| v as f32).collect();
    let mut buf = Buffer::new((width, height), data);
    mem_ds
        .rasterband(1)
        .map_err(|e| AcoliteError::Processing(format!("rasterband: {}", e)))?
        .write((0, 0), (width, height), &mut buf)
        .map_err(|e| AcoliteError::Processing(format!("write pixels: {}", e)))?;

    // CreateCopy to COG driver
    let cog_driver = DriverManager::get_driver_by_name("COG")
        .map_err(|e| AcoliteError::Processing(format!("COG driver: {}", e)))?;
    let options =
        gdal::raster::RasterCreationOptions::from_iter(["COMPRESS=DEFLATE", "BLOCKSIZE=512"]);
    mem_ds
        .create_copy(&cog_driver, path, &options)
        .map_err(|e| AcoliteError::Processing(format!("COG create_copy: {}", e)))?;

    log::info!("COG created: {}", path);
    Ok(())
}

/// Check if GDAL COG driver is available.
pub fn cog_available() -> bool {
    DriverManager::get_driver_by_name("COG").is_ok()
}
