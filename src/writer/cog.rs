//! Cloud Optimized GeoTIFF (COG) writer

use crate::{Result, AcoliteError};
use crate::core::{BandData, Metadata};
use std::process::Command;

/// Write bands as Cloud Optimized GeoTIFF
pub fn write_cog(
    output_path: &str,
    bands: &[BandData<f64>],
    _metadata: &Metadata,
) -> Result<()> {
    if bands.is_empty() {
        return Err(AcoliteError::Processing("No bands to write".into()));
    }

    let temp_path = format!("{}.tmp.tif", output_path);
    write_temp_geotiff(&temp_path, bands)?;

    log::info!("Converting to COG: {}", output_path);
    let result = Command::new("gdal_translate")
        .args(["-of", "COG", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=2",
               "-co", "BLOCKSIZE=512", &temp_path, output_path])
        .output();

    let _ = std::fs::remove_file(&temp_path);

    match result {
        Ok(out) if out.status.success() => {
            log::info!("COG created: {}", output_path);
            Ok(())
        }
        _ => {
            // Fallback: write was already done as temp, but it's deleted.
            // Re-write as plain tiff at output_path
            write_temp_geotiff(output_path, bands)?;
            log::warn!("COG unavailable, wrote regular GeoTIFF");
            Ok(())
        }
    }
}

fn write_temp_geotiff(path: &str, bands: &[BandData<f64>]) -> Result<()> {
    use std::fs::File;
    use std::io::BufWriter;
    use tiff::encoder::{TiffEncoder, colortype};

    let (height, width) = bands[0].data.dim();
    let file = File::create(path).map_err(AcoliteError::Io)?;
    let mut enc = TiffEncoder::new(BufWriter::new(file))
        .map_err(|e| AcoliteError::Processing(format!("TIFF: {}", e)))?;

    let data: Vec<f32> = bands[0].data.iter().map(|&v| v as f32).collect();
    enc.write_image::<colortype::Gray32Float>(width as u32, height as u32, &data)
        .map_err(|e| AcoliteError::Processing(format!("TIFF write: {}", e)))?;
    Ok(())
}

/// Check if COG tools are available
pub fn cog_available() -> bool {
    Command::new("gdal_translate").arg("--version").output()
        .map(|o| o.status.success()).unwrap_or(false)
}
