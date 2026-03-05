//! Cloud Optimized GeoTIFF (COG) writer

use crate::{Result, AcoliteError};
use crate::core::{BandData, Metadata};
use std::process::Command;

/// Write bands as Cloud Optimized GeoTIFF
pub fn write_cog(
    output_path: &str,
    bands: &[BandData<f64>],
    metadata: &Metadata,
) -> Result<()> {
    // First write as regular GeoTIFF
    let temp_path = format!("{}.tmp.tif", output_path);
    crate::io::write_geotiff_multiband(&temp_path, bands, metadata)?;
    
    // Convert to COG using gdal_translate
    log::info!("Converting to COG: {}", output_path);
    
    let output = Command::new("gdal_translate")
        .args(&[
            "-of", "COG",
            "-co", "COMPRESS=DEFLATE",
            "-co", "PREDICTOR=2",
            "-co", "BLOCKSIZE=512",
            &temp_path,
            output_path,
        ])
        .output();
    
    // Clean up temp file
    let _ = std::fs::remove_file(&temp_path);
    
    match output {
        Ok(output) if output.status.success() => {
            log::info!("COG created: {}", output_path);
            Ok(())
        }
        Ok(output) => {
            let stderr = String::from_utf8_lossy(&output.stderr);
            Err(AcoliteError::Processing(format!("gdal_translate failed: {}", stderr)))
        }
        Err(e) => {
            log::warn!("gdal_translate not available: {}", e);
            // Fallback: just rename temp file
            std::fs::rename(&temp_path, output_path)
                .map_err(|e| AcoliteError::Io(e))?;
            log::info!("Created regular GeoTIFF (COG conversion unavailable)");
            Ok(())
        }
    }
}

/// Check if COG tools are available
pub fn cog_available() -> bool {
    Command::new("gdal_translate")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Projection, GeoTransform};
    use ndarray::Array2;
    use chrono::Utc;
    use tempfile::TempDir;
    
    #[test]
    fn test_cog_available() {
        let available = cog_available();
        println!("COG tools available: {}", available);
    }
    
    #[test]
    fn test_write_cog() {
        let temp = TempDir::new().unwrap();
        let output = temp.path().join("test_cog.tif");
        
        let proj = Projection::from_epsg(32610);
        let geotrans = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);
        
        let bands: Vec<_> = (0..3)
            .map(|i| {
                BandData::new(
                    Array2::from_elem((100, 100), i as f64),
                    500.0 + i as f64 * 100.0,
                    50.0,
                    format!("B{}", i + 1),
                    proj.clone(),
                    geotrans.clone(),
                )
            })
            .collect();
        
        let metadata = Metadata::new("TEST".to_string(), Utc::now());
        
        let result = write_cog(output.to_str().unwrap(), &bands, &metadata);
        assert!(result.is_ok());
        assert!(output.exists());
    }
}
