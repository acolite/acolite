//! Zarr I/O for intermediate storage

use crate::{Result, AcoliteError, core::BandData};
use ndarray::Array2;
use std::path::{Path, PathBuf};
use std::fs;
use serde_json;

/// Zarr writer for intermediate band storage
pub struct ZarrWriter {
    base_path: PathBuf,
}

impl ZarrWriter {
    pub fn new<P: AsRef<Path>>(base_path: P) -> Self {
        Self {
            base_path: base_path.as_ref().to_path_buf(),
        }
    }
    
    /// Write band to Zarr array
    pub fn write_band(&self, band: &BandData<f64>) -> Result<()> {
        let band_path = self.base_path.join(&band.name);
        fs::create_dir_all(&band_path)
            .map_err(|e| AcoliteError::Io(e))?;
        
        // Write .zarray metadata
        let zarray = ZarrayMetadata {
            chunks: vec![band.data.nrows(), band.data.ncols()],
            compressor: None,
            dtype: "<f8".to_string(),
            fill_value: 0.0,
            filters: None,
            order: "C".to_string(),
            shape: vec![band.data.nrows(), band.data.ncols()],
            zarr_format: 2,
        };
        
        let zarray_json = serde_json::to_string_pretty(&zarray)
            .map_err(|e| AcoliteError::Processing(format!("JSON error: {}", e)))?;
        
        fs::write(band_path.join(".zarray"), zarray_json)
            .map_err(|e| AcoliteError::Io(e))?;
        
        // Write data chunk (simplified - single chunk)
        let data_bytes: Vec<u8> = band.data.iter()
            .flat_map(|&v| v.to_le_bytes())
            .collect();
        
        fs::write(band_path.join("0.0"), data_bytes)
            .map_err(|e| AcoliteError::Io(e))?;
        
        log::info!("Wrote band {} to Zarr: {:?}", band.name, band_path);
        Ok(())
    }
    
    /// Write multiple bands
    pub fn write_bands(&self, bands: &[BandData<f64>]) -> Result<()> {
        for band in bands {
            self.write_band(band)?;
        }
        Ok(())
    }
}

/// Zarr reader for intermediate storage
pub struct ZarrReader {
    base_path: PathBuf,
}

impl ZarrReader {
    pub fn new<P: AsRef<Path>>(base_path: P) -> Self {
        Self {
            base_path: base_path.as_ref().to_path_buf(),
        }
    }
    
    /// Read band from Zarr array
    pub fn read_band(&self, band_name: &str) -> Result<Array2<f64>> {
        let band_path = self.base_path.join(band_name);
        
        // Read .zarray metadata
        let zarray_json = fs::read_to_string(band_path.join(".zarray"))
            .map_err(|e| AcoliteError::Io(e))?;
        
        let zarray: ZarrayMetadata = serde_json::from_str(&zarray_json)
            .map_err(|e| AcoliteError::Processing(format!("JSON error: {}", e)))?;
        
        // Read data chunk
        let data_bytes = fs::read(band_path.join("0.0"))
            .map_err(|e| AcoliteError::Io(e))?;
        
        // Convert bytes to f64 array
        let values: Vec<f64> = data_bytes.chunks_exact(8)
            .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
            .collect();
        
        let shape = (zarray.shape[0], zarray.shape[1]);
        Array2::from_shape_vec(shape, values)
            .map_err(|e| AcoliteError::Processing(format!("Array shape error: {}", e)))
    }
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct ZarrayMetadata {
    chunks: Vec<usize>,
    compressor: Option<serde_json::Value>,
    dtype: String,
    fill_value: f64,
    filters: Option<serde_json::Value>,
    order: String,
    shape: Vec<usize>,
    zarr_format: u32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Projection, GeoTransform};
    use tempfile::TempDir;
    
    #[test]
    fn test_zarr_write_read() {
        let temp_dir = TempDir::new().unwrap();
        let zarr_path = temp_dir.path().join("test.zarr");
        
        // Write
        let writer = ZarrWriter::new(&zarr_path);
        let proj = Projection::from_epsg(32610);
        let geotrans = GeoTransform::new(0.0, 30.0, 0.0, -30.0);
        
        let data = Array2::from_shape_fn((10, 10), |(i, j)| (i * 10 + j) as f64);
        let band = BandData::new(
            data.clone(),
            443.0,
            16.0,
            "B1".to_string(),
            proj,
            geotrans,
        );
        
        writer.write_band(&band).unwrap();
        
        // Read
        let reader = ZarrReader::new(&zarr_path);
        let read_data = reader.read_band("B1").unwrap();
        
        assert_eq!(read_data.shape(), data.shape());
        assert_eq!(read_data[[0, 0]], data[[0, 0]]);
        assert_eq!(read_data[[9, 9]], data[[9, 9]]);
    }
}
