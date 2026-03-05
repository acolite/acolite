//! GeoTIFF I/O using tiff crate

use crate::{Result, AcoliteError};
use crate::core::{BandData, Projection, GeoTransform, Metadata};
use ndarray::Array2;
use std::fs::File;
use std::io::BufWriter;
use tiff::decoder::{Decoder, DecodingResult};
use tiff::encoder::{TiffEncoder, colortype};

/// Read a TIFF band (basic implementation)
pub fn read_geotiff_band(path: &str, _band_index: isize) -> Result<BandData<f64>> {
    let file = File::open(path)
        .map_err(|e| AcoliteError::Io(e))?;
    
    let mut decoder = Decoder::new(file)
        .map_err(|e| AcoliteError::Processing(format!("TIFF decode error: {}", e)))?;
    
    let (width, height) = decoder.dimensions()
        .map_err(|e| AcoliteError::Processing(format!("Get dimensions failed: {}", e)))?;
    
    let image = decoder.read_image()
        .map_err(|e| AcoliteError::Processing(format!("Read image failed: {}", e)))?;
    
    let data = match image {
        DecodingResult::U8(buf) => {
            Array2::from_shape_vec(
                (height as usize, width as usize),
                buf.into_iter().map(|v| v as f64).collect()
            ).map_err(|e| AcoliteError::Processing(format!("Shape error: {}", e)))?
        },
        DecodingResult::U16(buf) => {
            Array2::from_shape_vec(
                (height as usize, width as usize),
                buf.into_iter().map(|v| v as f64).collect()
            ).map_err(|e| AcoliteError::Processing(format!("Shape error: {}", e)))?
        },
        DecodingResult::F32(buf) => {
            Array2::from_shape_vec(
                (height as usize, width as usize),
                buf.into_iter().map(|v| v as f64).collect()
            ).map_err(|e| AcoliteError::Processing(format!("Shape error: {}", e)))?
        },
        _ => return Err(AcoliteError::Processing("Unsupported TIFF format".to_string())),
    };
    
    // Default geotransform and projection
    let geotrans = GeoTransform::new(0.0, 1.0, 0.0, -1.0);
    let proj = Projection::from_epsg(4326);
    
    Ok(BandData::new(
        data,
        0.0,
        0.0,
        "Band_1".to_string(),
        proj,
        geotrans,
    ))
}

/// Write a GeoTIFF band (basic implementation)
pub fn write_geotiff_band(
    path: &str,
    band: &BandData<f64>,
    _metadata: &Metadata,
) -> Result<()> {
    let file = File::create(path)
        .map_err(|e| AcoliteError::Io(e))?;
    
    let mut encoder = TiffEncoder::new(BufWriter::new(file))
        .map_err(|e| AcoliteError::Processing(format!("TIFF encode error: {}", e)))?;
    
    let (height, width) = band.data.dim();
    
    // Convert to f32 for TIFF
    let data_vec: Vec<f32> = band.data.iter().map(|&v| v as f32).collect();
    
    encoder.write_image::<colortype::Gray32Float>(width as u32, height as u32, &data_vec)
        .map_err(|e| AcoliteError::Processing(format!("Write failed: {}", e)))?;
    
    Ok(())
}

/// Write multiple bands (writes first band only for now)
pub fn write_geotiff_multiband(
    path: &str,
    bands: &[BandData<f64>],
    metadata: &Metadata,
) -> Result<()> {
    if bands.is_empty() {
        return Err(AcoliteError::Processing("No bands to write".to_string()));
    }
    
    // For now, just write the first band
    write_geotiff_band(path, &bands[0], metadata)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;
    use chrono::Utc;
    
    #[test]
    fn test_geotiff_write_read() {
        let temp = TempDir::new().unwrap();
        let path = temp.path().join("test.tif");
        
        let proj = Projection::from_epsg(32610);
        let geotrans = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);
        let data = Array2::from_shape_fn((100, 100), |(i, j)| (i + j) as f64);
        
        let band = BandData::new(
            data.clone(),
            665.0,
            30.0,
            "Red".to_string(),
            proj,
            geotrans,
        );
        
        let metadata = Metadata::new("TEST".to_string(), Utc::now());
        
        // Write
        write_geotiff_band(path.to_str().unwrap(), &band, &metadata).unwrap();
        
        // Read
        let read_band = read_geotiff_band(path.to_str().unwrap(), 1).unwrap();
        
        assert_eq!(read_band.data.dim(), (100, 100));
        assert_eq!(read_band.data[[0, 0]], 0.0);
        assert_eq!(read_band.data[[50, 50]], 100.0);
    }
}
