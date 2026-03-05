//! GeoZarr writer for hyperspectral data

use crate::{Result, AcoliteError};
use crate::core::{BandData, Metadata};
use std::fs;
use std::path::Path;

/// GeoZarr metadata following the spec
#[derive(serde::Serialize)]
struct GeoZarrMetadata {
    multiscales: Vec<Multiscale>,
}

#[derive(serde::Serialize)]
struct Multiscale {
    version: String,
    name: String,
    datasets: Vec<Dataset>,
    axes: Vec<Axis>,
}

#[derive(serde::Serialize)]
struct Dataset {
    path: String,
}

#[derive(serde::Serialize)]
struct Axis {
    name: String,
    #[serde(rename = "type")]
    axis_type: String,
    unit: Option<String>,
}

/// Zarray metadata for each array
#[derive(serde::Serialize)]
struct ZarrayMetadata {
    zarr_format: u8,
    shape: Vec<usize>,
    chunks: Vec<usize>,
    dtype: String,
    compressor: Option<serde_json::Value>,
    fill_value: f64,
    order: String,
    filters: Option<Vec<serde_json::Value>>,
}

/// Write hyperspectral bands as GeoZarr
pub fn write_geozarr(
    output_path: &str,
    bands: &[BandData<f64>],
    metadata: &Metadata,
) -> Result<()> {
    let root = Path::new(output_path);
    fs::create_dir_all(root)
        .map_err(|e| AcoliteError::Io(e))?;
    
    if bands.is_empty() {
        return Err(AcoliteError::Processing("No bands to write".to_string()));
    }
    
    let (height, width) = bands[0].data.dim();
    let nbands = bands.len();
    
    log::info!("Writing GeoZarr: {} bands ({}×{})", nbands, height, width);
    
    // Create data array directory
    let data_dir = root.join("data");
    fs::create_dir_all(&data_dir)
        .map_err(|e| AcoliteError::Io(e))?;
    
    // Write .zarray metadata
    let zarray = ZarrayMetadata {
        zarr_format: 2,
        shape: vec![nbands, height, width],
        chunks: vec![1, height.min(512), width.min(512)],
        dtype: "<f8".to_string(),
        compressor: Some(serde_json::json!({
            "id": "zlib",
            "level": 5
        })),
        fill_value: f64::NAN,
        order: "C".to_string(),
        filters: None,
    };
    
    let zarray_path = data_dir.join(".zarray");
    let zarray_json = serde_json::to_string_pretty(&zarray)
        .map_err(|e| AcoliteError::Processing(format!("JSON error: {}", e)))?;
    fs::write(zarray_path, zarray_json)
        .map_err(|e| AcoliteError::Io(e))?;
    
    // Write each band as a chunk
    for (i, band) in bands.iter().enumerate() {
        let chunk_name = format!("{}.0.0", i);
        let chunk_path = data_dir.join(&chunk_name);
        
        // Write raw f64 data (little-endian)
        let data_vec: Vec<u8> = band.data.iter()
            .flat_map(|&v| v.to_le_bytes())
            .collect();
        
        fs::write(&chunk_path, &data_vec)
            .map_err(|e| AcoliteError::Io(e))?;
    }
    
    // Write GeoZarr metadata
    let geozarr_meta = GeoZarrMetadata {
        multiscales: vec![Multiscale {
            version: "0.1".to_string(),
            name: metadata.sensor.clone(),
            datasets: vec![Dataset {
                path: "data".to_string(),
            }],
            axes: vec![
                Axis {
                    name: "band".to_string(),
                    axis_type: "spectral".to_string(),
                    unit: Some("nm".to_string()),
                },
                Axis {
                    name: "y".to_string(),
                    axis_type: "space".to_string(),
                    unit: Some("meter".to_string()),
                },
                Axis {
                    name: "x".to_string(),
                    axis_type: "space".to_string(),
                    unit: Some("meter".to_string()),
                },
            ],
        }],
    };
    
    let geozarr_path = root.join(".zattrs");
    let geozarr_json = serde_json::to_string_pretty(&geozarr_meta)
        .map_err(|e| AcoliteError::Processing(format!("JSON error: {}", e)))?;
    fs::write(geozarr_path, geozarr_json)
        .map_err(|e| AcoliteError::Io(e))?;
    
    // Write wavelengths as separate array
    let wl_dir = root.join("wavelengths");
    fs::create_dir_all(&wl_dir)
        .map_err(|e| AcoliteError::Io(e))?;
    
    let wl_zarray = ZarrayMetadata {
        zarr_format: 2,
        shape: vec![nbands],
        chunks: vec![nbands],
        dtype: "<f8".to_string(),
        compressor: None,
        fill_value: 0.0,
        order: "C".to_string(),
        filters: None,
    };
    
    let wl_zarray_path = wl_dir.join(".zarray");
    let wl_zarray_json = serde_json::to_string_pretty(&wl_zarray)
        .map_err(|e| AcoliteError::Processing(format!("JSON error: {}", e)))?;
    fs::write(wl_zarray_path, wl_zarray_json)
        .map_err(|e| AcoliteError::Io(e))?;
    
    let wavelengths: Vec<u8> = bands.iter()
        .flat_map(|b| b.wavelength.to_le_bytes())
        .collect();
    
    fs::write(wl_dir.join("0"), wavelengths)
        .map_err(|e| AcoliteError::Io(e))?;
    
    log::info!("GeoZarr written: {}", output_path);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Projection, GeoTransform};
    use ndarray::Array2;
    use chrono::Utc;
    use tempfile::TempDir;
    
    #[test]
    fn test_write_geozarr() {
        let temp = TempDir::new().unwrap();
        let output = temp.path().join("test.zarr");
        
        let proj = Projection::from_epsg(32610);
        let geotrans = GeoTransform::new(500000.0, 30.0, 4000000.0, -30.0);
        
        // Create hyperspectral-like data (many bands)
        let bands: Vec<_> = (0..50)
            .map(|i| {
                BandData::new(
                    Array2::from_elem((100, 100), i as f64 * 0.01),
                    400.0 + i as f64 * 10.0, // 400-890 nm
                    5.0,
                    format!("Band_{}", i + 1),
                    proj.clone(),
                    geotrans.clone(),
                )
            })
            .collect();
        
        let metadata = Metadata::new("HYPERSPECTRAL".to_string(), Utc::now());
        
        let result = write_geozarr(output.to_str().unwrap(), &bands, &metadata);
        assert!(result.is_ok());
        assert!(output.exists());
        assert!(output.join("data").exists());
        assert!(output.join("wavelengths").exists());
    }
}
