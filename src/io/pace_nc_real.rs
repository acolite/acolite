//! Real PACE OCI NetCDF reader using Python

use crate::{Result, AcoliteError};
use ndarray::Array2;
use std::process::Command;

/// Read PACE band data from NetCDF (blue bands)
pub fn read_pace_band(nc_path: &str, band_idx: usize) -> Result<Array2<f32>> {
    let output = Command::new("python3")
        .arg("-c")
        .arg(format!(
            r#"
import netCDF4 as nc
import sys

f = nc.Dataset('{}', 'r')
obs = f.groups['observation_data']
data = obs.variables['rhot_blue']

arr = data[{}, :, :]
print(f'{{arr.shape[0]}},{{arr.shape[1]}}')

for row in arr:
    print(','.join(map(str, row.flatten())))

f.close()
"#,
            nc_path, band_idx
        ))
        .output()
        .map_err(|e| AcoliteError::Processing(format!("Python failed: {}", e)))?;

    if !output.status.success() {
        let err = String::from_utf8_lossy(&output.stderr);
        return Err(AcoliteError::Processing(format!("NetCDF read failed: {}", err)));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.trim().lines().collect();

    let dims: Vec<usize> = lines[0].split(',').map(|s| s.parse().unwrap()).collect();
    let (rows, cols) = (dims[0], dims[1]);

    let mut data = Vec::with_capacity(rows * cols);
    for line in &lines[1..] {
        for val in line.split(',') {
            data.push(val.parse::<f32>().unwrap_or(f32::NAN));
        }
    }

    Array2::from_shape_vec((rows, cols), data)
        .map_err(|e| AcoliteError::Processing(format!("Array error: {}", e)))
}

/// Read PACE geolocation
pub fn read_pace_geolocation(nc_path: &str) -> Result<(Array2<f32>, Array2<f32>)> {
    let output = Command::new("python3")
        .arg("-c")
        .arg(format!(
            r#"
import netCDF4 as nc

f = nc.Dataset('{}', 'r')
geo = f.groups['geolocation_data']

lat = geo.variables['latitude'][:]
lon = geo.variables['longitude'][:]

print(f'{{lat.shape[0]}},{{lat.shape[1]}}')

for row in lat:
    print(','.join(map(str, row.flatten())))

print('LONSTART')

for row in lon:
    print(','.join(map(str, row.flatten())))

f.close()
"#,
            nc_path
        ))
        .output()
        .map_err(|e| AcoliteError::Processing(format!("Python failed: {}", e)))?;

    if !output.status.success() {
        return Err(AcoliteError::Processing("Geolocation read failed".to_string()));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.trim().lines().collect();

    let dims: Vec<usize> = lines[0].split(',').map(|s| s.parse().unwrap()).collect();
    let (rows, cols) = (dims[0], dims[1]);

    let lon_start = lines.iter().position(|&l| l == "LONSTART").unwrap();

    let mut lat_data = Vec::with_capacity(rows * cols);
    for line in &lines[1..lon_start] {
        for val in line.split(',') {
            lat_data.push(val.parse::<f32>().unwrap_or(f32::NAN));
        }
    }

    let mut lon_data = Vec::with_capacity(rows * cols);
    for line in &lines[lon_start + 1..] {
        for val in line.split(',') {
            lon_data.push(val.parse::<f32>().unwrap_or(f32::NAN));
        }
    }

    let lat = Array2::from_shape_vec((rows, cols), lat_data)
        .map_err(|e| AcoliteError::Processing(format!("Lat error: {}", e)))?;
    let lon = Array2::from_shape_vec((rows, cols), lon_data)
        .map_err(|e| AcoliteError::Processing(format!("Lon error: {}", e)))?;

    Ok((lat, lon))
}

/// Get PACE dimensions
pub fn get_pace_dimensions(nc_path: &str) -> Result<(usize, usize, usize)> {
    let output = Command::new("python3")
        .arg("-c")
        .arg(format!(
            r#"
import netCDF4 as nc
f = nc.Dataset('{}', 'r')
print(f'{{f.dimensions["blue_bands"].size}},{{f.dimensions["scans"].size}},{{f.dimensions["pixels"].size}}')
f.close()
"#,
            nc_path
        ))
        .output()
        .map_err(|e| AcoliteError::Processing(format!("Python failed: {}", e)))?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    let dims: Vec<usize> = stdout.trim().split(',').map(|s| s.parse().unwrap()).collect();

    Ok((dims[0], dims[1], dims[2]))
}
