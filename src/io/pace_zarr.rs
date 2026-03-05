//! Direct PACE to GeoZarr conversion using Python zarr

use crate::{Result, AcoliteError};
use std::process::Command;

/// Convert PACE NetCDF directly to GeoZarr (bypasses slow band-by-band reading)
pub fn pace_to_geozarr(
    nc_path: &str,
    zarr_path: &str,
    band_indices: &[usize],
) -> Result<()> {
    log::info!("Converting PACE to GeoZarr: {} bands", band_indices.len());
    
    let indices_str = band_indices.iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",");
    
    let output = Command::new("python3")
        .arg("-c")
        .arg(format!(
            r#"
import netCDF4 as nc
import zarr
import numpy as np
import json
from pathlib import Path

# Open NetCDF
f = nc.Dataset('{}', 'r')
obs = f.groups['observation_data']
geo = f.groups['geolocation_data']
params = f.groups['sensor_band_parameters']

# Get data
rhot = obs.variables['rhot_blue']
lat = geo.variables['latitude'][:]
lon = geo.variables['longitude'][:]
wavelengths = params.variables['blue_wavelength'][:]

# Band indices to extract
band_idx = [{}]

print(f'Reading {{len(band_idx)}} bands from {{rhot.shape}}')

# Create Zarr store
store = zarr.DirectoryStore('{}')
root = zarr.group(store=store, overwrite=True)

# Write each band
for i, idx in enumerate(band_idx):
    band_data = rhot[idx, :, :]
    wl = wavelengths[idx]
    
    band_name = f'band_{{int(wl)}}'
    z = root.create_dataset(
        band_name,
        data=band_data,
        chunks=(256, 256),
        compressor=zarr.Blosc(cname='zstd', clevel=3),
        dtype='f4'
    )
    z.attrs['wavelength'] = float(wl)
    z.attrs['units'] = 'reflectance'
    
    print(f'  Band {{i+1}}/{{len(band_idx)}}: {{band_name}} ({{wl:.1f}} nm)')

# Write geolocation
root.create_dataset('latitude', data=lat, chunks=(256, 256), compressor=zarr.Blosc(cname='zstd', clevel=3))
root.create_dataset('longitude', data=lon, chunks=(256, 256), compressor=zarr.Blosc(cname='zstd', clevel=3))

# Write metadata
root.attrs['sensor'] = 'PACE_OCI'
root.attrs['num_bands'] = len(band_idx)
root.attrs['shape'] = list(rhot.shape[1:])
root.attrs['crs'] = 'EPSG:4326'

f.close()

# Get size
import os
total_size = sum(os.path.getsize(os.path.join(dp, f)) 
                 for dp, dn, fn in os.walk('{}') 
                 for f in fn)
print(f'\\nGeoZarr written: {{total_size / 1e6:.1f}} MB')
"#,
            nc_path, indices_str, zarr_path, zarr_path
        ))
        .output()
        .map_err(|e| AcoliteError::Processing(format!("Python failed: {}", e)))?;

    if !output.status.success() {
        let err = String::from_utf8_lossy(&output.stderr);
        return Err(AcoliteError::Processing(format!("Conversion failed: {}", err)));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    log::info!("{}", stdout.trim());

    Ok(())
}

/// Convert PACE NetCDF to GeoZarr with atmospheric correction
pub fn pace_to_geozarr_corrected(
    nc_path: &str,
    zarr_path: &str,
    band_indices: &[usize],
) -> Result<()> {
    log::info!("Converting PACE to GeoZarr with AC: {} bands", band_indices.len());
    
    let indices_str = band_indices.iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",");
    
    let output = Command::new("python3")
        .arg("-c")
        .arg(format!(
            r#"
import netCDF4 as nc
import zarr
import numpy as np

f = nc.Dataset('{}', 'r')
obs = f.groups['observation_data']
geo = f.groups['geolocation_data']
params = f.groups['sensor_band_parameters']

rhot = obs.variables['rhot_blue']
lat = geo.variables['latitude'][:]
lon = geo.variables['longitude'][:]
wavelengths = params.variables['blue_wavelength'][:]

band_idx = [{}]

print(f'Processing {{len(band_idx)}} bands with atmospheric correction')

store = zarr.DirectoryStore('{}')
root = zarr.group(store=store, overwrite=True)

# Simple atmospheric correction (Rayleigh approximation)
for i, idx in enumerate(band_idx):
    band_data = rhot[idx, :, :]
    wl = wavelengths[idx]
    
    # Simple Rayleigh correction (placeholder)
    rayleigh = 0.05 * (550.0 / wl) ** 4
    corrected = np.maximum(band_data - rayleigh, 0.0)
    
    band_name = f'rhos_{{int(wl)}}'
    z = root.create_dataset(
        band_name,
        data=corrected,
        chunks=(256, 256),
        compressor=zarr.Blosc(cname='zstd', clevel=3),
        dtype='f4'
    )
    z.attrs['wavelength'] = float(wl)
    z.attrs['units'] = 'surface_reflectance'
    z.attrs['correction'] = 'rayleigh'
    
    print(f'  Band {{i+1}}/{{len(band_idx)}}: {{band_name}} ({{wl:.1f}} nm)')

root.create_dataset('latitude', data=lat, chunks=(256, 256), compressor=zarr.Blosc(cname='zstd', clevel=3))
root.create_dataset('longitude', data=lon, chunks=(256, 256), compressor=zarr.Blosc(cname='zstd', clevel=3))

root.attrs['sensor'] = 'PACE_OCI'
root.attrs['num_bands'] = len(band_idx)
root.attrs['shape'] = list(rhot.shape[1:])
root.attrs['crs'] = 'EPSG:4326'
root.attrs['processing'] = 'atmospheric_correction'

f.close()

import os
total_size = sum(os.path.getsize(os.path.join(dp, f)) 
                 for dp, dn, fn in os.walk('{}') 
                 for f in fn)
print(f'\\nGeoZarr written: {{total_size / 1e6:.1f}} MB')
"#,
            nc_path, indices_str, zarr_path, zarr_path
        ))
        .output()
        .map_err(|e| AcoliteError::Processing(format!("Python failed: {}", e)))?;

    if !output.status.success() {
        let err = String::from_utf8_lossy(&output.stderr);
        return Err(AcoliteError::Processing(format!("Conversion failed: {}", err)));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    log::info!("{}", stdout.trim());

    Ok(())
}
