//! PACE OCI L1B NetCDF reader
//!
//! Reads PACE_OCI.*.L1B.*.nc files following the Python ACOLITE `pace/l1_convert.py` structure.
//! Groups: geolocation_data (lat/lon/geometry), sensor_band_parameters (wavelengths/F0),
//!         observation_data (rhot_blue, rhot_red, rhot_SWIR)

use crate::core::{BandData, GeoTransform, Metadata, Projection};
use crate::{AcoliteError, Result};
use chrono::Utc;
use ndarray::Array2;
use std::path::Path;

/// A loaded PACE OCI scene
pub struct PaceScene {
    pub metadata: Metadata,
    pub bands: Vec<BandData<f32>>,
    pub lat: Array2<f64>,
    pub lon: Array2<f64>,
    pub sza: Array2<f32>,
    pub saa: Array2<f32>,
    pub vza: Array2<f32>,
    pub vaa: Array2<f32>,
}

/// Subset specification: (row_offset, col_offset, nrows, ncols)
pub type SubsetSpec = (usize, usize, usize, usize);

/// Load a PACE OCI L1B NetCDF file
///
/// Reads rhot (TOA reflectance) from blue, red, and SWIR detectors,
/// along with geolocation and geometry data.
/// Optionally subset with `limit` as [south, west, north, east].
#[cfg(feature = "netcdf")]
pub fn load_pace_l1b(path: &Path, limit: Option<&[f64; 4]>) -> Result<PaceScene> {
    let file = netcdf::open(path)
        .map_err(|e| AcoliteError::NetCdf(format!("Open {:?}: {}", path, e)))?;

    // Read global attributes for sensor/time
    let platform = read_attr_string(&file, "platform").unwrap_or_else(|| "PACE".into());
    let instrument = read_attr_string(&file, "instrument").unwrap_or_else(|| "OCI".into());
    let sensor = format!("{}_{}", platform, instrument);
    let isodate = read_attr_string(&file, "time_coverage_start")
        .unwrap_or_else(|| Utc::now().to_rfc3339());
    let datetime = chrono::DateTime::parse_from_rfc3339(&isodate)
        .map(|dt| dt.with_timezone(&Utc))
        .unwrap_or_else(|_| Utc::now());

    let processing_level = read_attr_string(&file, "processing_level").unwrap_or_default();
    if processing_level != "L1B" {
        return Err(AcoliteError::Processing(
            format!("Expected L1B, got '{}'", processing_level),
        ));
    }

    // Read geolocation group
    let geo = file.group("geolocation_data")
        .map_err(|e| AcoliteError::NetCdf(format!("geolocation_data: {}", e)))?
        .ok_or_else(|| AcoliteError::NetCdf("No geolocation_data group".into()))?;

    let lat_full = read_2d_f64(&geo, "latitude")?;
    let lon_full = read_2d_f64(&geo, "longitude")?;

    // Determine subset
    let sub = if let Some(lim) = limit {
        find_subset(&lat_full, &lon_full, lim)
    } else {
        None
    };

    let (lat, lon) = if let Some(s) = &sub {
        (subset_2d(&lat_full, s), subset_2d(&lon_full, s))
    } else {
        (lat_full, lon_full)
    };

    let lat = flip_ud(&lat);
    let lon = flip_ud(&lon);
    let (nrows, ncols) = (lat.nrows(), lat.ncols());

    // Read geometry
    let sza = flip_ud_f32(&read_geo_f32(&geo, "solar_zenith", &sub)?);
    let saa = flip_ud_f32(&read_geo_f32(&geo, "solar_azimuth", &sub)?);
    let vza = flip_ud_f32(&read_geo_f32(&geo, "sensor_zenith", &sub)?);
    let vaa = flip_ud_f32(&read_geo_f32(&geo, "sensor_azimuth", &sub)?);

    // Build metadata from center pixel
    let mut metadata = Metadata::new(sensor.clone(), datetime);
    let (cy, cx) = (nrows / 2, ncols / 2);
    metadata.set_geometry(sza[[cy, cx]] as f64, saa[[cy, cx]] as f64);
    metadata.view_zenith = Some(vza[[cy, cx]] as f64);
    metadata.view_azimuth = Some(vaa[[cy, cx]] as f64);
    metadata.add_attribute("processing_level".into(), "L1B".into());
    metadata.add_attribute("isodate".into(), isodate);

    // Approximate geotransform from lat/lon corners
    let proj = Projection::from_epsg(4326);
    let geotrans = GeoTransform::new(
        lon[[0, 0]],
        if ncols > 1 { (lon[[0, ncols - 1]] - lon[[0, 0]]) / ncols as f64 } else { 0.01 },
        lat[[0, 0]],
        if nrows > 1 { (lat[[nrows - 1, 0]] - lat[[0, 0]]) / nrows as f64 } else { -0.01 },
    );

    // Read band data from all three detectors
    let obs = file.group("observation_data")
        .map_err(|e| AcoliteError::NetCdf(format!("observation_data: {}", e)))?
        .ok_or_else(|| AcoliteError::NetCdf("No observation_data group".into()))?;
    let bp = file.group("sensor_band_parameters")
        .map_err(|e| AcoliteError::NetCdf(format!("sensor_band_parameters: {}", e)))?
        .ok_or_else(|| AcoliteError::NetCdf("No sensor_band_parameters group".into()))?;

    let mut bands = Vec::new();

    for det in &["blue", "red", "SWIR"] {
        let wv_name = format!("{}_wavelength", det);
        let rhot_name = format!("rhot_{}", det);

        let wavelengths = read_1d_f32(&bp, &wv_name)?;
        let bandwidths = if *det == "SWIR" {
            read_1d_f32(&bp, "SWIR_bandpass").unwrap_or_else(|_| vec![20.0; wavelengths.len()])
        } else {
            vec![5.0; wavelengths.len()]
        };

        // rhot is 3D: (nbands, nrows_full, ncols_full)
        let rhot_var = obs.variable(&rhot_name)
            .ok_or_else(|| AcoliteError::NetCdf(format!("Missing {}", rhot_name)))?;
        let dims = rhot_var.dimensions();
        if dims.len() != 3 {
            return Err(AcoliteError::NetCdf(format!("{} not 3D", rhot_name)));
        }
        let nbands_det = dims[0].len();
        let full_rows = dims[1].len();
        let full_cols = dims[2].len();

        for wi in 0..nbands_det {
            if wi >= wavelengths.len() { break; }
            let wave = wavelengths[wi];
            if !wave.is_finite() { continue; }
            let bw = if wi < bandwidths.len() { bandwidths[wi] } else { 5.0 };

            // Read single band slice [wi, :, :] or [wi, sub_rows, sub_cols]
            let (r_start, r_end, c_start, c_end) = if let Some(s) = &sub {
                (s.0, s.0 + s.2, s.1, s.1 + s.3)
            } else {
                (0, full_rows, 0, full_cols)
            };

            let band_data: Vec<f32> = rhot_var
                .get_values((wi..wi + 1, r_start..r_end, c_start..c_end))
                .map_err(|e| AcoliteError::NetCdf(format!("Read {}[{}]: {}", rhot_name, wi, e)))?;

            let data = Array2::from_shape_vec((nrows, ncols), band_data)
                .map_err(|e| AcoliteError::Processing(format!("Shape: {}", e)))?;
            let data = flip_ud_f32(&data);

            let name = format!("rhot_{}_{:.0}", det, wave);
            bands.push(BandData::new(data, wave as f64, bw as f64, name, proj.clone(), geotrans.clone()));
        }
        log::info!("Read {} bands from {} detector", nbands_det, det);
    }

    log::info!("Loaded PACE OCI: {} bands, {}×{} pixels", bands.len(), nrows, ncols);

    Ok(PaceScene { metadata, bands, lat, lon, sza, saa, vza, vaa })
}

// --- Helper functions ---

#[cfg(feature = "netcdf")]
fn read_attr_string(file: &netcdf::File, name: &str) -> Option<String> {
    file.attribute(name)
        .and_then(|a| a.value().ok())
        .and_then(|v| match v {
            netcdf::AttributeValue::Str(s) => Some(s),
            _ => None,
        })
}

#[cfg(feature = "netcdf")]
fn read_2d_f64(group: &netcdf::Group, name: &str) -> Result<Array2<f64>> {
    let var = group.variable(name)
        .ok_or_else(|| AcoliteError::NetCdf(format!("Missing {}", name)))?;
    let dims = var.dimensions();
    let (nrows, ncols) = (dims[0].len(), dims[1].len());
    let data: Vec<f64> = var.get_values(..)
        .map_err(|e| AcoliteError::NetCdf(format!("Read {}: {}", name, e)))?;
    Array2::from_shape_vec((nrows, ncols), data)
        .map_err(|e| AcoliteError::Processing(format!("Shape {}: {}", name, e)))
}

#[cfg(feature = "netcdf")]
fn read_geo_f32(group: &netcdf::Group, name: &str, sub: &Option<SubsetSpec>) -> Result<Array2<f32>> {
    let var = group.variable(name)
        .ok_or_else(|| AcoliteError::NetCdf(format!("Missing {}", name)))?;
    let dims = var.dimensions();
    let (full_rows, full_cols) = (dims[0].len(), dims[1].len());

    let (data, nrows, ncols) = if let Some(s) = sub {
        let d: Vec<f32> = var.get_values((s.0..s.0 + s.2, s.1..s.1 + s.3))
            .map_err(|e| AcoliteError::NetCdf(format!("Read {}: {}", name, e)))?;
        (d, s.2, s.3)
    } else {
        let d: Vec<f32> = var.get_values((0..full_rows, 0..full_cols))
            .map_err(|e| AcoliteError::NetCdf(format!("Read {}: {}", name, e)))?;
        (d, full_rows, full_cols)
    };

    Array2::from_shape_vec((nrows, ncols), data)
        .map_err(|e| AcoliteError::Processing(format!("Shape {}: {}", name, e)))
}

#[cfg(feature = "netcdf")]
fn read_1d_f32(group: &netcdf::Group, name: &str) -> Result<Vec<f32>> {
    let var = group.variable(name)
        .ok_or_else(|| AcoliteError::NetCdf(format!("Missing {}", name)))?;
    var.get_values(..)
        .map_err(|e| AcoliteError::NetCdf(format!("Read {}: {}", name, e)))
}

fn find_subset(lat: &Array2<f64>, lon: &Array2<f64>, limit: &[f64; 4]) -> Option<SubsetSpec> {
    let (south, west, north, east) = (limit[0], limit[1], limit[2], limit[3]);
    let (nrows, ncols) = (lat.nrows(), lat.ncols());
    let (mut rmin, mut rmax, mut cmin, mut cmax) = (nrows, 0usize, ncols, 0usize);

    for r in 0..nrows {
        for c in 0..ncols {
            let (la, lo) = (lat[[r, c]], lon[[r, c]]);
            if la >= south && la <= north && lo >= west && lo <= east {
                rmin = rmin.min(r);
                rmax = rmax.max(r);
                cmin = cmin.min(c);
                cmax = cmax.max(c);
            }
        }
    }

    if rmax >= rmin && cmax >= cmin {
        Some((rmin, cmin, rmax - rmin + 1, cmax - cmin + 1))
    } else {
        None
    }
}

fn subset_2d(arr: &Array2<f64>, sub: &SubsetSpec) -> Array2<f64> {
    arr.slice(ndarray::s![sub.0..sub.0 + sub.2, sub.1..sub.1 + sub.3]).to_owned()
}

fn flip_ud(arr: &Array2<f64>) -> Array2<f64> {
    let n = arr.nrows();
    let mut out = arr.clone();
    for i in 0..n / 2 {
        for j in 0..out.ncols() { out.swap([i, j], [n - 1 - i, j]); }
    }
    out
}

fn flip_ud_f32(arr: &Array2<f32>) -> Array2<f32> {
    let n = arr.nrows();
    let mut out = arr.clone();
    for i in 0..n / 2 {
        for j in 0..out.ncols() { out.swap([i, j], [n - 1 - i, j]); }
    }
    out
}
