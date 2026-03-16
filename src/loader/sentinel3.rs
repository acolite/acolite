//! Sentinel-3 OLCI scene loader
//!
//! Reads .SEN3 bundles: radiance NetCDFs, tie-point grids, instrument data,
//! and xfdumanifest.xml metadata. Implements smile correction per the
//! OLCI L2 ATBD Instrumental Correction specification.

use crate::core::{BandData, GeoTransform, Metadata, Projection};
use crate::error::{AcoliteError, Result};
use ndarray::Array2;
use std::collections::HashMap;
use std::path::Path;

/// OLCI band definitions: (name, wavelength_nm, bandwidth_nm, E0 W/m²/µm)
/// E0 values from OLCI L1 product handbook Table 4.3
pub const OLCI_BANDS: [(&str, f64, f64, f64); 21] = [
    ("Oa01", 400.0, 15.0, 1714.9),
    ("Oa02", 412.5, 10.0, 1872.4),
    ("Oa03", 442.5, 10.0, 1926.6),
    ("Oa04", 490.0, 10.0, 1930.2),
    ("Oa05", 510.0, 10.0, 1804.2),
    ("Oa06", 560.0, 10.0, 1651.5),
    ("Oa07", 620.0, 10.0, 1531.4),
    ("Oa08", 665.0, 10.0, 1475.8),
    ("Oa09", 673.75, 7.5, 1408.9),
    ("Oa10", 681.25, 7.5, 1368.1),
    ("Oa11", 708.75, 10.0, 1247.2),
    ("Oa12", 753.75, 7.5, 1141.1),
    ("Oa13", 761.25, 2.5, 1118.0),
    ("Oa14", 764.375, 3.75, 1099.7),
    ("Oa15", 767.5, 2.5, 1083.4),
    ("Oa16", 778.75, 15.0, 1047.1),
    ("Oa17", 865.0, 20.0, 879.1),
    ("Oa18", 885.0, 10.0, 838.5),
    ("Oa19", 900.0, 10.0, 802.8),
    ("Oa20", 940.0, 20.0, 714.6),
    ("Oa21", 1020.0, 40.0, 603.3),
];

/// Smile correction bounding bands (lower_water, upper_water, switch_water)
/// From OLCI band_info table used in Python ACOLITE
const SMILE_BOUNDS: [(usize, usize, bool); 21] = [
    (1, 2, true),   // Oa01
    (1, 3, true),   // Oa02
    (2, 4, true),   // Oa03
    (3, 5, true),   // Oa04
    (4, 6, true),   // Oa05
    (5, 7, true),   // Oa06
    (6, 8, true),   // Oa07
    (7, 9, true),   // Oa08
    (8, 10, true),  // Oa09
    (9, 11, true),  // Oa10
    (10, 12, true), // Oa11
    (11, 16, true), // Oa12
    (12, 14, false), // Oa13 - O2A absorption, no water correction
    (12, 15, false), // Oa14
    (14, 16, false), // Oa15
    (12, 17, true), // Oa16
    (16, 18, true), // Oa17
    (17, 19, true), // Oa18
    (18, 20, true), // Oa19
    (19, 21, true), // Oa20
    (20, 21, true), // Oa21
];

/// Per-detector instrument data needed for smile correction
#[derive(Debug, Clone)]
pub struct OlciInstrumentData {
    /// Per-detector central wavelength [band][detector] in nm
    pub lambda0: Vec<Vec<f64>>,
    /// Per-detector solar flux [band][detector] in mW/m²/nm
    pub solar_flux: Vec<Vec<f64>>,
    /// Per-detector FWHM [band][detector] in nm
    pub fwhm: Vec<Vec<f64>>,
    /// Detector index per pixel [rows × cols]
    pub detector_index: Array2<u32>,
}

impl OlciInstrumentData {
    /// Subset detector_index to a pixel window (lambda0/solar_flux/fwhm are per-detector, not spatial)
    pub fn subset(&self, r0: usize, c0: usize, nr: usize, nc: usize) -> Self {
        use ndarray::s;
        Self {
            lambda0: self.lambda0.clone(),
            solar_flux: self.solar_flux.clone(),
            fwhm: self.fwhm.clone(),
            detector_index: self.detector_index.slice(s![r0..r0+nr, c0..c0+nc]).to_owned(),
        }
    }
}

/// Tie-point grids interpolated to full resolution
#[derive(Debug, Clone)]
pub struct OlciTpg {
    pub sza: Array2<f64>,
    pub oza: Array2<f64>,
    pub saa: Array2<f64>,
    pub oaa: Array2<f64>,
    pub latitude: Array2<f64>,
    pub longitude: Array2<f64>,
    pub total_ozone: Option<Array2<f64>>,
    pub total_columnar_water_vapour: Option<Array2<f64>>,
    pub sea_level_pressure: Option<Array2<f64>>,
}

impl OlciTpg {
    /// Subset all TPG arrays to a pixel window
    pub fn subset(&self, r0: usize, c0: usize, nr: usize, nc: usize) -> Self {
        use ndarray::s;
        let sub = |a: &Array2<f64>| a.slice(s![r0..r0+nr, c0..c0+nc]).to_owned();
        Self {
            sza: sub(&self.sza), oza: sub(&self.oza),
            saa: sub(&self.saa), oaa: sub(&self.oaa),
            latitude: sub(&self.latitude), longitude: sub(&self.longitude),
            total_ozone: self.total_ozone.as_ref().map(|a| sub(a)),
            total_columnar_water_vapour: self.total_columnar_water_vapour.as_ref().map(|a| sub(a)),
            sea_level_pressure: self.sea_level_pressure.as_ref().map(|a| sub(a)),
        }
    }

    /// Compute relative azimuth angle
    pub fn raa(&self) -> Array2<f64> {
        let mut raa = &self.saa - &self.oaa;
        raa.mapv_inplace(|v| {
            let a = v.abs();
            if a > 180.0 { (360.0 - a).abs() } else { a }
        });
        raa
    }

    /// Mean geometry for scene-level processing
    pub fn mean_geometry(&self) -> (f64, f64, f64) {
        let sza = nanmean(&self.sza);
        let vza = nanmean(&self.oza);
        let raa = nanmean(&self.raa());
        (sza, vza, raa)
    }

    /// Center coordinates
    pub fn center_lonlat(&self) -> (f64, f64) {
        (nanmean(&self.longitude), nanmean(&self.latitude))
    }
}

/// Full OLCI scene
#[derive(Debug)]
pub struct OlciScene {
    pub sensor: String,
    pub metadata: Metadata,
    /// Radiance data per band name
    pub radiance: HashMap<String, Array2<f64>>,
    pub instrument: OlciInstrumentData,
    pub tpg: OlciTpg,
    pub data_shape: (usize, usize),
    pub product_type: Option<String>, // "FR" or "RR"
}

impl OlciScene {
    /// Convert radiance to TOA reflectance: ρ_t = π·L / (F0·cos(θ_s))
    /// Uses per-detector F0 if smile-corrected=false, nominal E0 if smile-corrected=true
    pub fn to_toa_reflectance(&self, smile_corrected: bool) -> HashMap<String, Array2<f64>> {
        let mu = self.tpg.sza.mapv(|v| (v * std::f64::consts::PI / 180.0).cos());
        let mut result = HashMap::new();

        for (i, (name, _wl, _bw, e0)) in OLCI_BANDS.iter().enumerate() {
            let key = name.to_string();
            if let Some(rad) = self.radiance.get(&key) {
                let rhot = if smile_corrected {
                    // Use nominal E0 after smile correction
                    rad.mapv(|v| std::f64::consts::PI * v / e0) / &mu
                } else {
                    // Use per-detector solar flux
                    let mut out = Array2::zeros(rad.dim());
                    let di = &self.instrument.detector_index;
                    for ((r, c), &l) in rad.indexed_iter() {
                        let d = di[[r, c]] as usize;
                        let f0 = self.instrument.solar_flux[i][d.min(self.instrument.solar_flux[i].len() - 1)];
                        out[[r, c]] = std::f64::consts::PI * l / (f0 * mu[[r, c]]);
                    }
                    out
                };
                result.insert(key, rhot);
            }
        }
        result
    }

    /// Apply smile correction (OLCI L2 ATBD Instrumental Correction)
    ///
    /// For each band:
    ///   1. Compute reflectance using per-detector F0
    ///   2. Recompute radiance using nominal E0 → difference = smile effect
    ///   3. For water bands: add bounding-band interpolation correction
    ///
    /// Optionally applies gas transmittance before/after (matching Python behavior)
    pub fn apply_smile_correction(&mut self, tt_gas: Option<&HashMap<String, f64>>) {
        let di = &self.instrument.detector_index;
        let nbands = OLCI_BANDS.len();

        // Step 1: Remove gas transmittance if requested
        if let Some(tg) = tt_gas {
            for (name, _, _, _) in &OLCI_BANDS {
                if let (Some(rad), Some(&t)) = (self.radiance.get_mut(&name.to_string()), tg.get(*name)) {
                    rad.mapv_inplace(|v| v / t);
                }
            }
        }

        // Step 2: Compute smile correction per band
        let mut smile: Vec<Array2<f64>> = Vec::with_capacity(nbands);
        for (i, (name, _wl, _bw, e0)) in OLCI_BANDS.iter().enumerate() {
            let rad = match self.radiance.get(&name.to_string()) {
                Some(r) => r,
                None => { smile.push(Array2::zeros((1, 1))); continue; }
            };
            let (rows, cols) = rad.dim();
            let mut s = Array2::zeros((rows, cols));

            for r in 0..rows {
                for c in 0..cols {
                    let d = di[[r, c]] as usize;
                    let sf = &self.instrument.solar_flux[i];
                    let f0_det = sf[d.min(sf.len() - 1)];
                    // Reflectance with per-detector F0
                    let refl = rad[[r, c]] / f0_det;
                    // Radiance at nominal E0
                    let rad_nominal = refl * e0;
                    s[[r, c]] = rad_nominal - rad[[r, c]];
                }
            }

            // Bounding band correction for water bands
            let (b1_idx, b2_idx, do_water) = SMILE_BOUNDS[i];
            if do_water && b1_idx >= 1 && b2_idx >= 1 && b1_idx <= nbands && b2_idx <= nbands {
                let b1 = b1_idx - 1;
                let b2 = b2_idx - 1;
                if let (Some(rad1), Some(rad2)) = (
                    self.radiance.get(&OLCI_BANDS[b1].0.to_string()),
                    self.radiance.get(&OLCI_BANDS[b2].0.to_string()),
                ) {
                    for r in 0..rows {
                        for c in 0..cols {
                            let d = di[[r, c]] as usize;
                            let sf1 = &self.instrument.solar_flux[b1];
                            let sf2 = &self.instrument.solar_flux[b2];
                            let f0_1 = sf1[d.min(sf1.len() - 1)];
                            let f0_2 = sf2[d.min(sf2.len() - 1)];
                            let r21_diff = rad2[[r, c]] / f0_2 - rad1[[r, c]] / f0_1;

                            let l0 = &self.instrument.lambda0;
                            let l_i = l0[i][d.min(l0[i].len() - 1)];
                            let l_b1 = l0[b1][d.min(l0[b1].len() - 1)];
                            let l_b2 = l0[b2][d.min(l0[b2].len() - 1)];
                            let denom = l_b2 - l_b1;
                            if denom.abs() > 0.001 {
                                let wdiff_ratio = (OLCI_BANDS[i].1 - l_i) / denom;
                                let sf_i = &self.instrument.solar_flux[i];
                                s[[r, c]] += r21_diff * wdiff_ratio * sf_i[d.min(sf_i.len() - 1)];
                            }
                        }
                    }
                }
            }
            smile.push(s);
        }

        // Step 3: Apply smile correction to radiance
        for (i, (name, _, _, _)) in OLCI_BANDS.iter().enumerate() {
            if let Some(rad) = self.radiance.get_mut(&name.to_string()) {
                if smile[i].dim() == rad.dim() {
                    *rad = &*rad + &smile[i];
                }
            }
        }

        // Step 4: Restore gas transmittance
        if let Some(tg) = tt_gas {
            for (name, _, _, _) in &OLCI_BANDS {
                if let (Some(rad), Some(&t)) = (self.radiance.get_mut(&name.to_string()), tg.get(*name)) {
                    rad.mapv_inplace(|v| v * t);
                }
            }
        }
    }

    /// Convert to BandData vector for pipeline processing
    pub fn to_band_data(&self, rhot: &HashMap<String, Array2<f64>>) -> Vec<BandData<f64>> {
        let proj = Projection::from_epsg(4326); // geographic for OLCI
        let (clon, clat) = self.tpg.center_lonlat();
        let geotrans = GeoTransform::new(clon, 0.003, clat, -0.003); // ~300m approx

        OLCI_BANDS
            .iter()
            .filter_map(|(name, wl, bw, _)| {
                rhot.get(&name.to_string()).map(|data| {
                    BandData::new(
                        data.clone(),
                        *wl,
                        *bw,
                        name.to_string(),
                        proj.clone(),
                        geotrans.clone(),
                    )
                })
            })
            .collect()
    }
}

/// Parse xfdumanifest.xml to determine sensor identity
pub fn parse_s3_manifest(manifest_path: &Path) -> Result<(String, HashMap<String, String>)> {
    let content = std::fs::read_to_string(manifest_path)
        .map_err(|e| AcoliteError::Io(e))?;

    let mut attrs = HashMap::new();
    let sensor = if content.contains("Ocean Land Colour Instrument") {
        let number = extract_xml_text(&content, "sentinel-safe:number").unwrap_or_default();
        attrs.insert("instrument".into(), "OLCI".into());
        format!("S3{}_OLCI", number)
    } else if content.contains("MEdium Resolution Imaging Spectrometer") {
        attrs.insert("instrument".into(), "MERIS".into());
        "EN1_MERIS".into()
    } else if content.contains("Sea and Land Surface Temperature Radiometer") {
        let number = extract_xml_text(&content, "sentinel-safe:number").unwrap_or_default();
        attrs.insert("instrument".into(), "SLSTR".into());
        format!("S3{}_SLSTR", number)
    } else {
        return Err(AcoliteError::UnsupportedSensor("Unknown S3 instrument".into()));
    };

    // Detect FR/RR
    for tag in &["sentinel3:productType", "envisat:productType"] {
        if let Some(val) = extract_xml_text(&content, tag) {
            attrs.insert(tag.to_string(), val.clone());
            if val.contains("FR") {
                attrs.insert("product_type".into(), "FR".into());
            } else if val.contains("RR") {
                attrs.insert("product_type".into(), "RR".into());
            }
        }
    }

    Ok((sensor, attrs))
}

/// Load a full OLCI scene from a .SEN3 directory
#[cfg(feature = "netcdf")]
pub fn load_olci_scene(
    sen3_dir: &Path,
    limit: Option<&[f64; 4]>,
) -> Result<OlciScene> {
    // Find manifest
    let manifest = sen3_dir.join("xfdumanifest.xml");
    if !manifest.exists() {
        return Err(AcoliteError::InvalidMetadata(
            format!("No xfdumanifest.xml in {}", sen3_dir.display()),
        ));
    }
    let (sensor, attrs) = parse_s3_manifest(&manifest)?;

    // Load tie-point grids
    let tpg = load_tpg(sen3_dir, limit)?;
    let (sza, vza, raa) = tpg.mean_geometry();
    let (clon, clat) = tpg.center_lonlat();

    // Build metadata
    let mut metadata = Metadata::new(sensor.clone(), chrono::Utc::now());
    metadata.set_geometry(sza, 0.0);
    metadata.view_zenith = Some(vza);
    for (k, v) in &attrs {
        metadata.add_attribute(k.clone(), v.clone());
    }

    // Load instrument data
    let instrument = load_instrument_data(sen3_dir)?;
    let data_shape = (instrument.detector_index.nrows(), instrument.detector_index.ncols());

    // Determine subset if limit provided
    let sub = limit.and_then(|lim| compute_subset(sen3_dir, &tpg, lim));

    // Subset TPGs and instrument data if limit applied
    let (tpg, instrument, data_shape) = if let Some((r0, c0, nr, nc)) = sub {
        let sub_tpg = tpg.subset(r0, c0, nr, nc);
        let sub_inst = instrument.subset(r0, c0, nr, nc);
        let shape = (nr, nc);
        (sub_tpg, sub_inst, shape)
    } else {
        (tpg, instrument, data_shape)
    };

    // Load radiance bands — read only the subset region directly from NetCDF
    let mut radiance = HashMap::new();
    for (_i, (name, _, _, _)) in OLCI_BANDS.iter().enumerate() {
        let nc_name = format!("{}_radiance.nc", name);
        let nc_path = sen3_dir.join(&nc_name);
        if !nc_path.exists() { continue; }

        let ds_name = format!("{}_radiance", name);
        if let Some((r0, c0, nr, nc)) = sub {
            match read_nc_scaled_subset(&nc_path, &ds_name, r0, c0, nr, nc) {
                Ok(data) => { radiance.insert(name.to_string(), data); }
                Err(e) => log::warn!("Failed to read {}: {}", nc_name, e),
            }
        } else {
            match read_nc_scaled(&nc_path, &ds_name) {
                Ok(data) => { radiance.insert(name.to_string(), data); }
                Err(e) => log::warn!("Failed to read {}: {}", nc_name, e),
            }
        }
    }

    Ok(OlciScene {
        sensor,
        metadata,
        radiance,
        instrument,
        tpg,
        data_shape,
        product_type: attrs.get("product_type").cloned(),
    })
}

#[cfg(feature = "netcdf")]
fn load_tpg(sen3_dir: &Path, _limit: Option<&[f64; 4]>) -> Result<OlciTpg> {
    let geo_file = sen3_dir.join("tie_geo_coordinates.nc");
    let met_file = sen3_dir.join("tie_meteo.nc");
    let geom_file = sen3_dir.join("tie_geometries.nc");

    // Read tie-point grids (subsampled)
    let lat_tp = read_nc_2d_f64(&geo_file, "latitude")?;
    let lon_tp = read_nc_2d_f64(&geo_file, "longitude")?;
    let sza_tp = read_nc_2d_f64(&geom_file, "SZA")?;
    let oza_tp = read_nc_2d_f64(&geom_file, "OZA")?;
    let saa_tp = read_nc_2d_f64(&geom_file, "SAA")?;
    let oaa_tp = read_nc_2d_f64(&geom_file, "OAA")?;

    // Determine full resolution from instrument_data detector_index
    let inst_file = sen3_dir.join("instrument_data.nc");
    let nc_inst = netcdf::open(&inst_file)
        .map_err(|e| AcoliteError::NetCdf(format!("{}: {}", inst_file.display(), e)))?;
    let di_var = nc_inst.variable("detector_index")
        .ok_or_else(|| AcoliteError::NetCdf("No detector_index".into()))?;
    let full_rows = di_var.dimensions()[0].len();
    let full_cols = di_var.dimensions()[1].len();
    drop(nc_inst);

    // Read subsampling factors from a radiance file (OLCI standard: ac=64, al=1)
    let (ac_sub, al_sub) = read_subsampling_factors(sen3_dir).unwrap_or((64, 1));

    // Interpolate TPGs to full resolution using bilinear interpolation
    // matching Python's RegularGridInterpolator with proper subsampling factors
    // and half-pixel offset (pixel centres)
    let sza = interp_tpg_to_full(&sza_tp, full_rows, full_cols, ac_sub, al_sub);
    let oza = interp_tpg_to_full(&oza_tp, full_rows, full_cols, ac_sub, al_sub);
    let saa = interp_tpg_to_full(&saa_tp, full_rows, full_cols, ac_sub, al_sub);
    let oaa = interp_tpg_to_full(&oaa_tp, full_rows, full_cols, ac_sub, al_sub);
    let latitude = interp_tpg_to_full(&lat_tp, full_rows, full_cols, ac_sub, al_sub);
    let longitude = interp_tpg_to_full(&lon_tp, full_rows, full_cols, ac_sub, al_sub);

    let total_ozone = read_nc_2d_f64(&met_file, "total_ozone")
        .ok().map(|tp| interp_tpg_to_full(&tp, full_rows, full_cols, ac_sub, al_sub));
    let total_columnar_water_vapour = read_nc_2d_f64(&met_file, "total_columnar_water_vapour")
        .ok().map(|tp| interp_tpg_to_full(&tp, full_rows, full_cols, ac_sub, al_sub));
    let sea_level_pressure = read_nc_2d_f64(&met_file, "sea_level_pressure")
        .ok().map(|tp| interp_tpg_to_full(&tp, full_rows, full_cols, ac_sub, al_sub));

    Ok(OlciTpg {
        sza, oza, saa, oaa, latitude, longitude,
        total_ozone, total_columnar_water_vapour, sea_level_pressure,
    })
}

/// Read ac/al subsampling factors from a radiance NetCDF file.
#[cfg(feature = "netcdf")]
fn read_subsampling_factors(sen3_dir: &Path) -> Option<(usize, usize)> {
    let p = sen3_dir.join("Oa01_radiance.nc");
    let nc = netcdf::open(&p).ok()?;
    let ac = nc.attribute("ac_subsampling_factor")
        .and_then(|a| a.value().ok())
        .and_then(|v| match v {
            netcdf::AttributeValue::Int(i) => Some(i as usize),
            netcdf::AttributeValue::Short(i) => Some(i as usize),
            _ => None,
        })?;
    let al = nc.attribute("al_subsampling_factor")
        .and_then(|a| a.value().ok())
        .and_then(|v| match v {
            netcdf::AttributeValue::Int(i) => Some(i as usize),
            netcdf::AttributeValue::Short(i) => Some(i as usize),
            _ => None,
        })?;
    Some((ac, al))
}

/// Bilinear interpolation of a tie-point grid to full resolution.
///
/// Matches Python ACOLITE's RegularGridInterpolator approach:
///   tpx = arange(tp_cols) * ac_sub
///   tpy = arange(tp_rows) * al_sub
///   evaluate at (pixel + 0.5) for each full-res pixel
fn interp_tpg_to_full(
    tpg: &Array2<f64>, full_rows: usize, full_cols: usize,
    ac_sub: usize, al_sub: usize,
) -> Array2<f64> {
    let (tp_rows, tp_cols) = tpg.dim();
    let mut out = Array2::zeros((full_rows, full_cols));
    let ac_sub_f = ac_sub as f64;
    let al_sub_f = al_sub as f64;
    let tp_rows_m1 = (tp_rows - 1) as f64;
    let tp_cols_m1 = (tp_cols - 1) as f64;

    for r in 0..full_rows {
        // Python: suby = r + 0.5, tpy = arange(tp_rows) * al_sub
        // Normalised: ry = (r + 0.5) / al_sub
        let ry = (r as f64 + 0.5) / al_sub_f;
        let ry = ry.clamp(0.0, tp_rows_m1);
        let ry0 = (ry as usize).min(tp_rows - 2);
        let fy = ry - ry0 as f64;

        for c in 0..full_cols {
            let cx = (c as f64 + 0.5) / ac_sub_f;
            let cx = cx.clamp(0.0, tp_cols_m1);
            let cx0 = (cx as usize).min(tp_cols - 2);
            let fx = cx - cx0 as f64;

            let v00 = tpg[[ry0, cx0]];
            let v01 = tpg[[ry0, cx0 + 1]];
            let v10 = tpg[[ry0 + 1, cx0]];
            let v11 = tpg[[ry0 + 1, cx0 + 1]];

            out[[r, c]] = v00 * (1.0 - fx) * (1.0 - fy)
                        + v01 * fx * (1.0 - fy)
                        + v10 * (1.0 - fx) * fy
                        + v11 * fx * fy;
        }
    }
    out
}

#[cfg(feature = "netcdf")]
fn load_instrument_data(sen3_dir: &Path) -> Result<OlciInstrumentData> {
    let inst_file = sen3_dir.join("instrument_data.nc");
    let nc = netcdf::open(&inst_file)
        .map_err(|e| AcoliteError::NetCdf(format!("{}: {}", inst_file.display(), e)))?;

    let lambda0 = read_nc_2d_f32_to_vecs(&nc, "lambda0")?;
    let solar_flux = read_nc_2d_f32_to_vecs(&nc, "solar_flux")?;
    let fwhm = read_nc_2d_f32_to_vecs(&nc, "FWHM")?;

    let di_var = nc.variable("detector_index")
        .ok_or_else(|| AcoliteError::NetCdf("No detector_index".into()))?;
    let di_raw: Vec<i16> = di_var.get_values(..)
        .map_err(|e| AcoliteError::NetCdf(e.to_string()))?;
    let dims = di_var.dimensions();
    let (rows, cols) = (dims[0].len(), dims[1].len());
    let di_u32: Vec<u32> = di_raw.iter().map(|&v| v.max(0) as u32).collect();
    let detector_index = Array2::from_shape_vec((rows, cols), di_u32)
        .map_err(|e| AcoliteError::Processing(e.to_string()))?;

    Ok(OlciInstrumentData { lambda0, solar_flux, fwhm, detector_index })
}

#[cfg(feature = "netcdf")]
fn read_nc_2d_f64(path: &Path, var_name: &str) -> Result<Array2<f64>> {
    let nc = netcdf::open(path)
        .map_err(|e| AcoliteError::NetCdf(format!("{}: {}", path.display(), e)))?;
    let var = nc.variable(var_name)
        .ok_or_else(|| AcoliteError::NetCdf(format!("No variable {} in {}", var_name, path.display())))?;
    let dims = var.dimensions();
    let (rows, cols) = (dims[0].len(), dims[1].len());

    // Check for CF-convention scale_factor/add_offset → data is packed integer
    let scale: f64 = var.attribute("scale_factor")
        .and_then(|a| a.value().ok())
        .map(|v| match v { netcdf::AttributeValue::Float(f) => f as f64, netcdf::AttributeValue::Double(d) => d, _ => 1.0 })
        .unwrap_or(1.0);
    let offset: f64 = var.attribute("add_offset")
        .and_then(|a| a.value().ok())
        .map(|v| match v { netcdf::AttributeValue::Float(f) => f as f64, netcdf::AttributeValue::Double(d) => d, _ => 0.0 })
        .unwrap_or(0.0);

    if scale != 1.0 || offset != 0.0 {
        // Packed integer data — try integer types
        if let Ok(data) = var.get_values::<u32, _>(..) {
            let out: Vec<f64> = data.iter().map(|&v| v as f64 * scale + offset).collect();
            return Array2::from_shape_vec((rows, cols), out).map_err(|e| AcoliteError::Processing(e.to_string()));
        }
        if let Ok(data) = var.get_values::<i32, _>(..) {
            let out: Vec<f64> = data.iter().map(|&v| v as f64 * scale + offset).collect();
            return Array2::from_shape_vec((rows, cols), out).map_err(|e| AcoliteError::Processing(e.to_string()));
        }
        if let Ok(data) = var.get_values::<u16, _>(..) {
            let out: Vec<f64> = data.iter().map(|&v| if v == 65535 { f64::NAN } else { v as f64 * scale + offset }).collect();
            return Array2::from_shape_vec((rows, cols), out).map_err(|e| AcoliteError::Processing(e.to_string()));
        }
        if let Ok(data) = var.get_values::<i16, _>(..) {
            let out: Vec<f64> = data.iter().map(|&v| v as f64 * scale + offset).collect();
            return Array2::from_shape_vec((rows, cols), out).map_err(|e| AcoliteError::Processing(e.to_string()));
        }
    }

    // Native float data (no scaling)
    if let Ok(data) = var.get_values::<f32, _>(..) {
        let out: Vec<f64> = data.iter().map(|&v| v as f64).collect();
        return Array2::from_shape_vec((rows, cols), out).map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    if let Ok(data) = var.get_values::<f64, _>(..) {
        return Array2::from_shape_vec((rows, cols), data).map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    Err(AcoliteError::NetCdf(format!("Cannot read {} in {}", var_name, path.display())))
}

#[cfg(feature = "netcdf")]
fn read_nc_2d_f32_to_vecs(nc: &netcdf::File, var_name: &str) -> Result<Vec<Vec<f64>>> {
    let var = nc.variable(var_name)
        .ok_or_else(|| AcoliteError::NetCdf(format!("No variable {}", var_name)))?;
    let data: Vec<f32> = var.get_values(..)
        .map_err(|e| AcoliteError::NetCdf(e.to_string()))?;
    let dims = var.dimensions();
    let (nbands, ndet) = (dims[0].len(), dims[1].len());
    Ok((0..nbands)
        .map(|b| data[b * ndet..(b + 1) * ndet].iter().map(|&v| v as f64).collect())
        .collect())
}

/// Read a 2D NetCDF variable, always applying scale_factor and add_offset.
/// Reads as native type (u16/i16) to avoid the netcdf crate silently casting
/// without applying CF-convention scaling.
#[cfg(feature = "netcdf")]
fn read_nc_scaled(path: &Path, var_name: &str) -> Result<Array2<f64>> {
    let nc = netcdf::open(path)
        .map_err(|e| AcoliteError::NetCdf(format!("{}: {}", path.display(), e)))?;
    let var = nc.variable(var_name)
        .ok_or_else(|| AcoliteError::NetCdf(format!("No variable {} in {}", var_name, path.display())))?;
    let dims = var.dimensions();
    let (rows, cols) = (dims[0].len(), dims[1].len());

    let scale: f64 = var.attribute("scale_factor")
        .and_then(|a| a.value().ok())
        .map(|v| match v {
            netcdf::AttributeValue::Float(f) => f as f64,
            netcdf::AttributeValue::Double(d) => d,
            _ => 1.0,
        })
        .unwrap_or(1.0);
    let offset: f64 = var.attribute("add_offset")
        .and_then(|a| a.value().ok())
        .map(|v| match v {
            netcdf::AttributeValue::Float(f) => f as f64,
            netcdf::AttributeValue::Double(d) => d,
            _ => 0.0,
        })
        .unwrap_or(0.0);

    // Try u16 first (OLCI radiance), then i16, then f32, then f64
    if let Ok(data) = var.get_values::<u16, _>(..) {
        let fill: u16 = 65535; // OLCI standard fill
        let f64_data: Vec<f64> = data.iter().map(|&v| {
            if v == fill { f64::NAN } else { v as f64 * scale + offset }
        }).collect();
        return Array2::from_shape_vec((rows, cols), f64_data)
            .map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    if let Ok(data) = var.get_values::<i16, _>(..) {
        let f64_data: Vec<f64> = data.iter().map(|&v| v as f64 * scale + offset).collect();
        return Array2::from_shape_vec((rows, cols), f64_data)
            .map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    if let Ok(data) = var.get_values::<f32, _>(..) {
        let f64_data: Vec<f64> = data.iter().map(|&v| {
            if v.is_nan() { f64::NAN } else { v as f64 * scale + offset }
        }).collect();
        return Array2::from_shape_vec((rows, cols), f64_data)
            .map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    // f64 — scale should already be 1.0 if native f64
    let data: Vec<f64> = var.get_values(..)
        .map_err(|e| AcoliteError::NetCdf(format!("Cannot read {}: {}", var_name, e)))?;
    Array2::from_shape_vec((rows, cols), data)
        .map_err(|e| AcoliteError::Processing(e.to_string()))
}

/// Read a 2D subset from a NetCDF variable, applying scale_factor and add_offset.
/// Reads only the requested region [r0..r0+nr, c0..c0+nc] from disk.
#[cfg(feature = "netcdf")]
fn read_nc_scaled_subset(
    path: &Path, var_name: &str,
    r0: usize, c0: usize, nr: usize, nc_count: usize,
) -> Result<Array2<f64>> {
    let nc = netcdf::open(path)
        .map_err(|e| AcoliteError::NetCdf(format!("{}: {}", path.display(), e)))?;
    let var = nc.variable(var_name)
        .ok_or_else(|| AcoliteError::NetCdf(format!("No variable {} in {}", var_name, path.display())))?;

    let scale: f64 = var.attribute("scale_factor")
        .and_then(|a| a.value().ok())
        .map(|v| match v {
            netcdf::AttributeValue::Float(f) => f as f64,
            netcdf::AttributeValue::Double(d) => d,
            _ => 1.0,
        })
        .unwrap_or(1.0);
    let offset: f64 = var.attribute("add_offset")
        .and_then(|a| a.value().ok())
        .map(|v| match v {
            netcdf::AttributeValue::Float(f) => f as f64,
            netcdf::AttributeValue::Double(d) => d,
            _ => 0.0,
        })
        .unwrap_or(0.0);

    let extent = ndarray::s![r0..r0+nr, c0..c0+nc_count];

    if let Ok(data) = var.get_values::<u16, _>(extent) {
        let fill: u16 = 65535;
        let f64_data: Vec<f64> = data.iter().map(|&v| {
            if v == fill { f64::NAN } else { v as f64 * scale + offset }
        }).collect();
        return Array2::from_shape_vec((nr, nc_count), f64_data)
            .map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    if let Ok(data) = var.get_values::<i16, _>(extent) {
        let f64_data: Vec<f64> = data.iter().map(|&v| v as f64 * scale + offset).collect();
        return Array2::from_shape_vec((nr, nc_count), f64_data)
            .map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    if let Ok(data) = var.get_values::<f32, _>(extent) {
        let f64_data: Vec<f64> = data.iter().map(|&v| {
            if v.is_nan() { f64::NAN } else { v as f64 * scale + offset }
        }).collect();
        return Array2::from_shape_vec((nr, nc_count), f64_data)
            .map_err(|e| AcoliteError::Processing(e.to_string()));
    }
    let data: Vec<f64> = var.get_values(extent)
        .map_err(|e| AcoliteError::NetCdf(format!("Cannot read {}: {}", var_name, e)))?;
    Array2::from_shape_vec((nr, nc_count), data)
        .map_err(|e| AcoliteError::Processing(e.to_string()))
}

fn compute_subset(sen3_dir: &Path, tpg: &OlciTpg, limit: &[f64; 4]) -> Option<(usize, usize, usize, usize)> {
    let (south, west, north, east) = (limit[0], limit[1], limit[2], limit[3]);

    // Try full-resolution geo_coordinates.nc first (matches Python use_tpg=False default)
    let geo_path = sen3_dir.join("geo_coordinates.nc");
    if let Ok(nc) = netcdf::open(&geo_path) {
        if let (Some(lat_var), Some(lon_var)) = (nc.variable("latitude"), nc.variable("longitude")) {
            // Read as i32 and apply scale_factor (stored as scaled integers)
            let scale = |var: &netcdf::Variable| -> Option<(f64, ndarray::ArrayD<i32>)> {
                let sf = var.attribute("scale_factor")
                    .and_then(|a| a.value().ok())
                    .and_then(|v| match v {
                        netcdf::AttributeValue::Double(d) => Some(d),
                        netcdf::AttributeValue::Float(f) => Some(f as f64),
                        _ => None,
                    })
                    .unwrap_or(1.0);
                var.get::<i32, _>(..).ok().map(|arr| (sf, arr))
            };
            if let (Some((lat_sf, lat_raw)), Some((lon_sf, lon_raw))) = (scale(&lat_var), scale(&lon_var)) {
                let (rows, cols) = (lat_raw.shape()[0], lat_raw.shape()[1]);
                let mut r_min = rows;
                let mut r_max = 0usize;
                let mut c_min = cols;
                let mut c_max = 0usize;
                for r in 0..rows {
                    for c in 0..cols {
                        let la = lat_raw[[r, c]] as f64 * lat_sf;
                        let lo = lon_raw[[r, c]] as f64 * lon_sf;
                        if la >= south && la <= north && lo >= west && lo <= east {
                            r_min = r_min.min(r);
                            r_max = r_max.max(r);
                            c_min = c_min.min(c);
                            c_max = c_max.max(c);
                        }
                    }
                }
                if r_max > r_min && c_max > c_min {
                    // Match Python geolocation_sub: returns (max - min), not (max - min + 1)
                    return Some((r_min, c_min, r_max - r_min, c_max - c_min));
                }
                return None;
            }
        }
    }

    // Fallback: use interpolated TPG
    let lat = &tpg.latitude;
    let lon = &tpg.longitude;
    let (rows, cols) = lat.dim();
    let mut r_min = rows;
    let mut r_max = 0usize;
    let mut c_min = cols;
    let mut c_max = 0usize;
    for r in 0..rows {
        for c in 0..cols {
            let la = lat[[r, c]];
            let lo = lon[[r, c]];
            if la >= south && la <= north && lo >= west && lo <= east {
                r_min = r_min.min(r);
                r_max = r_max.max(r);
                c_min = c_min.min(c);
                c_max = c_max.max(c);
            }
        }
    }
    if r_max > r_min && c_max > c_min {
        Some((r_min, c_min, r_max - r_min, c_max - c_min))
    } else {
        None
    }
}

fn extract_xml_text(content: &str, tag: &str) -> Option<String> {
    let open = format!("<{}", tag);
    let close = format!("</{}>", tag);
    if let Some(start) = content.find(&open) {
        let after_open = &content[start..];
        if let Some(gt) = after_open.find('>') {
            let after_gt = &after_open[gt + 1..];
            if let Some(end) = after_gt.find(&close) {
                return Some(after_gt[..end].trim().to_string());
            }
        }
    }
    None
}

fn nanmean(arr: &Array2<f64>) -> f64 {
    let (sum, count) = arr.iter().fold((0.0, 0u64), |(s, c), &v| {
        if v.is_finite() { (s + v, c + 1) } else { (s, c) }
    });
    if count > 0 { sum / count as f64 } else { f64::NAN }
}
