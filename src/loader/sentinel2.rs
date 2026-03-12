//! Sentinel-2 SAFE format loader
//!
//! Reads JP2 band files from a .SAFE directory and parses XML metadata.

use crate::core::{BandData, GeoTransform, Metadata, Projection};
use crate::error::{AcoliteError, Result};
use gdal::Dataset;
use ndarray::Array2;
use std::collections::HashMap;
use std::path::Path;

/// S2 band definitions: (band_name, lut_band, wavelength_nm, bandwidth_nm, native_res_m)
const S2_BANDS: &[(&str, &str, f64, f64, u32)] = &[
    ("B01", "1",  442.7,  21.0, 60),
    ("B02", "2",  492.4,  66.0, 10),
    ("B03", "3",  559.8,  36.0, 10),
    ("B04", "4",  664.6,  31.0, 10),
    ("B05", "5",  704.1,  15.0, 20),
    ("B06", "6",  740.5,  15.0, 20),
    ("B07", "7",  782.8,  20.0, 20),
    ("B08", "8",  832.8, 106.0, 10),
    ("B8A", "8A", 864.7,  21.0, 20),
    ("B09", "9",  945.1,  20.0, 60),
    ("B10", "10", 1373.5, 31.0, 60),
    ("B11", "11", 1613.7, 91.0, 20),
    ("B12", "12", 2202.4, 175.0, 20),
];

/// Bands used for DSF atmospheric correction (skip B09=water vapour, B10=cirrus)
pub const S2_AC_BANDS: &[&str] = &["B01","B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12"];

pub struct S2Scene {
    pub bands: Vec<BandData<u16>>,
    pub metadata: Metadata,
    pub sensor_lut: String,
    pub lut_band_map: Vec<(String, String)>,
    /// Quantification value (typically 10000)
    pub quantification_value: f64,
    /// Per-band radiometric add offset (processing baseline ≥ 4.0), keyed by band_id "0".."12"
    pub radio_add_offset: HashMap<String, f64>,
}

/// Load a Sentinel-2 SAFE directory, resampling all bands to `target_res` metres.
pub fn load_sentinel2_scene(safe_dir: &Path, target_res: u32) -> Result<S2Scene> {
    let safe_name = safe_dir
        .file_name()
        .ok_or_else(|| AcoliteError::InvalidMetadata("SAFE path has no filename".into()))?
        .to_string_lossy();

    // Detect sensor
    let sensor_lut = if safe_name.starts_with("S2A") {
        "S2A_MSI"
    } else if safe_name.starts_with("S2B") {
        "S2B_MSI"
    } else if safe_name.starts_with("S2C") {
        "S2C_MSI"
    } else {
        return Err(AcoliteError::Processing(format!("Unknown S2 sensor: {}", safe_name)));
    };

    // Find granule directory
    let granule_dir = safe_dir.join("GRANULE");
    let granule = std::fs::read_dir(&granule_dir)
        .map_err(|e| AcoliteError::Processing(format!("Read GRANULE: {}", e)))?
        .filter_map(|e| e.ok())
        .find(|e| e.file_type().map(|t| t.is_dir()).unwrap_or(false))
        .ok_or_else(|| AcoliteError::Processing("No granule found".into()))?;
    let img_data = granule.path().join("IMG_DATA");

    // Parse granule metadata for geometry
    let metadata = parse_s2_metadata(safe_dir, &granule.path(), sensor_lut)?;

    // Load bands, resampling to target resolution
    let mut bands = Vec::new();
    let mut lut_band_map = Vec::new();

    for &(bname, lut_bn, wl, bw, _native_res) in S2_BANDS {
        // Find JP2 file
        let jp2 = find_band_jp2(&img_data, bname);
        let jp2 = match jp2 {
            Some(p) => p,
            None => {
                log::warn!("Band {} not found in {:?}", bname, img_data);
                continue;
            }
        };

        let ds = Dataset::open(&jp2)
            .map_err(|e| AcoliteError::Gdal(format!("Open {}: {}", bname, e)))?;
        let gt = ds.geo_transform()
            .map_err(|e| AcoliteError::Gdal(format!("GT {}: {}", bname, e)))?;
        let (src_w, src_h) = ds.raster_size();

        // Compute target dimensions
        let src_res = gt[1].abs();
        let scale = src_res / target_res as f64;
        let tgt_w = (src_w as f64 * scale).round() as usize;
        let tgt_h = (src_h as f64 * scale).round() as usize;

        let rb = ds.rasterband(1)
            .map_err(|e| AcoliteError::Gdal(format!("Band {}: {}", bname, e)))?;

        let mut data = Array2::<u16>::zeros((tgt_h, tgt_w));
        rb.read_into_slice(
            (0, 0), (src_w, src_h), (tgt_w, tgt_h),
            data.as_slice_mut()
                .ok_or_else(|| AcoliteError::Processing(format!("Non-contiguous array for {}", bname)))?,
            None,
        ).map_err(|e| AcoliteError::Gdal(format!("Read {}: {}", bname, e)))?;

        let proj = Projection::from_wkt(ds.projection());
        let geotrans = GeoTransform::new(gt[0], target_res as f64, gt[3], -(target_res as f64));

        bands.push(BandData::new(data, wl, bw, bname.to_string(), proj, geotrans));
        lut_band_map.push((bname.to_string(), lut_bn.to_string()));
    }

    // Parse radiometric calibration from product-level metadata
    let (quant, offsets) = parse_s2_radiometric(safe_dir)?;

    Ok(S2Scene {
        bands,
        metadata,
        sensor_lut: sensor_lut.to_string(),
        lut_band_map,
        quantification_value: quant,
        radio_add_offset: offsets,
    })
}

/// Parse QUANTIFICATION_VALUE and RADIO_ADD_OFFSET from MTD_MSIL1C.xml
fn parse_s2_radiometric(safe_dir: &Path) -> Result<(f64, HashMap<String, f64>)> {
    let mtd_path = safe_dir.join("MTD_MSIL1C.xml");
    let content = std::fs::read_to_string(&mtd_path)
        .map_err(|e| AcoliteError::Processing(format!("Read MTD_MSIL1C.xml: {}", e)))?;

    // Parse QUANTIFICATION_VALUE
    let quant = extract_simple_tag(&content, "QUANTIFICATION_VALUE").unwrap_or(10000.0);

    // Parse RADIO_ADD_OFFSET elements: <RADIO_ADD_OFFSET band_id="0">-1000</RADIO_ADD_OFFSET>
    let mut offsets = HashMap::new();
    for part in content.split("<RADIO_ADD_OFFSET").skip(1) {
        if let (Some(bid_start), Some(val_end)) = (part.find("band_id=\""), part.find("</RADIO_ADD_OFFSET>")) {
            let bid_s = bid_start + 9; // len of 'band_id="'
            if let Some(bid_e) = part[bid_s..].find('"') {
                let band_id = &part[bid_s..bid_s + bid_e];
                // Value is between > and </
                if let Some(gt) = part[..val_end].find('>') {
                    if let Ok(v) = part[gt + 1..val_end].trim().parse::<f64>() {
                        offsets.insert(band_id.to_string(), v);
                    }
                }
            }
        }
    }

    Ok((quant, offsets))
}

fn find_band_jp2(img_data: &Path, band_name: &str) -> Option<std::path::PathBuf> {
    let entries = std::fs::read_dir(img_data).ok()?;
    for entry in entries.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if name.ends_with(".jp2") && name.contains(&format!("_{}", band_name)) {
            // Ensure exact band match (B8A vs B8)
            let suffix = format!("_{}.jp2", band_name);
            if name.ends_with(&suffix) {
                return Some(entry.path());
            }
        }
    }
    None
}

fn parse_s2_metadata(safe_dir: &Path, granule_path: &Path, sensor: &str) -> Result<Metadata> {
    use chrono::{DateTime, Utc};

    let mtd_tl = granule_path.join("MTD_TL.xml");
    let content = std::fs::read_to_string(&mtd_tl)
        .map_err(|e| AcoliteError::Processing(format!("Read MTD_TL: {}", e)))?;

    // Parse sun angles
    let sza = extract_xml_value(&content, "ZENITH_ANGLE", "Sun_Angles_Grid")
        .or_else(|| extract_xml_value(&content, "Mean_Sun_Zenith_Angle", ""))
        .or_else(|| extract_simple_tag(&content, "ZENITH_ANGLE"))
        .unwrap_or(30.0);
    let saa = extract_xml_value(&content, "AZIMUTH_ANGLE", "Sun_Angles_Grid")
        .or_else(|| extract_xml_value(&content, "Mean_Sun_Azimuth_Angle", ""))
        .or_else(|| extract_simple_tag(&content, "AZIMUTH_ANGLE"))
        .unwrap_or(150.0);

    // Parse mean viewing angles
    let vza = extract_mean_viewing_angle(&content, "ZENITH_ANGLE").unwrap_or(5.0);
    let vaa = extract_mean_viewing_angle(&content, "AZIMUTH_ANGLE").unwrap_or(100.0);

    // Parse sensing time
    let sensing_time = extract_tag_text(&content, "SENSING_TIME")
        .and_then(|s| DateTime::parse_from_rfc3339(&s).ok())
        .map(|dt| dt.with_timezone(&Utc))
        .unwrap_or_else(Utc::now);

    let mut metadata = Metadata::new(sensor.to_string(), sensing_time);
    metadata.set_geometry(sza, saa);
    metadata.view_zenith = Some(vza);
    metadata.view_azimuth = Some(vaa);

    // Parse MGRS tile from granule dir name
    let granule_name = granule_path.file_name().unwrap_or_default().to_string_lossy();
    if let Some(pos) = granule_name.find("_T") {
        let tile = &granule_name[pos+1..pos+7.min(granule_name.len())];
        metadata.add_attribute("MGRS_TILE".into(), tile.to_string());
    }

    Ok(metadata)
}

fn extract_tag_text(xml: &str, tag: &str) -> Option<String> {
    let open = format!("<{}", tag);
    let close = format!("</{}>", tag);
    if let Some(start) = xml.find(&open) {
        let after_open = &xml[start..];
        if let Some(gt) = after_open.find('>') {
            let content_start = start + gt + 1;
            if let Some(end) = xml[content_start..].find(&close) {
                return Some(xml[content_start..content_start + end].trim().to_string());
            }
        }
    }
    None
}

fn extract_simple_tag(xml: &str, tag: &str) -> Option<f64> {
    extract_tag_text(xml, tag).and_then(|s| s.parse().ok())
}

fn extract_xml_value(xml: &str, tag: &str, context: &str) -> Option<f64> {
    // Find tag within context section
    let search_in = if context.is_empty() {
        xml
    } else if let Some(pos) = xml.find(context) {
        &xml[pos..]
    } else {
        xml
    };
    extract_simple_tag(search_in, tag)
}

fn extract_mean_viewing_angle(xml: &str, angle_type: &str) -> Option<f64> {
    // Look for Mean_Viewing_Incidence_Angle_List section
    let section = "Mean_Viewing_Incidence_Angle_List";
    let start = xml.find(section)?;
    let section_xml = &xml[start..];
    // Collect all angle values and average them
    let tag = angle_type;
    let open = format!("<{}", tag);
    let close = format!("</{}>", tag);
    let mut values = Vec::new();
    let mut pos = 0;
    while let Some(idx) = section_xml[pos..].find(&open) {
        let abs_pos = pos + idx;
        if let Some(gt) = section_xml[abs_pos..].find('>') {
            let content_start = abs_pos + gt + 1;
            if let Some(end) = section_xml[content_start..].find(&close) {
                if let Ok(v) = section_xml[content_start..content_start + end].trim().parse::<f64>() {
                    values.push(v);
                }
            }
        }
        pos = abs_pos + 1;
        // Stop after reasonable number
        if values.len() > 100 { break; }
    }
    if values.is_empty() { return None; }
    Some(values.iter().sum::<f64>() / values.len() as f64)
}
