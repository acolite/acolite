//! CMR (Common Metadata Repository) search — single implementation

use crate::{AcoliteError, Result};

/// A granule found via CMR search
#[derive(Debug, Clone)]
pub struct CmrGranule {
    pub id: String,
    pub urls: Vec<String>,
}

/// Search NASA CMR for granules by short_name + bbox + temporal
pub fn search_cmr(
    collection: &str,
    bbox: &[f64; 4],
    start_date: &str,
    end_date: &str,
) -> Result<Vec<CmrGranule>> {
    let url = "https://cmr.earthdata.nasa.gov/search/granules.json";
    let params = [
        ("short_name", collection.to_string()),
        (
            "bounding_box",
            format!("{},{},{},{}", bbox[0], bbox[1], bbox[2], bbox[3]),
        ),
        ("temporal", format!("{},{}", start_date, end_date)),
        ("page_size", "50".into()),
    ];
    log::info!(
        "CMR search: {} [{} to {}]",
        collection,
        start_date,
        end_date
    );
    parse_cmr_response(do_cmr_query(url, &params)?)
}

/// Search NASA CMR for granules by echo_collection_id (OBDAAC style)
pub fn search_cmr_collection_id(
    collection_id: &str,
    lat: f64,
    lon: f64,
    start_date: &str,
    end_date: &str,
) -> Result<Vec<CmrGranule>> {
    let url = "https://cmr.earthdata.nasa.gov/search/granules.json";
    // Ensure dates have time component for CMR temporal format
    let start_iso = ensure_iso(start_date);
    let end_iso = ensure_iso(end_date);
    let temporal = format!("{},{}", start_iso, end_iso);
    let params = [
        ("echo_collection_id", collection_id.to_string()),
        ("point", format!("{},{}", lon, lat)),
        ("temporal", temporal),
        ("page_size", "50".into()),
    ];
    log::info!(
        "CMR collection search: {} at ({},{})",
        collection_id,
        lat,
        lon
    );
    parse_cmr_response(do_cmr_query(url, &params)?)
}

/// Search by collection_id + scene name (producerGranuleId)
pub fn search_cmr_scene(collection_id: &str, scene: &str) -> Result<Vec<CmrGranule>> {
    let url = "https://cmr.earthdata.nasa.gov/search/granules.json";
    let params = [
        ("echo_collection_id", collection_id.to_string()),
        ("producerGranuleId", scene.to_string()),
    ];
    log::info!("CMR scene search: {} in {}", scene, collection_id);
    parse_cmr_response(do_cmr_query(url, &params)?)
}

/// Ensure date string has ISO 8601 time component for CMR
fn ensure_iso(date: &str) -> String {
    if date.contains('T') {
        if date.ends_with('Z') {
            date.to_string()
        } else {
            format!("{}Z", date)
        }
    } else {
        format!("{}T00:00:00Z", date)
    }
}

fn do_cmr_query(url: &str, params: &[(&str, String)]) -> Result<serde_json::Value> {
    let client = reqwest::blocking::Client::new();
    let response = client
        .get(url)
        .query(params)
        .send()
        .map_err(|e| AcoliteError::Processing(format!("CMR request failed: {}", e)))?;
    if !response.status().is_success() {
        return Err(AcoliteError::Processing(format!(
            "CMR HTTP {}",
            response.status()
        )));
    }
    response
        .json()
        .map_err(|e| AcoliteError::Processing(format!("CMR parse: {}", e)))
}

fn parse_cmr_response(json: serde_json::Value) -> Result<Vec<CmrGranule>> {
    let entries = json["feed"]["entry"]
        .as_array()
        .ok_or_else(|| AcoliteError::Processing("No entries in CMR response".into()))?;

    let granules: Vec<CmrGranule> = entries
        .iter()
        .filter_map(|entry| {
            let id = entry["producer_granule_id"]
                .as_str()
                .or_else(|| entry["title"].as_str())?
                .to_string();
            let urls: Vec<String> = entry["links"]
                .as_array()?
                .iter()
                .filter_map(|link| {
                    let href = link["href"].as_str()?;
                    if href.starts_with("http")
                        && (href.ends_with(".nc")
                            || href.ends_with(".h5")
                            || href.ends_with(".zip"))
                    {
                        Some(href.to_string())
                    } else {
                        None
                    }
                })
                .collect();
            if urls.is_empty() {
                None
            } else {
                Some(CmrGranule { id, urls })
            }
        })
        .collect();

    log::info!("Found {} granules", granules.len());
    Ok(granules)
}

/// Convenience: search Landsat Collection 2
pub fn search_landsat(bbox: &[f64; 4], start: &str, end: &str) -> Result<Vec<CmrGranule>> {
    search_cmr("LANDSAT_OT_C2_L1", bbox, start, end)
}

/// PACE OCI collection IDs (from data/API/pace_oci_collection_id.json)
pub mod pace_collections {
    /// v3.0 L1B
    pub const L1B_V3: &str = "C3392966952-OB_CLOUD";
    /// v3.0 L2 AOP
    pub const L2_AOP_V3: &str = "C3385049983-OB_CLOUD";
}

/// Search PACE OCI L1B via OBDAAC
pub fn search_pace_l1b(lat: f64, lon: f64, start: &str, end: &str) -> Result<Vec<CmrGranule>> {
    search_cmr_collection_id(pace_collections::L1B_V3, lat, lon, start, end)
}

/// Search PACE OCI L1B by scene name
pub fn search_pace_scene(scene: &str) -> Result<Vec<CmrGranule>> {
    search_cmr_scene(pace_collections::L1B_V3, scene)
}
