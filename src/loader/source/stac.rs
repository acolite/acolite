//! STAC API client — single implementation

use crate::{AcoliteError, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct StacItem {
    pub id: String,
    pub bbox: Vec<f64>,
    pub properties: HashMap<String, serde_json::Value>,
    pub assets: HashMap<String, StacAsset>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct StacAsset {
    pub href: String,
    #[serde(rename = "type")]
    pub media_type: Option<String>,
    pub title: Option<String>,
}

#[derive(Deserialize)]
struct SearchResponse {
    features: Vec<StacItem>,
}

pub struct StacClient {
    base_url: String,
    client: reqwest::blocking::Client,
}

impl StacClient {
    pub fn new(base_url: &str) -> Self {
        Self {
            base_url: base_url.to_string(),
            client: reqwest::blocking::Client::new(),
        }
    }

    pub fn earth_search() -> Self {
        Self::new("https://earth-search.aws.element84.com/v1")
    }

    pub fn usgs() -> Self {
        Self::new("https://landsatlook.usgs.gov/stac-server")
    }

    pub fn search(
        &self,
        collection: &str,
        bbox: &[f64; 4],
        datetime: &str,
        limit: usize,
    ) -> Result<Vec<StacItem>> {
        let body = serde_json::json!({
            "collections": [collection],
            "bbox": bbox,
            "datetime": datetime,
            "limit": limit,
        });

        let response = self
            .client
            .post(format!("{}/search", self.base_url))
            .json(&body)
            .send()
            .map_err(|e| AcoliteError::Processing(format!("STAC search failed: {}", e)))?;

        if !response.status().is_success() {
            return Err(AcoliteError::Processing(format!(
                "STAC HTTP {}",
                response.status()
            )));
        }

        let result: SearchResponse = response
            .json()
            .map_err(|e| AcoliteError::Processing(format!("STAC parse failed: {}", e)))?;

        log::info!(
            "STAC: found {} items in {}",
            result.features.len(),
            collection
        );
        Ok(result.features)
    }

    /// Search Landsat Collection 2 Level-1 via USGS STAC
    pub fn search_landsat_c2(
        &self,
        bbox: &[f64; 4],
        datetime: &str,
        limit: usize,
    ) -> Result<Vec<StacItem>> {
        self.search("landsat-c2l1", bbox, datetime, limit)
    }
}

/// Download a Landsat scene's band assets from a STAC item to a local directory.
/// Returns the local directory path containing the downloaded files.
pub fn download_stac_landsat(
    item: &StacItem,
    output_dir: &std::path::Path,
    band_keys: &[&str],
) -> Result<std::path::PathBuf> {
    let scene_dir = output_dir.join(&item.id);
    std::fs::create_dir_all(&scene_dir).map_err(AcoliteError::Io)?;

    let client = reqwest::blocking::Client::builder()
        .redirect(reqwest::redirect::Policy::limited(10))
        .build()
        .map_err(|e| AcoliteError::Processing(format!("HTTP client: {}", e)))?;

    for key in band_keys {
        let asset = match item.assets.get(*key) {
            Some(a) => a,
            None => {
                log::warn!("Asset '{}' not found in STAC item {}", key, item.id);
                continue;
            }
        };

        let filename = asset.href.rsplit('/').next().unwrap_or(key);
        let local_path = scene_dir.join(filename);

        if local_path.exists() {
            log::info!("Cached: {:?}", local_path);
            continue;
        }

        log::info!("Downloading {} → {:?}", key, local_path);
        let response = client
            .get(&asset.href)
            .send()
            .map_err(|e| AcoliteError::Processing(format!("Download {}: {}", key, e)))?;

        if !response.status().is_success() {
            return Err(AcoliteError::Processing(format!(
                "HTTP {} for asset {}",
                response.status(),
                key
            )));
        }

        let bytes = response
            .bytes()
            .map_err(|e| AcoliteError::Processing(format!("Read {}: {}", key, e)))?;
        std::fs::write(&local_path, bytes).map_err(AcoliteError::Io)?;
    }

    Ok(scene_dir)
}
