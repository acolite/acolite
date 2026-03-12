//! STAC API client — single implementation

use crate::{Result, AcoliteError};
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
            return Err(AcoliteError::Processing(format!("STAC HTTP {}", response.status())));
        }

        let result: SearchResponse = response
            .json()
            .map_err(|e| AcoliteError::Processing(format!("STAC parse failed: {}", e)))?;

        log::info!("STAC: found {} items in {}", result.features.len(), collection);
        Ok(result.features)
    }
}
