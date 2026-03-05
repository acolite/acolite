//! STAC API client for satellite data discovery and download

use crate::{Result, AcoliteError};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// STAC Item representing a satellite scene
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct StacItem {
    pub id: String,
    pub bbox: Vec<f64>,
    pub geometry: serde_json::Value,
    pub properties: HashMap<String, serde_json::Value>,
    pub assets: HashMap<String, StacAsset>,
}

/// STAC Asset representing a downloadable file
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct StacAsset {
    pub href: String,
    #[serde(rename = "type")]
    pub media_type: Option<String>,
    pub title: Option<String>,
    pub roles: Option<Vec<String>>,
}

/// STAC search response
#[derive(Debug, Deserialize)]
struct StacSearchResponse {
    features: Vec<StacItem>,
}

/// STAC API client
pub struct StacClient {
    base_url: String,
    client: reqwest::blocking::Client,
}

impl StacClient {
    /// Create new STAC client
    pub fn new(base_url: String) -> Self {
        Self {
            base_url,
            client: reqwest::blocking::Client::new(),
        }
    }
    
    /// Create client for Microsoft Planetary Computer
    pub fn planetary_computer() -> Self {
        Self::new("https://planetarycomputer.microsoft.com/api/stac/v1".to_string())
    }
    
    /// Create client for Earth Search (AWS)
    pub fn earth_search() -> Self {
        Self::new("https://earth-search.aws.element84.com/v1".to_string())
    }
    
    /// Search for items
    pub fn search(
        &self,
        collection: &str,
        bbox: &[f64; 4],
        datetime: &str,
        limit: usize,
    ) -> Result<Vec<StacItem>> {
        let search_url = format!("{}/search", self.base_url);
        
        let body = serde_json::json!({
            "collections": [collection],
            "bbox": bbox,
            "datetime": datetime,
            "limit": limit,
        });
        
        log::info!("STAC search: {} in {}", collection, datetime);
        
        let response = self.client
            .post(&search_url)
            .json(&body)
            .send()
            .map_err(|e| AcoliteError::Processing(format!("STAC search failed: {}", e)))?;
        
        if !response.status().is_success() {
            return Err(AcoliteError::Processing(format!("HTTP {}", response.status())));
        }
        
        let search_response: StacSearchResponse = response.json()
            .map_err(|e| AcoliteError::Processing(format!("Parse failed: {}", e)))?;
        
        log::info!("Found {} items", search_response.features.len());
        Ok(search_response.features)
    }
    
    /// Download asset to file
    pub fn download_asset(&self, asset: &StacAsset, output_path: &str) -> Result<()> {
        log::info!("Downloading: {}", asset.href);
        
        let response = self.client
            .get(&asset.href)
            .send()
            .map_err(|e| AcoliteError::Processing(format!("Download failed: {}", e)))?;
        
        if !response.status().is_success() {
            return Err(AcoliteError::Processing(format!("HTTP {}", response.status())));
        }
        
        let bytes = response.bytes()
            .map_err(|e| AcoliteError::Processing(format!("Read failed: {}", e)))?;
        
        std::fs::write(output_path, bytes)
            .map_err(|e| AcoliteError::Io(e))?;
        
        log::info!("Downloaded {} bytes", std::fs::metadata(output_path).unwrap().len());
        Ok(())
    }
}

/// Helper to find Landsat 8/9 items
pub fn search_landsat(
    bbox: &[f64; 4],
    datetime: &str,
    max_cloud: f64,
) -> Result<Vec<StacItem>> {
    let client = StacClient::earth_search();
    let items = client.search("landsat-c2-l1", bbox, datetime, 10)?;
    
    // Filter by cloud cover
    Ok(items.into_iter()
        .filter(|item| {
            item.properties.get("eo:cloud_cover")
                .and_then(|v| v.as_f64())
                .map(|cc| cc <= max_cloud)
                .unwrap_or(false)
        })
        .collect())
}

/// Helper to find Sentinel-2 items
pub fn search_sentinel2(
    bbox: &[f64; 4],
    datetime: &str,
    max_cloud: f64,
) -> Result<Vec<StacItem>> {
    let client = StacClient::earth_search();
    let items = client.search("sentinel-2-l1c", bbox, datetime, 10)?;
    
    // Filter by cloud cover
    Ok(items.into_iter()
        .filter(|item| {
            item.properties.get("eo:cloud_cover")
                .and_then(|v| v.as_f64())
                .map(|cc| cc <= max_cloud)
                .unwrap_or(false)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_stac_client_creation() {
        let client = StacClient::planetary_computer();
        assert!(client.base_url.contains("planetarycomputer"));
        
        let client = StacClient::earth_search();
        assert!(client.base_url.contains("earth-search"));
    }
    
    #[test]
    #[ignore] // Requires network
    fn test_stac_search() {
        let client = StacClient::earth_search();
        let bbox = [-122.5, 37.5, -122.0, 38.0]; // San Francisco Bay
        let items = client.search("landsat-c2-l1", &bbox, "2023-01-01/2023-12-31", 5);
        
        if let Ok(items) = items {
            assert!(!items.is_empty());
            println!("Found {} Landsat scenes", items.len());
        }
    }
}
