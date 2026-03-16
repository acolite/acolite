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

    /// Search Sentinel-2 L1C via Element84 Earth Search
    pub fn search_sentinel2_l1c(
        &self,
        bbox: &[f64; 4],
        datetime: &str,
        limit: usize,
    ) -> Result<Vec<StacItem>> {
        self.search("sentinel-2-l1c", bbox, datetime, limit)
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

/// Download a Sentinel-2 SAFE directory from CDSE using OAuth2 bearer token.
/// `url` is the OData download URL (ending in `/$value`).
/// Returns the path to the extracted .SAFE directory.
pub fn download_cdse_sentinel2(
    url: &str,
    output_dir: &std::path::Path,
    creds: &crate::auth::CdseCredentials,
) -> Result<std::path::PathBuf> {
    let token = creds.access_token()?;

    let client = reqwest::blocking::Client::builder()
        .redirect(reqwest::redirect::Policy::limited(10))
        .build()
        .map_err(|e| AcoliteError::Processing(format!("HTTP client: {}", e)))?;

    // Discover scene name from OData metadata
    let meta_url = url.trim_end_matches("/$value").trim_end_matches("$value");
    let meta_resp = client
        .get(meta_url)
        .send()
        .ok()
        .and_then(|r| r.json::<serde_json::Value>().ok());
    let scene_name = meta_resp
        .as_ref()
        .and_then(|j| j["Name"].as_str().or_else(|| j["value"][0]["Name"].as_str()))
        .unwrap_or("S2_scene")
        .to_string();

    let safe_dir = output_dir.join(&scene_name);
    if safe_dir.exists() {
        log::info!("CDSE scene cached: {:?}", safe_dir);
        return Ok(safe_dir);
    }

    let zip_path = output_dir.join(format!(
        "{}.zip",
        scene_name.trim_end_matches(".SAFE")
    ));
    std::fs::create_dir_all(output_dir).map_err(AcoliteError::Io)?;

    log::info!("Downloading CDSE scene: {}", scene_name);
    let resp = client
        .get(url)
        .bearer_auth(&token)
        .send()
        .map_err(|e| AcoliteError::Processing(format!("CDSE download: {}", e)))?;

    if !resp.status().is_success() {
        return Err(AcoliteError::Processing(format!(
            "CDSE HTTP {}",
            resp.status()
        )));
    }

    let bytes = resp
        .bytes()
        .map_err(|e| AcoliteError::Processing(format!("CDSE read: {}", e)))?;
    std::fs::write(&zip_path, &bytes).map_err(AcoliteError::Io)?;
    log::info!("Downloaded {:.1} MB", bytes.len() as f64 / 1e6);

    // Extract zip
    let file = std::fs::File::open(&zip_path).map_err(AcoliteError::Io)?;
    let mut archive =
        zip::ZipArchive::new(file).map_err(|e| AcoliteError::Processing(format!("Zip: {}", e)))?;
    for i in 0..archive.len() {
        let mut entry = archive
            .by_index(i)
            .map_err(|e| AcoliteError::Processing(format!("Zip entry: {}", e)))?;
        let out_path = output_dir.join(entry.mangled_name());
        if entry.is_dir() {
            std::fs::create_dir_all(&out_path).map_err(AcoliteError::Io)?;
        } else {
            if let Some(p) = out_path.parent() {
                std::fs::create_dir_all(p).map_err(AcoliteError::Io)?;
            }
            let mut outfile = std::fs::File::create(&out_path).map_err(AcoliteError::Io)?;
            std::io::copy(&mut entry, &mut outfile).map_err(AcoliteError::Io)?;
        }
    }
    let _ = std::fs::remove_file(&zip_path);
    log::info!("Extracted → {:?}", safe_dir);

    Ok(safe_dir)
}

/// Download a Sentinel-2 scene from a STAC item (Element84 Earth Search).
/// Downloads all JP2 band files + metadata XML into a reconstructed SAFE-like directory.
pub fn download_stac_sentinel2(
    item: &StacItem,
    output_dir: &std::path::Path,
) -> Result<std::path::PathBuf> {
    let scene_dir = output_dir.join(&item.id);
    let img_dir = scene_dir.join("IMG_DATA");
    std::fs::create_dir_all(&img_dir).map_err(AcoliteError::Io)?;

    let client = reqwest::blocking::Client::builder()
        .redirect(reqwest::redirect::Policy::limited(10))
        .build()
        .map_err(|e| AcoliteError::Processing(format!("HTTP client: {}", e)))?;

    // Download band assets (keys like "B01", "B02", ..., "B8A", etc.) and metadata
    let band_keys: Vec<&str> = vec![
        "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11",
        "B12",
    ];
    let meta_keys: Vec<&str> = vec!["granule_metadata", "tileinfo_metadata", "product_metadata"];

    for key in band_keys.iter().chain(meta_keys.iter()) {
        let asset = match item.assets.get(*key) {
            Some(a) => a,
            None => continue,
        };
        let filename = asset.href.rsplit('/').next().unwrap_or(key);
        let dest = if filename.ends_with(".jp2") {
            img_dir.join(filename)
        } else {
            scene_dir.join(filename)
        };
        if dest.exists() {
            continue;
        }
        log::info!("Downloading S2 asset {} → {:?}", key, dest);
        let resp = client
            .get(&asset.href)
            .send()
            .map_err(|e| AcoliteError::Processing(format!("Download {}: {}", key, e)))?;
        if !resp.status().is_success() {
            log::warn!("HTTP {} for S2 asset {}", resp.status(), key);
            continue;
        }
        let bytes = resp
            .bytes()
            .map_err(|e| AcoliteError::Processing(format!("Read {}: {}", key, e)))?;
        std::fs::write(&dest, bytes).map_err(AcoliteError::Io)?;
    }

    Ok(scene_dir)
}
