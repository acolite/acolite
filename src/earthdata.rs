//! NASA Earthdata authentication and data access

use crate::{Result, AcoliteError};
use std::fs;
use std::path::PathBuf;
use std::io::Write;
use futures_util::StreamExt;

/// Earthdata credentials
#[derive(Debug, Clone)]
pub struct EarthdataAuth {
    pub username: String,
    pub password: String,
    pub token: Option<String>,
}

impl EarthdataAuth {
    /// Load credentials from ~/.easi-workflows-auth.conf
    pub fn load() -> Result<Self> {
        let home = dirs::home_dir()
            .ok_or_else(|| AcoliteError::Processing("Home directory not found".to_string()))?;
        
        let auth_path = home.join(".easi-workflows-auth.conf");
        
        if !auth_path.exists() {
            return Err(AcoliteError::Processing(
                format!("Auth file not found: {:?}\nCreate with:\n[earthdata]\n    user: YOUR_USERNAME\n    password: YOUR_PASSWORD", auth_path)
            ));
        }
        
        let content = fs::read_to_string(&auth_path)
            .map_err(|e| AcoliteError::Io(e))?;
        
        let mut username = None;
        let mut password = None;
        let mut token = None;
        let mut in_earthdata_section = false;
        
        for line in content.lines() {
            let line = line.trim();
            
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            
            // Check for [earthdata] section
            if line == "[earthdata]" {
                in_earthdata_section = true;
                continue;
            }
            
            // Check for other sections
            if line.starts_with('[') && line.ends_with(']') {
                in_earthdata_section = false;
                continue;
            }
            
            // Parse key: value pairs in earthdata section
            if in_earthdata_section {
                if let Some((key, value)) = line.split_once(':') {
                    let key = key.trim();
                    let value = value.trim();
                    
                    match key {
                        "user" => username = Some(value.to_string()),
                        "password" => password = Some(value.to_string()),
                        "token" => token = Some(value.to_string()),
                        _ => {}
                    }
                }
            }
        }
        
        let username = username.ok_or_else(|| 
            AcoliteError::Processing("user not found in [earthdata] section".to_string()))?;
        let password = password.ok_or_else(|| 
            AcoliteError::Processing("password not found in [earthdata] section".to_string()))?;
        
        Ok(Self { username, password, token })
    }
}

/// Search PACE data from NASA CMR
pub fn search_pace_data(
    bbox: &[f64; 4],
    start_date: &str,
    end_date: &str,
) -> Result<Vec<String>> {
    let client = reqwest::blocking::Client::new();
    
    let cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.json";
    
    let params = [
        ("short_name", "PACE_OCI_L1B_SCI"),
        ("bounding_box", &format!("{},{},{},{}", bbox[0], bbox[1], bbox[2], bbox[3])),
        ("temporal", &format!("{},{}", start_date, end_date)),
        ("page_size", "50"),
    ];
    
    log::info!("Searching PACE data in CMR...");
    log::debug!("CMR URL: {}", cmr_url);
    log::debug!("Params: {:?}", params);
    log::debug!("Bbox: {:?}", bbox);
    
    let response = client
        .get(cmr_url)
        .query(&params)
        .send()
        .map_err(|e| AcoliteError::Processing(format!("CMR search failed: {}", e)))?;
    
    if !response.status().is_success() {
        return Err(AcoliteError::Processing(format!("HTTP {}", response.status())));
    }
    
    let json: serde_json::Value = response.json()
        .map_err(|e| AcoliteError::Processing(format!("Parse failed: {}", e)))?;
    
    let mut urls = Vec::new();
    
    if let Some(entries) = json["feed"]["entry"].as_array() {
        for entry in entries {
            if let Some(links) = entry["links"].as_array() {
                for link in links {
                    if link["rel"].as_str() == Some("http://esipfed.org/ns/fedsearch/1.1/data#") {
                        if let Some(href) = link["href"].as_str() {
                            if href.ends_with(".nc") {
                                urls.push(href.to_string());
                            }
                        }
                    }
                }
            }
        }
    }
    
    log::info!("Found {} PACE granules", urls.len());
    Ok(urls)
}

/// Get temporary S3 credentials from TEA endpoint
pub fn get_s3_credentials(auth: &EarthdataAuth) -> Result<S3Credentials> {
    let client = reqwest::blocking::Client::new();
    let tea_url = "https://obdaac-tea.earthdatacloud.nasa.gov/s3credentials";
    
    log::info!("Getting S3 credentials from TEA...");
    
    let response = if let Some(token) = &auth.token {
        client
            .get(tea_url)
            .bearer_auth(token)
            .send()
            .map_err(|e| AcoliteError::Processing(format!("TEA request failed: {}", e)))?
    } else {
        return Err(AcoliteError::Processing("Token required for S3 access".to_string()));
    };
    
    if !response.status().is_success() {
        return Err(AcoliteError::Processing(format!("TEA HTTP {}", response.status())));
    }
    
    let creds: S3Credentials = response.json()
        .map_err(|e| AcoliteError::Processing(format!("Parse S3 creds failed: {}", e)))?;
    
    log::info!("Got S3 credentials (expires: {})", creds.expiration);
    Ok(creds)
}

/// S3 temporary credentials from TEA
#[derive(Debug, Clone, serde::Deserialize)]
pub struct S3Credentials {
    #[serde(rename = "accessKeyId")]
    pub access_key_id: String,
    #[serde(rename = "secretAccessKey")]
    pub secret_access_key: String,
    #[serde(rename = "sessionToken")]
    pub session_token: String,
    pub expiration: String,
}

/// Download file from S3 via TEA with streaming (async, no timeout)
pub async fn download_s3_streaming(
    s3_url: &str,
    output_path: &str,
    auth: &EarthdataAuth,
) -> Result<()> {
    let s3_path = s3_url.strip_prefix("s3://")
        .ok_or_else(|| AcoliteError::Processing("Invalid S3 URL".to_string()))?;
    
    let https_url = format!("https://obdaac-tea.earthdatacloud.nasa.gov/{}", s3_path);
    
    let token = auth.token.as_ref()
        .ok_or_else(|| AcoliteError::Processing("Token required".to_string()))?;
    
    log::info!("Streaming download: {}", s3_path);
    
    let client = reqwest::Client::builder()
        .build()
        .map_err(|e| AcoliteError::Processing(format!("Client failed: {}", e)))?;
    
    let response = client
        .get(&https_url)
        .bearer_auth(token)
        .send()
        .await
        .map_err(|e| AcoliteError::Processing(format!("Request failed: {}", e)))?;
    
    if !response.status().is_success() {
        return Err(AcoliteError::Processing(format!("HTTP {}", response.status())));
    }
    
    let total_size = response.content_length().unwrap_or(0);
    log::info!("File size: {} MB", total_size / 1_000_000);
    
    let mut file = fs::File::create(output_path)
        .map_err(|e| AcoliteError::Io(e))?;
    
    let mut stream = response.bytes_stream();
    let mut downloaded = 0u64;
    let mut last_progress = 0;
    
    while let Some(chunk) = stream.next().await {
        let chunk = chunk.map_err(|e| AcoliteError::Processing(format!("Stream error: {}", e)))?;
        file.write_all(&chunk).map_err(|e| AcoliteError::Io(e))?;
        
        downloaded += chunk.len() as u64;
        
        if total_size > 0 {
            let progress = (downloaded * 100 / total_size) as u32;
            if progress >= last_progress + 10 {
                log::info!("Progress: {}% ({} MB)", progress, downloaded / 1_000_000);
                last_progress = progress;
            }
        }
    }
    
    log::info!("Download complete: {} MB", downloaded / 1_000_000);
    Ok(())
}

/// Blocking wrapper for streaming download
pub fn download_from_s3_blocking(
    s3_url: &str,
    output_path: &str,
    auth: &EarthdataAuth,
) -> Result<()> {
    tokio::runtime::Runtime::new()
        .map_err(|e| AcoliteError::Processing(format!("Tokio runtime failed: {}", e)))?
        .block_on(download_s3_streaming(s3_url, output_path, auth))
}

/// Download file from S3 using TEA proxy (handles auth automatically)
pub fn download_from_s3_with_token(
    s3_url: &str,
    output_path: &str,
    auth: &EarthdataAuth,
) -> Result<()> {
    download_from_s3_blocking(s3_url, output_path, auth)
}

/// Download file from S3 using temporary credentials
pub fn download_from_s3(
    s3_url: &str,
    output_path: &str,
    creds: &S3Credentials,
) -> Result<()> {
    // Deprecated - use download_from_s3_with_token instead
    Err(AcoliteError::Processing("Use download_from_s3_with_token".to_string()))
}

/// Download file from NASA DAAC with authentication
pub fn download_pace_file(
    url: &str,
    output_path: &str,
    auth: &EarthdataAuth,
) -> Result<()> {
    log::info!("Downloading: {}", url);
    
    let client = reqwest::blocking::Client::builder()
        .redirect(reqwest::redirect::Policy::limited(10))
        .build()
        .map_err(|e| AcoliteError::Processing(format!("Client build failed: {}", e)))?;
    
    // Use token if available, otherwise basic auth
    let response = if let Some(token) = &auth.token {
        log::info!("Using bearer token authentication");
        client
            .get(url)
            .bearer_auth(token)
            .send()
            .map_err(|e| AcoliteError::Processing(format!("Download failed: {}", e)))?
    } else {
        log::info!("Using basic authentication");
        client
            .get(url)
            .basic_auth(&auth.username, Some(&auth.password))
            .send()
            .map_err(|e| AcoliteError::Processing(format!("Download failed: {}", e)))?
    };
    
    if !response.status().is_success() {
        return Err(AcoliteError::Processing(
            format!("HTTP {} - Check credentials and token", response.status())
        ));
    }
    
    let bytes = response.bytes()
        .map_err(|e| AcoliteError::Processing(format!("Read failed: {}", e)))?;
    
    fs::write(output_path, bytes)
        .map_err(|e| AcoliteError::Io(e))?;
    
    let size = fs::metadata(output_path)
        .map_err(|e| AcoliteError::Io(e))?
        .len();
    
    log::info!("Downloaded {} MB", size / 1_000_000);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_auth_load() {
        // This will fail if auth file doesn't exist, which is expected
        match EarthdataAuth::load() {
            Ok(auth) => {
                assert!(!auth.username.is_empty());
                assert!(!auth.password.is_empty());
            }
            Err(_) => {
                // Expected if no auth file
            }
        }
    }
    
    #[test]
    #[ignore] // Requires network
    fn test_search_pace() {
        let bbox = [-80.0, 30.0, -70.0, 40.0]; // Atlantic
        let result = search_pace_data(&bbox, "2024-03-01", "2024-03-31");
        
        if let Ok(urls) = result {
            println!("Found {} PACE files", urls.len());
        }
    }
}
