//! File download — single implementation supporting EarthData and AWS

use crate::auth::Credentials;
use crate::{AcoliteError, Result};
use std::path::Path;

/// Download a file with EarthData authentication
pub fn download_file(url: &str, output: &Path, creds: &Credentials) -> Result<()> {
    if output.exists() {
        log::info!("Cached: {:?}", output);
        return Ok(());
    }

    log::info!("Downloading: {}", url);

    let client = reqwest::blocking::Client::builder()
        .redirect(reqwest::redirect::Policy::limited(10))
        .build()
        .map_err(|e| AcoliteError::Processing(format!("HTTP client: {}", e)))?;

    let request = if let Some(token) = creds.token() {
        client.get(url).bearer_auth(token)
    } else {
        client
            .get(url)
            .basic_auth(&creds.username, Some(creds.password()))
    };

    let response = request
        .send()
        .map_err(|e| AcoliteError::Processing(format!("Download failed: {}", e)))?;

    if !response.status().is_success() {
        return Err(AcoliteError::Processing(format!(
            "HTTP {}",
            response.status()
        )));
    }

    let bytes = response
        .bytes()
        .map_err(|e| AcoliteError::Processing(format!("Read failed: {}", e)))?;

    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(AcoliteError::Io)?;
    }

    std::fs::write(output, bytes).map_err(AcoliteError::Io)?;

    let size = std::fs::metadata(output).map_err(AcoliteError::Io)?.len();
    log::info!("Downloaded {:.1} MB → {:?}", size as f64 / 1e6, output);
    Ok(())
}

/// Download from AWS S3 using CLI (for Requester Pays buckets)
pub fn download_s3_requester_pays(s3_url: &str, output: &Path, profile: &str) -> Result<()> {
    if output.exists() {
        log::info!("Cached: {:?}", output);
        return Ok(());
    }

    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(AcoliteError::Io)?;
    }

    let status = std::process::Command::new("aws")
        .args([
            "s3",
            "cp",
            s3_url,
            output.to_str().unwrap_or_default(),
            "--request-payer",
            "requester",
            "--profile",
            profile,
        ])
        .status()
        .map_err(|e| AcoliteError::Processing(format!("aws cli: {}", e)))?;

    if !status.success() {
        return Err(AcoliteError::Processing(format!(
            "aws s3 cp failed for {}",
            s3_url
        )));
    }

    let size = std::fs::metadata(output).map_err(AcoliteError::Io)?.len();
    log::info!("Downloaded {:.1} MB → {:?}", size as f64 / 1e6, output);
    Ok(())
}
