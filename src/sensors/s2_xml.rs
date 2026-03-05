//! Sentinel-2 XML metadata parser

use crate::{Result, AcoliteError, core::Metadata};
use chrono::{DateTime, Utc, NaiveDateTime};
use std::fs;

/// Parse Sentinel-2 MTD_MSIL1C.xml metadata
pub fn parse_s2_metadata(path: &str) -> Result<Metadata> {
    let content = fs::read_to_string(path)
        .map_err(|e| AcoliteError::Io(e))?;
    
    // Simple XML parsing (in production, use proper XML parser)
    let sensor = extract_tag(&content, "SPACECRAFT_NAME")
        .unwrap_or_else(|| "S2A_MSI".to_string());
    
    let datetime_str = extract_tag(&content, "PRODUCT_START_TIME")
        .or_else(|| extract_tag(&content, "DATATAKE_SENSING_START"))
        .ok_or_else(|| AcoliteError::InvalidMetadata("Missing datetime".to_string()))?;
    
    let datetime = parse_s2_datetime(&datetime_str)?;
    
    let mut metadata = Metadata::new(sensor, datetime);
    
    // Extract sun angles
    if let Some(zenith) = extract_tag(&content, "ZENITH_ANGLE") {
        if let Ok(z) = zenith.parse::<f64>() {
            if let Some(azimuth) = extract_tag(&content, "AZIMUTH_ANGLE") {
                if let Ok(a) = azimuth.parse::<f64>() {
                    metadata.set_geometry(z, a);
                }
            }
        }
    }
    
    Ok(metadata)
}

fn extract_tag(content: &str, tag: &str) -> Option<String> {
    let start_tag = format!("<{}>", tag);
    let end_tag = format!("</{}>", tag);
    
    if let Some(start) = content.find(&start_tag) {
        let start_pos = start + start_tag.len();
        if let Some(end) = content[start_pos..].find(&end_tag) {
            return Some(content[start_pos..start_pos + end].trim().to_string());
        }
    }
    None
}

fn parse_s2_datetime(s: &str) -> Result<DateTime<Utc>> {
    // S2 format: 2023-06-15T10:30:45.123Z
    let s_clean = s.trim_end_matches('Z');
    
    NaiveDateTime::parse_from_str(s_clean, "%Y-%m-%dT%H:%M:%S%.f")
        .or_else(|_| NaiveDateTime::parse_from_str(s_clean, "%Y-%m-%dT%H:%M:%S"))
        .map(|dt| DateTime::from_naive_utc_and_offset(dt, Utc))
        .map_err(|_| AcoliteError::InvalidMetadata("Invalid datetime format".to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse_s2_datetime() {
        let dt = parse_s2_datetime("2023-06-15T10:30:45.123Z").unwrap();
        assert_eq!(dt.format("%Y-%m-%d").to_string(), "2023-06-15");
    }
    
    #[test]
    fn test_extract_tag() {
        let xml = "<SPACECRAFT_NAME>Sentinel-2A</SPACECRAFT_NAME>";
        assert_eq!(extract_tag(xml, "SPACECRAFT_NAME"), Some("Sentinel-2A".to_string()));
    }
}
