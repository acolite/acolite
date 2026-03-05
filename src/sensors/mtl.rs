//! Landsat MTL (Metadata) file parser

use crate::{Result, AcoliteError, core::Metadata};
use chrono::{DateTime, Utc, NaiveDateTime};
use std::collections::HashMap;
use std::fs;

/// Parse Landsat MTL file
pub fn parse_mtl(path: &str) -> Result<Metadata> {
    let content = fs::read_to_string(path)
        .map_err(|e| AcoliteError::Io(e))?;
    
    let mut attrs = parse_mtl_content(&content)?;
    
    // Extract sensor name
    let sensor = attrs.get("SPACECRAFT_ID")
        .and_then(|s| attrs.get("SENSOR_ID").map(|sen| format!("{}_{}", s, sen)))
        .ok_or_else(|| AcoliteError::InvalidMetadata("Missing sensor ID".to_string()))?;
    
    // Parse datetime
    let date = attrs.get("DATE_ACQUIRED")
        .ok_or_else(|| AcoliteError::InvalidMetadata("Missing DATE_ACQUIRED".to_string()))?;
    let time = attrs.get("SCENE_CENTER_TIME")
        .ok_or_else(|| AcoliteError::InvalidMetadata("Missing SCENE_CENTER_TIME".to_string()))?;
    
    let datetime = parse_datetime(date, time)?;
    
    let mut metadata = Metadata::new(sensor, datetime);
    
    // Parse geometry
    if let (Some(sza), Some(saa)) = (attrs.get("SUN_ELEVATION"), attrs.get("SUN_AZIMUTH")) {
        let sun_zenith = 90.0 - sza.parse::<f64>()
            .map_err(|_| AcoliteError::InvalidMetadata("Invalid SUN_ELEVATION".to_string()))?;
        let sun_azimuth = saa.parse::<f64>()
            .map_err(|_| AcoliteError::InvalidMetadata("Invalid SUN_AZIMUTH".to_string()))?;
        metadata.set_geometry(sun_zenith, sun_azimuth);
    }
    
    metadata.attributes = attrs;
    Ok(metadata)
}

fn parse_mtl_content(content: &str) -> Result<HashMap<String, String>> {
    let mut attrs = HashMap::new();
    
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with("GROUP") || line.starts_with("END") {
            continue;
        }
        
        if let Some((key, value)) = line.split_once('=') {
            let key = key.trim().to_string();
            let value = value.trim().trim_matches('"').to_string();
            attrs.insert(key, value);
        }
    }
    
    Ok(attrs)
}

fn parse_datetime(date: &str, time: &str) -> Result<DateTime<Utc>> {
    let time_clean = time.trim_end_matches('Z');
    let datetime_str = format!("{} {}", date, time_clean);
    
    NaiveDateTime::parse_from_str(&datetime_str, "%Y-%m-%d %H:%M:%S%.f")
        .map(|dt| DateTime::from_naive_utc_and_offset(dt, Utc))
        .map_err(|_| AcoliteError::InvalidMetadata("Invalid datetime format".to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse_datetime() {
        let dt = parse_datetime("2023-06-15", "10:30:45.123Z").unwrap();
        assert_eq!(dt.format("%Y-%m-%d").to_string(), "2023-06-15");
    }
    
    #[test]
    fn test_parse_mtl_content() {
        let content = r#"
            GROUP = L1_METADATA_FILE
              SPACECRAFT_ID = "LANDSAT_8"
              SENSOR_ID = "OLI_TIRS"
              DATE_ACQUIRED = "2023-06-15"
            END_GROUP = L1_METADATA_FILE
        "#;
        
        let attrs = parse_mtl_content(content).unwrap();
        assert_eq!(attrs.get("SPACECRAFT_ID").unwrap(), "LANDSAT_8");
        assert_eq!(attrs.get("SENSOR_ID").unwrap(), "OLI_TIRS");
    }
}
