//! Sentinel-3 OLCI sensor implementation

use crate::core::Metadata;
use crate::sensors::Sensor;
use crate::error::Result;

pub struct Sentinel3Sensor;

impl Sensor for Sentinel3Sensor {
    fn name(&self) -> &str {
        "S3_OLCI"
    }
    
    fn parse_metadata(&self, _path: &str) -> Result<Metadata> {
        // Placeholder - would parse xfdumanifest.xml
        Ok(Metadata::new(self.name().to_string(), chrono::Utc::now()))
    }
    
    fn band_names(&self) -> Vec<String> {
        (1..=21).map(|i| format!("Oa{:02}", i)).collect()
    }
    
    fn wavelength(&self, band: &str) -> Option<f64> {
        let wavelengths = [
            400.0, 412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 665.0, 673.75, 681.25,
            708.75, 753.75, 761.25, 764.375, 767.5, 778.75, 865.0, 885.0, 900.0, 940.0, 1020.0
        ];
        
        if let Some(idx) = band.strip_prefix("Oa") {
            if let Ok(n) = idx.parse::<usize>() {
                if n >= 1 && n <= 21 {
                    return Some(wavelengths[n - 1]);
                }
            }
        }
        None
    }
    
    fn bandwidth(&self, band: &str) -> Option<f64> {
        let bandwidths = [
            15.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 7.5, 7.5,
            10.0, 7.5, 2.5, 3.75, 2.5, 15.0, 20.0, 10.0, 10.0, 20.0, 40.0
        ];
        
        if let Some(idx) = band.strip_prefix("Oa") {
            if let Ok(n) = idx.parse::<usize>() {
                if n >= 1 && n <= 21 {
                    return Some(bandwidths[n - 1]);
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_sentinel3_bands() {
        let sensor = Sentinel3Sensor;
        assert_eq!(sensor.name(), "S3_OLCI");
        assert_eq!(sensor.band_names().len(), 21);
        assert_eq!(sensor.wavelength("Oa01"), Some(400.0));
        assert_eq!(sensor.wavelength("Oa21"), Some(1020.0));
        assert_eq!(sensor.bandwidth("Oa01"), Some(15.0));
    }
}
