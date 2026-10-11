//! Sentinel-2 sensor implementation

use crate::{core::Metadata, sensors::Sensor, AcoliteError, Result};
use std::collections::HashMap;

pub struct Sentinel2Sensor {
    name: String,
    bands: HashMap<String, (f64, f64, u32)>, // wavelength, bandwidth, resolution
}

impl Sentinel2Sensor {
    pub fn new_s2a() -> Self {
        let mut bands = HashMap::new();
        // 10m bands
        bands.insert("B02".to_string(), (492.4, 66.0, 10)); // Blue
        bands.insert("B03".to_string(), (559.8, 36.0, 10)); // Green
        bands.insert("B04".to_string(), (664.6, 31.0, 10)); // Red
        bands.insert("B08".to_string(), (832.8, 106.0, 10)); // NIR

        // 20m bands
        bands.insert("B05".to_string(), (704.1, 15.0, 20)); // Red Edge 1
        bands.insert("B06".to_string(), (740.5, 15.0, 20)); // Red Edge 2
        bands.insert("B07".to_string(), (782.8, 20.0, 20)); // Red Edge 3
        bands.insert("B8A".to_string(), (864.7, 21.0, 20)); // NIR narrow
        bands.insert("B11".to_string(), (1613.7, 91.0, 20)); // SWIR 1
        bands.insert("B12".to_string(), (2202.4, 175.0, 20)); // SWIR 2

        // 60m bands
        bands.insert("B01".to_string(), (442.7, 21.0, 60)); // Coastal
        bands.insert("B09".to_string(), (945.1, 20.0, 60)); // Water vapor
        bands.insert("B10".to_string(), (1373.5, 31.0, 60)); // Cirrus

        Self {
            name: "S2A_MSI".to_string(),
            bands,
        }
    }

    pub fn new_s2b() -> Self {
        let mut sensor = Self::new_s2a();
        sensor.name = "S2B_MSI".to_string();
        // S2B has slightly different wavelengths but we'll use S2A for now
        sensor
    }

    /// Get band resolution
    pub fn resolution(&self, band: &str) -> Option<u32> {
        self.bands.get(band).map(|(_, _, r)| *r)
    }

    /// Get bands at specific resolution
    pub fn bands_at_resolution(&self, resolution: u32) -> Vec<String> {
        self.bands
            .iter()
            .filter(|(_, (_, _, r))| *r == resolution)
            .map(|(name, _)| name.clone())
            .collect()
    }
}

impl Sensor for Sentinel2Sensor {
    fn name(&self) -> &str {
        &self.name
    }

    fn parse_metadata(&self, path: &str) -> Result<Metadata> {
        crate::sensors::s2_xml::parse_s2_metadata(path)
    }

    fn band_names(&self) -> Vec<String> {
        self.bands.keys().cloned().collect()
    }

    fn wavelength(&self, band: &str) -> Option<f64> {
        self.bands.get(band).map(|(w, _, _)| *w)
    }

    fn bandwidth(&self, band: &str) -> Option<f64> {
        self.bands.get(band).map(|(_, b, _)| *b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_s2a_bands() {
        let s2a = Sentinel2Sensor::new_s2a();
        assert_eq!(s2a.name(), "S2A_MSI");

        let bands = s2a.band_names();
        assert!(bands.len() >= 13);

        // Check 10m bands
        let bands_10m = s2a.bands_at_resolution(10);
        assert_eq!(bands_10m.len(), 4);

        // Check wavelengths
        assert_eq!(s2a.wavelength("B04"), Some(664.6));
        assert_eq!(s2a.resolution("B04"), Some(10));
    }
}
