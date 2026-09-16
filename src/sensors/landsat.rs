//! Landsat sensor implementation

use crate::{
    core::Metadata,
    sensors::{parse_mtl, Sensor},
    AcoliteError, Result,
};
use std::collections::HashMap;

pub struct LandsatSensor {
    name: String,
    bands: HashMap<String, (f64, f64)>, // wavelength, bandwidth
}

impl LandsatSensor {
    pub fn new_l8() -> Self {
        let mut bands = HashMap::new();
        bands.insert("B1".to_string(), (443.0, 16.0));
        bands.insert("B2".to_string(), (482.0, 60.0));
        bands.insert("B3".to_string(), (561.0, 57.0));
        bands.insert("B4".to_string(), (655.0, 37.0));
        bands.insert("B5".to_string(), (865.0, 28.0));
        bands.insert("B6".to_string(), (1609.0, 85.0));
        bands.insert("B7".to_string(), (2201.0, 187.0));

        Self {
            name: "L8_OLI".to_string(),
            bands,
        }
    }

    pub fn new_l9() -> Self {
        let mut sensor = Self::new_l8();
        sensor.name = "L9_OLI".to_string();
        sensor
    }
}

impl Sensor for LandsatSensor {
    fn name(&self) -> &str {
        &self.name
    }

    fn parse_metadata(&self, path: &str) -> Result<Metadata> {
        parse_mtl(path)
    }

    fn band_names(&self) -> Vec<String> {
        self.bands.keys().cloned().collect()
    }

    fn wavelength(&self, band: &str) -> Option<f64> {
        self.bands.get(band).map(|(w, _)| *w)
    }

    fn bandwidth(&self, band: &str) -> Option<f64> {
        self.bands.get(band).map(|(_, b)| *b)
    }
}
