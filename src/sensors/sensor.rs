//! Sensor trait definition

use crate::{core::Metadata, Result};

/// Trait for all sensor implementations
pub trait Sensor {
    /// Get sensor name
    fn name(&self) -> &str;

    /// Parse metadata from sensor-specific files
    fn parse_metadata(&self, path: &str) -> Result<Metadata>;

    /// Get band names for this sensor
    fn band_names(&self) -> Vec<String>;

    /// Get wavelength for a band
    fn wavelength(&self, band: &str) -> Option<f64>;

    /// Get bandwidth for a band
    fn bandwidth(&self, band: &str) -> Option<f64>;
}
