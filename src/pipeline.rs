//! Processing pipeline

use crate::{Result, core::{BandData, Metadata}};
use ndarray::Array2;

/// Processing configuration
#[derive(Debug, Clone)]
pub struct ProcessingConfig {
    pub apply_rayleigh: bool,
    pub output_reflectance: bool,
    pub parallel: bool,
}

impl Default for ProcessingConfig {
    fn default() -> Self {
        Self {
            apply_rayleigh: true,
            output_reflectance: true,
            parallel: true,
        }
    }
}

/// Processing pipeline for a single scene
pub struct Pipeline {
    config: ProcessingConfig,
    metadata: Metadata,
}

impl Pipeline {
    pub fn new(metadata: Metadata, config: ProcessingConfig) -> Self {
        Self { config, metadata }
    }
    
    /// Process a single band from DN to surface reflectance
    pub fn process_band(&self, band: BandData<u16>) -> Result<BandData<f64>> {
        // Step 1: DN to TOA reflectance (placeholder - needs calibration coefficients)
        let toa = band.data.mapv(|v| v as f64 / 10000.0);
        
        // Step 2: Rayleigh correction (if enabled)
        let corrected = if self.config.apply_rayleigh {
            // Placeholder - will use full Rayleigh correction
            toa
        } else {
            toa
        };
        
        // Create output band
        Ok(BandData::new(
            corrected,
            band.wavelength,
            band.bandwidth,
            band.name,
            band.projection,
            band.geotransform,
        ))
    }
    
    /// Get metadata
    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Projection, GeoTransform};
    use ndarray::Array2;
    use chrono::Utc;
    
    #[test]
    fn test_pipeline_creation() {
        let metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
        let config = ProcessingConfig::default();
        let pipeline = Pipeline::new(metadata, config);
        assert_eq!(pipeline.metadata().sensor, "L8_OLI");
    }
}
