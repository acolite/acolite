//! Processing pipeline

use crate::ac::ancillary::Ancillary;
use crate::ac::{dsf_correction_simple, gas_correction, rayleigh_correction};
use crate::{
    core::{BandData, Metadata},
    Result,
};

/// Processing configuration
#[derive(Debug, Clone)]
pub struct ProcessingConfig {
    pub apply_rayleigh: bool,
    pub apply_gas: bool,
    pub apply_aerosol: bool,
    pub output_reflectance: bool,
    pub parallel: bool,
    pub ozone: f64,
    pub water_vapor: f64,
    pub pressure: f64,
}

impl Default for ProcessingConfig {
    fn default() -> Self {
        Self {
            apply_rayleigh: true,
            apply_gas: true,
            apply_aerosol: true,
            output_reflectance: true,
            parallel: true,
            ozone: 0.3,        // cm-atm
            water_vapor: 1.5,  // g/cm²
            pressure: 1013.25, // hPa
        }
    }
}

impl ProcessingConfig {
    /// Create config from downloaded ancillary data
    pub fn from_ancillary(anc: &Ancillary) -> Self {
        Self {
            ozone: anc.uoz,
            water_vapor: anc.uwv,
            pressure: anc.pressure,
            ..Self::default()
        }
    }
}

/// Processing pipeline for a single scene
pub struct Pipeline {
    config: ProcessingConfig,
    metadata: Metadata,
    aot: Option<f64>,
}

impl Pipeline {
    pub fn new(metadata: Metadata, config: ProcessingConfig) -> Self {
        Self {
            config,
            metadata,
            aot: None,
        }
    }

    /// Set AOT (from DSF optimization)
    pub fn set_aot(&mut self, aot: f64) {
        self.aot = Some(aot);
    }

    /// Process a single band from DN to surface reflectance
    pub fn process_band(&self, band: BandData<u16>) -> Result<BandData<f64>> {
        // Step 1: DN to TOA reflectance (placeholder - needs calibration coefficients)
        let mut toa = band.data.mapv(|v| v as f64 / 10000.0);

        // Step 2: Gas correction
        if self.config.apply_gas {
            toa = gas_correction(
                &toa,
                band.wavelength,
                self.config.ozone,
                self.config.water_vapor,
                self.metadata.sun_zenith,
                self.metadata.view_zenith.unwrap_or(0.0),
            );
        }

        // Step 3: Rayleigh correction
        if self.config.apply_rayleigh {
            toa = rayleigh_correction(
                &toa,
                band.wavelength,
                self.metadata.sun_zenith,
                self.metadata.view_zenith.unwrap_or(0.0),
                self.config.pressure,
            )?;
        }

        // Step 4: Aerosol correction (DSF)
        let corrected = if let (true, Some(aot)) = (self.config.apply_aerosol, self.aot) {
            dsf_correction_simple(
                &toa,
                aot,
                band.wavelength,
                self.metadata.sun_zenith,
                self.metadata.view_zenith.unwrap_or(0.0),
            )
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

    /// Process a single band that is already TOA reflectance (f32, e.g. PACE OCI rhot)
    pub fn process_band_f32(&self, band: BandData<f32>) -> Result<BandData<f64>> {
        let mut toa = band.data.mapv(|v| v as f64);

        if self.config.apply_gas {
            toa = gas_correction(
                &toa,
                band.wavelength,
                self.config.ozone,
                self.config.water_vapor,
                self.metadata.sun_zenith,
                self.metadata.view_zenith.unwrap_or(0.0),
            );
        }
        if self.config.apply_rayleigh {
            toa = rayleigh_correction(
                &toa,
                band.wavelength,
                self.metadata.sun_zenith,
                self.metadata.view_zenith.unwrap_or(0.0),
                self.config.pressure,
            )?;
        }
        let corrected = if let (true, Some(aot)) = (self.config.apply_aerosol, self.aot) {
            dsf_correction_simple(
                &toa,
                aot,
                band.wavelength,
                self.metadata.sun_zenith,
                self.metadata.view_zenith.unwrap_or(0.0),
            )
        } else {
            toa
        };

        Ok(BandData::new(
            corrected,
            band.wavelength,
            band.bandwidth,
            band.name,
            band.projection,
            band.geotransform,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{GeoTransform, Projection};
    use chrono::Utc;
    use ndarray::Array2;

    #[test]
    fn test_pipeline_creation() {
        let metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
        let config = ProcessingConfig::default();
        let pipeline = Pipeline::new(metadata, config);
        assert_eq!(pipeline.metadata().sensor, "L8_OLI");
    }

    #[test]
    fn test_full_pipeline() {
        let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
        metadata.set_geometry(30.0, 135.0);

        let config = ProcessingConfig::default();
        let mut pipeline = Pipeline::new(metadata, config);
        pipeline.set_aot(0.1);

        let proj = Projection::from_epsg(32610);
        let geotrans = GeoTransform::new(0.0, 30.0, 0.0, -30.0);

        let band = BandData::new(
            Array2::from_elem((10, 10), 1000u16),
            443.0,
            16.0,
            "B1".to_string(),
            proj,
            geotrans,
        );

        let result = pipeline.process_band(band);
        assert!(result.is_ok());
    }
}
