//! Parallel processing utilities

use crate::{core::BandData, pipeline::Pipeline, Result};
use rayon::prelude::*;

/// Process multiple bands in parallel
pub fn process_bands_parallel(
    pipeline: &Pipeline,
    bands: Vec<BandData<u16>>,
) -> Result<Vec<BandData<f64>>> {
    bands
        .into_par_iter()
        .map(|band| pipeline.process_band(band))
        .collect()
}

/// Process f32 bands in parallel (for sensors like PACE where data is already TOA reflectance)
pub fn process_bands_parallel_f32(
    pipeline: &Pipeline,
    bands: Vec<BandData<f32>>,
) -> Result<Vec<BandData<f64>>> {
    bands
        .into_par_iter()
        .map(|band| pipeline.process_band_f32(band))
        .collect()
}

/// Process bands sequentially (for comparison/debugging)
pub fn process_bands_sequential(
    pipeline: &Pipeline,
    bands: Vec<BandData<u16>>,
) -> Result<Vec<BandData<f64>>> {
    bands
        .into_iter()
        .map(|band| pipeline.process_band(band))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{GeoTransform, Metadata, Projection};
    use crate::pipeline::ProcessingConfig;
    use chrono::Utc;
    use ndarray::Array2;

    #[test]
    fn test_parallel_processing() {
        let metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
        let config = ProcessingConfig::default();
        let pipeline = Pipeline::new(metadata, config);

        let proj = Projection::from_epsg(32610);
        let geotrans = GeoTransform::new(0.0, 30.0, 0.0, -30.0);

        let bands = vec![BandData::new(
            Array2::zeros((10, 10)),
            443.0,
            16.0,
            "B1".to_string(),
            proj.clone(),
            geotrans.clone(),
        )];

        let result = process_bands_parallel(&pipeline, bands);
        assert!(result.is_ok());
    }
}
