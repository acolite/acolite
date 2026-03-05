//! Dark Spectrum Fitting (DSF) atmospheric correction

use ndarray::Array2;
use crate::Result;

/// DSF configuration
#[derive(Debug, Clone)]
pub struct DsfConfig {
    pub dark_bands: Vec<String>,
    pub min_percentile: f64,
    pub max_iterations: usize,
}

impl Default for DsfConfig {
    fn default() -> Self {
        Self {
            dark_bands: vec!["B5".to_string(), "B6".to_string(), "B7".to_string()],
            min_percentile: 5.0,
            max_iterations: 10,
        }
    }
}

/// Estimate dark spectrum from NIR/SWIR bands
pub fn estimate_dark_spectrum(
    bands: &[Array2<f64>],
    percentile: f64,
) -> Vec<f64> {
    bands.iter().map(|band| {
        let mut values: Vec<f64> = band.iter().copied().collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let idx = ((values.len() as f64) * percentile / 100.0) as usize;
        values[idx.min(values.len() - 1)]
    }).collect()
}

/// Optimize AOT using dark spectrum
pub fn optimize_aot(
    dark_spectrum: &[f64],
    wavelengths: &[f64],
    _sun_zenith: f64,
    _view_zenith: f64,
) -> f64 {
    // Simplified AOT estimation
    // Full implementation would use LUT interpolation
    let mean_dark = dark_spectrum.iter().sum::<f64>() / dark_spectrum.len() as f64;
    
    // Simple linear relationship (placeholder)
    (mean_dark * 10.0).min(1.0).max(0.01)
}

/// Apply aerosol correction using DSF
pub fn dsf_correction(
    toa: &Array2<f64>,
    aot: f64,
    wavelength: f64,
    _sun_zenith: f64,
    _view_zenith: f64,
) -> Array2<f64> {
    // Simplified aerosol correction
    // Full implementation would use aerosol LUT
    let path_reflectance = aot * 0.1 * (550.0 / wavelength).powf(1.3);
    
    toa.mapv(|v| (v - path_reflectance).max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;
    
    #[test]
    fn test_estimate_dark_spectrum() {
        let band1 = arr2(&[[0.1, 0.2], [0.3, 0.4]]);
        let band2 = arr2(&[[0.05, 0.15], [0.25, 0.35]]);
        
        let dark = estimate_dark_spectrum(&[band1, band2], 25.0);
        assert_eq!(dark.len(), 2);
        assert!(dark[0] > 0.0);
        assert!(dark[1] > 0.0);
    }
    
    #[test]
    fn test_optimize_aot() {
        let dark = vec![0.01, 0.02, 0.03];
        let wl = vec![865.0, 1609.0, 2201.0];
        
        let aot = optimize_aot(&dark, &wl, 30.0, 0.0);
        assert!(aot > 0.0 && aot < 1.0);
    }
}
