//! Rayleigh scattering correction

use crate::Result;
use ndarray::Array2;

/// Apply Rayleigh correction to TOA reflectance
pub fn rayleigh_correction(
    toa_reflectance: &Array2<f64>,
    _wavelength: f64,
    _sun_zenith: f64,
    _view_zenith: f64,
    _pressure: f64,
) -> Result<Array2<f64>> {
    // Placeholder implementation
    // Will implement full Rayleigh LUT interpolation
    let corrected = toa_reflectance.clone();
    Ok(corrected)
}

/// Calculate Rayleigh optical thickness
pub fn rayleigh_optical_thickness(wavelength: f64, pressure: f64) -> f64 {
    // Hansen and Travis (1974) formula
    let p0 = 1013.25; // Standard pressure
    let lambda_um = wavelength / 1000.0;

    let tau_r0 = 0.008569
        * lambda_um.powf(-4.0)
        * (1.0 + 0.0113 * lambda_um.powf(-2.0) + 0.00013 * lambda_um.powf(-4.0));

    tau_r0 * (pressure / p0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rayleigh_optical_thickness() {
        let tau = rayleigh_optical_thickness(443.0, 1013.25);
        assert!(tau > 0.0);
        assert!(tau < 1.0);
    }
}
