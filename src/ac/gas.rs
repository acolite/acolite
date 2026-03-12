//! Gas correction (ozone and water vapor)

use ndarray::Array2;
use crate::Result;
use crate::ac::lut::interp_lut_1d;

/// Ozone absorption coefficient (k_o3) from Anderson et al.
const K_O3_WAVELENGTHS: [f64; 7] = [400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0];
const K_O3_VALUES: [f64; 7] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // Placeholder

/// Calculate ozone transmittance
pub fn ozone_transmittance(wavelength: f64, ozone: f64, airmass: f64) -> f64 {
    let k_o3 = interp_lut_1d(wavelength, &K_O3_WAVELENGTHS, &K_O3_VALUES);
    (-k_o3 * ozone * airmass).exp()
}

/// Calculate water vapor transmittance (simplified)
pub fn water_vapor_transmittance(wavelength: f64, wv: f64, airmass: f64) -> f64 {
    // Simplified model - will be replaced with LUT
    if wavelength > 900.0 {
        let tau = 0.0001 * wv * airmass;
        (-tau).exp()
    } else {
        1.0
    }
}

/// Apply gas correction to TOA reflectance
pub fn gas_correction(
    toa: &Array2<f64>,
    wavelength: f64,
    ozone: f64,
    water_vapor: f64,
    sun_zenith: f64,
    view_zenith: f64,
) -> Array2<f64> {
    let airmass = 1.0 / sun_zenith.to_radians().cos() + 1.0 / view_zenith.to_radians().cos();
    
    let t_o3 = ozone_transmittance(wavelength, ozone, airmass);
    let t_wv = water_vapor_transmittance(wavelength, water_vapor, airmass);
    let t_gas = t_o3 * t_wv;
    
    toa.mapv(|v| v / t_gas)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_ozone_transmittance() {
        let t = ozone_transmittance(500.0, 0.3, 2.0);
        assert!(t > 0.0 && t <= 1.0);
    }
    
    #[test]
    fn test_water_vapor_transmittance() {
        let t = water_vapor_transmittance(940.0, 2.0, 2.0);
        assert!(t > 0.0 && t <= 1.0);
    }
}
