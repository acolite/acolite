//! Radiometric calibration functions

use ndarray::Array2;
use crate::Result;

/// Convert DN to radiance
pub fn dn_to_radiance(
    dn: &Array2<u16>,
    mult: f64,
    add: f64,
) -> Array2<f64> {
    dn.mapv(|v| (v as f64) * mult + add)
}

/// Convert DN to TOA reflectance
pub fn dn_to_reflectance(
    dn: &Array2<u16>,
    mult: f64,
    add: f64,
    sun_elevation: f64,
) -> Array2<f64> {
    let sin_sun_elev = sun_elevation.to_radians().sin();
    dn.mapv(|v| ((v as f64) * mult + add) / sin_sun_elev)
}

/// Convert radiance to reflectance
pub fn radiance_to_reflectance(
    radiance: &Array2<f64>,
    sun_zenith: f64,
    distance: f64,
    f0: f64,
) -> Array2<f64> {
    let cos_sza = sun_zenith.to_radians().cos();
    let factor = (std::f64::consts::PI * distance * distance) / (f0 * cos_sza);
    radiance.mapv(|v| v * factor)
}

/// Calculate Earth-Sun distance correction factor
pub fn earth_sun_distance(doy: u16) -> f64 {
    let d = doy as f64;
    let theta = 2.0 * std::f64::consts::PI * (d - 1.0) / 365.0;
    1.0 - 0.01672 * theta.cos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;
    
    #[test]
    fn test_dn_to_radiance() {
        let dn = arr2(&[[1000, 2000], [3000, 4000]]);
        let rad = dn_to_radiance(&dn, 0.01, 0.0);
        assert_eq!(rad[[0, 0]], 10.0);
        assert_eq!(rad[[1, 1]], 40.0);
    }
    
    #[test]
    fn test_earth_sun_distance() {
        let d = earth_sun_distance(1);
        assert!(d > 0.98 && d < 1.02);
        
        let d_mid = earth_sun_distance(182);
        assert!(d_mid > 0.98 && d_mid < 1.02);
    }
}
