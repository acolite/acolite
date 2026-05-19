//! Sun glint correction using Cox-Munk wave slope model and Fresnel reflectance.
//!
//! Computes sky reflectance (rsky) from wind speed and solar/view geometry,
//! then subtracts it from surface reflectance in-place.
//! Reference: Cox & Munk (1954), Mobley (1994) for Fresnel/Snell.

use crate::Result;
use ndarray::Array2;

/// Refractive index of water (WOPP, T=27°C, S=0 PSU) at common wavelengths.
/// Loaded lazily from data file in Python; here we interpolate from a small LUT
/// covering the sensor range 400–2200 nm (subset of data/Shared/WOPP).
const REFRI_WAVE: [f64; 10] = [
    400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 2200.0,
];
const REFRI_N: [f64; 10] = [
    1.3434, 1.3374, 1.3324, 1.3289, 1.3274, 1.3230, 1.3176, 1.3123, 1.3068, 1.2780,
];

/// Linearly interpolate refractive index for a given wavelength (nm).
fn refri_interp(wave_nm: f64) -> f64 {
    if wave_nm <= REFRI_WAVE[0] {
        return REFRI_N[0];
    }
    if wave_nm >= *REFRI_WAVE.last().unwrap() {
        return *REFRI_N.last().unwrap();
    }
    for i in 0..REFRI_WAVE.len() - 1 {
        if wave_nm >= REFRI_WAVE[i] && wave_nm <= REFRI_WAVE[i + 1] {
            let f = (wave_nm - REFRI_WAVE[i]) / (REFRI_WAVE[i + 1] - REFRI_WAVE[i]);
            return REFRI_N[i] + f * (REFRI_N[i + 1] - REFRI_N[i]);
        }
    }
    1.34 // fallback
}

/// Fresnel reflectance at incidence angle `theta` (radians) for air→water.
/// Matches Python `ac.rayleigh.sky_refl(theta, n_w)`.
fn fresnel(theta: f64, n_w: f64) -> f64 {
    let sin_t = theta.sin();
    let theta_t = (sin_t / n_w).asin(); // Snell's law
    let diff = theta - theta_t;
    let sum = theta + theta_t;
    if sum.abs() < 1e-12 {
        return 0.0;
    }
    0.5 * ((diff.sin() / sum.sin()).powi(2) + (diff.tan() / sum.tan()).powi(2))
}

/// Cox-Munk mean-square slope: σ² = 0.003 + 0.00512 * wind_speed
fn cox_munk_sigma2(wind: f32) -> f64 {
    0.003 + 0.00512 * wind as f64
}

/// Compute the specular (glint) half-angle ω from geometry.
/// ω = arccos(cos_2ω) / 2 where cos_2ω = μs·μv + sin(sza)·sin(vza)·cos(raa)
fn glint_half_angle(sza_rad: f64, vza_rad: f64, raa_rad: f64) -> f64 {
    let cos2omega = sza_rad.cos() * vza_rad.cos() + sza_rad.sin() * vza_rad.sin() * raa_rad.cos();
    cos2omega.clamp(-1.0, 1.0).acos() / 2.0
}

/// Compute scalar sky glint reflectance (rsky) for given geometry and wind.
///
/// This is the simplified Cox-Munk + Fresnel model:
/// 1. Compute specular half-angle ω from sun/view geometry
/// 2. Fresnel reflectance Rf at ω with wavelength-dependent refractive index
/// 3. Cox-Munk slope probability P(ω) ∝ exp(-tan²ω / σ²) / (σ² cos⁴ω)
/// 4. rsky = Rf · P / (4 · μs · μv)
pub fn compute_rsky(wind: f32, sza_deg: f32, vza_deg: f32, raa_deg: f32, wave_nm: u32) -> f64 {
    let sza = (sza_deg as f64).to_radians();
    let vza = (vza_deg as f64).to_radians();
    let raa = (raa_deg as f64).to_radians();

    let omega = glint_half_angle(sza, vza, raa);
    let n_w = refri_interp(wave_nm as f64);
    let rf = fresnel(omega, n_w);

    let sigma2 = cox_munk_sigma2(wind);
    let cos_omega = omega.cos();
    let tan_omega = omega.tan();

    // Cox-Munk slope probability density
    let p = (-tan_omega * tan_omega / sigma2).exp() / (sigma2 * cos_omega.powi(4));

    let mu_s = sza.cos();
    let mu_v = vza.cos();
    let denom = 4.0 * mu_s * mu_v;
    if denom.abs() < 1e-12 {
        return 0.0;
    }

    rf * p / denom
}

/// Subtract sun glint reflectance from surface reflectance in-place.
///
/// Computes rsky via Cox-Munk + Fresnel and subtracts from every pixel of `rhos`.
pub fn glint_correct(
    rhos: &mut Array2<f32>,
    wind: f32,
    sza: f32,
    vza: f32,
    raa: f32,
    wave_nm: u32,
) -> Result<()> {
    let rsky = compute_rsky(wind, sza, vza, raa, wave_nm) as f32;
    if rsky.abs() < 1e-12 {
        return Ok(());
    }
    rhos.mapv_inplace(|v| v - rsky);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_zero_wind_minimal_correction() {
        // Zero wind → very narrow slope distribution, but at non-specular geometry
        // the glint probability is near zero
        let rsky = compute_rsky(0.0, 30.0, 10.0, 90.0, 550);
        // Off-specular with zero wind: essentially no glint
        assert!(
            rsky < 1e-6,
            "rsky={rsky} should be near zero for off-specular + zero wind"
        );
    }

    #[test]
    fn test_high_wind_near_specular_reduces_rhos() {
        let mut rhos = Array2::from_elem((10, 10), 0.05_f32);
        let original = rhos[[5, 5]];
        // Near-specular geometry (sza≈vza, raa≈180°) + high wind
        glint_correct(&mut rhos, 10.0, 30.0, 30.0, 170.0, 550).unwrap();
        assert!(
            rhos[[5, 5]] < original,
            "rhos should decrease after glint correction"
        );
    }

    #[test]
    fn test_fresnel_normal_incidence() {
        // At normal incidence (θ=0), Fresnel → ((n-1)/(n+1))² for unpolarized
        let rf = fresnel(0.001, 1.34); // near-normal
        let expected = ((1.34_f64 - 1.0) / (1.34_f64 + 1.0)).powi(2);
        assert!((rf - expected).abs() < 0.005, "rf={rf} expected≈{expected}");
    }

    #[test]
    fn test_cox_munk_sigma2() {
        assert!((cox_munk_sigma2(0.0) - 0.003).abs() < 1e-9);
        assert!((cox_munk_sigma2(5.0) - 0.0286).abs() < 1e-4);
    }

    #[test]
    fn test_refri_interp_boundaries() {
        let n400 = refri_interp(400.0);
        assert!((n400 - 1.3434).abs() < 1e-4);
        let n550 = refri_interp(550.0);
        assert!(n550 > 1.33 && n550 < 1.35);
    }
}
