//! Atmospheric correction algorithms

pub mod rayleigh;
pub mod calibration;

pub use rayleigh::{rayleigh_correction, rayleigh_optical_thickness};
pub use calibration::{dn_to_radiance, dn_to_reflectance, earth_sun_distance};
