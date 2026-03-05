//! Atmospheric correction algorithms

pub mod rayleigh;
pub mod calibration;
pub mod lut;
pub mod gas;
pub mod dsf;

pub use rayleigh::{rayleigh_correction, rayleigh_optical_thickness};
pub use calibration::{dn_to_radiance, dn_to_reflectance, earth_sun_distance};
pub use lut::{LutManager, interp_lut_1d, interp_lut_2d};
pub use gas::{gas_correction, ozone_transmittance, water_vapor_transmittance};
pub use dsf::{DsfConfig, estimate_dark_spectrum, optimize_aot, dsf_correction};
