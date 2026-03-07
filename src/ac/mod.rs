//! Atmospheric correction algorithms

pub mod rayleigh;
pub mod calibration;
pub mod lut;
pub mod gas;
pub mod dsf;
pub mod interp;
pub mod aerlut;
pub mod gas_lut;

pub use rayleigh::{rayleigh_correction, rayleigh_optical_thickness};
pub use calibration::{dn_to_radiance, dn_to_reflectance, earth_sun_distance};
pub use lut::{LutManager, interp_lut_1d, interp_lut_2d};
pub use gas::{gas_correction, ozone_transmittance, water_vapor_transmittance};
pub use dsf::{DsfConfig, DsfMode, DarkSpectrumMethod, AotCompute, estimate_dark_spectrum};
pub use interp::RegularGridInterpolator;

// LUT-based exports (require full-io feature)
#[cfg(feature = "full-io")]
pub use aerlut::{AerosolLut, load_acolite_luts, load_sensor_lut, GenericAerosolLut, load_generic_luts, load_generic_lut, rsr_convolve_gauss};
#[cfg(feature = "full-io")]
pub use dsf::{optimize_aot, dsf_correct_band, DsfResult, optimize_aot_tiled, dsf_correct_band_tiled, TiledDsfResult, optimize_aot_fixed};
#[cfg(feature = "full-io")]
pub use dsf::{optimize_aot_generic, optimize_aot_fixed_generic, optimize_aot_tiled_generic, dsf_correct_band_generic, dsf_correct_band_tiled_generic};
pub use dsf::{optimize_aot_simple, dsf_correction_simple};
pub use gas_lut::{compute_gas_transmittance, read_rsr, read_ko3, GasTransmittance, BandRsr};
pub use gas_lut::{compute_wv_transmittance as compute_gas_transmittance_wv, compute_other_gas_transmittance as compute_gas_transmittance_other};
pub use gas_lut::{compute_gas_transmittance_hyper, HyperGasTransmittance};
#[cfg(feature = "full-io")]
pub use gas_lut::{load_wv_lut, load_gas_lut, WvLut, GasLut};
