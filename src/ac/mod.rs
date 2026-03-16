//! Atmospheric correction algorithms

pub mod aerlut;
pub mod ancillary;
pub mod calibration;
pub mod dsf;
pub mod gas;
pub mod gas_lut;
pub mod interp;
pub mod lut;
pub mod rayleigh;

pub use ancillary::Ancillary;
pub use calibration::{dn_to_radiance, dn_to_reflectance, earth_sun_distance};
pub use dsf::{estimate_dark_spectrum, AotCompute, DarkSpectrumMethod, DsfConfig, DsfMode};
pub use gas::{gas_correction, ozone_transmittance, water_vapor_transmittance};
pub use interp::RegularGridInterpolator;
pub use lut::{interp_lut_1d, interp_lut_2d, LutManager};
pub use rayleigh::{rayleigh_correction, rayleigh_optical_thickness};

// LUT-based exports (require full-io feature)
#[cfg(feature = "full-io")]
pub use aerlut::{
    load_acolite_luts, load_generic_lut, load_generic_luts, load_sensor_lut, rsr_convolve_gauss,
    AerosolLut, GenericAerosolLut,
};
#[cfg(feature = "full-io")]
pub use dsf::{
    dsf_correct_band, dsf_correct_band_tiled, optimize_aot, optimize_aot_fixed, optimize_aot_tiled,
    DsfResult, TiledDsfResult,
};
#[cfg(feature = "full-io")]
pub use dsf::{
    dsf_correct_band_generic, dsf_correct_band_tiled_generic, optimize_aot_fixed_generic,
    optimize_aot_generic, optimize_aot_tiled_generic,
};
pub use dsf::{dsf_correction_simple, optimize_aot_simple};
pub use gas_lut::{compute_gas_transmittance, read_ko3, read_rsr, BandRsr, GasTransmittance};
#[cfg(feature = "full-io")]
pub use gas_lut::{compute_gas_transmittance_hyper, HyperGasTransmittance};
pub use gas_lut::{
    compute_other_gas_transmittance as compute_gas_transmittance_other,
    compute_wv_transmittance as compute_gas_transmittance_wv,
};
#[cfg(feature = "full-io")]
pub use gas_lut::{load_gas_lut, load_wv_lut, GasLut, WvLut};
