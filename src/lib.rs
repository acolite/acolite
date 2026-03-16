//! ACOLITE-RS: High-performance atmospheric correction for aquatic remote sensing
//!
//! Architecture: Loader → AC (Processor) → Writer

pub mod ac;
pub mod auth;
pub mod core;
pub mod error;
pub mod loader;
pub mod parallel;
pub mod pipeline;
pub mod resample;
pub mod sensors;
pub mod simd;
pub mod writer;

pub use auth::{aws_profile, CdseCredentials, Credentials};
pub use core::{BandData, GeoTransform, Metadata, Projection};
pub use error::{AcoliteError, Result};
pub use loader::landsat::{find_mtl, parse_reflectance_coeffs, ReflectanceCoeffs};
#[cfg(feature = "gdal-support")]
pub use loader::load_sentinel2_scene;
pub use loader::source::cmr::{search_landsat, search_pace_l1b, search_pace_scene, CmrGranule};
pub use loader::source::stac::{download_cdse_sentinel2, download_stac_sentinel2};
pub use loader::{load_landsat_bands, load_landsat_scene};
#[cfg(feature = "netcdf")]
pub use loader::{load_pace_l1b, PaceScene};
pub use loader::{S2Scene, S2_AC_BANDS};
pub use pipeline::{Pipeline, ProcessingConfig};
pub use resample::{resample, ResampleMethod};
pub use sensors::{LandsatSensor, PaceOciSensor, Sentinel2Sensor, Sentinel3Sensor};
pub use writer::{cog_available, write_auto, write_cog};
#[cfg(feature = "zarr")]
pub use writer::write_geozarr;

// Ancillary data
pub use ac::ancillary::{
    download as ancillary_download, get as ancillary_get, list_files as ancillary_list_files,
    Ancillary,
};

// LUT-based atmospheric correction
pub use ac::{compute_gas_transmittance, GasTransmittance};
#[cfg(feature = "full-io")]
pub use ac::{compute_gas_transmittance_hyper, HyperGasTransmittance};
#[cfg(feature = "full-io")]
pub use ac::{
    dsf_correct_band, dsf_correct_band_tiled, load_acolite_luts, optimize_aot, optimize_aot_tiled,
    AerosolLut, DsfResult, TiledDsfResult,
};
#[cfg(feature = "full-io")]
pub use ac::{
    dsf_correct_band_generic, dsf_correct_band_tiled_generic, load_generic_luts,
    optimize_aot_fixed_generic, optimize_aot_generic, optimize_aot_tiled_generic,
    rsr_convolve_gauss, GenericAerosolLut,
};
pub use ac::{
    estimate_dark_spectrum, AotCompute, DarkSpectrumMethod, DsfConfig, RegularGridInterpolator,
};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
