//! ACOLITE-RS: High-performance atmospheric correction for aquatic remote sensing
//!
//! This is a Rust port of ACOLITE, focusing on performance through parallelization
//! and efficient memory management.

pub mod core;
pub mod io;
pub mod sensors;
pub mod ac;
pub mod error;
pub mod pipeline;
pub mod parallel;
pub mod resample;
pub mod simd;
pub mod stac;
pub mod earthdata;

pub use error::{AcoliteError, Result};
pub use pipeline::{Pipeline, ProcessingConfig};
pub use resample::{resample, ResampleMethod};
pub use io::{NetCdfWriter, NetCdfReader, ZarrWriter, ZarrReader, 
             read_geotiff_band, write_geotiff_band, write_geotiff_multiband,
             write_cog, cog_available, write_geozarr};
pub use sensors::{LandsatSensor, Sentinel2Sensor, Sentinel3Sensor, PaceOciSensor};
pub use stac::{StacClient, StacItem, StacAsset, search_landsat, search_sentinel2};
pub use earthdata::{EarthdataAuth, search_pace_data, download_pace_file};

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
