//! ACOLITE-RS: High-performance atmospheric correction for aquatic remote sensing
//!
//! Architecture: Loader → AC (Processor) → Writer

pub mod core;
pub mod error;
pub mod sensors;
pub mod auth;
pub mod loader;
pub mod writer;
pub mod ac;
pub mod pipeline;
pub mod parallel;
pub mod resample;
pub mod simd;

pub use error::{AcoliteError, Result};
pub use pipeline::{Pipeline, ProcessingConfig};
pub use core::{BandData, GeoTransform, Metadata, Projection};
pub use auth::{Credentials, aws_profile};
pub use loader::{load_landsat_scene, load_landsat_bands};
#[cfg(feature = "netcdf")]
pub use loader::{load_pace_l1b, PaceScene};
pub use writer::{write_cog, write_geozarr, write_auto, cog_available};
pub use sensors::{LandsatSensor, Sentinel2Sensor, Sentinel3Sensor, PaceOciSensor};
pub use resample::{resample, ResampleMethod};
pub use loader::source::cmr::{search_pace_l1b, search_pace_scene, CmrGranule};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
