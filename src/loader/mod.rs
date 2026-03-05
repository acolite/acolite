//! Scene loading — reads satellite data from disk or remote sources

pub mod geotiff;
pub mod landsat;
pub mod pace;
pub mod source;

pub use geotiff::read_geotiff_band;
pub use landsat::{load_landsat_scene, load_landsat_bands};
#[cfg(feature = "netcdf")]
pub use pace::{load_pace_l1b, PaceScene};
