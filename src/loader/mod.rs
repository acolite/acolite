//! Scene loading — reads satellite data from disk or remote sources

pub mod geotiff;
pub mod landsat;
pub mod sentinel2;
pub mod pace;
pub mod source;

pub use geotiff::read_geotiff_band;
pub use landsat::{load_landsat_scene, load_landsat_bands};
pub use sentinel2::{load_sentinel2_scene, S2Scene, S2_AC_BANDS};
#[cfg(feature = "netcdf")]
pub use pace::{load_pace_l1b, PaceScene};
