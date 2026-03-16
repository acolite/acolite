//! Scene loading — reads satellite data from disk or remote sources

pub mod geotiff;
pub mod landsat;
pub mod pace;
pub mod sentinel2;
pub mod sentinel3;
pub mod sentinel3_l2r;
pub mod source;

pub use geotiff::read_geotiff_band;
pub use landsat::{load_landsat_bands, load_landsat_scene};
pub use landsat::load_landsat_scene_limit;
#[cfg(feature = "netcdf")]
pub use pace::{load_pace_l1b, PaceScene};
#[cfg(feature = "gdal-support")]
pub use sentinel2::load_sentinel2_scene;
#[cfg(feature = "gdal-support")]
pub use sentinel2::load_sentinel2_scene_limit;
pub use sentinel2::{S2Scene, S2_AC_BANDS};
pub use sentinel3::{OlciScene, OLCI_BANDS, parse_s3_manifest};
#[cfg(feature = "netcdf")]
pub use sentinel3::load_olci_scene;
