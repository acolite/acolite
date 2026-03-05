//! I/O module for reading and writing various formats

pub mod netcdf;
pub mod zarr;
pub mod geotiff;

pub use self::netcdf::{NetCdfReader, NetCdfWriter};
pub use self::zarr::{ZarrReader, ZarrWriter};
pub use self::geotiff::{read_geotiff_band, write_geotiff_band, write_geotiff_multiband};
