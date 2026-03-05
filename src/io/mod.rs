//! I/O module for reading and writing various formats

pub mod netcdf;
pub mod zarr;
pub mod geotiff;
pub mod cog;
pub mod geozarr;
pub mod pace_nc;
pub mod pace_nc_real;
pub mod pace_zarr;

pub use self::netcdf::{NetCdfReader, NetCdfWriter};
pub use self::zarr::{ZarrReader, ZarrWriter};
pub use self::geotiff::{read_geotiff_band, write_geotiff_band, write_geotiff_multiband};
pub use self::cog::{write_cog, cog_available};
pub use self::geozarr::write_geozarr;
pub use self::pace_nc::read_pace_band_simple;
