//! I/O module for reading and writing various formats

pub mod netcdf;
pub mod zarr;

pub use self::netcdf::{NetCdfReader, NetCdfWriter};
pub use self::zarr::{ZarrReader, ZarrWriter};
