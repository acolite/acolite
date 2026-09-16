//! Error types for ACOLITE-RS

use thiserror::Error;

#[derive(Error, Debug)]
pub enum AcoliteError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("NetCDF error: {0}")]
    NetCdf(String),

    #[error("GDAL error: {0}")]
    Gdal(String),

    #[error("Invalid metadata: {0}")]
    InvalidMetadata(String),

    #[error("Unsupported sensor: {0}")]
    UnsupportedSensor(String),

    #[error("Processing error: {0}")]
    Processing(String),

    #[error("Configuration error: {0}")]
    Config(String),
}

pub type Result<T> = std::result::Result<T, AcoliteError>;
