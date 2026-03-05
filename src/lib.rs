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

pub use error::{AcoliteError, Result};
pub use pipeline::{Pipeline, ProcessingConfig};
pub use resample::{resample, ResampleMethod};
pub use io::{NetCdfWriter, NetCdfReader, ZarrWriter, ZarrReader};

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
