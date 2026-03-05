//! ACOLITE-RS: High-performance atmospheric correction for aquatic remote sensing
//!
//! This is a Rust port of ACOLITE, focusing on performance through parallelization
//! and efficient memory management.

pub mod core;
pub mod io;
pub mod sensors;
pub mod ac;
pub mod error;

pub use error::{AcoliteError, Result};

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
