//! Core data structures and types

pub mod projection;
pub mod metadata;
pub mod band;

pub use projection::Projection;
pub use metadata::Metadata;
pub use band::{BandData, GeoTransform};
