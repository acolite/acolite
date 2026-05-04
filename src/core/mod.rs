//! Core data structures and types

pub mod band;
pub mod metadata;
pub mod projection;

pub use band::{BandData, GeoTransform};
pub use metadata::Metadata;
pub use projection::Projection;
