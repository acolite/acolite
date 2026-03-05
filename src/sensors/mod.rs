//! Sensor abstraction and implementations

pub mod sensor;
pub mod landsat;
pub mod mtl;

pub use sensor::Sensor;
pub use mtl::parse_mtl;
