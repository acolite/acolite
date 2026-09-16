//! Sensor abstraction and implementations

pub mod landsat;
pub mod mtl;
pub mod pace;
pub mod s2_xml;
pub mod sensor;
pub mod sentinel2;
pub mod sentinel3;

pub use landsat::LandsatSensor;
pub use mtl::parse_mtl;
pub use pace::PaceOciSensor;
pub use s2_xml::parse_s2_metadata;
pub use sensor::Sensor;
pub use sentinel2::Sentinel2Sensor;
pub use sentinel3::Sentinel3Sensor;
