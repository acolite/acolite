//! Sensor abstraction and implementations

pub mod sensor;
pub mod landsat;
pub mod sentinel2;
pub mod sentinel3;
pub mod pace;
pub mod mtl;
pub mod s2_xml;

pub use sensor::Sensor;
pub use mtl::parse_mtl;
pub use s2_xml::parse_s2_metadata;
pub use landsat::LandsatSensor;
pub use sentinel2::Sentinel2Sensor;
pub use sentinel3::Sentinel3Sensor;
pub use pace::PaceOciSensor;
