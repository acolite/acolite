//! Metadata handling

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    pub sensor: String,
    pub datetime: DateTime<Utc>,
    pub sun_zenith: f64,
    pub sun_azimuth: f64,
    pub view_zenith: Option<f64>,
    pub view_azimuth: Option<f64>,
    pub attributes: HashMap<String, String>,
}

impl Metadata {
    pub fn new(sensor: String, datetime: DateTime<Utc>) -> Self {
        Self {
            sensor,
            datetime,
            sun_zenith: 0.0,
            sun_azimuth: 0.0,
            view_zenith: None,
            view_azimuth: None,
            attributes: HashMap::new(),
        }
    }

    pub fn set_geometry(&mut self, sun_zenith: f64, sun_azimuth: f64) {
        self.sun_zenith = sun_zenith;
        self.sun_azimuth = sun_azimuth;
    }

    pub fn add_attribute(&mut self, key: String, value: String) {
        self.attributes.insert(key, value);
    }
}
