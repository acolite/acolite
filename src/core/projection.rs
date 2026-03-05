//! Geospatial projection handling

use crate::Result;

#[derive(Debug, Clone)]
pub struct Projection {
    pub epsg: Option<i32>,
    pub wkt: Option<String>,
    pub proj4: Option<String>,
}

impl Projection {
    pub fn from_epsg(epsg: i32) -> Self {
        Self {
            epsg: Some(epsg),
            wkt: None,
            proj4: None,
        }
    }
    
    pub fn from_wkt(wkt: String) -> Self {
        Self {
            epsg: None,
            wkt: Some(wkt),
            proj4: None,
        }
    }
    
    pub fn is_geographic(&self) -> bool {
        self.epsg.map(|e| e == 4326).unwrap_or(false)
    }
}
