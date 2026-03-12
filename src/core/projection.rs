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
    
    pub fn to_wkt(&self) -> String {
        if let Some(ref wkt) = self.wkt {
            return wkt.clone();
        }
        
        if let Some(epsg) = self.epsg {
            // Simple WKT for common EPSG codes
            match epsg {
                4326 => r#"GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]"#.to_string(),
                _ => format!(r#"PROJCS["EPSG:{}"]"#, epsg),
            }
        } else {
            String::new()
        }
    }
    
    pub fn is_geographic(&self) -> bool {
        self.epsg.map(|e| e == 4326).unwrap_or(false)
    }
}
