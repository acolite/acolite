//! PACE OCI sensor implementation

use crate::core::Metadata;
use crate::sensors::Sensor;
use crate::error::Result;

pub struct PaceOciSensor;

impl Sensor for PaceOciSensor {
    fn name(&self) -> &str {
        "PACE_OCI"
    }
    
    fn parse_metadata(&self, _path: &str) -> Result<Metadata> {
        // Placeholder - would parse NetCDF-4 metadata
        Ok(Metadata::new(self.name().to_string(), chrono::Utc::now()))
    }
    
    fn band_names(&self) -> Vec<String> {
        // 286 hyperspectral bands
        // UV: 315-400 nm (sparse)
        // Continuous: 340-890 nm (5nm sampling)
        // SWIR: 940-2260 nm (discrete)
        
        let mut names = Vec::new();
        
        // UV bands (315-400 nm, ~10 bands)
        for wl in [315, 325, 340, 350, 360, 370, 380, 390, 400] {
            names.push(format!("OCI_{}", wl));
        }
        
        // Continuous hyperspectral (340-890 nm, 5nm sampling, ~110 bands)
        let mut wl = 340;
        while wl <= 890 {
            names.push(format!("OCI_{}", wl));
            wl += 5;
        }
        
        // SWIR bands (940-2260 nm, ~20 bands)
        for wl in [940, 1038, 1250, 1378, 1615, 2130, 2260] {
            names.push(format!("OCI_{}", wl));
        }
        
        names
    }
    
    fn wavelength(&self, band: &str) -> Option<f64> {
        if let Some(wl_str) = band.strip_prefix("OCI_") {
            wl_str.parse::<f64>().ok()
        } else {
            None
        }
    }
    
    fn bandwidth(&self, band: &str) -> Option<f64> {
        let wl = self.wavelength(band)?;
        
        // Bandwidth varies by wavelength range
        if wl < 400.0 {
            Some(10.0) // UV bands
        } else if wl <= 890.0 {
            Some(5.0) // Continuous hyperspectral
        } else {
            Some(20.0) // SWIR bands
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_pace_oci_bands() {
        let sensor = PaceOciSensor;
        assert_eq!(sensor.name(), "PACE_OCI");
        
        let bands = sensor.band_names();
        assert!(bands.len() > 100); // Should have 100+ bands
        
        // Check specific wavelengths
        assert_eq!(sensor.wavelength("OCI_340"), Some(340.0));
        assert_eq!(sensor.wavelength("OCI_443"), Some(443.0));
        assert_eq!(sensor.wavelength("OCI_555"), Some(555.0));
        assert_eq!(sensor.wavelength("OCI_2260"), Some(2260.0));
        
        // Check bandwidths (5nm for continuous spectrum)
        assert_eq!(sensor.bandwidth("OCI_555"), Some(5.0));
        assert_eq!(sensor.bandwidth("OCI_2260"), Some(20.0));
    }
}
