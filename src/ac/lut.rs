//! LUT (Look-Up Table) management

use crate::{Result, AcoliteError};
use std::path::{Path, PathBuf};
use std::fs;
use ndarray::{Array1, Array2, Array3};
use serde::{Deserialize, Serialize};

/// LUT manager for atmospheric correction
pub struct LutManager {
    lut_dir: PathBuf,
}

impl LutManager {
    pub fn new(lut_dir: Option<PathBuf>) -> Self {
        let lut_dir = lut_dir.unwrap_or_else(|| {
            dirs::home_dir()
                .unwrap_or_else(|| PathBuf::from("."))
                .join(".acolite")
                .join("luts")
        });
        
        Self { lut_dir }
    }
    
    /// Ensure LUT directory exists
    pub fn ensure_dir(&self) -> Result<()> {
        fs::create_dir_all(&self.lut_dir)
            .map_err(|e| AcoliteError::Io(e))
    }
    
    /// Get path to a specific LUT
    pub fn lut_path(&self, name: &str) -> PathBuf {
        self.lut_dir.join(name)
    }
    
    /// Check if LUT exists
    pub fn has_lut(&self, name: &str) -> bool {
        self.lut_path(name).exists()
    }
    
    /// Download LUT from GitHub if not present
    pub fn download_lut(&self, name: &str, url: &str) -> Result<()> {
        if self.has_lut(name) {
            return Ok(());
        }
        
        self.ensure_dir()?;
        
        log::info!("Downloading LUT: {}", name);
        let response = reqwest::blocking::get(url)
            .map_err(|e| AcoliteError::Processing(format!("Download failed: {}", e)))?;
        
        if !response.status().is_success() {
            return Err(AcoliteError::Processing(format!("HTTP {}", response.status())));
        }
        
        let bytes = response.bytes()
            .map_err(|e| AcoliteError::Processing(format!("Read failed: {}", e)))?;
        
        let size = bytes.len();
        fs::write(self.lut_path(name), bytes)
            .map_err(|e| AcoliteError::Io(e))?;
        
        log::info!("Downloaded {} ({} bytes)", name, size);
        Ok(())
    }
    
    /// Load Rayleigh LUT (simplified JSON format)
    pub fn load_rayleigh_lut(&self) -> Result<RayleighLut> {
        let path = self.lut_path("rayleigh_lut.json");
        
        if !path.exists() {
            // Create default LUT
            return Ok(RayleighLut::default());
        }
        
        let content = fs::read_to_string(path)
            .map_err(|e| AcoliteError::Io(e))?;
        
        serde_json::from_str(&content)
            .map_err(|e| AcoliteError::Processing(format!("JSON parse error: {}", e)))
    }
}

/// Rayleigh LUT structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RayleighLut {
    pub wavelengths: Vec<f64>,
    pub sun_zenith: Vec<f64>,
    pub view_zenith: Vec<f64>,
    pub relative_azimuth: Vec<f64>,
    pub pressure: Vec<f64>,
    // Simplified: store as flat vec, reshape on use
    pub reflectance: Vec<f64>,
}

impl Default for RayleighLut {
    fn default() -> Self {
        // Minimal default LUT
        Self {
            wavelengths: vec![400.0, 500.0, 600.0, 700.0, 800.0],
            sun_zenith: vec![0.0, 30.0, 60.0],
            view_zenith: vec![0.0, 30.0, 60.0],
            relative_azimuth: vec![0.0, 90.0, 180.0],
            pressure: vec![1013.25],
            reflectance: vec![0.0; 5 * 3 * 3 * 3 * 1],
        }
    }
}

impl RayleighLut {
    /// Interpolate Rayleigh reflectance
    pub fn interpolate(
        &self,
        wavelength: f64,
        sun_zenith: f64,
        view_zenith: f64,
        relative_azimuth: f64,
        pressure: f64,
    ) -> f64 {
        // Simplified 1D interpolation on wavelength
        // Full implementation would do 5D interpolation
        interp_lut_1d(wavelength, &self.wavelengths, &self.reflectance[0..self.wavelengths.len()])
    }
}

/// Multi-dimensional LUT interpolation
pub fn interp_lut_1d(x: f64, x_vals: &[f64], y_vals: &[f64]) -> f64 {
    if x <= x_vals[0] {
        return y_vals[0];
    }
    if x >= x_vals[x_vals.len() - 1] {
        return y_vals[y_vals.len() - 1];
    }
    
    for i in 0..x_vals.len() - 1 {
        if x >= x_vals[i] && x <= x_vals[i + 1] {
            let t = (x - x_vals[i]) / (x_vals[i + 1] - x_vals[i]);
            return y_vals[i] * (1.0 - t) + y_vals[i + 1] * t;
        }
    }
    
    y_vals[0]
}

/// 2D bilinear interpolation
pub fn interp_lut_2d(
    x: f64, y: f64,
    x_vals: &[f64], y_vals: &[f64],
    z_vals: &[Vec<f64>],
) -> f64 {
    let x_interp: Vec<f64> = y_vals.iter().enumerate()
        .map(|(i, _)| interp_lut_1d(x, x_vals, &z_vals[i]))
        .collect();
    
    interp_lut_1d(y, y_vals, &x_interp)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_interp_1d() {
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 10.0, 20.0];
        
        assert_eq!(interp_lut_1d(0.5, &x, &y), 5.0);
        assert_eq!(interp_lut_1d(1.5, &x, &y), 15.0);
        assert_eq!(interp_lut_1d(-1.0, &x, &y), 0.0);
        assert_eq!(interp_lut_1d(3.0, &x, &y), 20.0);
    }
    
    #[test]
    fn test_rayleigh_lut_default() {
        let lut = RayleighLut::default();
        assert_eq!(lut.wavelengths.len(), 5);
        
        let r = lut.interpolate(500.0, 30.0, 0.0, 90.0, 1013.25);
        assert!(r >= 0.0);
    }
    
    #[test]
    #[ignore] // Requires network
    fn test_lut_download() {
        use tempfile::TempDir;
        
        let temp = TempDir::new().unwrap();
        let manager = LutManager::new(Some(temp.path().to_path_buf()));
        
        // Download README from acolite_luts repo
        let url = "https://raw.githubusercontent.com/acolite/acolite_luts/main/README.md";
        let result = manager.download_lut("test_readme.md", url);
        
        assert!(result.is_ok());
        assert!(manager.has_lut("test_readme.md"));
    }
}
