//! LUT (Look-Up Table) management

use crate::{Result, AcoliteError};
use std::path::{Path, PathBuf};
use std::fs;
use ndarray::Array;

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
}
