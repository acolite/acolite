//! Band data structures

use ndarray::Array2;
use crate::core::Projection;

#[derive(Debug, Clone)]
pub struct GeoTransform {
    pub x_origin: f64,
    pub pixel_width: f64,
    pub x_rotation: f64,
    pub y_origin: f64,
    pub y_rotation: f64,
    pub pixel_height: f64,
}

impl GeoTransform {
    pub fn new(x_origin: f64, pixel_width: f64, y_origin: f64, pixel_height: f64) -> Self {
        Self {
            x_origin,
            pixel_width,
            x_rotation: 0.0,
            y_origin,
            y_rotation: 0.0,
            pixel_height,
        }
    }
    
    pub fn pixel_to_geo(&self, x: usize, y: usize) -> (f64, f64) {
        let geo_x = self.x_origin + (x as f64) * self.pixel_width;
        let geo_y = self.y_origin + (y as f64) * self.pixel_height;
        (geo_x, geo_y)
    }
}

#[derive(Debug, Clone)]
pub struct BandData<T> {
    pub data: Array2<T>,
    pub wavelength: f64,
    pub bandwidth: f64,
    pub name: String,
    pub projection: Projection,
    pub geotransform: GeoTransform,
}

impl<T: Clone> BandData<T> {
    pub fn new(
        data: Array2<T>,
        wavelength: f64,
        bandwidth: f64,
        name: String,
        projection: Projection,
        geotransform: GeoTransform,
    ) -> Self {
        Self {
            data,
            wavelength,
            bandwidth,
            name,
            projection,
            geotransform,
        }
    }
    
    pub fn shape(&self) -> (usize, usize) {
        (self.data.nrows(), self.data.ncols())
    }
    
    pub fn width(&self) -> usize {
        self.data.ncols()
    }
    
    pub fn height(&self) -> usize {
        self.data.nrows()
    }
}
