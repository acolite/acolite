//! Multi-resolution resampling for Sentinel-2

use ndarray::Array2;
use crate::Result;

/// Resample method
#[derive(Debug, Clone, Copy)]
pub enum ResampleMethod {
    NearestNeighbor,
    Bilinear,
    Average,
}

/// Resample array to target resolution
pub fn resample(
    data: &Array2<f64>,
    from_res: u32,
    to_res: u32,
    method: ResampleMethod,
) -> Result<Array2<f64>> {
    if from_res == to_res {
        return Ok(data.clone());
    }
    
    let scale = to_res as f64 / from_res as f64;
    let (h, w) = data.dim();
    let new_h = (h as f64 / scale).round() as usize;
    let new_w = (w as f64 / scale).round() as usize;
    
    match method {
        ResampleMethod::NearestNeighbor => resample_nearest(data, new_h, new_w, scale),
        ResampleMethod::Bilinear => resample_bilinear(data, new_h, new_w, scale),
        ResampleMethod::Average => resample_average(data, new_h, new_w, scale),
    }
}

fn resample_nearest(data: &Array2<f64>, new_h: usize, new_w: usize, scale: f64) -> Result<Array2<f64>> {
    let mut result = Array2::zeros((new_h, new_w));
    
    for i in 0..new_h {
        for j in 0..new_w {
            let src_i = (i as f64 * scale).round() as usize;
            let src_j = (j as f64 * scale).round() as usize;
            result[[i, j]] = data[[src_i.min(data.nrows() - 1), src_j.min(data.ncols() - 1)]];
        }
    }
    
    Ok(result)
}

fn resample_bilinear(data: &Array2<f64>, new_h: usize, new_w: usize, scale: f64) -> Result<Array2<f64>> {
    let mut result = Array2::zeros((new_h, new_w));
    
    for i in 0..new_h {
        for j in 0..new_w {
            let src_i = i as f64 * scale;
            let src_j = j as f64 * scale;
            
            let i0 = src_i.floor() as usize;
            let i1 = (i0 + 1).min(data.nrows() - 1);
            let j0 = src_j.floor() as usize;
            let j1 = (j0 + 1).min(data.ncols() - 1);
            
            let di = src_i - i0 as f64;
            let dj = src_j - j0 as f64;
            
            let v00 = data[[i0, j0]];
            let v01 = data[[i0, j1]];
            let v10 = data[[i1, j0]];
            let v11 = data[[i1, j1]];
            
            result[[i, j]] = v00 * (1.0 - di) * (1.0 - dj)
                           + v01 * (1.0 - di) * dj
                           + v10 * di * (1.0 - dj)
                           + v11 * di * dj;
        }
    }
    
    Ok(result)
}

fn resample_average(data: &Array2<f64>, new_h: usize, new_w: usize, scale: f64) -> Result<Array2<f64>> {
    let mut result = Array2::zeros((new_h, new_w));
    
    for i in 0..new_h {
        for j in 0..new_w {
            let src_i_start = (i as f64 * scale) as usize;
            let src_i_end = ((i + 1) as f64 * scale).ceil() as usize;
            let src_j_start = (j as f64 * scale) as usize;
            let src_j_end = ((j + 1) as f64 * scale).ceil() as usize;
            
            let mut sum = 0.0;
            let mut count = 0;
            
            for si in src_i_start..src_i_end.min(data.nrows()) {
                for sj in src_j_start..src_j_end.min(data.ncols()) {
                    sum += data[[si, sj]];
                    count += 1;
                }
            }
            
            result[[i, j]] = if count > 0 { sum / count as f64 } else { 0.0 };
        }
    }
    
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;
    
    #[test]
    fn test_resample_nearest() {
        let data = arr2(&[[1.0, 2.0], [3.0, 4.0]]);
        // Upsampling from 20m to 10m: 2x2 -> 4x4
        let result = resample(&data, 20, 10, ResampleMethod::NearestNeighbor).unwrap();
        assert_eq!(result.dim(), (4, 4));
    }
    
    #[test]
    fn test_resample_same_resolution() {
        let data = arr2(&[[1.0, 2.0], [3.0, 4.0]]);
        let result = resample(&data, 10, 10, ResampleMethod::NearestNeighbor).unwrap();
        assert_eq!(result, data);
    }
}
