//! N-dimensional regular grid interpolation (equivalent to scipy.interpolate.RegularGridInterpolator)

/// N-dimensional regular grid interpolator.
/// Stores axes and a flat data array in C-order (last axis varies fastest).
pub struct RegularGridInterpolator {
    axes: Vec<Vec<f64>>,
    data: Vec<f32>,
    strides: Vec<usize>,
}

impl RegularGridInterpolator {
    pub fn new(axes: Vec<Vec<f64>>, data: Vec<f32>) -> Self {
        let ndim = axes.len();
        let mut strides = vec![1usize; ndim];
        for i in (0..ndim - 1).rev() {
            strides[i] = strides[i + 1] * axes[i + 1].len();
        }
        let expected = axes.iter().map(|a| a.len()).product::<usize>();
        assert_eq!(
            data.len(),
            expected,
            "data len {} != expected {}",
            data.len(),
            expected
        );
        Self {
            axes,
            data,
            strides,
        }
    }

    /// Interpolate at a single point. `point` must have same length as number of axes.
    pub fn interpolate(&self, point: &[f64]) -> f64 {
        let ndim = self.axes.len();
        assert_eq!(point.len(), ndim);

        // Find bracketing indices and fractional positions for each axis
        let mut lo = vec![0usize; ndim];
        let mut frac = vec![0.0f64; ndim];
        for d in 0..ndim {
            let ax = &self.axes[d];
            let v = point[d].clamp(ax[0], ax[ax.len() - 1]);
            // binary search
            let mut i = match ax
                .binary_search_by(|a| a.partial_cmp(&v).unwrap_or(std::cmp::Ordering::Equal))
            {
                Ok(i) => i,
                Err(i) => {
                    if i == 0 {
                        0
                    } else {
                        i - 1
                    }
                }
            };
            if i >= ax.len() - 1 {
                i = ax.len() - 2;
            }
            lo[d] = i;
            let span = ax[i + 1] - ax[i];
            frac[d] = if span > 0.0 { (v - ax[i]) / span } else { 0.0 };
        }

        // Multilinear interpolation: iterate over 2^ndim corners
        let ncorners = 1usize << ndim;
        let mut result = 0.0f64;
        for corner in 0..ncorners {
            let mut weight = 1.0f64;
            let mut flat_idx = 0usize;
            for d in 0..ndim {
                let hi = (corner >> (ndim - 1 - d)) & 1;
                let idx = lo[d] + hi;
                flat_idx += idx * self.strides[d];
                weight *= if hi == 1 { frac[d] } else { 1.0 - frac[d] };
            }
            result += weight * self.data[flat_idx] as f64;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1d_interp() {
        let axes = vec![vec![0.0, 1.0, 2.0]];
        let data = vec![0.0f32, 10.0, 20.0];
        let rgi = RegularGridInterpolator::new(axes, data);
        assert!((rgi.interpolate(&[0.5]) - 5.0).abs() < 1e-6);
        assert!((rgi.interpolate(&[1.5]) - 15.0).abs() < 1e-6);
    }

    #[test]
    fn test_2d_interp() {
        // 2x3 grid: axes [0,1] x [0,1,2], data = x + y
        let axes = vec![vec![0.0, 1.0], vec![0.0, 1.0, 2.0]];
        let data = vec![0.0f32, 1.0, 2.0, 1.0, 2.0, 3.0]; // row-major
        let rgi = RegularGridInterpolator::new(axes, data);
        assert!((rgi.interpolate(&[0.5, 1.0]) - 1.5).abs() < 1e-6);
        assert!((rgi.interpolate(&[0.0, 0.5]) - 0.5).abs() < 1e-6);
    }
}
