//! SIMD-optimized operations for performance

use ndarray::Array2;

/// SIMD-optimized array operations
pub trait SimdOps {
    /// Element-wise multiplication with scalar (SIMD)
    fn mul_scalar_simd(&self, scalar: f64) -> Self;
    
    /// Element-wise addition with scalar (SIMD)
    fn add_scalar_simd(&self, scalar: f64) -> Self;
}

impl SimdOps for Array2<f64> {
    fn mul_scalar_simd(&self, scalar: f64) -> Self {
        // Use ndarray's built-in SIMD when available
        self.mapv(|v| v * scalar)
    }
    
    fn add_scalar_simd(&self, scalar: f64) -> Self {
        self.mapv(|v| v + scalar)
    }
}

/// Fast parallel reduction operations
pub fn parallel_sum(data: &Array2<f64>) -> f64 {
    use rayon::prelude::*;
    
    let rows: Vec<_> = data.axis_iter(ndarray::Axis(0)).collect();
    rows.par_iter()
        .map(|row| row.sum())
        .sum()
}

/// Fast parallel min/max
pub fn parallel_minmax(data: &Array2<f64>) -> (f64, f64) {
    use rayon::prelude::*;
    
    let rows: Vec<_> = data.axis_iter(ndarray::Axis(0)).collect();
    rows.par_iter()
        .map(|row| {
            let min = row.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let max = row.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            (min, max)
        })
        .reduce(
            || (f64::INFINITY, f64::NEG_INFINITY),
            |(min1, max1), (min2, max2)| (min1.min(min2), max1.max(max2))
        )
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;
    
    #[test]
    fn test_simd_ops() {
        let data = arr2(&[[1.0, 2.0], [3.0, 4.0]]);
        
        let result = data.mul_scalar_simd(2.0);
        assert_eq!(result[[0, 0]], 2.0);
        assert_eq!(result[[1, 1]], 8.0);
        
        let result = data.add_scalar_simd(10.0);
        assert_eq!(result[[0, 0]], 11.0);
        assert_eq!(result[[1, 1]], 14.0);
    }
    
    #[test]
    fn test_parallel_sum() {
        let data = arr2(&[[1.0, 2.0], [3.0, 4.0]]);
        let sum = parallel_sum(&data);
        assert_eq!(sum, 10.0);
    }
    
    #[test]
    fn test_parallel_minmax() {
        let data = arr2(&[[1.0, 5.0], [2.0, 4.0]]);
        let (min, max) = parallel_minmax(&data);
        assert_eq!(min, 1.0);
        assert_eq!(max, 5.0);
    }
}
