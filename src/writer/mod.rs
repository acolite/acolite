//! Scene writing — output processed data to various formats
//!
//! Dispatch strategy based on band count:
//! - Hyperspectral (>50 bands): GeoZarr — 3D chunked array, efficient for 100s of bands
//! - Multispectral/superspectral (≤50 bands): per-band COG — widely compatible

pub mod cog;
#[cfg(feature = "zarr")]
pub mod geozarr;

pub use cog::{cog_available, write_cog};
#[cfg(feature = "zarr")]
pub use geozarr::write_geozarr;

use crate::core::{BandData, Metadata};
use crate::Result;

/// Band count threshold: above this, use GeoZarr; at or below, use COG
pub const HYPERSPECTRAL_THRESHOLD: usize = 50;

/// Automatically choose output format based on band count
///
/// - >50 bands → GeoZarr (.zarr directory) [requires `zarr` feature]
/// - ≤50 bands → Cloud Optimized GeoTIFF (.tif)
pub fn write_auto(output_path: &str, bands: &[BandData<f64>], metadata: &Metadata) -> Result<()> {
    #[cfg(feature = "zarr")]
    if bands.len() > HYPERSPECTRAL_THRESHOLD {
        let zarr_path = if output_path.ends_with(".zarr") {
            output_path.to_string()
        } else {
            format!("{}.zarr", output_path.trim_end_matches(".tif"))
        };
        log::info!("{} bands → GeoZarr format", bands.len());
        return write_geozarr(&zarr_path, bands, metadata);
    }

    let tif_path = if output_path.ends_with(".tif") {
        output_path.to_string()
    } else {
        format!("{}.tif", output_path.trim_end_matches(".zarr"))
    };
    log::info!("{} bands → COG format", bands.len());
    write_cog(&tif_path, bands, metadata)
}
