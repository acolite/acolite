//! Scene loading — reads satellite data from disk or remote sources

pub mod geotiff;
pub mod landsat;
pub mod pace;
pub mod sentinel2;
pub mod sentinel3;
pub mod sentinel3_l2r;
pub mod source;

pub use geotiff::read_geotiff_band;
pub use landsat::{load_landsat_bands, load_landsat_scene};
pub use landsat::load_landsat_scene_limit;
#[cfg(feature = "netcdf")]
pub use pace::{load_pace_l1b, PaceScene};
#[cfg(feature = "gdal-support")]
pub use sentinel2::load_sentinel2_scene;
#[cfg(feature = "gdal-support")]
pub use sentinel2::load_sentinel2_scene_limit;
pub use sentinel2::{S2Scene, S2_AC_BANDS};
pub use sentinel3::{OlciScene, OLCI_BANDS, parse_s3_manifest};
#[cfg(feature = "netcdf")]
pub use sentinel3::load_olci_scene;

/// Convert a geographic bounding box `[south, west, north, east]` (degrees) to
/// pixel row/col bounds within a UTM-projected raster described by `geotransform`
/// `(x_origin, pixel_width, y_origin, pixel_height)` and image `(rows, cols)`.
///
/// Uses a pure-Rust WGS-84 UTM approximation accurate to ~1 m, sufficient for
/// subsetting. Returns `(r_min, c_min, n_rows, n_cols)` or `None` if the limit
/// does not intersect the image.
///
/// The UTM zone is inferred from the raster origin's easting (x_origin).
pub fn latlon_limit_to_pixel_subset(
    x_origin: f64,
    pixel_width: f64,
    y_origin: f64,
    pixel_height: f64,
    rows: usize,
    cols: usize,
    limit: &[f64; 4],
    wkt: Option<&str>,
) -> Option<(usize, usize, usize, usize)> {
    let (south, west, north, east) = (limit[0], limit[1], limit[2], limit[3]);

    // Detect if raster is geographic (pixel_width < 1 degree ≈ 111 km)
    if pixel_width.abs() < 1.0 && pixel_height.abs() < 1.0 {
        // Geographic CRS — direct pixel computation
        let c_min = ((west - x_origin) / pixel_width).floor().max(0.0) as usize;
        let c_max = ((east - x_origin) / pixel_width).ceil().min(cols as f64) as usize;
        let r_min = ((north - y_origin) / pixel_height).floor().max(0.0) as usize;
        let r_max = ((south - y_origin) / pixel_height).ceil().min(rows as f64) as usize;
        if r_max > r_min && c_max > c_min {
            return Some((r_min, c_min, r_max - r_min, c_max - c_min));
        }
        return None;
    }

    // Projected CRS (UTM or similar) — infer UTM zone from WKT or raster origin
    let utm_zone = infer_utm_zone_from_wkt(wkt).or_else(|| infer_utm_zone_from_easting(x_origin))?;
    let is_north = y_origin > 0.0; // northern hemisphere if northing > 0

    // Convert all four corners of the limit box to UTM
    let corners = [
        (south, west), (south, east), (north, west), (north, east),
    ];
    let utm_corners: Vec<(f64, f64)> = corners.iter()
        .map(|&(lat, lon)| latlon_to_utm(lat, lon, utm_zone, is_north))
        .collect();

    let utm_e_min = utm_corners.iter().map(|c| c.0).fold(f64::INFINITY, f64::min);
    let utm_e_max = utm_corners.iter().map(|c| c.0).fold(f64::NEG_INFINITY, f64::max);
    let utm_n_min = utm_corners.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
    let utm_n_max = utm_corners.iter().map(|c| c.1).fold(f64::NEG_INFINITY, f64::max);

    // Add a 1-pixel margin to avoid off-by-one clipping
    let margin_e = pixel_width.abs();
    let margin_n = pixel_height.abs();

    let c_min = ((utm_e_min - margin_e - x_origin) / pixel_width).floor().max(0.0) as usize;
    let c_max = ((utm_e_max + margin_e - x_origin) / pixel_width).ceil().min(cols as f64) as usize;
    // pixel_height is negative for north-up rasters
    let r_min = ((utm_n_max + margin_n - y_origin) / pixel_height).floor().max(0.0) as usize;
    let r_max = ((utm_n_min - margin_n - y_origin) / pixel_height).ceil().min(rows as f64) as usize;

    if r_max > r_min && c_max > c_min {
        Some((r_min, c_min, r_max - r_min, c_max - c_min))
    } else {
        None
    }
}

/// Infer UTM zone number (1–60) from a WKT projection string.
fn infer_utm_zone_from_wkt(wkt: Option<&str>) -> Option<u8> {
    let wkt = wkt?;
    // Look for "zone N" or "UTM zone N" patterns
    let lower = wkt.to_lowercase();
    if let Some(pos) = lower.find("zone ") {
        let rest = &lower[pos + 5..];
        let zone_str: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        return zone_str.parse().ok();
    }
    // Look for EPSG codes 32601–32660 (N) or 32701–32760 (S)
    if let Some(pos) = lower.find("epsg\",") {
        let rest = &lower[pos + 6..];
        let epsg_str: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        if let Ok(epsg) = epsg_str.parse::<u32>() {
            if (32601..=32660).contains(&epsg) { return Some((epsg - 32600) as u8); }
            if (32701..=32760).contains(&epsg) { return Some((epsg - 32700) as u8); }
        }
    }
    None
}

/// Infer UTM zone from easting value (UTM eastings are 100,000–900,000 m).
fn infer_utm_zone_from_easting(_easting: f64) -> Option<u8> {
    // Without knowing the zone, we can't infer it from easting alone.
    // Return None — caller must use WKT or fall back to no subsetting.
    None
}

/// Convert WGS-84 lat/lon (degrees) to UTM easting/northing (metres).
/// Uses the standard Karney series approximation (accurate to ~1 mm).
fn latlon_to_utm(lat_deg: f64, lon_deg: f64, zone: u8, _is_north: bool) -> (f64, f64) {
    use std::f64::consts::PI;
    let a = 6_378_137.0_f64;          // WGS-84 semi-major axis
    let f = 1.0 / 298.257_223_563;    // WGS-84 flattening
    let k0 = 0.9996;                   // UTM scale factor
    let e0 = 500_000.0;               // false easting
    let n0 = if lat_deg < 0.0 { 10_000_000.0 } else { 0.0 }; // false northing

    let lon0 = ((zone as f64 - 1.0) * 6.0 - 180.0 + 3.0).to_radians(); // central meridian
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();

    let n = f / (2.0 - f);
    let a_bar = a / (1.0 + n) * (1.0 + n*n/4.0 + n*n*n*n/64.0);

    let e = (2.0*f - f*f).sqrt();
    let e2 = e*e;
    let ep2 = e2 / (1.0 - e2);

    let t = lat.tan();
    let sigma = (e * (e * t / (1.0 + t*t).sqrt()).atanh()).sinh();
    let t_prime = t * (1.0 + sigma*sigma).sqrt() - sigma * (1.0 + t*t).sqrt();
    let xi_prime = (t_prime / (lon - lon0).cos()).atan();
    let eta_prime = ((lon - lon0).sin() / (1.0 + t_prime*t_prime).sqrt()).asinh();

    // Series coefficients (3rd order)
    let alpha = [
        n/2.0 - 2.0*n*n/3.0 + 5.0*n*n*n/16.0,
        13.0*n*n/48.0 - 3.0*n*n*n/5.0,
        61.0*n*n*n/240.0,
    ];

    let mut xi = xi_prime;
    let mut eta = eta_prime;
    for (j, &a_j) in alpha.iter().enumerate() {
        let j1 = (j + 1) as f64;
        xi += a_j * (2.0*j1*xi_prime).sin() * (2.0*j1*eta_prime).cosh();
        eta += a_j * (2.0*j1*xi_prime).cos() * (2.0*j1*eta_prime).sinh();
    }

    let easting = e0 + k0 * a_bar * eta;
    let northing = n0 + k0 * a_bar * xi;
    (easting, northing)
}
