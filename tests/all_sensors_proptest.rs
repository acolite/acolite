//! Property-based tests for all sensors: Landsat, Sentinel-2, Sentinel-3, PACE
//!
//! Verifies physical invariants hold across random inputs for the shared
//! atmospheric correction pipeline used by all sensors.

use acolite_rs::ac::dsf::{estimate_dark_spectrum, DarkSpectrumMethod};
use acolite_rs::ac::gas::{gas_correction, ozone_transmittance, water_vapor_transmittance};
use acolite_rs::ac::rayleigh::rayleigh_optical_thickness;
use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::pipeline::{Pipeline, ProcessingConfig};
use chrono::Utc;
use ndarray::Array2;
use proptest::prelude::*;

// ─── Strategies ──────────────────────────────────────────────────────────

fn any_wavelength() -> impl Strategy<Value = f64> {
    prop_oneof![
        // Landsat OLI
        Just(443.0), Just(482.0), Just(561.0), Just(655.0), Just(865.0), Just(1609.0), Just(2201.0),
        // Sentinel-2 MSI
        Just(442.7), Just(492.4), Just(559.8), Just(664.6), Just(704.1), Just(740.5),
        Just(782.8), Just(832.8), Just(864.7), Just(1613.7), Just(2202.4),
        // Sentinel-3 OLCI
        Just(400.0), Just(412.5), Just(442.5), Just(490.0), Just(510.0), Just(560.0),
        Just(620.0), Just(665.0), Just(708.75), Just(753.75), Just(865.0), Just(1020.0),
        // PACE OCI (sample)
        Just(350.0), Just(380.0), Just(410.0), Just(440.0), Just(500.0), Just(550.0),
        Just(600.0), Just(700.0), Just(800.0), Just(1000.0), Just(1250.0), Just(2100.0),
    ]
}

fn valid_sza() -> impl Strategy<Value = f64> { 1.0..85.0f64 }
fn valid_vza() -> impl Strategy<Value = f64> { 0.0..60.0f64 }
fn valid_pressure() -> impl Strategy<Value = f64> { 500.0..1100.0f64 }
fn valid_ozone() -> impl Strategy<Value = f64> { 0.1..0.6f64 }
fn valid_wv() -> impl Strategy<Value = f64> { 0.1..5.0f64 }
fn valid_aot() -> impl Strategy<Value = f64> { 0.01..2.0f64 }
fn valid_toa() -> impl Strategy<Value = f64> { 0.001..0.8f64 }

// ─── Rayleigh properties (all sensors) ───────────────────────────────────

proptest! {
    #[test]
    fn prop_rayleigh_tau_positive(wl in any_wavelength(), p in valid_pressure()) {
        let tau = rayleigh_optical_thickness(wl, p);
        prop_assert!(tau > 0.0, "τ_ray must be positive: got {} at {}nm", tau, wl);
    }

    #[test]
    fn prop_rayleigh_tau_bounded(wl in any_wavelength(), p in valid_pressure()) {
        let tau = rayleigh_optical_thickness(wl, p);
        prop_assert!(tau < 2.0, "τ_ray must be < 2.0: got {} at {}nm", tau, wl);
    }

    #[test]
    fn prop_rayleigh_tau_pressure_proportional(
        wl in any_wavelength(),
        p1 in 500.0..800.0f64,
        p2 in 800.0..1100.0f64,
    ) {
        let tau1 = rayleigh_optical_thickness(wl, p1);
        let tau2 = rayleigh_optical_thickness(wl, p2);
        prop_assert!(tau2 > tau1, "Higher pressure → higher τ_ray");
    }

    #[test]
    fn prop_rayleigh_tau_wavelength_inverse(
        p in valid_pressure(),
    ) {
        // Blue light scatters more than red (λ⁻⁴ dependence)
        let tau_blue = rayleigh_optical_thickness(443.0, p);
        let tau_red = rayleigh_optical_thickness(865.0, p);
        prop_assert!(tau_blue > tau_red, "Blue τ > Red τ");
    }
}

// ─── Gas transmittance properties (all sensors) ──────────────────────────

proptest! {
    #[test]
    fn prop_ozone_transmittance_unit_interval(
        wl in any_wavelength(),
        uoz in valid_ozone(),
        sza in valid_sza(),
    ) {
        let airmass = 1.0 / sza.to_radians().cos();
        let t = ozone_transmittance(wl, uoz, airmass);
        prop_assert!(t >= 0.0 && t <= 1.0, "O3 transmittance must be [0,1]: got {}", t);
    }

    #[test]
    fn prop_wv_transmittance_unit_interval(
        wl in any_wavelength(),
        uwv in valid_wv(),
        sza in valid_sza(),
    ) {
        let airmass = 1.0 / sza.to_radians().cos();
        let t = water_vapor_transmittance(wl, uwv, airmass);
        prop_assert!(t >= 0.0 && t <= 1.0, "WV transmittance must be [0,1]: got {}", t);
    }

    #[test]
    fn prop_gas_correction_increases_or_preserves(
        wl in any_wavelength(),
        uoz in valid_ozone(),
        uwv in valid_wv(),
        sza in valid_sza(),
        vza in valid_vza(),
    ) {
        let toa = Array2::from_elem((10, 10), 0.1);
        let corrected = gas_correction(&toa, wl, uoz, uwv, sza, vza);
        let mean_corr = corrected.iter().sum::<f64>() / corrected.len() as f64;
        // Gas correction divides by transmittance, so result ≥ input
        prop_assert!(mean_corr >= 0.099, "Gas correction should increase signal: got {}", mean_corr);
    }
}

// ─── Dark spectrum properties (all sensors) ──────────────────────────────

proptest! {
    #[test]
    fn prop_dark_spectrum_bounded_by_data(
        base in 0.01..0.5f64,
    ) {
        let bands: Vec<Array2<f64>> = (0..7).map(|i| {
            Array2::from_elem((100, 100), base + i as f64 * 0.01)
        }).collect();
        let dark = estimate_dark_spectrum(&bands, &DarkSpectrumMethod::Percentile(1.0));
        for (i, &d) in dark.iter().enumerate() {
            let band_min = base + i as f64 * 0.01;
            let band_max = band_min + 0.001; // uniform data, percentile ≈ min
            prop_assert!(d >= band_min - 0.01 && d <= band_max + 0.01,
                "Dark spectrum band {} = {} outside [{}, {}]", i, d, band_min, band_max);
        }
    }

    #[test]
    fn prop_dark_spectrum_intercept_stable(
        base in 0.01..0.3f64,
    ) {
        let bands: Vec<Array2<f64>> = (0..7).map(|_| {
            Array2::from_elem((100, 100), base)
        }).collect();
        let dark = estimate_dark_spectrum(&bands, &DarkSpectrumMethod::Intercept(50));
        for &d in &dark {
            prop_assert!((d - base).abs() < 0.02,
                "Intercept dark spectrum {} far from uniform value {}", d, base);
        }
    }
}

// ─── Pipeline properties (all sensors) ───────────────────────────────────

proptest! {
    #[test]
    fn prop_pipeline_output_physical(
        sza in valid_sza(),
        vza in valid_vza(),
        wl in any_wavelength(),
        dn in 500u16..8000u16,
    ) {
        let mut meta = Metadata::new("TEST_SENSOR".into(), Utc::now());
        meta.set_geometry(sza, 0.0);
        meta.view_zenith = Some(vza);

        let config = ProcessingConfig {
            apply_rayleigh: true,
            apply_gas: true,
            apply_aerosol: true,
            pressure: 1013.25,
            ozone: 0.3,
            water_vapor: 1.5,
            ..ProcessingConfig::default()
        };

        let mut pipeline = Pipeline::new(meta, config);
        pipeline.set_aot(0.1);

        let band = BandData::new(
            Array2::from_elem((10, 10), dn),
            wl, 10.0, "test".into(),
            Projection::from_epsg(4326),
            GeoTransform::new(0.0, 0.01, 0.0, -0.01),
        );

        let result = pipeline.process_band(band).unwrap();
        for &v in result.data.iter() {
            if v.is_finite() {
                prop_assert!(v >= -0.1 && v <= 1.5,
                    "ρs={} out of physical range at {}nm", v, wl);
            }
        }
    }

    #[test]
    fn prop_pipeline_deterministic(
        sza in valid_sza(),
        wl in any_wavelength(),
        dn in 1000u16..5000u16,
    ) {
        let make = || {
            let mut meta = Metadata::new("TEST".into(), Utc::now());
            meta.set_geometry(sza, 0.0);
            meta.view_zenith = Some(5.0);
            let config = ProcessingConfig::default();
            let mut p = Pipeline::new(meta, config);
            p.set_aot(0.1);
            let band = BandData::new(
                Array2::from_elem((10, 10), dn),
                wl, 10.0, "test".into(),
                Projection::from_epsg(4326),
                GeoTransform::new(0.0, 0.01, 0.0, -0.01),
            );
            p.process_band(band).unwrap()
        };
        let r1 = make();
        let r2 = make();
        prop_assert_eq!(r1.data, r2.data, "Pipeline must be deterministic");
    }
}

// ─── Subset processing properties ────────────────────────────────────────

proptest! {
    /// Geographic CRS subset: pixel bounds must be strictly inside the image.
    #[test]
    fn prop_geographic_subset_within_image(
        x_origin in -180.0..170.0f64,
        y_origin in -80.0..80.0f64,
        pixel_size in 0.001..0.1f64,
        rows in 100usize..2000,
        cols in 100usize..2000,
        // limit is a sub-box strictly inside the image
        r_frac_min in 0.1..0.4f64,
        r_frac_max in 0.6..0.9f64,
        c_frac_min in 0.1..0.4f64,
        c_frac_max in 0.6..0.9f64,
    ) {
        let pixel_height = -pixel_size; // north-up
        let south = y_origin + pixel_height * (r_frac_max * rows as f64);
        let north = y_origin + pixel_height * (r_frac_min * rows as f64);
        let west  = x_origin + pixel_size  * (c_frac_min * cols as f64);
        let east  = x_origin + pixel_size  * (c_frac_max * cols as f64);
        let limit = [south, west, north, east];

        let result = acolite_rs::loader::latlon_limit_to_pixel_subset(
            x_origin, pixel_size, y_origin, pixel_height,
            rows, cols, &limit, None,
        );
        let (r0, c0, nr, nc) = result.expect("Should find subset inside image");
        prop_assert!(r0 + nr <= rows, "Row subset exceeds image: {}+{} > {}", r0, nr, rows);
        prop_assert!(c0 + nc <= cols, "Col subset exceeds image: {}+{} > {}", c0, nc, cols);
        prop_assert!(nr > 0 && nc > 0, "Subset must be non-empty");
    }

    /// Subset of full image returns the full image (or close to it).
    #[test]
    fn prop_full_image_limit_returns_full(
        x_origin in -180.0..170.0f64,
        y_origin in -80.0..80.0f64,
        pixel_size in 0.001..0.1f64,
        rows in 50usize..500,
        cols in 50usize..500,
    ) {
        let pixel_height = -pixel_size;
        let south = y_origin + pixel_height * rows as f64 - pixel_size;
        let north = y_origin + pixel_size;
        let west  = x_origin - pixel_size;
        let east  = x_origin + pixel_size * cols as f64 + pixel_size;
        let limit = [south, west, north, east];

        let result = acolite_rs::loader::latlon_limit_to_pixel_subset(
            x_origin, pixel_size, y_origin, pixel_height,
            rows, cols, &limit, None,
        );
        let (_, _, nr, nc) = result.expect("Full-image limit should always find subset");
        prop_assert!(nr > 0 && nc > 0);
        // Should cover most of the image (within 2 pixels of edge)
        prop_assert!(nr >= rows.saturating_sub(2), "nr={} should be ≈ rows={}", nr, rows);
        prop_assert!(nc >= cols.saturating_sub(2), "nc={} should be ≈ cols={}", nc, cols);
    }

    /// Limit outside image returns None.
    #[test]
    fn prop_disjoint_limit_returns_none(
        x_origin in 0.0..100.0f64,
        y_origin in 0.0..50.0f64,
        pixel_size in 0.01..0.1f64,
        rows in 100usize..500,
        cols in 100usize..500,
    ) {
        let pixel_height = -pixel_size;
        // Limit entirely to the right of the image
        let image_east = x_origin + pixel_size * cols as f64;
        let limit = [y_origin - 10.0, image_east + 1.0, y_origin - 5.0, image_east + 5.0];

        let result = acolite_rs::loader::latlon_limit_to_pixel_subset(
            x_origin, pixel_size, y_origin, pixel_height,
            rows, cols, &limit, None,
        );
        prop_assert!(result.is_none(), "Disjoint limit should return None");
    }

    /// UTM subset: corners of a 0.5°×0.5° box in UTM zone 54S (South Australia)
    /// should produce a non-empty subset when the image covers that area.
    #[test]
    fn prop_utm_subset_sa_region(
        // Simulate a UTM zone 54S image covering Gulf St Vincent area
        // UTM 54S: easting ~300000–500000, northing ~6100000–6300000
        e_origin in 300_000.0..400_000.0f64,
        n_origin in 6_200_000.0..6_300_000.0f64,
        pixel_m in 10.0..100.0f64,
        rows in 500usize..2000,
        cols in 500usize..2000,
    ) {
        let pixel_height = -pixel_m;
        // Limit: Gulf St Vincent area [-35.0, 138.0, -34.5, 138.7]
        let limit = [-35.0_f64, 138.0, -34.5, 138.7];
        // WKT hinting UTM zone 54
        let wkt = r#"PROJCS["WGS 84 / UTM zone 54S",GEOGCS["WGS 84"],PROJECTION["Transverse_Mercator"],PARAMETER["central_meridian",141],AUTHORITY["EPSG","32754"]]"#;

        let result = acolite_rs::loader::latlon_limit_to_pixel_subset(
            e_origin, pixel_m, n_origin, pixel_height,
            rows, cols, &limit, Some(wkt),
        );
        // If the image covers the limit area, we should get a subset
        // (we can't assert non-None since the synthetic image may not overlap)
        if let Some((r0, c0, nr, nc)) = result {
            prop_assert!(r0 + nr <= rows, "Row subset out of bounds");
            prop_assert!(c0 + nc <= cols, "Col subset out of bounds");
            prop_assert!(nr > 0 && nc > 0);
        }
    }
}

// ─── Sensor-specific wavelength coverage ─────────────────────────────────

#[test]
fn test_landsat_bands_cover_vis_nir_swir() {
    use acolite_rs::sensors::{LandsatSensor, Sensor};
    let s = LandsatSensor::new_l9();
    let names = s.band_names();
    let wls: Vec<f64> = names.iter().filter_map(|n| s.wavelength(n)).collect();
    assert!(wls.iter().any(|&w| w < 500.0), "Landsat needs coastal/blue band");
    assert!(wls.iter().any(|&w| w > 800.0 && w < 1000.0), "Landsat needs NIR band");
    assert!(wls.iter().any(|&w| w > 1500.0), "Landsat needs SWIR band");
}

#[test]
fn test_sentinel2_bands_cover_vis_nir_swir() {
    use acolite_rs::sensors::{Sentinel2Sensor, Sensor};
    let s = Sentinel2Sensor::new_s2a();
    let names = s.band_names();
    let wls: Vec<f64> = names.iter().filter_map(|n| s.wavelength(n)).collect();
    assert!(wls.iter().any(|&w| w < 500.0), "S2 needs blue band");
    assert!(wls.iter().any(|&w| w > 700.0 && w < 900.0), "S2 needs red-edge/NIR");
    assert!(wls.iter().any(|&w| w > 1500.0), "S2 needs SWIR band");
    assert_eq!(names.len(), 13, "S2 should have 13 bands");
}

#[test]
fn test_sentinel3_olci_21_bands() {
    use acolite_rs::loader::sentinel3::OLCI_BANDS;
    assert_eq!(OLCI_BANDS.len(), 21);
    for i in 1..OLCI_BANDS.len() {
        assert!(OLCI_BANDS[i].1 > OLCI_BANDS[i-1].1,
            "OLCI bands not monotonic: {} vs {}", OLCI_BANDS[i-1].1, OLCI_BANDS[i].1);
    }
}

#[test]
fn test_pace_oci_hyperspectral() {
    use acolite_rs::sensors::{PaceOciSensor, Sensor};
    let s = PaceOciSensor;
    let names = s.band_names();
    assert!(names.len() > 100, "PACE OCI should have >100 bands, got {}", names.len());
    let first_wl = s.wavelength(&names[0]).unwrap();
    let last_wl = s.wavelength(&names[names.len() - 1]).unwrap();
    assert!(first_wl < 400.0, "PACE should start in UV, got {}", first_wl);
    assert!(last_wl > 2000.0, "PACE should extend to SWIR, got {}", last_wl);
}
