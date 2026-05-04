//! Property-based tests for Sentinel-3 OLCI atmospheric correction
//!
//! Uses proptest to verify physical invariants hold across random inputs.

use acolite_rs::ac::dsf::{estimate_dark_spectrum, DarkSpectrumMethod};
use acolite_rs::ac::gas::{gas_correction, ozone_transmittance, water_vapor_transmittance};
use acolite_rs::ac::rayleigh::rayleigh_optical_thickness;
use acolite_rs::core::{BandData, GeoTransform, Metadata, Projection};
use acolite_rs::loader::sentinel3::OLCI_BANDS;
use acolite_rs::pipeline::{Pipeline, ProcessingConfig};
use chrono::Utc;
use ndarray::Array2;
use proptest::prelude::*;

// ─── Strategies ──────────────────────────────────────────────────────────

fn olci_wavelength() -> impl Strategy<Value = f64> {
    prop_oneof![
        Just(400.0), Just(412.5), Just(442.5), Just(490.0), Just(510.0),
        Just(560.0), Just(620.0), Just(665.0), Just(673.75), Just(681.25),
        Just(708.75), Just(753.75), Just(761.25), Just(764.375), Just(767.5),
        Just(778.75), Just(865.0), Just(885.0), Just(900.0), Just(940.0),
        Just(1020.0),
    ]
}

fn valid_sza() -> impl Strategy<Value = f64> {
    1.0..85.0f64
}

fn valid_vza() -> impl Strategy<Value = f64> {
    0.0..60.0f64
}

fn valid_pressure() -> impl Strategy<Value = f64> {
    500.0..1100.0f64
}

fn valid_ozone() -> impl Strategy<Value = f64> {
    0.1..0.6f64
}

fn valid_wv() -> impl Strategy<Value = f64> {
    0.1..5.0f64
}

fn valid_aot() -> impl Strategy<Value = f64> {
    0.01..2.0f64
}

fn valid_toa() -> impl Strategy<Value = f64> {
    0.001..0.5f64
}

// ─── Rayleigh Properties ─────────────────────────────────────────────────

proptest! {
    #[test]
    fn prop_rayleigh_tau_positive(wl in olci_wavelength(), p in valid_pressure()) {
        let tau = rayleigh_optical_thickness(wl, p);
        prop_assert!(tau > 0.0, "τ_ray must be positive: {} at {} nm, {} hPa", tau, wl, p);
    }

    #[test]
    fn prop_rayleigh_tau_bounded(wl in olci_wavelength(), p in valid_pressure()) {
        let tau = rayleigh_optical_thickness(wl, p);
        prop_assert!(tau < 1.0, "τ_ray must be < 1: {} at {} nm", tau, wl);
    }

    #[test]
    fn prop_rayleigh_tau_pressure_proportional(
        wl in olci_wavelength(),
        p1 in 500.0..800.0f64,
        p2 in 800.0..1100.0f64,
    ) {
        let tau1 = rayleigh_optical_thickness(wl, p1);
        let tau2 = rayleigh_optical_thickness(wl, p2);
        prop_assert!(tau2 > tau1, "Higher pressure should give higher τ_ray");
        // Check linearity: τ ∝ P
        let ratio = (tau2 / tau1) / (p2 / p1);
        prop_assert!((ratio - 1.0).abs() < 0.01, "τ_ray should scale linearly with P");
    }

    #[test]
    fn prop_rayleigh_tau_wavelength_inverse(
        p in valid_pressure(),
    ) {
        // τ_ray ∝ λ^-4 approximately
        let tau_blue = rayleigh_optical_thickness(400.0, p);
        let tau_nir = rayleigh_optical_thickness(865.0, p);
        prop_assert!(tau_blue > tau_nir * 10.0, "Blue τ should be >> NIR τ");
    }
}

// ─── Gas Transmittance Properties ────────────────────────────────────────

proptest! {
    #[test]
    fn prop_ozone_transmittance_unit_interval(
        wl in olci_wavelength(),
        uoz in valid_ozone(),
        sza in valid_sza(),
    ) {
        let airmass = 1.0 / sza.to_radians().cos();
        let t = ozone_transmittance(wl, uoz, airmass);
        prop_assert!(t > 0.0 && t <= 1.0, "O3 transmittance out of [0,1]: {}", t);
    }

    #[test]
    fn prop_wv_transmittance_unit_interval(
        wl in olci_wavelength(),
        uwv in valid_wv(),
        sza in valid_sza(),
    ) {
        let airmass = 1.0 / sza.to_radians().cos();
        let t = water_vapor_transmittance(wl, uwv, airmass);
        prop_assert!(t > 0.0 && t <= 1.0, "WV transmittance out of [0,1]: {}", t);
    }

    #[test]
    fn prop_gas_correction_increases_or_preserves(
        wl in olci_wavelength(),
        uoz in valid_ozone(),
        uwv in valid_wv(),
        sza in valid_sza(),
        vza in valid_vza(),
        toa_val in valid_toa(),
    ) {
        let toa = Array2::from_elem((5, 5), toa_val);
        let corrected = gas_correction(&toa, wl, uoz, uwv, sza, vza);
        for (&orig, &corr) in toa.iter().zip(corrected.iter()) {
            prop_assert!(
                corr >= orig - 1e-12,
                "Gas correction decreased signal: {} → {} at {} nm",
                orig, corr, wl
            );
        }
    }
}

// ─── Dark Spectrum Properties ────────────────────────────────────────────

proptest! {
    #[test]
    fn prop_dark_spectrum_bounded_by_data(
        base_val in 0.01..0.3f64,
        noise in 0.0..0.05f64,
    ) {
        let bands: Vec<Array2<f64>> = (0..5).map(|i| {
            Array2::from_elem((50, 50), base_val + i as f64 * 0.01 + noise)
        }).collect();

        let dark = estimate_dark_spectrum(&bands, &DarkSpectrumMethod::Percentile(1.0));
        for (i, &v) in dark.iter().enumerate() {
            let expected = base_val + i as f64 * 0.01 + noise;
            prop_assert!(
                v <= expected + 1e-6,
                "Dark spectrum {} exceeds data {} for band {}",
                v, expected, i
            );
            prop_assert!(v > 0.0, "Dark spectrum must be positive");
        }
    }

    #[test]
    fn prop_dark_spectrum_intercept_stable(
        base_val in 0.01..0.2f64,
    ) {
        // Uniform data → intercept = value
        let bands: Vec<Array2<f64>> = (0..3).map(|_| {
            Array2::from_elem((100, 100), base_val)
        }).collect();

        let dark = estimate_dark_spectrum(&bands, &DarkSpectrumMethod::Intercept(200));
        for &v in &dark {
            prop_assert!(
                (v - base_val).abs() < 0.001,
                "Intercept of uniform data should equal value: {} vs {}",
                v, base_val
            );
        }
    }
}

// ─── Pipeline Properties ─────────────────────────────────────────────────

proptest! {
    #[test]
    fn prop_pipeline_deterministic(
        sza in valid_sza(),
        aot in valid_aot(),
        dn in 100u16..2000u16,
    ) {
        let run = |_| {
            let mut metadata = Metadata::new("S3A_OLCI".into(), Utc::now());
            metadata.set_geometry(sza, 140.0);
            metadata.view_zenith = Some(15.0);
            let config = ProcessingConfig::default();
            let mut pipeline = Pipeline::new(metadata, config);
            pipeline.set_aot(aot);

            let proj = Projection::from_epsg(4326);
            let geo = GeoTransform::new(3.0, 0.003, 51.0, -0.003);
            let band = BandData::new(
                Array2::from_elem((10, 10), dn),
                490.0, 10.0, "Oa04".into(),
                proj, geo,
            );
            pipeline.process_band(band).unwrap().data[[5, 5]]
        };

        let v1 = run(0);
        let v2 = run(1);
        prop_assert_eq!(v1, v2, "Pipeline must be deterministic");
    }

    #[test]
    fn prop_pipeline_output_finite(
        sza in 10.0..70.0f64,
        dn in 200u16..1500u16,
        wl_idx in 0..21usize,
    ) {
        let (name, wl, bw, _) = OLCI_BANDS[wl_idx];
        let mut metadata = Metadata::new("S3A_OLCI".into(), Utc::now());
        metadata.set_geometry(sza, 140.0);
        metadata.view_zenith = Some(15.0);
        let config = ProcessingConfig::default();
        let mut pipeline = Pipeline::new(metadata, config);
        pipeline.set_aot(0.1);

        let proj = Projection::from_epsg(4326);
        let geo = GeoTransform::new(3.0, 0.003, 51.0, -0.003);
        let band = BandData::new(
            Array2::from_elem((10, 10), dn),
            wl, bw, name.to_string(),
            proj, geo,
        );
        let result = pipeline.process_band(band).unwrap();
        let val = result.data[[5, 5]];
        prop_assert!(val.is_finite(), "Output must be finite for {} at SZA={}", name, sza);
    }
}

// ─── E0 / Solar Irradiance Properties ────────────────────────────────────

proptest! {
    #[test]
    fn prop_e0_positive_all_bands(idx in 0..21usize) {
        let (_, _, _, e0) = OLCI_BANDS[idx];
        prop_assert!(e0 > 0.0, "E0 must be positive");
        prop_assert!(e0 < 3000.0, "E0 unreasonably large: {}", e0);
    }
}
