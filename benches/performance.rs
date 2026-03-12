use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use acolite_rs::ac::rayleigh::rayleigh_optical_thickness;
use acolite_rs::ac::calibration::earth_sun_distance;
use acolite_rs::{Pipeline, ProcessingConfig};
use acolite_rs::core::{Metadata, BandData, Projection, GeoTransform};
use acolite_rs::parallel::{process_bands_parallel, process_bands_sequential};
use ndarray::Array2;
use chrono::Utc;

fn benchmark_rayleigh_tau(c: &mut Criterion) {
    c.bench_function("rayleigh_optical_thickness", |b| {
        b.iter(|| {
            rayleigh_optical_thickness(black_box(443.0), black_box(1013.25))
        })
    });
}

fn benchmark_earth_sun_distance(c: &mut Criterion) {
    c.bench_function("earth_sun_distance", |b| {
        b.iter(|| {
            earth_sun_distance(black_box(180))
        })
    });
}

fn make_landsat_bands(size: usize, nbands: usize) -> Vec<BandData<u16>> {
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(0.0, 30.0, 0.0, -30.0);
    let wavelengths = [443.0, 482.0, 561.0, 655.0, 865.0, 1609.0, 2201.0];

    (0..nbands).map(|i| {
        BandData::new(
            Array2::from_elem((size, size), 8000u16 + (i as u16) * 500),
            wavelengths[i % wavelengths.len()],
            50.0,
            format!("B{}", i + 1),
            proj.clone(),
            geotrans.clone(),
        )
    }).collect()
}

fn benchmark_parallel_processing(c: &mut Criterion) {
    let mut group = c.benchmark_group("landsat_band_processing");

    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(35.0, 140.0);
    let config = ProcessingConfig::default();
    let pipeline = Pipeline::new(metadata, config);

    for size in [100, 500, 1000, 2000].iter() {
        let bands = make_landsat_bands(*size, 7);

        group.bench_with_input(
            BenchmarkId::new("parallel_7band", size),
            size,
            |b, _| b.iter(|| process_bands_parallel(&pipeline, bands.clone()))
        );

        group.bench_with_input(
            BenchmarkId::new("sequential_7band", size),
            size,
            |b, _| b.iter(|| process_bands_sequential(&pipeline, bands.clone()))
        );
    }

    group.finish();
}

fn benchmark_landsat_full_scene(c: &mut Criterion) {
    // Simulate a realistic Landsat scene: 7 bands × 7000×7000 is too large for bench,
    // so use 2000×2000 as proxy
    let mut metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    metadata.set_geometry(35.0, 140.0);

    let config = ProcessingConfig {
        apply_rayleigh: true,
        apply_gas: true,
        apply_aerosol: true,
        output_reflectance: true,
        parallel: true,
        ozone: 0.3,
        water_vapor: 1.5,
    };
    let mut pipeline = Pipeline::new(metadata, config);
    pipeline.set_aot(0.1);

    let bands = make_landsat_bands(2000, 7);

    c.bench_function("landsat_full_scene_2k", |b| {
        b.iter(|| process_bands_parallel(&pipeline, bands.clone()))
    });
}

fn benchmark_s2_multi_resolution(c: &mut Criterion) {
    let mut group = c.benchmark_group("s2_multi_resolution");

    let mut metadata = Metadata::new("S2A_MSI".to_string(), Utc::now());
    metadata.set_geometry(28.0, 135.0);
    let config = ProcessingConfig::default();
    let pipeline = Pipeline::new(metadata, config);

    let proj = Projection::from_epsg(32632);

    // 10m bands (4 bands)
    for size in [500, 1000, 2000].iter() {
        let gt = GeoTransform::new(600000.0, 10.0, 5000000.0, -10.0);
        let bands: Vec<_> = [(492.4, "B02"), (559.8, "B03"), (664.6, "B04"), (832.8, "B08")]
            .iter()
            .map(|&(wl, name)| BandData::new(
                Array2::from_elem((*size, *size), 1500u16),
                wl, 30.0, name.to_string(), proj.clone(), gt.clone(),
            ))
            .collect();

        group.bench_with_input(
            BenchmarkId::new("10m_4band_par", size),
            size,
            |b, _| b.iter(|| process_bands_parallel(&pipeline, bands.clone())),
        );
    }

    // All 13 bands at 500×500
    {
        let gt = GeoTransform::new(600000.0, 10.0, 5000000.0, -10.0);
        let wls = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 945.1, 1373.5, 1613.7, 2202.4];
        let bands: Vec<_> = wls.iter().enumerate()
            .map(|(i, &wl)| BandData::new(
                Array2::from_elem((500, 500), 1200u16 + (i as u16) * 50),
                wl, 20.0, format!("B{}", i + 1), proj.clone(), gt.clone(),
            ))
            .collect();

        group.bench_function("13band_500_par", |b| {
            b.iter(|| process_bands_parallel(&pipeline, bands.clone()))
        });
        group.bench_function("13band_500_seq", |b| {
            b.iter(|| process_bands_sequential(&pipeline, bands.clone()))
        });
    }

    group.finish();
}

fn benchmark_resample(c: &mut Criterion) {
    use acolite_rs::{resample, ResampleMethod};

    let mut group = c.benchmark_group("resample");

    for size in [500, 1000, 2000].iter() {
        let data = Array2::from_elem((*size, *size), 0.15f64);

        group.bench_with_input(
            BenchmarkId::new("20m_to_10m_bilinear", size),
            size,
            |b, _| b.iter(|| resample(&data, 20, 10, ResampleMethod::Bilinear)),
        );
        group.bench_with_input(
            BenchmarkId::new("10m_to_20m_average", size),
            size,
            |b, _| b.iter(|| resample(&data, 10, 20, ResampleMethod::Average)),
        );
    }

    group.finish();
}

criterion_group!(benches,
    benchmark_rayleigh_tau,
    benchmark_earth_sun_distance,
    benchmark_parallel_processing,
    benchmark_landsat_full_scene,
    benchmark_s2_multi_resolution,
    benchmark_resample,
);
criterion_main!(benches);
