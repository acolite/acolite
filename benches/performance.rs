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

fn benchmark_parallel_processing(c: &mut Criterion) {
    let mut group = c.benchmark_group("band_processing");
    
    let metadata = Metadata::new("L8_OLI".to_string(), Utc::now());
    let config = ProcessingConfig::default();
    let pipeline = Pipeline::new(metadata, config);
    
    let proj = Projection::from_epsg(32610);
    let geotrans = GeoTransform::new(0.0, 30.0, 0.0, -30.0);
    
    for size in [100, 500, 1000].iter() {
        let bands: Vec<_> = (0..7).map(|i| {
            BandData::new(
                Array2::zeros((*size, *size)),
                400.0 + i as f64 * 100.0,
                50.0,
                format!("B{}", i + 1),
                proj.clone(),
                geotrans.clone(),
            )
        }).collect();
        
        group.bench_with_input(
            BenchmarkId::new("parallel", size),
            size,
            |b, _| b.iter(|| process_bands_parallel(&pipeline, bands.clone()))
        );
        
        group.bench_with_input(
            BenchmarkId::new("sequential", size),
            size,
            |b, _| b.iter(|| process_bands_sequential(&pipeline, bands.clone()))
        );
    }
    
    group.finish();
}

criterion_group!(benches, benchmark_rayleigh_tau, benchmark_earth_sun_distance, benchmark_parallel_processing);
criterion_main!(benches);
