use criterion::{black_box, criterion_group, criterion_main, Criterion};
use acolite_rs::ac::rayleigh::rayleigh_optical_thickness;

fn benchmark_rayleigh_tau(c: &mut Criterion) {
    c.bench_function("rayleigh_optical_thickness", |b| {
        b.iter(|| {
            rayleigh_optical_thickness(black_box(443.0), black_box(1013.25))
        })
    });
}

criterion_group!(benches, benchmark_rayleigh_tau);
criterion_main!(benches);
