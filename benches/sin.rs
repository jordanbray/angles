use criterion::{black_box, criterion_group, criterion_main, Criterion};
use integer_angles::Angle;
use nalgebra::RealField;

fn bench_sin_angle(c: &mut Criterion) {
    c.bench_function("sin_angle", |b| {
        b.iter(|| black_box(Angle::pi()).sin::<f64>())
    });
}

fn bench_sin_builtin(c: &mut Criterion) {
    c.bench_function("sin_buitin", |b| b.iter(|| black_box(f64::pi()).sin()));
}

criterion_group!(benches, bench_sin_angle, bench_sin_builtin);
criterion_main!(benches);
