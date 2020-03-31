use angles::Angle;
use nalgebra::RealField;
use criterion::{black_box, criterion_group, criterion_main, Criterion};


fn bench_sin_angle(c: &mut Criterion) {
    c.bench_function("sin_angle", |b| b.iter(|| black_box(Angle::pi()).sin()));
}

fn bench_sin_builtin(c: &mut Criterion) {
    c.bench_function("sin_buitin", |b| b.iter(|| black_box(f64::pi()).sin()));
}

criterion_group!(benches, bench_sin_angle, bench_sin_builtin);
criterion_main!(benches);
