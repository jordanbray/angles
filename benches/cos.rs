use angles::Angle;
use nalgebra::RealField;
use criterion::{black_box, criterion_group, criterion_main, Criterion};


fn bench_cos_angle(c: &mut Criterion) {
    c.bench_function("cos_angle", |b| b.iter(|| black_box(Angle::pi()).cos::<f64>()));
}

fn bench_cos_builtin(c: &mut Criterion) {
    c.bench_function("cos_buitin", |b| b.iter(|| black_box(f64::pi()).cos()));
}

criterion_group!(benches, bench_cos_angle, bench_cos_builtin);
criterion_main!(benches);
