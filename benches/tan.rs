use angles::Angle;
use nalgebra::RealField;
use criterion::{black_box, criterion_group, criterion_main, Criterion};


fn bench_tan_angle(c: &mut Criterion) {
    c.bench_function("tan_angle", |b| b.iter(|| black_box(Angle::pi()).tan()));
}

fn bench_tan_builtin(c: &mut Criterion) {
    c.bench_function("tan_buitin", |b| b.iter(|| black_box(f64::pi()).tan()));
}

criterion_group!(benches, bench_tan_angle, bench_tan_builtin);
criterion_main!(benches);
