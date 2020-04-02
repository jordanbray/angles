use integer_angles::Angle;
use criterion::{black_box, criterion_group, criterion_main, Criterion};


fn bench_atan_angle(c: &mut Criterion) {
    c.bench_function("atan_angle", |b| b.iter(|| Angle::atan(black_box(1.1))));
}

fn bench_atan_builtin(c: &mut Criterion) {
    c.bench_function("atan_buitin", |b| b.iter(|| f64::atan(black_box(1.1))));
}

criterion_group!(benches, bench_atan_angle, bench_atan_builtin);
criterion_main!(benches);
