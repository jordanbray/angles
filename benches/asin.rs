use angles::Angle;
use criterion::{black_box, criterion_group, criterion_main, Criterion};


fn bench_asin_angle(c: &mut Criterion) {
    c.bench_function("asin_angle", |b| b.iter(|| Angle::asin(black_box(1.1))));
}

fn bench_asin_builtin(c: &mut Criterion) {
    c.bench_function("asin_buitin", |b| b.iter(|| f64::asin(black_box(1.1))));
}

criterion_group!(benches, bench_asin_angle, bench_asin_builtin);
criterion_main!(benches);
