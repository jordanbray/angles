use angles::Angle;
use criterion::{black_box, criterion_group, criterion_main, Criterion};


fn bench_acos_angle(c: &mut Criterion) {
    c.bench_function("acos_angle", |b| b.iter(|| Angle::acos(black_box(1.1))));
}

fn bench_acos_builtin(c: &mut Criterion) {
    c.bench_function("acos_buitin", |b| b.iter(|| f64::acos(black_box(1.1))));
}

criterion_group!(benches, bench_acos_angle, bench_acos_builtin);
criterion_main!(benches);
