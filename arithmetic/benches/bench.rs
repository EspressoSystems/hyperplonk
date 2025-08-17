// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

#[macro_use]
extern crate criterion;
use ark_poly::Polynomial;
use arithmetic::fix_variables;
use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{ops::Range, test_rng};
use criterion::{black_box, BenchmarkId, Criterion};

const NUM_VARIABLES_RANGE: Range<usize> = 10..21;

fn evaluation_op_bench<F: Field>(c: &mut Criterion) {
    let mut rng = test_rng();
    let mut group = c.benchmark_group("Evaluate");
    for nv in NUM_VARIABLES_RANGE {
        group.bench_with_input(BenchmarkId::new("evaluate native", nv), &nv, |b, &nv| {
            let poly = DenseMultilinearExtension::<F>::rand(nv, &mut rng);
            let point: Vec<_> = (0..nv).map(|_| F::rand(&mut rng)).collect();
            b.iter(|| black_box(poly.evaluate(&point)))
        });

        group.bench_with_input(BenchmarkId::new("evaluate optimized", nv), &nv, |b, &nv| {
            let poly = DenseMultilinearExtension::<F>::rand(nv, &mut rng);
            let point: Vec<_> = (0..nv).map(|_| F::rand(&mut rng)).collect();
            b.iter(|| black_box(fix_variables(&poly, &point)))
        });
    }
    group.finish();
}

fn bench_bls_381(c: &mut Criterion) {
    evaluation_op_bench::<Fr>(c);
}

criterion_group!(benches, bench_bls_381);
criterion_main!(benches);
