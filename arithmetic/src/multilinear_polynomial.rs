use ark_ff::{Field, PrimeField};
use ark_std::{end_timer, rand::RngCore, start_timer};
#[cfg(feature = "parallel")]
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use std::rc::Rc;

pub use ark_poly::DenseMultilinearExtension;

// Build a randomize list of mle-s whose sum is zero.
pub fn random_zero_mle_list<F: PrimeField, R: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut R,
) -> Vec<Rc<DenseMultilinearExtension<F>>> {
    let start = start_timer!(|| "sample random zero mle list");

    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv))
    }
    for _ in 0..(1 << nv) {
        multiplicands[0].push(F::zero());
        for e in multiplicands.iter_mut().skip(1) {
            e.push(F::rand(rng));
        }
    }

    let list = multiplicands
        .into_iter()
        .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    list
}

pub fn fix_variables<F: Field>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &[F],
) -> DenseMultilinearExtension<F> {
    assert!(
        partial_point.len() <= poly.num_vars,
        "invalid size of partial point"
    );
    let nv = poly.num_vars;
    let mut poly = poly.evaluations.to_vec();
    let dim = partial_point.len();
    // evaluate single variable of partial point from left to right
    for i in 0..dim {
        poly = fix_one_variable_helper(&poly, nv - i, &partial_point[i]);
    }

    DenseMultilinearExtension::<F>::from_evaluations_slice(nv - dim, &poly[..(1 << (nv - dim))])
}

pub fn fix_first_variable<F: Field>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &F,
) -> DenseMultilinearExtension<F> {
    assert!(poly.num_vars != 0, "invalid size of partial point");

    let nv = poly.num_vars;
    let res = fix_one_variable_helper(&poly.evaluations, nv, partial_point);
    DenseMultilinearExtension::<F>::from_evaluations_slice(nv - 1, &res)
}

fn fix_one_variable_helper<F: Field>(data: &[F], nv: usize, point: &F) -> Vec<F> {
    let mut res = vec![F::zero(); 1 << (nv - 1)];
    let one_minus_p = F::one() - point;

    // evaluate single variable of partial point from left to right
    #[cfg(not(feature = "parallel"))]
    for b in 0..(1 << (nv - 1)) {
        res[b] = data[b << 1] * one_minus_p + data[(b << 1) + 1] * point;
    }

    #[cfg(feature = "parallel")]
    if nv >= 13 {
        // on my computer we parallelization doesn't help till nv >= 13
        res.par_iter_mut().enumerate().for_each(|(i, x)| {
            *x = data[i << 1] * one_minus_p + data[(i << 1) + 1] * point;
        });
    } else {
        for b in 0..(1 << (nv - 1)) {
            res[b] = data[b << 1] * one_minus_p + data[(b << 1) + 1] * point;
        }
    }

    res
}
