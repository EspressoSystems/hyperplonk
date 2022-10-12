use crate::{util::get_batched_nv, ArithErrors};
use ark_ff::{Field, PrimeField};
use ark_poly::MultilinearExtension;
use ark_std::{end_timer, rand::RngCore, start_timer};
#[cfg(feature = "parallel")]
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use std::rc::Rc;

pub use ark_poly::DenseMultilinearExtension;

/// Sample a random list of multilinear polynomials.
/// Returns
/// - the list of polynomials,
/// - its sum of polynomial evaluations over the boolean hypercube.
pub fn random_mle_list<F: PrimeField, R: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut R,
) -> (Vec<Rc<DenseMultilinearExtension<F>>>, F) {
    let start = start_timer!(|| "sample random mle list");
    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv))
    }
    let mut sum = F::zero();

    for _ in 0..(1 << nv) {
        let mut product = F::one();

        for e in multiplicands.iter_mut() {
            let val = F::rand(rng);
            e.push(val);
            product *= val;
        }
        sum += product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    (list, sum)
}

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

/// An MLE that represent an identity permutation: `f(index) \mapto index`
pub fn identity_permutation_mle<F: PrimeField>(
    num_vars: usize,
) -> Rc<DenseMultilinearExtension<F>> {
    let s_id_vec = (0..1u64 << num_vars).map(F::from).collect();
    Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars, s_id_vec,
    ))
}

/// An MLE that represent a random permutation
pub fn random_permutation_mle<F: PrimeField, R: RngCore>(
    num_vars: usize,
    rng: &mut R,
) -> Rc<DenseMultilinearExtension<F>> {
    let len = 1u64 << num_vars;
    let mut s_id_vec: Vec<F> = (0..len).map(F::from).collect();
    let mut s_perm_vec = vec![];
    for _ in 0..len {
        let index = rng.next_u64() as usize % s_id_vec.len();
        s_perm_vec.push(s_id_vec.remove(index));
    }
    Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars, s_perm_vec,
    ))
}

pub fn evaluate_opt<F: Field>(poly: &DenseMultilinearExtension<F>, point: &[F]) -> F {
    assert_eq!(poly.num_vars, point.len());
    fix_variables(poly, point).evaluations[0]
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
    for (i, point) in partial_point.iter().enumerate().take(dim) {
        poly = fix_one_variable_helper(&poly, nv - i, point);
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

pub fn evaluate_no_par<F: Field>(poly: &DenseMultilinearExtension<F>, point: &[F]) -> F {
    assert_eq!(poly.num_vars, point.len());
    fix_variables_no_par(poly, point).evaluations[0]
}

fn fix_variables_no_par<F: Field>(
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
    for i in 1..dim + 1 {
        let r = partial_point[i - 1];
        let one_minus_r = F::one() - r;
        for b in 0..(1 << (nv - i)) {
            poly[b] = poly[b << 1] * one_minus_r + poly[(b << 1) + 1] * r;
        }
    }
    DenseMultilinearExtension::from_evaluations_slice(nv - dim, &poly[..(1 << (nv - dim))])
}

/// TODO
pub fn pad_and_merge_polynomials<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
    target_nv: usize,
) -> Result<Rc<DenseMultilinearExtension<F>>, ArithErrors> {
    let nv = polynomials[0].num_vars();
    for poly in polynomials.iter() {
        if nv != poly.num_vars() {
            return Err(ArithErrors::InvalidParameters(
                "num_vars do not match for polynomials".to_string(),
            ));
        }
    }

    let merged_nv = get_batched_nv(nv, polynomials.len());
    if merged_nv > target_nv {
        return Err(ArithErrors::InvalidParameters(
            "target is smaller than merged nv".to_string(),
        ));
    }

    let mut scalars = vec![];
    for poly in polynomials.iter() {
        scalars.extend_from_slice(poly.to_evaluations().as_slice());
    }
    scalars.extend_from_slice(vec![F::zero(); (1 << target_nv) - scalars.len()].as_ref());
    Ok(Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        target_nv, scalars,
    )))
}

/// merge a set of polynomials. Returns an error if the
/// polynomials do not share a same number of nvs.
pub fn merge_polynomials<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
) -> Result<Rc<DenseMultilinearExtension<F>>, ArithErrors> {
    let nv = polynomials[0].num_vars();
    for poly in polynomials.iter() {
        if nv != poly.num_vars() {
            return Err(ArithErrors::InvalidParameters(
                "num_vars do not match for polynomials".to_string(),
            ));
        }
    }

    let merged_nv = get_batched_nv(nv, polynomials.len());
    let mut scalars = vec![];
    for poly in polynomials.iter() {
        scalars.extend_from_slice(poly.to_evaluations().as_slice());
    }
    scalars.extend_from_slice(vec![F::zero(); (1 << merged_nv) - scalars.len()].as_ref());
    Ok(Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        merged_nv, scalars,
    )))
}

pub fn fix_last_variables_no_par<F: PrimeField>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &[F],
) -> DenseMultilinearExtension<F> {
    let mut res = fix_last_variable_no_par(poly, partial_point.last().unwrap());
    for p in partial_point.iter().rev().skip(1) {
        res = fix_last_variable_no_par(&res, p);
    }
    res
}

pub fn fix_last_variable_no_par<F: PrimeField>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &F,
) -> DenseMultilinearExtension<F> {
    let nv = poly.num_vars();
    let half_len = 1 << (nv - 1);
    let mut res = vec![F::zero(); half_len];
    let one_minus_p = F::one() - partial_point;
    for (i, e) in res.iter_mut().enumerate().take(half_len) {
        *e = one_minus_p * poly.evaluations[i] + *partial_point * poly.evaluations[i + half_len];
    }
    DenseMultilinearExtension::from_evaluations_vec(nv - 1, res)
}
pub fn fix_last_variables<F: PrimeField>(
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
    for (i, point) in partial_point.iter().rev().enumerate().take(dim) {
        poly = fix_last_variable_helper(&poly, nv - i, point);
    }

    DenseMultilinearExtension::<F>::from_evaluations_slice(nv - dim, &poly[..(1 << (nv - dim))])
}

pub fn fix_last_variable<F: PrimeField>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &F,
) -> DenseMultilinearExtension<F> {
    assert!(poly.num_vars != 0, "invalid size of partial point");

    let nv = poly.num_vars;
    let res = fix_last_variable_helper(&poly.evaluations, nv, partial_point);
    DenseMultilinearExtension::<F>::from_evaluations_slice(nv - 1, &res)
}

fn fix_last_variable_helper<F: Field>(data: &[F], nv: usize, point: &F) -> Vec<F> {
    let one_minus_p = F::one() - point;
    let half_len = 1 << (nv - 1);
    let mut res = vec![F::zero(); half_len];

    // evaluate single variable of partial point from left to right
    #[cfg(not(feature = "parallel"))]
    for b in 0..half_len {
        res[b] = data[b] * one_minus_p + data[b + half_len] * point;
    }

    #[cfg(feature = "parallel")]
    if nv >= 13 {
        // on my computer we parallelization doesn't help till nv >= 13
        res.par_iter_mut().enumerate().for_each(|(i, x)| {
            *x = data[i] * one_minus_p + data[i + half_len] * point;
        });
    } else {
        for b in 0..(1 << (nv - 1)) {
            res[b] = data[b] * one_minus_p + data[b + half_len] * point;
        }
    }

    res
}
