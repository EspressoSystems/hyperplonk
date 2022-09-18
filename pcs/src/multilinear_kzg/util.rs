// Copyright (c) 2022 Espresso Systems (espressosys.com)
// This file is part of the Jellyfish library.

// You should have received a copy of the MIT License
// along with the Jellyfish library. If not, see <https://mit-license.org/>.

//! Useful utilities for KZG PCS
use crate::prelude::PCSError;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseMultilinearExtension, EvaluationDomain, Evaluations,
    MultilinearExtension, Polynomial, Radix2EvaluationDomain,
};
use ark_std::{end_timer, format, log2, rc::Rc, start_timer, string::ToString, vec, vec::Vec};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

/// Decompose an integer into a binary vector in little endian.
#[allow(dead_code)]
pub(crate) fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push(i & 1 == 1);
        i >>= 1;
    }
    res
}

/// given the evaluation input `point` of the `index`-th polynomial,
/// obtain the evaluation point in the merged polynomial
pub(crate) fn gen_eval_point<F: PrimeField>(index: usize, index_len: usize, point: &[F]) -> Vec<F> {
    let mut index_vec: Vec<F> = bit_decompose(index as u64, index_len)
        .into_iter()
        .map(|x| F::from(x))
        .collect();
    index_vec.reverse();
    [point, &index_vec].concat()
}

/// For an MLE w with `mle_num_vars` variables, and `point_len` number of
/// points, compute the degree of the univariate polynomial `q(x):= w(l(x))`
/// where l(x) is a list of polynomials that go through all points.
// uni_degree is computed as `mle_num_vars * point_len`:
// - each l(x) is of degree `point_len`
// - mle has degree one
// - worst case is `\prod_{i=0}^{mle_num_vars-1} l_i(x) < point_len * mle_num_vars`
#[inline]
pub fn compute_qx_degree(mle_num_vars: usize, point_len: usize) -> usize {
    mle_num_vars * (1 << log2(point_len))
}

/// get the domain for the univariate polynomial
#[inline]
pub(crate) fn get_uni_domain<F: PrimeField>(
    uni_poly_degree: usize,
) -> Result<Radix2EvaluationDomain<F>, PCSError> {
    let domain = match Radix2EvaluationDomain::<F>::new(uni_poly_degree) {
        Some(p) => p,
        None => {
            return Err(PCSError::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };
    Ok(domain)
}

/// Compute W \circ l.
///
/// Given an MLE W, and a list of univariate polynomials l, generate the
/// univariate polynomial that composes W with l.
///
/// Returns an error if l's length does not matches number of variables in W.
pub(crate) fn compute_w_circ_l_with_prefix<F: PrimeField>(
    w: &DenseMultilinearExtension<F>,
    l: &[DensePolynomial<F>],
    num_points: usize,
) -> Result<DensePolynomial<F>, PCSError> {
    let timer = start_timer!(|| "compute W \\circ l");

    if w.num_vars != l.len() {
        return Err(PCSError::InvalidParameters(format!(
            "l's length ({}) does not match num_variables ({})",
            l.len(),
            w.num_vars(),
        )));
    }

    let uni_degree = (w.num_vars + log2(num_points) as usize) * num_points;

    let domain = match Radix2EvaluationDomain::<F>::new(uni_degree) {
        Some(p) => p,
        None => {
            return Err(PCSError::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };
    let res_eval = (0..domain.size())
        .into_par_iter()
        .map(|i| {
            let l_eval: Vec<F> = l.iter().map(|x| x.evaluate(&domain.element(i))).collect();
            w.evaluate(l_eval.as_ref()).unwrap()
        })
        .collect();
    let evaluation = Evaluations::from_vec_and_domain(res_eval, domain);
    let res = evaluation.interpolate();

    end_timer!(timer);
    Ok(res)
}

/// Compute W \circ l.
///
/// Given an MLE W, and a list of univariate polynomials l, generate the
/// univariate polynomial that composes W with l.
///
/// Returns an error if l's length does not matches number of variables in W.
pub(crate) fn compute_w_circ_l<F: PrimeField>(
    w: &DenseMultilinearExtension<F>,
    l: &[DensePolynomial<F>],
    num_points: usize,
) -> Result<DensePolynomial<F>, PCSError> {
    let timer = start_timer!(|| "compute W \\circ l");

    if w.num_vars != l.len() {
        return Err(PCSError::InvalidParameters(format!(
            "l's length ({}) does not match num_variables ({})",
            l.len(),
            w.num_vars(),
        )));
    }

    let uni_degree = compute_qx_degree(w.num_vars(), num_points);

    let domain = match Radix2EvaluationDomain::<F>::new(uni_degree) {
        Some(p) => p,
        None => {
            return Err(PCSError::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };
    let step = start_timer!(|| "compute eval");

    let res_eval = (0..domain.size())
        .into_par_iter()
        .map(|i| {
            let l_eval: Vec<F> = l.iter().map(|x| x.evaluate(&domain.element(i))).collect();
            w.evaluate(l_eval.as_ref()).unwrap()
        })
        .collect();

    end_timer!(step);
    let step = start_timer!(|| "interpolation");

    let evaluation = Evaluations::from_vec_and_domain(res_eval, domain);
    let res = evaluation.interpolate();
    end_timer!(step);
    end_timer!(timer);
    Ok(res)
}

/// Return the number of variables that one need for an MLE to
/// batch the list of MLEs
#[inline]
pub fn get_batched_nv(num_var: usize, polynomials_len: usize) -> usize {
    num_var + log2(polynomials_len) as usize
}

/// merge a set of polynomials. Returns an error if the
/// polynomials do not share a same number of nvs.
pub fn merge_polynomials<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
) -> Result<Rc<DenseMultilinearExtension<F>>, PCSError> {
    let nv = polynomials[0].num_vars();
    for poly in polynomials.iter() {
        if nv != poly.num_vars() {
            return Err(PCSError::InvalidParameters(
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

/// Given a list of points, build `l(points)` which is a list of univariate
/// polynomials that goes through the points; extend the dimension of the points
/// by `log(points.len())`.
pub(crate) fn build_l_with_prefix<F: PrimeField>(
    points: &[Vec<F>],
    domain: &Radix2EvaluationDomain<F>,
) -> Result<Vec<DensePolynomial<F>>, PCSError> {
    let prefix_len = log2(points.len()) as usize;
    let mut uni_polys = Vec::new();

    println!("with prefix domain size {}", domain.size());
    // 1.1 build the indexes and the univariate polys that go through the indexes
    let indexes: Vec<Vec<bool>> = (0..points.len())
        .map(|x| bit_decompose(x as u64, prefix_len))
        .collect();
    for i in 0..prefix_len {
        let eval: Vec<F> = indexes
            .iter()
            .map(|x| F::from(x[prefix_len - i - 1]))
            .collect();

        uni_polys.push(Evaluations::from_vec_and_domain(eval, *domain).interpolate());
    }

    // 1.2 build the actual univariate polys that go through the points
    uni_polys.extend_from_slice(build_l(points, domain)?.as_slice());

    Ok(uni_polys)
}

/// Given a list of points, build `l(points)` which is a list of univariate
/// polynomials that goes through the points.
pub(crate) fn build_l<F: PrimeField>(
    points: &[Vec<F>],
    domain: &Radix2EvaluationDomain<F>,
) -> Result<Vec<DensePolynomial<F>>, PCSError> {
    let mut uni_polys = Vec::new();
    let num_var = points[0].len();
    // build the actual univariate polys that go through the points
    for i in 0..num_var {
        let mut eval: Vec<F> = points.iter().map(|x| x[i]).collect();
        eval.extend_from_slice(vec![F::zero(); domain.size as usize - eval.len()].as_slice());
        uni_polys.push(Evaluations::from_vec_and_domain(eval, *domain).interpolate())
    }
    Ok(uni_polys)
}

/// Input a list of multilinear polynomials and a list of points,
/// generate a list of evaluations.
// Note that this function is only used for testing verifications.
// In practice verifier does not see polynomials, and the `mle_values`
// are included in the `batch_proof`.
#[cfg(test)]
pub(crate) fn generate_evaluations<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
    points: &[Vec<F>],
) -> Result<Vec<F>, PCSError> {
    if polynomials.len() != points.len() {
        return Err(PCSError::InvalidParameters(
            "polynomial length does not match point length".to_string(),
        ));
    }
    let uni_poly_degree = points.len();
    let merge_poly = merge_polynomials(polynomials)?;

    let domain = get_uni_domain::<F>(uni_poly_degree)?;
    let uni_polys = build_l_with_prefix(points, &domain)?;
    let mut mle_values = vec![];

    for i in 0..uni_poly_degree {
        let point: Vec<F> = uni_polys
            .iter()
            .map(|poly| poly.evaluate(&domain.element(i)))
            .collect();

        let mle_value = merge_poly.evaluate(&point).unwrap();
        mle_values.push(mle_value)
    }
    Ok(mle_values)
}

/// Input a list of multilinear polynomials and a list of points,
/// generate a list of evaluations.
// Note that this function is only used for testing verifications.
// In practice verifier does not see polynomials, and the `mle_values`
// are included in the `batch_proof`.
#[cfg(test)]
pub(crate) fn generate_evaluations_single_poly<F: PrimeField>(
    polynomial: &Rc<DenseMultilinearExtension<F>>,
    points: &[Vec<F>],
) -> Result<Vec<F>, PCSError> {
    let uni_poly_degree = points.len();

    let domain = get_uni_domain::<F>(uni_poly_degree)?;
    let uni_polys = build_l(points, &domain)?;
    let mut mle_values = vec![];

    for i in 0..uni_poly_degree {
        let point: Vec<F> = uni_polys
            .iter()
            .map(|poly| poly.evaluate(&domain.element(i)))
            .collect();

        let mle_value = polynomial.evaluate(&point).unwrap();
        mle_values.push(mle_value)
    }
    Ok(mle_values)
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_ff::field_new;
    use ark_poly::UVPolynomial;
    use ark_std::{One, Zero};

    #[test]
    fn test_w_circ_l() -> Result<(), PCSError> {
        test_w_circ_l_helper::<Fr>()
    }

    fn test_w_circ_l_helper<F: PrimeField>() -> Result<(), PCSError> {
        {
            // Example from page 53:
            // W = 3x1x2 + 2x2 whose evaluations are
            // 0, 0 |-> 0
            // 1, 0 |-> 0
            // 0, 1 |-> 2
            // 1, 1 |-> 5
            let w_eval = vec![F::zero(), F::zero(), F::from(2u64), F::from(5u64)];
            let w = DenseMultilinearExtension::from_evaluations_vec(2, w_eval);

            // l0 =   t + 2
            // l1 = -2t + 4
            let l0 = DensePolynomial::from_coefficients_vec(vec![F::from(2u64), F::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![F::from(4u64), -F::from(2u64)]);

            // res = -6t^2 - 4t + 32
            let res = compute_w_circ_l(&w, [l0, l1].as_ref(), 4)?;
            let res_rec = DensePolynomial::from_coefficients_vec(vec![
                F::from(32u64),
                -F::from(4u64),
                -F::from(6u64),
            ]);
            assert_eq!(res, res_rec);
        }
        {
            // A random example
            // W = x1x2x3 - 2x1x2 + 3x2x3 - 4x1x3 + 5x1 - 6x2 + 7x3
            // 0, 0, 0 |->  0
            // 1, 0, 0 |->  5
            // 0, 1, 0 |-> -6
            // 1, 1, 0 |-> -3
            // 0, 0, 1 |->  7
            // 1, 0, 1 |->  8
            // 0, 1, 1 |->  4
            // 1, 1, 1 |->  4
            let w_eval = vec![
                F::zero(),
                F::from(5u64),
                -F::from(6u64),
                -F::from(3u64),
                F::from(7u64),
                F::from(8u64),
                F::from(4u64),
                F::from(4u64),
            ];
            let w = DenseMultilinearExtension::from_evaluations_vec(3, w_eval);

            // l0 =   t + 2
            // l1 =  3t - 4
            // l2 = -5t + 6
            let l0 = DensePolynomial::from_coefficients_vec(vec![F::from(2u64), F::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![-F::from(4u64), F::from(3u64)]);
            let l2 = DensePolynomial::from_coefficients_vec(vec![F::from(6u64), -F::from(5u64)]);
            let res = compute_w_circ_l(&w, [l0, l1, l2].as_ref(), 8)?;

            // res = -15t^3 - 23t^2 + 130t - 76
            let res_rec = DensePolynomial::from_coefficients_vec(vec![
                -F::from(76u64),
                F::from(130u64),
                -F::from(23u64),
                -F::from(15u64),
            ]);

            assert_eq!(res, res_rec);
        }
        Ok(())
    }

    #[test]
    fn test_w_circ_l_with_prefix() -> Result<(), PCSError> {
        test_w_circ_l_with_prefix_helper::<Fr>()
    }

    fn test_w_circ_l_with_prefix_helper<F: PrimeField>() -> Result<(), PCSError> {
        {
            // Example from page 53:
            // W = 3x1x2 + 2x2 whose evaluations are
            // 0, 0 |-> 0
            // 1, 0 |-> 0
            // 0, 1 |-> 2
            // 1, 1 |-> 5
            let w_eval = vec![F::zero(), F::zero(), F::from(2u64), F::from(5u64)];
            let w = DenseMultilinearExtension::from_evaluations_vec(2, w_eval);

            // l0 =   t + 2
            // l1 = -2t + 4
            let l0 = DensePolynomial::from_coefficients_vec(vec![F::from(2u64), F::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![F::from(4u64), -F::from(2u64)]);

            // res = -6t^2 - 4t + 32
            let res = compute_w_circ_l_with_prefix(&w, [l0, l1].as_ref(), 4)?;
            let res_rec = DensePolynomial::from_coefficients_vec(vec![
                F::from(32u64),
                -F::from(4u64),
                -F::from(6u64),
            ]);
            assert_eq!(res, res_rec);
        }
        {
            // A random example
            // W = x1x2x3 - 2x1x2 + 3x2x3 - 4x1x3 + 5x1 - 6x2 + 7x3
            // 0, 0, 0 |->  0
            // 1, 0, 0 |->  5
            // 0, 1, 0 |-> -6
            // 1, 1, 0 |-> -3
            // 0, 0, 1 |->  7
            // 1, 0, 1 |->  8
            // 0, 1, 1 |->  4
            // 1, 1, 1 |->  4
            let w_eval = vec![
                F::zero(),
                F::from(5u64),
                -F::from(6u64),
                -F::from(3u64),
                F::from(7u64),
                F::from(8u64),
                F::from(4u64),
                F::from(4u64),
            ];
            let w = DenseMultilinearExtension::from_evaluations_vec(3, w_eval);

            // l0 =   t + 2
            // l1 =  3t - 4
            // l2 = -5t + 6
            let l0 = DensePolynomial::from_coefficients_vec(vec![F::from(2u64), F::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![-F::from(4u64), F::from(3u64)]);
            let l2 = DensePolynomial::from_coefficients_vec(vec![F::from(6u64), -F::from(5u64)]);
            let res = compute_w_circ_l_with_prefix(&w, [l0, l1, l2].as_ref(), 8)?;

            // res = -15t^3 - 23t^2 + 130t - 76
            let res_rec = DensePolynomial::from_coefficients_vec(vec![
                -F::from(76u64),
                F::from(130u64),
                -F::from(23u64),
                -F::from(15u64),
            ]);

            assert_eq!(res, res_rec);
        }
        Ok(())
    }

    #[test]
    fn test_build_l_with_prefix() -> Result<(), PCSError> {
        test_build_l_with_prefix_helper::<Fr>()
    }

    fn test_build_l_with_prefix_helper<F: PrimeField>() -> Result<(), PCSError> {
        // point 1 is [1, 2]
        let point1 = vec![Fr::from(1u64), Fr::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![Fr::from(3u64), Fr::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![Fr::from(5u64), Fr::from(6u64)];

        {
            let domain = get_uni_domain::<Fr>(2)?;
            let l = build_l_with_prefix(&[point1.clone(), point2.clone()], &domain)?;

            // roots: [1, -1]
            // l0 = -1/2 * x + 1/2
            // l1 = -x + 2
            // l2 = -x + 3
            let l0 = DensePolynomial::from_coefficients_vec(vec![
                Fr::one() / Fr::from(2u64),
                -Fr::one() / Fr::from(2u64),
            ]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![Fr::from(2u64), -Fr::one()]);
            let l2 = DensePolynomial::from_coefficients_vec(vec![Fr::from(3u64), -Fr::one()]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
            assert_eq!(l2, l[2], "l2 not equal");
        }

        {
            let domain = get_uni_domain::<Fr>(3)?;
            let l = build_l_with_prefix(&[point1, point2, point3], &domain)?;

            // sage: q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
            // sage: P.<x> = PolynomialRing(Zmod(q))
            // sage: root1 = 1
            // sage: root2 = 0x8D51CCCE760304D0EC030002760300000001000000000000
            // sage: root3 = -1
            // sage: root4 = -root2
            // Arkwork's code is a bit wired: it also interpolate (root4, 0)
            // which returns a degree 3 polynomial, instead of degree 2

            // ========================
            // l0: [0, 0, 1]
            // ========================
            // sage: points = [(root1, 0), (root2, 0), (root3, 1), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 13108968793781547619861935127046491459422638125131909455650914674984645296128*x^3 +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888385*x^2 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296128*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888385
            let l0 = DensePolynomial::from_coefficients_vec(vec![
                field_new!(
                    Fr,
                    "39326906381344642859585805381139474378267914375395728366952744024953935888385"
                ),
                field_new!(
                    Fr,
                    "13108968793781547619861935127046491459422638125131909455650914674984645296128"
                ),
                field_new!(
                    Fr,
                    "39326906381344642859585805381139474378267914375395728366952744024953935888385"
                ),
                field_new!(
                    Fr,
                    "13108968793781547619861935127046491459422638125131909455650914674984645296128"
                ),
            ]);

            // ========================
            // l1: [0, 1, 0]
            // ========================
            // sage: points = [(root1, 0), (root2, 1), (root3, 0), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 866286206518413079694067382671935694567563117191340490752*x^3 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296128*x^2 +
            // 52435875175126190478581454301667552757996485117855702128036095582747240693761*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888385
            let l1 = DensePolynomial::from_coefficients_vec(vec![
                field_new!(
                    Fr,
                    "39326906381344642859585805381139474378267914375395728366952744024953935888385"
                ),
                field_new!(
                    Fr,
                    "52435875175126190478581454301667552757996485117855702128036095582747240693761"
                ),
                field_new!(
                    Fr,
                    "13108968793781547619861935127046491459422638125131909455650914674984645296128"
                ),
                field_new!(
                    Fr,
                    "866286206518413079694067382671935694567563117191340490752"
                ),
            ]);

            // ========================
            // l2: [1, 3, 5]
            // ========================
            // sage: points = [(root1, 1), (root2, 3), (root3, 5), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 2598858619555239239082202148015807083702689351574021472255*x^3 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296129*x^2 +
            // 52435875175126190476848881888630726598608350352511830738900969348364559712256*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888387
            let l2 = DensePolynomial::from_coefficients_vec(vec![
                field_new!(
                    Fr,
                    "39326906381344642859585805381139474378267914375395728366952744024953935888387"
                ),
                field_new!(
                    Fr,
                    "52435875175126190476848881888630726598608350352511830738900969348364559712256"
                ),
                field_new!(
                    Fr,
                    "13108968793781547619861935127046491459422638125131909455650914674984645296129"
                ),
                field_new!(
                    Fr,
                    "2598858619555239239082202148015807083702689351574021472255"
                ),
            ]);

            // ========================
            // l3: [2, 4, 6]
            // ========================
            // sage: points = [(root1, 2), (root2, 4), (root3, 6), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 3465144826073652318776269530687742778270252468765361963007*x^3 +
            // x^2 +
            // 52435875175126190475982595682112313518914282969839895044333406231173219221504*x +
            // 3
            let l3 = DensePolynomial::from_coefficients_vec(vec![
                Fr::from(3u64),
                field_new!(
                    Fr,
                    "52435875175126190475982595682112313518914282969839895044333406231173219221504"
                ),
                Fr::one(),
                field_new!(
                    Fr,
                    "3465144826073652318776269530687742778270252468765361963007"
                ),
            ]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
            assert_eq!(l2, l[2], "l2 not equal");
            assert_eq!(l3, l[3], "l3 not equal");
        }
        Ok(())
    }

    #[test]
    fn test_build_l() -> Result<(), PCSError> {
        test_build_l_helper::<Fr>()
    }

    fn test_build_l_helper<F: PrimeField>() -> Result<(), PCSError> {
        // point 1 is [1, 2]
        let point1 = vec![Fr::from(1u64), Fr::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![Fr::from(3u64), Fr::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![Fr::from(5u64), Fr::from(6u64)];

        {
            let domain = get_uni_domain::<Fr>(2)?;
            let l = build_l(&[point1.clone(), point2.clone()], &domain)?;

            // roots: [1, -1]
            // l0 = -x + 2
            // l1 = -x + 3
            let l0 = DensePolynomial::from_coefficients_vec(vec![Fr::from(2u64), -Fr::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![Fr::from(3u64), -Fr::one()]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
        }

        {
            let domain = get_uni_domain::<Fr>(3)?;
            let l = build_l(&[point1, point2, point3], &domain)?;

            // sage: q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
            // sage: P.<x> = PolynomialRing(Zmod(q))
            // sage: root1 = 1
            // sage: root2 = 0x8D51CCCE760304D0EC030002760300000001000000000000
            // sage: root3 = -1
            // sage: root4 = -root2
            // Arkwork's code is a bit wired: it also interpolate (root4, 0)
            // which returns a degree 3 polynomial, instead of degree 2

            // ========================
            // l0: [1, 3, 5]
            // ========================
            // sage: points = [(root1, 1), (root2, 3), (root3, 5), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 2598858619555239239082202148015807083702689351574021472255*x^3 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296129*x^2 +
            // 52435875175126190476848881888630726598608350352511830738900969348364559712256*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888387
            let l0 = DensePolynomial::from_coefficients_vec(vec![
                field_new!(
                    Fr,
                    "39326906381344642859585805381139474378267914375395728366952744024953935888387"
                ),
                field_new!(
                    Fr,
                    "52435875175126190476848881888630726598608350352511830738900969348364559712256"
                ),
                field_new!(
                    Fr,
                    "13108968793781547619861935127046491459422638125131909455650914674984645296129"
                ),
                field_new!(
                    Fr,
                    "2598858619555239239082202148015807083702689351574021472255"
                ),
            ]);

            // ========================
            // l1: [2, 4, 6]
            // ========================
            // sage: points = [(root1, 2), (root2, 4), (root3, 6), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 3465144826073652318776269530687742778270252468765361963007*x^3 +
            // x^2 +
            // 52435875175126190475982595682112313518914282969839895044333406231173219221504*x +
            // 3
            let l1 = DensePolynomial::from_coefficients_vec(vec![
                Fr::from(3u64),
                field_new!(
                    Fr,
                    "52435875175126190475982595682112313518914282969839895044333406231173219221504"
                ),
                Fr::one(),
                field_new!(
                    Fr,
                    "3465144826073652318776269530687742778270252468765361963007"
                ),
            ]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
        }
        Ok(())
    }

    #[test]
    fn test_qx() -> Result<(), PCSError> {
        // Example from page 53:
        // W1 = 3x1x2 + 2x2
        let w_eval = vec![Fr::zero(), Fr::from(2u64), Fr::zero(), Fr::from(5u64)];
        let w = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        let r = Fr::from(42u64);

        // point 1 is [1, 2]
        let point1 = vec![Fr::from(1u64), Fr::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![Fr::from(3u64), Fr::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![Fr::from(5u64), Fr::from(6u64)];

        {
            let domain = get_uni_domain::<Fr>(2)?;
            let l = build_l(&[point1.clone(), point2.clone()], &domain)?;

            let q_x = compute_w_circ_l(&w, &l, 2)?;

            let point: Vec<Fr> = l.iter().map(|poly| poly.evaluate(&r)).collect();

            assert_eq!(
                q_x.evaluate(&r),
                w.evaluate(&point).unwrap(),
                "q(r) != w(l(r))"
            );
        }

        {
            let domain = get_uni_domain::<Fr>(3)?;

            let l = build_l(&[point1, point2, point3], &domain)?;
            let q_x = compute_w_circ_l(&w, &l, 3)?;

            let point: Vec<Fr> = vec![l[0].evaluate(&r), l[1].evaluate(&r)];

            assert_eq!(
                q_x.evaluate(&r),
                w.evaluate(&point).unwrap(),
                "q(r) != w(l(r))"
            );
        }
        Ok(())
    }

    #[test]
    fn test_qx_with_prefix() -> Result<(), PCSError> {
        // Example from page 53:
        // W1 = 3x1x2 + 2x2
        let w_eval = vec![Fr::zero(), Fr::from(2u64), Fr::zero(), Fr::from(5u64)];
        let w1 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        // W2 = x1x2 + x1
        let w_eval = vec![Fr::zero(), Fr::zero(), Fr::from(1u64), Fr::from(2u64)];
        let w2 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        // W3 = x1 + x2
        let w_eval = vec![Fr::zero(), Fr::one(), Fr::from(1u64), Fr::from(2u64)];
        let w3 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        let r = Fr::from(42u64);

        // point 1 is [1, 2]
        let point1 = vec![Fr::from(1u64), Fr::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![Fr::from(3u64), Fr::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![Fr::from(5u64), Fr::from(6u64)];

        {
            let domain = get_uni_domain::<Fr>(2)?;
            // w = (3x1x2 + 2x2)(1-x0) + (x1x2 + x1)x0
            // with evaluations: [0,2,0,5,0,0,1,2]
            let w = merge_polynomials(&[w1.clone(), w2.clone()])?;

            let l = build_l_with_prefix(&[point1.clone(), point2.clone()], &domain)?;

            // sage: P.<x> = PolynomialRing(ZZ)
            // sage: l0 = -1/2 * x + 1/2
            // sage: l1 = -x + 2
            // sage: l2 = -x + 3
            // sage: w = (3 * l1 * l2 + 2 * l2) * (1-l0) + (l1 * l2 + l1) * l0
            // sage: w
            // x^3 - 7/2*x^2 - 7/2*x + 16
            //
            // q(x) = x^3 - 7/2*x^2 - 7/2*x + 16
            let q_x = compute_w_circ_l_with_prefix(&w, &l, 2)?;

            let point: Vec<Fr> = l.iter().map(|poly| poly.evaluate(&r)).collect();

            assert_eq!(
                q_x.evaluate(&r),
                w.evaluate(&point).unwrap(),
                "q(r) != w(l(r))"
            );
        }

        {
            let domain = get_uni_domain::<Fr>(3)?;
            let w = merge_polynomials(&[w1, w2, w3])?;

            let l = build_l_with_prefix(&[point1, point2, point3], &domain)?;
            let q_x = compute_w_circ_l_with_prefix(&w, &l, 3)?;

            let point: Vec<Fr> = vec![
                l[0].evaluate(&r),
                l[1].evaluate(&r),
                l[2].evaluate(&r),
                l[3].evaluate(&r),
            ];

            assert_eq!(
                q_x.evaluate(&r),
                w.evaluate(&point).unwrap(),
                "q(r) != w(l(r))"
            );
        }
        Ok(())
    }
}
