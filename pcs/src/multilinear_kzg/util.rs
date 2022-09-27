// Copyright (c) 2022 Espresso Systems (espressosys.com)
// This file is part of the Jellyfish library.

// You should have received a copy of the MIT License
// along with the Jellyfish library. If not, see <https://mit-license.org/>.

//! Useful utilities for KZG PCS
use crate::prelude::PCSError;
use arithmetic::evaluate_no_par;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseMultilinearExtension, EvaluationDomain, Evaluations,
    MultilinearExtension, Polynomial, Radix2EvaluationDomain,
};
use ark_std::{end_timer, format, log2, start_timer, string::ToString, vec::Vec};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

/// For an MLE w with `mle_num_vars` variables, and `point_len` number of
/// points, compute the degree of the univariate polynomial `q(x):= w(l(x))`
/// where l(x) is a list of polynomials that go through all points.
// uni_degree is computed as `mle_num_vars * point_len`:
// - each l(x) is of degree `point_len`
// - mle has degree one
// - worst case is `\prod_{i=0}^{mle_num_vars-1} l_i(x) < point_len * mle_num_vars`
#[inline]
pub fn compute_qx_degree(mle_num_vars: usize, point_len: usize) -> usize {
    mle_num_vars * ((1 << log2(point_len)) - 1) + 1
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
    with_suffix: bool,
) -> Result<DensePolynomial<F>, PCSError> {
    let timer = start_timer!(|| "compute W \\circ l");

    if w.num_vars != l.len() {
        return Err(PCSError::InvalidParameters(format!(
            "l's length ({}) does not match num_variables ({})",
            l.len(),
            w.num_vars(),
        )));
    }
    let uni_degree = if with_suffix {
        compute_qx_degree(w.num_vars() + log2(num_points) as usize, num_points)
    } else {
        compute_qx_degree(w.num_vars(), num_points)
    };

    let domain = match Radix2EvaluationDomain::<F>::new(uni_degree) {
        Some(p) => p,
        None => {
            return Err(PCSError::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };

    let step = start_timer!(|| format!("compute eval {}-dim domain", domain.size()));
    let res_eval = (0..domain.size())
        .into_par_iter()
        .map(|i| {
            let l_eval: Vec<F> = l.iter().map(|x| x.evaluate(&domain.element(i))).collect();
            evaluate_no_par(w, &l_eval)
        })
        .collect();
    end_timer!(step);

    let evaluation = Evaluations::from_vec_and_domain(res_eval, domain);
    let res = evaluation.interpolate();

    end_timer!(timer);
    Ok(res)
}

/// Input a list of multilinear polynomials and a list of points,
/// generate a list of evaluations.
// Note that this function is only used for testing verifications.
// In practice verifier does not see polynomials, and the `mle_values`
// are included in the `batch_proof`.
#[cfg(test)]
pub(crate) fn generate_evaluations_multi_poly<F: PrimeField>(
    polynomials: &[std::rc::Rc<DenseMultilinearExtension<F>>],
    points: &[Vec<F>],
) -> Result<Vec<F>, PCSError> {
    use arithmetic::{build_l, get_uni_domain, merge_polynomials};

    if polynomials.len() != points.len() {
        return Err(PCSError::InvalidParameters(
            "polynomial length does not match point length".to_string(),
        ));
    }
    let uni_poly_degree = points.len();
    let merge_poly = merge_polynomials(polynomials)?;

    let domain = get_uni_domain::<F>(uni_poly_degree)?;
    let uni_polys = build_l(points, &domain, true)?;
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
    polynomial: &std::rc::Rc<DenseMultilinearExtension<F>>,
    points: &[Vec<F>],
) -> Result<Vec<F>, PCSError> {
    use arithmetic::{build_l, get_uni_domain};

    let uni_poly_degree = points.len();

    let domain = get_uni_domain::<F>(uni_poly_degree)?;
    let uni_polys = build_l(points, &domain, false)?;
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
    use arithmetic::{build_l, get_uni_domain, merge_polynomials};
    use ark_bls12_381::Fr;
    use ark_poly::UVPolynomial;
    use ark_std::{One, Zero};
    use std::rc::Rc;

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
            let res = compute_w_circ_l(&w, [l0, l1].as_ref(), 4, false)?;
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
            let res = compute_w_circ_l(&w, [l0, l1, l2].as_ref(), 8, false)?;

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
            let res = compute_w_circ_l(&w, [l0, l1].as_ref(), 4, true)?;
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
            let res = compute_w_circ_l(&w, [l0, l1, l2].as_ref(), 8, true)?;

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
            let l = build_l(&[point1.clone(), point2.clone()], &domain, false)?;

            let q_x = compute_w_circ_l(&w, &l, 2, false)?;

            let point: Vec<Fr> = l.iter().map(|poly| poly.evaluate(&r)).collect();

            assert_eq!(
                q_x.evaluate(&r),
                w.evaluate(&point).unwrap(),
                "q(r) != w(l(r))"
            );
        }

        {
            let domain = get_uni_domain::<Fr>(3)?;

            let l = build_l(&[point1, point2, point3], &domain, false)?;
            let q_x = compute_w_circ_l(&w, &l, 3, false)?;

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

            let l = build_l(&[point1.clone(), point2.clone()], &domain, true)?;

            // sage: P.<x> = PolynomialRing(ZZ)
            // sage: l0 = -1/2 * x + 1/2
            // sage: l1 = -x + 2
            // sage: l2 = -x + 3
            // sage: w = (3 * l1 * l2 + 2 * l2) * (1-l0) + (l1 * l2 + l1) * l0
            // sage: w
            // x^3 - 7/2*x^2 - 7/2*x + 16
            //
            // q(x) = x^3 - 7/2*x^2 - 7/2*x + 16
            let q_x = compute_w_circ_l(&w, &l, 2, true)?;

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

            let l = build_l(&[point1, point2, point3], &domain, true)?;
            let q_x = compute_w_circ_l(&w, &l, 3, true)?;

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
