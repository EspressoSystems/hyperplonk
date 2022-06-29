//! Useful utilities for KZG PCS

use crate::PCSErrors;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseMultilinearExtension, EvaluationDomain, Evaluations,
    MultilinearExtension, Polynomial, Radix2EvaluationDomain,
};
use ark_std::{end_timer, start_timer};

/// Compute W \circ l.
///
/// Given an MLE W, and a list of univariate polynomials l, evaluate W at l.
///
/// Returns an error if l's length does not matches number of variables in W.
pub(crate) fn compute_w_circ_l<F: PrimeField>(
    w: &DenseMultilinearExtension<F>,
    l: &[DensePolynomial<F>],
) -> Result<DensePolynomial<F>, PCSErrors> {
    let timer = start_timer!(|| "compute W \\circ l");

    if w.num_vars != l.len() {
        return Err(PCSErrors::InvalidParameters(format!(
            "l's length ({}) does not match num_variables ({})",
            l.len(),
            w.num_vars(),
        )));
    }

    let mut res_eval: Vec<F> = vec![];
    let uni_degree = l.len() + 1;

    let domain = match Radix2EvaluationDomain::<F>::new(uni_degree) {
        Some(p) => p,
        None => {
            return Err(PCSErrors::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };

    for point in domain.elements() {
        let l_eval: Vec<F> = l.iter().map(|x| x.evaluate(&point)).collect();
        res_eval.push(w.evaluate(l_eval.as_ref()).unwrap())
    }
    let evaluation = Evaluations::from_vec_and_domain(res_eval, domain);
    let res = evaluation.interpolate();

    end_timer!(timer);
    Ok(res)
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_poly::UVPolynomial;

    #[test]
    fn test_w_circ_l() -> Result<(), PCSErrors> {
        test_w_circ_l_helper::<Fr>()
    }

    fn test_w_circ_l_helper<F: PrimeField>() -> Result<(), PCSErrors> {
        {
            // Example from page 53:
            // W = 3x1x2 + 2x2 whose evaluations are
            // 0, 0 |-> 0
            // 0, 1 |-> 2
            // 1, 0 |-> 0
            // 1, 1 |-> 5
            let w_eval = vec![F::zero(), F::from(2u64), F::zero(), F::from(5u64)];
            let w = DenseMultilinearExtension::from_evaluations_vec(2, w_eval);

            // l0 =   t + 2
            // l1 = -2t + 4
            let l0 = DensePolynomial::from_coefficients_vec(vec![F::from(2u64), F::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![F::from(4u64), -F::from(2u64)]);

            // res = -6t^2 - 4t + 32
            let res = compute_w_circ_l(&w, [l1, l0].as_ref())?;
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
            // 0, 0, 1 |->  7
            // 0, 1, 0 |-> -6
            // 0, 1, 1 |->  4
            // 1, 0, 0 |->  5
            // 1, 0, 1 |->  8
            // 1, 1, 0 |-> -3
            // 1, 1, 1 |->  4
            let w_eval = vec![
                F::zero(),
                F::from(7u64),
                -F::from(6u64),
                F::from(4u64),
                F::from(5u64),
                F::from(8u64),
                -F::from(3u64),
                F::from(4u64),
            ];
            let w = DenseMultilinearExtension::from_evaluations_vec(3, w_eval);

            // l0 =   t + 2
            // l1 =  3t - 4
            // l2 = -5t + 6
            let l0 = DensePolynomial::from_coefficients_vec(vec![F::from(2u64), F::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![-F::from(4u64), F::from(3u64)]);
            let l2 = DensePolynomial::from_coefficients_vec(vec![F::from(6u64), -F::from(5u64)]);
            let res = compute_w_circ_l(&w, [l2, l1, l0].as_ref())?;

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
}
