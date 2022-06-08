//! This module implements useful functions for the permutation check protocol.

use crate::{utils::get_index, PolyIOPErrors, VirtualPolynomial};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use std::rc::Rc;

/// Compute `Q(x)` which a virtual polynomial of degree 2 defined as
///
/// Q(x) := prod(1,x) - prod(x, 0) - prod(x, 1)
///       + alpha * (
///             (w(x) + beta * s_perm(x) + gamma) * prod(0, x)
///           - (w(x) + beta * s_id(x)   + gamma))
///
/// Returns an error when the num_vars in w/s_id/s_perm does not match
///
/// Cost: linear in N.
#[allow(dead_code)]
// TODO: remove
pub(super) fn build_q_x<F: PrimeField>(
    alpha: &F,
    beta: &F,
    gamma: &F,
    w: &DenseMultilinearExtension<F>,
    s_id: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
) -> Result<VirtualPolynomial<F>, PolyIOPErrors> {
    let start = start_timer!(|| "compute Q(x)");

    let num_vars = w.num_vars;
    if num_vars != s_id.num_vars || num_vars != s_perm.num_vars {
        return Err(PolyIOPErrors::InvalidParameters(
            "num of variables do not match".to_string(),
        ));
    }

    let prods = compute_products(beta, gamma, w, s_id, s_perm)?;
    // prods consists of the following:
    // - prod(x)
    // - prod(0, x)
    // - prod(1, x)
    // - prod(x, 0)
    // - prod(x, 1)
    // - numerator
    // - denominator
    let prod_0x = Rc::new(prods[1].clone());
    let numerator = Rc::new(prods[5].clone());
    let denominator = Rc::new(prods[6].clone());

    // compute (w(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
    // which is prods[6] * prod[1] * alpha
    let mut res = VirtualPolynomial::new_from_mle(denominator, *alpha);
    res.mul_by_mle(prod_0x, F::one())?;

    //   (w(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
    // - (w(x) + beta * s_id(x)   + gamma) * alpha
    res.add_mle_list([numerator], -*alpha)?;

    // Q(x) := prod(1,x) - prod(x, 0) - prod(x, 1)
    //       + alpha * (
    //             (w(x) + beta * s_perm(x) + gamma) * prod(0, x)
    //           - (w(x) + beta * s_id(x)   + gamma))
    let tmp = &(&prods[2] - &prods[3]) - &prods[4];
    let tmp = Rc::new(tmp);
    res.add_mle_list([tmp], F::one())?;

    end_timer!(start);
    Ok(res)
}

/// Compute the following 5 polynomials
/// - prod(x)
/// - prod(0, x)
/// - prod(1, x)
/// - prod(x, 0)
/// - prod(x, 1)
/// - numerator
/// - denominator
///
/// where
/// - `prod(0,x) := prod(0, x0, x1, …, xn)` which is the MLE over the
/// evaluations of the following polynomial on the boolean hypercube {0,1}^n:
///
///  (w(x) + \beta s_id(x) + \gamma)/(w(x) + \beta s_perm(x) + \gamma)
///
///   where
///   - beta and gamma are challenges
///   - w(x), s_id(x), s_perm(x) are mle-s
///
/// - `prod(1,x) := prod(x, 0) * prod(x, 1)`
/// - numerator is the MLE for `w(x) + \beta s_id(x) + \gamma`
/// - denominator is the MLE for `w(x) + \beta s_perm(x) + \gamma`
///
/// The caller needs to check num_vars matches in w/s_id/s_perm
/// Cost: linear in N.
#[allow(dead_code)]
// TODO: remove
pub(super) fn compute_products<F: PrimeField>(
    beta: &F,
    gamma: &F,
    w: &DenseMultilinearExtension<F>,
    s_id: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
) -> Result<[DenseMultilinearExtension<F>; 7], PolyIOPErrors> {
    let start = start_timer!(|| "compute all prod polynomial");

    let num_vars = w.num_vars;

    // ===================================
    // prod(0, x)
    // ===================================
    let (prod_0x, numerator, denominator) = compute_prod_0(beta, gamma, w, s_id, s_perm)?;

    // ===================================
    // prod(1, x)
    // ===================================
    //
    // `prod(1, x)` can be computed via recursing the following formula for 2^n-1
    // times
    //
    // `prod(1, x_1, ..., x_n) :=
    //      prod(x_1, x_2, ..., x_n, 0) * prod(x_1, x_2, ..., x_n, 1)`
    //
    // At any given step, the right hand side of the equation
    // is available via either eval_0x or the current view of eval_1x
    let eval_0x = &prod_0x.evaluations;
    let mut eval_1x = vec![];
    for x in 0..(1 << num_vars) - 1 {
        // sign will decide if the evaluation should be looked up from eval_0x or
        // eval_1x; x_zero_index is the index for the evaluation (x_2, ..., x_n,
        // 0); x_one_index is the index for the evaluation (x_2, ..., x_n, 1);
        let (x_zero_index, x_one_index, sign) = get_index(x, num_vars);
        if !sign {
            eval_1x.push(eval_0x[x_zero_index] * eval_0x[x_one_index]);
        } else {
            // sanity check: if we are trying to look up from the eval_1x table,
            // then the target index must already exist
            if x_zero_index >= eval_1x.len() || x_one_index >= eval_1x.len() {
                return Err(PolyIOPErrors::ShouldNotArrive);
            }
            eval_1x.push(eval_1x[x_zero_index] * eval_1x[x_one_index]);
        }
    }
    // prod(1, 1, ..., 1) := 0
    eval_1x.push(F::zero());

    // ===================================
    // prod(x)
    // ===================================
    // prod(x)'s evaluation is indeed `e := [eval_0x[..], eval_1x[..]].concat()`
    let eval = [eval_0x.as_slice(), eval_1x.as_slice()].concat();

    // ===================================
    // prod(x, 0) and prod(x, 1)
    // ===================================
    //
    // now we compute eval_x0 and eval_x1
    // eval_0x will be the odd coefficients of eval
    // and eval_1x will be the even coefficients of eval
    let mut eval_x0 = vec![];
    let mut eval_x1 = vec![];
    for (x, &prod_x) in eval.iter().enumerate() {
        if x & 1 == 0 {
            eval_x0.push(prod_x);
        } else {
            eval_x1.push(prod_x);
        }
    }

    let prod_1x = DenseMultilinearExtension::from_evaluations_vec(num_vars, eval_1x);
    let prod_x0 = DenseMultilinearExtension::from_evaluations_vec(num_vars, eval_x0);
    let prod_x1 = DenseMultilinearExtension::from_evaluations_vec(num_vars, eval_x1);
    let prod = DenseMultilinearExtension::from_evaluations_vec(num_vars + 1, eval);

    end_timer!(start);
    Ok([
        prod,
        prod_0x,
        prod_1x,
        prod_x0,
        prod_x1,
        numerator,
        denominator,
    ])
}

/// Returns three MLEs:
/// - prod(0,x)
/// - numerator
/// - denominator
///
/// where
/// - `prod(0,x) := prod(0, x1, …, xn)` which is the MLE over the
/// evaluations of the following polynomial on the boolean hypercube {0,1}^n:
///
///  (w(x) + \beta s_id(x) + \gamma)/(w(x) + \beta s_perm(x) + \gamma)
///
///  where
///  - beta and gamma are challenges
///  - w(x), s_id(x), s_perm(x) are mle-s
///
/// - numerator is the MLE for `w(x) + \beta s_id(x) + \gamma`
/// - denominator is the MLE for `w(x) + \beta s_perm(x) + \gamma`
///
/// The caller needs to check num_vars matches in w/s_id/s_perm
/// Cost: linear in N.
#[allow(dead_code)]
#[allow(clippy::type_complexity)]
// TODO: remove
pub(super) fn compute_prod_0<F: PrimeField>(
    beta: &F,
    gamma: &F,
    w: &DenseMultilinearExtension<F>,
    s_id: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
) -> Result<
    (
        DenseMultilinearExtension<F>,
        DenseMultilinearExtension<F>,
        DenseMultilinearExtension<F>,
    ),
    PolyIOPErrors,
> {
    let start = start_timer!(|| "compute prod(1,x)");

    let num_vars = w.num_vars;
    let mut prod_0x_evals = vec![];
    let mut numerator_evals = vec![];
    let mut denominator_evals = vec![];

    for (wi, (s_id_i, s_perm_i)) in w.iter().zip(s_id.iter().zip(s_perm.iter())) {
        let tmp = *wi + *gamma;
        let numerator = tmp + *beta * *s_id_i;
        let denominator = tmp + *beta * *s_perm_i;

        prod_0x_evals.push(numerator / denominator);
        numerator_evals.push(numerator);
        denominator_evals.push(denominator);
    }

    let prod_0x = DenseMultilinearExtension::from_evaluations_vec(num_vars, prod_0x_evals);
    let numerator = DenseMultilinearExtension::from_evaluations_vec(num_vars, numerator_evals);
    let denominator = DenseMultilinearExtension::from_evaluations_vec(num_vars, denominator_evals);

    end_timer!(start);
    Ok((prod_0x, numerator, denominator))
}

/// An MLE that represent an identity permutation: `f(index) \mapto index`
pub(super) fn identity_permutation_mle<F: PrimeField>(
    num_vars: usize,
) -> DenseMultilinearExtension<F> {
    let s_id_vec = (0..1u64 << num_vars).map(|index| F::from(index)).collect();
    DenseMultilinearExtension::from_evaluations_vec(num_vars, s_id_vec)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::bit_decompose;
    use ark_bls12_381::Fr;
    use ark_ff::UniformRand;
    use ark_poly::MultilinearExtension;
    use ark_std::test_rng;

    #[test]
    fn test_compute_prod_0() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        for num_vars in 2..6 {
            let w_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let w = DenseMultilinearExtension::from_evaluations_vec(num_vars, w_vec);

            let s_id = identity_permutation_mle::<Fr>(num_vars);

            let s_perm_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let s_perm = DenseMultilinearExtension::from_evaluations_vec(num_vars, s_perm_vec);

            let beta = Fr::rand(&mut rng);
            let gamma = Fr::rand(&mut rng);

            let prod_0 = compute_prod_0(&beta, &gamma, &w, &s_id, &s_perm)?;

            for i in 0..1 << num_vars {
                let r: Vec<Fr> = bit_decompose(i, num_vars)
                    .iter()
                    .map(|&x| Fr::from(x))
                    .collect();

                let eval = prod_0.0.evaluate(&r).unwrap();

                let w_eval = w.evaluate(&r).unwrap();
                let s_id_eval = s_id.evaluate(&r).unwrap();
                let s_perm_eval = s_perm.evaluate(&r).unwrap();
                let eval_rec =
                    (w_eval + beta * s_id_eval + gamma) / (w_eval + beta * s_perm_eval + gamma);

                assert_eq!(eval, eval_rec);
            }
        }
        Ok(())
    }

    #[test]
    fn test_compute_prod() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        for num_vars in 2..6 {
            let w_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let w = DenseMultilinearExtension::from_evaluations_vec(num_vars, w_vec);

            let s_id = identity_permutation_mle::<Fr>(num_vars);

            let s_perm_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let s_perm = DenseMultilinearExtension::from_evaluations_vec(num_vars, s_perm_vec);

            let beta = Fr::rand(&mut rng);
            let gamma = Fr::rand(&mut rng);

            // TODO: also test other 4 polynomials
            let res = compute_products(&beta, &gamma, &w, &s_id, &s_perm)?;

            for i in 0..1 << num_vars {
                let r: Vec<Fr> = bit_decompose(i, num_vars)
                    .iter()
                    .map(|&x| Fr::from(x))
                    .collect();

                let eval = res[1].evaluate(&r).unwrap();

                let w_eval = w.evaluate(&r).unwrap();
                let s_id_eval = s_id.evaluate(&r).unwrap();
                let s_perm_eval = s_perm.evaluate(&r).unwrap();
                let eval_rec =
                    (w_eval + beta * s_id_eval + gamma) / (w_eval + beta * s_perm_eval + gamma);

                assert_eq!(eval, eval_rec);
            }
        }
        Ok(())
    }

    #[test]
    fn test_compute_qx() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        for num_vars in 2..6 {
            let w_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let w = DenseMultilinearExtension::from_evaluations_vec(num_vars, w_vec);

            let s_id = identity_permutation_mle::<Fr>(num_vars);

            let s_perm_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let s_perm = DenseMultilinearExtension::from_evaluations_vec(num_vars, s_perm_vec);

            let alpha = Fr::rand(&mut rng);
            let beta = Fr::rand(&mut rng);
            let gamma = Fr::rand(&mut rng);

            // TODO: test the correctness of res
            let _res = build_q_x(&alpha, &beta, &gamma, &w, &s_id, &s_perm)?;
        }
        Ok(())
    }
}
