//! This module implements useful functions for the permutation check protocol.

use crate::{utils::get_index, PolyIOPErrors, VirtualPolynomial};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, rand::RngCore, start_timer};
use std::rc::Rc;

/// Compute `Q(x)` which a virtual polynomial of degree 2 defined as
///
/// Q(x) := prod(1,x)
///       - prod(x, 0) * prod(x, 1)
///       + alpha * (
///             (g(x) + beta * s_perm(x) + gamma) * prod(0, x)
///           - (f(x) + beta * s_id(x)   + gamma))
///
/// The caller needs to check num_vars matches in f/g/s_id/s_perm
///
/// Cost: linear in N.
pub(super) fn build_q_x<F: PrimeField>(
    alpha: &F,
    beta: &F,
    gamma: &F,
    fx: &DenseMultilinearExtension<F>,
    gx: &DenseMultilinearExtension<F>,
    s_id: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
) -> Result<VirtualPolynomial<F>, PolyIOPErrors> {
    let start = start_timer!(|| "compute Q(x)");

    let prods = compute_products(beta, gamma, fx, gx, s_id, s_perm)?;
    // prods consists of the following:
    // - prod(x)
    // - prod(0, x)
    // - prod(1, x)
    // - prod(x, 0)
    // - prod(x, 1)
    // - numerator
    // - denominator
    let prod_0x = Rc::new(prods[1].clone());
    let prod_1x = Rc::new(prods[2].clone());
    let prod_x1 = Rc::new(prods[3].clone());
    let prod_x0 = Rc::new(prods[4].clone());
    let numerator = Rc::new(prods[5].clone());
    let denominator = Rc::new(prods[6].clone());

    // compute (g(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
    // which is prods[6] * prod[1] * alpha
    let mut res = VirtualPolynomial::new_from_mle(denominator, F::one());
    res.mul_by_mle(prod_0x, *alpha)?;

    //   (g(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
    // - (f(x) + beta * s_id(x)   + gamma) * alpha
    res.add_mle_list([numerator], -*alpha)?;

    // Q(x) := prod(1,x) - prod(x, 0) * prod(x, 1)
    //       + alpha * (
    //             (g(x) + beta * s_perm(x) + gamma) * prod(0, x)
    //           - (f(x) + beta * s_id(x)   + gamma))
    res.add_mle_list([prod_x0, prod_x1], -F::one())?;
    res.add_mle_list([prod_1x], F::one())?;

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
///  (f(x) + \beta s_id(x) + \gamma)/(g(x) + \beta s_perm(x) + \gamma)
///
///   where
///   - beta and gamma are challenges
///   - f(x), g(x), s_id(x), s_perm(x) are mle-s
///
/// - `prod(1,x) := prod(x, 0) * prod(x, 1)`
/// - numerator is the MLE for `f(x) + \beta s_id(x) + \gamma`
/// - denominator is the MLE for `g(x) + \beta s_perm(x) + \gamma`
///
/// The caller needs to check num_vars matches in f/g/s_id/s_perm
/// Cost: linear in N.
pub(super) fn compute_products<F: PrimeField>(
    beta: &F,
    gamma: &F,
    fx: &DenseMultilinearExtension<F>,
    gx: &DenseMultilinearExtension<F>,
    s_id: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
) -> Result<[DenseMultilinearExtension<F>; 7], PolyIOPErrors> {
    let start = start_timer!(|| "compute all prod polynomial");

    let num_vars = fx.num_vars;

    // ===================================
    // prod(0, x)
    // ===================================
    let (prod_0x, numerator, denominator) = compute_prod_0(beta, gamma, fx, gx, s_id, s_perm)?;

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
///  (f(x) + \beta s_id(x) + \gamma)/(g(x) + \beta s_perm(x) + \gamma)
///
///  where
///  - beta and gamma are challenges
///  - f(x), g(x), s_id(x), s_perm(x) are mle-s
///
/// - numerator is the MLE for `f(x) + \beta s_id(x) + \gamma`
/// - denominator is the MLE for `g(x) + \beta s_perm(x) + \gamma`
///
/// The caller needs to check num_vars matches in f/g/s_id/s_perm
/// Cost: linear in N.
#[allow(clippy::type_complexity)]
pub(super) fn compute_prod_0<F: PrimeField>(
    beta: &F,
    gamma: &F,
    fx: &DenseMultilinearExtension<F>,
    gx: &DenseMultilinearExtension<F>,
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

    let num_vars = fx.num_vars;
    let mut prod_0x_evals = vec![];
    let mut numerator_evals = vec![];
    let mut denominator_evals = vec![];

    for (&fi, (&gi, (&s_id_i, &s_perm_i))) in
        fx.iter().zip(gx.iter().zip(s_id.iter().zip(s_perm.iter())))
    {
        let numerator = fi + *beta * s_id_i + gamma;
        let denominator = gi + *beta * s_perm_i + gamma;

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
pub fn identity_permutation_mle<F: PrimeField>(num_vars: usize) -> DenseMultilinearExtension<F> {
    let s_id_vec = (0..1u64 << num_vars).map(F::from).collect();
    DenseMultilinearExtension::from_evaluations_vec(num_vars, s_id_vec)
}

/// An MLE that represent a random permutation
pub fn random_permutation_mle<F: PrimeField, R: RngCore>(
    num_vars: usize,
    rng: &mut R,
) -> DenseMultilinearExtension<F> {
    let len = 1u64 << num_vars;
    let mut s_id_vec: Vec<F> = (0..len).map(F::from).collect();
    let mut s_perm_vec = vec![];
    for _ in 0..len {
        let index = rng.next_u64() as usize % s_id_vec.len();
        s_perm_vec.push(s_id_vec.remove(index));
    }
    DenseMultilinearExtension::from_evaluations_vec(num_vars, s_perm_vec)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::bit_decompose;
    use ark_bls12_381::Fr;
    use ark_ff::{UniformRand, Zero};
    use ark_poly::MultilinearExtension;
    use ark_std::test_rng;

    #[test]
    fn test_compute_prod_0() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        for num_vars in 2..6 {
            let f = DenseMultilinearExtension::rand(num_vars, &mut rng);
            let g = DenseMultilinearExtension::rand(num_vars, &mut rng);

            let s_id = identity_permutation_mle::<Fr>(num_vars);
            let s_perm = random_permutation_mle(num_vars, &mut rng);

            let beta = Fr::rand(&mut rng);
            let gamma = Fr::rand(&mut rng);

            let (prod_0, numerator, denominator) =
                compute_prod_0(&beta, &gamma, &f, &g, &s_id, &s_perm)?;

            for i in 0..1 << num_vars {
                let r: Vec<Fr> = bit_decompose(i, num_vars)
                    .iter()
                    .map(|&x| Fr::from(x))
                    .collect();

                let prod_0_eval = prod_0.evaluate(&r).unwrap();
                let numerator_eval = numerator.evaluate(&r).unwrap();
                let denominator_eval = denominator.evaluate(&r).unwrap();

                let f_eval = f.evaluate(&r).unwrap();
                let g_eval = g.evaluate(&r).unwrap();
                let s_id_eval = s_id.evaluate(&r).unwrap();
                let s_perm_eval = s_perm.evaluate(&r).unwrap();

                let numerator_eval_rec = f_eval + beta * s_id_eval + gamma;
                let denominator_eval_rec = g_eval + beta * s_perm_eval + gamma;
                let prod_0_eval_rec = numerator_eval_rec / denominator_eval_rec;

                assert_eq!(numerator_eval, numerator_eval_rec);
                assert_eq!(denominator_eval, denominator_eval_rec);
                assert_eq!(prod_0_eval, prod_0_eval_rec);
            }
        }
        Ok(())
    }

    #[test]
    fn test_compute_prod() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        for num_vars in 2..6 {
            let f = DenseMultilinearExtension::rand(num_vars, &mut rng);
            let g = DenseMultilinearExtension::rand(num_vars, &mut rng);

            let s_id = identity_permutation_mle::<Fr>(num_vars);
            let s_perm = random_permutation_mle(num_vars, &mut rng);

            let beta = Fr::rand(&mut rng);
            let gamma = Fr::rand(&mut rng);

            let res = compute_products(&beta, &gamma, &f, &g, &s_id, &s_perm)?;

            for i in 0..1 << num_vars {
                let r: Vec<Fr> = bit_decompose(i, num_vars)
                    .iter()
                    .map(|&x| Fr::from(x))
                    .collect();

                let eval = res[1].evaluate(&r).unwrap();

                let f_eval = f.evaluate(&r).unwrap();
                let g_eval = g.evaluate(&r).unwrap();
                let s_id_eval = s_id.evaluate(&r).unwrap();
                let s_perm_eval = s_perm.evaluate(&r).unwrap();
                let eval_rec =
                    (f_eval + beta * s_id_eval + gamma) / (g_eval + beta * s_perm_eval + gamma);

                assert_eq!(eval, eval_rec);
            }
        }
        Ok(())
    }

    #[test]
    fn test_compute_qx() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        for num_vars in 2..6 {
            let f = DenseMultilinearExtension::rand(num_vars, &mut rng);
            let g = DenseMultilinearExtension::rand(num_vars, &mut rng);

            let s_id = identity_permutation_mle::<Fr>(num_vars);
            let s_perm = random_permutation_mle(num_vars, &mut rng);

            let alpha = Fr::rand(&mut rng);
            let beta = Fr::rand(&mut rng);
            let gamma = Fr::rand(&mut rng);

            let qx = build_q_x(&alpha, &beta, &gamma, &f, &g, &s_id, &s_perm)?;

            // test q_x is a 0 over boolean hypercube
            for i in 0..1 << num_vars {
                let bit_sequence = bit_decompose(i, num_vars);
                let eval: Vec<Fr> = bit_sequence.iter().map(|x| Fr::from(*x as u64)).collect();
                let res = qx.evaluate(&eval)?;
                assert!(res.is_zero())
            }
        }
        Ok(())
    }
}
