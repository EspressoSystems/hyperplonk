//! This module implements useful functions for the permutation check protocol.

use crate::PolyIOPErrors;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, rand::RngCore, start_timer};

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

    let s_id = identity_permutation_mle::<F>(num_vars);

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
    use ark_ff::UniformRand;
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

            let (prod_0, numerator, denominator) = compute_prod_0(&beta, &gamma, &f, &g, &s_perm)?;

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
}
