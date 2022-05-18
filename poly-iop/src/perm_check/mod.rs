//! This module implements the permutation check protocol.

use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;

use crate::PolyIOPErrors;

/// Compute `prod(0,x)` which is
///
///  (w(x) + \beta s_id(x) + \gamma)/(w(x) + \beta s_perm(x) + \gamma)
///
/// where
/// - beta and gamma are challenges
/// - w(x), s_id(x), s_perm(x) are mle-s
#[allow(dead_code)]
// TODO: remove
fn compute_prod_0<F: PrimeField>(
    beta: &F,
    gamma: &F,
    w: &DenseMultilinearExtension<F>,
    s_id: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
) -> Result<DenseMultilinearExtension<F>, PolyIOPErrors> {
    let num_vars = w.num_vars;

    if num_vars != s_id.num_vars || num_vars != s_perm.num_vars {
        return Err(PolyIOPErrors::InvalidParameters(
            "num of variables do not match".to_string(),
        ));
    }
    let eval: Vec<F> = w
        .iter()
        .zip(s_id.iter().zip(s_perm.iter()))
        .map(|(wi, (s_id_i, s_perm_i))| {
            let tmp = *wi + *gamma;
            (tmp + *beta * *s_id_i) / (tmp + *beta * *s_perm_i)
        })
        .collect();

    let res = DenseMultilinearExtension::from_evaluations_vec(num_vars, eval);

    Ok(res)
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

            let s_id_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let s_id = DenseMultilinearExtension::from_evaluations_vec(num_vars, s_id_vec);

            let s_perm_vec: Vec<Fr> = (0..(1 << num_vars)).map(|_| Fr::rand(&mut rng)).collect();
            let s_perm = DenseMultilinearExtension::from_evaluations_vec(num_vars, s_perm_vec);

            let beta = Fr::rand(&mut rng);
            let gamma = Fr::rand(&mut rng);
            for i in 0..1 << num_vars {
                let r: Vec<Fr> = bit_decompose(i, num_vars)
                    .iter()
                    .map(|&x| Fr::from(x))
                    .collect();

                let prod_0 = compute_prod_0(&beta, &gamma, &w, &s_id, &s_perm)?;

                let eval = prod_0.evaluate(&r).unwrap();

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
}
