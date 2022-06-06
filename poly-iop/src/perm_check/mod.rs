//! This module implements the permutation check protocol.

use crate::{utils::bit_decompose, PolyIOPErrors};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};

/// Compute `prod(0,x) := prod(0, x1, …, xn)` which is the MLE over the
/// evaluations of the following polynomial on the boolean hypercube {0,1}^n:
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

/// TODO: description
/// Input prod(0, x)
/// Return [prod(1,x):= prod(x,0) * prod(x,1), prod(x,0), prod(x,1)]
#[allow(dead_code)]
// TODO: remove
fn compute_prod_1<F: PrimeField>(
    prod_0: &DenseMultilinearExtension<F>,
) -> Result<[DenseMultilinearExtension<F>; 3], PolyIOPErrors> {
    let num_vars = prod_0.num_vars - 1;

    let mut prod_x0_evals = vec![];
    let mut prod_x1_evals = vec![];
    let mut prod_1x_evals = vec![];

    for i in 0..1u64 << num_vars {
        let sequence = bit_decompose(i, num_vars);
        let sequence_f: Vec<F> = sequence.iter().map(|&x| F::from(x)).collect();

        // prod(x, 0)
        let eval = prod_0.evaluate(&[sequence_f.as_slice(), &[F::zero()]].concat());
        let eval_x_0 = match eval {
            Some(p) => p,
            None => {
                return Err(PolyIOPErrors::InvalidParameters(
                    "Evaluation failed".to_string(),
                ))
            },
        };
        prod_x0_evals.push(eval_x_0);

        // prod(x, 1)
        let eval = prod_0.evaluate(&[sequence_f.as_slice(), &[F::one()]].concat());
        let eval_x_1 = match eval {
            Some(p) => p,
            None => {
                return Err(PolyIOPErrors::InvalidParameters(
                    "Evaluation failed".to_string(),
                ))
            },
        };
        prod_x1_evals.push(eval_x_1);

        // prod(1, x)
        prod_1x_evals.push(eval_x_0 * eval_x_1);
    }

    let prod_x_0 = DenseMultilinearExtension::from_evaluations_vec(num_vars, prod_x0_evals);
    let prod_x_1 = DenseMultilinearExtension::from_evaluations_vec(num_vars, prod_x1_evals);
    let prod_1_x = DenseMultilinearExtension::from_evaluations_vec(num_vars, prod_1x_evals);

    Ok([prod_1_x, prod_x_0, prod_x_1])
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

    #[test]
    fn test_compute_prod_1() -> Result<(), PolyIOPErrors> {
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
                // println!("{:?}", prod_0);
                let res = compute_prod_1(&prod_0)?;
                let _prod_x0 = &res[0];
                let _prod_x1 = &res[1];
                let _prod_1x = &res[2];

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
