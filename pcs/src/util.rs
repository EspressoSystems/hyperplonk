use ark_ec::{msm::VariableBaseMSM, AffineCurve};
use ark_ff::PrimeField;
use rayon::prelude::*;

pub fn multi_scalar_mul<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    #[cfg(feature = "parallel")]
    return parallel_msm(bases, scalars);
    #[cfg(not(feature = "parallel"))]
    VariableBaseMSM::multi_scalar_mul(bases, scalars)
}

fn parallel_msm<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    let num_thread = rayon::current_num_threads();
    let block_size = (bases.len() + num_thread - 1) / num_thread;

    let base_blocks = bases.par_chunks(block_size);
    let scalar_blocks = scalars.par_chunks(block_size);
    let results = base_blocks
        .into_par_iter()
        .zip(scalar_blocks.into_par_iter())
        .map(|(x, g)| VariableBaseMSM::multi_scalar_mul(x, g))
        .collect::<Vec<_>>();

    let mut res = results[0];
    for e in results.iter().skip(1) {
        res += e
    }

    res
}

#[cfg(test)]
mod test {
    use ark_bls12_381::{Fr, G1Projective};
    use ark_ec::{msm::VariableBaseMSM, ProjectiveCurve};
    use ark_ff::{PrimeField, UniformRand};
    use ark_std::test_rng;

    use super::multi_scalar_mul;

    #[test]
    fn test_para_msm() {
        let mut rng = test_rng();
        for i in 10usize..15 {
            let size = 1 << i;
            let bases: Vec<_> = (0..size)
                .map(|_| G1Projective::rand(&mut rng).into_affine())
                .collect();
            let scalars: Vec<_> = (0..size).map(|_| Fr::rand(&mut rng).into_repr()).collect();
            let res1 = multi_scalar_mul(&bases, &scalars).into_affine();
            let res2 = VariableBaseMSM::multi_scalar_mul(&bases, &scalars).into_affine();
            assert_eq!(res1, res2)
        }
    }
}
