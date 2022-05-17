use ark_ec::{
    msm::{FixedBaseMSM, VariableBaseMSM},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::PrimeField;
use ark_poly::MultilinearExtension;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{rand::RngCore, vec::Vec, One, Zero};

use crate::{
    KZGMultilinearPC, MultilinearCommitmentScheme, PCSErrors, ProverParam, UniversalParams,
    VerifierParam,
};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// commitment
pub struct Commitment<E: PairingEngine> {
    /// number of variables
    pub nv: usize,
    /// product of g as described by the vRAM paper
    pub g_product: E::G1Affine,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// proof of opening
pub struct Proof<E: PairingEngine> {
    /// Evaluation of quotients
    pub proofs: Vec<E::G2Affine>,
}

impl<E: PairingEngine> MultilinearCommitmentScheme<E> for KZGMultilinearPC<E> {
    type ProverParam = ProverParam<E>;
    type VerifierParam = VerifierParam<E>;
    type SRS = UniversalParams<E>;
    type Commitment = Commitment<E>;
    type Proof = Proof<E>;

    /// Generate SRS from RNG.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn setup<R: RngCore>(rng: &mut R, num_vars: usize) -> Result<Self::SRS, PCSErrors> {
        Self::SRS::gen_srs_for_testing(rng, num_vars)
    }

    /// Generate a commitment for a polynomial
    fn commit(
        prover_param: &Self::ProverParam,
        poly: &impl MultilinearExtension<E::Fr>,
    ) -> Result<Self::Commitment, PCSErrors> {
        let nv = poly.num_vars();
        let scalars: Vec<_> = poly
            .to_evaluations()
            .into_iter()
            .map(|x| x.into_repr())
            .collect();
        let g_product = VariableBaseMSM::multi_scalar_mul(
            &prover_param.powers_of_g[0].evals,
            scalars.as_slice(),
        )
        .into_affine();
        Ok(Commitment { nv, g_product })
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: &Self::ProverParam,
        polynomial: &impl MultilinearExtension<E::Fr>,
        point: &[E::Fr],
    ) -> Result<Self::Proof, PCSErrors> {
        assert_eq!(
            polynomial.num_vars(),
            prover_param.num_vars,
            "Invalid size of polynomial"
        );
        let nv = polynomial.num_vars();
        let mut r: Vec<Vec<E::Fr>> = (0..nv + 1).map(|_| Vec::new()).collect();
        let mut q: Vec<Vec<E::Fr>> = (0..nv + 1).map(|_| Vec::new()).collect();

        r[nv] = polynomial.to_evaluations();

        let mut proofs = Vec::new();

        // todo: refactor the following code for better readability
        for (i, (&p, hi)) in point
            .iter()
            .zip(prover_param.powers_of_h.iter())
            .take(nv)
            .enumerate()
        {
            let k = nv - i;
            let point_at_k = p;
            q[k] = (0..(1 << (k - 1))).map(|_| E::Fr::zero()).collect();
            r[k - 1] = (0..(1 << (k - 1))).map(|_| E::Fr::zero()).collect();
            for b in 0..(1 << (k - 1)) {
                q[k][b] = r[k][(b << 1) + 1] - r[k][b << 1];
                r[k - 1][b] =
                    r[k][b << 1] * (E::Fr::one() - point_at_k) + (r[k][(b << 1) + 1] * point_at_k);
            }
            let scalars: Vec<_> = (0..(1 << k))
                .map(|x| q[k][x >> 1].into_repr()) // fine
                .collect();

            let pi_h = VariableBaseMSM::multi_scalar_mul(&hi.evals, &scalars).into_affine(); // no need to move outside and partition
            proofs.push(pi_h);
        }

        Ok(Proof { proofs })
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &[E::Fr],
        value: E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSErrors> {
        let scalar_size = E::Fr::size_in_bits();
        let window_size = FixedBaseMSM::get_mul_window_size(verifier_param.num_vars);

        let g_table = FixedBaseMSM::get_window_table(
            scalar_size,
            window_size,
            verifier_param.g.into_projective(),
        );
        let g_mul: Vec<E::G1Projective> =
            FixedBaseMSM::multi_scalar_mul(scalar_size, window_size, &g_table, point);

        let mut g1_vec: Vec<_> = (0..verifier_param.num_vars)
            .map(|i| verifier_param.g_mask[i].into_projective() - g_mul[i])
            .collect();
        g1_vec.push(verifier_param.g.mul(value) - commitment.g_product.into_projective());

        let g1_vec: Vec<E::G1Affine> = E::G1Projective::batch_normalization_into_affine(&g1_vec);
        let tmp = g1_vec[verifier_param.num_vars];

        let mut pairings: Vec<_> = g1_vec
            .into_iter()
            .take(verifier_param.num_vars)
            .map(E::G1Prepared::from)
            .zip(proof.proofs.iter().map(|&x| E::G2Prepared::from(x)))
            .collect();

        pairings.push((
            E::G1Prepared::from(tmp),
            E::G2Prepared::from(verifier_param.h),
        ));

        Ok(E::product_of_pairings(pairings.iter()) == E::Fqk::one())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
    use ark_std::{rand::RngCore, test_rng, vec::Vec, UniformRand};
    type E = Bls12_381;
    type Fr = <E as PairingEngine>::Fr;

    fn test_kzg_mlpc_helper<R: RngCore>(
        uni_params: &UniversalParams<E>,
        poly: &impl MultilinearExtension<Fr>,
        rng: &mut R,
    ) -> Result<(), PCSErrors> {
        let nv = poly.num_vars();
        assert_ne!(nv, 0);
        let (ck, vk) = uni_params.trim(nv)?;
        let point: Vec<_> = (0..nv).map(|_| Fr::rand(rng)).collect();
        let com = KZGMultilinearPC::commit(&ck, poly)?;
        let proof = KZGMultilinearPC::open(&ck, poly, &point)?;

        let value = poly.evaluate(&point).unwrap();
        assert!(KZGMultilinearPC::verify(&vk, &com, &point, value, &proof)?);

        Ok(())
    }

    #[test]
    fn setup_commit_verify_correct_polynomials() -> Result<(), PCSErrors> {
        let mut rng = test_rng();

        // normal polynomials
        let uni_params = KZGMultilinearPC::<E>::setup(&mut rng, 10)?;

        let poly1 = DenseMultilinearExtension::rand(8, &mut rng);
        test_kzg_mlpc_helper(&uni_params, &poly1, &mut rng)?;

        let poly2 = SparseMultilinearExtension::rand_with_config(9, 1 << 5, &mut rng);
        test_kzg_mlpc_helper(&uni_params, &poly2, &mut rng)?;

        // single-variate polynomials

        let poly3 = DenseMultilinearExtension::rand(1, &mut rng);
        test_kzg_mlpc_helper(&uni_params, &poly3, &mut rng)?;

        let poly4 = SparseMultilinearExtension::rand_with_config(1, 1 << 1, &mut rng);
        test_kzg_mlpc_helper(&uni_params, &poly4, &mut rng)?;
        Ok(())
    }

    #[test]
    fn setup_commit_verify_constant_polynomial() {
        let mut rng = test_rng();

        // normal polynomials
        assert!(KZGMultilinearPC::<E>::setup(&mut rng, 0).is_err());
    }

    #[test]
    fn setup_commit_verify_incorrect_polynomial_should_return_false() -> Result<(), PCSErrors> {
        let mut rng = test_rng();
        let nv = 8;
        let uni_params = KZGMultilinearPC::<E>::setup(&mut rng, nv)?;
        let poly = DenseMultilinearExtension::rand(nv, &mut rng);
        let nv = uni_params.prover_param.num_vars;
        let (ck, vk) = uni_params.trim(nv)?;
        let point: Vec<_> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();
        let com = KZGMultilinearPC::commit(&ck, &poly)?;
        let proof = KZGMultilinearPC::open(&ck, &poly, &point)?;

        let value = poly.evaluate(&point).unwrap();
        assert!(!KZGMultilinearPC::verify(
            &vk,
            &com,
            &point,
            value + &(1u16.into()),
            &proof
        )?);
        Ok(())
    }
}
