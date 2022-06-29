use ark_ec::{
    msm::{FixedBaseMSM, VariableBaseMSM},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::PrimeField;
use ark_poly::MultilinearExtension;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, rand::RngCore, start_timer, vec::Vec, One, Zero};

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
    pub proofs: Vec<E::G1Affine>,
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
        let setup_timer = start_timer!(|| format!("SRS setup for dim {}", num_vars));
        let res = Self::SRS::gen_srs_for_testing(rng, num_vars);
        end_timer!(setup_timer);
        res
    }

    /// Generate a commitment for a polynomial.
    ///
    /// This function takes `2^num_vars` number of scalar multiplications over
    /// G1.
    fn commit(
        prover_param: &Self::ProverParam,
        poly: &impl MultilinearExtension<E::Fr>,
    ) -> Result<Self::Commitment, PCSErrors> {
        let commit_timer = start_timer!(|| "commit");

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

        end_timer!(commit_timer);
        Ok(Commitment { nv, g_product })
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same. This function does not need to take the evaluation value as an
    /// input.
    ///
    /// This function takes 2^{num_var +1} number of scalar multiplications over
    /// G2:
    /// - it proceeds with `num_var` number of rounds,
    /// - at round i, we compute an MSM for `2^{num_var - i + 1}` number of G2
    ///   elements.
    fn open(
        prover_param: &Self::ProverParam,
        polynomial: &impl MultilinearExtension<E::Fr>,
        point: &[E::Fr],
    ) -> Result<Self::Proof, PCSErrors> {
        let open_timer = start_timer!(|| "open");

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

        for (i, (&point_at_k, gi)) in point
            .iter()
            .zip(prover_param.powers_of_g.iter())
            .take(nv)
            .enumerate()
        {
            let ith_round = start_timer!(|| format!("{}-th round", i));

            let k = nv - i;
            let cur_dim = 1 << (k - 1);
            let mut cur_q = vec![E::Fr::zero(); cur_dim];
            let mut cur_r = vec![E::Fr::zero(); cur_dim];

            for b in 0..(1 << (k - 1)) {
                // q_b = pre_r [2^b + 1] - pre_r [2^b]
                cur_q[b] = r[k][(b << 1) + 1] - r[k][b << 1];

                // r_b = pre_r [2^b]*(1-p) + pre_r [2^b + 1] * p
                cur_r[b] =
                    r[k][b << 1] * (E::Fr::one() - point_at_k) + (r[k][(b << 1) + 1] * point_at_k);
            }

            let scalars: Vec<_> = (0..(1 << k)).map(|x| cur_q[x >> 1].into_repr()).collect();

            q[k] = cur_q;
            r[k - 1] = cur_r;

            // this is a MSM over G1 and is likely to be the bottleneck
            proofs.push(VariableBaseMSM::multi_scalar_mul(&gi.evals, &scalars).into_affine());
            end_timer!(ith_round);
        }

        end_timer!(open_timer);
        Ok(Proof { proofs })
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    ///
    /// This function takes
    /// - num_var number of pairing product.
    /// - num_var number of MSM
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &[E::Fr],
        value: E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSErrors> {
        let verify_timer = start_timer!(|| "verify");
        let prepare_inputs_timer = start_timer!(|| "prepare pairing inputs");

        let scalar_size = E::Fr::size_in_bits();
        let window_size = FixedBaseMSM::get_mul_window_size(verifier_param.num_vars);

        let h_table = FixedBaseMSM::get_window_table(
            scalar_size,
            window_size,
            verifier_param.h.into_projective(),
        );
        let h_mul: Vec<E::G2Projective> =
            FixedBaseMSM::multi_scalar_mul(scalar_size, window_size, &h_table, point);

        let h_vec: Vec<_> = (0..verifier_param.num_vars)
            .map(|i| verifier_param.h_mask[i].into_projective() - h_mul[i])
            .collect();
        let h_vec: Vec<E::G2Affine> = E::G2Projective::batch_normalization_into_affine(&h_vec);
        end_timer!(prepare_inputs_timer);

        let pairing_product_timer = start_timer!(|| "pairing product");

        let mut pairings: Vec<_> = proof
            .proofs
            .iter()
            .map(|&x| E::G1Prepared::from(x))
            .zip(
                h_vec
                    .into_iter()
                    .take(verifier_param.num_vars)
                    .map(E::G2Prepared::from),
            )
            .collect();

        pairings.push((
            E::G1Prepared::from(
                (verifier_param.g.mul(value) - commitment.g_product.into_projective())
                    .into_affine(),
            ),
            E::G2Prepared::from(verifier_param.h),
        ));

        let res = E::product_of_pairings(pairings.iter()) == E::Fqk::one();

        end_timer!(pairing_product_timer);
        end_timer!(verify_timer);
        Ok(res)
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

        let value = Fr::rand(rng);
        assert!(!KZGMultilinearPC::verify(&vk, &com, &point, value, &proof)?);

        Ok(())
    }

    #[test]
    fn setup_commit_verify_correct_polynomials() -> Result<(), PCSErrors> {
        let mut rng = test_rng();

        let uni_params = KZGMultilinearPC::<E>::setup(&mut rng, 10)?;

        // normal polynomials
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
}
