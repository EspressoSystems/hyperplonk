//! Main module for univariate KZG commitment scheme

use crate::StructuredReferenceString;
use crate::{
    prelude::{Commitment, PCSErrors},
    PolynomialCommitmentScheme,
};
use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, rand::RngCore, start_timer, One, UniformRand, Zero};
use srs::{UnivariateProverParam, UnivariateUniversalParams, UnivariateVerifierParam};
use std::marker::PhantomData;

pub(crate) mod srs;

/// KZG Polynomial Commitment Scheme on univariate polynomial.
pub struct KZGUnivariatePCS<E: PairingEngine> {
    #[doc(hidden)]
    phantom: PhantomData<E>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// proof of opening
pub struct KZGUnivariateOpening<E: PairingEngine> {
    /// Evaluation of quotients
    pub proof: E::G1Affine,
}

impl<E: PairingEngine> PolynomialCommitmentScheme<E> for KZGUnivariatePCS<E> {
    type ProverParam = UnivariateProverParam<E::G1Affine>;
    type VerifierParam = UnivariateVerifierParam<E>;
    type SRS = UnivariateUniversalParams<E>;
    type Polynomial = DensePolynomial<E::Fr>;
    type Point = E::Fr;
    type Commitment = Commitment<E>;
    type Proof = KZGUnivariateOpening<E>;
    type BatchProof = Vec<Self::Proof>;
    type BatchCommitment = Vec<Self::Commitment>;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `log_size` is the log of maximum degree.
    /// - For multilinear polynomials, `log_size` is the number of variables.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore>(
        rng: &mut R,
        log_size: usize,
    ) -> Result<Self::SRS, PCSErrors> {
        Self::SRS::gen_srs_for_testing(rng, log_size)
    }

    /// Generate a commitment for a polynomial
    /// Note that the scheme is not hidding
    fn commit(
        prover_param: &Self::ProverParam,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, PCSErrors> {
        let commit_time =
            start_timer!(|| format!("Committing to polynomial of degree {} ", poly.degree()));

        if poly.degree() > prover_param.powers_of_g.len() {
            return Err(PCSErrors::InvalidParameters(format!(
                "poly degree {} is larger than allowed {}",
                poly.degree(),
                prover_param.powers_of_g.len()
            )));
        }

        let (num_leading_zeros, plain_coeffs) = skip_leading_zeros_and_convert_to_bigints(poly);

        let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let commitment = VariableBaseMSM::multi_scalar_mul(
            &prover_param.powers_of_g[num_leading_zeros..],
            &plain_coeffs,
        )
        .into_affine();
        end_timer!(msm_time);

        end_timer!(commit_time);
        Ok(Commitment { commitment })
    }

    /// Generate a commitment for a list of polynomials
    fn multi_commit(
        prover_param: &Self::ProverParam,
        polys: &[Self::Polynomial],
    ) -> Result<Self::BatchCommitment, PCSErrors> {
        let commit_time =
            start_timer!(|| format!("batch commit {} polynomials", polynomials.degree()));
        let res = polys
            .iter()
            .map(|poly| Self::commit(prover_param, poly))
            .collect::<Result<Vec<Self::Commitment>, PCSErrors>>()?;

        end_timer!(commit_time);
        Ok(res)
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: &Self::ProverParam,
        polynomial: &Self::Polynomial,
        point: &Self::Point,
    ) -> Result<Self::Proof, PCSErrors> {
        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));
        let divisor = Self::Polynomial::from_coefficients_vec(vec![-*point, E::Fr::one()]);

        let witness_time = start_timer!(|| "Computing witness polynomial");
        let witness_polynomial = polynomial / &divisor;
        end_timer!(witness_time);

        let (num_leading_zeros, witness_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(&witness_polynomial);

        let proof = VariableBaseMSM::multi_scalar_mul(
            &prover_param.powers_of_g[num_leading_zeros..],
            &witness_coeffs,
        )
        .into_affine();

        end_timer!(open_time);
        Ok(Self::Proof { proof })
    }

    /// Input a list of polynomials, and a same number of points,
    /// compute a multi-opening for all the polynomials.
    fn multi_open(
        prover_param: &Self::ProverParam,
        _multi_commitment: &Self::Commitment,
        polynomials: &[Self::Polynomial],
        points: &[Self::Point],
    ) -> Result<Self::BatchProof, PCSErrors> {
        let open_time =
            start_timer!(|| format!("batch opening {} polynomials", polynomials.degree()));
        if polynomials.len() != points.len() {
            return Err(PCSErrors::InvalidParameters(format!(
                "poly length {} is different from points length {}",
                polynomials.len(),
                points.len()
            )));
        }
        let batch_proof = polynomials
            .iter()
            .zip(points.iter())
            .map(|(poly, point)| Self::open(prover_param, poly, point))
            .collect::<Result<Self::BatchProof, PCSErrors>>()?;
        end_timer!(open_time);
        Ok(batch_proof)
    }
    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &Self::Point,
        value: &E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSErrors> {
        let check_time = start_timer!(|| "Checking evaluation");
        let pairing_inputs: Vec<(E::G1Prepared, E::G2Prepared)> = vec![
            (
                (commitment.commitment.into_projective() - verifier_param.g.mul(value.into_repr()))
                    .into_affine()
                    .into(),
                verifier_param.h.into(),
            ),
            (
                proof.proof.into(),
                (verifier_param.h.mul(point.into_repr()) - verifier_param.beta_h.into_projective())
                    .into_affine()
                    .into(),
            ),
        ];

        let res = E::product_of_pairings(pairing_inputs.iter());

        end_timer!(check_time, || format!("Result: {}", res));
        Ok(res.is_one())
    }

    /// Verifies that `value_i` is the evaluation at `x_i` of the polynomial
    /// `poly_i` committed inside `comm`.
    fn batch_verify<R: RngCore>(
        verifier_param: &Self::VerifierParam,
        multi_commitment: &Self::BatchCommitment,
        points: &[Self::Point],
        values: &[E::Fr],
        batch_proof: &Self::BatchProof,
        rng: &mut R,
    ) -> Result<bool, PCSErrors> {
        let check_time =
            start_timer!(|| format!("Checking {} evaluation proofs", commitments.len()));

        let mut total_c = <E::G1Projective>::zero();
        let mut total_w = <E::G1Projective>::zero();

        let combination_time = start_timer!(|| "Combining commitments and proofs");
        let mut randomizer = E::Fr::one();
        // Instead of multiplying g and gamma_g in each turn, we simply accumulate
        // their coefficients and perform a final multiplication at the end.
        let mut g_multiplier = E::Fr::zero();
        for (((c, z), v), proof) in multi_commitment
            .iter()
            .zip(points)
            .zip(values)
            .zip(batch_proof)
        {
            let w = proof.proof;
            let mut temp = w.mul(*z);
            temp.add_assign_mixed(&c.commitment);
            let c = temp;
            g_multiplier += &(randomizer * v);
            total_c += &c.mul(randomizer.into_repr());
            total_w += &w.mul(randomizer.into_repr());
            // We don't need to sample randomizers from the full field,
            // only from 128-bit strings.
            randomizer = u128::rand(rng).into();
        }
        total_c -= &verifier_param.g.mul(g_multiplier);
        end_timer!(combination_time);

        let to_affine_time = start_timer!(|| "Converting results to affine for pairing");
        let affine_points = E::G1Projective::batch_normalization_into_affine(&[-total_w, total_c]);
        let (total_w, total_c) = (affine_points[0], affine_points[1]);
        end_timer!(to_affine_time);

        let pairing_time = start_timer!(|| "Performing product of pairings");
        let result = E::product_of_pairings(&[
            (total_w.into(), verifier_param.beta_h.into()),
            (total_c.into(), verifier_param.h.into()),
        ])
        .is_one();
        end_timer!(pairing_time);
        end_timer!(check_time, || format!("Result: {}", result));
        Ok(result)
    }
}

fn skip_leading_zeros_and_convert_to_bigints<F: PrimeField, P: UVPolynomial<F>>(
    p: &P,
) -> (usize, Vec<F::BigInt>) {
    let mut num_leading_zeros = 0;
    while num_leading_zeros < p.coeffs().len() && p.coeffs()[num_leading_zeros].is_zero() {
        num_leading_zeros += 1;
    }
    let coeffs = convert_to_bigints(&p.coeffs()[num_leading_zeros..]);
    (num_leading_zeros, coeffs)
}

fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    let to_bigint_time = start_timer!(|| "Converting polynomial coeffs to bigints");
    let coeffs = p.iter().map(|s| s.into_repr()).collect::<Vec<_>>();
    end_timer!(to_bigint_time);
    coeffs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StructuredReferenceString;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_poly::univariate::DensePolynomial;
    use ark_std::{log2, test_rng, UniformRand};

    fn end_to_end_test_template<E>() -> Result<(), PCSErrors>
    where
        E: PairingEngine,
    {
        let rng = &mut test_rng();
        for _ in 0..100 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let log_degree = log2(degree) as usize;
            let pp = KZGUnivariatePCS::<E>::gen_srs_for_testing(rng, log_degree)?;
            let (ck, vk) = pp.trim(log_degree)?;
            let p = <DensePolynomial<E::Fr> as UVPolynomial<E::Fr>>::rand(degree, rng);
            let comm = KZGUnivariatePCS::<E>::commit(&ck, &p)?;
            let point = E::Fr::rand(rng);
            let value = p.evaluate(&point);
            let proof = KZGUnivariatePCS::<E>::open(&ck, &p, &point)?;
            assert!(
                KZGUnivariatePCS::<E>::verify(&vk, &comm, &point, &value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}",
                degree,
                p.degree(),
            );
        }
        Ok(())
    }

    fn linear_polynomial_test_template<E>() -> Result<(), PCSErrors>
    where
        E: PairingEngine,
    {
        let rng = &mut test_rng();
        for _ in 0..100 {
            let degree = 50;
            let log_degree = log2(degree) as usize;

            let pp = KZGUnivariatePCS::<E>::gen_srs_for_testing(rng, log_degree)?;
            let (ck, vk) = pp.trim(log_degree)?;
            let p = <DensePolynomial<E::Fr> as UVPolynomial<E::Fr>>::rand(degree, rng);
            let comm = KZGUnivariatePCS::<E>::commit(&ck, &p)?;
            let point = E::Fr::rand(rng);
            let value = p.evaluate(&point);
            let proof = KZGUnivariatePCS::<E>::open(&ck, &p, &point)?;
            assert!(
                KZGUnivariatePCS::<E>::verify(&vk, &comm, &point, &value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}",
                degree,
                p.degree(),
            );
        }
        Ok(())
    }

    fn batch_check_test_template<E>() -> Result<(), PCSErrors>
    where
        E: PairingEngine,
    {
        let rng = &mut test_rng();
        for _ in 0..10 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let log_degree = log2(degree) as usize;
            let pp = KZGUnivariatePCS::<E>::gen_srs_for_testing(rng, log_degree)?;
            let (ck, vk) = pp.trim(log_degree)?;
            let mut comms = Vec::new();
            let mut values = Vec::new();
            let mut points = Vec::new();
            let mut proofs = Vec::new();
            for _ in 0..10 {
                let p = <DensePolynomial<E::Fr> as UVPolynomial<E::Fr>>::rand(degree, rng);
                let comm = KZGUnivariatePCS::<E>::commit(&ck, &p)?;
                let point = E::Fr::rand(rng);
                let value = p.evaluate(&point);
                let proof = KZGUnivariatePCS::<E>::open(&ck, &p, &point)?;

                assert!(KZGUnivariatePCS::<E>::verify(
                    &vk, &comm, &point, &value, &proof
                )?);
                comms.push(comm);
                values.push(value);
                points.push(point);
                proofs.push(proof);
            }
            // assert!(KZGUnivariatePCS::<E>::batch_verify(
            //     &vk, &comms, &points, &values, &proofs, rng
            // )?);
        }
        Ok(())
    }

    #[test]
    fn end_to_end_test() {
        end_to_end_test_template::<Bls12_381>().expect("test failed for bls12-381");
    }

    #[test]
    fn linear_polynomial_test() {
        linear_polynomial_test_template::<Bls12_381>().expect("test failed for bls12-381");
    }
    #[test]
    fn batch_check_test() {
        batch_check_test_template::<Bls12_381>().expect("test failed for bls12-381");
    }
}
