// Copyright (c) 2022 Espresso Systems (espressosys.com)
// This file is part of the Jellyfish library.

// You should have received a copy of the MIT License
// along with the Jellyfish library. If not, see <https://mit-license.org/>.

//! Main module for univariate KZG commitment scheme

use crate::{prelude::Commitment, PCSError, PolynomialCommitmentScheme, StructuredReferenceString};
use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{
    borrow::Borrow,
    end_timer, format,
    marker::PhantomData,
    rand::{CryptoRng, RngCore},
    start_timer,
    string::ToString,
    vec,
    vec::Vec,
    One,
};
use srs::{UnivariateProverParam, UnivariateUniversalParams, UnivariateVerifierParam};

pub(crate) mod srs;

/// KZG Polynomial Commitment Scheme on univariate polynomial.
pub struct UnivariateKzgPCS<E: PairingEngine> {
    #[doc(hidden)]
    phantom: PhantomData<E>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
/// proof of opening
pub struct UnivariateKzgProof<E: PairingEngine> {
    /// Evaluation of quotients
    pub proof: E::G1Affine,
}
/// batch proof
pub type UnivariateKzgBatchProof<E> = Vec<UnivariateKzgProof<E>>;

impl<E: PairingEngine> PolynomialCommitmentScheme<E> for UnivariateKzgPCS<E> {
    // Parameters
    type ProverParam = UnivariateProverParam<E::G1Affine>;
    type VerifierParam = UnivariateVerifierParam<E>;
    type SRS = UnivariateUniversalParams<E>;
    // Polynomial and its associated types
    type Polynomial = DensePolynomial<E::Fr>;
    type Point = E::Fr;
    type Evaluation = E::Fr;
    // Polynomial and its associated types
    type Commitment = Commitment<E>;
    type Proof = UnivariateKzgProof<E>;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `supported_size` is the maximum degree.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore + CryptoRng>(
        rng: &mut R,
        supported_size: usize,
    ) -> Result<Self::SRS, PCSError> {
        Self::SRS::gen_srs_for_testing(rng, supported_size)
    }

    /// Trim the universal parameters to specialize the public parameters.
    /// Input `max_degree` for univariate.
    /// `supported_num_vars` must be None or an error is returned.
    fn trim(
        srs: impl Borrow<Self::SRS>,
        supported_degree: Option<usize>,
        supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        assert!(supported_num_vars.is_none());
        if supported_num_vars.is_some() {
            return Err(PCSError::InvalidParameters(
                "univariate should not receive a num_var param".to_string(),
            ));
        }
        srs.borrow().trim(supported_degree.unwrap())
    }

    /// Generate a commitment for a polynomial
    /// Note that the scheme is not hidding
    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, PCSError> {
        let prover_param = prover_param.borrow();
        let commit_time =
            start_timer!(|| format!("Committing to polynomial of degree {} ", poly.degree()));

        if poly.degree() >= prover_param.powers_of_g.len() {
            return Err(PCSError::InvalidParameters(format!(
                "uni poly degree {} is larger than allowed {}",
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
        Ok(Commitment(commitment))
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomial: &Self::Polynomial,
        point: &Self::Point,
    ) -> Result<(Self::Proof, Self::Evaluation), PCSError> {
        let open_time =
            start_timer!(|| format!("Opening polynomial of degree {}", polynomial.degree()));
        let divisor = Self::Polynomial::from_coefficients_vec(vec![-*point, E::Fr::one()]);

        let witness_time = start_timer!(|| "Computing witness polynomial");
        let witness_polynomial = polynomial / &divisor;
        end_timer!(witness_time);

        let (num_leading_zeros, witness_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(&witness_polynomial);

        let proof = VariableBaseMSM::multi_scalar_mul(
            &prover_param.borrow().powers_of_g[num_leading_zeros..],
            &witness_coeffs,
        )
        .into_affine();

        let eval = polynomial.evaluate(point);

        end_timer!(open_time);
        Ok((Self::Proof { proof }, eval))
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &Self::Point,
        value: &E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSError> {
        let check_time = start_timer!(|| "Checking evaluation");
        let pairing_inputs: Vec<(E::G1Prepared, E::G2Prepared)> = vec![
            (
                (verifier_param.g.mul(value.into_repr())
                    - proof.proof.mul(point.into_repr())
                    - commitment.0.into_projective())
                .into_affine()
                .into(),
                verifier_param.h.into(),
            ),
            (proof.proof.into(), verifier_param.beta_h.into()),
        ];

        let res = E::product_of_pairings(pairing_inputs.iter()).is_one();

        end_timer!(check_time, || format!("Result: {}", res));
        Ok(res)
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
    use ark_std::{test_rng, UniformRand};

    fn end_to_end_test_template<E>() -> Result<(), PCSError>
    where
        E: PairingEngine,
    {
        let rng = &mut test_rng();
        for _ in 0..100 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let pp = UnivariateKzgPCS::<E>::gen_srs_for_testing(rng, degree)?;
            let (ck, vk) = pp.trim(degree)?;
            let p = <DensePolynomial<E::Fr> as UVPolynomial<E::Fr>>::rand(degree, rng);
            let comm = UnivariateKzgPCS::<E>::commit(&ck, &p)?;
            let point = E::Fr::rand(rng);
            let (proof, value) = UnivariateKzgPCS::<E>::open(&ck, &p, &point)?;
            assert!(
                UnivariateKzgPCS::<E>::verify(&vk, &comm, &point, &value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}",
                degree,
                p.degree(),
            );
        }
        Ok(())
    }

    fn linear_polynomial_test_template<E>() -> Result<(), PCSError>
    where
        E: PairingEngine,
    {
        let rng = &mut test_rng();
        for _ in 0..100 {
            let degree = 50;

            let pp = UnivariateKzgPCS::<E>::gen_srs_for_testing(rng, degree)?;
            let (ck, vk) = pp.trim(degree)?;
            let p = <DensePolynomial<E::Fr> as UVPolynomial<E::Fr>>::rand(degree, rng);
            let comm = UnivariateKzgPCS::<E>::commit(&ck, &p)?;
            let point = E::Fr::rand(rng);
            let (proof, value) = UnivariateKzgPCS::<E>::open(&ck, &p, &point)?;
            assert!(
                UnivariateKzgPCS::<E>::verify(&vk, &comm, &point, &value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}",
                degree,
                p.degree(),
            );
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
}
