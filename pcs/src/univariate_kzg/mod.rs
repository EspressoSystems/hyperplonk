//! Main module for univariate KZG commitment scheme

use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::{MultilinearExtension, UVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{rand::RngCore, start_timer, end_timer};
use poly_iop::IOPTranscript;
use std::marker::PhantomData;
use ark_ec::msm::VariableBaseMSM;
use crate::{
    prelude::{Commitment, PCSErrors},
    PCSScheme,
};

use self::srs::{ProverParam, VerifierParam};

mod srs;

/// KZG Polynomial Commitment Scheme on univariate polynomial.
pub struct KZGUnivariatePCS<E: PairingEngine> {
    #[doc(hidden)]
    phantom: PhantomData<E>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// proof of opening
pub struct KZGUnivariateOpening<E: PairingEngine> {
    /// Evaluation of quotients
    pub proofs: E::G1Affine,
}

impl<E: PairingEngine> PCSScheme<E> for KZGUnivariatePCS<E> {
    type Polynomial;
    type ProverParam = ProverParam<E::G1Affine>;
    type VerifierParam = VerifierParam<E>;
    type SRS = ();
    type Commitment = Commitment<E, Self>;
    type Proof = KZGUnivariateOpening<E>;
    type BatchProof = ();
    type Transcript = ();

    /// Generate SRS from RNG.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn setup<R: RngCore>(rng: &mut R, num_vars: usize) -> Result<Self::SRS, PCSErrors> {
        todo!()
    }

    /// Generate a commitment for a polynomial
    fn commit(
        prover_param: &Self::ProverParam,
        poly: &impl MultilinearExtension<E::Fr>,
    ) -> Result<Self::Commitment, PCSErrors> {
        // Self::check_degree_is_too_large(polynomial.degree(), powers.size())?;

        let commit_time = start_timer!(|| format!(
            "Committing to polynomial of degree {} ",
            poly.degree(),
        ));

        let (num_leading_zeros, plain_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(poly);

        let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let mut commitment =
            VariableBaseMSM::multi_scalar_mul(&powers.powers_of_g[num_leading_zeros..], &plain_coeffs);
        end_timer!(msm_time);

        let mut randomness = Randomness::<E::Fr, P>::empty();
        if let Some(hiding_degree) = hiding_bound {
            let mut rng = rng.ok_or(Error::MissingRng)?;
            let sample_random_poly_time = start_timer!(|| format!(
                "Sampling a random polynomial of degree {}",
                hiding_degree
            ));

            randomness = Randomness::rand(hiding_degree, false, None, &mut rng);
            Self::check_hiding_bound(
                randomness.blinding_polynomial.degree(),
                powers.powers_of_gamma_g.len(),
            )?;
            end_timer!(sample_random_poly_time);
        }

        let random_ints = convert_to_bigints(&randomness.blinding_polynomial.coeffs());
        let msm_time = start_timer!(|| "MSM to compute commitment to random poly");
        let random_commitment =
            VariableBaseMSM::multi_scalar_mul(&prover_param.g_powers, random_ints.as_slice()).into_affine();
        end_timer!(msm_time);

        commitment.add_assign_mixed(&random_commitment);

        end_timer!(commit_time);
        Ok(Commitment(commitment.into()))
    }

    /// Generate a commitment for a list of polynomials
    fn multi_commit(
        prover_param: &Self::ProverParam,
        polys: &[impl MultilinearExtension<E::Fr>],
    ) -> Result<Self::Commitment, PCSErrors> {
        todo!()
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: &Self::ProverParam,
        polynomial: &impl MultilinearExtension<E::Fr>,
        point: &[E::Fr],
    ) -> Result<Self::Proof, PCSErrors> {
        todo!()
    }

    /// Input a list of MLEs, and a same number of points, and a transcript,
    /// compute a multi-opening for all the polynomials.
    #[allow(clippy::type_complexity)]
    // TODO: remove after we KZG-commit q(x)
    fn multi_open(
        prover_param: &Self::ProverParam,
        polynomials: &[impl MultilinearExtension<E::Fr>],
        point: &[&[E::Fr]],
        transcript: &mut Self::Transcript,
    ) -> Result<Self::BatchProof, PCSErrors> {
        todo!()
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &[E::Fr],
        value: &E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSErrors> {
        todo!()
    }

    /// Verifies that `value_i` is the evaluation at `x_i` of the polynomial
    /// `poly_i` committed inside `comm`.
    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        multi_commitment: &Self::Commitment,
        points: &[&[E::Fr]],
        batch_proof: &Self::BatchProof,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<bool, PCSErrors> {
        todo!()
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
    let coeffs = p.iter()
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    end_timer!(to_bigint_time);
    coeffs
}
