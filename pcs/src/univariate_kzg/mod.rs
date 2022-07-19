//! Main module for univariate KZG commitment scheme

use crate::{
    prelude::{Commitment, PCSErrors},
    PolynomialCommitmentScheme,
};
use ark_ec::{msm::VariableBaseMSM, PairingEngine, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, start_timer};
use poly_iop::IOPTranscript;
use srs::{ProverParam, UniversalParams, VerifierParam};
use std::marker::PhantomData;

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

impl<E: PairingEngine> PolynomialCommitmentScheme<E> for KZGUnivariatePCS<E> {
    type ProverParam = ProverParam<E::G1Affine>;
    type VerifierParam = VerifierParam<E>;
    type SRS = UniversalParams<E>;
    type Polynomial = DensePolynomial<E::Fr>;
    type Point = E::Fr;
    type Commitment = Commitment<E, Self>;
    type Proof = KZGUnivariateOpening<E>;
    type BatchProof = ();
    type Transcript = ();

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
        Ok(Commitment {
            commitment,
            phantom: PhantomData::default(),
        })
    }

    /// Generate a commitment for a list of polynomials
    fn multi_commit(
        _prover_param: &Self::ProverParam,
        _polys: &[Self::Polynomial],
    ) -> Result<Self::Commitment, PCSErrors> {
        todo!()
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: &Self::ProverParam,
        polynomial: &Self::Polynomial,
        point: &Self::Point,
    ) -> Result<Self::Proof, PCSErrors> {
        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));
        let (num_leading_zeros, witness_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(polynomial);

        let w = VariableBaseMSM::multi_scalar_mul(
            &prover_param.powers_of_g[num_leading_zeros..],
            &witness_coeffs,
        );

        end_timer!(open_time);
        Ok(Self::Proof {
            proofs: w.into_affine(),
        })
    }

    /// Input a list of MLEs, and a same number of points, and a transcript,
    /// compute a multi-opening for all the polynomials.
    fn multi_open(
        _prover_param: &Self::ProverParam,
        _polynomials: &[Self::Polynomial],
        _point: &[&Self::Point],
        _transcript: &mut Self::Transcript,
    ) -> Result<Self::BatchProof, PCSErrors> {
        todo!()
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
        todo!()
    }

    /// Verifies that `value_i` is the evaluation at `x_i` of the polynomial
    /// `poly_i` committed inside `comm`.
    fn batch_verify(
        _verifier_param: &Self::VerifierParam,
        _multi_commitment: &Self::Commitment,
        _points: &[&Self::Point],
        _batch_proof: &Self::BatchProof,
        _transcript: &mut IOPTranscript<E::Fr>,
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
    let coeffs = p.iter().map(|s| s.into_repr()).collect::<Vec<_>>();
    end_timer!(to_bigint_time);
    coeffs
}
