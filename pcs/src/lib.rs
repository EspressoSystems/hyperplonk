mod errors;
mod multilinear_kzg;
mod structs;

pub mod prelude;

use ark_ec::PairingEngine;
use ark_poly::MultilinearExtension;
use ark_std::rand::RngCore;
use errors::PCSErrors;
use poly_iop::IOPTranscript;

pub trait PCSScheme<E: PairingEngine> {
    type ProverParam;
    type VerifierParam;
    type SRS;
    type Commitment;
    type Proof;
    type BatchProof;
    type Transcript;

    /// Generate SRS from RNG.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn setup<R: RngCore>(rng: &mut R, num_vars: usize) -> Result<Self::SRS, PCSErrors>;

    /// Generate a commitment for a polynomial
    fn commit(
        prover_param: &Self::ProverParam,
        poly: &impl MultilinearExtension<E::Fr>,
    ) -> Result<Self::Commitment, PCSErrors>;

    /// Generate a commitment for a list of polynomials
    fn multi_commit(
        prover_param: &Self::ProverParam,
        polys: &[impl MultilinearExtension<E::Fr>],
    ) -> Result<Self::Commitment, PCSErrors>;

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: &Self::ProverParam,
        polynomial: &impl MultilinearExtension<E::Fr>,
        point: &[E::Fr],
    ) -> Result<Self::Proof, PCSErrors>;

    /// Input a list of MLEs, and a same number of points, and a transcript,
    /// compute a multi-opening for all the polynomials.
    #[allow(clippy::type_complexity)]
    // TODO: remove after we KZG-commit q(x)
    fn multi_open(
        prover_param: &Self::ProverParam,
        polynomials: &[impl MultilinearExtension<E::Fr>],
        point: &[&[E::Fr]],
        transcript: &mut Self::Transcript,
    ) -> Result<Self::BatchProof, PCSErrors>;

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &[E::Fr],
        value: &E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSErrors>;

    /// Verifies that `value_i` is the evaluation at `x_i` of the polynomial
    /// `poly_i` committed inside `comm`.
    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        multi_commitment: &Self::Commitment,
        points: &[&[E::Fr]],
        batch_proof: &Self::BatchProof,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<bool, PCSErrors>;
}

/// API definitions for structured reference string
pub trait StructuredReferenceString<E: PairingEngine>: Sized {
    type ProverParam;
    type VerifierParam;

    /// Extract the prover parameters from the public parameters.
    fn extract_prover_param(&self) -> Self::ProverParam;
    /// Extract the verifier parameters from the public parameters.
    fn extract_verifier_param(&self) -> Self::VerifierParam;

    /// Trim the universal parameters to specialize the public parameters
    /// for multilinear polynomials to the given `supported_num_vars`, and
    /// returns committer key and verifier key. `supported_num_vars` should
    /// be in range `1..=params.num_vars`
    fn trim(
        &self,
        supported_num_vars: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSErrors>;

    /// Build SRS for testing.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore>(rng: &mut R, num_vars: usize) -> Result<Self, PCSErrors>;
}
