mod commit;
mod errors;
mod param;
mod util;

use ark_ec::PairingEngine;
use ark_poly::MultilinearExtension;
use ark_std::rand::RngCore;
use poly_iop::IOPTranscript;
use std::marker::PhantomData;

pub use errors::PCSErrors;
pub use param::{ProverParam, UniversalParams, VerifierParam};

/// KZG Polynomial Commitment Scheme on multilinear extensions.
pub struct KZGMultilinearPC<E: PairingEngine> {
    #[doc(hidden)]
    phantom: PhantomData<E>,
}

pub trait MultilinearCommitmentScheme<E: PairingEngine> {
    type ProverParam;
    type VerifierParam;
    type SRS;
    type Commitment;
    type Proof;

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

    /// On input a list of polynomials and a point `point`, outputs a proof for
    /// the same.
    fn multi_open(
        prover_param: &Self::ProverParam,
        polynomials: &[impl MultilinearExtension<E::Fr>],
        point: &[&[E::Fr]],
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::Proof, Vec<E::Fr>), PCSErrors>;

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &[E::Fr],
        value: E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSErrors>;
}
