// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

mod errors;
mod multilinear_kzg;
mod structs;
mod univariate_kzg;

pub mod prelude;

use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use errors::PCSError;
use std::{borrow::Borrow, fmt::Debug, hash::Hash};
use transcript::IOPTranscript;

/// This trait defines APIs for polynomial commitment schemes.
/// Note that for our usage of PCS, we do not require the hiding property.
pub trait PolynomialCommitmentScheme<E: Pairing> {
    /// Prover parameters
    type ProverParam: Clone + Sync;
    /// Verifier parameters
    type VerifierParam: Clone + CanonicalSerialize + CanonicalDeserialize;
    /// Structured reference string
    type SRS: Clone + Debug;
    /// Polynomial and its associated types
    type Polynomial: Clone + Debug + Hash + PartialEq + Eq;
    /// Polynomial input domain
    type Point: Clone + Ord + Debug + Sync + Hash + PartialEq + Eq;
    /// Polynomial Evaluation
    type Evaluation: Field;
    /// Commitments
    type Commitment: Clone + CanonicalSerialize + CanonicalDeserialize + Debug + PartialEq + Eq;
    /// Proofs
    type Proof: Clone + CanonicalSerialize + CanonicalDeserialize + Debug + PartialEq + Eq;
    /// Batch proofs
    type BatchProof;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `supported_size` is the maximum degree.
    /// - For multilinear polynomials, `supported_size` is the number of
    ///   variables.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: Rng>(
        rng: &mut R,
        supported_size: usize,
    ) -> Result<Self::SRS, PCSError>;

    /// Trim the universal parameters to specialize the public parameters.
    /// Input both `supported_degree` for univariate and
    /// `supported_num_vars` for multilinear.
    /// ## Note on function signature
    /// Usually, data structure like SRS and ProverParam are huge and users
    /// might wish to keep them in heap using different kinds of smart pointers
    /// (instead of only in stack) therefore our `impl Borrow<_>` interface
    /// allows for passing in any pointer type, e.g.: `trim(srs: &Self::SRS,
    /// ..)` or `trim(srs: Box<Self::SRS>, ..)` or `trim(srs: Arc<Self::SRS>,
    /// ..)` etc.
    fn trim(
        srs: impl Borrow<Self::SRS>,
        supported_degree: Option<usize>,
        supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError>;

    /// Generate a commitment for a polynomial
    /// ## Note on function signature
    /// Usually, data structure like SRS and ProverParam are huge and users
    /// might wish to keep them in heap using different kinds of smart pointers
    /// (instead of only in stack) therefore our `impl Borrow<_>` interface
    /// allows for passing in any pointer type, e.g.: `commit(prover_param:
    /// &Self::ProverParam, ..)` or `commit(prover_param:
    /// Box<Self::ProverParam>, ..)` or `commit(prover_param:
    /// Arc<Self::ProverParam>, ..)` etc.
    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, PCSError>;

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same.
    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomial: &Self::Polynomial,
        point: &Self::Point,
    ) -> Result<(Self::Proof, Self::Evaluation), PCSError>;

    /// Input a list of multilinear extensions, and a same number of points, and
    /// a transcript, compute a multi-opening for all the polynomials.
    fn multi_open(
        _prover_param: impl Borrow<Self::ProverParam>,
        _polynomials: &[Self::Polynomial],
        _points: &[Self::Point],
        _evals: &[Self::Evaluation],
        _transcript: &mut IOPTranscript<E::ScalarField>,
    ) -> Result<Self::BatchProof, PCSError> {
        // the reason we use unimplemented!() is to enable developers to implement the
        // trait without always implementing the batching APIs.
        unimplemented!()
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &Self::Point,
        value: &E::ScalarField,
        proof: &Self::Proof,
    ) -> Result<bool, PCSError>;

    /// Verifies that `value_i` is the evaluation at `x_i` of the polynomial
    /// `poly_i` committed inside `comm`.
    fn batch_verify(
        _verifier_param: &Self::VerifierParam,
        _commitments: &[Self::Commitment],
        _points: &[Self::Point],
        _batch_proof: &Self::BatchProof,
        _transcript: &mut IOPTranscript<E::ScalarField>,
    ) -> Result<bool, PCSError> {
        // the reason we use unimplemented!() is to enable developers to implement the
        // trait without always implementing the batching APIs.
        unimplemented!()
    }
}

/// API definitions for structured reference string
pub trait StructuredReferenceString<E: Pairing>: Sized {
    /// Prover parameters
    type ProverParam;
    /// Verifier parameters
    type VerifierParam;

    /// Extract the prover parameters from the public parameters.
    fn extract_prover_param(&self, supported_size: usize) -> Self::ProverParam;
    /// Extract the verifier parameters from the public parameters.
    fn extract_verifier_param(&self, supported_size: usize) -> Self::VerifierParam;

    /// Trim the universal parameters to specialize the public parameters
    /// for polynomials to the given `supported_size`, and
    /// returns committer key and verifier key.
    ///
    /// - For univariate polynomials, `supported_size` is the maximum degree.
    /// - For multilinear polynomials, `supported_size` is 2 to the number of
    ///   variables.
    ///
    /// `supported_log_size` should be in range `1..=params.log_size`
    fn trim(
        &self,
        supported_size: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError>;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `supported_size` is the maximum degree.
    /// - For multilinear polynomials, `supported_size` is the number of
    ///   variables.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: Rng>(rng: &mut R, supported_size: usize) -> Result<Self, PCSError>;
}
