//! Structs for polynomials and extensions.

use ark_ff::PrimeField;
use std::marker::PhantomData;

#[derive(Clone, Debug, Default, PartialEq)]
/// Auxiliary information about the multilinear polynomial
pub struct DomainInfo<F: PrimeField> {
    /// max number of multiplicands in each product
    pub max_degree: usize,
    /// number of variables of the polynomial
    pub num_variables: usize,
    /// Associated field
    #[doc(hidden)]
    pub(crate) phantom: PhantomData<F>,
}

/// Subclaim when verifier is convinced
pub struct SubClaim<F: PrimeField> {
    /// the multi-dimensional point that this multilinear extension is evaluated
    /// to
    pub point: Vec<F>,
    /// the expected evaluation
    pub expected_evaluation: F,
}

/// A SumCheck proof is a list of messages from prover to verifier
/// through the interactive protocol.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IOPProof<F: PrimeField> {
    pub proofs: Vec<IOPProverMessage<F>>,
}

/// A message from the prover to the verifier at a given round
/// is a list of evaluations.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IOPProverMessage<F: PrimeField> {
    pub(crate) evaluations: Vec<F>,
}
