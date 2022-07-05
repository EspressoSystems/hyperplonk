//! Error module.

use ark_serialize::SerializationError;
use ark_std::string::String;
use displaydoc::Display;
use pcs::prelude::PCSErrors;
use poly_iop::prelude::PolyIOPErrors;

/// A `enum` specifying the possible failure modes of hyperplonk.
#[allow(dead_code)]
// todo: REMOVE
#[derive(Display, Debug)]
pub enum HyperPlonkErrors {
    /// Invalid Prover: {0}
    InvalidProver(String),
    /// Invalid Verifier: {0}
    InvalidVerifier(String),
    /// Invalid Proof: {0}
    InvalidProof(String),
    /// Invalid parameters: {0}
    InvalidParameters(String),
    /// An error during (de)serialization: {0}
    SerializationError(SerializationError),
    /// PolyIOP error {0}
    PolyIOPErrors(PolyIOPErrors),
    /// PCS error {0}
    PCSErrors(PCSErrors),
}

impl From<SerializationError> for HyperPlonkErrors {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationError(e)
    }
}

impl From<PolyIOPErrors> for HyperPlonkErrors {
    fn from(e: PolyIOPErrors) -> Self {
        Self::PolyIOPErrors(e)
    }
}

impl From<PCSErrors> for HyperPlonkErrors {
    fn from(e: PCSErrors) -> Self {
        Self::PCSErrors(e)
    }
}
