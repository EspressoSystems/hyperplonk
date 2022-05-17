//! Error module.

use ark_std::string::String;
use displaydoc::Display;
use poly_iop::PolyIOPErrors;

/// A `enum` specifying the possible failure modes of the PCS.
#[derive(Display, Debug)]
pub enum PCSErrors {
    /// Invalid Prover: {0}
    InvalidProver(String),
    /// Invalid Verifier: {0}
    InvalidVerifier(String),
    /// Invalid Proof: {0}
    InvalidProof(String),
    /// Invalid parameters: {0}
    InvalidParameters(String),
    /// An error during (de)serialization: {0}
    SerializationError(ark_serialize::SerializationError),
    /// An error from PolyIop
    PolyIOPError(PolyIOPErrors),
}

impl From<ark_serialize::SerializationError> for PCSErrors {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationError(e)
    }
}

impl From<PolyIOPErrors> for PCSErrors {
    fn from(e: PolyIOPErrors) -> Self {
        Self::PolyIOPError(e)
    }
}
