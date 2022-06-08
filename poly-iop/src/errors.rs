//! Error module.

use ark_std::string::String;
use displaydoc::Display;

/// A `enum` specifying the possible failure modes of the PolyIOP.
#[derive(Display, Debug)]
pub enum PolyIOPErrors {
    /// Invalid Prover: {0}
    InvalidProver(String),
    /// Invalid Verifier: {0}
    InvalidVerifier(String),
    /// Invalid Proof: {0}
    InvalidProof(String),
    /// Invalid parameters: {0}
    InvalidParameters(String),
    /// Invalid Transcript: {0}
    InvalidTranscript(String),
    /// Should not arrive to this point
    ShouldNotArrive,
    /// An error during (de)serialization: {0}
    SerializationError(ark_serialize::SerializationError),
}

impl From<ark_serialize::SerializationError> for PolyIOPErrors {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationError(e)
    }
}
