//! Error module.

use ark_std::string::String;
use displaydoc::Display;
use transcript::TranscriptErrors;

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
    /// Invalid challenge: {0}
    InvalidChallenge(String),
    /// Should not arrive to this point
    ShouldNotArrive,
    /// An error during (de)serialization: {0}
    SerializationErrors(ark_serialize::SerializationError),
    /// Transcript Error: {0}
    TranscriptErrors(TranscriptErrors),
}

impl From<ark_serialize::SerializationError> for PolyIOPErrors {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationErrors(e)
    }
}

impl From<TranscriptErrors> for PolyIOPErrors {
    fn from(e: TranscriptErrors) -> Self {
        Self::TranscriptErrors(e)
    }
}
