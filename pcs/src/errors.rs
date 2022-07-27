//! Error module.

use ark_serialize::SerializationError;
use ark_std::string::String;
use displaydoc::Display;
use transcript::TranscriptErrors;

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
    SerializationErrors(SerializationError),
    /// Transcript error {0}
    TranscriptErrors(TranscriptErrors),
}

impl From<SerializationError> for PCSErrors {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationErrors(e)
    }
}

impl From<TranscriptErrors> for PCSErrors {
    fn from(e: TranscriptErrors) -> Self {
        Self::TranscriptErrors(e)
    }
}
