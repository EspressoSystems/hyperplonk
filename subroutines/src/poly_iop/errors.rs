// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! Error module.

use crate::pcs::prelude::PCSError;
use arithmetic::ArithErrors;
use ark_std::string::String;
use displaydoc::Display;
use transcript::TranscriptError;

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
    TranscriptErrors(TranscriptError),
    /// Arithmetic Error: {0}
    ArithmeticErrors(ArithErrors),
    /// PCS error {0}
    PCSErrors(PCSError),
}

impl From<ark_serialize::SerializationError> for PolyIOPErrors {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationErrors(e)
    }
}

impl From<TranscriptError> for PolyIOPErrors {
    fn from(e: TranscriptError) -> Self {
        Self::TranscriptErrors(e)
    }
}

impl From<ArithErrors> for PolyIOPErrors {
    fn from(e: ArithErrors) -> Self {
        Self::ArithmeticErrors(e)
    }
}

impl From<PCSError> for PolyIOPErrors {
    fn from(e: PCSError) -> Self {
        Self::PCSErrors(e)
    }
}
