// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! Error module.

use ark_std::string::String;
use displaydoc::Display;

/// A `enum` specifying the possible failure modes of the Transcript.
#[derive(Display, Debug)]
pub enum TranscriptError {
    /// Invalid Transcript: {0}
    InvalidTranscript(String),
    /// An error during (de)serialization: {0}
    SerializationError(ark_serialize::SerializationError),
}

impl From<ark_serialize::SerializationError> for TranscriptError {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationError(e)
    }
}
