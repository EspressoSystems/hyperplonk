// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! Error module.

use ark_std::string::String;
use displaydoc::Display;

/// A `enum` specifying the possible failure modes of the arithmetics.
#[derive(Display, Debug)]
pub enum ArithErrors {
    /// Invalid parameters: {0}
    InvalidParameters(String),
    /// Should not arrive to this point
    ShouldNotArrive,
    /// An error during (de)serialization: {0}
    SerializationErrors(ark_serialize::SerializationError),
}

impl From<ark_serialize::SerializationError> for ArithErrors {
    fn from(e: ark_serialize::SerializationError) -> Self {
        Self::SerializationErrors(e)
    }
}
