use ark_ec::PairingEngine;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, Default, PartialEq, Eq)]
/// A commitment is an Affine point.
pub struct Commitment<E: PairingEngine> {
    /// the actual commitment is an affine point.
    pub commitment: E::G1Affine,
}
