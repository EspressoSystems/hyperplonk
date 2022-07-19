use std::marker::PhantomData;

use crate::PolynomialCommitmentScheme;
use ark_ec::PairingEngine;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// A commitment is an Affine point.
pub struct Commitment<E: PairingEngine, PCS: PolynomialCommitmentScheme<E>> {
    /// the actual commitment is an affine point.
    pub commitment: E::G1Affine,
    /// polynomial commitment scheme
    pub(crate) phantom: PhantomData<PCS>,
}
