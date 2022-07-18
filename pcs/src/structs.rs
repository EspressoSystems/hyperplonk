use std::marker::PhantomData;

use ark_ec::PairingEngine;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};

use crate::PCSScheme;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// A commitment is an Affine point.
/// Additionally, if it is a multilinear commitment, it also
/// has a field for number of variables.
pub struct Commitment<E: PairingEngine, PCS: PCSScheme<E>> {
    /// number of variables
    pub num_vars: Option<usize>,
    /// product of g as described by the vRAM paper
    pub g_product: E::G1Affine,
    /// polynomial commitment scheme
    pub(crate) phantom: PhantomData<PCS>,
}
