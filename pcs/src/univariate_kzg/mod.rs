//! Main module for univariate KZG commitment scheme


/// KZG Polynomial Commitment Scheme on univariate polynomial.
pub struct KZGUnivariatePC<E: PairingEngine> {
    #[doc(hidden)]
    phantom: PhantomData<E>,
}

// #[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
// /// A commitment is a group element
// pub struct Commitment<E: PairingEngine> {
//     /// product of g as described by the vRAM paper
//     pub g_product: E::G1Affine,
// }