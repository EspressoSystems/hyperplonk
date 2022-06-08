//! Main module for the Permutation Check protocol

use crate::errors::PolyIOPErrors;
use ark_ff::PrimeField;

mod util;

pub trait PermutationCheck<F: PrimeField> {
    type Proof;
    type VirtualPolynomial;
    type VPAuxInfo;
    type SubClaim;
    type Transcript;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a ZeroCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// ZeroCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// initialize the prover to argue for the sum of polynomial over
    /// {0,1}^`num_vars` is zero.
    fn prove(
        poly: &Self::VirtualPolynomial,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors>;

    /// verify the claimed sum using the proof
    fn verify(
        proof: &Self::Proof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, PolyIOPErrors>;
}
