//! Main module for the HyperPlonk PolyIOP.

use ark_ff::PrimeField;
use errors::HyperPlonkErrors;
use poly_iop::prelude::{PermutationCheck, SumCheck, ZeroCheck};

mod errors;
mod structs;

/// A trait for HyperPlonk Poly-IOPs.
/// A HyperPlonk is derived from SumChecks, ZeroChecks and PermutationChecks.
pub trait HyperPlonkPIOP<F: PrimeField>: SumCheck<F> + ZeroCheck<F> + PermutationCheck<F> {
    type Parameters;
    type ProvingKey;
    type Proof;
    type SubClaim;

    /// Generate the preprocessed polynomials output by the indexer.
    ///
    /// Inputs:
    /// - `params`: HyperPlonk instance parameters
    /// - `permutation`: the permutation for the copy constraints
    /// - `selectors`: the list of selector vectors for custom gates
    /// Outputs:
    /// - The HyperPlonk proving key, which includes the preprocessed
    ///   polynomials.
    fn preprocess(
        params: &Self::Parameters,
        permutation: &[F],
        selectors: &[&[F]],
    ) -> Result<Self::ProvingKey, HyperPlonkErrors>;

    /// Generate HyperPlonk PIOP proof.
    ///
    /// Inputs:
    /// - `pk`: circuit proving key
    /// - `pub_input`: online public input
    /// - `witness`: witness assignement
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - The HyperPlonk PIOP proof.
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[F],
        witness: &[&[F]],
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, HyperPlonkErrors>;

    /// Verify the HyperPlonk proof and generate the evaluation subclaims to be
    /// checked later by the SNARK verifier.
    ///
    /// Inputs:
    /// - `params`: instance parameter
    /// - `pub_input`: online public input
    /// - `proof`: HyperPlonk PIOP proof
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - Return error if the verification fails, otherwise return the
    ///   evaluation subclaim
    fn verify(
        params: &Self::Parameters,
        pub_input: &[F],
        proof: &Self::Proof,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, HyperPlonkErrors>;
}
