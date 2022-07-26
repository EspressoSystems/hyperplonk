//! Main module for the HyperPlonk PolyIOP.

use crate::{errors::PolyIOPErrors, perm_check::PermutationCheck, zero_check::ZeroCheck};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use std::rc::Rc;
use transcript::IOPTranscript;

/// A trait for HyperPlonk Poly-IOPs
pub trait HyperPlonkPIOP<F: PrimeField> {
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
    ) -> Result<Self::ProvingKey, PolyIOPErrors>;

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
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PolyIOPErrors>;

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
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::SubClaim, PolyIOPErrors>;
}

/// The sub-claim for the HyperPlonk PolyIOP, consists of the following:
///   - the SubClaim for the zero-check PIOP
///   - the SubClaim for the permutation-check PIOP
///   - the SubClaim for public input consistency
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkSubClaim<F: PrimeField, ZC: ZeroCheck<F>, PC: PermutationCheck<F>> {
    /// the SubClaim for the custom gate zerocheck
    pub zero_check_sub_claim: ZC::ZeroCheckSubClaim,
    /// the SubClaim for the permutation check
    pub perm_check_sub_claim: PC::PermutationCheckSubClaim,
    /// the public input consistency check
    pub pub_input_sub_claim: (Vec<F>, F), // (point, expected_eval)
}

/// The proof for the HyperPlonk PolyIOP, consists of the following:
///   - the zero-check proof for checking custom gate-satisfiability
///   - the permutation-check proof for checking the copy constraints
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkProof<F: PrimeField, ZC: ZeroCheck<F>, PC: PermutationCheck<F>> {
    /// the custom gate zerocheck proof
    pub zero_check_proof: ZC::Proof,
    /// the permutation check proof for copy constraints
    pub perm_check_proof: PC::Proof,
}

/// The HyperPlonk instance parameters, consists of the following:
///   - the number of variables in the poly-IOP
///   - binary log of the number of public input variables
///   - binary log of the number of selectors
///   - binary log of the number of witness wires
///   - the customized gate function
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkParams {
    /// the number of variables in polys
    pub nv: usize,
    /// binary log of the public input length
    pub log_pub_input_len: usize,
    // binary log of the number of selectors
    pub log_n_selectors: usize,
    /// binary log of the number of witness wires
    pub log_n_wires: usize,
    /// customized gate function
    // TODO: define a struct for it.
    pub gate_func: Vec<Vec<usize>>,
}

/// The HyperPlonk proving key, consists of the following:
///   - the hyperplonk instance parameters
///   - the preprocessed polynomials output by the indexer
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkProvingKey<F: PrimeField> {
    /// hyperplonk instance parameters
    pub params: HyperPlonkParams,
    /// the preprocessed index polynomials
    pub index_oracles: Vec<Rc<DenseMultilinearExtension<F>>>,
}

#[cfg(test)]
mod test {}
