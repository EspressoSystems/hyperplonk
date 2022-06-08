//! Main module for the Permutation Check protocol

use crate::{
    errors::PolyIOPErrors,
    perm_check::util::{build_q_x, identity_permutation_mle},
    structs::{IOPProof, IOPProverMessage, SubClaim},
    transcript::IOPTranscript,
    virtual_poly::{VPAuxInfo, VirtualPolynomial},
    PolyIOP, ZeroCheck,
};
use ark_ff::PrimeField;
use ark_std::{end_timer, start_timer};

mod util;

pub trait PermutationCheck<F: PrimeField> {
    type Proof;
    type VirtualPolynomial;
    type VPAuxInfo;
    type SubClaim;
    type Transcript;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
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

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PermutationProof<F: PrimeField> {
    // A zero check proof to argue that Q(x) is 0.
    // See `build_qx` for definition of Q(x)
    iop_proof: IOPProof<F>,

    // A final proof that `prod(1,...,1,0) = 1`
    final_proof: IOPProverMessage<F>,
}

impl<F: PrimeField> PermutationCheck<F> for PolyIOP<F> {
    type Proof = PermutationProof<F>;
    type VirtualPolynomial = VirtualPolynomial<F>;
    type VPAuxInfo = VPAuxInfo<F>;
    type SubClaim = SubClaim<F>;
    type Transcript = IOPTranscript<F>;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<F>::new(b"Initializing PermutationCheck transcript")
    }

    /// Initialize the prover to argue for the sum of polynomial f(x) over
    /// {0,1}^`num_vars` is zero.
    ///
    /// f(x) is zero if \hat f(x) := f(x) * eq(x,r) is also a zero polynomial
    /// for a random r sampled from transcript.
    ///
    /// This function will build the \hat f(x) and then invoke the sumcheck
    /// protocol to generate a proof for which the sum of \hat f(x) is zero
    fn prove(
        poly: &Self::VirtualPolynomial,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check prove");

        let num_vars = poly.aux_info.num_variables;

        // sample challenges := [alpha, beta, gamma]
        let challenges = transcript.get_and_append_challenge_vectors(b"vector r", 3)?;
        let s_id = identity_permutation_mle::<F>(num_vars);

        // compute q(x)
        let q_x = build_q_x(
            &challenges[0],
            &challenges[1],
            &challenges[2],
            &poly.flattened_ml_extensions[0],
            &s_id,
            &poly.flattened_ml_extensions[1],
        )?;

        let iop_proof = <Self as ZeroCheck<F>>::prove(&q_x, transcript)?;

        // // Additional query on prod(1, ..., 1, 0) = 1
        // let r = transcript.get_and_append_challenge_vectors(b"vector r", num_vars)?;

        end_timer!(start);
        Ok(Self::Proof {
            iop_proof,
            final_proof: IOPProverMessage::default(),
        })
    }

    /// Verify that the polynomial's sum is zero using the proof.
    /// Return a Self::Subclaim that consists of the
    ///
    /// - a Subclaim that the sum is zero
    /// - the initial challenge `r` that is used to build `eq(x, r)`
    ///
    /// This function will check that \hat f(x)'s sum is zero. It does not check
    /// `\hat f(x)` is build correctly. The caller needs to makes sure that
    /// `\hat f(x) = f(x) * eq(x, r)`
    fn verify(
        proof: &Self::Proof,
        fx_aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check verify");
        let _subclaim = <Self as ZeroCheck<F>>::verify(&proof.iop_proof, &fx_aux_info, transcript)?;

        end_timer!(start);

        todo!()
    }
}
