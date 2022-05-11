// TODO: some of the struct is generic for Sum Checks and Zero Checks.
// If so move them to src/structs.rs

use super::SumCheckVerifier;
use crate::{
    errors::PolyIOPErrors,
    structs::{DomainInfo, IOPProverMessage, SubClaim},
    transcript::IOPTranscript,
    PolyIOP,
};
use ark_ff::PrimeField;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

/// Verifier State
pub struct VerifierState<F: PrimeField> {
    round: usize,
    num_vars: usize,
    max_degree: usize,
    finished: bool,
    /// a list storing the univariate polynomial in evaluation form sent by the
    /// prover at each round
    polynomials_received: Vec<Vec<F>>,
    /// a list storing the randomness sampled by the verifier at each round
    challenges: Vec<F>,
}

impl<F: PrimeField> SumCheckVerifier<F> for PolyIOP<F> {
    type DomainInfo = DomainInfo<F>;
    type VerifierState = VerifierState<F>;
    type ProverMessage = IOPProverMessage<F>;
    type Challenge = F;
    type Transcript = IOPTranscript<F>;
    type SubClaim = SubClaim<F>;

    /// initialize the verifier
    fn verifier_init(index_info: &Self::DomainInfo) -> Self::VerifierState {
        VerifierState {
            round: 1,
            num_vars: index_info.num_variables,
            max_degree: index_info.max_degree,
            finished: false,
            polynomials_received: Vec::with_capacity(index_info.num_variables),
            challenges: Vec::with_capacity(index_info.num_variables),
        }
    }

    /// Run verifier at current round, given prover message
    ///
    /// Normally, this function should perform actual verification. Instead,
    /// `verify_round` only samples and stores randomness and perform
    /// verifications altogether in `check_and_generate_subclaim` at
    /// the last step.
    fn verify_round_and_update_state(
        verifier_state: &mut VerifierState<F>,
        prover_msg: &Self::ProverMessage,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Challenge, PolyIOPErrors> {
        if verifier_state.finished {
            return Err(PolyIOPErrors::InvalidVerifier(
                "Incorrect verifier state: Verifier is already finished.".to_string(),
            ));
        }

        // Now, verifier should check if the received P(0) + P(1) = expected. The check
        // is moved to `check_and_generate_subclaim`, and will be done after the
        // last round.

        let challenge = transcript.get_and_append_challenge(b"Internal round")?;
        verifier_state.challenges.push(challenge);
        verifier_state
            .polynomials_received
            .push(prover_msg.evaluations.to_vec());

        // Now, verifier should set `expected` to P(r).
        // This operation is also moved to `check_and_generate_subclaim`,
        // and will be done after the last round.

        if verifier_state.round == verifier_state.num_vars {
            // accept and close
            verifier_state.finished = true;
        } else {
            verifier_state.round += 1;
        }
        Ok(challenge)
    }

    /// verify the sumcheck phase, and generate the subclaim
    ///
    /// If the asserted sum is correct, then the multilinear polynomial
    /// evaluated at `subclaim.point` is `subclaim.expected_evaluation`.
    /// Otherwise, it is highly unlikely that those two will be equal.
    /// Larger field size guarantees smaller soundness error.
    fn check_and_generate_subclaim(
        verifier_state: &Self::VerifierState,
        asserted_sum: &F,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        if !verifier_state.finished {
            return Err(PolyIOPErrors::InvalidVerifier(
                "Incorrect verifier state: Verifier has not finished.".to_string(),
            ));
        }

        if verifier_state.polynomials_received.len() != verifier_state.num_vars {
            return Err(PolyIOPErrors::InvalidVerifier(
                "insufficient rounds".to_string(),
            ));
        }
        let mut expected_vec = verifier_state
            .polynomials_received
            .clone()
            .into_par_iter()
            .zip(verifier_state.challenges.clone().into_par_iter())
            .map(|(evaluations, challenge)| {
                if evaluations.len() != verifier_state.max_degree + 1 {
                    return Err(PolyIOPErrors::InvalidVerifier(format!(
                        "incorrect number of evaluations: {} vs {}",
                        evaluations.len(),
                        verifier_state.max_degree + 1
                    )));
                }
                Ok(interpolate_uni_poly::<F>(&evaluations, challenge))
            })
            .collect::<Result<Vec<_>, PolyIOPErrors>>()?;
        // insert the asserted_sum to the first position of the expected vector
        expected_vec.insert(0, *asserted_sum);

        for (evaluations, &expected) in verifier_state
            .polynomials_received
            .iter()
            .zip(expected_vec.iter())
            .take(verifier_state.num_vars)
        {
            if evaluations[0] + evaluations[1] != expected {
                return Err(PolyIOPErrors::InvalidProof(
                    "Prover message is not consistent with the claim.".to_string(),
                ));
            }
        }

        Ok(SubClaim {
            point: verifier_state.challenges.to_vec(),
            // the last expected value (unchecked) will be included in the subclaim
            expected_evaluation: expected_vec[verifier_state.num_vars],
        })
    }
}

/// interpolate a uni-variate degree-`p_i.len()-1` polynomial and evaluate this
/// polynomial at `eval_at`.
pub(crate) fn interpolate_uni_poly<F: PrimeField>(p_i: &[F], eval_at: F) -> F {
    let mut result = F::zero();
    let mut i = F::zero();
    for term in p_i.iter() {
        let mut term = *term;
        let mut j = F::zero();
        for _ in 0..p_i.len() {
            if j != i {
                term = term * (eval_at - j) / (i - j)
            }
            j += F::one();
        }
        i += F::one();
        result += term;
    }

    result
}
