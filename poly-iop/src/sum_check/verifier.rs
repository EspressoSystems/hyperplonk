//! Verifier subroutines for a SumCheck protocol.

use super::SumCheckVerifier;
use crate::{
    errors::PolyIOPErrors,
    structs::{IOPProverMessage, IOPVerifierState, SubClaim},
    transcript::IOPTranscript,
    virtual_poly::VPAuxInfo,
};
use ark_ff::PrimeField;
use ark_std::{end_timer, start_timer};

#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

impl<F: PrimeField> SumCheckVerifier<F> for IOPVerifierState<F> {
    type VPAuxInfo = VPAuxInfo<F>;
    type ProverMessage = IOPProverMessage<F>;
    type Challenge = F;
    type Transcript = IOPTranscript<F>;
    type SubClaim = SubClaim<F>;

    /// Initialize the verifier's state.
    fn verifier_init(index_info: &Self::VPAuxInfo) -> Self {
        let start = start_timer!(|| "sum check verifier init");
        let res = Self {
            round: 1,
            num_vars: index_info.num_variables,
            max_degree: index_info.max_degree,
            finished: false,
            polynomials_received: Vec::with_capacity(index_info.num_variables),
            challenges: Vec::with_capacity(index_info.num_variables),
        };
        end_timer!(start);
        res
    }

    /// Run verifier for the current round, given a prover message.
    ///
    /// Note that `verify_round_and_update_state` only samples and stores
    /// challenges; and update the verifier's state accordingly. The actual
    /// verifications are deferred (in batch) to `check_and_generate_subclaim`
    /// at the last step.
    fn verify_round_and_update_state(
        &mut self,
        prover_msg: &Self::ProverMessage,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Challenge, PolyIOPErrors> {
        let start =
            start_timer!(|| format!("sum check verify {}-th round and update state", self.round));

        if self.finished {
            return Err(PolyIOPErrors::InvalidVerifier(
                "Incorrect verifier state: Verifier is already finished.".to_string(),
            ));
        }

        // In an interactive protocol, the verifier should
        //
        // 1. check if the received 'P(0) + P(1) = expected`.
        // 2. set `expected` to P(r)`
        //
        // When we turn the protocol to a non-interactive one, it is sufficient to defer
        // such checks to `check_and_generate_subclaim` after the last round.

        let challenge = transcript.get_and_append_challenge(b"Internal round")?;
        self.challenges.push(challenge);
        self.polynomials_received
            .push(prover_msg.evaluations.to_vec());

        if self.round == self.num_vars {
            // accept and close
            self.finished = true;
        } else {
            // proceed to the next round
            self.round += 1;
        }

        end_timer!(start);
        Ok(challenge)
    }

    /// This function verifies the deferred checks in the interactive version of
    /// the protocol; and generate the subclaim. Returns an error if the
    /// proof failed to verify.
    ///
    /// If the asserted sum is correct, then the multilinear polynomial
    /// evaluated at `subclaim.point` will be `subclaim.expected_evaluation`.
    /// Otherwise, it is highly unlikely that those two will be equal.
    /// Larger field size guarantees smaller soundness error.
    fn check_and_generate_subclaim(
        &self,
        asserted_sum: &F,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "sum check check and generate subclaim");
        if !self.finished {
            return Err(PolyIOPErrors::InvalidVerifier(
                "Incorrect verifier state: Verifier has not finished.".to_string(),
            ));
        }

        if self.polynomials_received.len() != self.num_vars {
            return Err(PolyIOPErrors::InvalidVerifier(
                "insufficient rounds".to_string(),
            ));
        }

        // the deferred check during the interactive phase:
        // 2. set `expected` to P(r)`
        #[cfg(feature = "parallel")]
        let mut expected_vec = self
            .polynomials_received
            .clone()
            .into_par_iter()
            .zip(self.challenges.clone().into_par_iter())
            .map(|(evaluations, challenge)| {
                if evaluations.len() != self.max_degree + 1 {
                    return Err(PolyIOPErrors::InvalidVerifier(format!(
                        "incorrect number of evaluations: {} vs {}",
                        evaluations.len(),
                        self.max_degree + 1
                    )));
                }
                interpolate_uni_poly::<F>(&evaluations, challenge)
            })
            .collect::<Result<Vec<_>, PolyIOPErrors>>()?;

        #[cfg(not(feature = "parallel"))]
        let mut expected_vec = self
            .polynomials_received
            .clone()
            .into_iter()
            .zip(self.challenges.clone().into_iter())
            .map(|(evaluations, challenge)| {
                if evaluations.len() != self.max_degree + 1 {
                    return Err(PolyIOPErrors::InvalidVerifier(format!(
                        "incorrect number of evaluations: {} vs {}",
                        evaluations.len(),
                        self.max_degree + 1
                    )));
                }
                interpolate_uni_poly::<F>(&evaluations, challenge)
            })
            .collect::<Result<Vec<_>, PolyIOPErrors>>()?;

        // insert the asserted_sum to the first position of the expected vector
        expected_vec.insert(0, *asserted_sum);

        for (evaluations, &expected) in self
            .polynomials_received
            .iter()
            .zip(expected_vec.iter())
            .take(self.num_vars)
        {
            // the deferred check during the interactive phase:
            // 1. check if the received 'P(0) + P(1) = expected`.
            if evaluations[0] + evaluations[1] != expected {
                return Err(PolyIOPErrors::InvalidProof(
                    "Prover message is not consistent with the claim.".to_string(),
                ));
            }
        }
        end_timer!(start);
        Ok(SubClaim {
            point: self.challenges.clone(),
            // the last expected value (not checked within this function) will be included in the
            // subclaim
            expected_evaluation: expected_vec[self.num_vars],
        })
    }
}

/// Interpolate a uni-variate degree-`p_i.len()-1` polynomial and evaluate this
/// polynomial at `eval_at`:
///
///   \sum_{i=0}^len p_i * (\prod_{j!=i} (eval_at - j)/(i-j) )
///
/// This implementation is linear in number of inputs in terms of field
/// operations. It also has a quadratic term in primitive operations which is
/// negligible compared to field operations.
fn interpolate_uni_poly<F: PrimeField>(p_i: &[F], eval_at: F) -> Result<F, PolyIOPErrors> {
    let start = start_timer!(|| "sum check interpolate uni poly opt");

    let len = p_i.len();
    let mut evals = vec![];
    let mut prod = eval_at;
    evals.push(eval_at);

    // `prod = \prod_{j} (eval_at - j)`
    for e in 1..len {
        let tmp = eval_at - F::from(e as u64);
        evals.push(tmp);
        prod *= tmp;
    }
    let mut res = F::zero();
    // we want to compute \prod (j!=i) (i-j) for a given i
    //
    // we start from the last step, which is
    //  denom[len-1] = (len-1) * (len-2) *... * 2 * 1
    // the step before that is
    //  denom[len-2] = (len-2) * (len-3) * ... * 2 * 1 * -1
    // and the step before that is
    //  denom[len-3] = (len-3) * (len-4) * ... * 2 * 1 * -1 * -2
    //
    // that is, we only need to store the current denom[i] (as a fraction number),
    // and the one before this will be derived from
    //  denom[i-1] = denom[i] * (len-i) / i
    //

    // We know
    //  - 2^57 < (factorial(12))^2 < 2^58
    //  - 2^122 < (factorial(20))^2 < 2^123
    //  - 2^122 < factorial(33) < 2^123
    // so we will be able to compute the denom
    //  - for len <= 12 with i64
    //  - for len <= 20 with i128
    //  - for len <= 33 with i128 in quadratic time
    //  - for len >  33 with BigInt
    if p_i.len() <= 12 {
        let mut denom_up = u64_factorial(len - 1) as i64;
        let mut denom_down = 1u64;

        for i in (0..len).rev() {
            let demon_up_f = if denom_up < 0 {
                -F::from((-denom_up) as u64)
            } else {
                F::from(denom_up as u64)
            };

            res += p_i[i] * prod * F::from(denom_down) / (demon_up_f * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                denom_up *= -(len as i64 - i as i64);
                denom_down *= i as u64;
            }
        }
    } else if p_i.len() <= 20 {
        let mut denom_up = u128_factorial(len - 1) as i128;
        let mut denom_down = 1u128;

        for i in (0..len).rev() {
            let demon_up_f = if denom_up < 0 {
                -F::from((-denom_up) as u128)
            } else {
                F::from(denom_up as u128)
            };

            res += p_i[i] * prod * F::from(denom_down) / (demon_up_f * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                denom_up *= -(len as i128 - i as i128);
                denom_down *= i as u128;
            }
        }
    } else if p_i.len() <= 33 {
        for i in 0..len {
            // res += p_i * prod / (divisor * (eval_at - j))
            let divisor = get_divisor(i, len)?;
            let divisor_f = {
                if divisor < 0 {
                    -F::from((-divisor) as u128)
                } else {
                    F::from(divisor as u128)
                }
            };
            res += p_i[i] * prod / (divisor_f * evals[i]);
        }
    } else {
        let mut denom_up = field_factorial::<F>(len - 1);
        let mut denom_down = F::one();

        for i in (0..len).rev() {
            res += p_i[i] * prod * denom_down / (denom_up * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                denom_up *= -F::from((len - i) as u64);
                denom_down *= F::from(i as u64);
            }
        }
    }
    end_timer!(start);
    Ok(res)
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn field_factorial<F: PrimeField>(a: usize) -> F {
    let mut res = 1u64;
    for i in 1..=a {
        res *= i as u64;
    }
    F::from(res)
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u128_factorial(a: usize) -> u128 {
    let mut res = 1u128;
    for i in 1..=a {
        res *= i as u128;
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u64_factorial(a: usize) -> u64 {
    let mut res = 1u64;
    for i in 1..=a {
        res *= i as u64;
    }
    res
}

/// Compute \prod_{j!=i)^len (i-j). This function takes O(n^2) number of
/// primitive operations which is negligible compared to field operations.
// We know
//  - factorial(20) ~ 2^61
//  - factorial(33) ~ 2^123
// so we will be able to store the result for len<=20 with i64;
// for len<=33 with i128; and we do not currently support len>33.
#[inline]
fn get_divisor(i: usize, len: usize) -> Result<i128, PolyIOPErrors> {
    if len <= 20 {
        let mut res = 1i64;
        for j in 0..len {
            if j != i {
                res *= i as i64 - j as i64;
            }
        }
        Ok(res as i128)
    } else if len <= 33 {
        let mut res = 1i128;
        for j in 0..len {
            if j != i {
                res *= i as i128 - j as i128;
            }
        }
        Ok(res)
    } else {
        Err(PolyIOPErrors::InvalidParameters(
            "Do not support number variable > 33".to_string(),
        ))
    }
}
