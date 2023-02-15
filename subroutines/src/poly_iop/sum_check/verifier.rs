// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! Verifier subroutines for a SumCheck protocol.

use super::{SumCheckSubClaim, SumCheckVerifier};
use crate::poly_iop::{
    errors::PolyIOPErrors,
    structs::{IOPProverMessage, IOPVerifierState},
};
use arithmetic::VPAuxInfo;
use ark_ff::PrimeField;
use ark_std::{end_timer, start_timer};
use transcript::IOPTranscript;

#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

impl<F: PrimeField> SumCheckVerifier<F> for IOPVerifierState<F> {
    type VPAuxInfo = VPAuxInfo<F>;
    type ProverMessage = IOPProverMessage<F>;
    type Challenge = F;
    type Transcript = IOPTranscript<F>;
    type SumCheckSubClaim = SumCheckSubClaim<F>;

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
    ) -> Result<Self::SumCheckSubClaim, PolyIOPErrors> {
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
        Ok(SumCheckSubClaim {
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
/// TODO: The quadratic term can be removed by precomputing the lagrange
/// coefficients.
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
    // i.e., for any i, the one before this will be derived from
    //  denom[i-1] = denom[i] * (len-i) / i
    //
    // that is, we only need to store
    // - the last denom for i = len-1, and
    // - the ratio between current step and fhe last step, which is the product of
    //   (len-i) / i from all previous steps and we store this product as a fraction
    //   number to reduce field divisions.

    // We know
    //  - 2^61 < factorial(20) < 2^62
    //  - 2^122 < factorial(33) < 2^123
    // so we will be able to compute the ratio
    //  - for len <= 20 with i64
    //  - for len <= 33 with i128
    //  - for len >  33 with BigInt
    if p_i.len() <= 20 {
        let last_denominator = F::from(u64_factorial(len - 1));
        let mut ratio_numerator = 1i64;
        let mut ratio_denominator = 1u64;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                -F::from((-ratio_numerator) as u64)
            } else {
                F::from(ratio_numerator as u64)
            };

            res += p_i[i] * prod * F::from(ratio_denominator)
                / (last_denominator * ratio_numerator_f * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i64 - i as i64);
                ratio_denominator *= i as u64;
            }
        }
    } else if p_i.len() <= 33 {
        let last_denominator = F::from(u128_factorial(len - 1));
        let mut ratio_numerator = 1i128;
        let mut ratio_denominator = 1u128;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                -F::from((-ratio_numerator) as u128)
            } else {
                F::from(ratio_numerator as u128)
            };

            res += p_i[i] * prod * F::from(ratio_denominator)
                / (last_denominator * ratio_numerator_f * evals[i]);

            // compute denom for the next step is current_denom * (len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i128 - i as i128);
                ratio_denominator *= i as u128;
            }
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
    let mut res = F::one();
    for i in 2..=a {
        res *= F::from(i as u64);
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u128_factorial(a: usize) -> u128 {
    let mut res = 1u128;
    for i in 2..=a {
        res *= i as u128;
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u64_factorial(a: usize) -> u64 {
    let mut res = 1u64;
    for i in 2..=a {
        res *= i as u64;
    }
    res
}

#[cfg(test)]
mod test {
    use super::interpolate_uni_poly;
    use crate::poly_iop::errors::PolyIOPErrors;
    use ark_bls12_381::Fr;
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{vec::Vec, UniformRand};

    #[test]
    fn test_interpolation() -> Result<(), PolyIOPErrors> {
        let mut prng = ark_std::test_rng();

        // test a polynomial with 20 known points, i.e., with degree 19
        let poly = DensePolynomial::<Fr>::rand(20 - 1, &mut prng);
        let evals = (0..20)
            .map(|i| poly.evaluate(&Fr::from(i)))
            .collect::<Vec<Fr>>();
        let query = Fr::rand(&mut prng);

        assert_eq!(poly.evaluate(&query), interpolate_uni_poly(&evals, query)?);

        // test a polynomial with 33 known points, i.e., with degree 32
        let poly = DensePolynomial::<Fr>::rand(33 - 1, &mut prng);
        let evals = (0..33)
            .map(|i| poly.evaluate(&Fr::from(i)))
            .collect::<Vec<Fr>>();
        let query = Fr::rand(&mut prng);

        assert_eq!(poly.evaluate(&query), interpolate_uni_poly(&evals, query)?);

        // test a polynomial with 64 known points, i.e., with degree 63
        let poly = DensePolynomial::<Fr>::rand(64 - 1, &mut prng);
        let evals = (0..64)
            .map(|i| poly.evaluate(&Fr::from(i)))
            .collect::<Vec<Fr>>();
        let query = Fr::rand(&mut prng);

        assert_eq!(poly.evaluate(&query), interpolate_uni_poly(&evals, query)?);

        Ok(())
    }
}
