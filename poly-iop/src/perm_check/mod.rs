//! Main module for the Permutation Check protocol

use crate::{
    errors::PolyIOPErrors,
    perm_check::util::{build_q_x, compute_prod_0, identity_permutation_mle},
    structs::IOPProof,
    transcript::IOPTranscript,
    utils::get_index,
    PolyIOP, VirtualPolynomial, ZeroCheck,
};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};

pub mod util;

/// A PermutationCheck is derived from ZeroCheck.
pub trait PermutationCheck<F: PrimeField>: ZeroCheck<F> {
    type PermutationCheckSubClaim;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Compute the following 5 polynomials
    /// - prod(x)
    /// - prod(0, x)
    /// - prod(1, x)
    /// - prod(x, 0)
    /// - prod(x, 1)
    /// - numerator
    /// - denominator
    ///
    /// where
    /// - `prod(0,x) := prod(0, x0, x1, …, xn)` which is the MLE over the
    /// evaluations of the following polynomial on the boolean hypercube
    /// {0,1}^n:
    ///
    ///  (f(x) + \beta s_id(x) + \gamma)/(g(x) + \beta s_perm(x) + \gamma)
    ///
    ///   where
    ///   - beta and gamma are challenges
    ///   - f(x), g(x), s_id(x), s_perm(x) are mle-s
    ///
    /// - `prod(1,x) := prod(x, 0) * prod(x, 1)`
    /// - numerator is the MLE for `f(x) + \beta s_id(x) + \gamma`
    /// - denominator is the MLE for `g(x) + \beta s_perm(x) + \gamma`
    ///
    /// The caller needs to check num_vars matches in f/g/s_id/s_perm
    /// Cost: linear in N.
    fn compute_products(
        beta: &F,
        gamma: &F,
        fx: &DenseMultilinearExtension<F>,
        gx: &DenseMultilinearExtension<F>,
        s_id: &DenseMultilinearExtension<F>,
        s_perm: &DenseMultilinearExtension<F>,
    ) -> Result<[DenseMultilinearExtension<F>; 7], PolyIOPErrors>;

    /// Initialize the prover to argue that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    /// Inputs:
    /// - fx
    /// - gx
    /// - permutation defined by a polynomial
    /// - alpha: a challenge generated after fixing the PC commitment of
    ///   \tilde{prod}.
    /// - transcript: a transcript that is used to generate the challenges beta
    ///   and gamma
    fn prove(
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        s_perm: &Self::MultilinearExtension,
        alpha: &F,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors>;

    /// Verify that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::Proof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors>;
}

/// A permutation subclaim consists of
/// - A zero check IOP subclaim for Q(x) is 0, consists of the following:
///  (See `build_qx` for definition of Q(x).)
///   - the SubClaim from the SumCheck
///   - the initial challenge r which is used to build eq(x, r) in ZeroCheck
/// - A final query for `prod(1, ..., 1, 0) = 1`.
// Note that this final query is in fact a constant that
// is independent from the proof. So we should avoid
// (de)serialize it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PermutationCheckSubClaim<F: PrimeField, ZC: ZeroCheck<F>> {
    // the SubClaim from the ZeroCheck
    zero_check_sub_claim: ZC::ZeroCheckSubClaim,
    // final query which consists of
    // - the vector `(1, ..., 1, 0)`
    // - the evaluation `1`
    final_query: (Vec<F>, F),
}

impl<F: PrimeField> PermutationCheck<F> for PolyIOP<F> {
    /// A Permutation SubClaim is indeed a ZeroCheck SubClaim that consists of
    /// - the SubClaim from the SumCheck
    /// - the initial challenge r which is used to build eq(x, r)
    type PermutationCheckSubClaim = PermutationCheckSubClaim<F, Self>;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<F>::new(b"Initializing PermutationCheck transcript")
    }

    /// Compute the following 5 polynomials
    /// - prod(x)
    /// - prod(0, x)
    /// - prod(1, x)
    /// - prod(x, 0)
    /// - prod(x, 1)
    /// - numerator
    /// - denominator
    ///
    /// where
    /// - `prod(0,x) := prod(0, x0, x1, …, xn)` which is the MLE over the
    /// evaluations of the following polynomial on the boolean hypercube
    /// {0,1}^n:
    ///
    ///  (f(x) + \beta s_id(x) + \gamma)/(g(x) + \beta s_perm(x) + \gamma)
    ///
    ///   where
    ///   - beta and gamma are challenges
    ///   - f(x), g(x), s_id(x), s_perm(x) are mle-s
    ///
    /// - `prod(1,x) := prod(x, 0) * prod(x, 1)`
    /// - numerator is the MLE for `f(x) + \beta s_id(x) + \gamma`
    /// - denominator is the MLE for `g(x) + \beta s_perm(x) + \gamma`
    ///
    /// The caller needs to check num_vars matches in f/g/s_id/s_perm
    /// Cost: linear in N.
    fn compute_products(
        beta: &F,
        gamma: &F,
        fx: &DenseMultilinearExtension<F>,
        gx: &DenseMultilinearExtension<F>,
        s_id: &DenseMultilinearExtension<F>,
        s_perm: &DenseMultilinearExtension<F>,
    ) -> Result<[DenseMultilinearExtension<F>; 7], PolyIOPErrors> {
        let start = start_timer!(|| "compute all prod polynomial");

        let num_vars = fx.num_vars;

        // ===================================
        // prod(0, x)
        // ===================================
        let (prod_0x, numerator, denominator) = compute_prod_0(beta, gamma, fx, gx, s_id, s_perm)?;

        // ===================================
        // prod(1, x)
        // ===================================
        //
        // `prod(1, x)` can be computed via recursing the following formula for 2^n-1
        // times
        //
        // `prod(1, x_1, ..., x_n) :=
        //      prod(x_1, x_2, ..., x_n, 0) * prod(x_1, x_2, ..., x_n, 1)`
        //
        // At any given step, the right hand side of the equation
        // is available via either eval_0x or the current view of eval_1x
        let eval_0x = &prod_0x.evaluations;
        let mut eval_1x = vec![];
        for x in 0..(1 << num_vars) - 1 {
            // sign will decide if the evaluation should be looked up from eval_0x or
            // eval_1x; x_zero_index is the index for the evaluation (x_2, ..., x_n,
            // 0); x_one_index is the index for the evaluation (x_2, ..., x_n, 1);
            let (x_zero_index, x_one_index, sign) = get_index(x, num_vars);
            if !sign {
                eval_1x.push(eval_0x[x_zero_index] * eval_0x[x_one_index]);
            } else {
                // sanity check: if we are trying to look up from the eval_1x table,
                // then the target index must already exist
                if x_zero_index >= eval_1x.len() || x_one_index >= eval_1x.len() {
                    return Err(PolyIOPErrors::ShouldNotArrive);
                }
                eval_1x.push(eval_1x[x_zero_index] * eval_1x[x_one_index]);
            }
        }
        // prod(1, 1, ..., 1) := 0
        eval_1x.push(F::zero());

        // ===================================
        // prod(x)
        // ===================================
        // prod(x)'s evaluation is indeed `e := [eval_0x[..], eval_1x[..]].concat()`
        let eval = [eval_0x.as_slice(), eval_1x.as_slice()].concat();

        // ===================================
        // prod(x, 0) and prod(x, 1)
        // ===================================
        //
        // now we compute eval_x0 and eval_x1
        // eval_0x will be the odd coefficients of eval
        // and eval_1x will be the even coefficients of eval
        let mut eval_x0 = vec![];
        let mut eval_x1 = vec![];
        for (x, &prod_x) in eval.iter().enumerate() {
            if x & 1 == 0 {
                eval_x0.push(prod_x);
            } else {
                eval_x1.push(prod_x);
            }
        }

        let prod_1x = DenseMultilinearExtension::from_evaluations_vec(num_vars, eval_1x);
        let prod_x0 = DenseMultilinearExtension::from_evaluations_vec(num_vars, eval_x0);
        let prod_x1 = DenseMultilinearExtension::from_evaluations_vec(num_vars, eval_x1);
        let prod = DenseMultilinearExtension::from_evaluations_vec(num_vars + 1, eval);

        end_timer!(start);
        Ok([
            prod,
            prod_0x,
            prod_1x,
            prod_x0,
            prod_x1,
            numerator,
            denominator,
        ])
    }

    /// Initialize the prover to argue that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    /// Inputs:
    /// - fx
    /// - gx
    /// - permutation defined by a polynomial
    /// - alpha: a challenge generated after fixing the PC commitment of
    ///   \tilde{prod}.
    /// - transcript: a transcript that is used to generate the challenges beta
    ///   and gamma
    ///
    /// Cost: O(N)
    fn prove(
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        s_perm: &Self::MultilinearExtension,
        alpha: &F,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors> {
        let res = prove_internal(fx, gx, s_perm, alpha, transcript)?;
        Ok(res.0)
    }

    /// Verify that an MLE g(x) is a permutation of an
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::Proof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check verify");

        // invoke the zero check on the iop_proof
        let zero_check_sub_claim = <Self as ZeroCheck<F>>::verify(proof, aux_info, transcript)?;

        let mut final_query = vec![F::one(); aux_info.num_variables];
        final_query[aux_info.num_variables - 1] = F::zero();
        let final_eval = F::one();

        end_timer!(start);

        Ok(PermutationCheckSubClaim {
            zero_check_sub_claim,
            final_query: (final_query, final_eval),
        })
    }
}

/// Initialize the prover to argue that an MLE g(x) is a permutation of
/// MLE f(x) over a permutation given by s_perm.
/// Inputs:
/// - fx
/// - gx
/// - permutation defined by a polynomial
/// - alpha: a challenge generated after fixing the PC commitment of
///   \tilde{prod}.
/// - transcript: a transcript that is used to generate the challenges beta and
///   gamma
///
/// Cost: O(N)
fn prove_internal<F: PrimeField>(
    fx: &DenseMultilinearExtension<F>,
    gx: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
    alpha: &F,
    transcript: &mut IOPTranscript<F>,
) -> Result<(IOPProof<F>, VirtualPolynomial<F>), PolyIOPErrors> {
    let start = start_timer!(|| "Permutation check prove");

    let num_vars = fx.num_vars;
    if num_vars != gx.num_vars || num_vars != s_perm.num_vars {
        return Err(PolyIOPErrors::InvalidParameters(
            "num of variables do not match".to_string(),
        ));
    }

    // identity permutation
    let s_id = identity_permutation_mle::<F>(num_vars);

    // compute q(x)
    let beta = transcript.get_and_append_challenge(b"beta")?;
    let gamma = transcript.get_and_append_challenge(b"gamma")?;

    let q_x = build_q_x(alpha, &beta, &gamma, fx, gx, &s_id, s_perm)?;

    let iop_proof = <PolyIOP<F> as ZeroCheck<F>>::prove(&q_x, transcript)?;

    end_timer!(start);
    Ok((iop_proof, q_x))
}

#[cfg(test)]
mod test {

    use super::PermutationCheck;
    use crate::{
        errors::PolyIOPErrors,
        perm_check::{prove_internal, util::identity_permutation_mle},
        random_permutation_mle,
        virtual_poly::VPAuxInfo,
        PolyIOP,
    };
    use ark_bls12_381::Fr;
    use ark_ff::UniformRand;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use std::marker::PhantomData;

    fn test_permutation_check(nv: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        {
            // good path: w is a permutation of w itself under the identify map
            let w = DenseMultilinearExtension::rand(nv, &mut rng);

            // s_perm is the identity map
            let s_perm = identity_permutation_mle(nv);

            let alpha = Fr::rand(&mut rng);

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let (proof, q_x) = prove_internal(&w, &w, &s_perm, &alpha, &mut transcript)?;

            let poly_info = VPAuxInfo {
                max_degree: 2,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let subclaim =
                <PolyIOP<Fr> as PermutationCheck<Fr>>::verify(&proof, &poly_info, &mut transcript)?
                    .zero_check_sub_claim;
            assert_eq!(
                q_x.evaluate(&subclaim.sum_check_sub_claim.point).unwrap(),
                subclaim.sum_check_sub_claim.expected_evaluation,
                "wrong subclaim"
            );
        }

        {
            // bad path 1: w is a not permutation of w itself under a random map
            let w = DenseMultilinearExtension::rand(nv, &mut rng);

            // s_perm is a random map
            let s_perm = random_permutation_mle(nv, &mut rng);

            let alpha = Fr::rand(&mut rng);

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let (proof, q_x) = prove_internal(&w, &w, &s_perm, &alpha, &mut transcript)?;

            let poly_info = VPAuxInfo {
                max_degree: 2,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            if nv != 1 {
                assert!(<PolyIOP<Fr> as PermutationCheck<Fr>>::verify(
                    &proof,
                    &poly_info,
                    &mut transcript
                )
                .is_err());
            } else {
                // a trivial poly is always a permutation of itself, so the zero check should
                // pass
                let subclaim = <PolyIOP<Fr> as PermutationCheck<Fr>>::verify(
                    &proof,
                    &poly_info,
                    &mut transcript,
                )?
                .zero_check_sub_claim;

                // the evaluation should fail because a different s_perm is used for proof and
                // for w |-> w mapping
                assert_ne!(
                    q_x.evaluate(&subclaim.sum_check_sub_claim.point).unwrap(),
                    subclaim.sum_check_sub_claim.expected_evaluation,
                    "wrong subclaim"
                );
            }
        }

        {
            // bad path 2: f is a not permutation of g under a identity map
            let f = DenseMultilinearExtension::rand(nv, &mut rng);
            let g = DenseMultilinearExtension::rand(nv, &mut rng);

            // s_perm is the identity map
            let s_perm = identity_permutation_mle(nv);

            let alpha = Fr::rand(&mut rng);

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let (proof, q_x) = prove_internal(&f, &g, &s_perm, &alpha, &mut transcript)?;

            let poly_info = VPAuxInfo {
                max_degree: 2,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            if nv != 1 {
                assert!(<PolyIOP<Fr> as PermutationCheck<Fr>>::verify(
                    &proof,
                    &poly_info,
                    &mut transcript
                )
                .is_err());
            } else {
                // a trivial poly is always a permutation of itself, so the zero check should
                // pass
                let subclaim = <PolyIOP<Fr> as PermutationCheck<Fr>>::verify(
                    &proof,
                    &poly_info,
                    &mut transcript,
                )?
                .zero_check_sub_claim;

                // the evaluation should fail because a different s_perm is used for proof and
                // for f |-> g mapping
                assert_ne!(
                    q_x.evaluate(&subclaim.sum_check_sub_claim.point).unwrap(),
                    subclaim.sum_check_sub_claim.expected_evaluation,
                    "wrong subclaim"
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_trivial_polynomial() -> Result<(), PolyIOPErrors> {
        test_permutation_check(1)
    }
    #[test]
    fn test_normal_polynomial() -> Result<(), PolyIOPErrors> {
        test_permutation_check(5)
    }

    #[test]
    fn zero_polynomial_should_error() -> Result<(), PolyIOPErrors> {
        assert!(test_permutation_check(0).is_err());
        Ok(())
    }
}
