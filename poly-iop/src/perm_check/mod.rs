//! Main module for the Permutation Check protocol

use crate::{
    errors::PolyIOPErrors,
    perm_check::util::{build_prod_partial_eval, compute_prod_0},
    structs::IOPProof,
    utils::get_index,
    PolyIOP, ZeroCheck,
};
use arithmetic::VirtualPolynomial;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_serialize::CanonicalSerialize;
use ark_std::{end_timer, start_timer};
use std::rc::Rc;
use transcript::IOPTranscript;

pub mod util;

/// A PermutationCheck is derived from ZeroCheck.
///
/// A Permutation Check IOP takes the following steps:
///
/// Inputs:
/// - f(x)
/// - g(x)
/// - permutation s_perm(x)
///
/// Steps:
/// 1. `generate_challenge` from current transcript (generate beta, gamma)
/// 2. `compute_product` to build `prod(x)` etc. from f, g and s_perm
/// 3. push a commitment of `prod(x)` to the transcript (done by the snark
/// caller)
/// 4. `update_challenge` with the updated transcript (generate alpha)
/// 5. `prove` to generate the proof
pub trait PermutationCheck<F: PrimeField>: ZeroCheck<F> {
    type PermutationCheckSubClaim;
    type PermutationChallenge;
    type PermutationProof;

    /// Generate the preprocessed polynomial for the permutation check.
    ///
    /// The algorithm takes as input a permutation and outputs a merged
    /// multilinear polynomial s(X0, X1, ..., Xn) such that
    /// - s(0, X1, ..., Xn) = s_id(X1, ..., Xn) (identity permutation
    ///   polynomial)
    /// - s(1, X1, ..., Xn) = s_perm(X1, ..., Xn) (permutation polynomial)
    fn preprocess(
        permutation: &[F],
        aux_info: &Self::VPAuxInfo,
    ) -> Result<DenseMultilinearExtension<F>, PolyIOPErrors>;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Step 1 of the IOP.
    /// Generate challenge beta and gamma from a transcript.
    fn generate_challenge(
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationChallenge, PolyIOPErrors>;

    /// Step 2 of the IOP.
    ///
    /// Input:
    /// - f(x), g(x), s_perm(x) are mle-s
    /// - challenges, consists of beta and gamma
    ///
    /// Output: the evaluations for the following 3 polynomials
    /// - prod(x)
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
    ///   where s_id(x) is an identity permutation
    ///
    /// - numerator is the MLE for `f(x) + \beta s_id(x) + \gamma`
    /// - denominator is the MLE for `g(x) + \beta s_perm(x) + \gamma`
    ///
    /// The caller needs to check num_vars matches in f/g/s_id/s_perm
    /// Cost: linear in N.
    ///
    /// TODO: replace argument `s_perm` with the merged polynomial `s`.
    fn compute_prod_evals(
        challenge: &Self::PermutationChallenge,
        fx: &DenseMultilinearExtension<F>,
        gx: &DenseMultilinearExtension<F>,
        s_perm: &DenseMultilinearExtension<F>,
    ) -> Result<[DenseMultilinearExtension<F>; 3], PolyIOPErrors>;

    /// Step 3 of the IOP.
    /// push a commitment of `prod(x)` to the transcript
    /// IMPORTANT: this step is done by the snark caller
    fn commit_prod_x() {
        unimplemented!()
    }

    /// Step 4 of the IOP.
    /// Update the challenge with alpha; returns an error if
    /// alpha already exists.
    fn update_challenge(
        challenge: &mut Self::PermutationChallenge,
        transcript: &mut Self::Transcript,
        prod_x_binding: &impl CanonicalSerialize,
    ) -> Result<(), PolyIOPErrors>;

    /// Step 5 of the IOP.
    ///
    /// Initialize the prover to argue that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    /// Inputs:
    /// - 3 MLEs from the second step
    /// - challenge: `Self::Challenge` that has been updated
    /// - transcript: a transcript that is used to generate the challenges beta
    ///   and gamma
    /// Cost: O(N)
    fn prove(
        prod_x_and_aux_info: &[DenseMultilinearExtension<F>; 3],
        challenge: &Self::PermutationChallenge,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<(Self::PermutationProof, VirtualPolynomial<F>), PolyIOPErrors>;

    /// Verify that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::PermutationProof,
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
    pub zero_check_sub_claim: ZC::ZeroCheckSubClaim,
    // final query which consists of
    // - the vector `(1, ..., 1, 0)`
    // - the evaluation `1`
    final_query: (Vec<F>, F),
}

#[derive(Debug, Clone)]
pub struct PermutationChallenge<F: PrimeField> {
    pub alpha: Option<F>,
    pub beta: F,
    pub gamma: F,
}

/// A PermutationCheck is derived from ZeroCheck.
///
/// A Permutation Check IOP takes the following steps:
///
/// Inputs:
/// - f(x)
/// - g(x)
/// - permutation s_perm(x)
///
/// Steps:
/// 1. `generate_challenge` from current transcript (generate beta, gamma)
/// 2. `compute_product` to build `prod(x)` etc. from f, g and s_perm
/// 3. push a commitment of `prod(x)` to the transcript (done by the snark
/// caller)
/// 4. `update_challenge` with the updated transcript (generate alpha)
/// 5. `prove` to generate the proof
impl<F: PrimeField> PermutationCheck<F> for PolyIOP<F> {
    /// A Permutation SubClaim is indeed a ZeroCheck SubClaim that consists of
    /// - the SubClaim from the SumCheck
    /// - the initial challenge r which is used to build eq(x, r)
    type PermutationCheckSubClaim = PermutationCheckSubClaim<F, Self>;

    type PermutationProof = Self::SumCheckProof;
    type PermutationChallenge = PermutationChallenge<F>;

    /// Generate the preprocessed polynomial for the permutation check.
    ///
    /// The algorithm takes as input a permutation and outputs a merged
    /// multilinear polynomial s(X0, X1, ..., Xn) such that
    /// - s(0, X1, ..., Xn) = s_id(X1, ..., Xn) (identity permutation
    ///   polynomial)
    /// - s(1, X1, ..., Xn) = s_perm(X1, ..., Xn) (permutation polynomial)
    fn preprocess(
        _permutation: &[F],
        _aux_info: &Self::VPAuxInfo,
    ) -> Result<DenseMultilinearExtension<F>, PolyIOPErrors> {
        unimplemented!();
    }

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<F>::new(b"Initializing PermutationCheck transcript")
    }

    /// Step 1 of the IOP.
    /// Generate challenge beta and gamma from a transcript.
    fn generate_challenge(
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationChallenge, PolyIOPErrors> {
        Ok(Self::PermutationChallenge {
            beta: transcript.get_and_append_challenge(b"beta")?,
            gamma: transcript.get_and_append_challenge(b"gamma")?,
            alpha: None,
        })
    }

    /// Step 2 of the IOP.
    ///
    /// Input:
    /// - f(x), g(x), s_perm(x) are mle-s
    /// - challenges, consists of beta and gamma
    ///
    /// Output: the evaluations for the following 3 polynomials
    /// - prod(x)
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
    ///   where s_id(x) is an identity permutation
    ///
    /// - numerator is the MLE for `f(x) + \beta s_id(x) + \gamma`
    /// - denominator is the MLE for `g(x) + \beta s_perm(x) + \gamma`
    ///
    /// The caller needs to check num_vars matches in f/g/s_id/s_perm
    /// Cost: linear in N.
    ///
    /// TODO: replace argument `s_perm` with the merged polynomial `s`.
    fn compute_prod_evals(
        challenge: &Self::PermutationChallenge,
        fx: &DenseMultilinearExtension<F>,
        gx: &DenseMultilinearExtension<F>,
        s_perm: &DenseMultilinearExtension<F>,
    ) -> Result<[DenseMultilinearExtension<F>; 3], PolyIOPErrors> {
        let start = start_timer!(|| "compute evaluations of prod polynomial");

        if challenge.alpha.is_some() {
            return Err(PolyIOPErrors::InvalidChallenge(
                "alpha is already sampled".to_string(),
            ));
        }

        let num_vars = fx.num_vars;

        // ===================================
        // prod(0, x)
        // ===================================
        let (prod_0x_eval, numerator_eval, denominator_eval) =
            compute_prod_0(&challenge.beta, &challenge.gamma, fx, gx, s_perm)?;

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
        let mut prod_1x_eval = vec![];
        for x in 0..(1 << num_vars) - 1 {
            // sign will decide if the evaluation should be looked up from eval_0x or
            // eval_1x; x_zero_index is the index for the evaluation (x_2, ..., x_n,
            // 0); x_one_index is the index for the evaluation (x_2, ..., x_n, 1);
            let (x_zero_index, x_one_index, sign) = get_index(x, num_vars);
            if !sign {
                prod_1x_eval.push(prod_0x_eval[x_zero_index] * prod_0x_eval[x_one_index]);
            } else {
                // sanity check: if we are trying to look up from the eval_1x table,
                // then the target index must already exist
                if x_zero_index >= prod_1x_eval.len() || x_one_index >= prod_1x_eval.len() {
                    return Err(PolyIOPErrors::ShouldNotArrive);
                }
                prod_1x_eval.push(prod_1x_eval[x_zero_index] * prod_1x_eval[x_one_index]);
            }
        }

        // prod(1, 1, ..., 1) := 0
        prod_1x_eval.push(F::zero());

        // ===================================
        // prod(x)
        // ===================================
        // prod(x)'s evaluation is indeed `e := [eval_0x[..], eval_1x[..]].concat()`
        let eval = [prod_0x_eval.as_slice(), prod_1x_eval.as_slice()].concat();

        let fx = DenseMultilinearExtension::from_evaluations_vec(num_vars + 1, eval);
        let numerator = DenseMultilinearExtension::from_evaluations_vec(num_vars, numerator_eval);
        let denominator =
            DenseMultilinearExtension::from_evaluations_vec(num_vars, denominator_eval);

        end_timer!(start);
        Ok([fx, numerator, denominator])
    }

    /// Step 4 of the IOP.
    /// Update the challenge with alpha; returns an error if
    /// alpha already exists.
    fn update_challenge(
        challenge: &mut Self::PermutationChallenge,
        transcript: &mut Self::Transcript,
        prod_x_binding: &impl CanonicalSerialize,
    ) -> Result<(), PolyIOPErrors> {
        if challenge.alpha.is_some() {
            return Err(PolyIOPErrors::InvalidChallenge(
                "alpha should not be sampled at the current stage".to_string(),
            ));
        }
        transcript.append_serializable_element(b"prod(x)", prod_x_binding)?;
        challenge.alpha = Some(transcript.get_and_append_challenge(b"alpha")?);
        Ok(())
    }

    /// Step 5 of the IOP.
    ///
    /// Initialize the prover to argue that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    /// Inputs:
    /// - 3 MLEs from the second step
    /// - challenge: `Self::Challenge` that has been updated
    /// - transcript: a transcript that is used to generate the challenges beta
    ///   and gamma
    /// Cost: O(N)
    fn prove(
        prod_x_and_aux_info: &[DenseMultilinearExtension<F>; 3],
        challenge: &Self::PermutationChallenge,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<(Self::PermutationProof, VirtualPolynomial<F>), PolyIOPErrors> {
        let alpha = match challenge.alpha {
            Some(p) => p,
            None => {
                return Err(PolyIOPErrors::InvalidChallenge(
                    "alpha is not sampled yet".to_string(),
                ))
            },
        };

        let (proof, _q_x) = prove_internal(prod_x_and_aux_info, &alpha, transcript)?;
        Ok((proof, _q_x))
    }

    /// Verify that an MLE g(x) is a permutation of an
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::PermutationProof,
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

/// Step 5 of the IOP.
///
/// Generate a proof to argue that an MLE g(x) is a permutation of
/// MLE f(x) over a permutation given by s_perm.
/// Inputs:
/// - 3 MLEs from the second step
/// - challenge alpha
/// - transcript: a transcript that is used to generate the challenges beta and
///   gamma
///
/// Returns proof and Q(x) for testing purpose.
///
/// Cost: O(N)
fn prove_internal<F: PrimeField>(
    prod_x_and_aux_info: &[DenseMultilinearExtension<F>; 3],
    alpha: &F,
    transcript: &mut IOPTranscript<F>,
) -> Result<(IOPProof<F>, VirtualPolynomial<F>), PolyIOPErrors> {
    let start = start_timer!(|| "Permutation check prove");

    let prod_partial_evals = build_prod_partial_eval(&prod_x_and_aux_info[0])?;
    let prod_0x = Rc::new(prod_partial_evals[0].clone());
    let prod_1x = Rc::new(prod_partial_evals[1].clone());
    let prod_x0 = Rc::new(prod_partial_evals[2].clone());
    let prod_x1 = Rc::new(prod_partial_evals[3].clone());
    let numerator = Rc::new(prod_x_and_aux_info[1].clone());
    let denominator = Rc::new(prod_x_and_aux_info[2].clone());

    // compute (g(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
    // which is prods[6] * prod[1] * alpha
    let mut q_x = VirtualPolynomial::new_from_mle(denominator, F::one());
    q_x.mul_by_mle(prod_0x, *alpha)?;

    //   (g(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
    // - (f(x) + beta * s_id(x)   + gamma) * alpha
    q_x.add_mle_list([numerator], -*alpha)?;

    // Q(x) := prod(1,x) - prod(x, 0) * prod(x, 1)
    //       + alpha * (
    //             (g(x) + beta * s_perm(x) + gamma) * prod(0, x)
    //           - (f(x) + beta * s_id(x)   + gamma))
    q_x.add_mle_list([prod_x0, prod_x1], -F::one())?;
    q_x.add_mle_list([prod_1x], F::one())?;

    let iop_proof = <PolyIOP<F> as ZeroCheck<F>>::prove(&q_x, transcript)?;

    end_timer!(start);
    Ok((iop_proof, q_x))
}

#[cfg(test)]
mod test {
    use super::{util::build_prod_partial_eval, PermutationCheck};
    use crate::{
        errors::PolyIOPErrors,
        perm_check::prove_internal,
        prelude::{identity_permutation_mle, random_permutation_mle},
        structs::IOPProof,
        utils::bit_decompose,
        PolyIOP,
    };
    use arithmetic::{VPAuxInfo, VirtualPolynomial};
    use ark_bls12_381::{Fr, G1Affine};
    use ark_ec::{AffineCurve, ProjectiveCurve};
    use ark_ff::{UniformRand, Zero};
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use std::marker::PhantomData;

    /// This is a mock function to generate some commitment element for testing.
    fn mock_commit<G: AffineCurve>(_f: &DenseMultilinearExtension<G::ScalarField>) -> G {
        let mut rng = test_rng();
        G::Projective::rand(&mut rng).into_affine()
    }

    fn test_permutation_check_helper(
        f: &DenseMultilinearExtension<Fr>,
        g: &DenseMultilinearExtension<Fr>,
        s_perm: &DenseMultilinearExtension<Fr>,
    ) -> Result<(IOPProof<Fr>, VirtualPolynomial<Fr>), PolyIOPErrors> {
        let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;

        let mut challenge =
            <PolyIOP<Fr> as PermutationCheck<Fr>>::generate_challenge(&mut transcript)?;

        let prod_x_and_aux =
            <PolyIOP<Fr> as PermutationCheck<Fr>>::compute_prod_evals(&challenge, f, g, s_perm)?;

        let prod_x_binding: G1Affine = mock_commit(&prod_x_and_aux[0]);

        <PolyIOP<Fr> as PermutationCheck<Fr>>::update_challenge(
            &mut challenge,
            &mut transcript,
            &prod_x_binding,
        )?;
        let alpha = challenge.alpha.unwrap();

        prove_internal(&prod_x_and_aux, &alpha, &mut transcript)
    }

    fn test_permutation_check(nv: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        let poly_info = VPAuxInfo {
            max_degree: 2,
            num_variables: nv,
            phantom: PhantomData::default(),
        };

        {
            // good path: w is a permutation of w itself under the identify map
            let w = DenseMultilinearExtension::rand(nv, &mut rng);

            // s_perm is the identity map
            let s_perm = identity_permutation_mle(nv);

            let (proof, q_x) = test_permutation_check_helper(&w, &w, &s_perm)?;

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

            // test q_x is a 0 over boolean hypercube
            for i in 0..1 << nv {
                let bit_sequence = bit_decompose(i, nv);
                let eval: Vec<Fr> = bit_sequence.iter().map(|x| Fr::from(*x as u64)).collect();
                let res = q_x.evaluate(&eval)?;
                assert!(res.is_zero())
            }
        }

        {
            // bad path 1: w is a not permutation of w itself under a random map
            let w = DenseMultilinearExtension::rand(nv, &mut rng);

            // s_perm is a random map
            let s_perm = random_permutation_mle(nv, &mut rng);

            let (proof, q_x) = test_permutation_check_helper(&w, &w, &s_perm)?;

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

            let (proof, q_x) = test_permutation_check_helper(&f, &g, &s_perm)?;

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

    #[test]
    fn test_compute_prod() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        for num_vars in 2..6 {
            let f = DenseMultilinearExtension::rand(num_vars, &mut rng);
            let g = DenseMultilinearExtension::rand(num_vars, &mut rng);

            let s_id = identity_permutation_mle::<Fr>(num_vars);
            let s_perm = random_permutation_mle(num_vars, &mut rng);

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let challenge =
                <PolyIOP<Fr> as PermutationCheck<Fr>>::generate_challenge(&mut transcript)?;

            let prod_and_aux = <PolyIOP<Fr> as PermutationCheck<Fr>>::compute_prod_evals(
                &challenge, &f, &g, &s_perm,
            )?;
            let prod_partial_eval = build_prod_partial_eval(&prod_and_aux[0])?;

            for i in 0..1 << num_vars {
                let r: Vec<Fr> = bit_decompose(i, num_vars)
                    .iter()
                    .map(|&x| Fr::from(x))
                    .collect();

                let eval = prod_partial_eval[0].evaluate(&r).unwrap();

                let f_eval = f.evaluate(&r).unwrap();
                let g_eval = g.evaluate(&r).unwrap();
                let s_id_eval = s_id.evaluate(&r).unwrap();
                let s_perm_eval = s_perm.evaluate(&r).unwrap();
                let eval_rec = (f_eval + challenge.beta * s_id_eval + challenge.gamma)
                    / (g_eval + challenge.beta * s_perm_eval + challenge.gamma);

                assert_eq!(eval, eval_rec);
            }
        }
        Ok(())
    }
}
