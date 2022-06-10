//! Main module for the Permutation Check protocol

use crate::{
    errors::PolyIOPErrors,
    perm_check::util::{build_q_x, identity_permutation_mle},
    transcript::IOPTranscript,
    PolyIOP, ZeroCheck,
};
use ark_ff::PrimeField;
use ark_std::{end_timer, start_timer};

mod util;

/// A ZeroCheck is derived from ZeroCheck.
pub trait PermutationCheck<F: PrimeField>: ZeroCheck<F> {
    type SubClaim;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Initialize the prover to argue that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    fn prove(
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        s_perm: &Self::MultilinearExtension,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors>;

    /// Verify that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::Proof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, PolyIOPErrors>;
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
    type SubClaim = PermutationCheckSubClaim<F, Self>;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<F>::new(b"Initializing PermutationCheck transcript")
    }

    /// Initialize the prover to argue that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm
    /// Cost: O(N)
    fn prove(
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        s_perm: &Self::MultilinearExtension,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check prove");

        let num_vars = fx.num_vars;
        if num_vars != gx.num_vars || num_vars != s_perm.num_vars {
            return Err(PolyIOPErrors::InvalidParameters(
                "num of variables do not match".to_string(),
            ));
        }

        // sample challenges := [alpha, beta, gamma]
        let challenges = transcript.get_and_append_challenge_vectors(b"q(x) challenge", 3)?;

        // identity permutation
        let s_id = identity_permutation_mle::<F>(num_vars);

        // compute q(x)
        let q_x = build_q_x(
            &challenges[0],
            &challenges[1],
            &challenges[2],
            fx,
            gx,
            &s_id,
            s_perm,
        )?;

        let iop_proof = <Self as ZeroCheck<F>>::prove(&q_x, transcript)?;

        end_timer!(start);
        Ok(iop_proof)
    }

    /// Verify that an MLE g(x) is a permutation of an
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::Proof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check verify");

        // invoke the zero check on the iop_proof
        let zero_check_sub_claim = <Self as ZeroCheck<F>>::verify(&proof, aux_info, transcript)?;

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

#[cfg(test)]
mod test {

    use super::PermutationCheck;
    use crate::{
        errors::PolyIOPErrors, perm_check::util::identity_permutation_mle, virtual_poly::VPAuxInfo,
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
            // good path: zero virtual poly
            let w_vec: Vec<Fr> = (0..(1 << nv)).map(|_| Fr::rand(&mut rng)).collect();
            let w = DenseMultilinearExtension::from_evaluations_vec(nv, w_vec);

            let s_id = identity_permutation_mle(nv);

            let mut transcript = <PolyIOP<Fr> as PermutationCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let proof =
                <PolyIOP<Fr> as PermutationCheck<Fr>>::prove(&w, &w, &s_id, &mut transcript)?;

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
                w.evaluate(&subclaim.sum_check_sub_claim.point).unwrap(),
                subclaim.sum_check_sub_claim.expected_evaluation,
                "wrong subclaim"
            );
        }

        // {
        //     // bad path: random virtual poly whose sum is not zero
        //     let (poly, _sum) =
        //         VirtualPolynomial::rand(nv, num_multiplicands_range, num_products,
        // &mut rng)?;

        //     let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
        //     transcript.append_message(b"testing", b"initializing transcript for
        // testing")?;     let proof = <PolyIOP<Fr> as
        // ZeroCheck<Fr>>::prove(&poly, &mut transcript)?;

        //     let poly_info = poly.aux_info.clone();
        //     let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
        //     transcript.append_message(b"testing", b"initializing transcript for
        // testing")?;

        //     assert!(
        //         <PolyIOP<Fr> as ZeroCheck<Fr>>::verify(&proof, &poly_info, &mut
        // transcript)             .is_err()
        //     );
        // }

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
