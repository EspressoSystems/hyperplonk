//! Main module for the Permutation Check protocol

use self::util::computer_num_and_denom;
use crate::{errors::PolyIOPErrors, prelude::ProductCheck, PolyIOP};
use ark_ec::PairingEngine;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use pcs::PolynomialCommitmentScheme;
use std::rc::Rc;
use transcript::IOPTranscript;

/// A permutation subclaim consists of
/// - the SubClaim from the ProductCheck
/// - Challenges beta and gamma
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PermutationCheckSubClaim<E, PCS, PC>
where
    E: PairingEngine,
    PC: ProductCheck<E, PCS>,
    PCS: PolynomialCommitmentScheme<E>,
{
    /// the SubClaim from the ProductCheck
    pub product_check_sub_claim: PC::ProductCheckSubClaim,
    /// Challenges beta and gamma
    pub challenges: (E::Fr, E::Fr),
}

pub mod util;

/// A PermutationCheck is derived from ZeroCheck.
///
/// A Permutation Check IOP takes the following steps:
///
/// Inputs:
/// - f(x)
/// - g(x)
/// - permutation s_perm(x)
pub trait PermutationCheck<E, PCS>: ProductCheck<E, PCS>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
{
    type PermutationCheckSubClaim;
    type PermutationProof;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Inputs:
    /// - f(x)
    /// - g(x)
    /// - permutation s_perm(x)
    /// Outputs:
    /// - a permutation check proof proving that g is a permutation of f under
    ///   s_perm
    /// - the Q(x) polynomial build during product check
    ///
    /// Cost: O(N)
    fn prove(
        pcs_param: &PCS::ProverParam,
        fx: &Self::Polynomial,
        gx: &Self::Polynomial,
        s_perm: &Self::Polynomial,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::PermutationProof, Self::Polynomial), PolyIOPErrors>;

    /// Verify that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::PermutationProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors>;
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
impl<E, PCS> PermutationCheck<E, PCS> for PolyIOP<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E, Polynomial = Rc<DenseMultilinearExtension<E::Fr>>>,
{
    /// A Permutation SubClaim is indeed a ZeroCheck SubClaim that consists of
    /// - the SubClaim from the SumCheck
    /// - the initial challenge r which is used to build eq(x, r)
    type PermutationCheckSubClaim = PermutationCheckSubClaim<E, PCS, Self>;
    type PermutationProof = Self::ProductCheckProof;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<E::Fr>::new(b"Initializing PermutationCheck transcript")
    }

    /// Inputs:
    /// - f(x)
    /// - g(x)
    /// - permutation s_perm(x)
    /// Outputs:
    /// - a permutation check proof proving that g is a permutation of f under
    ///   s_perm
    /// - the Q(x) polynomial build during product check
    ///
    /// Cost: O(N)
    fn prove(
        pcs_param: &PCS::ProverParam,
        fx: &Self::Polynomial,
        gx: &Self::Polynomial,
        s_perm: &Self::Polynomial,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::PermutationProof, Self::Polynomial), PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check prove");
        if fx.num_vars != gx.num_vars {
            return Err(PolyIOPErrors::InvalidParameters(
                "fx and gx have different number of variables".to_string(),
            ));
        }

        if fx.num_vars != s_perm.num_vars {
            return Err(PolyIOPErrors::InvalidParameters(
                "fx and s_perm have different number of variables".to_string(),
            ));
        }

        // generate challenge `beta` and `gamma` from current transcript
        let beta = transcript.get_and_append_challenge(b"beta")?;
        let gamma = transcript.get_and_append_challenge(b"gamma")?;
        let (numerator, denominator) = computer_num_and_denom(&beta, &gamma, fx, gx, s_perm)?;

        // invoke product check on numerator and denominator
        let (proof, poly) =
            <Self as ProductCheck<E, PCS>>::prove(pcs_param, &numerator, &denominator, transcript)?;

        end_timer!(start);
        Ok((proof, poly))
    }

    /// Verify that an MLE g(x) is a permutation of an
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::PermutationProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check verify");

        let beta = transcript.get_and_append_challenge(b"beta")?;
        let gamma = transcript.get_and_append_challenge(b"gamma")?;

        // invoke the zero check on the iop_proof
        let product_check_sub_claim =
            <Self as ProductCheck<E, PCS>>::verify(proof, aux_info.num_variables, transcript)?;

        end_timer!(start);
        Ok(PermutationCheckSubClaim {
            product_check_sub_claim,
            challenges: (beta, gamma),
        })
    }
}

#[cfg(test)]
mod test {
    use super::{util::build_prod_partial_eval, PermutationCheck};
    use crate::{
        errors::PolyIOPErrors,
        perm_check::util::computer_num_and_denom,
        prelude::{identity_permutation_mle, random_permutation_mle},
        utils::bit_decompose,
        PolyIOP,
    };
    use arithmetic::{VPAuxInfo, VirtualPolynomial};
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_ff::{One, Zero};
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use pcs::{prelude::KZGMultilinearPCS, PolynomialCommitmentScheme};
    use std::{marker::PhantomData, rc::Rc};

    type KZG = KZGMultilinearPCS<Bls12_381>;

    fn test_permutation_check_helper<E, PCS>(
        pcs_param: &PCS::ProverParam,
        fx: &Rc<DenseMultilinearExtension<E::Fr>>,
        gx: &Rc<DenseMultilinearExtension<E::Fr>>,
        s_perm: &Rc<DenseMultilinearExtension<E::Fr>>,
    ) -> Result<(), PolyIOPErrors>
    where
        E: PairingEngine,
        PCS: PolynomialCommitmentScheme<E, Polynomial = Rc<DenseMultilinearExtension<E::Fr>>>,
    {
        let nv = fx.num_vars;
        let poly_info = VPAuxInfo {
            max_degree: 2,
            num_variables: nv,
            phantom: PhantomData::default(),
        };

        // prover
        let mut transcript = <PolyIOP<E::Fr> as PermutationCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let (proof, q_x) = <PolyIOP<E::Fr> as PermutationCheck<E, PCS>>::prove(
            pcs_param,
            fx,
            gx,
            s_perm,
            &mut transcript,
        )?;

        // verifier
        let mut transcript = <PolyIOP<E::Fr> as PermutationCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let perm_check_sum_claim = <PolyIOP<E::Fr> as PermutationCheck<E, PCS>>::verify(
            &proof,
            &poly_info,
            &mut transcript,
        )?;

        let prod_partial_evals = build_prod_partial_eval(&q_x)?;
        let prod_0x = prod_partial_evals[0].clone();
        let prod_1x = prod_partial_evals[1].clone();
        let prod_x0 = prod_partial_evals[2].clone();
        let prod_x1 = prod_partial_evals[3].clone();

        let (numerator, denominator) = computer_num_and_denom(
            &perm_check_sum_claim.challenges.0,
            &perm_check_sum_claim.challenges.1,
            fx,
            gx,
            &s_perm,
        )?;

        // compute (g(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
        // which is prods[6] * prod[1] * alpha
        let mut q_x = VirtualPolynomial::new_from_mle(&denominator, E::Fr::one());
        q_x.mul_by_mle(
            prod_0x,
            perm_check_sum_claim.product_check_sub_claim.challenge,
        )?;

        //   (g(x) + beta * s_perm(x) + gamma) * prod(0, x) * alpha
        // - (f(x) + beta * s_id(x)   + gamma) * alpha
        q_x.add_mle_list(
            [numerator],
            -perm_check_sum_claim.product_check_sub_claim.challenge,
        )?;

        // Q(x) := prod(1,x) - prod(x, 0) * prod(x, 1)
        //       + alpha * (
        //             (g(x) + beta * s_perm(x) + gamma) * prod(0, x)
        //           - (f(x) + beta * s_id(x)   + gamma))
        q_x.add_mle_list([prod_x0, prod_x1], -E::Fr::one())?;
        q_x.add_mle_list([prod_1x], E::Fr::one())?;

        if q_x
            .evaluate(
                &perm_check_sum_claim
                    .product_check_sub_claim
                    .zero_check_sub_claim
                    .sum_check_sub_claim
                    .point,
            )
            .unwrap()
            != perm_check_sum_claim
                .product_check_sub_claim
                .zero_check_sub_claim
                .sum_check_sub_claim
                .expected_evaluation
        {
            return Err(PolyIOPErrors::InvalidVerifier("wrong subclaim".to_string()));
        };

        // test q_x is a 0 over boolean hypercube
        for i in 0..1 << nv {
            let bit_sequence = bit_decompose(i, nv);
            let eval: Vec<E::Fr> = bit_sequence
                .iter()
                .map(|x| E::Fr::from(*x as u64))
                .collect();
            let res = q_x.evaluate(&eval).unwrap();
            if !res.is_zero() {}
        }
        Ok(())
    }

    fn test_permutation_check(nv: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        let srs = KZGMultilinearPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = KZGMultilinearPCS::<Bls12_381>::trim(&srs, nv + 1, Some(nv + 1))?;

        {
            // good path: w is a permutation of w itself under the identify map
            let w = Rc::new(DenseMultilinearExtension::rand(nv, &mut rng));
            // s_perm is the identity map
            let s_perm = identity_permutation_mle(nv);
            test_permutation_check_helper::<Bls12_381, KZG>(&pcs_param, &w, &w, &s_perm)?;
        }

        {
            // bad path 1: w is a not permutation of w itself under a random map
            let w = Rc::new(DenseMultilinearExtension::rand(nv, &mut rng));
            // s_perm is a random map
            let s_perm = random_permutation_mle(nv, &mut rng);

            if nv == 1 {
                test_permutation_check_helper::<Bls12_381, KZG>(&pcs_param, &w, &w, &s_perm)?;
            } else {
                assert!(test_permutation_check_helper::<Bls12_381, KZG>(
                    &pcs_param, &w, &w, &s_perm
                )
                .is_err());
            }
        }

        {
            // bad path 2: f is a not permutation of g under a identity map
            let f = Rc::new(DenseMultilinearExtension::rand(nv, &mut rng));
            let g = Rc::new(DenseMultilinearExtension::rand(nv, &mut rng));
            // s_perm is the identity map
            let s_perm = identity_permutation_mle(nv);

            assert!(
                test_permutation_check_helper::<Bls12_381, KZG>(&pcs_param, &f, &g, &s_perm)
                    .is_err()
            );
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
