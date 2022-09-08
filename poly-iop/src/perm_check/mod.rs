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

/// A PermutationCheck w.r.t. `(f, g, perm)`
/// proves that g is a permutation of f under
/// permutation `perm`
/// It is derived from ProductCheck.
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
    /// - the product polynomial build during product check
    ///
    /// Cost: O(N)
    fn prove(
        pcs_param: &PCS::ProverParam,
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        s_perm: &Self::MultilinearExtension,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::PermutationProof, Self::MultilinearExtension), PolyIOPErrors>;

    /// Verify that an MLE g(x) is a permutation of
    /// MLE f(x) over a permutation given by s_perm.
    fn verify(
        proof: &Self::PermutationProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors>;
}

impl<E, PCS> PermutationCheck<E, PCS> for PolyIOP<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E, Polynomial = Rc<DenseMultilinearExtension<E::Fr>>>,
{
    type PermutationCheckSubClaim = PermutationCheckSubClaim<E, PCS, Self>;
    type PermutationProof = Self::ProductCheckProof;

    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<E::Fr>::new(b"Initializing PermutationCheck transcript")
    }

    fn prove(
        pcs_param: &PCS::ProverParam,
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        s_perm: &Self::MultilinearExtension,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::PermutationProof, Self::MultilinearExtension), PolyIOPErrors> {
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
        let (proof, prod_poly) =
            <Self as ProductCheck<E, PCS>>::prove(pcs_param, &numerator, &denominator, transcript)?;

        end_timer!(start);
        Ok((proof, prod_poly))
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
            <Self as ProductCheck<E, PCS>>::verify(proof, aux_info, transcript)?;

        end_timer!(start);
        Ok(PermutationCheckSubClaim {
            product_check_sub_claim,
            challenges: (beta, gamma),
        })
    }
}

#[cfg(test)]
mod test {
    use super::PermutationCheck;
    use crate::{errors::PolyIOPErrors, PolyIOP};
    use arithmetic::{evaluate_opt, identity_permutation_mle, random_permutation_mle, VPAuxInfo};
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use pcs::{prelude::MultilinearKzgPCS, PolynomialCommitmentScheme};
    use std::{marker::PhantomData, rc::Rc};

    type KZG = MultilinearKzgPCS<Bls12_381>;

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
        let (proof, prod_x) = <PolyIOP<E::Fr> as PermutationCheck<E, PCS>>::prove(
            pcs_param,
            fx,
            gx,
            s_perm,
            &mut transcript,
        )?;

        // verifier
        let mut transcript = <PolyIOP<E::Fr> as PermutationCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let perm_check_sub_claim = <PolyIOP<E::Fr> as PermutationCheck<E, PCS>>::verify(
            &proof,
            &poly_info,
            &mut transcript,
        )?;

        // check product subclaim
        if evaluate_opt(
            &prod_x,
            &perm_check_sub_claim.product_check_sub_claim.final_query.0,
        ) != perm_check_sub_claim.product_check_sub_claim.final_query.1
        {
            return Err(PolyIOPErrors::InvalidVerifier("wrong subclaim".to_string()));
        };

        Ok(())
    }

    fn test_permutation_check(nv: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        let srs = MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = MultilinearKzgPCS::<Bls12_381>::trim(&srs, nv + 1, Some(nv + 1))?;

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
