//! Main module for the Product Check protocol

use crate::{
    errors::PolyIOPErrors,
    prod_check::util::{compute_product_poly, prove_zero_check},
    PolyIOP, ZeroCheck,
};
use ark_ec::PairingEngine;
use ark_ff::{One, PrimeField, Zero};
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use pcs::prelude::PolynomialCommitmentScheme;
use std::rc::Rc;
use transcript::IOPTranscript;

mod util;

/// A product-check proves that two n-variate multilinear polynomials `f(x),
/// g(x)` satisfy:
/// \prod_{x \in {0,1}^n} f(x) = \prod_{x \in {0,1}^n} g(x)
///
/// A ProductCheck is derived from ZeroCheck.
///
/// Prover steps:
/// 0. assume that the commitments of `f` and `g` has been pushed to the
/// transcript before 1. build `prod(x0, ..., x_n)` from f and g,
///    such that `prod(0, x1, ..., xn)` equals `f/g` over domain {0,1}^n
/// 2. push commitments of `prod(x)` to the transcript,
///    and `generate_challenge` from current transcript (generate alpha)
/// 3. generate the zerocheck proof for the virtual polynomial
///    prod(1, x) - prod(x, 0) * prod(x, 1) + alpha * (f(x) - prod(0, x) * g(x))
///
/// Verifier steps:
/// 1. Extract commitments of `f(x)`, `g(x)`, `prod(x)` from the proof, push
/// them to the transcript (done by the snark caller)
/// 2. `generate_challenge` from current transcript (generate alpha)
/// 3. `verify` to verify the zerocheck proof and generate the subclaim for
/// polynomial evaluations
pub trait ProductCheck<E, PCS>: ZeroCheck<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
{
    type ProductCheckSubClaim;
    type ProductProof;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a ProductCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// ProductCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Generate a proof for product check.
    ///
    /// Inputs:
    /// - fx: the numerator multilinear polynomial
    /// - gx: the denominator multilinear polynomial
    /// - transcript: the IOP transcript
    ///
    /// Outputs
    /// - the product check proof
    ///
    /// Cost: O(N)
    fn prove(
        fx: &DenseMultilinearExtension<E::Fr>,
        gx: &DenseMultilinearExtension<E::Fr>,
        transcript: &mut IOPTranscript<E::Fr>,
        pk: &PCS::ProverParam,
    ) -> Result<Self::ProductProof, PolyIOPErrors>;

    /// Verify that for witness multilinear polynomials f(x), g(x)
    /// it holds that `\prod_{x \in {0,1}^n} f(x) = \prod_{x \in {0,1}^n} g(x)`
    fn verify(
        proof: &Self::ProductProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::ProductCheckSubClaim, PolyIOPErrors>;
}

/// A product check subclaim consists of
/// - A zero check IOP subclaim for
/// `Q(x) = prod(1, x) - prod(x, 0) * prod(x, 1) + alpha * (f(x) - prod(0, x) *
/// g(x))` is 0, consists of the following:
///   - the SubClaim from the SumCheck
///   - the initial challenge r which is used to build eq(x, r) in ZeroCheck
/// - A final query for `prod(1, ..., 1, 0) = claimed_product`.
// Note that this final query is in fact a constant that
// is independent from the proof. So we should avoid
// (de)serialize it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProductCheckSubClaim<F: PrimeField, ZC: ZeroCheck<F>> {
    // the SubClaim from the ZeroCheck
    zero_check_sub_claim: ZC::ZeroCheckSubClaim,
    // final query which consists of
    // - the vector `(1, ..., 1, 0)`
    // The expected final query evaluation is 1
    final_query: (Vec<F>, F),
}

/// A product check proof consists of
/// - a zerocheck proof
/// - a product polynomial commitment
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProductProof<E: PairingEngine, PCS: PolynomialCommitmentScheme<E>, ZC: ZeroCheck<E::Fr>>
{
    zero_check_proof: ZC::ZeroCheckProof,
    prod_x_comm: PCS::Commitment,
}

impl<E, PCS> ProductCheck<E, PCS> for PolyIOP<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E, Polynomial = Rc<DenseMultilinearExtension<E::Fr>>>,
{
    type ProductCheckSubClaim = ProductCheckSubClaim<E::Fr, Self>;
    type ProductProof = ProductProof<E, PCS, Self>;

    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<E::Fr>::new(b"Initializing ProductCheck transcript")
    }

    fn prove(
        fx: &DenseMultilinearExtension<E::Fr>,
        gx: &DenseMultilinearExtension<E::Fr>,
        transcript: &mut IOPTranscript<E::Fr>,
        pk: &PCS::ProverParam,
    ) -> Result<Self::ProductProof, PolyIOPErrors> {
        let start = start_timer!(|| "prod_check prove");

        if fx.num_vars != gx.num_vars {
            return Err(PolyIOPErrors::InvalidParameters(
                "fx and gx have different number of variables".to_string(),
            ));
        }

        // compute the product polynomial
        let prod_x = compute_product_poly(fx, gx)?;

        // generate challenge
        let prod_x_comm = PCS::commit(pk, &Rc::new(prod_x.clone()))?;
        transcript.append_serializable_element(b"prod(x)", &prod_x_comm)?;
        let alpha = transcript.get_and_append_challenge(b"alpha")?;

        // build the zero-check proof
        let (zero_check_proof, _) = prove_zero_check(fx, gx, &prod_x, &alpha, transcript)?;

        end_timer!(start);

        Ok(ProductProof {
            zero_check_proof,
            prod_x_comm,
        })
    }

    fn verify(
        proof: &Self::ProductProof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::ProductCheckSubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "prod_check verify");

        // update transcript
        transcript.append_serializable_element(b"prod(x)", &proof.prod_x_comm)?;
        // TODO: add alpha to `aux_info`
        let _alpha = transcript.get_and_append_challenge(b"alpha")?;

        // invoke the zero check on the iop_proof
        let zero_check_sub_claim =
            <Self as ZeroCheck<E::Fr>>::verify(&proof.zero_check_proof, aux_info, transcript)?;
        let mut final_query = vec![E::Fr::one(); aux_info.num_variables];
        final_query[aux_info.num_variables - 1] = E::Fr::zero();
        let final_eval = E::Fr::one();

        end_timer!(start);

        Ok(ProductCheckSubClaim {
            zero_check_sub_claim,
            final_query: (final_query, final_eval),
        })
    }
}

#[cfg(test)]
mod test {}
