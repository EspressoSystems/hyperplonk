//! Main module for the Product Check protocol

use crate::{
    errors::PolyIOPErrors,
    prod_check::util::{compute_product_poly, prove_zero_check},
    zero_check::ZeroCheck,
    PolyIOP,
};
use arithmetic::VPAuxInfo;
use ark_ec::PairingEngine;
use ark_ff::{One, PrimeField, Zero};
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use pcs::PolynomialCommitmentScheme;
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
/// 1. build `prod(x0, ..., x_n)` from f and g,
///    such that `prod(0, x1, ..., xn)` equals `f/g` over domain {0,1}^n
/// 2. push commitments of `prod(x)` to the transcript,
///    and `generate_challenge` from current transcript (generate alpha)
/// 3. generate the zerocheck proof for the virtual polynomial
///    prod(1, x) - prod(x, 0) * prod(x, 1) + alpha * (f(x) - prod(0, x) * g(x))
///
/// Verifier steps:
/// 1. Extract commitments of `prod(x)` from the proof, push
/// them to the transcript
/// 2. `generate_challenge` from current transcript (generate alpha)
/// 3. `verify` to verify the zerocheck proof and generate the subclaim for
/// polynomial evaluations
pub trait ProductCheck<E, PCS>: ZeroCheck<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
{
    type ProductCheckSubClaim;
    type ProductCheckProof;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a ProductCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// ProductCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Generate a proof for product check, showing that witness multilinear
    /// polynomials f(x), g(x) satisfy `\prod_{x \in {0,1}^n} f(x) =
    /// \prod_{x \in {0,1}^n} g(x)`
    ///
    /// Inputs:
    /// - fx: the numerator multilinear polynomial
    /// - gx: the denominator multilinear polynomial
    /// - transcript: the IOP transcript
    /// - pk: PCS committing key
    ///
    /// Outputs
    /// - the product check proof
    /// - the product polynomial (used for testing)
    ///
    /// Cost: O(N)
    fn prove(
        pcs_param: &PCS::ProverParam,
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::ProductCheckProof, Self::MultilinearExtension), PolyIOPErrors>;

    /// Verify that for witness multilinear polynomials f(x), g(x)
    /// it holds that `\prod_{x \in {0,1}^n} f(x) = \prod_{x \in {0,1}^n} g(x)`
    fn verify(
        proof: &Self::ProductCheckProof,
        aux_info: &VPAuxInfo<E::Fr>,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::ProductCheckSubClaim, PolyIOPErrors>;
}

/// A product check subclaim consists of
/// - A zero check IOP subclaim for
/// `Q(x) = prod(1, x) - prod(x, 0) * prod(x, 1) + alpha * (f(x) - prod(0,
/// x) * g(x)) = 0`
/// - The random challenge `alpha`
/// - A final query for `prod(1, ..., 1, 0) = 1`.
// Note that this final query is in fact a constant that
// is independent from the proof. So we should avoid
// (de)serialize it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProductCheckSubClaim<F: PrimeField, ZC: ZeroCheck<F>> {
    // the SubClaim from the ZeroCheck
    pub zero_check_sub_claim: ZC::ZeroCheckSubClaim,
    // final query which consists of
    // - the vector `(1, ..., 1, 0)` (needs to be reversed because Arkwork's MLE uses big-endian
    //   format for points)
    // The expected final query evaluation is 1
    pub final_query: (Vec<F>, F),
    pub alpha: F,
}

/// A product check proof consists of
/// - a zerocheck proof
/// - a product polynomial commitment
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProductCheckProof<
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
    ZC: ZeroCheck<E::Fr>,
> {
    pub zero_check_proof: ZC::ZeroCheckProof,
    pub prod_x_comm: PCS::Commitment,
}

impl<E, PCS> ProductCheck<E, PCS> for PolyIOP<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E, Polynomial = Rc<DenseMultilinearExtension<E::Fr>>>,
{
    type ProductCheckSubClaim = ProductCheckSubClaim<E::Fr, Self>;
    type ProductCheckProof = ProductCheckProof<E, PCS, Self>;

    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<E::Fr>::new(b"Initializing ProductCheck transcript")
    }

    fn prove(
        pcs_param: &PCS::ProverParam,
        fx: &Self::MultilinearExtension,
        gx: &Self::MultilinearExtension,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::ProductCheckProof, Self::MultilinearExtension), PolyIOPErrors> {
        let start = start_timer!(|| "prod_check prove");

        if fx.num_vars != gx.num_vars {
            return Err(PolyIOPErrors::InvalidParameters(
                "fx and gx have different number of variables".to_string(),
            ));
        }

        // compute the product polynomial
        let prod_x = compute_product_poly(fx, gx)?;

        // generate challenge
        let prod_x_comm = PCS::commit(pcs_param, &prod_x)?;
        transcript.append_serializable_element(b"prod(x)", &prod_x_comm)?;
        let alpha = transcript.get_and_append_challenge(b"alpha")?;

        // build the zero-check proof
        let (zero_check_proof, _) = prove_zero_check(fx, gx, &prod_x, &alpha, transcript)?;

        end_timer!(start);

        Ok((
            ProductCheckProof {
                zero_check_proof,
                prod_x_comm,
            },
            prod_x,
        ))
    }

    fn verify(
        proof: &Self::ProductCheckProof,
        aux_info: &VPAuxInfo<E::Fr>,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::ProductCheckSubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "prod_check verify");

        // update transcript and generate challenge
        transcript.append_serializable_element(b"prod(x)", &proof.prod_x_comm)?;
        let alpha = transcript.get_and_append_challenge(b"alpha")?;

        // invoke the zero check on the iop_proof
        // the virtual poly info for Q(x)
        let zero_check_sub_claim =
            <Self as ZeroCheck<E::Fr>>::verify(&proof.zero_check_proof, aux_info, transcript)?;

        // the final query is on prod_x, hence has length `num_vars` + 1
        let mut final_query = vec![E::Fr::one(); aux_info.num_variables + 1];
        // the point has to be reversed because Arkworks uses big-endian.
        final_query[0] = E::Fr::zero();
        let final_eval = E::Fr::one();

        end_timer!(start);

        Ok(ProductCheckSubClaim {
            zero_check_sub_claim,
            final_query: (final_query, final_eval),
            alpha,
        })
    }
}

#[cfg(test)]
mod test {
    use super::ProductCheck;
    use crate::{errors::PolyIOPErrors, PolyIOP};
    use arithmetic::VPAuxInfo;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ec::PairingEngine;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use pcs::{prelude::MultilinearKzgPCS, PolynomialCommitmentScheme};
    use std::{marker::PhantomData, rc::Rc};

    // f and g are guaranteed to have the same product
    fn test_product_check_helper<E, PCS>(
        f: &DenseMultilinearExtension<E::Fr>,
        g: &DenseMultilinearExtension<E::Fr>,
        pcs_param: &PCS::ProverParam,
    ) -> Result<(), PolyIOPErrors>
    where
        E: PairingEngine,
        PCS: PolynomialCommitmentScheme<E, Polynomial = Rc<DenseMultilinearExtension<E::Fr>>>,
    {
        let mut transcript = <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;

        let (proof, prod_x) = <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::prove(
            pcs_param,
            &Rc::new(f.clone()),
            &Rc::new(g.clone()),
            &mut transcript,
        )?;

        let mut transcript = <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;

        let aux_info = VPAuxInfo {
            max_degree: 2,
            num_variables: f.num_vars,
            phantom: PhantomData::default(),
        };
        let prod_subclaim =
            <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::verify(&proof, &aux_info, &mut transcript)?;
        assert_eq!(
            prod_x.evaluate(&prod_subclaim.final_query.0).unwrap(),
            prod_subclaim.final_query.1,
            "different product"
        );

        // bad path
        let mut transcript = <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;

        let h = f + g;
        let (bad_proof, prod_x_bad) = <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::prove(
            pcs_param,
            &Rc::new(f.clone()),
            &Rc::new(h),
            &mut transcript,
        )?;

        let mut transcript = <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;
        let bad_subclaim = <PolyIOP<E::Fr> as ProductCheck<E, PCS>>::verify(
            &bad_proof,
            &aux_info,
            &mut transcript,
        )?;
        assert_ne!(
            prod_x_bad.evaluate(&bad_subclaim.final_query.0).unwrap(),
            bad_subclaim.final_query.1,
            "can't detect wrong proof"
        );

        Ok(())
    }

    fn test_product_check(nv: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        let f: DenseMultilinearExtension<Fr> = DenseMultilinearExtension::rand(nv, &mut rng);
        let mut g = f.clone();
        g.evaluations.reverse();

        let srs = MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = MultilinearKzgPCS::<Bls12_381>::trim(&srs, None, Some(nv + 1))?;

        test_product_check_helper::<Bls12_381, MultilinearKzgPCS<Bls12_381>>(&f, &g, &pcs_param)?;

        Ok(())
    }

    #[test]
    fn test_trivial_polynomial() -> Result<(), PolyIOPErrors> {
        test_product_check(1)
    }
    #[test]
    fn test_normal_polynomial() -> Result<(), PolyIOPErrors> {
        test_product_check(10)
    }
}
