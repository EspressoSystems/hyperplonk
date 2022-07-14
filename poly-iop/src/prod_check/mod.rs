//! Main module for the Permutation Check protocol

use crate::{errors::PolyIOPErrors, transcript::IOPTranscript, VirtualPolynomial, ZeroCheck};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;

/// A ProductCheck is derived from ZeroCheck.
///
/// A ProductCheck IOP takes the following steps:
///
/// Inputs:
/// - f(x)
///
/// Prover steps:
/// 1. `compute_product_poly` to build `prod(x0, ..., x_n)` from virtual
/// polynomial f 2. push commitments of `f(x)`, `prod(x)` to the transcript
/// (done by the snark caller)
/// 3. `generate_challenge` from current transcript (generate alpha)
/// 4. `prove` to generate the zerocheck proof for the virtual polynomial
///     prod(1, x) - prod(x, 0) * prod(x, 1) + alpha * (f(x) - prod(0, x))
///
/// Verifier steps:
/// 1. Extract commitments of `f(x)`, `prod(x)` from the proof, push them to the
/// transcript (done by the snark caller)
/// 2. `generate_challenge` from current transcript (generate alpha)
/// 3. `verify` to verify the zerocheck proof and generate the subclaim for
/// polynomial evaluations
pub trait ProductCheck<F: PrimeField>: ZeroCheck<F> {
    type ProductCheckSubClaim;
    type ProductChallenge;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a ProductCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// ProductCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Generate random challenge `alpha` from a transcript.
    fn generate_challenge(
        transcript: &mut Self::Transcript,
    ) -> Result<Self::ProductChallenge, PolyIOPErrors>;

    /// Compute the product polynomial `prod(x)` where
    ///
    ///  - `prod(0,x) := prod(0, x1, …, xn)` is the MLE over the
    /// evaluations of `f(x)` on the boolean hypercube {0,1}^n
    ///
    /// - `prod(1,x)` is a MLE over the evaluations of `prod(x, 0) * prod(x, 1)`
    /// on the boolean hypercube {0,1}^n
    ///
    /// The caller needs to check num_vars matches in f
    /// Cost: linear in N.
    fn compute_product_poly(
        fx: &VirtualPolynomial<F>,
    ) -> Result<DenseMultilinearExtension<F>, PolyIOPErrors>;

    /// Initialize the prover to argue that for a virtual polynomial f(x),
    /// it holds that `s = \prod_{x \in {0,1}^n} f(x)`
    ///
    /// Inputs:
    /// - fx: the virtual polynomial
    /// - prod_x: the product polynomial
    /// - transcript: a transcript that is used to generate the challenges alpha
    /// - claimed_product: the claimed product value
    ///
    /// Cost: O(N)
    fn prove(
        fx: &VirtualPolynomial<F>,
        prod_x: &DenseMultilinearExtension<F>,
        transcript: &mut IOPTranscript<F>,
        claimed_product: F,
    ) -> Result<Self::Proof, PolyIOPErrors>;

    /// Verify that for a witness virtual polynomial f(x),
    /// it holds that `s = \prod_{x \in {0,1}^n} f(x)`
    fn verify(
        proof: &Self::Proof,
        aux_info: &Self::VPAuxInfo,
        transcript: &mut Self::Transcript,
        claimed_product: F,
    ) -> Result<Self::ProductCheckSubClaim, PolyIOPErrors>;
}

/// A product check subclaim consists of
/// - A zero check IOP subclaim for
/// `Q(x) = prod(1, x) - prod(x, 0) * prod(x, 1) + alpha * (f(x) - prod(0, x)`
/// is 0, consists of the following:
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
    // - the evaluation `claimed_product`
    final_query: (Vec<F>, F),
}

/// The random challenges in a product check protocol
#[allow(dead_code)]
pub struct ProductChallenge<F: PrimeField> {
    alpha: F,
}

#[cfg(test)]
mod test {}
