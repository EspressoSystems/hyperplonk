//! This module implements the sum check protocol.
//! Currently this is a simple wrapper of the sumcheck protocol
//! from Arkworks.

use crate::{
    errors::PolyIOPErrors,
    poly_list::PolynomialList,
    structs::{AuxInfo, IOPProof, SubClaim},
    transcript::IOPTranscript,
    PolyIOP,
};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;

mod prover;
mod verifier;

pub use prover::ProverState;
pub use verifier::VerifierState;

pub trait SumCheck<F: PrimeField> {
    type Proof;
    type Poly;
    type PolyList;
    type AuxInfo;
    type SubClaim;
    type ProverState;
    type VerifierState;

    /// extract sum from the proof
    fn extract_sum(proof: &Self::Proof) -> F;

    /// generate proof of the sum of polynomial over {0,1}^`num_vars`
    ///
    /// The polynomial is represented by a list of products of polynomials along
    /// with its coefficient that is meant to be added together.
    ///
    /// This data structure of the polynomial is a list of list of
    /// `(coefficient, DenseMultilinearExtension)`.
    /// * Number of products n = `polynomial.products.len()`,
    /// * Number of multiplicands of ith product m_i =
    ///   `polynomial.products[i].1.len()`,
    /// * Coefficient of ith product c_i = `polynomial.products[i].0`
    ///
    /// The resulting polynomial is
    ///
    /// $$\sum_{i=0}^{n}C_i\cdot\prod_{j=0}^{m_i}P_{ij}$$
    fn prove(poly: &Self::PolyList) -> Result<Self::Proof, PolyIOPErrors>;

    /// verify the claimed sum using the proof
    fn verify(
        sum: F,
        proof: &Self::Proof,
        aux_info: &Self::AuxInfo,
    ) -> Result<Self::SubClaim, PolyIOPErrors>;
}

pub trait SumCheckProver<F: PrimeField> {
    type PolyList;
    type ProverState;
    type ProverMessage;

    /// initialize the prover to argue for the sum of polynomial over
    /// {0,1}^`num_vars`
    ///
    /// The polynomial is represented by a list of products of polynomials along
    /// with its coefficient that is meant to be added together.
    ///
    /// This data structure of the polynomial is a list of list of
    /// `(coefficient, DenseMultilinearExtension)`.
    /// * Number of products n = `polynomial.products.len()`,
    /// * Number of multiplicands of ith product m_i =
    ///   `polynomial.products[i].1.len()`,
    /// * Coefficient of ith product c_i = `polynomial.products[i].0`
    ///
    /// The resulting polynomial is
    ///
    /// $$\sum_{i=0}^{n}C_i\cdot\prod_{j=0}^{m_i}P_{ij}$$
    fn prover_init(polynomial: &Self::PolyList) -> Self::ProverState;

    /// receive message from verifier, generate prover message, and proceed to
    /// next round
    ///
    /// Main algorithm used is from section 3.2 of [XZZPS19](https://eprint.iacr.org/2019/317.pdf#subsection.3.2).
    fn prove_round_and_update_state(
        prover_state: &mut Self::ProverState,
        challenge: &Option<F>,
    ) -> Self::ProverMessage;
}

pub trait SumCheckVerifier<F: PrimeField> {
    type AuxInfo;
    type VerifierState;
    type ProverMessage;
    type Challenge;
    type Transcript;
    type SubClaim;

    /// initialize the verifier
    fn verifier_init(index_info: &Self::AuxInfo) -> Self::VerifierState;

    /// Run verifier at current round, given prover message
    ///
    /// Normally, this function should perform actual verification. Instead,
    /// `verify_round` only samples and stores randomness and perform
    /// verifications altogether in `check_and_generate_subclaim` at
    /// the last step.
    fn verify_round_and_update_state(
        verifier_state: &mut VerifierState<F>,
        prover_msg: &Self::ProverMessage,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Challenge, PolyIOPErrors>;

    /// verify the sumcheck phase, and generate the subclaim
    ///
    /// If the asserted sum is correct, then the multilinear polynomial
    /// evaluated at `subclaim.point` is `subclaim.expected_evaluation`.
    /// Otherwise, it is highly unlikely that those two will be equal.
    /// Larger field size guarantees smaller soundness error.
    fn check_and_generate_subclaim(
        verifier_state: &Self::VerifierState,
        asserted_sum: &F,
    ) -> Result<Self::SubClaim, PolyIOPErrors>;
}

impl<F: PrimeField> SumCheck<F> for PolyIOP<F> {
    type Proof = IOPProof<F>;

    type Poly = DenseMultilinearExtension<F>;

    type PolyList = PolynomialList<F>;

    type ProverState = ProverState<F>;

    type AuxInfo = AuxInfo<F>;

    type SubClaim = SubClaim<F>;

    type VerifierState = VerifierState<F>;

    fn extract_sum(proof: &Self::Proof) -> F {
        proof.proofs[0].evaluations[0] + proof.proofs[0].evaluations[1]
    }

    /// generate proof of the sum of polynomial over {0,1}^`num_vars`
    ///
    /// The polynomial is represented by a list of products of polynomials along
    /// with its coefficient that is meant to be added together.
    ///
    /// This data structure of the polynomial is a list of list of
    /// `(coefficient, DenseMultilinearExtension)`.
    /// * Number of products n = `polynomial.products.len()`,
    /// * Number of multiplicands of ith product m_i =
    ///   `polynomial.products[i].1.len()`,
    /// * Coefficient of ith product c_i = `polynomial.products[i].0`
    ///
    /// The resulting polynomial is
    ///
    /// $$\sum_{i=0}^{n}C_i\cdot\prod_{j=0}^{m_i}P_{ij}$$
    fn prove(poly: &Self::PolyList) -> Result<Self::Proof, PolyIOPErrors> {
        let mut fs_rng = IOPTranscript::<F>::new(b"Initializing transcript");
        fs_rng.append_aux_info(&poly.aux_info)?;

        let mut prover_state = Self::prover_init(&poly);
        let mut challenge = None;
        let mut prover_msgs = Vec::with_capacity(poly.aux_info.num_variables);
        for _ in 0..poly.aux_info.num_variables {
            let prover_msg = Self::prove_round_and_update_state(&mut prover_state, &challenge);
            fs_rng.append_prover_message(&prover_msg)?;
            prover_msgs.push(prover_msg);
            challenge = Some(fs_rng.get_and_append_challenge(b"Internal round")?);
        }

        Ok(IOPProof {
            proofs: prover_msgs,
        })
    }

    /// verify the claimed sum using the proof
    fn verify(
        claimed_sum: F,
        proof: &Self::Proof,
        polynomial_info: &Self::AuxInfo,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        let mut fs_rng = IOPTranscript::<F>::new(b"Initializing transcript");
        fs_rng.append_aux_info(&polynomial_info)?;
        let mut verifier_state = Self::verifier_init(polynomial_info);
        for i in 0..polynomial_info.num_variables {
            let prover_msg = proof.proofs.get(i).expect("proof is incomplete");
            fs_rng.append_prover_message(&prover_msg)?;

            Self::verify_round_and_update_state(&mut verifier_state, prover_msg, &mut fs_rng)?;
        }

        Ok(Self::check_and_generate_subclaim(
            &verifier_state,
            &claimed_sum,
        )?)
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use ark_bls12_381::Fr;
    use ark_ff::UniformRand;
    use ark_poly::MultilinearExtension;
    use ark_std::{
        rand::{Rng, RngCore},
        test_rng,
    };
    use std::rc::Rc;

    fn random_product<F: PrimeField, R: RngCore>(
        nv: usize,
        num_multiplicands: usize,
        rng: &mut R,
    ) -> (Vec<Rc<DenseMultilinearExtension<F>>>, F) {
        let mut multiplicands = Vec::with_capacity(num_multiplicands);
        for _ in 0..num_multiplicands {
            multiplicands.push(Vec::with_capacity(1 << nv))
        }
        let mut sum = F::zero();

        for _ in 0..(1 << nv) {
            let mut product = F::one();
            for i in 0..num_multiplicands {
                let val = F::rand(rng);
                multiplicands[i].push(val);
                product *= val;
            }
            sum += product;
        }

        return (
            multiplicands
                .into_iter()
                .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
                .collect(),
            sum,
        );
    }

    fn random_list_of_products<F: PrimeField, R: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut R,
    ) -> (PolynomialList<F>, F) {
        let mut sum = F::zero();
        let mut poly = PolynomialList::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let (product, product_sum) = random_product(nv, num_multiplicands, rng);
            let coefficient = F::rand(rng);
            poly.add_product(product.into_iter(), coefficient);
            sum += product_sum * coefficient;
        }

        (poly, sum)
    }

    fn test_polynomial(nv: usize, num_multiplicands_range: (usize, usize), num_products: usize) {
        let mut rng = test_rng();
        let (poly, asserted_sum) =
            random_list_of_products::<Fr, _>(nv, num_multiplicands_range, num_products, &mut rng);
        let proof = PolyIOP::prove(&poly).expect("fail to prove");
        let poly_info = poly.aux_info.clone();
        let subclaim = PolyIOP::verify(asserted_sum, &proof, &poly_info).expect("fail to verify");
        assert!(
            poly.evaluate(&subclaim.point) == subclaim.expected_evaluation,
            "wrong subclaim"
        );
    }

    fn test_protocol(nv: usize, num_multiplicands_range: (usize, usize), num_products: usize) {
        let mut rng = test_rng();
        let (poly, asserted_sum) =
            random_list_of_products::<Fr, _>(nv, num_multiplicands_range, num_products, &mut rng);
        let poly_info = poly.aux_info.clone();
        let mut prover_state = PolyIOP::prover_init(&poly);
        let mut verifier_state = PolyIOP::verifier_init(&poly_info);
        let mut challenge = None;
        let mut transcript = IOPTranscript::new(b"a test transcript");
        for _ in 0..poly.aux_info.num_variables {
            let prover_message =
                PolyIOP::prove_round_and_update_state(&mut prover_state, &challenge);

            challenge = Some(
                PolyIOP::verify_round_and_update_state(
                    &mut verifier_state,
                    &prover_message,
                    &mut transcript,
                )
                .unwrap(),
            );
        }
        let subclaim = PolyIOP::check_and_generate_subclaim(&verifier_state, &asserted_sum)
            .expect("fail to generate subclaim");
        assert!(
            poly.evaluate(&subclaim.point) == subclaim.expected_evaluation,
            "wrong subclaim"
        );
    }

    #[test]
    fn test_trivial_polynomial() {
        let nv = 1;
        let num_multiplicands_range = (4, 13);
        let num_products = 5;

        test_polynomial(nv, num_multiplicands_range, num_products);
        test_protocol(nv, num_multiplicands_range, num_products);
    }
    #[test]
    fn test_normal_polynomial() {
        let nv = 12;
        let num_multiplicands_range = (4, 9);
        let num_products = 5;

        test_polynomial(nv, num_multiplicands_range, num_products);
        test_protocol(nv, num_multiplicands_range, num_products);
    }
    #[test]
    #[should_panic]
    fn zero_polynomial_should_error() {
        let nv = 0;
        let num_multiplicands_range = (4, 13);
        let num_products = 5;

        test_polynomial(nv, num_multiplicands_range, num_products);
        test_protocol(nv, num_multiplicands_range, num_products);
    }

    #[test]
    fn test_extract_sum() {
        let mut rng = test_rng();
        let (poly, asserted_sum) = random_list_of_products::<Fr, _>(8, (3, 4), 3, &mut rng);

        let proof = PolyIOP::prove(&poly).expect("fail to prove");
        assert_eq!(PolyIOP::extract_sum(&proof), asserted_sum);
    }

    #[test]
    /// Test that the memory usage of shared-reference is linear to number of
    /// unique MLExtensions instead of total number of multiplicands.
    fn test_shared_reference() {
        let mut rng = test_rng();
        let ml_extensions: Vec<_> = (0..5)
            .map(|_| Rc::new(DenseMultilinearExtension::<Fr>::rand(8, &mut rng)))
            .collect();
        let mut poly = PolynomialList::new(8);
        poly.add_product(
            vec![
                ml_extensions[2].clone(),
                ml_extensions[3].clone(),
                ml_extensions[0].clone(),
            ],
            Fr::rand(&mut rng),
        );
        poly.add_product(
            vec![
                ml_extensions[1].clone(),
                ml_extensions[4].clone(),
                ml_extensions[4].clone(),
            ],
            Fr::rand(&mut rng),
        );
        poly.add_product(
            vec![
                ml_extensions[3].clone(),
                ml_extensions[2].clone(),
                ml_extensions[1].clone(),
            ],
            Fr::rand(&mut rng),
        );
        poly.add_product(
            vec![ml_extensions[0].clone(), ml_extensions[0].clone()],
            Fr::rand(&mut rng),
        );
        poly.add_product(vec![ml_extensions[4].clone()], Fr::rand(&mut rng));

        assert_eq!(poly.flattened_ml_extensions.len(), 5);

        // test memory usage for prover
        let prover = PolyIOP::prover_init(&poly);
        assert_eq!(prover.flattened_ml_extensions.len(), 5);
        drop(prover);

        let poly_info = poly.aux_info.clone();
        let proof = PolyIOP::prove(&poly).expect("fail to prove");
        let asserted_sum = PolyIOP::extract_sum(&proof);
        let subclaim = PolyIOP::verify(asserted_sum, &proof, &poly_info).expect("fail to verify");
        assert!(
            poly.evaluate(&subclaim.point) == subclaim.expected_evaluation,
            "wrong subclaim"
        );
    }
}
