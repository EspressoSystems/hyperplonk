mod prover;
mod verifier;

use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
pub use prover::ProverState;
use std::rc::Rc;
pub use verifier::VerifierState;

use crate::{
    errors::PolyIOPErrors,
    structs::{DomainInfo, IOPProof, SubClaim},
    sum_check::SumCheck,
    transcript::IOPTranscript,
    virtual_poly::VirtualPolynomial,
    PolyIOP,
};

pub trait ZeroCheck<F: PrimeField> {
    type Proof;
    type PolyList;
    type DomainInfo;
    type SubClaim;
    type Transcript;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a ZeroCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// ZeroCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    fn prove(
        poly: &Self::PolyList,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors>;

    /// verify the claimed sum using the proof
    fn verify(
        proof: &Self::Proof,
        domain_info: &Self::DomainInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, PolyIOPErrors>;
}

impl<F: PrimeField> ZeroCheck<F> for PolyIOP<F> {
    type Proof = IOPProof<F>;
    type PolyList = VirtualPolynomial<F>;
    type DomainInfo = DomainInfo<F>;
    type SubClaim = SubClaim<F>;
    type Transcript = IOPTranscript<F>;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a ZeroCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// ZeroCheck prover/verifier.
    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<F>::new(b"Initializing ZeroCheck transcript")
    }

    fn prove(
        poly: &Self::PolyList,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors> {
        let length = poly.domain_info.num_variables;
        let r = transcript.get_and_append_challenge_vectors(b"vector r", length)?;

        let f_hat = build_f_hat(poly, r.as_ref());

        <Self as SumCheck<F>>::prove(&f_hat, transcript)
    }

    /// verify the claimed sum using the proof
    fn verify(
        proof: &Self::Proof,
        domain_info: &Self::DomainInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        <Self as SumCheck<F>>::verify(F::zero(), proof, domain_info, transcript)
    }
}

// Input poly f(x) and a random vector r, output
//      \hat f(x) = \sum_{x_i \in eval_x} f(x_i) eq(x, r)
// where
//      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
fn build_f_hat<F: PrimeField>(poly: &VirtualPolynomial<F>, r: &[F]) -> VirtualPolynomial<F> {
    assert_eq!(poly.domain_info.num_variables, r.len());
    let mut res = poly.clone();
    let eq_x_r = build_eq_x_r(r);
    res.add_product(eq_x_r, F::one());

    res
}

// Evaluate
//      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
// over r, which is
//      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
fn build_eq_x_r<F: PrimeField>(r: &[F]) -> Vec<Rc<DenseMultilinearExtension<F>>> {
    let num_var = r.len();

    let mut res = vec![];
    for (i, &ri) in r.iter().enumerate() {
        // we want to build a polynomial for (x_i * r_i + (1-x_i)*(1-r_i)).
        // we compute the i-th evaluation, i.e., x = (0,...0, 1, 0..,0) where x[i] = 1
        // as:
        //  1 - r_1,
        //  1 - r_2,
        //  ...
        //  1 - r_{i-1}
        //  r_i
        //  1 - r_{i+1}
        //  ...
        //  r_{num_var}
        let one_minus_ri = F::one() - ri;
        let mut current_eval = vec![];
        for j in 0..num_var {
            if i == j {
                current_eval.push(ri);
            } else {
                current_eval.push(one_minus_ri);
            }
        }
        res.push(Rc::new(DenseMultilinearExtension::from_evaluations_vec(
            num_var,
            current_eval,
        )))
    }

    res
}

#[cfg(test)]
mod test {

    use super::ZeroCheck;
    use crate::{virtual_poly::test::random_list_of_products, PolyIOP, VirtualPolynomial};
    use ark_bls12_381::Fr;
    use ark_ff::UniformRand;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use std::rc::Rc;

    fn test_polynomial(nv: usize, num_multiplicands_range: (usize, usize), num_products: usize) {
        let mut rng = test_rng();
        let mut transcript = PolyIOP::init_transcript();

        let (poly, _asserted_sum) =
            random_list_of_products::<Fr, _>(nv, num_multiplicands_range, num_products, &mut rng);
        let proof = PolyIOP::prove(&poly, &mut transcript).expect("fail to prove");
        let poly_info = poly.domain_info.clone();
        let mut transcript = PolyIOP::init_transcript();
        let subclaim =
            PolyIOP::verify(&proof, &poly_info, &mut transcript).expect("fail to verify");
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
    }
    #[test]
    fn test_normal_polynomial() {
        let nv = 12;
        let num_multiplicands_range = (4, 9);
        let num_products = 5;

        test_polynomial(nv, num_multiplicands_range, num_products);
    }
    #[test]
    #[should_panic]
    fn zero_polynomial_should_error() {
        let nv = 0;
        let num_multiplicands_range = (4, 13);
        let num_products = 5;

        test_polynomial(nv, num_multiplicands_range, num_products);
    }

    #[test]
    /// Test that the memory usage of shared-reference is linear to number of
    /// unique MLExtensions instead of total number of multiplicands.
    fn test_shared_reference() {
        let mut rng = test_rng();
        let ml_extensions: Vec<_> = (0..5)
            .map(|_| Rc::new(DenseMultilinearExtension::<Fr>::rand(8, &mut rng)))
            .collect();
        let mut poly = VirtualPolynomial::new(8);
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

        let mut transcript = PolyIOP::init_transcript();
        let poly_info = poly.domain_info.clone();
        let proof = PolyIOP::prove(&poly, &mut transcript).expect("fail to prove");

        let mut transcript = PolyIOP::init_transcript();
        let subclaim =
            PolyIOP::verify(&proof, &poly_info, &mut transcript).expect("fail to verify");
        assert!(
            poly.evaluate(&subclaim.point) == subclaim.expected_evaluation,
            "wrong subclaim"
        );
    }
}
