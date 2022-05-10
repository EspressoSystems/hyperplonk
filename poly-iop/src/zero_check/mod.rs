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
        println!(
            "sum: {}",
            proof.proofs[0].evaluations[0] + proof.proofs[0].evaluations[1]
        );

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

    res.add_product([eq_x_r; 1], F::one());
    // let num_var = r.len();
    // for i in 0..1 << 2 {
    //     let bit_sequence: Vec<F> = bit_decompose(i, num_var)
    //         .iter()
    //         .map(|&x| F::from(x as u64))
    //         .collect();
    //     println!("i {}, eval {}", i, res.evaluate(&bit_sequence))
    // }

    res
}

// Evaluate
//      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
// over r, which is
//      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
fn build_eq_x_r<F: PrimeField>(r: &[F]) -> Rc<DenseMultilinearExtension<F>> {
    // we build eq(x,r) from its evaluations
    // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
    // for example, with num_vars = 4, x is a binary vector of 4, then
    //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
    //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
    //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
    //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
    //  ....
    //  1 1 1 1 -> r0       * r1        * r2        * r3
    // we will need 2^num_var evaluations

    // First, we build array for {1 - r_i}
    let one_minus_r: Vec<F> = r.iter().map(|ri| F::one() - ri).collect();

    let num_var = r.len();
    let mut eval = vec![];

    // TODO: optimize the following code
    // currently, a naive implementation requires num_var * 2^num_var
    // field multiplications.
    for i in 0..1 << num_var {
        let mut current_eval = F::one();
        let bit_sequence = bit_decompose(i, num_var);

        for (&bit, (ri, one_minus_ri)) in bit_sequence.iter().zip(r.iter().zip(one_minus_r.iter()))
        {
            if bit {
                current_eval *= *ri;
            } else {
                current_eval *= *one_minus_ri;
            }
        }
        eval.push(current_eval);
    }
    let res = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_var, eval,
    ));

    res
}

fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push(i & 1 == 1);
        i >>= 1;
    }
    res
}

#[cfg(test)]
mod test {

    use super::ZeroCheck;
    use crate::{virtual_poly::test::random_zero_list_of_products, PolyIOP, VirtualPolynomial};
    use ark_bls12_381::Fr;
    use ark_ff::UniformRand;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use std::rc::Rc;

    fn test_polynomial(nv: usize, num_multiplicands_range: (usize, usize), num_products: usize) {
        let mut rng = test_rng();
        let mut transcript = PolyIOP::init_transcript();

        let poly = random_zero_list_of_products::<Fr, _>(
            nv,
            num_multiplicands_range,
            num_products,
            &mut rng,
        );
        println!("{:?}", poly);

        let proof = PolyIOP::prove(&poly, &mut transcript).expect("fail to prove");
        println!(
            "{:?}",
            proof.proofs[0].evaluations[0] + proof.proofs[0].evaluations[1]
        );

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
        let nv = 16;
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
