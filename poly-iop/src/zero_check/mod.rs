use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use std::rc::Rc;

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

    /// initialize the prover to argue for the sum of polynomial over
    /// {0,1}^`num_vars` is zero.
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

    /// A ZeroCheck SubClaim consists of
    /// - the SubClaim from the ZeroCheck
    /// - the initial challenge r which is used to build eq(x, r)
    type SubClaim = (SubClaim<F>, Vec<F>);
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

    /// initialize the prover to argue for the sum of polynomial over
    /// {0,1}^`num_vars` is zero.
    fn prove(
        poly: &Self::PolyList,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, PolyIOPErrors> {
        let start = start_timer!(|| "zero check prove");

        let length = poly.domain_info.num_variables;
        let r = transcript.get_and_append_challenge_vectors(b"vector r", length)?;
        let f_hat = build_f_hat(poly, r.as_ref())?;
        let res = <Self as SumCheck<F>>::prove(&f_hat, transcript);

        end_timer!(start);
        res
    }

    /// Verify the claimed sum using the proof.
    /// the initial challenge `r` is also returned.
    /// The caller needs to makes sure that `\hat f = f * eq(x, r)`
    fn verify(
        proof: &Self::Proof,
        fx_domain_info: &Self::DomainInfo,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "zero check verify");

        // check that the sum is zero
        if proof.proofs[0].evaluations[0] + proof.proofs[0].evaluations[1] != F::zero() {
            return Err(PolyIOPErrors::InvalidProof(format!(
                "zero check: sum {} is not zero",
                proof.proofs[0].evaluations[0] + proof.proofs[0].evaluations[1]
            )));
        }

        // check the correctness of r (To be completed)
        let length = fx_domain_info.num_variables;
        let r = transcript.get_and_append_challenge_vectors(b"vector r", length)?;

        // hat_fx's max degree is increased by eq(x, r).degree() which is 1
        let mut hat_fx_domain_info = fx_domain_info.clone();
        hat_fx_domain_info.max_degree += 1;
        let subclaim =
            <Self as SumCheck<F>>::verify(F::zero(), proof, &hat_fx_domain_info, transcript)?;

        end_timer!(start);
        Ok((subclaim, r))
    }
}

// Input poly f(x) and a random vector r, output
//      \hat f(x) = \sum_{x_i \in eval_x} f(x_i) eq(x, r)
// where
//      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
fn build_f_hat<F: PrimeField>(
    poly: &VirtualPolynomial<F>,
    r: &[F],
) -> Result<VirtualPolynomial<F>, PolyIOPErrors> {
    let start = start_timer!(|| "zero check build hat f");

    assert_eq!(poly.domain_info.num_variables, r.len());

    let eq_x_r = build_eq_x_r(r);
    let mut res = poly.clone();
    res.mul_by_mle(eq_x_r, F::one())?;

    end_timer!(start);
    Ok(res)
}

// Evaluate
//      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
// over r, which is
//      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
fn build_eq_x_r<F: PrimeField>(r: &[F]) -> Rc<DenseMultilinearExtension<F>> {
    let start = start_timer!(|| "zero check build build eq_x_r");

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
            current_eval *= if bit { *ri } else { *one_minus_ri };
        }
        eval.push(current_eval);
    }
    let mle = DenseMultilinearExtension::from_evaluations_vec(num_var, eval);

    // println!("eq(x,r): {:?}, {}", mle, mle.evaluate(r).unwrap());
    let res = Rc::new(mle);
    end_timer!(start);
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
    use crate::{errors::PolyIOPErrors, PolyIOP, VirtualPolynomial};
    use ark_bls12_381::Fr;
    use ark_std::test_rng;

    fn test_zerocheck(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
    ) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        {
            // good path: zero virtual poly
            let poly =
                VirtualPolynomial::rand_zero(nv, num_multiplicands_range, num_products, &mut rng)?;

            let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let proof = <PolyIOP<Fr> as ZeroCheck<Fr>>::prove(&poly, &mut transcript)?;

            let poly_info = poly.domain_info.clone();
            let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let subclaim =
                <PolyIOP<Fr> as ZeroCheck<Fr>>::verify(&proof, &poly_info, &mut transcript)?.0;
            assert!(
                poly.evaluate(&subclaim.point)? == subclaim.expected_evaluation,
                "wrong subclaim"
            );
        }

        {
            // bad path: random virtual poly
            let (poly, _sum) =
                VirtualPolynomial::rand(nv, num_multiplicands_range, num_products, &mut rng)?;

            let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let proof = <PolyIOP<Fr> as ZeroCheck<Fr>>::prove(&poly, &mut transcript)?;

            let poly_info = poly.domain_info.clone();
            let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;

            assert!(
                <PolyIOP<Fr> as ZeroCheck<Fr>>::verify(&proof, &poly_info, &mut transcript)
                    .is_err()
            );
        }

        Ok(())
    }

    #[test]
    fn test_trivial_polynomial() -> Result<(), PolyIOPErrors> {
        let nv = 1;
        let num_multiplicands_range = (4, 5);
        let num_products = 1;

        test_zerocheck(nv, num_multiplicands_range, num_products)
    }
    #[test]
    fn test_normal_polynomial() -> Result<(), PolyIOPErrors> {
        let nv = 5;
        let num_multiplicands_range = (4, 9);
        let num_products = 5;

        test_zerocheck(nv, num_multiplicands_range, num_products)
    }
    #[test]

    fn zero_polynomial_should_error() -> Result<(), PolyIOPErrors> {
        let nv = 0;
        let num_multiplicands_range = (4, 13);
        let num_products = 5;

        assert!(test_zerocheck(nv, num_multiplicands_range, num_products).is_err());
        Ok(())
    }
}
