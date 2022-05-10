mod prover;
mod verifier;

use ark_ff::PrimeField;
use ark_std::test_rng;
pub use prover::ProverState;
pub use verifier::VerifierState;

use crate::{
    errors::PolyIOPErrors,
    poly_list::PolynomialList,
    structs::{DomainInfo, IOPProof, SubClaim},
    sum_check::SumCheck,
    PolyIOP,
};

pub trait ZeroCheck<F: PrimeField> {
    type Proof;
    type PolyList;
    type DomainInfo;
    type SubClaim;

    fn prove(poly: &Self::PolyList) -> Result<Self::Proof, PolyIOPErrors>;

    /// verify the claimed sum using the proof
    fn verify(
        hat_f: F,
        proof: &Self::Proof,
        domain_info: &Self::DomainInfo,
    ) -> Result<Self::SubClaim, PolyIOPErrors>;
}

impl<F: PrimeField> ZeroCheck<F> for PolyIOP<F> {
    type Proof = IOPProof<F>;

    type PolyList = PolynomialList<F>;

    type DomainInfo = DomainInfo<F>;

    type SubClaim = SubClaim<F>;

    fn prove(poly: &Self::PolyList) -> Result<Self::Proof, PolyIOPErrors> {
        // TODO: sample eval_x from Transcript
        let mut rng = test_rng();
        let length = poly.domain_info.num_variables;
        let r: Vec<F> = (0..length).map(|_| F::rand(&mut rng)).collect();

        let g_r = build_g_r(&poly, r.as_ref());

        <Self as SumCheck<F>>::prove(&g_r)
    }

    /// verify the claimed sum using the proof
    fn verify(
        hat_f: F,
        proof: &Self::Proof,
        domain_info: &Self::DomainInfo,
    ) -> Result<Self::SubClaim, PolyIOPErrors> {
        <Self as SumCheck<F>>::verify(hat_f, proof, domain_info)
    }
}

// Input poly f(x) and a random vector r, output
//      g(r) = \sum_{x_i \in eval_x} f(x_i) eq(x, r)
// where
//      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
fn build_g_r<F: PrimeField>(poly: &PolynomialList<F>, r: &[F]) -> PolynomialList<F> {
    assert_eq!(poly.domain_info.num_variables, r.len());
    let num_var = r.len();

    let res = PolynomialList::new(num_var);
    for _ in 0..num_var {}

    res
}
