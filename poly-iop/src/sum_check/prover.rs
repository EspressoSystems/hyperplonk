//! Prover
// TODO: some of the struct is generic for Sum Checks and Zero Checks.
// If so move them to src/structs.rs
use super::SumCheckProver;
use crate::{poly_list::PolynomialList, structs::IOPProverMessage, PolyIOP};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::vec::Vec;

/// Prover State
pub struct ProverState<F: PrimeField> {
    /// sampled randomness given by the verifier
    pub challenges: Vec<F>,
    /// Stores the list of products that is meant to be added together. Each
    /// multiplicand is represented by the index in flattened_ml_extensions
    pub list_of_products: Vec<(F, Vec<usize>)>,
    /// Stores a list of multilinear extensions in which `self.list_of_products`
    /// points to
    pub flattened_ml_extensions: Vec<DenseMultilinearExtension<F>>,
    pub(crate) num_vars: usize,
    pub(crate) max_multiplicands: usize,
    pub(crate) round: usize,
}

impl<F: PrimeField> SumCheckProver<F> for PolyIOP<F> {
    type PolyList = PolynomialList<F>;
    type ProverState = ProverState<F>;
    type ProverMessage = IOPProverMessage<F>;

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
    fn prover_init(polynomial: &Self::PolyList) -> Self::ProverState {
        if polynomial.aux_info.num_variables == 0 {
            panic!("Attempt to prove a constant.")
        }

        // create a deep copy of all unique MLExtensions
        let flattened_ml_extensions = polynomial
            .flattened_ml_extensions
            .iter()
            .map(|x| x.as_ref().clone())
            .collect();

        ProverState {
            challenges: Vec::with_capacity(polynomial.aux_info.num_variables),
            list_of_products: polynomial.products.clone(),
            flattened_ml_extensions,
            num_vars: polynomial.aux_info.num_variables,
            max_multiplicands: polynomial.aux_info.max_multiplicands,
            round: 0,
        }
    }

    /// receive message from verifier, generate prover message, and proceed to
    /// next round
    ///
    /// Main algorithm used is from section 3.2 of [XZZPS19](https://eprint.iacr.org/2019/317.pdf#subsection.3.2).
    fn prove_round_and_update_state(
        prover_state: &mut Self::ProverState,
        challenge: &Option<F>,
    ) -> Self::ProverMessage {
        if let Some(chal) = challenge {
            if prover_state.round == 0 {
                panic!("first round should be prover first.");
            }
            prover_state.challenges.push(*chal);

            // fix argument
            let i = prover_state.round;
            let r = prover_state.challenges[i - 1];
            for multiplicand in prover_state.flattened_ml_extensions.iter_mut() {
                *multiplicand = multiplicand.fix_variables(&[r]);
            }
        } else if prover_state.round > 0 {
            panic!("verifier message is empty");
        }

        prover_state.round += 1;

        if prover_state.round > prover_state.num_vars {
            panic!("Prover is not active");
        }

        let i = prover_state.round;
        let nv = prover_state.num_vars;
        let degree = prover_state.max_multiplicands; // the degree of univariate polynomial sent by prover at this round

        let mut products_sum = Vec::with_capacity(degree + 1);
        products_sum.resize(degree + 1, F::zero());

        // generate sum
        for b in 0..1 << (nv - i) {
            let mut t_as_field = F::zero();
            for e in products_sum.iter_mut().take(degree + 1) {
                // evaluate P_round(t)
                for (coefficient, products) in &prover_state.list_of_products {
                    let num_multiplicands = products.len();
                    let mut product = *coefficient;
                    for &f in products.iter().take(num_multiplicands) {
                        let table = &prover_state.flattened_ml_extensions[f]; // j's range is checked in init
                        product *= table[b << 1] * (F::one() - t_as_field)
                            + table[(b << 1) + 1] * t_as_field;
                    }
                    *e += product;
                }
                t_as_field += F::one();
            }
        }

        IOPProverMessage {
            evaluations: products_sum,
        }
    }
}
