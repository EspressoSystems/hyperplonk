//! Prover
// TODO: some of the struct is generic for Sum Checks and Zero Checks.
// If so move them to src/structs.rs
use super::SumCheckProver;
use crate::{
    errors::PolyIOPErrors, structs::IOPProverMessage, virtual_poly::VirtualPolynomial, PolyIOP,
};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, start_timer, vec::Vec};

#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

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
    pub(crate) max_degree: usize,
    pub(crate) round: usize,
}

impl<F: PrimeField> SumCheckProver<F> for PolyIOP<F> {
    type PolyList = VirtualPolynomial<F>;
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
    fn prover_init(polynomial: &Self::PolyList) -> Result<Self::ProverState, PolyIOPErrors> {
        let start = start_timer!(|| "prover init");
        if polynomial.domain_info.num_variables == 0 {
            return Err(PolyIOPErrors::InvalidParameters(
                "Attempt to prove a constant.".to_string(),
            ));
        }

        // create a deep copy of all unique MLExtensions
        let flattened_ml_extensions = polynomial
            .flattened_ml_extensions
            .iter()
            .map(|x| x.as_ref().clone())
            .collect();

        end_timer!(start);

        Ok(ProverState {
            challenges: Vec::with_capacity(polynomial.domain_info.num_variables),
            list_of_products: polynomial.products.clone(),
            flattened_ml_extensions,
            num_vars: polynomial.domain_info.num_variables,
            max_degree: polynomial.domain_info.max_degree,
            round: 0,
        })
    }

    /// receive message from verifier, generate prover message, and proceed to
    /// next round
    ///
    /// Main algorithm used is from section 3.2 of [XZZPS19](https://eprint.iacr.org/2019/317.pdf#subsection.3.2).
    fn prove_round_and_update_state(
        prover_state: &mut Self::ProverState,
        challenge: &Option<F>,
    ) -> Result<Self::ProverMessage, PolyIOPErrors> {
        let start =
            start_timer!(|| format!("prove {}-th round and update state", prover_state.round));

        let fix_argument = start_timer!(|| "fix argument");
        if let Some(chal) = challenge {
            if prover_state.round == 0 {
                return Err(PolyIOPErrors::InvalidProver(
                    "first round should be prover first.".to_string(),
                ));
            }
            prover_state.challenges.push(*chal);

            // fix argument
            let i = prover_state.round;
            let r = prover_state.challenges[i - 1];
            #[cfg(feature = "parallel")]
            prover_state
                .flattened_ml_extensions
                .par_iter_mut()
                .for_each(|multiplicand| *multiplicand = multiplicand.fix_variables(&[r]));

            #[cfg(not(feature = "parallel"))]
            prover_state
                .flattened_ml_extensions
                .iter_mut()
                .for_each(|multiplicand| *multiplicand = multiplicand.fix_variables(&[r]));
        } else if prover_state.round > 0 {
            return Err(PolyIOPErrors::InvalidProver(
                "verifier message is empty".to_string(),
            ));
        }
        end_timer!(fix_argument);

        prover_state.round += 1;

        if prover_state.round > prover_state.num_vars {
            return Err(PolyIOPErrors::InvalidProver(
                "Prover is not active".to_string(),
            ));
        }

        let i = prover_state.round;
        let nv = prover_state.num_vars;
        let degree = prover_state.max_degree; // the degree of univariate polynomial sent by prover at this round

        let mut products_sum = Vec::with_capacity(degree + 1);
        products_sum.resize(degree + 1, F::zero());

        let compute_sum = start_timer!(|| "compute sum");
        // generate sum
        for b in 0..1 << (nv - i) {
            #[cfg(feature = "parallel")]
            products_sum
                .par_iter_mut()
                .take(degree + 1)
                .enumerate()
                .for_each(|(i, e)| {
                    // evaluate P_round(t)
                    for (coefficient, products) in &prover_state.list_of_products {
                        let num_multiplicands = products.len();
                        let mut product = *coefficient;
                        for &f in products.iter().take(num_multiplicands) {
                            let table = &prover_state.flattened_ml_extensions[f]; // j's range is checked in init
                            product *= table[b << 1] * (F::one() - F::from(i as u64))
                                + table[(b << 1) + 1] * F::from(i as u64);
                        }
                        *e += product;
                    }
                });
            #[cfg(not(feature = "parallel"))]
            products_sum
                .iter_mut()
                .take(degree + 1)
                .enumerate()
                .for_each(|(i, e)| {
                    // evaluate P_round(t)
                    for (coefficient, products) in &prover_state.list_of_products {
                        let num_multiplicands = products.len();
                        let mut product = *coefficient;
                        for &f in products.iter().take(num_multiplicands) {
                            let table = &prover_state.flattened_ml_extensions[f]; // j's range is checked in init
                            product *= table[b << 1] * (F::one() - F::from(i as u64))
                                + table[(b << 1) + 1] * F::from(i as u64);
                        }
                        *e += product;
                    }
                });
        }
        end_timer!(compute_sum);
        end_timer!(start);
        Ok(IOPProverMessage {
            evaluations: products_sum,
        })
    }
}
