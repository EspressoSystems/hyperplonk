//! Prover subroutines for a SumCheck protocol.

use super::SumCheckProver;
use crate::{
    errors::PolyIOPErrors,
    structs::{IOPProverMessage, IOPProverState},
    virtual_poly::VirtualPolynomial,
};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, start_timer, vec::Vec};
use std::rc::Rc;

#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

impl<F: PrimeField> SumCheckProver<F> for IOPProverState<F> {
    type VirtualPolynomial = VirtualPolynomial<F>;
    type ProverMessage = IOPProverMessage<F>;

    /// Initialize the prover state to argue for the sum of the input polynomial
    /// over {0,1}^`num_vars`.
    fn prover_init(polynomial: &Self::VirtualPolynomial) -> Result<Self, PolyIOPErrors> {
        let start = start_timer!(|| "sum check prover init");
        if polynomial.aux_info.num_variables == 0 {
            return Err(PolyIOPErrors::InvalidParameters(
                "Attempt to prove a constant.".to_string(),
            ));
        }
        end_timer!(start);

        Ok(Self {
            challenges: Vec::with_capacity(polynomial.aux_info.num_variables),
            round: 0,
            poly: polynomial.clone(),
        })
    }

    /// Receive message from verifier, generate prover message, and proceed to
    /// next round.
    ///
    /// Main algorithm used is from section 3.2 of [XZZPS19](https://eprint.iacr.org/2019/317.pdf#subsection.3.2).
    fn prove_round_and_update_state(
        &mut self,
        challenge: &Option<F>,
    ) -> Result<Self::ProverMessage, PolyIOPErrors> {
        let start =
            start_timer!(|| format!("sum check prove {}-th round and update state", self.round));

        let fix_argument = start_timer!(|| "fix argument");

        let mut flattened_ml_extensions: Vec<DenseMultilinearExtension<F>> = self
            .poly
            .flattened_ml_extensions
            .iter()
            .map(|x| x.as_ref().clone())
            .collect();

        if let Some(chal) = challenge {
            if self.round == 0 {
                return Err(PolyIOPErrors::InvalidProver(
                    "first round should be prover first.".to_string(),
                ));
            }
            self.challenges.push(*chal);

            // fix argument and evaluate f(x) over r for r.len() <= x.len()
            // i.e.:
            // for each mle g(x_0, ... x_n) within the flattened_mle
            // eval g(x_0,...x_n) over [r_0, ..., r_m] with m <= n
            // and mutate g to g(r_0, ...r_m, x_{m+1},... x_n)
            //
            // TODO(ZZ): have a question on this. Let's chat @binyi.
            let r = self.challenges[self.round - 1];
            #[cfg(feature = "parallel")]
            flattened_ml_extensions
                .par_iter_mut()
                .for_each(|mle| *mle = mle.fix_variables(&[r]));

            #[cfg(not(feature = "parallel"))]
            flattened_ml_extensions
                .iter_mut()
                .for_each(|mle| *mle = mle.fix_variables(&[r]));
        } else if self.round > 0 {
            return Err(PolyIOPErrors::InvalidProver(
                "verifier message is empty".to_string(),
            ));
        }
        end_timer!(fix_argument);

        self.round += 1;

        if self.round > self.poly.aux_info.num_variables {
            return Err(PolyIOPErrors::InvalidProver(
                "Prover is not active".to_string(),
            ));
        }

        let products_list = self.poly.products.clone();
        let i = self.round;
        let nv = self.poly.aux_info.num_variables;
        let degree = self.poly.aux_info.max_degree; // the degree of univariate polynomial sent by prover at this round

        let mut products_sum = Vec::with_capacity(degree + 1);
        products_sum.resize(degree + 1, F::zero());

        let compute_sum = start_timer!(|| "compute sum");
        // generate sum
        #[cfg(feature = "parallel")]
        products_sum.par_iter_mut().enumerate().for_each(|(t, e)| {
            for b in 0..1 << (nv - i) {
                // evaluate P_round(t)
                for (coefficient, products) in products_list.iter() {
                    let num_multiplicands = products.len();
                    let mut product = *coefficient;
                    for &f in products.iter().take(num_multiplicands) {
                        let table = &flattened_ml_extensions[f]; // f's range is checked in init
                        product *= table[b << 1] * (F::one() - F::from(t as u64))
                            + table[(b << 1) + 1] * F::from(t as u64);
                    }
                    *e += product;
                }
            }
        });

        #[cfg(not(feature = "parallel"))]
        for b in 0..1 << (nv - i) {
            products_sum
                .iter_mut()
                .take(degree + 1)
                .enumerate()
                .for_each(|(t, e)| {
                    // evaluate P_round(t)
                    for (coefficient, products) in products_list.iter() {
                        let num_multiplicands = products.len();
                        let mut product = *coefficient;
                        for &f in products.iter().take(num_multiplicands) {
                            let table = &flattened_ml_extensions[f]; // f's range is checked in init
                            product *= table[b << 1] * (F::one() - F::from(t as u64))
                                + table[(b << 1) + 1] * F::from(t as u64);
                        }
                        *e += product;
                    }
                });
        }

        self.poly.flattened_ml_extensions = flattened_ml_extensions
            .iter()
            .map(|x| Rc::new(x.clone()))
            .collect();

        end_timer!(compute_sum);
        end_timer!(start);
        Ok(IOPProverMessage {
            evaluations: products_sum,
        })
    }
}
