//! Prover subroutines for a SumCheck protocol.

use super::SumCheckProver;
use crate::poly_iop::{
    errors::PolyIOPErrors,
    structs::{IOPProverMessage, IOPProverState},
};
use arithmetic::{fix_first_variable, VirtualPolynomial};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer, vec::Vec};
use rayon::prelude::IntoParallelIterator;
use std::rc::Rc;

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

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
        // let start =
        //     start_timer!(|| format!("sum check prove {}-th round and update state",
        // self.round));

        if self.round >= self.poly.aux_info.num_variables {
            return Err(PolyIOPErrors::InvalidProver(
                "Prover is not active".to_string(),
            ));
        }

        // let fix_argument = start_timer!(|| "fix argument");

        // Step 1:
        // fix argument and evaluate f(x) over x_m = r; where r is the challenge
        // for the current round, and m is the round number, indexed from 1
        //
        // i.e.:
        // at round m <= n, for each mle g(x_1, ... x_n) within the flattened_mle
        // which has already been evaluated to
        //
        //    g(r_1, ..., r_{m-1}, x_m ... x_n)
        //
        // eval g over r_m, and mutate g to g(r_1, ... r_m,, x_{m+1}... x_n)
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

            let r = self.challenges[self.round - 1];
            #[cfg(feature = "parallel")]
            flattened_ml_extensions
                .par_iter_mut()
                .for_each(|mle| *mle = fix_first_variable(mle, &r));
            #[cfg(not(feature = "parallel"))]
            flattened_ml_extensions
                .iter_mut()
                .for_each(|mle| *mle = fix_variables(mle, &[r]));
        } else if self.round > 0 {
            return Err(PolyIOPErrors::InvalidProver(
                "verifier message is empty".to_string(),
            ));
        }
        // end_timer!(fix_argument);

        self.round += 1;

        let products_list = self.poly.products.clone();
        let mut products_sum = vec![F::zero(); self.poly.aux_info.max_degree + 1];

        // let compute_sum = start_timer!(|| "compute sum");

        // Step 2: generate sum for the partial evaluated polynomial:
        // f(r_1, ... r_m,, x_{m+1}... x_n)

        #[cfg(feature = "parallel")]
        let flag = (self.poly.aux_info.max_degree == 2)
            && (products_list.len() == 1)
            && (products_list[0].0 == F::one());
        if flag {
            for (t, e) in products_sum.iter_mut().enumerate() {
                let evals = (0..1 << (self.poly.aux_info.num_variables - self.round))
                    .into_par_iter()
                    .map(|b| {
                        // evaluate P_round(t)
                        let table0 = &flattened_ml_extensions[products_list[0].1[0]];
                        let table1 = &flattened_ml_extensions[products_list[0].1[1]];
                        let val = if t == 0 {
                            table0[b << 1] * table1[b << 1]
                        } else if t == 1 {
                            table0[(b << 1) + 1] * table1[(b << 1) + 1]
                        } else {
                            (table0[(b << 1) + 1] + table0[(b << 1) + 1] - table0[b << 1])
                                * (table1[(b << 1) + 1] + table1[(b << 1) + 1] - table1[b << 1])
                        };
                        val
                    })
                    .collect::<Vec<F>>();
                for val in evals.iter() {
                    *e += val
                }
            }
        } else {
            for (t, e) in products_sum.iter_mut().enumerate() {
                let t = F::from(t as u64);
                let one_minus_t = F::one() - t;
                let products = (0..1 << (self.poly.aux_info.num_variables - self.round))
                    .into_par_iter()
                    .map(|b| {
                        // evaluate P_round(t)
                        let mut tmp = F::zero();
                        products_list.iter().for_each(|(coefficient, products)| {
                            let num_mles = products.len();
                            let mut product = *coefficient;
                            for &f in products.iter().take(num_mles) {
                                let table = &flattened_ml_extensions[f]; // f's range is checked in init
                                product *= table[b << 1] * one_minus_t + table[(b << 1) + 1] * t;
                            }
                            tmp += product;
                        });

                        tmp
                    })
                    .collect::<Vec<F>>();

                for i in products.iter() {
                    *e += i
                }
            }
        }

        #[cfg(not(feature = "parallel"))]
        products_sum.iter_mut().enumerate().for_each(|(t, e)| {
            let t = F::from(t as u64);
            let one_minus_t = F::one() - t;

            for b in 0..1 << (self.poly.aux_info.num_variables - self.round) {
                // evaluate P_round(t)
                for (coefficient, products) in products_list.iter() {
                    let num_mles = products.len();
                    let mut product = *coefficient;
                    for &f in products.iter().take(num_mles) {
                        let table = &flattened_ml_extensions[f]; // f's range is checked in init
                        product *= table[b << 1] * one_minus_t + table[(b << 1) + 1] * t;
                    }
                    *e += product;
                }
            }
        });

        // update prover's state to the partial evaluated polynomial
        self.poly.flattened_ml_extensions = flattened_ml_extensions
            .iter()
            .map(|x| Rc::new(x.clone()))
            .collect();

        // end_timer!(compute_sum);
        // end_timer!(start);
        Ok(IOPProverMessage {
            evaluations: products_sum,
        })
    }
}
