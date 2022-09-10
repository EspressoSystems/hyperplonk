// todo: remove
#![allow(dead_code)]

use ark_ff::PrimeField;
use ark_poly::MultilinearExtension;
use ark_std::{log2, test_rng};
use poly_iop::prelude::identity_permutation_mle;

use crate::{
    custom_gate::CustomizedGates,
    selectors::SelectorColumn,
    structs::{HyperPlonkIndex, HyperPlonkParams},
    witness::WitnessColumn,
};

pub(crate) struct MockCircuit<F: PrimeField> {
    witnesses: Vec<WitnessColumn<F>>,
    index: HyperPlonkIndex<F>,
}

impl<F: PrimeField> MockCircuit<F> {
    /// Generate a mock plonk circuit for the input constraint size.
    pub(crate) fn mock_circuit(num_constraints: usize, gate: &CustomizedGates) -> MockCircuit<F> {
        let mut rng = test_rng();
        let nv = log2(num_constraints);
        let num_selectors = gate.num_selectors();
        let num_witnesses = gate.num_witnesses();
        let log_n_wires = log2(num_witnesses);
        let merged_nv = nv + log_n_wires;

        let mut selectors: Vec<SelectorColumn<F>> = vec![SelectorColumn::default(); num_selectors];
        let mut witnesses: Vec<WitnessColumn<F>> = vec![WitnessColumn::default(); num_witnesses];

        for _cs_counter in 0..1 << nv {
            let mut cur_selectors: Vec<F> = (0..(num_selectors - 1))
                .map(|_| F::rand(&mut rng))
                .collect();
            let cur_witness: Vec<F> = (0..num_witnesses).map(|_| F::rand(&mut rng)).collect();
            let mut last_selector = F::zero();
            for (index, (coeff, q, wit)) in gate.gates.iter().enumerate() {
                if index != num_selectors - 1 {
                    let mut cur_monomial = if *coeff < 0 {
                        -F::from((-coeff) as u64)
                    } else {
                        F::from(*coeff as u64)
                    };
                    cur_monomial = match q {
                        Some(p) => cur_monomial * cur_selectors[*p],
                        None => cur_monomial,
                    };
                    for wit_index in wit.iter() {
                        cur_monomial *= cur_witness[*wit_index];
                    }
                    last_selector += cur_monomial;
                } else {
                    let mut cur_monomial = if *coeff < 0 {
                        -F::from((-coeff) as u64)
                    } else {
                        F::from(*coeff as u64)
                    };
                    for wit_index in wit.iter() {
                        cur_monomial *= cur_witness[*wit_index];
                    }
                    last_selector /= -cur_monomial;
                }
            }
            cur_selectors.push(last_selector);
            for i in 0..num_selectors {
                selectors[i].append(cur_selectors[i]);
            }
            for i in 0..num_witnesses {
                witnesses[i].append(cur_witness[i]);
            }
        }

        let params = HyperPlonkParams {
            nv: nv as usize,
            log_pub_input_len: 0,
            log_n_selectors: log2(num_selectors) as usize,
            log_n_wires: log2(num_witnesses) as usize,
            gate_func: gate.clone(),
        };

        let permutation = identity_permutation_mle(merged_nv as usize).to_evaluations();
        let index = HyperPlonkIndex {
            params,
            permutation,
            selectors,
        };

        Self { witnesses, index }
    }

    pub(crate) fn is_satisfied(&self) -> bool {
        for current_row in 0..1 << self.index.params.nv {
            let mut cur = F::zero();
            for (coeff, q, wit) in self.index.params.gate_func.gates.iter() {
                let mut cur_monomial = if *coeff < 0 {
                    -F::from((-coeff) as u64)
                } else {
                    F::from(*coeff as u64)
                };
                cur_monomial = match q {
                    Some(p) => cur_monomial * self.index.selectors[*p].0[current_row],
                    None => cur_monomial,
                };
                for wit_index in wit.iter() {
                    cur_monomial *= self.witnesses[*wit_index].0[current_row];
                }
                cur += cur_monomial;
            }
            if !cur.is_zero() {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{errors::HyperPlonkErrors, HyperPlonkSNARK};
    use ark_bls12_381::{Bls12_381, Fr};
    use pcs::{prelude::MultilinearKzgPCS, PolynomialCommitmentScheme};
    use poly_iop::PolyIOP;

    #[test]
    fn test_mock_circuit_sat() {
        for i in 1..10 {
            let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
            let circuit = MockCircuit::<Fr>::mock_circuit(1 << i, &vanilla_gate);
            assert!(circuit.is_satisfied());

            let jf_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
            let circuit = MockCircuit::<Fr>::mock_circuit(1 << i, &jf_gate);
            assert!(circuit.is_satisfied());
        }
    }

    #[test]
    fn test_mock_circuit_zkp() -> Result<(), HyperPlonkErrors> {
        for nv in 1..10 {
            println!("nv: {}", nv);
            let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
            let circuit = MockCircuit::<Fr>::mock_circuit(1 << nv, &vanilla_gate);
            assert!(circuit.is_satisfied());

            let index = circuit.index;

            let mut rng = test_rng();
            let pcs_srs = MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, 16)?;

            // generate pk and vks
            let (pk, vk) = <PolyIOP<Fr> as HyperPlonkSNARK<
                Bls12_381,
                MultilinearKzgPCS<Bls12_381>,
            >>::preprocess(&index, &pcs_srs)?;

            // generate a proof and verify
            let proof = <PolyIOP<Fr> as HyperPlonkSNARK<
                Bls12_381,
                MultilinearKzgPCS<Bls12_381>,
            >>::prove(
                &pk, &[circuit.witnesses[0].0[0]], &circuit.witnesses
            )?;

            let _verify = <PolyIOP<Fr> as HyperPlonkSNARK<
                Bls12_381,
                MultilinearKzgPCS<Bls12_381>,
            >>::verify(&vk, &[circuit.witnesses[0].0[0]], &proof)?;
        }
        Ok(())
    }
}
