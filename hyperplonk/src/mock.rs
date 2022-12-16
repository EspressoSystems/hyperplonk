use arithmetic::identity_permutation;
use ark_ff::PrimeField;
use ark_std::{log2, test_rng};

use crate::{
    custom_gate::CustomizedGates,
    selectors::SelectorColumn,
    structs::{HyperPlonkIndex, HyperPlonkParams},
    witness::WitnessColumn,
};

pub struct MockCircuit<F: PrimeField> {
    pub public_inputs: Vec<F>,
    pub witnesses: Vec<WitnessColumn<F>>,
    pub index: HyperPlonkIndex<F>,
}

impl<F: PrimeField> MockCircuit<F> {
    /// Number of variables in a multilinear system
    pub fn num_variables(&self) -> usize {
        self.index.num_variables()
    }

    /// number of selector columns
    pub fn num_selector_columns(&self) -> usize {
        self.index.num_selector_columns()
    }

    /// number of witness columns
    pub fn num_witness_columns(&self) -> usize {
        self.index.num_witness_columns()
    }
}

impl<F: PrimeField> MockCircuit<F> {
    /// Generate a mock plonk circuit for the input constraint size.
    pub fn new(num_constraints: usize, gate: &CustomizedGates) -> MockCircuit<F> {
        let mut rng = test_rng();
        let nv = log2(num_constraints);
        let num_selectors = gate.num_selector_columns();
        let num_witnesses = gate.num_witness_columns();
        let log_n_wires = log2(num_witnesses);
        let merged_nv = nv + log_n_wires;

        let mut selectors: Vec<SelectorColumn<F>> = vec![SelectorColumn::default(); num_selectors];
        let mut witnesses: Vec<WitnessColumn<F>> = vec![WitnessColumn::default(); num_witnesses];

        for _cs_counter in 0..num_constraints {
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
        let pub_input_len = ark_std::cmp::min(4, num_constraints);
        let public_inputs = witnesses[0].0[0..pub_input_len].to_vec();

        let params = HyperPlonkParams {
            num_constraints,
            num_pub_input: public_inputs.len(),
            gate_func: gate.clone(),
        };

        let permutation = identity_permutation(merged_nv as usize, 1);
        let index = HyperPlonkIndex {
            params,
            permutation,
            selectors,
        };

        Self {
            public_inputs,
            witnesses,
            index,
        }
    }

    pub fn is_satisfied(&self) -> bool {
        for current_row in 0..self.num_variables() {
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
    use ark_bw6_761::{Fr as BwFr, BW6_761};
    use subroutines::{
        pcs::{
            prelude::{MultilinearKzgPCS, MultilinearUniversalParams},
            PolynomialCommitmentScheme,
        },
        poly_iop::PolyIOP,
    };

    const SUPPORTED_SIZE: usize = 20;
    const MIN_NUM_VARS: usize = 8;
    const MAX_NUM_VARS: usize = 15;
    const CUSTOM_DEGREE: [usize; 6] = [1, 2, 4, 8, 16, 32];

    fn test_mock_circuit_zkp_helper_bw_curve(
        nv: usize,
        gate: &CustomizedGates,
        pcs_srs: &MultilinearUniversalParams<BW6_761>,
    ) -> Result<(), HyperPlonkErrors> {
        let circuit = MockCircuit::<BwFr>::new(1 << nv, gate);
        assert!(circuit.is_satisfied());

        let index = circuit.index;

        // generate pk and vks
        let (pk, _vk) =
            <PolyIOP<BwFr> as HyperPlonkSNARK<BW6_761, MultilinearKzgPCS<BW6_761>>>::preprocess(
                &index, &pcs_srs,
            )?;
        // generate a proof and verify
        let _proof =
            <PolyIOP<BwFr> as HyperPlonkSNARK<BW6_761, MultilinearKzgPCS<BW6_761>>>::prove(
                &pk,
                &circuit.witnesses[0].0[..],
                &circuit.witnesses[..],
            )?;

        // let verify =
        //     <PolyIOP<BwFr> as HyperPlonkSNARK<BW6_761,
        // MultilinearKzgPCS<BW6_761>>>::verify(         &vk,
        //         &circuit.witnesses[0].0[..],
        //         &proof,
        //     )?;
        // assert!(verify);
        Ok(())
    }

    #[test]
    fn test_mock_circuit_e2e_bw_curve() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs = MultilinearKzgPCS::<BW6_761>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;

        let turboplonk_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        test_mock_circuit_zkp_helper_bw_curve(17, &turboplonk_gate, &pcs_srs)?;
        test_mock_circuit_zkp_helper_bw_curve(20, &turboplonk_gate, &pcs_srs)?;

        Ok(())
    }

    #[test]
    fn test_mock_circuit_sat() {
        for i in 1..10 {
            let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
            let circuit = MockCircuit::<Fr>::new(1 << i, &vanilla_gate);
            assert!(circuit.is_satisfied());

            let jf_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
            let circuit = MockCircuit::<Fr>::new(1 << i, &jf_gate);
            assert!(circuit.is_satisfied());

            for num_witness in 2..10 {
                for degree in CUSTOM_DEGREE {
                    let mock_gate = CustomizedGates::mock_gate(num_witness, degree);
                    let circuit = MockCircuit::<Fr>::new(1 << i, &mock_gate);
                    assert!(circuit.is_satisfied());
                }
            }
        }
    }

    fn test_mock_circuit_zkp_helper(
        nv: usize,
        gate: &CustomizedGates,
        pcs_srs: &MultilinearUniversalParams<Bls12_381>,
    ) -> Result<(), HyperPlonkErrors> {
        let circuit = MockCircuit::<Fr>::new(1 << nv, gate);
        assert!(circuit.is_satisfied());

        let index = circuit.index;
        // generate pk and vks
        let (pk, vk) =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::preprocess(
                &index, &pcs_srs,
            )?;
        // generate a proof and verify
        let proof =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::prove(
                &pk,
                &circuit.public_inputs,
                &circuit.witnesses,
            )?;

        let verify =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::verify(
                &vk,
                &circuit.public_inputs,
                &proof,
            )?;
        assert!(verify);
        Ok(())
    }

    #[test]
    fn test_mock_circuit_zkp() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs =
            MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
        for nv in MIN_NUM_VARS..MAX_NUM_VARS {
            let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
            test_mock_circuit_zkp_helper(nv, &vanilla_gate, &pcs_srs)?;
        }
        for nv in MIN_NUM_VARS..MAX_NUM_VARS {
            let tubro_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
            test_mock_circuit_zkp_helper(nv, &tubro_gate, &pcs_srs)?;
        }
        let nv = 5;
        for num_witness in 2..10 {
            for degree in CUSTOM_DEGREE {
                let mock_gate = CustomizedGates::mock_gate(num_witness, degree);
                test_mock_circuit_zkp_helper(nv, &mock_gate, &pcs_srs)?;
            }
        }

        Ok(())
    }

    #[test]
    fn test_mock_circuit_e2e() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs =
            MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
        let nv = MAX_NUM_VARS;

        let turboplonk_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        test_mock_circuit_zkp_helper(nv, &turboplonk_gate, &pcs_srs)?;

        Ok(())
    }

    #[test]
    fn test_mock_long_selector_e2e() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs =
            MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
        let nv = MAX_NUM_VARS;

        let long_selector_gate = CustomizedGates::super_long_selector_gate();
        test_mock_circuit_zkp_helper(nv, &long_selector_gate, &pcs_srs)?;

        Ok(())
    }
}
