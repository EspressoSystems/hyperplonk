// todo: remove
#![allow(dead_code)]

use ark_ff::PrimeField;

use crate::{custom_gate::CustomizedGates, structs::HyperPlonkIndex, witness::WitnessColumn};

pub(crate) struct MockCircuit<F: PrimeField> {
    witnesses: Vec<WitnessColumn<F>>,
    index: HyperPlonkIndex<F>,
    num_constraints: usize,
}

impl<F: PrimeField> MockCircuit<F> {
    /// Generate a mock plonk circuit for the input constraint size.
    pub(crate) fn mock_circuit(_num_constraints: usize, _gate: &CustomizedGates) -> MockCircuit<F> {
        todo!()
    }
}
