// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

use crate::{build_mle, errors::HyperPlonkErrors};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::log2;
use std::sync::Arc;

/// A row of witnesses of width `#wires`
#[derive(Debug, Clone)]
pub struct WitnessRow<F: PrimeField>(pub(crate) Vec<F>);

/// A column of witnesses of length `#constraints`
#[derive(Debug, Clone, Default)]
pub struct WitnessColumn<F: PrimeField>(pub(crate) Vec<F>);

impl<F: PrimeField> WitnessColumn<F> {
    /// the number of variables of the multilinear polynomial that presents a
    /// column.
    pub fn get_nv(&self) -> usize {
        log2(self.0.len()) as usize
    }

    /// Append a new element to the witness column
    pub fn append(&mut self, new_element: F) {
        self.0.push(new_element)
    }

    /// Build witness columns from rows
    pub fn from_witness_rows(
        witness_rows: &[WitnessRow<F>],
    ) -> Result<Vec<Self>, HyperPlonkErrors> {
        if witness_rows.is_empty() {
            return Err(HyperPlonkErrors::InvalidParameters(
                "empty witness rows".to_string(),
            ));
        }

        let mut res = Vec::with_capacity(witness_rows.len());
        let num_columns = witness_rows[0].0.len();

        for i in 0..num_columns {
            let mut cur_column = Vec::new();
            for row in witness_rows.iter() {
                cur_column.push(row.0[i])
            }
            res.push(Self(cur_column))
        }

        Ok(res)
    }

    pub fn coeff_ref(&self) -> &[F] {
        self.0.as_ref()
    }
}

impl<F: PrimeField> From<&WitnessColumn<F>> for DenseMultilinearExtension<F> {
    fn from(witness: &WitnessColumn<F>) -> Self {
        let nv = witness.get_nv();
        Self::from_evaluations_slice(nv, witness.0.as_ref())
    }
}

impl<F: PrimeField> WitnessRow<F> {
    /// Build MLE from matrix of witnesses.
    ///
    /// Given a matrix := [row1, row2, ...] where
    /// row1:= (a1, a2, ...)
    /// row2:= (b1, b2, ...)
    /// row3:= (c1, c2, ...)
    ///
    /// output mle(a1,b1,c1, ...), mle(a2,b2,c2, ...), ...
    pub fn build_mles(
        matrix: &[Self],
    ) -> Result<Vec<Arc<DenseMultilinearExtension<F>>>, HyperPlonkErrors> {
        build_mle!(matrix)
    }
}
