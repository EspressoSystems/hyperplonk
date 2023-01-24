// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

use crate::{build_mle, errors::HyperPlonkErrors};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::log2;
use std::sync::Arc;

/// A row of selector of width `#selectors`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectorRow<F: PrimeField>(pub(crate) Vec<F>);

/// A column of selectors of length `#constraints`
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SelectorColumn<F: PrimeField>(pub(crate) Vec<F>);

impl<F: PrimeField> SelectorColumn<F> {
    /// the number of variables of the multilinear polynomial that presents a
    /// column.
    pub fn get_nv(&self) -> usize {
        log2(self.0.len()) as usize
    }

    /// Append a new element to the selector column
    pub fn append(&mut self, new_element: F) {
        self.0.push(new_element)
    }

    /// Build selector columns from rows
    pub fn from_selector_rows(
        selector_rows: &[SelectorRow<F>],
    ) -> Result<Vec<Self>, HyperPlonkErrors> {
        if selector_rows.is_empty() {
            return Err(HyperPlonkErrors::InvalidParameters(
                "empty witness rows".to_string(),
            ));
        }

        let mut res = Vec::with_capacity(selector_rows.len());
        let num_colnumns = selector_rows[0].0.len();

        for i in 0..num_colnumns {
            let mut cur_column = Vec::new();
            for row in selector_rows.iter() {
                cur_column.push(row.0[i])
            }
            res.push(Self(cur_column))
        }

        Ok(res)
    }
}

impl<F: PrimeField> From<&SelectorColumn<F>> for DenseMultilinearExtension<F> {
    fn from(witness: &SelectorColumn<F>) -> Self {
        let nv = witness.get_nv();
        Self::from_evaluations_slice(nv, witness.0.as_ref())
    }
}

impl<F: PrimeField> SelectorRow<F> {
    /// Build MLE from matrix of selectors.
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
