use crate::errors::HyperPlonkErrors;
use arithmetic::prelude::DenseMultilinearExtension;
use ark_ff::PrimeField;
use ark_std::log2;

/// A column of selectors of length `#constraints`
#[derive(Debug, Clone)]
pub struct SelectorColumn<F: PrimeField>(pub(crate) Vec<F>);

impl<F: PrimeField> SelectorColumn<F> {
    /// the number of variables for MLE to present a column.
    pub fn get_nv(&self) -> usize {
        log2(self.0.len()) as usize
    }

    /// Build selector columns from rows
    pub fn from_selector_rows(
        selector_rows: &[SelectorColumn<F>],
    ) -> Result<Vec<Self>, HyperPlonkErrors> {
        if selector_rows.is_empty() {
            return Err(HyperPlonkErrors::InvalidParameters(
                "empty witness rows".to_string(),
            ));
        }

        let mut res = Vec::with_capacity(selector_rows.len());
        let num_wires = selector_rows[0].0.len();

        for i in 0..num_wires {
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
