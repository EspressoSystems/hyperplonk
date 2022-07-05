use crate::{build_mle, errors::HyperPlonkErrors};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::log2;
use std::rc::Rc;

/// A row of witnesses of width `SelectorWidth`
#[derive(Debug, Clone)]
pub struct WitnessRow<F: PrimeField>(pub(crate) Vec<F>);

impl<F: PrimeField> WitnessRow<F> {
    /// Build MLE from rows of witnesses.
    pub fn build_mles(
        rows: &[Self],
    ) -> Result<Vec<Rc<DenseMultilinearExtension<F>>>, HyperPlonkErrors> {
        build_mle!(rows)
    }
}
