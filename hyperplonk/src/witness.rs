use crate::{build_mle, errors::HyperPlonkErrors};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::log2;
use std::rc::Rc;

/// A row of witnesses of width `SelectorWidth`
#[derive(Debug, Clone)]
pub struct WitnessRow<F: PrimeField>(pub(crate) Vec<F>);

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
    ) -> Result<Vec<Rc<DenseMultilinearExtension<F>>>, HyperPlonkErrors> {
        build_mle!(matrix)
    }
}
