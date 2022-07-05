use crate::{build_mle, errors::HyperPlonkErrors};
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::log2;
use std::rc::Rc;

#[cfg(test)]
use ark_std::rand::RngCore;

/// A row of selectors of width `SelectorWidth`
#[derive(Debug, Clone)]
pub struct SelectorRow<F: PrimeField>(Vec<F>);

impl<F: PrimeField> SelectorRow<F> {
    /// Build MLE from rows of witnesses.
    pub fn build_mles(
        rows: &[Self],
    ) -> Result<Vec<Rc<DenseMultilinearExtension<F>>>, HyperPlonkErrors> {
        build_mle!(rows)
    }
}

#[cfg(test)]
impl<F: PrimeField> SelectorRow<F> {
    /// sample a row of random selector
    pub(crate) fn rand<R: RngCore>(rng: &mut R, num_vars: usize) -> Self {
        Self((0..1 << num_vars).map(|_| F::rand(rng)).collect())
    }
    /// sample a set of random selectors
    pub(crate) fn rand_selectors<R: RngCore>(
        rng: &mut R,
        num_vars: usize,
        num_selectors: usize,
    ) -> Vec<Self> {
        (0..num_selectors)
            .map(|_| Self((0..1 << num_vars).map(|_| F::rand(rng)).collect::<Vec<F>>()))
            .collect()
    }
}
