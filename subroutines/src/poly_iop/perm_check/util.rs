//! This module implements useful functions for the permutation check protocol.

use crate::poly_iop::errors::PolyIOPErrors;
use arithmetic::identity_permutation_mle;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use std::rc::Rc;

/// Returns the evaluations of two MLEs:
/// - numerator
/// - denominator
///
///  where
///  - beta and gamma are challenges
///  - f(x), g(x), s_id(x), s_perm(x) are mle-s
///
/// - numerator is the MLE for `f(x) + \beta s_id(x) + \gamma`
/// - denominator is the MLE for `g(x) + \beta s_perm(x) + \gamma`
#[allow(clippy::type_complexity)]
pub(super) fn computer_num_and_denom<F: PrimeField>(
    beta: &F,
    gamma: &F,
    fx: &DenseMultilinearExtension<F>,
    gx: &DenseMultilinearExtension<F>,
    s_perm: &DenseMultilinearExtension<F>,
) -> Result<
    (
        Rc<DenseMultilinearExtension<F>>,
        Rc<DenseMultilinearExtension<F>>,
    ),
    PolyIOPErrors,
> {
    let start = start_timer!(|| "compute numerator and denominator");

    let num_vars = fx.num_vars;
    let mut numerator_evals = vec![];
    let mut denominator_evals = vec![];
    let s_id = identity_permutation_mle::<F>(num_vars);

    for (&fi, (&gi, (&s_id_i, &s_perm_i))) in
        fx.iter().zip(gx.iter().zip(s_id.iter().zip(s_perm.iter())))
    {
        let numerator = fi + *beta * s_id_i + gamma;
        let denominator = gi + *beta * s_perm_i + gamma;

        numerator_evals.push(numerator);
        denominator_evals.push(denominator);
    }
    let numerator = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars,
        numerator_evals,
    ));
    let denominator = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars,
        denominator_evals,
    ));

    end_timer!(start);
    Ok((numerator, denominator))
}
