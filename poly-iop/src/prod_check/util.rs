//! This module implements useful functions for the product check protocol.

use crate::{
    errors::PolyIOPErrors, structs::IOPProof, utils::get_index, zero_check::ZeroCheck, PolyIOP,
};
use arithmetic::VirtualPolynomial;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use std::rc::Rc;
use transcript::IOPTranscript;

/// Compute the product polynomial `prod(x)` where
///
///  - `prod(0,x) := prod(0, x1, …, xn)` is the MLE over the
/// evaluations of `f(x)/g(x)` on the boolean hypercube {0,1}^n
///
/// - `prod(1,x)` is a MLE over the evaluations of `prod(x, 0) * prod(x, 1)`
/// on the boolean hypercube {0,1}^n
///
/// The caller needs to check num_vars matches in f and g
/// Cost: linear in N.
pub(super) fn compute_product_poly<F: PrimeField>(
    fx: &Rc<DenseMultilinearExtension<F>>,
    gx: &Rc<DenseMultilinearExtension<F>>,
) -> Result<Rc<DenseMultilinearExtension<F>>, PolyIOPErrors> {
    let start = start_timer!(|| "compute evaluations of prod polynomial");
    let num_vars = fx.num_vars;

    // ===================================
    // prod(0, x)
    // ===================================
    let prod_0x_eval = compute_prod_0(fx, gx)?;

    // ===================================
    // prod(1, x)
    // ===================================
    //
    // `prod(1, x)` can be computed via recursing the following formula for 2^n-1
    // times
    //
    // `prod(1, x_1, ..., x_n) :=
    //      prod(x_1, x_2, ..., x_n, 0) * prod(x_1, x_2, ..., x_n, 1)`
    //
    // At any given step, the right hand side of the equation
    // is available via either eval_0x or the current view of eval_1x
    let mut prod_1x_eval = vec![];
    for x in 0..(1 << num_vars) - 1 {
        // sign will decide if the evaluation should be looked up from eval_0x or
        // eval_1x; x_zero_index is the index for the evaluation (x_2, ..., x_n,
        // 0); x_one_index is the index for the evaluation (x_2, ..., x_n, 1);
        let (x_zero_index, x_one_index, sign) = get_index(x, num_vars);
        if !sign {
            prod_1x_eval.push(prod_0x_eval[x_zero_index] * prod_0x_eval[x_one_index]);
        } else {
            // sanity check: if we are trying to look up from the eval_1x table,
            // then the target index must already exist
            if x_zero_index >= prod_1x_eval.len() || x_one_index >= prod_1x_eval.len() {
                return Err(PolyIOPErrors::ShouldNotArrive);
            }
            prod_1x_eval.push(prod_1x_eval[x_zero_index] * prod_1x_eval[x_one_index]);
        }
    }

    // prod(1, 1, ..., 1) := 0
    prod_1x_eval.push(F::zero());

    // ===================================
    // prod(x)
    // ===================================
    // prod(x)'s evaluation is indeed `e := [eval_0x[..], eval_1x[..]].concat()`
    let eval = [prod_0x_eval.as_slice(), prod_1x_eval.as_slice()].concat();

    let prod_x = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars + 1,
        eval,
    ));

    end_timer!(start);
    Ok(prod_x)
}

/// generate the zerocheck proof for the virtual polynomial
///    prod(1, x) - prod(x, 0) * prod(x, 1) + alpha * (f(x) - prod(0, x) * g(x))
///
/// Returns proof and Q(x) for testing purpose.
///
/// Cost: O(N)
pub(super) fn prove_zero_check<F: PrimeField>(
    fx: &Rc<DenseMultilinearExtension<F>>,
    gx: &Rc<DenseMultilinearExtension<F>>,
    prod_x: &Rc<DenseMultilinearExtension<F>>,
    alpha: &F,
    transcript: &mut IOPTranscript<F>,
) -> Result<(IOPProof<F>, VirtualPolynomial<F>), PolyIOPErrors> {
    let start = start_timer!(|| "zerocheck in product check");

    let prod_partial_evals = build_prod_partial_eval(prod_x)?;
    let prod_0x = prod_partial_evals[0].clone();
    let prod_1x = prod_partial_evals[1].clone();
    let prod_x0 = prod_partial_evals[2].clone();
    let prod_x1 = prod_partial_evals[3].clone();

    // compute g(x) * prod(0, x) * alpha
    let mut q_x = VirtualPolynomial::new_from_mle(gx, F::one());
    q_x.mul_by_mle(prod_0x, *alpha)?;

    //   g(x) * prod(0, x) * alpha
    // - f(x) * alpha
    q_x.add_mle_list([fx.clone()], -*alpha)?;

    // Q(x) := prod(1,x) - prod(x, 0) * prod(x, 1)
    //       + alpha * (
    //             g(x) * prod(0, x)
    //           - f(x))
    q_x.add_mle_list([prod_x0, prod_x1], -F::one())?;
    q_x.add_mle_list([prod_1x], F::one())?;

    let iop_proof = <PolyIOP<F> as ZeroCheck<F>>::prove(&q_x, transcript)?;

    end_timer!(start);
    Ok((iop_proof, q_x))
}

/// Helper function of the IOP.
///
/// Input:
/// - prod(x)
///
/// Output: the following 4 polynomials
/// - prod(0, x)
/// - prod(1, x)
/// - prod(x, 0)
/// - prod(x, 1)
fn build_prod_partial_eval<F: PrimeField>(
    prod_x: &Rc<DenseMultilinearExtension<F>>,
) -> Result<[Rc<DenseMultilinearExtension<F>>; 4], PolyIOPErrors> {
    let start = start_timer!(|| "build partial prod polynomial");

    let prod_x_eval = &prod_x.evaluations;
    let num_vars = prod_x.num_vars - 1;

    // prod(0, x)
    let prod_0_x = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
        num_vars,
        &prod_x_eval[0..1 << num_vars],
    ));
    // prod(1, x)
    let prod_1_x = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
        num_vars,
        &prod_x_eval[1 << num_vars..1 << (num_vars + 1)],
    ));

    // ===================================
    // prod(x, 0) and prod(x, 1)
    // ===================================
    //
    // now we compute eval_x0 and eval_x1
    // eval_0x will be the odd coefficients of eval
    // and eval_1x will be the even coefficients of eval
    let mut eval_x0 = vec![];
    let mut eval_x1 = vec![];
    for (x, &prod_x) in prod_x_eval.iter().enumerate() {
        if x & 1 == 0 {
            eval_x0.push(prod_x);
        } else {
            eval_x1.push(prod_x);
        }
    }
    let prod_x_0 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars, eval_x0,
    ));
    let prod_x_1 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars, eval_x1,
    ));

    end_timer!(start);

    Ok([prod_0_x, prod_1_x, prod_x_0, prod_x_1])
}

/// Returns the evaluations of
/// - `prod(0,x) := prod(0, x1, …, xn)` which is the MLE over the
/// evaluations of f(x)/g(x) on the boolean hypercube {0,1}^n:
///
/// The caller needs to check num_vars matches in f/g
/// Cost: linear in N.
fn compute_prod_0<F: PrimeField>(
    fx: &DenseMultilinearExtension<F>,
    gx: &DenseMultilinearExtension<F>,
) -> Result<Vec<F>, PolyIOPErrors> {
    let start = start_timer!(|| "compute prod(0,x)");

    // let mut prod_0x_evals = vec![];
    // for (&fi, &gi) in fx.iter().zip(gx.iter()) {
    //     prod_0x_evals.push(fi / gi);
    // }
    let input = fx
        .iter()
        .zip(gx.iter())
        .map(|(&fi, &gi)| (fi, gi))
        .collect::<Vec<_>>();
    let prod_0x_evals = input.par_iter().map(|(x, y)| *x / *y).collect::<Vec<_>>();

    end_timer!(start);
    Ok(prod_0x_evals)
}
