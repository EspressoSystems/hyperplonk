// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! This module implements useful functions for the product check protocol.

use crate::poly_iop::{errors::PolyIOPErrors, structs::IOPProof, zero_check::ZeroCheck, PolyIOP};
use arithmetic::{get_index, VirtualPolynomial};
use ark_ff::{batch_inversion, PrimeField};
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, start_timer};
use std::sync::Arc;
use transcript::IOPTranscript;

/// Compute multilinear fractional polynomial s.t. frac(x) = f1(x) * ... * fk(x)
/// / (g1(x) * ... * gk(x)) for all x \in {0,1}^n
///
/// The caller needs to sanity-check that the number of polynomials and
/// variables match in fxs and gxs; and gi(x) has no zero entries.
pub(super) fn compute_frac_poly<F: PrimeField>(
    fxs: &[Arc<DenseMultilinearExtension<F>>],
    gxs: &[Arc<DenseMultilinearExtension<F>>],
) -> Result<Arc<DenseMultilinearExtension<F>>, PolyIOPErrors> {
    let start = start_timer!(|| "compute frac(x)");

    let mut f_evals = vec![F::one(); 1 << fxs[0].num_vars];
    for fx in fxs.iter() {
        for (f_eval, fi) in f_evals.iter_mut().zip(fx.iter()) {
            *f_eval *= fi;
        }
    }
    let mut g_evals = vec![F::one(); 1 << gxs[0].num_vars];
    for gx in gxs.iter() {
        for (g_eval, gi) in g_evals.iter_mut().zip(gx.iter()) {
            *g_eval *= gi;
        }
    }
    batch_inversion(&mut g_evals[..]);

    for (f_eval, g_eval) in f_evals.iter_mut().zip(g_evals.iter()) {
        if *g_eval == F::zero() {
            return Err(PolyIOPErrors::InvalidParameters(
                "gxs has zero entries in the boolean hypercube".to_string(),
            ));
        }
        *f_eval *= g_eval;
    }

    end_timer!(start);
    Ok(Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        fxs[0].num_vars,
        f_evals,
    )))
}

/// Compute the product polynomial `prod(x)` such that
/// `prod(x) = [(1-x1)*frac(x2, ..., xn, 0) + x1*prod(x2, ..., xn, 0)] *
/// [(1-x1)*frac(x2, ..., xn, 1) + x1*prod(x2, ..., xn, 1)]` on the boolean
/// hypercube {0,1}^n
///
/// The caller needs to check num_vars matches in f and g
/// Cost: linear in N.
pub(super) fn compute_product_poly<F: PrimeField>(
    frac_poly: &Arc<DenseMultilinearExtension<F>>,
) -> Result<Arc<DenseMultilinearExtension<F>>, PolyIOPErrors> {
    let start = start_timer!(|| "compute evaluations of prod polynomial");
    let num_vars = frac_poly.num_vars;
    let frac_evals = &frac_poly.evaluations;

    // ===================================
    // prod(x)
    // ===================================
    //
    // `prod(x)` can be computed via recursing the following formula for 2^n-1
    // times
    //
    // `prod(x_1, ..., x_n) :=
    //      [(1-x1)*frac(x2, ..., xn, 0) + x1*prod(x2, ..., xn, 0)] *
    //      [(1-x1)*frac(x2, ..., xn, 1) + x1*prod(x2, ..., xn, 1)]`
    //
    // At any given step, the right hand side of the equation
    // is available via either frac_x or the current view of prod_x
    let mut prod_x_evals = vec![];
    for x in 0..(1 << num_vars) - 1 {
        // sign will decide if the evaluation should be looked up from frac_x or
        // prod_x; x_zero_index is the index for the evaluation (x_2, ..., x_n,
        // 0); x_one_index is the index for the evaluation (x_2, ..., x_n, 1);
        let (x_zero_index, x_one_index, sign) = get_index(x, num_vars);
        if !sign {
            prod_x_evals.push(frac_evals[x_zero_index] * frac_evals[x_one_index]);
        } else {
            // sanity check: if we are trying to look up from the prod_x_evals table,
            // then the target index must already exist
            if x_zero_index >= prod_x_evals.len() || x_one_index >= prod_x_evals.len() {
                return Err(PolyIOPErrors::ShouldNotArrive);
            }
            prod_x_evals.push(prod_x_evals[x_zero_index] * prod_x_evals[x_one_index]);
        }
    }

    // prod(1, 1, ..., 1) := 0
    prod_x_evals.push(F::zero());
    end_timer!(start);

    Ok(Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars,
        prod_x_evals,
    )))
}

/// generate the zerocheck proof for the virtual polynomial
///    prod(x) - p1(x) * p2(x) + alpha * [frac(x) * g1(x) * ... * gk(x) - f1(x)
/// * ... * fk(x)] where p1(x) = (1-x1) * frac(x2, ..., xn, 0) + x1 * prod(x2,
///   ..., xn, 0), p2(x) = (1-x1) * frac(x2, ..., xn, 1) + x1 * prod(x2, ...,
///   xn, 1)
/// Returns proof.
///
/// Cost: O(N)
pub(super) fn prove_zero_check<F: PrimeField>(
    fxs: &[Arc<DenseMultilinearExtension<F>>],
    gxs: &[Arc<DenseMultilinearExtension<F>>],
    frac_poly: &Arc<DenseMultilinearExtension<F>>,
    prod_x: &Arc<DenseMultilinearExtension<F>>,
    alpha: &F,
    transcript: &mut IOPTranscript<F>,
) -> Result<(IOPProof<F>, VirtualPolynomial<F>), PolyIOPErrors> {
    let start = start_timer!(|| "zerocheck in product check");
    let num_vars = frac_poly.num_vars;

    // compute p1(x) = (1-x1) * frac(x2, ..., xn, 0) + x1 * prod(x2, ..., xn, 0)
    // compute p2(x) = (1-x1) * frac(x2, ..., xn, 1) + x1 * prod(x2, ..., xn, 1)
    let mut p1_evals = vec![F::zero(); 1 << num_vars];
    let mut p2_evals = vec![F::zero(); 1 << num_vars];
    for x in 0..1 << num_vars {
        let (x0, x1, sign) = get_index(x, num_vars);
        if !sign {
            p1_evals[x] = frac_poly.evaluations[x0];
            p2_evals[x] = frac_poly.evaluations[x1];
        } else {
            p1_evals[x] = prod_x.evaluations[x0];
            p2_evals[x] = prod_x.evaluations[x1];
        }
    }
    let p1 = Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars, p1_evals,
    ));
    let p2 = Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_vars, p2_evals,
    ));

    // compute Q(x)
    // prod(x)
    let mut q_x = VirtualPolynomial::new_from_mle(prod_x, F::one());

    //   prod(x)
    // - p1(x) * p2(x)
    q_x.add_mle_list([p1, p2], -F::one())?;

    //   prod(x)
    // - p1(x) * p2(x)
    // + alpha * frac(x) * g1(x) * ... * gk(x)
    let mut mle_list = gxs.to_vec();
    mle_list.push(frac_poly.clone());
    q_x.add_mle_list(mle_list, *alpha)?;

    //   prod(x)
    // - p1(x) * p2(x)
    // + alpha * frac(x) * g1(x) * ... * gk(x)
    // - alpha * f1(x) * ... * fk(x)]
    q_x.add_mle_list(fxs.to_vec(), -*alpha)?;

    let iop_proof = <PolyIOP<F> as ZeroCheck<F>>::prove(&q_x, transcript)?;

    end_timer!(start);
    Ok((iop_proof, q_x))
}
