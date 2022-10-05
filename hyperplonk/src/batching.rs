//! Sumcheck based batch opening and verify commitment.
// TODO: refactoring this code to somewhere else
// currently IOP depends on PCS because perm check requires commitment.
// The sumcheck based batch opening therefore cannot stay in the PCS repo --
// which creates a cyclic dependency.

use std::rc::Rc;

use arithmetic::{
    bit_decompose, build_eq_x_r, gen_eval_point, DenseMultilinearExtension, VirtualPolynomial,
};
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_std::{end_timer, log2, start_timer, One, Zero};
use pcs::{
    prelude::{
        Commitment, MultilinearKzgBatchProof, MultilinearKzgPCS, MultilinearProverParam,
        UnivariateProverParam,
    },
    PolynomialCommitmentScheme,
};
use poly_iop::{
    prelude::{IOPProof, SumCheck},
    PolyIOP,
};
use transcript::IOPTranscript;

use crate::prelude::HyperPlonkErrors;

pub(crate) struct NewBatchProof<E: PairingEngine> {
    /// commitment to tilde g MLE
    pub(crate) tilde_g_commit: Commitment<E>,
    /// A sum check proof proving tilde g's sum
    pub(crate) sum_check_proof: IOPProof<E::Fr>,
}

/// Steps:
/// 1. todo...
pub(crate) fn multi_open_internal<E: PairingEngine>(
    uni_prover_param: UnivariateProverParam<E::G1Affine>,
    ml_prover_param: MultilinearProverParam<E>,
    polynomials: &[Rc<DenseMultilinearExtension<E::Fr>>],
    points: &[Vec<E::Fr>],
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<NewBatchProof<E>, HyperPlonkErrors> {
    let open_timer = start_timer!(|| "multi open");

    // TODO: sanity checks

    let num_var = polynomials[0].num_vars;
    let k = polynomials.len();
    let ell = log2(k) as usize;

    // challenge point t
    let t = transcript.get_and_append_challenge_vectors("t".as_ref(), ell)?;

    // eq(t, i) for i in [0..k]
    let eq_t_i_list = (0..k)
        .map(|index| get_eq_t_i(t.as_ref(), index, ell))
        .collect::<Vec<_>>();

    // \tilde g(i, b) = eq(t, i) * f_i(b)
    let mut tilde_g_eval = vec![E::Fr::zero(); 1 << (ell + num_var)];
    for (index, f_i) in polynomials.iter().enumerate() {
        for &f_i_eval in f_i.iter() {
            tilde_g_eval.push(f_i_eval * eq_t_i_list[index])
        }
    }
    let tilde_g = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_var + ell,
        tilde_g_eval,
    ));

    // commit to tilde g
    let tilde_g_commit = MultilinearKzgPCS::commit(&(ml_prover_param, uni_prover_param), &tilde_g)?;

    // \tilde eq is the concat of all points, padded with 0s
    let mut tilde_eq_eval = points.iter().flatten().copied().collect::<Vec<E::Fr>>();
    tilde_eq_eval.resize(1 << (num_var + ell), E::Fr::zero());
    let tilde_eq = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        num_var + ell,
        tilde_eq_eval,
    ));

    // built the virtual polynomial for SumCheck
    let mut sum_check_vp = VirtualPolynomial::new(num_var + ell);
    sum_check_vp.add_mle_list([tilde_g, tilde_eq], E::Fr::one())?;

    let proof = <PolyIOP<E::Fr> as SumCheck<E::Fr>>::prove(&sum_check_vp, transcript)?;

    end_timer!(open_timer);
    Ok(NewBatchProof {
        tilde_g_commit,
        sum_check_proof: proof,
    })
}

// generate an eval which is eq(\vec t, <i>)
fn get_eq_t_i<F: PrimeField>(t: &[F], index: usize, index_len: usize) -> F {
    todo!()
}
