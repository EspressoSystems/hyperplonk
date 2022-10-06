//! Sumcheck based batch opening and verify commitment.
// TODO: refactoring this code to somewhere else
// currently IOP depends on PCS because perm check requires commitment.
// The sumcheck based batch opening therefore cannot stay in the PCS repo --
// which creates a cyclic dependency.

use arithmetic::{DenseMultilinearExtension, VirtualPolynomial};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_poly::MultilinearExtension;
use ark_std::{end_timer, log2, start_timer, One, Zero};
use pcs::{
    prelude::{
        Commitment, MultilinearKzgPCS, MultilinearProverParam, MultilinearVerifierParam,
        UnivariateProverParam, UnivariateVerifierParam,
    },
    PolynomialCommitmentScheme,
};
use poly_iop::{
    prelude::{IOPProof, SumCheck},
    PolyIOP,
};
use std::rc::Rc;
use transcript::IOPTranscript;

use crate::prelude::HyperPlonkErrors;

pub(crate) struct NewBatchProof<E: PairingEngine> {
    /// A sum check proof proving tilde g's sum
    pub(crate) sum_check_proof: IOPProof<E::Fr>,
    /// \tilde g(a1, a2)
    pub(crate) tilde_g_eval: E::Fr,
    /// f_i(a2)
    pub(crate) f_i_eval: Vec<E::Fr>,
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
    let merged_num_var = num_var + ell;

    // challenge point t
    let t = transcript.get_and_append_challenge_vectors("t".as_ref(), ell)?;

    // eq(t, i) for i in [0..k]
    let eq_t_i_list = get_eq_x(t.as_ref(), ell);

    // \tilde g(i, b) = eq(t, i) * f_i(b)
    let mut tilde_g_eval = vec![];
    for (index, f_i) in polynomials.iter().enumerate() {
        for &f_i_eval in f_i.iter() {
            tilde_g_eval.push(f_i_eval * eq_t_i_list[index])
        }
    }
    tilde_g_eval.resize(1 << (ell + num_var), E::Fr::zero());
    let tilde_g = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        merged_num_var,
        tilde_g_eval,
    ));

    // evaluate eq(b, z_i) at boolean hypercube
    // merge all evals into a nv + ell mle
    let mut tilde_eq_eval = vec![];
    for point in points.iter() {
        let eq_b_zi = get_eq_x(&point, num_var);
        tilde_eq_eval.extend_from_slice(eq_b_zi.as_slice());
    }
    tilde_eq_eval.resize(1 << (ell + num_var), E::Fr::zero());
    let tilde_eq = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        merged_num_var,
        tilde_eq_eval,
    ));

    // built the virtual polynomial for SumCheck
    let mut sum_check_vp = VirtualPolynomial::new(num_var + ell);
    sum_check_vp.add_mle_list([tilde_g.clone(), tilde_eq], E::Fr::one())?;

    let proof = <PolyIOP<E::Fr> as SumCheck<E::Fr>>::prove(&sum_check_vp, transcript)?;
    let tilde_g_eval = tilde_g.evaluate(&proof.point).unwrap();

    let a2 = &proof.point[ell..];
    let f_i_eval = polynomials
        .iter()
        .map(|p| p.evaluate(a2).unwrap())
        .collect::<Vec<E::Fr>>();

    end_timer!(open_timer);
    Ok(NewBatchProof {
        sum_check_proof: proof,
        tilde_g_eval,
        f_i_eval,
    })
}

/// Steps:
/// 1. todo...
pub(crate) fn batch_internal<E: PairingEngine>(
    uni_prover_param: UnivariateVerifierParam<E>,
    ml_prover_param: MultilinearVerifierParam<E>,
    num_vars: usize,
    f_i_commitments: &[Commitment<E>],

    // f_i(a_2)
    f_i_eval: &[E::Fr],
    proof: &NewBatchProof<E>,
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<bool, HyperPlonkErrors> {
    let open_timer = start_timer!(|| "batch verification");

    // TODO: sanity checks

    let k = f_i_commitments.len();
    let ell = log2(k) as usize;
    let merged_num_var = num_vars + ell;

    // challenge point t
    let t = transcript.get_and_append_challenge_vectors("t".as_ref(), ell)?;

    // sum check point (a1, a2)
    let a1 = &proof.sum_check_proof.point[0..ell];
    let a2 = &proof.sum_check_proof.point[ell..];

    // build g' commitment
    let eq_a1_list = get_eq_x(a1, ell);
    let eq_t_list = get_eq_x(t.as_ref(), ell);

    let mut g_prime_eval = E::Fr::zero();
    let mut g_prime_commit = E::G1Affine::zero().into_projective();
    for i in 0..k {
        let tmp = eq_a1_list[i] * eq_t_list[i];
        g_prime_eval += tmp * f_i_eval[i];
        g_prime_commit += &f_i_commitments[i].0.mul(tmp);
    }

    // ensure g'(a_2) == \tilde g(a1, a2)
    if proof.tilde_g_eval != g_prime_eval {
        return Ok(false);
    }

    // verify commitment
    let res = MultilinearKzgPCS::verify(
        &(ml_prover_param, uni_prover_param),
        &Commitment(g_prime_commit.into_affine()),
        a2.to_vec().as_ref(),
        &g_prime_eval,
        // what is this proof?
        proof,
    )?;

    end_timer!(open_timer);
    Ok(res)
}

// generate a list of evals which are eq(\vec t, <i>) for i in 0..2^index_len
fn get_eq_x<F: PrimeField>(t: &[F], index_len: usize) -> Vec<F> {
    todo!()
}
