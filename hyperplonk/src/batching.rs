//! Sumcheck based batch opening and verify commitment.
// TODO: refactoring this code to somewhere else
// currently IOP depends on PCS because perm check requires commitment.
// The sumcheck based batch opening therefore cannot stay in the PCS repo --
// which creates a cyclic dependency.

use arithmetic::{
    build_eq_x_r_vec, fix_last_variables, DenseMultilinearExtension, VPAuxInfo, VirtualPolynomial,
};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_poly::MultilinearExtension;
use ark_std::{end_timer, log2, start_timer, One, Zero};
use pcs::{prelude::Commitment, PolynomialCommitmentScheme};
use poly_iop::{
    prelude::{IOPProof, SumCheck},
    PolyIOP,
};
use std::{marker::PhantomData, rc::Rc};
use transcript::IOPTranscript;

use crate::prelude::HyperPlonkErrors;

#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct NewBatchProof<E, PCS>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
{
    /// A sum check proof proving tilde g's sum
    pub(crate) sum_check_proof: IOPProof<E::Fr>,
    /// \tilde g(a1, a2)
    pub(crate) tilde_g_eval: E::Fr,
    /// f_i(a2)
    pub(crate) f_i_eval_at_a2: Vec<E::Fr>,
    /// f_i(point_i)
    pub(crate) f_i_eval_at_point_i: Vec<E::Fr>,
    /// proof for g'(a_2)
    pub(crate) g_prime_proof: PCS::Proof,
}

/// Steps:
/// 1. todo...
pub(crate) fn multi_open_internal<E, PCS>(
    prover_param: &PCS::ProverParam,
    polynomials: &[PCS::Polynomial],
    points: &[PCS::Point],
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<NewBatchProof<E, PCS>, HyperPlonkErrors>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<
        E,
        Polynomial = Rc<DenseMultilinearExtension<E::Fr>>,
        Point = Vec<E::Fr>,
        Evaluation = E::Fr,
    >,
{
    let open_timer = start_timer!(|| format!("multi open {} points", points.len()));

    // TODO: sanity checks

    let num_var = polynomials[0].num_vars;
    let k = polynomials.len();
    let ell = log2(k) as usize;
    let merged_num_var = num_var + ell;

    // println!("ell {}, num_var {}", ell, num_var);

    // challenge point t
    let t = transcript.get_and_append_challenge_vectors("t".as_ref(), ell)?;

    // eq(t, i) for i in [0..k]
    let eq_t_i_list = build_eq_x_r_vec(t.as_ref())?;

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
        let eq_b_zi = build_eq_x_r_vec(&point)?;
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

    // (a1, a2) := sumcheck's point
    let step = start_timer!(|| "open at a2");
    let a1 = &proof.point[num_var..];
    let a2 = &proof.point[..num_var];
    let f_i_eval_at_a2 = polynomials
        .iter()
        .map(|p| p.evaluate(a2).unwrap())
        .collect::<Vec<E::Fr>>();
    end_timer!(step);

    // build g'(a2)
    let step = start_timer!(|| "evaluate at a2");
    let g_prime = Rc::new(fix_last_variables(&tilde_g, a1));
    end_timer!(step);

    let (g_prime_proof, g_prime_eval) = PCS::open(prover_param, &g_prime, a2.to_vec().as_ref())?;
    assert_eq!(g_prime_eval, tilde_g_eval);

    let step = start_timer!(|| "evaluate fi(pi)");
    let f_i_eval_at_point_i = polynomials
        .iter()
        .zip(points.iter())
        .map(|(f, p)| f.evaluate(p).unwrap())
        .collect();
    end_timer!(step);
    end_timer!(open_timer);
    Ok(NewBatchProof {
        sum_check_proof: proof,
        tilde_g_eval,
        f_i_eval_at_a2,
        f_i_eval_at_point_i,
        g_prime_proof,
    })
}

/// Steps:
/// 1. todo...
pub(crate) fn batch_verify_internal<E, PCS>(
    verifier_param: &PCS::VerifierParam,
    f_i_commitments: &[Commitment<E>],
    proof: &NewBatchProof<E, PCS>,
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<bool, HyperPlonkErrors>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<
        E,
        Polynomial = Rc<DenseMultilinearExtension<E::Fr>>,
        Point = Vec<E::Fr>,
        Evaluation = E::Fr,
        Commitment = Commitment<E>,
    >,
{
    let open_timer = start_timer!(|| "batch verification");

    // TODO: sanity checks

    let k = f_i_commitments.len();
    let ell = log2(k) as usize;
    let num_var = proof.sum_check_proof.point.len() - ell;
    // println!("ell {}, num_var {}", ell, num_var);

    // challenge point t
    let t = transcript.get_and_append_challenge_vectors("t".as_ref(), ell)?;

    // sum check point (a1, a2)
    let a1 = &proof.sum_check_proof.point[num_var..];
    let a2 = &proof.sum_check_proof.point[..num_var];

    // build g' commitment
    let eq_a1_list = build_eq_x_r_vec(a1)?;
    let eq_t_list = build_eq_x_r_vec(t.as_ref())?;

    let mut g_prime_eval = E::Fr::zero();
    let mut g_prime_commit = E::G1Affine::zero().into_projective();
    for i in 0..k {
        let tmp = eq_a1_list[i] * eq_t_list[i];
        g_prime_eval += tmp * proof.f_i_eval_at_a2[i];
        g_prime_commit += &f_i_commitments[i].0.mul(tmp);
    }

    // ensure g'(a_2) == \tilde g(a1, a2)
    if proof.tilde_g_eval != g_prime_eval {
        // println!("eval not match");
        return Ok(false);
    }

    // ensure \sum_i eq(t, <i>) * f_i_evals matches the sum via SumCheck
    // verification
    let mut sum = E::Fr::zero();
    for i in 0..k {
        sum += eq_t_list[i] * proof.f_i_eval_at_point_i[i];
    }
    let aux_info = VPAuxInfo {
        max_degree: 2,
        num_variables: num_var + ell,
        phantom: PhantomData,
    };
    let _subclaim = <PolyIOP<E::Fr> as SumCheck<E::Fr>>::verify(
        sum,
        &proof.sum_check_proof,
        &aux_info,
        transcript,
    )?;

    // verify commitment
    let res = PCS::verify(
        verifier_param,
        &Commitment(g_prime_commit.into_affine()),
        a2.to_vec().as_ref(),
        &g_prime_eval,
        &proof.g_prime_proof,
    )?;

    // println!("res {}", res);
    end_timer!(open_timer);
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use arithmetic::get_batched_nv;
    use ark_bls12_381::Bls12_381 as E;
    use ark_ec::PairingEngine;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::{
        log2,
        rand::{CryptoRng, RngCore},
        test_rng,
        vec::Vec,
        UniformRand,
    };
    use pcs::{
        prelude::{
            compute_qx_degree, MultilinearKzgPCS, MultilinearUniversalParams,
            UnivariateUniversalParams,
        },
        StructuredReferenceString,
    };

    type Fr = <E as PairingEngine>::Fr;

    fn test_multi_open_helper<R: RngCore + CryptoRng>(
        uni_params: &UnivariateUniversalParams<E>,
        ml_params: &MultilinearUniversalParams<E>,
        polys: &[Rc<DenseMultilinearExtension<Fr>>],
        rng: &mut R,
    ) -> Result<(), HyperPlonkErrors> {
        let merged_nv = get_batched_nv(polys[0].num_vars(), polys.len());
        let qx_degree = compute_qx_degree(merged_nv, polys.len());
        let padded_qx_degree = 1usize << log2(qx_degree);

        let (uni_ck, uni_vk) = uni_params.trim(padded_qx_degree)?;
        let (ml_ck, ml_vk) = ml_params.trim(merged_nv)?;

        let mut points = Vec::new();
        for poly in polys.iter() {
            let point = (0..poly.num_vars())
                .map(|_| Fr::rand(rng))
                .collect::<Vec<Fr>>();
            points.push(point);
        }

        // let evals = polys
        //     .iter()
        //     .zip(points.iter())
        //     .map(|(f, p)| f.evaluate(p).unwrap())
        //     .collect::<Vec<_>>();

        let commitments = polys
            .iter()
            .map(|poly| MultilinearKzgPCS::commit(&(ml_ck.clone(), uni_ck.clone()), poly).unwrap())
            .collect::<Vec<_>>();

        let mut transcript = IOPTranscript::new("test transcript".as_ref());
        transcript.append_field_element("init".as_ref(), &Fr::zero())?;

        println!("prove");
        let batch_proof = multi_open_internal::<E, MultilinearKzgPCS<E>>(
            &(ml_ck.clone(), uni_ck.clone()),
            polys,
            &points,
            &mut transcript,
        )?;

        // good path
        println!("verify");
        let mut transcript = IOPTranscript::new("test transcript".as_ref());
        transcript.append_field_element("init".as_ref(), &Fr::zero())?;
        assert!(batch_verify_internal::<E, MultilinearKzgPCS<E>>(
            &(ml_vk, uni_vk),
            &commitments,
            &batch_proof,
            &mut transcript
        )?);

        Ok(())
    }

    #[test]
    fn test_multi_open_internal() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();

        let uni_params =
            UnivariateUniversalParams::<E>::gen_srs_for_testing(&mut rng, 1usize << 15)?;
        let ml_params = MultilinearUniversalParams::<E>::gen_srs_for_testing(&mut rng, 15)?;
        for num_poly in 2..10 {
            for nv in 1..5 {
                let polys1: Vec<_> = (0..num_poly)
                    .map(|_| Rc::new(DenseMultilinearExtension::rand(nv, &mut rng)))
                    .collect();
                test_multi_open_helper(&uni_params, &ml_params, &polys1, &mut rng)?;
            }
        }

        Ok(())
    }
}
