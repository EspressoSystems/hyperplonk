use crate::{
    errors::HyperPlonkErrors,
    selectors::SelectorColumn,
    structs::{HyperPlonkParams, HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkVerifyingKey},
    utils::{build_f, eval_f},
    witness::WitnessColumn,
};
use arithmetic::{DenseMultilinearExtension, VPAuxInfo};
use ark_ec::PairingEngine;
use ark_poly::MultilinearExtension;
use ark_std::{end_timer, log2, start_timer, One, Zero};
use pcs::prelude::{compute_qx_degree, merge_polynomials, PCSErrors, PolynomialCommitmentScheme};
use poly_iop::{
    identity_permutation_mle,
    prelude::{IOPProof, PermutationCheck, PolyIOP, SumCheck, ZeroCheck},
};
use std::{marker::PhantomData, rc::Rc};
use transcript::IOPTranscript;

/// Internal function to generate
/// - zero check proof
/// - PCS openings and evaluations for witness and selector polynomials at zero
///   check point
#[allow(clippy::type_complexity)]
pub(crate) fn zero_check_prover_subroutine<E, PCS>(
    pk: &HyperPlonkProvingKey<E, PCS>,
    witness_polys: &[PCS::Polynomial],
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<
    (
        IOPProof<E::Fr>,
        Vec<PCS::Proof>,
        Vec<PCS::Evaluation>,
        Vec<PCS::Proof>,
        Vec<PCS::Evaluation>,
    ),
    HyperPlonkErrors,
>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<
        E,
        Polynomial = Rc<DenseMultilinearExtension<E::Fr>>,
        Point = Vec<E::Fr>,
        Evaluation = E::Fr,
    >,
{
    let step = start_timer!(|| "ZeroCheck on f");

    let fx = build_f(
        &pk.params.gate_func,
        pk.params.nv,
        &pk.selector_oracles,
        &witness_polys,
    )?;

    let zero_check_proof = <PolyIOP<E::Fr> as ZeroCheck<E::Fr>>::prove(&fx, transcript)?;

    let mut witness_zero_check_evals = vec![];
    let mut witness_zero_check_openings = vec![];
    // 4.2 open zero check proof
    // TODO: batch opening
    for wire_poly in witness_polys {
        // Open zero check proof
        let (zero_proof, zero_eval) =
            PCS::open(&pk.pcs_param, &wire_poly, &zero_check_proof.point)?;
        {
            let eval = wire_poly.evaluate(&zero_check_proof.point).ok_or_else(|| {
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
            })?;
            if eval != zero_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }
        witness_zero_check_evals.push(zero_eval);
        witness_zero_check_openings.push(zero_proof);
    }

    // Open selector polynomial at zero_check_point
    let mut selector_oracle_openings = vec![];
    let mut selector_oracle_evals = vec![];

    // TODO: parallelization
    for selector_poly in pk.selector_oracles.iter() {
        // Open zero check proof
        // during verification, use this eval against subclaim
        let (zero_proof, zero_eval) =
            PCS::open(&pk.pcs_param, selector_poly, &zero_check_proof.point)?;

        #[cfg(feature = "extensive_sanity_checks")]
        {
            let eval = selector_poly
                .evaluate(&zero_check_proof.point)
                .ok_or_else(|| {
                    HyperPlonkErrors::InvalidParameters(
                        "evaluation dimension does not match".to_string(),
                    )
                })?;
            if eval != zero_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }
        selector_oracle_openings.push(zero_proof);
        selector_oracle_evals.push(zero_eval);
    }

    end_timer!(step);

    Ok((
        zero_check_proof,
        witness_zero_check_openings,
        witness_zero_check_evals,
        selector_oracle_openings,
        selector_oracle_evals,
    ))
}

#[allow(clippy::type_complexity)]
pub(crate) fn zero_check_verifier_subroutine<E, PCS, ZC, PC>(
    vk: &HyperPlonkVerifyingKey<E, PCS>,
    proof: &HyperPlonkProof<E, PCS, ZC, PC>,
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<bool, HyperPlonkErrors>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<
        E,
        Polynomial = Rc<DenseMultilinearExtension<E::Fr>>,
        Point = Vec<E::Fr>,
        Evaluation = E::Fr,
    >,
    ZC: ZeroCheck<E::Fr, ZeroCheckProof = IOPProof<E::Fr>>,
    PC: PermutationCheck<E::Fr>,
{
    // =======================================================================
    // 1. Verify zero_check_proof on
    //     `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`
    //
    // where `f` is the constraint polynomial i.e.,
    //
    //     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
    //     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
    //
    // =======================================================================
    let step = start_timer!(|| "verify zero check");

    let num_var = vk.params.nv;

    // Zero check and sum check have different AuxInfo because `w_merged` and
    // `Prod(x)` have degree and num_vars
    let zero_check_aux_info = VPAuxInfo::<E::Fr> {
        // TODO: get the real max degree from gate_func
        // Here we use 6 is because the test has q[0] * w[0]^5 which is
        // degree 6
        max_degree: 6,
        num_variables: num_var,
        phantom: PhantomData::default(),
    };

    // push witness to transcript
    transcript.append_serializable_element(b"w", &proof.w_merged_com)?;

    let zero_check_sub_claim = <PolyIOP<E::Fr> as ZeroCheck<E::Fr>>::verify(
        &proof.zero_check_proof,
        &zero_check_aux_info,
        transcript,
    )?;

    let zero_check_point = &zero_check_sub_claim.sum_check_sub_claim.point;

    // check zero check subclaim
    let f_eval = eval_f(
        &vk.params.gate_func,
        &proof.selector_oracle_evals,
        &proof.witness_zero_check_evals,
    )?;
    if f_eval != zero_check_sub_claim.expected_evaluation {
        return Err(HyperPlonkErrors::InvalidProof(
            "zero check evaluation failed".to_string(),
        ));
    }

    // =======================================================================
    // 3.2 check zero check evaluations
    // =======================================================================
    // witness for zero check
    // TODO: batch verification
    for (commitment, (opening, eval)) in proof.witness_commits.iter().zip(
        proof
            .witness_zero_check_openings
            .iter()
            .zip(proof.witness_zero_check_evals.iter()),
    ) {
        if !PCS::verify(&vk.pcs_param, commitment, zero_check_point, eval, opening)? {
            return Err(HyperPlonkErrors::InvalidProof(
                "pcs verification failed".to_string(),
            ));
        }
    }

    // selector for zero check
    // TODO: for now we only support a single selector polynomial
    for (opening, eval) in proof
        .selector_oracle_openings
        .iter()
        .zip(proof.selector_oracle_evals.iter())
    {
        if !PCS::verify(
            &vk.pcs_param,
            &vk.selector_com[0],
            zero_check_point,
            eval,
            opening,
        )? {
            return Err(HyperPlonkErrors::InvalidProof(
                "pcs verification failed".to_string(),
            ));
        }
    }

    end_timer!(step);
    Ok(true)
}
