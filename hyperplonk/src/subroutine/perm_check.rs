use crate::{
    errors::HyperPlonkErrors,
    structs::{HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkVerifyingKey},
};
use arithmetic::VPAuxInfo;
use ark_ec::PairingEngine;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, start_timer, One, Zero};
use pcs::prelude::{merge_polynomials, PolynomialCommitmentScheme};
use poly_iop::{
    identity_permutation_mle,
    prelude::{build_prod_partial_eval, IOPProof, PermutationCheck, ZeroCheck},
    PolyIOP,
};
use std::{marker::PhantomData, rc::Rc};
use transcript::IOPTranscript;

/// Estimate the PCS parameter sizes for permutation check
/// Returns
/// - degree of univariate polynomial
/// - number over vars in multilinear polynomial
pub(crate) fn estimate_perm_check_param_size(
    num_vars: usize,
    log_n_wires: usize,
) -> (usize, usize) {
    // we first get the num_var for m_merged that is to be used for perm check
    let merged_nv = num_vars + log_n_wires;
    // we merge all the product into a single MLE
    // whose number variable = merged_nv + 2
    let merged_prod_nv = merged_nv + 2;
    // to batch open its commitment we will need a univariate q(x)
    // whose degree is merged_nv * 4
    let uni_degree = merged_nv * 4;
    (merged_prod_nv, uni_degree)
}

/// Internal function to generate
/// - permutation check proof
/// - PCS openings and evaluations for witness and permutation polynomials at
///   zero check point
pub(crate) fn perm_check_prover_subroutine<E, PCS>(
    pk: &HyperPlonkProvingKey<E, PCS>,
    witness_polys: &[PCS::Polynomial],
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<
    (
        IOPProof<E::Fr>,
        PCS::Proof,
        PCS::Evaluation,
        PCS::Proof,
        PCS::Evaluation,
        PCS::Commitment,
        PCS::BatchProof,
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
    // =======================================================================
    // 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracles`, and
    // obtain a PermCheckSubClaim.
    //
    // 3.1. `generate_challenge` from current transcript (generate beta, gamma)
    // 3.2. `compute_product` to build `prod(x)` etc. from f, g and s_perm
    // 3.3. push a commitment of `prod(x)` to the transcript
    // 3.4. `update_challenge` with the updated transcript
    // 3.5. `prove` to generate the proof
    // 3.6. open `prod(0,x)`, `prod(1, x)`, `prod(x, 0)`, `prod(x, 1)` at
    //      zero_check.point
    // =======================================================================
    let step = start_timer!(|| "Permutation check on w_i(x)");
    let w_merged = merge_polynomials(&witness_polys)?;

    // 3.1 `generate_challenge` from current transcript (generate beta, gamma)
    let mut permutation_challenge =
        <PolyIOP<E::Fr> as PermutationCheck<E::Fr>>::generate_challenge(transcript)?;

    // 3.2. `compute_product` to build `prod(x)` etc. from f, g and s_perm

    // This function returns 3 MLEs:
    // - prod(x)
    // - numerator
    // - denominator
    // See function signature for details.
    let prod_x_and_aux_info = <PolyIOP<E::Fr> as PermutationCheck<E::Fr>>::compute_prod_evals(
        &permutation_challenge,
        &w_merged,
        &w_merged,
        &pk.permutation_oracles,
    )?;
    let prod_x = prod_x_and_aux_info[0].clone();

    // 3.3 push a commitment of `prod(x)` to the transcript
    let prod_com = PCS::commit(&pk.pcs_param, &prod_x)?;

    // 3.4. `update_challenge` with the updated transcript
    <PolyIOP<E::Fr> as PermutationCheck<E::Fr>>::update_challenge(
        &mut permutation_challenge,
        transcript,
        &prod_com,
    )?;

    // 3.5. `prove` to generate the proof
    let perm_check_proof = <PolyIOP<E::Fr> as PermutationCheck<E::Fr>>::prove(
        &prod_x_and_aux_info,
        &permutation_challenge,
        transcript,
    )?;

    // 3.6 open prod(0,x), prod(1, x), prod(x, 0), prod(x, 1) at zero_check.point

    let point_0_x = [perm_check_proof.point.as_slice(), &[E::Fr::zero()]].concat();
    let point_1_x = [perm_check_proof.point.as_slice(), &[E::Fr::one()]].concat();
    let point_x_0 = [&[E::Fr::zero()], perm_check_proof.point.as_slice()].concat();
    let point_x_1 = [&[E::Fr::one()], perm_check_proof.point.as_slice()].concat();

    // prod(0, x)
    let (prod_0_x_opening, prod_0_x_eval) = PCS::open(&pk.pcs_param, &prod_x, &point_0_x)?;
    #[cfg(feature = "extensive_sanity_checks")]
    {
        // sanity check
        let eval = prod_x.evaluate(&point_0_x).ok_or_else(|| {
            HyperPlonkErrors::InvalidParameters("evaluation dimension does not match".to_string())
        })?;
        if eval != prod_0_x_eval {
            return Err(HyperPlonkErrors::InvalidProver(
                "Evaluation is different from PCS opening".to_string(),
            ));
        }
    }
    // prod(1, x)
    let (prod_1_x_opening, prod_1_x_eval) = PCS::open(&pk.pcs_param, &prod_x, &point_1_x)?;
    #[cfg(feature = "extensive_sanity_checks")]
    {
        // sanity check
        let eval = prod_x.evaluate(&point_1_x).ok_or_else(|| {
            HyperPlonkErrors::InvalidParameters("evaluation dimension does not match".to_string())
        })?;
        if eval != prod_1_x_eval {
            return Err(HyperPlonkErrors::InvalidProver(
                "Evaluation is different from PCS opening".to_string(),
            ));
        }
    }
    // prod(x, 0)
    let tmp_point = [&[E::Fr::zero()], perm_check_proof.point.as_slice()].concat();
    let (prod_x_0_opening, prod_x_0_eval) = PCS::open(&pk.pcs_param, &prod_x, &tmp_point)?;
    #[cfg(feature = "extensive_sanity_checks")]
    {
        // sanity check
        let eval = prod_x.evaluate(&tmp_point).ok_or_else(|| {
            HyperPlonkErrors::InvalidParameters("evaluation dimension does not match".to_string())
        })?;

        if eval != prod_x_0_eval {
            return Err(HyperPlonkErrors::InvalidProver(
                "Evaluation is different from PCS opening".to_string(),
            ));
        }
    }
    // prod(x, 1)
    let tmp_point = [&[E::Fr::one()], perm_check_proof.point.as_slice()].concat();
    let (prod_x_1_opening, prod_x_1_eval) = PCS::open(&pk.pcs_param, &prod_x, &tmp_point)?;
    #[cfg(feature = "extensive_sanity_checks")]
    {
        // sanity check
        let eval = prod_x.evaluate(&tmp_point).ok_or_else(|| {
            HyperPlonkErrors::InvalidParameters("evaluation dimension does not match".to_string())
        })?;
        if eval != prod_x_1_eval {
            return Err(HyperPlonkErrors::InvalidProver(
                "Evaluation is different from PCS opening".to_string(),
            ));
        }
    }
    let prod_openings = vec![
        prod_0_x_opening,
        prod_1_x_opening,
        prod_x_0_opening,
        prod_x_1_opening,
    ];
    let prod_evals = vec![prod_0_x_eval, prod_1_x_eval, prod_x_0_eval, prod_x_1_eval];
    println!("prod evals: {:?}\n", prod_evals);



    let (prod_opening, _prod_evals) = 
    PCS::multi_open(
        &pk.pcs_param,
        &prod_com,
        &[
            prod_x.clone(),
            prod_x.clone(),
            prod_x.clone(),
            prod_x.clone(),
        ],
        &[point_0_x, point_1_x,
        point_x_0,point_x_1
        ],
    )?;

    println!("here");
    println!("prod evals: {:?}\n", prod_evals);



    println!("here");

    let prod_partial_mles = build_prod_partial_eval(&prod_x)?;
    let (prod_opening, _prod_evals) = PCS::multi_open(
        &pk.pcs_param,
        &prod_com,
        &prod_partial_mles,
        // &[
        //     prod_x.clone(),
        //     prod_x.clone(),
        //     prod_x.clone(),
        //     prod_x.clone(),
        // ],
        &[perm_check_proof.point.clone(),
        perm_check_proof.point.clone(),perm_check_proof.point.clone(),perm_check_proof.point.clone(),
        
        ],
    )?;

    println!("here");
    println!("prod evals: {:?}\n", prod_evals);

    // =======================================================================
    // 4. Generate evaluations and corresponding proofs
    // - permutation check evaluations and proofs
    //   - wi_poly(r_perm_check) where r_perm_check is from perm_check_proof
    //   - selector_poly(r_perm_check)
    //
    // =======================================================================
    // 4.1 permutation check

    // open permutation check proof
    let (witness_perm_check_opening, witness_perm_check_eval) =
        PCS::open(&pk.pcs_param, &w_merged, &perm_check_proof.point)?;

    #[cfg(feature = "extensive_sanity_checks")]
    {
        // sanity checks
        let eval = w_merged.evaluate(&perm_check_proof.point).ok_or_else(|| {
            HyperPlonkErrors::InvalidParameters("evaluation dimension does not match".to_string())
        })?;
        if eval != witness_perm_check_eval {
            return Err(HyperPlonkErrors::InvalidProver(
                "Evaluation is different from PCS opening".to_string(),
            ));
        }
    }

    // Open permutation polynomial at perm_check_point
    let (s_perm_opening, s_perm_eval) = PCS::open(
        &pk.pcs_param,
        &pk.permutation_oracles,
        &perm_check_proof.point,
    )?;

    #[cfg(feature = "extensive_sanity_checks")]
    {
        // sanity check
        let eval = pk
            .permutation_oracles
            .evaluate(&perm_check_proof.point)
            .ok_or_else(|| {
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
            })?;
        if eval != s_perm_eval {
            return Err(HyperPlonkErrors::InvalidProver(
                "Evaluation is different from PCS opening".to_string(),
            ));
        }
    }

    end_timer!(step);
    Ok((
        perm_check_proof,
        witness_perm_check_opening,
        witness_perm_check_eval,
        s_perm_opening,
        s_perm_eval,
        prod_com,
        prod_opening,
        prod_evals,
    ))
}

/// Internal function to verify the zero check component
/// is correct in the proof.
pub(crate) fn perm_check_verifier_subroutine<E, PCS, ZC, PC>(
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
    ZC: ZeroCheck<E::Fr>,
    PC: PermutationCheck<E::Fr, PermutationProof = IOPProof<E::Fr>>,
{
    let num_var = vk.params.nv;

    let log_num_witness_polys = vk.params.log_n_wires;
    // number of variables in merged polynomial for Multilinear-KZG
    let merged_nv = num_var + log_num_witness_polys;

    // =======================================================================
    // 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracles`
    // =======================================================================
    let step = start_timer!(|| "verify permutation check");
    // Zero check and sum check have different AuxInfo because `w_merged` and
    // `Prod(x)` have degree and num_vars
    let perm_check_aux_info = VPAuxInfo::<E::Fr> {
        // Prod(x) has a max degree of 2
        max_degree: 2,
        // degree of merged poly
        num_variables: merged_nv,
        phantom: PhantomData::default(),
    };
    let mut challenge =
        <PolyIOP<E::Fr> as PermutationCheck<E::Fr>>::generate_challenge(transcript)?;
    <PolyIOP<E::Fr> as PermutationCheck<E::Fr>>::update_challenge(
        &mut challenge,
        transcript,
        &proof.prod_commit,
    )?;

    let perm_check_sub_claim = <PolyIOP<E::Fr> as PermutationCheck<E::Fr>>::verify(
        &proof.perm_check_proof,
        &perm_check_aux_info,
        transcript,
    )?;
    let perm_check_point = &perm_check_sub_claim
        .zero_check_sub_claim
        .sum_check_sub_claim
        .point;

    // check perm check subclaim:
    // proof.witness_perm_check_eval ?= perm_check_sub_claim.expected_eval
    //
    // Q(x) := prod(1,x) - prod(x, 0) * prod(x, 1)
    //       + alpha * (
    //             (g(x) + beta * s_perm(x) + gamma) * prod(0, x)
    //           - (f(x) + beta * s_id(x)   + gamma))
    // where
    // - Q(x) is perm_check_sub_claim.zero_check.exp_eval
    // - prod(1, x) ... from prod(x) evaluated over (1, zero_point)
    // - g(x), f(x) are both w_merged over (zero_point)
    // - s_perm(x) and s_id(x) from vk_param.perm_oracle
    // - alpha, beta, gamma from challenge
    let alpha = challenge
        .alpha
        .ok_or_else(|| HyperPlonkErrors::InvalidVerifier("alpha is not set".to_string()))?;

    let s_id = identity_permutation_mle::<E::Fr>(perm_check_point.len());
    let s_id_eval = s_id.evaluate(perm_check_point).ok_or_else(|| {
        HyperPlonkErrors::InvalidVerifier("unable to evaluate s_id(x)".to_string())
    })?;

    let q_x_rec = proof.prod_evals[1] - proof.prod_evals[2] * proof.prod_evals[3]
        + alpha
            * ((proof.witness_perm_check_eval
                + challenge.beta * proof.perm_oracle_eval
                + challenge.gamma)
                * proof.prod_evals[0]
                - (proof.witness_perm_check_eval + challenge.beta * s_id_eval + challenge.gamma));

    if q_x_rec
        != perm_check_sub_claim
            .zero_check_sub_claim
            .expected_evaluation
    {
        return Err(HyperPlonkErrors::InvalidVerifier(
            "evaluation failed".to_string(),
        ));
    }

    end_timer!(step);
    // =======================================================================
    // 3.1 check permutation check evaluations
    // =======================================================================
    // witness for permutation check
    if !PCS::verify(
        &vk.pcs_param,
        &proof.w_merged_com,
        perm_check_point,
        &proof.witness_perm_check_eval,
        &proof.witness_perm_check_opening,
    )? {
        return Err(HyperPlonkErrors::InvalidProof(
            "pcs verification failed".to_string(),
        ));
    }

    if !PCS::verify(
        &vk.pcs_param,
        &vk.perm_com,
        perm_check_point,
        &proof.perm_oracle_eval,
        &proof.perm_oracle_opening,
    )? {
        return Err(HyperPlonkErrors::InvalidProof(
            "pcs verification failed".to_string(),
        ));
    }

    // prod(x) for permutation check
    // TODO: batch verification

    // // prod(0, x)
    // if !PCS::verify(
    //     &vk.pcs_param,
    //     &proof.prod_commit,
    //     &[perm_check_point.as_slice(), &[E::Fr::zero()]].concat(),
    //     &proof.prod_evals[0],
    //     &proof.prod_openings[0],
    // )? {
    //     return Err(HyperPlonkErrors::InvalidProof(
    //         "pcs verification failed".to_string(),
    //     ));
    // }
    // // prod(1, x)
    // if !PCS::verify(
    //     &vk.pcs_param,
    //     &proof.prod_commit,
    //     &[perm_check_point.as_slice(), &[E::Fr::one()]].concat(),
    //     &proof.prod_evals[1],
    //     &proof.prod_openings[1],
    // )? {
    //     return Err(HyperPlonkErrors::InvalidProof(
    //         "pcs verification failed".to_string(),
    //     ));
    // }
    // // prod(x, 0)
    // if !PCS::verify(
    //     &vk.pcs_param,
    //     &proof.prod_commit,
    //     &[&[E::Fr::zero()], perm_check_point.as_slice()].concat(),
    //     &proof.prod_evals[2],
    //     &proof.prod_openings[2],
    // )? {
    //     return Err(HyperPlonkErrors::InvalidProof(
    //         "pcs verification failed".to_string(),
    //     ));
    // }
    // // prod(x, 1)
    // if !PCS::verify(
    //     &vk.pcs_param,
    //     &proof.prod_commit,
    //     &[&[E::Fr::one()], perm_check_point.as_slice()].concat(),
    //     &proof.prod_evals[3],
    //     &proof.prod_openings[3],
    // )? {
    //     return Err(HyperPlonkErrors::InvalidProof(
    //         "pcs verification failed".to_string(),
    //     ));
    // }
    Ok(true)
}
