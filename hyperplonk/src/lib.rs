//! Main module for the HyperPlonk PolyIOP.

use ark_ec::PairingEngine;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, log2, start_timer, One, Zero};
use errors::HyperPlonkErrors;
use pcs::prelude::{
    compute_qx_degree, get_batched_nv, merge_polynomials, PCSErrors, PolynomialCommitmentScheme,
};
use poly_iop::{
    prelude::{PermutationCheck, SumCheck, VPAuxInfo, ZeroCheck},
    PolyIOP,
};
use selectors::SelectorRow;
use std::{marker::PhantomData, rc::Rc};
use structs::{
    HyperPlonkParams, HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkSubClaim,
    HyperPlonkVerifyingKey,
};
use utils::{build_f, eval_f};
use witness::WitnessRow;

mod errors;
mod selectors;
mod structs;
mod utils;
mod witness;

/// A trait for HyperPlonk Poly-IOPs.
/// A HyperPlonk is derived from SumChecks, ZeroChecks and PermutationChecks.
pub trait HyperPlonkPIOP<E, PCS>:
    SumCheck<E::Fr> + ZeroCheck<E::Fr> + PermutationCheck<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
{
    type Parameters;
    type ProvingKey;
    type VerifyingKey;
    type Proof;
    type SubClaim;

    /// Generate the preprocessed polynomials output by the indexer.
    ///
    /// Inputs:
    /// - `params`: HyperPlonk instance parameters
    /// - `permutation`: the permutation for the copy constraints
    /// - `selectors`: the list of selector vectors for custom gates
    /// Outputs:
    /// - The HyperPlonk proving key, which includes the preprocessed
    ///   polynomials.
    fn preprocess(
        params: &Self::Parameters,
        pcs_srs: &PCS::SRS,
        permutation: &[E::Fr],
        selectors: &[SelectorRow<E::Fr>],
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), HyperPlonkErrors>;

    /// Generate HyperPlonk PIOP proof.
    ///
    /// Inputs:
    /// - `pk`: circuit proving key
    /// - `pub_input`: online public input
    /// - `witness`: witness assignment
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - The HyperPlonk PIOP proof.
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[E::Fr],
        witness: &[WitnessRow<E::Fr>],
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, HyperPlonkErrors>;

    /// Verify the HyperPlonk proof and generate the evaluation subclaims to be
    /// checked later by the SNARK verifier.
    ///
    /// Inputs:
    /// - `params`: instance parameter
    /// - `pub_input`: online public input
    /// - `proof`: HyperPlonk PIOP proof
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - Return error if the verification fails, otherwise return the
    ///   evaluation subclaim
    fn verify(
        params: &Self::VerifyingKey,
        pub_input: &[E::Fr],
        proof: &Self::Proof,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, HyperPlonkErrors>;
}

impl<E, PCS> HyperPlonkPIOP<E, PCS> for PolyIOP<E::Fr>
where
    E: PairingEngine,
    // Ideally we want to access polynomial as PCS::Polynomial, instead of instantiating it here.
    // But since PCS::Polynomial can be both univariate or multivariate in our implementation
    // we cannot bound PCS::Polynomial with a property trait bound.
    PCS: PolynomialCommitmentScheme<
        E,
        Polynomial = Rc<DenseMultilinearExtension<E::Fr>>,
        Point = Vec<E::Fr>,
        Evaluation = E::Fr,
    >,
{
    type Parameters = HyperPlonkParams;
    type ProvingKey = HyperPlonkProvingKey<E, PCS>;
    type VerifyingKey = HyperPlonkVerifyingKey<E, PCS>;
    type Proof = HyperPlonkProof<E, PCS, Self, Self>;
    type SubClaim = HyperPlonkSubClaim<E::Fr, Self, Self>;

    /// Generate the preprocessed polynomials output by the indexer.
    ///
    /// Inputs:
    /// - `params`: HyperPlonk instance parameters
    /// - `permutation`: the permutation for the copy constraints
    /// - `selectors`: the list of selector vectors for custom gates
    /// Outputs:
    /// - The HyperPlonk proving key, which includes the preprocessed
    ///   polynomials.
    fn preprocess(
        params: &Self::Parameters,
        pcs_srs: &PCS::SRS,
        permutation: &[E::Fr],
        selectors: &[SelectorRow<E::Fr>],
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), HyperPlonkErrors> {
        let num_vars = params.nv;
        let log_num_polys = params.log_n_wires;

        // number of variables in merged polynomial for Multilinear-KZG
        let merged_nv = get_batched_nv(num_vars, 1 << log_num_polys);
        // degree of q(x) for Univariate-KZG
        let supported_uni_degree = compute_qx_degree(num_vars, 1 << log_num_polys);

        // extract PCS prover and verifier keys from SRS
        let (pcs_prover_param, pcs_verifier_param) = PCS::trim(
            pcs_srs,
            log2(supported_uni_degree) as usize,
            Some(merged_nv),
        )?;

        // build permutation oracles
        let permutation_oracles = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
            params.nv,
            permutation,
        ));

        // build selector oracles and commit to it
        let selector_oracles = SelectorRow::build_mles(selectors)?;
        let selector_com = selector_oracles
            .iter()
            .map(|poly| PCS::commit(&pcs_prover_param, poly))
            .collect::<Result<Vec<PCS::Commitment>, PCSErrors>>()?;

        Ok((
            Self::ProvingKey {
                params: params.clone(),
                permutation_oracles,
                selector_oracles,
                pcs_param: pcs_prover_param,
            },
            Self::VerifyingKey {
                params: params.clone(),
                pcs_param: pcs_verifier_param,
                selector_com,
            },
        ))
    }

    /// Generate HyperPlonk PIOP proof.
    ///
    /// Inputs:
    /// - `pk`: circuit proving key
    /// - `pub_input`: online public input of length 2^\ell
    /// - `witness`: witness assignment of length 2^n
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - The HyperPlonk PIOP proof.
    ///
    /// Steps:
    ///
    /// 1. Commit Witness polynomials `w_i(x)` and append commitment to
    /// transcript
    ///
    /// 2. Run ZeroCheck on
    ///
    ///     `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`  
    ///
    /// where `f` is the constraint polynomial i.e.,
    /// ```ignore
    ///     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
    ///     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
    /// ```
    /// in vanilla plonk, and obtain a ZeroCheckSubClaim
    ///
    /// 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracles`, and
    /// obtain a PermCheckSubClaim.
    ///
    /// 4. Generate evaluations and corresponding proofs
    /// - permutation check evaluations and proofs
    /// - zero check evaluations and proofs
    /// - public input consistency checks
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[E::Fr],
        witnesses: &[WitnessRow<E::Fr>],
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, HyperPlonkErrors> {
        let start = start_timer!("hyperplonk proving");

        // witness assignment of length 2^n
        let num_vars = pk.params.log_n_wires;
        //  online public input of length 2^\ell
        let ell = pk.params.log_pub_input_len;

        let witness_polys = WitnessRow::<E::Fr>::build_mles(witnesses)?;
        let pi_poly = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
            ell as usize,
            pub_input,
        ));

        // =======================================================================
        // 0. sanity checks
        // =======================================================================
        // public input length
        if pub_input.len() != 1 << ell {
            return Err(HyperPlonkErrors::InvalidProver(format!(
                "Public input length is not correct: got {}, expect {}",
                pub_input.len(),
                1 << ell
            )));
        }
        // witnesses length
        for (i, w) in witnesses.iter().enumerate() {
            if w.0.len() != 1 << num_vars {
                return Err(HyperPlonkErrors::InvalidProver(format!(
                    "{}-th witness length is not correct: got {}, expect {}",
                    i,
                    pub_input.len(),
                    1 << ell
                )));
            }
        }
        // check public input matches witness[0]'s first 2^ell elements
        let pi_in_w0 =
            Rc::new(witness_polys[0].fix_variables(vec![E::Fr::zero(); num_vars - ell].as_ref()));

        if pi_in_w0 != pi_poly {
            return Err(HyperPlonkErrors::InvalidProver(format!(
                "Public input {:?} does not match witness[0] {:?}",
                pi_poly, pi_in_w0,
            )));
        }

        // =======================================================================
        // 1. Commit Witness polynomials `w_i(x)` and append commitment to
        // transcript
        // =======================================================================
        let step = start_timer!("commit witnesses");
        let mut witness_commits = vec![];
        // TODO: batch commit
        for wi_poly in witness_polys.iter() {
            let wi_com = PCS::commit(&pk.pcs_param, wi_poly)?;
            transcript.append_serializable_element(b"w", &wi_com)?;
            witness_commits.push(wi_com);
        }

        end_timer!(step);

        // =======================================================================
        // 2 Run ZeroCheck on
        //
        //     `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`
        //
        // where `f` is the constraint polynomial i.e.,
        //
        //     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
        //     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
        //
        // in vanilla plonk, and obtain a ZeroCheckSubClaim
        // =======================================================================
        let step = start_timer!("ZeroCheck on f");

        let fx = build_f(
            &pk.params.gate_func,
            pk.params.nv,
            &pk.selector_oracles,
            &witness_polys,
        )?;
        let zero_check_proof = <Self as ZeroCheck<E::Fr>>::prove(&fx, transcript)?;
        end_timer!(step);

        // =======================================================================
        // 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracles`, and
        // obtain a PermCheckSubClaim.
        //
        // 3.1. `generate_challenge` from current transcript (generate beta, gamma)
        // 3.2. `compute_product` to build `prod(x)` etc. from f, g and s_perm
        // 3.3. push a commitment of `prod(x)` to the transcript
        // 3.4. `update_challenge` with the updated transcript
        // 3.5. `prove` to generate the proof
        // =======================================================================
        let step = start_timer!("Permutation check on w_i(x)");

        // 3.1 `generate_challenge` from current transcript (generate beta, gamma)
        let mut permutation_challenge = Self::generate_challenge(transcript)?;

        // 3.2. `compute_product` to build `prod(x)` etc. from f, g and s_perm
        // s_perm is the second half of permutation oracle
        let s_perm = pk.permutation_oracles.fix_variables(&[E::Fr::one()]);
        let w_merged = merge_polynomials(&witness_polys)?;

        // This function returns 3 MLEs:
        // - prod(x)
        // - numerator
        // - denominator
        // See function signature for details.
        let prod_x_and_aux_info =
            Self::compute_prod_evals(&permutation_challenge, &w_merged, &w_merged, &s_perm)?;

        // 3.3 push a commitment of `prod(x)` to the transcript
        let prod_com = PCS::commit(&pk.pcs_param, &Rc::new(prod_x_and_aux_info[0].clone()))?;

        // 3.4. `update_challenge` with the updated transcript
        Self::update_challenge(&mut permutation_challenge, transcript, &prod_com)?;

        // 3.5. `prove` to generate the proof
        let perm_check_proof = <Self as PermutationCheck<E::Fr>>::prove(
            &prod_x_and_aux_info,
            &permutation_challenge,
            transcript,
        )?;
        end_timer!(step);

        // =======================================================================
        // 4. Generate evaluations and corresponding proofs
        // - permutation check evaluations and proofs
        //   - wi_poly(r_perm_check) where r_perm_check is from perm_check_proof
        //   - selector_poly(r_perm_check)
        //
        // - zero check evaluations and proofs
        //   - wi_poly(r_zero_check) where r_zero_check is from zero_check_proof
        //   - selector_poly(r_zero_check)
        //
        // - public input consistency checks
        //   - pi_poly(r_pi) where r_pi is sampled from transcript
        // =======================================================================
        let step = start_timer!("opening and evaluations");

        let mut witness_perm_check_evals = vec![];
        let mut witness_zero_check_evals = vec![];
        let mut witness_perm_check_openings = vec![];
        let mut witness_zero_check_openings = vec![];
        // TODO: parallelization
        // TODO: Batch opening
        for wire_poly in witness_polys {
            // Open permutation check proof
            let (perm_proof, perm_eval) =
                PCS::open(&pk.pcs_param, &wire_poly, &perm_check_proof.point)?;

            // Open zero check proof
            let (zero_proof, zero_eval) =
                PCS::open(&pk.pcs_param, &wire_poly, &zero_check_proof.point)?;

            {
                // sanity checks
                if wire_poly.evaluate(&perm_check_proof.point).ok_or(
                    HyperPlonkErrors::InvalidParameters(
                        "evaluation dimension does not match".to_string(),
                    ),
                )? != perm_eval
                {
                    return Err(HyperPlonkErrors::InvalidProver(
                        "Evaluation is different from PCS opening".to_string(),
                    ));
                }
                if wire_poly.evaluate(&zero_check_proof.point).ok_or(
                    HyperPlonkErrors::InvalidParameters(
                        "evaluation dimension does not match".to_string(),
                    ),
                )? != zero_eval
                {
                    return Err(HyperPlonkErrors::InvalidProver(
                        "Evaluation is different from PCS opening".to_string(),
                    ));
                }
            }
            witness_perm_check_evals.push(perm_eval);
            witness_zero_check_evals.push(zero_eval);
            witness_perm_check_openings.push(perm_proof);
            witness_zero_check_openings.push(zero_proof);
        }

        let mut selector_perm_check_evals = vec![];
        let mut selector_zero_check_evals = vec![];
        let mut selector_perm_check_openings = vec![];
        let mut selector_zero_check_openings = vec![];
        // TODO: parallelization
        for selector_poly in pk.selector_oracles.iter() {
            // Open permutation check proof
            let (perm_proof, perm_eval) =
                PCS::open(&pk.pcs_param, &selector_poly, &perm_check_proof.point)?;

            // Open zero check proof
            let (zero_proof, zero_eval) =
                PCS::open(&pk.pcs_param, &selector_poly, &zero_check_proof.point)?;

            {
                // sanity check
                if selector_poly.evaluate(&perm_check_proof.point).ok_or(
                    HyperPlonkErrors::InvalidParameters(
                        "evaluation dimension does not match".to_string(),
                    ),
                )? != perm_eval
                {
                    return Err(HyperPlonkErrors::InvalidProver(
                        "Evaluation is different from PCS opening".to_string(),
                    ));
                }
                if selector_poly.evaluate(&zero_check_proof.point).ok_or(
                    HyperPlonkErrors::InvalidParameters(
                        "evaluation dimension does not match".to_string(),
                    ),
                )? != zero_eval
                {
                    return Err(HyperPlonkErrors::InvalidProver(
                        "Evaluation is different from PCS opening".to_string(),
                    ));
                }
            }

            selector_perm_check_openings.push(perm_proof);
            selector_perm_check_evals.push(perm_eval);
            selector_zero_check_openings.push(zero_proof);
            selector_zero_check_evals.push(zero_eval);
        }

        // - public input consistency checks
        let r_pi = transcript.get_and_append_challenge_vectors(b"r_pi", ell)?;
        let (pi_opening, pi_eval) = PCS::open(&pk.pcs_param, &pi_poly, &zero_check_proof.point)?;
        {
            // sanity check
            if pi_poly
                .evaluate(&r_pi)
                .ok_or(HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                ))?
                != pi_eval
            {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }
        end_timer!(step);

        end_timer!(start);

        Ok(HyperPlonkProof {
            // PCS components
            witness_commits,
            witness_perm_check_openings,
            witness_zero_check_openings,
            witness_perm_check_evals,
            witness_zero_check_evals,
            selector_perm_check_openings,
            selector_zero_check_openings,
            selector_perm_check_evals,
            selector_zero_check_evals,
            pi_eval,
            pi_opening,
            // IOP components
            zero_check_proof,
            perm_check_proof,
        })
    }

    /// Verify the HyperPlonk proof and generate the evaluation subclaims to be
    /// checked later by the SNARK verifier.
    ///
    /// Inputs:
    /// - `params`: instance parameter
    /// - `pub_input`: online public input
    /// - `proof`: HyperPlonk PIOP proof
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - Return error if the verification fails, otherwise return the
    ///   evaluation subclaim
    ///
    /// 1. Verify zero_check_proof on
    ///
    ///     `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`
    ///
    /// where `f` is the constraint polynomial i.e.,
    /// ```ignore
    ///     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
    ///     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
    /// ```
    /// in vanilla plonk, and obtain a ZeroCheckSubClaim
    ///
    /// 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracles`
    ///
    /// 3. Verify the opening against the commitment:
    /// - check permutation check evaluations
    /// - check zero check evaluations
    /// - public input consistency checks
    fn verify(
        vk: &Self::VerifyingKey,
        pub_input: &[E::Fr],
        proof: &Self::Proof,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::SubClaim, HyperPlonkErrors> {
        let start = start_timer!("hyperplonk verification");

        // witness assignment of length 2^n
        let num_var = vk.params.log_n_wires;
        //  online public input of length 2^\ell
        let ell = vk.params.log_pub_input_len;

        let aux_info = VPAuxInfo::<E::Fr> {
            // TODO: get the real max degree from gate_func
            max_degree: 2,
            num_variables: num_var,
            phantom: PhantomData::default(),
        };
        let pi_poly = DenseMultilinearExtension::from_evaluations_slice(ell as usize, pub_input);

        // =======================================================================
        // 0. sanity checks
        // =======================================================================
        // public input length
        if pub_input.len() != 1 << ell {
            return Err(HyperPlonkErrors::InvalidProver(format!(
                "Public input length is not correct: got {}, expect {}",
                pub_input.len(),
                1 << ell
            )));
        }

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
        let step = start_timer!("verify zero check");

        // push witness to transcript
        for wi_com in proof.witness_commits.iter() {
            transcript.append_serializable_element(b"w", wi_com)?;
        }
        let zero_check_sub_claim =
            <Self as ZeroCheck<E::Fr>>::verify(&proof.zero_check_proof, &aux_info, transcript)?;
        end_timer!(step);

        // =======================================================================
        // 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracles`
        // =======================================================================
        let step = start_timer!("verify permutation check");
        let perm_check_sub_claim = <Self as PermutationCheck<E::Fr>>::verify(
            &proof.perm_check_proof,
            &aux_info,
            transcript,
        )?;
        end_timer!(step);

        // =======================================================================
        // 3. Verify the opening against the commitment
        // =======================================================================

        let step = start_timer!("verify commitments");
        let perm_point = &perm_check_sub_claim
            .zero_check_sub_claim
            .sum_check_sub_claim
            .point;
        let zero_point = &zero_check_sub_claim.sum_check_sub_claim.point;

        // =======================================================================
        // 3.1 check permutation check evaluations
        // =======================================================================
        // witness for permutation check
        for (commitment, (opening, eval)) in proof.witness_commits.iter().zip(
            proof
                .witness_perm_check_openings
                .iter()
                .zip(proof.witness_perm_check_evals.iter()),
        ) {
            if !PCS::verify(&vk.pcs_param, commitment, &perm_point, eval, opening)? {
                return Err(HyperPlonkErrors::InvalidProof(
                    "pcs verification failed".to_string(),
                ));
            }
        }

        // selector for permutation check
        for (commitment, (opening, eval)) in proof.witness_commits.iter().zip(
            proof
                .witness_perm_check_openings
                .iter()
                .zip(proof.witness_perm_check_evals.iter()),
        ) {
            if !PCS::verify(&vk.pcs_param, commitment, &perm_point, eval, opening)? {
                return Err(HyperPlonkErrors::InvalidProof(
                    "pcs verification failed".to_string(),
                ));
            }
        }

        // =======================================================================
        // 3.2 check zero check evaluations
        // =======================================================================
        // witness for zero check
        for (commitment, (opening, eval)) in proof.witness_commits.iter().zip(
            proof
                .witness_zero_check_openings
                .iter()
                .zip(proof.witness_zero_check_evals.iter()),
        ) {
            if !PCS::verify(&vk.pcs_param, commitment, &zero_point, eval, opening)? {
                return Err(HyperPlonkErrors::InvalidProof(
                    "pcs verification failed".to_string(),
                ));
            }
        }

        // selector for zero check
        for (commitment, (opening, eval)) in proof.witness_commits.iter().zip(
            proof
                .witness_zero_check_openings
                .iter()
                .zip(proof.witness_zero_check_evals.iter()),
        ) {
            if !PCS::verify(&vk.pcs_param, commitment, &zero_point, eval, opening)? {
                return Err(HyperPlonkErrors::InvalidProof(
                    "pcs verification failed".to_string(),
                ));
            }
        }

        let f_eval = eval_f(
            &vk.params.gate_func,
            &proof.selector_zero_check_evals,
            &proof.witness_zero_check_evals,
        )?;

        if f_eval != zero_check_sub_claim.sum_check_sub_claim.expected_evaluation {
            return Err(HyperPlonkErrors::InvalidProof(
                "zero check evaluation failed".to_string(),
            ));
        }

        // =======================================================================
        // 3.3 public input consistency checks
        // =======================================================================
        let mut r_pi = transcript.get_and_append_challenge_vectors(b"r_pi", ell)?;
        let pi_eval = pi_poly
            .evaluate(&r_pi)
            .ok_or(HyperPlonkErrors::InvalidParameters(
                "evaluation dimension does not match".to_string(),
            ))?;
        r_pi = [vec![E::Fr::zero(); num_var - ell], r_pi].concat();
        if !PCS::verify(
            &vk.pcs_param,
            &proof.witness_commits[0],
            &r_pi,
            &pi_eval,
            &proof.pi_opening,
        )? {
            return Err(HyperPlonkErrors::InvalidProof(
                "pcs verification failed".to_string(),
            ));
        }

        end_timer!(step);

        end_timer!(start);

        Ok(HyperPlonkSubClaim {
            zero_check_sub_claim,
            perm_check_sub_claim,
            pub_input_sub_claim: (vec![], E::Fr::default()), // FIXME
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::CustomizedGates;
    use ark_bls12_381::Bls12_381;
    use ark_std::{test_rng, UniformRand};
    use pcs::prelude::KZGMultilinearPCS;
    use transcript::IOPTranscript;

    #[test]
    fn test_hyperplonk_e2e() -> Result<(), HyperPlonkErrors> {
        // Example:
        //     q_L(X) * W_1(X)^5 - W_2(X) = 0
        // is represented as
        // vec![
        //     ( 1,    Some(id_qL),    vec![id_W1, id_W1, id_W1, id_W1, id_W1]),
        //     (-1,    None,           vec![id_W2])
        // ]
        //
        // 1 public input
        // 1 selector,
        // 2 witnesses,
        // 2 variables for MLE,
        // 4 wires,
        let gates = CustomizedGates {
            gates: vec![(1, Some(0), vec![0, 0, 0, 0, 0]), (-1, None, vec![1])],
        };
        test_hyperplonk_helper::<Bls12_381>(2, 0, 0, 1, gates)
    }

    fn test_hyperplonk_helper<E: PairingEngine>(
        nv: usize,
        log_pub_input_len: usize,
        log_n_selectors: usize,
        log_n_wires: usize,
        gate_func: CustomizedGates,
    ) -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        // system parameters
        let params = HyperPlonkParams {
            nv,
            log_pub_input_len,
            log_n_selectors,
            log_n_wires,
            gate_func,
        };
        let pcs_srs = KZGMultilinearPCS::<E>::gen_srs_for_testing(&mut rng, 15)?;

        let permutation: Vec<E::Fr> = (0..1 << nv).map(|_| E::Fr::rand(&mut rng)).collect();
        let selectors = SelectorRow::<E::Fr>::rand_selectors(&mut rng, nv, 1);
        let w1 = WitnessRow(vec![
            E::Fr::one(),
            E::Fr::one(),
            E::Fr::zero(),
            E::Fr::zero(),
        ]);
        let w2 = WitnessRow(vec![
            E::Fr::from(2u64),
            E::Fr::from(32u64),
            E::Fr::zero(),
            E::Fr::zero(),
        ]);

        // generate pk and vks
        let (pk, vk) = <PolyIOP<E::Fr> as HyperPlonkPIOP<E, KZGMultilinearPCS<E>>>::preprocess(
            &params,
            &pcs_srs,
            &permutation,
            &selectors,
        )?;

        // generate a proof and verify
        let mut transcript = IOPTranscript::<E::Fr>::new(b"test hyperplonk");
        let proof = <PolyIOP<E::Fr> as HyperPlonkPIOP<E, KZGMultilinearPCS<E>>>::prove(
            &pk,
            &[],
            &[w1.clone(), w2.clone(), w1, w2],
            &mut transcript,
        )?;
        let mut transcript = IOPTranscript::<E::Fr>::new(b"test hyperplonk");
        let _sub_claim = <PolyIOP<E::Fr> as HyperPlonkPIOP<E, KZGMultilinearPCS<E>>>::verify(
            &vk,
            &[],
            &proof,
            &mut transcript,
        )?;

        Ok(())
    }
}
