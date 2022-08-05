//! Main module for the HyperPlonk PolyIOP.

use crate::utils::eval_f;
use arithmetic::VPAuxInfo;
use ark_ec::PairingEngine;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, log2, start_timer, One, Zero};
use errors::HyperPlonkErrors;
use pcs::prelude::{compute_qx_degree, merge_polynomials, PCSErrors, PolynomialCommitmentScheme};
use poly_iop::{
    identity_permutation_mle,
    prelude::{PermutationCheck, SumCheck, ZeroCheck},
    PolyIOP,
};
use selectors::SelectorColumn;
use std::{marker::PhantomData, rc::Rc};
use structs::{HyperPlonkParams, HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkVerifyingKey};
use utils::build_f;
use witness::WitnessColumn;

mod errors;
mod selectors;
mod structs;
mod utils;
mod witness;

/// A trait for HyperPlonk Poly-IOPs.
/// A HyperPlonk is derived from SumChecks, ZeroChecks and PermutationChecks.
pub trait HyperPlonkSNARK<E, PCS>:
    SumCheck<E::Fr> + ZeroCheck<E::Fr> + PermutationCheck<E::Fr>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
{
    type Parameters;
    type ProvingKey;
    type VerifyingKey;
    type Proof;

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
        selectors: &[SelectorColumn<E::Fr>],
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), HyperPlonkErrors>;

    /// Generate HyperPlonk SNARK proof.
    ///
    /// Inputs:
    /// - `pk`: circuit proving key
    /// - `pub_input`: online public input
    /// - `witness`: witness assignment
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - The HyperPlonk SNARK proof.
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[E::Fr],
        witnesses: &[WitnessColumn<E::Fr>],
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, HyperPlonkErrors>;

    /// Verify the HyperPlonk proof and generate the evaluation subclaims to be
    /// checked later by the SNARK verifier.
    ///
    /// Inputs:
    /// - `params`: instance parameter
    /// - `pub_input`: online public input
    /// - `proof`: HyperPlonk SNARK proof
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
    ) -> Result<bool, HyperPlonkErrors>;
}

impl<E, PCS> HyperPlonkSNARK<E, PCS> for PolyIOP<E::Fr>
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
        selectors: &[SelectorColumn<E::Fr>],
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), HyperPlonkErrors> {
        let num_vars = params.nv;
        let log_num_witness_polys = params.log_n_wires;

        // number of variables in merged polynomial for Multilinear-KZG
        let merged_nv = num_vars + log_num_witness_polys;
        // degree of q(x) for Univariate-KZG
        let supported_uni_degree = compute_qx_degree(num_vars, 1 << log_num_witness_polys);

        // extract PCS prover and verifier keys from SRS
        let (pcs_prover_param, pcs_verifier_param) = PCS::trim(
            pcs_srs,
            log2(supported_uni_degree) as usize,
            Some(merged_nv + 1),
        )?;

        // build permutation oracles
        let permutation_oracles = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
            merged_nv,
            permutation,
        ));
        let perm_com = PCS::commit(&pcs_prover_param, &permutation_oracles)?;

        // build selector oracles and commit to it
        let selector_oracles: Vec<Rc<DenseMultilinearExtension<E::Fr>>> = selectors
            .iter()
            .map(|s| Rc::new(DenseMultilinearExtension::from(s)))
            .collect();

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
                perm_com,
            },
        ))
    }

    /// Generate HyperPlonk SNARK proof.
    ///
    /// Inputs:
    /// - `pk`: circuit proving key
    /// - `pub_input`: online public input of length 2^\ell
    /// - `witness`: witness assignment of length 2^n
    /// - `transcript`: the transcript used for generating pseudorandom
    ///   challenges
    /// Outputs:
    /// - The HyperPlonk SNARK proof.
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
    ///
    /// TODO: this function is gigantic -- refactor it to smaller ones
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[E::Fr],
        witnesses: &[WitnessColumn<E::Fr>],
        transcript: &mut Self::Transcript,
    ) -> Result<Self::Proof, HyperPlonkErrors> {
        let start = start_timer!(|| "hyperplonk proving");

        // witness assignment of length 2^n
        let num_vars = pk.params.nv;
        let log_num_witness_polys = pk.params.log_n_wires;
        // number of variables in merged polynomial for Multilinear-KZG
        let merged_nv = num_vars + log_num_witness_polys;
        // degree of q(x) for Univariate-KZG
        let _supported_uni_degree = compute_qx_degree(num_vars, 1 << log_num_witness_polys);
        //  online public input of length 2^\ell
        let ell = pk.params.log_pub_input_len;

        let witness_polys: Vec<Rc<DenseMultilinearExtension<E::Fr>>> = witnesses
            .iter()
            .map(|w| Rc::new(DenseMultilinearExtension::from(w)))
            .collect();
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
        let step = start_timer!(|| "commit witnesses");
        let mut witness_commits = vec![];
        // TODO: batch commit
        for wi_poly in witness_polys.iter() {
            let wi_com = PCS::commit(&pk.pcs_param, wi_poly)?;
            witness_commits.push(wi_com);
        }

        let w_merged = merge_polynomials(&witness_polys)?;
        if w_merged.num_vars != merged_nv {
            return Err(HyperPlonkErrors::InvalidParameters(format!(
                "merged witness poly has a different num_var ({}) from expected ({})",
                w_merged.num_vars, merged_nv
            )));
        }
        let w_merged_com = PCS::commit(&pk.pcs_param, &Rc::new(w_merged.clone()))?;

        transcript.append_serializable_element(b"w", &w_merged_com)?;
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
        let step = start_timer!(|| "ZeroCheck on f");

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
        // 3.6. open `prod(0,x)`, `prod(1, x)`, `prod(x, 0)`, `prod(x, 1)` at
        //      zero_check.point
        // =======================================================================
        let step = start_timer!(|| "Permutation check on w_i(x)");

        // 3.1 `generate_challenge` from current transcript (generate beta, gamma)
        let mut permutation_challenge = Self::generate_challenge(transcript)?;

        // 3.2. `compute_product` to build `prod(x)` etc. from f, g and s_perm

        // This function returns 3 MLEs:
        // - prod(x)
        // - numerator
        // - denominator
        // See function signature for details.
        let prod_x_and_aux_info = Self::compute_prod_evals(
            &permutation_challenge,
            &w_merged,
            &w_merged,
            &pk.permutation_oracles,
        )?;
        let prod_x = Rc::new(prod_x_and_aux_info[0].clone());

        // 3.3 push a commitment of `prod(x)` to the transcript
        let prod_com = PCS::commit(&pk.pcs_param, &prod_x)?;

        // 3.4. `update_challenge` with the updated transcript
        Self::update_challenge(&mut permutation_challenge, transcript, &prod_com)?;

        // 3.5. `prove` to generate the proof
        let perm_check_proof = <Self as PermutationCheck<E::Fr>>::prove(
            &prod_x_and_aux_info,
            &permutation_challenge,
            transcript,
        )?;

        // 3.6 open prod(0,x), prod(1, x), prod(x, 0), prod(x, 1) at zero_check.point
        // prod(0, x)
        let tmp_point = [perm_check_proof.point.as_slice(), &[E::Fr::zero()]].concat();
        let (prod_0_x_opening, prod_0_x_eval) = PCS::open(&pk.pcs_param, &prod_x, &tmp_point)?;
        #[cfg(feature = "extensive_sanity_checks")]
        {
            // sanity check
            let eval = prod_x.evaluate(&tmp_point).ok_or_else(|| {
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
            })?;
            if eval != prod_0_x_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }
        // prod(1, x)
        let tmp_point = [perm_check_proof.point.as_slice(), &[E::Fr::one()]].concat();
        let (prod_1_x_opening, prod_1_x_eval) = PCS::open(&pk.pcs_param, &prod_x, &tmp_point)?;
        #[cfg(feature = "extensive_sanity_checks")]
        {
            // sanity check
            let eval = prod_x.evaluate(&tmp_point).ok_or_else(|| {
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
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
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
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
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
            })?;
            if eval != prod_x_1_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }
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
        let step = start_timer!(|| "opening and evaluations");

        // 4.1 permutation check
        let mut witness_zero_check_evals = vec![];
        let mut witness_zero_check_openings = vec![];
        // TODO: parallelization
        // TODO: Batch opening

        // open permutation check proof
        let (witness_perm_check_opening, witness_perm_check_eval) = PCS::open(
            &pk.pcs_param,
            &Rc::new(w_merged.clone()),
            &perm_check_proof.point,
        )?;

        #[cfg(feature = "extensive_sanity_checks")]
        {
            // sanity checks
            let eval = w_merged.evaluate(&perm_check_proof.point).ok_or_else(|| {
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
            })?;
            if eval != witness_perm_check_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }

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

        // 4.3 public input consistency checks
        let r_pi = transcript.get_and_append_challenge_vectors(b"r_pi", ell)?;

        let (pi_opening, pi_eval) = PCS::open(&pk.pcs_param, &pi_in_w0, &r_pi)?;

        #[cfg(feature = "extensive_sanity_checks")]
        {
            // sanity check
            let eval = pi_poly.evaluate(&r_pi).ok_or_else(|| {
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
            })?;
            if eval != pi_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }

        end_timer!(step);
        end_timer!(start);

        Ok(HyperPlonkProof {
            // =======================================================================
            // PCS components: common
            // =======================================================================
            witness_commits,
            w_merged_com,
            // =======================================================================
            // PCS components: permutation check
            // =======================================================================
            // We do not validate prod(x), this is checked by subclaim
            prod_commit: prod_com,
            prod_evals: vec![prod_0_x_eval, prod_1_x_eval, prod_x_0_eval, prod_x_1_eval],
            prod_openings: vec![
                prod_0_x_opening,
                prod_1_x_opening,
                prod_x_0_opening,
                prod_x_1_opening,
            ],
            witness_perm_check_opening,
            witness_perm_check_eval,
            perm_oracle_opening: s_perm_opening,
            perm_oracle_eval: s_perm_eval,
            // =======================================================================
            // PCS components: zero check
            // =======================================================================
            witness_zero_check_openings,
            witness_zero_check_evals,
            selector_oracle_openings,
            selector_oracle_evals,
            // =======================================================================
            // PCS components: public inputs
            // =======================================================================
            pi_eval,
            pi_opening,
            // =======================================================================
            // IOP components
            // =======================================================================
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
    /// - `proof`: HyperPlonk SNARK proof
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
    /// 4. check subclaim validity // todo
    fn verify(
        vk: &Self::VerifyingKey,
        pub_input: &[E::Fr],
        proof: &Self::Proof,
        transcript: &mut Self::Transcript,
    ) -> Result<bool, HyperPlonkErrors> {
        let start = start_timer!(|| "hyperplonk verification");

        // witness assignment of length 2^n
        let num_var = vk.params.nv;
        let log_num_witness_polys = vk.params.log_n_wires;
        // number of variables in merged polynomial for Multilinear-KZG
        let merged_nv = num_var + log_num_witness_polys;

        //  online public input of length 2^\ell
        let ell = vk.params.log_pub_input_len;

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
        let step = start_timer!(|| "verify zero check");
        // Zero check and sum check have different AuxInfo because `w_merged` and
        // `Prod(x)` have degree and num_vars
        let zero_check_aux_info = VPAuxInfo::<E::Fr> {
            // TODO: get the real max degree from gate_func
            // Here we use 6 is because the test has q[0] * w[0]^5 which is degree 6
            max_degree: 6,
            num_variables: num_var,
            phantom: PhantomData::default(),
        };

        // push witness to transcript
        transcript.append_serializable_element(b"w", &proof.w_merged_com)?;

        let zero_check_sub_claim = <Self as ZeroCheck<E::Fr>>::verify(
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

        end_timer!(step);
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
        let mut challenge = <Self as PermutationCheck<E::Fr>>::generate_challenge(transcript)?;
        <Self as PermutationCheck<E::Fr>>::update_challenge(
            &mut challenge,
            transcript,
            &proof.prod_commit,
        )?;

        let perm_check_sub_claim = <Self as PermutationCheck<E::Fr>>::verify(
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
                    - (proof.witness_perm_check_eval
                        + challenge.beta * s_id_eval
                        + challenge.gamma));

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
        // 3. Verify the opening against the commitment
        // =======================================================================
        let step = start_timer!(|| "verify commitments");

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

        // prod(0, x)
        if !PCS::verify(
            &vk.pcs_param,
            &proof.prod_commit,
            &[perm_check_point.as_slice(), &[E::Fr::zero()]].concat(),
            &proof.prod_evals[0],
            &proof.prod_openings[0],
        )? {
            return Err(HyperPlonkErrors::InvalidProof(
                "pcs verification failed".to_string(),
            ));
        }
        // prod(1, x)
        if !PCS::verify(
            &vk.pcs_param,
            &proof.prod_commit,
            &[perm_check_point.as_slice(), &[E::Fr::one()]].concat(),
            &proof.prod_evals[1],
            &proof.prod_openings[1],
        )? {
            return Err(HyperPlonkErrors::InvalidProof(
                "pcs verification failed".to_string(),
            ));
        }
        // prod(x, 0)
        if !PCS::verify(
            &vk.pcs_param,
            &proof.prod_commit,
            &[&[E::Fr::zero()], perm_check_point.as_slice()].concat(),
            &proof.prod_evals[2],
            &proof.prod_openings[2],
        )? {
            return Err(HyperPlonkErrors::InvalidProof(
                "pcs verification failed".to_string(),
            ));
        }
        // prod(x, 1)
        if !PCS::verify(
            &vk.pcs_param,
            &proof.prod_commit,
            &[&[E::Fr::one()], perm_check_point.as_slice()].concat(),
            &proof.prod_evals[3],
            &proof.prod_openings[3],
        )? {
            return Err(HyperPlonkErrors::InvalidProof(
                "pcs verification failed".to_string(),
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
                perm_check_point,
                eval,
                opening,
            )? {
                return Err(HyperPlonkErrors::InvalidProof(
                    "pcs verification failed".to_string(),
                ));
            }
        }

        // =======================================================================
        // 3.3 public input consistency checks
        // =======================================================================
        let mut r_pi = transcript.get_and_append_challenge_vectors(b"r_pi", ell)?;
        let pi_eval = pi_poly.evaluate(&r_pi).ok_or_else(|| {
            HyperPlonkErrors::InvalidParameters("evaluation dimension does not match".to_string())
        })?;
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
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{selectors::SelectorColumn, structs::CustomizedGates, witness::WitnessColumn};
    use ark_bls12_381::Bls12_381;
    use ark_std::test_rng;
    use pcs::prelude::KZGMultilinearPCS;
    use poly_iop::random_permutation_mle;
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
        // 4 public input
        // 1 selector,
        // 2 witnesses,
        // 2 variables for MLE,
        // 4 wires,
        let gates = CustomizedGates {
            gates: vec![(1, Some(0), vec![0, 0, 0, 0, 0]), (-1, None, vec![1])],
        };
        test_hyperplonk_helper::<Bls12_381>(2, 2, 0, 1, gates)
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
        let merged_nv = nv + log_n_wires;

        let s_perm = random_permutation_mle(merged_nv, &mut rng);

        let q1 = SelectorColumn(vec![E::Fr::one(), E::Fr::one(), E::Fr::one(), E::Fr::one()]);
        // w1 := [0, 1, 2, 3]
        let w1 = WitnessColumn(vec![
            E::Fr::zero(),
            E::Fr::one(),
            E::Fr::from(2u64),
            E::Fr::from(3u64),
        ]);
        // w2 := [0^5, 1^5, 2^5, 3^5]
        let w2 = WitnessColumn(vec![
            E::Fr::zero(),
            E::Fr::one(),
            E::Fr::from(32u64),
            E::Fr::from(243u64),
        ]);
        // public input = w1
        let pi = w1.clone();

        // generate pk and vks
        let (pk, vk) = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, KZGMultilinearPCS<E>>>::preprocess(
            &params,
            &pcs_srs,
            // &perm,
            &s_perm.evaluations,
            &[q1],
        )?;

        // generate a proof and verify
        let mut transcript = IOPTranscript::<E::Fr>::new(b"test hyperplonk");
        let proof = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, KZGMultilinearPCS<E>>>::prove(
            &pk,
            &pi.0,
            &[w1, w2],
            &mut transcript,
        )?;

        let mut transcript = IOPTranscript::<E::Fr>::new(b"test hyperplonk");
        let _sub_claim = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, KZGMultilinearPCS<E>>>::verify(
            &vk,
            &pi.0,
            &proof,
            &mut transcript,
        )?;

        Ok(())
    }
}
