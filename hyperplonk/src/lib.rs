//! Main module for the HyperPlonk PolyIOP.

use crate::subroutine::{
    perm_check::{
        estimate_perm_check_param_size, perm_check_prover_subroutine,
        perm_check_verifier_subroutine,
    },
    zero_check::{
        estimate_zero_check_param_size, zero_check_prover_subroutine,
        zero_check_verifier_subroutine,
    },
};
use ark_ec::PairingEngine;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, log2, start_timer, Zero};
use errors::HyperPlonkErrors;
use pcs::prelude::{compute_qx_degree, merge_polynomials, PCSErrors, PolynomialCommitmentScheme};
use poly_iop::{
    prelude::{PermutationCheck, SumCheck, ZeroCheck},
    PolyIOP,
};
use selectors::SelectorColumn;
use std::{cmp::max, rc::Rc};
use structs::{HyperPlonkParams, HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkVerifyingKey};
use transcript::IOPTranscript;
use witness::WitnessColumn;

mod errors;
mod selectors;
mod structs;
mod subroutine;
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
    /// Outputs:
    /// - The HyperPlonk SNARK proof.
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[E::Fr],
        witnesses: &[WitnessColumn<E::Fr>],
    ) -> Result<Self::Proof, HyperPlonkErrors>;

    /// Verify the HyperPlonk proof.
    ///
    /// Inputs:
    /// - `params`: instance parameter
    /// - `pub_input`: online public input
    /// - `proof`: HyperPlonk SNARK proof challenges
    /// Outputs:
    /// - Return a boolean on whether the verification is successful
    fn verify(
        params: &Self::VerifyingKey,
        pub_input: &[E::Fr],
        proof: &Self::Proof,
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
        let (zero_check_nv, zero_check_uni_degree) =
            estimate_zero_check_param_size(num_vars, params.log_n_wires);
        println!(
            "zero check nv: {}, uni {}",
            zero_check_nv, zero_check_uni_degree
        );
        let (perm_check_nv, perm_check_uni_degree) =
            estimate_perm_check_param_size(num_vars, params.log_n_wires);
        let nv = max(zero_check_nv, perm_check_nv);
        let uni_degree = max(zero_check_uni_degree, perm_check_uni_degree);

        // // number of variables in merged polynomial for Multilinear-KZG
        // let merged_nv = num_vars + log_num_witness_polys;
        // // degree of q(x) for Univariate-KZG
        // let q_x_degree = compute_qx_degree(num_vars, 1 << log_num_witness_polys);
        // let prod_x_degree = merged_nv * 4;
        // let supported_uni_degree = max(prod_x_degree, q_x_degree);
        // println!("q(x) degree: {}\nnum_vars: {}\nmerged_vars {}\nprod x {}",
        // q_x_degree, num_vars, merged_nv, prod_x_degree); extract PCS prover
        // and verifier keys from SRS
        let (pcs_prover_param, pcs_verifier_param) =
            PCS::trim(pcs_srs, log2(uni_degree) as usize, Some(nv + 1))?;

        println!("here");
        // build permutation oracles
        let permutation_oracles = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
            num_vars + log_num_witness_polys,
            permutation,
        ));
        println!("here");
        let perm_com = PCS::commit(&pcs_prover_param, &permutation_oracles)?;

        println!("here");
        // build selector oracles and commit to it
        let selector_oracles: Vec<Rc<DenseMultilinearExtension<E::Fr>>> = selectors
            .iter()
            .map(|s| Rc::new(DenseMultilinearExtension::from(s)))
            .collect();

        println!("here");
        let selector_com = selector_oracles
            .iter()
            .map(|poly| PCS::commit(&pcs_prover_param, poly))
            .collect::<Result<Vec<PCS::Commitment>, PCSErrors>>()?;

        println!("here");
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
    ) -> Result<Self::Proof, HyperPlonkErrors> {
        let start = start_timer!(|| "hyperplonk proving");
        let mut transcript = IOPTranscript::<E::Fr>::new(b"hyperplonk");

        // witness assignment of length 2^n
        let num_vars = pk.params.nv;
        let log_num_witness_polys = pk.params.log_n_wires;
        // number of variables in merged polynomial for Multilinear-KZG
        let merged_nv = num_vars + log_num_witness_polys;
        // degree of q(x) for Univariate-KZG
        let _supported_uni_degree = compute_qx_degree(num_vars, 1 << log_num_witness_polys);
        //  online public input of length 2^\ell
        let ell = pk.params.log_pub_input_len;

        println!("here");
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
        println!("here");
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
        let w_merged_com = PCS::commit(&pk.pcs_param, &w_merged)?;

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
        let (
            zero_check_proof,
            witness_zero_check_openings,
            witness_zero_check_evals,
            selector_oracle_openings,
            selector_oracle_evals,
        ) = zero_check_prover_subroutine(pk, &witness_polys, &mut transcript)?;

        // =======================================================================
        // 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracles`, and
        // obtain a PermCheckSubClaim.
        // =======================================================================
        let (
            perm_check_proof,
            witness_perm_check_opening,
            witness_perm_check_eval,
            perm_oracle_opening,
            perm_oracle_eval,
            prod_com,
            prod_opening,
            prod_evals,
        ) = perm_check_prover_subroutine(pk, &witness_polys, &mut transcript)?;

        // =======================================================================
        // 4. Generate evaluations and corresponding proofs
        //
        // - public input consistency checks
        //   - pi_poly(r_pi) where r_pi is sampled from transcript
        // =======================================================================
        let step = start_timer!(|| "opening and evaluations");

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
            prod_evals,
            prod_opening,
            witness_perm_check_opening,
            witness_perm_check_eval,
            perm_oracle_opening,
            perm_oracle_eval,
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

    /// Verify the HyperPlonk proof.
    ///
    /// Inputs:
    /// - `vk`: verification key
    /// - `pub_input`: online public input
    /// - `proof`: HyperPlonk SNARK proof
    /// Outputs:
    /// - Return a boolean on whether the verification is successful
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
    ) -> Result<bool, HyperPlonkErrors> {
        let start = start_timer!(|| "hyperplonk verification");

        let mut transcript = IOPTranscript::<E::Fr>::new(b"hyperplonk");
        let num_var = vk.params.nv;

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
        if !zero_check_verifier_subroutine(vk, proof, &mut transcript)? {
            return Ok(false);
        }

        // =======================================================================
        // 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracles`
        // =======================================================================
        if !perm_check_verifier_subroutine(vk, proof, &mut transcript)? {
            return Ok(false);
        }

        // =======================================================================
        // 3. Verify the opening against the commitment
        // =======================================================================
        let step = start_timer!(|| "verify commitments");

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
    use ark_std::{test_rng, One};
    use pcs::prelude::KZGMultilinearPCS;
    use poly_iop::random_permutation_mle;

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
            &s_perm.evaluations,
            &[q1],
        )?;

        // generate a proof and verify
        let proof = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, KZGMultilinearPCS<E>>>::prove(
            &pk,
            &pi.0,
            &[w1, w2],
        )?;

        let _sub_claim = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, KZGMultilinearPCS<E>>>::verify(
            &vk, &pi.0, &proof,
        )?;

        Ok(())
    }
}
