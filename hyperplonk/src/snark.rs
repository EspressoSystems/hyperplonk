use crate::{
    batching::batch_verify_internal,
    errors::HyperPlonkErrors,
    structs::{HyperPlonkIndex, HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkVerifyingKey},
    utils::{build_f, eval_f, prover_sanity_check, PcsAccumulator},
    witness::WitnessColumn,
    HyperPlonkSNARK,
};
use arithmetic::{
    evaluate_opt, gen_eval_point, identity_permutation_mle, merge_polynomials, VPAuxInfo,
};
use ark_ec::PairingEngine;
use ark_poly::DenseMultilinearExtension;
use ark_std::{end_timer, log2, start_timer, One, Zero};
use pcs::prelude::{compute_qx_degree, Commitment, PolynomialCommitmentScheme};
use poly_iop::{
    prelude::{PermutationCheck, ZeroCheck},
    PolyIOP,
};
use std::{cmp::max, marker::PhantomData, rc::Rc};
use transcript::IOPTranscript;

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
        Commitment = Commitment<E>,
    >,
{
    type Index = HyperPlonkIndex<E::Fr>;
    type ProvingKey = HyperPlonkProvingKey<E, PCS>;
    type VerifyingKey = HyperPlonkVerifyingKey<E, PCS>;
    type Proof = HyperPlonkProof<E, Self, PCS>;

    fn preprocess(
        index: &Self::Index,
        pcs_srs: &PCS::SRS,
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), HyperPlonkErrors> {
        let num_vars = index.num_variables();

        let log_num_witness_polys = log2(index.num_witness_columns()) as usize;
        let log_num_selector_polys = log2(index.num_selector_columns()) as usize;

        let witness_merged_nv = num_vars + log_num_witness_polys;
        let selector_merged_nv = num_vars + log_num_selector_polys;

        let max_nv = max(witness_merged_nv + 1, selector_merged_nv);
        let max_points = max(
            // prod(x) has 5 points
            5,
            max(
                // selector points
                index.num_selector_columns(),
                // witness points + public input point + perm point
                index.num_witness_columns() + 2,
            ),
        );

        let supported_uni_degree = compute_qx_degree(max_nv, max_points);
        let supported_ml_degree = max_nv;

        // extract PCS prover and verifier keys from SRS
        let (pcs_prover_param, pcs_verifier_param) =
            PCS::trim(pcs_srs, supported_uni_degree, Some(supported_ml_degree))?;

        // build permutation oracles
        let permutation_oracle = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
            witness_merged_nv,
            &index.permutation,
        ));
        let perm_com = PCS::commit(&pcs_prover_param, &permutation_oracle)?;

        // build selector oracles and commit to it
        let selector_oracles: Vec<Rc<DenseMultilinearExtension<E::Fr>>> = index
            .selectors
            .iter()
            .map(|s| Rc::new(DenseMultilinearExtension::from(s)))
            .collect();

        let selector_merged = merge_polynomials(&selector_oracles)?;
        let selector_com = PCS::commit(&pcs_prover_param, &selector_merged)?;

        Ok((
            Self::ProvingKey {
                params: index.params.clone(),
                permutation_oracle: permutation_oracle.clone(),
                selector_oracles,
                selector_com: selector_com.clone(),
                pcs_param: pcs_prover_param,
            },
            Self::VerifyingKey {
                params: index.params.clone(),
                permutation_oracle,
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
    /// 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracle`, and
    /// obtain a PermCheckSubClaim.
    ///
    /// 4. Generate evaluations and corresponding proofs
    /// - 4.1. (deferred) batch opening prod(x) at
    ///   - [0, perm_check_point]
    ///   - [1, perm_check_point]
    ///   - [perm_check_point, 0]
    ///   - [perm_check_point, 1]
    ///   - [1,...1, 0]
    ///
    /// - 4.2. permutation check evaluations and proofs
    ///   - 4.2.1. (deferred) wi_poly(perm_check_point)
    ///
    /// - 4.3. zero check evaluations and proofs
    ///   - 4.3.1. (deferred) wi_poly(zero_check_point)
    ///   - 4.3.2. (deferred) selector_poly(zero_check_point)
    ///
    /// - 4.4. public input consistency checks
    ///   - pi_poly(r_pi) where r_pi is sampled from transcript
    ///
    /// - 5. deferred batch opening
    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[E::Fr],
        witnesses: &[WitnessColumn<E::Fr>],
    ) -> Result<Self::Proof, HyperPlonkErrors> {
        let start = start_timer!(|| "hyperplonk proving");
        let mut transcript = IOPTranscript::<E::Fr>::new(b"hyperplonk");

        prover_sanity_check(&pk.params, pub_input, witnesses)?;

        // witness assignment of length 2^n
        let num_vars = pk.params.num_variables();
        let log_num_witness_polys = log2(pk.params.num_witness_columns()) as usize;
        let log_num_selector_polys = log2(pk.params.num_selector_columns()) as usize;
        // number of variables in merged polynomial for Multilinear-KZG
        let merged_nv = num_vars + log_num_witness_polys;
        // online public input of length 2^\ell
        let ell = log2(pk.params.num_pub_input) as usize;

        // We use accumulators to store the polynomials and their eval points.
        // They are batch opened at a later stage.
        // This includes
        // - witnesses
        // - prod(x)
        // - selectors
        //
        // Accumulator for w_merged and its points
        let mut w_merged_pcs_acc = PcsAccumulator::<E, PCS>::new();
        // Accumulator for prod(x) and its points
        let mut prod_pcs_acc = PcsAccumulator::<E, PCS>::new();
        // Accumulator for prod(x) and its points
        let mut selector_pcs_acc = PcsAccumulator::<E, PCS>::new();

        let witness_polys: Vec<Rc<DenseMultilinearExtension<E::Fr>>> = witnesses
            .iter()
            .map(|w| Rc::new(DenseMultilinearExtension::from(w)))
            .collect();

        // =======================================================================
        // 1. Commit Witness polynomials `w_i(x)` and append commitment to
        // transcript
        // =======================================================================
        let step = start_timer!(|| "commit witnesses");
        let w_merged = merge_polynomials(&witness_polys)?;
        if w_merged.num_vars != merged_nv {
            return Err(HyperPlonkErrors::InvalidParameters(format!(
                "merged witness poly has a different num_vars ({}) from expected ({})",
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
        let step = start_timer!(|| "ZeroCheck on f");

        let fx = build_f(
            &pk.params.gate_func,
            pk.params.num_variables(),
            &pk.selector_oracles,
            &witness_polys,
        )?;

        let zero_check_proof = <Self as ZeroCheck<E::Fr>>::prove(&fx, &mut transcript)?;
        end_timer!(step);
        // =======================================================================
        // 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracle`, and
        // obtain a PermCheckSubClaim.
        // =======================================================================
        let step = start_timer!(|| "Permutation check on w_i(x)");

        let (perm_check_proof, prod_x) = <Self as PermutationCheck<E, PCS>>::prove(
            &pk.pcs_param,
            &w_merged,
            &w_merged,
            &pk.permutation_oracle,
            &mut transcript,
        )?;
        let perm_check_point = &perm_check_proof.zero_check_proof.point;

        end_timer!(step);
        // =======================================================================
        // 4. Generate evaluations and corresponding proofs
        // - 4.1. (deferred) batch opening prod(x) at
        //   - [0, perm_check_point]
        //   - [1, perm_check_point]
        //   - [perm_check_point, 0]
        //   - [perm_check_point, 1]
        //   - [1,...1, 0]
        //
        // - 4.2. permutation check evaluations and proofs
        //   - 4.2.1. (deferred) wi_poly(perm_check_point)
        //
        // - 4.3. zero check evaluations and proofs
        //   - 4.3.1. (deferred) wi_poly(zero_check_point)
        //   - 4.3.2. (deferred) selector_poly(zero_check_point)
        //
        // - 4.4. (deferred) public input consistency checks
        //   - pi_poly(r_pi) where r_pi is sampled from transcript
        // =======================================================================
        let step = start_timer!(|| "opening and evaluations");

        // 4.1 (deferred) open prod(0,x), prod(1, x), prod(x, 0), prod(x, 1)
        // perm_check_point
        // prod(0, x)
        let tmp_point1 = [perm_check_point.as_slice(), &[E::Fr::zero()]].concat();
        // prod(1, x)
        let tmp_point2 = [perm_check_point.as_slice(), &[E::Fr::one()]].concat();
        // prod(x, 0)
        let tmp_point3 = [&[E::Fr::zero()], perm_check_point.as_slice()].concat();
        // prod(x, 1)
        let tmp_point4 = [&[E::Fr::one()], perm_check_point.as_slice()].concat();
        // prod(1, ..., 1, 0)
        let tmp_point5 = [vec![E::Fr::zero()], vec![E::Fr::one(); merged_nv]].concat();

        prod_pcs_acc.insert_poly_and_points(&prod_x, &perm_check_proof.prod_x_comm, &tmp_point1);
        prod_pcs_acc.insert_poly_and_points(&prod_x, &perm_check_proof.prod_x_comm, &tmp_point2);
        prod_pcs_acc.insert_poly_and_points(&prod_x, &perm_check_proof.prod_x_comm, &tmp_point3);
        prod_pcs_acc.insert_poly_and_points(&prod_x, &perm_check_proof.prod_x_comm, &tmp_point4);
        prod_pcs_acc.insert_poly_and_points(&prod_x, &perm_check_proof.prod_x_comm, &tmp_point5);

        // 4.2  permutation check
        //   - 4.2.1. (deferred) wi_poly(perm_check_point)
        w_merged_pcs_acc.insert_poly_and_points(&w_merged, &w_merged_com, perm_check_point);

        // - 4.3. zero check evaluations and proofs
        //   - 4.3.1 (deferred) wi_poly(zero_check_point)
        for i in 0..witness_polys.len() {
            let tmp_point = gen_eval_point(i, log_num_witness_polys, &zero_check_proof.point);
            // Deferred opening zero check proof
            // w_merged_pcs_acc.insert_point(&tmp_point);

            w_merged_pcs_acc.insert_poly_and_points(&w_merged, &w_merged_com, &tmp_point);
        }

        //   - 4.3.2. (deferred) selector_poly(zero_check_point)
        let selector_merged = merge_polynomials(&pk.selector_oracles)?;
        // selector_pcs_acc.init_poly(selector_merged, pk.selector_com.clone())?;
        for i in 0..pk.selector_oracles.len() {
            let tmp_point = gen_eval_point(i, log_num_selector_polys, &zero_check_proof.point);
            // Deferred opening zero check proof
            selector_pcs_acc.insert_poly_and_points(&selector_merged, &pk.selector_com, &tmp_point);
        }

        // - 4.4. public input consistency checks
        //   - pi_poly(r_pi) where r_pi is sampled from transcript
        let r_pi = transcript.get_and_append_challenge_vectors(b"r_pi", ell)?;
        let tmp_point = [
            vec![E::Fr::zero(); num_vars - ell],
            r_pi,
            vec![E::Fr::zero(); log_num_witness_polys],
        ]
        .concat();

        w_merged_pcs_acc.insert_poly_and_points(&w_merged, &w_merged_com, &tmp_point);
        end_timer!(step);

        // =======================================================================
        // 5. deferred batch opening
        // =======================================================================
        let step = start_timer!(|| "deferred batch openings");
        let sub_step = start_timer!(|| "open witness");
        let w_merged_batch_opening = w_merged_pcs_acc.multi_open(&pk.pcs_param, &mut transcript)?;
        end_timer!(sub_step);

        let sub_step = start_timer!(|| "open prod(x)");
        // let (prod_batch_openings, prod_batch_evals) =
        // prod_pcs_acc.batch_open(&pk.pcs_param)?;
        let prod_batch_openings = prod_pcs_acc.multi_open(&pk.pcs_param, &mut transcript)?;
        end_timer!(sub_step);

        let sub_step = start_timer!(|| "open selector");
        let selector_batch_openings =
            selector_pcs_acc.multi_open(&pk.pcs_param, &mut transcript)?;
        // let (selector_batch_opening, selector_batch_evals) = selector_pcs_acc
        //     .batch_open_overlapped_points(&pk.pcs_param,
        // zero_check_proof.point.len())?;
        end_timer!(sub_step);
        end_timer!(step);
        end_timer!(start);

        Ok(HyperPlonkProof {
            // =======================================================================
            // witness related
            // =======================================================================
            /// PCS commit for witnesses
            w_merged_batch_opening,
            w_merged_com,
            // Batch opening for witness commitment
            // - PermCheck eval: 1 point
            // - ZeroCheck evals: #witness points
            // - public input eval: 1 point
            // w_merged_batch_opening,
            // Evaluations of Witness
            // - PermCheck eval: 1 point
            // - ZeroCheck evals: #witness points
            // - public input eval: 1 point
            // w_merged_batch_evals,
            // =======================================================================
            // prod(x) related
            // =======================================================================
            // prod(x)'s openings
            // - prod(0, x),
            // - prod(1, x),
            // - prod(x, 0),
            // - prod(x, 1),
            // - prod(1, ..., 1,0)
            prod_batch_openings,
            // prod(x)'s evaluations
            // - prod(0, x),
            // - prod(1, x),
            // - prod(x, 0),
            // - prod(x, 1),
            // - prod(1, ..., 1,0)
            // prod_batch_evals,
            // =======================================================================
            // selectors related
            // =======================================================================
            // PCS openings for selectors on zero check point
            selector_batch_openings,
            // Evaluates of selectors on zero check point
            // selector_batch_evals,
            // =======================================================================
            // IOP proofs
            // =======================================================================
            // the custom gate zerocheck proof
            zero_check_proof,
            // the permutation check proof for copy constraints
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
    /// 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracle`
    ///
    /// 3. check subclaim validity
    ///
    /// 4. Verify the opening against the commitment:
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
        // witness assignment of length 2^n
        let num_vars = vk.params.num_variables();
        let log_num_witness_polys = log2(vk.params.num_witness_columns()) as usize;
        // number of variables in merged polynomial for Multilinear-KZG
        let merged_nv = num_vars + log_num_witness_polys;

        //  online public input of length 2^\ell
        let ell = log2(vk.params.num_pub_input) as usize;

        let pi_poly = DenseMultilinearExtension::from_evaluations_slice(ell as usize, pub_input);

        // =======================================================================
        // 0. sanity checks
        // =======================================================================
        // public input length
        if pub_input.len() != vk.params.num_pub_input {
            return Err(HyperPlonkErrors::InvalidProver(format!(
                "Public input length is not correct: got {}, expect {}",
                pub_input.len(),
                1 << ell
            )));
        }
        if proof.selector_batch_openings.f_i_eval_at_point_i.len()
            != vk.params.num_selector_columns()
        {
            return Err(HyperPlonkErrors::InvalidVerifier(format!(
                "Selector length is not correct: got {}, expect {}",
                proof.selector_batch_openings.f_i_eval_at_point_i.len(),
                1 << vk.params.num_selector_columns()
            )));
        }
        if proof.w_merged_batch_opening.f_i_eval_at_point_i.len()
            != vk.params.num_witness_columns() + 2
        {
            return Err(HyperPlonkErrors::InvalidVerifier(format!(
                "Witness length is not correct: got {}, expect {}",
                proof.w_merged_batch_opening.f_i_eval_at_point_i.len() - 2,
                vk.params.num_witness_columns()
            )));
        }
        if proof.prod_batch_openings.f_i_eval_at_point_i.len() != 5 {
            return Err(HyperPlonkErrors::InvalidVerifier(format!(
                "the number of product polynomial evaluations is not correct: got {}, expect {}",
                proof.prod_batch_openings.f_i_eval_at_point_i.len(),
                5
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
        // Zero check and perm check have different AuxInfo
        let zero_check_aux_info = VPAuxInfo::<E::Fr> {
            max_degree: vk.params.gate_func.degree(),
            num_variables: num_vars,
            phantom: PhantomData::default(),
        };
        // push witness to transcript
        transcript.append_serializable_element(b"w", &proof.w_merged_com)?;

        let zero_check_sub_claim = <Self as ZeroCheck<E::Fr>>::verify(
            &proof.zero_check_proof,
            &zero_check_aux_info,
            &mut transcript,
        )?;

        let zero_check_point = &zero_check_sub_claim.point;

        // check zero check subclaim
        let f_eval = eval_f(
            &vk.params.gate_func,
            &&proof.selector_batch_openings.f_i_eval_at_point_i[..vk.params.num_selector_columns()],
            &proof.w_merged_batch_opening.f_i_eval_at_point_i[1..],
        )?;
        if f_eval != zero_check_sub_claim.expected_evaluation {
            return Err(HyperPlonkErrors::InvalidProof(
                "zero check evaluation failed".to_string(),
            ));
        }

        end_timer!(step);
        // =======================================================================
        // 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracle`
        // =======================================================================
        let step = start_timer!(|| "verify permutation check");

        // Zero check and perm check have different AuxInfo
        let perm_check_aux_info = VPAuxInfo::<E::Fr> {
            // Prod(x) has a max degree of 2
            max_degree: 2,
            // degree of merged poly
            num_variables: merged_nv,
            phantom: PhantomData::default(),
        };
        let perm_check_sub_claim = <Self as PermutationCheck<E, PCS>>::verify(
            &proof.perm_check_proof,
            &perm_check_aux_info,
            &mut transcript,
        )?;

        let perm_check_point = &perm_check_sub_claim
            .product_check_sub_claim
            .zero_check_sub_claim
            .point;

        let alpha = perm_check_sub_claim.product_check_sub_claim.alpha;
        let (beta, gamma) = perm_check_sub_claim.challenges;

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

        // we evaluate MLE directly instead of using s_id/s_perm PCS verify
        // Verification takes n pairings while evaluate takes 2^n field ops.
        let s_id = identity_permutation_mle::<E::Fr>(perm_check_point.len());
        let s_id_eval = evaluate_opt(&s_id, perm_check_point);
        let s_perm_eval = evaluate_opt(&vk.permutation_oracle, perm_check_point);

        let q_x_rec = proof.prod_batch_openings.f_i_eval_at_point_i[1]
            - proof.prod_batch_openings.f_i_eval_at_point_i[2]
                * proof.prod_batch_openings.f_i_eval_at_point_i[3]
            + alpha
                * ((proof.w_merged_batch_opening.f_i_eval_at_point_i[0]
                    + beta * s_perm_eval
                    + gamma)
                    * proof.prod_batch_openings.f_i_eval_at_point_i[0]
                    - (proof.w_merged_batch_opening.f_i_eval_at_point_i[0]
                        + beta * s_id_eval
                        + gamma));

        if q_x_rec
            != perm_check_sub_claim
                .product_check_sub_claim
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
        // 3.1 open prod(x)' evaluations
        // =======================================================================
        let prod_final_query = perm_check_sub_claim.product_check_sub_claim.final_query;
        let prod_points = [
            [perm_check_point.as_slice(), &[E::Fr::zero()]].concat(),
            [perm_check_point.as_slice(), &[E::Fr::one()]].concat(),
            [&[E::Fr::zero()], perm_check_point.as_slice()].concat(),
            [&[E::Fr::one()], perm_check_point.as_slice()].concat(),
            prod_final_query.0,
        ];

        // if !PCS::batch_verify_single_poly(
        //     &vk.pcs_param,
        //     &proof.perm_check_proof.prod_x_comm,
        //     &points,
        //     &proof.prod_batch_evals,
        //     &proof.prod_batch_openings,
        // )? {
        //     return Err(HyperPlonkErrors::InvalidProof(
        //         "prod(0, x) pcs verification failed".to_string(),
        //     ));
        // }

        // =======================================================================
        // 3.2 open selectors' evaluations
        // =======================================================================
        let log_num_selector_polys = log2(vk.params.num_selector_columns()) as usize;
        let mut selector_points = vec![];
        for i in 0..vk.params.num_selector_columns() {
            let tmp_point =
                gen_eval_point(i, log_num_selector_polys, &proof.zero_check_proof.point);
            selector_points.push(tmp_point);
        }

        // if proof.se.len() == 2 {
        //     if !PCS::batch_verify_single_poly(
        //         &vk.pcs_param,
        //         &vk.selector_com,
        //         &points,
        //         &proof.selector_batch_evals,
        //         &proof.selector_batch_opening,
        //     )? {
        //         return Err(HyperPlonkErrors::InvalidProof(
        //             "selector pcs verification failed".to_string(),
        //         ));
        //     }
        // } else if !PCS::batch_verify_single_poly_overlap_points(
        //     &vk.pcs_param,
        //     &vk.selector_com,
        //     &points,
        //     &proof.selector_batch_evals,
        //     &proof.selector_batch_opening,
        //     proof.zero_check_proof.point.len(),
        // )? {
        //     return Err(HyperPlonkErrors::InvalidProof(
        //         "selector pcs verification failed".to_string(),
        //     ));
        // }

        // =======================================================================
        // 3.2 open witnesses' evaluations
        // =======================================================================
        let mut r_pi = transcript.get_and_append_challenge_vectors(b"r_pi", ell)?;
        let pi_eval = evaluate_opt(&pi_poly, &r_pi);
        assert_eq!(
            pi_eval,
            proof.w_merged_batch_opening.f_i_eval_at_point_i
                [proof.w_merged_batch_opening.f_i_eval_at_point_i.len() - 1]
        );

        r_pi = [
            vec![E::Fr::zero(); num_vars - ell],
            r_pi,
            vec![E::Fr::zero(); log_num_witness_polys],
        ]
        .concat();

        let mut w_merged_points = vec![perm_check_point.clone()];

        for i in 0..proof.w_merged_batch_opening.f_i_eval_at_point_i.len() - 2 {
            w_merged_points.push(gen_eval_point(i, log_num_witness_polys, zero_check_point))
        }
        w_merged_points.push(r_pi);
        // if !PCS::batch_verify_single_poly(
        //     &vk.pcs_param,
        //     &proof.w_merged_com,
        //     &points,
        //     &proof.w_merged_batch_evals,
        //     &proof.w_merged_batch_opening,
        // )? {
        //     return Err(HyperPlonkErrors::InvalidProof(
        //         "witness for permutation check pcs verification failed".to_string(),
        //     ));
        // }

        let res = batch_verify_internal(
            &vk.pcs_param,
            vec![proof.w_merged_com; vk.params.num_witness_columns() + 2].as_ref(),
            &w_merged_points,
            &proof.w_merged_batch_opening,
            &mut transcript,
        )?;
        assert!(res);
        let res = batch_verify_internal(
            &vk.pcs_param,
            vec![proof.perm_check_proof.prod_x_comm; 5].as_ref(),
            &prod_points,
            &proof.prod_batch_openings,
            &mut transcript,
        )?;
        assert!(res);
        let res = batch_verify_internal(
            &vk.pcs_param,
            vec![vk.selector_com; vk.params.num_selector_columns()].as_ref(),
            &selector_points,
            &proof.selector_batch_openings,
            &mut transcript,
        )?;

        assert!(res);
        end_timer!(step);
        end_timer!(start);
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        custom_gate::CustomizedGates, selectors::SelectorColumn, structs::HyperPlonkParams,
        witness::WitnessColumn,
    };
    use arithmetic::random_permutation_mle;
    use ark_bls12_381::Bls12_381;
    use ark_std::test_rng;
    use pcs::prelude::MultilinearKzgPCS;

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
        test_hyperplonk_helper::<Bls12_381>(gates)
    }

    fn test_hyperplonk_helper<E: PairingEngine>(
        gate_func: CustomizedGates,
    ) -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs = MultilinearKzgPCS::<E>::gen_srs_for_testing(&mut rng, 16)?;

        let num_constraints = 4;
        let num_pub_input = 4;
        let nv = log2(num_constraints) as usize;
        let merged_nv = nv + log2(gate_func.num_witness_columns()) as usize;

        // generate index
        let params = HyperPlonkParams {
            num_constraints,
            num_pub_input,
            gate_func,
        };
        let permutation = identity_permutation_mle(merged_nv).evaluations.clone();
        let q1 = SelectorColumn(vec![E::Fr::one(), E::Fr::one(), E::Fr::one(), E::Fr::one()]);
        let index = HyperPlonkIndex {
            params,
            permutation,
            selectors: vec![q1],
        };

        // generate pk and vks
        let (pk, vk) = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, MultilinearKzgPCS<E>>>::preprocess(
            &index, &pcs_srs,
        )?;

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

        // generate a proof and verify
        let proof = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, MultilinearKzgPCS<E>>>::prove(
            &pk,
            &pi.0,
            &[w1.clone(), w2.clone()],
        )?;

        let _verify = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, MultilinearKzgPCS<E>>>::verify(
            &vk, &pi.0, &proof,
        )?;

        // bad path 1: wrong permutation
        let rand_perm: Vec<E::Fr> = random_permutation_mle(merged_nv, &mut rng)
            .evaluations
            .clone();
        let mut bad_index = index.clone();
        bad_index.permutation = rand_perm;
        // generate pk and vks
        let (_, bad_vk) = <PolyIOP<E::Fr> as HyperPlonkSNARK<E, MultilinearKzgPCS<E>>>::preprocess(
            &bad_index, &pcs_srs,
        )?;
        assert!(
            <PolyIOP<E::Fr> as HyperPlonkSNARK<E, MultilinearKzgPCS<E>>>::verify(
                &bad_vk, &pi.0, &proof,
            )
            .is_err()
        );

        // bad path 2: wrong witness
        let mut w1_bad = w1.clone();
        w1_bad.0[0] = E::Fr::one();
        assert!(
            <PolyIOP<E::Fr> as HyperPlonkSNARK<E, MultilinearKzgPCS<E>>>::prove(
                &pk,
                &pi.0,
                &[w1_bad, w2],
            )
            .is_err()
        );

        Ok(())
    }
}
