// Copyright (c) 2022 Espresso Systems (espressosys.com)
// This file is part of the Jellyfish library.

// You should have received a copy of the MIT License
// along with the Jellyfish library. If not, see <https://mit-license.org/>.

use crate::{
    multilinear_kzg::{
        open_internal,
        srs::{MultilinearProverParam, MultilinearVerifierParam},
        util::{build_l, compute_w_circ_l, get_uni_domain},
        verify_internal, MultilinearKzgBatchProof,
    },
    prelude::{Commitment, UnivariateProverParam, UnivariateVerifierParam},
    univariate_kzg::UnivariateKzgPCS,
    PCSError, PolynomialCommitmentScheme,
};
use ark_ec::PairingEngine;
use ark_poly::{DenseMultilinearExtension, EvaluationDomain, MultilinearExtension, Polynomial};
use ark_std::{end_timer, format, rc::Rc, start_timer, string::ToString, vec, vec::Vec};
use transcript::IOPTranscript;

/// Input
/// - the prover parameters for univariate KZG,
/// - the prover parameters for multilinear KZG,
/// - a single MLE,
/// - a commitment to the MLE
/// - and a list of points,
/// compute a multi-opening for this polynomial.
///
/// For simplicity, this API requires each MLE to have only one point. If
/// the caller wish to use more than one points per MLE, it should be
/// handled at the caller layer.
///
///
/// Returns the proof, consists of
/// - the multilinear KZG opening
/// - the univariate KZG commitment to q(x)
/// - the openings and evaluations of q(x) at omega^i and r
///
/// Steps:
/// 1. build `l(points)` which is a list of univariate polynomials that goes
/// through the points
/// 3. build `q(x)` which is a univariate polynomial `W circ l`
/// 4. commit to q(x) and sample r from transcript
/// transcript contains: w commitment, points, q(x)'s commitment
/// 5. build q(omega^i) and their openings
/// 6. build q(r) and its opening
/// 7. get a point `p := l(r)`
/// 8. output an opening of `w` over point `p`
/// 9. output `w(p)`
pub(crate) fn multi_open_same_poly_internal<E: PairingEngine>(
    uni_prover_param: &UnivariateProverParam<E::G1Affine>,
    ml_prover_param: &MultilinearProverParam<E>,
    polynomial: &Rc<DenseMultilinearExtension<E::Fr>>,
    commitment: &Commitment<E>,
    points: &[Vec<E::Fr>],
) -> Result<(MultilinearKzgBatchProof<E>, Vec<E::Fr>), PCSError> {
    let open_timer = start_timer!(|| "multi open");

    // ===================================
    // Sanity checks on inputs
    // ===================================
    let points_len = points.len();
    if points_len == 0 {
        return Err(PCSError::InvalidParameters("points is empty".to_string()));
    }

    let num_var = polynomial.num_vars();
    for point in points.iter() {
        if point.len() != num_var {
            return Err(PCSError::InvalidParameters(
                "points do not have same num_vars".to_string(),
            ));
        }
    }

    let domain = get_uni_domain::<E::Fr>(points_len)?;

    // 1. build `l(points)` which is a list of univariate polynomials that goes
    // through the points
    let uni_polys = build_l(points, &domain)?;

    // 3. build `q(x)` which is a univariate polynomial `W circ l`
    let q_x = compute_w_circ_l(&polynomial, &uni_polys, points.len())?;

    // 4. commit to q(x) and sample r from transcript
    // transcript contains: w commitment, points, q(x)'s commitment
    let mut transcript = IOPTranscript::new(b"ml kzg");
    transcript.append_serializable_element(b"w", commitment)?;
    for point in points {
        // println!("points proving {:?}", points);
        transcript.append_serializable_element(b"w", point)?;
    }

    let q_x_commit = UnivariateKzgPCS::<E>::commit(uni_prover_param, &q_x)?;
    transcript.append_serializable_element(b"q(x)", &q_x_commit)?;
    let r = transcript.get_and_append_challenge(b"r")?;
    // 5. build q(omega^i) and their openings
    let mut q_x_opens = vec![];
    let mut q_x_evals = vec![];
    for i in 0..points_len {
        let (q_x_open, q_x_eval) =
            UnivariateKzgPCS::<E>::open(uni_prover_param, &q_x, &domain.element(i))?;
        q_x_opens.push(q_x_open);
        q_x_evals.push(q_x_eval);

        #[cfg(feature = "extensive_sanity_checks")]
        {
            // sanity check
            let point: Vec<E::Fr> = uni_polys
                .iter()
                .map(|poly| poly.evaluate(&domain.element(i)))
                .collect();
            let mle_eval = polynomial.evaluate(&point).unwrap();
            if mle_eval != q_x_eval {
                return Err(PCSError::InvalidProver(
                    "Q(omega) does not match W(l(omega))".to_string(),
                ));
            }
        }
    }

    // 6. build q(r) and its opening
    let (q_x_open, q_r_value) = UnivariateKzgPCS::<E>::open(uni_prover_param, &q_x, &r)?;
    q_x_opens.push(q_x_open);
    q_x_evals.push(q_r_value);

    // 7. get a point `p := l(r)`
    let point: Vec<E::Fr> = uni_polys.iter().map(|poly| poly.evaluate(&r)).collect();
    // 8. output an opening of `w` over point `p`
    let (mle_opening, mle_eval) = open_internal(ml_prover_param, &polynomial, &point)?;

    // 9. output value that is `w` evaluated at `p` (which should match `q(r)`)
    if mle_eval != q_r_value {
        return Err(PCSError::InvalidProver(
            "Q(r) does not match W(l(r))".to_string(),
        ));
    }
    end_timer!(open_timer);
    Ok((
        MultilinearKzgBatchProof {
            proof: mle_opening,
            q_x_commit,
            q_x_opens,
        },
        q_x_evals,
    ))
}

/// Verifies that the `multi_commitment` is a valid commitment
/// to a list of MLEs for the given openings and evaluations in
/// the batch_proof.
///
/// steps:
///
/// 1. push w, points and q_com into transcript
/// 2. sample `r` from transcript
/// 3. check `q(r) == batch_proof.q_x_value.last` and
/// `q(omega^i) == batch_proof.q_x_value[i]`
/// 4. build `l(points)` which is a list of univariate
/// polynomials that goes through the points
/// 5. get a point `p := l(r)`
/// 6. verifies `p` is valid against multilinear KZG proof
#[allow(dead_code)]
pub(crate) fn batch_verify_same_poly_internal<E: PairingEngine>(
    uni_verifier_param: &UnivariateVerifierParam<E>,
    ml_verifier_param: &MultilinearVerifierParam<E>,
    multi_commitment: &Commitment<E>,
    points: &[Vec<E::Fr>],
    values: &[E::Fr],
    batch_proof: &MultilinearKzgBatchProof<E>,
) -> Result<bool, PCSError> {
    let verify_timer = start_timer!(|| "batch verify");

    // ===================================
    // Sanity checks on inputs
    // ===================================
    let points_len = points.len();
    if points_len == 0 {
        return Err(PCSError::InvalidParameters("points is empty".to_string()));
    }

    // add one here because we also have q(r) and its opening
    if points_len + 1 != batch_proof.q_x_opens.len() {
        return Err(PCSError::InvalidParameters(format!(
            "openings length {} does not match point length {}",
            points_len + 1,
            batch_proof.q_x_opens.len()
        )));
    }

    if points_len + 1 != values.len() {
        return Err(PCSError::InvalidParameters(format!(
            "values length {} does not match point length {}",
            points_len + 1,
            values.len()
        )));
    }

    let num_var = points[0].len();
    for point in points.iter().skip(1) {
        if point.len() != num_var {
            return Err(PCSError::InvalidParameters(format!(
                "points do not have same num_vars ({} vs {})",
                point.len(),
                num_var,
            )));
        }
    }

    let domain = get_uni_domain::<E::Fr>(points_len)?;
    // 1. push w, points and q_com into transcript
    let mut transcript = IOPTranscript::new(b"ml kzg");
    transcript.append_serializable_element(b"w", multi_commitment)?;

    for point in points {
        transcript.append_serializable_element(b"w", point)?;
    }

    transcript.append_serializable_element(b"q(x)", &batch_proof.q_x_commit)?;
    // 2. sample `r` from transcript
    let r = transcript.get_and_append_challenge(b"r")?;
    // 3. check `q(r) == batch_proof.q_x_value.last` and `q(omega^i) =
    // batch_proof.q_x_value[i]`
    for (i, value) in values.iter().enumerate().take(points_len) {
        if !UnivariateKzgPCS::verify(
            uni_verifier_param,
            &batch_proof.q_x_commit,
            &domain.element(i),
            value,
            &batch_proof.q_x_opens[i],
        )? {
            #[cfg(debug_assertion)]
            println!("q(omega^{}) verification failed", i);
            return Ok(false);
        }
    }

    if !UnivariateKzgPCS::verify(
        uni_verifier_param,
        &batch_proof.q_x_commit,
        &r,
        &values[points_len],
        &batch_proof.q_x_opens[points_len],
    )? {
        #[cfg(debug_assertion)]
        println!("q(r) verification failed");
        return Ok(false);
    }
    // 4. build `l(points)` which is a list of univariate polynomials that goes
    // through the points
    let uni_polys = build_l(points, &domain)?;

    // 5. get a point `p := l(r)`
    let point: Vec<E::Fr> = uni_polys.iter().map(|x| x.evaluate(&r)).collect();
    // 6. verifies `p` is valid against multilinear KZG proof
    let res = verify_internal(
        ml_verifier_param,
        multi_commitment,
        &point,
        &values[points_len],
        &batch_proof.proof,
    )?;
    #[cfg(debug_assertion)]
    if !res {
        println!("multilinear KZG verification failed");
    }

    end_timer!(verify_timer);

    Ok(res)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        multilinear_kzg::{
            srs::MultilinearUniversalParams,
            util::{compute_qx_degree, generate_evaluations_single_poly},
            MultilinearKzgPCS,
        },
        prelude::UnivariateUniversalParams,
        StructuredReferenceString,
    };
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
    type Fr = <E as PairingEngine>::Fr;

    fn test_same_poly_multi_open_internal_helper<R: RngCore + CryptoRng>(
        uni_params: &UnivariateUniversalParams<E>,
        ml_params: &MultilinearUniversalParams<E>,
        poly: &Rc<DenseMultilinearExtension<Fr>>,
        point_len: usize,
        rng: &mut R,
    ) -> Result<(), PCSError> {
        let nv = poly.num_vars;
        let qx_degree = compute_qx_degree(nv, point_len);
        let padded_qx_degree = 1usize << log2(qx_degree);

        let (uni_ck, uni_vk) = uni_params.trim(padded_qx_degree)?;
        let (ml_ck, ml_vk) = ml_params.trim(nv)?;

        let mut points = Vec::new();
        let mut eval = Vec::new();
        for _ in 0..point_len {
            let point = (0..nv).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();
            eval.push(poly.evaluate(&point).unwrap());
            points.push(point);
        }

        let evals = generate_evaluations_single_poly(poly, &points)?;
        let com = MultilinearKzgPCS::commit(&(ml_ck.clone(), uni_ck.clone()), poly)?;
        let (batch_proof, evaluations) =
            multi_open_same_poly_internal(&uni_ck, &ml_ck, poly, &com, &points)?;

        for (a, b) in evals.iter().zip(evaluations.iter()) {
            assert_eq!(a, b)
        }

        // good path
        assert!(batch_verify_same_poly_internal(
            &uni_vk,
            &ml_vk,
            &com,
            &points,
            &evaluations,
            &batch_proof,
        )?);

        Ok(())
    }

    #[test]
    fn test_same_poly_multi_open_internal() -> Result<(), PCSError> {
        let mut rng = test_rng();

        let uni_params =
            UnivariateUniversalParams::<E>::gen_srs_for_testing(&mut rng, 1usize << 15)?;
        let ml_params = MultilinearUniversalParams::<E>::gen_srs_for_testing(&mut rng, 15)?;
        for nv in 1..10 {
            for point_len in 1..10 {
                // normal polynomials
                let polys1 = Rc::new(DenseMultilinearExtension::rand(nv, &mut rng));
                test_same_poly_multi_open_internal_helper(
                    &uni_params,
                    &ml_params,
                    &polys1,
                    point_len,
                    &mut rng,
                )?;
            }
        }
        Ok(())
    }
}
