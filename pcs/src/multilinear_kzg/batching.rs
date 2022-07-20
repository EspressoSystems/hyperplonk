use super::{
    srs::{MultilinearProverParam, MultilinearVerifierParam},
    util::{build_l, compute_w_circ_l, merge_polynomials},
    BatchProof,
};
use crate::{
    multilinear_kzg::util::get_uni_domain,
    prelude::{Commitment, KZGMultilinearPCS, UnivariateProverParam, UnivariateVerifierParam},
    univariate_kzg::KZGUnivariatePCS,
    PCSErrors, PolynomialCommitmentScheme,
};
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, EvaluationDomain, MultilinearExtension, Polynomial};
use ark_std::{end_timer, start_timer, vec::Vec};
use poly_iop::IOPTranscript;

/// Input
/// - the prover parameters for univariate KZG,
/// - the prover parameters for multilinear KZG,
/// - a list of MLEs,
/// - a commitment to all MLEs
/// - and a same number of points,
/// compute a multi-opening for all the polynomials.
///
/// For simplicity, this API requires each MLE to have only one point. If
/// the caller wish to use more than one points per MLE, it should be
/// handled at the caller layer.
///
/// Returns an error if the lengths do not match.
///
/// Returns:
/// - the proof,
/// - q(x), which is a univariate polynomial `w circ l` where `w` is the merged
///   MLE, and `l` is a list of polynomials that go through all the points.
///   TODO: change this field to a commitment to `q(x)`
/// - and a value which is `w` evaluated at `p:= l(r)` from some `r` from the
///   transcript.
///
/// Steps:
/// 1. build `l(points)` which is a list of univariate polynomials that goes
/// through the points
/// 2. build MLE `w` which is the merge of all MLEs.
/// 3. build `q(x)` which is a univariate polynomial `W circ l`
/// 4. output `q(x)`'
/// transcript contains: w commitment, points, q(x)'s commitment
/// 5. sample `r` from transcript
/// 6. get a point `p := l(r)`
/// 7. output an opening of `w` over point `p`
/// 8. output `w(p)`
fn multi_open_internal<E: PairingEngine>(
    uni_prover_param: &UnivariateProverParam<E::G1Affine>,
    ml_prover_param: &MultilinearProverParam<E>,
    polynomials: &[DenseMultilinearExtension<E::Fr>],
    multi_commitment: &Commitment<E>,
    points: &[Vec<E::Fr>],
) -> Result<BatchProof<E>, PCSErrors> {
    let open_timer = start_timer!(|| "multi open");

    // ===================================
    // Sanity checks on inputs
    // ===================================
    let uni_poly_degree = points.len();
    if uni_poly_degree == 0 {
        return Err(PCSErrors::InvalidParameters("points is empty".to_string()));
    }

    if uni_poly_degree != polynomials.len() {
        return Err(PCSErrors::InvalidParameters(
            "polynomial length does not match point length".to_string(),
        ));
    }

    let num_var = polynomials[0].num_vars();
    for poly in polynomials.iter().skip(1) {
        if poly.num_vars() != num_var {
            return Err(PCSErrors::InvalidParameters(
                "polynomials do not have same num_vars".to_string(),
            ));
        }
    }
    for point in points.iter() {
        if point.len() != num_var {
            return Err(PCSErrors::InvalidParameters(
                "points do not have same num_vars".to_string(),
            ));
        }
    }

    let domain = get_uni_domain::<E::Fr>(uni_poly_degree)?;

    let mut transcript = IOPTranscript::new(b"ml kzg");
    transcript.append_serializable_element(b"w", multi_commitment)?;
    for point in points {
        transcript.append_serializable_element(b"w", point)?;
    }

    // 1. build `l(points)` which is a list of univariate polynomials that goes
    // through the points
    let uni_polys = build_l(num_var, points, &domain)?;

    // 2. build MLE `w` which is the merge of all MLEs.
    let merge_poly = merge_polynomials(polynomials)?;

    // 3. build `q(x)` which is a univariate polynomial `W circ l`
    let q_x = compute_w_circ_l(&merge_poly, &uni_polys)?;

    // 4. output `q(x)`' and put it into transcript
    let mut q_x_opens = vec![];
    let mut q_x_evals = vec![];
    for i in 0..uni_poly_degree {
        let q_x_eval = q_x.evaluate(&domain.element(i));
        let q_x_open = KZGUnivariatePCS::<E>::open(uni_prover_param, &q_x, &domain.element(i))?;
        q_x_opens.push(q_x_open);
        q_x_evals.push(q_x_eval);

        // sanity check
        let point: Vec<E::Fr> = uni_polys
            .iter()
            .rev()
            .map(|poly| poly.evaluate(&domain.element(i)))
            .collect();
        let mle_eval = merge_poly.evaluate(&point).unwrap();
        if mle_eval != q_x_eval {
            return Err(PCSErrors::InvalidProver(
                "Q(omega) does not match W(l(omega))".to_string(),
            ));
        }
    }

    let q_x_commit = KZGUnivariatePCS::<E>::commit(uni_prover_param, &q_x)?;
    transcript.append_serializable_element(b"q(x)", &q_x_commit)?;
    // 5. sample `r` from transcript
    let r = transcript.get_and_append_challenge(b"r")?;
    let q_x_open = KZGUnivariatePCS::<E>::open(uni_prover_param, &q_x, &r)?;
    q_x_opens.push(q_x_open);

    // 6. get a point `p := l(r)`
    let point: Vec<E::Fr> = uni_polys
        .iter()
        .rev()
        .map(|poly| poly.evaluate(&r))
        .collect();

    // 7. output an opening of `w` over point `p`
    let opening = KZGMultilinearPCS::<E>::open(ml_prover_param, &merge_poly, &point)?;

    // 8. output value that is `w` evaluated at `p` (which should match `q(r)`)
    let value = merge_poly.evaluate(&point).unwrap();
    let value2 = q_x.evaluate(&r);

    if value != value2 {
        return Err(PCSErrors::InvalidProver(
            "Q(r) does not match W(l(r))".to_string(),
        ));
    }
    q_x_evals.push(value);
    end_timer!(open_timer);

    Ok(BatchProof {
        proof: opening,
        q_x_commit,
        q_x_opens,
        q_x_evals,
    })
}

/// Verifies that `value` is the evaluation at `x_i` of the polynomial
/// `poly_i` committed inside `comm`.
/// steps:
///
/// 1. put `q(x)`'s evaluations over `(1, omega,...)` into transcript
/// 2. sample `r` from transcript
/// 3. check `q(r) == value`
/// 4. build `l(points)` which is a list of univariate polynomials that goes
/// through the points
/// 5. get a point `p := l(r)`
/// 6. verifies `p` is verifies against proof
fn batch_verify_internal<E: PairingEngine>(
    uni_verifier_param: &UnivariateVerifierParam<E>,
    ml_verifier_param: &MultilinearVerifierParam<E>,
    multi_commitment: &Commitment<E>,
    points: &[Vec<E::Fr>],
    _values: &[E::Fr],
    batch_proof: &BatchProof<E>,
) -> Result<bool, PCSErrors> {
    let verify_timer = start_timer!(|| "batch verify");

    // ===================================
    // Sanity checks on inputs
    // ===================================
    let uni_poly_degree = points.len();
    if uni_poly_degree == 0 {
        return Err(PCSErrors::InvalidParameters("points is empty".to_string()));
    }

    if uni_poly_degree+1 != batch_proof.q_x_opens.len() {
        return Err(PCSErrors::InvalidParameters(
            "openings length does not match point length".to_string(),
        ));
    }

    if uni_poly_degree+1 != batch_proof.q_x_evals.len() {
        return Err(PCSErrors::InvalidParameters(
            "values length does not match point length".to_string(),
        ));
    }

    let domain = get_uni_domain::<E::Fr>(uni_poly_degree)?;

    let mut transcript = IOPTranscript::new(b"ml kzg");
    transcript.append_serializable_element(b"w", multi_commitment)?;
    for point in points {
        transcript.append_serializable_element(b"w", point)?;
    }

    let num_var = points[0].len();

    for point in points.iter().skip(1) {
        if point.len() != num_var {
            return Err(PCSErrors::InvalidParameters(format!(
                "points do not have same num_vars ({} vs {})",
                point.len(),
                num_var,
            )));
        }
    }

    // 1. put `q(x)`'s evaluations over `(1, omega,...)` into transcript
    transcript.append_serializable_element(b"q(x)", &batch_proof.q_x_commit)?;

    // 2. sample `r` from transcript

    let r = transcript.get_and_append_challenge(b"r")?;

    // 3. check `q(r) == proof.value` and `q(omega^i) = proof.values[i]`
    for i in 0..uni_poly_degree {
        if !KZGUnivariatePCS::verify(
            uni_verifier_param,
            &batch_proof.q_x_commit,
            &domain.element(i),
            &batch_proof.q_x_evals[i],
            &batch_proof.q_x_opens[i],
        )? {
            println!("{}-th verification failed", i);
            return Ok(false);
        }
    }

    if !KZGUnivariatePCS::verify(
        uni_verifier_param,
        &batch_proof.q_x_commit,
        &r,
        &batch_proof.q_x_evals[uni_poly_degree],
        &batch_proof.q_x_opens[uni_poly_degree],
    )? {
        return Ok(false);
    }

    // 4. build `l(points)` which is a list of univariate polynomials that goes
    // through the points
    let uni_polys = build_l(num_var, points, &domain)?;

    // 5. get a point `p := l(r)`
    let point: Vec<E::Fr> = uni_polys.iter().rev().map(|x| x.evaluate(&r)).collect();

    // 6. verifies `p` is verifies against proof
    let res = KZGMultilinearPCS::verify(
        ml_verifier_param,
        multi_commitment,
        &point,
        &batch_proof.q_x_evals[uni_poly_degree],
        &batch_proof.proof,
    );
    end_timer!(verify_timer);

    res
}

#[cfg(test)]
mod tests {
    use super::{
        super::{util::get_batched_nv, *},
        *,
    };
    use crate::{
        multilinear_kzg::util::{compute_uni_degree, generate_values},
        prelude::UnivariateUniversalParams,
        univariate_kzg::KZGUnivariatePCS,
        StructuredReferenceString,
    };
    use ark_bls12_381::Bls12_381 as E;
    use ark_ec::PairingEngine;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::{rand::RngCore, test_rng, vec::Vec, UniformRand};
    type Fr = <E as PairingEngine>::Fr;

    fn test_multi_commit_helper<R: RngCore>(
        uni_params: &UnivariateUniversalParams<E>,
        ml_params: &MultilinearUniversalParams<E>,
        polys: &[DenseMultilinearExtension<Fr>],
        rng: &mut R,
    ) -> Result<(), PCSErrors> {
        let nv = get_batched_nv(polys[0].num_vars(), polys.len());
        let uni_degree = compute_uni_degree(polys[0].num_vars(), polys.len());

        let (uni_ck, uni_vk) = uni_params.trim(uni_degree)?;
        let (ml_ck, ml_vk) = ml_params.trim(nv)?;

        let mut points = Vec::new();
        for poly in polys.iter() {
            let point = (0..poly.num_vars())
                .map(|_| Fr::rand(rng))
                .collect::<Vec<Fr>>();
            points.push(point);
        }

        let evals = generate_values(polys, &points)?;

        let com = KZGMultilinearPCS::multi_commit(&ml_ck, polys)?;
        let batch_proof = multi_open_internal(&uni_ck, &ml_ck, polys, &com, &points)?;

        for (a, b) in evals.iter().zip(batch_proof.q_x_evals.iter()) {
            assert_eq!(a, b)
        }

        // good path
        assert!(batch_verify_internal(
            &uni_vk,
            &ml_vk,
            &com,
            &points,
            &evals,
            &batch_proof,
        )?);

        // bad commitment
        assert!(!batch_verify_internal(
            &uni_vk,
            &ml_vk,
            &Commitment {
                commitment: <E as PairingEngine>::G1Affine::default()
            },
            &points,
            &evals,
            &batch_proof,
        )?);

        // bad points
        assert!(
            batch_verify_internal(&uni_vk, &ml_vk, &com, &points[1..], &evals, &batch_proof,)
                .is_err()
        );

        // bad proof
        assert!(batch_verify_internal(
            &uni_vk,
            &ml_vk,
            &com,
            &points,
            &evals,
            &BatchProof {
                proof: Proof { proofs: Vec::new() },
                q_x_commit: Commitment {
                    commitment: <E as PairingEngine>::G1Affine::default()
                },
                q_x_opens: vec![],
                q_x_evals: vec![],
            },
        )
        .is_err());

        println!("2");
        // bad value
        assert!(batch_verify_internal(&uni_vk, &ml_vk, &com, &points, &[], &batch_proof,).is_err());

        println!("3");
        // bad q(x) commit
        let mut wrong_proof = batch_proof.clone();
        wrong_proof.q_x_commit = Commitment {
            commitment: <E as PairingEngine>::G1Affine::default(),
        };
        assert!(!batch_verify_internal(
            &uni_vk,
            &ml_vk,
            &com,
            &points,
            &evals,
            &wrong_proof,
        )?);

        Ok(())
    }

    #[test]
    fn test_multi_commit_internal() -> Result<(), PCSErrors> {
        let mut rng = test_rng();

        let uni_params = KZGUnivariatePCS::<E>::gen_srs_for_testing(&mut rng, 15)?;
        let ml_params = KZGMultilinearPCS::<E>::gen_srs_for_testing(&mut rng, 15)?;

        // normal polynomials
        let polys1: Vec<_> = (0..2)
            .map(|_| DenseMultilinearExtension::rand(4, &mut rng))
            .collect();
        test_multi_commit_helper(&uni_params, &ml_params, &polys1, &mut rng)?;

        // single-variate polynomials
        let polys1: Vec<_> = (0..5)
            .map(|_| DenseMultilinearExtension::rand(1, &mut rng))
            .collect();
        test_multi_commit_helper(&uni_params, &ml_params, &polys1, &mut rng)?;

        Ok(())
    }
}
