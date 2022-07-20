use super::{
    srs::{MultilinearProverParam, MultilinearUniversalParams, MultilinearVerifierParam},
    util::{build_l, compute_w_circ_l, merge_polynomials},
    BatchProof,
};
use crate::{
    prelude::{Commitment, KZGMultilinearPCS, UnivariateProverParam, UnivariateVerifierParam},
    univariate_kzg::KZGUnivariatePCS,
    PCSErrors, PolynomialCommitmentScheme, multilinear_kzg::util::get_uni_domain,
};
use ark_ec::{
    msm::{FixedBaseMSM, VariableBaseMSM},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseMultilinearExtension, MultilinearExtension, Polynomial,
    UVPolynomial,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, rand::RngCore, start_timer, vec::Vec, One, Zero};
use poly_iop::IOPTranscript;
use std::marker::PhantomData;

/// Input
/// - the prover parameters,
/// - a list of MLEs,
/// - and a same number of points,
/// - and a transcript,
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
    multi_commitment: &Commitment<E>,
    polynomials: &[DenseMultilinearExtension<E::Fr>],
    points: &[Vec<E::Fr>],
) -> Result<BatchProof<E>, PCSErrors> {
    let open_timer = start_timer!(|| "multi open");
    let mut transcript = IOPTranscript::new(b"ml kzg");
    transcript.append_serializable_element(b"w", multi_commitment)?;
    for point in points {
        transcript.append_serializable_element(b"w", point)?;
    }

    if points.len() != polynomials.len() {
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

    let num_vars = points[0].len();
    let uni_poly_degree = points.len();
    let domain = get_uni_domain::<E::Fr>(
        num_vars, uni_poly_degree
    )?;

    // 1. build `l(points)` which is a list of univariate polynomials that goes
    // through the points
    let uni_polys = build_l(num_var, points)?;

    // 2. build MLE `w` which is the merge of all MLEs.
    let merge_poly = merge_polynomials(polynomials)?;

    // 3. build `q(x)` which is a univariate polynomial `W circ l`
    let q_x = compute_w_circ_l(&merge_poly, &uni_polys)?;

    // 4. output `q(x)`' and put it into transcript
    //
    // TODO: use KZG commit for q(x)
    // TODO: unwrap
    let mut q_x_opens = vec![];
    let q_x_eval = q_x.clone().evaluate_over_domain(domain);
    for i in 0..uni_poly_degree{
        // assert_eq!(q_x_eval[i], polynomials[i].evaluate(&points[i]).unwrap());


    }


    q_x.coeffs
        .iter()
        .for_each(|x| transcript.append_field_element(b"q(x)", x).unwrap());

    let q_x_commit = KZGUnivariatePCS::<E>::commit(uni_prover_param, &q_x)?;

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

    end_timer!(open_timer);

    Ok(BatchProof {
        proof: opening,
        q_x_com: q_x.coeffs,
        value,
        q_x_commit,
        q_x_opens,
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
    values: &[E::Fr],
    batch_proof: &BatchProof<E>,
) -> Result<bool, PCSErrors> {
    let verify_timer = start_timer!(|| "batch verify");

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

    // TODO: verify commitment of `q(x)` instead of receiving full `q(x)`

    // 1. put `q(x)`'s evaluations over `(1, omega,...)` into transcript
    // TODO: unwrap
    batch_proof
        .q_x_com
        .iter()
        .for_each(|x| transcript.append_field_element(b"q(x)", x).unwrap());

    // 2. sample `r` from transcript
    let r = transcript.get_and_append_challenge(b"r")?;

    // 3. check `q(r) == value`
    let q_x = DensePolynomial::from_coefficients_slice(&batch_proof.q_x_com);
    let q_r = q_x.evaluate(&r);
    if q_r != batch_proof.value {
        return Ok(false);
    }
    if !KZGUnivariatePCS::verify(
        uni_verifier_param,
        &batch_proof.q_x_commit,
        &r,
        &q_r,
        &batch_proof.q_x_opens[0],
    )? {
        return Ok(false);
    }

    // 4. build `l(points)` which is a list of univariate polynomials that goes
    // through the points
    let uni_polys = build_l(num_var, points)?;

    // 5. get a point `p := l(r)`
    let point: Vec<E::Fr> = uni_polys.iter().rev().map(|x| x.evaluate(&r)).collect();

    // 6. verifies `p` is verifies against proof
    let res = KZGMultilinearPCS::verify(
        ml_verifier_param,
        multi_commitment,
        &point,
        &batch_proof.value,
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
        multilinear_kzg::util::{compute_uni_degree, get_srs_size},
        prelude::UnivariateUniversalParams,
        univariate_kzg::KZGUnivariatePCS,
        StructuredReferenceString,
    };
    use ark_bls12_381::{Bls12_381 as E, G1Affine};
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
        // let points_ref: Vec<&Vec<Fr>> = points.iter().map(|x| x.as_ref()).collect();

        let com = KZGMultilinearPCS::multi_commit(&ml_ck, polys)?;
        let batch_proof = multi_open_internal(&uni_ck, &ml_ck, &com, polys, &points)?;

        // good path
        assert!(batch_verify_internal(
            &uni_vk,
            &ml_vk,
            &com,
            &points,
            &[],
            &batch_proof,
        )?);

        // // bad commitment
        // assert!(!batch_verify_internal(
        //     &uni_vk,
        //     &ml_vk,
        //     &Commitment {
        //         commitment: <E as PairingEngine>::G1Affine::default()
        //     },
        //     &points,
        //     &[],
        //     &batch_proof,
        // )?);

        // // bad points
        // assert!(!batch_verify_internal(
        //     &uni_vk,
        //     &ml_vk,
        //     &com,
        //     &points[1..],
        //     &[],
        //     &batch_proof,
        // )?);

        // // bad proof
        // assert!(!batch_verify_internal(
        //     &uni_vk,
        //     &ml_vk,
        //     &com,
        //     &points,
        //     &[],
        //     &BatchProof {
        //         proof: Proof { proofs: Vec::new() },
        //         value: batch_proof.value,
        //         q_x_com: batch_proof.q_x_com.clone(),
        //         q_x_commits: vec![],
        //         q_x_opens: vec![],
        //     },
        // )?);

        // // bad value
        // assert!(!batch_verify_internal(
        //     &uni_vk,
        //     &ml_vk,
        //     &com,
        //     &points,
        //     &[],
        //     &BatchProof {
        //         proof: batch_proof.proof.clone(),
        //         value: Fr::one(),
        //         q_x_com: batch_proof.q_x_com,
        //         q_x_commits: vec![],
        //         q_x_opens: vec![],
        //     },
        // )?);

        // // bad q(x) commit
        // assert!(!batch_verify_internal(
        //     &uni_vk,
        //     &ml_vk,
        //     &com,
        //     &points,
        //     &[],
        //     &BatchProof {
        //         proof: batch_proof.proof,
        //         value: batch_proof.value,
        //         q_x_com: Vec::new(),
        //         q_x_commits: vec![],
        //         q_x_opens: vec![],
        //     },
        // )?);

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
