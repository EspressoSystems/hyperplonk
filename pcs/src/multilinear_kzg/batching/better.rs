use crate::{
    multilinear_kzg::{
        open_internal,
        srs::{MultilinearProverParam, MultilinearVerifierParam},
        util::{bit_decompose, build_l, compute_w_circ_l, gen_eval_point, get_uni_domain},
        verify_internal, MultilinearKzgBatchProof,
    },
    prelude::{Commitment, UnivariateProverParam, UnivariateVerifierParam},
    univariate_kzg::UnivariateKzgPCS,
    PCSError, PolynomialCommitmentScheme,
};
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::{
    DenseMultilinearExtension, EvaluationDomain, MultilinearExtension, Polynomial,
    Radix2EvaluationDomain,
};
use ark_std::{end_timer, format, rc::Rc, start_timer, string::ToString, vec, vec::Vec};
use transcript::IOPTranscript;

#[derive(Debug)]
pub struct PolyAndPoints<E: PairingEngine> {
    polynomial: Rc<DenseMultilinearExtension<E::Fr>>,
    commitment: Commitment<E>,
    points: Vec<Vec<E::Fr>>,
}

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
pub(crate) fn multi_open_better<E: PairingEngine>(
    uni_prover_param: &UnivariateProverParam<E::G1Affine>,
    ml_prover_param: &MultilinearProverParam<E>,
    poly_and_points_vec: &[PolyAndPoints<E>],
) -> Result<(MultilinearKzgBatchProof<E>, Vec<E::Fr>), PCSError> {
    let open_timer = start_timer!(|| "multi open");

    let num_vars = poly_and_points_vec
        .iter()
        .map(|x| x.polynomial.num_vars())
        .max()
        .unwrap();

    let mut points = vec![];
    let mut index = 0;

    for poly_and_points in poly_and_points_vec.iter() {
        for point in poly_and_points.points.iter() {
            points.push(gen_eval_point(index, num_vars - point.len(), point));
            index += 1;
        }
    }
    let domain = get_uni_domain(points.len())?;
    let l = build_l(&points, &domain)?;

    end_timer!(open_timer);

    todo!()
}

#[cfg(test)]
mod test {

    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ff::UniformRand;
    use ark_std::test_rng;

    use crate::{
        multilinear_kzg::{srs::MultilinearUniversalParams, MultilinearKzgPCS},
        univariate_kzg::srs::UnivariateUniversalParams,
        StructuredReferenceString,
    };

    use super::*;

    type E = Bls12_381;
    #[test]
    fn test_batch_open_better() -> Result<(), PCSError> {
        let mut rng = test_rng();

        let uni_params =
            UnivariateUniversalParams::<E>::gen_srs_for_testing(&mut rng, 1usize << 15)?;
        let ml_params = MultilinearUniversalParams::<E>::gen_srs_for_testing(&mut rng, 15)?;

        let (uni_ck, uni_vk) = uni_params.trim(200)?;
        let (ml_ck, ml_vk) = ml_params.trim(5)?;

        // First poly and points
        let a_nv = 3;
        let a = Rc::new(DenseMultilinearExtension::rand(a_nv, &mut rng));
        let mut a_points = Vec::new();
        for _ in 0..3 {
            let point = (0..a.num_vars())
                .map(|_| Fr::rand(&mut rng))
                .collect::<Vec<Fr>>();
            a_points.push(point);
        }
        let a_com = MultilinearKzgPCS::<E>::commit(&(ml_ck.clone(), uni_ck.clone()), &a)?;
        let a_poly_and_point = PolyAndPoints {
            polynomial: a,
            commitment: a_com,
            points: a_points,
        };

        let b_nv = 2;
        let b = Rc::new(DenseMultilinearExtension::rand(b_nv, &mut rng));
        let mut b_points = Vec::new();
        for _ in 0..4 {
            let point = (0..b.num_vars())
                .map(|_| Fr::rand(&mut rng))
                .collect::<Vec<Fr>>();
            b_points.push(point);
        }

        let b_com = MultilinearKzgPCS::<E>::commit(&(ml_ck.clone(), uni_ck.clone()), &b)?;

        let b_poly_and_point = PolyAndPoints {
            polynomial: b,
            commitment: b_com,
            points: b_points,
        };

        let (proof, open) =
            multi_open_better(&uni_ck, &ml_ck, &[a_poly_and_point, b_poly_and_point])?;

        Ok(())
    }
}
