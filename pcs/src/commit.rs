use crate::{
    util::{build_l, compute_w_circ_l, merge_polynomials},
    KZGMultilinearPC, MultilinearCommitmentScheme, PCSErrors, ProverParam, UniversalParams,
    VerifierParam,
};
use ark_ec::{
    msm::{FixedBaseMSM, VariableBaseMSM},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, MultilinearExtension, Polynomial, UVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, log2, rand::RngCore, start_timer, vec::Vec, One, Zero};
use poly_iop::IOPTranscript;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// commitment
pub struct Commitment<E: PairingEngine> {
    /// number of variables
    pub nv: usize,
    /// product of g as described by the vRAM paper
    pub g_product: E::G1Affine,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// proof of opening
pub struct Proof<E: PairingEngine> {
    /// Evaluation of quotients
    pub proofs: Vec<E::G1Affine>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// proof of batch opening
pub struct BatchProof<E: PairingEngine> {
    /// The actual proof
    pub proof: Proof<E>,
    /// The value which is `w` evaluated at `p:= l(r)`, where
    /// - `w` is the merged MLE
    /// - `l` is the list of univariate polys that goes through all points
    /// - `r` is sampled from the transcript.
    pub value: E::Fr,
    /// Commitment to q(x)
    // This is currently set to the entire coefficient list of q(x)
    // TODO: replace me with a KZG commit
    pub q_x_com: Vec<E::Fr>,
}

impl<E: PairingEngine> MultilinearCommitmentScheme<E> for KZGMultilinearPC<E> {
    type ProverParam = ProverParam<E>;
    type VerifierParam = VerifierParam<E>;
    type SRS = UniversalParams<E>;
    type Commitment = Commitment<E>;
    type Proof = Proof<E>;
    type Transcript = IOPTranscript<E::Fr>;
    type BatchProof = BatchProof<E>;

    /// Generate SRS from RNG.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn setup<R: RngCore>(rng: &mut R, num_vars: usize) -> Result<Self::SRS, PCSErrors> {
        let setup_timer = start_timer!(|| format!("SRS setup for dim {}", num_vars));
        let res = Self::SRS::gen_srs_for_testing(rng, num_vars);
        end_timer!(setup_timer);
        res
    }

    /// Generate a commitment for a polynomial.
    ///
    /// This function takes `2^num_vars` number of scalar multiplications over
    /// G1.
    fn commit(
        prover_param: &Self::ProverParam,
        poly: &impl MultilinearExtension<E::Fr>,
    ) -> Result<Self::Commitment, PCSErrors> {
        let commit_timer = start_timer!(|| "commit");

        let nv = poly.num_vars();
        let scalars: Vec<_> = poly
            .to_evaluations()
            .into_iter()
            .map(|x| x.into_repr())
            .collect();
        let g_product = VariableBaseMSM::multi_scalar_mul(
            &prover_param.powers_of_g[0].evals,
            scalars.as_slice(),
        )
        .into_affine();

        end_timer!(commit_timer);
        Ok(Commitment { nv, g_product })
    }

    /// Generate a commitment for a list of polynomials.
    ///
    /// This function takes `2^(num_vars + log(polys.len())` number of scalar
    /// multiplications over G1.
    fn multi_commit(
        prover_param: &Self::ProverParam,
        polys: &[impl MultilinearExtension<E::Fr>],
    ) -> Result<Self::Commitment, PCSErrors> {
        let commit_timer = start_timer!(|| "multi commit");
        let poly = merge_polynomials(polys)?;

        let scalars: Vec<_> = poly
            .to_evaluations()
            .iter()
            .map(|x| x.into_repr())
            .collect();

        let g_product = VariableBaseMSM::multi_scalar_mul(
            &prover_param.powers_of_g[0].evals,
            scalars.as_slice(),
        )
        .into_affine();

        end_timer!(commit_timer);
        Ok(Commitment {
            nv: poly.num_vars,
            g_product,
        })
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same. This function does not need to take the evaluation value as an
    /// input.
    ///
    /// This function takes 2^{num_var +1} number of scalar multiplications over
    /// G1:
    /// - it proceeds with `num_var` number of rounds,
    /// - at round i, we compute an MSM for `2^{num_var - i + 1}` number of G2
    ///   elements.
    fn open(
        prover_param: &Self::ProverParam,
        polynomial: &impl MultilinearExtension<E::Fr>,
        point: &[E::Fr],
    ) -> Result<Self::Proof, PCSErrors> {
        let open_timer = start_timer!(|| "open");

        assert_eq!(
            polynomial.num_vars(),
            prover_param.num_vars,
            "Invalid size of polynomial"
        );
        let nv = polynomial.num_vars();
        let mut r: Vec<Vec<E::Fr>> = (0..nv + 1).map(|_| Vec::new()).collect();
        let mut q: Vec<Vec<E::Fr>> = (0..nv + 1).map(|_| Vec::new()).collect();

        r[nv] = polynomial.to_evaluations();

        let mut proofs = Vec::new();

        for (i, (&point_at_k, gi)) in point
            .iter()
            .zip(prover_param.powers_of_g.iter())
            .take(nv)
            .enumerate()
        {
            let ith_round = start_timer!(|| format!("{}-th round", i));

            let k = nv - i;
            let cur_dim = 1 << (k - 1);
            let mut cur_q = vec![E::Fr::zero(); cur_dim];
            let mut cur_r = vec![E::Fr::zero(); cur_dim];

            for b in 0..(1 << (k - 1)) {
                // q_b = pre_r [2^b + 1] - pre_r [2^b]
                cur_q[b] = r[k][(b << 1) + 1] - r[k][b << 1];

                // r_b = pre_r [2^b]*(1-p) + pre_r [2^b + 1] * p
                cur_r[b] =
                    r[k][b << 1] * (E::Fr::one() - point_at_k) + (r[k][(b << 1) + 1] * point_at_k);
            }

            let scalars: Vec<_> = (0..(1 << k)).map(|x| cur_q[x >> 1].into_repr()).collect();

            q[k] = cur_q;
            r[k - 1] = cur_r;

            // this is a MSM over G1 and is likely to be the bottleneck
            proofs.push(VariableBaseMSM::multi_scalar_mul(&gi.evals, &scalars).into_affine());
            end_timer!(ith_round);
        }

        end_timer!(open_timer);
        Ok(Proof { proofs })
    }

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
    /// - q(x), which is a univariate polynomial `w circ l` where `w` is the
    ///   merged MLE, and `l` is a list of polynomials that go through all the
    ///   points. TODO: change this field to a commitment to `q(x)`
    /// - and a value which is `w` evaluated at `p:= l(r)` from some `r` from
    ///   the transcript.
    ///
    /// Steps:
    /// 1. build `l(points)` which is a list of univariate polynomials that goes
    /// through the points
    /// 2. build MLE `w` which is the merge of all MLEs.
    /// 3. build `q(x)` which is a univariate polynomial `W circ l`
    /// 4. output `q(x)`' and put it into transcript
    /// 5. sample `r` from transcript
    /// 6. get a point `p := l(r)`
    /// 7. output an opening of `w` over point `p`
    /// 8. output `w(p)`
    fn multi_open(
        prover_param: &Self::ProverParam,
        polynomials: &[impl MultilinearExtension<E::Fr>],
        points: &[&[E::Fr]],
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<Self::BatchProof, PCSErrors> {
        let open_timer = start_timer!(|| "multi open");

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
        for &point in points.iter() {
            if point.len() != num_var {
                return Err(PCSErrors::InvalidParameters(
                    "points do not have same num_vars".to_string(),
                ));
            }
        }

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
        q_x.coeffs
            .iter()
            .for_each(|x| transcript.append_field_element(b"q(x)", x).unwrap());

        // 5. sample `r` from transcript
        let r = transcript.get_and_append_challenge(b"r")?;

        // 6. get a point `p := l(r)`
        let point: Vec<E::Fr> = uni_polys.iter().map(|poly| poly.evaluate(&r)).collect();

        // 7. output an opening of `w` over point `p`
        let opening = Self::open(prover_param, &merge_poly, &point)?;

        // 8. output value that is `w` evaluated at `p`
        let value = merge_poly.evaluate(&point).unwrap();
        end_timer!(open_timer);

        Ok(Self::BatchProof {
            proof: opening,
            q_x_com: q_x.coeffs,
            value,
        })
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    ///
    /// This function takes
    /// - num_var number of pairing product.
    /// - num_var number of MSM
    fn verify(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        point: &[E::Fr],
        value: &E::Fr,
        proof: &Self::Proof,
    ) -> Result<bool, PCSErrors> {
        let verify_timer = start_timer!(|| "verify");
        let prepare_inputs_timer = start_timer!(|| "prepare pairing inputs");

        let scalar_size = E::Fr::size_in_bits();
        let window_size = FixedBaseMSM::get_mul_window_size(verifier_param.num_vars);

        let h_table = FixedBaseMSM::get_window_table(
            scalar_size,
            window_size,
            verifier_param.h.into_projective(),
        );
        let h_mul: Vec<E::G2Projective> =
            FixedBaseMSM::multi_scalar_mul(scalar_size, window_size, &h_table, point);

        let h_vec: Vec<_> = (0..verifier_param.num_vars)
            .map(|i| verifier_param.h_mask[i].into_projective() - h_mul[i])
            .collect();
        let h_vec: Vec<E::G2Affine> = E::G2Projective::batch_normalization_into_affine(&h_vec);
        end_timer!(prepare_inputs_timer);

        let pairing_product_timer = start_timer!(|| "pairing product");

        let mut pairings: Vec<_> = proof
            .proofs
            .iter()
            .map(|&x| E::G1Prepared::from(x))
            .zip(
                h_vec
                    .into_iter()
                    .take(verifier_param.num_vars)
                    .map(E::G2Prepared::from),
            )
            .collect();

        pairings.push((
            E::G1Prepared::from(
                (verifier_param.g.mul(*value) - commitment.g_product.into_projective())
                    .into_affine(),
            ),
            E::G2Prepared::from(verifier_param.h),
        ));

        let res = E::product_of_pairings(pairings.iter()) == E::Fqk::one();

        end_timer!(pairing_product_timer);
        end_timer!(verify_timer);
        Ok(res)
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
    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        multi_commitment: &Self::Commitment,
        points: &[&[E::Fr]],
        batch_proof: &Self::BatchProof,
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<bool, PCSErrors> {
        let verify_timer = start_timer!(|| "batch verify");
        let num_var = points[0].len();

        for &point in points.iter().skip(1) {
            if point.len() != num_var {
                return Err(PCSErrors::InvalidParameters(format!(
                    "points do not have same num_vars ({} vs {})",
                    point.len(),
                    num_var,
                )));
            }
        }
        if num_var + log2(points.len()) as usize != multi_commitment.nv {
            return Err(PCSErrors::InvalidParameters(format!(
                "points and multi_commitment do not have same num_vars ({} vs {})",
                num_var + log2(points.len()) as usize,
                num_var,
            )));
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
            println!("univariate failed");
            return Ok(false);
        }

        // 4. build `l(points)` which is a list of univariate polynomials that goes
        // through the points
        let uni_polys = build_l(num_var, points)?;

        // 5. get a point `p := l(r)`
        let point: Vec<E::Fr> = uni_polys.iter().map(|x| x.evaluate(&r)).collect();

        // 6. verifies `p` is verifies against proof
        let res = Self::verify(
            verifier_param,
            multi_commitment,
            &point,
            &batch_proof.value,
            &batch_proof.proof,
        );
        end_timer!(verify_timer);

        res
    }
}

#[cfg(test)]
mod tests {
    use crate::util::get_batched_nv;

    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::{rand::RngCore, test_rng, vec::Vec, UniformRand};
    type E = Bls12_381;
    type Fr = <E as PairingEngine>::Fr;

    fn test_single_helper<R: RngCore>(
        uni_params: &UniversalParams<E>,
        poly: &impl MultilinearExtension<Fr>,
        rng: &mut R,
    ) -> Result<(), PCSErrors> {
        let nv = poly.num_vars();
        assert_ne!(nv, 0);
        let (ck, vk) = uni_params.trim(nv)?;
        let point: Vec<_> = (0..nv).map(|_| Fr::rand(rng)).collect();
        let com = KZGMultilinearPC::commit(&ck, poly)?;
        let proof = KZGMultilinearPC::open(&ck, poly, &point)?;

        let value = poly.evaluate(&point).unwrap();
        assert!(KZGMultilinearPC::verify(&vk, &com, &point, &value, &proof)?);

        let value = Fr::rand(rng);
        assert!(!KZGMultilinearPC::verify(
            &vk, &com, &point, &value, &proof
        )?);

        Ok(())
    }

    #[test]
    fn test_single_commit() -> Result<(), PCSErrors> {
        let mut rng = test_rng();

        let uni_params = KZGMultilinearPC::<E>::setup(&mut rng, 10)?;

        // normal polynomials
        let poly1 = DenseMultilinearExtension::rand(8, &mut rng);
        test_single_helper(&uni_params, &poly1, &mut rng)?;

        // single-variate polynomials
        let poly2 = DenseMultilinearExtension::rand(1, &mut rng);
        test_single_helper(&uni_params, &poly2, &mut rng)?;

        Ok(())
    }

    fn test_multi_commit_helper<R: RngCore>(
        uni_params: &UniversalParams<E>,
        polys: &[impl MultilinearExtension<Fr>],
        rng: &mut R,
    ) -> Result<(), PCSErrors> {
        let mut transcript = IOPTranscript::new(b"test");

        let nv = get_batched_nv(polys[0].num_vars(), polys.len());
        let (ck, vk) = uni_params.trim(nv)?;
        let mut points = Vec::new();

        for poly in polys.iter() {
            let point = (0..poly.num_vars())
                .map(|_| Fr::rand(rng))
                .collect::<Vec<Fr>>();
            points.push(point);
        }
        let points_ref: Vec<&[Fr]> = points.iter().map(|x| x.as_ref()).collect();

        let com = KZGMultilinearPC::multi_commit(&ck, polys)?;
        let batch_proof = KZGMultilinearPC::multi_open(&ck, polys, &points_ref, &mut transcript)?;

        let mut transcript = IOPTranscript::new(b"test");
        assert!(KZGMultilinearPC::batch_verify(
            &vk,
            &com,
            &points_ref,
            &batch_proof,
            &mut transcript
        )?);

        Ok(())
    }

    #[test]
    fn test_multi_commit() -> Result<(), PCSErrors> {
        let mut rng = test_rng();

        let uni_params = KZGMultilinearPC::<E>::setup(&mut rng, 15)?;

        // normal polynomials
        let polys1: Vec<_> = (0..5)
            .map(|_| DenseMultilinearExtension::rand(4, &mut rng))
            .collect();
        test_multi_commit_helper(&uni_params, &polys1, &mut rng)?;

        // single-variate polynomials
        let polys1: Vec<_> = (0..5)
            .map(|_| DenseMultilinearExtension::rand(1, &mut rng))
            .collect();
        test_multi_commit_helper(&uni_params, &polys1, &mut rng)?;

        Ok(())
    }

    #[test]
    fn setup_commit_verify_constant_polynomial() {
        let mut rng = test_rng();

        // normal polynomials
        assert!(KZGMultilinearPC::<E>::setup(&mut rng, 0).is_err());
    }
}
