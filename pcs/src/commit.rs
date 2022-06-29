use crate::{
    util::compute_w_circ_l, KZGMultilinearPC, MultilinearCommitmentScheme, PCSErrors, ProverParam,
    UniversalParams, VerifierParam,
};
use ark_ec::{
    msm::{FixedBaseMSM, VariableBaseMSM},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::PrimeField;
use ark_poly::{
    DenseMultilinearExtension, EvaluationDomain, Evaluations, MultilinearExtension, Polynomial,
    Radix2EvaluationDomain,
};
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

impl<E: PairingEngine> MultilinearCommitmentScheme<E> for KZGMultilinearPC<E> {
    type ProverParam = ProverParam<E>;
    type VerifierParam = VerifierParam<E>;
    type SRS = UniversalParams<E>;
    type Commitment = Commitment<E>;
    type Proof = Proof<E>;

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
        let commit_timer = start_timer!(|| "commit");
        let poly = merge_polynomials(polys);

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

    /// Steps:
    /// 1. build `l(points)` which is a list of univariate polynomials that goes
    /// through the points
    /// 2. build MLE `w` which is the merge of all MLEs.
    /// 3. build `q(x)` which is a univariate polynomial `W circ l`
    /// 4. output `q(x)`'s evaluations over `(1, omega,...)`, and put it into
    /// transcript
    /// 5. sample `r` from transcript
    /// 6. get a point `p := l(r)`
    /// 7. output an opening of `w` over point `p`
    fn multi_open(
        prover_param: &Self::ProverParam,
        polynomials: &[impl MultilinearExtension<E::Fr>],
        points: &[&[E::Fr]],
        transcript: &mut IOPTranscript<E::Fr>,
    ) -> Result<(Self::Proof, Vec<E::Fr>), PCSErrors> {
        let open_timer = start_timer!(|| "open");

        if points.len() != polynomials.len() {
            return Err(PCSErrors::InvalidParameters(
                "polynomial length does not match point length".to_string(),
            ));
        }

        let uni_degree = points.len() + 1;
        let small_domain = match Radix2EvaluationDomain::<E::Fr>::new(uni_degree) {
            Some(p) => p,
            None => {
                return Err(PCSErrors::InvalidParameters(
                    "failed to build radix 2 domain".to_string(),
                ))
            },
        };
        let num_var = polynomials[0].num_vars();
        let large_domain = match Radix2EvaluationDomain::<E::Fr>::new(1 << get_nv(polynomials)) {
            Some(p) => p,
            None => {
                return Err(PCSErrors::InvalidParameters(
                    "failed to build radix 2 domain".to_string(),
                ))
            },
        };

        // 1. build `l(points)` which is a list of univariate polynomials that goes
        // through the points
        let mut uni_polys = Vec::new();
        for i in 0..1 << num_var {
            let eval: Vec<E::Fr> = points.iter().map(|x| x[i]).collect();
            uni_polys.push(Evaluations::from_vec_and_domain(eval, small_domain).interpolate())
        }
        // 2. build MLE `w` which is the merge of all MLEs.
        let merge_poly = merge_polynomials(polynomials);

        // 3. build `q(x)` which is a univariate polynomial `W circ l`
        let q_x = compute_w_circ_l(&merge_poly, &uni_polys)?;

        // 4. output `q(x)`'s evaluations over `(1, omega,...)`, and put it into
        // transcript
        //
        // TODO: use KZG commit for q(x)
        let evals = q_x.evaluate_over_domain(large_domain);
        evals
            .evals
            .iter()
            .for_each(|x| transcript.append_field_element(b"q(x)", x).unwrap());

        // 5. sample `r` from transcript
        let r = transcript.get_and_append_challenge(b"r")?;

        // 6. get a point `p := l(r)`
        let point: Vec<E::Fr> = uni_polys.iter().map(|poly| poly.evaluate(&r)).collect();

        // 7. output an opening of `w` over point `p`
        let opening = Self::open(prover_param, &merge_poly, &point)?;

        end_timer!(open_timer);

        Ok((opening, evals.evals))
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
        value: E::Fr,
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
                (verifier_param.g.mul(value) - commitment.g_product.into_projective())
                    .into_affine(),
            ),
            E::G2Prepared::from(verifier_param.h),
        ));

        let res = E::product_of_pairings(pairings.iter()) == E::Fqk::one();

        end_timer!(pairing_product_timer);
        end_timer!(verify_timer);
        Ok(res)
    }
}

/// Return the number of variables that one need for an MLE to
/// batch the list of MLEs
#[inline]
fn get_nv<F: PrimeField>(polynomials: &[impl MultilinearExtension<F>]) -> usize {
    polynomials.iter().map(|x| x.num_vars()).max().unwrap() + log2(polynomials.len()) as usize
}

fn merge_polynomials<F: PrimeField>(
    polynomials: &[impl MultilinearExtension<F>],
) -> DenseMultilinearExtension<F> {
    let individual_max = polynomials.iter().map(|x| x.num_vars()).max()
    // safe unwrap since polynomials is not empty
    .unwrap();
    let new_nv = get_nv(polynomials);
    let mut scalars = vec![];
    for poly in polynomials.iter() {
        let mut scalar = poly.to_evaluations();

        if poly.num_vars() < individual_max {
            scalar
                .extend_from_slice(vec![F::zero(); (1 << individual_max) - scalar.len()].as_ref());
        }
        scalars.extend_from_slice(scalar.as_slice());
    }
    scalars.extend_from_slice(vec![F::zero(); (1 << new_nv) - scalars.len()].as_ref());
    DenseMultilinearExtension::from_evaluations_vec(new_nv, scalars)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::{
        rand::{Rng, RngCore},
        test_rng,
        vec::Vec,
        UniformRand,
    };
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
        assert!(KZGMultilinearPC::verify(&vk, &com, &point, value, &proof)?);

        let value = Fr::rand(rng);
        assert!(!KZGMultilinearPC::verify(&vk, &com, &point, value, &proof)?);

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

        let nv = get_nv(polys);
        let (ck, _vk) = uni_params.trim(nv)?;
        let mut points = Vec::new();

        for poly in polys.iter() {
            let point = (0..poly.num_vars())
                .map(|_| Fr::rand(rng))
                .collect::<Vec<Fr>>();
            points.push(point);
        }
        let points_ref: Vec<&[Fr]> = points.iter().map(|x| x.as_ref()).collect();

        let _com = KZGMultilinearPC::multi_commit(&ck, polys)?;
        let _proof = KZGMultilinearPC::multi_open(&ck, polys, &points_ref, &mut transcript)?;

        Ok(())
    }

    #[test]
    fn test_multi_commit() -> Result<(), PCSErrors> {
        let mut rng = test_rng();

        let uni_params = KZGMultilinearPC::<E>::setup(&mut rng, 20)?;

        // normal polynomials
        let polys1: Vec<_> = (0..5)
            .map(|_| DenseMultilinearExtension::rand(rng.gen_range(5..9), &mut rng))
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
