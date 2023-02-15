// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! Main module for multilinear KZG commitment scheme

pub(crate) mod batching;
pub(crate) mod srs;
pub(crate) mod util;

use crate::{
    pcs::{prelude::Commitment, PCSError, PolynomialCommitmentScheme, StructuredReferenceString},
    BatchProof,
};
use arithmetic::evaluate_opt;
use ark_ec::{
    pairing::Pairing,
    scalar_mul::{fixed_base::FixedBase, variable_base::VariableBaseMSM},
    AffineRepr, CurveGroup,
};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow, end_timer, format, marker::PhantomData, rand::Rng, start_timer,
    string::ToString, sync::Arc, vec, vec::Vec, One, Zero,
};
use std::ops::Mul;
// use batching::{batch_verify_internal, multi_open_internal};
use srs::{MultilinearProverParam, MultilinearUniversalParams, MultilinearVerifierParam};
use transcript::IOPTranscript;

use self::batching::{batch_verify_internal, multi_open_internal};

/// KZG Polynomial Commitment Scheme on multilinear polynomials.
pub struct MultilinearKzgPCS<E: Pairing> {
    #[doc(hidden)]
    phantom: PhantomData<E>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
/// proof of opening
pub struct MultilinearKzgProof<E: Pairing> {
    /// Evaluation of quotients
    pub proofs: Vec<E::G1Affine>,
}

impl<E: Pairing> PolynomialCommitmentScheme<E> for MultilinearKzgPCS<E> {
    // Parameters
    type ProverParam = MultilinearProverParam<E>;
    type VerifierParam = MultilinearVerifierParam<E>;
    type SRS = MultilinearUniversalParams<E>;
    // Polynomial and its associated types
    type Polynomial = Arc<DenseMultilinearExtension<E::ScalarField>>;
    type Point = Vec<E::ScalarField>;
    type Evaluation = E::ScalarField;
    // Commitments and proofs
    type Commitment = Commitment<E>;
    type Proof = MultilinearKzgProof<E>;
    type BatchProof = BatchProof<E, Self>;

    /// Build SRS for testing.
    ///
    /// - For univariate polynomials, `log_size` is the log of maximum degree.
    /// - For multilinear polynomials, `log_size` is the number of variables.
    ///
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: Rng>(rng: &mut R, log_size: usize) -> Result<Self::SRS, PCSError> {
        MultilinearUniversalParams::<E>::gen_srs_for_testing(rng, log_size)
    }

    /// Trim the universal parameters to specialize the public parameters.
    /// Input both `supported_log_degree` for univariate and
    /// `supported_num_vars` for multilinear.
    fn trim(
        srs: impl Borrow<Self::SRS>,
        supported_degree: Option<usize>,
        supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        assert!(supported_degree.is_none());

        let supported_num_vars = match supported_num_vars {
            Some(p) => p,
            None => {
                return Err(PCSError::InvalidParameters(
                    "multilinear should receive a num_var param".to_string(),
                ))
            },
        };
        let (ml_ck, ml_vk) = srs.borrow().trim(supported_num_vars)?;

        Ok((ml_ck, ml_vk))
    }

    /// Generate a commitment for a polynomial.
    ///
    /// This function takes `2^num_vars` number of scalar multiplications over
    /// G1.
    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, PCSError> {
        let prover_param = prover_param.borrow();
        let commit_timer = start_timer!(|| "commit");
        if prover_param.num_vars < poly.num_vars {
            return Err(PCSError::InvalidParameters(format!(
                "MlE length ({}) exceeds param limit ({})",
                poly.num_vars, prover_param.num_vars
            )));
        }
        let ignored = prover_param.num_vars - poly.num_vars;
        let scalars: Vec<_> = poly.to_evaluations();
        let msm_timer = start_timer!(|| format!(
            "msm of size {}",
            prover_param.powers_of_g[ignored].evals.len()
        ));
        let commitment =
            E::G1::msm_unchecked(&prover_param.powers_of_g[ignored].evals, scalars.as_slice())
                .into_affine();
        end_timer!(msm_timer);

        end_timer!(commit_timer);
        Ok(Commitment(commitment))
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the
    /// same. This function does not need to take the evaluation value as an
    /// input.
    ///
    /// This function takes 2^{num_var +1} number of scalar multiplications over
    /// G1:
    /// - it prodceeds with `num_var` number of rounds,
    /// - at round i, we compute an MSM for `2^{num_var - i + 1}` number of G2
    ///   elements.
    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomial: &Self::Polynomial,
        point: &Self::Point,
    ) -> Result<(Self::Proof, Self::Evaluation), PCSError> {
        open_internal(prover_param.borrow(), polynomial, point)
    }

    /// Input a list of multilinear extensions, and a same number of points, and
    /// a transcript, compute a multi-opening for all the polynomials.
    fn multi_open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomials: &[Self::Polynomial],
        points: &[Self::Point],
        evals: &[Self::Evaluation],
        transcript: &mut IOPTranscript<E::ScalarField>,
    ) -> Result<BatchProof<E, Self>, PCSError> {
        multi_open_internal(
            prover_param.borrow(),
            polynomials,
            points,
            evals,
            transcript,
        )
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
        point: &Self::Point,
        value: &E::ScalarField,
        proof: &Self::Proof,
    ) -> Result<bool, PCSError> {
        verify_internal(verifier_param, commitment, point, value, proof)
    }

    /// Verifies that `value_i` is the evaluation at `x_i` of the polynomial
    /// `poly_i` committed inside `comm`.
    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        commitments: &[Self::Commitment],
        points: &[Self::Point],
        batch_proof: &Self::BatchProof,
        transcript: &mut IOPTranscript<E::ScalarField>,
    ) -> Result<bool, PCSError> {
        batch_verify_internal(verifier_param, commitments, points, batch_proof, transcript)
    }
}

/// On input a polynomial `p` and a point `point`, outputs a proof for the
/// same. This function does not need to take the evaluation value as an
/// input.
///
/// This function takes 2^{num_var} number of scalar multiplications over
/// G1:
/// - it proceeds with `num_var` number of rounds,
/// - at round i, we compute an MSM for `2^{num_var - i}` number of G1 elements.
fn open_internal<E: Pairing>(
    prover_param: &MultilinearProverParam<E>,
    polynomial: &DenseMultilinearExtension<E::ScalarField>,
    point: &[E::ScalarField],
) -> Result<(MultilinearKzgProof<E>, E::ScalarField), PCSError> {
    let open_timer = start_timer!(|| format!("open mle with {} variable", polynomial.num_vars));

    if polynomial.num_vars() > prover_param.num_vars {
        return Err(PCSError::InvalidParameters(format!(
            "Polynomial num_vars {} exceed the limit {}",
            polynomial.num_vars, prover_param.num_vars
        )));
    }

    if polynomial.num_vars() != point.len() {
        return Err(PCSError::InvalidParameters(format!(
            "Polynomial num_vars {} does not match point len {}",
            polynomial.num_vars,
            point.len()
        )));
    }

    let nv = polynomial.num_vars();
    // the first `ignored` SRS vectors are unused for opening.
    let ignored = prover_param.num_vars - nv + 1;
    let mut f = polynomial.to_evaluations();

    let mut proofs = Vec::new();

    for (i, (&point_at_k, gi)) in point
        .iter()
        .zip(prover_param.powers_of_g[ignored..ignored + nv].iter())
        .enumerate()
    {
        let ith_round = start_timer!(|| format!("{}-th round", i));

        let k = nv - 1 - i;
        let cur_dim = 1 << k;
        let mut q = vec![E::ScalarField::zero(); cur_dim];
        let mut r = vec![E::ScalarField::zero(); cur_dim];

        let ith_round_eval = start_timer!(|| format!("{}-th round eval", i));
        for b in 0..(1 << k) {
            // q[b] = f[1, b] - f[0, b]
            q[b] = f[(b << 1) + 1] - f[b << 1];

            // r[b] = f[0, b] + q[b] * p
            r[b] = f[b << 1] + (q[b] * point_at_k);
        }
        f = r;
        end_timer!(ith_round_eval);

        // this is a MSM over G1 and is likely to be the bottleneck
        let msm_timer = start_timer!(|| format!("msm of size {} at round {}", gi.evals.len(), i));

        proofs.push(E::G1::msm_unchecked(&gi.evals, &q).into_affine());
        end_timer!(msm_timer);

        end_timer!(ith_round);
    }
    let eval = evaluate_opt(polynomial, point);
    end_timer!(open_timer);
    Ok((MultilinearKzgProof { proofs }, eval))
}

/// Verifies that `value` is the evaluation at `x` of the polynomial
/// committed inside `comm`.
///
/// This function takes
/// - num_var number of pairing product.
/// - num_var number of MSM
fn verify_internal<E: Pairing>(
    verifier_param: &MultilinearVerifierParam<E>,
    commitment: &Commitment<E>,
    point: &[E::ScalarField],
    value: &E::ScalarField,
    proof: &MultilinearKzgProof<E>,
) -> Result<bool, PCSError> {
    let verify_timer = start_timer!(|| "verify");
    let num_var = point.len();

    if num_var > verifier_param.num_vars {
        return Err(PCSError::InvalidParameters(format!(
            "point length ({}) exceeds param limit ({})",
            num_var, verifier_param.num_vars
        )));
    }

    let prepare_inputs_timer = start_timer!(|| "prepare pairing inputs");

    let scalar_size = E::ScalarField::MODULUS_BIT_SIZE as usize;
    let window_size = FixedBase::get_mul_window_size(num_var);

    let h_table =
        FixedBase::get_window_table(scalar_size, window_size, verifier_param.h.into_group());
    let h_mul: Vec<E::G2> = FixedBase::msm(scalar_size, window_size, &h_table, point);

    let ignored = verifier_param.num_vars - num_var;
    let h_vec: Vec<_> = (0..num_var)
        .map(|i| verifier_param.h_mask[ignored + i].into_group() - h_mul[i])
        .collect();
    let h_vec: Vec<E::G2Affine> = E::G2::normalize_batch(&h_vec);
    end_timer!(prepare_inputs_timer);

    let pairing_product_timer = start_timer!(|| "pairing product");

    let mut pairings: Vec<_> = proof
        .proofs
        .iter()
        .map(|&x| E::G1Prepared::from(x))
        .zip(h_vec.into_iter().take(num_var).map(E::G2Prepared::from))
        .collect();

    pairings.push((
        E::G1Prepared::from(
            (verifier_param.g.mul(*value) - commitment.0.into_group()).into_affine(),
        ),
        E::G2Prepared::from(verifier_param.h),
    ));

    let ps = pairings.iter().map(|(p, _)| p.clone());
    let hs = pairings.iter().map(|(_, h)| h.clone());

    let res = E::multi_pairing(ps, hs) == ark_ec::pairing::PairingOutput(E::TargetField::one());

    end_timer!(pairing_product_timer);
    end_timer!(verify_timer);
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::{test_rng, vec::Vec, UniformRand};

    type E = Bls12_381;
    type Fr = <E as Pairing>::ScalarField;

    fn test_single_helper<R: Rng>(
        params: &MultilinearUniversalParams<E>,
        poly: &Arc<DenseMultilinearExtension<Fr>>,
        rng: &mut R,
    ) -> Result<(), PCSError> {
        let nv = poly.num_vars();
        assert_ne!(nv, 0);
        let (ck, vk) = MultilinearKzgPCS::trim(params, None, Some(nv))?;
        let point: Vec<_> = (0..nv).map(|_| Fr::rand(rng)).collect();
        let com = MultilinearKzgPCS::commit(&ck, poly)?;
        let (proof, value) = MultilinearKzgPCS::open(&ck, poly, &point)?;

        assert!(MultilinearKzgPCS::verify(
            &vk, &com, &point, &value, &proof
        )?);

        let value = Fr::rand(rng);
        assert!(!MultilinearKzgPCS::verify(
            &vk, &com, &point, &value, &proof
        )?);

        Ok(())
    }

    #[test]
    fn test_single_commit() -> Result<(), PCSError> {
        let mut rng = test_rng();

        let params = MultilinearKzgPCS::<E>::gen_srs_for_testing(&mut rng, 10)?;

        // normal polynomials
        let poly1 = Arc::new(DenseMultilinearExtension::rand(8, &mut rng));
        test_single_helper(&params, &poly1, &mut rng)?;

        // single-variate polynomials
        let poly2 = Arc::new(DenseMultilinearExtension::rand(1, &mut rng));
        test_single_helper(&params, &poly2, &mut rng)?;

        Ok(())
    }

    #[test]
    fn setup_commit_verify_constant_polynomial() {
        let mut rng = test_rng();

        // normal polynomials
        assert!(MultilinearKzgPCS::<E>::gen_srs_for_testing(&mut rng, 0).is_err());
    }
}
