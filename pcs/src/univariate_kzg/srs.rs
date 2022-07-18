use crate::{prelude::PCSErrors, StructuredReferenceString};
use ark_ec::{msm::FixedBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, rand::RngCore, start_timer, One, UniformRand};
use std::collections::BTreeMap;

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
/// Adapted from
/// https://github.com/arkworks-rs/poly-commit/blob/master/src/kzg10/data_structures.rs#L20
#[derive(Clone, Debug)]
pub struct UniversalParams<E: PairingEngine> {
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to
    /// `degree`.
    pub powers_of_g: Vec<E::G1Affine>,
    /// Group elements of the form `{ \beta^i \gamma G }`, where `i` ranges from
    /// 0 to `degree`.
    pub powers_of_gamma_g: BTreeMap<usize, E::G1Affine>,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
    // /// Group elements of the form `{ \beta^i G2 }`, where `i` ranges from `0`
    // /// to `-degree`.
    // pub neg_powers_of_h: BTreeMap<usize, E::G2Affine>,
    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_h: E::G2Prepared,
}

/// Prover Parameters
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct ProverParam<C: AffineCurve> {
    pub g_powers: Vec<C>,
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
/// https://github.com/arkworks-rs/poly-commit/blob/master/src/kzg10/data_structures.rs#L236
#[derive(Clone, Debug)]
pub struct VerifierParam<E: PairingEngine> {
    /// The generator of G1.
    pub g: E::G1Affine,
    /// The generator of G1 that is used for making a commitment hiding.
    pub gamma_g: E::G1Affine,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_h: E::G2Prepared,
}

impl<E: PairingEngine> StructuredReferenceString<E> for UniversalParams<E> {
    type ProverParam = ProverParam<E::G1Affine>;
    type VerifierParam = VerifierParam<E>;

    /// Extract the prover parameters from the public parameters.
    fn extract_prover_param(&self) -> Self::ProverParam {
        todo!()
    }

    /// Extract the verifier parameters from the public parameters.
    fn extract_verifier_param(&self) -> Self::VerifierParam {
        todo!()
    }

    /// Trim the universal parameters to specialize the public parameters
    /// for multilinear polynomials to the given `supported_num_vars`, and
    /// returns committer key and verifier key. `supported_num_vars` should
    /// be in range `1..=params.num_vars`
    fn trim(
        &self,
        supported_num_vars: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSErrors> {
        todo!()
    }

    /// Build SRS for testing.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore>(rng: &mut R, log_degree: usize) -> Result<Self, PCSErrors> {
        let max_degree = 1usize << log_degree;
        let setup_time = start_timer!(|| format!("KZG10::Setup with degree {}", max_degree));
        let beta = E::Fr::rand(rng);
        let g = E::G1Projective::rand(rng);
        let gamma_g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        let mut powers_of_beta = vec![E::Fr::one()];

        let mut cur = beta;
        for _ in 0..max_degree {
            powers_of_beta.push(cur);
            cur *= &beta;
        }

        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);

        let scalar_bits = E::Fr::size_in_bits();
        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);
        let powers_of_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers_of_beta,
        );
        end_timer!(g_time);
        let gamma_g_time = start_timer!(|| "Generating powers of gamma * G");
        let gamma_g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, gamma_g);
        let mut powers_of_gamma_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &gamma_g_table,
            &powers_of_beta,
        );
        // Add an additional power of gamma_g, because we want to be able to support
        // up to D queries.
        powers_of_gamma_g.push(powers_of_gamma_g.last().unwrap().mul(beta.into_repr()));
        end_timer!(gamma_g_time);

        let powers_of_g = E::G1Projective::batch_normalization_into_affine(&powers_of_g);
        let powers_of_gamma_g =
            E::G1Projective::batch_normalization_into_affine(&powers_of_gamma_g)
                .into_iter()
                .enumerate()
                .collect();

        let h = h.into_affine();
        let beta_h = h.mul(beta).into_affine();
        let prepared_h = h.into();
        let prepared_beta_h = beta_h.into();

        let pp = UniversalParams {
            powers_of_g,
            powers_of_gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        };
        end_timer!(setup_time);
        Ok(pp)
    }
}
