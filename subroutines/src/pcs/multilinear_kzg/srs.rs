// Copyright (c) 2022 Espresso Systems (espressosys.com)
// This file is part of the Jellyfish library.

// You should have received a copy of the MIT License
// along with the Jellyfish library. If not, see <https://mit-license.org/>.

//! Implementing Structured Reference Strings for multilinear polynomial KZG
use std::ops::Add;

use ark_ec::{AffineCurve, msm::FixedBaseMSM, PairingEngine, ProjectiveCurve};
use ark_ff::{PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{
    end_timer, format,
    rand::{CryptoRng, RngCore},
    start_timer,
    string::ToString,
    UniformRand,
    vec::Vec,
};
use rayon::iter::IntoParallelRefIterator;

use arithmetic::{build_eq_x_r_vec};

use crate::pcs::{
    prelude::PCSError, StructuredReferenceString,
};

/// Evaluations over {0,1}^n for G1 or G2
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct Evaluations<C: AffineCurve> {
    /// The evaluations.
    pub evals: Vec<C>,
}

/// Universal Parameter
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct MultilinearUniversalParams<E: PairingEngine> {
    /// prover parameters
    pub prover_param: MultilinearProverParam<E>,
    /// h^randomness: h^t1, h^t2, ..., **h^{t_nv}**
    pub h_mask: Vec<E::G2Affine>,
}

/// Prover Parameters
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct MultilinearProverParam<E: PairingEngine> {
    /// number of variables
    pub num_vars: usize,

    /// `pp_{0}`, `pp_{1}`, ...,pp_{nu_vars} defined
    /// by XZZPD19 where pp_0=g and pp_{i}=g^{eq((t_1,..t_i),(X_1,..X_i))}
    pub powers_of_g: Vec<Evaluations<E::G1Affine>>,
    /// generator for G1
    pub g: E::G1Affine,
    /// generator for G2
    pub h: E::G2Affine,
}

/// Verifier Parameters
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct MultilinearVerifierParam<E: PairingEngine> {
    /// number of variables
    pub num_vars: usize,
    /// generator of G1
    pub g: E::G1Affine,
    /// generator of G2
    pub h: E::G2Affine,
    /// h^randomness: h^t1, h^t2, ..., **h^{t_nv}**
    pub h_mask: Vec<E::G2Affine>,
}

impl<E: PairingEngine> StructuredReferenceString<E> for MultilinearUniversalParams<E> {
    type ProverParam = MultilinearProverParam<E>;
    type VerifierParam = MultilinearVerifierParam<E>;

    /// Extract the prover parameters from the public parameters.
    fn extract_prover_param(&self, supported_num_vars: usize) -> Self::ProverParam {
        Self::ProverParam {
            powers_of_g: self.prover_param.powers_of_g[..supported_num_vars + 1].to_vec(),
            g: self.prover_param.g,
            h: self.prover_param.h,
            num_vars: supported_num_vars,
        }
    }

    /// Extract the verifier parameters from the public parameters.
    fn extract_verifier_param(&self, supported_num_vars: usize) -> Self::VerifierParam {
        let to_reduce = self.prover_param.num_vars - supported_num_vars;
        Self::VerifierParam {
            num_vars: supported_num_vars,
            g: self.prover_param.g,
            h: self.prover_param.h,
            h_mask: self.h_mask[to_reduce..].to_vec(),
        }
    }

    /// Trim the universal parameters to specialize the public parameters
    /// for multilinear polynomials to the given `supported_num_vars`, and
    /// returns committer key and verifier key. `supported_num_vars` should
    /// be in range `1..=params.num_vars`
    fn trim(
        &self,
        supported_num_vars: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        if supported_num_vars > self.prover_param.num_vars {
            return Err(PCSError::InvalidParameters(format!(
                "SRS does not support target number of vars {}",
                supported_num_vars
            )));
        }

        let ck = Self::ProverParam {
            powers_of_g: self.prover_param.powers_of_g[self.prover_param.num_vars-supported_num_vars..] .to_vec(),
            g: self.prover_param.g,
            h: self.prover_param.h,
            num_vars: supported_num_vars,
        };
        let vk = Self::VerifierParam {
            num_vars: supported_num_vars,
            g: self.prover_param.g,
            h: self.prover_param.h,
            h_mask: self.h_mask[self.prover_param.num_vars-supported_num_vars..].to_vec(),
        };
        Ok((ck, vk))
    }

    /// Build SRS for testing.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore + CryptoRng>(
        rng: &mut R,
        num_vars: usize,
    ) -> Result<Self, PCSError> {
        if num_vars == 0 {
            return Err(PCSError::InvalidParameters(
                "constant polynomial not supported".to_string(),
            ));
        }

        let total_timer = start_timer!(|| "SRS generation");

        let pp_generation_timer = start_timer!(|| "Prover Param generation");

        let g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);
        let mut powers_of_g = vec![Evaluations { evals: vec![] }; num_vars+1];

        let t: Vec<_> = (0..num_vars).map(|_| E::Fr::rand(rng)).collect();
        let scalar_bits = E::Fr::size_in_bits();
        let t_evals = build_eq_x_r_vec(&t)?;
        let scalars = t_evals.iter().map(|t_eval| t_eval.into_repr());
        powers_of_g[0] = Evaluations { evals: scalars.map(|s| g.mul(s).into_affine()).collect() };
        
        for i in 1..num_vars+1 {
            let lastvec = &powers_of_g[i-1].evals;
            powers_of_g[i] = Evaluations {
                evals: lastvec.iter().step_by(2).zip(lastvec.iter().skip(1).step_by(2)).map(|(a, b)| a.add(*b)).collect()
            };
        }
        let pp = Self::ProverParam {
            num_vars,
            powers_of_g,
            g: g.into_affine(),
            h: h.into_affine(),
        };

        end_timer!(pp_generation_timer);

        let vp_generation_timer = start_timer!(|| "VP generation");
        let h_mask = {
            let window_size = FixedBaseMSM::get_mul_window_size(num_vars);
            let h_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, h);
            E::G2Projective::batch_normalization_into_affine(&FixedBaseMSM::multi_scalar_mul(
                scalar_bits,
                window_size,
                &h_table,
                &t,
            ))
        };
        end_timer!(vp_generation_timer);
        end_timer!(total_timer);
        Ok(Self {
            prover_param: pp,
            h_mask,
        })
    }
}


#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_std::test_rng;

    use super::*;

    type E = Bls12_381;

    #[test]
    fn test_srs_gen() -> Result<(), PCSError> {
        let mut rng = test_rng();
        for nv in 4..11 {
            let params = MultilinearUniversalParams::<E>::gen_srs_for_testing(&mut rng, nv)?;
            assert_eq!(params.prover_param.g, params.prover_param.powers_of_g[nv].evals[0]);
        }

        Ok(())
    }
}
