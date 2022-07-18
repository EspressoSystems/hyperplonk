use crate::{prelude::PCSErrors, StructuredReferenceString};
use ark_ec::{msm::FixedBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField};
use ark_poly::DenseMultilinearExtension;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{end_timer, rand::RngCore, start_timer, vec::Vec, UniformRand};
use std::collections::LinkedList;

/// Evaluations over {0,1}^n for G1 or G2
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct Evaluations<C: AffineCurve> {
    /// The evaluations.
    pub evals: Vec<C>,
}

/// Universal Parameter
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct UniversalParams<E: PairingEngine> {
    /// prover parameters
    pub prover_param: ProverParam<E>,
    /// h^randomness: h^t1, h^t2, ..., **h^{t_nv}**
    pub h_mask: Vec<E::G2Affine>,
}

/// Prover Parameters
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct ProverParam<E: PairingEngine> {
    /// number of variables
    pub num_vars: usize,
    /// `pp_{num_vars}`, `pp_{num_vars - 1}`, `pp_{num_vars - 2}`, ..., defined
    /// by XZZPD19
    pub powers_of_g: Vec<Evaluations<E::G1Affine>>,
    /// generator for G1
    pub g: E::G1Affine,
    /// generator for G2
    pub h: E::G2Affine,
}

/// Verifier Parameters
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct VerifierParam<E: PairingEngine> {
    /// number of variables
    pub num_vars: usize,
    /// generator of G1
    pub g: E::G1Affine,
    /// generator of G2
    pub h: E::G2Affine,
    /// h^randomness: h^t1, h^t2, ..., **h^{t_nv}**
    pub h_mask: Vec<E::G2Affine>,
}

impl<E: PairingEngine> StructuredReferenceString<E> for UniversalParams<E> {
    type ProverParam = ProverParam<E>;
    type VerifierParam = VerifierParam<E>;

    /// Extract the prover parameters from the public parameters.
    fn extract_prover_param(&self) -> ProverParam<E> {
        self.prover_param.clone()
    }

    /// Extract the verifier parameters from the public parameters.
    fn extract_verifier_param(&self) -> VerifierParam<E> {
        VerifierParam {
            num_vars: self.prover_param.num_vars,
            g: self.prover_param.g,
            h: self.prover_param.h,
            h_mask: self.h_mask.clone(),
        }
    }

    /// Trim the universal parameters to specialize the public parameters
    /// for multilinear polynomials to the given `supported_num_vars`, and
    /// returns committer key and verifier key. `supported_num_vars` should
    /// be in range `1..=params.num_vars`
    fn trim(
        &self,
        supported_num_vars: usize,
    ) -> Result<(ProverParam<E>, VerifierParam<E>), PCSErrors> {
        if supported_num_vars > self.prover_param.num_vars {
            return Err(PCSErrors::InvalidParameters(format!(
                "SRS does not support target number of vars {}",
                supported_num_vars
            )));
        }

        let to_reduce = self.prover_param.num_vars - supported_num_vars;
        let ck = ProverParam {
            powers_of_g: self.prover_param.powers_of_g[to_reduce..].to_vec(),
            g: self.prover_param.g,
            h: self.prover_param.h,
            num_vars: supported_num_vars,
        };
        let vk = VerifierParam {
            num_vars: supported_num_vars,
            g: self.prover_param.g,
            h: self.prover_param.h,
            h_mask: self.h_mask[to_reduce..].to_vec(),
        };
        Ok((ck, vk))
    }

    /// Build SRS for testing.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing<R: RngCore>(rng: &mut R, num_vars: usize) -> Result<Self, PCSErrors> {
        if num_vars == 0 {
            return Err(PCSErrors::InvalidParameters(
                "constant polynomial not supported".to_string(),
            ));
        }

        let total_timer = start_timer!(|| "SRS generation");

        let pp_generation_timer = start_timer!(|| "Prover Param generation");

        let g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        let mut powers_of_g = Vec::new();

        let t: Vec<_> = (0..num_vars).map(|_| E::Fr::rand(rng)).collect();
        let scalar_bits = E::Fr::size_in_bits();

        let mut eq: LinkedList<DenseMultilinearExtension<E::Fr>> =
            LinkedList::from_iter(eq_extension(&t).into_iter());
        let mut eq_arr = LinkedList::new();
        let mut base = eq.pop_back().unwrap().evaluations;

        for i in (0..num_vars).rev() {
            eq_arr.push_front(remove_dummy_variable(&base, i)?);
            if i != 0 {
                let mul = eq.pop_back().unwrap().evaluations;
                base = base
                    .into_iter()
                    .zip(mul.into_iter())
                    .map(|(a, b)| a * b)
                    .collect();
            }
        }

        let mut pp_powers = Vec::new();
        let mut total_scalars = 0;
        for i in 0..num_vars {
            let eq = eq_arr.pop_front().unwrap();
            let pp_k_powers = (0..(1 << (num_vars - i))).map(|x| eq[x]);
            pp_powers.extend(pp_k_powers);
            total_scalars += 1 << (num_vars - i);
        }
        let window_size = FixedBaseMSM::get_mul_window_size(total_scalars);
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);

        let pp_g = E::G1Projective::batch_normalization_into_affine(
            &FixedBaseMSM::multi_scalar_mul(scalar_bits, window_size, &g_table, &pp_powers),
        );

        let mut start = 0;
        for i in 0..num_vars {
            let size = 1 << (num_vars - i);
            let pp_k_g = Evaluations {
                evals: pp_g[start..(start + size)].to_vec(),
            };
            powers_of_g.push(pp_k_g);
            start += size;
        }

        let pp = ProverParam {
            num_vars,
            g: g.into_affine(),
            h: h.into_affine(),
            powers_of_g,
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

/// fix first `pad` variables of `poly` represented in evaluation form to zero
fn remove_dummy_variable<F: Field>(poly: &[F], pad: usize) -> Result<Vec<F>, PCSErrors> {
    if pad == 0 {
        return Ok(poly.to_vec());
    }
    if !poly.len().is_power_of_two() {
        return Err(PCSErrors::InvalidParameters(
            "Size of polynomial should be power of two.".to_string(),
        ));
    }
    let nv = ark_std::log2(poly.len()) as usize - pad;

    Ok((0..(1 << nv)).map(|x| poly[x << pad]).collect())
}

/// Generate eq(t,x), a product of multilinear polynomials with fixed t.
/// eq(a,b) is takes extensions of a,b in {0,1}^num_vars such that if a and b in
/// {0,1}^num_vars are equal then this polynomial evaluates to 1.
fn eq_extension<F: PrimeField>(t: &[F]) -> Vec<DenseMultilinearExtension<F>> {
    let start = start_timer!(|| "eq extension");

    let dim = t.len();
    let mut result = Vec::new();
    for (i, &ti) in t.iter().enumerate().take(dim) {
        let mut poly = Vec::with_capacity(1 << dim);
        for x in 0..(1 << dim) {
            let xi = if x >> i & 1 == 1 { F::one() } else { F::zero() };
            let ti_xi = ti * xi;
            poly.push(ti_xi + ti_xi - xi - ti + F::one());
        }
        result.push(DenseMultilinearExtension::from_evaluations_vec(dim, poly));
    }

    end_timer!(start);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_std::test_rng;
    type E = Bls12_381;

    #[test]
    fn test_srs_gen() -> Result<(), PCSErrors> {
        let mut rng = test_rng();
        for nv in 4..10 {
            let _ = UniversalParams::<E>::gen_srs_for_testing(&mut rng, nv)?;
        }

        Ok(())
    }
}