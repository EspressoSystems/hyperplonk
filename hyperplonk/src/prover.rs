use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::{univariate::DenseOrSparsePolynomial, DenseMultilinearExtension};
use pcs::MultilinearCommitmentScheme;
use poly_iop::prelude::VirtualPolynomial;

use crate::{errors::HyperPlonkErrors, structs::HyperPlonkProvingKey, witness::WitnessRow};

/// Round 1:
/// 1. Compute and commit wire witness polynomials.
/// 2. Compute public input polynomial.
/// Return the wire witness polynomials and their commitments,
/// also return the public input polynomial.
#[allow(clippy::type_complexity)]
pub(crate) fn run_1st_round<E: PairingEngine, PCS: MultilinearCommitmentScheme<E>>(
    pk: &HyperPlonkProvingKey<E, PCS>,
    pub_input: &[E::Fr],
    num_vars: usize,
    witness: &[WitnessRow<E::Fr>],
) -> Result<
    (
        Vec<DenseMultilinearExtension<E::Fr>>,
        PCS::Commitment,
        DenseMultilinearExtension<E::Fr>,
    ),
    HyperPlonkErrors,
> {
    let witness_polys = WitnessRow::build_mles(witness)?;
    let witnesses_commitment = PCS::multi_commit(&pk.pcs_param, &witness_polys)?;
    let pi = DenseMultilinearExtension::<E::Fr>::from_evaluations_slice(num_vars, pub_input);
    Ok((witness_polys, witnesses_commitment, pi))
}

/// Build a virtual polynomial `f := f(w_0(x),...w_d(x))`  where `f` is
/// the constraint polynomial, i.e., `f(a, b, c) = q_l a(x) + q_r b(x) +
/// q_m a(x)b(x) - q_o c(x)` in vanilla plonk
pub(crate) fn build_f<F: PrimeField>(
    gate_func: Vec<Vec<usize>>,
    w_list: &[DenseOrSparsePolynomial<F>],
) -> Result<VirtualPolynomial<F>, HyperPlonkErrors> {
    todo!()
}
