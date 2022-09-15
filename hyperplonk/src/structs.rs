//! Main module for the HyperPlonk PolyIOP.

use crate::{custom_gate::CustomizedGates, selectors::SelectorColumn};
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::log2;
use pcs::PolynomialCommitmentScheme;
use poly_iop::prelude::{PermutationCheck, ZeroCheck};
use std::rc::Rc;

/// The proof for the HyperPlonk PolyIOP, consists of the following:
///   - a batch commitment to all the witness MLEs
///   - a batch opening to all the MLEs at certain index
///   - the zero-check proof for checking custom gate-satisfiability
///   - the permutation-check proof for checking the copy constraints
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkProof<E, PC, PCS>
where
    E: PairingEngine,
    PC: PermutationCheck<E, PCS>,
    PCS: PolynomialCommitmentScheme<E>,
{
    // =======================================================================
    // PCS components: common
    // =======================================================================
    /// PCS commit for witnesses
    pub w_merged_com: PCS::Commitment,
    pub w_merged_batch_opening: PCS::BatchProof,
    pub w_merged_batch_evals: Vec<E::Fr>,
    // =======================================================================
    // PCS components: permutation check
    // =======================================================================
    /// prod(x)'s evaluations
    /// sequence: prod(0,x), prod(1, x), prod(x, 0), prod(x, 1), prod(1, ..., 1,
    /// 0)
    pub prod_batch_evals: Vec<E::Fr>,
    /// prod(x)'s openings
    /// sequence: prod(0,x), prod(1, x), prod(x, 0), prod(x, 1), prod(1, ..., 1,
    /// 0)
    pub prod_batch_openings: PCS::BatchProof,
    /// PCS openings for witness on permutation check point
    // // TODO: replace me with a batch opening
    // pub witness_perm_check_opening: PCS::Proof,
    // /// Evaluates of witnesses on permutation check point
    // pub witness_perm_check_eval: E::Fr,
    /// PCS openings for selectors on permutation check point
    // TODO: replace me with a batch opening
    pub perm_oracle_opening: PCS::Proof,
    /// Evaluates of selectors on permutation check point
    pub perm_oracle_eval: E::Fr,
    // =======================================================================
    // PCS components: zero check
    // =======================================================================
    // /// PCS openings for witness on zero check point
    // // TODO: replace me with a batch opening
    // pub witness_zero_check_openings: Vec<PCS::Proof>,
    // /// Evaluates of witnesses on zero check point
    // pub witness_zero_check_evals: Vec<E::Fr>,
    /// PCS openings for selectors on zero check point
    // TODO: replace me with a batch opening
    pub selector_oracle_openings: Vec<PCS::Proof>,
    /// Evaluates of selectors on zero check point
    pub selector_oracle_evals: Vec<E::Fr>,
    // =======================================================================
    // PCS components: public inputs
    // =======================================================================
    /// Evaluates of public inputs on r_pi from transcript
    pub pi_eval: E::Fr,
    /// Opening of public inputs on r_pi from transcript
    pub pi_opening: PCS::Proof,
    // =======================================================================
    // IOP components
    // =======================================================================
    /// the custom gate zerocheck proof
    pub zero_check_proof: <PC as ZeroCheck<E::Fr>>::ZeroCheckProof,
    /// the permutation check proof for copy constraints
    pub perm_check_proof: PC::PermutationProof,
}

/// The HyperPlonk instance parameters, consists of the following:
///   - the number of constraints
///   - number of public input columns
///   - number of selector columns
///   - number of witness columns
///   - the customized gate function
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HyperPlonkParams {
    /// the number of constraints
    pub num_constraints: usize,
    /// number of public input
    // public input is only 1 column and is implicitly the first witness column.
    // this size must not exceed number of constraints.
    pub num_pub_input: usize,
    /// customized gate function
    pub gate_func: CustomizedGates,
}

impl HyperPlonkParams {
    /// Number of variables in a multilinear system
    pub fn num_variables(&self) -> usize {
        log2(self.num_constraints) as usize
    }

    /// number of selector columns
    pub fn num_selector_columns(&self) -> usize {
        self.gate_func.num_selector_columns()
    }

    /// number of witness columns
    pub fn num_witness_columns(&self) -> usize {
        self.gate_func.num_witness_columns()
    }
}

/// The HyperPlonk index, consists of the following:
///   - HyperPlonk parameters
///   - the wire permutation
///   - the selector vectors
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HyperPlonkIndex<F: PrimeField> {
    pub params: HyperPlonkParams,
    pub permutation: Vec<F>,
    pub selectors: Vec<SelectorColumn<F>>,
}

impl<F: PrimeField> HyperPlonkIndex<F> {
    /// Number of variables in a multilinear system
    pub fn num_variables(&self) -> usize {
        self.params.num_variables()
    }

    /// number of selector columns
    pub fn num_selector_columns(&self) -> usize {
        self.params.num_selector_columns()
    }

    /// number of witness columns
    pub fn num_witness_columns(&self) -> usize {
        self.params.num_witness_columns()
    }
}

/// The HyperPlonk proving key, consists of the following:
///   - the hyperplonk instance parameters
///   - the preprocessed polynomials output by the indexer
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkProvingKey<E: PairingEngine, PCS: PolynomialCommitmentScheme<E>> {
    /// hyperplonk instance parameters
    pub params: HyperPlonkParams,
    /// the preprocessed permutation polynomials
    pub permutation_oracle: Rc<DenseMultilinearExtension<E::Fr>>,
    /// the preprocessed selector polynomials
    // TODO: merge the list into a single MLE
    pub selector_oracles: Vec<Rc<DenseMultilinearExtension<E::Fr>>>,
    /// the parameters for PCS commitment
    pub pcs_param: PCS::ProverParam,
}

/// The HyperPlonk verifying key, consists of the following:
///   - the hyperplonk instance parameters
///   - the preprocessed polynomials output by the indexer
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkVerifyingKey<E: PairingEngine, PCS: PolynomialCommitmentScheme<E>> {
    /// hyperplonk instance parameters
    pub params: HyperPlonkParams,
    /// the parameters for PCS commitment
    pub pcs_param: PCS::VerifierParam,
    /// Selector's commitment
    // TODO: replace me with a batch commitment
    pub selector_com: Vec<PCS::Commitment>,
    /// Permutation oracle's commitment
    pub perm_com: PCS::Commitment,
}
