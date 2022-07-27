//! Main module for the HyperPlonk PolyIOP.

use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use pcs::PolynomialCommitmentScheme;
use poly_iop::prelude::{PermutationCheck, ZeroCheck};
use std::rc::Rc;

/// The sub-claim for the HyperPlonk PolyIOP, consists of the following:
///   - the SubClaim for the zero-check PIOP
///   - the SubClaim for the permutation-check PIOP
///   - the SubClaim for public input consistency
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkSubClaim<F: PrimeField, ZC: ZeroCheck<F>, PC: PermutationCheck<F>> {
    /// the SubClaim for the custom gate zerocheck
    pub zero_check_sub_claim: ZC::ZeroCheckSubClaim,
    /// the SubClaim for the permutation check
    pub perm_check_sub_claim: PC::PermutationCheckSubClaim,
    /// the public input consistency check
    pub pub_input_sub_claim: (Vec<F>, F), // (point, expected_eval)
}

/// The proof for the HyperPlonk PolyIOP, consists of the following:
///   - a batch commitment to all the witness MLEs
///   - a batch opening to all the MLEs at certain index
///   - the zero-check proof for checking custom gate-satisfiability
///   - the permutation-check proof for checking the copy constraints
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkProof<
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<E>,
    ZC: ZeroCheck<E::Fr>,
    PC: PermutationCheck<E::Fr>,
> {
    // =======================================================================
    // PCS components
    // =======================================================================
    /// PCS commit for witnesses
    // TODO: replace me with a batch commitment
    pub witness_commits: Vec<PCS::Commitment>,
    /// PCS openings for witness on permutation check point
    // TODO: replace me with a batch opening
    pub witness_perm_check_openings: Vec<PCS::Proof>,
    /// PCS openings for witness on zero check point
    // TODO: replace me with a batch opening
    pub witness_zero_check_openings: Vec<PCS::Proof>,
    /// Evaluates of witnesses on permutation check point
    pub witness_perm_check_evals: Vec<E::Fr>,
    /// Evaluates of witnesses on zero check point
    pub witness_zero_check_evals: Vec<E::Fr>,
    /// PCS openings for selectors on permutation check point
    // TODO: replace me with a batch opening
    pub selector_perm_check_openings: Vec<PCS::Proof>,
    /// PCS openings for selectors on zero check point
    // TODO: replace me with a batch opening
    pub selector_zero_check_openings: Vec<PCS::Proof>,
    /// Evaluates of selectors on permutation check point
    pub selector_perm_check_evals: Vec<E::Fr>,
    /// Evaluates of selectors on zero check point
    pub selector_zero_check_evals: Vec<E::Fr>,
    /// Evaluates of public inputs on r_pi from transcript
    pub pi_eval: E::Fr,
    /// Opening of public inputs on r_pi from transcript
    pub pi_opening: PCS::Proof,
    // =======================================================================
    // IOP components
    // =======================================================================
    /// the custom gate zerocheck proof
    pub zero_check_proof: ZC::ZeroCheckProof,
    /// the permutation check proof for copy constraints
    pub perm_check_proof: PC::PermutationProof,
}

/// The HyperPlonk instance parameters, consists of the following:
///   - the number of variables in the poly-IOP
///   - binary log of the number of public input variables
///   - binary log of the number of selectors
///   - binary log of the number of witness wires
///   - the customized gate function
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkParams {
    /// the number of variables in polys
    pub nv: usize,
    /// binary log of the public input length
    pub log_pub_input_len: usize,
    // binary log of the number of selectors
    pub log_n_selectors: usize,
    /// binary log of the number of witness wires
    pub log_n_wires: usize,
    /// customized gate function
    pub gate_func: CustomizedGates,
}

/// The HyperPlonk proving key, consists of the following:
///   - the hyperplonk instance parameters
///   - the preprocessed polynomials output by the indexer
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HyperPlonkProvingKey<E: PairingEngine, PCS: PolynomialCommitmentScheme<E>> {
    /// hyperplonk instance parameters
    pub params: HyperPlonkParams,
    /// the preprocessed permutation polynomials
    pub permutation_oracles: Rc<DenseMultilinearExtension<E::Fr>>,
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
}

/// Customized gate is a list of tuples of
///     (coefficient, selector_index, wire_indices)
///
/// Example:
///     q_L(X) * W_1(X)^5 - W_2(X)
/// is represented as
/// vec![
///     ( 1,    Some(id_qL),    vec![id_W1, id_W1, id_W1, id_W1, id_W1]),
///     (-1,    None,           vec![id_W2])
/// ]
///
/// NOTE: here coeff is a signed integer, instead of a field element
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CustomizedGates {
    pub(crate) gates: Vec<(i64, Option<usize>, Vec<usize>)>,
}
